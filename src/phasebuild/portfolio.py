"""Bounded, opt-in exploration of fixed bundled P+ protocols.

Rfree is displayed after ranking and is never an input to the ranking key.
The recommendation is a heuristic for review, not a correctness certificate.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading

from .data import fingerprint


PROTOCOLS = {
    'pplus_w14': dict(trust_wang_radius=14., trust_stability=False, trust_work_weight=False),
    'pplus_w8': dict(trust_wang_radius=8., trust_stability=False, trust_work_weight=False),
    'stable_w14': dict(trust_wang_radius=14., trust_stability=True, trust_work_weight=False),
    'pplus_w8_weight': dict(trust_wang_radius=8., trust_stability=False, trust_work_weight=True),
}
# Broad review gates for final models, deliberately distinct from the much
# stricter work-weight pilot policy. These are declared heuristics, not a full
# geometry validation or universal protein-structure acceptance standard.
BOND_GATE = .030
ANGLE_GATE = 3.0


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def plan(config, output, names=None, *, inspect=True):
    from .sad_pipeline import DEFAULTS, validate
    allowed = set(DEFAULTS) | {'sad', 'native', 'sequence', 'heavy', 'fpp', 'output', 'runtime', 'reference'}
    if set(config) - allowed:
        raise ValueError('Unknown pipeline configuration keys: ' + ', '.join(sorted(set(config) - allowed)))
    names = list(PROTOCOLS) if names is None else list(names)
    if not names or len(set(names)) != len(names) or set(names) - set(PROTOCOLS):
        raise ValueError('Select distinct predefined portfolio protocols')
    out = Path(output).expanduser().resolve()
    if out.exists():
        raise ValueError('Portfolio output exists; choose a new directory')
    if config.get('method') != 'pplus_bundled' or not config.get('build_model', True):
        raise ValueError('Portfolio requires bundled P+ with model building')
    branches = []
    for name in names:
        cfg = validate({**config, **PROTOCOLS[name], 'output': str(out / name),
                        'trust_cleanup': True, 'trust_partial_cleanup': True}, inspect=inspect)
        branches.append({'name': name, 'config': cfg})
    return {'schema': 1, 'output': str(out), 'branches': branches,
            'ranking_geometry_gates': {'bond_rmsd_angstrom': BOND_GATE, 'angle_rmsd_degrees': ANGLE_GATE},
            'selection_uses_Rfree': False, 'automatic_phase_or_model_merge': False,
            'protocol': 'Fixed common input, flags, seed and build budget; no score-directed early stopping. '
                        'All branches enable sequence cleanup; only declared Wang radius, stability and post-cleanup weight differ.',
            'limitation': 'Experimental fixed portfolio developed using known cases. Heuristic selection does not establish '
                          'an independent end-to-end validation score or guarantee a target Rfree.'}


def sequence_audit(model, sequence, copies, anchor=12):
    """Conservative exact-sequence coverage; repeated sequence positions capped.

    Only unique sequence anchors are counted. Uncertain or unknown assignments
    reduce support instead of inflating completeness. No coordinate reference.
    """
    import gemmi
    from .oasis_options import read_sequences
    multiplicity = Counter(seq for _, seq in read_sequences(sequence))
    sequences = sorted(multiplicity)
    kmers = defaultdict(list)
    for si, seq in enumerate(sequences):
        for pos in range(len(seq) - anchor + 1):
            kmers[seq[pos:pos + anchor]].append((si, pos))
    structure = gemmi.read_structure(str(model))
    if len(structure) != 1:
        raise ValueError('Portfolio audit requires a single coordinate model')
    assigned = Counter()
    protein = 0
    ambiguous = 0
    for chain in structure[0]:
        letters = []
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info.is_amino_acid() and residue.find_atom('CA', '*'):
                protein += 1
                letter = info.one_letter_code.upper()
                letters.append(letter if letter in 'ACDEFGHIKLMNPQRSTVWY' else '?')
            else:
                letters.append('?')
        text = ''.join(letters)
        mappings = defaultdict(set)
        for start in range(len(text) - anchor + 1):
            targets = kmers.get(text[start:start + anchor], [])
            if len(targets) != 1:
                continue
            si, pos = targets[0]
            for offset in range(anchor):
                mappings[start + offset].add((si, pos + offset))
        for targets in mappings.values():
            if len(targets) == 1:
                assigned[next(iter(targets))] += 1
            else:
                ambiguous += 1
    expected = copies * sum(len(seq) * count for seq, count in multiplicity.items())
    covered = sum(min(count, copies * multiplicity[sequences[si]]) for (si, _), count in assigned.items())
    duplicates = sum(max(0, count - copies * multiplicity[sequences[si]]) for (si, _), count in assigned.items())
    return {'expected_residues': expected, 'protein_CA': protein, 'anchored_residues': sum(assigned.values()),
            'unique_covered_positions': covered, 'coverage': covered / expected,
            'duplicate_positions': duplicates, 'duplicate_fraction': duplicates / expected,
            'unanchored_residues': protein - sum(assigned.values()),
            'unanchored_fraction': (protein - sum(assigned.values())) / max(protein, 1),
            'ambiguous_residues': ambiguous, 'anchor_length': anchor,
            'scope': 'Exact unique sequence support only; does not establish coordinate correctness.'}


def rank_candidates(candidates):
    """Geometry gate, then Pareto support set. Rfree is not read here."""
    rows = []
    for candidate in candidates:
        row = dict(candidate)
        metric = row.get('primary_metric') or {}
        audit = row.get('sequence_audit') or {}
        values = [metric.get('r_work'), metric.get('bond_rmsd_angstrom'), metric.get('angle_rmsd_degrees'),
                  audit.get('coverage'), audit.get('duplicate_fraction'), audit.get('unanchored_fraction')]
        valid = all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in values)
        eligible = (row.get('status') == 'completed' and valid and 0 < values[0] < 1
                    and 0 <= values[1] <= BOND_GATE and 0 <= values[2] <= ANGLE_GATE
                    and 0 < values[3] <= 1 and all(0 <= v <= 1 for v in values[4:]))
        row['eligible'] = eligible
        row['eligibility_reason'] = ('Geometry and sequence gates passed' if eligible else
                                     'Failed, missing audit, or outside declared geometry/sequence gates')
        # Explicit tie-break for display only: coverage, duplication, unanchored fraction, work R.
        row['_key'] = (-values[3], values[4], values[5], values[0]) if eligible else (math.inf,) * 4
        rows.append(row)
    for row in rows:
        row['pareto'] = row['eligible'] and not any(
            other['eligible'] and all(a <= b for a, b in zip(other['_key'], row['_key']))
            and any(a < b for a, b in zip(other['_key'], row['_key'])) for other in rows)
    rows.sort(key=lambda r: (not r['eligible'], not r['pareto'], r['_key'], r['name']))
    for row in rows:
        del row['_key']
    frontier = [r['name'] for r in rows if r['pareto']]
    return {'candidates': rows, 'review_candidates': [r['name'] for r in rows if r['eligible']][:3],
            'recommended': frontier[0] if len(frontier) == 1 else None,
            'pareto_candidates': frontier, 'selection_uses_Rfree': False,
            'reason': ('One candidate dominates on the declared work/sequence criteria; provisional recommendation.'
                       if len(frontier) == 1 else 'No unique dominant candidate; retain alternatives for review.'),
            'policy': 'Broad heuristic review gates: bond RMSD<=0.030A, angle RMSD<=3.0deg and positive anchored coverage; '
                      'Pareto maximize exact sequence coverage and minimize duplicate fraction, unanchored fraction, work R. '
                      'Display ties by those same criteria, then protocol name. No Rfree ranking; no merges.',
            'limitation': 'Geometry RMSDs are restraint agreement, not full stereochemical validation. '
                          'Work R and sequence support can favor an incorrect model; recommendation remains heuristic.'}


def collect(name, directory, config):
    directory = Path(directory)
    state = json.loads((directory / 'pipeline_report.json').read_text())
    if state['status'] != 'completed':
        raise ValueError('Branch did not complete')
    if state['config'] != config:
        raise ValueError('Branch configuration differs from declared protocol')
    for stage in state['completed']:
        for artifact in stage['artifacts']:
            if fingerprint(artifact['path']) != artifact:
                raise ValueError('Branch stage artifact changed: ' + artifact['path'])
    model = state['primary_model']
    if fingerprint(model['path']) != model:
        raise ValueError('Branch model hash changed')
    for key in ['sad', 'native', 'sequence', 'heavy']:
        if state['config'][key] != config[key] or state['inputs'][key] != fingerprint(config[key] or config['sad']):
            raise ValueError('Branch input identity differs: ' + key)
    if state['config']['free_value'] != config['free_value']:
        raise ValueError('Branch free partition differs')
    return {'name': name, 'status': 'completed', 'directory': str(directory), 'model': model,
            'primary_metric': state['primary_metric'],
            'sequence_audit': sequence_audit(model['path'], config['sequence'], config['copies'])}


def launch_branch(config_file, directory, stop):
    root = Path(__file__).resolve().parents[2]
    command = [sys.executable, '-B', '-u', str(root / 'run_sad_portfolio.py'), '--child-config', str(config_file)]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    with (Path(config_file).with_suffix('.log')).open('x') as log:
        process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            while process.poll() is None:
                if stop.wait(.25):
                    process.terminate()  # Pipeline SIGTERM handler also cancels its stage process group.
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    raise InterruptedError('Portfolio cancelled; branch artifacts retained')
            if process.returncode:
                raise RuntimeError('Branch exited with code ' + str(process.returncode))
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=15)


def run(config, output, *, workers=2, names=None, runner=None):
    if isinstance(workers, bool) or workers not in (1, 2):
        raise ValueError('Portfolio workers must be 1 or 2')
    declaration = plan(config, output, names)
    out = Path(declaration['output'])
    out.mkdir(parents=True, exist_ok=False)
    write_json(out / 'PLAN.json', declaration)
    from .sad_pipeline import source_identity
    write_json(out / 'SOURCES.json', {'pipeline': source_identity(), 'portfolio_entry':
               fingerprint(Path(__file__).resolve().parents[2] / 'run_sad_portfolio.py')})
    inputs = {key: fingerprint(declaration['branches'][0]['config'][key] or declaration['branches'][0]['config']['sad']) for key in ['sad', 'native', 'sequence', 'heavy']}
    write_json(out / 'INPUTS.json', inputs)
    stop = threading.Event()
    handlers = {}
    if threading.current_thread() is threading.main_thread():
        for sig in (signal.SIGINT, signal.SIGTERM):
            handlers[sig] = signal.signal(sig, lambda *_: stop.set())
    rows = []
    runner = runner or launch_branch

    def execute(branch):
        name, cfg = branch['name'], branch['config']
        if stop.is_set():
            return {'name': name, 'status': 'cancelled'}
        cfg_file = out / (name + '.json')
        write_json(cfg_file, cfg)
        try:
            runner(cfg_file, cfg['output'], stop)
            return collect(name, cfg['output'], cfg)
        except Exception as error:
            return {'name': name, 'status': 'cancelled' if stop.is_set() else 'failed', 'error': str(error),
                    'directory': cfg['output']}

    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(execute, branch) for branch in declaration['branches']]
            for future in as_completed(futures):
                rows.append(future.result())
                write_json(out / 'PROGRESS.json', {'completed': rows, 'planned': len(futures)})
        ranking = rank_candidates(rows)
        unchanged = all(fingerprint(item['path']) == item for item in inputs.values())
        if not unchanged:
            ranking.update(recommended=None, review_candidates=[], input_integrity_error=True)
        ranking.update(status=('cancelled' if stop.is_set() else 'completed' if all(r['status'] == 'completed' for r in rows)
                               else 'partial' if any(r['status'] == 'completed' for r in rows) else 'failed'),
                       inputs_unchanged=unchanged, workers=workers, automatic_phase_or_model_merge=False)
        if not unchanged:
            ranking['status'] = 'failed'
        write_json(out / 'PORTFOLIO_REPORT.json', ranking)
        lines = ['# P+ fixed-protocol exploration', '', 'Status: ' + ranking['status'], '',
                 ranking['reason'], '', 'Provisional recommendation: ' + str(ranking['recommended']), '',
                 '| Protocol | Status | Coverage | Work R | Rfree (display only) | Geometry gate |',
                 '|---|---|---:|---:|---:|---|']
        for row in ranking['candidates']:
            metric = row.get('primary_metric') or {}
            audit = row.get('sequence_audit') or {}
            def percent(value):
                return '—' if value is None else f'{100 * value:.2f}%'
            lines.append('| ' + ' | '.join([row['name'], row['status'], percent(audit.get('coverage')),
                         percent(metric.get('r_work')), percent(metric.get('r_free')), str(row['eligible'])]) + ' |')
        lines += ['', ranking['policy'], '', ranking['limitation'], '', declaration['limitation'], '',
                  'All original branch outputs remain in their individual directories. No phase/model averaging.']
        (out / 'REPORT.md').write_text('\n'.join(lines) + '\n')
        return ranking
    finally:
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
