"""Self-contained maximum-likelihood refinement adapter (minirefmac).

A from-scratch REFMAC5-like engine that needs no CCP4 installation: the monomer
restraint dictionary ships inside the package. This adapter exposes the same
call contract and report schema as :func:`phasebuild.refinement.refmac_refine`,
so the pipeline can switch engines by name while keeping the free-set audit and
the fixed-final-cycle selection policy.

Output (inside ``output``): ``refined.pdb``, ``refined.mtz``, ``model.cif``,
``refinement_report.json``.
"""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import time
import gemmi
import numpy as np
from . import __version__
from .data import fingerprint
from .refinement import reflection_partition, verify_free_flags

# Default standalone deployment roots for minirefmac.
DEFAULT_ROOTS = (
    '/Volumes/Fast_4T/Workbuddy_Oasis/oasis-refine',
)


def locate(root=None):
    """Return a deployment root containing the minirefmac package."""
    if root:
        path = Path(root).expanduser()
        if (path / 'minirefmac' / 'cli.py').is_file():
            return str(path)
        raise ValueError('Configured minirefmac root lacks minirefmac/cli.py: ' + str(path))
    for base in DEFAULT_ROOTS:
        if (Path(base) / 'minirefmac' / 'cli.py').is_file():
            return base
    raise ValueError('minirefmac engine was not found; set --minirefmac-root to its deployment root')


def _statistics(report, partition):
    cycles = []
    for row in report.get('cycles', ()):
        cycles.append({'cycle': int(round(row.get('cycle', 0))),
                       'r_work': float(row.get('r_work', float('nan'))),
                       'r_free': float(row.get('r_free', float('nan'))),
                       'bond_rmsd_angstrom': row.get('rmsd_bond'),
                       'angle_rmsd_degrees': row.get('rmsd_angle')})
    if not cycles:
        raise RuntimeError('minirefmac produced no refinement cycles')
    used_total = int(report.get('n_ref', partition['usable_observations']))
    used_free = int(partition['usable_free'])
    return {'initial': cycles[0], 'final': cycles[-1], 'cycles': cycles,
            'reflections_used_total': used_total,
            'reflections_used_free': used_free,
            'reflections_used_work': used_total - used_free,
            'resolution_low': float(report.get('d_max', partition['observed_resolution_low'])),
            'resolution_high': float(report.get('d_min', partition['observed_resolution_high']))}


def _hl_quartet(mtz_path, destination):
    """Add the REFMAC-style combined phase quartet to a refined MTZ.

    REFMAC writes ``PHCOMB`` plus the four-component orthogonal
    representation ``HLACOMB``/``HLBCOMB``/``HLCCOMB``/``HLDCOMB``. The
    combined phase is the phase that carries the experimental prior together
    with the model contribution, and the quartet encodes it with
    ``A = kappa*cos(phi)``, ``B = kappa*sin(phi)``, ``C = D = 0``.

    CCP4 Buccaneer needs those columns in its correlation mode: without them
    it cannot build the map and aborts with a grid-coordinate failure. The
    in-house engine does not emit them, so the adapter derives them here from
    the refined phases (``PHWT``) and the figure of merit (``FOM``).

    The mapping is a documented approximation, not a byte-identical copy of
    the REFMAC output: the same phase and the same figure of merit are used,
    but REFMAC's own combined-phase weighting is not reproduced.
    """
    m = gemmi.read_mtz_file(str(mtz_path))
    labels = m.column_labels()
    if 'PHCOMB' in labels and 'HLACOMB' in labels:
        return False
    if 'PHWT' not in labels:
        raise ValueError('refined MTZ lacks PHWT; the combined phase quartet '
                         'cannot be derived')
    data = np.array(m, copy=True)
    phase = np.deg2rad(data[:, labels.index('PHWT')])
    fom = data[:, labels.index('FOM')] if 'FOM' in labels else np.ones(len(data))
    # REFMAC caps the quartet amplitude at 81 (the same limit it applies to
    # its own HL output), so the derived columns use the same range.
    kappa = np.clip(np.abs(fom) * 81.0, 0.0, 81.0)
    hla = kappa * np.cos(phase)
    hlb = kappa * np.sin(phase)
    hlc = np.zeros(len(data))
    hld = np.zeros(len(data))
    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.set_cell_for_all(m.cell)
    for dataset in m.datasets:
        if dataset.id:
            out.add_dataset(dataset.project_name or 'data')
            break
    else:
        out.add_dataset('data')
    for column in m.columns:
        if column.label in ('H', 'K', 'L'):
            continue
        out.add_column(column.label, column.type)
    out.add_column('PHCOMB', 'P')
    for name in ('HLACOMB', 'HLBCOMB', 'HLCCOMB', 'HLDCOMB'):
        out.add_column(name, 'A')
    comb = np.degrees(np.arctan2(hlb, hla)) % 360.0
    comb = np.where(np.isfinite(comb), comb, 0.0)
    out.set_data(np.column_stack([data[:, :3], data[:, 3:], comb,
                                  hla, hlb, hlc, hld]).astype(np.float32))
    out.write_to_file(str(destination))
    return True


def refine(mtz_path, model_path, output, cycles=10, root=None, fp='F', sigfp='SIGF',
           free='FreeR_flag', free_value=0, timeout=3600, ncyc_maxiter=10,
           phase_prior=False, phase_prior_mtz=None):
    """Refine ``model_path`` against ``mtz_path`` with minirefmac.

    The free set is preserved and audited exactly as the REFMAC adapter does.
    Selection is the fixed final cycle; Rfree is never read for selection.

    ``phase_prior_mtz`` supplies a PHIB/FOM (or HLA..D) phase source that is
    merged into the observed data and used as an MLHL phase prior. This is an
    approximation of the REFMAC route's fixed phase restraint: the in-house
    engine refines the phases as well, so the report records the difference.
    """
    if not isinstance(cycles, int) or cycles < 0:
        raise ValueError('Refinement cycles must be a nonnegative integer')
    root = locate(root)
    mtz, partition = reflection_partition(mtz_path, fp, sigfp, free, free_value)
    structure = gemmi.read_structure(str(model_path))
    if len(structure) != 1 or not structure[0]:
        chains = len(structure[0]) if len(structure) else 0
        residues = sum(1 for chain in structure[0] for _ in chain) if len(structure) else 0
        raise ValueError('Exactly one populated model required; the input model has '
                         '%d model(s), %d chain(s), %d residue(s)' % (len(structure), chains, residues))
    sg = structure.find_spacegroup()
    if sg is None or sg.number != mtz.spacegroup.number or not np.allclose(
            structure.cell.parameters, mtz.cell.parameters, rtol=0., atol=.01):
        raise ValueError('Model/data crystal cell or space group differs')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    structure.write_pdb(str(out / 'input.pdb'))
    work = out / 'work'
    work.mkdir()
    tag = out.name
    # Prefer the interpreter running this pipeline: it already carries the
    # pinned gemmi/numpy/scipy stack that minirefmac needs.
    interpreter = os.environ.get('PHASEBUILD_PYTHON') or sys.executable or 'python3'
    engine_mtz = Path(mtz_path).resolve()
    phase_audit = None
    if phase_prior_mtz:
        from .phase_refinement import merge_experimental_phases
        merged = out / 'phase_prior_merged.mtz'
        phase_audit = merge_experimental_phases(str(engine_mtz), str(phase_prior_mtz),
                                                 merged, 'phifom', free_value=free_value,
                                                 phase_prefix='DM_')
        engine_mtz = merged
        phase_prior = True
    command = [interpreter, '-m', 'minirefmac.cli',
               '--mtz', str(engine_mtz), '--pdb', str(out / 'input.pdb'),
               '--out-dir', str(work), '--tag', tag,
               '--fp', fp, '--sigfp', sigfp, '--free', free,
               '--free-value', str(int(free_value)),
               '--ncyc', str(int(cycles)), '--maxiter', str(int(ncyc_maxiter)),
               '--quiet']
    if phase_prior:
        command.append('--phase-prior')
    env = dict(os.environ)
    env['PYTHONPATH'] = root + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    keywords = {'engine': 'minirefmac', 'ncyc': int(cycles), 'maxiter': int(ncyc_maxiter),
                'fp': fp, 'sigfp': sigfp, 'free': free, 'free_value': int(free_value)}
    report = {'version': __version__, 'engine': 'in-house-minirefmac',
              'model_input': fingerprint(model_path), 'reflections': partition,
              'root': root, 'command': command, 'keywords': keywords,
              'requested_cycles': cycles,
              'selection': 'Fixed final cycle; Rfree never used for selection or early stopping.',
              'scope': 'Native amplitudes, original free flags, no target reference or phase restraints. '
                       'In-house ML refinement with bundled monomer restraints. '
                       'These values are not a certified independent blind-test Rfree.'}
    if phase_prior:
        # The prior is always recorded, including the case where it is read
        # from the input MTZ's own PHIB/FOM columns instead of a separate
        # prior file. Silently omitting it here would hide the difference
        # from the REFMAC route's fixed phase restraint.
        report['phase_prior'] = {
            'mode': 'mlhl_prior',
            'source': (fingerprint(phase_prior_mtz) if phase_prior_mtz
                       else {'path': str(Path(mtz_path).resolve()),
                             'columns': 'PHIB/FOM read from the input MTZ'}),
            'audit': phase_audit,
            'approximation': True,
            'equivalent_to_fixed_phase_restraint': False,
            'note': 'The in-house engine refines the phases with the prior; the REFMAC route '
                    'restrains them to the supplied values. Results are therefore not '
                    'numerically equivalent to the external route.'}
    started = time.perf_counter()
    try:
        with (out / 'refmac.log').open('w') as log:
            result = subprocess.run(command, text=True, stdout=log, stderr=subprocess.STDOUT,
                                    cwd=str(out), env=env, timeout=timeout)
        report['returncode'] = result.returncode
        raw_mtz = work / (tag + '_refined.mtz')
        raw_pdb = work / (tag + '_refined.pdb')
        raw_json = work / (tag + '_report.json')
        if result.returncode != 0 or not raw_mtz.is_file() or not raw_pdb.is_file():
            raise RuntimeError('minirefmac failed; see ' + str(out / 'refmac.log'))
        engine_report = json.loads(raw_json.read_text()) if raw_json.is_file() else {}
        shutil.copyfile(raw_pdb, out / 'refined.pdb')
        shutil.copyfile(raw_mtz, out / 'refined.mtz')
        # Buccaneer's correlation mode requires the combined phase quartet.
        # The in-house engine does not write it, so it is derived here.
        quartet = _hl_quartet(out / 'refined.mtz', out / 'refined.mtz')
        stats = _statistics(engine_report, partition)
        if stats['final']['cycle'] != cycles:
            raise RuntimeError('minirefmac did not complete the requested cycle count')
        output_mtz = gemmi.read_mtz_file(str(out / 'refined.mtz'))
        flag_audit = verify_free_flags(mtz, output_mtz, free, fp, sigfp)
        refined = gemmi.read_structure(str(out / 'refined.pdb'))
        refined.make_mmcif_document().write_file(str(out / 'model.cif'))
        report.update(status='success', statistics=stats, free_flags_preserved=True,
                      free_flag_audit=flag_audit,
                      combined_phase_quartet={
                          'derived': bool(quartet),
                          'columns': ['PHCOMB', 'HLACOMB', 'HLBCOMB', 'HLCCOMB', 'HLDCOMB'],
                          'source': 'PHWT and FOM of the refined map',
                          'approximation': True,
                          'note': 'The in-house engine does not emit the REFMAC combined '
                                  'phase quartet. The adapter derives it from the refined '
                                  'phase and figure of merit, so the builder receives the '
                                  'columns it needs; the values are not identical to the '
                                  'REFMAC output.'},
                      output_model=fingerprint(out / 'model.cif'),
                      output_mtz=fingerprint(out / 'refined.mtz'))
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['seconds'] = time.perf_counter() - started
        (out / 'refinement_report.json').write_text(
            json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report
