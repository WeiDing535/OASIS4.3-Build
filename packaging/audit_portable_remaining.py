"""Read-only endpoint checks for the portable remaining-SAD cohort.

``audit_completed(case, state, run)`` reads construction/refinement evidence,
never a target reference, and never writes a file.  A failure is returned as an
audit issue so execution failures, scientific failures and missing evaluations
can retain their own denominators in the caller's report.
"""
import hashlib
import json
import math
from pathlib import Path
import re
import traceback
import xml.etree.ElementTree as ET

import gemmi
import numpy as np

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if (_HERE.parent / 'engine').is_dir() else _HERE.parents[1]
BUNDLE = ROOT / 'dist/local_macos_x86_64'
PRIOR_LABELS = ('HLA', 'HLB', 'HLC', 'HLD')
AUDITOR_REVISION = 'portable-remaining-fixedHL-exactHKL-v1'
SANDBOX_PROFILE = '(version 1)(allow default)(deny file-read* (subpath "/Applications/ccp4-9"))'


def _read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _fp(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': _sha(path)}


def _check_fp(item):
    actual = _fp(item['path'])
    _require(actual['path'] == item['path'] and actual['sha256'] == item['sha256'],
             'File fingerprint differs: ' + str(item['path']))
    return actual


def _contained(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def _flag(command, name):
    _require(command.count(name) == 1, 'Missing/duplicate command flag: ' + name)
    return command[command.index(name) + 1]


def _equal(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def _declared(case, path_key, hash_key):
    path = Path(case[path_key]).resolve()
    digest = case.get(hash_key)
    if digest is None:
        # Legacy preparation manifests store paths relative to their run root.
        matches = [v for k, v in case.get('prepared_sha256', {}).items()
                   if path.as_posix().endswith('/' + k)]
        _require(len(matches) == 1, 'Missing frozen input hash: ' + path_key)
        digest = matches[0]
    item = {'path': str(path), 'sha256': digest}
    _check_fp(item)
    return item


def _parse_refmac(log_path, xml_path):
    """Independently reconcile the four-decimal raw table with REFMAC XML."""
    text = Path(log_path).read_text()
    root = ET.parse(str(xml_path)).getroot()
    _require(int(root.findtext('err_level', '999')) == 0, 'REFMAC XML reports an error')
    parts = text.split('Ncyc    Rfact    Rfree')
    _require(len(parts) > 1, 'REFMAC raw R-factor table missing')
    cycles = []
    for line in parts[-1].split('$TEXT:')[0].splitlines():
        fields = line.split()
        if len(fields) < 11 or not fields[0].isdigit():
            continue
        try:
            values = [float(value) for value in fields[1:11]]
        except ValueError:
            continue
        _require(np.isfinite(values).all() and min(values[:2]) >= 0.,
                 'Invalid finite/nonnegative REFMAC numeric cycle table')
        cycles.append({'cycle': int(fields[0]), 'r_work': values[0], 'r_free': values[1],
                       'bond_rmsd_angstrom': values[5], 'angle_rmsd_degrees': values[7]})
    xml_cycles = root.findall('./Overall_stats/stats_vs_cycle/new_cycle')
    _require([r['cycle'] for r in cycles] == list(range(11)), 'REFMAC cycles must be 0..10')
    _require(len(xml_cycles) == 11, 'REFMAC XML cycle count is not eleven including cycle zero')
    for row, xml in zip(cycles, xml_cycles):
        _require(row['cycle'] == int(xml.findtext('cycle')), 'REFMAC log/XML cycle mismatch')
        for key, name in [('r_work', 'r_factor'), ('r_free', 'r_free')]:
            _require(abs(row[key] - float(xml.findtext(name))) <= .0011,
                     'REFMAC log/XML R-factor mismatch')
    total = int(root.findtext('./Overall_stats/n_reflections_all'))
    free = int(root.findtext('./Overall_stats/n_reflections_free'))
    _require(0 < free < total, 'REFMAC working/test reflection counts invalid')
    return {'initial': cycles[0], 'final': cycles[-1], 'cycles': cycles,
            'reflections_used_total': total, 'reflections_used_free': free,
            'reflections_used_work': total - free,
            'resolution_low': float(root.findtext('./Overall_stats/resolution_low')),
            'resolution_high': float(root.findtext('./Overall_stats/resolution_high')),
            'r_factor_precision_decimals': 4}


def _hkl_rows(mtz):
    data = np.array(mtz, copy=False)
    _require(np.isfinite(data[:, :3]).all() and np.all(data[:, :3] == np.rint(data[:, :3])),
             'Noninteger/missing Miller indices')
    index = {tuple(row[:3].astype(int)): i for i, row in enumerate(data)}
    _require(len(index) == len(data), 'Duplicate Miller indices')
    return data, index


def _free_flags(native, refined):
    before, bi = _hkl_rows(native)
    after, ai = _hkl_rows(refined)
    bl, al = native.column_labels(), refined.column_labels()
    for label in ['FreeR_flag', 'F', 'SIGF']:
        _require(bl.count(label) == al.count(label) == 1, 'Required unique output column: ' + label)
    _require(not (ai.keys() - bi.keys()), 'Refinement added unexpected HKLs')
    dropped = bi.keys() - ai.keys()
    for key in dropped:
        row = before[bi[key]]
        f, sigma = row[bl.index('F')], row[bl.index('SIGF')]
        _require(not (np.isfinite(f) and np.isfinite(sigma) and f >= 0 and sigma > 0),
                 'Refinement removed an observed reflection')
    retained = np.array([bi[k] for k in ai], dtype=int)
    a = before[retained, bl.index('FreeR_flag')]
    b = after[:, al.index('FreeR_flag')]
    _require(np.array_equal(a, b, equal_nan=True), 'Original free-set flags changed by exact HKL')
    original_obs = before[retained][:, [bl.index('F'), bl.index('SIGF')]]
    output_obs = after[:, [al.index('F'), al.index('SIGF')]]
    same = (original_obs == output_obs) | (np.isnan(original_obs) & np.isnan(output_obs))
    return {'all_retained_flags_unchanged': True, 'retained_reflections': len(after),
            'removed_unobserved_reflections': len(dropped),
            'output_amplitude_or_sigma_changed_rows': int((~same.all(axis=1)).sum()),
            'output_scaling_policy': 'Normal REFMAC output scaling allowed; unchanged native input and exact-HKL original free flags mandatory'}


def _check_build(build, command, case, previous_map, previous_model, binary, library, number):
    _require(command[0] == binary['path'], 'Construction binary path differs')
    _require(Path(_flag(command, '-mtzin')).resolve() == Path(previous_map).resolve(),
             'Construction uses another map instead of its own previous output')
    _require(Path(_flag(command, '-pdbin')).resolve() == Path(previous_model).resolve(),
             'Construction uses another model instead of its own previous output')
    _require(_flag(command, '-seqin') == str(Path(case['sequence']).resolve()), 'Construction sequence changed')
    _require(_flag(command, '-mtzin-ref') == library[0]['path'] and
             _flag(command, '-pdbin-ref') == library[1]['path'], 'Training resources changed')
    exact_flags = {'-colin-ref-fo': '/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
                   '-colin-ref-hl': '/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
                   '-colin-fo': '/*/*/[F,SIGF]', '-colin-phifom': '/*/*/[PHWT,FOM]',
                   '-colin-fc': '/*/*/[FWT,PHWT]', '-cycles': '3' if number == 1 else '2', '-jobs': '1',
                   '-pdbout': str(build / 'model.pdb'), '-xmlout': str(build / 'buccaneer.xml')}
    for name, value in exact_flags.items():
        _require(_flag(command, name) == value, 'Construction scientific flag changed: ' + name)
    boolean_flags = ['-find', '-grow', '-join', '-link', '-sequence', '-correct', '-filter',
                     '-ncsbuild', '-prune', '-rebuild', '-model-filter',
                     '-no-correlation-mode' if number == 1 else '-correlation-mode']
    if case.get('build_semet', False):
        boolean_flags.append('-build-semet')
    for name in boolean_flags:
        _require(command.count(name) == 1, 'Missing/duplicate construction stage: ' + name)
    _require(('-build-semet' in command) == bool(case.get('build_semet', False)), 'SeMet flag changed')
    _require(('-correlation-mode' in command) == (number > 1) and
             ('-no-correlation-mode' in command) == (number == 1), 'Likelihood target changed')
    value_flags = ['-mtzin', '-pdbin', '-seqin', '-mtzin-ref', '-pdbin-ref'] + list(exact_flags)
    _require(len(command) == 1 + 2 * len(value_flags) + len(boolean_flags),
             'Unexpected extra construction command argument')
    xml = ET.parse(str(build / 'buccaneer.xml')).getroot()
    _require([int(c.findtext('Number')) for c in xml.findall('./Cycles/Cycle')] ==
             list(range(1, 4 if number == 1 else 3)), 'Buccaneer raw XML cycles changed')
    _require('cbuccaneer: Normal termination' in (build / 'buccaneer.log').read_text(),
             'Buccaneer did not report normal termination')
    for name in ['model.pdb', 'model.cif']:
        _require((build / name).is_file(), 'Completed construction model missing: ' + name)


def _sources(report, run):
    _require(report.get('source_changed_during_run') is False, 'Scientific sources changed while running')
    engine = Path(run).resolve() / 'engine'
    _require(report.get('sources'), 'Missing recorded scientific sources')
    for item in report['sources']:
        _require(_contained(item['path'], engine / 'src'), 'Construction imported a source outside frozen engine')
        _check_fp(item)
    for item in [report['binary']] + report['training_library']:
        _require(_contained(item['path'], BUNDLE / 'runtime'), 'Construction resource outside bundled runtime')
        _check_fp(item)


def _audit_free(case, state, run):
    folder = Path(run).resolve() / 'cases' / case['id']
    primary = Path(state['primary_model']).resolve()
    pipeline = primary.parent
    _require(_contained(pipeline, folder), 'Primary output escapes case directory')
    report = _read(pipeline / 'hybrid_report.json')
    _require(state['status'] == report['status'] == 'completed', 'Incomplete fixed endpoint')
    _require(state.get('completed_rounds', 5) == 5, 'State completed round count differs')
    _require(report['parameters'] == {'rounds': 5, 'refine_cycles': 10,
             'build_semet': bool(case.get('build_semet', False)), 'model_filter': True},
             'Fixed scientific parameters differ')
    _sources(report, run)
    inputs = {'initial_map': _declared(case, 'initial_mtz', 'initial_sha256'),
              'native': _declared(case, 'native', 'native_sha256'),
              'sequence': _declared(case, 'sequence', 'sequence_sha256')}
    phase = _declared(case, 'phase_source', 'phase_sha256')
    seed_map = _declared(case, 'seed_map', 'seed_map_sha256')
    _require(len(report['inputs']) == 4, 'Unexpected/missing construction inputs')
    for key, item in inputs.items():
        _require(report['inputs'][key] == item, 'Construction input substituted: ' + key)
    seed = report['inputs']['phasebuild_seed_model']
    _check_fp(seed)
    _require(_contained(seed['path'], folder), 'Native seed is outside its declared case')
    _require(state['seed_model_fingerprint'] == seed and
             str(Path(state['seed_model']).resolve()) == seed['path'],
             'Construction seed differs from frozen state input')
    if case.get('seed_model'):
        _require(seed == _declared(case, 'seed_model', 'seed_sha256'), 'Frozen native seed differs')
    _require(report.get('resumed_first_construction') is None, 'Undeclared first-construction recovery')
    _require(len(report['stages']) == 5, 'Not five construction/refinement rounds')
    _require(_equal(report['primary_metric'], report['stages'][-1]['statistics']), 'Primary metric is not final round')
    native, source = [gemmi.read_mtz_file(str(p)) for p in [case['native'], case['phase_source']]]
    nd, ni = _hkl_rows(native)
    sd, si = _hkl_rows(source)
    _require(native.spacegroup.hm == source.spacegroup.hm and
             np.allclose(native.cell.parameters, source.cell.parameters, rtol=0., atol=.01),
             'Native and original experimental prior crystal settings differ')
    nl, sl = native.column_labels(), source.column_labels()
    for label, kind in [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')]:
        _require(nl.count(label) == 1 and native.column_with_label(label).type == kind,
                 'Native required column/type differs: ' + label)
    for label in PRIOR_LABELS:
        _require(sl.count(label) == 1 and source.column_with_label(label).type == 'A',
                 'Original phase required column/type differs: ' + label)
    expected = np.full((len(nd), 4), np.nan, dtype=np.float32)
    source_cols = [sl.index(s) for s in PRIOR_LABELS]
    for hkl, index in ni.items():
        if hkl in si:
            values = sd[si[hkl], source_cols]
            if np.isfinite(values).all():
                expected[index] = values
    missing = ~np.isfinite(expected).all(axis=1)
    _require((~missing).any(), 'No finite exact-HKL original phase priors')
    observed = (np.isfinite(nd[:, nl.index('F')]) & np.isfinite(nd[:, nl.index('SIGF')]) &
                (nd[:, nl.index('F')] >= 0.) & (nd[:, nl.index('SIGF')] > 0.))
    flags = nd[:, nl.index('FreeR_flag')]
    _require(np.isfinite(flags[observed]).all() and np.all(flags[observed] == np.rint(flags[observed])),
             'Usable native observations lack valid original flags')
    _require(np.any(observed & (flags == 0)) and np.any(observed & (flags != 0)),
             'Both original working and test sets must exist')
    distances = native.make_d_array()[observed]
    keywords = ('LABIN FP=F SIGFP=SIGF FREE=FreeR_flag HLA=SAD_HLA HLB=SAD_HLB HLC=SAD_HLC HLD=SAD_HLD\n'
                'FREE 0\nREFI TYPE REST RESO %.6f %.6f\n' %
                (float(distances.max()) * 1.0001, float(distances.min()) * .9999) +
                'REFI RESI MLKF\nREFI BREF ISOT\nMAKE HYDR N HOUT N BUIL N\n'
                'SCAL TYPE BULK\nSCAL MLSC WORK\nWEIGHT AUTO\nNCYC 10\nEND\n')
    previous_map, previous_model = inputs['initial_map']['path'], seed['path']
    rounds = []
    for number, stage in enumerate(report['stages'], 1):
        build = pipeline / ('round_%02d_build' % number)
        refine = pipeline / ('round_%02d_refine' % number)
        bc = _read(build / 'command.json')
        rr = _read(refine / 'refinement_report.json')
        _require(stage['round'] == number and stage['buccaneer_cycles'] == (3 if number == 1 else 2)
                 and stage['refmac_cycles'] == 10 and stage['construction'] == 'external-Buccaneer',
                 'Completed round scientific schedule differs')
        _check_build(build, bc, case, previous_map, previous_model,
                     report['binary'], report['training_library'], number)
        _require(rr['status'] == 'success' and rr['returncode'] == 0 and rr['requested_cycles'] == 10,
                 'REFMAC did not complete ten cycles')
        _require(rr['phase_target_verified'] and rr['phase_modified_likelihood_active'] and
                 rr['use_experimental_prior'] and rr['prior_representation'] == 'hl',
                 'Experimental-HL phase likelihood not active')
        _require(re.search(r'Residual\s*:\s*Phase modified Maximum Likelihood for Fs',
                           (refine / 'refmac.log').read_text()) is not None,
                 'REFMAC raw log lacks active phase likelihood')
        _require(rr['original_native_input'] == inputs['native'] and
                 rr['fixed_original_experimental_prior'] == phase,
                 'Native input or fixed original HL changed between rounds')
        _require(rr['keywords'] == (refine / 'keywords.txt').read_text() == keywords,
                 'REFMAC exact scientific keywords changed')
        _require(rr['command'][1:] == ['HKLIN', 'input.mtz', 'XYZIN', 'input.pdb',
                 'HKLOUT', 'refined.mtz', 'XYZOUT', 'refined.pdb', 'XMLOUT', 'refmac.xml'],
                 'REFMAC staged file arguments changed')
        _require(rr['command'][0] == rr['binary']['path'] and
                 _contained(rr['binary']['path'], BUNDLE / 'runtime'), 'REFMAC outside bundled runtime')
        for key in ['binary', 'source', 'wrapper_source', 'model_input', 'input', 'staged_input', 'output_model', 'output_mtz']:
            _check_fp(rr[key])
        for key in ['source', 'wrapper_source']:
            _require(_contained(rr[key]['path'], Path(run) / 'engine/src'), 'REFMAC source outside frozen engine')
        _require(rr['model_input'] == _fp(build / 'model.cif'), 'Refinement uses another construction model')
        _require(rr['staged_input'] == _fp(refine / 'input.mtz') and
                 rr['input']['sha256'] == rr['staged_input']['sha256'], 'Staged MTZ is not byte-identical')
        staged = gemmi.read_mtz_file(str(refine / 'input.mtz'))
        arr = np.array(staged, copy=False)
        _require(staged.column_labels() == nl + ['SAD_' + s for s in PRIOR_LABELS],
                 'Native columns reordered, renamed or unexpected columns added')
        _require(np.array_equal(arr[:, :nd.shape[1]], nd, equal_nan=True),
                 'Original native columns/row order changed')
        _require(np.array_equal(arr[:, nd.shape[1]:], expected, equal_nan=True),
                 'Exact-HKL original HL join/missing NaN policy differs')
        merge = rr['phase_merge_audit']
        _require(merge['native'] == inputs['native'] and merge['phase_source'] == phase and
                 merge['representation'] == 'hl' and merge['output'] == rr['input'], 'Prior merge audit provenance differs')
        _require(merge['matched_prior_rows'] == int((~missing).sum()) and
                 merge['missing_prior_rows'] == int(missing.sum()), 'Prior finite/missing row counts differ')
        parsed = _parse_refmac(refine / 'refmac.log', refine / 'refmac.xml')
        _require(all(_equal(value, rr['statistics'].get(key)) for key, value in parsed.items()),
                 'REFMAC recorded statistics differ from raw log/XML')
        _require(_equal(stage['statistics'], parsed['final']), 'Round R factors differ from raw log')
        _require(rr['free_flags_preserved'] and stage['free_flags_preserved'], 'Missing original free-set preservation')
        free_audit = _free_flags(native, gemmi.read_mtz_file(str(refine / 'refined.mtz')))
        rounds.append({'round': number, 'buccaneer_cycles': stage['buccaneer_cycles'],
                       'refmac_cycles': 10, 'statistics': parsed['final'],
                       'raw_native_columns_and_order_unchanged': True,
                       'same_original_HL_exactHKL_missing_NaN': True,
                       'free_flag_audit': free_audit, 'own_previous_map': _fp(previous_map),
                       'own_previous_model': _fp(previous_model),
                       'construction_command': _fp(build / 'command.json'),
                       'refinement_report': _fp(refine / 'refinement_report.json')})
        previous_map, previous_model = refine / 'refined.mtz', refine / 'model.cif'
    _require(primary.name == 'primary_model.cif' and _sha(primary) == _sha(previous_model),
             'Primary model is not fixed last-round output')
    _require(report['primary_model'] == _fp(primary), 'Primary model recorded fingerprint differs')
    _require(_equal(state['primary_metric'], report['primary_metric']), 'State final R factors differ')
    return {'id': case['id'], 'applicable': True, 'passed': True, 'audit_passed': True,
            'auditor_revision': AUDITOR_REVISION, 'issues': [], 'rounds': rounds,
            'primary_model': _fp(primary), 'primary_metric': report['primary_metric'],
            'input_fingerprints': dict(inputs, phase_source=phase, seed_map=seed_map, seed_model=seed),
            'prior_rows_finite': int((~missing).sum()), 'prior_rows_missing_NaN': int(missing.sum()),
            'fixed_last_round_verified': True, 'reference_read': False,
            'scientific_sources_unchanged': True, 'total_refmac_cycles': 50}


def _audit_no_free(case, state, run):
    folder = Path(run).resolve() / 'cases' / case['id']
    _require(state['status'] == 'completed_no_free' and not case.get('has_free'),
             'No-free branch has inconsistent status')
    _require(state.get('primary_metric') is None and state.get('r_free') is None,
             'No-free case cannot have an Rfree metric')
    _require(not any(folder.rglob('refinement_report.json')) and not any(folder.rglob('refmac.xml')),
             'No-free case unexpectedly contains refinement evidence')
    build = state.get('build_directory')
    if build is None and state.get('build_report'):
        candidate = state['build_report']
        build = Path(candidate['path'] if isinstance(candidate, dict) else candidate).parent
    if build is None:
        options = [folder / 'build', folder / 'pipeline/round_01_build', folder / 'one_pass_build']
        available = [p for p in options if (p / 'command.json').is_file()]
        _require(len(available) == 1, 'Missing/ambiguous one-pass construction directory')
        build = available[0]
    build = Path(build).resolve()
    _require(_contained(build, folder), 'No-free construction escapes its case')
    command = _read(build / 'command.json')
    initial = _declared(case, 'initial_mtz', 'initial_sha256')
    _declared(case, 'sequence', 'sequence_sha256')
    _declared(case, 'phase_source', 'phase_sha256')
    _declared(case, 'seed_map', 'seed_map_sha256')
    seed = _flag(command, '-pdbin')
    _require(_contained(seed, folder), 'No-free initial seed outside case')
    _require(state['seed_model_fingerprint'] == _fp(seed) and
             str(Path(state['seed_model']).resolve()) == str(Path(seed).resolve()),
             'No-free construction seed differs from frozen state')
    binary = _fp(command[0])
    library = [_fp(_flag(command, name)) for name in ['-mtzin-ref', '-pdbin-ref']]
    for item in [binary] + library:
        _require(_contained(item['path'], BUNDLE / 'runtime'), 'No-free construction outside bundled runtime')
    _check_build(build, command, case, initial['path'], seed, binary, library, 1)
    other = list(folder.rglob('buccaneer.xml'))
    _require(other == [build / 'buccaneer.xml'] or set(other) == {build / 'buccaneer.xml'},
             'No-free case must have exactly one three-cycle construction')
    primary = Path(state['primary_model']).resolve()
    _require(_contained(primary, folder) and _sha(primary) == _sha(build / 'model.cif'),
             'No-free primary model is not its only construction output')
    return {'id': case['id'], 'applicable': False, 'passed': True, 'audit_passed': True,
            'auditor_revision': AUDITOR_REVISION, 'issues': [], 'reference_read': False,
            'classification': 'completed_no_original_free_set', 'primary_model': _fp(primary),
            'buccaneer_builds': 1, 'buccaneer_cycles': 3, 'total_refmac_cycles': 0,
            'primary_metric': None, 'r_free': None,
            'strict_success_eligible': False,
            'limitation': 'No supplied original free set; one-pass construction only; no regenerated free flags or invented Rfree'}


def audit_completed(case, state, run):
    """Return read-only protocol evidence; never read a reference or edit state.

    The caller freezes/verifies the common runtime inventory once.  This
    function additionally checks only binaries/resources actually used by the
    endpoint and all recorded scientific source files.
    """
    applicable = bool(case.get('has_free', True))
    try:
        return _audit_free(case, state, run) if applicable else _audit_no_free(case, state, run)
    except Exception as error:
        return {'id': case.get('id'), 'applicable': applicable, 'passed': False, 'audit_passed': False,
                'auditor_revision': AUDITOR_REVISION, 'reference_read': False,
                'issues': [str(error)], 'traceback': traceback.format_exc()}
