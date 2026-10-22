"""Read-only scientific audit of the fixed four-case portable regression.

This script does not launch programs or construct models.
Until all four cases have reached a terminal state it exits 2 without writing
results. Once terminal, it saves immutable package-manifest snapshots and writes
audit.json, comparison.json and REPORT.md inside the validation directory.
The explicit --finalize-delivery-metadata option corrects only four known app
manifest display lengths after successful scientific auditing.
Failed executions retain their place
in the denominator and are distinguished from scientifically unsuccessful
completed models. Run using a Python with Gemmi/NumPy/SciPy, preferably the
bundled Python with -B. Target references are read only after a completed fixed endpoint.
"""
import sys
sys.dont_write_bytecode = True

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import traceback
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
IDS = ('5hhk', '7bll', '7tem', '7twc')
TERMINAL = {'completed', 'failed', 'timeout', 'cancelled'}
PRIOR_LABELS = ('HLA', 'HLB', 'HLC', 'HLD')
CRITERIA = {'r_free_max': .35, 'n_ca_c_cutoff_angstrom': 1.,
            'n_ca_c_completeness_min': .8, 'n_ca_c_precision_min': .9,
            'joint_sequence_backbone_coverage_min': .8}
KNOWN_APP_SIZE_METADATA = {
    'src/phasebuild/hybrid_recycle.py': (8705, 8727),
    'src/phasebuild/phase_refinement.py': (12671, 12754),
    'src/phasebuild/refinement.py': (10579, 10657),
    'src/phasebuild/runtime.py': (2855, 3428),
}


def read(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': sha256(path)}


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def checked_fingerprint(item):
    require(fingerprint(item['path']) == item, 'File fingerprint changed: ' + item['path'])
    return item


def check_nested_fingerprints(value):
    """Include any additional execution-failure evidence in the recovery audit."""
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            checked_fingerprint(value)
            return 1
        return sum(check_nested_fingerprints(item) for item in value.values())
    if isinstance(value, list):
        return sum(check_nested_fingerprints(item) for item in value)
    return 0


def contained(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def numeric_equal(a, b):
    """Metric JSON may encode NaN; all finite values must agree exactly."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(numeric_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(numeric_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def flag(command, name):
    require(command.count(name) == 1, 'Missing or duplicate command flag: ' + name)
    return command[command.index(name) + 1]


def scientific_construction_command(command):
    """Compare all science flags while permitting the declared relocated paths."""
    result = list(command)
    result[0] = '<binary>'
    for name in ['-mtzin', '-mtzin-ref', '-pdbin-ref', '-pdbin', '-seqin', '-pdbout', '-xmlout']:
        flag(result, name)
        result[result.index(name) + 1] = '<' + name + '>'
    return result


def additional_execution_history(history):
    """Verify later startup failures without treating them as model failures."""
    records = []
    for attempt in history.get('additional_failed_execution_attempts', []):
        directory = Path(attempt['directory']).resolve()
        for key in ['state', 'pipeline_report', 'refinement_report', 'refinement_log',
                    'app_manifest_snapshot', 'driver_snapshot']:
            checked_fingerprint(attempt[key])
        state, pipeline, refinement = [read(attempt[key]['path']) for key in
                                       ['state', 'pipeline_report', 'refinement_report']]
        require(state['status'] == pipeline['status'] == refinement['status'] == 'failed',
                'Additional historical execution is not a failed startup')
        require(pipeline['stages'] == [] and pipeline['current_round'] == 1
                and pipeline['current'] == 'refmac', 'Additional failure had completed scientific stages')
        require(refinement['returncode'] != 0 and refinement['requested_cycles'] == 10
                and not refinement.get('statistics', {}).get('cycles'), 'Additional failure contains completed REFMAC cycles')
        log = Path(attempt['refinement_log']['path']).read_text()
        require('Cannot get CCP4_SCR from environment' in log and 'Ncyc    Rfact    Rfree' not in log,
                'Additional failure is not the documented pre-cycle scratch-environment failure')
        xml_path = Path(attempt['refinement_report']['path']).parent / 'refmac.xml'
        if xml_path.exists():
            require(not ET.parse(str(xml_path)).getroot().findall('./Overall_stats/stats_vs_cycle/new_cycle'),
                    'Additional failed execution XML contains a REFMAC cycle')
        old_manifest = read(directory.parent / 'manifest.json')
        app_manifest = read(attempt['app_manifest_snapshot']['path'])
        app_record = next(row for row in old_manifest['manifests'] if Path(row['path']).parent.name == 'app')
        require(attempt['app_manifest_snapshot']['sha256'] == app_record['sha256']
                and attempt['driver_snapshot']['sha256'] == old_manifest['driver']['sha256'],
                'Additional failed execution snapshot differs from its recorded package/driver')
        app_root = Path(attempt['app_manifest_snapshot']['path']).parent
        for row in app_manifest['files']:
            require(sha256(app_root / row['relative_path']) == row['sha256'],
                    'Additional old app snapshot changed: ' + row['relative_path'])
        middle_path = directory.parent / 'before_scratch_repair_runtime_manifest.json'
        middle_record = next(row for row in old_manifest['manifests'] if Path(row['path']).parent.name == 'runtime')
        require(middle_path.is_file() and sha256(middle_path) == middle_record['sha256'],
                'Intermediate pre-scratch runtime manifest was not preserved exactly')
        records.append({'directory': str(directory), 'state': attempt['state'],
                        'pipeline_report': attempt['pipeline_report'],
                        'refinement_report': attempt['refinement_report'],
                        'refinement_log': attempt['refinement_log'],
                        'app_manifest_snapshot': attempt['app_manifest_snapshot'],
                        'driver_snapshot': attempt['driver_snapshot'],
                        'app_snapshot_files_checked': len(app_manifest['files']),
                        'before_scratch_runtime_manifest': fingerprint(middle_path),
                        'refmac_cycles_executed': 0, 'classification': 'execution_startup_failure',
                        'not_counted_as_scientific_failure': True})
    return records


def recovery_chain(cid, output, bundle, trial_manifest, state, report):
    """Verify an execution-only recovery without rerunning a completed build."""
    resume_root = trial_manifest.get('execution_recovery_from')
    recovered = report.get('resumed_first_construction')
    if resume_root is None:
        require(recovered is None and '--existing-first-build' not in state['command'],
                'Undeclared reuse of an existing first construction')
        return None
    previous_output = Path(resume_root).resolve()
    source = previous_output / cid / 'pipeline/round_01_build'
    destination = output / cid / 'pipeline/round_01_build'
    history = read(output / cid / 'failure_history.json')
    history_fingerprint_count = check_nested_fingerprints(history)
    additional_failures = additional_execution_history(history)
    require(history['previous_failed_output'] == str(previous_output / cid), 'Wrong recovery failure directory')
    for key in ['previous_failed_state', 'previous_failed_pipeline_report', 'previous_failed_refinement_log',
                'previous_failed_refinement_report', 'first_construction_model', 'first_construction_command',
                'old_runtime_manifest', 'old_app_manifest', 'old_driver']:
        checked_fingerprint(history[key])
    previous_state = read(history['previous_failed_state']['path'])
    previous_report = read(history['previous_failed_pipeline_report']['path'])
    previous_refine = read(history['previous_failed_refinement_report']['path'])
    require(previous_state['status'] == previous_report['status'] == previous_refine['status'] == 'failed',
            'Recovery source is not a preserved failed execution')
    require(previous_report['stages'] == [] and previous_report['current_round'] == 1
            and previous_report['current'] == 'refmac', 'Recovery would repeat a completed refinement')
    require(previous_report['inputs'] == report['inputs'] and previous_state['inputs'] == state['inputs']
            and previous_report['parameters'] == report['parameters'], 'Recovery changed original construction inputs/protocol')
    require(previous_report['binary'] == report['binary']
            and previous_report['training_library'] == report['training_library'], 'Recovery changed construction engine/data')
    require('Cannot open environ.def' in Path(history['previous_failed_refinement_log']['path']).read_text(),
            'Preserved previous execution has an unexpected failure cause')
    require(previous_refine['returncode'] != 0 and previous_refine['requested_cycles'] == 10,
            'Previous refinement was already completed or uses different cycles')
    require(flag(state['command'], '--existing-first-build') == str(source), 'Launcher recovery source differs')
    require(recovered and recovered['construction_executed_again'] is False
            and recovered['source'] == str(source), 'First construction was not recorded as reused exactly once')
    require(recovered['source_report'] == history['previous_failed_pipeline_report'], 'Recovery source report differs')
    require(recovered['original_model'] == history['first_construction_model']
            and recovered['copied_model'] == fingerprint(destination / 'model.cif'), 'Recovered model provenance differs')
    require(recovered['original_model']['sha256'] == recovered['copied_model']['sha256'], 'Recovered model bytes changed')
    expected_files = {str(p.relative_to(source)) for p in source.rglob('*') if p.is_file()}
    require(expected_files == {str(p.relative_to(destination)) for p in destination.rglob('*') if p.is_file()}
            == {row['relative_path'] for row in recovered['files']}, 'Recovery file inventories differ')
    for row in recovered['files']:
        relative = row['relative_path']
        require(contained(source / relative, source) and contained(destination / relative, destination),
                'Recovered construction file escapes case directory')
        require(row['original'] == fingerprint(source / relative)
                and row['copied'] == fingerprint(destination / relative)
                and row['original']['sha256'] == row['copied']['sha256'], 'Recovered construction file changed: ' + relative)
    require(fingerprint(destination / 'command.json')['sha256'] == history['first_construction_command']['sha256'],
            'Recovered command was rewritten instead of preserved')
    # The old manifest points to the then-current package. Check preserved
    # snapshots against those recorded digests, not against the repaired files.
    old_trial = read(previous_output / 'manifest.json')
    manifest_hashes = {Path(item['path']).parent.name: item['sha256'] for item in old_trial['manifests']}
    require(history['old_runtime_manifest']['sha256'] == manifest_hashes['runtime']
            and history['old_app_manifest']['sha256'] == manifest_hashes['app']
            and history['old_driver']['sha256'] == old_trial['driver']['sha256'], 'Previous package/driver snapshots do not match old run manifest')
    snapshot_root = Path(history['old_driver']['path']).parent
    old_app = read(history['old_app_manifest']['path'])
    for row in old_app['files']:
        require(sha256(snapshot_root / row['relative_path']) == row['sha256'],
                'Old app source snapshot changed: ' + row['relative_path'])
    old_runtime = read(history['old_runtime_manifest']['path'])
    new_runtime = read(bundle / 'runtime/runtime_manifest.json')
    old_members = {row['relative_path']: row for row in old_runtime['files']}
    new_members = {row['relative_path']: row for row in new_runtime['files']}
    require(not (old_members.keys() - new_members.keys()), 'Recovery removed existing runtime files')
    changed = sorted(k for k in old_members if old_members[k]['sha256'] != new_members[k]['sha256'])
    require(set(changed) <= {'README_LOCAL_EVALUATION.txt'}, 'Recovery altered an existing scientific runtime file')
    additions = sorted(new_members.keys() - old_members.keys())
    expected_additions = {'include/default.def', 'include/environ.def'} | {
        'share/startup_checks/' + repair + '/' + program + '/' + filename
        for repair in ['include_definition_repair_20261020', 'scratch_environment_repair_20261020']
        for program in ['cbuccaneer', 'refmac5'] for filename in ['result.json', 'stdout.txt', 'stderr.txt']}
    if additional_failures:
        require(set(additions) == expected_additions,
                'Unexpected runtime additions beyond include/scratch repair evidence')
    else:
        include_additions = {x for x in expected_additions if 'scratch_environment_repair' not in x}
        require(set(additions) in (include_additions, expected_additions),
                'Unexpected runtime additions beyond include repair and startup evidence')
    return {'source': str(source), 'first_construction_executed_again': False,
            'copied_files_checked': len(expected_files), 'failure_history': fingerprint(output / cid / 'failure_history.json'),
            'old_runtime_manifest': history['old_runtime_manifest'], 'old_app_manifest': history['old_app_manifest'],
            'old_driver': history['old_driver'], 'old_app_snapshot_files_checked': len(old_app['files']),
            'failure_history_fingerprints_checked': history_fingerprint_count,
            'additional_failed_execution_attempts': additional_failures,
            'runtime_added_files': additions, 'runtime_changed_documentation': changed,
            'old_run_manifest': fingerprint(previous_output / 'manifest.json'),
            'old_failed_outputs_preserved': True, 'recovery_scientific_protocol_unchanged': True}


def only_known_size_metadata_corrections(original, delivered):
    expected = json.loads(json.dumps(original))
    corrected = set()
    for row in expected['files']:
        if row['relative_path'] in KNOWN_APP_SIZE_METADATA:
            old, new = KNOWN_APP_SIZE_METADATA[row['relative_path']]
            require(row['bytes'] == old, 'Unexpected original manifest display size')
            row['bytes'] = new
            corrected.add(row['relative_path'])
    require(corrected == set(KNOWN_APP_SIZE_METADATA), 'Incomplete known manifest metadata correction')
    return expected == delivered


def freeze_manifests(output, bundle, trial_manifest):
    """Preserve run-time bytes before any delivery-only metadata correction."""
    frozen_root = output / 'frozen_package_metadata'
    records = {}
    for item in trial_manifest['manifests']:
        relative = Path(item['path']).resolve().relative_to(bundle)
        destination = frozen_root / relative
        if destination.exists():
            require(sha256(destination) == item['sha256'], 'Frozen run manifest changed: ' + str(destination))
        else:
            require(sha256(item['path']) == item['sha256'], 'Cannot snapshot changed run manifest')
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open('xb') as stream:
                stream.write(Path(item['path']).read_bytes())
            require(sha256(destination) == item['sha256'], 'Manifest snapshot did not preserve original bytes')
        records[item['path']] = fingerprint(destination)
    return records


def freeze_launch_metadata(output, bundle):
    """Preserve launchers and release bookkeeping without changing their bytes."""
    root = output / 'frozen_package_metadata'
    path = root / 'launch_metadata_snapshot.json'
    if path.exists():
        snapshot = read(path)
        for row in snapshot['files']:
            checked_fingerprint(row['snapshot'])
            if row['relative_path'] != 'bundle.json':
                checked_fingerprint(row['original'])
        return fingerprint(path)
    records = []
    for name in ['phasebuild', 'PhaseBuild.command', 'bundle.json']:
        original = fingerprint(bundle / name)
        destination = root / ('launchers/' + name if name != 'bundle.json' else name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            require(sha256(destination) == original['sha256'], 'Existing launch metadata snapshot differs')
        else:
            with destination.open('xb') as stream:
                stream.write((bundle / name).read_bytes())
        require(sha256(destination) == original['sha256'], 'Launch metadata bytes changed during snapshot')
        records.append({'relative_path': name, 'original': original, 'snapshot': fingerprint(destination)})
    snapshot = {'captured_unix': time.time(), 'files': records,
                'scope': 'Captured during the frozen validation run; bundle.json is release bookkeeping and may be updated after audit, while launchers must remain unchanged'}
    with path.open('x') as stream:
        stream.write(json.dumps(snapshot, indent=2) + '\n')
    return fingerprint(path)


def inventory(bundle, trial_manifest, frozen):
    """Verify every recorded app, Python, binary, library, and data file."""
    expected = {str((bundle / name).resolve()): name for name in
                ('runtime/runtime_manifest.json', 'python/python_manifest.json', 'app/app_manifest.json')}
    recorded = trial_manifest['manifests']
    require({x['path'] for x in recorded} == set(expected), 'Wrong or incomplete frozen package manifests')
    records = []
    for item in recorded:
        frozen_record = frozen[item['path']]
        checked_fingerprint(frozen_record)
        manifest = read(frozen_record['path'])
        current = fingerprint(item['path'])
        current_metadata_corrected = False
        if current['sha256'] != item['sha256']:
            require(Path(item['path']).parent.name == 'app'
                    and only_known_size_metadata_corrections(manifest, read(item['path'])),
                    'Run manifest changed beyond documented delivery-size metadata correction')
            current_metadata_corrected = True
        root = Path(item['path']).parent
        seen = set()
        total = 0
        metadata_warnings = []
        for row in manifest['files']:
            rel = row['relative_path']
            require(rel not in seen, 'Duplicate inventory member: ' + rel)
            seen.add(rel)
            path = root / rel
            require(contained(path, root), 'Inventory member escapes package: ' + rel)
            require(sha256(path) == row['sha256'], 'Package member changed: ' + str(path))
            size = row.get('bytes', row.get('size'))
            if size is not None and path.stat().st_size != size:
                actual = path.stat().st_size
                require(root.name == 'app' and rel in KNOWN_APP_SIZE_METADATA
                        and (size, actual) == KNOWN_APP_SIZE_METADATA[rel],
                        'Package member size differs beyond four documented display-only fields: ' + str(path))
                metadata_warnings.append({'relative_path': rel, 'frozen_declared_bytes': size,
                                          'actual_bytes': actual, 'sha256_verified': True,
                                          'classification': 'manifest_display_length_only'})
            total += path.stat().st_size
        for row in manifest.get('symlinks', []):
            path = root / row['relative_path']
            require(path.is_symlink() and os.readlink(path) == row['target'],
                    'Package symlink changed: ' + str(path))
            require(contained(path, root), 'Package symlink escapes package: ' + str(path))
        records.append({'manifest': item, 'frozen_manifest': frozen_record, 'current_manifest': current,
                        'files_checked': len(seen), 'bytes_checked': total,
                        'symlinks_checked': len(manifest.get('symlinks', [])), 'passed': True,
                        'scientific_file_hashes_passed': True, 'metadata_warnings': metadata_warnings,
                        'current_manifest_has_only_documented_size_corrections': current_metadata_corrected})
    return records


def structure_metrics(evaluation, metric):
    n = evaluation['reference_residues']
    joint = evaluation['sequence_evaluation']['correct_type_and_backbone'] / n
    success = (metric['r_free'] <= CRITERIA['r_free_max']
               and evaluation['n_ca_c_completeness'] >= CRITERIA['n_ca_c_completeness_min']
               and evaluation['n_ca_c_precision'] >= CRITERIA['n_ca_c_precision_min']
               and joint >= CRITERIA['joint_sequence_backbone_coverage_min'])
    return {'r_work': metric['r_work'], 'r_free': metric['r_free'],
            'built_residues': evaluation['built_residues'], 'reference_residues': n,
            'n_ca_c_correct_residues': evaluation['n_ca_c_correct_residues'],
            'n_ca_c_completeness': evaluation['n_ca_c_completeness'],
            'n_ca_c_precision': evaluation['n_ca_c_precision'],
            'correct_type_and_backbone': evaluation['sequence_evaluation']['correct_type_and_backbone'],
            'joint_sequence_backbone_coverage': joint, 'strict_success': bool(success)}


def audit_case(cid, output, bundle, case, state, trial_manifest, np, gemmi, parse_statistics, verify_free_flags, evaluate):
    """Audit a completed case. Any audit failure is separate from scientific failure."""
    folder = output / cid
    pipeline = folder / 'pipeline'
    old_folder = ROOT / 'runs/sad_research_v2_20261019/full_hybrid_hl_recovered' / cid
    old = read(old_folder / 'pipeline/hybrid_report.json')
    report = read(pipeline / 'hybrid_report.json')
    require(state['status'] == 'completed' and report['status'] == 'completed', 'Incomplete fixed endpoint')
    require(report['parameters'] == old['parameters'] == {
        'rounds': 5, 'refine_cycles': 10, 'build_semet': bool(case['build_semet']), 'model_filter': True},
        'Scientific parameters differ from declared originalHL control')
    require(report['inputs'] == old['inputs'], 'Construction inputs differ from previous fixedHL workflow')
    for item in report['inputs'].values():
        checked_fingerprint(item)
    declared = {'map': (case['initial_mtz'], case['initial_sha256']),
                'native': (case['native'], case['native_sha256']),
                'sequence': (case['sequence'], case['sequence_sha256']),
                'phase_prior': (case['phase_source'], case['phase_sha256']),
                'model': (old['inputs']['phasebuild_seed_model']['path'],
                          old['inputs']['phasebuild_seed_model']['sha256'])}
    require(set(state['inputs']) == set(declared), 'Missing/extra validation construction inputs')
    for key, (path, digest) in declared.items():
        require(state['inputs'][key] == {'path': str(Path(path).resolve()), 'sha256': digest},
                'Declared input differs: ' + key)
        checked_fingerprint(state['inputs'][key])
    command = state['command']
    require(command[:2] == ['/usr/bin/sandbox-exec', '-p'], 'Sandbox launcher is missing')
    require(command[2] == '(version 1)(allow default)(deny file-read* (subpath "/Applications/ccp4-9"))',
            'Original CCP4 installation is not wholly denied')
    require(command[3:5] == [str(bundle / 'phasebuild'), 'hybrid-hl'], 'Wrong portable entrypoint')
    require(state['original_ccp4_read_denied'] is True, 'Missing installation-denial declaration')
    for name, key in [('--map', 'map'), ('--native', 'native'), ('--sequence', 'sequence'),
                      ('--model', 'model'), ('--phase-prior', 'phase_prior')]:
        require(flag(command, name) == state['inputs'][key]['path'], 'Launcher input substitution: ' + name)
    require(flag(command, '--rounds') == '5' and flag(command, '--cycles') == '10', 'Launcher cycle changes')
    require(flag(command, '--phase-prior-sha256') == case['phase_sha256'], 'Launcher prior hash changed')
    require(flag(command, '--output') == str(pipeline), 'Launcher output path changed')
    require(('--build-semet' in command) == bool(case['build_semet']), 'Experiment SeMet metadata changed')
    recovery = recovery_chain(cid, output, bundle, trial_manifest, state, report)
    require(not report['source_changed_during_run'], 'Construction sources changed during execution')
    require(len(report['stages']) == 5, 'Not exactly five completed construction/refinement stages')
    require(report['primary_metric'] == report['stages'][-1]['statistics'], 'Selected metric is not final round')
    for item in report['sources']:
        require(contained(item['path'], bundle / 'app/src'), 'Imported construction source outside package')
        checked_fingerprint(item)
    require(contained(report['binary']['path'], bundle / 'runtime'), 'Buccaneer is outside portable runtime')
    checked_fingerprint(report['binary'])
    require(report['binary']['sha256'] == old['binary']['sha256'], 'Buccaneer binary differs from fixedHL control')
    require(len(report['training_library']) == 2, 'Wrong training resources')
    for item, old_item in zip(report['training_library'], old['training_library']):
        require(contained(item['path'], bundle / 'runtime'), 'Training data outside portable runtime')
        checked_fingerprint(item)
        require(item['sha256'] == old_item['sha256'], 'Training resource differs from fixedHL control')

    native = gemmi.read_mtz_file(case['native'])
    source = gemmi.read_mtz_file(case['phase_source'])
    native_data = np.array(native, copy=False)
    native_row_index = {tuple(row[:3].astype(int)): i for i, row in enumerate(native_data)}
    require(len(native_row_index) == len(native_data), 'Duplicate native HKLs')
    source_data = np.array(source, copy=False)
    cols = [source.column_labels().index(x) for x in PRIOR_LABELS]
    lookup = {tuple(row[:3].astype(int)): row[cols] for row in source_data}
    require(len(lookup) == len(source_data), 'Duplicate source phase HKLs')
    require(np.all(source_data[:, :3] == np.rint(source_data[:, :3])), 'Noninteger source phase HKLs')
    expected = np.full((len(native_data), 4), np.nan, dtype=np.float32)
    for i, row in enumerate(native_data):
        prior = lookup.get(tuple(row[:3].astype(int)))
        if prior is not None and np.isfinite(prior).all():
            expected[i] = prior
    missing = ~np.isfinite(expected).all(axis=1)
    require(missing.any(), 'Expected missing high-resolution experimental priors were not found')
    previous_map = Path(case['initial_mtz']).resolve()
    previous_model = Path(old['inputs']['phasebuild_seed_model']['path']).resolve()
    rounds = []
    for number, stage in enumerate(report['stages'], 1):
        build = pipeline / ('round_%02d_build' % number)
        refine = pipeline / ('round_%02d_refine' % number)
        bc = read(build / 'command.json')
        old_bc = read(old_folder / 'pipeline' / ('round_%02d_build' % number) / 'command.json')
        rr = read(refine / 'refinement_report.json')
        old_rr = read(old_folder / 'pipeline' / ('round_%02d_refine' % number) / 'refinement_report.json')
        require(stage['round'] == number and stage['buccaneer_cycles'] == (3 if number == 1 else 2)
                and stage['refmac_cycles'] == 10 and stage['construction'] == 'external-Buccaneer',
                'Wrong stage/cycle schedule in round ' + str(number))
        require(bc[0] == report['binary']['path'], 'Construction engine path substitution')
        require(scientific_construction_command(bc) == scientific_construction_command(old_bc),
                'Buccaneer science flags differ from fixedHL control')
        require(Path(flag(bc, '-mtzin')).resolve() == previous_map, 'Construction did not use its own preceding map')
        require(Path(flag(bc, '-pdbin')).resolve() == previous_model, 'Construction did not use its own preceding model')
        require(flag(bc, '-mtzin-ref') == report['training_library'][0]['path']
                and flag(bc, '-pdbin-ref') == report['training_library'][1]['path']
                and flag(bc, '-seqin') == report['inputs']['sequence']['path'], 'Construction training/sequence path changed')
        output_origin = Path(recovery['source']) if number == 1 and recovery else build
        require(flag(bc, '-pdbout') == str(output_origin / 'model.pdb')
                and flag(bc, '-xmlout') == str(output_origin / 'buccaneer.xml'), 'Construction output paths do not match executed or preserved command')
        require(flag(bc, '-cycles') == ('3' if number == 1 else '2'), 'Buccaneer cycle flag changed')
        require(flag(bc, '-jobs') == '1', 'Buccaneer jobs changed')
        require(('-no-correlation-mode' in bc) == (number == 1)
                and ('-correlation-mode' in bc) == (number > 1), 'Wrong Buccaneer likelihood target')
        require(('-build-semet' in bc) == bool(case['build_semet']), 'SeMet flag changed')
        for name in ['find', 'grow', 'join', 'link', 'sequence', 'correct', 'filter', 'ncsbuild', 'prune', 'rebuild', 'model-filter']:
            require('-' + name in bc, 'Missing construction stage: ' + name)
        xml = ET.parse(str(build / 'buccaneer.xml')).getroot()
        require([int(c.findtext('Number')) for c in xml.findall('./Cycles/Cycle')]
                == list(range(1, 4 if number == 1 else 3)), 'Buccaneer XML cycle count differs')
        require('cbuccaneer: Normal termination' in (build / 'buccaneer.log').read_text(),
                'Buccaneer did not report normal termination')
        require(rr['status'] == 'success' and rr['returncode'] == 0 and rr['requested_cycles'] == 10,
                'Refinement did not finish requested ten cycles')
        require(rr['phase_target_verified'] and rr['phase_modified_likelihood_active']
                and rr['use_experimental_prior'] and rr['prior_representation'] == 'hl', 'Wrong phase likelihood')
        require(rr['fixed_original_experimental_prior'] == state['inputs']['phase_prior'], 'Prior changed between rounds')
        require(rr['original_native_input'] == state['inputs']['native'], 'Native input changed between rounds')
        require(rr['keywords'] == old_rr['keywords'] == (refine / 'keywords.txt').read_text(),
                'REFMAC scientific keywords differ from fixedHL control')
        require(rr['command'][1:] == ['HKLIN', 'input.mtz', 'XYZIN', 'input.pdb',
                'HKLOUT', 'refined.mtz', 'XYZOUT', 'refined.pdb', 'XMLOUT', 'refmac.xml'], 'Wrong staged REFMAC arguments')
        require(contained(rr['binary']['path'], bundle / 'runtime'), 'REFMAC outside portable runtime')
        checked_fingerprint(rr['binary'])
        require(rr['binary']['sha256'] == old_rr['binary']['sha256'], 'REFMAC binary differs from fixedHL control')
        for name in ['source', 'wrapper_source', 'model_input', 'input', 'staged_input', 'output_model', 'output_mtz']:
            checked_fingerprint(rr[name])
        require(rr['model_input'] == fingerprint(build / 'model.cif'), 'Refinement used another model')
        require(rr['staged_input'] == fingerprint(refine / 'input.mtz')
                and rr['input']['sha256'] == rr['staged_input']['sha256'], 'Staged MTZ bytes changed')
        staged = gemmi.read_mtz_file(str(refine / 'input.mtz'))
        arr = np.array(staged, copy=False)
        require(staged.column_labels()[:len(native.column_labels())] == native.column_labels(), 'Native columns renamed/reordered')
        require(np.array_equal(arr[:, :native_data.shape[1]], native_data, equal_nan=True),
                'Native observations/HKL order/free flags changed')
        require(np.array_equal(arr[:, [staged.column_labels().index('SAD_' + s) for s in PRIOR_LABELS]],
                               expected, equal_nan=True), 'Experimental HL exact-HKL join or missing NaNs differ')
        parsed = parse_statistics(refine / 'refmac.log', refine / 'refmac.xml')
        require([x['cycle'] for x in parsed['cycles']] == list(range(11)), 'REFMAC raw cycle schedule differs')
        require(all(numeric_equal(value, rr['statistics'].get(key)) for key, value in parsed.items()),
                'REFMAC raw statistics differ from report')
        require(stage['statistics'] == parsed['final'], 'Stage R factors differ from raw log')
        refined = gemmi.read_mtz_file(str(refine / 'refined.mtz'))
        free_audit = verify_free_flags(native, refined, 'FreeR_flag', 'F', 'SIGF')
        require(free_audit['all_retained_flags_unchanged'] and rr['free_flags_preserved'], 'Original free flags changed')
        refined_data = np.array(refined, copy=False)
        retained_indices = [native_row_index[tuple(row[:3].astype(int))] for row in refined_data]
        refined_observations = refined_data[:, [refined.column_labels().index(x) for x in ['F', 'SIGF']]]
        original_observations = native_data[retained_indices][:, [native.column_labels().index(x) for x in ['F', 'SIGF']]]
        # REFMAC may apply its normal scaling to output F/SIGF. Input staging
        # above remains bitwise unchanged, and the original free partition is
        # checked independently by exact HKL. Output scaling is recorded only.
        same_observations = ((refined_observations == original_observations)
                             | (np.isnan(refined_observations) & np.isnan(original_observations)))
        changed_observation_rows = int((~same_observations.all(axis=1)).sum())
        rounds.append({'round': number, 'buccaneer_cycles': stage['buccaneer_cycles'],
                       'refmac_cycles': 10, 'statistics': parsed['final'], 'same_original_HL': True,
                       'raw_native_columns_unchanged': True, 'exact_HKL_HL_with_missing_NaN': True,
                       'output_amplitude_or_sigma_changed_rows': changed_observation_rows,
                       'output_scaling_policy': 'REFMAC output F/SIGF scaling allowed; exact original input columns and raw-HKL free flags mandatory',
                       'free_flags_preserved': True, 'previous_own_map': fingerprint(previous_map),
                       'construction_command': fingerprint(build / 'command.json'),
                       'refinement_report': fingerprint(refine / 'refinement_report.json'),
                       'input_mtz': fingerprint(refine / 'input.mtz'),
                       'output_mtz': fingerprint(refine / 'refined.mtz')})
        previous_map, previous_model = refine / 'refined.mtz', refine / 'model.cif'
    require(fingerprint(pipeline / 'primary_model.cif')['sha256'] == fingerprint(previous_model)['sha256'],
            'Primary model is not fixed final round')
    require(report['primary_model'] == fingerprint(pipeline / 'primary_model.cif'), 'Primary model fingerprint differs')
    require(state['primary_metric'] == report['primary_metric'], 'State R factors differ')
    saved = read(folder / 'evaluation.json')
    reference = ROOT / 'runs/buccaneer_sad100_20261019/evaluation' / cid / 'protein_reference.cif'
    require(saved['model'] == report['primary_model'] and saved['reference'] == fingerprint(reference),
            'Evaluation model/reference differs')
    require(saved['cutoff_angstrom'] == 1., 'Strict evaluation distance threshold changed')
    # Recompute posterior evaluation only; construction is already completed.
    recomputed = evaluate(pipeline / 'primary_model.cif', reference)
    fields = ['built_residues', 'reference_residues', 'n_ca_c_correct_residues',
              'n_ca_c_completeness', 'n_ca_c_precision', 'sequence_evaluation']
    require(all(numeric_equal(saved[k], recomputed[k]) for k in fields), 'Stored structure metrics fail recomputation')
    current_metrics = structure_metrics(saved, report['primary_metric'])
    require(state['strict_success'] == current_metrics['strict_success'], 'Stored success differs from fixed criteria')
    old_evaluation = read(old_folder / 'evaluation.json')
    checked_fingerprint(old['primary_model'])
    require(old_evaluation['reference'] == fingerprint(reference), 'Paired control reference differs')
    previous_metrics = structure_metrics(old_evaluation, old['primary_metric'])
    return {'id': cid, 'execution_status': 'completed', 'classification':
            'strict_success' if current_metrics['strict_success'] else 'scientific_not_solved',
            'portable': current_metrics, 'previous_fixed_HL': previous_metrics,
            'r_free_delta': current_metrics['r_free'] - previous_metrics['r_free'],
            'r_work_delta': current_metrics['r_work'] - previous_metrics['r_work'],
            'reference': fingerprint(reference), 'primary_model': report['primary_model'],
            'original_ccp4_read_denied': True, 'rounds': rounds,
            'execution_recovery': recovery,
            'prior_rows_finite': int((~missing).sum()), 'prior_rows_missing_NaN': int(missing.sum()),
            'scientific_parameters_match_control': True, 'audit_passed': True}


def write_results(output, audit, rows):
    comparison = {'expected_cases': 4, 'cases': rows, 'criteria': CRITERIA,
                  'strict_successes': sum(r.get('portable', {}).get('strict_success', False) and r.get('audit_passed', False) for r in rows),
                  'execution_completed': sum(r['execution_status'] == 'completed' for r in rows),
                  'execution_failed': sum(r['classification'] == 'execution_failure' for r in rows),
                  'scientific_not_solved': sum(r['classification'] == 'scientific_not_solved' for r in rows),
                  'audit_failed': sum(r['classification'] == 'audit_failure' for r in rows),
                  'primary_endpoint': 'Fixed last round; Buccaneer 3/2/2/2/2, REFMAC 5x10; never minimum Rfree',
                  'method': 'PhaseBuild seed plus external Buccaneer and REFMAC; not native PhaseBuild',
                  'scientific_consistency_passed': audit['scientific_consistency_passed'],
                  'manifest_display_metadata_warnings': audit['metadata_warnings'],
                  'public_redistribution_ready': False}
    text = ['# PhaseBuild 独立运行包四例回归审计', '',
            '固定终点严格成功：**%d/4**；完成执行 %d/4；执行失败 %d；科学未达标 %d；审计失败 %d。' % (
                comparison['strict_successes'], comparison['execution_completed'], comparison['execution_failed'],
                comparison['scientific_not_solved'], comparison['audit_failed']), '',
            '每例从预先固定的 PhaseBuild 种子开始，Buccaneer 为 3/2/2/2/2 循环，每轮 REFMAC 10 周期，共 50 周期。'
            '所有末轮结果保留在分母中；未按最低 Rfree 或参考结构选择模型。', '',
            '| 案例 | 执行/审计 | Rwork | Rfree | 原混合 Rfree | ΔRfree（百分点） | 严格主链覆盖 | 精确率 | 序列与主链联合覆盖 | 达标 |',
            '|---|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for row in rows:
        m = row.get('portable')
        if m:
            text.append('| %s | %s | %.2f%% | %.2f%% | %.2f%% | %+.2f | %.2f%% | %.2f%% | %.2f%% | %s |' % (
                row['id'], '通过' if row.get('audit_passed') else '审计失败',
                m['r_work'] * 100, m['r_free'] * 100, row['previous_fixed_HL']['r_free'] * 100,
                row['r_free_delta'] * 100, m['n_ca_c_completeness'] * 100, m['n_ca_c_precision'] * 100,
                m['joint_sequence_backbone_coverage'] * 100, '是' if m['strict_success'] else '否'))
        else:
            text.append('| %s | %s | — | — | — | — | — | — | — | 未证实 |' % (row['id'], row['classification']))
    text.extend(['', '严格标准：Rfree ≤35%，N/CA/C 同一对称像在 1 Å 内正确覆盖 ≥80%、精确率 ≥90%，序列与主链联合正确覆盖 ≥80%。', '',
                 '## 可移植性与科学协议检查', '',
                 '建模子进程的 sandbox-exec 配置拒绝读取整个 `/Applications/ccp4-9`；启动环境仅有系统 PATH 和 LANG。'
                 '包内 Python、Buccaneer、REFMAC、动态库、训练数据和字典按运行前清单逐文件核验。'
                 '逐轮直接读取 MTZ，核验原 native 列与行顺序、原自由集、固定原始 HL 精确 HKL 左连接及缺失高分辨率 NaN。'
                 '每轮建模仅使用本例上一轮输出地图和模型。REFMAC 关键词与原固定 HL 对照逐字一致。', '',
                 '若恢复运行，首轮已完成的 3 个 Buccaneer 循环按文件原样复制，仅恢复首轮失败的 REF10，'
                 '随后完成 4×(Buccaneer2+REF10)。旧命令、XML、日志、模型、失败记录和旧包/驱动快照均保留并按哈希核对。'
                 '输入 native 的 F/SIGF 和自由集保持不变；REFMAC 输出 F/SIGF 可有正常缩放，输出自由集仍须按精确 HKL 不变。', '',
                 '本包包含外部 Buccaneer 和 REFMAC，其贡献保留署名。本结果说明本机在拒绝访问原 CCP4 安装目录时可运行，'
                 '不能替代全新机器/不同系统的安装验证。公开再分发许可和对应源码条件尚未完成，因此这是本地评估包。', '',
                 '四例是已知成功开发集，不代表 100 套盲测成功率。重原子及上游组成已给定；上游定相使用过自由反射，'
                 '因此 Rfree 不是独立的端到端盲测指标。目标参考仅在构建完成后用于评价；建模未读取目标参考。', '',
                 '初次预运行中尚未开始任何建模的作用域/启动错误属于执行故障，保留原记录，未计作科学失败。'
                 '本表只计当前预先声明的四例回归；当前执行失败仍占四例分母。', '',
                 '运行时 app 清单含 %d 项已知 bytes 显示长度旧值，实际源码 SHA 全部正确。'
                 '这与科学一致性分开记录；运行原清单按原始 SHA 留存在 frozen_package_metadata。'
                 '交付前仅修正这些长度字段，补充审计见 delivery_metadata_audit.json（若已生成）。' % len(audit['metadata_warnings']), '',
                 '审计问题：' + ('无。' if not audit['issues'] else '见 audit.json；存在问题时不得将结果宣称为通过。'), ''])
    # Do not leave half-written JSON on interruption. Existing raw outputs are never changed.
    for filename, value in [('audit.json', audit), ('comparison.json', comparison)]:
        temporary = output / (filename + '.tmp')
        temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
        temporary.replace(output / filename)
    temporary = output / 'REPORT.md.tmp'
    temporary.write_text('\n'.join(text))
    temporary.replace(output / 'REPORT.md')


def finalize_delivery_metadata(output, bundle):
    """Correct four display lengths after the completed, preserved run audit.

    No scientific source, model, input, binary, library or data bytes change.
    The original run manifests and scientific audit remain unchanged.
    """
    audit = read(output / 'audit.json')
    require(audit['status'] == 'passed' and audit['scientific_consistency_passed']
            and all(row['execution_status'] == 'completed' and row['audit_passed'] for row in audit['cases']),
            'Delivery metadata may be finalized only after all four completed scientific audits pass')
    manifest = read(output / 'manifest.json')
    frozen = freeze_manifests(output, bundle, manifest)
    app_path = bundle / 'app/app_manifest.json'
    original_record = next(x for x in manifest['manifests'] if Path(x['path']).parent.name == 'app')
    original = read(frozen[original_record['path']]['path'])
    delivered = json.loads(json.dumps(original))
    changes = []
    for row in delivered['files']:
        rel = row['relative_path']
        path = bundle / 'app' / rel
        require(sha256(path) == row['sha256'], 'App file changed before metadata finalization: ' + rel)
        if rel in KNOWN_APP_SIZE_METADATA:
            old, actual = KNOWN_APP_SIZE_METADATA[rel]
            require(row['bytes'] == old and path.stat().st_size == actual, 'Unknown size correction requested')
            changes.append({'relative_path': rel, 'old_bytes': old, 'new_bytes': actual,
                            'sha256_before_and_after': row['sha256']})
            row['bytes'] = actual
        else:
            require(path.stat().st_size == row['bytes'], 'Another app size mismatch discovered')
    require({x['relative_path'] for x in changes} == set(KNOWN_APP_SIZE_METADATA)
            and only_known_size_metadata_corrections(original, delivered), 'Metadata correction exceeds declared four fields')
    current = read(app_path)
    require(current == original or current == delivered, 'Current app manifest differs beyond allowed display metadata')
    for item in manifest['manifests']:
        if item['path'] != str(app_path):
            checked_fingerprint(item)
    already_corrected = current == delivered
    if not already_corrected:
        temporary = app_path.with_name('app_manifest.json.metadata.tmp')
        require(not temporary.exists(), 'Unexpected unfinished metadata edit exists')
        with temporary.open('x') as stream:
            stream.write(json.dumps(delivered, indent=2, ensure_ascii=False) + '\n')
        temporary.replace(app_path)
    require(read(app_path) == delivered, 'Corrected app manifest did not persist exactly')
    # Re-read every member after correction; neither source SHA nor inventory
    # membership may have changed. This is a metadata supplement, not modelling.
    for row in read(app_path)['files']:
        path = bundle / 'app' / row['relative_path']
        require(sha256(path) == row['sha256'] and path.stat().st_size == row['bytes'],
                'Delivered app inventory verification failed: ' + row['relative_path'])
    supplement = {'status': 'passed', 'updated_unix': time.time(),
                  'scope': 'Delivery-only correction of four manifest display byte counts; no scientific rerun',
                  'scientific_audit': fingerprint(output / 'audit.json'),
                  'run_manifest': fingerprint(output / 'manifest.json'),
                  'frozen_package_metadata': list(frozen.values()),
                  'frozen_app_manifest': frozen[original_record['path']],
                  'delivered_app_manifest': fingerprint(app_path), 'changes': changes,
                  'only_four_bytes_fields_changed': True,
                  'all_other_manifest_fields_and_file_hashes_unchanged': True,
                  'app_files_checked': len(delivered['files']),
                  'original_run_manifest_and_scientific_audit_modified': False,
                  'scientific_outputs_sources_binaries_and_data_modified': False,
                  'model_construction_or_refinement_rerun': False,
                  'correction_was_already_present': already_corrected}
    path = output / 'delivery_metadata_audit.json'
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(supplement, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)
    return supplement


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'runs/portable_validation_20261019_retry3')
    parser.add_argument('--bundle', type=Path, default=ROOT / 'dist/local_macos_x86_64')
    parser.add_argument('--finalize-delivery-metadata', action='store_true',
                        help='After all completed audits pass, correct only four known app-manifest bytes display fields')
    args = parser.parse_args(argv)
    output, bundle = args.output.resolve(), args.bundle.resolve()
    states, pending = {}, []
    for cid in IDS:
        path = output / cid / 'state.json'
        if not path.is_file():
            pending.append(cid + ': missing state')
            continue
        states[cid] = read(path)
        if states[cid].get('status') not in TERMINAL:
            pending.append(cid + ': ' + str(states[cid].get('status')))
    if pending:
        print(json.dumps({'status': 'pending', 'no_results_written': True, 'cases': pending}, ensure_ascii=False))
        return 2
    # Lazy numerical/scientific imports: a pending run does not need its engine.
    sys.path.insert(0, str(bundle / 'app/src'))
    import gemmi
    import numpy as np
    from phasebuild.phase_refinement import parse_statistics
    from phasebuild.refinement import verify_free_flags
    from phasebuild.evaluate import evaluate
    issues, rows, package_checks = [], [], []
    launch_metadata = None
    manifest = read(output / 'manifest.json')
    require(manifest['cases'] == list(IDS), 'Declared denominator differs from four fixed cases')
    try:
        frozen = freeze_manifests(output, bundle, manifest)
        launch_metadata = freeze_launch_metadata(output, bundle)
        package_checks = inventory(bundle, manifest, frozen)
        checked_fingerprint(manifest['driver'])
        checked_fingerprint(manifest['doctor'])
        doctor = read(manifest['doctor']['path'])
        require(doctor['status'] == 'success' and all(x['passed'] for x in doctor['verification'].values()),
                'Pre-run portable runtime verification did not pass')
    except Exception as error:
        issues.append({'scope': 'package_or_driver', 'error': str(error), 'traceback': traceback.format_exc()})
    cases = {c['id']: c for c in read(ROOT / 'runs/buccaneer_sad100_20261019/manifest.json')['cases']}
    for cid in IDS:
        state = states[cid]
        if state['status'] != 'completed':
            failure = {'id': cid, 'execution_status': state['status'], 'classification': 'execution_failure',
                       'error': state.get('error'), 'strict_success': False, 'audit_passed': False,
                       'completed_refinement_stages': len(read(output / cid / 'pipeline/hybrid_report.json').get('stages', []))
                       if (output / cid / 'pipeline/hybrid_report.json').exists() else 0}
            if manifest.get('execution_recovery_from'):
                try:
                    failure['execution_recovery'] = recovery_chain(cid, output, bundle, manifest, state,
                        read(output / cid / 'pipeline/hybrid_report.json'))
                except Exception as error:
                    failure['recovery_audit_error'] = str(error)
                    issues.append({'scope': cid, 'error': str(error), 'traceback': traceback.format_exc()})
            rows.append(failure)
            continue
        try:
            row = audit_case(cid, output, bundle, cases[cid], state, manifest, np, gemmi,
                             parse_statistics, verify_free_flags, evaluate)
            rows.append(row)
        except Exception as error:
            issues.append({'scope': cid, 'error': str(error), 'traceback': traceback.format_exc()})
            rows.append({'id': cid, 'execution_status': state['status'], 'classification': 'audit_failure',
                         'error': str(error), 'strict_success': False, 'audit_passed': False})
    if any(x['scope'] == 'package_or_driver' for x in issues):
        for row in rows:
            if row.get('audit_passed'):
                row.update(audit_passed=False, classification='audit_failure')
    audit = {'updated_unix': time.time(), 'status': 'passed' if not issues else 'failed',
             'scientific_consistency_passed': not issues and all(row['execution_status'] == 'completed'
                                                               and row['audit_passed'] for row in rows),
             'metadata_warnings': [warning for package in package_checks for warning in package['metadata_warnings']],
             'expected_cases': 4, 'terminal_cases': 4, 'issues': issues,
             'package_checks': package_checks, 'audit_script': fingerprint(__file__),
             'launch_metadata_snapshot': launch_metadata,
             'run_manifest': fingerprint(output / 'manifest.json'), 'cases': rows,
             'writes': ['audit.json', 'comparison.json', 'REPORT.md', 'frozen_package_metadata/*'],
             'no_structural_experiment_launched': True, 'package_and_raw_outputs_modified': False,
             'reference_policy': 'Only posterior evaluation after fixed endpoint; no target reference in construction',
             'scope': 'Four known development cases, supplied heavy atoms/upstream composition; upstream phases used free reflections'}
    write_results(output, audit, rows)
    if args.finalize_delivery_metadata and not issues:
        finalize_delivery_metadata(output, bundle)
    print(json.dumps({'status': audit['status'], 'issues': len(issues), 'output': str(output),
                      'scientific_consistency_passed': audit['scientific_consistency_passed'],
                      'manifest_display_metadata_warnings': len(audit['metadata_warnings']),
                      'strict_successes': sum(r.get('portable', {}).get('strict_success', False)
                                              and r.get('audit_passed', False) for r in rows)}, ensure_ascii=False))
    return 1 if issues else 0


if __name__ == '__main__':
    raise SystemExit(main())
