"""Fixed four-case packaging regression; references are only read after building."""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import signal
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / 'dist/local_macos_x86_64'
OUT = Path(os.environ.get('PHASEBUILD_VALIDATION_OUTPUT', str(ROOT / 'runs/portable_validation_20261019'))).resolve()
DOCTOR = Path(os.environ.get('PHASEBUILD_DOCTOR_REPORT', str(ROOT / 'runs/standalone_feasibility_20261019/doctor_denied_ccp4.json'))).resolve()
RESUME_ROOT = Path(os.environ['PHASEBUILD_RESUME_ROOT']).resolve() if os.environ.get('PHASEBUILD_RESUME_ROOT') else None
IDS = ['5hhk', '7bll', '7tem', '7twc']
PROFILE = '(version 1)(allow default)(deny file-read* (subpath "/Applications/ccp4-9"))'
sys.path.insert(0, str(BUNDLE / 'app/src'))


def read(path):
    return json.loads(Path(path).read_text())


def fingerprint(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


def run(case):
    from phasebuild.evaluate import evaluate
    cid = case['id']
    folder = OUT / cid
    folder.mkdir(exist_ok=False)
    prior = ROOT / 'runs/sad_research_v2_20261019/full_hybrid_hl_recovered' / cid / 'pipeline/hybrid_report.json'
    previous = read(prior)
    seed = previous['inputs']['phasebuild_seed_model']['path']
    inputs = {'map': fingerprint(case['initial_mtz']), 'native': fingerprint(case['native']),
              'sequence': fingerprint(case['sequence']), 'model': fingerprint(seed),
              'phase_prior': fingerprint(case['phase_source'])}
    assert inputs['map']['sha256'] == case['initial_sha256']
    assert inputs['native']['sha256'] == case['native_sha256']
    assert inputs['sequence']['sha256'] == case['sequence_sha256']
    assert inputs['phase_prior']['sha256'] == case['phase_sha256']
    assert inputs['model'] == previous['inputs']['phasebuild_seed_model']
    command = ['/usr/bin/sandbox-exec', '-p', PROFILE, str(BUNDLE / 'phasebuild'), 'hybrid-hl',
               '--map', inputs['map']['path'], '--native', inputs['native']['path'],
               '--sequence', inputs['sequence']['path'], '--model', seed,
               '--phase-prior', inputs['phase_prior']['path'],
               '--phase-prior-sha256', inputs['phase_prior']['sha256'],
               '--rounds', '5', '--cycles', '10', '--output', str(folder / 'pipeline')]
    if case['build_semet']:
        command += ['--build-semet']
    if RESUME_ROOT is not None:
        source = RESUME_ROOT / cid / 'pipeline/round_01_build'
        command += ['--existing-first-build', str(source)]
        history = {'previous_failed_output': str(RESUME_ROOT / cid),
                   'previous_failed_state': fingerprint(RESUME_ROOT / cid / 'state.json'),
                   'previous_failed_pipeline_report': fingerprint(RESUME_ROOT / cid / 'pipeline/hybrid_report.json'),
                   'previous_failed_refinement_log': fingerprint(RESUME_ROOT / cid / 'pipeline/round_01_refine/refmac.log'),
                   'previous_failed_refinement_report': fingerprint(RESUME_ROOT / cid / 'pipeline/round_01_refine/refinement_report.json'),
                   'first_construction_model': fingerprint(source / 'model.cif'),
                   'first_construction_command': fingerprint(source / 'command.json'),
                   'old_runtime_manifest': fingerprint(RESUME_ROOT / 'before_include_repair_runtime_manifest.json'),
                   'old_app_manifest': fingerprint(RESUME_ROOT / 'failed_package_metadata/app/app_manifest.json'),
                   'old_driver': fingerprint(RESUME_ROOT / 'failed_app_snapshot/validation_driver.py'),
                   'execution_repair': 'Add original include/environ.def/default.def and CINCL locator; reuse completed3cyclefirstconstruction; unchanged science protocol'}
        other_failed = ROOT / 'runs/portable_validation_20261019_retry2'
        if other_failed.exists() and other_failed != RESUME_ROOT:
            history['additional_failed_execution_attempts'] = [{
                'directory': str(other_failed / cid),
                'state': fingerprint(other_failed / cid / 'state.json'),
                'pipeline_report': fingerprint(other_failed / cid / 'pipeline/hybrid_report.json'),
                'refinement_report': fingerprint(other_failed / cid / 'pipeline/round_01_refine/refinement_report.json'),
                'refinement_log': fingerprint(other_failed / cid / 'pipeline/round_01_refine/refmac.log'),
                'app_manifest_snapshot': fingerprint(other_failed / 'failed_app_snapshot/app_manifest.json'),
                'driver_snapshot': fingerprint(other_failed / 'validation_driver.py'),
                'repair': 'Set independent per-job CCP4_SCR; previous execution stopped before cycle0'}]
        history['execution_repair'] += '; set per-job writable CCP4_SCR without inherited installation settings'
        save(folder / 'failure_history.json', history)
    state = {'id': cid, 'status': 'running', 'started': time.time(), 'pid': os.getpid(),
             'inputs': inputs, 'command': command, 'original_ccp4_read_denied': True,
             'protocol': 'Same initial PhaseBuild coordinates and fixed originalHL hybrid5x10; no scientific retuning'}
    save(folder / 'state.json', state)
    try:
        with (folder / 'runner.log').open('w') as log:
            process = subprocess.Popen(command, cwd=BUNDLE, env={'PATH':'/usr/bin:/bin','LANG':'C'},
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            state['child_pid'] = process.pid
            save(folder / 'state.json', state)
            try:
                code = process.wait(timeout=21600)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                try: process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL); process.wait(timeout=3)
                raise RuntimeError('Portable case exceeded6hour execution timeout')
            if code:
                raise RuntimeError('Portable pipeline exit code ' + str(code))
        report = read(folder / 'pipeline/hybrid_report.json')
        assert report['status'] == 'completed' and len(report['stages']) == 5
        # The target reference is first opened here, after a completed fixed50 output.
        reference = ROOT / 'runs/buccaneer_sad100_20261019/evaluation' / cid / 'protein_reference.cif'
        evaluation = evaluate(folder / 'pipeline/primary_model.cif', reference)
        save(folder / 'evaluation.json', evaluation)
        joint = evaluation['sequence_evaluation']['correct_type_and_backbone'] / evaluation['reference_residues']
        metric = report['primary_metric']
        passed = (metric['r_free'] <= .35 and evaluation['n_ca_c_completeness'] >= .8
                  and evaluation['n_ca_c_precision'] >= .9 and joint >= .8)
        state.update(status='completed', primary_metric=metric,
                     built=evaluation['built_residues'], correct=evaluation['n_ca_c_correct_residues'],
                     coverage=evaluation['n_ca_c_completeness'], precision=evaluation['n_ca_c_precision'],
                     joint_coverage=joint, strict_success=bool(passed),
                     rfree_delta_vs_previous=metric['r_free'] - previous['primary_metric']['r_free'])
    except BaseException as error:
        state.update(status='failed', error=str(error), traceback=traceback.format_exc())
    state.pop('child_pid', None)
    state['finished'] = time.time()
    save(folder / 'state.json', state)
    print(cid, state['status'], state.get('primary_metric'), state.get('error'), flush=True)
    return state


def main():
    OUT.mkdir(exist_ok=False)
    cases = {case['id']:case for case in read(ROOT / 'runs/buccaneer_sad100_20261019/manifest.json')['cases']}
    doctor = read(DOCTOR)
    assert doctor['status'] == 'success' and all(row['passed'] for row in doctor['verification'].values())
    save(OUT / 'manifest.json', {'cases': IDS, 'sandbox_profile': PROFILE,
        'doctor': fingerprint(DOCTOR),
        'manifests': [fingerprint(BUNDLE / name) for name in
            ['runtime/runtime_manifest.json','python/python_manifest.json','app/app_manifest.json']],
        'driver': fingerprint(__file__), 'primary_endpoint': 'Fixed last round,5x10REF; not minimumRfree',
        'scope': 'Four known development cases; initial native seed reused solely as declared construction input; no target reference in modelling',
        'public_redistribution_ready': False,
        'execution_recovery_from': str(RESUME_ROOT) if RESUME_ROOT is not None else None})
    save(OUT / 'process.json', {'pid':os.getpid(),'status':'running','started':time.time()})
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(guarded_run, cases[cid]) for cid in IDS]
        for future in concurrent.futures.as_completed(futures):
            rows.append(future.result())
            save(OUT / 'summary.json', rows)
    save(OUT / 'process.json', {'pid':os.getpid(),'status':'finished','finished':time.time()})


def guarded_run(case):
    try:
        return run(case)
    except BaseException as error:
        folder = OUT / case['id']; folder.mkdir(exist_ok=True)
        state = {'id':case['id'],'status':'failed','failure_class':'execution_setup',
                 'error':str(error),'traceback':traceback.format_exc(),'finished':time.time()}
        save(folder / 'state.json', state)
        print(case['id'], 'setup_failed', str(error), flush=True)
        return state


if __name__ == '__main__':
    main()
