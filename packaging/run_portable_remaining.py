"""Frozen portable hybrid-HL continuation on the 96 remaining SAD cases.

The four completed development cases are never scheduled again. References
are accessed only by the post-build evaluator, after the declared endpoint.
"""
import argparse
import concurrent.futures
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
import traceback

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if (_HERE.parent / 'engine').is_dir() else _HERE.parents[1]
BUNDLE = ROOT / 'dist/local_macos_x86_64'
PREVIOUS = ROOT / 'runs/buccaneer_sad100_20261019'
NATIVE = ROOT / 'runs/sad100_20261019'
FOUR = ROOT / 'runs/portable_validation_20261019_retry3'
EXCLUDED = ('5hhk', '7bll', '7tem', '7twc')
PROFILE = '(version 1)(allow default)(deny file-read* (subpath "/Applications/ccp4-9"))'
STOP = threading.Event()
ACTIVE = {}
ACTIVE_LOCK = threading.Lock()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4194304), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fp(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': sha(path)}


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.%d.%d.tmp' % (os.getpid(), threading.get_ident()))
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def copy_checked(source, target, expected):
    if sha(source) != expected:
        raise ValueError('Source changed: ' + str(source))
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if sha(target) != expected:
        raise ValueError('Copy differs: ' + str(target))


def prepare(out, workers, doctor):
    if out.exists():
        raise ValueError('Output already exists; never overwrite: ' + str(out))
    verification = read(doctor)
    if verification['status'] != 'success' or not all(v['passed'] for v in verification['verification'].values()):
        raise ValueError('Portable doctor must pass before freezing')
    for kind, directory, name in [('engines', 'runtime', 'runtime_manifest.json'),
                                  ('python', 'python', 'python_manifest.json'),
                                  ('app', 'app', 'app_manifest.json')]:
        if verification['verification'][kind]['manifest_sha256'] != sha(BUNDLE / directory / name):
            raise ValueError('Doctor is stale for ' + kind)
    recipe = read(ROOT / 'research/remaining96_seed_recipe_20261020.json')
    previous = read(PREVIOUS / 'manifest.json')
    native = {c['id']: c for c in read(NATIVE / 'manifest.json')['cases']}
    out.mkdir(parents=True)
    for directory, name in [('runtime', 'runtime_manifest.json'), ('python', 'python_manifest.json'),
                             ('app', 'app_manifest.json')]:
        copy_checked(BUNDLE / directory / name, out / 'frozen_package_metadata' / directory / name,
                     sha(BUNDLE / directory / name))
    shutil.copytree(BUNDLE / 'app', out / 'engine',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '.gui', 'gui_runs', 'gui_inbox'))
    for row in read(BUNDLE / 'app/app_manifest.json')['files']:
        if sha(out / 'engine' / row['relative_path']) != row['sha256']:
            raise ValueError('Frozen pipeline engine differs')
    for relative, item in recipe['sources'].items():
        copy_checked(item['path'], out / 'seed_engine' / relative, item['sha256'])
    for name, item in recipe['assets'].items():
        copy_checked(item['path'], out / 'seed_engine/assets' / (name + '.npz'), item['sha256'])
    save(out / 'seed_recipe.json', recipe)
    # Freeze orchestration separately from both scientific engines.
    for name in ['run_portable_remaining.py', 'summarize_portable_remaining.py', 'audit_portable_remaining.py']:
        shutil.copy2(ROOT / 'packaging' / name, out / name)
    cases = []
    for original in previous['cases']:
        if original['id'] in EXCLUDED:
            continue
        case = dict(original)
        if case['supported']:
            old = native[case['id']]
            case['seed_map'] = old['map']
            case['seed_map_sha256'] = old['prepared_sha256'][str(Path(old['map']).relative_to(NATIVE))]
            specs = [('initial_mtz', 'initial_sha256', 'initial.mtz'),
                     ('phase_source', 'phase_sha256', 'original_hl.mtz'),
                     ('sequence', 'sequence_sha256', 'sequence.fasta'),
                     ('seed_map', 'seed_map_sha256', 'seed_map.mtz')]
            if case['has_free']:
                specs.append(('native', 'native_sha256', 'native.mtz'))
            case['original_inputs'] = {}
            for key, digest, filename in specs:
                source = Path(case[key]); target = out / 'inputs' / case['id'] / filename
                copy_checked(source, target, case[digest])
                case['original_inputs'][key] = {'path': str(source), 'sha256': case[digest]}
                case[key] = str(target)
        cases.append(case)
    if len(cases) != 96 or sum(c['supported'] for c in cases) != 86:
        raise ValueError('Unexpected cohort')
    rules = dict(previous['rules'])
    rules.update(phase_schedule='Every refinement uses the same original experimental HL; subsequent build maps/models are own previous output',
                 seed_recipe='512 contextual seeds -> bounded weak trace -> original overlap assembly -> fixed 0.05 discriminant pruning',
                 seed_capacity='int(1.25*sequence_residues), identical to four-case prefix',
                 endpoint='Fixed round5 with50REF cycles; no minimumRfree selection',
                 seed_timeout_seconds=[21600,21600,7200,1800], pipeline_timeout_seconds=172800,
                 missing_free_policy='One PhaseBuild-seeded3cycleBuccaneer construction, no refinement/Rfree',
                 workers=workers)
    manifests = [fp(BUNDLE / d / n) for d, n in [('runtime','runtime_manifest.json'),
                 ('python','python_manifest.json'),('app','app_manifest.json')]]
    engines = {str(p.relative_to(out)): sha(p) for d in ['engine','seed_engine']
               for p in sorted((out/d).rglob('*')) if p.is_file()}
    orchestration = {name: sha(out/name) for name in
                     ['run_portable_remaining.py','summarize_portable_remaining.py','audit_portable_remaining.py']}
    save(out / 'manifest.json', {'created': time.time(), 'cases': cases,
         'excluded_development_cases': list(EXCLUDED), 'previous_four_run': str(FOUR),
         'previous_buccaneer_run': str(PREVIOUS), 'previous_native_run': str(NATIVE),
         'original_corpus_manifest': str(PREVIOUS/'manifest.json'),
         'buccaneer_baseline': str(PREVIOUS), 'native_baseline': str(NATIVE),
         'reference_evaluation_root': str(PREVIOUS/'evaluation'),
         'rules': rules, 'workers': workers, 'bundle': str(BUNDLE), 'bundle_manifests': manifests,
         'doctor': fp(doctor), 'engine_sha256': engines, 'orchestration_sha256': orchestration,
         'sandbox_profile': PROFILE, 'public_redistribution_ready': False,
         'scope': 'Follow-up96, four development cases separately retained; no new phasing, reference/model-based selection or reference input'})
    for case in cases:
        save(out / 'states' / (case['id']+'.json'), {'id': case['id'], 'sample_class': case['sample_class'],
             'status': 'queued' if case['supported'] else 'unsupported_nonprotein',
             'completed_rounds': 0, 'current': '', 'primary_metric': None})
    save(out / 'batch_process.json', {'status': 'prepared', 'workers': workers, 'prepared': time.time()})


def load_helper(out, name):
    spec = importlib.util.spec_from_file_location(name, out / (name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def environment():
    result = {'PATH': '/usr/bin:/bin', 'LANG': 'C', 'PYTHONUTF8': '1',
              'PYTHONDONTWRITEBYTECODE':'1', 'PYTHONNOUSERSITE':'1',
              'PHASEBUILD_STANDALONE':'1', 'PHASEBUILD_RUNTIME':str(BUNDLE/'runtime')}
    result.update({name: '1' for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
                  'MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']})
    return result


def terminate(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=4)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=3)
    except ProcessLookupError:
        pass


def execute(command, log_path, timeout, state, state_path, cwd):
    if STOP.is_set():
        raise InterruptedError('Batch stop requested')
    wrapped = ['/usr/bin/sandbox-exec','-p',PROFILE] + [str(v) for v in command]
    state['command'] = wrapped
    save(state_path, state)
    log_path.parent.mkdir(parents=True,exist_ok=True)
    with log_path.open('x',encoding='utf-8') as log:
        process = subprocess.Popen(wrapped,cwd=str(cwd),env=environment(),stdout=log,
                    stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
        with ACTIVE_LOCK:
            ACTIVE[state['id']] = process
        state['child_pid'] = process.pid
        save(state_path, state)
        try:
            code = process.wait(timeout=timeout)
            if code:
                raise RuntimeError('Program exit%d; inspect%s' % (code,log_path))
        except BaseException:
            terminate(process)
            raise
        finally:
            with ACTIVE_LOCK:
                ACTIVE.pop(state['id'],None)
            state.pop('child_pid',None)


def no_free_build(case, state, out, state_path):
    import gemmi
    sys.path.insert(0,str(out/'engine/src'))
    from phasebuild.runtime import runtime_environment
    folder = out/'cases'/case['id']/'pipeline/round_01_build'
    folder.mkdir(parents=True,exist_ok=False)
    runtime = BUNDLE/'runtime'; library=runtime/'lib/data/reference_structures/reference-1tqw'
    command=[str(runtime/'bin/cbuccaneer'),'-mtzin',case['initial_mtz'],
        '-mtzin-ref',str(library.with_suffix('.mtz')),'-pdbin-ref',str(library.with_suffix('.pdb')),
        '-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
        '-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
        '-colin-fo','/*/*/[F,SIGF]','-colin-phifom','/*/*/[PHWT,FOM]',
        '-colin-fc','/*/*/[FWT,PHWT]','-pdbin',state['seed_model'],'-model-filter',
        '-seqin',case['sequence'],'-pdbout',str(folder/'model.pdb'),'-xmlout',str(folder/'buccaneer.xml')]
    command+=['-'+s for s in ['find','grow','join','link','sequence','correct','filter','ncsbuild','prune','rebuild']]
    command+=['-cycles','3','-jobs','1','-no-correlation-mode']
    if case['build_semet']:
        command+=['-build-semet']
    save(folder/'command.json',command)
    # Direct binary invocation needs the same self-contained runtime environment.
    launcher=folder/'execute_build.py'
    launcher.write_text('import os,sys\nsys.path.insert(0,'+repr(str(out/'engine/src'))+')\n'
        'from phasebuild.runtime import runtime_environment\n'
        'os.execve('+repr(command[0])+','+repr(command)+',runtime_environment('+repr(str(runtime))+',scratch='+repr(str(folder/'scratch'))+'))\n',encoding='utf-8')
    execute([BUNDLE/'python/bin/python3','-u',launcher],folder/'buccaneer.log',21600,state,state_path,folder)
    model=gemmi.read_structure(str(folder/'model.pdb'))
    model.make_mmcif_document().write_file(str(folder/'model.cif'))
    save(folder/'build_report.json',{'status':'success','engine':'external-Buccaneer-seeded-single-build',
         'command':command,'cycles':3,'model':fp(folder/'model.cif'),'no_refinement':True,'r_free':None})
    state.update(status='completed_no_free',completed_rounds=1,primary_metric=None,
                 primary_model=str(folder/'model.cif'),build_directory=str(folder),build_report=str(folder/'build_report.json'))


def run_case(case, out):
    path=out/'states'/(case['id']+'.json');state=read(path)
    if state['status']!='queued':
        return state
    if STOP.is_set():
        return state
    caseout=out/'cases'/case['id']
    state.update(status='running',started=time.time(),pid=os.getpid(),current='preflight')
    save(path,state)
    try:
        caseout.mkdir(parents=True,exist_ok=False)
        for key,digest in [('initial_mtz','initial_sha256'),('phase_source','phase_sha256'),
                           ('sequence','sequence_sha256'),('seed_map','seed_map_sha256')]:
            if sha(case[key])!=case[digest]:
                raise ValueError('Frozen input changed:'+key)
        if case['has_free'] and sha(case['native'])!=case['native_sha256']:
            raise ValueError('Frozen native changed')
        py=BUNDLE/'python/bin/python3';engine=out/'seed_engine';profile=engine/'assets/backbone_profile.npz'
        capacity=int(1.25*case['sequence_residues'])
        directories=[caseout/name for name in ['01_context','02_weak','03_join','04_seed']]
        commands=[
            [py,'-u',engine/'run_context_seeds.py',case['seed_map'],'--profile',profile,'--count','512','--output',directories[0]],
            [py,'-u',engine/'run_weak_trace.py',case['seed_map'],'--profile',profile,'--candidates',directories[0]/'candidates.npz',
             '--capacity',str(capacity),'--max-length','30','--output',directories[1]],
            [py,'-u',engine/'run_fragment_join.py',case['seed_map'],'--profile',profile,'--candidate-archive',directories[1]/'candidates.npz',
             '--capacity',str(capacity),'--maximum-length','40','--beam-width','3','--output',directories[2]],
            [py,'-u',engine/'run_discriminant_prune.py',case['seed_map'],'--model',directories[2]/'backbone.cif',
             '--discriminant',engine/'assets/discriminant.npz','--calibration-decoy-rate','0.05','--minimum-length','4','--bridge','1','--output',directories[3]]]
        for number,(command,timeout,directory) in enumerate(zip(commands,[21600,21600,7200,1800],directories),1):
            state.update(current='seed_stage_%d'%number);save(path,state)
            execute(command,caseout/('seed_stage_%d.log'%number),timeout,state,path,engine)
            save(directory/'stage_completion.json',{'status':'completed','command':[str(v) for v in command],
                 'report':fp(directory/'build_report.json'),'finished':time.time()})
        state.update(seed_model=str(directories[-1]/'backbone.cif'),seed_model_fingerprint=fp(directories[-1]/'backbone.cif'),current='hybrid_hl')
        save(path,state)
        if not case['has_free']:
            no_free_build(case,state,out,path)
        else:
            command=[py,'-u',out/'engine/run_phasebuild_hybrid_hl.py','--map',case['initial_mtz'],
                '--native',case['native'],'--model',state['seed_model'],'--sequence',case['sequence'],
                '--phase-prior',case['phase_source'],'--phase-prior-sha256',case['phase_sha256'],
                '--runtime',BUNDLE/'runtime','--rounds','5','--cycles','10','--output',caseout/'pipeline']
            if case['build_semet']:
                command+=['--build-semet']
            execute(command,caseout/'pipeline.log',172800,state,path,out/'engine')
            report=read(caseout/'pipeline/hybrid_report.json')
            if report['status']!='completed' or len(report['stages'])!=5:
                raise RuntimeError('Full fixed endpoint missing')
            state.update(status='completed',completed_rounds=5,primary_model=str(caseout/'pipeline/primary_model.cif'),
                         pipeline_report=str(caseout/'pipeline/hybrid_report.json'),primary_metric=report['primary_metric'])
    except BaseException as error:
        state.update(status='cancelled' if STOP.is_set() else 'timeout' if isinstance(error,subprocess.TimeoutExpired) else 'seed_failed' if state['current'].startswith('seed_stage_') else 'failed',
                     failure_stage=state['current'],error=str(error),traceback=traceback.format_exc())
    state.update(current='',finished=time.time(),seconds=time.time()-state['started']);state.pop('child_pid',None)
    save(path,state)
    return state


def postprocess(case, state, out, evaluator, auditor):
    if state['status'] not in ('completed','completed_no_free'):
        return state
    try:
        state['endpoint_audit']=auditor.audit_completed(case,state,out)
        if not state['endpoint_audit'].get('passed',False):
            raise RuntimeError('Endpoint audit failed: '+str(state['endpoint_audit'].get('issues')))
        state=evaluator.evaluate_completed(case,state,out)
        state.pop('evaluation_error',None);state.pop('evaluation_traceback',None)
    except Exception as error:
        state.update(evaluation_error=str(error),evaluation_traceback=traceback.format_exc())
    save(out/'states'/(case['id']+'.json'),state)
    return state


def run(out):
    lock=(out/'batch.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=read(out/'manifest.json')
    for relative,digest in {**manifest['engine_sha256'],**manifest['orchestration_sha256']}.items():
        if sha(out/relative)!=digest:
            raise ValueError('Frozen source changed: '+relative)
    for record in manifest['bundle_manifests']:
        if sha(record['path'])!=record['sha256']:
            raise ValueError('Shared portable manifest changed')
    sys.path.insert(0,str(out/'engine/src'))
    evaluator=load_helper(out,'summarize_portable_remaining');auditor=load_helper(out,'audit_portable_remaining')
    workers=manifest['workers']; started=time.time()
    save(out/'batch_process.json',{'pid':os.getpid(),'status':'running','started':started,'workers':workers})
    def stop(signum,frame):
        STOP.set()
        with ACTIVE_LOCK:
            children=list(ACTIVE.values())
        for process in children:
            terminate(process)
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures={pool.submit(run_case,case,out):case for case in sorted(manifest['cases'],key=lambda c:c.get('sequence_residues',0)) if case['supported'] and read(out/'states'/(case['id']+'.json'))['status']=='queued'}
            while futures:
                done,_=concurrent.futures.wait(futures,timeout=30,return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    case=futures.pop(future)
                    try:
                        state=postprocess(case,future.result(),out,evaluator,auditor)
                    except Exception as error:
                        state=read(out/'states'/(case['id']+'.json'));state.update(status='failed',error=str(error),traceback=traceback.format_exc());save(out/'states'/(case['id']+'.json'),state)
                    print(case['id'],state['status'],state.get('primary_metric'),state.get('evaluation_error',''),flush=True)
                evaluator.summarize(out)
        evaluator.summarize(out)
        save(out/'batch_process.json',{'pid':os.getpid(),'status':'stopped' if STOP.is_set() else 'finished','started':started,'finished':time.time(),'workers':workers})
    except BaseException as error:
        STOP.set()
        with ACTIVE_LOCK:
            children=list(ACTIVE.values())
        for process in children:
            terminate(process)
        save(out/'batch_process.json',{'pid':os.getpid(),'status':'failed','error':str(error),'started':started,'finished':time.time()})
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','run','status','postprocess'])
    parser.add_argument('--root',type=Path,default=ROOT/'runs/portable_remaining96_20261020')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--doctor',type=Path,default=FOUR/'doctor_release.json')
    args=parser.parse_args();out=args.root.resolve()
    if not 1<=args.workers<=8:
        parser.error('workers must be1..8')
    if args.action=='prepare':
        prepare(out,args.workers,args.doctor)
    elif args.action=='run':
        run(out)
    else:
        sys.path.insert(0,str(out/'engine/src'))
        evaluator=load_helper(out,'summarize_portable_remaining')
        if args.action=='postprocess':
            auditor=load_helper(out,'audit_portable_remaining')
            for case in read(out/'manifest.json')['cases']:
                state=read(out/'states'/(case['id']+'.json'))
                if state['status'] in ('completed','completed_no_free') and (not state.get('evaluation') or state.get('evaluation_error')):
                    postprocess(case,state,out,evaluator,auditor)
        summary=evaluator.summarize(out)
        print(json.dumps({k:v for k,v in summary.items() if k!='cases'},indent=2,ensure_ascii=True))


if __name__=='__main__':
    main()
