"""Checkpointed SAD -> DM -> model-building orchestration, no new numerical target."""
import json
import hashlib
from contextlib import contextmanager
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid
import shutil
from .data import fingerprint
from .runtime import resolve_runtime, runtime_info
from .sad_bridge import inspect_inputs, prepare_sad, adapt_maps, FLAG

DEFAULTS = dict(copies=1, phase_cycles=5, dm_cycles=5, rounds=5, cycles=10,
                free_value=0, seed=20261023, labels='', build_semet=False)


def validate(config, inspect=True):
    cfg = {**DEFAULTS, **config}
    for key,low,high in [('copies',1,100),('phase_cycles',1,100),('dm_cycles',1,100),
                         ('rounds',1,20),('cycles',1,200),('free_value',0,999),('seed',0,2147483647)]:
        value=cfg[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or int(value)!=value or not low<=value<=high:
            raise ValueError('Invalid integer parameter: '+key)
        cfg[key]=int(value)
    if not isinstance(cfg['build_semet'],bool):raise ValueError('build_semet must be boolean')
    if cfg['build_semet']:raise ValueError('Unified entry does not yet support automatic SeMet composition; use the independent workflows')
    for key in ['sad','native','sequence','heavy','reference']:
        path=Path(cfg.get(key,'')).expanduser().resolve()
        if not path.is_file():raise ValueError('Missing input file: '+key)
        cfg[key]=str(path)
    cfg['output']=str(Path(cfg['output']).expanduser().resolve())
    cfg['runtime']=str(resolve_runtime(cfg.get('runtime') or None,required=('cbuccaneer','refmac5')))
    for suffix in ['.pdb','.mtz']:
        if not (Path(cfg['runtime'])/'lib/data/reference_structures'/('reference-1tqw'+suffix)).is_file():
            raise ValueError('External construction reference unavailable')
    if not cfg.get('fpp') or not isinstance(cfg['fpp'],str):raise ValueError('Specify element=f-double-prime')
    if inspect:
        if sys.version_info<(3,12):raise ValueError('SAD 全流程需要 Python 3.12；请使用整合版 conda 环境启动 GUI')
        from oasis_py.cli import pairs
        from oasis_py.data import read_mtz, read_atoms
        from oasis_py.sequence import sequence_asu
        from oasis_py.reference_data import read_reference
        from .sequencing import read_sequences
        sad,native,flags,partition=inspect_inputs(cfg['sad'],cfg['native'],cfg['free_value'])
        data,_=read_mtz(cfg['sad'],'sad',pairs(cfg['labels']))
        heavy=read_atoms(cfg['heavy'],data)
        import gemmi
        elements={}
        for name,value in pairs(cfg['fpp']).items():
            element=gemmi.Element(name)
            number=float(value)
            if not element.atomic_number or not math.isfinite(number) or number<=0:raise ValueError('Invalid f-double-prime')
            elements[element.name]=number
        if set(heavy.elements)-set(elements):raise ValueError('Missing f-double-prime for a heavy-atom element')
        # Both engines must agree on the same explicit sequence/copy count.
        read_sequences(cfg['sequence'])
        sequence_asu(cfg['sequence'],cfg['copies'],data.cell,data.spacegroup)
        read_reference(cfg['reference'])
    return cfg


def stage_plan(cfg, python=None):
    root=Path(__file__).resolve().parents[2];out=Path(cfg['output']);python=python or sys.executable
    prepare=out/'01_inputs';phase=out/'02_sad';dm=out/'03_dm';bridge=out/'04_bridge';build=out/'05_build'
    common=['--sad',str(prepare/'sad.mtz'),'--sequence',cfg['sequence'],'--heavy',cfg['heavy'],
            '--fpp',cfg['fpp'],'--copies',str(cfg['copies']),'--labels',cfg['labels'],
            '--free-value',str(cfg['free_value']),'--seed',str(cfg['seed'])]
    phase_cmd=[python,'-u',str(root/'run_oasis_stage.py'),'phase']+common+['--cycles',str(cfg['phase_cycles']),'--output',str(phase/'experimental.mtz')]
    dm_cmd=[python,'-u',str(root/'run_oasis_stage.py'),'dm']+common+['--cycles',str(cfg['dm_cycles']),
        '--initial',str(phase/'experimental.mtz'),'--reference',cfg['reference'],'--output',str(dm/'map_phases.mtz')]
    build_cmd=[python,'-u',str(root/'run_build_from_phases.py'),'--map',str(bridge/'building.mtz'),
        '--native',cfg['native'],'--sequence',str(prepare/'building.fasta'),'--phase-prior',str(phase/'experimental.mtz'),
        '--runtime',cfg['runtime'],'--free-value',str(cfg['free_value']),'--rounds',str(cfg['rounds']),
        '--cycles',str(cfg['cycles']),'--output',str(build)]
    if cfg['build_semet']:build_cmd.append('--build-semet')
    return [dict(id='prepare',name='输入与原自由集核验',directory=str(prepare),command=None,
                 outputs=[str(prepare/'sad.mtz'),str(prepare/'sad.json'),str(prepare/'building.fasta')]),
            dict(id='sad',name='SAD 定相：保存密度修饰前实验 HL',directory=str(phase),command=phase_cmd,
                 outputs=[str(phase/'experimental.mtz'),str(phase/'experimental.json')]),
            dict(id='dm',name='密度修饰：生成建模相位',directory=str(dm),command=dm_cmd,
                 outputs=[str(dm/'map_phases.mtz'),str(dm/'map_phases.json')]),
            dict(id='bridge',name='按 HKL 生成工作集建模图',directory=str(bridge),command=None,
                 outputs=[str(bridge/'building.mtz'),str(bridge/'building.audit.json')]),
            dict(id='build',name='Buccaneer 构建与固定原 HL 精修',directory=str(build),command=build_cmd,
                 outputs=[str(build/'primary_model.cif'),str(build/'hybrid_report.json')])]


def source_identity():
    root=Path(__file__).resolve().parents[2]
    files=list((root/'src').rglob('*.py'))+[root/'run_oasis_stage.py',root/'run_sad_to_model.py',root/'run_build_from_phases.py']
    return {str(p.relative_to(root)):fingerprint(p)['sha256'] for p in sorted(files)}


def run_command(command, log, environment):
    process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
    try:
        code=process.wait()
        if code:raise RuntimeError('Stage process exited with code '+str(code))
    finally:
        if process.poll() is None:
            try:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=1.5)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait(timeout=1)
            except ProcessLookupError:pass


@contextmanager
def output_lock(out):
    import fcntl
    with (Path(out)/'.pipeline.lock').open('a+') as handle:
        try:fcntl.flock(handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Pipeline is already running in this output directory')
        try:
            handle.seek(0);handle.truncate();handle.write(str(os.getpid())+'\n');handle.flush()
            yield
        finally:fcntl.flock(handle.fileno(),fcntl.LOCK_UN)


def runtime_identity(runtime):
    runtime=Path(runtime);manifest=runtime/'runtime_manifest.json'
    if manifest.is_file():
        rows=json.loads(manifest.read_text())['files']
        for row in rows:
            if fingerprint(runtime/row['relative_path'])['sha256']!=row['sha256']:
                raise ValueError('Runtime file differs from manifest: '+row['relative_path'])
        return {'manifest':fingerprint(manifest),'verified_files':len(rows)}
    paths=list((runtime/'lib/data').rglob('*'))
    paths += list((runtime/'lib').glob('*.dylib'))+list((runtime/'lib').glob('*.so*'))
    paths += [runtime/'include/environ.def',runtime/'include/default.def']
    digest=hashlib.sha256();count=0
    for path in sorted(set(paths)):
        if not path.is_file():continue
        digest.update(str(path.relative_to(runtime)).encode()+b'\0'+fingerprint(path)['sha256'].encode()+b'\n');count+=1
    return {'resources_sha256':digest.hexdigest(),'files':count}


def run_pipeline(config, resume=False, runner=None, progress=print):
    cfg=validate(config);out=Path(cfg['output'])
    if not resume:out.mkdir(parents=True,exist_ok=False)
    elif not out.is_dir():raise ValueError('Resume output directory is missing')
    with output_lock(out):
        return _run_pipeline(cfg,resume,runner,progress)


def _run_pipeline(config, resume=False, runner=None, progress=print):
    cfg=config
    out=Path(cfg['output']);state_path=out/'pipeline_report.json';stages=stage_plan(cfg)
    inputs={key:fingerprint(cfg[key]) for key in ['sad','native','sequence','heavy','reference']}
    runtime=Path(cfg['runtime']);inputs.update({name:fingerprint(runtime/'bin'/name) for name in ['cbuccaneer','refmac5']})
    # Lock exact resources too; replacement reference/dictionaries cannot silently
    # enter a resumed experiment. The complete runtime manifest is preferred.
    manifest=runtime/'runtime_manifest.json'
    if manifest.is_file():inputs['runtime_manifest']=fingerprint(manifest)
    for suffix in ['.mtz','.pdb']:
        inputs['builder_reference'+suffix]=fingerprint(runtime/'lib/data/reference_structures'/('reference-1tqw'+suffix))
    inputs['runtime_resources']=runtime_identity(runtime)
    import gemmi,numpy,scipy
    environment={'python':sys.version,'executable':sys.executable,'gemmi':gemmi.__version__,
                 'numpy':numpy.__version__,'scipy':scipy.__version__}
    sources=source_identity()
    if resume:
        state=json.loads(state_path.read_text())
        if state['config']!=cfg or state['inputs']!=inputs or state['sources']!=sources or state.get('environment')!=environment:
            raise ValueError('Resume refused: input, protocol, runtime or source identity changed')
        if [s['id'] for s in state['completed']] != [s['id'] for s in stages[:len(state['completed'])]]:
            raise ValueError('Checkpoint stage order differs from the protocol')
        for saved in state['completed']:
            for item in saved['artifacts']:
                if fingerprint(item['path'])!=item:raise ValueError('Checkpoint artifact changed: '+item['path'])
    else:
        state=dict(engine='oasis-sad-to-model',status='running',config=cfg,inputs=inputs,sources=sources,
                   completed=[],created=time.time(),runtime=runtime_info(runtime),environment=environment,
                   scope='SAD and DM use working reflections only; original flags unchanged. Later model-derived maps are not certified work-only; Rfree is not an independent end-to-end blind score.')
    def save():
        state['updated']=time.time();temp=state_path.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');temp.replace(state_path)
    state.update(status='running',error=None);save()
    old=None
    if __import__('threading').current_thread() is __import__('threading').main_thread():
        def stop(signum,frame):raise InterruptedError('Unified SAD pipeline cancelled')
        old=signal.signal(signal.SIGTERM,stop)
    runner=runner or run_command
    try:
        for index,stage in enumerate(stages):
            if index<len(state['completed']):continue
            folder=Path(stage['directory'])
            if folder.exists():
                archived=folder.with_name(folder.name+'.interrupted-'+uuid.uuid4().hex[:8]);folder.rename(archived)
                state.setdefault('preserved_interrupted_directories',[]).append(str(archived))
            state.update(current=stage['id'],current_name=stage['name']);save();progress(stage['name'])
            if stage['id']=='prepare':
                folder.mkdir()
                prepare_sad(cfg['sad'],cfg['native'],folder/'sad.mtz',cfg['free_value'])
                from .sequencing import read_sequences
                lines=[]
                for copy in range(cfg['copies']):
                    for i,(name,seq) in enumerate(read_sequences(cfg['sequence'])):
                        lines.extend(['>copy{}_chain{}'.format(copy+1,i+1),seq])
                (folder/'building.fasta').write_text('\n'.join(lines)+'\n')
            elif stage['id']=='bridge':
                folder.mkdir();adapt_maps(cfg['native'],out/'02_sad/experimental.mtz',out/'03_dm/map_phases.mtz',folder/'building.mtz',cfg['free_value'])
            else:
                if stage['id']!='build':folder.mkdir()
                env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',
                         OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
                with (out/(stage['id']+'.log')).open('w') as log:runner(stage['command'],log,env)
            artifacts=[fingerprint(path) for path in stage['outputs']]
            state['completed'].append(dict(id=stage['id'],name=stage['name'],artifacts=artifacts));save()
        build_report=json.loads((out/'05_build/hybrid_report.json').read_text())
        if build_report.get('status')!='completed':raise RuntimeError('Construction/refinement endpoint is incomplete')
        model=out/'05_build/primary_model.cif';expected=build_report['primary_model']['sha256']
        if fingerprint(model)['sha256']!=expected:raise RuntimeError('Final endpoint hash differs from engine report')
        primary=out/'primary_model.cif'
        if primary.exists():
            if fingerprint(primary)['sha256']!=expected:raise ValueError('Existing primary model differs; refusing overwrite')
        else:shutil.copyfile(model,primary)
        state.update(status='completed',current='complete',current_name='全部阶段完成',primary_model=fingerprint(primary),primary_metric=build_report.get('primary_metric'),
                     experimental_prior=fingerprint(out/'02_sad/experimental.mtz'),dm_phases=fingerprint(out/'03_dm/map_phases.mtz'))
        save();return state
    except BaseException as error:
        state.update(status='cancelled' if isinstance(error,(InterruptedError,KeyboardInterrupt)) else 'failed',error=str(error));save();raise
    finally:
        if old is not None:signal.signal(signal.SIGTERM,old)
