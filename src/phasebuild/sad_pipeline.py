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
                free_value=0, seed=20261023, labels='', build_semet=False,
                method='likelihood', solvent_fraction=None, aoe=.5, rank_fit=False,
                triplet_weight=.1, kmin=.03, max_neighbors=600, max_per_target=96, build_model=True,
                observation_prior='', composition_sequence='', phenix_bin='', trust_kl=.0005, trust_cleanup=True, trust_stability=False, trust_partial_cleanup=False, trust_work_weight=False, trust_wang_radius=14., wavelength=None,
                elements='', auto_fpp=False, native_labels='', input_columns=False, auto_hand=False, building_recipe='legacy',refinement_prior='original_sad',phase_blur_scale=None,phase_blur_b=None,canonical_dm_seed=False,per_round_cleanup=False,dm_engine='resolve',dm_solvent=None,dm_mask_radius=3.0,dm_anisotropy=True,feedback_cycles=0,trace_resolution=2.0,dm_executable='')


def validate(config, inspect=True):
    cfg = {**DEFAULTS, **config}
    cfg.setdefault('native','')
    if not isinstance(cfg['auto_hand'],bool):raise ValueError('auto_hand must be boolean')
    if cfg['auto_hand'] and cfg['method']=='pplus_trust':raise ValueError('External comparison does not support automatic hands')
    if cfg['auto_hand'] and not cfg['build_model']:raise ValueError('Automatic hands require model building')
    if cfg['refinement_prior'] not in ('original_sad','dm_phifom'):raise ValueError('Unknown refinement phase source')
    if cfg['refinement_prior']=='dm_phifom' and cfg['method']!='pplus_bundled':raise ValueError('DM phase restraints currently require bundled P+')
    if not isinstance(cfg['canonical_dm_seed'],bool):raise ValueError('Canonical DM seed must be boolean')
    if cfg['canonical_dm_seed'] and cfg['method'] not in ('pplus_trust','pplus_bundled'):raise ValueError('Canonical DM seed requires original-likelihood P+ density modification')
    if not isinstance(cfg['per_round_cleanup'],bool):raise ValueError('Per-round cleanup must be boolean')
    if cfg['per_round_cleanup'] and (cfg['method']!='pplus_bundled' or cfg['refinement_prior']!='dm_phifom' or not cfg['build_model'] or cfg['copies']!=1):raise ValueError('Per-round cleanup requires single-copy bundled P+ building with DM phase restraints')
    scale,blur_b=cfg['phase_blur_scale'],cfg['phase_blur_b']
    if (scale is None)!=(blur_b is None):raise ValueError('Phase blur scale and B must be supplied together')
    if scale is not None:
        from .phase_refinement import validate_phase_blur
        cfg['phase_blur_scale'],cfg['phase_blur_b']=validate_phase_blur((scale,blur_b))
        if cfg['method']!='pplus_bundled' or cfg['refinement_prior']!='dm_phifom' or not cfg['build_model']:raise ValueError('Phase confidence correction requires bundled P+ building with DM phase restraints')
    if cfg['building_recipe'] not in ('legacy','ccp4','pipeline'):raise ValueError('Unknown building recipe')
    if not isinstance(cfg['input_columns'],bool):raise ValueError('input_columns must be boolean')
    if cfg['method']=='pplus_bundled' and (cfg.get('labels') or cfg.get('native_labels') or not cfg.get('native')):cfg['input_columns']=True
    if not isinstance(cfg['auto_fpp'],bool):raise ValueError('auto_fpp must be boolean')
    if not isinstance(cfg['elements'],str):raise ValueError('元素必须为文本，如 Se 或 S')
    if cfg['elements'].strip():cfg['auto_fpp']=True
    if cfg['auto_fpp']:
        from .anomalous_scattering import calculate
        scattering=calculate(cfg['elements'],cfg['wavelength'])
        cfg['elements']=','.join(scattering['elements']);cfg['fpp']=scattering['fpp_argument']
    if not isinstance(cfg['build_model'],bool):raise ValueError('build_model must be boolean')
    bundled=cfg.get('method')=='pplus_bundled'
    if not cfg['build_model']:
        cfg.update(rounds=DEFAULTS['rounds'],cycles=DEFAULTS['cycles'])
        if not bundled:cfg['runtime']=None
    if cfg['method'] not in ('likelihood','pplus','pplus_trust','pplus_internal','pplus_bundled'):raise ValueError('Unknown phase method')
    trusted=cfg['method']=='pplus_trust'
    optimized=cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled')
    if os.environ.get('OASIS_INSTALLATION_FREE')=='1' and trusted:raise ValueError('免安装版禁止调用 Phenix；请选择内部 P+、连续 SAD 或旧 P+')
    if cfg['method']=='pplus_internal':
        cfg.update(phenix_bin='',phase_cycles=100)
        if cfg.get('observation_prior'):raise ValueError('内部 P+ 不接受外部 Phaser 先验')
    if bundled:
        cfg.update(phenix_bin='',phase_cycles=100,dm_cycles=40)
        if cfg.get('observation_prior'):raise ValueError('自带模式不接受外部相位文件')
        cfg['runtime']=str(resolve_runtime(cfg.get('runtime') or None))
        from .bundled_pplus import engine_layout
        engine_layout(cfg['runtime'])
    if trusted:
        from .trust_workflow import phenix_directory
        cfg['phenix_bin']=phenix_directory(cfg['phenix_bin'])
        cfg.update(phase_cycles=100,dm_cycles=40)
    for flag in ['trust_stability','trust_partial_cleanup','trust_work_weight']:
        if not isinstance(cfg[flag],bool):raise ValueError(flag+' must be boolean')
        if cfg[flag] and not (trusted or bundled):raise ValueError(flag+' requires original-likelihood P+')
    radius=cfg['trust_wang_radius']
    if isinstance(radius,bool) or not isinstance(radius,(int,float)) or not math.isfinite(radius) or not 2<=radius<=30:raise ValueError('Wang radius must be 2–30 Angstrom')
    if radius!=14 and not (trusted or bundled):raise ValueError('Wang radius is active only for original-likelihood P+')
    if cfg['trust_work_weight'] and (not cfg['build_model'] or not cfg['trust_cleanup']):raise ValueError('Work-only weight calibration requires model building and cleanup enabled')
    if cfg['trust_partial_cleanup'] and (not cfg['build_model'] or not cfg['trust_cleanup']):raise ValueError('Partial cleanup requires model building and cleanup enabled')
    if cfg['dm_engine'] not in ('resolve','parrot'):raise ValueError('Unknown density-modification engine')
    if cfg['dm_engine']=='parrot' and cfg['method'] not in ('pplus_trust','pplus_internal','pplus_bundled'):
        raise ValueError('The Parrot engine currently requires a P+ phasing method')
    if cfg['dm_solvent'] is not None:
        value=cfg['dm_solvent']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<value<1:
            raise ValueError('dm_solvent must be between 0 and 1, or blank for sequence estimate')
    if isinstance(cfg['dm_mask_radius'],bool) or not isinstance(cfg['dm_mask_radius'],(int,float)) or not math.isfinite(cfg['dm_mask_radius']) or not 1<=cfg['dm_mask_radius']<=10:
        raise ValueError('dm_mask_radius must be 1-10 Angstrom')
    if not isinstance(cfg['dm_anisotropy'],bool):raise ValueError('dm_anisotropy must be boolean')
    if isinstance(cfg['feedback_cycles'],bool) or not isinstance(cfg['feedback_cycles'],int) or not 0<=cfg['feedback_cycles']<=8:
        raise ValueError('feedback_cycles must be 0-8')
    if cfg['feedback_cycles'] and cfg['dm_engine']!='parrot':
        raise ValueError('Trace feedback requires the Parrot density-modification engine')
    if config.get('trace_resolution') is not None:
        value=cfg['trace_resolution']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not 0.5<=value<=10:
            raise ValueError('trace_resolution must be 0.5-10 Angstrom')
    if cfg['dm_engine']=='parrot':
        from . import parrot_dm
        cfg['dm_executable']=parrot_dm.locate(cfg.get('dm_executable') or None,cfg.get('runtime') or None)
    if optimized:
        if isinstance(cfg['trust_kl'],bool) or not isinstance(cfg['trust_kl'],(int,float)) or not math.isfinite(cfg['trust_kl']) or not 0<cfg['trust_kl']<=.05:
            raise ValueError('P+ KL 上限须大于 0 且不超过 0.05')
        if not isinstance(cfg['trust_cleanup'],bool):raise ValueError('trust_cleanup must be boolean')
        for key in ['observation_prior','composition_sequence']:
            if key=='composition_sequence' and not cfg[key]:continue
            if key=='observation_prior' and not trusted:continue
            path=Path(cfg[key]).expanduser().resolve()
            if not path.is_file():raise ValueError('Missing input file: '+key)
            cfg[key]=str(path)
    if not isinstance(cfg['rank_fit'],bool):raise ValueError('rank_fit must be boolean')
    if cfg['method']!='pplus' and cfg['rank_fit']:raise ValueError('Rank fitting is available only for P+')
    for key,low,high in [('aoe',0,1),('triplet_weight',0,1),('kmin',0,1)]:
        value=cfg[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not low<=value<=high:
            raise ValueError('Invalid numeric parameter: '+key)
    if not 0<cfg['aoe']<1:raise ValueError('AOE must lie strictly between 0 and 1')
    if cfg['solvent_fraction'] is not None:
        value=cfg['solvent_fraction']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<value<1:
            raise ValueError('Solvent fraction must be between 0 and 1, or blank for sequence estimate')
    for key,low,high in [('copies',1,100),('phase_cycles',1,100),('dm_cycles',1,100),
                         ('rounds',1,20),('cycles',1,200),('free_value',0,999),('seed',0,2147483647),
                         ('max_neighbors',0,10000),('max_per_target',0,4096)]:
        value=cfg[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or int(value)!=value or not low<=value<=high:
            raise ValueError('Invalid integer parameter: '+key)
        cfg[key]=int(value)
    if not isinstance(cfg['build_semet'],bool):raise ValueError('build_semet must be boolean')
    if cfg['build_semet'] and not optimized:raise ValueError('Unified entry does not yet support automatic SeMet composition; use the independent workflows')
    for key in ['sad','native','sequence','heavy']+([] if trusted or bundled else ['reference']):
        if key=='native' and cfg['input_columns'] and not cfg.get(key):continue
        path=Path(cfg.get(key,'')).expanduser().resolve()
        if not path.is_file():raise ValueError('Missing input file: '+key)
        cfg[key]=str(path)
    if cfg['per_round_cleanup']:
        from .oasis_options import read_sequences
        if len(read_sequences(cfg['sequence']))!=1:raise ValueError('Per-round cleanup requires one sequence record')
    cfg['output']=str(Path(cfg['output']).expanduser().resolve())
    if cfg['build_model'] or cfg['auto_hand']:
        cfg['runtime']=str(resolve_runtime(cfg.get('runtime') or None,required=('cbuccaneer','refmac5')))
        for suffix in ['.pdb','.mtz']:
            if not (Path(cfg['runtime'])/'lib/data/reference_structures'/('reference-1tqw'+suffix)).is_file():
                raise ValueError('External construction reference unavailable')
    if cfg['auto_hand'] and cfg['method']!='pplus_bundled':
        from .bundled_pplus import engine_layout
        engine_layout(cfg['runtime'])
        if cfg.get('wavelength') is None:raise ValueError('Automatic hands require experimental SAD wavelength')
    if not cfg.get('fpp') or not isinstance(cfg['fpp'],str):raise ValueError('Specify element=f-double-prime')
    if inspect and cfg['input_columns']:
        from .mtz_inputs import checked_config
        with checked_config(cfg) as (selected,selection):validate(selected,inspect=True)
        return cfg
    if inspect:
        if sys.version_info<(3,12):raise ValueError('SAD 全流程需要 Python 3.12；请使用整合版 conda 环境启动 GUI')
        from .oasis_options import pairs
        from oasis_py.data import read_mtz, read_atoms
        from oasis_py.sequence import sequence_asu
        from oasis_py.reference_data import read_reference
        from .oasis_options import read_sequences
        sad,native,flags,partition=inspect_inputs(cfg['sad'],cfg['native'],cfg['free_value'])
        forbidden=set(native.column_labels()) & {'FWT','PHWT','PHCOMB','FOM','SAD_HLA','SAD_HLB','SAD_HLC','SAD_HLD'}
        if forbidden:raise ValueError('Native must contain original observations without downstream map/prior columns: '+', '.join(sorted(forbidden)))
        if set(pairs(cfg['labels']))-{'f','sigf','dano','sigdano','fplus','sigfplus','fminus','sigfminus'}:
            raise ValueError('Unknown SAD column mapping name')
        data,_=read_mtz(cfg['sad'],'sad',pairs(cfg['labels']))
        if bundled:
            from .heavy_scattering import read_heavy
            heavy=read_heavy(cfg['heavy'],data)
        else:heavy=read_atoms(cfg['heavy'],data)
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
        if not trusted and not bundled:read_reference(cfg['reference'])
        if optimized:
            if trusted:from .trust_workflow import preflight
            elif bundled:from .bundled_pplus import preflight
            else:from .internal_pplus import preflight
            preflight(cfg)
    return cfg


def stage_plan(cfg, python=None):
    root=Path(__file__).resolve().parents[2];out=Path(cfg['output']);python=python or sys.executable
    prepare=out/'01_inputs';phase=out/'02_sad';dm=out/'03_dm';bridge=out/'04_bridge';build=out/'05_build';hand=out/'02_hand'
    scout=cfg.get('auto_hand') and cfg['method']!='pplus_bundled'
    common=['--sad',str(hand/'selected_sad.mtz') if scout else str(prepare/'sad.mtz'),'--sequence',cfg['sequence'],'--heavy',str(hand/'refined_heavy.pdb') if scout else cfg['heavy'],
            '--fpp',cfg['fpp'],'--copies',str(cfg['copies']),'--labels',cfg['labels'],
            '--free-value',str(cfg['free_value']),'--seed',str(cfg['seed']),'--method',cfg['method']]
    phase_cmd=[python,'-u',str(root/'run_oasis_stage.py'),'phase']+common+['--cycles',str(cfg['phase_cycles']),'--output',str(phase/'experimental.mtz')]
    for key in ['aoe','triplet_weight','kmin','max_neighbors','max_per_target']:
        phase_cmd.extend(['--'+key.replace('_','-'),str(cfg[key])])
    if cfg['rank_fit']:phase_cmd.append('--rank-fit')
    phase_cmd.extend(['--observation-output',str(phase/'observations.mtz')])
    dm_cmd=[python,'-u',str(root/'run_oasis_stage.py'),'dm']+common+['--cycles',str(cfg['dm_cycles']),
        '--initial',str(phase/'experimental.mtz'),'--reference',cfg['reference'],'--output',str(dm/'map_phases.mtz')]
    if cfg['solvent_fraction'] is not None:dm_cmd.extend(['--solvent-fraction',str(cfg['solvent_fraction'])])
    build_cmd=[python,'-u',str(root/'run_build_from_phases.py'),'--map',str(bridge/'building.mtz'),
        '--native',str(hand/'selected_native.mtz') if scout else (str(phase/'selected_native.mtz') if cfg.get('auto_hand') else (str(prepare/'native.mtz') if cfg.get('input_columns') else cfg['native'])),'--sequence',str(prepare/'building.fasta'),'--phase-prior',str(dm/'map_phases.mtz') if cfg.get('refinement_prior')=='dm_phifom' else str(phase/'observations.mtz'),
        '--runtime',cfg['runtime'],'--free-value',str(cfg['free_value']),'--rounds',str(cfg['rounds']),
        '--cycles',str(cfg['cycles']),'--building-recipe',cfg.get('building_recipe','legacy'),'--refinement-prior',cfg.get('refinement_prior','original_sad'),'--output',str(build)]
    if cfg.get('per_round_cleanup'):build_cmd.extend(['--per-round-cleanup','--cleanup-map',str(dm/'map_phases.mtz'),'--asu-copies',str(cfg['copies'])])
    if cfg.get('method')=='pplus_bundled':build_cmd.append('--work-round-guard')
    if cfg.get('phase_blur_scale') is not None:build_cmd.extend(['--phase-blur-scale',str(cfg['phase_blur_scale']),'--phase-blur-b',str(cfg['phase_blur_b'])])
    if cfg['build_semet']:build_cmd.append('--build-semet')
    stages=[dict(id='prepare',name='输入与原自由集核验',directory=str(prepare),command=None,
                 outputs=[str(prepare/'sad.mtz'),str(prepare/'sad.json')]+([str(prepare/'building.fasta')] if cfg.get('build_model',True) else [])),
            dict(id='sad',name=('连续 SAD' if cfg['method']=='likelihood' else '直接法 P+')+' 定相：分别保存定相结果与观测 HL',directory=str(phase),command=phase_cmd,
                 outputs=[str(phase/'experimental.mtz'),str(phase/'experimental.json'),
                          str(phase/'observations.mtz'),str(phase/'observations.json')]),
            dict(id='dm',name=('MLHL-best' if cfg['method']=='likelihood' else 'Legacy')+' 密度修饰：生成建模相位',directory=str(dm),command=dm_cmd,
                 outputs=[str(dm/'map_phases.mtz'),str(dm/'map_phases.json')]),
            dict(id='bridge',name='按 HKL 生成工作集建模图',directory=str(bridge),command=None,
                 outputs=[str(bridge/'building.mtz'),str(bridge/'building.audit.json')]),
            dict(id='build',name='Buccaneer 构建与固定观测 HL 精修',directory=str(build),command=build_cmd,
                 outputs=[str(build/'primary_model.cif'),str(build/'hybrid_report.json')])]
    if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled'):
        def trusted_command(action,folder):
            return [python,'-u',str(root/'run_pplus_trust_stage.py'),action,'--config',str(out/'pipeline_report.json'),'--output',str(folder)]
        stages[1].update(name='内部 SAD 似然 → P+ 信任域更新' if cfg['method']=='pplus_internal' else 'P+ 优化：原始 Phaser SAD HL → KL 信任域更新',command=trusted_command('phase',phase))
        stages[2].update(name='内部 MLHL 密度修饰（无需 Phenix）' if cfg['method']=='pplus_internal' else 'RESOLVE Wang%g → classic40 密度修饰'%cfg.get('trust_wang_radius',14.),command=trusted_command('dm',dm))
        stages[3].update(name='实际 DM 地图按 native 尺度回接；排除原自由集',command=trusted_command('bridge',bridge))
        finish=out/'06_finish'
        stages.append(dict(id='finish',name=('序列锚定去重与固定后精修' if cfg.get('trust_partial_cleanup') else '单拷贝完整链核验与固定后精修') if cfg['trust_cleanup'] else '保留固定末轮模型',directory=str(finish),
            command=trusted_command('finish',finish),outputs=[str(finish/'primary_model.cif'),str(finish/'endpoint.json')]))
    if cfg['method']=='pplus_bundled':
        stages[1]['name']='自带 Phaser SAD → P+ 信任域更新'
        stages[1]['outputs'] += [str(phase/'frontend.json'),str(phase/'refined_heavy.pdb')]
        stages[2]['name']='自带 RESOLVE/Wang%g → CCTBX classic40'%cfg.get('trust_wang_radius',14.)
    if cfg.get('auto_hand') and cfg['method']=='pplus_bundled':
        stages[1]['name']='双手 P+ / Wang8 密度修饰与工作集测序判别'
        stages[1]['outputs'] += [str(phase/'hand_selection.json'),str(phase/'selected_native.mtz'),str(phase/'selected_sad.mtz')]
        stages[2]['name']='保存所选手的密度修饰结果（校验缓存）'
    if cfg.get('input_columns'):
        stages[0]['outputs'] += [str(prepare/'native.mtz'),str(prepare/'selected_sad.mtz'),str(prepare/'columns.json')]
        stages[0]['name']='输入列选择与工作集 / 自由集核验'
    if scout:
        stages.insert(1,dict(id='hand',name='双手工作集测序判别',directory=str(hand),
            command=[python,'-u',str(root/'run_hand_preselection.py'),'--config',str(out/'pipeline_report.json'),'--output',str(hand)],
            outputs=[str(hand/name) for name in ('hand_selection.json','selected_native.mtz','selected_sad.mtz','refined_heavy.pdb')]))
    if cfg.get('building_recipe')=='pipeline' and cfg.get('build_model',True):
        stages[4].update(name='CCP4 buccaneer_pipeline（与对照流程同一构建器）',
            command=[python,'-u',str(root/'run_build_pipeline.py'),'--map',str(dm/'map_phases.mtz'),
                     '--native',str(hand/'selected_native.mtz') if scout else (str(phase/'selected_native.mtz') if cfg.get('auto_hand') else (str(prepare/'native.mtz') if cfg.get('input_columns') else cfg['native'])),
                     '--sequence',str(prepare/'building.fasta'),'--runtime',cfg['runtime'],
                     '--rounds',str(cfg['rounds']),'--free-value',str(cfg['free_value']),'--output',str(build)],
            outputs=[str(build/'primary_model.cif'),str(build/'hybrid_report.json')])
    if cfg.get('dm_engine')=='parrot':
        stages[2].update(name=('Parrot 各向异性校正/溶剂掩膜密度修饰（工作集）' if not cfg.get('feedback_cycles')
                               else 'Parrot 密度修饰 + %d 轮工作集 trace/Fcalc 反馈'%cfg['feedback_cycles']),
            command=[python,'-u',str(root/'run_trace_feedback.py'),'--config',str(out/'pipeline_report.json'),'--output',str(dm)],
            outputs=[str(dm/'map_phases.mtz'),str(dm/'map_phases.json')])
    return stages if cfg.get('build_model',True) else stages[:3]


def source_identity():
    root=Path(__file__).resolve().parents[2]
    files=list((root/'src').rglob('*.py'))+[root/'run_oasis_stage.py',root/'run_sad_to_model.py',root/'run_build_from_phases.py',root/'run_pplus_trust_stage.py',root/'run_hand_preselection.py']
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
    inputs={key:fingerprint(cfg[key] or cfg['sad']) for key in ['sad','native','sequence','heavy']+([] if cfg['method'] in ('pplus_trust','pplus_bundled') else ['reference'])}
    if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled'):
        for key in ['observation_prior','composition_sequence']:
            if cfg[key]:inputs[key]=fingerprint(cfg[key])
    if cfg['method']=='pplus_trust':
        for name in ['phenix.python','mmtbx.density_modification']:
            inputs[name]=fingerprint(Path(cfg['phenix_bin'])/name)
    runtime=Path(cfg['runtime']) if cfg['build_model'] or cfg['method']=='pplus_bundled' else None
    if runtime and (cfg['build_model'] or cfg['auto_hand']):inputs.update({name:fingerprint(runtime/'bin'/name) for name in ['cbuccaneer','refmac5']})
    if cfg['method']=='pplus_bundled':inputs['native_phaser']=fingerprint(runtime/'bin/phaser')
    # Lock exact resources too; replacement reference/dictionaries cannot silently
    # enter a resumed experiment. The complete runtime manifest is preferred.
    if runtime:
        manifest=runtime/'runtime_manifest.json'
        if manifest.is_file():inputs['runtime_manifest']=fingerprint(manifest)
        if cfg['build_model'] or cfg['auto_hand']:
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
                   completed=[],created=time.time(),runtime=runtime_info(runtime) if runtime else None,environment=environment,
                   scope='SAD and DM use working reflections only; original flags unchanged. Later model-derived maps are not certified work-only; Rfree is not an independent end-to-end blind score.')
    def save():
        state['updated']=time.time();temp=state_path.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');temp.replace(state_path)
    if cfg.get('auto_fpp'):
        from .anomalous_scattering import calculate
        state['anomalous_scattering']=calculate(cfg['elements'],cfg['wavelength'])
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
                if cfg.get('input_columns'):
                    from .mtz_inputs import prepare
                    selection=prepare(cfg,folder)
                    state['column_selection']=selection['column_selection'];save()
                else:prepare_sad(cfg['sad'],cfg['native'],folder/'sad.mtz',cfg['free_value'])
                if cfg['build_model']:
                    from .oasis_options import read_sequences
                    lines=[]
                    for copy in range(cfg['copies']):
                        for i,(name,seq) in enumerate(read_sequences(cfg['sequence'])):
                            lines.extend(['>copy{}_chain{}'.format(copy+1,i+1),seq])
                    (folder/'building.fasta').write_text('\n'.join(lines)+'\n')
            elif stage['id']=='bridge' and cfg['method'] not in ('pplus_trust','pplus_internal','pplus_bundled'):
                folder.mkdir();adapt_maps(str(out/'02_hand/selected_native.mtz') if cfg.get('auto_hand') and cfg['method']!='pplus_bundled' else (str(out/'01_inputs/native.mtz') if cfg.get('input_columns') else cfg['native']),out/'02_sad/observations.mtz',out/'03_dm/map_phases.mtz',folder/'building.mtz',cfg['free_value'])
            else:
                if stage['id']!='build':folder.mkdir()
                env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',
                         OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
                log_path=out/(stage['id']+'.log')
                if log_path.exists():
                    preserved=out/(stage['id']+'.interrupted-'+uuid.uuid4().hex+'.log')
                    log_path.rename(preserved)
                    state.setdefault('preserved_stage_logs',[]).append(
                        dict(stage=stage['id'],artifact=fingerprint(preserved)))
                    save()
                with log_path.open('x') as log:runner(stage['command'],log,env)
            artifacts=[fingerprint(path) for path in stage['outputs']]
            state['completed'].append(dict(id=stage['id'],name=stage['name'],artifacts=artifacts));save()
        if not cfg['build_model']:
            state.update(status='completed',current='complete',current_name='定相与密度修饰完成（未启用建模）',
                         primary_model=None,primary_metric=None,experimental_prior=fingerprint(out/'02_sad/observations.mtz'),
                         pre_dm_result=fingerprint(out/'02_sad/experimental.mtz'),
                         dm_phases=fingerprint(out/'03_dm/map_phases.mtz'))
            save();return state
        build_report=json.loads((out/'05_build/hybrid_report.json').read_text())
        if build_report.get('status')!='completed':raise RuntimeError('Construction/refinement endpoint is incomplete')
        endpoint=build_report
        if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled'):endpoint=json.loads((out/'06_finish/endpoint.json').read_text())
        model=out/('06_finish/primary_model.cif' if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled') else '05_build/primary_model.cif');expected=endpoint['primary_model']['sha256']
        if fingerprint(model)['sha256']!=expected:raise RuntimeError('Final endpoint hash differs from engine report')
        primary=out/'primary_model.cif'
        if primary.exists():
            if fingerprint(primary)['sha256']!=expected:raise ValueError('Existing primary model differs; refusing overwrite')
        else:shutil.copyfile(model,primary)
        state.update(status='completed',current='complete',current_name='全部阶段完成',primary_model=fingerprint(primary),primary_metric=endpoint.get('primary_metric'),
                     experimental_prior=fingerprint(out/'02_sad/observations.mtz'),
                     pre_dm_result=fingerprint(out/'02_sad/experimental.mtz'),dm_phases=fingerprint(out/'03_dm/map_phases.mtz'))
        save();return state
    except BaseException as error:
        state.update(status='cancelled' if isinstance(error,(InterruptedError,KeyboardInterrupt)) else 'failed',error=str(error));save();raise
    finally:
        if old is not None:signal.signal(signal.SIGTERM,old)
