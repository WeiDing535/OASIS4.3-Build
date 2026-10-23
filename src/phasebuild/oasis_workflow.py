"""Focused OASIS4.3 workflow: two phase methods, paired DM, Buccaneer/REFMAC."""
import json
import math
from pathlib import Path
import sys
from .sad_pipeline import DEFAULTS, validate, stage_plan
from .runtime import default_runtime
from .data import fingerprint

METHODS={'likelihood':('连续 SAD','MLHL-best'),'pplus':('直接法 P+','Legacy')}


def defaults(project):
    project=Path(project)
    candidates=[project/'runtime',project/'dist/OASIS4.3-Build-macOS-Intel-0.8.3-local/runtime']
    runtime=next((p for p in candidates if (p/'bin/cbuccaneer').is_file() and (p/'bin/refmac5').is_file()),None)
    runtime=str(runtime) if runtime else (default_runtime() or '')
    cfg={k:v for k,v in DEFAULTS.items() if k!='build_semet'}
    cfg.update(sad='',native='',sequence='',heavy='',fpp='',runtime=runtime,
               reference=str(Path(runtime)/'lib/data/reference_structures/reference-1tqw.mtz') if runtime else '',
               output=str(project/'runs/oasis_run'),mode='sad_to_model')
    return cfg


def make_plan(config, project, python=None, inspect=True):
    if not isinstance(config,dict):raise ValueError('配置必须是 JSON 对象')
    cfg=defaults(project)
    unknown=set(config)-set(cfg)
    if unknown:raise ValueError('此版本不支持配置项：'+', '.join(sorted(unknown)))
    cfg.update(config)
    if cfg.pop('mode')!='sad_to_model':raise ValueError('OASIS4.3 仅支持 SAD → 密度修饰 → Buccaneer → REFMAC')
    if not isinstance(cfg['output'],str) or not cfg['output'].strip():raise ValueError('请填写新的输出目录')
    outer=Path(cfg['output']).expanduser().resolve()
    if outer.exists():raise ValueError('输出目录已存在，请选择新目录')
    cfg['output']=str(outer/'00_sad_to_model')
    checked=validate(cfg,inspect=inspect)
    command=[python or sys.executable,'-u',str(Path(project).resolve()/'run_sad_to_model.py')]
    for key,value in checked.items():
        if key=='build_semet':continue
        if key=='rank_fit':
            if value:command.append('--rank-fit')
        elif value is not None:command.extend(['--'+key.replace('_','-'),str(value)])
    checked.pop('build_semet')
    checked.update(output=str(outer),mode='sad_to_model')
    return checked,[dict(name='OASIS4.3：定相 → 密度修饰 → Buccaneer → REFMAC',kind='sad_to_model',
                         output=str(outer/'00_sad_to_model'),command=command)]


def describe_plan(cfg):
    values={**DEFAULTS,**cfg,'output':str(Path(cfg['output'])/'00_sad_to_model')}
    phase,dm=METHODS[cfg['method']]
    active=dict(method=phase,dm_engine=dm,solvent_fraction=cfg['solvent_fraction'],
                solvent_source='按序列与 ASU 份数估算' if cfg['solvent_fraction'] is None else '用户指定',
                copies=cfg['copies'],phase_cycles=cfg['phase_cycles'],dm_cycles=cfg['dm_cycles'],
                rounds=cfg['rounds'],refmac_cycles=cfg['cycles'],seed=cfg['seed'],
                free_value=cfg['free_value'],kmin=cfg['kmin'],max_neighbors=cfg['max_neighbors'],
                max_per_target=cfg['max_per_target'])
    active.update({'triplet_weight':cfg['triplet_weight']} if cfg['method']=='likelihood' else {'aoe':cfg['aoe'],'rank_fit':cfg['rank_fit']})
    import gemmi
    from oasis_py.sequence import sequence_asu
    from .refinement import reflection_partition
    data=gemmi.read_mtz_file(cfg['sad'])
    sample=sequence_asu(cfg['sequence'],cfg['copies'],data.cell,data.spacegroup)
    _,partition=reflection_partition(cfg['native'],free_value=cfg['free_value'])
    if active['solvent_fraction'] is None:active['solvent_fraction']=sample['solvent_fraction']
    active.update(asu_residues=sample['asu_residues'],asu_chains=sample['asu_chains'],
                  native_work=partition['usable_work'],native_free=partition['usable_free'])
    return dict(stages=stage_plan(values),effective_options=active,sample_warnings=sample['warnings'],
                fixed_protocol='Buccaneer 首轮 3 周期、后续每轮 2 周期；单线程；REFMAC 固定原始 SAD HL；固定最后一轮结果。')


def read_json(path):
    try:return json.loads(Path(path).read_text())
    except (OSError,ValueError):return {}


def pipeline_root(directory):
    root=Path(directory).resolve()
    return root if (root/'pipeline_report.json').is_file() else root/'00_sad_to_model'


def validate_import(directory):
    state=read_json(pipeline_root(directory)/'pipeline_report.json')
    if state.get('engine')!='oasis-sad-to-model' or state.get('config',{}).get('method','likelihood') not in METHODS:
        raise ValueError('请选择含 OASIS4.3 pipeline_report.json 的流程结果目录')


def phase_statistics(path):
    if not path.is_file():return None
    import gemmi,numpy as np
    m=gemmi.read_mtz_file(str(path));a=np.asarray(m.column_with_label('FOM'),dtype=float)
    valid=a[np.isfinite(a)&(a>=0)&(a<=1)]
    return dict(reflections=m.nreflections,mean_fom=float(valid.mean()) if len(valid) else None,
                resolution=m.resolution_high(),spacegroup=m.spacegroup.xhm())


def collect_results(directory):
    root=Path(directory).resolve()
    if not root.is_dir():raise ValueError('结果目录不存在')
    pipe=pipeline_root(root);state=read_json(pipe/'pipeline_report.json')
    cfg=state.get('config',{});completed={item['id'] for item in state.get('completed',[])}
    report=read_json(pipe/'05_build/hybrid_report.json')
    labels=[('prepare','输入核验'),('sad','SAD 定相'),('dm','密度修饰'),('bridge','建模地图'),('build','Buccaneer / REFMAC')]
    stages=[]
    for key,label in labels:
        status='completed' if key in completed else state.get('status','waiting') if state.get('current')==key else 'waiting'
        if status=='completed' and key=='build' and state.get('status')!='completed':status='running'
        stages.append(dict(id=key,name=label,status=status))
    metrics=[]
    for number in range(1,cfg.get('rounds',0)+1):
        folder=pipe/'05_build'/('round_%02d_refine'%number)
        ref=read_json(folder/'refinement_report.json')
        if ref.get('status')!='success':continue
        stats=ref.get('statistics',{}).get('final',{})
        metrics.append(dict(stage='第 %d 轮'%number,round=number,r_work=stats.get('r_work'),r_free=stats.get('r_free'),
            cycles=ref.get('requested_cycles'),free_preserved=ref.get('free_flags_preserved'),
            engine='REFMAC · 原始 SAD HL',model=str((folder/'model.cif').relative_to(root))))
    complete=state.get('status')=='completed' and report.get('status')=='completed'
    model=pipe/'primary_model.cif' if complete else root/metrics[-1]['model'] if metrics else None
    if complete:
        expected=state.get('primary_model',{}).get('sha256')
        if not model.is_file() or not expected or fingerprint(model)['sha256']!=expected:
            raise ValueError('最终模型与流程记录的哈希不一致')
        if len(metrics)!=cfg['rounds']:raise ValueError('最终模型缺少完整的逐轮精修报告')
    preview=[]
    if model and model.is_file():
        import gemmi
        structure=gemmi.read_structure(str(model))
        if len(structure):
            for chain in structure[0]:
                points=[]
                for residue in chain:
                    ca=residue.find_atom('CA','*')
                    if ca:points.append([ca.pos.x,ca.pos.y,ca.pos.z,str(residue.seqid),residue.name])
                if points:preview.append(dict(chain=chain.name,points=points))
    artifacts=[]
    def add(path,role,group):
        if path.is_file():artifacts.append(dict(path=str(path.relative_to(root)),role=role,group=group))
    add(pipe/'primary_model.cif','最终模型（固定末轮）','主要结果') if complete else None
    for key,name,role in [('sad','02_sad/experimental.mtz','原始 SAD 相位与 HL'),
                         ('dm','03_dm/map_phases.mtz','密度修饰后相位与 FOM'),
                         ('bridge','04_bridge/building.mtz','初始 Buccaneer 建模地图')]:
        if key in completed:add(pipe/name,role,'相位与地图')
    if complete:add(pipe/'05_build'/('round_%02d_refine'%cfg['rounds'])/'refined.mtz','最后一轮精修地图','主要结果')
    for name in ['pipeline_report.json','01_inputs/sad.json','02_sad/experimental.json','03_dm/map_phases.json',
                 '04_bridge/building.audit.json','05_build/hybrid_report.json']:
        add(pipe/name,Path(name).name,'运行记录')
    for name in ['gui.log']:
        add(root/name,'流程日志','日志')
    for name in ['sad.log','dm.log','build.log']:
        add(pipe/name,name,'日志')
    for number in range(1,cfg.get('rounds',0)+1):
        for suffix,file,role in [('build','buccaneer.log','Buccaneer 日志'),('refine','refmac.log','REFMAC 日志'),
                                 ('refine','refinement_report.json','精修统计'),('build','buccaneer.xml','Buccaneer 统计')]:
            add(pipe/'05_build'/('round_%02d_%s'%(number,suffix))/file,'第 %d 轮 '%number+role,'逐轮记录')
    raw=phase_statistics(pipe/'02_sad/experimental.mtz') if 'sad' in completed else None
    dm=phase_statistics(pipe/'03_dm/map_phases.mtz') if 'dm' in completed else None
    dm_meta=read_json(pipe/'03_dm/map_phases.json') if 'dm' in completed else {}
    raw_meta=read_json(pipe/'02_sad/experimental.json') if 'sad' in completed else {}
    history=dm_meta.get('density_modification',{}).get('cycles',[])
    partition=read_json(pipe/'01_inputs/sad.json').get('native_partition',{})
    method=cfg.get('method','likelihood')
    return dict(root=str(root),status=state.get('status','waiting'),method=method,
                method_name=METHODS[method][0],dm_engine=METHODS[method][1],config=cfg,
                stages=stages,sad=raw,dm=dm,dm_history=history,solvent_fraction=dm_meta.get('solvent_fraction_used'),
                solvent_source=dm_meta.get('solvent_source'),partition=partition,metrics=metrics,
                sample=raw_meta.get('sequence_asu',{}),
                primary_metric=metrics[-1] if metrics else None,primary_complete=complete,
                final_model=str(model.relative_to(root)) if model and model.is_file() else None,
                preview=preview,residues=sum(len(c['points']) for c in preview),chains=len(preview),
                artifacts=artifacts,error=state.get('error'),selection_policy='fixed_final_round')


def save_selected_result(directory):
    result=collect_results(directory)
    if not result['primary_complete']:raise ValueError('流程尚未产生完整的最终模型')
    return result


def enrich_job(job):
    pipe=pipeline_root(job['output']);state=read_json(pipe/'pipeline_report.json')
    if not state:return job
    current=state.get('current','');build=read_json(pipe/'05_build/hybrid_report.json')
    round_number=len(build.get('stages',[]))+1
    label=state.get('current_name','')
    logs=[Path(job['output'])/'gui.log']
    if current in ('sad','dm','build'):logs.append(pipe/(current+'.log'))
    if current=='build' and build.get('status')=='running':
        sub=build.get('current','buccaneer')
        label='第 %d / %d 轮 · %s'%(round_number,state['config']['rounds'],'REFMAC 精修' if sub=='refmac' else 'Buccaneer 构建')
        folder=pipe/'05_build'/('round_%02d_%s'%(round_number,'refine' if sub=='refmac' else 'build'))
        logs+=list(folder.glob('*.log'))
    job['pipeline_progress']=dict(completed=len(state.get('completed',[])),total=5,current=label,
                                  revision=[state.get('updated'),build.get('current'),len(build.get('stages',[]))])
    parts=[]
    for path in logs:
        if not path.is_file():continue
        with path.open('rb') as handle:
            handle.seek(max(0,path.stat().st_size-24000));tail=handle.read().decode(errors='replace')
        parts.append('── '+str(path.relative_to(Path(job['output'])))+' ──\n'+tail)
    job['log']='\n\n'.join(parts)[-60000:]
    return job
