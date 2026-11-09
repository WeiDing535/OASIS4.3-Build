"""Focused OASIS4.3 workflow: explicit phase methods, paired DM, Buccaneer/REFMAC."""
import json
import math
from pathlib import Path
import sys
from .sad_pipeline import DEFAULTS, validate, stage_plan
from .runtime import default_runtime
from .data import fingerprint

METHODS={'likelihood':('连续 SAD','MLHL-best'),'pplus':('旧版直接法 P+','Legacy'),'pplus_trust':('外部对照（Phaser SAD 前端）','Wang14 + classic40'),'pplus_internal':('内部 P+（实验性）','内部 MLHL-best'),'pplus_bundled':('P+ 自带引擎版','自带 RESOLVE + CCTBX')}


def defaults(project):
    project=Path(project)
    candidates=[project/'runtime',project/'dist/OASIS4.3-dev12-portable local/runtime',project/'dist/OASIS4.3-dev11-portable local/runtime',project/'dist/OASIS4.3-dev10-portable local/runtime',project/'dist/OASIS4.3-dev9-portable local/runtime',project/'dist/OASIS4.3-dev8-portable local/runtime',project/'dist/OASIS4.3-dev7-portable local/runtime',project/'dist/OASIS4.3-dev6-portable local/runtime',project/'dist/OASIS4.3-dev5d-portable local/runtime',project/'dist/OASIS4.3-dev5-portable local/runtime',project/'dist/OASIS4.3-Build-macOS-Intel-0.8.3-local/runtime']
    runtime=next((p for p in candidates if (p/'bin/cbuccaneer').is_file() and (p/'bin/refmac5').is_file()),None)
    runtime=str(runtime) if runtime else (default_runtime() or '')
    cfg=dict(DEFAULTS)
    cfg['phenix_bin']=''
    cfg['method']='pplus_bundled'
    cfg['dm_cycles']=40
    cfg.update(sad='',native='',sequence='',heavy='',fpp='',runtime=runtime,
               reference=str(Path(runtime)/'lib/data/reference_structures/reference-1tqw.mtz') if runtime else '',
               output=str(project/'runs/oasis_run'),mode='sad_to_model',
               explore=False,explore_workers=2,protocol='pplus_w8_weight',elements='',auto_fpp=True,input_columns=True,native_labels='',auto_hand=True,building_recipe='ccp4',refinement_prior='dm_phifom',
               trust_wang_radius=8.,trust_partial_cleanup=True,trust_work_weight=True,trust_stability=False)
    return cfg


def make_plan(config, project, python=None, inspect=True):
    if not isinstance(config,dict):raise ValueError('配置必须是 JSON 对象')
    cfg=defaults(project)
    unknown=set(config)-set(cfg)
    if unknown:raise ValueError('此版本不支持配置项：'+', '.join(sorted(unknown)))
    cfg.update(config)
    from .portfolio import PROTOCOLS
    protocol=cfg.pop('protocol') if 'protocol' in config else None
    if protocol is None:cfg.pop('protocol',None)
    elif protocol not in PROTOCOLS:raise ValueError('未知计算模式：'+str(protocol))
    else:
        if config.get('method','pplus_bundled')!='pplus_bundled':raise ValueError('P+ 参数方案仅适用于自带 P+；其他定相方法请勿附带 protocol')
        cfg.update(PROTOCOLS[protocol])
    if 'auto_hand' not in config:cfg['auto_hand']=False
    if 'building_recipe' not in config:cfg['building_recipe']='legacy'
    if 'refinement_prior' not in config:cfg['refinement_prior']='original_sad'
    if 'auto_fpp' not in config and 'elements' not in config:
        cfg['auto_fpp']=False
        if 'input_columns' not in config:cfg['input_columns']=False
        for key in ['trust_wang_radius','trust_partial_cleanup','trust_work_weight','trust_stability']:
            if key not in config and protocol is None:cfg[key]=DEFAULTS[key]
    if protocol is not None and not cfg['build_model']:cfg.update(trust_work_weight=False,trust_partial_cleanup=False)
    if 'method' not in config:cfg['method']='pplus_bundled' if protocol is not None else DEFAULTS['method'] # old saved configs retain their historical route
    if 'explore' not in config:cfg['explore']=False # never silently migrate saved single-run configs
    explore=cfg.pop('explore');workers=cfg.pop('explore_workers')
    if not isinstance(explore,bool):raise ValueError('多流程探索选项必须为布尔值')
    if isinstance(workers,bool) or workers not in (1,2):raise ValueError('多流程并发数只能为 1 或 2')
    if explore and (cfg['method']!='pplus_bundled' or not cfg['build_model']):raise ValueError('多流程探索需要 P+ 自带引擎版并启用建模')
    if cfg.pop('mode')!='sad_to_model':raise ValueError('OASIS4.3 仅支持 SAD → 密度修饰 → Buccaneer → REFMAC')
    if not isinstance(cfg['output'],str) or not cfg['output'].strip():raise ValueError('请填写新的输出目录')
    outer=Path(cfg['output']).expanduser().resolve()
    if outer.exists():raise ValueError('输出目录已存在，请选择新目录')
    cfg['output']=str(outer/'00_sad_to_model')
    checked=validate(cfg,inspect=inspect)
    if explore:
        from .portfolio import plan
        plan(checked,outer/'00_portfolio',inspect=inspect)
        command=[python or sys.executable,'-B','-u',str(Path(project).resolve()/'run_sad_portfolio.py'),
                 '--gui-config',str(outer/'gui_job.json'),'--output',str(outer/'00_portfolio'),'--workers',str(workers)]
        checked.update(output=str(outer),mode='sad_to_model',explore=True,explore_workers=workers)
        if protocol is not None:checked['protocol']=protocol
        return checked,[dict(name='P+ 多流程探索与候选推荐',kind='sad_portfolio',output=str(outer/'00_portfolio'),command=command)]
    command=[python or sys.executable,'-u',str(Path(project).resolve()/'run_sad_to_model.py')]
    for key,value in checked.items():
        if key=='build_semet':
            if value:command.append('--build-semet')
        elif key in ('trust_cleanup','trust_stability','trust_partial_cleanup','trust_work_weight','auto_fpp','input_columns','auto_hand','canonical_dm_seed','per_round_cleanup','dm_anisotropy'):command.append('--'+key.replace('_','-') if value else '--no-'+key.replace('_','-'))
        elif key=='build_model':command.append('--build-model' if value else '--no-build-model')
        elif key=='rank_fit':
            if value:command.append('--rank-fit')
        elif value is not None:command.extend(['--'+key.replace('_','-'),str(value)])
    checked.update(output=str(outer),mode='sad_to_model',explore=False,explore_workers=workers)
    if protocol is not None:checked['protocol']=protocol
    return checked,[dict(name='OASIS4.3：定相 → 密度修饰'+(' → Buccaneer → REFMAC' if checked['build_model'] else '（未启用建模）'),kind='sad_to_model',
                         output=str(outer/'00_sad_to_model'),command=command)]


def describe_plan(cfg):
    if cfg.get('explore'):
        from .portfolio_gui import describe
        return describe(cfg)
    if cfg.get('input_columns'):
        from .mtz_inputs import checked_config
        with checked_config(cfg) as (selected,selection):result=describe_plan(selected)
        values={**DEFAULTS,**cfg,'output':str(Path(cfg['output'])/'00_sad_to_model')}
        result['stages']=stage_plan(values)
        result['column_selection']=selection
        selected_lines=['SAD '+k+' = '+v['label']+' ['+v['dataset']+'; 列 '+str(v['index'])+']' for k,v in selection['sad_columns'].items()]
        selected_lines += ['观测 '+k+' = '+v['label']+' ['+v['dataset']+'; 列 '+str(v['index'])+']' for k,v in selection['native_columns'].items()]
        if not selection['native_columns']:selected_lines.append('平均观测：由所选 SAD Friedel 振幅对生成')
        selected_lines.append('观测文件：'+('使用 SAD 文件（未另选观测文件）' if selection['native_omitted'] else selection['native']['path']))
        selected_lines.append('自由集：'+('缺失，将新建可重复的 5% 自由集；Rfree 不可直接与原自由集比较' if selection['flags_generated'] else '保留 '+selection['free_column']['label']+'，自由值 '+str(selection['free_value'])))
        result['effective_options']['selected_columns']='\n'.join(selected_lines)
        if selection['flags_generated']:result['sample_warnings'].append('未提供已有自由集：将生成可重复的 5% 自由集；Rfree 不可与原有精修结果直接比较。')
        return result
    if cfg.get('explore'):
        from .portfolio_gui import describe
        return describe(cfg)
    values={**DEFAULTS,**cfg,'output':str(Path(cfg['output'])/'00_sad_to_model')}
    phase,dm=METHODS[cfg['method']]
    active=dict(method=phase,dm_engine=dm,solvent_fraction=cfg['solvent_fraction'],
                solvent_source='按序列与 ASU 份数估算' if cfg['solvent_fraction'] is None else '用户指定',
                copies=cfg['copies'],phase_cycles=cfg['phase_cycles'],dm_cycles=cfg['dm_cycles'],
                rounds=cfg['rounds'],refmac_cycles=cfg['cycles'],seed=cfg['seed'],
                free_value=cfg['free_value'],kmin=cfg['kmin'],max_neighbors=cfg['max_neighbors'],
                max_per_target=cfg['max_per_target'])
    if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled'):
        active.update(trust_work_weight=cfg.get('trust_work_weight',False),trust_stability=cfg.get('trust_stability',False),trust_partial_cleanup=cfg.get('trust_partial_cleanup',False),trust_kl=cfg['trust_kl'],trust_cleanup=cfg['trust_cleanup'],frontend='OASIS 内部异常差似然（不读取 Phaser 文件）' if cfg['method']=='pplus_internal' else '导入原始 Phaser SAD 纯观测似然（不重新运行 Phaser）',
                      build_semet=cfg['build_semet'],phenix_bin=cfg['phenix_bin'])
    else:active.update({'triplet_weight':cfg['triplet_weight']} if cfg['method']=='likelihood' else {'aoe':cfg['aoe'],'rank_fit':cfg['rank_fit']})
    if cfg['method']=='pplus_bundled':active.update(frontend='自带原生 Phaser：从工作集计算 SAD 似然',wavelength=cfg['wavelength'],engine_origin='第三方内核随包提供；无需额外安装套件')
    if cfg.get('auto_fpp'):
        from .anomalous_scattering import calculate
        scattering=calculate(cfg['elements'],cfg['wavelength'])
        active.update(elements=cfg['elements'],fpp=cfg['fpp'],fpp_source=scattering['source'],
                      fpp_kind='由必填实验波长计算的理论值，不是吸收边实测值')
    if cfg['method']=='pplus_bundled' and cfg.get('trust_wang_radius')==8 and cfg.get('trust_work_weight') and not cfg.get('trust_stability'):
        active.update(method='P+ Wang8 + 精修权重搜索',dm_engine='Wang8 + classic40')
    if cfg['method'] in ('pplus_trust','pplus_bundled'):active['trust_wang_radius']=cfg.get('trust_wang_radius',14.)
    if cfg['method']!='pplus_trust':active.pop('phenix_bin',None)
    from .portfolio import PROTOCOLS
    from .portfolio_gui import NAMES
    active_profile=cfg.get('protocol') or next((name for name,options in PROTOCOLS.items() if all(cfg.get(k)==v for k,v in options.items())),None)
    if active_profile and cfg['method']=='pplus_bundled':
        active['protocol']=active_profile
        active['method']=NAMES[active_profile]
        active['dm_engine']='Wang%g + classic40'%cfg['trust_wang_radius']
    active['per_round_cleanup']=cfg.get('per_round_cleanup',False)
    active['canonical_dm_seed']=cfg.get('canonical_dm_seed',False)
    active['auto_hand']=cfg.get('auto_hand',False)
    active['refinement_prior']=cfg.get('refinement_prior','original_sad')
    active['building_recipe']=cfg.get('building_recipe','legacy')
    if cfg.get('phase_blur_scale') is not None:active['phase_blur']={'SCBLUR':cfg['phase_blur_scale'],'BBLUR':cfg['phase_blur_b'],'policy':'预先声明的 DM 相位置信度修正，不按 Rfree 选择'}
    active['build_model']=cfg.get('build_model',True)
    if not active['build_model']:
        active.pop('rounds');active.pop('refmac_cycles')
    import gemmi
    from oasis_py.sequence import sequence_asu
    from .refinement import reflection_partition
    data=gemmi.read_mtz_file(cfg['sad'])
    sample=sequence_asu(cfg['sequence'],cfg['copies'],data.cell,data.spacegroup)
    if active.get('copies')=='auto':active['copies']=sample['copies']
    if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled'):
        from .trust_workflow import sample_info
        sample=sample_info(cfg,data)
    _,partition=reflection_partition(cfg['native'],free_value=cfg['free_value'])
    if active['solvent_fraction'] is None:active['solvent_fraction']=sample['solvent_fraction']
    active.update(asu_residues=sample['asu_residues'],asu_chains=sample['asu_chains'],
                  native_work=partition['usable_work'],native_free=partition['usable_free'])
    protocol_extra=' P+ 优化另执行原始 HL 条件分支积分（1440 点），最多 100 次更新；默认 KL=0.0005。清理启用时，仅唯一完整链可进入固定 30 周期、权重 0.005 的后精修；未满足条件时保留原末轮结果。' if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled') else ''
    if cfg.get('auto_hand'):protocol_extra+=' 同一流程内部计算双手，各作工作集三周期试建，以唯一测序残基及明确优势判别；证据不足时保留两手并报告未判定，不按 Rfree 选手。'
    if cfg.get('refinement_prior')=='dm_phifom':protocol_extra+=' 精修使用冻结的工作集密度修饰相位/FOM近似约束，标记为 DM 来源；原始 SAD 似然另存，不冒称同一种先验。'
    if cfg.get('per_round_cleanup'):protocol_extra+=' 每轮精修前按固定工作集DM地图与唯一序列锚清理片段；不满足门控时保留当轮原模型，不按Rfree决定。'
    if cfg.get('phase_blur_scale') is not None:protocol_extra+=' DM 相位约束使用预先声明的 SCBLUR=%g、BBLUR=%g；原始相位文件不变，相位角不旋转。'%(cfg['phase_blur_scale'],cfg['phase_blur_b'])
    if cfg.get('trust_partial_cleanup'):protocol_extra=protocol_extra.replace('仅唯一完整链可进入固定 30 周期、','唯一完整链或满足部分模型门控的片段可进入固定 30 周期、')
    if cfg.get('trust_work_weight'):protocol_extra=protocol_extra.replace('权重 0.005','由工作集几何门控选择的固定权重（0.005–0.04）')
    if cfg.get('trust_stability'):protocol_extra+=' 关系稳定性收缩开启；这不是独立观测交叉验证。'
    if cfg['method']=='pplus_internal':protocol_extra+=' 内部异常差似然与内部 MLHL DM 均为近似实现，非 Phaser/RESOLVE 等价替代；此前外部流程的 27.42% 不适用于此模式。'
    return dict(stages=stage_plan(values),effective_options=active,sample_warnings=sample['warnings'],
                fixed_protocol=('Buccaneer 首轮 3 周期、后续每轮 2 周期；单线程；REFMAC '+('冻结工作集 DM 相位/FOM 近似约束' if cfg.get('refinement_prior')=='dm_phifom' else '固定原始 SAD 观测 HL')+'；固定最后一轮结果。'+protocol_extra) if active['build_model'] else
                '未启用 Buccaneer / REFMAC：流程在密度修饰后结束，不检查或调用建模精修引擎；分别保留定相结果、纯观测 HL、最终相位、FOM 与日志。')


def read_json(path):
    try:return json.loads(Path(path).read_text())
    except (OSError,ValueError):return {}


def pipeline_root(directory):
    root=Path(directory).resolve()
    return root if (root/'pipeline_report.json').is_file() else root/'00_sad_to_model'


def validate_import(directory):
    from .portfolio_gui import locate
    if locate(directory) is not None:
        from .portfolio_gui import validate_import as validate_portfolio
        return validate_portfolio(directory)
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
    from .portfolio_gui import locate,collect_results as collect_portfolio
    if locate(root) is not None or read_json(root/'gui_job.json').get('config',{}).get('explore'):
        return collect_portfolio(root)
    pipe=pipeline_root(root);state=read_json(pipe/'pipeline_report.json')
    cfg=state.get('config') or read_json(root/'gui_job.json').get('config',{})
    completed={item['id'] for item in state.get('completed',[])}
    prior_label='观测似然 HL' if any(Path(a['path']).name=='observations.mtz'
        for item in state.get('completed',[]) for a in item.get('artifacts',[])) else 'DM 前 HL（旧协议）'
    if cfg.get('refinement_prior')=='dm_phifom':prior_label='DM 工作集相位/FOM 约束'
    building=cfg.get('build_model',True)
    report=read_json(pipe/'05_build/hybrid_report.json') if building else {}
    labels=[('prepare','输入核验'),('sad','P+ 优化' if cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else 'SAD 定相'),('dm','密度修饰'),('bridge','建模地图'),('build','Buccaneer / REFMAC')]
    if cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled'):labels.append(('finish','序列锚定去重 / 后精修' if cfg.get('trust_partial_cleanup') else '完整链清理 / 后精修'))
    if not building:labels=labels[:3]
    stages=[]
    for key,label in labels:
        status='completed' if key in completed else state.get('status','waiting') if state.get('current')==key else 'waiting'
        if status=='completed' and key=='build' and state.get('status')!='completed':status='running'
        stages.append(dict(id=key,name=label,status=status))
    metrics=[]
    for number in range(1,(cfg.get('rounds',0) if building else 0)+1):
        folder=pipe/'05_build'/('round_%02d_refine'%number)
        ref=read_json(folder/'refinement_report.json')
        if ref.get('status')!='success':continue
        stats=ref.get('statistics',{}).get('final',{})
        metrics.append(dict(stage='第 %d 轮'%number,round=number,r_work=stats.get('r_work'),r_free=stats.get('r_free'),
            cycles=ref.get('requested_cycles'),free_preserved=ref.get('free_flags_preserved'),
            engine='REFMAC · '+prior_label,model=str((folder/'model.cif').relative_to(root))))
    endpoint=read_json(pipe/'06_finish/endpoint.json') if cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else {}
    post_label='部分模型去重后精修' if endpoint.get('cleanup_kind')=='explicit-single-copy-sequence-anchored-partial-model-reduction' else '完整链后精修'
    if endpoint.get('cleanup_status')=='completed' and 'finish' in completed:
        ref=read_json(pipe/'06_finish/refinement_report.json');stats=ref.get('statistics',{}).get('final',{})
        metrics.append(dict(stage=post_label,round=cfg['rounds']+1,r_work=stats.get('r_work'),r_free=stats.get('r_free'),
            cycles=30,free_preserved=ref.get('free_flags_preserved'),engine='REFMAC · 振幅似然 / 固定权重',model=str((pipe/'06_finish/model.cif').relative_to(root))))
    pipeline_complete=state.get('status')=='completed'
    if pipeline_complete and not building:
        expected=state.get('dm_phases',{}).get('sha256')
        phase_file=pipe/'03_dm/map_phases.mtz'
        if not phase_file.is_file() or not expected or fingerprint(phase_file)['sha256']!=expected:
            raise ValueError('最终密度修饰相位与流程记录的哈希不一致')
    complete=pipeline_complete and building and report.get('status')=='completed'
    model=pipe/'primary_model.cif' if complete else root/metrics[-1]['model'] if metrics else None
    if complete:
        expected=state.get('primary_model',{}).get('sha256')
        if not model.is_file() or not expected or fingerprint(model)['sha256']!=expected:
            raise ValueError('最终模型与流程记录的哈希不一致')
        if len(metrics)!=cfg['rounds']+(endpoint.get('cleanup_status')=='completed'):raise ValueError('最终模型缺少完整的逐轮精修报告')
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
    add(pipe/'primary_model.cif',('最终模型（'+post_label+'）') if endpoint.get('cleanup_status')=='completed' else '最终模型（固定末轮）','主要结果') if complete else None
    for key,name,role in [('sad','02_sad/experimental.mtz','DM 前 SAD 相位与 HL'),
                         ('sad','02_sad/observations.mtz','SAD 纯观测似然 HL（用于 P+ 与密度修饰）' if cfg.get('refinement_prior')=='dm_phifom' else 'SAD 纯观测似然 HL（固定精修先验）'),
                         ('dm','03_dm/map_phases.mtz','密度修饰后相位与 FOM'),
                         ('bridge','04_bridge/building.mtz','初始 Buccaneer 建模地图')]:
        if key in completed:add(pipe/name,'最终密度修饰相位与 FOM' if key=='dm' and not building and pipeline_complete else 'P+ 更新相位 / FOM / 概率（不含 HL）' if key=='sad' and name.endswith('experimental.mtz') and cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else role,
                              '主要结果' if key=='dm' and not building and pipeline_complete else '相位与地图')
    if complete:add(pipe/'05_build'/('round_%02d_refine'%cfg['rounds'])/'refined.mtz','最后一轮精修地图','主要结果')
    if cfg.get('method') in ('pplus_trust','pplus_bundled'):
        for name in ['resolve.log','classic.log','resolve.inp','seedmap.mtz']:
            add(pipe/'03_dm'/name,name,'日志' if name.endswith('.log') else '相位与地图')
    if cfg.get('auto_hand') and cfg.get('method') in ('pplus_internal','likelihood','pplus'):
        add(pipe/'02_hand/hand_selection.json','双手工作集试建判别记录','运行记录')
        add(pipe/'02_hand/refined_heavy.pdb','选中手的重原子亚结构','相位与地图')
        for hand in ('input_hand','opposite_hand'):
            folder=pipe/'02_hand/hands'/hand
            add(folder/'candidate.json',hand+' 候选记录','运行记录')
            add(folder/'pilot/buccaneer.log',hand+' 工作集试建日志','日志')
    if cfg.get('method')=='pplus_bundled':
        add(pipe/'02_sad/frontend.json','原生 Phaser 前端谱系','运行记录')
        add(pipe/'02_sad/refined_heavy.pdb','工作集精修后的重原子亚结构','相位与地图')
        add(pipe/'02_sad/phaser/phaser.log','原生 Phaser 日志','日志')
        add(pipe/'02_sad/hand_selection.json','双手试建证据与判别记录（不使用 Rfree）','运行记录')
        for hand in ('input_hand','opposite_hand'):
            folder=pipe/'02_sad/hands'/hand
            add(folder/'candidate.json',hand+' 候选记录','运行记录')
            add(folder/'pilot/buccaneer.log',hand+' 工作集试建日志','日志')
            add(folder/'03_dm/map_phases.mtz',hand+' 密度修饰相位','相位与地图')
    if endpoint:
        add(pipe/'06_finish/endpoint.json','末轮与清理结果记录','运行记录')
        add(pipe/'06_finish/refined.mtz',post_label+'地图','主要结果')
        add(pipe/'06_finish/refmac.log',post_label+'日志','日志')
    for name in ['pipeline_report.json','01_inputs/sad.json','02_sad/experimental.json','02_sad/observations.json','03_dm/map_phases.json',
                 '04_bridge/building.audit.json','05_build/hybrid_report.json']:
        add(pipe/name,Path(name).name,'运行记录')
    for name in ['gui.log']:
        add(root/name,'流程日志','日志')
    for name in ['sad.log','dm.log','build.log']:
        add(pipe/name,name,'日志')
    for stage in ['sad','dm','build']:
        for path in sorted(pipe.glob(stage+'.interrupted-*.log')):
            add(path,'保留的中断日志：'+path.name,'日志')
    for number in range(1,(cfg.get('rounds',0) if building else 0)+1):
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
                hand_selection=read_json(pipe/'02_sad/hand_selection.json') or read_json(pipe/'02_hand/hand_selection.json'),
                method_name='P+ Wang8 + 精修权重搜索' if method=='pplus_bundled' and cfg.get('trust_wang_radius')==8 and cfg.get('trust_work_weight') and not cfg.get('trust_stability') else METHODS[method][0],dm_engine=METHODS[method][1],config=cfg,
                stages=stages,sad=raw,dm=dm,dm_history=history,solvent_fraction=dm_meta.get('solvent_fraction_used'),
                solvent_source=dm_meta.get('solvent_source'),partition=partition,metrics=metrics,
                sample=raw_meta.get('sequence_asu',{}),anomalous_scattering=state.get('anomalous_scattering'),
                primary_metric=metrics[-1] if metrics else None,primary_complete=complete,
                build_model=building,pipeline_complete=pipeline_complete,
                final_phases=str((pipe/'03_dm/map_phases.mtz').relative_to(root)) if pipeline_complete and not building else None,
                final_model=str(model.relative_to(root)) if model and model.is_file() else None,
                preview=preview,residues=sum(len(c['points']) for c in preview),chains=len(preview),
                artifacts=artifacts,error=state.get('error'),selection_policy=endpoint.get('selection','fixed_final_round'),cleanup_status=endpoint.get('cleanup_status'),cleanup_reason=endpoint.get('cleanup_reason'),cleanup_kind=endpoint.get('cleanup_kind'))


def save_selected_result(directory):
    result=collect_results(directory)
    if result.get('portfolio'):
        if not result['pipeline_complete']:raise ValueError('多流程探索尚未完整完成；已有候选保留')
        return result
    if not result['pipeline_complete'] or (result['build_model'] and not result['primary_complete']):
        raise ValueError('流程尚未产生完整的最终结果')
    return result


def enrich_job(job):
    from .portfolio_gui import locate,enrich_job as enrich_portfolio
    if job.get('config',{}).get('explore') or locate(job['output']) is not None:
        return enrich_portfolio(job)
    pipe=pipeline_root(job['output']);state=read_json(pipe/'pipeline_report.json')
    if not state:return job
    if state.get('status')=='failed':
        detail=read_json(pipe/'05_build/hybrid_report.json').get('error') if state.get('current')=='build' else None
        job['error']=detail or state.get('error') or job.get('error')
        if detail and 'empty or invalid model' in detail:job['error']='Buccaneer 第 %s 轮没有生成有效模型；已有相位与中间结果保留。请检查波长、所选数据列、序列链划分、重原子模型及其手性。详细原因：'%read_json(pipe/'05_build/hybrid_report.json').get('current_round','?')+detail
    current=state.get('current','');build=read_json(pipe/'05_build/hybrid_report.json')
    round_number=len(build.get('stages',[]))+1
    label=state.get('current_name','')
    logs=[Path(job['output'])/'gui.log']
    if current in ('sad','dm','build','bridge','finish'):logs.append(pipe/(current+'.log'))
    if current=='sad':logs.append(pipe/'02_sad/phaser/phaser.log')
    if current=='dm':logs.extend([pipe/'03_dm/resolve.log',pipe/'03_dm/classic.log'])
    if current=='finish':logs.append(pipe/'06_finish/refmac.log')
    if current=='build' and build.get('status')=='running':
        sub=build.get('current','buccaneer')
        label='第 %d / %d 轮 · %s'%(round_number,state['config']['rounds'],'REFMAC 精修' if sub=='refmac' else 'Buccaneer 构建')
        folder=pipe/'05_build'/('round_%02d_%s'%(round_number,'refine' if sub=='refmac' else 'build'))
        logs+=list(folder.glob('*.log'))
    job['pipeline_progress']=dict(completed=len(state.get('completed',[])),total=(6 if state['config'].get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else 5) if state['config'].get('build_model',True) else 3,current=label,
                                  revision=[state.get('updated'),build.get('current'),len(build.get('stages',[]))])
    parts=[]
    for path in logs:
        if not path.is_file():continue
        with path.open('rb') as handle:
            handle.seek(max(0,path.stat().st_size-24000));tail=handle.read().decode(errors='replace')
        parts.append('── '+str(path.relative_to(Path(job['output'])))+' ──\n'+tail)
    job['log']='\n\n'.join(parts)[-60000:]
    return job
