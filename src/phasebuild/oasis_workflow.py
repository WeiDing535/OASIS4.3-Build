"""Focused OASIS4.3 workflow: explicit phase methods, paired DM, Buccaneer/REFMAC."""
import json
import math
from pathlib import Path
import sys
from .sad_pipeline import DEFAULTS, validate, stage_plan
from .runtime import default_runtime
from .data import fingerprint

METHODS={'likelihood':('continuous SAD','MLHL-best'),'pplus':('legacy direct-method P+','Legacy'),'pplus_trust':('external reference (Phaser SAD frontend)','external Parrot'),'pplus_internal':('internal P+ (experimental)','internal MLHL-best'),'pplus_bundled':('P+ external engine version','external CCP4 Parrot')}


_CLI_CACHE={}


def cli_options(project):
    """Every --flag run_sad_to_model.py actually accepts.

    The GUI turns the whole validated config into argv, so any bookkeeping key
    added to the config (copies_source, dm_solvent_source, ...) would otherwise be
    passed as an unknown flag and abort the run.
    """
    script=Path(project).resolve()/'run_sad_to_model.py'
    if not script.is_file():script=Path(__file__).resolve().parents[2]/'run_sad_to_model.py'
    cached=_CLI_CACHE.get(str(script))
    if cached is not None:return cached
    import importlib.util
    spec=importlib.util.spec_from_file_location('oasis_cli_interface',script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    options=set(module.build_parser()._option_string_actions)
    _CLI_CACHE[str(script)]=options
    return options


def defaults(project):
    project=Path(project)
    # No in-package runtime: discover only a user-pinned or standard host CCP4.
    runtime=default_runtime() or ''
    cfg=dict(DEFAULTS)
    cfg['method']='pplus_bundled'
    # Calibration recipe: external Parrot density modification + 2 feedback cycles + 40 modification cycles.
    cfg.update(dm_engine='parrot',dm_cycles=40,feedback_cycles=2,max_per_target=48,dm_solvent=0.5)
    # Endpoint combination is left empty by default: equivalent to a single endpoint, determined by --feedback-cycles; when filled, the endpoint is selected by work-set R.
    cfg.update(feedback_endpoints='',dm_solvents='')
    cfg.update(sad='',native='',sequence='',heavy='',fpp='',runtime=runtime,
               reference=str(Path(runtime)/'lib/data/reference_structures/reference-1tqw.mtz') if runtime else '',
               output=str(project/'runs/oasis_run'),mode='sad_to_model',
               explore=False,explore_workers=2,protocol='pplus_w8_weight',elements='',auto_fpp=True,input_columns=True,native_labels='',auto_hand=True,building_recipe='pipeline',refinement_prior='dm_phifom',
               trust_wang_radius=8.,trust_partial_cleanup=True,trust_work_weight=True,trust_stability=False)
    return cfg


def make_plan(config, project, python=None, inspect=True):
    if not isinstance(config,dict):raise ValueError('configuration must be a JSON object')
    cfg=defaults(project)
    unknown=set(config)-set(cfg)
    if unknown:raise ValueError('unsupported configuration items in this version: '+', '.join(sorted(unknown)))
    cfg.update(config)
    from .portfolio import PROTOCOLS
    protocol=cfg.pop('protocol') if 'protocol' in config else None
    if protocol is None:cfg.pop('protocol',None)
    elif protocol not in PROTOCOLS:raise ValueError('unknown computation mode: '+str(protocol))
    else:
        if config.get('method','pplus_bundled')!='pplus_bundled':raise ValueError('P+ parameter scheme only applies to external P+; do not include protocol with other phasing methods')
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
    if cfg.get('dm_engine') in ('parrot','handdm') and cfg.get('method') not in ('pplus_trust','pplus_internal','pplus_bundled'):
        # The external DM engines only matter for the P+ routes; likelihood/pplus
        # run their own legacy DM engine, so the choice is nominal and feedback is disabled.
        cfg['feedback_cycles']=0
    if 'explore' not in config:cfg['explore']=False # never silently migrate saved single-run configs
    explore=cfg.pop('explore');workers=cfg.pop('explore_workers')
    if not isinstance(explore,bool):raise ValueError('multi-workflow exploration option must be a boolean')
    if isinstance(workers,bool) or workers not in (1,2):raise ValueError('multi-workflow concurrency can only be 1 or 2')
    if explore and (cfg['method']!='pplus_bundled' or not cfg['build_model']):raise ValueError('multi-workflow exploration requires the P+ external engine version and model building enabled')
    if cfg.pop('mode')!='sad_to_model':raise ValueError('OASIS4.3 only supports SAD → density modification → Buccaneer → REFMAC')
    if not isinstance(cfg['output'],str) or not cfg['output'].strip():raise ValueError('please provide a new output directory')
    outer=Path(cfg['output']).expanduser().resolve()
    if outer.exists():raise ValueError('output directory already exists, please select a new directory')
    cfg['output']=str(outer/'00_sad_to_model')
    checked=validate(cfg,inspect=inspect)
    if explore:
        from .portfolio import plan
        plan(checked,outer/'00_portfolio',inspect=inspect)
        command=[python or sys.executable,'-B','-u',str(Path(project).resolve()/'run_sad_portfolio.py'),
                 '--gui-config',str(outer/'gui_job.json'),'--output',str(outer/'00_portfolio'),'--workers',str(workers)]
        checked.update(output=str(outer),mode='sad_to_model',explore=True,explore_workers=workers)
        if protocol is not None:checked['protocol']=protocol
        return checked,[dict(name='P+ multi-workflow exploration and candidate recommendation',kind='sad_portfolio',output=str(outer/'00_portfolio'),command=command)]
    command=[python or sys.executable,'-u',str(Path(project).resolve()/'run_sad_to_model.py')]
    supported=cli_options(project)
    for key,value in checked.items():
        if '--'+key.replace('_','-') not in supported:continue
        if key=='build_semet':
            if value:command.append('--build-semet')
        elif key in ('trust_cleanup','trust_stability','trust_partial_cleanup','trust_work_weight','auto_fpp','input_columns','auto_hand','canonical_dm_seed','per_round_cleanup','dm_anisotropy'):command.append('--'+key.replace('_','-') if value else '--no-'+key.replace('_','-'))
        elif key=='build_model':command.append('--build-model' if value else '--no-build-model')
        elif key=='rank_fit':
            if value:command.append('--rank-fit')
        elif value is not None:command.extend(['--'+key.replace('_','-'),str(value)])
    checked.update(output=str(outer),mode='sad_to_model',explore=False,explore_workers=workers)
    if protocol is not None:checked['protocol']=protocol
    return checked,[dict(name='OASIS4.3: phasing → density modification'+(' → Buccaneer → REFMAC' if checked['build_model'] else ' (model building not enabled)'),kind='sad_to_model',
                         output=str(outer/'00_sad_to_model'),command=command)]


def _dm_engine_label(cfg):
    """Human-readable label for the density-modification engine in use."""
    return 'external CCP4 Parrot' if cfg.get('dm_engine')!='handdm' else 'handdm (in-house)'


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
        selected_lines=['SAD '+k+' = '+v['label']+' ['+v['dataset']+'; column '+str(v['index'])+']' for k,v in selection['sad_columns'].items()]
        selected_lines += ['observed '+k+' = '+v['label']+' ['+v['dataset']+'; column '+str(v['index'])+']' for k,v in selection['native_columns'].items()]
        if not selection['native_columns']:selected_lines.append('mean observations: generated from the selected SAD Friedel amplitude pairs')
        selected_lines.append('observation file: '+('use the SAD file (no separate observation file selected)' if selection['native_omitted'] else selection['native']['path']))
        selected_lines.append('free set: '+('missing; a reproducible 5% free set will be created; Rfree cannot be directly compared with the original free set' if selection['flags_generated'] else 'keep '+selection['free_column']['label']+', free value '+str(selection['free_value'])))
        result['effective_options']['selected_columns']='\n'.join(selected_lines)
        if selection['flags_generated']:result['sample_warnings'].append('no existing free set provided: a reproducible 5% free set will be generated; Rfree cannot be directly compared with previous refinement results.')
        return result
    if cfg.get('explore'):
        from .portfolio_gui import describe
        return describe(cfg)
    values={**DEFAULTS,**cfg,'output':str(Path(cfg['output'])/'00_sad_to_model')}
    phase,dm=METHODS[cfg['method']]
    active=dict(method=phase,dm_engine=dm,solvent_fraction=cfg['solvent_fraction'],
                solvent_source='estimated from the sequence and number of ASU copies' if cfg['solvent_fraction'] is None else 'user specified',
                copies=cfg['copies'],phase_cycles=cfg['phase_cycles'],dm_cycles=cfg['dm_cycles'],
                rounds=cfg['rounds'],refmac_cycles=cfg['cycles'],seed=cfg['seed'],
                free_value=cfg['free_value'],kmin=cfg['kmin'],max_neighbors=cfg['max_neighbors'],
                max_per_target=cfg['max_per_target'])
    if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled'):
        active.update(trust_work_weight=cfg.get('trust_work_weight',False),trust_stability=cfg.get('trust_stability',False),trust_partial_cleanup=cfg.get('trust_partial_cleanup',False),trust_kl=cfg['trust_kl'],trust_cleanup=cfg['trust_cleanup'],frontend='OASIS internal anomalous-difference likelihood (does not read Phaser files)' if cfg['method']=='pplus_internal' else 'import raw Phaser SAD pure-observation likelihood (does not re-run Phaser)',
                      build_semet=cfg['build_semet'])
    else:active.update({'triplet_weight':cfg['triplet_weight']} if cfg['method']=='likelihood' else {'aoe':cfg['aoe'],'rank_fit':cfg['rank_fit']})
    if cfg['method']=='pplus_bundled':active.update(frontend='external native Phaser: compute SAD likelihood from the work set',wavelength=cfg['wavelength'],engine_origin='third-party CCP4 engine is provided by the user, with the path specified via --ccp4-root')
    if cfg.get('auto_fpp'):
        from .anomalous_scattering import calculate
        scattering=calculate(cfg['elements'],cfg['wavelength'])
        active.update(elements=cfg['elements'],fpp=cfg['fpp'],fpp_source=scattering['source'],
                      fpp_kind='theoretical value computed from the required experimental wavelength, not an absorption-edge measurement')
    if cfg['method']=='pplus_bundled' and cfg.get('trust_wang_radius')==8 and cfg.get('trust_work_weight') and not cfg.get('trust_stability'):
        active.update(method='P+ Wang8 + refinement weight search',dm_engine=_dm_engine_label(cfg))
    if cfg['method'] in ('pplus_trust','pplus_bundled'):active['trust_wang_radius']=cfg.get('trust_wang_radius',14.)
    from .portfolio import PROTOCOLS
    from .portfolio_gui import NAMES
    active_profile=cfg.get('protocol') or next((name for name,options in PROTOCOLS.items() if all(cfg.get(k)==v for k,v in options.items())),None)
    if active_profile and cfg['method']=='pplus_bundled':
        active['protocol']=active_profile
        active['method']=NAMES[active_profile]
        active['dm_engine']=_dm_engine_label(cfg)
    active['per_round_cleanup']=cfg.get('per_round_cleanup',False)
    active['canonical_dm_seed']=cfg.get('canonical_dm_seed',False)
    active['auto_hand']=cfg.get('auto_hand',False)
    active['refinement_prior']=cfg.get('refinement_prior','original_sad')
    active['building_recipe']=cfg.get('building_recipe','legacy')
    if cfg.get('phase_blur_scale') is not None:active['phase_blur']={'SCBLUR':cfg['phase_blur_scale'],'BBLUR':cfg['phase_blur_b'],'policy':'pre-declared DM phase-confidence correction, not selected by Rfree'}
    active['build_model']=cfg.get('build_model',True)
    if cfg.get('feedback_endpoints'):
        active['feedback_endpoints']=cfg['feedback_endpoints']
        if cfg.get('dm_solvents'):active['dm_solvents']=cfg['dm_solvents']
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
    protocol_extra=' P+ optimization additionally performs raw HL conditional-branch integration (1440 points), up to 100 updates; default KL=0.0005. When cleanup is enabled, only a unique complete chain may enter fixed 30 cycles, weight 0.005 post-refinement; when conditions are not met, the original final-cycle result is kept.' if cfg['method'] in ('pplus_trust','pplus_internal','pplus_bundled') else ''
    if cfg.get('auto_hand'):protocol_extra+=' within the same workflow, both hands are computed, each with a three-cycle work-set pilot build, and decided by unique sequencing residues and a clear advantage; when evidence is insufficient, both hands are kept and reported as undetermined, without selecting by Rfree.'
    if cfg.get('refinement_prior')=='dm_phifom':protocol_extra+=' refinement uses frozen work-set density-modification phase/FOM approximate restraints, labeled as DM origin; raw SAD likelihood is saved separately and not presented as the same prior.'
    if cfg.get('per_round_cleanup'):protocol_extra+=' before each refinement cycle, fragments are cleaned according to the fixed work-set DM map and unique sequence anchor; when gating conditions are not met, the original model for that cycle is kept, without deciding by Rfree.'
    if cfg.get('phase_blur_scale') is not None:protocol_extra+=' DM phase restraints use pre-declared SCBLUR=%g, BBLUR=%g; the original phase file is unchanged and phase angles are not rotated.'%(cfg['phase_blur_scale'],cfg['phase_blur_b'])
    if cfg.get('trust_partial_cleanup'):protocol_extra=protocol_extra.replace('only a unique complete chain may enter fixed 30 cycles, ','a unique complete chain or fragments satisfying partial-model gating may enter fixed 30 cycles, ')
    if cfg.get('trust_work_weight'):protocol_extra=protocol_extra.replace('weight 0.005','fixed weight selected by work-set geometry gating (0.005–0.04)')
    if cfg.get('trust_stability'):protocol_extra+=' relationship-stability shrinkage is enabled; this is not independent-observation cross-validation.'
    if cfg['method']=='pplus_internal':protocol_extra+=' internal anomalous-difference likelihood and internal MLHL DM are approximate implementations, not equivalent replacements for Phaser/RESOLVE; the 27.42% from the previous external workflow does not apply to this mode.'
    return dict(stages=stage_plan(values),effective_options=active,sample_warnings=sample['warnings'],
                fixed_protocol=('Buccaneer: 3 cycles in the first round, 2 cycles in each subsequent round; single-threaded; REFMAC '+('frozen work-set DM phase/FOM approximate restraints' if cfg.get('refinement_prior')=='dm_phifom' else 'fixed raw SAD observed HL')+'; final-cycle result fixed.'+protocol_extra) if active['build_model'] else
                'Buccaneer / REFMAC not enabled: the workflow ends after density modification, without checking or calling the model-building refinement engine; phasing results, pure-observation HL, final phases, FOM, and logs are kept separately.')


def read_json(path):
    try:return json.loads(Path(path).read_text())
    except (OSError,ValueError):return {}


def pipeline_root(directory):
    root=Path(directory).resolve()
    return root if (root/'pipeline_report.json').is_file() else root/'00_sad_to_model'


def selected_endpoint(directory):
    """Resolve a fixed-endpoint grid to its retained endpoint directory.

    ``--feedback-endpoints`` runs one full pipeline per predeclared endpoint and
    keeps them all, so the case root only carries ``endpoint_selection.json``.
    Returns ``(endpoint_directory, selection)`` or ``(None, {})``.
    """
    root=Path(directory).resolve()
    for base in (root,root/'00_sad_to_model'):
        selection=read_json(base/'endpoint_selection.json')
        endpoints=selection.get('endpoints') or []
        if not endpoints:continue
        chosen=None
        for item in endpoints:
            if (item.get('feedback')==selection.get('selected_feedback')
                    and item.get('solvent')==selection.get('selected_solvent')):
                chosen=item;break
        if chosen is None:
            chosen=min(endpoints,key=lambda item:(item.get('status')!='completed',
                        (item.get('primary_metric') or {}).get('r_work',float('inf'))))
        path=Path(chosen.get('output') or '')
        if not path.is_absolute():path=base/path
        if not path.is_dir():path=base/Path(str(chosen.get('output') or '')).name
        if path.is_dir():return path,selection
    return None,{}


def _running_endpoint(directory):
    """Newest endpoint directory of an in-progress grid, for progress reporting."""
    root=Path(directory).resolve()
    candidates=[path for path in (root/'00_sad_to_model').glob('endpoint_*') if path.is_dir()]
    if not candidates:return None
    def stamp(path):
        for name in ('pipeline_report.json','gui.log'):
            if (path/name).is_file():return (path/name).stat().st_mtime
        return path.stat().st_mtime
    return max(candidates,key=stamp)


def validate_import(directory):
    from .portfolio_gui import locate
    if locate(directory) is not None:
        from .portfolio_gui import validate_import as validate_portfolio
        return validate_portfolio(directory)
    endpoint,_=selected_endpoint(directory)
    if endpoint is not None:
        return validate_import(endpoint)
    state=read_json(pipeline_root(directory)/'pipeline_report.json')
    if state.get('engine')!='oasis-sad-to-model' or state.get('config',{}).get('method','likelihood') not in METHODS:
        raise ValueError('please select a workflow result directory containing OASIS4.3 pipeline_report.json')


def phase_statistics(path):
    if not path.is_file():return None
    import gemmi,numpy as np
    m=gemmi.read_mtz_file(str(path));a=np.asarray(m.column_with_label('FOM'),dtype=float)
    valid=a[np.isfinite(a)&(a>=0)&(a<=1)]
    return dict(reflections=m.nreflections,mean_fom=float(valid.mean()) if len(valid) else None,
                resolution=m.resolution_high(),spacegroup=m.spacegroup.xhm())


def _engine_label(method,cfg,dm_meta):
 """Report the density-modification engine that actually ran, not the method default."""
 engine=(dm_meta.get('parameters') or {}).get('engine') or cfg.get('dm_engine')
 if engine=='parrot':return 'CCP4 Parrot + %d cycles of trace feedback'%cfg['feedback_cycles'] if cfg.get('feedback_cycles') else 'CCP4 Parrot'
 return METHODS[method][1]


def collect_results(directory):
    root=Path(directory).resolve()
    if not root.is_dir():raise ValueError('result directory does not exist')
    from .portfolio_gui import locate,collect_results as collect_portfolio
    if locate(root) is not None or read_json(root/'gui_job.json').get('config',{}).get('explore'):
        return collect_portfolio(root)
    endpoint,selection=selected_endpoint(root)
    if endpoint is not None:
        result=collect_results(endpoint)
        prefix=endpoint.relative_to(root)
        for item in result.get('artifacts') or []:
            item['path']=str(prefix/item['path'])
        for key in ('final_model','final_phases'):
            if result.get(key):result[key]=str(prefix/result[key])
        result['root']=str(root);result['output']=str(root)
        result['endpoint_root']=str(endpoint)
        result['endpoints']=selection.get('endpoints') or []
        result['endpoint_selection']=selection
        result['selection_policy']=selection.get('rule') or result.get('selection_policy')
        return result
    pipe=pipeline_root(root);state=read_json(pipe/'pipeline_report.json')
    cfg=state.get('config') or read_json(root/'gui_job.json').get('config',{})
    completed={item['id'] for item in state.get('completed',[])}
    prior_label='observation likelihood HL' if any(Path(a['path']).name=='observations.mtz'
        for item in state.get('completed',[]) for a in item.get('artifacts',[])) else 'pre-DM HL (legacy protocol)'
    if cfg.get('refinement_prior')=='dm_phifom':prior_label='DM work-set phase/FOM restraints'
    building=cfg.get('build_model',True)
    report=read_json(pipe/'05_build/hybrid_report.json') if building else {}
    labels=[('prepare','input verification'),('sad','P+ optimization' if cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else 'SAD phasing'),('dm','density modification'),('bridge','model-building map'),('build','Buccaneer / REFMAC')]
    if cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled'):labels.append(('finish','sequence-anchored deduplication / post-refinement' if cfg.get('trust_partial_cleanup') else 'complete-chain cleanup / post-refinement'))
    if not building:labels=labels[:3]
    stages=[]
    for key,label in labels:
        status='completed' if key in completed else state.get('status','waiting') if state.get('current')==key else 'waiting'
        if status=='completed' and key=='build' and state.get('status')!='completed':status='running'
        stages.append(dict(id=key,name=label,status=status))
    metrics=[]
    # Two builder layouts exist: the in-package pipeline writes hybrid_report.json
    # with one entry per Buccaneer/REFMAC cycle, the older route writes one folder
    # per refinement round.
    build_report=read_json(pipe/'05_build/hybrid_report.json') if building else {}
    if build_report.get('status')=='completed' and build_report.get('cycles'):
        for row in build_report['cycles']:
            number=row.get('cycle')
            metrics.append(dict(stage='cycle %s'%number,round=number,r_work=row.get('r_work'),r_free=row.get('r_free'),
                cycles=row.get('refmac_cycles'),free_preserved=build_report.get('free_flags_preserved'),
                engine='Buccaneer + REFMAC',model=str((pipe/'05_build/primary_model.cif').relative_to(root))))
    else:
        for number in range(1,(cfg.get('rounds',0) if building else 0)+1):
            folder=pipe/'05_build'/('round_%02d_refine'%number)
            ref=read_json(folder/'refinement_report.json')
            if ref.get('status')!='success':continue
            stats=ref.get('statistics',{}).get('final',{})
            metrics.append(dict(stage='cycle %d'%number,round=number,r_work=stats.get('r_work'),r_free=stats.get('r_free'),
                cycles=ref.get('requested_cycles'),free_preserved=ref.get('free_flags_preserved'),
                engine='REFMAC · '+prior_label,model=str((folder/'model.cif').relative_to(root))))
    endpoint=read_json(pipe/'06_finish/endpoint.json') if cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else {}
    post_label='partial-model deduplication post-refinement' if endpoint.get('cleanup_kind')=='explicit-single-copy-sequence-anchored-partial-model-reduction' else 'complete-chain post-refinement'
    if endpoint.get('cleanup_status')=='completed' and 'finish' in completed:
        ref=read_json(pipe/'06_finish/refinement_report.json');stats=ref.get('statistics',{}).get('final',{})
        metrics.append(dict(stage=post_label,round=cfg['rounds']+1,r_work=stats.get('r_work'),r_free=stats.get('r_free'),
            cycles=30,free_preserved=ref.get('free_flags_preserved'),engine='REFMAC · amplitude likelihood / fixed weight',model=str((pipe/'06_finish/model.cif').relative_to(root))))
    pipeline_complete=state.get('status')=='completed'
    if pipeline_complete and not building:
        expected=state.get('dm_phases',{}).get('sha256')
        phase_file=pipe/'03_dm/map_phases.mtz'
        if not phase_file.is_file() or not expected or fingerprint(phase_file)['sha256']!=expected:
            raise ValueError('final density-modified phases have a hash mismatch with the workflow record')
    complete=pipeline_complete and building and report.get('status')=='completed'
    model=pipe/'primary_model.cif' if complete else root/metrics[-1]['model'] if metrics else None
    if complete:
        expected=state.get('primary_model',{}).get('sha256')
        if not model.is_file() or not expected or fingerprint(model)['sha256']!=expected:
            raise ValueError('final model has a hash mismatch with the workflow record')
        if len(metrics)!=cfg['rounds']+(endpoint.get('cleanup_status')=='completed'):raise ValueError('final model lacks a complete per-cycle refinement report')
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
    add(pipe/'primary_model.cif',('final model ('+post_label+')') if endpoint.get('cleanup_status')=='completed' else 'final model (fixed final cycle)','primary results') if complete else None
    for key,name,role in [('sad','02_sad/experimental.mtz','pre-DM SAD phases and HL'),
                         ('sad','02_sad/observations.mtz','SAD pure-observation likelihood HL (for P+ and density modification)' if cfg.get('refinement_prior')=='dm_phifom' else 'SAD pure-observation likelihood HL (fixed refinement prior)'),
                         ('dm','03_dm/map_phases.mtz','post-density-modification phases and FOM'),
                         ('bridge','04_bridge/building.mtz','initial Buccaneer model-building map')]:
        if key in completed:add(pipe/name,'final density-modified phases and FOM' if key=='dm' and not building and pipeline_complete else 'P+-updated phases / FOM / probabilities (without HL)' if key=='sad' and name.endswith('experimental.mtz') and cfg.get('method') in ('pplus_trust','pplus_internal','pplus_bundled') else role,
                              'primary results' if key=='dm' and not building and pipeline_complete else 'phases and maps')
    if complete:
        last=pipe/'05_build'/(('cycle_%02d'%cfg['rounds']) if build_report.get('status')=='completed' else ('round_%02d_refine'%cfg['rounds']))
        add(last/(('refine.mtz') if build_report.get('status')=='completed' else 'refined.mtz'),'final-cycle refinement map','primary results')
    if cfg.get('method') in ('pplus_trust','pplus_bundled'):
        for name in ['resolve.log','classic.log','resolve.inp','seedmap.mtz']:
            add(pipe/'03_dm'/name,name,'log' if name.endswith('.log') else 'phases and maps')
    if cfg.get('auto_hand') and cfg.get('method') in ('pplus_internal','likelihood','pplus'):
        add(pipe/'02_hand/hand_selection.json','two-hand work-set pilot-build discrimination record','run record')
        add(pipe/'02_hand/refined_heavy.pdb','heavy-atom substructure of the selected hand','phases and maps')
        for hand in ('input_hand','opposite_hand'):
            folder=pipe/'02_hand/hands'/hand
            add(folder/'candidate.json',hand+' candidate record','run record')
            add(folder/'pilot/buccaneer.log',hand+' work-set pilot-build log','log')
    if cfg.get('method')=='pplus_bundled':
        add(pipe/'02_sad/frontend.json','native Phaser frontend lineage','run record')
        add(pipe/'02_sad/refined_heavy.pdb','work-set refined heavy-atom substructure','phases and maps')
        add(pipe/'02_sad/phaser/phaser.log','native Phaser log','log')
        add(pipe/'02_sad/hand_selection.json','two-hand pilot-build evidence and discrimination record (Rfree not used)','run record')
        for hand in ('input_hand','opposite_hand'):
            folder=pipe/'02_sad/hands'/hand
            add(folder/'candidate.json',hand+' candidate record','run record')
            add(folder/'pilot/buccaneer.log',hand+' work-set pilot-build log','log')
            add(folder/'03_dm/map_phases.mtz',hand+' density-modified phases','phases and maps')
    if endpoint:
        add(pipe/'06_finish/endpoint.json','final-cycle and cleanup result record','run record')
        add(pipe/'06_finish/refined.mtz',post_label+'map','primary results')
        add(pipe/'06_finish/refmac.log',post_label+'log','log')
    for name in ['pipeline_report.json','01_inputs/sad.json','02_sad/experimental.json','02_sad/observations.json','03_dm/map_phases.json',
                 '04_bridge/building.audit.json','05_build/hybrid_report.json']:
        add(pipe/name,Path(name).name,'run record')
    for name in ['gui.log']:
        add(root/name,'workflow log','log')
    for name in ['sad.log','dm.log','build.log']:
        add(pipe/name,name,'log')
    for stage in ['sad','dm','build']:
        for path in sorted(pipe.glob(stage+'.interrupted-*.log')):
            add(path,'retained interrupted log: '+path.name,'log')
    for number in range(1,(cfg.get('rounds',0) if building else 0)+1):
        for suffix,file,role in [('build','buccaneer.log','Buccaneer log'),('refine','refmac.log','REFMAC log'),
                                 ('refine','refinement_report.json','refinement statistics'),('build','buccaneer.xml','Buccaneer statistics')]:
            add(pipe/'05_build'/('round_%02d_%s'%(number,suffix))/file,'cycle %d '%number+role,'per-cycle record')
    raw=phase_statistics(pipe/'02_sad/experimental.mtz') if 'sad' in completed else None
    dm=phase_statistics(pipe/'03_dm/map_phases.mtz') if 'dm' in completed else None
    dm_meta=read_json(pipe/'03_dm/map_phases.json') if 'dm' in completed else {}
    raw_meta=read_json(pipe/'02_sad/experimental.json') if 'sad' in completed else {}
    history=dm_meta.get('density_modification',{}).get('cycles',[])
    partition=read_json(pipe/'01_inputs/sad.json').get('native_partition',{})
    method=cfg.get('method','likelihood')
    return dict(root=str(root),status=state.get('status','waiting'),method=method,
                hand_selection=read_json(pipe/'02_sad/hand_selection.json') or read_json(pipe/'02_hand/hand_selection.json'),
                method_name='P+ Wang8 + refinement weight search' if method=='pplus_bundled' and cfg.get('trust_wang_radius')==8 and cfg.get('trust_work_weight') and not cfg.get('trust_stability') else METHODS[method][0],dm_engine=_engine_label(method,cfg,dm_meta),config=cfg,
                stages=stages,sad=raw,dm=dm,dm_history=history,solvent_fraction=dm_meta.get('solvent_fraction_used',(dm_meta.get('parameters') or {}).get('solvent_fraction')),
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
        if not result['pipeline_complete']:raise ValueError('multi-workflow exploration has not completed fully; existing candidates are kept')
        return result
    if not result['pipeline_complete'] or (result['build_model'] and not result['primary_complete']):
        raise ValueError('workflow has not yet produced complete final results')
    return result


def enrich_job(job):
    from .portfolio_gui import locate,enrich_job as enrich_portfolio
    if job.get('config',{}).get('explore') or locate(job['output']) is not None:
        return enrich_portfolio(job)
    pipe=pipeline_root(job['output']);state=read_json(pipe/'pipeline_report.json')
    if not state:
        # A fixed-endpoint grid has one pipeline per endpoint; report the
        # endpoint that is running (or the retained one) instead of nothing.
        endpoint,_=selected_endpoint(job['output'])
        running=endpoint or _running_endpoint(job['output'])
        if running is None:return job
        sub=dict(job);sub['output']=str(running)
        enriched=enrich_job(sub)
        enriched['output']=job['output'];enriched['endpoint']=running.name
        progress=enriched.setdefault('pipeline_progress',{})
        label=progress.get('current') or ''
        progress['current']=('endpoint %s · %s'%(running.name,label) if label else 'endpoint '+running.name)
        top=Path(job['output'])/'gui.log'
        if top.is_file():
            with top.open('rb') as handle:
                handle.seek(max(0,top.stat().st_size-24000))
                enriched['log']=(enriched.get('log') or '')+'\n\n── gui.log ──\n'+handle.read().decode(errors='replace')
        return enriched
    if state.get('status')=='failed':
        detail=read_json(pipe/'05_build/hybrid_report.json').get('error') if state.get('current')=='build' else None
        job['error']=detail or state.get('error') or job.get('error')
        if detail and 'empty or invalid model' in detail:job['error']='Buccaneer produced no valid model in cycle %s; existing phases and intermediate results are kept. Please check the wavelength, selected data columns, sequence chain partitioning, heavy-atom model, and its chirality. Details: '%read_json(pipe/'05_build/hybrid_report.json').get('current_round','?')+detail
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
        label='cycle %d / %d · %s'%(round_number,state['config']['rounds'],'REFMAC refinement' if sub=='refmac' else 'Buccaneer building')
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
