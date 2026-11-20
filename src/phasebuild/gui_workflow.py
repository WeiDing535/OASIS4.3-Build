"""Validated GUI plans; subprocess argument lists share the existing CLI engines."""
import json
import math
import re
import sys
from pathlib import Path
from .runtime import default_runtime, resolve_runtime

INTEGRATED_MODES = ('hybrid_hl', 'build_hl', 'sad_to_model')
INTEGRATED_ENGINES = ('hybrid-phasebuild-seeded-external-buccaneer-refmac',
                      'phasebuild-integrated-external-buccaneer-refmac')


def defaults(project):
    project = Path(project)
    return dict(sad_input='', heavy='', fpp='', copies='auto', phase_cycles=5, dm_cycles=5, sad_labels='',
                dm_reference=str(Path(default_runtime())/'lib/data/reference_structures/reference-1tqw.mtz') if default_runtime() else '',
                mode='fresh', map='', native='', sequence='', model='', sequence_report='',
                phase_prior='', build_semet=False,
                output=str(project/'runs'/'gui_run'), profile=str(project/'models/1tqw_1p9_symmetric.npz'),
                residue_profile=str(project/'models/residue_1tqw_1p9.npz'),
                rotamers=str(project/'models/rotamers_1tqw_v08.json'),
                ccp4=default_runtime(), refine=True, rounds=5, cycles=10,
                map_mode='rebuild', map_amplitude='FWT',map_phase='PHWT',rotations=512, refine_count=96, max_length=30,
                max_residues=65, beam_width=4, selection_cutoff=2.5, min_score=2.0,
                amplitude='F', sigma='SIGF', free_column='FreeR_flag', free_value=0)


def _inspect_hybrid_inputs(cfg):
    """Validate the fixed hybrid schema in memory, without creating an MTZ."""
    import gemmi
    import numpy as np
    from .refinement import reflection_partition
    from .sequencing import read_sequences
    native, _ = reflection_partition(cfg['native'], free_value=cfg['free_value'])
    maps = gemmi.read_mtz_file(cfg['map'])
    prior = gemmi.read_mtz_file(cfg['phase_prior'])
    required = [('map',maps,{'F':'F','SIGF':'Q','FWT':'F','PHWT':'P','FOM':'W'}),
                ('original experimental phases',prior,{label:'A' for label in ('HLA','HLB','HLC','HLD')})]
    for title, mtz, columns in required:
        if (mtz.spacegroup.hm != native.spacegroup.hm or
                not np.allclose(mtz.cell.parameters,native.cell.parameters,rtol=0.,atol=.01)):
            raise ValueError(title+' has a unit cell or space group setting inconsistent with native')
        for label, kind in columns.items():
            if mtz.column_labels().count(label) != 1 or mtz.column_with_label(label).type != kind:
                raise ValueError('{} must contain a unique {} column ({} type)'.format(title,label,kind))
    if any(label in native.column_labels() for label in ('SAD_HLA','SAD_HLB','SAD_HLC','SAD_HLD')):
        raise ValueError('native already contains SAD_HL prior columns; please provide raw native data without merged phases')
    native_data, prior_data = np.array(native,copy=False), np.array(prior,copy=False)
    for title, values in [('native',native_data),('original experimental phases',prior_data)]:
        if not np.isfinite(values[:,:3]).all() or not np.all(values[:,:3] == np.rint(values[:,:3])):
            raise ValueError(title+' contains invalid integer HKL')
        if len(set(map(tuple,values[:,:3]))) != len(values):
            raise ValueError(title+' contains duplicate HKL')
    native_hkl = set(map(tuple,native_data[:,:3]))
    hl_columns = [prior.column_labels().index(label) for label in ('HLA','HLB','HLC','HLD')]
    if not any(tuple(row[:3]) in native_hkl and np.isfinite(row[hl_columns]).all() for row in prior_data):
        raise ValueError('original experimental HL and native have no exact HKL intersection with finite values')
    if cfg['mode'] == 'hybrid_hl':
        structure = gemmi.read_structure(cfg['model'])
        group = structure.find_spacegroup()
        if (len(structure) != 1 or not any(r.find_atom('CA','*') for c in structure[0] for r in c)):
            raise ValueError('input seed model must contain one model and Cα atoms')
        if (group is None or group.hm != native.spacegroup.hm or
                not np.allclose(structure.cell.parameters,native.cell.parameters,rtol=0.,atol=.01)):
            raise ValueError('seed model has a unit cell or space group setting inconsistent with native')
    if not read_sequences(cfg['sequence']):
        raise ValueError('sequence file is empty')



def make_sad_plan(cfg, project, python=None, inspect=True):
    from .sad_pipeline import validate
    out = Path(cfg['output']).expanduser().resolve()
    if out.exists():
        raise ValueError('output directory already exists, please choose a new name: '+str(out))
    values = dict(sad=cfg['sad_input'], native=cfg['native'], sequence=cfg['sequence'],
                  heavy=cfg['heavy'], reference=cfg['dm_reference'], fpp=cfg['fpp'],
                  copies=cfg['copies'], phase_cycles=cfg['phase_cycles'], dm_cycles=cfg['dm_cycles'],
                  rounds=cfg['rounds'], cycles=cfg['cycles'], free_value=cfg['free_value'],
                  labels=cfg['sad_labels'], runtime=cfg['ccp4'], output=str(out/'00_sad_to_model'),
                  build_semet=False)
    checked=validate(values, inspect=inspect)
    cfg.update(output=str(out),refine=True,model='',map='',phase_prior='',ccp4=checked['runtime'])
    command=[python or sys.executable,'-u',str(project/'run_sad_to_model.py')]
    for key,value in checked.items():
        if key=='build_model':command.append('--build-model' if value else '--no-build-model')
        elif key in ('build_semet','rank_fit'):
            if value:command.append('--'+key.replace('_','-'))
        elif key in ('dm_anisotropy',):
            command.append('--'+key.replace('_','-') if value else '--no-'+key.replace('_','-'))
        elif value is not None:command.extend(['--'+key.replace('_','-'),str(value)])
    return cfg,[dict(name='SAD → density modification → model building and refinement (resumable)',kind='sad_to_model',
                     output=checked['output'],command=command)]


def make_plan(config, project, python=None, inspect=True):
    """Read/validate inputs before any output is created. Never overwrite a run."""
    project = Path(project).resolve()
    cfg = defaults(project)
    if not isinstance(config, dict):
        raise ValueError('configuration must be a JSON object')
    if set(config)-set(cfg):
        raise ValueError('unknown configuration item: '+', '.join(sorted(set(config)-set(cfg))))
    cfg.update(config)
    if cfg['mode'] == 'sad_to_model':
        return make_sad_plan(cfg, project, python=python, inspect=inspect)
    if cfg['mode'] not in ('fresh', 'extend', 'refine','evaluate') + INTEGRATED_MODES:
        raise ValueError('please select a valid job type')
    if not isinstance(cfg['refine'], bool):
        raise ValueError('refinement switch must be a boolean')
    if not isinstance(cfg['build_semet'], bool):
        raise ValueError('SeMet switch must be a boolean')
    if cfg['map_mode'] not in ('rebuild', 'updated', 'fixed'):
        raise ValueError('please select a valid map strategy')
    for key, low, high in [('rounds',1,50),('cycles',1,200),('rotations',8,8192),
                           ('refine_count',1,2048),('max_length',3,1000),
                           ('max_residues',3,20000),('beam_width',1,32),('free_value',0,999)]:
        try:
            value = float(cfg[key])
            if not math.isfinite(value) or not value.is_integer() or not low <= value <= high:
                raise ValueError()
            cfg[key] = int(value)
        except (TypeError, ValueError):
            raise ValueError('{} must be an integer between {} and {}'.format(key,low,high))
    for key in ('selection_cutoff','min_score'):
        try:
            cfg[key] = float(cfg[key])
            if not math.isfinite(cfg[key]) or not 0 < cfg[key] <= 20:
                raise ValueError()
        except (TypeError, ValueError):
            raise ValueError(key+' must be between 0 and 20')
    hybrid = cfg['mode'] in INTEGRATED_MODES
    if hybrid and cfg['rounds'] > 20:
        raise ValueError('number of build cycles for Buccaneer integration mode must be between 1 and 20')
    if hybrid:
        cfg.update(refine=True,map_amplitude='FWT',map_phase='PHWT',amplitude='F',
                   sigma='SIGF',free_column='FreeR_flag',map_mode='updated')
        if cfg['mode'] == 'build_hl':
            cfg['model'] = ''
    refining = hybrid or cfg['mode']=='refine' or (cfg['mode']!='evaluate' and cfg['refine'])
    evaluating = cfg['mode']=='evaluate'
    if not refining and not evaluating and cfg['rounds'] != 1:
        raise ValueError('only one build is supported when refinement is disabled; please set the number of cycles to 1')
    required = ['map','native','sequence','phase_prior'] if hybrid else ['native','model'] if cfg['mode'] in ('refine','evaluate') else [
        'map','sequence','profile','residue_profile','rotamers']
    if cfg['mode'] == 'hybrid_hl':
        required.append('model')
    if cfg['mode']=='extend':
        required += ['model','sequence_report']
    if refining:
        required += ['native']
        if hybrid:
            runtime = resolve_runtime(cfg['ccp4'] or None, required=('cbuccaneer','refmac5'))
            cfg['ccp4'] = str(runtime)
            for suffix in ('.mtz','.pdb'):
                library = runtime/'lib/data/reference_structures'/('reference-1tqw'+suffix)
                if not library.is_file():
                    raise ValueError('engine training resource does not exist: '+str(library))
        else:
            binary = Path(str(cfg['ccp4'])).expanduser()/'bin/refmac5'
            if not binary.is_file():
                raise ValueError('could not find REFMAC: '+str(binary))
    for key in set(required):
        if not isinstance(cfg[key],str) or not cfg[key].strip():
            raise ValueError('please provide input: '+key)
        path = Path(cfg[key]).expanduser().resolve()
        if not path.is_file():
            raise ValueError('file does not exist: '+str(path))
        cfg[key] = str(path)
    if not isinstance(cfg['output'],str) or not cfg['output'].strip():
        raise ValueError('please provide a new output directory')
    out = Path(cfg['output']).expanduser().resolve()
    if out.exists():
        raise ValueError('output directory already exists, please choose a new name: '+str(out))
    parent = out.parent
    while not parent.exists():
        parent = parent.parent
    if not parent.is_dir():
        raise ValueError('parent path of the output directory is not a folder')
    cfg['output'] = str(out)
    if inspect:
        from .data import Density
        from .profile import DensityProfile
        from .residue_profile import ResidueProfile
        from .rotamers import load_library
        from .sequencing import read_sequences
        import gemmi
        if hybrid:
            _inspect_hybrid_inputs(cfg)
        elif cfg['mode'] not in ('refine','evaluate'):
            density = Density(cfg['map'],cfg['map_amplitude'],cfg['map_phase'])
            if not read_sequences(cfg['sequence']):
                raise ValueError('sequence file is empty')
            DensityProfile.load(cfg['profile'])
            ResidueProfile.load(cfg['residue_profile'])
            library = load_library(cfg['rotamers'])
            if any(not item['component']['bonds'] for item in library['entries'].values()):
                raise ValueError('rotamer library contains an empty bond graph; please use the repaired v08 chemical library')
        if cfg['mode'] in ('extend','refine','evaluate'):
            model = gemmi.read_structure(cfg['model'])
            if not len(model) or not any(r.find_atom('CA','*') for c in model[0] for r in c):
                raise ValueError('input model has no Cα atoms')
            if cfg['mode']=='extend':
                from .data import fingerprint
                from .loops import sequence_bindings
                previous = json.loads(Path(cfg['sequence_report']).read_text())
                if previous['sequence_input']['sha256'] != fingerprint(cfg['sequence'])['sha256']:
                    raise ValueError('sequence assignment report is inconsistent with the selected sequence')
                bindings = sequence_bindings(model,previous['sequence_assignment'],read_sequences(cfg['sequence']))
                if not any(binding is not None for binding in bindings):
                    raise ValueError('model has no sequence anchor available for extension')
                if not model.cell.is_similar(density.cell,.001,.1) or gemmi.find_spacegroup_by_name(model.spacegroup_hm).number != density.sg.number:
                    raise ValueError('model and phase map have inconsistent unit cells or space groups')
        if (refining or evaluating) and not hybrid:
            from .refinement import reflection_partition
            native, _ = reflection_partition(cfg['native'],cfg['amplitude'],cfg['sigma'],cfg['free_column'],cfg['free_value'])
            comparison = model if cfg['mode'] in ('refine','evaluate') else density.mtz
            comparison_sg = gemmi.find_spacegroup_by_name(comparison.spacegroup_hm) if cfg['mode'] in ('refine','evaluate') else comparison.spacegroup
            if comparison_sg is None or comparison_sg.number != native.spacegroup.number or not comparison.cell.is_similar(native.cell,.001,.1):
                raise ValueError('native data and model/phase map have inconsistent unit cells or space groups')
    steps = []
    python = python or sys.executable

    def add(name, script, args, directory, kind):
        steps.append(dict(name=name,command=[python,'-u',str(project/script)]+[str(v) for v in args],
                          output=str(directory),kind=kind))

    def refine(model, directory, number):
        add('REFMAC refinement cycle {}'.format(number),'run_refine.py',[
            cfg['native'],'--model',model,'--output',directory,'--cycles',cfg['cycles'],
            '--ccp4',cfg['ccp4'],'--amplitude',cfg['amplitude'],'--sigma',cfg['sigma'],
            '--free-column',cfg['free_column'],'--free-value',cfg['free_value']],directory,'refine')
        return directory/'model.cif', directory/'refined.mtz'

    common = ['--sequence',cfg['sequence'],'--profile',cfg['profile'],
              '--residue-profile',cfg['residue_profile'],'--rotamers',cfg['rotamers']]
    if hybrid:
        unseeded = cfg['mode'] == 'build_hl'
        directory = out/('01_build_hl' if unseeded else '01_hybrid_hl')
        args = ['--map',cfg['map'],'--native',cfg['native'],
                '--sequence',cfg['sequence'],'--phase-prior',cfg['phase_prior'],
                '--runtime',cfg['ccp4'],'--output',directory,
                '--rounds',cfg['rounds'],'--cycles',cfg['cycles'],
                '--free-value',cfg['free_value']]
        if not unseeded:
            args += ['--model',cfg['model']]
        if cfg['build_semet']:
            args.append('--build-semet')
        add(('no initial model: ' if unseeded else 'seeded model: ') +
            'Buccaneer building + original HL restraint refinement (fixed final cycle)',
            'run_build_from_phases.py' if unseeded else 'run_phasebuild_hybrid_hl.py',
            args,directory,cfg['mode'])
    elif evaluating:
        directory=out/'01_evaluate'
        add('compute Rwork/Rfree independently (coordinates not refined)','run_rfactor.py',[cfg['native'],'--model',cfg['model'],
            '--output',directory,'--amplitude',cfg['amplitude'],'--sigma',cfg['sigma'],
            '--free-column',cfg['free_column'],'--free-value',cfg['free_value']],directory,'evaluate')
    elif cfg['mode']=='refine':
        refine(cfg['model'],out/'01_refine',1)
    else:
        current_map = cfg['map']
        current_amplitude,current_phase=cfg['map_amplitude'],cfg['map_phase']
        def map_args():
            return ['--map-amplitude',current_amplitude,'--map-phase',current_phase]
        if cfg['mode']=='fresh':
            trace = out/'01_trace'
            add('whole-map backbone search','run_phasebuild.py',['trace',current_map,'--amplitude',current_amplitude,'--phase',current_phase,'--profile',cfg['profile'],
                '--output',trace,'--rotations',cfg['rotations'],'--refine-count',cfg['refine_count'],
                '--max-length',cfg['max_length'],'--max-residues',cfg['max_residues'],
                '--beam-width',cfg['beam_width'],'--selection-cutoff',cfg['selection_cutoff'],
                '--min-score',cfg['min_score']],trace,'trace')
            improve = out/'02_sequence'
            add('fragment completion and sequence assignment','run_improve.py',[current_map,'--model',trace/'backbone.cif',
                '--profile',cfg['profile'],'--candidates',trace/'candidates.npz','--sequence',cfg['sequence'],
                '--residue-profile',cfg['residue_profile']]+([] if cfg['map_mode']=='rebuild' else ['--polish-terminals'])+map_args()+['--output',improve],improve,'build')
            complete = out/'03_sidechains'
            add('loop completion and side-chain building','run_complete.py',[current_map,'--model',improve/'model.cif',
                '--sequence-report',improve/'build_report.json']+common+map_args()+['--output',complete],complete,'build')
            model, report = complete/'model.cif', complete/'build_report.json'
            if refining:
                model, updated = refine(model,out/'round_01_refine',1)
                if cfg['map_mode'] in ('updated','rebuild'):
                    current_map = updated
                    current_amplitude,current_phase='FWT','PHWT'
            first_extension = 2
        else:
            model, report = Path(cfg['model']),Path(cfg['sequence_report'])
            first_extension = 1
        for number in range(first_extension,cfg['rounds']+1):
            # Reproduce the validated schedule: refine -> original-map rebuild
            # -> refine -> original-map extension -> updated-map extensions.
            rebuilding = cfg['mode']=='fresh' and cfg['map_mode']=='rebuild' and number==2
            if cfg['map_mode']=='rebuild':
                if rebuilding or number==3 and cfg['mode']=='fresh':
                    current_map = cfg['map']
                    current_amplitude,current_phase=cfg['map_amplitude'],cfg['map_phase']
            extension = out/('round_{:02d}_{}'.format(number,'rebuild' if rebuilding else 'extend'))
            add(('cycle {} cleanup, loop completion, and sequence reassignment after refinement' if rebuilding else 'cycle {} bidirectional terminal extension').format(number),
                'run_rebuild.py' if rebuilding else 'run_extend.py',[
                current_map,'--model',model,'--sequence-report',report]+common+map_args()+['--output',extension],extension,'build')
            model, report = extension/'model.cif',extension/'build_report.json'
            if refining:
                model, updated = refine(model,out/('round_{:02d}_refine'.format(number)),number)
                if cfg['map_mode'] in ('updated','rebuild'):
                    current_map = updated
                    current_amplitude,current_phase='FWT','PHWT'
    return cfg, steps


def collect_results(directory):
    """Read saved reports; never infer or fabricate R factors from coordinates."""
    root = Path(directory).resolve()
    if not root.is_dir():
        raise ValueError('result directory does not exist')
    rows, artifacts, hybrid_reports = [], [], []
    for path in sorted(root.rglob('*report.json')):
        if any(part in ('source_snapshot','runner_snapshot') for part in path.relative_to(root).parts):
            continue
        if not path.is_file() or path.stat().st_size > 20_000_000:
            continue
        if len(artifacts) >= 1000:
            break
        artifacts.append(str(path.relative_to(root)))
        if path.name=='hybrid_report.json':
            hybrid_data = json.loads(path.read_text())
            if hybrid_data.get('engine') in INTEGRATED_ENGINES:
                hybrid_reports.append((path,hybrid_data))
        if path.name not in ('refinement_report.json','r_factor_report.json'):
            continue
        data = json.loads(path.read_text())
        if data.get('status')!='success':
            continue
        final = data.get('statistics',{}).get('final',{})
        rows.append(dict(stage=str(path.parent.relative_to(root)),r_work=final.get('r_work'),
                         r_free=final.get('r_free'),cycles=data.get('requested_cycles'),engine=data.get('engine','external-REFMAC5'),
                         free_preserved=data.get('free_flags_preserved'),
                         model=str((path.parent/'model.cif').relative_to(root))))
    # Plans retain chronological order, independent of directory naming conventions.
    manifest = root/'gui_job.json'
    data = json.loads(manifest.read_text()) if manifest.is_file() else {}
    order = {str(Path(s['output']).resolve().relative_to(root)):i for i,s in enumerate(data.get('steps',[]))}
    def stage_order(row):
        stage=row['stage']
        if order:
            return (0,order.get(stage,len(order)),stage)
        # Legacy experiments nest cycles_40/50 under continuation/, alongside
        # cycles_30. Lexical directory order would display cycle 30 as final.
        numbered=re.search(r'(?:cycles_|round_)(\d+)',stage)
        if numbered:
            return (0,int(numbered.group(1)),stage)
        if stage=='.':
            return (0,-1,stage)
        return (1,max(p.stat().st_mtime for p in (root/stage).glob('*report.json')),stage)
    rows.sort(key=stage_order)
    if manifest.is_file():
        artifacts.append('gui_job.json')
    if (root/'selection.json').is_file():
        artifacts.append('selection.json')
    for suffix in ('*.cif','*.pdb','*.mtz','*.ccp4','*.log','*.txt','*.csv'):
        for path in sorted(root.rglob(suffix)):
            if path.is_file() and len(artifacts)<1000 and 'source_snapshot' not in path.parts:
                artifacts.append(str(path.relative_to(root)))
    for state in root.rglob('pipeline_report.json'):
        for relative in ['02_sad/experimental.json','03_dm/map_phases.json','04_bridge/building.audit.json']:
            item=state.parent/relative
            if item.is_file():artifacts.append(str(item.relative_to(root)))
    candidates = []
    for step in data.get('steps',[]):
        folder = Path(step['output'])
        for name in ('model.cif','backbone.cif'):
            if (folder/name).is_file():
                candidates.append(folder/name)
                break
    if not candidates:
        candidates = sorted(root.rglob('model.cif'),key=lambda p:p.stat().st_mtime)
    fixed_endpoint = bool(hybrid_reports) or data.get('config',{}).get('mode') in INTEGRATED_MODES
    refmac = [r for r in rows if not fixed_endpoint and r['engine']=='external-REFMAC5' and
              isinstance(r['r_free'],(int,float)) and math.isfinite(r['r_free']) and
              r['free_preserved'] is True and (root/r['model']).is_file()]
    best = min(refmac,key=lambda r:r['r_free']) if refmac else None
    last_model = str(candidates[-1].relative_to(root)) if candidates else None
    primary_metric, primary_complete = None, False
    if fixed_endpoint:
        # Hybrid protocols predeclare the last round. The lowest Rfree is never
        # used, including while viewing an interrupted or relocated experiment.
        candidates = [root/rows[-1]['model']] if rows and (root/rows[-1]['model']).is_file() else candidates
        primary_metric = rows[-1] if rows else None
        for report_path, report in sorted(hybrid_reports,key=lambda item:item[0].stat().st_mtime):
            primary = report_path.parent/'primary_model.cif'
            if report.get('status')!='completed' or not primary.is_file():
                continue
            expected = report.get('primary_model',{}).get('sha256')
            if expected:
                from .data import fingerprint
                if fingerprint(primary)['sha256'] != expected:
                    raise ValueError('final model (fixed final cycle) and hybrid_report.json hash mismatch')
            candidates = [primary]
            last_model = str(primary.relative_to(root))
            final_round = report['stages'][-1]['round']
            final_stage = str((report_path.parent/('round_%02d_refine'%final_round)).relative_to(root))
            primary_metric = next((row for row in rows if row['stage']==final_stage),None)
            if primary_metric is None:
                raise ValueError('final model (fixed final cycle) is missing the corresponding refinement report')
            primary_complete = True
    if best:
        candidates.append(root/best['model'])
    preview = []
    if candidates:
        import gemmi
        model = gemmi.read_structure(str(candidates[-1]))
        if len(model):
            for chain in model[0]:
                points = []
                for residue in chain:
                    ca = residue.find_atom('CA','*')
                    if ca:
                        points.append([ca.pos.x,ca.pos.y,ca.pos.z,str(residue.seqid),residue.name])
                if points:
                    preview.append(dict(chain=chain.name,points=points))
    return dict(root=str(root),best_metric=best,primary_metric=primary_metric,
                selection_policy='fixed_final_round' if fixed_endpoint else 'exploratory_lowest_refmac_rfree',
                primary_complete=primary_complete,last_model=last_model,metrics=rows,artifacts=sorted(set(artifacts)),preview=preview,
                final_model=str(candidates[-1].relative_to(root)) if candidates else None,
                residues=sum(len(c['points']) for c in preview))


def save_selected_result(directory):
    """Materialize a reproducible selection; retain every round and the last model."""
    import shutil
    root = Path(directory)
    results = collect_results(root)
    if results['selection_policy']=='fixed_final_round':
        if not results['primary_complete']:
            raise ValueError('integrated workflow did not generate a complete final model (fixed final cycle); check hybrid_report.json and the log')
        selection = dict(criterion='fixed_final_round',selected=results['primary_metric'],
                         last_model=results['last_model'],model=results['final_model'],
                         scope='Predeclared final round; no Rfree-based model selection. External Buccaneer performs construction.')
        (root/'selection.json').write_text(json.dumps(selection,indent=2,allow_nan=False)+'\n')
        return results
    best = results['best_metric']
    if best:
        shutil.copyfile(str(root/best['model']),str(root/'selected_model.cif'))
        selection = dict(criterion='lowest_REFMAC_Rfree_with_preserved_flags',
                         selected=best,last_model=results['last_model'],
                         model='selected_model.cif',
                         scope='Exploratory selection among a fixed planned set of rounds; '
                               'Rfree was used for selection, not an untouched final test.')
        (root/'selection.json').write_text(json.dumps(selection,indent=2,allow_nan=False)+'\n')
    return results
