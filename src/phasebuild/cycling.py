"""Fixed-budget experiments: extra refinement versus repeated native rebuilding."""
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
import gemmi
from .chemistry import BACKBONE
from .data import fingerprint
from .refinement import refmac_refine
from .sequencing import read_sequences
from .loops import sequence_bindings


def validate_start(native,sequence,start,sequence_report):
    start=Path(start).resolve()
    report=json.loads((start/'refinement_report.json').read_text())
    if report.get('status')!='success' or not report.get('free_flags_preserved'):
        raise ValueError('A successful audited starting refinement is required')
    protocol=report['reflections']
    if protocol['input']!=fingerprint(native) or any(protocol.get(k)!=v for k,v in
            [('fp','F'),('sigfp','SIGF'),('free_column','FreeR_flag'),('free_value',0)]):
        raise ValueError('Starting refinement uses a different observation/free-set protocol')
    if report['output_model']!=fingerprint(start/'model.cif') or report['output_mtz']!=fingerprint(start/'refined.mtz'):
        raise ValueError('Starting output model or map has changed')
    assignment=json.loads(Path(sequence_report).read_text())
    if assignment['sequence_input']['sha256']!=fingerprint(sequence)['sha256']:
        raise ValueError('Sequence report belongs to another sequence')
    sequence_bindings(gemmi.read_structure(str(start/'model.cif')),assignment['sequence_assignment'],read_sequences(sequence))
    return report


def model_counts(path):
    structure=gemmi.read_structure(str(path))
    return {'residues':sum(len(c) for c in structure[0]),'fragments':len(structure[0]),
            'heavy_atoms':sum(a.element.name not in ['H','D'] for c in structure[0] for r in c for a in r),
            'sidechain_atoms':sum(a.name not in BACKBONE and a.element.name not in ['H','D']
                                  for c in structure[0] for r in c for a in r)}


def run_cycles(native,sequence,start,sequence_report,original_map,profile,residue_profile,rotamers,
               output,mode='updated-map',additional_rounds=3,initial_cycles=20,progress=print,rebuild_stage='local'):
    if mode not in ['refine-only','fixed-map','updated-map'] or additional_rounds<1 or initial_cycles<0:
        raise ValueError('Invalid mode or cycle budget')
    if rebuild_stage not in ['local','extension']:
        raise ValueError('Unknown native rebuilding stage')
    initial=validate_start(native,sequence,start,sequence_report)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    project=Path(__file__).resolve().parents[2]
    source_files=list(sorted((project/'src/phasebuild').glob('*.py')))+[project/'run_rebuild.py',project/'run_phasebuild_cycles.py']
    if rebuild_stage=='extension':source_files.append(project/'run_extend.py')
    sources=[fingerprint(f) for f in source_files]
    snapshots=out/'source_snapshot';snapshots.mkdir()
    for file in source_files:shutil.copy2(file,snapshots/file.name)
    first=Path(start).resolve();previous=first;assignment=Path(sequence_report).resolve()
    stages=[{'cumulative_cycles':initial_cycles,'refinement_report':fingerprint(first/'refinement_report.json'),
             'r_work':initial['statistics']['final']['r_work'],'r_free':initial['statistics']['final']['r_free'],
             'model':fingerprint(first/'model.cif'),'counts':model_counts(first/'model.cif')}]
    report={'engine':'phasebuild-fixed-cycle-experiment','mode':mode,'rebuild_stage':rebuild_stage,'native':fingerprint(native),
            'sequence':fingerprint(sequence),'starting_refinement':fingerprint(first/'refinement_report.json'),
            'starting_sequence_report':fingerprint(sequence_report),'original_map':fingerprint(original_map),
            'profile':fingerprint(profile),'residue_profile':fingerprint(residue_profile),'rotamers':fingerprint(rotamers),
            'initial_cumulative_cycles':initial_cycles,'additional_rounds':additional_rounds,
            'refmac_cycles_per_round':10,'stages':stages,'sources':sources,
            'scope':'The starting cumulative budget is explicitly declared from prior audited runs. '
                    'Selected native stage fixed across rounds; no target reference or Buccaneer coordinates. '
                    'Every arm uses a fixed final round, not minimum Rfree selection. '
                    'Updated REFMAC maps may contain free-reflection information; no independent blind-test claim.'}
    started=time.perf_counter()
    try:
        for index in range(1,additional_rounds+1):
            cumulative=initial_cycles+index*10
            folder=out/('cycles_%d'%cumulative);folder.mkdir()
            progress('%s: cumulative refinement cycles %d'%(mode,cumulative))
            build_report=None;model=previous/'model.cif'
            if mode!='refine-only':
                map_path=Path(original_map).resolve() if mode=='fixed-map' else previous/'refined.mtz'
                runner=project/('run_extend.py' if rebuild_stage=='extension' else 'run_rebuild.py')
                command=[sys.executable,str(runner),str(map_path),
                         '--model',str(model),'--sequence-report',str(assignment),
                         '--sequence',str(Path(sequence).resolve()),'--profile',str(Path(profile).resolve()),
                         '--residue-profile',str(Path(residue_profile).resolve()),
                         '--rotamers',str(Path(rotamers).resolve()),'--output',str(folder/'build')]
                command+=(['--beam-width','6','--maximum-extension','12'] if rebuild_stage=='extension'
                          else ['--loop-starts','48','--seed','0'])
                with (folder/'rebuild.log').open('w') as log:
                    result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,cwd=str(project),timeout=300)
                if result.returncode:
                    raise RuntimeError('Native rebuilding failed; see '+str(folder/'rebuild.log'))
                assignment=folder/'build/build_report.json';model=folder/'build/model.cif'
                build_report=fingerprint(assignment)
            refined=refmac_refine(native,model,folder/'refine',cycles=10)
            if refined['keywords']!=initial['keywords'] or any(refined['statistics'][k]!=initial['statistics'][k]
                    for k in ['reflections_used_work','reflections_used_free']):
                raise RuntimeError('Refinement settings or R-factor reflection counts changed')
            previous=folder/'refine'
            stage={'cumulative_cycles':cumulative,'refinement_report':fingerprint(previous/'refinement_report.json'),
                   'build_report':build_report,'r_work':refined['statistics']['final']['r_work'],
                   'r_free':refined['statistics']['final']['r_free'],
                   'model':fingerprint(previous/'model.cif'),'counts':model_counts(previous/'model.cif')}
            stages.append(stage);progress(json.dumps(stage))
            (out/'progress.json').write_text(json.dumps(report,indent=2)+'\n')
        report.update(status='success',final_model=fingerprint(previous/'model.cif'),
                      final_refinement_report=fingerprint(previous/'refinement_report.json'))
    except Exception as error:
        report.update(status='failed',error=str(error));raise
    finally:
        report['seconds']=time.perf_counter()-started
        report['source_changed_during_run']=any(fingerprint(s['path'])!=s for s in sources)
        (out/'experiment_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
