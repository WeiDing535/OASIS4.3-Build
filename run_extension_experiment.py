"""Fixed three-stage anchored-growth experiment from an audited 20-cycle model."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import fingerprint
from phasebuild.refinement import refmac_refine
from phasebuild.cycling import validate_start,run_cycles


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',required=True);p.add_argument('--sequence',required=True)
    p.add_argument('--start-refinement',required=True);p.add_argument('--sequence-report',required=True)
    p.add_argument('--initial-map',required=True);p.add_argument('--original-map',required=True)
    p.add_argument('--profile',required=True);p.add_argument('--residue-profile',required=True)
    p.add_argument('--rotamers',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();validate_start(a.native,a.sequence,a.start_refinement,a.sequence_report)
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    project=Path(__file__).resolve().parent
    report={'engine':'native-extension-fixed-50-cycle-experiment','initial_map':fingerprint(a.initial_map),
            'starting_refinement':fingerprint(Path(a.start_refinement)/'refinement_report.json'),
            'scope':'Start at prior audited 20 cycles, then exactly three growth/10-cycle refinement blocks. '
                    'First map explicitly supplied; subsequent maps from own refinement. Fixed final model, no Rfree selection.'}
    try:
        first=out/'cycles_30'
        first.mkdir()
        command=[sys.executable,str(project/'run_extend.py'),str(Path(a.initial_map).resolve()),
                 '--model',str((Path(a.start_refinement)/'model.cif').resolve()),
                 '--sequence-report',str(Path(a.sequence_report).resolve()),
                 '--sequence',str(Path(a.sequence).resolve()),'--profile',str(Path(a.profile).resolve()),
                 '--residue-profile',str(Path(a.residue_profile).resolve()),'--rotamers',str(Path(a.rotamers).resolve()),
                 '--output',str(first/'build')]
        with (first/'rebuild.log').open('w') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,cwd=str(project),timeout=300)
        if result.returncode:raise RuntimeError('Initial growth failed; see '+str(first/'rebuild.log'))
        refmac_refine(a.native,first/'build/model.cif',first/'refine',cycles=10)
        print('Initial growth/refinement completed at cumulative 30 cycles',flush=True)
        continuation=run_cycles(a.native,a.sequence,first/'refine',first/'build/build_report.json',a.original_map,
                               a.profile,a.residue_profile,a.rotamers,out/'continuation',mode='updated-map',
                               initial_cycles=30,additional_rounds=2,rebuild_stage='extension',
                               progress=lambda message:print(message,flush=True))
        report.update(status='success',continuation=fingerprint(out/'continuation/experiment_report.json'),
                      first_build=fingerprint(first/'build/build_report.json'),
                      first_refinement=fingerprint(first/'refine/refinement_report.json'),
                      final_model=continuation['final_model'],final_refinement_report=continuation['final_refinement_report'])
    except Exception as error:
        report.update(status='failed',error=str(error));raise
    finally:
        (out/'pipeline_report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
