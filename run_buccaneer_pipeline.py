"""Fixed five-round Buccaneer/REFMAC control, reusing an audited first round."""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.benchmark import buccaneer_full
from phasebuild.refinement import refmac_refine
from phasebuild.data import fingerprint


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',required=True);p.add_argument('--sequence',required=True)
    p.add_argument('--first-build',required=True);p.add_argument('--first-refinement',required=True)
    p.add_argument('--output',required=True);p.add_argument('--rounds',type=int,default=5)
    args=p.parse_args()
    if args.rounds<2:raise ValueError('At least two rounds required')
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    first=Path(args.first_refinement).resolve()
    build=json.loads((Path(args.first_build)/'build_report.json').read_text())
    refine=json.loads((first/'refinement_report.json').read_text())
    if build['sequence']['sha256']!=fingerprint(args.sequence)['sha256'] or refine['reflections']['input']!=fingerprint(args.native):
        raise ValueError('First-round input provenance differs')
    if refine['model_input']['sha256']!=build['output']['sha256'] or refine['requested_cycles']!=10 or build['cycles']!=3:
        raise ValueError('Expected audited three-cycle build plus ten-cycle first refinement')
    if any(refine['reflections'].get(k)!=v for k,v in
           [('fp','F'),('sigfp','SIGF'),('free_column','FreeR_flag'),('free_value',0)]):
        raise ValueError('First round uses a different observation/free-set protocol')
    if not build.get('correlation_mode',False):
        raise ValueError('This pipeline requires correlation-mode first building')
    stages=[{'round':1,'build_report':fingerprint(Path(args.first_build)/'build_report.json'),
             'refinement_report':fingerprint(first/'refinement_report.json'),
             'r_work':refine['statistics']['final']['r_work'],'r_free':refine['statistics']['final']['r_free']}]
    report={'engine':'external-Buccaneer-REFMAC-fixed-pipeline','native':fingerprint(args.native),
            'sequence':fingerprint(args.sequence),'requested_rounds':args.rounds,'stages':stages,
            'scope':'Three initial and two subsequent Buccaneer cycles, ten REFMAC cycles per round. '
                    'Subsequent input model and map are only previous Buccaneer/REFMAC outputs. '
                    'Correlation mode, no waters, final fixed round; not lowest-Rfree selection. '
                    'Higher compute/refinement budget than PhaseBuild 10+10; no independent Rfree claim.'}
    started=time.perf_counter()
    try:
        previous=first
        for index in range(2,args.rounds+1):
            print('Buccaneer pipeline round %d/%d'%(index,args.rounds),flush=True)
            folder=out/('round_%d'%index)
            buccaneer_full(previous/'refined.mtz',args.sequence,folder/'build',cycles=2,
                           correlation=True,starting_model=previous/'refined.pdb')
            result=refmac_refine(args.native,folder/'build/model.cif',folder/'refine',cycles=10)
            stages.append({'round':index,'build_report':fingerprint(folder/'build/build_report.json'),
                           'refinement_report':fingerprint(folder/'refine/refinement_report.json'),
                           'r_work':result['statistics']['final']['r_work'],
                           'r_free':result['statistics']['final']['r_free']})
            previous=folder/'refine'
            print(json.dumps(stages[-1]),flush=True)
        report.update(status='success',final_model=fingerprint(previous/'model.cif'),
                      final_refinement_report=fingerprint(previous/'refinement_report.json'))
    except Exception as error:
        report.update(status='failed',error=str(error));raise
    finally:
        report['additional_seconds']=time.perf_counter()-started
        (out/'pipeline_report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
