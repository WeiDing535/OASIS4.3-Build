"""Sequence weak-map backbone candidates without target/reference inputs."""
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import Density,fingerprint
from phasebuild.model import read_backbone
from phasebuild.residue_profile import ResidueProfile
from phasebuild.sequencing import read_sequences,annotate_structure
from phasebuild.fast_sequencing import assign_sequence

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mtz');p.add_argument('--model',required=True);p.add_argument('--sequence',required=True);p.add_argument('--residue-profile',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=False);start=time.time();d=Density(a.mtz);records=read_sequences(a.sequence);chains=read_backbone(a.model,d);profile=ResidueProfile.load(a.residue_profile);result=assign_sequence(chains,d,profile,records)
 annotate_structure(a.model,out/'model.cif',result,records)
 report=dict(engine='phasebuild-batched-sequence-placement',parameters=vars(a),input=d.source,sequence_input=fingerprint(a.sequence),sequence_assignment=result,seconds=time.time()-start)
 (out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(assigned=result['assigned_residues'],shuffles=result['shuffles'],seconds=report['seconds'])),flush=True)
if __name__=='__main__':main()
