"""Preserve residue chemistry/labels while removing classifier-unsupported runs."""
import argparse,json,sys,time
from pathlib import Path
import gemmi
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import Density,fingerprint
from phasebuild.backbone_discriminant import BackboneDiscriminant
from phasebuild.structure_pruning import prune_structure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mtz');parser.add_argument('--model',required=True)
    parser.add_argument('--discriminant',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    start=time.time();density=Density(args.mtz);classifier=BackboneDiscriminant.load(args.discriminant)
    sources=[fingerprint(Path(__file__).parent/'src/phasebuild'/name) for name in
             ['structure_pruning.py','backbone_discriminant.py','profile.py','data.py']]
    model=gemmi.read_structure(args.model)
    retained,report=prune_structure(model,density,classifier)
    retained.make_mmcif_document().write_file(str(out/'model.cif'))
    report.update(engine='phasebuild-structure-preserving-discriminant-prune',parameters=vars(args),
                  scoring_map=density.source,model_input=fingerprint(args.model),
                  discriminant=fingerprint(args.discriminant),sources=sources,
                  source_changed_during_run=any(fingerprint(s['path'])!=s for s in sources),
                  model_output=fingerprint(out/'model.cif'),seconds=time.time()-start)
    (out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['input_residues','output_residues','output_fragments','seconds']}),flush=True)


if __name__=='__main__':main()
