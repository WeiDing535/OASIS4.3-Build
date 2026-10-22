"""Unrelated held-out residue identity and sequence-window validation."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
from phasebuild.data import Density, fingerprint
from phasebuild.residue_profile import ResidueProfile, training_groups, ALPHABET
from phasebuild.sequencing import windows, placement_scores


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mtz',required=True)
    p.add_argument('--structure',required=True)
    p.add_argument('--profile',required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    density=Density(args.mtz,'FC_BEST.F_phi.F','FC_BEST.F_phi.phi',resolution=1.9)
    target=ResidueProfile.load(args.profile)
    xyz,labels,_=training_groups(args.structure)
    scores=target.scores(density,xyz)
    predicted=np.argmax(scores,axis=1)
    confusion=np.zeros((20,20),dtype=int)
    np.add.at(confusion,(labels,predicted),1)
    sequence=''.join(ALPHABET[i] for i in labels)
    choices=windows([('validation_catalog',sequence)],12)
    successes=0;count=0
    for start in range(0,len(sequence)-11,12):
        best=int(np.argmax(placement_scores(scores[start:start+12],choices)))
        predicted_start=choices[best][1]
        successes+=int(sequence[predicted_start:predicted_start+12]==sequence[start:start+12])
        count+=1
    report={'validation_map':fingerprint(args.mtz),'validation_structure':fingerprint(args.structure),
            'profile':fingerprint(args.profile),'residues':len(labels),
            'top1_identity_accuracy':float(np.mean(labels==predicted)),
            'top3_identity_accuracy':float(np.mean(np.any(np.argsort(scores,axis=1)[:,-3:]==labels[:,None],axis=1))),
            'majority_class_baseline':float(np.bincount(labels,minlength=20).max()/len(labels)),
            'window_length':12,'windows_tested':count,'sequence_exact_windows':successes,
            'sequence_window_accuracy':successes/count if count else None,
            'confusion_true_by_predicted':confusion.tolist(),'alphabet':ALPHABET,
            'scope':'Known backbone coordinates on unrelated 1AJR, not a de novo building benchmark. '
                    'Window validation uses a concatenated residue catalog and accepts equivalent repeated sequences.'}
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if not k.startswith('confusion')},indent=2))


if __name__=='__main__':
    main()
