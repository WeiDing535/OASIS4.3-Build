"""Experimental local sequence anchoring with a search-matched permutation null."""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import Density, fingerprint
from phasebuild.local_sequencing import assign_local_sequence
from phasebuild.model import read_backbone
from phasebuild.residue_profile import ResidueProfile
from phasebuild.sequencing import read_sequences, annotate_structure
from phasebuild.trace import write_backbone


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mtz')
    parser.add_argument('--model', required=True)
    parser.add_argument('--sequence', required=True)
    parser.add_argument('--residue-profile', required=True)
    parser.add_argument('--target-kind', choices=['density', 'logistic'], default='density')
    parser.add_argument('--output', required=True)
    parser.add_argument('--amplitude', default='FWT')
    parser.add_argument('--phase', default='PHWT')
    parser.add_argument('--minimum-length', type=int, default=8)
    parser.add_argument('--residue-penalty', type=float, default=.5)
    parser.add_argument('--shuffles', type=int, default=499)
    parser.add_argument('--alpha', type=float, default=.05)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    output = Path(args.output); output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    density = Density(args.mtz, args.amplitude, args.phase)
    records = read_sequences(args.sequence)
    chains = read_backbone(args.model, density)
    if args.target_kind == 'logistic':
        from phasebuild.residue_logistic import ResidueLogisticProfile
        target = ResidueLogisticProfile.load(args.residue_profile)
    else:
        target = ResidueProfile.load(args.residue_profile)
    split, report = assign_local_sequence(chains, density, target, records,
                                          minimum_length=args.minimum_length, shuffles=args.shuffles,
                                          alpha=args.alpha, seed=args.seed, residue_penalty=args.residue_penalty)
    write_backbone(output/'backbone.cif', split, density)
    annotate_structure(output/'backbone.cif', output/'model.cif', report, records)
    result = {'engine': 'phasebuild-local-sequence-anchors', 'parameters': vars(args),
              'input': density.source, 'model_input': fingerprint(args.model),
              'sequence_input': fingerprint(args.sequence), 'residue_profile': fingerprint(args.residue_profile),
              'sequence_assignment': report, 'seconds': time.time()-started,
              'scope': 'Experimental reassignment from backbone density; previous residue labels are not input evidence.'}
    (output/'build_report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'assigned': report['assigned_residues'], 'shuffles': report['shuffles'],
                      'input_fragments': len(chains), 'output_fragments': len(split),
                      'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    main()
