"""Stable model building from experimental phases with no initial protein model.

PhaseBuild integrates external Buccaneer construction and REFMAC refinement.
The first experimental HL distribution stays fixed throughout refinement.
"""
import argparse
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild import hybrid_recycle
from phasebuild.fixed_phase_prior import make_fixed_prior_refiner


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('map', 'native', 'sequence', 'phase-prior', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--phase-prior-sha256')
    parser.add_argument('--rounds', type=int, default=5)
    parser.add_argument('--cycles', type=int, default=10)
    parser.add_argument('--free-value', type=int, default=0)
    parser.add_argument('--build-semet', action='store_true')
    parser.add_argument('--runtime')
    args = parser.parse_args(argv)
    hybrid_recycle.refmac_refine = make_fixed_prior_refiner(
        args.phase_prior, expected_sha256=args.phase_prior_sha256, free_value=args.free_value)
    def progress(label):
        match = re.search(r'round\s+(\d+)', label)
        number = int(match.group(1)) if match else None
        stage = 'refine' if 'REFMAC' in label else 'build'
        print(json.dumps(dict(event='progress', stage=stage, label=label, round=number,
                              completed=(number - 1) if number else 0, total=args.rounds)), flush=True)
    report = hybrid_recycle.recycle(args.map, args.native, None, args.sequence, args.output,
        args.rounds, args.cycles, ccp4=args.runtime, build_semet=args.build_semet, progress=progress)
    print(json.dumps(dict(event='building_complete', stage='complete',
                          primary_metric=report['primary_metric'],
                          primary_model=report['primary_model'])), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
