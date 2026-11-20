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
    parser.add_argument('--refinement-prior',choices=['original_sad','dm_phifom'],default='original_sad')
    parser.add_argument('--per-round-cleanup',action=argparse.BooleanOptionalAction,default=False)
    parser.add_argument('--cleanup-map')
    parser.add_argument('--asu-copies',type=int,default=1)
    parser.add_argument('--phase-blur-scale',type=float)
    parser.add_argument('--phase-blur-b',type=float)
    parser.add_argument('--rounds', type=int, default=5)
    parser.add_argument('--cycles', type=int, default=10)
    parser.add_argument('--free-value', type=int, default=0)
    parser.add_argument('--build-semet', action='store_true')
    parser.add_argument('--runtime')
    parser.add_argument('--building-recipe',choices=['legacy','ccp4'],default='legacy')
    parser.add_argument('--work-round-guard',action='store_true',help='Keep final round unless work R deteriorates by more than 0.05 absolute')
    args = parser.parse_args(argv)
    if (args.phase_blur_scale is None)!=(args.phase_blur_b is None):parser.error('Phase blur scale and B must be supplied together')
    blur=(args.phase_blur_scale,args.phase_blur_b) if args.phase_blur_scale is not None else None
    if blur is not None and args.refinement_prior!='dm_phifom':parser.error('DM phase confidence correction requires --refinement-prior dm_phifom')
    hybrid_recycle.refmac_refine = make_fixed_prior_refiner(
        args.phase_prior, expected_sha256=args.phase_prior_sha256, free_value=args.free_value,representation='phifom' if args.refinement_prior=='dm_phifom' else 'hl',source_kind='density_modified_work' if args.refinement_prior=='dm_phifom' else 'original_sad',**({'phase_blur':blur} if blur is not None else {}))
    if args.per_round_cleanup:
        if args.refinement_prior!='dm_phifom':parser.error('Per-round cleanup requires DM phase restraints')
        from phasebuild.iterative_cleanup import wrap_refiner
        hybrid_recycle.refmac_refine=wrap_refiner(hybrid_recycle.refmac_refine,args.cleanup_map or args.phase_prior,args.native,args.sequence,asu_copies=args.asu_copies,free_value=args.free_value)
    def progress(label):
        match = re.search(r'round\s+(\d+)', label)
        number = int(match.group(1)) if match else None
        stage = 'refine' if 'REFMAC' in label else 'build'
        print(json.dumps(dict(event='progress', stage=stage, label=label, round=number,
                              completed=(number - 1) if number else 0, total=args.rounds)), flush=True)
    report = hybrid_recycle.recycle(args.map, args.native, None, args.sequence, args.output,
        args.rounds, args.cycles, ccp4=args.runtime, build_semet=args.build_semet, progress=progress,building_recipe=args.building_recipe,free_value=args.free_value,work_round_guard=args.work_round_guard)
    print(json.dumps(dict(event='building_complete', stage='complete',
                          primary_metric=report['primary_metric'],
                          primary_model=report['primary_model'])), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
