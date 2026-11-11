"""Refine/evaluate a PhaseBuild model using native F/SIGF and existing free flags.

Two engines are available:

``refmac``   the bundled REFMAC5, the default and the one the ten-case protocol
             was calibrated against;
``internal`` the in-tree NumPy/SciPy refinement engine, which needs no external
             executable at all.  It is offered for licence-free experiments and
             for measuring how far a self-written refinement gets; see
             ``docs/NO_CCP4_REPLACEMENT.md`` for the measured gap.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild.refinement import refmac_refine


def internal_engine(args):
    from phasebuild.internal_refine import refine_structure
    from phasebuild.runtime import resolve_runtime
    monomer_root = args.monomer_root
    if monomer_root is None:
        if not args.ccp4:
            raise SystemExit('--monomer-root or --runtime is required with --engine internal')
        monomer_root = str(resolve_runtime(args.ccp4) / 'lib/data/monomers')
    report = refine_structure(args.mtz, args.model, args.output, monomer_root,
                              cycles=args.cycles, free_value=args.free_value,
                              resolution=args.resolution)
    statistics = dict(report['final'])
    statistics['engine'] = 'internal'
    statistics['input'] = {'r_work': report['initial']['r_work'],
                           'r_free': report['initial']['r_free']}
    return statistics


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mtz')
    p.add_argument('--model', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--cycles', type=int, default=10)
    p.add_argument('--runtime', '--ccp4', dest='ccp4',
                   help='Bundled runtime or existing CCP4 directory')
    p.add_argument('--engine', choices=('refmac', 'internal'), default='refmac')
    p.add_argument('--monomer-root', dest='monomer_root',
                   help='Monomer library root for --engine internal '
                        '(default: <runtime>/lib/data/monomers)')
    p.add_argument('--resolution', type=float,
                   help='High-resolution cutoff for --engine internal')
    p.add_argument('--amplitude', default='F')
    p.add_argument('--sigma', default='SIGF')
    p.add_argument('--free-column', default='FreeR_flag')
    p.add_argument('--free-value', type=int, default=0)
    args = p.parse_args()
    if args.engine == 'internal':
        statistics = internal_engine(args)
    else:
        report = refmac_refine(args.mtz, args.model, args.output, args.cycles, args.ccp4,
                               args.amplitude, args.sigma, args.free_column, args.free_value)
        statistics = report['statistics']
    print(json.dumps(statistics, indent=2))


if __name__ == '__main__':
    main()
