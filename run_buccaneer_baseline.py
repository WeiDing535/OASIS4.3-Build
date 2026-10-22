"""A separate established-tool control; never reads the target reference PDB.

The CCP4 library reference is an unrelated training structure, not INS.
No sequence assignment or refinement is requested, for a backbone-only control.
"""
import argparse
import json
import subprocess
import time
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mtz', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--ccp4', default='/Applications/ccp4-9')
    args = p.parse_args()
    root, out = Path(args.ccp4), Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    library = root / 'lib/data/reference_structures/reference-1tqw'
    cmd = [str(root / 'bin/cbuccaneer'), '-mtzin', str(Path(args.mtz).resolve()),
           '-mtzin-ref', str(library.with_suffix('.mtz')), '-pdbin-ref', str(library.with_suffix('.pdb')),
           '-colin-ref-fo', '/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
           '-colin-ref-hl', '/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
           '-colin-fo', '/*/*/[F,SIGF]', '-colin-phifom', '/*/*/[PHWT,FOM]',
           '-colin-fc', '/*/*/[FWT,PHWT]', '-pdbout', str(out / 'backbone.pdb'),
           '-xmlout', str(out / 'buccaneer.xml'), '-find', '-grow', '-join', '-link', '-prune',
           '-cycles', '1', '-jobs', '1']
    started = time.perf_counter()
    with (out / 'buccaneer.log').open('w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    report = {'command': cmd, 'returncode': result.returncode,
              'elapsed_seconds': time.perf_counter() - started,
              'scope': 'One-cycle Buccaneer backbone-only control with unrelated CCP4 library 1tqw. '
                       'No target reference, sequence assignment, or refinement. '
                       'Initial map uses FWT/PHWT; subsequent internal weighting is tool-specific.'}
    (out / 'baseline_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if result.returncode:
        raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
