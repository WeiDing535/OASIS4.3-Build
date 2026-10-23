"""Separate pre-DM SAD and post-DM files using the preserved OASIS engine."""
import argparse
import json
import sys
from pathlib import Path
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_thread_key]='1'
sys.path.insert(0, str(Path(__file__).resolve().parent/'src'))


def main(argv=None):
    if sys.version_info < (3,12):raise RuntimeError('Unified SAD stages require Python 3.12 or later')
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['phase','dm']);p.add_argument('--sad',required=True)
    p.add_argument('--sequence',required=True);p.add_argument('--heavy',required=True)
    p.add_argument('--fpp',required=True);p.add_argument('--labels',default='')
    p.add_argument('--copies',type=int,default=1);p.add_argument('--cycles',type=int,default=5)
    p.add_argument('--free-value',type=int,default=0);p.add_argument('--seed',type=int,default=20261023)
    p.add_argument('--reference');p.add_argument('--initial');p.add_argument('--output',required=True)
    a=p.parse_args(argv)
    from oasis_py.cli import main as oasis_main, pairs
    from oasis_py.data import read_mtz,write_mtz
    from oasis_py.sequence import sequence_asu
    from oasis_py.reference_data import read_reference
    from phasebuild.sad_bridge import FLAG
    if a.stage=='phase':
        args=['phase',a.sad,'--mode','sad','--sad-method','likelihood','--heavy',a.heavy,
              '--sequence',a.sequence,'--copies',str(a.copies),'--fpp',a.fpp,
              '--cycles',str(a.cycles),'--seed',str(a.seed),'--free-label',FLAG,
              '--free-value',str(a.free_value),'--max-neighbors','600','--max-per-target','96','--out',a.output]
        if a.labels:args+=['--labels',a.labels]
        code=oasis_main(args)
        if code:raise RuntimeError('OASIS experimental SAD stage failed')
        import gemmi
        m=gemmi.read_mtz_file(a.output)
        if not all(x in m.column_labels() for x in ['HLA','HLB','HLC','HLD']):
            raise RuntimeError('Experimental SAD stage did not produce HL; refusing posterior substitution')
    else:
        import gemmi
        import numpy as np
        from oasis_py.phasing import PhaseResult
        from oasis_py.density_mlhl import density_modify_mlhl
        if not a.initial or not a.reference:raise ValueError('DM requires original prior and reference')
        target=Path(a.output)
        if target.exists() or target.with_suffix('.json').exists():raise FileExistsError(target)
        read_reference(a.reference)
        data,_=read_mtz(a.sad,'sad',pairs(a.labels),free_label=FLAG,free_value=a.free_value)
        initial=gemmi.read_mtz_file(a.initial);arr=np.array(initial,copy=False);labels=initial.column_labels()
        lookup={tuple(map(int,row[:3])):row for row in arr}
        if len(lookup)!=len(arr):raise ValueError('Duplicate original SAD HKLs')
        if set(lookup)!=set(map(tuple,data.hkl)):raise ValueError('Raw SAD/data HKL sets differ')
        aligned=np.array([lookup[tuple(h)] for h in data.hkl])
        prior=PhaseResult(np.deg2rad(aligned[:,labels.index('PHIB')]),aligned[:,labels.index('FOM')],None,
            aligned[:,[labels.index(x) for x in ['HLA','HLB','HLC','HLD']]],{'source':'pre-DM experimental SAD'})
        seq=sequence_asu(a.sequence,a.copies,data.cell,data.spacegroup)
        result=density_modify_mlhl(data,prior,seq['solvent_fraction'],Path(a.reference),cycles=a.cycles,
                                   feedback='best',seed=a.seed)
        write_mtz(target,data,result.phase,result.fom,result.hl,result.probability)
        target.with_suffix('.json').write_text(json.dumps({'role':'post-DM map phases; not the fixed experimental prior',
            'initial':a.initial,'sequence_asu':seq,'density_modification':result.metadata},indent=2,allow_nan=False)+'\n')
    print(json.dumps({'stage':a.stage,'output':a.output}),flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
