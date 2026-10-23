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
    p.add_argument('--method',choices=['likelihood','pplus'],default='likelihood')
    p.add_argument('--solvent-fraction',type=float)
    p.add_argument('--aoe',type=float,default=.5);p.add_argument('--triplet-weight',type=float,default=.1)
    p.add_argument('--kmin',type=float,default=.03)
    p.add_argument('--max-neighbors',type=int,default=600);p.add_argument('--max-per-target',type=int,default=96)
    p.add_argument('--rank-fit',action='store_true')
    a=p.parse_args(argv)
    from phasebuild.oasis_options import pairs
    from oasis_py.data import read_mtz,write_mtz
    from oasis_py.sequence import sequence_asu
    from oasis_py.reference_data import read_reference
    from phasebuild.sad_bridge import FLAG
    if a.stage=='phase':
        target=Path(a.output)
        if target.exists() or target.with_suffix('.json').exists():raise FileExistsError(target)
        import gemmi,numpy,scipy
        from oasis_py.phasing import phase_native
        from oasis_py.sad import phase_sad_likelihood
        from oasis_py import __version__ as oasis_version
        data,input_meta=read_mtz(a.sad,'sad',pairs(a.labels),free_label=FLAG,free_value=a.free_value)
        from oasis_py.data import read_atoms
        heavy=read_atoms(a.heavy,data)
        seq=sequence_asu(a.sequence,a.copies,data.cell,data.spacegroup)
        fpp={gemmi.Element(el).name:float(value) for el,value in pairs(a.fpp).items()}
        options=dict(cycles=a.cycles,kmin=a.kmin,max_neighbors=a.max_neighbors,max_per_target=a.max_per_target)
        if a.method=='likelihood':
            result=phase_sad_likelihood(data,seq['asu_composition'],heavy,fpp,triplet_weight=a.triplet_weight,**options)
        else:
            result=phase_native(data,'sad',seq['asu_composition'],heavy=heavy,fpp=fpp,
                                seed=a.seed,aoe=a.aoe,fit_difference=a.rank_fit,**options)
        write_mtz(target,data,result.phase,result.fom,result.hl,result.probability)
        meta={**result.metadata,'input':input_meta,'sequence_asu':seq,'source':data.source,
              'phase_method':a.method,'dependencies':dict(oasis_py=oasis_version,python=sys.version,
                  numpy=numpy.__version__,scipy=scipy.__version__,gemmi=gemmi.__version__)}
        target.with_suffix('.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
        import gemmi
        m=gemmi.read_mtz_file(a.output)
        if not all(x in m.column_labels() for x in ['HLA','HLB','HLC','HLD']):
            raise RuntimeError('Experimental SAD stage did not produce HL; refusing posterior substitution')
    else:
        import gemmi
        import numpy as np
        from oasis_py.phasing import PhaseResult
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
        initial_meta=json.loads(Path(a.initial).with_suffix('.json').read_text())
        if (a.method=='likelihood') != (initial_meta.get('method')=='sad_likelihood'):
            raise ValueError('Selected method differs from original SAD result metadata')
        prior=PhaseResult(np.deg2rad(aligned[:,labels.index('PHIB')]),aligned[:,labels.index('FOM')],
            aligned[:,labels.index('PPLUS')] if 'PPLUS' in labels else None,
            aligned[:,[labels.index(x) for x in ['HLA','HLB','HLC','HLD']]],initial_meta)
        seq=sequence_asu(a.sequence,a.copies,data.cell,data.spacegroup)
        solvent=seq['solvent_fraction'] if a.solvent_fraction is None else a.solvent_fraction
        if a.method=='likelihood':
            from oasis_py.density_mlhl import density_modify_mlhl
            result=density_modify_mlhl(data,prior,solvent,Path(a.reference),cycles=a.cycles,feedback='best',seed=a.seed)
        else:
            from oasis_py.native_density import density_modify_native
            result=density_modify_native(data,prior,solvent,Path(a.reference),cycles=a.cycles,seed=a.seed)
        write_mtz(target,data,result.phase,result.fom,result.hl,result.probability)
        target.with_suffix('.json').write_text(json.dumps({'role':'post-DM map phases; not the fixed experimental prior',
            'initial':a.initial,'phase_method':a.method,'dm_engine':'mlhl-best' if a.method=='likelihood' else 'legacy',
            'solvent_fraction_used':solvent,'solvent_source':'sequence_estimate' if a.solvent_fraction is None else 'user',
            'sequence_asu':seq,'density_modification':result.metadata},indent=2,allow_nan=False)+'\n')
    print(json.dumps({'stage':a.stage,'output':a.output}),flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
