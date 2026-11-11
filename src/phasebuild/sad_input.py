"""Research SAD-pair reader retaining standard single-valued centric rows."""
from pathlib import Path
import numpy as np
import gemmi
from oasis_py.data import Reflections


def read_sad_pairs(path,free_label='FreeR_flag',free_value=0):
    """Read F(+/-), preserving centric singletons without inventing a pair.

    A positive sigma placeholder satisfies the shared Reflections schema, but
    singleton_mask MUST be passed to phase_pplus_trial: its centric rows have
    no anomalous observation likelihood and never enter error calibration.
    """
    m=gemmi.read_mtz_file(str(path));m.ensure_asu();a=np.array(m,copy=True)
    labels=m.column_labels()
    def column(name,types):
        if labels.count(name)!=1 or m.column_with_label(name).type not in types:
            raise ValueError('Missing or ambiguous SAD column/type: '+name)
        return a[:,labels.index(name)]
    fp,sp,fm,sm=[column(k,t) for k,t in [('F(+)','FG'),('SIGF(+)','QL'),('F(-)','FG'),('SIGF(-)','QL')]]
    hkl=m.make_miller_array();centric=m.spacegroup.operations().centric_flag_array(hkl)
    plus=np.isfinite(fp)&np.isfinite(sp)&(fp>0)&(sp>0)
    minus=np.isfinite(fm)&np.isfinite(sm)&(fm>0)&(sm>0)
    paired=plus&minus;single=centric&(plus^minus)
    keep=(paired|single)&~m.spacegroup.operations().systematic_absences(hkl)&np.any(hkl!=0,axis=1)
    if free_label is not None:
        flags=column(free_label,'I');keep&=np.isfinite(flags)&(flags!=free_value)
    f=(fp+fm)/2;sig=np.hypot(sp,sm)/2;dano=fp-fm;sd=np.hypot(sp,sm)
    f[single]=np.where(plus[single],fp[single],fm[single])
    sig[single]=np.where(plus[single],sp[single],sm[single])
    dano[single]=0.;sd[single]=sig[single]
    data=Reflections(hkl[keep],f[keep],sig[keep],m.cell,m.spacegroup,
        dano=dano[keep],sigdano=sd[keep],source=str(Path(path).resolve()))
    singleton_mask=single[keep]
    meta={'input_rows':len(a),'retained_rows':int(keep.sum()),'centric_singletons_retained':int(singleton_mask.sum()),
        'singletons_are_measured_bijvoet_pairs':False,
        'singleton_likelihood':'No anomalous likelihood; phase choices restricted by crystal symmetry; sigma placeholder not used for error calibration',
        'noncentric_missing_pairs_excluded':int(np.sum(~centric&~paired)),
        'free_value_excluded':free_value if free_label is not None else None}
    return data,singleton_mask,meta


def read_phase_prior(path,data,metadata=None):
    """Read a persisted research prior by HKL, independent of MTZ row order."""
    from oasis_py.phasing import PhaseResult
    m=gemmi.read_mtz_file(str(path))
    if m.spacegroup.operations()!=data.spacegroup.operations() or not np.allclose(m.cell.parameters,data.cell.parameters,rtol=1e-7,atol=1e-5):
        raise ValueError('Phase-prior cell or symmetry differs from observations')
    m.ensure_asu();a=np.array(m,copy=True);labels=m.column_labels()
    lookup={tuple(h):i for i,h in enumerate(m.make_miller_array())}
    if len(lookup)!=len(a):raise ValueError('Duplicate phase-prior HKLs')
    if set(lookup)!=set(map(tuple,data.hkl)):raise ValueError('Phase-prior/data HKL sets differ')
    a=a[[lookup[tuple(h)] for h in data.hkl]]
    required=['PHIB','FOM','HLA','HLB','HLC','HLD']
    if any(labels.count(k)!=1 for k in required):raise ValueError('Missing or ambiguous phase-prior columns')
    values=a[:,[labels.index(k) for k in required]].astype(float)
    if not np.isfinite(values).all() or np.any((values[:,1]<0)|(values[:,1]>1)):
        raise ValueError('Invalid phase-prior coefficients or FOM')
    probability=a[:,labels.index('PPLUS')].astype(float) if labels.count('PPLUS')==1 else None
    return PhaseResult(np.deg2rad(values[:,0]),values[:,1],probability,values[:,2:],dict(metadata or {}))
