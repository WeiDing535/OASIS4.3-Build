"""Fixed-coordinate R factors, independent of REFMAC and CCP4 installations.

Atomic scattering and mask FFT: Gemmi. Work-only scale/solvent fitting: SciPy.
This is an explicitly isotropic amplitude-LS evaluation, not coordinate refinement.
"""
import json
import time
from pathlib import Path
import gemmi
import numpy as np
from scipy.optimize import least_squares
from . import __version__
from .data import fingerprint
from .refinement import reflection_partition


def fit_amplitudes(fo,fc,stol2,work,fmask=None):
    fo=np.asarray(fo,dtype=float);fc=np.asarray(fc,dtype=complex)
    stol2=np.asarray(stol2,dtype=float);work=np.asarray(work,dtype=bool)
    if not (fo.shape==fc.shape==stol2.shape==work.shape) or fo.ndim!=1:
        raise ValueError('Reflection array shapes differ')
    if not np.isfinite(fo).all() or not np.isfinite(fc).all() or not np.isfinite(stol2).all() or np.any(fo<0) or np.any(stol2<0):
        raise ValueError('Invalid amplitudes or reciprocal spacings')
    solvent=fmask is not None
    mask=np.asarray(fmask,dtype=complex) if solvent else np.zeros_like(fc)
    if mask.shape!=fc.shape or not np.isfinite(mask).all():
        raise ValueError('Invalid solvent factors')
    if work.sum()<8 or not (~work).any() or fo[work].sum()<=0 or fo[~work].sum()<=0:
        raise ValueError('Need at least eight work reflections and a nonempty free set with positive amplitude sums')
    # The fitter never receives free-set amplitudes, even for initialization.
    obs,atomic,x,sol=fo[work],fc[work],stol2[work],mask[work]
    normalizer=max(float(np.mean(obs)),1e-6)
    def predicted(parameters,atomic,x,sol):
        logk,b=parameters[:2]
        contribution=atomic
        if solvent:
            contribution=atomic+parameters[2]*np.exp(-parameters[3]*x)*sol
        return np.exp(logk-b*x)*np.abs(contribution)
    start_k=max(float(np.dot(np.abs(atomic),obs)/max(np.vdot(atomic,atomic).real,1e-20)),1e-8)
    lower=[-20.,-100.]+([0.,0.] if solvent else [])
    upper=[20.,200.]+([2.,300.] if solvent else [])
    solutions=[]
    for ksol,bsol in ([(.25,30.),(.45,80.)] if solvent else [(0.,0.)]):
        start=[float(np.clip(np.log(start_k),-19,19)),0.]+([ksol,bsol] if solvent else [])
        result=least_squares(lambda p:(predicted(p,atomic,x,sol)-obs)/normalizer,
                             start,bounds=(lower,upper),max_nfev=400,ftol=1e-10,xtol=1e-10,gtol=1e-10)
        if result.success and np.isfinite(result.fun).all():
            solutions.append(result)
    if not solutions:
        raise RuntimeError('Work-set scale/solvent fitting did not converge')
    best=min(solutions,key=lambda r:float(np.dot(r.fun,r.fun)))
    calculated=predicted(best.x,fc,stol2,mask)
    params=dict(k_overall=float(np.exp(best.x[0])),b_overall=float(best.x[1]),
                k_sol=float(best.x[2]) if solvent else 0.,b_sol=float(best.x[3]) if solvent else 0.,
                fitted_on='work reflections only',objective='unweighted amplitude least squares',
                bound_hit=bool(np.any(np.isclose(best.x,lower,atol=1e-4)) or np.any(np.isclose(best.x,upper,atol=1e-4))))
    return calculated,params


def atomic_factors(structure,hkl,dmin,method='fft'):
    structure.setup_cell_images()
    if method=='direct':
        calculator=gemmi.StructureFactorCalculatorX(structure.cell)
        return np.array([calculator.calculate_sf_from_model(structure[0],index) for index in hkl])
    if method!='fft':
        raise ValueError('Unknown structure-factor method')
    if np.prod(np.ceil(np.array(structure.cell.parameters[:3])*4/dmin))>80_000_000:
        raise ValueError('FFT grid exceeds this prototype memory limit; use --method direct or lower resolution')
    calculator=gemmi.DensityCalculatorX()
    calculator.d_min=dmin;calculator.rate=2.;calculator.cutoff=1e-6
    calculator.grid.setup_from(structure)
    calculator.set_refmac_compatible_blur(structure[0])
    calculator.put_model_density_on_grid(structure[0])
    reciprocal=gemmi.transform_map_to_f_phi(calculator.grid,half_l=True)
    return reciprocal.get_value_by_hkl(hkl,unblur=calculator.blur).astype(complex)


def mask_factors(structure,hkl,dmin):
    spacing=min(.6,dmin/2.-1e-8)
    if np.prod(np.ceil(np.array(structure.cell.parameters[:3])/spacing))>80_000_000:
        raise ValueError('Solvent grid exceeds this prototype memory limit')
    grid=gemmi.FloatGrid()
    grid.setup_from(structure,spacing=spacing)
    # These radii are built into Gemmi; no external CCP4 dictionary is opened.
    masker=gemmi.SolventMasker(gemmi.AtomicRadiiSet.Refmac)
    masker.put_mask_on_float_grid(grid,structure[0])
    return gemmi.transform_map_to_f_phi(grid,half_l=True).get_value_by_hkl(hkl).astype(complex)


def calculate_rfactors(mtz_path,model_path,output,fp='F',sigfp='SIGF',free='FreeR_flag',
                      free_value=0,bulk_solvent=True,method='fft',progress=print):
    started=time.perf_counter()
    mtz,partition=reflection_partition(mtz_path,fp,sigfp,free,free_value)
    structure=gemmi.read_structure(str(model_path))
    if len(structure)!=1 or not structure[0]:
        raise ValueError('Exactly one populated model required')
    sg=structure.find_spacegroup()
    if sg is None or sg.number!=mtz.spacegroup.number or not np.allclose(
            structure.cell.parameters,mtz.cell.parameters,rtol=0,atol=.01):
        raise ValueError('Model/data crystal cell or space group differs')
    atoms=[cra.atom for cra in structure[0].all()]
    if not atoms or any(a.element.atomic_number==0 or a.element.it92 is None or
                        not np.isfinite([a.pos.x,a.pos.y,a.pos.z,a.occ,a.b_iso,*a.aniso.elements_voigt()]).all()
                        or a.occ<0 for a in atoms):
        raise ValueError('Model has empty, unknown-element or invalid atoms')
    labels=mtz.column_labels();raw=np.array(mtz,copy=True)
    fo=raw[:,labels.index(fp)];sigma=raw[:,labels.index(sigfp)];flags=raw[:,labels.index(free)]
    observed=np.isfinite(fo)&np.isfinite(sigma)&(fo>=0)&(sigma>0)
    raw=raw[observed];fo=fo[observed];sigma=sigma[observed];flags=flags[observed]
    hkl=np.ascontiguousarray(raw[:,:3],dtype=np.int32)
    ops=sg.operations()
    allowed=np.array([np.any(index) and not ops.is_systematically_absent(index) for index in hkl])
    excluded=hkl[~allowed].tolist()
    hkl=np.ascontiguousarray(hkl[allowed]);fo=fo[allowed];sigma=sigma[allowed];flags=flags[allowed]
    asu=gemmi.ReciprocalAsu(sg)
    keys=[tuple(asu.to_asu(index,ops)[0]) for index in hkl]
    if len(set(keys))!=len(keys):
        raise ValueError('Use merged non-anomalous data: equivalent/Friedel duplicate observations found')
    work=flags!=free_value
    inv_d2=np.array([mtz.cell.calculate_1_d2(index) for index in hkl])
    if not len(hkl):
        raise ValueError('No usable non-absent reflections')
    dmin=float(1/np.sqrt(inv_d2.max()))-1e-7
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    report=dict(version=__version__,engine='phasebuild-rfactor-gemmi',status='running',model_input=fingerprint(model_path),
                reflections=partition,external_refmac_used=False,coordinates_refined=False,
                scattering='Gemmi neutral-atom IT92, no anomalous addends',method=method,
                scale_model='isotropic overall k/B plus flat bulk solvent' if bulk_solvent else 'isotropic overall k/B',
                excluded_origin_or_systematic_absences=excluded,
                scope='Fixed-coordinate amplitude R factors, scale/solvent fitted on work only. '
                      'No coordinate/ADP refinement, anisotropic overall scale, twinning or anomalous scattering. '
                      'Not numerically interchangeable with REFMAC; prior model/phase free-set leakage is not removed.')
    try:
        progress('Calculating atomic structure factors ('+method+')')
        fc=atomic_factors(structure,hkl,dmin,method)
        fmask=mask_factors(structure,hkl,dmin) if bulk_solvent else None
        progress('Fitting overall scale and solvent using work reflections only')
        predicted,params=fit_amplitudes(fo,fc,inv_d2/4,work,fmask)
        r=lambda selection:float(np.abs(fo[selection]-predicted[selection]).sum()/fo[selection].sum())
        final=dict(r_work=r(work),r_free=r(~work))
        report.update(status='success',parameters=params,statistics=dict(final=final,
            reflections_used_total=int(len(fo)),reflections_used_work=int(work.sum()),
            reflections_used_free=int((~work).sum()),resolution_high=float(1/np.sqrt(inv_d2.max())),
            resolution_low=float(1/np.sqrt(inv_d2.min()))),free_flags_preserved=True,
            free_set_used_for_fitting=False,selection='Fixed model and protocol; no free-set parameter selection')
        structure.make_mmcif_document().write_file(str(out/'model.cif'))
        np.savetxt(str(out/'rfactors.csv'),np.column_stack([hkl,fo,sigma,flags,predicted]),delimiter=',',
                   header='H,K,L,Fobs,SIGF,original_free_flag,Fmodel',comments='',fmt=['%d']*3+['%.9g']*4)
        progress('Rwork {:.2%}, Rfree {:.2%}'.format(final['r_work'],final['r_free']))
    except Exception as error:
        report.update(status='failed',error=str(error));raise
    finally:
        report['seconds']=time.perf_counter()-started
        (out/'r_factor_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
