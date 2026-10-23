# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Controlled integration of the Clipper kernel into our existing Python DM.

Only the MLHL kernel and the 2mFo-DFc feedback formula follow LGPL Clipper.
Map modification remains our older approximation. This is NOT a Parrot port,
nor Clipper's full spline weighting. No external program is executed here.
"""
import numpy as np
import gemmi
from scipy.optimize import minimize
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from oasis_py.clipper_likelihood import hl_target,reflection_classes
from oasis_py.native_density import (initial_prior,working_grid_size,metric_filter,
    reference_histogram,protein_mask,histogram_transform,apply_transform)


def fit_shells(data,fc,prior,samples=72):
    """Keep legacy shell partition/optimizer/cap; replace the likelihood kernel.

    This deliberately does NOT claim Clipper's spline scaling/optimization.
    Centric observations now enter their proper discrete likelihood.
    """
    n=len(data.f);scale=np.sqrt(np.mean(data.f**2/data.epsilon))
    fo=data.f/scale;sig=data.sigf/scale;calc=fc/scale
    eps,centric,allowed=reflection_classes(data)
    field=np.zeros(n,complex);effective_d=np.zeros(n);history=[]
    for ix in np.array_split(np.argsort(data.stol2,kind='stable'),max(1,min(12,n//250))):
        f=fo[ix];c=calc[ix];a=abs(c)
        d0=np.clip(np.sqrt(np.mean(f*f)/max(np.mean(a*a),1e-9)),.01,100)
        def target(x):
            d,t=np.exp(x)
            return hl_target(f,sig[ix],c,prior[ix],d,t,eps[ix],centric[ix],allowed[ix],samples=samples)
        x0=np.log([d0,.5]);offset=float(target(x0).value.mean())
        def objective(x):
            v=target(x);d,t=np.exp(x)
            return float(v.value.mean()-offset),np.array([v.ds.mean()*d,v.dw.mean()*t])
        fit=minimize(objective,x0,jac=True,method='L-BFGS-B',bounds=[(-7,7),(-9,5)],options={'maxiter':60})
        d,t=np.exp(fit.x)
        variance=2*sig[ix]**2+eps[ix]*np.where(centric[ix],2.,1.)*t
        raw=2*d*fo[ix]*calc[ix]/variance
        limiter=np.minimum(1,20/np.maximum(abs(raw),1e-12))
        field[ix]=raw*limiter;effective_d[ix]=d*limiter
        history.append(dict(nfit=len(ix),ncentric=int(centric[ix].sum()),D=float(d),tau=float(t),
            converged=bool(fit.success),optimizer_message=str(fit.message),
            field_cap_count=int(np.sum(limiter<1)),iterations=int(fit.nit)))
    return field,effective_d,history


def density_modify_mlhl(data,initial,solvent_fraction,reference_mtz,cycles=5,
                   feedback='mfo',seed=20261016,samples=72,cycle_callback=None):
    if (not isinstance(cycles,int) or cycles<1 or not np.isfinite(solvent_fraction)
            or not 0<solvent_fraction<1):
        raise ValueError('Invalid MLHL density modification cycles or solvent fraction')
    if feedback not in {'mfo','best'}:raise ValueError('feedback must be mfo or best')
    prior,kind,mismatch=initial_prior(data,initial)
    phase,fom=posterior_moments(data,prior)
    coefficients=data.f*fom*np.exp(1j*phase)
    def asu(values):
        return gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),
                                   np.ascontiguousarray(values,dtype=np.complex64))
    shape,rate=working_grid_size(asu(coefficients),8_000_000)
    kernel=metric_filter(data.cell,shape,4.)
    # Unchanged legacy histogram protocol, including its reference mask fraction.
    target=reference_histogram(data,fom,reference_mtz)
    rng=np.random.default_rng(seed);history=[];hl=prior.copy()
    for cycle in range(cycles):
        grid=asu(coefficients).transform_f_phi_to_map(exact_size=shape)
        rho=np.array(grid.array,dtype=float);mask=protein_mask(grid,solvent_fraction,kernel)
        modified,mapping=histogram_transform(rho,mask,target)
        noise=rng.normal(size=shape);grid.array[:]=noise;grid.symmetrize_avg()
        noise=np.array(grid.array,dtype=float);noise-=noise.mean()
        noise*=max(rho.std(),1e-12)*.01/max(noise.std(),1e-12)
        perturbed=apply_transform(rho+noise,mask,mapping)
        gamma=float(np.clip(np.sum(noise*(perturbed-modified))/max(np.sum(noise*noise),1e-30),0,.95))
        grid.array[:]=modified-gamma*rho;grid.symmetrize_avg()
        fc=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(data.hkl).astype(complex)
        field,d,fit=fit_shells(data,fc,prior,samples=samples)
        hl=prior.copy();hl[:,0]+=field.real;hl[:,1]+=field.imag
        phase,fom=posterior_moments(data,hl)
        mfo=data.f*fom*np.exp(1j*phase)
        coefficients=mfo if feedback=='mfo' else 2*mfo-d*fc
        history.append(dict(cycle=cycle+1,gamma=gamma,mean_fom=float(fom.mean()),
            protein_volume_fraction=float(mask.mean()),shell_fit=fit,
            feedback_amplitude_rms=float(np.sqrt(np.mean(abs(coefficients)**2)))))
        if cycle_callback is not None:cycle_callback(cycle+1,phase.copy(),fom.copy())
    return PhaseResult(phase,fom,None,hl,dict(method='experimental_clipper_kernel_legacy_dm',
        backend='NumPy/SciPy/Gemmi; no Parrot or Clipper runtime execution',
        feedback=feedback,cycles=history,grid_shape=shape,grid_sample_rate=rate,
        prior_kind=kind,input_HL_moment_mismatch=mismatch,seed=seed,
        likelihood_grid=samples,output_moment_grid=720,
        target_complete_model_used=False,prior_repeatedly_accumulated=False,
        beta_correction='deleted',default_algorithm_changed=False,
        limitations=['Only per-reflection Clipper MLHL kernel and coefficient formula ported',
                     'Legacy shell scaling, L-BFGS optimizer and field cap 20 retained',
                     'Reference histogram, masks and perturbation correction are legacy approximations',
                     'No NCS or anisotropy; not a complete Parrot/Clipper SFweight_spline reproduction']))
