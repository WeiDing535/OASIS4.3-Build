# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Experimental symmetry-constrained anisotropic Wilson amplitude correction."""
from dataclasses import replace
import numpy as np
from scipy.optimize import minimize


def fit_anisotropy(data, min_snr=3., max_scale=3.):
    """Fit trace-free Cartesian B tensor with radial nuisance parameters.

    Exponential/Gamma Wilson model for positive F²/epsilon, using strong
    observations only. Not an intensity likelihood or ellipsoidal truncation.
    Correction is limited to [1/max_scale,max_scale] to control noise gain.
    """
    if min_snr<=0 or max_scale<1:raise ValueError('Invalid anisotropy settings')
    if not data.cell.is_compatible_with_spacegroup(data.spacegroup):raise ValueError('Cell incompatible with crystal symmetry')
    basis=[]
    for a,b in [(0,0),(1,1),(0,1),(0,2),(1,2)]:
        m=np.zeros((3,3));m[a,b]=1;m[b,a]=1
        if a==b:m[2,2]=-1
        basis.append(m)
    orth=np.array(data.cell.orth.mat);frac=np.array(data.cell.frac.mat)
    rotations=[orth@(np.array(op.rot,dtype=float)/op.DEN)@frac for op in data.spacegroup.operations()]
    projected=np.array([np.mean([r@m@r.T for r in rotations],axis=0) for m in basis])
    flat=projected.reshape(5,9);u,s,v=np.linalg.svd(flat,full_matrices=False)
    keep=s>1e-7
    if not keep.any():return np.ones(len(data.f)),{'status':'isotropic_by_crystal_symmetry','tensor_cartesian_B':np.zeros((3,3)).tolist(),'independent_parameters':0}
    tensors=v[keep].reshape(-1,3,3)
    q=data.hkl@frac
    features=.5*np.einsum('ni,kij,nj->nk',q,tensors,q)
    use=(data.f/data.sigf>=min_snr)&(~data.centric)
    if use.sum()<200:return np.ones(len(data.f)),{'status':'too_few_strong_acentric_observations','independent_parameters':len(tensors)}
    # Equal-population radial shells absorb isotropic falloff and scattering.
    edges=np.unique(np.quantile(data.stol2[use],np.linspace(0,1,13)))
    shell=np.clip(np.searchsorted(edges[1:-1],data.stol2,side='right'),0,len(edges)-2)
    y=data.f[use]**2/data.epsilon[use];x=features[use];groups=shell[use];n_shell=len(edges)-1
    median=max(float(np.median(y)),1e-20);y=y/median
    intercept=np.array([np.log(max(float(np.mean(y[groups==i])),1e-12)) if (groups==i).any() else 0 for i in range(n_shell)])
    # Dimensionless conditioning; output converted back to B in A².
    scale=np.sqrt(np.mean(x*x,axis=0));scale=np.maximum(scale,1e-12);x=x/scale
    def objective(parameters):
        mu=parameters[:n_shell][groups]-x@parameters[n_shell:]
        residual=y*np.exp(np.clip(-mu,-60,60));gradient=1-residual
        return float(np.mean(mu+residual)),np.r_[np.bincount(groups,weights=gradient,minlength=n_shell),-x.T@gradient]/len(y)
    result=minimize(objective,np.r_[intercept,np.zeros(len(tensors))],jac=True,method='L-BFGS-B',bounds=[(-40,40)]*n_shell+[(-10,10)]*len(tensors),options={'maxiter':150})
    tensor=np.einsum('k,kij->ij',result.x[n_shell:]/scale,tensors)
    logscale=.25*np.einsum('ni,ij,nj->n',q,tensor,q)
    factor=np.exp(np.clip(logscale,-np.log(max_scale),np.log(max_scale)))
    if not result.success:return np.ones(len(data.f)),{'status':'fit_not_converged_no_correction','optimizer_message':str(result.message)}
    return factor,{'status':'fitted','tensor_cartesian_B':tensor.tolist(),'eigenvalues_B':np.linalg.eigvalsh(tensor).tolist(),'independent_parameters':len(tensors),'strong_acentric_n':int(use.sum()),'min_snr':min_snr,'max_scale':max_scale,'clipped_fraction':float(np.mean(abs(logscale)>np.log(max_scale))),'likelihood':'positive_amplitude Wilson approximation with radial nuisance shells'}


def corrected_reflections(data,factors):
    factors=np.asarray(factors,dtype=float)
    if factors.shape!=data.f.shape or not np.all(np.isfinite(factors)) or np.any(factors<=0):raise ValueError('Invalid anisotropy factors')
    arrays={name:None if getattr(data,name) is None else getattr(data,name)*factors for name in ['dano','sigdano','fderivative','sigderivative']}
    return replace(data,f=data.f*factors,sigf=data.sigf*factors,**arrays)
