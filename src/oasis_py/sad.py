# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Experimental SAD posterior integration with an observation-only scale fit.

No reference protein coordinates or phases are accepted by this module.
The small-anomalous-signal approximation is explicit; this is not Phaser's
correlated joint Bijvoet likelihood or a claim of equivalent accuracy.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize_scalar
from .data import structure_factors,normalize
from .phasing import PhaseResult,sad_hl,hl_moments,_sim_field
from .relations import build_relations


def fit_anomalous_scale(data,imaginary_sf):
    """Fit DANO variance = 2 alpha² |Q|² + sigmaD² on acentric observations.

    This marginal second-moment approximation calibrates scale without
    knowing protein phases. Its fitted alpha also absorbs f'' uncertainty;
    it must not be interpreted as an independently measured f''.
    """
    a=np.abs(imaginary_sf);use=(~data.centric)&(a>1e-8)
    if use.sum()<20:raise ValueError('At least 20 acentric observations required for SAD scale fitting')
    signal=2*a[use]**2;noise=data.sigdano[use]**2;d2=data.dano[use]**2
    # A noise-subtracted moment may be negative even when the weighted
    # likelihood contains strong signal. Never use that difference to bound
    # the admissible scale: it previously excluded valid positive optima.
    center=np.sqrt(max(float(np.sum(d2)),1e-30)/max(float(np.sum(signal)),1e-30))
    def objective(logscale):
        variance=noise+signal*np.exp(2*logscale)
        return float(np.mean(np.log(variance)+d2/variance))
    null_objective=float(np.mean(np.log(noise)+d2/noise))
    grid=np.linspace(np.log(max(center,1e-30))-20,np.log(max(center,1e-30))+8,81)
    scores=np.array([objective(x) for x in grid])
    candidates=[]
    for i in range(1,len(grid)-1):
        if scores[i]<=scores[i-1] and scores[i]<=scores[i+1]:
            fitted=minimize_scalar(objective,bounds=(grid[i-1],grid[i+1]),method='bounded')
            if fitted.success:candidates.append((float(fitted.fun),float(fitted.x)))
    best=int(np.argmin(scores));candidates.append((float(scores[best]),float(grid[best])))
    value,logalpha=min(candidates)
    if value>=null_objective:alpha=0.;value=null_objective
    else:alpha=float(np.exp(logalpha))
    return alpha,{'anomalous_scale':alpha,'acentric_observations':int(use.sum()),
        'marginal_variance_objective':value,'scale_fit_uses_reference_phases':False,
        'anomalous_scale_fit_version':'wide_grid_null_v2',
        'null_marginal_variance_objective':null_objective,
        'loglikelihood_gain_vs_null':float(.5*use.sum()*(null_objective-value)),
        'scale_fit_grid_size':len(grid),'scale_fit_search_log_bounds':[float(grid[0]),float(grid[-1])]}


def posterior_moments(data,hl,samples=720):
    hl=np.asarray(hl,dtype=float)
    if hl.shape!=(len(data.f),4) or not np.all(np.isfinite(hl)) or samples<72:
        raise ValueError('Invalid SAD posterior coefficients or angular grid')
    phase,fom=hl_moments(hl,samples=samples)
    for i in np.flatnonzero(data.centric):
        h=data.hkl[i].tolist()
        op=next(op for op in data.spacegroup.operations() if op.apply_to_hkl(h)==[-v for v in h])
        allowed=-op.phase_shift(h)/2
        theta=np.array([allowed,allowed+np.pi]);basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
        logp=hl[i]@basis;w=np.exp(logp-np.max(logp));w/=w.sum()
        z=np.sum(w*np.exp(1j*theta));phase[i]=np.angle(z)%(2*np.pi);fom[i]=abs(z)
    return phase,np.clip(fom,0,1)


def phase_sad_likelihood(data,composition,heavy,fpp,cycles=5,kmin=.03,
                         max_neighbors=1000,max_per_target=128,samples=720,
                         damping=.5,model_error_fraction=.1,triplet_weight=.1):
    """Integrate continuous SAD likelihood plus a partial-structure/triplet prior.

    Each update combines the unchanged observation likelihood with the
    current triplet field (no repeated multiplication of measured data).
    Fixed phase seeds and protein models are deliberately not accepted.
    """
    if data.dano is None or data.sigdano is None:raise ValueError('SAD requires anomalous differences and sigmas')
    if cycles<1 or not 0<damping<=1 or not 0<=model_error_fraction<=1 or samples<72 or not 0<=triplet_weight<=1:
        raise ValueError('Invalid SAD iteration or integration parameters')
    if not fpp or any(e not in fpp or fpp[e]<=0 for e in heavy.elements):
        raise ValueError('Positive f-double-prime required for every heavy element')
    e,expected,kappa,normalization=normalize(data,composition)
    raw=structure_factors(data,heavy,fpp)
    alpha,fit=fit_anomalous_scale(data,raw);sf=alpha*raw
    sigma=np.hypot(data.sigdano,model_error_fraction*np.abs(sf))
    observations=sad_hl(np.angle(sf),np.abs(sf),data.dano,sigma)
    prior=_sim_field(data,heavy,composition,expected,normalization['scale'])
    graph=build_relations(data,e,kappa,kmin,max_neighbors,max_per_target)
    field=prior.copy();history=[]
    for cycle in range(cycles):
        posterior=observations.copy();posterior[:,0]+=field.real;posterior[:,1]+=field.imag
        phase,fom=posterior_moments(data,posterior,samples)
        history.append({'cycle':cycle+1,'mean_fom':float(fom.mean())})
        if cycle+1<cycles:
            # Empirical tempering reduces overconfidence from dependent
            # triplets; it is not an exact independence correction.
            updated=prior+triplet_weight*graph.field(phase,fom)
            field=(1-damping)*field+damping*updated
    # Branch mass is not represented by the continuous posterior's first moment.
    # Keep PhaseResult schema but explicitly report its definition.
    probability=None
    meta={'status':'research_SAD_continuous_posterior','method':'sad_likelihood',
          'normalization':normalization,'relations':graph.statistics,'cycles':history,
          'samples':samples,'damping':damping,'model_error_fraction':model_error_fraction,
          'triplet_weight':triplet_weight,
          'pplus_definition':'not_applicable_continuous_posterior_no_PPLUS_column',
          'hl_definition':'SAD_observation_likelihood_plus_final_triplet_and_partial_structure_field',**fit}
    return PhaseResult(phase,fom,probability,posterior,meta)
