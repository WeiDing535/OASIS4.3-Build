# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Native P+ research implementation with stable circular statistics."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
from .data import structure_factors,normalize,form_factor
from .relations import build_relations


@dataclass
class PhaseResult:
    phase: np.ndarray
    fom: np.ndarray
    probability: np.ndarray | None
    hl: np.ndarray | None
    metadata: dict


def error_weights(exponent,target=.5):
    a=np.asarray(exponent,dtype=float)
    if np.any(a<0) or not np.all(np.isfinite(a)) or not 0<target<1:
        raise ValueError('Invalid error exponents or target')
    if np.mean(a==0)>=target:
        raise ValueError('Requested error target cannot be reached with zero uncertainties')
    # Solve on a rescaled dimensionless axis so absolute root tolerances do
    # not fail when all exponents are extremely large or small.
    scale=float(np.median(a[a>0]))
    with np.errstate(over='ignore'):
        relative=a/scale
    def weights(n):
        with np.errstate(over='ignore'):
            return np.exp(-np.minimum(relative*n*n,745))
    def fun(n):return np.mean(weights(n))-target
    high=1.
    while fun(high)>0:
        high*=2
        if high>1e12:raise ValueError('Error weights cannot be calibrated')
    n=brentq(fun,0,high,xtol=1e-10)
    return weights(n),n/np.sqrt(scale)


def circular_doublet(reference,difference,probability,error):
    p=np.asarray(probability,dtype=float)
    d=np.asarray(difference,dtype=float)
    g=np.asarray(error,dtype=float)
    if np.any((p<0)|(p>1)) or np.any((g<0)|(g>1)):
        raise ValueError('Probabilities and error weights must lie in [0,1]')
    z=g*(np.cos(d)+1j*(2*p-1)*np.sin(d))
    # np.angle(0)=0 is an arbitrary phase with zero confidence.
    return (reference+np.angle(z))%(2*np.pi),np.minimum(np.abs(z),1.)


def iterate_pplus(reference,difference,error,sim_field,relations,cycles=2,
                  fixed_phase=None,fixed_fom=None,fixed_mask=None):
    ref=np.asarray(reference,dtype=float);d=np.asarray(difference,dtype=float)
    g=np.asarray(error,dtype=float);sim=np.asarray(sim_field,dtype=complex)
    n=len(ref)
    if cycles<1 or any(a.shape!=(n,) for a in [d,g,sim]):
        raise ValueError('Invalid phase iteration arrays or cycle count')
    if not all(np.all(np.isfinite(a)) for a in [ref,d,g,sim]):
        raise ValueError('Non-finite phase iteration inputs')
    mask=np.zeros(n,bool) if fixed_mask is None else np.asarray(fixed_mask,dtype=bool)
    if mask.shape!=(n,):raise ValueError('Invalid fixed-phase mask')
    if mask.any():
        if fixed_phase is None or fixed_fom is None:raise ValueError('Fixed phases require FOMs')
        fp=np.asarray(fixed_phase);fm=np.asarray(fixed_fom)
        if fp.shape!=(n,) or fm.shape!=(n,) or not np.all(np.isfinite(fp[mask])) or not np.all(np.isfinite(fm[mask])) or np.any((fm[mask]<0)|(fm[mask]>1)):
            raise ValueError('Invalid fixed phases/FOMs')
    probability=np.full(n,.5)
    phase,fom=circular_doublet(ref,d,probability,g)
    history=[]
    for cycle in range(cycles):
        if cycle>0:
            field=relations.field(phase,fom)+sim
            score=np.sin(d)*np.imag(field*np.exp(-1j*ref))
            probability=.5+.5*np.tanh(score)
            phase,fom=circular_doublet(ref,d,probability,g)
        if mask.any():phase[mask]=fp[mask]%(2*np.pi);fom[mask]=fm[mask]
        history.append({'cycle':cycle+1,'mean_fom':float(fom.mean()),
                        'mean_pplus':float(probability.mean())})
    return phase,fom,probability,history


def sad_hl(reference,amplitude,dano,sigma):
    """Exact harmonic coefficients of the small-anomalous-signal likelihood.

    log P = -(DANO-2*A*cos(phi-reference))^2/(2*sigma^2) + const.
    These are raw experimental constraints, not P+ posterior coefficients.
    """
    ref=np.asarray(reference);a=np.asarray(amplitude);sig=np.asarray(sigma)
    if np.any(sig<=0):raise ValueError('Likelihood uncertainties must be positive')
    first=2*np.asarray(dano)*a/sig**2
    second=-a*a/sig**2
    return np.column_stack([first*np.cos(ref),first*np.sin(ref),
                            second*np.cos(2*ref),second*np.sin(2*ref)])


def hl_moments(hl,samples=720,chunk=256):
    """Integrate periodic probability distributions with log-sum-exp."""
    theta=np.arange(samples)*2*np.pi/samples
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    moments=[]
    for start in range(0,len(hl),chunk):
        lp=np.asarray(hl[start:start+chunk])@basis
        weights=np.exp(lp-logsumexp(lp,axis=1)[:,None])
        moments.extend(weights@np.exp(1j*theta))
    z=np.asarray(moments)
    return np.angle(z)%(2*np.pi),np.abs(z)


def sir_hl(data,heavy_sf,sigma,samples=720):
    """Fourier projection of exact SIR amplitude likelihood; truncation recorded."""
    theta=np.arange(samples)*2*np.pi/samples
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    out=[];residual=[]
    for start in range(0,len(data.f),128):
        sl=slice(start,start+128)
        calc=np.abs(data.f[sl,None]*np.exp(1j*theta)+heavy_sf[sl,None])
        lp=-.5*((data.fderivative[sl,None]-calc)/sigma[sl,None])**2
        lp-=lp.mean(axis=1)[:,None]
        hl=2*lp@basis.T/samples
        out.append(hl);residual.extend(np.sqrt(np.mean((lp-hl@basis)**2,axis=1)))
    return np.vstack(out),float(np.mean(residual))


def _sim_field(data,atoms,composition,expected,scale):
    known=structure_factors(data,atoms)*np.sqrt(scale)
    totals=sum(n*form_factor(el,data.stol2)**2 for el,n in composition.items())
    unknown=np.zeros(len(data.f))
    for el,count in composition.items():
        used=sum(float(o)**2 for e,o in zip(atoms.elements,atoms.occupancy) if e==el)
        unknown+=max(count-used,0)*form_factor(el,data.stol2)**2
    fraction=np.clip(unknown/totals,1e-3,1)
    chi=np.minimum(2*data.f*np.abs(known)/(expected*fraction),1000.)
    return chi*np.exp(1j*np.angle(known))


def phase_native(data,mode,composition,model=None,heavy=None,fpp=None,cycles=2,
                 seed=1,pseudo_fraction=.05,aoe=.5,fit_difference=False,
                 kmin=.03,max_neighbors=1000,max_per_target=256,
                 fixed_phase=None,fixed_fom=None,fixed_mask=None):
    """Native prototype, not a numerical translation or validated replacement.

    Differences from legacy: float64 continuous phases, robust Wilson fit,
    symmetry hash relations, optional rank-fit, analytic SAD HL coefficients.
    SIR input datasets must already share a common amplitude scale.
    """
    if mode not in {'sad','sir','dmr'}:raise ValueError('Unsupported mode')
    e,expected,kappa,stats=normalize(data,composition)
    scale=stats['scale'];hl=None;extra={}
    if mode=='dmr':
        if model is None:raise ValueError('DMR requires a partial model')
        if not 0<pseudo_fraction<=1:raise ValueError('pseudo_fraction must lie in (0,1]')
        rng=np.random.default_rng(seed)
        count=max(1,round(len(model.elements)*pseudo_fraction))
        pseudo=model.subset(rng.choice(len(model.elements),count,replace=False))
        sf=structure_factors(data,pseudo,{el:2. for el in pseudo.elements})*np.sqrt(scale)
        known=structure_factors(data,model)*np.sqrt(scale)
        ref=np.angle(sf)
        d=np.abs(np.angle(known*np.exp(-1j*ref)))
        exponent=data.sigf**2/(2*np.maximum(np.abs(sf),1e-8)**2)
        extra={'pseudo_atoms':count,'seed':seed}
    else:
        if heavy is None:raise ValueError(f'{mode.upper()} requires an ASU heavy-atom model')
        if mode=='sad':
            if data.dano is None or data.sigdano is None:raise ValueError('SAD requires Bijvoet differences and uncertainties')
            if not fpp or any(el not in fpp or fpp[el]<=0 for el in heavy.elements):
                raise ValueError('Provide a positive f-double-prime for each heavy element')
            sf=structure_factors(data,heavy,fpp)*np.sqrt(scale)
            ref=np.angle(sf)
            ratio=data.dano/(2*np.maximum(np.abs(sf),1e-8))
            exponent=data.sigdano**2/(4*np.maximum(np.abs(sf),1e-8)**2)
        else:
            if data.fderivative is None:raise ValueError('SIR requires derivative amplitudes and uncertainties')
            sf=structure_factors(data,heavy)*np.sqrt(scale)
            ref=np.angle(sf)
            ratio=(data.fderivative**2-data.f**2-np.abs(sf)**2)/(2*data.f*np.maximum(np.abs(sf),1e-8))
            sd=np.hypot(data.sigf,data.sigderivative)
            exponent=sd**2*data.fderivative**2/(np.maximum(np.abs(sf),1e-8)**2*data.f**2)
        d=np.arccos(np.clip(ratio,-1,1))
        extra['clipped_doublet_fraction']=float(np.mean(np.abs(ratio)>1))
    if fit_difference:
        noncentric=np.flatnonzero(~data.centric)
        order=noncentric[np.argsort(d[noncentric],kind='stable')]
        d[order]=np.arange(1,len(order)+1)*np.pi/max(len(order),1)
    # Centric mean protein phases are restricted by crystallographic symmetry.
    # The SAD imaginary substructure is in quadrature with allowed phases.
    if mode in {'sad','dmr'}:d[data.centric]=np.pi/2
    else:
        # Convert SIR reference to the midpoint of the two allowed phases.
        ref[data.centric]+=np.pi/2;d[data.centric]=np.pi/2
    g,n=error_weights(exponent,aoe)
    if mode=='sad':
        hl=sad_hl(ref,np.abs(sf),data.dano,n*data.sigdano/np.sqrt(2))
    elif mode=='sir':
        hl,rms=sir_hl(data,sf,n*np.hypot(data.sigf,data.sigderivative)/np.sqrt(2))
        extra['sir_hl_mean_log_rms']=rms
    if model is None:
        sim=_sim_field(data,heavy,composition,expected,scale)
    else:
        # A model should include known heavy atoms when they are part of the
        # protein; do not add their contribution twice.
        sim=_sim_field(data,model,composition,expected,scale)
    graph=build_relations(data,e,kappa,kmin,max_neighbors,max_per_target)
    phase,fom,p,history=iterate_pplus(ref,d,g,sim,graph,cycles,
                                    fixed_phase,fixed_fom,fixed_mask)
    meta={'status':'research_prototype_not_legacy_equivalent','mode':mode,
          'normalization':stats,'relations':graph.statistics,'cycles':history,
          'error_scale':float(n),'mean_error_weight':float(g.mean()),
          'rank_fit_difference':fit_difference,**extra}
    return PhaseResult(phase,fom,p,hl,meta)
