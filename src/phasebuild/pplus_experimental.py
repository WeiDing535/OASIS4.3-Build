"""Opt-in P+ research trials; the preserved oasis_py implementation is unchanged.

The bimodal DM prior is a maximum-entropy, two-harmonic approximation matching
the first two moments of a weighted pair of von Mises kernels. It is not the
original observation likelihood, nor an exact representation of that mixture.
"""
from dataclasses import dataclass
from itertools import product
import hashlib
import json
import numpy as np
from scipy.special import i0e, i1e, ive, logsumexp
from scipy.stats import rankdata
from oasis_py.data import normalize, structure_factors
from oasis_py.phasing import (PhaseResult, error_weights, iterate_pplus,
                             sad_hl, _sim_field)
from oasis_py.relations import build_relations
from oasis_py.sad import posterior_moments


def rank_cosine(difference, centric):
    """Fit acentric cos(delta) quantiles, preserving ties and centric rows."""
    d=np.asarray(difference,dtype=float)
    centric=np.asarray(centric,dtype=bool)
    if d.ndim!=1 or centric.shape!=d.shape or not np.isfinite(d).all() or np.any((d<0)|(d>np.pi)):
        raise ValueError('Invalid doublet angles or centric mask')
    result=d.copy();use=~centric
    if use.any():
        ranks=rankdata(d[use],method='average')
        result[use]=np.arccos(1-2*(ranks-.5)/use.sum())
    return result


def kernel_moments(error):
    """First/second circular moments of VM kernels; cap explicitly at .9999."""
    g=np.asarray(error,dtype=float)
    if not np.isfinite(g).all() or np.any((g<0)|(g>1)):
        raise ValueError('Error attenuation must be finite and in [0,1]')
    r=np.minimum(g,.9999)
    k=np.where(r<.53,2*r+r**3+5*r**5/6,
        np.where(r<.85,-.4+1.39*r+.43/np.maximum(1-r,1e-12),
                 1/np.maximum(r**3-4*r*r+3*r,1e-12)))
    for _ in range(12):
        a=i1e(k)/i0e(k)
        derivative=np.where(k>1e-8,1-a/np.maximum(k,1e-30)-a*a,.5)
        k=np.maximum(k-(a-r)/np.maximum(derivative,1e-12),0)
    return r,ive(2,k)/i0e(k)


def maximum_entropy_hl(target, samples=720, tolerance=2e-8, max_iterations=400, chunk=128):
    """Solve the convex log-partition dual for four circular moment constraints.

    Analytic covariance Hessians and per-row backtracking are used. A failed
    moment fit raises; no silent first-moment-only fallback is permitted.
    """
    target=np.asarray(target,dtype=float)
    if target.ndim!=2 or target.shape[1]!=4 or not np.isfinite(target).all() or samples<180:
        raise ValueError('Invalid circular moment constraints or grid')
    theta=np.arange(samples)*2*np.pi/samples
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    products=(basis[:,None,:]*basis[None,:,:]).reshape(16,samples)
    outputs=[];worst=0.;iterations=0
    for start in range(0,len(target),chunk):
        wanted=target[start:start+chunk];h=np.zeros_like(wanted)
        def evaluate(coeff):
            lp=coeff@basis;z=logsumexp(lp,axis=1)
            w=np.exp(lp-z[:,None]);mean=w@basis.T
            return z-np.sum(coeff*wanted,axis=1),w,mean
        for iteration in range(max_iterations):
            objective,weights,mean=evaluate(h)
            delta=mean-wanted
            converged=np.max(np.abs(delta),axis=1)<tolerance
            if converged.all():break
            covariance=(weights@products.T).reshape(-1,4,4)-mean[:,:,None]*mean[:,None,:]
            step=np.linalg.solve(covariance+np.eye(4)*1e-11,delta)
            step[converged]=0
            # A fixed coefficient-step cap stalls close, sharp doublets whose
            # moment-matching coefficients legitimately exceed 1e6. Grow the
            # trust radius with the current iterate; Armijo still controls
            # every accepted update of the same convex objective.
            radius=np.maximum(1000.,.5*np.max(np.abs(h),axis=1))
            step*=np.minimum(1,radius/np.maximum(np.max(np.abs(step),axis=1),1e-30))[:,None]
            rate=np.ones(len(h));slope=np.sum(delta*step,axis=1)
            for _ in range(30):
                candidate=h-rate[:,None]*step
                trial,_,_=evaluate(candidate)
                accepted=trial<=objective-1e-4*rate*slope+1e-12
                if accepted.all():break
                rate[~accepted]*=.5
            h=candidate
        residual=float(np.max(np.abs(evaluate(h)[2]-wanted))) if len(h) else 0.
        if residual>tolerance*2:
            raise RuntimeError('Bimodal moment fit did not converge: %.3g'%residual)
        outputs.append(h);worst=max(worst,residual);iterations=max(iterations,iteration+1)
    return np.vstack(outputs) if outputs else np.empty((0,4)),{
        'maximum_moment_residual':worst,'maximum_iterations_used':iterations,
        'samples':samples,'tolerance':tolerance,
        'distribution':'maximum entropy matching first and second circular moments; higher mixture moments not retained'}


def bimodal_prior(data, reference, difference, probability, error):
    ref,d,p,g=[np.asarray(a,dtype=float) for a in [reference,difference,probability,error]]
    if any(a.shape!=(len(data.f),) for a in [ref,d,p,g]) or not all(np.isfinite(a).all() for a in [ref,d,p,g]) or np.any((p<0)|(p>1)) or np.any((d<0)|(d>np.pi)):
        raise ValueError('Invalid P+ branch arrays')
    first,second=kernel_moments(g)
    z1=first*np.exp(1j*ref)*(np.cos(d)+1j*(2*p-1)*np.sin(d))
    z2=second*np.exp(2j*ref)*(np.cos(2*d)+1j*(2*p-1)*np.sin(2*d))
    target=np.column_stack([z1.real,z1.imag,z2.real,z2.imag])
    hl=np.zeros((len(ref),4));use=~data.centric
    hl[use],diagnostics=maximum_entropy_hl(target[use])
    # Centric priors are integrated on allowed two-state phases downstream.
    hl[~use,0]=np.arctanh(np.minimum(np.abs(z1[~use]),.9999))*np.cos(np.angle(z1[~use]))
    hl[~use,1]=np.arctanh(np.minimum(np.abs(z1[~use]),.9999))*np.sin(np.angle(z1[~use]))
    phase,fom=posterior_moments(data,hl)
    diagnostics.update(method='pplus_bimodal_maxentropy_experimental',
        role='P+ posterior approximation for DM only; never an experimental observation prior',
        attenuation_cap=.9999,capped_attenuations=int(np.sum(g>.9999)),
        maximum_first_moment_difference=float(np.max(np.abs(fom*np.exp(1j*phase)-z1))))
    return PhaseResult(phase,fom,p.copy(),hl,diagnostics)


@dataclass
class PplusTrial:
    raw: PhaseResult
    dm_initial: PhaseResult


@dataclass(frozen=True)
class PplusNormalization:
    identity: str
    e: np.ndarray
    expected: np.ndarray
    kappa: float
    statistics: dict


def _normalization_identity(data,composition):
    digest=hashlib.sha256(data.hkl.tobytes()+data.f.tobytes())
    digest.update(json.dumps({'cell':list(data.cell.parameters),
        'symmetry':[op.triplet() for op in data.spacegroup.operations()],
        'composition':{k:float(v) for k,v in composition.items()}},sort_keys=True).encode())
    return digest.hexdigest()


def prepare_pplus_normalization(data,composition):
    """Fit once per fixed input, avoiding optimizer-rounding confounding."""
    e,expected,kappa,stats=normalize(data,composition)
    e.setflags(write=False);expected.setflags(write=False)
    return PplusNormalization(_normalization_identity(data,composition),e,expected,kappa,dict(stats))


def _four_moments(hl,samples=720,chunk=128):
    theta=np.arange(samples)*2*np.pi/samples
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    output=[]
    for start in range(0,len(hl),chunk):
        lp=np.asarray(hl[start:start+chunk],dtype=float)@basis
        weights=np.exp(lp-logsumexp(lp,axis=1)[:,None])
        output.append(weights@basis.T)
    return np.vstack(output) if output else np.empty((0,4))


def mtz_safe_bimodal_prior(data,prior,moment_tolerance=1e-5):
    """Quantize HL explicitly, preserving both moments within a fixed bound.

    Large cancelling HL coefficients can lose branch weights in float32.
    Search neighboring representable coefficients only when naive rounding
    exceeds the numerical bound. Never change the target distribution or
    observation likelihood. PHIB/FOM describe the actual quantized HL.
    """
    hl=np.asarray(prior.hl,dtype=float)
    if hl.shape!=(len(data.f),4) or not np.isfinite(hl).all() or not np.isfinite(moment_tolerance) or moment_tolerance<=0:
        raise ValueError('Invalid HL or MTZ moment tolerance')
    target=_four_moments(hl);quantized=hl.astype(np.float32)
    if not np.isfinite(quantized).all():raise ValueError('HL exceeds finite MTZ float32 range')
    def error(values,wanted):
        delta=_four_moments(values)-wanted
        return np.maximum(np.linalg.norm(delta[:,:2],axis=1),np.linalg.norm(delta[:,2:],axis=1))
    initial_error=error(quantized,target);repaired=0;rescaled=0
    for i in np.flatnonzero((initial_error>moment_tolerance)&~data.centric):
        center=quantized[i].copy();best=center;best_error=initial_error[i]
        for radius in [1,2,4,8]:
            levels=[]
            for c in center:
                down=[c];up=[c]
                for _ in range(radius):
                    down.append(np.nextafter(down[-1],np.float32(-np.inf)))
                    up.append(np.nextafter(up[-1],np.float32(np.inf)))
                levels.append(list(reversed(down[1:]))+[c]+up[1:])
            candidates=np.asarray(list(product(*levels)),dtype=np.float32)
            scores=error(candidates,target[i:i+1])
            j=int(np.argmin(scores))
            if scores[j]<best_error:best,best_error=candidates[j],float(scores[j])
            if best_error<=moment_tolerance:break
        if best_error>moment_tolerance:
            # Closely spaced sharp modes have large, cancelling coefficients.
            # Moving along this near-flat coefficient direction before rounding
            # can avoid a coarse float32 lattice. The target moments and their
            # hard error bound stay unchanged; this is not a confidence fit.
            scales=np.linspace(.9999,1.0001,257)
            candidates=(scales[:,None]*hl[i]).astype(np.float32)
            scores=error(candidates,target[i:i+1]);j=int(np.argmin(scores))
            if scores[j]<best_error:best,best_error=candidates[j],float(scores[j])
            if best_error<=moment_tolerance:rescaled+=1
        if best_error>moment_tolerance:
            raise RuntimeError('MTZ float32 HL cannot retain circular moments: %.3g'%best_error)
        quantized[i]=best;repaired+=1
    phase,fom=posterior_moments(data,quantized.astype(float))
    first_error=float(np.max(abs(fom*np.exp(1j*phase)-prior.fom*np.exp(1j*prior.phase))))
    final_error=error(quantized,target)
    second_error=float(np.max(final_error[~data.centric])) if np.any(~data.centric) else 0.
    if max(first_error,second_error)>moment_tolerance:
        raise RuntimeError('Quantized P+ prior exceeds fixed moment error tolerance')
    metadata={**prior.metadata,'mtz_hl_quantization':{
        'precision':'float32 HL, with PHIB/FOM recomputed from the stored coefficients',
        'moment_tolerance':moment_tolerance,'rows_with_neighbor_rounding':repaired,
        'rows_with_small_coefficient_rescaling':rescaled,'maximum_rescaling_search_fraction':.0001,
        'maximum_naive_four_moment_error':float(np.max(initial_error[~data.centric])) if np.any(~data.centric) else 0.,
        'maximum_final_four_moment_error':second_error,
        'maximum_first_moment_difference_from_unquantized':first_error,
        'maximum_ulp_radius':8,'selection':'Numerical circular-moment error only; no Rfree, model or map fit'}}
    return PhaseResult(phase,fom,None if prior.probability is None else prior.probability.copy(),quantized.astype(float),metadata)


@dataclass(frozen=True)
class PplusCalibration:
    identity: str
    reference: np.ndarray
    difference: np.ndarray
    error: np.ndarray
    observations: np.ndarray
    sim: np.ndarray
    graph: object
    metadata: dict


def _calibration_identity(data,composition,heavy,fpp,single,aoe,kmin,max_neighbors,max_per_target):
    digest=hashlib.sha256(_normalization_identity(data,composition).encode())
    for a in [data.sigf,data.dano,data.sigdano,heavy.fractional,heavy.occupancy,heavy.b_iso,single]:
        digest.update(np.asarray(a).tobytes())
    digest.update(json.dumps({'elements':heavy.elements,'fpp':fpp,'aoe':aoe,'kmin':kmin,
        'max_neighbors':max_neighbors,'max_per_target':max_per_target},sort_keys=True).encode())
    return digest.hexdigest()


def prepare_pplus_calibration(data,composition,heavy,fpp,*,aoe=.5,kmin=.03,max_neighbors=600,
                              max_per_target=96,centric_singletons=None,normalization_context=None):
    """Freeze all measured-data and structural fields, not only Wilson scale."""
    if data.dano is None or data.sigdano is None:
        raise ValueError('P+ trial requires SAD differences and uncertainties')
    if not fpp or any(el not in fpp or not np.isfinite(fpp[el]) or fpp[el]<=0 for el in heavy.elements):
        raise ValueError('Positive finite f-double-prime required')
    context=normalization_context or prepare_pplus_normalization(data,composition)
    if context.identity!=_normalization_identity(data,composition):
        raise ValueError('Frozen normalization does not match reflections/composition')
    e,expected,kappa,stats=context.e,context.expected,context.kappa,dict(context.statistics)
    sf=structure_factors(data,heavy,fpp)*np.sqrt(stats['scale']);ref=np.angle(sf)
    ratio=data.dano/(2*np.maximum(np.abs(sf),1e-8))
    d=np.arccos(np.clip(ratio,-1,1));d[data.centric]=np.pi/2
    exponent=data.sigdano**2/(4*np.maximum(np.abs(sf),1e-8)**2)
    single=np.zeros(len(data.f),bool) if centric_singletons is None else np.asarray(centric_singletons,dtype=bool)
    if single.shape!=(len(data.f),) or np.any(single&~data.centric) or np.all(single):
        raise ValueError('Invalid centric singleton mask or no measured anomalous pairs')
    g=np.ones(len(data.f))
    g[~single],n=error_weights(exponent[~single],aoe)
    observations=sad_hl(ref,np.abs(sf),data.dano,n*data.sigdano/np.sqrt(2))
    observations[single]=0.
    sim=_sim_field(data,heavy,composition,expected,stats['scale'])
    graph=build_relations(data,e,kappa,kmin,max_neighbors,max_per_target)
    meta=dict(method='pplus_trial_experimental',normalization=stats,error_scale=float(n),
        mean_error_weight=float(g.mean()),relations=graph.statistics,
        original_clipped_doublet_fraction=float(np.mean(np.abs(ratio)>1)),
        centric_singletons_without_anomalous_likelihood=int(single.sum()),
        calibrated_error_mean=float(g[~single].mean()),
        normalization_context_sha256=context.identity,
        hl_definition='Unmodified original observation likelihood; rank fitting is not applied to observations')
    for a in [ref,d,g,observations,sim]:a.setflags(write=False)
    identity=_calibration_identity(data,composition,heavy,fpp,single,aoe,kmin,max_neighbors,max_per_target)
    return PplusCalibration(identity,ref,d,g,observations,sim,graph,meta)


def phase_pplus_trial(data,composition,heavy,fpp,*,cosine_rank=False,preserve_bimodal=False,
                      cycles=5,aoe=.5,kmin=.03,max_neighbors=600,max_per_target=96,
                      centric_singletons=None,normalization_context=None,calibration_context=None):
    """Fixed two-factor trial: cosine-rank fitting and DM branch-moment retention."""
    if calibration_context is None:
        context=prepare_pplus_calibration(data,composition,heavy,fpp,aoe=aoe,kmin=kmin,
            max_neighbors=max_neighbors,max_per_target=max_per_target,centric_singletons=centric_singletons,
            normalization_context=normalization_context)
    else:
        if normalization_context is not None:raise ValueError('Use one complete calibration context')
        context=calibration_context
        single=np.zeros(len(data.f),bool) if centric_singletons is None else np.asarray(centric_singletons,dtype=bool)
        if context.identity!=_calibration_identity(data,composition,heavy,fpp,single,aoe,kmin,max_neighbors,max_per_target):
            raise ValueError('Frozen P+ calibration does not match input/settings')
    d=rank_cosine(context.difference,data.centric) if cosine_rank else context.difference
    phase,fom,p,history=iterate_pplus(context.reference,d,context.error,context.sim,context.graph,cycles)
    meta={**context.metadata,'cycles':history,'cosine_rank':cosine_rank,'preserve_bimodal':preserve_bimodal,
        'calibration_context_sha256':context.identity}
    raw=PhaseResult(phase,fom,p,context.observations.copy(),meta)
    initial=mtz_safe_bimodal_prior(data,bimodal_prior(data,context.reference,d,p,context.error)) if preserve_bimodal else raw
    return PplusTrial(raw,initial)
