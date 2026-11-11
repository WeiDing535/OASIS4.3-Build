"""Research map seed: reweight the two original SAD likelihood branches.

P+ supplies branch masses; observed SAD HL supplies conditional shapes.
Only first moments initialize a map. No reconstructed HL or additional
experimental likelihood is produced, and original observations stay fixed.
"""
import numpy as np
from scipy.special import logsumexp
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from oasis_py.clipper_likelihood import reflection_classes


def _moments(hl,reference,probability,samples):
    t=np.arange(samples)*2*np.pi/samples
    plus=np.zeros(samples);minus=np.zeros(samples)
    plus[1:samples//2]=1.;minus[samples//2+1:]=1.
    plus[[0,samples//2]]=.5;minus[[0,samples//2]]=.5
    with np.errstate(divide='ignore'):lpplus=np.log(plus);lpminus=np.log(minus)
    out=[]
    for start in range(0,len(hl),128):
        sl=slice(start,start+128);theta=reference[sl,None]+t
        lp=hl[sl,0,None]*np.cos(theta)+hl[sl,1,None]*np.sin(theta)+hl[sl,2,None]*np.cos(2*theta)+hl[sl,3,None]*np.sin(2*theta)
        pos=lp+lpplus;neg=lp+lpminus
        wp=np.exp(pos-logsumexp(pos,axis=1)[:,None]);wm=np.exp(neg-logsumexp(neg,axis=1)[:,None])
        e=np.exp(1j*theta);zp=np.sum(wp*e,axis=1);zm=np.sum(wm*e,axis=1)
        out.append(probability[sl]*zp+(1-probability[sl])*zm)
    return np.concatenate(out) if out else np.empty(0,complex)


def branch_weighted_seed(data,observation_hl,reference,probability,*,tolerance=1e-6,max_samples=11520):
    """Use q(theta)=P+(branch)*L(theta)/integral_branch L, not VM kernels.

    Half weight at branch-boundary quadrature nodes gives a trapezoid on each
    semicircle. Refinement checks grid convergence, not protein/map/Rfree fit.
    Centric reflections use their allowed discrete phases; an undefined
    branch reference leaves their original likelihood moment unchanged.
    """
    n=len(data.f);hl=np.asarray(observation_hl,dtype=float);ref=np.asarray(reference,dtype=float);p=np.asarray(probability,dtype=float)
    if hl.shape!=(n,4) or ref.shape!=(n,) or p.shape!=(n,) or not all(np.isfinite(a).all() for a in [hl,ref,p]) or np.any((p<0)|(p>1)) or tolerance<=0 or not np.isfinite(tolerance) or max_samples<1440:
        raise ValueError('Invalid SAD branch seed inputs')
    _,centric,allowed=reflection_classes(data);indices=np.flatnonzero(~centric)
    z=np.empty(n,complex);sizes=np.zeros(n,int);worst=0.
    current=_moments(hl[indices],ref[indices],p[indices],720);active=np.arange(len(indices));samples=720
    while len(active):
        samples*=2
        if samples>max_samples:raise RuntimeError('SAD branch moment quadrature did not converge')
        absolute=indices[active];refined=_moments(hl[absolute],ref[absolute],p[absolute],samples)
        error=abs(refined-current[active]);passed=error<=tolerance
        if passed.any():worst=max(worst,float(error[passed].max()))
        z[absolute[passed]]=refined[passed];sizes[absolute[passed]]=samples
        current[active]=refined;active=active[~passed]
    undefined=0
    if centric.any():
        j=np.flatnonzero(centric);s=np.sin(allowed[j]-ref[j]);valid=abs(s)>1e-7
        plus_phase=np.where(s>=0,allowed[j],allowed[j]+np.pi)
        z[j[valid]]=(2*p[j[valid]]-1)*np.exp(1j*plus_phase[valid])
        if (~valid).any():
            phase,fom=posterior_moments(data,hl);z[j[~valid]]=fom[j[~valid]]*np.exp(1j*phase[j[~valid]]);undefined=int(sum(~valid))
    if not np.isfinite(z).all() or np.any(abs(z)>1+1e-10):raise RuntimeError('Invalid conditional SAD moment')
    return PhaseResult(np.angle(z)%(2*np.pi),np.minimum(abs(z),1),p.copy(),None,{
        'method':'original_SAD_conditional_branch_map_seed_research',
        'role':'Initial map first moment only; original SAD HL remains the fixed DM likelihood',
        'branch_mass_source':'P+ probability, used without Rfree or true-phase tuning',
        'within_branch_shape':'Original observed SAD likelihood, not von Mises approximation',
        'quadrature_tolerance':tolerance,'maximum_accepted_grid_change':worst,
        'maximum_samples_used':int(sizes.max()) if n else 0,
        'centric_rows':int(sum(centric)),'centric_undefined_branch_reference_rows':undefined,
        'original_observation_modified':False,'posterior_HL_exported':False})
