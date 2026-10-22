"""Bounded weak-map growth: accumulated evidence may bridge short density troughs.

An opt-in research path. No target/reference or sequence labels are accepted.
"""
import numpy as np
from scipy.optimize import minimize
from .geometry import CONFORMATIONS,extend,prepend,oxygen


def growth_evidence(scores,baseline=2.,maximum_weak_steps=2,previous_weak_steps=0):
    """Allow short negative increments, never unbounded low-evidence runs."""
    gain=float(scores)-baseline
    weak=previous_weak_steps+1 if gain<0 else 0
    return gain,weak,weak<=maximum_weak_steps


def clash_free(residue,chain,density,used=(),backward=False):
    other=chain[1:] if backward else chain[:-1]
    if len(other):
        if np.min(np.linalg.norm(residue[:3,None,:]-other[:,:3].reshape(1,-1,3),axis=-1))<2.1:return False
        if any(density.distance(residue[1],ca)<2.8 for ca in other[:,1]):return False
    if any(density.distance(residue[1],ca)<2.5 for ca in used):return False
    return True


def grow_weak(density,xyz,maximum=30,beam_width=4,baseline=2.,maximum_weak_steps=2,backward=False,used=()):
    if maximum_weak_steps<0 or maximum_weak_steps>3:raise ValueError('Weak-step budget must be 0..3')
    if len(xyz)<2 or maximum<len(xyz):return xyz.copy()
    beams=[(0.,xyz.copy(),0)];best=beams[0]
    pairs=list(dict.fromkeys([(psi,phi) for phi,psi in CONFORMATIONS]+[(psi,phi) for phi in [-150.,-120.,-90.,-60.] for psi in [-60.,-30.,0.,120.,150.,180.]]))
    back_pairs=[(phi,psi) for phi in [-150.,-120.,-90.,-63.,60.] for psi in [-60.,-42.,0.,40.,120.,150.,180.]]
    for _ in range(maximum-len(xyz)):
        coarse=[]
        for bi,(cumulative,chain,weak) in enumerate(beams):
            variants=np.array([prepend(chain[0],*p) if backward else extend(chain[-1],*p) for p in (back_pairs if backward else pairs)])
            values=density.score(variants);atoms=density.sample(variants[:,:3]).min(axis=1)
            for pi in np.argsort(values)[::-1]:
                if values[pi]<.5 or atoms[pi]<-.5:continue
                if not clash_free(variants[pi],chain,density,used,backward):continue
                coarse.append((cumulative+float(values[pi])-baseline,bi,pi))
                # Several torsion alternatives per live beam, then a global budget.
                if sum(x[1]==bi for x in coarse)>=beam_width:break
        coarse.sort(reverse=True);branches=[]
        for _,bi,pi in coarse[:beam_width*2]:
            cumulative,chain,weak=beams[bi]
            if backward:
                start=back_pairs[pi]
                def proposal(params):return prepend(chain[0],*params)
                def loss(params):return -float(density.score(proposal(params)))
                fit=minimize(loss,start,method='L-BFGS-B',bounds=[(v-25.,v+25.) for v in start],options={'maxiter':20,'eps':.15,'ftol':1e-5})
                residue=proposal(fit.x);value=-float(fit.fun)
                new=np.concatenate([residue[None],chain])
            else:
                psi,phi=pairs[pi]
                def proposal(params):
                    residue=extend(chain[-1],params[0],params[1]);residue[3]=oxygen(*residue[:3],params[2]);previous=chain[-1].copy();previous[3]=oxygen(*previous[:3],params[0]);return residue,previous
                def loss(params):
                    a,b=proposal(params);return -float(.8*density.score(a)+.2*density.score(b))
                fit=minimize(loss,[psi,phi,-40.],method='L-BFGS-B',bounds=[(psi-20.,psi+20.),(phi-20.,phi+20.),(-180.,180.)],options={'maxiter':18,'eps':.15,'ftol':1e-5})
                residue,previous=proposal(fit.x);value=-float(fit.fun)
                new=np.concatenate([chain[:-1],previous[None],residue[None]])
            if np.min(density.sample(residue[:3]))<-.5 or value<.5:continue
            if not clash_free(residue,chain,density,used,backward):continue
            gain,nweak,ok=growth_evidence(value,baseline,maximum_weak_steps,weak)
            if not ok:continue
            branches.append((cumulative+gain,new,nweak))
        if not branches:break
        branches.sort(key=lambda x:x[0],reverse=True);beams=branches[:beam_width]
        if beams[0][0]>best[0]:best=beams[0]
    return best[1]


def supported_intervals(scores,cutoff=2.5,weak_steps=2):
    """Bridge only bounded weak islands with positive total map evidence."""
    scores=np.asarray(scores);strong=scores>=cutoff
    good=strong.copy()
    i=0
    while i<len(scores):
        if strong[i]:i+=1;continue
        start=i
        while i<len(scores) and not strong[i]:i+=1
        if start>0 and i<len(scores) and i-start<=weak_steps and scores[start:i].min()>=cutoff-1.:
            good[start:i]=True
    starts=np.flatnonzero(good&~np.r_[False,good[:-1]]);ends=np.flatnonzero(good&~np.r_[good[1:],False])+1
    return [(int(a),int(b)) for a,b in zip(starts,ends) if b-a>=3 and np.sum(scores[a:b]-cutoff)>0]
