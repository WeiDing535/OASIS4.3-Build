# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Experimental NumPy port of Clipper's per-reflection MLHL target.

Derived from SFweight_spline<T>::TargetHL in Clipper 2.1 sfweight.cpp.
Copyright (C) 2000–2006 Kevin Cowtan and University of York.
Python implementation changes (2026): array operations, log-sum-exp stability,
explicit input checks and exact centric phase integration on the same grid.

This file is LGPL-2.1-or-later, including Clipper's additional linking
permission. See THIRD_PARTY_NOTICES.md and licenses/clipper-LGPL.txt.
The complete corresponding upstream source URL and SHA are recorded in
docs/THIRD_PARTY.md and resources/clipper-source.json. No Parrot application code is used.

This is the likelihood kernel only: it does not reproduce Clipper's spline
scaling/optimization or Parrot's map modification, masks and map simulation.
No target complete coordinates or reference phases are accepted.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import logsumexp


@dataclass
class HLTarget:
    value: np.ndarray
    ds: np.ndarray
    dw: np.ndarray
    dss: np.ndarray
    dww: np.ndarray
    hl: np.ndarray
    phase: np.ndarray
    fom: np.ndarray


def reflection_classes(data):
    """Gemmi epsilon, centricity and allowed phase in Clipper's convention."""
    centric=data.centric
    allowed=np.zeros(len(data.f))
    for i in np.flatnonzero(centric):
        h=data.hkl[i].tolist()
        op=next(op for op in data.spacegroup.operations()
                if op.apply_to_hkl(h)==[-v for v in h])
        allowed[i]=(-op.phase_shift(h)/2)%np.pi
    return data.epsilon,centric,allowed


def hl_target(fo,sigfo,fc,prior_hl,scale,error,epsilon,centric,allowed_phase,
              samples=24,chunk=512):
    """Per-reflection negative log likelihood, derivatives and phase moments.

    d = 2 sigmaFo² + epsilon_c * error; epsilon_c=epsilon*(2 if centric).
    logP(phi) = prior_HL(phi) + 2 Fo scale |Fc| cos(phi-argFc)/d.
    NLL = cf*log(d) + (Fo²+scale²|Fc|²)/d - log(mean_phi(exp(logP))),
    where cf=1 acentric, 1/2 centric. The latter average uses the two allowed
    phases. The prior is unnormalized, matching Clipper's TargetHL; its
    normalization is independent of scale/error and hence their optimization.

    samples must be a multiple of 24, so allowed crystallographic centric
    phases lie exactly on Clipper's angular grid. No cross derivative is
    returned: upstream TargetHL leaves dsw undefined.
    """
    fo=np.asarray(fo,dtype=float);n=len(fo)
    if fo.ndim!=1 or n==0:raise ValueError('Expected a nonempty amplitude vector')
    def vector(value,dtype=float):
        try:return np.broadcast_to(np.asarray(value,dtype=dtype),(n,))
        except ValueError as exc:raise ValueError('Inconsistent likelihood array shapes') from exc
    sig=vector(sigfo);fc=vector(fc,complex);s=vector(scale);w=vector(error)
    eps=vector(epsilon);cen=vector(centric,bool);allowed=vector(allowed_phase)
    prior=np.asarray(prior_hl,dtype=float)
    if prior.shape!=(n,4):raise ValueError('HL coefficients must have shape (n,4)')
    if (any(not np.all(np.isfinite(x)) for x in [fo,sig,fc,s,w,eps,allowed,prior])
            or np.any(fo<0) or np.any(sig<=0) or np.any(s<0) or np.any(w<0)
            or np.any(eps<1)):
        raise ValueError('Invalid likelihood inputs, sigmas, scale or error variance')
    if not isinstance(samples,int) or samples<24 or samples%24 or not isinstance(chunk,int) or chunk<1:
        raise ValueError('Use an angular grid multiple of 24 and a positive chunk size')
    acentric_theta=np.arange(samples)*2*np.pi/samples
    # Upstream rounds centric phase to its grid. Reject inconsistent supplied
    # restrictions rather than silently changing a caller's space-group setting.
    if np.any(abs(np.sin(samples*allowed[cen]/2))>1e-7):
        raise ValueError('Centric phase restriction is not on the angular grid')
    epc=eps*np.where(cen,2.,1.);cf=np.where(cen,.5,1.)
    variance=2*sig**2+epc*w
    amplitude=abs(fc);angle=np.angle(fc)
    xs=2*fo*amplitude/variance
    field=s*xs*np.exp(1j*angle)
    hl=prior.copy();hl[:,0]+=field.real;hl[:,1]+=field.imag
    value=np.empty(n);ds=np.empty(n);dw=np.empty(n)
    dss=np.empty(n);dww=np.empty(n);phase=np.empty(n);fom=np.empty(n)
    # Separate classes avoid allocating n x samples data for centric rows.
    for flag in [False,True]:
        indices=np.flatnonzero(cen==flag)
        for start in range(0,len(indices),chunk):
            ix=indices[start:start+chunk]
            theta=(allowed[ix,None]+np.array([0.,np.pi])[None,:]
                   if flag else np.broadcast_to(acentric_theta,(len(ix),samples)))
            c=np.cos(theta);sn=np.sin(theta)
            logp=hl[ix,0,None]*c+hl[ix,1,None]*sn+hl[ix,2,None]*np.cos(2*theta)+hl[ix,3,None]*np.sin(2*theta)
            norm=logsumexp(logp,axis=1)
            prob=np.exp(logp-norm[:,None])
            q1=xs[ix,None]*np.cos(theta-angle[ix,None]);q2=s[ix,None]*q1
            mean1=np.sum(prob*q1,axis=1);mean2=np.sum(prob*q2,axis=1)
            v=variance[ix];ec=epc[ix];sf2=fo[ix]**2+s[ix]**2*amplitude[ix]**2
            value[ix]=cf[ix]*np.log(v)+sf2/v-norm+np.log(theta.shape[1])
            ds[ix]=2*s[ix]*amplitude[ix]**2/v-mean1
            dw[ix]=ec*(cf[ix]/v-sf2/v**2)+ec/v*mean2
            dss[ix]=2*amplitude[ix]**2/v-np.sum(prob*q1**2,axis=1)+mean1**2
            dww[ix]=ec**2*(-cf[ix]/v**2+2*sf2/v**3)-(ec/v)**2*(np.sum(prob*(2+q2)*q2,axis=1)-mean2**2)
            z=np.sum(prob*(c+1j*sn),axis=1)
            phase[ix]=np.angle(z)%(2*np.pi);fom[ix]=np.clip(abs(z),0,1)
    return HLTarget(value,ds,dw,dss,dww,hl,phase,fom)
