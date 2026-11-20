"""Experimental observation-only missing-substructure variance model.

This is a marginal DANO model, not a full correlated Bijvoet likelihood.
Held-out anomalous observations gate its use; no target phases are accepted.
"""
import numpy as np
from scipy.optimize import minimize
from oasis_py.sad import fit_anomalous_scale


def _fit(data, q, indices, baseline_alpha):
    noise_scale=float(np.median(data.sigdano[indices]**2))
    # Normalize the scale coordinate, including when the null wins globally.
    unit=max(baseline_alpha, np.sqrt(np.mean(data.dano[indices]**2)
             /max(np.mean(2*abs(q[indices])**2),1e-30)), 1e-15)
    signal=2*abs(unit*q[indices])**2/noise_scale
    noise=data.sigdano[indices]**2/noise_scale
    observed=data.dano[indices]**2/noise_scale
    eps=data.epsilon[indices]
    def objective(x):
        variance=noise+x[0]*signal+x[1]*eps
        residual=(1-observed/variance)/variance
        return float(np.mean(np.log(variance)+observed/variance)),np.array([
            np.mean(residual*signal),np.mean(residual*eps)])
    extra=max(float(np.mean(observed-noise))/float(eps.mean()),0.)
    fits=[minimize(objective,x,jac=True,bounds=[(0,None),(0,None)],method='L-BFGS-B',
                   options={'maxiter':200,'ftol':1e-12,'gtol':1e-8})
          for x in [[(baseline_alpha/unit)**2,0.], [0.,extra], [1.,extra]]]
    fit=min((f for f in fits if f.success and np.isfinite(f.fun)),key=lambda f:f.fun,default=None)
    if fit is None:raise ArithmeticError('Anomalous variance optimization failed')
    return unit*np.sqrt(fit.x[0]),float(fit.x[1]*noise_scale)


def missing_substructure_variance(data,q,max_shells=8):
    """Return scaled anomalous SF and effective sigma with held-out selection.

    V(DANO)=sigmaD² + 2 alpha(shell)² |Q|² + epsilon*tau(shell).
    The final term covers unresolved anomalous scattering/error and cannot
    supply phases of missing sites. Select only if held-out marginal log
    likelihood gains exceed two standard errors. Fixed HKL split, not a phase
    quality oracle. Final selected parameters are refitted on all acentric rows.
    """
    if data.dano is None:raise ValueError('SAD observations required')
    q=np.asarray(q,dtype=complex)
    if q.shape!=data.f.shape or not np.all(np.isfinite(q)):raise ValueError('Invalid anomalous SF')
    use=np.flatnonzero(~data.centric)
    if len(use)<1000:raise ValueError('At least 1000 acentric observations for shell-error validation')
    h=data.hkl.astype(np.int64)
    validation=((h[:,0]*73856093 ^ h[:,1]*19349663 ^ h[:,2]*83492791)%5)==0
    train=use[~validation[use]]
    alpha_train,_=fit_anomalous_scale(data.subset(train),q[train])
    alpha_all,base=fit_anomalous_scale(data,q)
    shells=np.array_split(use[np.argsort(data.stol2[use],kind='stable')],
                          max(1,min(max_shells,len(use)//1000)))
    gains=[];cv=[]
    for ix in shells:
        tr=ix[~validation[ix]];va=ix[validation[ix]]
        if len(tr)<100 or len(va)<30:raise ValueError('Insufficient observations in validation shell')
        alpha,tau=_fit(data,q,tr,alpha_train)
        v0=data.sigdano[va]**2+2*abs(alpha_train*q[va])**2
        v1=data.sigdano[va]**2+2*abs(alpha*q[va])**2+data.epsilon[va]*tau
        gain=.5*(np.log(v0/v1)+data.dano[va]**2*(1/v0-1/v1))
        gains.extend(gain.tolist());cv.append(dict(train=len(tr),validation=len(va),
                        alpha=float(alpha),extra_variance=tau,loglikelihood_gain=float(gain.sum())))
    gains=np.asarray(gains);se=float(gains.std(ddof=1)*np.sqrt(len(gains)))
    accepted=bool(gains.sum()>2*se)
    scaled=alpha_all*q;extra=np.zeros(len(q));full=[]
    if accepted:
        # Use resolution intervals to apply fitted errors to centric rows too.
        edges=[-np.inf]+[(data.stol2[a[-1]]+data.stol2[b[0]])/2 for a,b in zip(shells[:-1],shells[1:])]+[np.inf]
        for i,ix in enumerate(shells):
            alpha,tau=_fit(data,q,ix,alpha_all)
            rows=(data.stol2>=edges[i])&(data.stol2<edges[i+1])
            scaled[rows]=alpha*q[rows];extra[rows]=tau*data.epsilon[rows]
            full.append(dict(acentric_rows=len(ix),alpha=float(alpha),extra_variance=tau,
                             stol2_lower=None if i==0 else float(edges[i]),
                             stol2_upper=None if i==len(shells)-1 else float(edges[i+1])))
    return scaled,np.sqrt(data.sigdano**2+extra),dict(
        model='heldout_shell_missing_substructure_variance_v1',accepted=accepted,
        validation_rule='fixed HKL hash one fifth; gain > two standard errors',
        validation_loglikelihood_gain=float(gains.sum()),validation_gain_standard_error=se,
        validation_shells=cv,full_shells=full,baseline_scale_fit=base,
        target_reference_used=False,
        limitation='Marginal variance diagnostic; unresolved sites are downweighted, not recovered')
