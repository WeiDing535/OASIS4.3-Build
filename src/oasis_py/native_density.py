# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Independent experimental classical density modification, inspired by Cowtan 2010.

Uses NumPy/SciPy and Gemmi FFT/symmetry; never executes Parrot. This is NOT
an implementation-equivalent Parrot port: histogram simulation, shell MLHL
model and perturbation correction are approximations. Optional anisotropic
Wilson correction and explicit-operator NCS averaging require validation.
"""
from pathlib import Path
import numpy as np
import gemmi
from scipy.fft import fftn, ifftn, fftfreq
from scipy.optimize import minimize, brentq
from scipy.special import logsumexp, expit, i0e, i1e
from .phasing import PhaseResult
from .sad import posterior_moments


def metric_filter(cell, shape, radius):
    """Periodic Gaussian convolution in Cartesian Angstroms, including skew cells."""
    axes=[fftfreq(n)*n for n in shape]
    h,k,l=np.meshgrid(*axes,indexing='ij',sparse=True)
    b=np.array(cell.frac.mat);g=b@b.T
    s2=g[0,0]*h*h+g[1,1]*k*k+g[2,2]*l*l+2*g[0,1]*h*k+2*g[0,2]*h*l+2*g[1,2]*k*l
    return np.exp(-2*np.pi**2*radius**2*s2)


def smooth(values, kernel):
    return ifftn(fftn(values)*kernel).real


def protein_mask(grid, solvent_fraction, kernel):
    rho=np.array(grid.array,dtype=float)
    variance=np.maximum(smooth(rho*rho,kernel)-smooth(rho,kernel)**2,0)
    # Explicit crystallographic averaging protects screw axes and skew grids.
    grid.array[:]=variance;grid.symmetrize_avg();variance=np.array(grid.array,dtype=float)
    width=max(float(variance.std())*.04,1e-12)
    lo,hi=float(variance.min()-30*width),float(variance.max()+30*width)
    for _ in range(45):
        threshold=(lo+hi)/2
        if expit((variance-threshold)/width).mean()>1-solvent_fraction:lo=threshold
        else:hi=threshold
    grid.array[:]=rho
    return expit((variance-(lo+hi)/2)/width)


def histogram_transform(rho, mask, target_quantiles):
    """Smooth monotonic empirical histogram map plus solvent flattening."""
    ordered=rho[mask>.5]
    if ordered.size<32:raise ValueError('Too few protein grid points for histogram matching')
    q=np.linspace(0,1,len(target_quantiles))
    source=np.quantile(ordered,q)
    # Standardize target independently; preserve working-map protein variance.
    target=np.asarray(target_quantiles)
    target=(target-target.mean())/max(target.std(),1e-12)
    target=target*ordered.std()+ordered.mean()
    source,ix=np.unique(source,return_index=True)
    matched=np.interp(rho,source,target[ix])
    # Zero Fourier origin is arbitrary: flatten to the current solvent mean.
    solvent=float(np.sum((1-mask)*rho)/max(np.sum(1-mask),1e-12))
    return mask*matched+(1-mask)*solvent, (source,target[ix],solvent)


def apply_transform(rho, mask, mapping):
    source,target,solvent=mapping
    return mask*np.interp(rho,source,target)+(1-mask)*solvent


def reference_histogram(data, fom, reference_mtz, solvent_fraction=.5, max_grid_points=8_000_000):
    """Generic reference radial amplitude/FOM matching; no target model accepted.

    The reference ordered region is estimated by local variance rather than
    an atomic mask. Bin medians replace Parrot's spline scaling/simulation.
    """
    from .reference_data import read_reference
    m=read_reference(reference_mtz)
    h=m.make_miller_array();s=m.cell.calculate_1_d2_array(h)
    f=np.array(m.column_with_label('FC_BEST.F_phi.F').array,dtype=float)
    phi=np.deg2rad(np.array(m.column_with_label('FC_BEST.F_phi.phi').array,dtype=float))
    work_s=data.stol2*4;keep=(s>0)&(s<=work_s.max())&np.isfinite(f)&(f>0)&np.isfinite(phi)
    h,s,f,phi=h[keep],s[keep],f[keep],phi[keep]
    if len(f)<100:raise ValueError('Reference has fewer than 100 valid coefficients in the working resolution range')
    ops=m.spacegroup.operations();eps=ops.epsilon_factor_array(h)
    edges=np.unique(np.quantile(work_s,np.linspace(0,1,13)))
    coeff=f.copy()
    for lower,upper in zip(edges[:-1],edges[1:]):
        w=(work_s>=lower)&(work_s<=upper);r=(s>=lower)&(s<=upper)
        if not w.any() or not r.any():continue
        target=np.sqrt(np.mean(data.f[w]**2/data.epsilon[w]))
        observed=np.sqrt(np.mean(f[r]**2/eps[r]))
        coeff[r]*=target/max(observed,1e-12)*float(np.mean(fom[w]))
    asu=gemmi.ComplexAsuData(m.cell,m.spacegroup,np.ascontiguousarray(h,dtype=np.int32),np.ascontiguousarray(coeff*np.exp(1j*phi),dtype=np.complex64))
    shape=asu.get_size_for_hkl(sample_rate=3.)
    if np.prod(shape)>max_grid_points:raise ValueError('Reference map exceeds grid memory limit')
    grid=asu.transform_f_phi_to_map(exact_size=shape)
    mask=protein_mask(grid,solvent_fraction,metric_filter(m.cell,shape,4.))
    return np.quantile(np.array(grid.array)[mask>.5],np.linspace(0,1,513))


def fit_shell_likelihood(data, fc, prior_hl, samples=72):
    """Approximate MLHL: integrate complex-Gaussian error against fixed prior.

    Sigma = epsilon*tau + sigF^2; D and tau fitted per resolution shell.
    Centric rows excluded from fit (their discrete moments handled separately).
    Quadrature approximates the acentric likelihood; no claim of calibrated FOM.
    """
    n=len(data.f);theta=np.arange(samples)*2*np.pi/samples
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    # Normalizing each reflection also conditions the optimization numerically.
    scale=np.sqrt(np.mean(data.f**2/data.epsilon));fo=data.f/scale;calc=fc/scale
    sig=data.sigf/scale;eps=data.epsilon
    field=np.zeros(n,complex);history=[]
    order=np.argsort(data.stol2,kind='stable')
    for ix in np.array_split(order,max(1,min(12,n//250))):
        use=ix[~data.centric[ix]]
        if len(use)<20:continue
        f=fo[use];c=calc[use];e=eps[use];sigma=sig[use]**2
        lp=prior_hl[use]@basis;lp-=logsumexp(lp,axis=1)[:,None]
        cos=np.cos(theta[None,:]-np.angle(c)[:,None]);a=np.abs(c)
        def objective(x):
            d,t=np.exp(x);v=e*t+sigma
            delta=f[:,None]**2+(d*a[:,None])**2-2*d*(f*a)[:,None]*cos
            prob=lp-delta/v[:,None];norm=logsumexp(prob,axis=1)
            weights=np.exp(prob-norm[:,None])
            value=np.mean(np.log(v)-norm)
            grad_d=np.mean(np.sum(weights*(2*d*a[:,None]**2-2*(f*a)[:,None]*cos),axis=1)/v)*d
            grad_t=np.mean(e/v-e*np.sum(weights*delta,axis=1)/v**2)*t
            return value,np.array([grad_d,grad_t])
        d0=np.clip(np.sqrt(np.mean(f*f)/max(np.mean(a*a),1e-9)),.01,100)
        fitted=minimize(objective,np.log([d0,.5]),jac=True,method='L-BFGS-B',bounds=[(-7,7),(-9,5)],options={'maxiter':60})
        d,t=np.exp(fitted.x);v=eps[ix]*t+sig[ix]**2
        field[ix]=2*d*fo[ix]*calc[ix]/v
        history.append({'nfit':len(use),'D':float(d),'tau':float(t),'converged':bool(fitted.success)})
    # Conservative cap is a research stabilizer, not Parrot's calibration.
    field*=np.minimum(1,20/np.maximum(np.abs(field),1e-12))
    return field,history


def initial_prior(data, initial):
    """Keep consistent HL; otherwise match the first moment by maximum entropy.

    Legacy P+ HL represents original observations, not final PHIB/FOM. A
    von Mises reconstruction preserves final phase/FOM but loses branch shape.
    Centric rows use a discrete two-state reconstruction instead.
    """
    phase=np.asarray(initial.phase,dtype=float);fom=np.asarray(initial.fom,dtype=float)
    if phase.shape!=(len(data.f),) or fom.shape!=phase.shape or not np.all(np.isfinite(phase)) or not np.all(np.isfinite(fom)) or np.any((fom<0)|(fom>1)):
        raise ValueError('Invalid initial phase/FOM arrays')
    if initial.hl is not None:
        hp,hf=posterior_moments(data,initial.hl)
        mismatch=float(np.max(np.abs(hf*np.exp(1j*hp)-fom*np.exp(1j*phase))))
        # MTZ stores single precision. Large HL amplify roundoff; accept a
        # documented absolute moment tolerance, without altering likelihood.
        if mismatch<1e-4:return np.array(initial.hl,dtype=float,copy=True),'original_consistent_HL',mismatch
        if initial.metadata.get('method')=='sad_likelihood':
            raise ValueError('Initial SAD HL and PHIB/FOM disagree beyond MTZ rounding tolerance')
    else:mismatch=None
    kappa=np.array([0. if m<1e-12 else brentq(lambda k:i1e(k)/i0e(k)-min(m,.9999),0,10000) for m in fom])
    kappa[data.centric]=np.arctanh(np.minimum(fom[data.centric],.9999))
    hl=np.zeros((len(phase),4));hl[:,0]=kappa*np.cos(phase);hl[:,1]=kappa*np.sin(phase)
    return hl,'moment_matched_von_Mises_branch_shape_not_preserved',mismatch


def working_grid_size(asu, max_grid_points):
    """Select only from oversampled, symmetry-compatible FFT grids."""
    for rate in (3.,2.5,2.1):
        shape=asu.get_size_for_hkl(sample_rate=rate)
        if np.prod(shape)<=max_grid_points:return shape,rate
    raise ValueError('Working map exceeds grid memory limit even at sample rate 2.1')


def _density_modify(data, initial, solvent_fraction, reference_mtz, cycles=5,
                    seed=20261016, max_grid_points=8_000_000, ncs_operators=None,
                    _alternate_initial=None, _pattern=None, _cycle_callback=None):
    """Independent dual-space phase improvement with a fixed experimental prior.

    Accepts reflections, a PhaseResult and a GENERIC reference MTZ only.
    No target complete coordinates/phases, external executable or subprocess.
    """
    if not 0<solvent_fraction<1 or not isinstance(cycles,int) or cycles<1:
        raise ValueError('Invalid native density modification settings')
    prior,prior_kind,input_mismatch=initial_prior(data,initial)
    # Experimental alternating engine: fixed priors are prepared ONCE. The
    # previous posterior supplies the next map, never the next prior.
    priors=None
    if _pattern is not None:
        if _pattern not in ('A','B','AB','BA') or _alternate_initial is None:
            raise ValueError('Alternation requires A/B priors and pattern A, B, AB or BA')
        other,other_kind,other_mismatch=initial_prior(data,_alternate_initial)
        priors={'A':prior,'B':other}
        # The caller always passes A as initial and B as alternate_initial.
        prior=priors[_pattern[0]]
    phase,fom=posterior_moments(data,prior)
    asu=gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),np.ascontiguousarray(data.f*np.exp(1j*phase),dtype=np.complex64))
    shape,sample_rate=working_grid_size(asu,max_grid_points)
    kernel=metric_filter(data.cell,shape,4.)
    rng=np.random.default_rng(seed);history=[];hl=prior.copy()
    target=reference_histogram(data,fom,reference_mtz,max_grid_points=max_grid_points)
    for cycle in range(cycles):
        if priors is not None:prior=priors[_pattern[cycle%len(_pattern)]]
        values=data.f*fom*np.exp(1j*phase)
        asu=gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),np.ascontiguousarray(values,dtype=np.complex64))
        grid=asu.transform_f_phi_to_map(exact_size=shape)
        rho=np.array(grid.array,dtype=float);mask=protein_mask(grid,solvent_fraction,kernel)
        if ncs_operators:
            from .ncs import prepare_ncs,apply_ncs
            plans,ncs_diagnostics=prepare_ncs(grid,mask,kernel,ncs_operators)
            averaged=apply_ncs(rho,data.cell,plans)
        else:plans=[];ncs_diagnostics=[];averaged=rho
        modified,mapping=histogram_transform(averaged,mask,target)
        # Real-space Hutchinson derivative estimate, fixed mask/mapping.
        # Symmetrize noise before estimating the response on the crystal subspace.
        noise=rng.normal(size=shape);grid.array[:]=noise;grid.symmetrize_avg()
        noise=np.array(grid.array,dtype=float);noise-=noise.mean()
        noise*=max(rho.std(),1e-12)*.01/max(noise.std(),1e-12)
        perturbation=rho+noise
        if plans:perturbation=apply_ncs(perturbation,data.cell,plans)
        perturbed=apply_transform(perturbation,mask,mapping)
        gamma=float(np.clip(np.sum(noise*(perturbed-modified))/max(np.sum(noise*noise),1e-30),0,.95))
        grid.array[:]=modified-gamma*rho;grid.symmetrize_avg()
        fc=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(data.hkl).astype(complex)
        field,fit=fit_shell_likelihood(data,fc,prior)
        updated=prior.copy();updated[:,0]+=field.real;updated[:,1]+=field.imag
        # Replace the DM evidence each cycle; never accumulate initial prior.
        hl=updated;phase,fom=posterior_moments(data,hl)
        history.append({'cycle':cycle+1,'gamma':gamma,'mean_fom':float(fom.mean()),'protein_volume_fraction':float(mask.mean()),'shell_fit':fit,'ncs':ncs_diagnostics})
        if priors is not None:history[-1]['active_prior']=_pattern[cycle%len(_pattern)]
        if _cycle_callback is not None:
            _cycle_callback(cycle+1,phase.copy(),fom.copy())
    result=PhaseResult(phase,fom,None,hl,{'method':'native_classical_dm_experimental','backend':'NumPy/SciPy/Gemmi; no Parrot execution','initial_method':initial.metadata.get('method','legacy_pplus'),'initial_scale_fit_version':initial.metadata.get('anomalous_scale_fit_version','legacy_bounded'),'initial_prior':prior_kind,'input_HL_moment_mismatch':input_mismatch,'cycles':history,'seed':seed,'solvent_fraction':solvent_fraction,'grid_shape':shape,'grid_sample_rate':sample_rate,'reference_mtz':str(Path(reference_mtz).resolve()),'target_complete_model_used':False,'prior_repeatedly_accumulated':False,'limitations':['NCS requires explicit operators; automatic discovery not implemented','Approximate radial histogram scaling and variance mask','Approximate shell MLHL and global perturbation correction; FOM not calibrated','Generic histogram simulated once using initial FOM, not every cycle','Not a Parrot reproduction; experimental validation required']})
    if priors is not None:
        result.metadata.update(method='alternating_fixed_prior_native_dm_experimental',
            pattern=_pattern, initial_method=initial.metadata.get('method') if _pattern[0]=='A' else _alternate_initial.metadata.get('method'),
            prior_A={'kind':prior_kind,'moment_mismatch':input_mismatch,'method':initial.metadata.get('method')},
            prior_B={'kind':other_kind,'moment_mismatch':other_mismatch,'method':_alternate_initial.metadata.get('method')},
            initial_scale_fit_version=_alternate_initial.metadata.get('anomalous_scale_fit_version','legacy_bounded'),
            alternation_scope='Fixed P+ / continuous SAD priors in shared-map density modification; no triplet re-estimation',
            histogram_initialization='First active prior; AB shares initialization with A, BA with B',
            final_active_prior=_pattern[(cycles-1)%len(_pattern)])
    return result


def density_modify_native(data, initial, solvent_fraction, reference_mtz, cycles=5,
                          seed=20261016, max_grid_points=8_000_000,
                          anisotropy=False, ncs_operators=None):
    """Native density modification with optional anisotropy and explicit NCS."""
    if anisotropy:
        from .anisotropy import fit_anisotropy,corrected_reflections
        factors,anisotropy_meta=fit_anisotropy(data)
        data=corrected_reflections(data,factors)
    else:anisotropy_meta={'status':'disabled'}
    if ncs_operators:
        from .ncs import validate_operators
        validate_operators(ncs_operators)
    result=_density_modify(data,initial,solvent_fraction,reference_mtz,
        cycles=cycles,seed=seed,max_grid_points=max_grid_points,ncs_operators=ncs_operators)
    result.metadata['anisotropy']=anisotropy_meta
    result.metadata['ncs_operator_count']=len(ncs_operators or [])
    return result
