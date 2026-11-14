"""Opt-in DM experiment separating an initial map from a fixed phase prior.

The numerical DM steps reuse preserved upstream primitives. P+ supplies only
initial map phase/FOM; original SAD HL supplies every cycle's likelihood.
No GUI/default changes, target protein model or automatic parameter search.
"""
from pathlib import Path
import hashlib
import numpy as np
import gemmi
from oasis_py import native_density as legacy
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from oasis_py.clipper_likelihood import reflection_classes


def _map_arrays(data,result):
    phase=np.array(result.phase,dtype=float,copy=True);fom=np.array(result.fom,dtype=float,copy=True)
    if phase.shape!=(len(data.f),) or fom.shape!=phase.shape or not np.isfinite(phase).all() or not np.isfinite(fom).all() or np.any((fom<0)|(fom>1)):
        raise ValueError('Invalid map phase/FOM')
    _,centric,allowed=reflection_classes(data)
    use=centric&(fom>1e-7)
    if np.any(abs(np.sin(phase[use]-allowed[use]))>1e-5):raise ValueError('Map seed violates centric allowed phases')
    return phase,fom


def density_modify_dual_channel(data,observation_prior,map_initial,solvent_fraction,reference_mtz,
                               *,cycles=5,seed=20261024,histogram='observation',
                               max_grid_points=8_000_000,_cycle_callback=None):
    """Keep original observation HL fixed while allowing a distinct map seed.

    histogram='observation' isolates the map seed from the reference histogram.
    histogram='map' is a declared ablation initialized from the seed FOM.
    map_initial.hl is deliberately ignored; it is not another likelihood.
    """
    if not 0<solvent_fraction<1 or type(cycles) is not int or cycles<1 or histogram not in ('observation','map'):
        raise ValueError('Invalid dual-channel DM settings')
    if observation_prior.hl is None:raise ValueError('Fixed observation HL required')
    prior=np.array(observation_prior.hl,dtype=float,copy=True)
    op,of=posterior_moments(data,prior)
    supplied_phase,supplied_fom=_map_arrays(data,observation_prior)
    mismatch=float(np.max(abs(of*np.exp(1j*op)-supplied_fom*np.exp(1j*supplied_phase))))
    if mismatch>1e-5:raise ValueError('Fixed prior HL and phase/FOM disagree')
    phase,fom=_map_arrays(data,map_initial)
    fixed_identity=hashlib.sha256(data.hkl.tobytes()+prior.tobytes()).hexdigest();prior.setflags(write=False)
    asu=gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),np.ascontiguousarray(data.f*np.exp(1j*phase),dtype=np.complex64))
    shape,rate=legacy.working_grid_size(asu,max_grid_points)
    kernel=legacy.metric_filter(data.cell,shape,4.)
    histogram_fom=of if histogram=='observation' else fom
    target=legacy.reference_histogram(data,histogram_fom,reference_mtz,max_grid_points=max_grid_points)
    rng=np.random.default_rng(seed);history=[];hl=prior.copy()
    for cycle in range(cycles):
        values=data.f*fom*np.exp(1j*phase)
        asu=gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),np.ascontiguousarray(values,dtype=np.complex64))
        grid=asu.transform_f_phi_to_map(exact_size=shape)
        rho=np.array(grid.array,dtype=float);mask=legacy.protein_mask(grid,solvent_fraction,kernel)
        modified,mapping=legacy.histogram_transform(rho,mask,target)
        noise=rng.normal(size=shape);grid.array[:]=noise;grid.symmetrize_avg()
        noise=np.array(grid.array,dtype=float);noise-=noise.mean()
        noise*=max(rho.std(),1e-12)*.01/max(noise.std(),1e-12)
        perturbed=legacy.apply_transform(rho+noise,mask,mapping)
        gamma=float(np.clip(np.sum(noise*(perturbed-modified))/max(np.sum(noise*noise),1e-30),0,.95))
        grid.array[:]=modified-gamma*rho;grid.symmetrize_avg()
        fc=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(data.hkl).astype(complex)
        field,fit=legacy.fit_shell_likelihood(data,fc,prior)
        hl=prior.copy();hl[:,0]+=field.real;hl[:,1]+=field.imag
        phase,fom=posterior_moments(data,hl)
        history.append({'cycle':cycle+1,'gamma':gamma,'mean_fom':float(fom.mean()),
                        'protein_volume_fraction':float(mask.mean()),'shell_fit':fit,
                        'fixed_observation_prior_sha256':fixed_identity})
        if _cycle_callback is not None:_cycle_callback(cycle+1,phase.copy(),fom.copy())
    return PhaseResult(phase,fom,None,hl,{
        'method':'dual_channel_native_dm_research','initial_map_method':map_initial.metadata.get('method'),
        'map_initial_HL_used_as_likelihood':False,'fixed_prior_method':observation_prior.metadata.get('method'),
        'fixed_prior_sha256':fixed_identity,'fixed_prior_moment_mismatch':mismatch,
        'prior_repeatedly_accumulated':False,'histogram_initialization':histogram,
        'cycles':history,'seed':seed,'solvent_fraction':solvent_fraction,'grid_shape':shape,
        'grid_sample_rate':rate,'reference_mtz':str(Path(reference_mtz).resolve()),
        'target_complete_model_used':False,'ncs_operator_count':0,'anisotropy':{'status':'disabled'},
        'limitations':['Experimental fork of approximate native DM primitives, not a reproduction of any external density modifier',
                       'Phase FOM is not calibrated; fixed generic histogram and approximate shell likelihood',
                       'No NCS or anisotropy in this controlled research comparison']})
