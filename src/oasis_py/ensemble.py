# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""First-moment model averaging, never multiplication of shared SAD evidence."""
import numpy as np
from .phasing import PhaseResult


def mix_phases(first,second,weight=.5):
    """Mixture p(phi)=weight*p1(phi)+(1-weight)*p2(phi).

    Components MUST describe identical HKL rows, cell and space-group setting.
    The exact mixture first moment is available without an HL approximation.
    Output omits HL and PPLUS; it is a final-phase/map ensemble, not independent
    likelihood evidence to multiply back into either original posterior.
    """
    if not np.isfinite(weight) or not 0<=weight<=1:raise ValueError('Mixture weight must lie in [0,1]')
    if first.phase.shape!=second.phase.shape or first.fom.shape!=first.phase.shape or second.fom.shape!=second.phase.shape:raise ValueError('Mixture component shapes disagree')
    for result in (first,second):
        if not np.isfinite(result.phase).all() or not np.isfinite(result.fom).all() or np.any((result.fom<0)|(result.fom>1)):raise ValueError('Invalid mixture phase/FOM')
    z=weight*first.fom*np.exp(1j*first.phase)+(1-weight)*second.fom*np.exp(1j*second.phase)
    return PhaseResult(np.angle(z)%(2*np.pi),np.clip(abs(z),0,1),None,None,
        {'method':'phase_mixture','weight_first':weight,'components':[first.metadata.get('method','legacy_pplus'),second.metadata.get('method','unknown')],'shared_SAD_likelihood_multiplied':False,'HL_omitted':True,'PPLUS_omitted':True,'component_FOM_calibrated':False,'requires':'same HKL order and crystallographic setting; separate validation before iterative feedback'})


def phase_sad_ensemble(data,composition,heavy,fpp,cycles=5,kmin=.03,
                       max_neighbors=600,max_per_target=96,triplet_weight=.1,
                       aoe=.5,density_options=None):
    """Run two SAD hypotheses; optionally modify each map before mixing.

    Equal mixing weight is fixed, not chosen using a target reference model.
    With density_options the sequence is P+ -> DM and continuous SAD -> DM,
    then model averaging. It is not DM of a moment-reconstructed raw mixture.
    """
    from .phasing import phase_native
    from .sad import phase_sad_likelihood
    first=phase_native(data,'sad',composition,heavy=heavy,fpp=fpp,cycles=cycles,
        aoe=aoe,kmin=kmin,max_neighbors=max_neighbors,max_per_target=max_per_target)
    second=phase_sad_likelihood(data,composition,heavy,fpp,cycles=cycles,
        kmin=kmin,max_neighbors=max_neighbors,max_per_target=max_per_target,
        triplet_weight=triplet_weight)
    original_metadata=[first.metadata,second.metadata]
    if density_options is not None:
        from .native_density import density_modify_native
        first=density_modify_native(data,first,**density_options)
        second=density_modify_native(data,second,**density_options)
    result=mix_phases(first,second)
    result.metadata.update(ensemble_stage='after_separate_native_dm' if density_options is not None else 'initial_phases',
        initial_component_metadata=original_metadata,
        final_component_metadata=[first.metadata,second.metadata],
        uses_target_reference_phases=False,
        initial_scale_fit_version=original_metadata[1]['anomalous_scale_fit_version'])
    return result


def mix_mtz_files(first_path,second_path,output_path):
    """Mix standard PHIB/FOM MTZ outputs only after matching HKL and data."""
    from pathlib import Path
    import json,gemmi
    from .data import read_mtz,write_mtz
    output_path=Path(output_path);sidecar=output_path.with_suffix('.json')
    if output_path.exists() or sidecar.exists():raise ValueError('Mixture output exists; choose a new path')
    first_data,_=read_mtz(first_path,'dmr',{'f':'FP','sigf':'SIGFP'})
    second_data,_=read_mtz(second_path,'dmr',{'f':'FP','sigf':'SIGFP'})
    if first_data.spacegroup.operations()!=second_data.spacegroup.operations() or not np.allclose(first_data.cell.parameters,second_data.cell.parameters,atol=1e-3,rtol=1e-5):raise ValueError('Mixture MTZ files must have matching cell and space-group setting')
    lookup={tuple(h):i for i,h in enumerate(second_data.hkl)}
    if len(lookup)!=len(first_data.hkl) or any(tuple(h) not in lookup for h in first_data.hkl):raise ValueError('Mixture MTZ files must contain the same HKL set')
    ix=np.array([lookup[tuple(h)] for h in first_data.hkl])
    if not np.allclose(first_data.f,second_data.f[ix],atol=1e-5,rtol=1e-5) or not np.allclose(first_data.sigf,second_data.sigf[ix],atol=1e-5,rtol=1e-5):raise ValueError('Mixture MTZ files must contain the same observed amplitudes and sigmas')
    def read_phase(path,data,hkl):
        m=gemmi.read_mtz_file(str(path));m.ensure_asu();indices={tuple(h):i for i,h in enumerate(m.make_miller_array())};order=[indices[tuple(h)] for h in hkl]
        phase_column=m.column_with_label('PHIB');fom_column=m.column_with_label('FOM')
        if phase_column is None or fom_column is None or phase_column.type!='P' or fom_column.type not in {'W','R'}:raise ValueError('Mixture requires standard PHIB phase and FOM weight columns')
        phase=np.deg2rad(phase_column.array[order].astype(float));fom=fom_column.array[order].astype(float)
        return PhaseResult(phase,fom,None,None,{'method':'supplied_mtz','source':str(Path(path).resolve())})
    first=read_phase(first_path,first_data,first_data.hkl);second=read_phase(second_path,second_data,first_data.hkl)
    result=mix_phases(first,second);result.metadata['sources']=[str(Path(p).resolve()) for p in [first_path,second_path]]
    output_path.parent.mkdir(parents=True,exist_ok=True);write_mtz(output_path,first_data,result.phase,result.fom)
    sidecar.write_text(json.dumps(result.metadata,indent=2,allow_nan=False)+'\n')
    return result.metadata
