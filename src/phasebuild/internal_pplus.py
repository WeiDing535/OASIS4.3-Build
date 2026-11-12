"""Installation-free phasing/DM adapter: NumPy, SciPy and Gemmi only.

This uses OASIS' approximate anomalous-difference likelihood, not Phaser's
complete SAD likelihood. It must be validated independently of the external
Phaser/RESOLVE development branch. Building uses separately bundled engines.
"""
from pathlib import Path
import json
import os
import gemmi
import numpy as np
from .trust_workflow import dump,sample_info
from .sad_bridge import FLAG,columns
from .data import fingerprint


def preflight(cfg):
    if cfg.get('observation_prior'):raise ValueError('内部 P+ 不接受外部 Phaser 先验；请清空原始似然文件')
    if cfg.get('labels'):raise ValueError('内部 P+ 当前要求标准 F(+)/SIGF(+)/F(-)/SIGF(-) 列名')
    sad=gemmi.read_mtz_file(cfg['sad'])
    columns(sad,{'F(+)':'G','SIGF(+)':'L','F(-)':'G','SIGF(-)':'L'})
    sample_info(cfg,sad)


def phase(cfg,output):
    from .sad_input import read_sad_pairs
    from .oasis_options import pairs
    from .pplus_experimental import prepare_pplus_normalization
    from .pplus_trust import phase_trusted_pplus
    from oasis_py.data import read_atoms,structure_factors,write_mtz
    from oasis_py.sad import fit_anomalous_scale,posterior_moments
    from oasis_py.phasing import sad_hl
    from oasis_py.relations import build_relations
    out=Path(output);data,single,input_meta=read_sad_pairs(cfg['prepared_sad'],free_label=FLAG,free_value=cfg['free_value'])
    seq=sample_info(cfg,data);heavy=read_atoms(cfg['heavy'],data,project_anisotropic=True)
    fpp={gemmi.Element(k).name:float(v) for k,v in pairs(cfg['fpp']).items()}
    raw=structure_factors(data,heavy,fpp);alpha,fit=fit_anomalous_scale(data,raw)
    sf=alpha*raw;error=.1;sigma=np.hypot(data.sigdano,error*abs(sf))
    hl=sad_hl(np.angle(sf),abs(sf),data.dano,sigma);hl[single]=0
    original_phase,original_fom=posterior_moments(data,hl)
    write_mtz(out/'observations.mtz',data,original_phase,original_fom,hl)
    original_meta={'method':'internal_anomalous_difference_observation_likelihood','frontend':'OASIS NumPy/SciPy/Gemmi; no Phaser',
        'approximation':'Anomalous-difference Gaussian likelihood with fitted scale; not Phaser complete SAD likelihood',
        'role':'Fixed original observation likelihood; no structural/triplet or protein-model terms',
        'anomalous_scale':float(alpha),'model_error_fraction':error,'fit':fit,'input':input_meta,
        'target_protein_coordinates_used':False,'external_phase_file_used':False}
    dump(out/'observations.json',original_meta)
    # Use precisely the stored float32 HL for branch integration and later DM.
    from .sad_input import read_phase_prior
    persisted=read_phase_prior(out/'observations.mtz',data)
    norm=prepare_pplus_normalization(data,seq['asu_composition'])
    graph=build_relations(data,norm.e,norm.kappa,cfg['kmin'],cfg['max_neighbors'],cfg['max_per_target'])
    result=phase_trusted_pplus(data,persisted.hl,np.angle(raw),graph,samples=1440,kl_budget=cfg['trust_kl'],cycles=100,damping=.25)
    if result.metadata['mean_abs_probability_change']<=1e-12:raise RuntimeError('Internal P+ has no nonzero branch update')
    seed=gemmi.Mtz(with_base=True);seed.spacegroup=data.spacegroup;seed.set_cell_for_all(data.cell);seed.add_dataset('internal_Pplus_seed_no_HL')
    for name,kind in [('PHIB','P'),('FOM','W'),('PPLUS','R')]:seed.add_column(name,kind)
    seed.set_data(np.column_stack([data.hkl,np.rad2deg(result.phase),result.fom,result.probability]).astype(np.float32));seed.write_to_file(str(out/'experimental.mtz'))
    dump(out/'experimental.json',{**result.metadata,'phase_method':'pplus_internal','frontend':original_meta['frontend'],
        'frontend_approximation':original_meta['approximation'],'sequence_asu':seq,'history':result.history,
        'input':input_meta,'normalization':norm.statistics,'relations':graph.statistics,'original_prior':fingerprint(out/'observations.mtz'),
        'Phenix_executed':False,'Phaser_input_used':False,'quality_claim':'Independent experimental route; no transfer of SEC17 external-branch Rfree claims'})


def seeded_dm(data,observation,seed_phase,seed_fom,solvent,reference,*,cycles,random_seed,trace_dir=None,
              combine='mlhl',mask_radius=4.0,ladder=None,gamma_scale=1.0,mask_kind='variance',
              wang_radius=None,wang_edge=1.0,map_weight=1.0,oracle=None,prior_weight=1.0):
    """Frozen original HL + separate initial P+ map, never posterior-as-data."""
    from oasis_py.density_mlhl import fit_shells
    from oasis_py.native_density import (working_grid_size,metric_filter,reference_histogram,
        protein_mask,wang_mask,histogram_transform,apply_transform)
    from oasis_py.sad import posterior_moments
    from oasis_py.phasing import PhaseResult
    if combine not in ('mlhl','map','normalized','fwt'):raise ValueError('Invalid DM combination '+repr(combine))
    if mask_kind not in ('variance','wang'):raise ValueError('Invalid DM mask kind '+repr(mask_kind))
    if ladder is not None and (not ladder or any(c is not None and float(c)<=0 for c in ladder)):
        raise ValueError('Invalid DM resolution ladder')
    prior=np.asarray(observation.hl,dtype=float).copy()
    if prior.shape!=(len(data.f),4) or not np.isfinite(prior).all():raise ValueError('Invalid fixed observation HL')
    phase=np.asarray(seed_phase,dtype=float);fom=np.asarray(seed_fom,dtype=float)
    if phase.shape!=(len(data.f),) or fom.shape!=phase.shape or not np.isfinite([phase,fom]).all() or np.any((fom<0)|(fom>1)):
        raise ValueError('Invalid independent P+ seed')
    coefficients=data.f*fom*np.exp(1j*phase)
    def asu(values):return gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),np.ascontiguousarray(values,dtype=np.complex64))
    shape,rate=working_grid_size(asu(coefficients),8_000_000);kernel=metric_filter(data.cell,shape,mask_radius)
    mask_override=(None if not os.environ.get('OASIS_INTERNAL_DM_MASK_NPY')
                   else np.load(os.environ['OASIS_INTERNAL_DM_MASK_NPY']).astype(float))
    if mask_override is not None and mask_override.shape!=tuple(shape):
        raise ValueError('External DM mask shape %s does not match the working grid %s'
                         %(mask_override.shape,tuple(shape)))
    target=reference_histogram(data,fom,reference);rng=np.random.default_rng(random_seed);history=[];hl=prior.copy()
    spacing=1.0/np.sqrt((data.hkl**2/np.array([data.cell.a,data.cell.b,data.cell.c])**2).sum(axis=1))
    steps=[None] if not ladder else list(ladder)
    schedule=[];remaining=cycles
    for index,cutoff in enumerate(steps):
        count=max(1,remaining) if index==len(steps)-1 else max(1,remaining//(len(steps)-index))
        schedule.extend([cutoff]*count);remaining-=count
    schedule=schedule[:cycles]
    for cycle,cutoff in enumerate(schedule):
        active=np.ones(len(data.f),bool) if cutoff is None else (spacing>=float(cutoff))
        if not active.any():continue
        current=np.where(active,coefficients,0).astype(complex)
        grid=asu(current).transform_f_phi_to_map(exact_size=shape);rho=np.array(grid.array,dtype=float)
        if mask_override is None:
            mask=(wang_mask(grid,solvent,kernel,data.cell,radius=wang_radius,edge=wang_edge)
                  if mask_kind=='wang' else protein_mask(grid,solvent,kernel))
        else:
            mask=mask_override.copy()
        modified,mapping=histogram_transform(rho,mask,target)
        noise=rng.normal(size=shape);grid.array[:]=noise;grid.symmetrize_avg();noise=np.array(grid.array,dtype=float);noise-=noise.mean()
        noise*=max(rho.std(),1e-12)*.01/max(noise.std(),1e-12)
        perturbed=apply_transform(rho+noise,mask,mapping)
        gamma=float(np.clip(np.sum(noise*(perturbed-modified))/max(np.sum(noise*noise),1e-30),0,.95))
        gamma=min(gamma*gamma_scale,.95)
        grid.array[:]=modified;grid.symmetrize_avg()
        fc_map=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(data.hkl).astype(complex)
        grid.array[:]=modified-gamma*rho;grid.symmetrize_avg()
        fc=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(data.hkl).astype(complex)
        if combine in ('map','normalized'):
            # Map iteration: the modified map itself carries the information
            # forward. 'normalized' additionally rescales each reflection to
            # the observed amplitude, which keeps |coefficients| a valid FOM.
            updated=fc_map.astype(complex) if combine=='map' else \
                data.f*np.clip(np.abs(fc_map)/np.maximum(data.f,1e-12),0,1)*np.exp(1j*np.angle(fc_map))
            fit={}
        elif combine=='fwt':
            # Standard FWT recurrence: the *modified map* is the partial model.
            # The experimental HL stays fixed, the map contributes its own
            # per-reflection likelihood (sigma_A weighted), and the model term
            # is subtracted from the 2mFo-DFc coefficients.
            field,d,fit=fit_shells(data,fc_map,prior,samples=72)
            weight=np.full(len(data.f),float(map_weight))
            if oracle is not None:
                # Diagnostic only: replace the missing per-reflection FOM with
                # the true map-phase agreement to measure what per-reflection
                # discrimination would buy.
                weight=weight*np.clip(np.cos(np.angle(fc_map)-oracle),0,1)
            hl=prior_weight*prior.copy();hl[:,0]+=weight*field.real;hl[:,1]+=weight*field.imag
            phase,fom=posterior_moments(data,hl)
            updated=2*data.f*fom*np.exp(1j*phase)-weight*d*fc_map
        else:
            field,d,fit=fit_shells(data,fc,prior,samples=72)
            hl=prior.copy();hl[:,0]+=field.real;hl[:,1]+=field.imag
            phase,fom=posterior_moments(data,hl)
            updated=2*data.f*fom*np.exp(1j*phase)-d*fc
        coefficients=np.where(active,updated,coefficients).astype(complex)
        if trace_dir is not None:
            _write_map_coefficients(Path(trace_dir)/f'cycle_{cycle+1:02d}.mtz',data,coefficients)
            if combine=='fwt':_write_map_coefficients(Path(trace_dir)/f'cycle_{cycle+1:02d}_map.mtz',data,fc_map)
            if cycle==len(schedule)-1:
                Path(trace_dir).mkdir(parents=True,exist_ok=True)
                np.save(Path(trace_dir)/'mask_last.npy',mask.astype(np.float32))
                np.save(Path(trace_dir)/'rho_last.npy',rho.astype(np.float32))
        map_fom=float(np.clip(np.abs(coefficients)/np.maximum(data.f,1e-12),0,1)[active].mean())
        row={'cycle':cycle+1,'mean_fom':float(fom.mean()),'map_fom':map_fom,'gamma':gamma,'shell_fit':fit,
             'resolution_cutoff':None if cutoff is None else float(cutoff)};history.append(row)
        print(json.dumps({'internal_DM_cycle':cycle+1,'mean_fom':row['mean_fom']}),flush=True)
    if combine!='mlhl':
        fom=np.clip(np.abs(coefficients)/np.maximum(data.f,1e-12),0,1);phase=np.angle(coefficients)
    return PhaseResult(phase,fom,None,hl,{'method':'internal_fixed_observation_MLHL_seeded_DM','backend':'NumPy/SciPy/Gemmi only; no external executable',
        'cycles':history,'grid_shape':shape,'grid_sample_rate':rate,'seed_map':'P+ first moment once','fixed_HL':'Original internal observation likelihood, never accumulated',
        'combination':combine,'map_weight':float(map_weight),'prior_weight':float(prior_weight),'external_mask':mask_override is not None,'oracle_weights':oracle is not None,
        'limitations':['Approximate mask, reference histogram, shell scale and field cap inherited from internal MLHL DM','No automatic NCS; not RESOLVE/Parrot/Phenix equivalent']}),coefficients


def _ladder_from_env(data):
    """Multi-resolution DM schedule: low-resolution shells first (Parrot-like)."""
    spec=os.environ.get('OASIS_INTERNAL_DM_LADDER','').strip()
    if not spec:return None
    d_min=float(1.0/np.sqrt(np.max((data.hkl**2/np.array([data.cell.a,data.cell.b,data.cell.c])**2).sum(axis=1))))
    if spec=='auto':
        cutoffs=[];value=8.0
        while value>d_min*1.25:cutoffs.append(round(value,2));value*=0.75
        cutoffs.append(None);return cutoffs
    values=[float(x) for x in spec.split(',') if x.strip()]
    return values+[None]


def _dm_options():
    """Ablation switches for the internal DM; all defaults reproduce shipping behaviour."""
    return {'combine':os.environ.get('OASIS_INTERNAL_DM_COMBINE','mlhl'),
        'gamma_scale':float(os.environ.get('OASIS_INTERNAL_DM_GAMMA','1.0')),
        'mask_kind':os.environ.get('OASIS_INTERNAL_DM_MASK','variance'),
        'wang_radius':(None if not os.environ.get('OASIS_INTERNAL_DM_WANG_RADIUS') else float(os.environ['OASIS_INTERNAL_DM_WANG_RADIUS'])),
        'wang_edge':float(os.environ.get('OASIS_INTERNAL_DM_WANG_EDGE','1.0')),
        'mask_radius':float(os.environ.get('OASIS_INTERNAL_DM_MASK_RADIUS','4.0')),
        'map_weight':float(os.environ.get('OASIS_INTERNAL_DM_MAP_WEIGHT','1.0'))}


def _write_map_coefficients(path,data,coefficients):
    """Debug aid: dump one DM cycle so the per-cycle map gain is measurable."""
    path.parent.mkdir(parents=True,exist_ok=True)
    m=gemmi.Mtz(with_base=True);m.spacegroup=data.spacegroup;m.set_cell_for_all(data.cell);m.add_dataset('internal_DM_cycle')
    for name,kind in [('FWT','F'),('PHWT','P')]:m.add_column(name,kind)
    m.set_data(np.column_stack([data.hkl,abs(coefficients),np.rad2deg(np.angle(coefficients))]).astype(np.float32))
    m.write_to_file(str(path))


def dm(cfg,output):
    from .sad_input import read_sad_pairs,read_phase_prior
    out=Path(output);parent=out.parent
    data,_,_=read_sad_pairs(cfg['prepared_sad'],free_label=FLAG,free_value=cfg['free_value'])
    original=read_phase_prior(parent/'02_sad/observations.mtz',data)
    seed=gemmi.read_mtz_file(str(parent/'02_sad/experimental.mtz'));lookup={tuple(h):i for i,h in enumerate(seed.make_miller_array())};ix=[lookup[tuple(h)] for h in data.hkl]
    phase=np.deg2rad(np.asarray(seed.column_with_label('PHIB'))[ix]);fom=np.asarray(seed.column_with_label('FOM'))[ix]
    seq=sample_info(cfg,data);solvent=seq['solvent_fraction'] if cfg['solvent_fraction'] is None else cfg['solvent_fraction']
    trace=os.environ.get('OASIS_INTERNAL_DM_TRACE')
    result,coeff=seeded_dm(data,original,phase,fom,solvent,cfg['reference'],cycles=cfg['dm_cycles'],random_seed=cfg['seed'],
        trace_dir=Path(trace)/Path(output).parent.name if trace else None,
        ladder=_ladder_from_env(data),**_dm_options())
    m=gemmi.Mtz(with_base=True);m.spacegroup=data.spacegroup;m.set_cell_for_all(data.cell);m.add_dataset('internal_DM')
    for name,kind in [('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHIB','P'),('FOM','W')]+[(x,'A') for x in ['HLA','HLB','HLC','HLD']]:m.add_column(name,kind)
    m.set_data(np.column_stack([data.hkl,data.f,data.sigf,abs(coeff),np.rad2deg(np.angle(coeff)),np.rad2deg(result.phase),result.fom,result.hl]).astype(np.float32));m.write_to_file(str(out/'map_phases.mtz'))
    dump(out/'map_phases.json',{'phase_method':'pplus_internal','dm_engine':'internal-MLHL-best','density_modification':result.metadata,
        'solvent_fraction_used':solvent,'solvent_source':'sequence_estimate' if cfg['solvent_fraction'] is None else 'user',
        'original_prior':fingerprint(parent/'02_sad/observations.mtz'),'Pplus_seed':fingerprint(parent/'02_sad/experimental.mtz'),
        'Phenix_executed':False,'external_programs':[]})
