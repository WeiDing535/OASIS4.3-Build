"""Installation-free phasing/DM adapter: NumPy, SciPy and Gemmi only.

This uses OASIS' approximate anomalous-difference likelihood, not Phaser's
complete SAD likelihood. It must be validated independently of the external
Phaser/RESOLVE development branch. Building uses separately bundled engines.
"""
from pathlib import Path
import json
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


def seeded_dm(data,observation,seed_phase,seed_fom,solvent,reference,*,cycles,random_seed):
    """Frozen original HL + separate initial P+ map, never posterior-as-data."""
    from oasis_py.density_mlhl import fit_shells
    from oasis_py.native_density import (working_grid_size,metric_filter,reference_histogram,
        protein_mask,histogram_transform,apply_transform)
    from oasis_py.sad import posterior_moments
    from oasis_py.phasing import PhaseResult
    prior=np.asarray(observation.hl,dtype=float).copy()
    if prior.shape!=(len(data.f),4) or not np.isfinite(prior).all():raise ValueError('Invalid fixed observation HL')
    phase=np.asarray(seed_phase,dtype=float);fom=np.asarray(seed_fom,dtype=float)
    if phase.shape!=(len(data.f),) or fom.shape!=phase.shape or not np.isfinite([phase,fom]).all() or np.any((fom<0)|(fom>1)):
        raise ValueError('Invalid independent P+ seed')
    coefficients=data.f*fom*np.exp(1j*phase)
    def asu(values):return gemmi.ComplexAsuData(data.cell,data.spacegroup,np.ascontiguousarray(data.hkl),np.ascontiguousarray(values,dtype=np.complex64))
    shape,rate=working_grid_size(asu(coefficients),8_000_000);kernel=metric_filter(data.cell,shape,4.)
    target=reference_histogram(data,fom,reference);rng=np.random.default_rng(random_seed);history=[];hl=prior.copy()
    for cycle in range(cycles):
        grid=asu(coefficients).transform_f_phi_to_map(exact_size=shape);rho=np.array(grid.array,dtype=float)
        mask=protein_mask(grid,solvent,kernel);modified,mapping=histogram_transform(rho,mask,target)
        noise=rng.normal(size=shape);grid.array[:]=noise;grid.symmetrize_avg();noise=np.array(grid.array,dtype=float);noise-=noise.mean()
        noise*=max(rho.std(),1e-12)*.01/max(noise.std(),1e-12)
        perturbed=apply_transform(rho+noise,mask,mapping)
        gamma=float(np.clip(np.sum(noise*(perturbed-modified))/max(np.sum(noise*noise),1e-30),0,.95))
        grid.array[:]=modified-gamma*rho;grid.symmetrize_avg()
        fc=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(data.hkl).astype(complex)
        field,d,fit=fit_shells(data,fc,prior,samples=72)
        hl=prior.copy();hl[:,0]+=field.real;hl[:,1]+=field.imag
        phase,fom=posterior_moments(data,hl)
        coefficients=2*data.f*fom*np.exp(1j*phase)-d*fc
        row={'cycle':cycle+1,'mean_fom':float(fom.mean()),'gamma':gamma,'shell_fit':fit};history.append(row)
        print(json.dumps({'internal_DM_cycle':cycle+1,'mean_fom':row['mean_fom']}),flush=True)
    return PhaseResult(phase,fom,None,hl,{'method':'internal_fixed_observation_MLHL_seeded_DM','backend':'NumPy/SciPy/Gemmi only; no external executable',
        'cycles':history,'grid_shape':shape,'grid_sample_rate':rate,'seed_map':'P+ first moment once','fixed_HL':'Original internal observation likelihood, never accumulated',
        'limitations':['Approximate mask, reference histogram, shell scale and field cap inherited from internal MLHL DM','No automatic NCS; not RESOLVE/Parrot/Phenix equivalent']}),coefficients


def dm(cfg,output):
    from .sad_input import read_sad_pairs,read_phase_prior
    out=Path(output);parent=out.parent
    data,_,_=read_sad_pairs(cfg['prepared_sad'],free_label=FLAG,free_value=cfg['free_value'])
    original=read_phase_prior(parent/'02_sad/observations.mtz',data)
    seed=gemmi.read_mtz_file(str(parent/'02_sad/experimental.mtz'));lookup={tuple(h):i for i,h in enumerate(seed.make_miller_array())};ix=[lookup[tuple(h)] for h in data.hkl]
    phase=np.deg2rad(np.asarray(seed.column_with_label('PHIB'))[ix]);fom=np.asarray(seed.column_with_label('FOM'))[ix]
    seq=sample_info(cfg,data);solvent=seq['solvent_fraction'] if cfg['solvent_fraction'] is None else cfg['solvent_fraction']
    result,coeff=seeded_dm(data,original,phase,fom,solvent,cfg['reference'],cycles=cfg['dm_cycles'],random_seed=cfg['seed'])
    m=gemmi.Mtz(with_base=True);m.spacegroup=data.spacegroup;m.set_cell_for_all(data.cell);m.add_dataset('internal_DM')
    for name,kind in [('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHIB','P'),('FOM','W')]+[(x,'A') for x in ['HLA','HLB','HLC','HLD']]:m.add_column(name,kind)
    m.set_data(np.column_stack([data.hkl,data.f,data.sigf,abs(coeff),np.rad2deg(np.angle(coeff)),np.rad2deg(result.phase),result.fom,result.hl]).astype(np.float32));m.write_to_file(str(out/'map_phases.mtz'))
    dump(out/'map_phases.json',{'phase_method':'pplus_internal','dm_engine':'internal-MLHL-best','density_modification':result.metadata,
        'solvent_fraction_used':solvent,'solvent_source':'sequence_estimate' if cfg['solvent_fraction'] is None else 'user',
        'original_prior':fingerprint(parent/'02_sad/observations.mtz'),'Pplus_seed':fingerprint(parent/'02_sad/experimental.mtz'),
        'Phenix_executed':False,'external_programs':[]})
