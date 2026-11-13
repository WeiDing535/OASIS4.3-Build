"""Audited GUI adapter for the demonstrated original-likelihood P+ protocol.

The frontend is an explicitly supplied original Phaser SAD probability MTZ,
not an AutoSol/DM/model posterior. No target protein coordinates are accepted.
"""
import json
import os
import signal
from pathlib import Path
import shutil
import subprocess
import gemmi
import numpy as np
from .data import fingerprint
from .sad_bridge import FLAG, columns, keys, crystal


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def sample_info(cfg, data):
    from oasis_py.sequence import sequence_asu
    from .oasis_options import read_sequences
    path=cfg.get('composition_sequence') or cfg['sequence']
    info=sequence_asu(path,cfg['copies'],data.cell,data.spacegroup)
    # Composition may explicitly encode SeMet; construction uses ordinary M.
    import re
    records=gemmi.read_pir_or_fasta(Path(path).read_text())
    normalized=[''.join('M' if t=='(MSE)' else t for t in re.findall(r'\([A-Z0-9]+\)|[A-Z]',r.seq.upper())) for r in records]
    if normalized!=[s for _,s in read_sequences(cfg['sequence'])]:
        raise ValueError('元素组成序列必须与建模序列逐链一致，仅允许 M 与 (MSE) 的显式替换')
    if bool(info['asu_composition'].get('Se',0)) != bool(cfg['build_semet']):
        raise ValueError('SeMet 建模选项须与显式 (MSE) 元素组成序列一致')
    return info


def read_inputs(cfg, sad_path):
    from .sad_input import read_sad_pairs, read_phase_prior
    from .oasis_options import pairs
    from oasis_py.data import read_atoms
    data,single,meta=read_sad_pairs(sad_path,free_label=FLAG,free_value=cfg['free_value'])
    prior=read_phase_prior(cfg['observation_prior'],data)
    # Require the exact demonstrated probability schema; map columns indicate a
    # downstream product and are never silently accepted as original likelihood.
    original=gemmi.read_mtz_file(cfg['observation_prior'])
    columns(original,{'FP':'F','SIGFP':'Q','PHIB':'P','FOM':'W',**{x:'A' for x in ('HLA','HLB','HLC','HLD')}})
    if set(original.column_labels()) & {'FWT','PHWT','PHCOMB','FC','PHIC'}:
        raise ValueError('请提供原始纯观测 Phaser HL，不能使用带地图或模型列的下游结果')
    if not np.isfinite(np.column_stack([np.asarray(original.column_with_label('FP')),np.asarray(original.column_with_label('SIGFP'))])).all() or np.any(np.asarray(original.column_with_label('SIGFP'))<=0):
        raise ValueError('Invalid original Phaser amplitudes/sigmas')
    seq=sample_info(cfg,data)
    if cfg.get('method')=='pplus_bundled':
        from .heavy_scattering import read_heavy
        heavy=read_heavy(cfg['heavy'],data)
    else:heavy=read_atoms(cfg['heavy'],data,project_anisotropic=True)
    fpp={gemmi.Element(k).name:float(v) for k,v in pairs(cfg['fpp']).items()}
    return data,prior,seq,heavy,fpp,meta


def preflight(cfg):
    # Work flags are joined in memory: checking inputs must not create a file.
    from .sad_bridge import inspect_inputs
    from .sad_input import read_sad_pairs,read_phase_prior
    sad,_,flags,_=inspect_inputs(cfg['sad'],cfg['native'],cfg['free_value'])
    if cfg['labels']:raise ValueError('P+ 优化目前要求标准 F(+)/SIGF(+)/F(-)/SIGF(-) 列名')
    columns(sad,{'F(+)':'G','SIGF(+)':'L','F(-)':'G','SIGF(-)':'L'})
    # The persisted prior defines precisely the retained work/singleton set.
    # Use the same reader without staging files by validating its retained set
    # against the original flags after reading the unflagged SAD input.
    data,_,_=read_sad_pairs(cfg['sad'],free_label=None)
    by_hkl=dict(zip(keys(sad),flags));keep=np.array([np.isfinite(by_hkl[tuple(h)]) and by_hkl[tuple(h)]!=cfg['free_value'] for h in data.hkl])
    from oasis_py.data import Reflections
    work=Reflections(data.hkl[keep],data.f[keep],data.sigf[keep],data.cell,data.spacegroup,
                     dano=data.dano[keep],sigdano=data.sigdano[keep],source=data.source)
    read_phase_prior(cfg['observation_prior'],work)
    original=gemmi.read_mtz_file(cfg['observation_prior'])
    columns(original,{'FP':'F','SIGFP':'Q','PHIB':'P','FOM':'W',**{x:'A' for x in ('HLA','HLB','HLC','HLD')}})
    if set(original.column_labels()) & {'FWT','PHWT','PHCOMB','FC','PHIC'}:
        raise ValueError('原始观测先验不能含下游地图或模型列')
    original_keys=keys(original);asu=gemmi.ReciprocalAsu(original.spacegroup);ops=original.spacegroup.operations()
    if any(tuple(asu.to_asu(list(k),ops)[0])!=k for k in original_keys):
        raise ValueError('Original Phaser prior requires reciprocal-ASU HKLs matching native; reindex explicitly')
    values=np.column_stack([np.asarray(original.column_with_label('FP')),np.asarray(original.column_with_label('SIGFP'))])
    if not np.isfinite(values).all() or np.any(values[:,1]<=0):raise ValueError('Invalid original Phaser amplitudes/sigmas')
    sample_info(cfg,work)


def phase(cfg, output):
    from .pplus_experimental import prepare_pplus_normalization
    from .pplus_trust import phase_trusted_pplus
    from oasis_py.data import structure_factors
    from oasis_py.relations import build_relations
    out=Path(output);data,prior,seq,heavy,fpp,meta=read_inputs(cfg,cfg['prepared_sad'])
    context=prepare_pplus_normalization(data,seq['asu_composition'])
    graph=build_relations(data,context.e,context.kappa,cfg['kmin'],cfg['max_neighbors'],cfg['max_per_target'])
    from .heavy_scattering import anomalous_factors,scattering_metadata
    reference=np.angle(anomalous_factors(data,heavy,fpp))
    algorithm=phase_trusted_pplus
    if cfg.get('trust_stability'):
        from .pplus_stability import phase_stability_pplus
        algorithm=phase_stability_pplus
    result=algorithm(data,prior.hl,reference,graph,samples=1440,
                              kl_budget=cfg['trust_kl'],cycles=100,damping=.25)
    if result.metadata['mean_abs_probability_change']<=1e-12 or not np.any(abs(result.first-prior.fom*np.exp(1j*prior.phase))>1e-7):
        raise RuntimeError('P+ did not produce a nonzero phase-moment update; refusing to label the control as P+')
    seed=gemmi.Mtz(with_base=True);seed.spacegroup=data.spacegroup;seed.set_cell_for_all(data.cell);seed.add_dataset('Pplus_map_seed_no_HL')
    for label,kind in [('PHIB','P'),('FOM','W'),('PPLUS','R')]:seed.add_column(label,kind)
    seed.set_data(np.column_stack([data.hkl,np.rad2deg(result.phase),result.fom,result.probability]).astype(np.float32));seed.write_to_file(str(out/'experimental.mtz'))
    shutil.copyfile(cfg['observation_prior'],out/'observations.mtz')
    if fingerprint(out/'observations.mtz')['sha256']!=fingerprint(cfg['observation_prior'])['sha256']:raise RuntimeError('Original probability copy differs')
    dump(out/'experimental.json',{**result.metadata,'history':result.history,'input':meta,'sequence_asu':seq,'phase_method':'pplus_trust',
        'frontend':'Supplied original Phaser SAD likelihood; no new Phaser run', 'original_observations':fingerprint(cfg['observation_prior']),
        'normalization':context.statistics,'relations':graph.statistics,'conditional_samples':1440,'heavy_scattering':scattering_metadata(heavy)})
    dump(out/'observations.json',{'role':'Supplied original Phaser SAD observation likelihood; fixed for DM and building refinement',
                                'input':fingerprint(cfg['observation_prior']),'output':fingerprint(out/'observations.mtz'),'unchanged':True})


def launch(command, out, name, stdin=None, environment=None, timeout=1800):
    dump(out/(name+'_command.json'),command)
    with (out/(name+'.log')).open('x') as log:
        process=subprocess.Popen(command,stdin=subprocess.PIPE if stdin is not None else subprocess.DEVNULL,
            text=True,cwd=out,stdout=log,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
        try:
            process.communicate(input=stdin,timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait(timeout=5)
            raise
    if process.returncode:raise RuntimeError(name+' failed; see '+str(out/(name+'.log')))


def bridge(cfg, output):
    from .refinement import reflection_partition
    out=Path(output);native,_=reflection_partition(cfg['native'],free_value=cfg['free_value'])
    dmfile=out.parent/'03_dm/map_phases.mtz';m=gemmi.read_mtz_file(str(dmfile));m.ensure_asu();crystal(native,m)
    labels=columns(m,{'F':'F','SIGF':'Q','FWT':'F','PHWT':'P','PHIB':'P','FOM':'W'})
    a=np.array(native,copy=True);original=a.copy();nl=native.column_labels();b=np.array(m,copy=True);lookup={k:i for i,k in enumerate(keys(m))}
    fi,si,flag=[nl.index(k) for k in ['F','SIGF','FreeR_flag']];maps=np.full((len(a),4),np.nan,np.float32);matched=0
    for i,k in enumerate(keys(native)):
        j=lookup.get(k)
        if a[i,flag]==cfg['free_value'] or j is None:a[i,[fi,si]]=np.nan;continue
        f,s,fwt,phwt,phib,fom=b[j,[labels.index(x) for x in ['F','SIGF','FWT','PHWT','PHIB','FOM']]]
        if not np.isfinite([f,s,fwt,phwt,phib,fom,a[i,fi],a[i,si]]).all() or min(f,s,a[i,fi],a[i,si])<=0:
            a[i,[fi,si]]=np.nan;continue
        if fwt<0 or not 0<=fom<=1:raise ValueError('Invalid DM map/FOM')
        maps[i]=[fwt*(a[i,fi]/f),phwt,phib,fom];matched+=1
    if not matched:raise ValueError('No usable native/DM map matches')
    for label,kind in [('FWT','F'),('PHWT','P'),('PHCOMB','P'),('FOM','W')]:
        if label in nl:raise ValueError('Native contains downstream map columns')
        native.add_column(label,kind)
    native.set_data(np.column_stack([a,maps]));native.write_to_file(str(out/'building.mtz'))
    free=original[:,flag]==cfg['free_value']
    if np.isfinite(a[free][:,[fi,si]]).any() or np.isfinite(maps[free]).any():raise RuntimeError('Free observations leaked to initial construction')
    dump(out/'building.audit.json',{'native':fingerprint(cfg['native']),'dm':fingerprint(dmfile),'matched_work_rows':matched,
        'original_flags_exact':bool(np.array_equal(original[:,flag],a[:,flag],equal_nan=True)),
        'initial_free_F_SIGF_and_map_excluded':True,'refinement_native_unchanged':True,
        'map_formula':'actual DM FWT * native F / DM observed F, per HKL; no extra FOM multiplication',
        'probability_formula':'PHCOMB=DM PHIB; FOM=DM FOM'})


def accept_cleanup_by_work_fit(raw_metric, cleaned_metric, max_work_rise=.02):
    """Predeclared work-only guard; never inspect Rfree for model selection."""
    raw=float(raw_metric['r_work']);cleaned=float(cleaned_metric['r_work'])
    if not np.isfinite([raw,cleaned]).all():
        raise ValueError('Nonfinite work R in cleanup comparison')
    return cleaned <= raw + max_work_rise


def finish(cfg, output):
    from .model_cleanup import cleanup_single_copy
    from .oasis_options import read_sequences
    from .runtime import runtime_environment
    from .refinement import reflection_partition,verify_free_flags
    from .phase_refinement import parse_statistics
    out=Path(output);parent=out.parent;report=json.loads((parent/'05_build/hybrid_report.json').read_text());raw=parent/'05_build/primary_model.cif'
    model=raw;metric=report['primary_metric'];status='disabled';reason=None;post_weight=.005;calibration=None
    if cfg['trust_cleanup']:
        try:
            records=read_sequences(cfg['sequence'])
            if len(records)!=1:raise ValueError('Single-copy cleanup requires exactly one sequence record')
            try:
                cleaned,audit=cleanup_single_copy(gemmi.read_structure(str(raw)),records[0][1],asu_copies=cfg['copies'])
            except ValueError:
                if not cfg.get('trust_partial_cleanup'):raise
                from .partial_model_cleanup import work_density,cleanup_partial_single_copy
                grid,map_audit=work_density(parent/'03_dm/map_phases.mtz',cfg['native'],cfg['free_value'])
                cleaned,audit=cleanup_partial_single_copy(gemmi.read_structure(str(raw)),records[0][1],grid,asu_copies=cfg['copies'])
                audit['work_map_audit']=map_audit
        except ValueError as error:status='refused';reason=str(error)
        else:
            cleaned.make_mmcif_document().write_file(str(out/'cleaned.cif'));dump(out/'cleanup.json',audit)
            if cfg.get('trust_work_weight'):
                from .work_weight import calibrate_work_weight
                post_weight,calibration=calibrate_work_weight(cfg['native'],out/'cleaned.cif',cfg['runtime'],out/'work_weight_calibration',cfg['free_value'])
            shutil.copyfile(cfg['native'],out/'input.mtz');cleaned.write_pdb(str(out/'input.pdb'))
            inp,part=reflection_partition(cfg['native'],free_value=cfg['free_value'])
            keywords='LABIN FP=F SIGFP=SIGF FREE=FreeR_flag\nFREE %d\nREFI TYPE REST RESO %.6f %.6f\nREFI RESI MLKF\nREFI BREF ISOT\nMAKE HYDR N HOUT N BUIL N\nSCAL TYPE BULK\nSCAL MLSC WORK\nWEIGHT MATRIX %.6g\nNCYC 30\nEND\n'%(cfg['free_value'],part['observed_resolution_low']*1.0001,part['observed_resolution_high']*.9999,post_weight)
            (out/'keywords.txt').write_text(keywords)
            command=[str(Path(cfg['runtime'])/'bin/refmac5'),'HKLIN','input.mtz','XYZIN','input.pdb','HKLOUT','refined.mtz','XYZOUT','refined.pdb','XMLOUT','refmac.xml']
            env=runtime_environment(Path(cfg['runtime']),scratch=out/'scratch')
            # The post-cleanup refinement is an optional improvement. A numerical
            # failure here must not discard a completed phasing/building endpoint;
            # the fixed last build round is kept instead and the failure recorded.
            try:
                launch(command,out,'refmac',keywords,env)
                stats=parse_statistics(out/'refmac.log',out/'refmac.xml')
                if stats['final']['cycle']!=30:raise RuntimeError('Post-cleanup refinement did not finish')
                flags=verify_free_flags(inp,gemmi.read_mtz_file(str(out/'refined.mtz')),'FreeR_flag','F','SIGF')
                gemmi.read_structure(str(out/'refined.pdb')).make_mmcif_document().write_file(str(out/'model.cif'))
            except (RuntimeError,ValueError,OSError) as error:
                status='refused_refinement_failure'
                reason='Post-cleanup refinement failed (%s); retained fixed last build round'%str(error)[:200]
                dump(out/'refinement_report.json',{'status':'failed','error':str(error),
                    'requested_cycles':30,'native_input':fingerprint(cfg['native']),
                    'retained':'fixed last build round'})
            else:
                accepted=accept_cleanup_by_work_fit(report['primary_metric'],stats['final'])
                if accepted:
                    model=out/'model.cif';metric=stats['final'];status='completed'
                else:
                    status='rejected_work_fit'
                    reason='Post-cleanup work R rose by more than 0.02 absolute; retained fixed last build round'
                dump(out/'refinement_report.json',{'status':'success','statistics':stats,'requested_cycles':30,'free_flags_preserved':True,'free_flag_audit':flags,'native_input':fingerprint(cfg['native'])})
    shutil.copyfile(model,out/'primary_model.cif')
    dump(out/'endpoint.json',{'status':'completed','cleanup_status':status,'cleanup_reason':reason,
        'raw_model':fingerprint(raw),'raw_metric':report['primary_metric'],'primary_model':fingerprint(out/'primary_model.cif'),
        'primary_metric':metric,'post_refinement_weight':post_weight if status=='completed' else None,'work_weight_calibration':calibration,'cleanup_kind':audit['operation'] if status=='completed' else None,'selection':'fixed last build round; optional declared full/partial sequence cleanup then fixed30 amplitude-only cycles; accept cleanup only when work R rises by at most 0.02 absolute; never select by Rfree',
        'refined_mtz':fingerprint(out/'refined.mtz') if status=='completed' else None})
