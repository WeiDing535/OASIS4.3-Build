"""Self-contained engine adapter: native Phaser -> genuine P+ -> CCP4 Parrot.

Independence means no separately installed suites. Phaser and Parrot are part of
the CCP4 that also supplies cbuccaneer/refmac5, so one CCP4 install runs the whole
chain; the in-tree P+ likelihood refinement is OASIS's own code.
"""
from pathlib import Path
import json
import os
import sys
import shutil
import gemmi
import numpy as np
from .sad_bridge import FLAG,columns,keys,crystal
from .data import fingerprint
from .trust_workflow import dump,launch,sample_info


def engine_layout(runtime,dm_engine=None):
    """Check the bundled engine; density modification is CCP4 Parrot's job."""
    runtime=Path(runtime).expanduser().resolve()
    if dm_engine not in (None,'parrot','internal'):raise ValueError('未知密度修饰引擎：'+str(dm_engine))
    phaser=runtime/'bin/phaser'
    if not phaser.is_file() or not os.access(phaser,os.X_OK):
        raise ValueError('P+ 引擎缺失：'+str(phaser)+'。若“CCP4 位置”填的是外部安装，请确认它有 bin/phaser，或清空该栏改用包内运行时')
    return runtime


def preflight(cfg):
    from .oasis_options import read_sequences
    if cfg.get('observation_prior'):raise ValueError('自带引擎模式直接计算 SAD 似然，不接受外部相位文件')
    value=cfg.get('wavelength')
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or not .1<=value<=10:
        raise ValueError('请填写实验 SAD 波长（Å，0.1–10）；程序不猜测波长')
    engine_layout(cfg['runtime'],cfg.get('dm_engine'))
    if cfg.get('labels'):raise ValueError('自带引擎模式当前要求标准 F(+)/SIGF(+)/F(-)/SIGF(-) 列名')
    sad=gemmi.read_mtz_file(cfg['sad']);columns(sad,{'F(+)':'G','SIGF(+)':'L','F(-)':'G','SIGF(-)':'L'})
    read_sequences(cfg['sequence']);sample_info(cfg,sad)


def stage_work(cfg,output):
    """Exact HKL join, measured work data only; no positional concatenation."""
    sad=gemmi.read_mtz_file(cfg['prepared_sad']);native=gemmi.read_mtz_file(cfg['native']);crystal(sad,native)
    sk,nk=keys(sad),keys(native);lookup={h:i for i,h in enumerate(nk)}
    sl=columns(sad,{FLAG:'I','F(+)':'G','SIGF(+)':'L','F(-)':'G','SIGF(-)':'L'})
    nl=columns(native,{'F':'F','SIGF':'Q','FreeR_flag':'I'})
    sa=np.array(sad,copy=True);na=np.array(native,copy=True);rows=[]
    for i,h in enumerate(sk):
        flag=sa[i,sl.index(FLAG)]
        if not np.isfinite(flag) or flag==cfg['free_value']:continue
        j=lookup.get(h)
        if j is None:raise ValueError('Working SAD HKL is absent from native')
        if flag!=na[j,nl.index('FreeR_flag')]:raise ValueError('Work/free flags disagree')
        rows.append([*h,*sa[i,[sl.index(c) for c in ['F(+)','SIGF(+)','F(-)','SIGF(-)']]],
                     na[j,nl.index('F')],na[j,nl.index('SIGF')],flag])
    if not rows:raise ValueError('No shared measured working reflections')
    m=gemmi.Mtz(with_base=True);m.spacegroup=sad.spacegroup;m.set_cell_for_all(sad.cell);dataset=m.add_dataset('measured_SAD_work');dataset.wavelength=cfg['wavelength']
    for c,t in [('F(+)','G'),('SIGF(+)','L'),('F(-)','G'),('SIGF(-)','L'),('F','F'),('SIGF','Q'),('FreeR_flag','I')]:m.add_column(c,t)
    m.set_data(np.asarray(rows,np.float32));m.write_to_file(str(output))
    assert not np.any(m.column_with_label('FreeR_flag').array==cfg['free_value'])
    return {'work_rows':m.nreflections,'measured_free_rows':0,'matching':'Unique integer HKL, never row order; preserve incoming SAD row order to avoid unnecessary floating-point accumulation changes','output':fingerprint(output)}


def extract_prior(raw_path,data,output):
    """Integrate original HL; unavailable likelihoods remain uninformative."""
    from oasis_py.data import write_mtz
    from oasis_py.sad import posterior_moments
    raw=gemmi.read_mtz_file(str(raw_path));crystal(raw,data)
    raw.ensure_asu();rk=keys(raw);labels=columns(raw,{c:'A' for c in ['HLA','HLB','HLC','HLD']});lookup={h:i for i,h in enumerate(rk)}
    a=np.array(raw,copy=False);hl=np.zeros((len(data.f),4));missing=0
    for i,h in enumerate(data.hkl):
        j=lookup.get(tuple(h))
        if j is None:missing+=1;continue
        v=a[j,[labels.index(c) for c in ['HLA','HLB','HLC','HLD']]]
        if not np.isfinite(v).all():missing+=1;continue
        hl[i]=v
    phase,fom=posterior_moments(data,hl)
    write_mtz(output,data,phase,fom,hl)
    return {'source':fingerprint(raw_path),'output':fingerprint(output),'uninformative_missing_HL_rows':missing,
            'missing_HL_policy':'Zero HL = uniform/allowed-centric likelihood; measured amplitudes remain measured, no protein phase substitution',
            'PHIB_FOM':'Reintegrated from persisted original HL; not copied inconsistently from another distribution'}


def phase(cfg,output):
    from .sad_input import read_sad_pairs
    from .oasis_options import pairs
    from .runtime import runtime_environment
    from . import trust_workflow
    out=Path(output)
    # Everything here comes from the runtime: Phaser for the frontend, the
    # in-tree P+ likelihood refinement, and CCP4 Parrot for density
    # modification.  No separately copied Phenix tree is involved any more.
    runtime=engine_layout(cfg['runtime'],cfg.get('dm_engine'));front=out/'phaser';front.mkdir()
    work=stage_work(cfg,front/'work.mtz');shutil.copyfile(cfg['sequence'],front/'sequence.fasta')
    if Path(cfg['heavy']).suffix.lower() in ('.pdb','.ent'):shutil.copyfile(cfg['heavy'],front/'heavy.pdb')
    else:gemmi.read_structure(cfg['heavy']).write_pdb(str(front/'heavy.pdb'))
    data,_,_=read_sad_pairs(cfg['prepared_sad'],free_label=FLAG,free_value=cfg['free_value']);seq=sample_info(cfg,data)
    solvent=seq['solvent_fraction'] if cfg['solvent_fraction'] is None else cfg['solvent_fraction']
    recipe=['TITLE OASIS bundled work-only SAD frontend','MODE EP_AUTO','HKLIN work.mtz',
        'ATOM CRYSTAL sample PDB heavy.pdb','CRYSTAL sample DATASET peak LABIN F+ = F(+) SIGF(+) = SIGF(+) F- = F(-) SIGF(-) = SIGF(-)',
        'WAVELENGTH %.10g'%cfg['wavelength'],'COMPOSITION PROTEIN SEQ sequence.fasta NUMBER %d'%cfg['copies'],
        'HAND BOTH','ROOT phaser','LLGCOMPLETE COMPLETE ON']
    for element in pairs(cfg['fpp']):recipe.append('LLGCOMPLETE SCATTERER '+gemmi.Element(element).name)
    recipe+=['AVERAGE SOLVENT %.16g'%solvent,'JOBS 1','END'];text='\n'.join(recipe)+'\n';(front/'keywords.txt').write_text(text)
    env=runtime_environment(runtime,scratch=front/'scratch')
    launch([str(runtime/'bin/phaser')],front,'phaser',text,env,timeout=3600)
    if cfg.get('auto_hand'):
        from .hand_selection import resolve
        selected=resolve(cfg,out,front,work)
        # The hand decision is the run's first real result: say it out loud so
        # the stage log (and the GUI log pane) shows more than a stage header.
        for row in (json.loads((out/'hand_selection.json').read_text()).get('candidates') or []):
            print('hand %s: %s, %s uniquely sequenced residues'
                  %(row.get('hand'),row.get('status'),row.get('unique')),flush=True)
        print('hand selection: %s (%s)'%(selected.get('selected_hand'),
              'verified' if selected.get('hand_verified',True) else 'unverified'),flush=True)
        scattering={'mode':'legacy_explicit','fpp_argument':cfg['fpp']}
        if cfg.get('auto_fpp'):
            from .anomalous_scattering import calculate
            scattering=calculate(cfg['elements'],cfg['wavelength'])
        meta=json.loads((out/'experimental.json').read_text());meta.update(frontend='Bundled native Phaser EP_AUTO; both hands evaluated by work-only sequence tracing',external_phase_file_used=False);dump(out/'experimental.json',meta)
        policy='Both Phaser hands: P+ -> Wang/classic40 -> masked work-only three-cycle sequence tracing; no refinement, Rfree or reference selection'
        if selected.get('unverified_selection'):
            policy+='. The sequence evidence was insufficient, so the configured input-hand continuation was used; handedness is unverified for this run, see hand_selection.json'
            print('警告：手性证据不足，按“证据不足也继续”的设置沿用了选中的手；本次结果手性未验证（详见 hand_selection.json）',file=sys.stderr)
        elif selected.get('single_hand_selection'):
            policy+='. Only one hand could be computed, so handedness is unverified for this run; see hand_selection.json for the failed hand'
            print('警告：只有一手算出结果，本次手性未交叉验证（详见 hand_selection.json）',file=sys.stderr)
        dump(out/'frontend.json',{'engine':'Bundled native Phaser','binary':fingerprint(runtime/'bin/phaser'),'work':work,
            'prior':{'output':fingerprint(out/'native_observations.mtz')},'refined_substructure':fingerprint(out/'refined_heavy.pdb'),
            'pplus_scattering':scattering,'wavelength':cfg['wavelength'],'hand_policy':policy,'hand_verified':selected.get('hand_verified',True),
            'unverified_selection':bool(selected.get('unverified_selection',False)),
            'hand_selection':fingerprint(out/'hand_selection.json'),'selected_hand':selected['selected_hand'],'selected_spacegroup':selected['selected_spacegroup'],
            'inputs':{k:fingerprint(cfg[k]) for k in ['sad','native','heavy','sequence']},'third_party_engines':True,'requires_separate_suite_installation':False})
        return
    candidates=[front/'phaser.1.mtz',front/'phaser.1.hand.mtz']
    selected=next((p for p in candidates if p.is_file() and gemmi.read_mtz_file(str(p)).spacegroup.hall==data.spacegroup.hall),None)
    if selected is None:raise RuntimeError('Phaser produced no solution matching the supplied data space-group setting')
    model=selected.with_suffix('.pdb')
    if not model.is_file():raise RuntimeError('Phaser did not save its refined substructure')
    shutil.copyfile(model,out/'refined_heavy.pdb')
    prior=extract_prior(selected,data,out/'native_observations.mtz')
    local={**cfg,'observation_prior':str(out/'native_observations.mtz'),'heavy':str(out/'refined_heavy.pdb')}
    trust_workflow.phase(local,out)
    meta=json.loads((out/'experimental.json').read_text());meta.update(phase_method='pplus_bundled',frontend='Self-contained native Phaser EP_AUTO; original likelihood computed here from work data',external_phase_file_used=False,refined_substructure=fingerprint(out/'refined_heavy.pdb'))
    dump(out/'experimental.json',meta)
    scattering={'mode':'legacy_explicit','fpp_argument':cfg['fpp']}
    if cfg.get('auto_fpp'):
        from .anomalous_scattering import calculate
        scattering=calculate(cfg['elements'],cfg['wavelength'])
    dump(out/'frontend.json',{'engine':'Bundled native Phaser','binary':fingerprint(runtime/'bin/phaser'),'work':work,'prior':prior,'refined_substructure':fingerprint(out/'refined_heavy.pdb'),
        'pplus_scattering':scattering,
        'wavelength':cfg['wavelength'],'hand_policy':'Retain both outputs; first candidate matching supplied space-group setting, no Rfree or target-reference selection',
        'inputs':{k:fingerprint(cfg[k]) for k in ['sad','native','heavy','sequence']},'third_party_engines':True,'requires_separate_suite_installation':False})


def dm(cfg,output):
    if cfg.get('auto_hand'):
        from .hand_selection import copy_selected_dm
        return copy_selected_dm(cfg,output)
    raise ValueError('P+ 自带引擎的密度修饰统一由 CCP4 Parrot 完成（--dm-engine parrot）；'
                     '旧的包内 RESOLVE/CCTBX 后端已移除')
