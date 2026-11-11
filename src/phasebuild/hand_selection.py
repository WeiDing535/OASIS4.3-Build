"""Resolve SAD hands using matched work-only three-cycle sequence tracing.

No refinement, test-set residual, reference protein or SHELXE phase is used.
Ambiguous evidence is an explicit failure, never a silent first-file choice.
"""
import json
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET
import gemmi
import numpy as np
from .data import fingerprint,audit
from .sad_bridge import FLAG,keys
from .runtime import runtime_environment
from .trust_workflow import dump


def allowed_spacegroups(original):
    ops=original.operations();ops.change_basis_forward(original.change_of_hand_op())
    other=gemmi.find_spacegroup_by_ops(ops)
    if other is None:raise ValueError('Cannot identify the enantiomorphic space-group setting')
    return {original.hall,other.hall}


def hand_input(source,target,output):
    """Use Phaser's indexed hand setting; measured amplitudes/flags unchanged."""
    m=gemmi.read_mtz_file(str(source))
    if target.hall not in allowed_spacegroups(m.spacegroup):raise ValueError('Candidate is not the input space group or its allowed enantiomorph')
    before=np.array(m,copy=True);old=m.spacegroup.hm;m.spacegroup=target
    # Both groups have the same Laue rotations. No guessed sign or Bijvoet swap.
    asu=gemmi.ReciprocalAsu(target);ops=target.operations()
    if any(tuple(asu.to_asu(list(k),ops)[0])!=k for k in keys(m)):raise ValueError('Candidate requires an unsupported HKL reindexing; refusing silent relabelling')
    m.write_to_file(str(output));restored=gemmi.read_mtz_file(str(output))
    if not np.array_equal(before,np.array(restored),equal_nan=True):raise RuntimeError('Hand staging changed measured values or free flags')
    return dict(source=fingerprint(source),output=fingerprint(output),input_spacegroup=old,output_spacegroup=target.hm,measured_values_unchanged=True,HKLs_unchanged=True,free_flags_unchanged=True)


def pilot_metrics(model,xml):
    st=gemmi.read_structure(str(model));root=ET.parse(xml).getroot()
    cycles=root.findall('./Cycles/Cycle')
    if [int(c.findtext('Number')) for c in cycles]!=[1,2,3]:raise ValueError('Hand pilot did not finish exactly three cycles')
    final=root.find('Final')
    if final is None:raise ValueError('Hand pilot final XML missing')
    counts={k:int(final.findtext(tag,'-1')) for k,tag in [('unique','ResiduesUnique'),('sequenced','ResiduesSequenced'),('built','ResiduesBuilt'),('longest','ResiduesLongestFragment')]}
    if min(counts.values())<0:raise ValueError('Hand pilot XML counts invalid')
    counts['CA_residues']=sum('CA' in res for mod in st for chain in mod for res in chain)
    # Buccaneer reports unique assignment and sequencing using different
    # accounting; a few unique assignments can exceed the sequenced count.
    if counts['unique']>counts['built'] or counts['sequenced']>counts['built']:
        raise ValueError('Inconsistent pilot sequence counts')
    return counts


def choose(reports,expected):
    """Predeclared sequence-evidence gate; no R or reference-phase argument."""
    if not isinstance(expected,int) or expected<=0:raise ValueError('Positive ASU residue count required')
    if len(reports)!=2:raise ValueError('双手判别需要两个完整候选，不能默认选择第一个文件')
    if any(r.get('status') not in ('completed','empty') for r in reports):raise ValueError('至少一手计算未完整完成，无法可靠判别；请查看候选日志')
    ranked=sorted(reports,key=lambda r:(r.get('unique',0),r.get('sequenced',0),r.get('longest',0)),reverse=True)
    best,second=ranked
    threshold=max(3,min(15,int(np.ceil(.15*expected))))
    # The pilot is deliberately short and may sequence only a small fraction
    # of a large/multicopy ASU.  Cap the absolute separation requirement just
    # like the evidence threshold; otherwise four-copy targets can have clear
    # 27-versus-0 evidence rejected merely because the full ASU is large.
    margin=max(3,min(15,int(np.ceil(.05*expected))))
    if best.get('unique',0)==0 and second.get('unique',0)==0:
        # Very weak maps may build short poly-Ala fragments without assigning
        # sequence.  Accept only a large, predeclared work-only build-count
        # separation; close build counts remain explicitly ambiguous.
        best,second=sorted(reports,key=lambda r:(r.get('built',0),r.get('longest',0)),reverse=True)
        built_margin=max(20,int(np.ceil(.05*expected)))
        if (best.get('built',0)-second.get('built',0)>=built_margin and
                best.get('built',0)>=1.75*max(1,second.get('built',0))):
            return best['hand'],dict(policy='Maximum residues built after matched work-only three-cycle tracing when neither hand sequences; large absolute and ratio separation required',minimum_built_margin=built_margin,minimum_built_ratio=1.75,expected_ASU_residues=expected,selection_uses_Rfree=False,selection_uses_reference=False,selection_uses_refinement=False)
    if best.get('unique',0)<threshold or best.get('unique',0)>1.35*expected:raise ValueError('两手都缺少足够的唯一测序残基；手性未能判定，请核验 SHELXD 子结构与数据是否匹配')
    if best['unique']-second.get('unique',0)<margin:raise ValueError('两手测序证据过于接近；保留两手输出，拒绝静默选择')
    return best['hand'],dict(policy='Maximum uniquely sequenced residues after matched work-only three-cycle tracing; minimum sequence evidence and separation required',minimum_unique=threshold,minimum_margin=margin,expected_ASU_residues=expected,selection_uses_Rfree=False,selection_uses_reference=False,selection_uses_refinement=False)


def pilot(cfg,root):
    from .hybrid_recycle import construction_command,_construction_phase_columns
    from .work_weight import mask_test_observations
    root=Path(root);out=root/'pilot';out.mkdir()
    source=root/'04_bridge/building.mtz';masked=out/'work_only.mtz'
    mask=mask_test_observations(source,masked,cfg['free_value'])
    # The sequence is expanded exactly as in the final construction stage.
    from .oasis_options import read_sequences
    lines=[]
    for copy in range(cfg['copies']):
        for i,(_,seq) in enumerate(read_sequences(cfg['sequence'])):lines.extend(['>copy%d_chain%d'%(copy+1,i+1),seq])
    sequence=out/'building.fasta';sequence.write_text('\n'.join(lines)+'\n')
    runtime=Path(cfg['runtime']);library=runtime/'lib/data/reference_structures/reference-1tqw'
    command=construction_command(runtime/'bin/cbuccaneer',library,masked,sequence,out,1,None,_construction_phase_columns(audit(masked)),cfg['build_semet'],cfg.get('building_recipe','ccp4'),cfg['free_value'])
    dump(out/'command.json',command);env=runtime_environment(runtime,scratch=out/'scratch')
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    with (out/'buccaneer.log').open('x') as log:
        result=subprocess.run(command,cwd=out,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    if result.returncode:raise RuntimeError('Hand pilot Buccaneer failed; see '+str(out/'buccaneer.log'))
    counts=pilot_metrics(out/'model.pdb',out/'buccaneer.xml')
    return dict(status='completed' if counts['CA_residues'] else 'empty',**counts,mask_audit=mask,model=fingerprint(out/'model.pdb'),xml=fingerprint(out/'buccaneer.xml'))


def resolve(cfg,out,front,work):
    from . import bundled_pplus as b,trust_workflow as trust
    from .sad_input import read_sad_pairs
    out=Path(out);front=Path(front);original=gemmi.read_mtz_file(cfg['prepared_sad'])
    candidates=[('input_hand',front/'phaser.1.mtz'),('opposite_hand',front/'phaser.1.hand.mtz')]
    if any(not path.is_file() or not path.with_suffix('.pdb').is_file() for _,path in candidates):raise ValueError('Phaser 未输出完整双手 MTZ/PDB，无法自动判别')
    directory=out/'hands';directory.mkdir();reports=[];roots={};expected=None
    for name,raw in candidates:
        root=directory/name;root.mkdir();roots[name]=root;phase=root/'02_sad';phase.mkdir();dm=root/'03_dm';dm.mkdir();bridge=root/'04_bridge';bridge.mkdir()
        target=gemmi.read_mtz_file(str(raw)).spacegroup
        staged=root/'prepared_sad.mtz';native=root/'native.mtz'
        sad_audit=hand_input(cfg['prepared_sad'],target,staged);native_audit=hand_input(cfg['native'],target,native)
        local={**cfg,'prepared_sad':str(staged),'native':str(native),'auto_hand':False,'labels':''}
        record=dict(hand=name,spacegroup=target.hm,raw=fingerprint(raw),sad_staging=sad_audit,native_staging=native_audit)
        try:
            data,_,_=read_sad_pairs(staged,free_label=FLAG,free_value=cfg['free_value']);sample=trust.sample_info(local,data);expected=sample['asu_residues']
            model=phase/'refined_heavy.pdb';shutil.copy2(raw.with_suffix('.pdb'),model)
            prior=phase/'native_observations.mtz';b.extract_prior(raw,data,prior)
            local.update(heavy=str(model),observation_prior=str(prior),density_runtime=str(Path(cfg['runtime'])/'density'))
            trust.phase(local,phase);trust.dm(local,dm);trust.bridge(local,bridge)
            record.update(pilot(local,root))
            record['phase_artifacts']=[fingerprint(phase/name) for name in ('experimental.mtz','experimental.json','observations.mtz','observations.json','native_observations.mtz','refined_heavy.pdb')]
            record['dm_artifacts']=[fingerprint(p) for p in sorted(dm.iterdir()) if p.is_file()]
        except Exception as exc:
            record.update(status='failed',error=str(exc))
        reports.append(record);dump(root/'candidate.json',record)
    try:selected,policy=choose(reports,int(expected or 0))
    except ValueError as exc:
        dump(out/'hand_selection.json',dict(status='ambiguous',candidates=reports,error=str(exc),selection_uses_Rfree=False));raise
    chosen=next(r for r in reports if r['hand']==selected);root=roots[selected]
    for row in chosen['phase_artifacts']:
        if fingerprint(row['path'])!=row:raise RuntimeError('Chosen hand phase artifact changed')
        shutil.copy2(row['path'],out/Path(row['path']).name)
    shutil.copy2(root/'native.mtz',out/'selected_native.mtz');shutil.copy2(root/'prepared_sad.mtz',out/'selected_sad.mtz')
    dm_artifacts=chosen['dm_artifacts']
    report=dict(status='selected',selected_hand=selected,selected_spacegroup=chosen['spacegroup'],candidates=reports,policy=policy,
        selected_native=fingerprint(out/'selected_native.mtz'),selected_sad=fingerprint(out/'selected_sad.mtz'),selected_dm_artifacts=dm_artifacts,
        work=work,pilots_reused_for_final_model=False,original_input_files_modified=False)
    dump(out/'hand_selection.json',report)
    meta=json.loads((out/'experimental.json').read_text());meta.update(phase_method='pplus_bundled',auto_hand=report['selected_hand'],refined_substructure=fingerprint(out/'refined_heavy.pdb'));dump(out/'experimental.json',meta)
    return report


def copy_selected_dm(cfg,out):
    out=Path(out);phase=out.parent/'02_sad';report=json.loads((phase/'hand_selection.json').read_text())
    if report.get('status')!='selected':raise ValueError('Hand remains ambiguous')
    for key in ('selected_native','selected_sad'):
        if fingerprint(report[key]['path'])!=report[key]:raise RuntimeError('Selected hand input changed')
    for row in report['selected_dm_artifacts']:
        if fingerprint(row['path'])!=row:raise RuntimeError('Selected hand density artifact changed')
        shutil.copy2(row['path'],out/Path(row['path']).name)
    meta=json.loads((out/'map_phases.json').read_text());meta.update(auto_hand=report['selected_hand'],selected_spacegroup=report['selected_spacegroup'],hand_selection=fingerprint(phase/'hand_selection.json'),dm_computed_during_hand_selection=True,dm_engine='Bundled RESOLVE/Wang%g + CCTBX classic40'%cfg.get('trust_wang_radius',14.));dump(out/'map_phases.json',meta)
