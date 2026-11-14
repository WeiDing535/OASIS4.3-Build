"""Work-only hand decisions, true enantiomorphic settings and frozen caches."""
from pathlib import Path
import json,shutil
from types import SimpleNamespace
import gemmi,numpy as np,pytest
from phasebuild import hand_selection as hs,hybrid_recycle as h,bundled_pplus as b,trust_workflow as trust,sad_pipeline
from phasebuild.data import fingerprint
from phasebuild.sad_bridge import prepare_sad,FLAG
from phasebuild.sad_input import read_sad_pairs
from test_sad_pipeline import mtz,inputs
from test_trust_workflow import trusted
from test_bundled_pplus import bundled


def reports(a,b):
 return [dict(hand='input_hand',status='completed' if a else 'empty',unique=a,sequenced=a,built=max(a,4),longest=a),dict(hand='opposite_hand',status='completed' if b else 'empty',unique=b,sequenced=b,built=max(b,4),longest=b)]


def test_pilot_metrics_accepts_buccaneer_unique_count_above_sequenced(tmp_path):
 model=tmp_path/'model.pdb';model.write_text('END\n')
 xml=tmp_path/'pilot.xml'
 xml.write_text('<BuccaneerResult><Cycles><Cycle><Number>1</Number></Cycle><Cycle><Number>2</Number></Cycle><Cycle><Number>3</Number></Cycle></Cycles><Final><ResiduesUnique>390</ResiduesUnique><ResiduesSequenced>387</ResiduesSequenced><ResiduesBuilt>406</ResiduesBuilt><ResiduesLongestFragment>95</ResiduesLongestFragment></Final></BuccaneerResult>')
 assert hs.pilot_metrics(model,xml)['unique']==390


@pytest.mark.parametrize('a,b,winner',[(40,6,'input_hand'),(6,40,'opposite_hand'),(0,40,'opposite_hand')])
def test_choice_uses_sequence_and_not_hand_order_or_test_score(a,b,winner):
 rows=reports(a,b);rows[0]['r_free']=.1;rows[1]['r_free']=.9
 selected,policy=hs.choose(rows,51);assert selected==winner
 assert not policy['selection_uses_Rfree'] and not policy['selection_uses_reference'] and not policy['selection_uses_refinement']
 assert hs.choose(rows[::-1],51)[0]==winner


@pytest.mark.parametrize('rows',[reports(0,0),reports(4,2),reports(40,39),reports(80,5),[dict(hand='one',status='completed',unique=40)]])
def test_ambiguous_or_unphysical_evidence_fails_closed(rows):
 with pytest.raises(ValueError):hs.choose(rows,51)


def test_insufficient_evidence_refuses_without_the_explicit_opt_in():
    for rows in (reports(0,0),reports(4,2),reports(40,39)):
        with pytest.raises(ValueError):
            hs.choose(rows,51,allow_unverified=False)


def test_opt_in_carries_the_input_hand_and_records_the_refusal():
    rows=reports(4,2)
    selected,policy=hs.choose(rows,51,allow_unverified=True)
    assert selected=='input_hand'
    assert policy['forced_unverified'] is True and policy['hand_verified'] is False
    assert policy['single_hand'] is False
    assert '缺少足够的唯一测序残基' in policy['refusal_reason']
    assert not policy['selection_uses_Rfree'] and not policy['selection_uses_reference'] and not policy['selection_uses_refinement']


def test_opt_in_never_rescues_a_run_where_neither_hand_computed():
    rows=reports(0,0)
    for row in rows:row.update(status='failed',error='missing engine')
    with pytest.raises(ValueError,match='两手都没算出结果'):
        hs.choose(rows,51,allow_unverified=True)


def test_opt_in_does_not_change_a_decision_the_evidence_supports():
    selected,policy=hs.choose(reports(40,6),51,allow_unverified=True)
    assert selected=='input_hand' and policy['hand_verified'] is True
    assert not policy.get('forced_unverified')


def test_opt_in_carries_the_surviving_hand_when_the_other_failed_short_of_evidence():
    rows=reports(4,0);rows[1].update(status='failed',error='missing engine')
    selected,policy=hs.choose(rows,51,allow_unverified=True)
    assert selected=='input_hand'
    assert policy['forced_unverified'] is True and policy['single_hand'] is True
    assert policy['hand_verified'] is False and 'missing engine' in policy['other_hand_error']


def test_failed_opposite_hand_continues_only_on_its_own_positive_evidence():
 rows=reports(40,0);rows[1].update(status='failed',error='missing engine')
 selected,policy=hs.choose(rows,51)
 assert selected=='input_hand' and policy['single_hand'] and not policy['hand_verified']
 assert 'missing engine' in policy['other_hand_error'] and not policy['selection_uses_Rfree']
 weak=reports(4,0);weak[1].update(status='failed',error='missing engine')
 with pytest.raises(ValueError,match='未完整'):hs.choose(weak,51)


def test_large_multicopy_asu_accepts_clear_short_pilot_separation():
 rows=reports(27,0)
 selected,policy=hs.choose(rows,748)
 assert selected=='input_hand' and policy['minimum_margin']==15
 assert not policy['selection_uses_Rfree']


def test_zero_sequence_can_use_large_work_only_build_separation():
 rows=reports(0,0);rows[0].update(built=57,longest=9);rows[1].update(built=18,longest=10)
 selected,policy=hs.choose(rows,200)
 assert selected=='input_hand' and policy['minimum_built_ratio']==1.75
 assert not policy['selection_uses_Rfree'] and not policy['selection_uses_reference']


def test_zero_sequence_close_build_counts_remain_ambiguous():
 rows=reports(0,0);rows[0].update(built=48,longest=12);rows[1].update(built=48,longest=14)
 with pytest.raises(ValueError,match='缺少足够'):hs.choose(rows,200)


def test_enantiomorphic_spacegroup_propagates_without_changing_observations(tmp_path):
 cell=gemmi.UnitCell(20,20,30,90,90,90);sg=gemmi.SpaceGroup('P 41');hk=gemmi.make_miller_array(cell,sg,4.)
 source=tmp_path/'source.mtz';m=mtz(source,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],np.column_stack([hk,np.ones(len(hk))*20,np.ones(len(hk)),np.arange(len(hk))%10]));m.spacegroup=sg;m.set_cell_for_all(cell);m.write_to_file(str(source));before=np.array(m)
 target=tmp_path/'hand.mtz';meta=hs.hand_input(source,gemmi.SpaceGroup('P 43'),target)
 out=gemmi.read_mtz_file(str(target));np.testing.assert_array_equal(np.array(out),before)
 assert out.spacegroup.hm=='P 43' and gemmi.read_mtz_file(str(source)).spacegroup.hm=='P 41'
 assert meta['free_flags_unchanged'] and meta['HKLs_unchanged']
 with pytest.raises(ValueError,match='enantiomorph'):hs.hand_input(source,gemmi.SpaceGroup('P 31'),tmp_path/'bad.mtz')


def test_ccp4_recipe_channels_and_explicit_parameters(tmp_path):
 m=tmp_path/'data.mtz';mtz(m,[('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHCOMB','P'),('FOM','W')]+[(k,'A') for k in ('HLACOMB','HLBCOMB','HLCCOMB','HLDCOMB')],[[1,0,0,10,1,5,30,40,.5,1,0,0,0]])
 args=(Path('/bin/cbuccaneer'),Path('/reference/1tqw'),m,tmp_path/'seq',tmp_path)
 first=h.construction_command(*args,1,None,{'probability_phase':'PHCOMB'},recipe='ccp4')
 assert '-colin-fc' not in first and '-model-filter' not in first and '-fast' in first
 assert first[first.index('-sequence-reliability')+1]=='0.95'
 second=h.construction_command(*args,2,tmp_path/'old.cif',{'probability_phase':'PHCOMB'},recipe='ccp4')
 assert second[second.index('-colin-hl')+1]=='/*/*/[HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB]'
 assert '-colin-phifom' not in second and '-colin-fc' in second
 other_free=h.construction_command(*args,1,None,{'probability_phase':'PHCOMB'},recipe='ccp4',free_value=7)
 assert '-colin-free' not in other_free


def test_pipeline_recipe_maps_to_the_ccp4_construction(tmp_path):
 m=tmp_path/'data.mtz';mtz(m,[('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHCOMB','P'),('FOM','W')],[[1,0,0,10,1,5,30,40,.5]])
 args=(Path('/bin/cbuccaneer'),Path('/reference/1tqw'),m,tmp_path/'seq',tmp_path)
 # The hand pilot forwards the pipeline-level recipe; it must not raise.
 command=h.construction_command(*args,1,None,{'probability_phase':'PHCOMB'},recipe='pipeline')
 assert '-fast' in command and '-model-filter' not in command
 with pytest.raises(ValueError,match='Unknown construction recipe'):
  h.construction_command(*args,1,None,{'probability_phase':'PHCOMB'},recipe='invented')

def test_auto_hand_stage_freezes_selected_data_and_recipe(bundled):
 cfg=sad_pipeline.validate({**bundled,'auto_hand':True,'building_recipe':'ccp4'});plan=sad_pipeline.stage_plan(cfg)
 assert str(Path(cfg['output'])/'02_sad/hand_selection.json') in plan[1]['outputs']
 build=plan[4]['command'];assert build[build.index('--native')+1].endswith('02_sad/selected_native.mtz')
 assert build[build.index('--building-recipe')+1]=='ccp4'
 with pytest.raises(ValueError):sad_pipeline.validate({**bundled,'method':'pplus','auto_hand':True})


def test_resolver_promotes_opposite_hand_and_checks_cache(bundled,tmp_path,monkeypatch):
 cfg=sad_pipeline.validate(bundled);prepared=tmp_path/'prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg.update(prepared_sad=str(prepared),auto_hand=True,building_recipe='ccp4')
 data,_,_=read_sad_pairs(prepared,free_label=FLAG,free_value=cfg['free_value']);before=fingerprint(prepared)
 out=tmp_path/'02_sad';out.mkdir();front=out/'phaser';front.mkdir()
 for stem in ('phaser.1','phaser.1.hand'):
  mtz(front/(stem+'.mtz'),[(name,'A') for name in ('HLA','HLB','HLC','HLD')],np.column_stack([data.hkl,np.ones((len(data.hkl),4))]));shutil.copy2(cfg['heavy'],front/(stem+'.pdb'))
 def fake_phase(local,folder):
  for name in ('experimental.mtz','observations.mtz'):shutil.copy2(local['observation_prior'],folder/name)
  for name in ('experimental.json','observations.json'):(folder/name).write_text('{}')
 def fake_dm(local,phase,folder,native,solvent):
  (folder/'map_phases.mtz').write_bytes(b'mock density artifact');(folder/'map_phases.json').write_text('{}')
 def fake_bridge(local,folder):(folder/'building.mtz').write_bytes(b'mock map')
 def fake_pilot(local,root):
  count=12 if root.name=='opposite_hand' else 2
  return dict(status='completed',unique=count,sequenced=count,built=count,longest=count)
 monkeypatch.setattr(trust,'phase',fake_phase);monkeypatch.setattr(hs,'density_modify',fake_dm);monkeypatch.setattr(trust,'bridge',fake_bridge);monkeypatch.setattr(hs,'pilot',fake_pilot)
 report=hs.resolve(cfg,out,front,{'work_rows':len(data.hkl)})
 assert report['selected_hand']=='opposite_hand' and not report['pilots_reused_for_final_model']
 assert fingerprint(prepared)==before
 dm=tmp_path/'03_dm';dm.mkdir();hs.copy_selected_dm(cfg,dm)
 assert (dm/'map_phases.mtz').read_bytes()==b'mock density artifact'
 assert json.loads((dm/'map_phases.json').read_text())['auto_hand']=='opposite_hand'
 row=report['selected_dm_artifacts'][0];Path(row['path']).write_bytes(b'changed')
 other=tmp_path/'changed_dm';other.mkdir()
 with pytest.raises(RuntimeError,match='artifact changed'):hs.copy_selected_dm(cfg,other)


def test_dm_restraint_role_is_explicit_and_free_phase_values_are_masked(inputs,tmp_path):
 from phasebuild.phase_refinement import merge_experimental_phases
 native=gemmi.read_mtz_file(inputs['native']);hkl=native.make_miller_array();phase=tmp_path/'improved.mtz'
 mtz(phase,[('PHIB','P'),('FOM','W')],np.column_stack([hkl,np.full(len(hkl),30),np.full(len(hkl),.8)]))
 out=tmp_path/'merged.mtz';before=np.array(native,copy=True)
 meta=merge_experimental_phases(inputs['native'],phase,out,'phifom',free_value=7,phase_prefix='DM_')
 m=gemmi.read_mtz_file(str(out));flags=np.asarray(m.column_with_label('FreeR_flag'));values=np.array(m)
 np.testing.assert_array_equal(values[:,:before.shape[1]],before)
 assert meta['matched_free_prior_rows']==0 and meta['phase_source_kind']=='density_modified_work'
 assert 'DM_PHIB' in m.column_labels() and 'SAD_PHIB' not in m.column_labels()
 assert np.isnan(np.asarray(m.column_with_label('DM_PHIB'))[flags==7]).all()


def test_dm_prior_refiner_does_not_claim_original_sad_likelihood(inputs,tmp_path,monkeypatch):
 from phasebuild import fixed_phase_prior as fixed
 native=gemmi.read_mtz_file(inputs['native']);hkl=native.make_miller_array();phase=tmp_path/'DM.mtz'
 mtz(phase,[('PHIB','P'),('FOM','W')],np.column_stack([hkl,np.full(len(hkl),30),np.full(len(hkl),.8)]))
 def fake(merged,model,out,**kwargs):
  Path(out).mkdir();assert kwargs['phase_prefix']=='DM_';return dict(status='success')
 monkeypatch.setattr(fixed,'refmac_with_experimental_prior',fake)
 wrapper=fixed.make_fixed_prior_refiner(phase,representation='phifom',free_value=7,source_kind='density_modified_work')
 result=wrapper(inputs['native'],inputs['heavy'],tmp_path/'refine',free_value=7)
 assert 'fixed_original_experimental_prior' not in result and 'fixed_density_modified_work_prior' in result
 assert result['phase_merge_audit']['matched_free_prior_rows']==0


def test_the_hand_pilot_uses_the_in_tree_oasis_dm_engine(tmp_path,monkeypatch):
    """Both hands must be modified by the same engine the final run will use."""
    from phasebuild import oasis_dm,trace_feedback
    seen={}
    def fake_oasis_dm(work,prior,output,solvent,**options):
        seen['engine']='oasis-dm';seen['options']=options;seen['solvent']=solvent
        Path(output).write_bytes(b'x');return {'engine':'oasis-dm v2 (in-tree)'}
    monkeypatch.setattr(oasis_dm,'run',fake_oasis_dm)
    monkeypatch.setattr(trace_feedback,'merge_phase_map',lambda *a,**k:dict(work_only=True))
    phase=tmp_path/'02_sad';phase.mkdir();(phase/'observations.mtz').write_bytes(b'x')
    (tmp_path/'03_dm').mkdir()
    native=tmp_path/'native.mtz';native.write_bytes(b'x')
    cfg={'dm_engine':'oasis-dm','dm_cycles':12,'dm_cmbwt':0.4,'dm_holdout':0.05,'free_value':0,
         'runtime':str(tmp_path),'dm_solvent':0.5}
    info=hs.density_modify(cfg,phase,tmp_path/'03_dm',native,0.5)
    assert seen['engine']=='oasis-dm' and info['engine'].startswith('oasis-dm')
    assert seen['options']['cmbwt']==0.4 and seen['options']['holdout']==0.05
    assert 'mask_radius' not in seen['options'] and 'anisotropy' not in seen['options']
    record=json.loads((tmp_path/'03_dm/map_phases.json').read_text())
    assert record['dm_engine']=='oasis-dm v2'
    assert record['dm_computed_during_hand_selection'] is True
