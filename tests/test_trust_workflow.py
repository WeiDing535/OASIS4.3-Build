"""Executable GUI protocol contracts for original-likelihood P+ integration."""
import json
from pathlib import Path
import shutil
import gemmi
import numpy as np
import pytest
from test_sad_pipeline import inputs,mtz,write_phases,FakeStages
from phasebuild import sad_pipeline as pipeline,oasis_workflow as workflow,trust_workflow as trust
from phasebuild.sad_bridge import prepare_sad


def test_cleanup_work_fit_guard_ignores_rfree():
    raw={'r_work': .2917, 'r_free': .3198}
    harmful={'r_work': .4503, 'r_free': .1000}
    acceptable={'r_work': .3000, 'r_free': .9000}
    assert not trust.accept_cleanup_by_work_fit(raw, harmful)
    assert trust.accept_cleanup_by_work_fit(raw, acceptable)


@pytest.fixture
def trusted(inputs,tmp_path):
    cfg={**inputs,'method':'pplus_trust','reference':'','trust_kl':.0005}
    n=gemmi.read_mtz_file(cfg['native']);a=np.array(n);f,s,flags=a[:,3],a[:,4],a[:,5]
    d=np.sin(np.arange(len(a)))
    mtz(Path(cfg['sad']),[('F(+)','G'),('SIGF(+)','L'),('F(-)','G'),('SIGF(-)','L')],
        np.column_stack([a[:,:3],f+d/2,s*np.sqrt(2),f-d/2,s*np.sqrt(2)]))
    prior=tmp_path/'observations.mtz';write_phases(prior,cfg['native'])
    cfg['observation_prior']=str(prior)
    bins=tmp_path/'phenix';bins.mkdir()
    for name in ['phenix.python','mmtbx.density_modification']:
        p=bins/name;p.write_text('#!/bin/sh\nexit 0\n');p.chmod(0o755)
    cfg['phenix_bin']=str(bins)
    return cfg


def test_new_gui_plan_preserves_explicit_route_and_legacy_missing_method(trusted,inputs):
    cfg,steps=workflow.make_plan(trusted,Path(__file__).resolve().parents[1])
    command=steps[0]['command']
    assert command.count('--trust-cleanup')==1 and 'True' not in command and 'False' not in command
    assert '--build-semet' not in command
    assert cfg['phase_cycles']==100 and cfg['dm_cycles']==12
    plan=workflow.describe_plan(cfg)
    assert [s['id'] for s in plan['stages']]==['prepare','sad','dm','bridge','build','finish']
    assert plan['effective_options']['trust_kl']==.0005 and 'aoe' not in plan['effective_options']
    builder=plan['stages'][4]['command'];assert builder[builder.index('--phase-prior')+1].endswith('/02_sad/observations.mtz')
    assert not Path(trusted['output']).exists()
    # A pre-existing config without method must retain the old continuous route.
    old={**inputs,'sad':inputs['sad'],'method':'likelihood'}
    old.pop('method');cfg,_=workflow.make_plan(old,Path(__file__).resolve().parents[1],inspect=False)
    assert cfg['method']=='likelihood'


@pytest.mark.parametrize('change',[
    {'trust_kl':0},{'trust_kl':True},{'trust_kl':float('nan')},{'trust_kl':.1},
    {'observation_prior':''},{'trust_cleanup':'yes'},
    {'build_semet':True},{'labels':'fplus=other'}])
def test_new_route_invalid_inputs_fail_without_creating_output(trusted,change):
    with pytest.raises(ValueError):pipeline.validate({**trusted,**change})
    assert not Path(trusted['output']).exists()


def test_prior_wrong_hkl_or_downstream_columns_refused(trusted):
    path=Path(trusted['observation_prior']);m=gemmi.read_mtz_file(str(path));a=np.array(m)
    m.add_column('FWT','F');m.set_data(np.column_stack([a,np.ones(len(a))]));m.write_to_file(str(path))
    with pytest.raises(ValueError,match='下游'):pipeline.validate(trusted)
    m.remove_column(m.column_labels().index('FWT'));m.set_data(np.array(m)[:-1]);m.write_to_file(str(path))
    with pytest.raises(ValueError,match='HKL sets differ'):pipeline.validate(trusted)


def test_explicit_semet_composition_matches_actual_build_sequence(trusted,tmp_path):
    Path(trusted['sequence']).write_text('>protein\nMAAAAAAAAAAA\n')
    comp=tmp_path/'composition.fasta';comp.write_text('>protein\n(MSE)AAAAAAAAAAA\n')
    cfg=pipeline.validate({**trusted,'composition_sequence':str(comp),'build_semet':True})
    info=trust.sample_info(cfg,gemmi.read_mtz_file(cfg['sad']))
    assert info['asu_composition']['Se']==2 and info['asu_residues']==24
    comp.write_text('>protein\n(MSE)CAAAAAAAAAA\n')
    with pytest.raises(ValueError,match='逐链一致'):pipeline.validate({**trusted,'composition_sequence':str(comp),'build_semet':True})


def test_real_new_pplus_exports_moments_only_and_exact_original_hl(trusted,tmp_path):
    cfg=pipeline.validate(trusted);prepared=tmp_path/'prepared.mtz'
    prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg['prepared_sad']=str(prepared)
    out=tmp_path/'phase';out.mkdir();original=Path(cfg['observation_prior']).read_bytes()
    trust.phase(cfg,out)
    assert (out/'observations.mtz').read_bytes()==original
    seed=gemmi.read_mtz_file(str(out/'experimental.mtz'))
    assert seed.column_labels()==['H','K','L','PHIB','FOM','PPLUS']
    assert len(seed.array)>0 and np.isfinite(seed.array).all()
    meta=json.loads((out/'experimental.json').read_text())
    assert meta['posterior_HL_exported'] is False and meta['original_observation_HL_modified'] is False
    assert meta['mean_abs_probability_change']>0


def make_dm(path,cfg):
    n=gemmi.read_mtz_file(cfg['native']);a=np.array(n)
    # Deliberately different DM amplitude scale and map phase from centroid.
    mtz(path,[('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHIB','P'),('FOM','W')],
        np.column_stack([a[:,:3],a[:,3]*2,a[:,4]*2,a[:,3]*.8,np.full(len(a),60),np.full(len(a),90),np.full(len(a),.3)])[::-1])


def test_actual_dm_amplitude_backscale_and_complete_initial_free_mask(trusted,tmp_path):
    cfg=pipeline.validate(trusted);out=tmp_path/'04_bridge';out.mkdir();(tmp_path/'03_dm').mkdir()
    make_dm(tmp_path/'03_dm/map_phases.mtz',cfg);before=Path(cfg['native']).read_bytes()
    trust.bridge(cfg,out);m=gemmi.read_mtz_file(str(out/'building.mtz'));a=np.array(m);labels=m.column_labels();work=a[:,labels.index('FreeR_flag')]!=7
    for k in ['F','SIGF','FWT','PHWT','PHCOMB','FOM']:assert np.isnan(np.asarray(m.column_with_label(k))[~work]).all()
    np.testing.assert_allclose(np.asarray(m.column_with_label('FWT'))[work],a[work,labels.index('F')]*.4,rtol=1e-6)
    np.testing.assert_array_equal(np.asarray(m.column_with_label('PHWT'))[work],60)
    np.testing.assert_array_equal(np.asarray(m.column_with_label('PHCOMB'))[work],90)
    assert Path(cfg['native']).read_bytes()==before


def test_six_stage_pipeline_retains_raw_model_on_conservative_cleanup_refusal(trusted):
    cfg=pipeline.validate(trusted)
    fake=FakeStages(cfg)
    def runner(command,log,env):
        if command[2].endswith('run_build_from_phases.py'):action='build'
        elif command[2].endswith('run_trace_feedback.py'):action='dm'
        else:action=command[3]
        out=Path(command[command.index('--output')+1])
        if action=='phase':
            prepared=out.parent/'01_inputs/sad.mtz';trust.phase({**cfg,'prepared_sad':str(prepared)},out)
        elif action=='dm':
            make_dm(out/'map_phases.mtz',cfg);trust.dump(out/'map_phases.json',{})
        elif action=='bridge':trust.bridge(cfg,out)
        elif action=='finish':trust.finish(cfg,out)
        else:
            fake(command,log,env);r=out/'round_01_refine';r.mkdir();shutil.copyfile(out/'primary_model.cif',r/'model.cif')
            trust.dump(r/'refinement_report.json',{'status':'success','statistics':{'final':{'r_work':.2,'r_free':.3}},'free_flags_preserved':True,'requested_cycles':1})
    state=pipeline.run_pipeline(cfg,runner=runner,progress=lambda s:None)
    assert state['status']=='completed' and len(state['completed'])==6
    result=workflow.collect_results(cfg['output'])
    assert result['primary_complete'] and result['cleanup_status']=='refused'
    assert 'single-copy' in result['cleanup_reason'] and result['primary_metric']['r_free']==.3
    # Every source/prior identity remains part of the restart boundary.
    with pytest.raises(ValueError,match='Resume refused'):
        pipeline.run_pipeline({**cfg,'trust_kl':.001},resume=True,runner=runner)


def test_stage_cli_refuses_nonempty_output_without_touching_artifacts(tmp_path):
    import subprocess,sys
    config=tmp_path/'pipeline_report.json';config.write_text(json.dumps({'config':{}}))
    out=tmp_path/'stage';out.mkdir();marker=out/'original.mtz';marker.write_bytes(b'preserve this exact existing artifact')
    script=Path(__file__).resolve().parents[1]/'run_pplus_trust_stage.py'
    result=subprocess.run([sys.executable,'-B',str(script),'phase','--config',str(config),'--output',str(out)],capture_output=True,text=True)
    assert result.returncode!=0 and 'new, empty stage directory' in result.stderr
    assert marker.read_bytes()==b'preserve this exact existing artifact'


def test_stability_opt_in_runs_real_updates_without_changing_original_HL(trusted,tmp_path):
    cfg=pipeline.validate({**trusted,'trust_stability':True});prepared=tmp_path/'stable-prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg['prepared_sad']=str(prepared)
    output=tmp_path/'stable-phase';output.mkdir();before=Path(cfg['observation_prior']).read_bytes();trust.phase(cfg,output)
    meta=json.loads((output/'experimental.json').read_text())
    assert meta['method']=='disjoint-canonical-triplet-stability-shrunk-Pplus'
    assert meta['mean_abs_probability_change']>0 and meta['stable_update_rows']>0
    assert meta['max_KL_from_original']<=cfg['trust_kl']+1e-12
    assert (output/'observations.mtz').read_bytes()==before
    assert not meta['posterior_HL_exported']


def test_new_research_flags_are_explicit_boolean_CLI_options(trusted):
    cfg,steps=workflow.make_plan({**trusted,'trust_stability':True,'trust_partial_cleanup':True},Path(__file__).resolve().parents[1])
    command=steps[0]['command'];assert '--trust-stability' in command and '--trust-partial-cleanup' in command
    assert 'True' not in command and 'False' not in command
    assert cfg['trust_stability'] and cfg['trust_partial_cleanup']
    for invalid in [{'trust_stability':'yes'},{'trust_partial_cleanup':1},{'trust_partial_cleanup':True,'trust_cleanup':False},{'trust_partial_cleanup':True,'build_model':False}]:
        with pytest.raises(ValueError):pipeline.validate({**trusted,**invalid})
