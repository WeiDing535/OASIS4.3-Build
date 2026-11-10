"""Focused product contracts, both real phase/DM routes and scientific handoff."""
import json
from pathlib import Path
import sys
import numpy as np
import gemmi
import pytest
from test_sad_pipeline import inputs,FakeStages
from phasebuild import sad_pipeline as pipeline,oasis_workflow as workflow
from phasebuild.data import fingerprint


def root():return Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('method',['likelihood','pplus'])
def test_both_real_phase_and_dm_routes_match_preserved_upstream(inputs,method):
    cfg={**inputs,'method':method,'solvent_fraction':.45}
    runner=FakeStages(cfg,real_phase=True)
    state=pipeline.run_pipeline(cfg,runner=runner,progress=lambda x:None)
    out=Path(cfg['output']);raw=out/'02_sad/experimental.mtz';dm=out/'03_dm/map_phases.mtz'
    meta=json.loads(dm.with_suffix('.json').read_text())
    assert meta['dm_engine']==('legacy' if method=='pplus' else 'mlhl-best')
    assert meta['phase_method']==method and meta['solvent_fraction_used']==.45 and meta['solvent_source']=='user'
    assert state['pre_dm_result']['sha256']==fingerprint(raw)['sha256']
    observations=out/'02_sad/observations.mtz'
    assert state['experimental_prior']['sha256']==fingerprint(observations)['sha256']
    build_command=pipeline.stage_plan(state['config'])[-1]['command']
    assert build_command[build_command.index('--phase-prior')+1]==str(observations)
    assert raw.read_bytes()!=dm.read_bytes()
    from oasis_py.cli import main as upstream
    ref=out/'upstream.mtz'
    assert upstream(['phase',str(out/'01_inputs/sad.mtz'),'--mode','sad','--sad-method',method,
        '--heavy',cfg['heavy'],'--sequence',cfg['sequence'],'--copies','2','--fpp',cfg['fpp'],
        '--free-label','PIPELINE_FREE','--free-value','7','--seed','20261023','--cycles','1',
        '--max-neighbors','600','--max-per-target','96','--out',str(ref)])==0
    a,b=gemmi.read_mtz_file(str(raw)),gemmi.read_mtz_file(str(ref))
    assert a.column_labels()==b.column_labels()
    np.testing.assert_array_equal(np.array(a),np.array(b))
    # DM consumes the persisted initial result, using the same upstream engine.
    from oasis_py.data import read_mtz
    from oasis_py.phasing import PhaseResult
    from oasis_py.density_mlhl import density_modify_mlhl
    from oasis_py.native_density import density_modify_native
    data,_=read_mtz(out/'01_inputs/sad.mtz','sad',free_label='PIPELINE_FREE',free_value=7)
    values=np.array(a);labels=a.column_labels();lookup={tuple(row[:3].astype(int)):row for row in values}
    ordered=np.array([lookup[tuple(h)] for h in data.hkl])
    initial=PhaseResult(np.deg2rad(ordered[:,labels.index('PHIB')]),ordered[:,labels.index('FOM')],None,
        ordered[:,[labels.index(k) for k in ['HLA','HLB','HLC','HLD']]],json.loads(raw.with_suffix('.json').read_text()))
    expected=(density_modify_native(data,initial,.45,Path(cfg['reference']),cycles=1,seed=20261023) if method=='pplus' else
              density_modify_mlhl(data,initial,.45,Path(cfg['reference']),cycles=1,seed=20261023,feedback='best'))
    actual=gemmi.read_mtz_file(str(dm));arr=np.array(actual);lab=actual.column_labels()
    np.testing.assert_allclose(arr[:,lab.index('FOM')],expected.fom,atol=1e-7,rtol=0)
    delta=np.angle(np.exp(1j*(np.deg2rad(arr[:,lab.index('PHIB')])-expected.phase)))
    np.testing.assert_allclose(delta,0,atol=1e-6,rtol=0)
    if method=='pplus':
        assert 'PPLUS' in labels
        assert 'moment_matched' in meta['density_modification']['initial_prior']
    changed={**cfg,'method':'pplus' if method=='likelihood' else 'likelihood'}
    with pytest.raises(ValueError,match='Resume refused'):pipeline.run_pipeline(changed,resume=True,runner=runner)


@pytest.mark.parametrize('update',[{'method':'mixture'},{'method':'likelihood','rank_fit':True},
    {'aoe':0},{'aoe':1},{'triplet_weight':1.1},{'solvent_fraction':0},
    {'solvent_fraction':float('nan')},{'max_neighbors':-.5},{'seed':True}])
def test_invalid_method_options_fail_before_creating_output(inputs,update):
    with pytest.raises(ValueError):pipeline.validate({**inputs,**update})
    assert not Path(inputs['output']).exists()


def test_pure_plan_rejects_removed_modes_and_preserves_method_options(inputs):
    cfg={**inputs,'mode':'sad_to_model','method':'pplus','rank_fit':True,'aoe':.6,'solvent_fraction':.4}
    checked,steps=workflow.make_plan(cfg,root())
    command=steps[0]['command'];assert '--rank-fit' in command
    assert command[command.index('--method')+1]=='pplus'
    assert command[command.index('--solvent-fraction')+1]=='0.4'
    plan=workflow.describe_plan(checked)
    assert len(plan['stages'])==5 and plan['effective_options']['dm_engine']=='Legacy'
    assert 'triplet_weight' not in plan['effective_options'] and plan['effective_options']['aoe']==.6
    assert not Path(inputs['output']).exists()
    for bad in [{'mode':'extend'},{'map':'hidden.mtz'},{'arbitrary_command':'echo unsafe'}]:
        with pytest.raises(ValueError):workflow.make_plan({**cfg,**bad},root())


def test_desktop_controller_only_accepts_the_two_method_pipeline(inputs):
    from phasebuild.desktop import Controller
    controller=Controller(root(),workflow=workflow)
    default=controller.defaults();assert default['method']==workflow.defaults(root())['method'] and 'profile' not in default
    config={**inputs,'method':'pplus'}
    cfg,steps=controller.plan(config)
    assert len(steps)==1 and steps[0]['command'][2].endswith('run_sad_to_model.py')
    assert cfg['output']==str(Path(inputs['output']))
    plan=controller.describe(cfg)
    assert len(plan['stages'])==5 and plan['effective_options']['dm_engine']=='Legacy'
    with pytest.raises(ValueError):controller.plan({**config,'mode':'refine'})
    assert not Path(inputs['output']).exists()


def test_results_keep_intermediate_separate_and_verify_final_hash(inputs):
    runner=FakeStages(inputs);pipeline.run_pipeline(inputs,runner=runner,progress=lambda x:None)
    out=Path(inputs['output']);state_path=out/'pipeline_report.json';state=json.loads(state_path.read_text())
    # Replace fake builder's minimal report with a coherent saved refinement endpoint.
    folder=out/'05_build/round_01_refine';folder.mkdir()
    (folder/'model.cif').write_bytes((out/'primary_model.cif').read_bytes())
    (folder/'refinement_report.json').write_text(json.dumps(dict(status='success',requested_cycles=1,
        free_flags_preserved=True,statistics={'final':{'r_work':.2,'r_free':.3}})))
    final=workflow.collect_results(out)
    assert final['primary_complete'] and final['metrics'][0]['r_free']==.3
    assert final['sad']['mean_fom']==.5 and final['dm']['mean_fom']==.5
    assert any(f['role']=='最终模型（固定末轮）' for f in final['artifacts'])
    interrupted=out/'03_dm.interrupted-example';interrupted.mkdir()
    (interrupted/'model.cif').write_text('invalid preserved partial model')
    assert workflow.collect_results(out)['final_model']=='primary_model.cif'
    state.update(status='failed',current='build',error='interrupted')
    state['completed']=state['completed'][:-1];state_path.write_text(json.dumps(state))
    partial=workflow.collect_results(out)
    assert not partial['primary_complete'] and partial['final_model']=='05_build/round_01_refine/model.cif'
    assert all(f['role']!='最终模型（固定末轮）' for f in partial['artifacts'])
    state.update(status='completed');state_path.write_text(json.dumps(state))
    (out/'primary_model.cif').write_bytes(b'tampered')
    with pytest.raises(ValueError,match='哈希'):workflow.collect_results(out)


def test_focused_source_runs_without_legacy_entrypoints_or_training_models(inputs,tmp_path):
    import importlib.util,subprocess,tomllib,zipfile
    spec=importlib.util.spec_from_file_location('focus_packager',root()/'package_oasis.py')
    packager=importlib.util.module_from_spec(spec);spec.loader.exec_module(packager)
    target=tmp_path/'focused.zip';packager.main(['--output',str(target)])
    with zipfile.ZipFile(target) as z:
        names=z.namelist()
        prefix='OASIS4.3-focused-development/'
        project=tomllib.loads((root()/'pyproject.toml').read_text())
        packaged=tomllib.loads(z.read(prefix+'pyproject.toml').decode())
        manifest=json.loads(z.read(prefix+'PACKAGE_MANIFEST.json'))
        assert packaged['project']['version']==manifest['version']==project['tool']['oasis43']['version']
        for key in ['requires-python','dependencies']:
            assert packaged['project'][key]==project['project'][key]
        assert not any('/gui_assets/' in n or '/models/' in n or n.endswith('/gui_workflow.py') for n in names)
        assert any(n.endswith('/phasebuild/desktop.py') for n in names) and not any(n.endswith('/phasebuild/gui.py') for n in names)
        assert not any(n.endswith('/run_phasebuild.py') or n.endswith('/run_refine.py') or n.endswith('/oasis_py/cli.py') for n in names)
        assert any(n.endswith('/run_desktop.py') for n in names)
        assert any(n.endswith('/clipper-2.1.20201109.tar.gz') for n in names)
        z.extractall(tmp_path/'source')
    extracted=tmp_path/'source/OASIS4.3-focused-development'
    argv=[sys.executable,'-B',str(extracted/'run_sad_to_model.py'),'--check','--method','pplus']
    for key,value in inputs.items():argv+=['--'+key.replace('_','-'),str(value)]
    result=subprocess.run(argv,cwd=extracted,text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert json.loads(result.stdout)['config']['method']=='pplus'
    probe=subprocess.run([sys.executable,'-B','-c',
        "import sys;sys.path.insert(0,'src');from phasebuild.desktop import Controller;from phasebuild import oasis_workflow;c=Controller('.',workflow=oasis_workflow);assert 'gui_workflow' not in sys.modules;assert c.defaults()['method'];print('focused imports passed')"],
        cwd=extracted,text=True,capture_output=True)
    assert probe.returncode==0,probe.stderr
    assert not Path(inputs['output']).exists()


@pytest.mark.parametrize('method',['likelihood','pplus'])
def test_phases_only_completes_without_resolving_or_fingerprinting_external_engines(inputs,method,monkeypatch):
    cfg={**inputs,'method':method,'build_model':False,'runtime':'/missing/unused/runtime',
         'rounds':None,'cycles':None}
    def forbidden(*a,**kw):raise AssertionError('Phases-only must never resolve external engines')
    monkeypatch.setattr(pipeline,'resolve_runtime',forbidden)
    monkeypatch.setattr(pipeline,'runtime_info',forbidden)
    runner=FakeStages(cfg,real_phase=True)
    state=pipeline.run_pipeline(cfg,runner=runner,progress=lambda x:None)
    out=Path(cfg['output'])
    assert state['status']=='completed' and state['primary_model'] is None and state['runtime'] is None
    assert runner.calls==['phase','dm'] and len(state['completed'])==3
    assert 'cbuccaneer' not in state['inputs'] and 'runtime_resources' not in state['inputs']
    assert not (out/'04_bridge').exists() and not (out/'05_build').exists()
    assert not (out/'01_inputs/building.fasta').exists()
    result=workflow.save_selected_result(out)
    assert result['pipeline_complete'] and not result['primary_complete'] and not result['build_model']
    assert len(result['stages'])==3 and result['final_model'] is None and result['metrics']==[]
    assert result['final_phases']=='03_dm/map_phases.mtz'
    assert any(f['group']=='主要结果' and f['path']==result['final_phases'] for f in result['artifacts'])
    pipeline.run_pipeline({**cfg,'runtime':'/another/irrelevant/path'},resume=True,runner=runner)
    assert runner.calls==['phase','dm']
    (out/'03_dm/map_phases.mtz').write_bytes(b'tampered')
    with pytest.raises(ValueError,match='哈希'):workflow.collect_results(out)


def test_optional_build_flag_in_preflight_and_resume_protocol(inputs):
    assert pipeline.validate(inputs)['build_model'] is True
    cfg={**inputs,'build_model':False,'runtime':''}
    checked,steps=workflow.make_plan(cfg,root())
    assert '--no-build-model' in steps[0]['command'] and '--runtime' not in steps[0]['command']
    plan=workflow.describe_plan(checked)
    assert len(plan['stages'])==3 and plan['effective_options']['build_model'] is False
    assert 'rounds' not in plan['effective_options'] and 'refmac_cycles' not in plan['effective_options']
    assert not Path(inputs['output']).exists()
    pipeline.run_pipeline(cfg,runner=FakeStages(cfg),progress=lambda x:None)
    with pytest.raises(ValueError,match='Resume refused'):
        pipeline.run_pipeline({**inputs,'build_model':True},resume=True,runner=FakeStages(inputs))
    with pytest.raises(ValueError,match='boolean'):pipeline.validate({**inputs,'build_model':'false'})


def test_desktop_phases_only_job_is_successful_and_reports_three_stages(inputs,tmp_path):
    from phasebuild.jobs import Jobs
    project=tmp_path/'project';project.mkdir()
    (project/'run_sad_to_model.py').symlink_to(root()/'run_sad_to_model.py')
    jobs=Jobs(project,workflow=workflow)
    job=jobs.start({**inputs,'build_model':False,'runtime':''})
    jobs.worker.join(timeout=25)
    assert not jobs.worker.is_alive()
    saved=jobs.get(job['id'])
    assert saved['status']=='success',saved.get('log')
    assert saved['pipeline_progress']['total']==saved['pipeline_progress']['completed']==3
    assert len(saved['plan']['stages'])==3
    assert workflow.collect_results(job['output'])['final_phases']=='00_sad_to_model/03_dm/map_phases.mtz'


def test_fixed_endpoint_grid_reports_its_retained_endpoint(tmp_path):
    """--feedback-endpoints keeps one pipeline per endpoint; the front ends must
    resolve the selection instead of reporting an empty case root."""
    base=tmp_path/'run'/'00_sad_to_model'
    for tag,fb,rwork in [('endpoint_fb2_s0.50',2,0.2834),('endpoint_fb3_s0.50',3,0.2911)]:
        folder=base/tag
        folder.mkdir(parents=True)
        (folder/'pipeline_report.json').write_text(json.dumps(
            {'engine':'oasis-sad-to-model','status':'completed','config':{'method':'pplus_bundled','rounds':1},
             'completed':[],'primary_metric':{'r_work':rwork,'r_free':rwork+0.03}}))
    (base/'endpoint_selection.json').write_text(json.dumps(
        {'rule':'lowest final work-set R across the predeclared endpoint grid; R-free never read',
         'endpoints':[{'feedback':2,'solvent':0.5,'status':'completed','output':str(base/'endpoint_fb2_s0.50'),
                       'primary_metric':{'r_work':0.2834,'r_free':0.3147}},
                      {'feedback':3,'solvent':0.5,'status':'completed','output':str(base/'endpoint_fb3_s0.50'),
                       'primary_metric':{'r_work':0.2911,'r_free':0.3273}}],
         'selected_feedback':2,'selected_solvent':0.5,'selected_metric':{'r_work':0.2834,'r_free':0.3147}}))
    chosen,selection=workflow.selected_endpoint(tmp_path/'run')
    assert chosen.name=='endpoint_fb2_s0.50' and selection['selected_feedback']==2
    result=workflow.collect_results(tmp_path/'run')
    assert result['endpoint_root'].endswith('endpoint_fb2_s0.50')
    assert result['root']==str((tmp_path/'run').resolve())
    assert result['status']=='completed' and result['pipeline_complete'] is True
    assert len(result['endpoints'])==2
    assert 'work-set R' in result['selection_policy']
    # Without a selection file the newest endpoint is reported while running.
    (base/'endpoint_selection.json').unlink()
    assert workflow.selected_endpoint(tmp_path/'run')==(None,{})
    job=workflow.enrich_job({'output':str(tmp_path/'run'),'config':{}})
    assert job['endpoint'].startswith('endpoint_fb')
