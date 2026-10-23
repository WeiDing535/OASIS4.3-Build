"""Focused product contracts, both real phase/DM routes and scientific handoff."""
import json
from pathlib import Path
import sys
import threading
import urllib.request
import urllib.error
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
    assert state['experimental_prior']['sha256']==fingerprint(raw)['sha256']
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


def test_focused_server_only_accepts_the_two_method_pipeline(inputs):
    from phasebuild.gui import create_server
    server=create_server(root(),workflow=workflow)
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    def call(route,data=None):
        req=urllib.request.Request('http://'+server.address+'/api/'+route,
            data=None if data is None else json.dumps(data).encode(),headers={'X-PhaseBuild-Token':server.token,'Content-Type':'application/json'})
        with urllib.request.urlopen(req) as response:return json.load(response)
    try:
        default=call('defaults')['config'];assert default['method']=='likelihood' and 'profile' not in default
        result=call('plan',{**inputs,'method':'pplus'});assert len(result['stages'])==5
        with pytest.raises(urllib.error.HTTPError) as err:call('plan',{**inputs,'mode':'refine'})
        assert err.value.code==400 and not Path(inputs['output']).exists()
    finally:server.shutdown();server.server_close();thread.join(timeout=3)


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
    import importlib.util,subprocess,zipfile
    spec=importlib.util.spec_from_file_location('focus_packager',root()/'package_oasis.py')
    packager=importlib.util.module_from_spec(spec);spec.loader.exec_module(packager)
    target=tmp_path/'focused.zip';packager.main(['--output',str(target)])
    with zipfile.ZipFile(target) as z:
        names=z.namelist()
        assert not any('/gui_assets/' in n or '/models/' in n or n.endswith('/gui_workflow.py') for n in names)
        assert not any(n.endswith('/run_phasebuild.py') or n.endswith('/run_refine.py') or n.endswith('/oasis_py/cli.py') for n in names)
        assert any(n.endswith('/oasis_assets/app.js') for n in names)
        assert any(n.endswith('/clipper-2.1.20201109.tar.gz') for n in names)
        z.extractall(tmp_path/'source')
    extracted=tmp_path/'source/OASIS4.3-focused-development'
    argv=[sys.executable,'-B',str(extracted/'run_sad_to_model.py'),'--check','--method','pplus']
    for key,value in inputs.items():argv+=['--'+key.replace('_','-'),str(value)]
    result=subprocess.run(argv,cwd=extracted,text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert json.loads(result.stdout)['config']['method']=='pplus'
    probe=subprocess.run([sys.executable,'-B','-c',
        "import sys;sys.path.insert(0,'src');from phasebuild.gui import create_server;from phasebuild import oasis_workflow; s=create_server('.',workflow=oasis_workflow);assert 'gui_workflow' not in sys.modules;s.server_close();print('focused imports passed')"],
        cwd=extracted,text=True,capture_output=True)
    assert probe.returncode==0,probe.stderr
    assert not Path(inputs['output']).exists()
