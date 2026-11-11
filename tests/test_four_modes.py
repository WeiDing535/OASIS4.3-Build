import pytest
from phasebuild import oasis_workflow as workflow
from phasebuild.portfolio import PROTOCOLS

@pytest.mark.parametrize('name', list(PROTOCOLS))
def test_profile_serializes_actual_pipeline_flags(tmp_path, monkeypatch, name):
    monkeypatch.setattr(workflow, 'validate', lambda cfg, **kw: dict(cfg))
    cfg, steps = workflow.make_plan(dict(protocol=name, output=str(tmp_path/'run')),tmp_path)
    command=steps[0]['command']
    assert cfg['method']=='pplus_bundled'
    assert '--protocol' not in command
    assert command[command.index('--trust-wang-radius')+1]==str(PROTOCOLS[name]['trust_wang_radius'])
    for key in ['trust_stability','trust_work_weight']:
        assert ('--'+('' if PROTOCOLS[name][key] else 'no-')+key.replace('_','-')) in command
    assert not (tmp_path/'run').exists()

def test_phase_only_keeps_profile_without_post_build_weight_search(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow, 'validate', lambda cfg, **kw: dict(cfg))
    cfg,steps=workflow.make_plan(dict(protocol='pplus_w8_weight',build_model=False,output=str(tmp_path/'run')),tmp_path)
    assert cfg['protocol']=='pplus_w8_weight'
    assert '--no-trust-work-weight' in steps[0]['command']
    assert '--no-trust-partial-cleanup' in steps[0]['command']

def test_unknown_profile_fails_before_creating_output(tmp_path):
    with pytest.raises(ValueError,match='未知计算模式'):
        workflow.make_plan(dict(protocol='unknown',output=str(tmp_path/'run')),tmp_path)
    assert not (tmp_path/'run').exists()

def test_profile_cannot_silently_replace_continuous_sad(tmp_path):
    with pytest.raises(ValueError,match='仅适用于自带 P+'):
        workflow.make_plan(dict(method='likelihood',protocol='pplus_w8_weight',output=str(tmp_path/'run')),tmp_path)
    assert not (tmp_path/'run').exists()

@pytest.mark.parametrize('method', ['likelihood','pplus','pplus_internal'])
def test_original_method_routes_are_preserved(tmp_path,monkeypatch,method):
    monkeypatch.setattr(workflow,'validate',lambda cfg,**kw:dict(cfg))
    cfg,steps=workflow.make_plan(dict(method=method,output=str(tmp_path/'run'),auto_hand=False,refinement_prior='original_sad',trust_wang_radius=14,trust_partial_cleanup=False,trust_work_weight=False),tmp_path)
    command=steps[0]['command']
    assert command[command.index('--method')+1]==method
    assert '--protocol' not in command and not cfg['explore']
    assert cfg['refinement_prior']=='original_sad'

def test_scientific_corrections_serialize_without_changing_phase_method(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow,'validate',lambda cfg,**kw:dict(cfg))
    cfg,steps=workflow.make_plan(dict(method='pplus_bundled',protocol='pplus_w8_weight',canonical_dm_seed=True,phase_blur_scale=.7,phase_blur_b=20.,output=str(tmp_path/'run')),tmp_path)
    command=steps[0]['command']
    assert '--canonical-dm-seed' in command
    assert command[command.index('--phase-blur-scale')+1]=='0.7'
    assert command[command.index('--phase-blur-b')+1]=='20.0'
    assert command[command.index('--method')+1]=='pplus_bundled'
    assert 'True' not in command and 'False' not in command
