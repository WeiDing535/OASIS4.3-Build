import json
from pathlib import Path
import pytest

from phasebuild.anomalous_scattering import calculate
from phasebuild import sad_pipeline, oasis_workflow
from test_sad_pipeline import inputs
from test_trust_workflow import trusted
from test_bundled_pplus import bundled


def test_known_theoretical_values_and_energy_units():
    sulfur=calculate('S',1.54)
    assert sulfur['factors']['S']['fpp']==pytest.approx(.5562888691,abs=1e-9)
    assert sulfur['energy_eV']==pytest.approx(8050.9219,abs=.001)
    selenium=calculate('Se',.71073)
    assert selenium['factors']['Se']['fpp']==pytest.approx(2.23336,abs=1e-5)
    assert selenium['wavelength_source'].startswith('explicit user input')


def test_multi_element_normalization_has_no_manual_factors():
    result=calculate('se， S,Se',.9792)
    assert result['elements']==['Se','S']
    assert result['factors']['Se']['fpp']==pytest.approx(3.8415219482,abs=1e-9)
    assert result['mode']=='automatic_theoretical'
    assert result['source_url'].startswith('https://gemmi.')


@pytest.mark.parametrize('wave',[None,True,0,-1,float('nan'),float('inf'),.09,11,'1.54'])
def test_actual_wavelength_required(wave):
    with pytest.raises(ValueError,match='波长'):
        calculate('Se',wave)


@pytest.mark.parametrize('elements',['','X','Se=6','H','He','Cf','not-an-element'])
def test_invalid_or_unsupported_elements_fail(elements):
    with pytest.raises(ValueError):
        calculate(elements,1.)


def test_table_validity_range_is_enforced():
    with pytest.raises(ValueError,match='短波长'):
        calculate('Au',.1)


def test_new_input_ignores_old_manual_value_and_records_exact_calculation(bundled):
    cfg=sad_pipeline.validate({**bundled,'elements':'S','auto_fpp':True,'wavelength':1.54,'fpp':'S=999'})
    assert cfg['fpp']==calculate('S',1.54)['fpp_argument']
    assert cfg['auto_fpp'] and cfg['elements']=='S'


def test_no_mtz_wavelength_fallback_even_when_data_exist(bundled):
    with pytest.raises(ValueError,match='波长'):
        sad_pipeline.validate({**bundled,'elements':'S','auto_fpp':True,'wavelength':None})
    assert not Path(bundled['output']).exists()


def test_gui_plan_transmits_auto_mode_without_manual_input(bundled):
    cfg={**bundled,'auto_fpp':True,'elements':'S','wavelength':1.54,
         'trust_wang_radius':8.,'trust_partial_cleanup':True,'trust_work_weight':True,'explore':False}
    checked,steps=oasis_workflow.make_plan(cfg,Path(__file__).resolve().parents[1])
    command=steps[0]['command']
    assert '--auto-fpp' in command and '--elements' in command
    assert 'True' not in command and 'False' not in command
    assert '--trust-work-weight' in command and '--trust-partial-cleanup' in command
    assert checked['auto_fpp'] and not checked['explore']
    view=oasis_workflow.describe_plan(checked)
    assert view['effective_options']['method']=='P+ Wang8 + 精修权重搜索'
    assert view['effective_options']['fpp']==checked['fpp']
    assert not Path(checked['output']).exists()


def test_new_cli_defaults_to_one_w8_weight_run(bundled,capsys):
    import importlib.util
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('fixed_cli',root/'run_sad_to_model.py')
    cli=importlib.util.module_from_spec(spec);spec.loader.exec_module(cli)
    args=['--check','--elements','S','--wavelength','1.54','--copies','2','--free-value','7']
    for key in ['sad','native','sequence','heavy','runtime','output']:
        args+=['--'+key,bundled[key]]
    assert cli.main(args)==0
    result=json.loads(capsys.readouterr().out)
    cfg=result['config']
    assert cfg['method']=='pplus_bundled' and cfg['trust_wang_radius']==8
    assert cfg['auto_fpp'] and cfg['trust_partial_cleanup'] and cfg['trust_work_weight']
    assert not cfg['trust_stability'] and not Path(cfg['output']).exists()
    assert len(result['stages'])==6


def test_explicit_legacy_cli_does_not_silently_change_old_controls(bundled,capsys):
    import importlib.util
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('legacy_cli',root/'run_sad_to_model.py')
    cli=importlib.util.module_from_spec(spec);spec.loader.exec_module(cli)
    args=['--check','--method','pplus_bundled','--wavelength','.9792','--copies','2','--free-value','7']
    for key in ['sad','native','sequence','heavy','runtime','output','fpp']:
        args+=['--'+key,bundled[key]]
    assert cli.main(args)==0
    cfg=json.loads(capsys.readouterr().out)['config']
    assert cfg['trust_wang_radius']==14 and not cfg['auto_fpp']
    assert not cfg['trust_partial_cleanup'] and not cfg['trust_work_weight']
