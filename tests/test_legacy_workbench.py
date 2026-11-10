"""Planning and result collection for the archived PhaseBuild workbench."""
import json
from pathlib import Path
import pytest
from phasebuild.gui_workflow import defaults,make_plan,collect_results


@pytest.fixture
def cfg(tmp_path):
    c=defaults(tmp_path)
    data=tmp_path/'input with spaces.mtz'
    data.write_text('fixture')
    for key in ('map','native','sequence','model','sequence_report','profile','residue_profile','rotamers'):
        c[key]=str(data)
    binary=tmp_path/'ccp4/bin/refmac5'
    binary.parent.mkdir(parents=True)
    binary.touch()
    c['ccp4']=str(binary.parents[1])
    return c


def test_fresh_plan_connects_reports_and_updated_maps(cfg,tmp_path):
    cfg.update(map_mode='updated',rounds=3)
    _,steps=make_plan(cfg,tmp_path,inspect=False)
    assert [s['kind'] for s in steps]==['trace','build','build','refine','build','refine','build','refine']
    assert str(Path(steps[3]['output'])/'refined.mtz') in steps[4]['command']
    assert str(Path(steps[2]['output'])/'build_report.json') in steps[4]['command']
    assert not Path(cfg['output']).exists()
    # Spaces remain one argument, with no shell quoting or interpolation.
    assert cfg['map'] in steps[0]['command']


def test_extension_fixed_map_uses_new_sequence_report(cfg,tmp_path):
    cfg.update(mode='extend',map_mode='fixed',rounds=2)
    _,steps=make_plan(cfg,tmp_path,inspect=False)
    assert cfg['map'] in steps[0]['command'] and cfg['map'] in steps[2]['command']
    assert str(Path(steps[0]['output'])/'build_report.json') in steps[2]['command']
    assert str(Path(steps[1]['output'])/'model.cif') in steps[2]['command']


def test_build_without_refine_has_no_native_requirement(cfg,tmp_path):
    cfg.update(mode='extend',refine=False,rounds=1,native='',ccp4='missing')
    _,steps=make_plan(cfg,tmp_path,inspect=False)
    assert len(steps)==1 and steps[0]['kind']=='build'


def test_refine_plan_does_not_need_maps_or_sequence(cfg,tmp_path):
    cfg.update(mode='refine',map='',sequence='',model=cfg['model'])
    _,steps=make_plan(cfg,tmp_path,inspect=False)
    assert len(steps)==1
    assert '--free-column' in steps[0]['command']


def test_evaluation_needs_no_ccp4(cfg,tmp_path):
    cfg.update(mode='evaluate',ccp4='/does/not/exist',map='',sequence='')
    _,steps=make_plan(cfg,tmp_path,inspect=False)
    assert len(steps)==1 and steps[0]['kind']=='evaluate'
    assert 'run_rfactor.py' in steps[0]['command'][2]


def test_selected_map_columns_propagate_then_switch_to_refined_labels(cfg,tmp_path):
    cfg.update(map_amplitude='MY_MAP',map_phase='MY_PHASE',rounds=2,map_mode='updated')
    _,steps=make_plan(cfg,tmp_path,inspect=False)
    assert 'MY_MAP' in steps[0]['command'] and 'MY_PHASE' in steps[0]['command']
    assert 'MY_MAP' in steps[1]['command'] and 'MY_PHASE' in steps[2]['command']
    assert 'FWT' in steps[4]['command'] and 'PHWT' in steps[4]['command']
    assert 'MY_MAP' not in steps[4]['command']


@pytest.mark.parametrize('update',[{'rounds':1.5},{'rounds':0},{'selection_cutoff':float('nan')},
                                  {'mode':'shell'},{'arbitrary_command':'echo no'},
                                  {'refine':False,'rounds':2}])
def test_invalid_config_rejected_before_output(cfg,tmp_path,update):
    cfg.update(update)
    with pytest.raises(ValueError):
        make_plan(cfg,tmp_path,inspect=False)
    assert not Path(cfg['output']).exists()


def test_existing_output_never_overwritten(cfg,tmp_path):
    output=Path(cfg['output']);output.mkdir(parents=True)
    marker=output/'keep';marker.write_text('preserved')
    with pytest.raises(ValueError,match='已存在'):
        make_plan(cfg,tmp_path,inspect=False)
    assert marker.read_text()=='preserved'


def test_complete_fresh_schedule_uses_original_then_updated_map(cfg,tmp_path):
    cfg.update(map_amplitude='ORIG_F',map_phase='ORIG_P')
    actual,steps=make_plan(cfg,tmp_path,inspect=False)
    assert actual['rounds']==5 and len(steps)==12
    assert '--polish-terminals' not in steps[1]['command']
    assert Path(steps[4]['command'][2]).name=='run_rebuild.py'
    assert cfg['map'] in steps[4]['command'] and cfg['map'] in steps[6]['command']
    assert 'ORIG_F' in steps[4]['command'] and 'ORIG_P' in steps[6]['command']
    assert str(Path(steps[2]['output'])/'build_report.json') in steps[4]['command']
    assert str(Path(steps[4]['output'])/'build_report.json') in steps[6]['command']
    assert str(Path(steps[3]['output'])/'model.cif') in steps[4]['command']
    assert str(Path(steps[7]['output'])/'refined.mtz') in steps[8]['command']
    assert 'FWT' in steps[8]['command'] and 'PHWT' in steps[8]['command']
    assert sum(s['kind']=='refine' for s in steps)==5
    assert all('reference' not in str(s['command']) for s in steps)


def test_selection_keeps_best_and_last_models_distinct(tmp_path):
    import gemmi
    from phasebuild.gui_workflow import save_selected_result
    steps=[]
    for n,rfree in [(1,.39),(2,.28),(3,.31)]:
        folder=tmp_path/('round_%02d_refine'%n);folder.mkdir()
        st=gemmi.Structure();st.cell=gemmi.UnitCell(30,30,30,90,90,90);st.spacegroup_hm='P 1'
        model=gemmi.Model('1');chain=gemmi.Chain('A');res=gemmi.Residue();res.name='GLY';res.seqid=gemmi.SeqId(1,' ')
        atom=gemmi.Atom();atom.name='CA';atom.element=gemmi.Element('C');atom.pos=gemmi.Position(n,2,3)
        res.add_atom(atom);chain.add_residue(res);model.add_chain(chain);st.add_model(model)
        st.make_mmcif_document().write_file(str(folder/'model.cif'))
        report=dict(status='success',engine='external-REFMAC5',free_flags_preserved=True,
                    statistics={'final':{'r_work':.2,'r_free':rfree}})
        (folder/'refinement_report.json').write_text(json.dumps(report))
        steps.append(dict(output=str(folder)))
    (tmp_path/'gui_job.json').write_text(json.dumps({'steps':steps}))
    result=save_selected_result(tmp_path)
    assert result['best_metric']['r_free']==.28
    assert result['final_model']=='round_02_refine/model.cif'
    assert result['last_model']=='round_03_refine/model.cif'
    assert result['preview'][0]['points'][0][0]==2
    assert (tmp_path/'selected_model.cif').read_bytes()==(tmp_path/result['final_model']).read_bytes()
    assert len(result['metrics'])==3
    assert collect_results(tmp_path)['final_model']==result['final_model']
