import json
from pathlib import Path
import pytest
import phasebuild.cycling as cycling


@pytest.mark.parametrize('change_free_count',[False,True])
def test_cycle_budget_is_fixed_and_free_set_changes_abort(tmp_path,monkeypatch,change_free_count):
    start=tmp_path/'start';start.mkdir()
    initial={'keywords':'fixed work-set protocol','statistics':{
        'final':{'r_work':.3099,'r_free':.367},'reflections_used_work':5847,'reflections_used_free':299}}
    (start/'refinement_report.json').write_text(json.dumps(initial))
    (start/'model.cif').write_text('starting model')
    sequence_report=tmp_path/'sequence_report.json';sequence_report.write_text('{}')
    inputs=[]
    for name in ['native','sequence','map','profile','residue_profile','rotamers']:
        p=tmp_path/name;p.write_text('input fingerprint');inputs.append(p)
    monkeypatch.setattr(cycling,'validate_start',lambda *a:initial)
    monkeypatch.setattr(cycling,'model_counts',lambda p:{'residues':38})
    calls=[]
    def fake_refine(native,model,output,cycles):
        assert cycles==10
        calls.append((model,output))
        folder=Path(output);folder.mkdir()
        # Later rounds deliberately have worse Rfree: output must still be last.
        result={'keywords':initial['keywords'],'statistics':{
            'final':{'r_work':.30,'r_free':[.35,.40,.42][len(calls)-1]},
            'reflections_used_work':5847,'reflections_used_free':300 if change_free_count else 299}}
        (folder/'model.cif').write_text('model %d'%len(calls))
        (folder/'refinement_report.json').write_text(json.dumps(result))
        return result
    monkeypatch.setattr(cycling,'refmac_refine',fake_refine)
    arguments=(inputs[0],inputs[1],start,sequence_report,*inputs[2:])
    out=tmp_path/'experiment'
    if change_free_count:
        with pytest.raises(RuntimeError,match='reflection counts changed'):
            cycling.run_cycles(*arguments,out,mode='refine-only',progress=lambda _:None)
        assert json.loads((out/'experiment_report.json').read_text())['status']=='failed'
    else:
        report=cycling.run_cycles(*arguments,out,mode='refine-only',progress=lambda _:None)
        assert [s['cumulative_cycles'] for s in report['stages']]==[20,30,40,50]
        assert report['stages'][-1]['r_free']==.42
        assert Path(report['final_model']['path']).read_text()=='model 3'
        assert report['status']=='success' and not report['source_changed_during_run']
        assert len(calls)==3


def test_extension_cycle_uses_its_own_previous_map_and_the_extension_runner(tmp_path,monkeypatch):
    from types import SimpleNamespace
    start=tmp_path/'start';start.mkdir()
    initial={'keywords':'same','statistics':{'final':{'r_work':.30,'r_free':.36},
             'reflections_used_work':5847,'reflections_used_free':299}}
    (start/'model.cif').write_text('starting model');(start/'refined.mtz').write_text('map 0')
    (start/'refinement_report.json').write_text(json.dumps(initial))
    assignment=tmp_path/'assignment.json';assignment.write_text('{}')
    inputs=[]
    for name in ['native','sequence','map','profile','residue_profile','rotamers']:
        file=tmp_path/name;file.write_text('input');inputs.append(file)
    monkeypatch.setattr(cycling,'validate_start',lambda *a:initial)
    monkeypatch.setattr(cycling,'model_counts',lambda p:{'residues':45})
    calls=[]
    def rebuild(command,**kwargs):
        assert Path(command[1]).name=='run_extend.py'
        assert '--loop-starts' not in command and '--maximum-extension' in command
        assert Path(command[2]).read_text()=='map %d'%len(calls)
        folder=Path(command[command.index('--output')+1]);folder.mkdir()
        (folder/'model.cif').write_text('grown');(folder/'build_report.json').write_text('{}')
        return SimpleNamespace(returncode=0)
    def refine(native,model,output,cycles):
        assert Path(model).parent.name=='build' and cycles==10
        calls.append(model);folder=Path(output);folder.mkdir()
        (folder/'model.cif').write_text('refined');(folder/'refined.mtz').write_text('map %d'%len(calls))
        (folder/'refinement_report.json').write_text(json.dumps(initial))
        return initial
    monkeypatch.setattr(cycling.subprocess,'run',rebuild);monkeypatch.setattr(cycling,'refmac_refine',refine)
    args=(inputs[0],inputs[1],start,assignment,*inputs[2:])
    result=cycling.run_cycles(*args,tmp_path/'out',mode='updated-map',additional_rounds=2,
                              initial_cycles=30,rebuild_stage='extension',progress=lambda _:None)
    assert result['status']=='success' and result['rebuild_stage']=='extension'
    assert [s['cumulative_cycles'] for s in result['stages']]==[30,40,50]
    assert len(calls)==2
