"""Hybrid GUI plans validate immutable inputs and keep the declared endpoint."""
import json
from pathlib import Path
from types import SimpleNamespace

import gemmi
import numpy as np
import pytest

from phasebuild.data import fingerprint
from phasebuild.geometry import template
from phasebuild.gui_workflow import defaults, make_plan, collect_results, save_selected_result
from phasebuild.trace import write_backbone


def mtz(path, columns, rows, length=40.):
    value = gemmi.Mtz(with_base=True)
    value.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    value.set_cell_for_all(gemmi.UnitCell(length,length,length,90.,90.,90.))
    value.add_dataset('test')
    for label, kind in columns:
        value.add_column(label,kind)
    value.set_data(np.asarray(rows,dtype=np.float32))
    value.update_reso()
    value.write_to_file(str(path))
    return value


HL = [(name,'A') for name in ('HLA','HLB','HLC','HLD')]


@pytest.fixture
def hybrid_cfg(tmp_path):
    runtime = tmp_path/'runtime'
    (runtime/'bin').mkdir(parents=True)
    for name in ('cbuccaneer','refmac5'):
        binary=runtime/'bin'/name
        binary.touch();binary.chmod(0o755)
    library=runtime/'lib/data/reference_structures'
    library.mkdir(parents=True)
    for suffix in ('.mtz','.pdb'):
        (library/('reference-1tqw'+suffix)).touch()
    native = tmp_path/'native.mtz'
    value = mtz(native,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],
                [[1,0,0,10,1,0],[2,0,0,20,1,1]])
    maps = tmp_path/'map.mtz'
    mtz(maps,[('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('FOM','W')],
        [[1,0,0,10,1,8,30,.6],[2,0,0,20,1,15,45,.5]])
    phase = tmp_path/'prior.mtz'
    mtz(phase,HL,[[1,0,0,2,1,0,0]])
    model = tmp_path/'seed.cif'
    write_backbone(model,[template(-63.,-42.,length=3)+[20.,20.,20.]],
                   SimpleNamespace(cell=value.cell,sg=value.spacegroup))
    sequence = tmp_path/'sequence.fasta';sequence.write_text('>protein\nACDEFGHIKLMN\n')
    cfg = defaults(tmp_path)
    cfg.update(mode='hybrid_hl',map=str(maps),native=str(native),model=str(model),
               sequence=str(sequence),phase_prior=str(phase),ccp4=str(runtime),
               profile='/missing/unrelated_profile',rotamers='',residue_profile='',
               sequence_report='',refine=False,build_semet=True)
    return cfg


def test_hybrid_plan_is_one_fixed_protocol_and_reads_inputs_without_writes(hybrid_cfg,tmp_path):
    before = {path:path.read_bytes() for path in tmp_path.rglob('*') if path.is_file()}
    cfg, steps = make_plan(hybrid_cfg,tmp_path,python='/python with space')
    assert len(steps)==1 and steps[0]['kind']=='hybrid_hl'
    command = steps[0]['command']
    assert command[:2]==['/python with space','-u']
    assert Path(command[2]).name=='run_phasebuild_hybrid_hl.py'
    for flag,key in [('--map','map'),('--native','native'),('--phase-prior','phase_prior'),
                     ('--model','model'),('--sequence','sequence'),('--runtime','ccp4')]:
        assert command[command.index(flag)+1]==cfg[key]
    assert command[command.index('--rounds')+1]=='5'
    assert command[command.index('--cycles')+1]=='10'
    assert '--build-semet' in command
    assert all('profile' not in item and 'rotamer' not in item for item in command)
    assert not Path(cfg['output']).exists()
    assert before=={path:path.read_bytes() for path in tmp_path.rglob('*') if path.is_file()}
    assert defaults(tmp_path)['mode']=='fresh'


@pytest.mark.parametrize('problem',['type','duplicate','no_overlap','different_cell','no_free'])
def test_hybrid_rejects_invalid_priors_or_missing_original_free_set(hybrid_cfg,tmp_path,problem):
    if problem=='no_free':
        mtz(hybrid_cfg['native'],[('F','F'),('SIGF','Q'),('FreeR_flag','I')],
            [[1,0,0,10,1,1],[2,0,0,20,1,1]])
    else:
        columns = [('HLA','F')]+HL[1:] if problem=='type' else HL
        rows = [[9 if problem=='no_overlap' else 1,0,0,2,1,0,0]]
        if problem=='duplicate':
            rows *= 2
        mtz(hybrid_cfg['phase_prior'],columns,rows,length=41. if problem=='different_cell' else 40.)
    with pytest.raises(ValueError):
        make_plan(hybrid_cfg,tmp_path)
    assert not Path(hybrid_cfg['output']).exists()


def write_hybrid_results(tmp_path,completed=True,mode='hybrid_hl'):
    pipeline=tmp_path/('01_build_hl' if mode=='build_hl' else '01_hybrid_hl');pipeline.mkdir()
    stages=[]
    for number,rfree in [(1,.30),(2,.25),(3,.29)]:
        folder=pipeline/('round_%02d_refine'%number);folder.mkdir()
        write_backbone(folder/'model.cif',[template(-63.,-42.,length=3)+[number,20.,20.]],
                       SimpleNamespace(cell=gemmi.UnitCell(40,40,40,90,90,90),
                                       sg=gemmi.find_spacegroup_by_name('P 1')))
        statistics=dict(r_work=.20,r_free=rfree)
        report=dict(status='success',engine='external-REFMAC-experimental-phase-probe',
                    free_flags_preserved=True,requested_cycles=10,statistics={'final':statistics})
        (folder/'refinement_report.json').write_text(json.dumps(report))
        stages.append(dict(round=number,statistics=statistics))
    report=dict(status='completed' if completed else 'running',
                engine=('phasebuild-integrated-external-buccaneer-refmac' if mode=='build_hl'
                        else 'hybrid-phasebuild-seeded-external-buccaneer-refmac'),stages=stages)
    if completed:
        primary=pipeline/'primary_model.cif'
        primary.write_bytes((folder/'model.cif').read_bytes())
        report.update(primary_model=fingerprint(primary),primary_metric=statistics)
    (pipeline/'hybrid_report.json').write_text(json.dumps(report))
    (tmp_path/'gui_job.json').write_text(json.dumps(dict(config={'mode':mode},
        steps=[dict(output=str(pipeline))])))
    return pipeline


def test_hybrid_results_keep_last_round_when_earlier_rfree_is_lower(tmp_path):
    pipeline = write_hybrid_results(tmp_path)
    result = save_selected_result(tmp_path)
    assert result['selection_policy']=='fixed_final_round'
    assert result['primary_complete'] is True
    assert result['best_metric'] is None
    assert result['primary_metric']['r_free']==.29
    assert result['final_model']=='01_hybrid_hl/primary_model.cif'
    last = gemmi.read_structure(str(pipeline/'round_03_refine/model.cif'))
    assert result['preview'][0]['points'][0][0]==last[0][0][0].find_atom('CA','*').pos.x
    assert not (tmp_path/'selected_model.cif').exists()
    selection=json.loads((tmp_path/'selection.json').read_text())
    assert selection['criterion']=='fixed_final_round'
    assert selection['model']==result['final_model']
    # Imported or relocated reports retain the primary result, without relying
    # on absolute fingerprint paths recorded by the original machine.
    import shutil
    relocated=tmp_path/'relocated';shutil.copytree(pipeline,relocated)
    assert collect_results(relocated)['final_model']=='primary_model.cif'


def test_incomplete_hybrid_is_not_promoted_to_fixed_endpoint(tmp_path):
    write_hybrid_results(tmp_path,completed=False)
    result=collect_results(tmp_path)
    assert result['best_metric'] is None and not result['primary_complete']
    assert result['primary_metric']['r_free']==.29
    assert result['final_model']=='01_hybrid_hl/round_03_refine/model.cif'
    with pytest.raises(ValueError,match='未生成完整固定末轮'):
        save_selected_result(tmp_path)
    assert not (tmp_path/'selection.json').exists()


def test_changed_primary_model_is_not_displayed_with_old_metrics(tmp_path):
    pipeline=write_hybrid_results(tmp_path)
    primary=pipeline/'primary_model.cif'
    primary.write_text(primary.read_text()+'\n# changed\n')
    with pytest.raises(ValueError,match='哈希不一致'):
        collect_results(tmp_path)


@pytest.mark.parametrize('mode',['hybrid_hl','build_hl'])
@pytest.mark.parametrize('free_value',[0,7])
def test_integrated_plan_preserves_original_free_value_and_model_contract(hybrid_cfg,tmp_path,mode,free_value):
    hybrid_cfg.update(mode=mode,free_value=free_value)
    mtz(hybrid_cfg['native'],[('F','F'),('SIGF','Q'),('FreeR_flag','I')],
        [[1,0,0,10,1,free_value],[2,0,0,20,1,1]])
    if mode=='build_hl':
        # A stale hidden seed from a previous task must not enter this mode.
        hybrid_cfg['model']='/missing/stale-seed.cif'
    before={p:p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    cfg,steps=make_plan(hybrid_cfg,tmp_path)
    command=steps[0]['command']
    assert cfg['free_value']==free_value
    assert command[command.index('--free-value')+1]==str(free_value)
    assert steps[0]['kind']==mode and len(steps)==1
    if mode=='build_hl':
        assert Path(command[2]).name=='run_build_from_phases.py'
        assert '--model' not in command and cfg['model']==''
    else:
        assert Path(command[2]).name=='run_phasebuild_hybrid_hl.py'
        assert command[command.index('--model')+1]==cfg['model']
    assert not Path(cfg['output']).exists()
    assert before=={p:p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}


@pytest.mark.parametrize('mode',['hybrid_hl','build_hl'])
def test_wrong_declared_free_value_fails_before_output(hybrid_cfg,tmp_path,mode):
    hybrid_cfg.update(mode=mode,free_value=7)
    with pytest.raises(ValueError):
        make_plan(hybrid_cfg,tmp_path)
    assert not Path(hybrid_cfg['output']).exists()


def test_seeded_mode_still_requires_an_initial_model(hybrid_cfg,tmp_path):
    hybrid_cfg['model']=''
    with pytest.raises(ValueError,match='model'):
        make_plan(hybrid_cfg,tmp_path)
    assert not Path(hybrid_cfg['output']).exists()


@pytest.mark.parametrize('completed',[False,True])
def test_unseeded_import_keeps_fixed_endpoint_without_gui_manifest(tmp_path,completed):
    pipeline=write_hybrid_results(tmp_path,completed=completed,mode='build_hl')
    (tmp_path/'gui_job.json').unlink()
    result=collect_results(pipeline)
    assert result['selection_policy']=='fixed_final_round'
    assert result['best_metric'] is None
    assert result['primary_metric']['r_free']==.29
    assert result['primary_complete'] is completed
    if completed:
        assert result['final_model']=='primary_model.cif'
        assert save_selected_result(pipeline)['final_model']=='primary_model.cif'
    else:
        with pytest.raises(ValueError,match='未生成完整固定末轮'):
            save_selected_result(pipeline)
        assert not (pipeline/'selection.json').exists()
