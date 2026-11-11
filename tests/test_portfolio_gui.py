import copy
import json
from pathlib import Path

import pytest

from phasebuild import oasis_workflow as workflow, portfolio_gui as adapter
from phasebuild.data import fingerprint
from phasebuild.portfolio import rank_candidates
from test_portfolio import candidate


def test_bundled_gui_defaults_to_fixed_w8_weight_with_auto_factors(tmp_path):
    for name in ['bin/cbuccaneer', 'bin/refmac5', 'bin/phaser', 'density/bin/oasis-density-python']:
        path = tmp_path / 'runtime' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    cfg = workflow.defaults(tmp_path)
    assert cfg['method'] == 'pplus_bundled'
    assert cfg['explore'] is False and cfg['auto_fpp'] is True
    assert cfg['trust_wang_radius'] == 8 and cfg['trust_partial_cleanup'] and cfg['trust_work_weight']
    assert not cfg['trust_stability'] and cfg['wavelength'] is None


def test_gui_plan_uses_saved_job_config_and_creates_no_files(tmp_path, monkeypatch):
    from phasebuild import portfolio
    monkeypatch.setattr(workflow, 'validate', lambda cfg, **kw: dict(cfg))
    seen = []
    monkeypatch.setattr(portfolio, 'plan', lambda cfg, output, **kw: seen.append((cfg, output, kw)))
    cfg = dict(method='pplus_bundled', explore=True, explore_workers=1, output=str(tmp_path / 'run'))
    checked, steps = workflow.make_plan(cfg, tmp_path, inspect=False)
    assert checked['explore'] and len(steps) == 1
    assert steps[0]['kind'] == 'sad_portfolio'
    assert '--gui-config' in steps[0]['command']
    assert str(tmp_path / 'run/gui_job.json') in steps[0]['command']
    assert not (tmp_path / 'run').exists()
    assert seen[0][2]['inspect'] is False


def test_historical_saved_config_keeps_single_run(tmp_path, monkeypatch):
    monkeypatch.setattr(workflow, 'validate', lambda cfg, **kw: dict(cfg))
    checked, steps = workflow.make_plan({'method': 'pplus_bundled', 'output': str(tmp_path / 'run')}, tmp_path)
    assert checked['explore'] is False
    assert steps[0]['kind'] == 'sad_to_model'
    assert '--explore' not in steps[0]['command'] and '--explore-workers' not in steps[0]['command']


@pytest.mark.parametrize('update', [{'method': 'likelihood'}, {'build_model': False}, {'explore_workers': True}, {'explore_workers': 3}, {'explore': 'true'}])
def test_invalid_gui_exploration_rejected_before_output(tmp_path, update):
    cfg = dict(method='pplus_bundled', explore=True, output=str(tmp_path / 'run'))
    with pytest.raises(ValueError):
        workflow.make_plan({**cfg, **update}, tmp_path)
    assert not (tmp_path / 'run').exists()


def make_completed(tmp_path, monkeypatch, *, ambiguous=False, identical=False):
    pipe = tmp_path / '00_portfolio'
    pipe.mkdir()
    rows, branches = [], []
    for name, score in [('pplus_w14', .25), ('pplus_w8', .28)]:
        path = pipe / name
        path.mkdir()
        model = path / 'primary_model.cif'
        model.write_text('placeholder model ' + ('identical' if identical else name))
        row = candidate(name, work=.25 if ambiguous or identical else score)
        row.update(directory=str(path), model=fingerprint(model))
        (path / 'pipeline_report.json').write_text(json.dumps({'status': 'completed', 'primary_metric': row['primary_metric']}))
        rows.append(row)
        branches.append({'name': name, 'config': {'output': str(path)}})
    report = rank_candidates(rows)
    report.update(status='completed', inputs_unchanged=True)
    (pipe / 'PORTFOLIO_REPORT.json').write_text(json.dumps(report))
    (pipe / 'PLAN.json').write_text(json.dumps({'branches': branches}))
    monkeypatch.setattr(adapter, 'collect', lambda name, *args: copy.deepcopy(next(r for r in rows if r['name'] == name)))
    single = dict(artifacts=[{'path': 'primary_model.cif', 'role': 'final', 'group': '主要结果'}],
                  final_model='primary_model.cif', final_phases=None, metrics=[{'model': '06_finish/model.cif'}],
                  primary_complete=True, build_model=True, config={}, preview=[], residues=0, chains=0)
    monkeypatch.setattr(workflow, 'collect_results', lambda root: copy.deepcopy(single))
    return pipe, report


def test_selected_model_and_downloads_stay_in_gui_result_root(tmp_path, monkeypatch):
    _, report = make_completed(tmp_path, monkeypatch)
    result = adapter.collect_results(tmp_path)
    assert result['pipeline_complete']
    assert result['portfolio']['recommended'] == report['recommended'] == 'pplus_w14'
    assert result['final_model'] == '00_portfolio/pplus_w14/primary_model.cif'
    assert result['metrics'][0]['model'] == '00_portfolio/pplus_w14/06_finish/model.cif'
    assert all(not Path(row['path']).is_absolute() for row in result['artifacts'])
    assert len([r for r in result['artifacts'] if r['group'] == '各流程候选']) == 2


def test_ambiguous_candidates_do_not_silently_export_one_model(tmp_path, monkeypatch):
    make_completed(tmp_path, monkeypatch, ambiguous=True)
    result = adapter.collect_results(tmp_path)
    assert result['pipeline_complete'] and not result['primary_complete']
    assert result['final_model'] is None
    assert result['portfolio']['recommended'] is None
    assert len([row for row in result['portfolio']['candidates'] if row['model']]) == 2


def test_tampered_recommendation_is_rejected(tmp_path, monkeypatch):
    pipe, report = make_completed(tmp_path, monkeypatch)
    report['recommended'] = 'pplus_w8'
    (pipe / 'PORTFOLIO_REPORT.json').write_text(json.dumps(report))
    with pytest.raises(ValueError, match='排序规则'):
        adapter.collect_results(tmp_path)


def test_identical_frontier_files_offer_explicit_representative_without_changing_ranking(tmp_path, monkeypatch):
    _, report = make_completed(tmp_path, monkeypatch, identical=True)
    result = adapter.collect_results(tmp_path)
    assert report['recommended'] is None
    assert result['portfolio']['recommended'] is None
    assert result['portfolio']['representative'] == 'pplus_w14'
    assert set(result['portfolio']['equivalent_candidates']) == {'pplus_w14', 'pplus_w8'}
    assert result['primary_complete']
    assert result['final_model'] == '00_portfolio/pplus_w14/primary_model.cif'


def test_representative_choice_does_not_follow_free_score(tmp_path, monkeypatch):
    pipe, report = make_completed(tmp_path, monkeypatch, identical=True)
    for row in report['candidates']:
        row['primary_metric']['r_free'] = .9 if row['name'] == 'pplus_w14' else .01
        path = pipe / row['name'] / 'pipeline_report.json'
        path.write_text(json.dumps({'status': 'completed', 'primary_metric': row['primary_metric']}))
    (pipe / 'PORTFOLIO_REPORT.json').write_text(json.dumps(report))
    result = adapter.collect_results(tmp_path)
    assert result['portfolio']['representative'] == 'pplus_w14'


def test_waiting_gui_job_has_complete_empty_result_schema(tmp_path):
    (tmp_path / 'gui_job.json').write_text(json.dumps({'config': {'explore': True}}))
    result = workflow.collect_results(tmp_path)
    assert result['portfolio'] and len(result['stages']) == 4
    assert result['final_model'] is None and not result['pipeline_complete']
    assert result['metrics'] == [] and result['artifacts'] == []
