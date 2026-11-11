import copy
import json
from pathlib import Path
import threading
import time

import gemmi
import pytest

from phasebuild import portfolio


def candidate(name, coverage=.9, work=.25, duplicate=0., unanchored=0., angle=1.1):
    return {'name': name, 'status': 'completed',
            'primary_metric': {'r_work': work, 'r_free': .3, 'bond_rmsd_angstrom': .004,
                               'angle_rmsd_degrees': angle},
            'sequence_audit': {'coverage': coverage, 'duplicate_fraction': duplicate,
                               'unanchored_fraction': unanchored}}


def test_free_scores_do_not_change_ranking_or_recommendation():
    rows = [candidate('a'), candidate('b', work=.28), candidate('c', coverage=.85, work=.24)]
    first = portfolio.rank_candidates(rows)
    altered = copy.deepcopy(rows)
    for i, row in enumerate(altered):
        row['primary_metric']['r_free'] = [999, -1, float('nan')][i]
    second = portfolio.rank_candidates(altered)
    for key in ['recommended', 'review_candidates', 'pareto_candidates', 'reason']:
        assert first[key] == second[key]
    assert first['recommended'] is None
    assert set(first['pareto_candidates']) == {'a', 'c'}
    assert rows[0].get('eligible') is None


def test_geometry_gate_rejects_lower_work_score_and_keeps_all_candidates():
    result = portfolio.rank_candidates([candidate('good'), candidate('distorted', work=.1, angle=4.1),
                                        {'name': 'failed', 'status': 'failed'}])
    assert result['recommended'] == 'good'
    assert result['review_candidates'] == ['good']
    assert len(result['candidates']) == 3


def test_no_eligible_candidate_means_no_recommendation():
    result = portfolio.rank_candidates([candidate('bad', coverage=0)])
    assert result['recommended'] is None
    assert result['review_candidates'] == []


def model_file(tmp_path, chains):
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model(1))
    for name, sequence in chains:
        chain = gemmi.Chain(name)
        for i, letter in enumerate(sequence):
            residue = gemmi.Residue()
            residue.name = gemmi.expand_one_letter_sequence(letter, gemmi.ResidueKind.AA)[0]
            residue.seqid = gemmi.SeqId(i + 1, ' ')
            atom = gemmi.Atom()
            atom.name = 'CA'
            atom.element = gemmi.Element('C')
            atom.pos = gemmi.Position(i * 3.8, 0, 0)
            residue.add_atom(atom)
            chain.add_residue(residue)
        structure[0].add_chain(chain)
    path = tmp_path / 'model.cif'
    structure.make_mmcif_document().write_file(str(path))
    return path


def test_duplicate_fragments_do_not_inflate_coverage(tmp_path):
    seq = 'ACDEFGHIKLMNPQRSTVWY'
    fasta = tmp_path / 'seq.fa'
    fasta.write_text('>test\n' + seq + '\n')
    model = model_file(tmp_path, [('A', seq), ('B', seq)])
    one = portfolio.sequence_audit(model, fasta, 1)
    assert one['coverage'] == 1 and one['duplicate_positions'] == len(seq)
    two = portfolio.sequence_audit(model, fasta, 2)
    assert two['coverage'] == 1 and two['duplicate_positions'] == 0


def test_repeated_sequence_and_short_fragments_are_not_false_coverage(tmp_path):
    fasta = tmp_path / 'seq.fa'
    fasta.write_text('>test\n' + 'A' * 30 + '\n')
    model = model_file(tmp_path, [('A', 'A' * 20)])
    result = portfolio.sequence_audit(model, fasta, 1)
    assert result['coverage'] == 0 and result['unanchored_residues'] == 20


def setup_config(tmp_path, monkeypatch):
    from phasebuild import sad_pipeline
    cfg = dict(sad_pipeline.DEFAULTS, method='pplus_bundled')
    for key in ['sad', 'native', 'sequence', 'heavy']:
        path = tmp_path / (key + '.dat')
        path.write_text(key)
        cfg[key] = str(path)
    monkeypatch.setattr(sad_pipeline, 'validate', lambda x, **kwargs: dict(x))
    return cfg


def test_plan_is_read_only_and_keeps_common_flags_and_budget(tmp_path, monkeypatch):
    cfg = setup_config(tmp_path, monkeypatch)
    out = tmp_path / 'outputs'
    result = portfolio.plan(cfg, out)
    assert not out.exists()
    assert len(result['branches']) == 4
    for branch in result['branches']:
        for key in ['sad', 'native', 'free_value', 'seed', 'rounds', 'cycles']:
            assert branch['config'][key] == cfg[key]
    with pytest.raises(ValueError, match='distinct'):
        portfolio.plan(cfg, out, ['pplus_w14', 'pplus_w14'])
    with pytest.raises(ValueError, match='Unknown'):
        portfolio.plan({**cfg, 'typo': True}, out)
    out.mkdir()
    with pytest.raises(ValueError, match='exists'):
        portfolio.plan(cfg, out)


def test_executor_bounds_concurrency_isolates_failure_and_preserves_results(tmp_path, monkeypatch):
    cfg = setup_config(tmp_path, monkeypatch)
    active = 0
    peak = 0
    lock = threading.Lock()

    def runner(path, directory, stop):
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
        try:
            time.sleep(.03)
            if path.stem == 'stable_w14':
                raise RuntimeError('engine failure')
        finally:
            with lock:
                active -= 1

    monkeypatch.setattr(portfolio, 'collect', lambda name, directory, cfg: candidate(name))
    out = tmp_path / 'portfolio'
    report = portfolio.run(cfg, out, workers=2, runner=runner)
    assert peak == 2
    assert report['status'] == 'partial' and report['inputs_unchanged']
    assert len(report['candidates']) == 4
    assert [r for r in report['candidates'] if r['name'] == 'stable_w14'][0]['status'] == 'failed'
    assert json.loads((out / 'PORTFOLIO_REPORT.json').read_text()) == report
    assert cfg['method'] == 'pplus_bundled'
    with pytest.raises(ValueError, match='exists'):
        portfolio.run(cfg, out, runner=runner)


def test_changed_inputs_suppress_recommendation(tmp_path, monkeypatch):
    cfg = setup_config(tmp_path, monkeypatch)
    def runner(path, directory, stop):
        Path(cfg['native']).write_text('changed by external process')
    monkeypatch.setattr(portfolio, 'collect', lambda name, directory, cfg: candidate(name))
    result = portfolio.run(cfg, tmp_path / 'out', names=['pplus_w14'], runner=runner)
    assert not result['inputs_unchanged']
    assert result['recommended'] is None and result['review_candidates'] == []
    assert result['status'] == 'failed'


def test_final_model_review_gates_are_distinct_from_strict_weight_pilot_gates():
    row = candidate('ordinary_refinement', angle=2.386)
    row['primary_metric']['bond_rmsd_angstrom'] = .0119
    assert portfolio.rank_candidates([row])['recommended'] == 'ordinary_refinement'


def test_cancellation_does_not_start_queued_branch(tmp_path, monkeypatch):
    cfg = setup_config(tmp_path, monkeypatch)
    called = []
    def runner(path, directory, stop):
        called.append(path.stem)
        stop.set()
        raise InterruptedError('cancelled')
    result = portfolio.run(cfg, tmp_path / 'out', workers=1, runner=runner)
    assert len(called) == 1
    assert result['status'] == 'cancelled'
    assert all(row['status'] == 'cancelled' for row in result['candidates'])
