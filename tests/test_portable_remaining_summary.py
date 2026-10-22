# -*- coding: utf-8 -*-
"""Accounting must not hide unfinished/invalid cases or fabricate successes."""
import importlib.util
import json
from pathlib import Path

import pytest


SPEC = importlib.util.spec_from_file_location(
    'portable_remaining_summary', Path(__file__).resolve().parents[1] / 'packaging/summarize_portable_remaining.py')
summary_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summary_module)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding='utf-8')


def completed_state(case, success=True):
    return {'id': case['id'], 'status': 'completed', 'completed_rounds': 5,
            'primary_metric': {'r_work': .25, 'r_free': .30},
            'endpoint_audit': {'passed': True},
            'evaluation': {'built': 100, 'reference_protein_residues': 100,
                           'strict_completeness': .95, 'strict_precision': .96,
                           'joint_completeness': .94, 'protein_success': success,
                           'full_structure_success': success and case['sample_class'] == 'protein'}}


@pytest.fixture
def corpus(tmp_path):
    root = tmp_path / 'run'
    cases = [{'id': 'p%03d' % index, 'sample_class': 'protein', 'supported': True, 'has_free': True}
             for index in range(79)]
    cases += [{'id': 'caoy_6e9c', 'sample_class': 'protein', 'supported': True, 'has_free': False}]
    cases += [{'id': 'c%03d' % index, 'sample_class': 'protein_nucleic_complex', 'supported': True, 'has_free': True}
              for index in range(6)]
    cases += [{'id': 'n%03d' % index, 'sample_class': 'nucleic_or_small_nonstandard', 'supported': False}
              for index in range(10)]
    development = [{'id': cid, 'sample_class': 'protein', 'supported': True, 'has_free': True}
                   for cid in sorted(summary_module.DEVELOPMENT_IDS)]
    original = tmp_path / 'original.json'
    write(original, {'cases': cases + development})
    four = tmp_path / 'four'
    write(four / 'comparison.json', {
        'scientific_consistency_passed': True,
        'cases': [{'id': case['id'], 'audit_passed': True, 'execution_status': 'completed',
                   'portable': {'r_work': .25, 'r_free': .3, 'built_residues': 100,
                                'reference_residues': 100, 'n_ca_c_completeness': .95,
                                'n_ca_c_precision': .96, 'joint_sequence_backbone_coverage': .94,
                                'strict_success': True}} for case in development]})
    write(root / 'manifest.json', {'cases': cases, 'excluded_development_cases': sorted(summary_module.DEVELOPMENT_IDS),
                                  'original_corpus_manifest': str(original), 'previous_four_run': str(four),
                                  'buccaneer_baseline': str(tmp_path / 'empty_bucc'),
                                  'native_baseline': str(tmp_path / 'empty_native'),
                                  'rules': {'thresholds': summary_module.CRITERIA}})
    for case in cases:
        if not case['supported']:
            state = {'id': case['id'], 'status': 'unsupported_nonprotein'}
        elif not case['has_free']:
            state = completed_state(case)
            state.update(status='completed_no_free', completed_rounds=1, primary_metric=None)
            state['evaluation'].update(protein_success=None, full_structure_success=None)
        else:
            state = completed_state(case)
        write(root / 'states' / (case['id'] + '.json'), state)
    return root, cases


def test_registered_denominators_and_development_cohort_remain_distinct(corpus):
    root, cases = corpus
    result = summary_module.summarize(root, final=True)
    assert result['complete'] and len(result['cases']) == 100
    assert result['remaining_96']['pure_protein']['registered'] == 80
    assert result['remaining_96']['pure_protein']['evaluated_fixed_endpoint'] == 79
    assert result['remaining_96']['pure_protein']['strict_success'] == 79
    assert result['combined_100']['pure_protein']['registered'] == 84
    assert result['combined_100']['pure_protein']['strict_success'] == 83
    assert result['remaining_96']['complex_protein_component']['strict_success'] == 6
    assert all(row['full_structure_success'] is False for row in result['cases']
               if row['sample_class'] == 'protein_nucleic_complex')
    assert sum(row['cohort'] == 'development_4' for row in result['cases']) == 4


def test_missing_case_is_retained_and_final_refused(corpus):
    root, cases = corpus
    (root / 'states/p000.json').unlink()
    result = summary_module.summarize(root)
    assert not result['complete'] and len(result['cases']) == 100
    assert result['remaining_96']['pure_protein']['registered'] == 80
    assert result['remaining_96']['classification_counts']['incomplete'] == 1
    with pytest.raises(ValueError, match='Final result is not ready'):
        summary_module.summarize(root, final=True)


def test_unscored_and_unaudited_completions_are_never_successes(corpus):
    root, cases = corpus
    state = completed_state(cases[0])
    state.pop('endpoint_audit')
    write(root / 'states/p000.json', state)
    other = completed_state(cases[1])
    other.pop('evaluation')
    write(root / 'states/p001.json', other)
    result = summary_module.summarize(root)
    rows = {row['id']: row for row in result['cases']}
    assert rows['p000']['classification'] == 'unaudited_completed'
    assert rows['p001']['classification'] == 'unscored_completed'
    assert rows['p000']['success'] is rows['p001']['success'] is None
    assert not result['complete']


def test_fabricated_no_free_rfactor_is_an_integrity_error(corpus):
    root, _ = corpus
    path = root / 'states/caoy_6e9c.json'
    state = json.loads(path.read_text(encoding='utf-8'))
    state['primary_metric'] = {'r_free': .3}
    write(path, state)
    result = summary_module.summarize(root)
    assert not result['complete']
    assert any('No-free record' in error for error in result['integrity_errors'])
    row = next(row for row in result['cases'] if row['id'] == 'caoy_6e9c')
    assert row['success'] is None


def test_low_rfree_with_poor_precision_is_a_scientific_failure():
    case = {'id': 'weak', 'sample_class': 'protein', 'supported': True, 'has_free': True}
    state = completed_state(case, success=False)
    state['evaluation']['strict_precision'] = .89
    row = summary_module.row_from_state(case, state, summary_module.CRITERIA)
    assert row['classification'] == 'scientific_not_solved'
    assert row['success'] is False


def test_stale_saved_success_cannot_override_frozen_criteria():
    case = {'id': 'weak', 'sample_class': 'protein', 'supported': True, 'has_free': True}
    state = completed_state(case)
    state['evaluation']['strict_precision'] = .89
    row = summary_module.row_from_state(case, state, summary_module.CRITERIA)
    assert row['classification'] == 'evaluation_failure'
    assert row['success'] is None and row['integrity_errors']


def test_reference_gate_requires_terminal_endpoint_and_independent_audit(tmp_path):
    case = {'id': 'p', 'sample_class': 'protein', 'supported': True, 'has_free': True}
    state = {'id': 'p', 'status': 'running'}
    with pytest.raises(ValueError, match='completed construction'):
        summary_module.evaluate_completed(case, state, tmp_path)
    state.update(status='completed', primary_model='/nonexistent/model.cif', completed_rounds=5)
    with pytest.raises(ValueError, match='independent scientific audit'):
        summary_module.evaluate_completed(case, state, tmp_path)


def test_endpoint_gate_rejects_cherry_picked_round_before_reference(tmp_path):
    case = {'id': 'p', 'sample_class': 'protein', 'supported': True, 'has_free': True}
    pipeline = tmp_path / 'pipeline'
    pipeline.mkdir()
    model = pipeline / 'primary_model.cif'
    model.write_text('not parsed before endpoint gate', encoding='utf-8')
    final = pipeline / 'round_05_refine/model.cif'
    final.parent.mkdir()
    final.write_bytes(model.read_bytes())
    stages = [{'round': n, 'refmac_cycles': 10, 'buccaneer_cycles': 3 if n == 1 else 2,
               'statistics': {'r_work': .3, 'r_free': .35 + n / 100}}
              for n in range(1, 6)]
    picked = stages[0]['statistics']
    report = {'status': 'completed', 'source_changed_during_run': False,
              'parameters': {'rounds': 5, 'refine_cycles': 10}, 'stages': stages,
              'primary_model': summary_module.fingerprint(model), 'primary_metric': picked}
    write(pipeline / 'hybrid_report.json', report)
    state = {'id': 'p', 'status': 'completed', 'completed_rounds': 5,
             'endpoint_audit': {'passed': True}, 'primary_model': str(model), 'primary_metric': picked}
    with pytest.raises(ValueError, match='fixed final round'):
        summary_module.evaluate_completed(case, state, tmp_path)
