# -*- coding: utf-8 -*-
"""Post-build evaluation and live/final accounting for the portable remaining 96.

The summarizer reads saved metrics only.  The separate evaluate_completed()
entry point may open a target reference only after the declared construction
endpoint has completed.  It returns the updated state; the batch driver owns
the atomic state write, so an evaluation cannot overwrite live worker state.
No program is launched and no construction input or frozen result is changed.
"""
import argparse
import collections
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import sys
import time

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2] if (_HERE.parent / 'engine').is_dir() else _HERE.parents[1]
BUCCANEER_BASELINE = ROOT / 'runs/buccaneer_sad100_20261019'
NATIVE_BASELINE = ROOT / 'runs/sad100_20261019'
DEVELOPMENT_IDS = {'5hhk', '7bll', '7tem', '7twc'}
COMPLETED = {'completed', 'completed_no_free'}
TERMINAL = COMPLETED | {'unsupported_nonprotein', 'seed_failed', 'failed', 'timeout', 'cancelled'}
CRITERIA = {'r_free_max': .35, 'backbone_completeness_min': .8,
            'backbone_precision_min': .9,
            'joint_sequence_backbone_completeness_min': .8,
            'matching_cutoff_angstrom': 1.}
LIMITATIONS = ('四个已知成功案例是开发集；剩余96例是预先固定流程的后续测试，'
               '合并100例不能称为100例独立盲测。重原子位置及上游组成/溶剂元数据已给定；'
               '上游OASIS定相使用过自由反射，Rfree是模型诊断，不是无泄漏的端到端盲测统计。'
               '该方法包含外部Buccaneer和REFMAC组件，不能归为原生PhaseBuild算法成果。')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(text, encoding='utf-8')
    temporary.replace(path)


def save(path, value):
    atomic_text(path, json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def numeric(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def thresholds(manifest):
    result = manifest.get('rules', {}).get('thresholds', CRITERIA)
    # This batch follows the already fixed success rule.  An accidental
    # manifest edit must never silently redefine the yield.
    for key, value in CRITERIA.items():
        if result.get(key) != value:
            raise ValueError('Frozen criterion differs: ' + key)
    return dict(CRITERIA)


def strict_success(metric, evaluation, rules=CRITERIA):
    fields = [(metric.get('r_free'), 'r_free'),
              (evaluation.get('strict_completeness'), 'strict_completeness'),
              (evaluation.get('strict_precision'), 'strict_precision'),
              (evaluation.get('joint_completeness'), 'joint_completeness')]
    if any(not numeric(value) for value, _ in fields):
        raise ValueError('Missing/nonfinite fixed-endpoint metric: ' +
                         ', '.join(name for value, name in fields if not numeric(value)))
    if metric['r_free'] < 0 or any(not 0 <= value <= 1 for value, _ in fields[1:]):
        raise ValueError('Invalid structural metric or negative Rfree')
    return bool(metric['r_free'] <= rules['r_free_max'] and
                evaluation['strict_completeness'] >= rules['backbone_completeness_min'] and
                evaluation['strict_precision'] >= rules['backbone_precision_min'] and
                evaluation['joint_completeness'] >= rules['joint_sequence_backbone_completeness_min'])


def path_from(value, run):
    if isinstance(value, dict):
        value = value['path']
    path = Path(value)
    return path.resolve() if path.is_absolute() else (Path(run) / path).resolve()


def completed_endpoint(case, state, run):
    """Check endpoint bookkeeping before any reference can be opened.

    Detailed per-reflection/native-HL and engine audits may be attached as
    state.endpoint_audit by the batch auditor.  These checks reject incomplete
    or best-Rfree substitutions even when a worker accidentally says completed.
    """
    if state.get('status') not in COMPLETED:
        raise ValueError('Evaluation only accepts a completed construction endpoint')
    if state.get('id') != case['id']:
        raise ValueError('Evaluation case/state identity differs')
    if (state.get('endpoint_audit') or {}).get('passed') is not True:
        raise ValueError('Completed endpoint must pass the independent scientific audit before reference evaluation')
    model = path_from(state['primary_model'], run)
    if not model.is_file():
        raise ValueError('Completed primary model is missing')
    if state['status'] == 'completed_no_free':
        if case.get('has_free') is not False or state.get('primary_metric'):
            raise ValueError('No-free endpoint cannot contain an invented R factor')
        if state.get('completed_rounds') != 1:
            raise ValueError('No-free case must contain exactly the one initial construction')
        return model, {'endpoint': 'one_initial_build_without_existing_free_set',
                       'refmac_cycles': 0, 'has_validated_rfree': False}
    if case.get('has_free') is not True or state.get('completed_rounds') != 5:
        raise ValueError('Completed endpoint requires original free set and five rounds')
    report_path = state.get('pipeline_report') or model.parent / 'hybrid_report.json'
    report_path = path_from(report_path, run)
    report = read(report_path)
    stages = report.get('stages', [])
    if report.get('status') != 'completed' or len(stages) != 5:
        raise ValueError('Pipeline has not completed all five frozen stages')
    if report.get('source_changed_during_run') is not False:
        raise ValueError('Pipeline source changed or source-integrity record is missing')
    for number, stage in enumerate(stages, 1):
        if (stage.get('round') != number or stage.get('refmac_cycles') != 10 or
                stage.get('buccaneer_cycles') != (3 if number == 1 else 2)):
            raise ValueError('Construction/refinement stage schedule differs')
    parameters = report.get('parameters', {})
    if parameters.get('rounds') != 5 or parameters.get('refine_cycles') != 10:
        raise ValueError('Pipeline scientific budget differs')
    if report.get('primary_metric') != stages[-1].get('statistics') or state.get('primary_metric') != report['primary_metric']:
        raise ValueError('Primary metric is not the fixed final round')
    actual = fingerprint(model)
    if report.get('primary_model') != actual:
        raise ValueError('Primary model differs from completed pipeline fingerprint')
    final_model = report_path.parent / 'round_05_refine/model.cif'
    if fingerprint(final_model)['sha256'] != actual['sha256']:
        raise ValueError('Primary model is not the fixed final model')
    return model, {'endpoint': 'fixed_final_round_5', 'refmac_cycles': 50,
                   'has_validated_rfree': True, 'pipeline_report': fingerprint(report_path)}


def evaluate_completed(case, state, run):
    """Return state with posterior structural evaluation; caller saves the state.

    Import the frozen bundled evaluator before calling, or put bundle/app/src on
    PYTHONPATH.  The target reference is never supplied to any builder here.
    Known caoy b/2 gauge correction is applied only to a new evaluation copy.
    """
    model, endpoint = completed_endpoint(case, state, run)
    run = Path(run).resolve()
    manifest = read(run / 'manifest.json')
    rules = thresholds(manifest)
    source_root = Path(manifest.get('reference_evaluation_root', BUCCANEER_BASELINE / 'evaluation')).resolve()
    reference = source_root / case['id'] / 'protein_reference.cif'
    if not reference.is_file() or reference.name != 'protein_reference.cif':
        raise ValueError('Residue-aware original reference is missing')
    # No native/target model can be used as a fallback reference.
    reference_digest = fingerprint(reference)
    declared = case.get('evaluation_reference')
    if declared and declared != reference_digest:
        raise ValueError('Declared residue-aware reference fingerprint changed')
    output = run / 'evaluation' / case['id']
    output.mkdir(parents=True, exist_ok=True)
    shift = [0., 0., 0.]
    if case['id'] == 'caoy_6e9c':
        previous = read(source_root / case['id'] / 'accuracy.json')
        shift = previous.get('evaluation_only_origin_shift')
        if shift != [0, .5, 0]:
            raise ValueError('Documented caoy evaluation-only origin correction differs')
    evaluation_model = model
    if any(shift):
        import gemmi
        structure = gemmi.read_structure(str(model))
        translation = structure.cell.orthogonalize(gemmi.Fractional(*shift))
        for chain in structure[0]:
            for residue in chain:
                for atom in residue:
                    atom.pos = gemmi.Position(atom.pos.x + translation.x,
                                              atom.pos.y + translation.y,
                                              atom.pos.z + translation.z)
        evaluation_model = output / 'evaluation_only_shift.cif'
        structure.make_mmcif_document().write_file(str(evaluation_model))
    from phasebuild.evaluate import evaluate
    accuracy = evaluate(evaluation_model, reference, cutoff=rules['matching_cutoff_angstrom'])
    accuracy.update(evaluation_only_origin_shift=shift, original_primary_model=fingerprint(model),
                    reference_usage='Posterior evaluation only; no reference used for construction or seed selection',
                    endpoint=endpoint)
    n = accuracy['reference_residues']
    if n <= 0:
        raise ValueError('Residue-aware reference has no protein residues')
    result = {'evaluator_revision': 'residue-aware-cached-original-reference-v3-portable',
              'built': accuracy['built_residues'], 'reference_protein_residues': n,
              'strict_backbone_correct': accuracy['n_ca_c_correct_residues'],
              'strict_completeness': accuracy['n_ca_c_completeness'],
              'strict_precision': accuracy['n_ca_c_precision'],
              'joint_sequence_backbone_correct': accuracy['sequence_evaluation']['correct_type_and_backbone'],
              'joint_completeness': accuracy['sequence_evaluation']['correct_type_and_backbone'] / n,
              'reference': reference_digest, 'original_primary_model': fingerprint(model),
              'evaluation_only_origin_shift': shift,
              'endpoint': endpoint}
    passed = strict_success(state.get('primary_metric') or {}, result, rules) if state['status'] == 'completed' else None
    result.update(protein_success=passed,
                  full_structure_success=bool(passed and case['sample_class'] == 'protein') if passed is not None else None)
    save(output / 'accuracy.json', accuracy)
    result['accuracy'] = fingerprint(output / 'accuracy.json')
    updated = dict(state, evaluation=result)
    updated.pop('evaluation_error', None)
    return updated


def row_from_state(case, state, rules):
    status = state.get('status', 'queued')
    metric = state.get('primary_metric') or {}
    evaluation = state.get('evaluation') or {}
    row = {'id': case['id'], 'sample_class': case['sample_class'], 'cohort': 'remaining_96',
           'status': status, 'current': state.get('current', ''),
           'completed_rounds': state.get('completed_rounds', 0),
           'r_work': metric.get('r_work'), 'r_free': metric.get('r_free'),
           'built': evaluation.get('built'), 'reference_protein_residues': evaluation.get('reference_protein_residues'),
           'strict_completeness': evaluation.get('strict_completeness'),
           'strict_precision': evaluation.get('strict_precision'),
           'joint_completeness': evaluation.get('joint_completeness'),
           'success': None, 'full_structure_success': None,
           'error': state.get('error', ''), 'evaluation_error': state.get('evaluation_error', ''),
           'endpoint_audit_passed': (state.get('endpoint_audit') or {}).get('passed'),
           'integrity_errors': []}
    if not case.get('supported', True) or status == 'unsupported_nonprotein':
        row['classification'] = 'unsupported_nonprotein'
        if status != 'unsupported_nonprotein':
            row['integrity_errors'].append('Unsupported nonprotein lacks explicit terminal record')
    elif status == 'completed_no_free':
        row['classification'] = 'no_existing_free_set'
        if case.get('has_free') is not False or metric:
            row['integrity_errors'].append('No-free record has a free-set/R-factor contradiction')
        if not evaluation or row['evaluation_error']:
            row['evaluation_pending'] = True
    elif status == 'completed':
        if row['evaluation_error']:
            row['classification'] = 'evaluation_failure'
        elif not evaluation:
            row['classification'] = 'unscored_completed'
        elif row['endpoint_audit_passed'] is False:
            row['classification'] = 'audit_failure'
        elif row['endpoint_audit_passed'] is not True:
            row['classification'] = 'unaudited_completed'
        else:
            try:
                if state.get('completed_rounds') != 5:
                    raise ValueError('Completed state has not reached five rounds')
                passed = strict_success(metric, evaluation, rules)
                if evaluation.get('protein_success') is not passed:
                    raise ValueError('Stored success does not agree with frozen thresholds')
                full = bool(passed and case['sample_class'] == 'protein')
                if evaluation.get('full_structure_success') is not full:
                    raise ValueError('Full-structure success improperly includes a complex')
                row.update(success=passed, full_structure_success=full,
                           classification='strict_success' if passed else 'scientific_not_solved')
            except ValueError as error:
                row['classification'] = 'evaluation_failure'
                row['integrity_errors'].append(str(error))
    elif status == 'seed_failed':
        row['classification'] = 'seed_failure'
    elif status in {'failed', 'cancelled'}:
        row['classification'] = 'execution_failure'
    elif status == 'timeout':
        row['classification'] = 'timeout'
    else:
        row['classification'] = 'incomplete'
        if status not in {'queued', 'pending', 'running', 'seed_running', 'evaluating'}:
            row['integrity_errors'].append('Unknown state status: ' + str(status))
    # Nonfinite R values cannot be serialized as invented JSON NaN statistics.
    for key in ['r_work', 'r_free', 'strict_completeness', 'strict_precision', 'joint_completeness']:
        if row[key] is not None and not numeric(row[key]):
            row['integrity_errors'].append('Nonfinite/non-numeric ' + key)
            row[key] = None
    return row


def previous_four_rows(manifest, original_cases):
    value = manifest.get('previous_four_run', ROOT / 'runs/portable_validation_20261019_retry3')
    root = Path(value['path'] if isinstance(value, dict) else value).resolve()
    comparison_path = root / 'comparison.json'
    comparison = read(comparison_path)
    if (comparison.get('scientific_consistency_passed') is not True or
            {item['id'] for item in comparison['cases']} != DEVELOPMENT_IDS):
        raise ValueError('Previous four-case audit is incomplete or identities differ')
    rows = []
    for item in comparison['cases']:
        result = item['portable']
        if item.get('audit_passed') is not True or item.get('execution_status') != 'completed':
            raise ValueError('Previous development case was not scientifically audited')
        case = original_cases[item['id']]
        evaluation = {'built': result['built_residues'], 'reference_protein_residues': result['reference_residues'],
                      'strict_completeness': result['n_ca_c_completeness'],
                      'strict_precision': result['n_ca_c_precision'],
                      'joint_completeness': result['joint_sequence_backbone_coverage'],
                      'protein_success': result['strict_success'],
                      'full_structure_success': result['strict_success']}
        state = {'status': 'completed', 'completed_rounds': 5,
                 'primary_metric': {'r_work': result['r_work'], 'r_free': result['r_free']},
                 'evaluation': evaluation, 'endpoint_audit': {'passed': True}}
        row = row_from_state(case, state, CRITERIA)
        row.update(cohort='development_4', source_comparison=fingerprint(comparison_path))
        rows.append(row)
    return rows


def baseline_rows(root):
    """Read saved baseline states only; never read a target structure here."""
    return {path.stem: read(path) for path in (Path(root) / 'states').glob('*.json')}


def baseline_metrics(state, prefix):
    metric, evaluation = state.get('primary_metric') or {}, state.get('evaluation') or {}
    success = evaluation.get('protein_success', evaluation.get('protein_success_by_frozen_rule'))
    if state.get('status') != 'completed':
        success = None
    return {prefix + '_status': state.get('status'), prefix + '_r_work': metric.get('r_work'),
            prefix + '_r_free': metric.get('r_free'),
            prefix + '_strict_completeness': evaluation.get('strict_completeness'),
            prefix + '_strict_precision': evaluation.get('strict_precision'),
            prefix + '_joint_completeness': evaluation.get('joint_completeness'),
            prefix + '_success': success}


def statistics_for(rows):
    evaluated = [row for row in rows if row['success'] is not None]
    valid_r = [row['r_free'] for row in evaluated if numeric(row['r_free'])]
    passed = sum(row['success'] is True for row in rows)
    return {'registered': len(rows), 'strict_success': passed,
            'yield_all_registered': passed / len(rows) if rows else None,
            'evaluated_fixed_endpoint': len(evaluated),
            'success_fraction_evaluated': passed / len(evaluated) if evaluated else None,
            'classification_counts': dict(collections.Counter(row['classification'] for row in rows)),
            'r_free': {'min': min(valid_r), 'median': statistics.median(valid_r), 'max': max(valid_r)} if valid_r else None}


def paired_statistics(rows, prefix):
    paired = [row for row in rows if row['success'] is not None and
              row.get(prefix + '_success') is not None and numeric(row.get(prefix + '_r_free'))]
    deltas = [row['r_free'] - row[prefix + '_r_free'] for row in paired]
    return {'n': len(paired), 'portable_success': sum(row['success'] is True for row in paired),
            'baseline_success': sum(row[prefix + '_success'] is True for row in paired),
            'r_free_delta_portable_minus_baseline': {
                'mean': statistics.mean(deltas), 'median': statistics.median(deltas),
                'portable_lower': sum(value < 0 for value in deltas),
                'portable_higher': sum(value > 0 for value in deltas),
                'ties': sum(value == 0 for value in deltas)} if deltas else None,
            'case_ids': [row['id'] for row in paired]}


def percent(value):
    return '—' if value is None else '{:.2%}'.format(value)


def render_report(summary):
    new = summary['remaining_96']
    total = summary['combined_100']
    pure = new['pure_protein']
    text = ['# PhaseBuild便携流程：剩余96例SAD测试', '',
            '**当前状态：{}。**'.format('全部完成' if summary['complete'] else '运行中，以下不是最终成功数'), '',
            '当前流程包含PhaseBuild种子、Buccaneer构建及原始实验HL约束的REFMAC精修。'
            '主终点固定第5轮累计50个REFMAC周期，不挑最低Rfree或参考最相似轮次。', '',
            '| 范围 | 注册数 | 完整固定终点评价 | 严格达标 |', '|---|---:|---:|---:|',
            '| 剩余纯蛋白 | {} | {} | {} |'.format(pure['registered'], pure['evaluated_fixed_endpoint'], pure['strict_success']),
            '| 剩余复合物（只评蛋白） | 6 | {} | {} |'.format(new['complex_protein_component']['evaluated_fixed_endpoint'], new['complex_protein_component']['strict_success']),
            '| 非蛋白、不支持 | 10 | — | — |',
            '| 已完成开发集纯蛋白 | 4 | 4 | {} |'.format(summary['development_4']['strict_success']),
            '| 合并纯蛋白（开发＋后续） | 84 | {} | {} |'.format(total['pure_protein']['evaluated_fixed_endpoint'], total['pure_protein']['strict_success']), '',
            '剩余纯蛋白的注册分母固定80；合并纯蛋白固定84。复合物不能记为完整结构解析成功。'
            '缺原自由集的caoy_6e9c仅单次构建，保留结构指标，不生成自由集或Rfree。', '',
            '预设判据：Rfree≤35%，N/CA/C在CA选择的同一晶体对称像下均<1 Å，'
            '严格主链覆盖率≥80%、精确率≥90%、序列与严格主链联合覆盖率≥80%。', '',
            '## 执行与评价分类', '', '| 分类 | 剩余96例数量 |', '|---|---:|']
    for label, count in sorted(new['classification_counts'].items()):
        text.append('| {} | {} |'.format(label, count))
    text += ['', '未完成、种子失败、程序失败、超时、评价失败和未评分保留在注册分母，'
             '不会被当作已解或从分母删除。失败与尚未完成并非同一结论。', '',
             '## 限制', '', LIMITATIONS, '',
             '当前只声明无需安装整套CCP4，包内仍有Buccaneer/REFMAC等组件；许可适用范围与本地评估包一致。', '',
             '## 全部100例逐例结果', '',
             '| 案例 | 分组 | 类型 | 状态/分类 | Rwork | Rfree | 严格覆盖 | 严格精确率 | 序列联合覆盖 | 达标 | 旧Buccaneer Rfree |',
             '|---|---|---|---|---:|---:|---:|---:|---:|---|---:|']
    for row in summary['cases']:
        text.append('| {} | {} | {} | {} / {} | {} | {} | {} | {} | {} | {} | {} |'.format(
            row['id'], row['cohort'], row['sample_class'], row['status'], row['classification'],
            percent(row['r_work']), percent(row['r_free']), percent(row['strict_completeness']),
            percent(row['strict_precision']), percent(row['joint_completeness']),
            '—' if row['success'] is None else '是' if row['success'] else '否', percent(row.get('buccaneer_r_free'))))
    if summary['integrity_errors']:
        text += ['', '## 待处理完整性问题', '']
        text += ['- ' + error for error in summary['integrity_errors']]
    text += ['', '机器可读结果：summary.json / comparison.json / comparison.csv；'
             '状态在states/，逐残基后验评价在evaluation/<id>/accuracy.json。', '']
    return '\n'.join(text)


def summarize(run, final=False):
    run = Path(run).resolve()
    manifest = read(run / 'manifest.json')
    rules = thresholds(manifest)
    original_manifest = Path(manifest.get('original_corpus_manifest', BUCCANEER_BASELINE / 'manifest.json'))
    original = read(original_manifest)
    original_cases = {case['id']: case for case in original['cases']}
    cases = manifest['cases']
    case_ids = [case['id'] for case in cases]
    excluded = set(manifest.get('excluded_development_cases', []))
    if (len(original_cases) != 100 or len(cases) != 96 or len(set(case_ids)) != 96 or
            excluded != DEVELOPMENT_IDS or set(case_ids) != set(original_cases) - excluded):
        raise ValueError('Remaining cohort must be exactly original100 minus the four declared development cases')
    if collections.Counter(case['sample_class'] for case in cases) != {
            'protein': 80, 'protein_nucleic_complex': 6, 'nucleic_or_small_nonstandard': 10}:
        raise ValueError('Remaining96 sample-class denominators differ from80/6/10')
    rows, errors = [], []
    known_paths = {case['id'] + '.json' for case in cases}
    extra = sorted(path.name for path in (run / 'states').glob('*.json') if path.name not in known_paths)
    if extra:
        errors.append('Unexpected state files: ' + ', '.join(extra))
    for case in cases:
        if case['sample_class'] != original_cases[case['id']]['sample_class']:
            raise ValueError('Original sample class changed: ' + case['id'])
        path = run / 'states' / (case['id'] + '.json')
        try:
            state = read(path) if path.exists() else {'status': 'queued'}
            if state.get('id', case['id']) != case['id']:
                raise ValueError('State identity does not match filename')
            row = row_from_state(case, state, rules)
        except (ValueError, OSError, KeyError, TypeError) as error:
            row = row_from_state(case, {'status': 'queued'}, rules)
            row['classification'] = 'state_read_failure'
            row['integrity_errors'].append(str(error))
        rows.append(row)
    four = previous_four_rows(manifest, original_cases)
    combined = sorted(rows + four, key=lambda row: row['id'])
    bucc = baseline_rows(manifest.get('buccaneer_baseline', BUCCANEER_BASELINE))
    native = baseline_rows(manifest.get('native_baseline', NATIVE_BASELINE))
    for row in combined:
        row.update(baseline_metrics(bucc.get(row['id'], {}), 'buccaneer'))
        row.update(baseline_metrics(native.get(row['id'], {}), 'native_phasebuild'))
        row['r_free_delta_vs_buccaneer'] = (row['r_free'] - row['buccaneer_r_free']
                                           if numeric(row['r_free']) and numeric(row['buccaneer_r_free']) else None)
        errors.extend(row['id'] + ': ' + error for error in row['integrity_errors'])
    terminal = all(row['status'] in TERMINAL for row in rows)
    unscored = [row['id'] for row in rows if row['classification'] in {'unscored_completed', 'unaudited_completed', 'evaluation_failure', 'audit_failure', 'state_read_failure'}
                or row.get('evaluation_pending')]
    pure = [row for row in rows if row['sample_class'] == 'protein']
    allpure = [row for row in combined if row['sample_class'] == 'protein']
    complex_rows = [row for row in rows if row['sample_class'] == 'protein_nucleic_complex']
    summary = {'updated': time.time(), 'run': str(run), 'status': 'complete' if terminal and not unscored and not errors else 'running_or_attention_required',
               'complete': terminal and not unscored and not errors, 'all_construction_terminal': terminal,
               'registered': 100, 'remaining_registered': 96, 'development_registered': 4,
               'criteria': rules, 'primary_endpoint': 'Fixed final round5,50REFMACcycles; no minimum-Rfree selection',
               'method': 'PhaseBuild seed plus integrated external Buccaneer and REFMAC with fixed original experimental HL',
               'remaining_96': {'status_counts': dict(collections.Counter(row['status'] for row in rows)),
                                'classification_counts': dict(collections.Counter(row['classification'] for row in rows)),
                                'pure_protein': statistics_for(pure),
                                'complex_protein_component': statistics_for(complex_rows),
                                'unsupported_nonprotein_registered': 10},
               'development_4': statistics_for(four),
               'combined_100': {'pure_protein': statistics_for(allpure),
                                'complex_protein_component': statistics_for(complex_rows),
                                'unsupported_nonprotein_registered': 10},
               'paired_remaining_pure_vs_buccaneer': paired_statistics(pure, 'buccaneer'),
               'paired_remaining_pure_vs_native': paired_statistics(pure, 'native_phasebuild'),
               'paired_combined_pure_vs_buccaneer': paired_statistics(allpure, 'buccaneer'),
               'integrity_errors': errors, 'completed_but_unscored_or_evaluation_failed': unscored,
               'successful_remaining_pure_cases': [row['id'] for row in pure if row['success'] is True],
               'limitations': LIMITATIONS, 'cases': combined}
    save(run / 'summary.json', summary)
    save(run / 'comparison.json', summary)
    csv_rows = [{key: value for key, value in row.items() if key not in {'integrity_errors', 'source_comparison'}} for row in combined]
    columns = list(dict.fromkeys(key for row in csv_rows for key in row))
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fieldnames=columns)
    writer.writeheader()
    writer.writerows(csv_rows)
    atomic_text(run / 'comparison.csv', buffer.getvalue())
    atomic_text(run / 'REPORT.md', render_report(summary))
    if final and not summary['complete']:
        raise ValueError('Final result is not ready: pending construction, missing evaluation or integrity error')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--final', action='store_true', help='Require all96 terminal with no unscored/invalid records')
    args = parser.parse_args(argv)
    result = summarize(args.run, final=args.final)
    print(json.dumps({key: result[key] for key in ['status', 'remaining_96', 'development_4', 'combined_100',
                                                 'completed_but_unscored_or_evaluation_failed', 'integrity_errors']},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
