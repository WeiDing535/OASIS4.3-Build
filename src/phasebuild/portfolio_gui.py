"""GUI adapter for the fixed portfolio; original numerical workflows unchanged."""
from pathlib import Path

from .portfolio import PROTOCOLS, collect, rank_candidates
from .data import fingerprint

NAMES = {'pplus_w14': 'P+ / Wang14', 'pplus_w8': 'P+ / Wang8',
         'stable_w14': 'stability shrinkage P+ / Wang14', 'pplus_w8_weight': 'P+ Wang8 + refinement weight search'}


def locate(directory):
    from .oasis_workflow import read_json
    root = Path(directory).resolve()
    for path in [root, root / '00_portfolio']:
        if (path / 'PLAN.json').is_file() and isinstance(read_json(path / 'PLAN.json').get('branches'), list):
            return path
    return None


def describe(cfg):
    from .oasis_workflow import describe_plan
    result = describe_plan({**cfg, 'explore': False})
    result['effective_options'].update(method='P+ four-mode comparison', dm_engine='Wang14 / Wang8 + classic40',
                                       explore=True, explore_workers=cfg.get('explore_workers', 2),
                                       selection_uses_Rfree=False, automatic_merge=False)
    for key in ['trust_stability', 'trust_work_weight', 'trust_wang_radius', 'trust_partial_cleanup', 'trust_cleanup']:
        result['effective_options'].pop(key, None)
    result['stages'] = [{'id': name, 'name': NAMES[name],
                         'directory': str(Path(cfg['output']) / '00_portfolio' / name)} for name in PROTOCOLS]
    result['fixed_protocol'] = ('four fixed pipelines use the same inputs, free set, seeds, and model-building budget; at most %d run concurrently. '
        'each pipeline applies full-chain/partial-model cleanup when its gate is met and refinement after a fixed 30 cycles; the weight search branch has four additional work-set trial refinements. '
        'a tentative recommendation is made by sequence support, repeat fragments, and work-set fit; Rfree is displayed only and not used for ranking. '
        'all results are retained; phases and models are not merged. '
        'when there is no single superior candidate, ties are retained for manual review.' % cfg.get('explore_workers', 2))
    return result


def validate_import(directory):
    from .oasis_workflow import read_json
    path = locate(directory)
    if path is None:
        raise ValueError('please select a result directory containing PLAN.json for the multi-pipeline workflow')
    plan = read_json(path / 'PLAN.json')
    names = [row.get('name') for row in plan.get('branches', [])]
    if not names or len(set(names)) != len(names) or set(names) - set(PROTOCOLS):
        raise ValueError('multi-pipeline record contains unknown or duplicate branches')


def collect_results(directory):
    from .oasis_workflow import read_json, collect_results as collect_single
    root = Path(directory).resolve()
    pipe = locate(root) or root / '00_portfolio'
    if locate(root):
        validate_import(root)
    declaration = read_json(pipe / 'PLAN.json')
    report = read_json(pipe / 'PORTFOLIO_REPORT.json')
    cfg = read_json(root / 'gui_job.json').get('config', {})
    stages, candidates, artifacts = [], [], []
    declared = declaration.get('branches', [{'name': name} for name in PROTOCOLS])
    final_by_name = {row['name']: row for row in report.get('candidates', [])}
    verified = []
    for branch in declared:
        name = branch['name']
        state = read_json(pipe / name / 'pipeline_report.json')
        final = final_by_name.get(name)
        status = final['status'] if final else state.get('status', 'waiting')
        metric = state.get('primary_metric') or {}
        model = pipe / name / 'primary_model.cif'
        if final and status == 'completed':
            checked = collect(name, pipe / name, branch['config'])
            for key in ['primary_metric', 'sequence_audit']:
                if checked[key] != final[key]:
                    raise ValueError('candidate statistics do not match the verified results: ' + name)
            if fingerprint(model)['sha256'] != final['model']['sha256']:
                raise ValueError('candidate model hash mismatch: ' + name)
            verified.append(checked)
        elif final:
            verified.append(final)
        row = dict(name=name, title=NAMES[name], status=status, r_work=metric.get('r_work'),
                   r_free=metric.get('r_free'), coverage=(final or {}).get('sequence_audit', {}).get('coverage'),
                   recommended=report.get('recommended') == name,
                   model=str(model.relative_to(root)) if final and status == 'completed' else None,
                   error=(final or {}).get('error') or state.get('error'))
        candidates.append(row)
        stages.append(dict(id=name, name=NAMES[name], status=status))
        if row['model']:
            artifacts.append(dict(path=row['model'], role=NAMES[name] + ' · candidate model', group='per-pipeline candidates'))
    if report:
        if set(final_by_name) != {b['name'] for b in declared}:
            raise ValueError('multi-pipeline result branches are inconsistent with the plan')
        reranked = rank_candidates(verified)
        if report.get('input_integrity_error') or report.get('inputs_unchanged') is not True:
            raise ValueError('multi-pipeline input integrity check failed')
        for key in ['recommended', 'review_candidates', 'pareto_candidates']:
            if reranked[key] != report.get(key):
                raise ValueError('recommendation is inconsistent with the fixed ranking rule')
    recommended = report.get('recommended')
    representative = None
    equivalent = []
    frontier = [row for row in verified if row['name'] in report.get('pareto_candidates', [])]
    if not recommended and report.get('status') == 'completed' and frontier:
        hashes = {row.get('model', {}).get('sha256') for row in frontier}
        if None not in hashes and len(hashes) == 1:
            equivalent = [row['name'] for row in frontier]
            representative = min(equivalent, key=list(PROTOCOLS).index)
    for row in candidates:
        row['representative'] = row['name'] == representative
    selected = recommended or representative
    if selected:
        selected_root = pipe / selected
        result = collect_single(selected_root)
        prefix = selected_root.relative_to(root)
        for file in result['artifacts']:
            file['path'] = str(prefix / file['path'])
        for key in ['final_model', 'final_phases']:
            if result.get(key):
                result[key] = str(prefix / result[key])
        for metric in result['metrics']:
            if metric.get('model'):
                metric['model'] = str(prefix / metric['model'])
        result['artifacts'].extend(artifacts)
    else:
        result = dict(config=cfg, sad=None, dm=None, dm_history=[], metrics=[], artifacts=artifacts,
                      primary_metric=None, primary_complete=False, build_model=True, final_model=None,
                      final_phases=None, preview=[], residues=0, chains=0, sample={}, partition={},
                      solvent_fraction=None, solvent_source=None, error=None)
    for name, role in [('PORTFOLIO_REPORT.json', 'multi-pipeline statistics and recommendation record'), ('REPORT.md', 'multi-pipeline report'),
                       ('PLAN.json', 'pre-declared pipeline plan'), ('SOURCES.json', 'source code identity'), ('INPUTS.json', 'input validation')]:
        path = pipe / name
        if path.is_file():
            result['artifacts'].append(dict(path=str(path.relative_to(root)), role=role, group='run record'))
    if (root / 'gui.log').is_file():
        result['artifacts'].append(dict(path='gui.log', role='multi-pipeline main log', group='log'))
    result.update(root=str(root), status=report.get('status', 'running' if declaration else 'waiting'),
                  stages=stages, method='pplus_bundled', method_name='P+ multi-pipeline exploration', dm_engine='external CCP4 Parrot · candidate recommendation',
                  pipeline_complete=report.get('status') == 'completed',
                  portfolio=dict(candidates=candidates, recommended=recommended,
                                 representative=representative, equivalent_candidates=equivalent,
                                 review_candidates=report.get('review_candidates', []),
                                 reason=report.get('reason'), selection_uses_Rfree=False,
                                 policy=report.get('policy')))
    return result


def enrich_job(job):
    from .oasis_workflow import read_json
    root = Path(job['output'])
    pipe = locate(root) or root / '00_portfolio'
    plan = read_json(pipe / 'PLAN.json')
    final = read_json(pipe / 'PORTFOLIO_REPORT.json')
    revisions, labels, completed = [], [], 0
    logs = [root / 'gui.log']
    for branch in plan.get('branches', []):
        name = branch['name']
        state = read_json(pipe / name / 'pipeline_report.json')
        revisions.append([name, state.get('updated')])
        completed += state.get('status') == 'completed'
        if state.get('status') == 'running':
            labels.append(NAMES[name] + ': ' + state.get('current_name', 'starting'))
            logs.append(pipe / (name + '.log'))
    job['pipeline_progress'] = dict(completed=completed, total=len(plan.get('branches', [])) or 4,
                                    current='; '.join(labels) or ('multi-pipeline completed' if final else 'multi-pipeline precheck / queued'),
                                    revision=[revisions, final.get('status')])
    parts = []
    for path in logs:
        if path.is_file():
            with path.open('rb') as handle:
                handle.seek(max(0, path.stat().st_size - 18000))
                parts.append(str(path.relative_to(root)) + '\n' + handle.read().decode(errors='replace'))
    job['log'] = '\n\n'.join(parts)[-60000:]
    return job
