"""Factory for explicit, fixed original-phase REFMAC recycling adapters.

Use in a separate frozen experiment: hybrid_recycle.refmac_refine =
make_fixed_prior_refiner(original_phase_mtz, expected_sha256=known_hash).
The returned call signature matches refinement.refmac_refine. Every call uses
the same original prior; posterior model phases are never promoted to priors.
"""
import json
from pathlib import Path
from .data import fingerprint
from .phase_refinement import merge_experimental_phases, refmac_with_experimental_prior, validate_phase_blur


def make_fixed_prior_refiner(original_phase_mtz, expected_sha256=None, representation='hl', free_value=0, source_kind='original_sad', phase_blur=None):
    phase_blur=validate_phase_blur(phase_blur)
    if isinstance(free_value, bool) or not isinstance(free_value, int):
        raise ValueError('Free-set test value must be an integer')
    if source_kind not in ('original_sad','density_modified_work'):raise ValueError('Invalid phase source kind')
    phase_prefix='DM_' if source_kind=='density_modified_work' else 'SAD_'
    declared_free_value = free_value
    fixed = fingerprint(original_phase_mtz)
    if expected_sha256 is not None and fixed['sha256'] != expected_sha256:
        raise ValueError('Original experimental phase source hash does not match declared protocol')
    if representation not in ['hl', 'phifom']:
        raise ValueError('Invalid fixed experimental prior representation')

    def refine(mtz_path, model_path, output, cycles=10, ccp4=None,
               fp='F', sigfp='SIGF', free='FreeR_flag', free_value=declared_free_value, timeout=3600):
        if (fp, sigfp, free, free_value) != ('F', 'SIGF', 'FreeR_flag', declared_free_value):
            raise ValueError('This audited fixed-prior branch requires F/SIGF/FreeR_flag with the declared original test value')
        if fingerprint(fixed['path']) != fixed:
            raise RuntimeError('Fixed original experimental phase source changed during recycling')
        native = fingerprint(mtz_path)
        out = Path(output).resolve()
        cache = out.parent / 'fixed_experimental_prior_inputs'
        cache.mkdir(parents=True, exist_ok=True)
        # Full native/prior hashes remain mandatory in the sidecar audit. Short
        # names avoid the legacy CCP4 pathname buffer; collisions are rejected
        # by the exact provenance check below rather than silently reused.
        stem = native['sha256'][:16] + '_' + fixed['sha256'][:16] + '_' + representation
        if source_kind!='original_sad':stem+='_'+source_kind
        if declared_free_value != 0:
            stem += '_free' + str(declared_free_value)
        merged = cache / (stem + '.mtz')
        audit_path = merged.with_suffix('.audit.json')
        if not merged.exists():
            audit = merge_experimental_phases(mtz_path, fixed['path'], merged, representation, **({'phase_prefix':phase_prefix} if phase_prefix!='SAD_' else {}), **({'free_value': declared_free_value} if declared_free_value != 0 else {}))
        else:
            if not audit_path.exists():
                raise RuntimeError('Existing merged prior cache has no audit')
            audit = json.loads(audit_path.read_text())
            if (audit['native'] != native or audit['phase_source'] != fixed or
                    audit['representation'] != representation or
                    audit['output'] != fingerprint(merged)):
                raise RuntimeError('Merged experimental-prior cache provenance changed')
        report = refmac_with_experimental_prior(merged, model_path, out, use_prior=True,
            cycles=cycles, representation=representation, timeout=timeout, ccp4=ccp4,
            **({'free_value': declared_free_value} if declared_free_value != 0 else {}),**({'phase_prefix':phase_prefix} if phase_prefix!='SAD_' else {}),**({'phase_blur':phase_blur} if phase_blur is not None else {}))
        report.update(wrapper_source=fingerprint(__file__), original_native_input=native,
            fixed_original_experimental_prior=fixed, phase_merge_audit=audit,
            recycling_phase_policy='Same fixed original experimental prior every round; no posterior model-derived prior substitution',
            scientific_protocol='Experimental-phase-aware likelihood; differs from amplitude-only frozen benchmark')
        if source_kind=='density_modified_work':
            if audit['matched_free_prior_rows']!=0:raise RuntimeError('DM restraint supplied on original free rows')
            report['fixed_density_modified_work_prior']=report.pop('fixed_original_experimental_prior')
            report.update(phase_source_kind=source_kind,recycling_phase_policy='Same frozen work-only DM phase/FOM restraint each round; not original SAD likelihood, never model posterior substitution',scientific_protocol='Explicit approximate phase restraint after work-only density modification; dependence on measured amplitudes is not an independent likelihood assumption')
        (out / 'refinement_report.json').write_text(json.dumps(report, indent=2) + '\n')
        return report

    refine.fixed_original_experimental_prior = fixed
    refine.representation = representation
    refine.free_value = declared_free_value
    return refine
