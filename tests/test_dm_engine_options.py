"""Configuration and stage-plan tests for the opt-in density-modification upgrade."""
import pytest

from phasebuild.sad_pipeline import DEFAULTS, validate, stage_plan


BASE=dict(method='pplus_trust', build_model=False, fpp='Se=3.8',
          sad='/nonexistent/sad.mtz', native='/nonexistent/native.mtz',
          sequence='/nonexistent/seq.fasta', heavy='/nonexistent/heavy.pdb',
          output='/nonexistent/out')


def _validate(**overrides):
    # Skip filesystem preflight and the Parrot binary probe: these tests only
    # exercise option validation, not engine discovery.
    import phasebuild.parrot_dm as pd
    cfg={**BASE, **overrides}
    original=pd.locate
    pd.locate=lambda *a, **k: '/usr/bin/true'
    try:
        return validate(cfg, inspect=False)
    finally:
        pd.locate=original


def test_defaults_route_through_external_parrot():
    assert DEFAULTS['dm_engine']=='parrot'
    assert DEFAULTS['feedback_cycles']==0
    assert DEFAULTS['dm_solvent'] is None
    assert DEFAULTS['dm_mask_radius']==3.0
    assert DEFAULTS['dm_anisotropy'] is True


def test_unknown_engine_is_rejected():
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot2')


def test_legacy_methods_tolerate_nominal_engine_choice():
    # likelihood/pplus run the legacy reference-based DM; dm_engine is nominal
    # there. BASE uses placeholder paths, so validation proceeds past the
    # engine checks and fails on the missing input file instead.
    with pytest.raises(ValueError,match='Missing input file'):
        _validate(method='likelihood', dm_engine='parrot')


@pytest.mark.parametrize('solvent',[0.0,1.0,-0.1,1.5,float('nan')])
def test_solvent_override_must_be_a_fraction(solvent):
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot', dm_solvent=solvent)


@pytest.mark.parametrize('radius',[0.5,11.0,float('inf')])
def test_mask_radius_bounds(radius):
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot', dm_mask_radius=radius)


@pytest.mark.parametrize('cycles',[-1,9,True])
def test_feedback_cycle_bounds(cycles):
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot', feedback_cycles=cycles)


def test_feedback_requires_parrot_engine():
    with pytest.raises(ValueError):
        _validate(dm_engine='resolve', feedback_cycles=1)


def test_unknown_building_recipe_is_rejected():
    with pytest.raises(ValueError):
        _validate(building_recipe='pipeline2')


def _stage_plan(tmp_path, **overrides):
    # Real placeholder files so the stage planner can fingerprint its inputs.
    for name in ('sad','native','sequence','heavy','reference','observation_prior'):
        (tmp_path/name).write_text('placeholder')
    cfg=_validate(sad=str(tmp_path/'sad'),
                  native=str(tmp_path/'native'), sequence=str(tmp_path/'sequence'),
                  heavy=str(tmp_path/'heavy'), reference=str(tmp_path/'reference'),
                  observation_prior=str(tmp_path/'observation_prior'),
                  output=str(tmp_path/'out'), **overrides)
    return cfg, stage_plan(cfg)


def test_stage_plan_default_engine_routes_to_trace_feedback(tmp_path):
    # The only DM engine is external CCP4 Parrot, always routed through the
    # trace-feedback stage script.
    cfg,stages=_stage_plan(tmp_path)
    dm=[s for s in stages if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert dm['outputs']==[cfg['output']+'/03_dm/map_phases.mtz',
                          cfg['output']+'/03_dm/map_phases.json']


def test_stage_plan_routes_parrot_engine_to_trace_feedback(tmp_path):
    cfg,stages=_stage_plan(tmp_path,feedback_cycles=2,dm_solvent=0.5)
    dm=[s for s in stages if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert 'trace/Fcalc' in dm['name']
