"""Configuration and stage-plan tests for the opt-in density-modification upgrade."""
import pytest

from phasebuild.sad_pipeline import DEFAULTS, validate, stage_plan


BASE=dict(method='pplus_trust', build_model=False, fpp='Se=3.8',
          sad='/nonexistent/sad.mtz', native='/nonexistent/native.mtz',
          sequence='/nonexistent/seq.fasta', heavy='/nonexistent/heavy.pdb',
          output='/nonexistent/out')


def _validate(monkeypatch, **overrides):
    # Skip filesystem preflight: these tests only exercise option validation.
    cfg={**BASE, **overrides}
    return validate(cfg, inspect=False)


def test_defaults_preserve_legacy_behaviour():
    assert DEFAULTS['dm_engine']=='resolve'
    assert DEFAULTS['feedback_cycles']==0
    assert DEFAULTS['dm_solvent'] is None
    assert DEFAULTS['dm_mask_radius']==3.0
    assert DEFAULTS['dm_anisotropy'] is True


def test_unknown_engine_is_rejected():
    with pytest.raises(ValueError):
        _validate(None, dm_engine='parrot2')


def test_parrot_requires_p_plus_method():
    with pytest.raises(ValueError):
        _validate(None, method='likelihood', dm_engine='parrot')


@pytest.mark.parametrize('solvent',[0.0,1.0,-0.1,1.5,float('nan')])
def test_solvent_override_must_be_a_fraction(solvent):
    with pytest.raises(ValueError):
        _validate(None, dm_engine='parrot', dm_solvent=solvent)


@pytest.mark.parametrize('radius',[0.5,11.0,float('inf')])
def test_mask_radius_bounds(radius):
    with pytest.raises(ValueError):
        _validate(None, dm_engine='parrot', dm_mask_radius=radius)


@pytest.mark.parametrize('cycles',[-1,9,True])
def test_feedback_cycle_bounds(cycles):
    with pytest.raises(ValueError):
        _validate(None, dm_engine='parrot', feedback_cycles=cycles)


def test_feedback_requires_parrot_engine():
    with pytest.raises(ValueError):
        _validate(None, dm_engine='resolve', feedback_cycles=1)


def test_unknown_building_recipe_is_rejected():
    with pytest.raises(ValueError):
        _validate(None, building_recipe='pipeline2')


def _stage_plan(monkeypatch, tmp_path, **overrides):
    # Real placeholder files so the stage planner can fingerprint its inputs.
    for name in ('sad','native','sequence','heavy','reference'):
        (tmp_path/name).write_text('placeholder')
    # A runnable placeholder keeps the Parrot route independent of any CCP4
    # installed on the test host; without it validate() falls back to the
    # internal engine and the planned command silently changes.
    parrot=tmp_path/'cparrot'
    parrot.write_text('#!/bin/sh\nexit 1\n')
    parrot.chmod(0o755)
    cfg=_validate(None, method='pplus_internal', dm_executable=str(parrot),
                  sad=str(tmp_path/'sad'),
                  native=str(tmp_path/'native'), sequence=str(tmp_path/'sequence'),
                  heavy=str(tmp_path/'heavy'), reference=str(tmp_path/'reference'),
                  output=str(tmp_path/'out'), **overrides)
    return cfg, stage_plan(cfg)


def test_stage_plan_default_dm_command_is_untouched(tmp_path):
    cfg,stages=_stage_plan(None,tmp_path)
    dm=[s for s in stages if s['id']=='dm'][0]
    assert 'run_pplus_trust_stage.py' in ' '.join(dm['command'])
    assert dm['outputs']==[cfg['output']+'/03_dm/map_phases.mtz',
                          cfg['output']+'/03_dm/map_phases.json']


def test_stage_plan_routes_parrot_engine_to_trace_feedback(tmp_path):
    cfg,stages=_stage_plan(None,tmp_path,dm_engine='parrot',feedback_cycles=2,dm_solvent=0.5)
    dm=[s for s in stages if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert 'trace/Fcalc' in dm['name']
