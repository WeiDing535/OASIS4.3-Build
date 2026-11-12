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


def test_missing_parrot_is_recorded_instead_of_hidden(monkeypatch,tmp_path):
    """A run without CCP4 must say which engine really ran."""
    from phasebuild import parrot_dm
    from phasebuild.sad_pipeline import fallback_notice
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('Parrot (cparrot) was not found')))
    for name in ('sad','native','sequence','heavy','reference'):
        (tmp_path/name).write_text('placeholder')
    cfg=_validate(None, method='pplus_internal', dm_engine='parrot',
                  feedback_cycles=2, dm_cycles=40, dm_executable='/nonexistent/cparrot',
                  sad=str(tmp_path/'sad'), native=str(tmp_path/'native'),
                  sequence=str(tmp_path/'sequence'), heavy=str(tmp_path/'heavy'),
                  reference=str(tmp_path/'reference'), output=str(tmp_path/'out'))
    assert cfg['dm_engine']=='resolve'
    assert cfg['dm_engine_requested']=='parrot'
    assert 'cparrot' in cfg['dm_engine_fallback']
    assert cfg['feedback_cycles']==0
    notice=fallback_notice(cfg)
    assert '已回退' in notice and 'Parrot' in notice and 'RESOLVE' in notice
    assert fallback_notice({})==''
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' not in ' '.join(dm['command'])


def test_engine_label_marks_the_fallback():
    from phasebuild.oasis_workflow import _engine_label
    label=_engine_label('pplus_bundled',{'dm_engine':'resolve','dm_engine_requested':'parrot',
                                         'dm_engine_fallback':'missing'},{})
    assert '已回退' in label
    assert '已回退' not in _engine_label('pplus_bundled',{'dm_engine':'resolve'},{})


def test_planned_config_shows_the_internal_stages(monkeypatch,tmp_path):
    """The stage list must match the engine the run will really use."""
    from phasebuild import parrot_dm
    from phasebuild.oasis_workflow import planned_config
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('Parrot (cparrot) was not found')))
    for name in ('sad','native','sequence','heavy','reference'):
        (tmp_path/name).write_text('placeholder')
    base=dict(_validate(None, method='pplus_internal', dm_engine='parrot',
                        feedback_cycles=2, dm_executable='/nonexistent/cparrot',
                        sad=str(tmp_path/'sad'), native=str(tmp_path/'native'),
                        sequence=str(tmp_path/'sequence'), heavy=str(tmp_path/'heavy'),
                        reference=str(tmp_path/'reference'), output=str(tmp_path/'out')))
    # A raw (not yet validated) plan request is what the window sends.
    raw={k:v for k,v in base.items() if k not in ('dm_engine_fallback','dm_engine_requested')}
    raw.update(dm_engine='parrot',dm_executable='/nonexistent/cparrot')
    values=planned_config(raw)
    assert values['dm_engine']=='resolve' and values['feedback_cycles']==0
    dm=[s for s in stage_plan(values) if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' not in ' '.join(dm['command'])


def test_plan_probe_detects_parrot_before_the_run(monkeypatch):
    """The plan must warn up front, not only once the run has started."""
    from phasebuild import parrot_dm
    from phasebuild.oasis_workflow import engine_fallback
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('Parrot (cparrot) was not found')))
    assert 'cparrot' in engine_fallback({'dm_engine':'parrot'})
    # A recorded fallback needs no probe at all.
    assert engine_fallback({'dm_engine':'resolve','dm_engine_fallback':'recorded'})=='recorded'
    # A reachable Parrot means no warning.
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k: '/opt/ccp4/bin/cparrot')
    assert engine_fallback({'dm_engine':'parrot'})==''
    assert engine_fallback({'dm_engine':'resolve'})==''


def test_parrot_layout_does_not_require_the_resolve_density_tree(tmp_path, monkeypatch):
    """A CCP4-derived runtime has phaser + cparrot but no RESOLVE/CCTBX tree."""
    from phasebuild.bundled_pplus import engine_layout
    runtime = tmp_path / 'runtime'
    for name in ('phaser', 'cparrot'):
        binary = runtime / 'bin' / name
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_text('#!/bin/sh\n')
        binary.chmod(0o755)
    assert engine_layout(runtime, 'parrot') == (runtime.resolve(), None)
    with pytest.raises(ValueError, match='P\\+ 引擎缺失'):
        engine_layout(runtime, 'resolve')
    with pytest.raises(ValueError, match='P\\+ 引擎缺失'):
        engine_layout(runtime)


def test_bundled_runtime_wins_over_a_host_ccp4_installation(tmp_path, monkeypatch):
    from phasebuild import parrot_dm
    bundled = tmp_path / 'bundle/runtime/bin/cparrot'
    bundled.parent.mkdir(parents=True)
    bundled.write_text('#!/bin/sh\n')
    bundled.chmod(0o755)
    host = tmp_path / 'host/bin/cparrot'
    host.parent.mkdir(parents=True)
    host.write_text('#!/bin/sh\n')
    host.chmod(0o755)
    monkeypatch.setenv('PHASEBUILD_RUNTIME', str(bundled.parents[1]))
    monkeypatch.setattr(parrot_dm.shutil, 'which', lambda name: str(host))
    assert parrot_dm.locate() == str(bundled)
    monkeypatch.delenv('PHASEBUILD_RUNTIME')
    assert parrot_dm.locate(runtime=bundled.parents[1]) == str(bundled)
