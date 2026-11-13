"""Configuration and stage-plan tests for the two density-modification engines.

There are only two engines now: CCP4 Parrot (the calibrated protocol, and the
only one P+ may use) and the in-tree MLHL kernel used by 连续 SAD / 旧版 P+ /
内部 P+.  The copied RESOLVE/CCTBX tree is gone, so no test here may assume a
density/ directory.
"""
from pathlib import Path

import pytest

from phasebuild.sad_pipeline import DEFAULTS, fallback_notice, validate, stage_plan


BASE=dict(method='pplus_trust', build_model=False, fpp='Se=3.8',
          sad='/nonexistent/sad.mtz', native='/nonexistent/native.mtz',
          sequence='/nonexistent/seq.fasta', heavy='/nonexistent/heavy.pdb',
          output='/nonexistent/out')


def _validate(monkeypatch=None, **overrides):
    # Skip filesystem preflight: these tests only exercise option validation.
    return validate({**BASE, **overrides}, inspect=False)


def _ccp4(tmp_path, names=('phaser','cparrot','cbuccaneer','refmac5')):
    runtime=tmp_path/'ccp4'
    for name in names:
        binary=runtime/'bin'/name
        binary.parent.mkdir(parents=True,exist_ok=True)
        binary.write_text('#!/bin/sh\nexit 1\n');binary.chmod(0o755)
    lib=runtime/'lib/data/reference_structures';lib.mkdir(parents=True,exist_ok=True)
    for suffix in ('.mtz','.pdb'):(lib/('reference-1tqw'+suffix)).write_text('placeholder')
    return runtime


def _inputs(tmp_path):
    for name in ('sad','native','sequence','heavy','reference'):
        (tmp_path/name).write_text('placeholder')
    return dict(sad=str(tmp_path/'sad'),native=str(tmp_path/'native'),
                sequence=str(tmp_path/'sequence'),heavy=str(tmp_path/'heavy'),
                reference=str(tmp_path/'reference'),output=str(tmp_path/'out'))


def test_defaults_use_ccp4_parrot():
    assert DEFAULTS['dm_engine']=='parrot'
    assert DEFAULTS['feedback_cycles']==0
    assert DEFAULTS['dm_solvent'] is None
    assert DEFAULTS['dm_mask_radius']==3.0
    assert DEFAULTS['dm_anisotropy'] is True


def test_removed_resolve_engine_is_rejected():
    with pytest.raises(ValueError,match='已移除'):
        _validate(dm_engine='resolve')
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot2')


@pytest.mark.parametrize('method',['pplus_trust','pplus_bundled'])
def test_the_two_p_plus_routes_refuse_the_internal_engine(method):
    with pytest.raises(ValueError,match='CCP4 Parrot'):
        _validate(method=method,dm_engine='internal')


@pytest.mark.parametrize('method',['likelihood','pplus','pplus_internal','pplus_trust','pplus_bundled'])
def test_parrot_is_the_engine_for_every_phasing_method(method,tmp_path,monkeypatch):
    """Both SAD front ends write an observation likelihood, so Parrot fits all."""
    from phasebuild import bundled_pplus, parrot_dm
    runtime=_ccp4(tmp_path)
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k:str(runtime/'bin/cparrot'))
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda value,dm_engine=None:Path(value))
    if method=='pplus_trust':(tmp_path/'prior').write_text('placeholder')
    extra={'observation_prior':str(tmp_path/'prior')} if method=='pplus_trust' else {}
    cfg=_validate(method=method,dm_engine='parrot',dm_executable=str(runtime/'bin/cparrot'),
                  runtime=str(runtime),**extra,**_inputs(tmp_path))
    assert cfg['dm_engine']=='parrot'
    assert 'dm_engine_fallback' not in cfg and 'dm_engine_requested' not in cfg


def test_missing_parrot_falls_back_for_the_continuous_routes(tmp_path,monkeypatch):
    """No CCP4 must not stop 连续 SAD, but the fallback has to be visible."""
    from phasebuild import parrot_dm
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('Parrot (cparrot) was not found')))
    cfg=_validate(method='likelihood',dm_engine='parrot',dm_executable='/nonexistent/cparrot',
                  **_inputs(tmp_path))
    assert cfg['dm_engine']=='internal' and cfg['dm_engine_requested']=='parrot'
    assert cfg['feedback_cycles']==0
    notice=fallback_notice(cfg)
    assert '内置 MLHL' in notice and 'cparrot' in notice
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_oasis_stage.py' in ' '.join(dm['command'])


def test_missing_parrot_is_fatal_for_the_p_plus_routes(tmp_path,monkeypatch):
    from phasebuild import parrot_dm
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('Parrot (cparrot) was not found')))
    with pytest.raises(ValueError,match='cparrot'):
        _validate(method='pplus_trust',dm_engine='parrot',dm_executable='/nonexistent/cparrot',
                  observation_prior=str(tmp_path/'prior'),**_inputs(tmp_path))


@pytest.mark.parametrize('solvent',[0.0,1.0,-0.1,1.5,float('nan')])
def test_solvent_override_must_be_a_fraction(solvent):
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot',dm_solvent=solvent)


@pytest.mark.parametrize('radius',[0.5,11.0,float('inf')])
def test_mask_radius_bounds(radius):
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot',dm_mask_radius=radius)


@pytest.mark.parametrize('cycles',[-1,9,True])
def test_feedback_cycle_bounds(cycles):
    with pytest.raises(ValueError):
        _validate(dm_engine='parrot',feedback_cycles=cycles)


def test_feedback_requires_parrot_engine():
    with pytest.raises(ValueError):
        _validate(method='likelihood',dm_engine='internal',feedback_cycles=1)


def test_unknown_building_recipe_is_rejected():
    with pytest.raises(ValueError):
        _validate(building_recipe='pipeline2')


def test_parrot_dm_stage_is_the_trace_feedback_driver(tmp_path,monkeypatch):
    from phasebuild import bundled_pplus, parrot_dm
    runtime=_ccp4(tmp_path)
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k:str(runtime/'bin/cparrot'))
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda value,dm_engine=None:Path(value))
    cfg=_validate(method='pplus_internal',build_model=True,dm_engine='parrot',
                  feedback_cycles=2,dm_solvent=.5,runtime=str(runtime),
                  dm_executable=str(runtime/'bin/cparrot'),**_inputs(tmp_path))
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert 'trace/Fcalc' in dm['name']


def test_internal_mlhl_stage_stays_on_the_oasis_driver(tmp_path):
    cfg=_validate(method='likelihood',build_model=True,dm_engine='internal',**_inputs(tmp_path))
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_oasis_stage.py' in ' '.join(dm['command']) and 'run_trace_feedback.py' not in ' '.join(dm['command'])


def test_phases_only_keeps_parrot_but_drops_the_feedback(tmp_path,monkeypatch):
    """No tracer means no Fcalc feedback; the DM map itself is still Parrot's."""
    from phasebuild import bundled_pplus, parrot_dm
    runtime=_ccp4(tmp_path)
    monkeypatch.setattr(parrot_dm,'locate',lambda *a,**k:str(runtime/'bin/cparrot'))
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda value,dm_engine=None:Path(value))
    cfg=_validate(method='likelihood',build_model=False,dm_engine='parrot',feedback_cycles=2,
                  runtime=str(runtime),dm_executable=str(runtime/'bin/cparrot'),**_inputs(tmp_path))
    assert cfg['dm_engine']=='parrot' and cfg['feedback_cycles']==0
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert 'Parrot' in dm['name']


def test_engine_layout_needs_phaser_but_no_density_tree(tmp_path):
    from phasebuild import bundled_pplus
    runtime=_ccp4(tmp_path,names=('phaser',))
    assert bundled_pplus.engine_layout(runtime,'parrot')==runtime.resolve()
    assert bundled_pplus.engine_layout(runtime) == runtime.resolve()
    with pytest.raises(ValueError,match='未知密度修饰引擎'):
        bundled_pplus.engine_layout(runtime,'resolve')
    with pytest.raises(ValueError,match='P\\+ 引擎缺失'):
        bundled_pplus.engine_layout(tmp_path/'empty','parrot')


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
