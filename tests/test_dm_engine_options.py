"""Configuration and stage-plan tests for the density-modification engines.

There are two engines and both are in-tree: oasis-dm v2 (the default, and the only
one the P+ routes use) and the MLHL kernel used by 连续 SAD / 旧版 P+ / 内部 P+.
Every external density modifier -- the copied RESOLVE/CCTBX tree and CCP4
The external density-modification engine -- has been removed from the program,
so no test here may assume a density/ directory or an external DM binary.
"""
import sys
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


def _ccp4(tmp_path, names=('phaser','cbuccaneer','refmac5')):
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


def test_defaults_use_the_in_tree_oasis_dm_engine():
    assert DEFAULTS['dm_engine']=='oasis-dm'
    assert DEFAULTS['dm_cmbwt']==0.4 and DEFAULTS['dm_holdout']==0.05
    assert DEFAULTS['dm_hist_lib']==''
    assert DEFAULTS['feedback_cycles']==0
    assert DEFAULTS['dm_solvent'] is None
    # The old external-engine knobs are gone, not merely defaulted away.
    assert 'dm_mask_radius' not in DEFAULTS and 'dm_anisotropy' not in DEFAULTS


@pytest.mark.parametrize('engine',['resolve','oasis_dm_v2',''])
def test_only_the_two_in_tree_engines_are_accepted(engine):
    with pytest.raises(ValueError,match='Unknown density-modification engine'):
        _validate(dm_engine=engine)


@pytest.mark.parametrize('weight',[-0.1,0.0,2.5,float('nan'),True,'strong'])
def test_cmbwt_must_be_a_positive_weight_or_auto(weight):
    with pytest.raises(ValueError,match='dm_cmbwt'):
        _validate(dm_engine='oasis-dm',dm_cmbwt=weight)


def test_cmbwt_auto_is_accepted(tmp_path):
    cfg=_validate(method='likelihood',dm_engine='oasis-dm',dm_cmbwt='auto',
                  **_inputs(tmp_path))
    assert cfg['dm_cmbwt']=='auto'


@pytest.mark.parametrize('fraction',[-0.1,0.6,True,'five'])
def test_holdout_must_be_a_small_fraction(fraction):
    with pytest.raises(ValueError,match='dm_holdout'):
        _validate(dm_engine='oasis-dm',dm_holdout=fraction)


@pytest.mark.parametrize('method',['pplus_trust','pplus_bundled'])
def test_the_two_p_plus_routes_refuse_the_internal_engine(method):
    with pytest.raises(ValueError,match='只能用 oasis-dm'):
        _validate(method=method,dm_engine='internal')


@pytest.mark.parametrize('engine',['oasis-dm','internal'])
@pytest.mark.parametrize('method',['likelihood','pplus','pplus_internal'])
def test_each_dm_engine_fits_the_continuous_routes(method,engine,tmp_path,monkeypatch):
    """Both in-tree engines serve every route that writes an observation likelihood."""
    from phasebuild import oasis_dm
    monkeypatch.setattr(oasis_dm,'locate',lambda *a,**k:(sys.executable,'/vendor/oasis_dm_v2.py'))
    cfg=_validate(method=method,dm_engine=engine,**_inputs(tmp_path))
    assert cfg['dm_engine']==engine
    assert 'dm_engine_fallback' not in cfg and 'dm_engine_requested' not in cfg


def test_oasis_dm_serves_both_p_plus_routes(tmp_path,monkeypatch):
    from phasebuild import bundled_pplus, oasis_dm
    runtime=_ccp4(tmp_path)
    monkeypatch.setattr(oasis_dm,'locate',lambda *a,**k:(sys.executable,'/vendor/oasis_dm_v2.py'))
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda value,dm_engine=None:Path(value))
    (tmp_path/'prior').write_text('placeholder')
    for method in ('pplus_trust','pplus_bundled'):
        extra={'observation_prior':str(tmp_path/'prior')} if method=='pplus_trust' else {}
        cfg=_validate(method=method,dm_engine='oasis-dm',runtime=str(runtime),
                      **extra,**_inputs(tmp_path))
        assert cfg['dm_engine']=='oasis-dm' and 'dm_engine_fallback' not in cfg


def test_missing_oasis_dm_falls_back_for_the_continuous_routes(tmp_path,monkeypatch):
    """A runtime without numpy/gemmi must not stop 连续 SAD, but must be visible."""
    from phasebuild import oasis_dm
    monkeypatch.setattr(oasis_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('oasis-dm v2 needs a Python with numpy, scipy and gemmi')))
    cfg=_validate(method='likelihood',dm_engine='oasis-dm',**_inputs(tmp_path))
    assert cfg['dm_engine']=='internal' and cfg['dm_engine_requested']=='oasis-dm'
    assert cfg['feedback_cycles']==0
    notice=fallback_notice(cfg)
    assert 'oasis-dm' in notice
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_oasis_stage.py' in ' '.join(dm['command'])


def test_missing_oasis_dm_is_fatal_for_the_p_plus_routes(tmp_path,monkeypatch):
    from phasebuild import oasis_dm
    monkeypatch.setattr(oasis_dm,'locate',lambda *a,**k: (_ for _ in ()).throw(
        ValueError('oasis-dm v2 needs a Python with numpy, scipy and gemmi')))
    with pytest.raises(ValueError,match='oasis-dm'):
        _validate(method='pplus_trust',dm_engine='oasis-dm',
                  observation_prior=str(tmp_path/'prior'),**_inputs(tmp_path))


@pytest.mark.parametrize('solvent',[0.0,1.0,-0.1,1.5,float('nan')])
def test_solvent_override_must_be_a_fraction(solvent):
    with pytest.raises(ValueError):
        _validate(dm_engine='oasis-dm',dm_solvent=solvent)


@pytest.mark.parametrize('cycles',[-1,9,True])
def test_feedback_cycle_bounds(cycles):
    with pytest.raises(ValueError):
        _validate(dm_engine='oasis-dm',feedback_cycles=cycles)


def test_feedback_requires_oasis_dm_engine():
    with pytest.raises(ValueError,match='oasis-dm density-modification engine'):
        _validate(method='likelihood',dm_engine='internal',feedback_cycles=1)


def test_unknown_building_recipe_is_rejected():
    with pytest.raises(ValueError):
        _validate(building_recipe='pipeline2')


def test_oasis_dm_stage_runs_the_trace_feedback_driver(tmp_path,monkeypatch):
    from phasebuild import bundled_pplus, oasis_dm
    runtime=_ccp4(tmp_path)
    monkeypatch.setattr(oasis_dm,'locate',lambda *a,**k:(sys.executable,'/vendor/oasis_dm_v2.py'))
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda value,dm_engine=None:Path(value))
    cfg=_validate(method='pplus_bundled',build_model=True,dm_engine='oasis-dm',
                  feedback_cycles=2,dm_solvent=.5,dm_cmbwt='auto',runtime=str(runtime),
                  **_inputs(tmp_path))
    assert cfg['dm_cycles']==12
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert 'oasis-dm v2' in dm['name'] and 'trace/Fcalc' in dm['name']


def test_internal_mlhl_stage_stays_on_the_oasis_driver(tmp_path):
    cfg=_validate(method='likelihood',build_model=True,dm_engine='internal',**_inputs(tmp_path))
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_oasis_stage.py' in ' '.join(dm['command']) and 'run_trace_feedback.py' not in ' '.join(dm['command'])


def test_phases_only_keeps_oasis_dm_but_drops_the_feedback(tmp_path,monkeypatch):
    """No tracer means no Fcalc feedback; the DM map itself is still oasis-dm's."""
    from phasebuild import bundled_pplus, oasis_dm
    runtime=_ccp4(tmp_path)
    monkeypatch.setattr(oasis_dm,'locate',lambda *a,**k:(sys.executable,'/vendor/oasis_dm_v2.py'))
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda value,dm_engine=None:Path(value))
    cfg=_validate(method='likelihood',build_model=False,dm_engine='oasis-dm',feedback_cycles=2,
                  runtime=str(runtime),**_inputs(tmp_path))
    assert cfg['dm_engine']=='oasis-dm' and cfg['feedback_cycles']==0
    dm=[s for s in stage_plan(cfg) if s['id']=='dm'][0]
    assert 'run_trace_feedback.py' in ' '.join(dm['command'])
    assert 'oasis-dm v2' in dm['name']


def test_engine_layout_needs_phaser_but_no_density_tree(tmp_path):
    from phasebuild import bundled_pplus
    runtime=_ccp4(tmp_path,names=('phaser',))
    assert bundled_pplus.engine_layout(runtime,'oasis-dm')==runtime.resolve()
    assert bundled_pplus.engine_layout(runtime,'internal')==runtime.resolve()
    assert bundled_pplus.engine_layout(runtime) == runtime.resolve()
    for gone in ('resolve',):
        with pytest.raises(ValueError,match='未知密度修饰引擎'):
            bundled_pplus.engine_layout(runtime,gone)
    with pytest.raises(ValueError,match='P\\+ 引擎缺失'):
        bundled_pplus.engine_layout(tmp_path/'empty','oasis-dm')


def test_feedback_model_hl_is_merged_into_the_experimental_quartet(tmp_path, monkeypatch):
    """The second feedback cycle feeds a model back; that merge must import cleanly."""
    import numpy as np
    import gemmi
    from phasebuild import oasis_dm

    rng = np.random.default_rng(5)
    hkl = np.array([(h, k, l) for h in range(6) for k in range(6) for l in range(6)][:200],
                   dtype=float)
    fp = 100.0 + 10.0 * rng.random(len(hkl))
    phib = 360.0 * rng.random(len(hkl))
    staged = tmp_path / 'staged.mtz'
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    m.set_cell_for_all(gemmi.UnitCell(50, 50, 50, 90, 90, 90))
    m.add_dataset('work')
    for label, kind in [('FP', 'F'), ('SIGFP', 'Q'), ('PHIB', 'P'), ('FOM', 'W'),
                        ('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A'),
                        ('DM_HOLDOUT', 'I')]:
        m.add_column(label, kind)
    amp = 0.5 * fp
    cols = [fp, np.full(len(hkl), 2.0), phib, np.full(len(hkl), 0.5),
            amp * np.cos(np.radians(phib)), amp * np.sin(np.radians(phib)),
            np.zeros(len(hkl)), np.zeros(len(hkl)), np.zeros(len(hkl))]
    m.set_data(np.column_stack([hkl] + cols).astype(np.float32))
    m.write_to_file(str(staged))

    def fake_fc(source, model_pdb, destination, phenix_bin=None, runtime=None, fcalc='fc_all'):
        src = gemmi.read_mtz_file(str(source))
        data = np.array(src, copy=True)
        labels = src.column_labels()
        extra = [data[:, labels.index('FP')], data[:, labels.index('PHIB')]]
        src.add_column('FC', 'F')
        src.add_column('PHIC', 'P')
        src.set_data(np.column_stack([data, np.column_stack(extra)]).astype(np.float32))
        src.write_to_file(str(destination))

    import phasebuild.fcalc as fcalc_module
    monkeypatch.setattr(fcalc_module, 'attach_fc', fake_fc)
    model = tmp_path / 'model.pdb'
    model.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
    out = tmp_path / 'combined.mtz'
    info = oasis_dm.combine_model_hl(staged, model, out, phenix_bin=None, runtime=None,
                                      fcalc='fc_all')
    # The model term is normalised to the experimental quartet, so a complete
    # model (FC == FP here) may at most double |HL| -- it must never swamp it.
    exp_amp = float(np.hypot(amp * np.cos(np.radians(phib)),
                             amp * np.sin(np.radians(phib))).mean())
    assert 0 < info['mean_added_hl'] <= 1.05 * exp_amp
    result = gemmi.read_mtz_file(str(out))
    labels = result.column_labels()
    assert 'PHIB' in labels and 'FOM' in labels and 'DM_HOLDOUT' in labels
    assert result.nreflections == len(hkl)
