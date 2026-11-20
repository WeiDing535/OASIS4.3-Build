"""In-house engines: handdm density modification and minirefmac refinement.

These cover the configuration surface and the MTZ column contract of the
in-house adapters. The engines themselves are exercised by their own
deployments; the pipeline-facing contract is what is asserted here.
"""
import json
from pathlib import Path

import gemmi
import numpy as np
import pytest

from phasebuild import handdm_dm, minirefmac_refine, refinement
from phasebuild.sad_pipeline import DEFAULTS, validate, stage_plan


BASE = dict(method='pplus_trust', build_model=False, fpp='Se=3.8',
            sad='/nonexistent/sad.mtz', native='/nonexistent/native.mtz',
            sequence='/nonexistent/seq.fasta', heavy='/nonexistent/heavy.pdb',
            output='/nonexistent/out')


def _validate(tmp_path, include_prior=True, **overrides):
    """Validate a configuration with every engine probe stubbed out.

    ``include_prior=False`` leaves the external phase file out of the
    configuration, which the bundled method rejects by design.
    """
    import phasebuild.parrot_dm as pd
    import phasebuild.handdm_dm as hd
    import phasebuild.minirefmac_refine as mr
    inputs = {}
    names = ['sad', 'native', 'sequence', 'heavy', 'reference']
    if include_prior:
        names.append('observation_prior')
    for name in names:
        path = tmp_path / name
        path.write_text('placeholder')
        inputs[name] = str(path)
    cfg = {**BASE, **inputs, **overrides}
    originals = (pd.locate, hd.locate, mr.locate)
    pd.locate = lambda *a, **k: '/usr/bin/true'
    hd.locate = lambda *a, **k: '/usr/bin/true'
    mr.locate = lambda *a, **k: '/opt/minirefmac'
    try:
        return validate(cfg, inspect=False)
    finally:
        pd.locate, hd.locate, mr.locate = originals


# --------------------------------------------------------------- configuration

def test_defaults_still_use_the_external_engines():
    assert DEFAULTS['dm_engine'] == 'parrot'
    assert DEFAULTS['refine_engine'] == 'refmac'


def test_handdm_is_an_accepted_density_modification_engine(tmp_path):
    cfg = _validate(tmp_path, dm_engine='handdm')
    assert cfg['dm_engine'] == 'handdm'


def test_minirefmac_is_an_accepted_refinement_engine(tmp_path):
    cfg = _validate(tmp_path, refine_engine='minirefmac')
    assert cfg['refine_engine'] == 'minirefmac'


@pytest.mark.parametrize('engine', ['parrot2', 'resolve', 'cparrot', ''])
def test_unknown_density_modification_engine_is_rejected(tmp_path, engine):
    with pytest.raises(ValueError):
        _validate(tmp_path, dm_engine=engine)


@pytest.mark.parametrize('engine', ['refmac5', 'refine', ''])
def test_unknown_refinement_engine_is_rejected(tmp_path, engine):
    with pytest.raises(ValueError):
        _validate(tmp_path, refine_engine=engine)


def test_inhouse_engines_are_selected_explicitly_not_by_default(tmp_path):
    """The in-house engines must never take over a run silently."""
    cfg = _validate(tmp_path)
    assert cfg['dm_engine'] == 'parrot' and cfg['refine_engine'] == 'refmac'


# ------------------------------------------------------------------ dispatch

def test_refine_model_rejects_an_unknown_engine():
    with pytest.raises(ValueError, match='Unknown refinement engine'):
        refinement.refine_model('/nonexistent.mtz', '/nonexistent.pdb', '/tmp/x',
                                engine='nonsense')


def test_refine_model_routes_minirefmac_to_the_inhouse_adapter(monkeypatch):
    seen = {}

    def fake_refine(mtz_path, model_path, output, cycles=10, root=None, **kwargs):
        seen.update(mtz=mtz_path, model=model_path, output=output, cycles=cycles, root=root)
        return {'status': 'success'}

    monkeypatch.setattr(minirefmac_refine, 'refine', fake_refine)
    result = refinement.refine_model('/in.mtz', '/in.pdb', '/out', cycles=3,
                                     engine='minirefmac', minirefmac_root='/opt/mr')
    assert result['status'] == 'success'
    assert seen == dict(mtz='/in.mtz', model='/in.pdb', output='/out', cycles=3,
                        root='/opt/mr')


# ------------------------------------------------------------ column contract

def _write_mtz(path, columns, rows):
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = gemmi.SpaceGroup('P 1 21 1')
    m.set_cell_for_all(gemmi.UnitCell(30, 40, 50, 90, 90, 90))
    m.add_dataset('test')
    for label, kind in columns:
        m.add_column(label, kind)
    m.set_data(np.asarray(rows, np.float32))
    m.write_to_file(str(path))
    return m


def test_handdm_staging_merges_amplitudes_with_the_hl_quartet(tmp_path):
    work = tmp_path / 'work.mtz'
    prior = tmp_path / 'prior.mtz'
    rows = [[1, 0, 0, 100.0, 10.0, 0], [2, 0, 0, 200.0, 20.0, 0]]
    _write_mtz(work, [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')], rows)
    hl_rows = [[1, 0, 0, 10.0, 0.0, -1.0, 0.0], [2, 0, 0, 5.0, 0.0, -1.0, 0.0]]
    _write_mtz(prior, [('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A')], hl_rows)
    staged = tmp_path / 'staged.mtz'
    count, quartet = handdm_dm._stage_input(work, prior, staged)
    assert count == 2 and quartet == ['HLA', 'HLB', 'HLC', 'HLD']
    m = gemmi.read_mtz_file(str(staged))
    labels = m.column_labels()
    for label in ('FOBS', 'SIGFOBS', 'PHIB', 'FOM', 'HLA', 'HLB', 'HLC', 'HLD', 'FreeR'):
        assert label in labels
    data = np.array(m)
    assert data[0, labels.index('FOBS')] == 100.0
    assert data[0, labels.index('PHIB')] == pytest.approx(0.0, abs=1e-6)


def test_handdm_remap_writes_the_pipeline_schema(tmp_path):
    raw = tmp_path / 'raw.mtz'
    # H K L | FOBS SIGFOBS FreeR PHIB FOM PHIC | HLA HLB HLC HLD
    #        | DM_FDM DM_PHIDM DM_FOMDM | DM_HLADM..DM_HLDDM
    #        | DM_FWT DM_PHWT | DM_DELFWT DM_PHDELFWT | DM_FCDM DM_PHICDM
    rows = [[1, 0, 0,
             100.0, 10.0, 0, 30.0, 0.8, 12.0,
             11.0, 10.0, 9.0, 0.0,
             55.0, 45.0, 0.75,
             3.0, 2.0, 1.0, 0.0,
             60.0, 40.0,
             1.0, 50.0,
             35.0, 2.0]]
    columns = [('FOBS', 'F'), ('SIGFOBS', 'Q'), ('FreeR', 'I'), ('PHIB', 'P'),
               ('FOM', 'W'), ('PHIC', 'P'),
               ('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A'),
               ('DM_FDM', 'F'), ('DM_PHIDM', 'P'), ('DM_FOMDM', 'W'),
               ('DM_HLADM', 'A'), ('DM_HLBDM', 'A'), ('DM_HLCDM', 'A'), ('DM_HLDDM', 'A'),
               ('DM_FWT', 'F'), ('DM_PHWT', 'P'), ('DM_DELFWT', 'F'), ('DM_PHDELFWT', 'P'),
               ('DM_FCDM', 'F'), ('DM_PHICDM', 'P')]
    _write_mtz(raw, columns, rows)
    out = tmp_path / 'dm.mtz'
    info = handdm_dm._remap(raw, out)
    m = gemmi.read_mtz_file(str(out))
    labels = m.column_labels()
    for label in ('F', 'SIGF', 'FWT', 'PHWT', 'PHIB', 'FOM',
                  'HLA', 'HLB', 'HLC', 'HLD'):
        assert label in labels
    data = np.array(m)
    assert data[0, labels.index('F')] == pytest.approx(100.0)
    assert data[0, labels.index('FWT')] == pytest.approx(60.0)
    assert data[0, labels.index('PHWT')] == pytest.approx(40.0)
    assert data[0, labels.index('PHIB')] == pytest.approx(45.0)
    assert data[0, labels.index('FOM')] == pytest.approx(0.75)
    assert data[0, labels.index('HLA')] == pytest.approx(3.0)
    assert info['mean_fom'] == pytest.approx(0.75)


def test_handdm_reports_a_missing_engine_clearly():
    with pytest.raises(ValueError, match='handdm'):
        handdm_dm.locate('/nonexistent/handdm_v2.py')


def test_minirefmac_reports_a_missing_deployment_clearly():
    with pytest.raises(ValueError, match='minirefmac'):
        minirefmac_refine.locate('/nonexistent/root')


def test_stage_plan_passes_the_refinement_engine_to_the_build_stage(tmp_path):
    import shutil
    (tmp_path / 'reference').write_text('placeholder')
    (tmp_path / 'heavy').write_text('placeholder')
    runtime = tmp_path / 'runtime'
    (runtime / 'bin').mkdir(parents=True)
    for name in ('phaser', 'cparrot', 'cbuccaneer', 'refmac5'):
        path = runtime / 'bin' / name
        path.write_text('#!/bin/sh\nexit 0\n')
        path.chmod(0o755)
    library = runtime / 'lib/data/reference_structures'
    library.mkdir(parents=True)
    shutil.copyfile(tmp_path / 'reference', library / 'reference-1tqw.mtz')
    shutil.copyfile(tmp_path / 'heavy', library / 'reference-1tqw.pdb')
    cfg = _validate(tmp_path, method='pplus_bundled', build_model=True,
                    include_prior=False, runtime=str(runtime), input_columns=False,
                    output=str(tmp_path / 'out'), refine_engine='minirefmac')
    stages = stage_plan(cfg)
    build = [s for s in stages if s['id'] == 'build'][0]
    command = ' '.join(build['command'])
    assert '--refine-engine minirefmac' in command


# ------------------------------------------------- protocol-difference guards

def test_handdm_refuses_the_trace_feedback_loop(tmp_path):
    """The in-house engine cannot combine a partial model with its map."""
    with pytest.raises(ValueError, match='trace feedback'):
        _validate(tmp_path, dm_engine='handdm', feedback_cycles=2)


def test_parrot_keeps_the_trace_feedback_loop(tmp_path):
    cfg = _validate(tmp_path, feedback_cycles=2)
    assert cfg['feedback_cycles'] == 2


def test_build_entry_refuses_per_round_cleanup_with_the_inhouse_refiner():
    import run_build_from_phases as entry
    with pytest.raises(SystemExit):
        entry.main(['--map', 'm.mtz', '--native', 'n.mtz', '--sequence', 's.fasta',
                    '--phase-prior', 'p.mtz', '--output', 'o',
                    '--refine-engine', 'minirefmac',
                    '--refinement-prior', 'dm_phifom', '--per-round-cleanup'])


# ------------------------------------------------------- building recipe path

def test_pipeline_recipe_passes_the_refinement_engine_to_the_builder(tmp_path):
    import shutil
    from phasebuild import ccp4_pipeline
    (tmp_path / 'reference').write_text('placeholder')
    (tmp_path / 'heavy').write_text('placeholder')
    runtime = tmp_path / 'runtime'
    (runtime / 'bin').mkdir(parents=True)
    for name in ('phaser', 'cparrot', 'cbuccaneer', 'refmac5'):
        path = runtime / 'bin' / name
        path.write_text('#!/bin/sh\nexit 0\n')
        path.chmod(0o755)
    library = runtime / 'lib/data/reference_structures'
    library.mkdir(parents=True)
    shutil.copyfile(tmp_path / 'reference', library / 'reference-1tqw.mtz')
    shutil.copyfile(tmp_path / 'heavy', library / 'reference-1tqw.pdb')
    cfg = _validate(tmp_path, method='pplus_bundled', build_model=True,
                    include_prior=False, runtime=str(runtime), input_columns=False,
                    building_recipe='pipeline', output=str(tmp_path / 'out'),
                    refine_engine='minirefmac')
    stages = stage_plan(cfg)
    build = [s for s in stages if s['id'] == 'build'][0]
    command = ' '.join(build['command'])
    assert '--refine-engine minirefmac' in command
    assert 'run_build_pipeline.py' in command


def test_pipeline_recipe_keeps_the_external_refiner_by_default(tmp_path):
    import shutil
    (tmp_path / 'reference').write_text('placeholder')
    (tmp_path / 'heavy').write_text('placeholder')
    runtime = tmp_path / 'runtime'
    (runtime / 'bin').mkdir(parents=True)
    for name in ('phaser', 'cparrot', 'cbuccaneer', 'refmac5'):
        path = runtime / 'bin' / name
        path.write_text('#!/bin/sh\nexit 0\n')
        path.chmod(0o755)
    library = runtime / 'lib/data/reference_structures'
    library.mkdir(parents=True)
    shutil.copyfile(tmp_path / 'reference', library / 'reference-1tqw.mtz')
    shutil.copyfile(tmp_path / 'heavy', library / 'reference-1tqw.pdb')
    cfg = _validate(tmp_path, method='pplus_bundled', build_model=True,
                    include_prior=False, runtime=str(runtime), input_columns=False,
                    building_recipe='pipeline', output=str(tmp_path / 'out'))
    stages = stage_plan(cfg)
    build = [s for s in stages if s['id'] == 'build'][0]
    assert '--refine-engine' not in ' '.join(build['command'])


def test_inhouse_builder_keeps_the_combined_hl_columns_the_adapter_derives():
    """Buccaneer's correlation mode needs the combined phase quartet.

    The in-house engine does not emit PHCOMB/HLACOMB..HLDCOMB itself, so the
    refinement adapter derives them before the next build round. The builder
    command therefore keeps the same ``colin-hl`` line as the REFMAC route;
    dropping it makes Buccaneer abort while building the map grid.
    """
    from phasebuild import ccp4_pipeline
    flags = {'title': 't', 'pdbin_ref': 'a.pdb', 'mtzin_ref': 'a.mtz',
             'seqin': 's.fasta', 'resolution': '2.000'}
    text = ccp4_pipeline.BUCCANEER_NTH.format(
        mtzin='m.mtz', pdbin='p.pdb', pdbout='o.pdb', **flags)
    assert 'colin-hl HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB' in text
    assert 'colin-fc FWT,PHWT' in text


def test_inhouse_refinement_derives_the_combined_phase_quartet(tmp_path):
    """The adapter adds PHCOMB and the HL quartet to the refined MTZ."""
    import numpy as np
    from phasebuild import minirefmac_refine
    source = tmp_path / 'refined.mtz'
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = gemmi.find_spacegroup_by_name('P 21 21 21')
    m.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
    m.add_dataset('data')
    for label, kind in [('FP', 'F'), ('SIGFP', 'Q'), ('FreeR_flag', 'I'),
                        ('PHWT', 'P'), ('FOM', 'W')]:
        m.add_column(label, kind)
    rows = np.zeros((3, 8), dtype=np.float32)
    rows[:, 3] = [100.0, 200.0, 300.0]
    rows[:, 4] = [10.0, 20.0, 30.0]
    rows[:, 5] = [0, 1, 0]
    rows[:, 6] = [45.0, 200.0, 310.0]
    rows[:, 7] = [0.8, 0.5, 0.2]
    m.set_data(rows)
    m.write_to_file(str(source))
    destination = tmp_path / 'with_quartet.mtz'
    assert minirefmac_refine._hl_quartet(source, destination) is True
    out = gemmi.read_mtz_file(str(destination))
    labels = out.column_labels()
    for name in ('PHCOMB', 'HLACOMB', 'HLBCOMB', 'HLCCOMB', 'HLDCOMB'):
        assert name in labels
    data = np.array(out)
    kappa = np.hypot(data[:, labels.index('HLACOMB')],
                     data[:, labels.index('HLBCOMB')])
    # The amplitude is the figure of merit scaled to the REFMAC range.
    assert np.allclose(kappa, np.clip(data[:, labels.index('FOM')] * 81.0, 0, 81))
    # The quartet phase is the refined phase.
    phase = np.degrees(np.arctan2(data[:, labels.index('HLBCOMB')],
                                   data[:, labels.index('HLACOMB')])) % 360.0
    expected = data[:, labels.index('PHWT')]
    delta = np.abs(((phase - expected + 180) % 360) - 180)
    assert np.allclose(delta, 0, atol=1e-3)
    # A second call is a no-op: the columns are already present.
    assert minirefmac_refine._hl_quartet(destination, tmp_path / 'again.mtz') is False


def _capture_handdm_command(tmp_path, monkeypatch, free_rows):
    """Run the adapter with the subprocess call stubbed out and return the command."""
    import phasebuild.handdm_dm as hd
    work = tmp_path / 'work.mtz'
    prior = tmp_path / 'prior.mtz'
    rows = [[1, 0, 0, 100.0, 10.0, 0], [2, 0, 0, 200.0, 20.0, 0]]
    _write_mtz(work, [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')], rows)
    hl_rows = [[1, 0, 0, 10.0, 0.0, -1.0, 0.0], [2, 0, 0, 5.0, 0.0, -1.0, 0.0]]
    _write_mtz(prior, [('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A')], hl_rows)
    captured = {}

    class _Completed:
        returncode = 0

    def fake_run(command, **kwargs):
        captured['command'] = list(command)
        raw = tmp_path / 'handdm_raw.mtz'
        _write_mtz(raw, [('FOBS', 'F'), ('SIGFOBS', 'Q'), ('DM_FWT', 'F'),
                          ('DM_PHWT', 'P'), ('DM_PHIDM', 'P'), ('DM_FOMDM', 'W'),
                          ('DM_HLADM', 'A'), ('DM_HLBDM', 'A'),
                          ('DM_HLCDM', 'A'), ('DM_HLDDM', 'A')],
                   [[1, 0, 0, 100.0, 10.0, 60.0, 40.0, 45.0, 0.75, 3.0, 2.0, 1.0, 0.0]])
        return _Completed()

    monkeypatch.setattr(hd, 'locate', lambda *a, **k: '/opt/handdm_v2.py')
    monkeypatch.setattr(hd, '_hist_lib', lambda *a, **k: '/opt/hist.lib')
    monkeypatch.setattr(hd.subprocess, 'run', fake_run)
    monkeypatch.setattr(hd, '_free_rows', lambda *a, **k: free_rows)
    hd.run(work, prior, tmp_path / 'out.mtz', 0.5, cycles=40)
    return captured['command']


def test_work_only_input_disables_the_engine_early_stop(tmp_path, monkeypatch):
    """Without free reflections the early-stop criterion cannot be evaluated."""
    command = _capture_handdm_command(tmp_path, monkeypatch, free_rows=0)
    assert '--no-auto-stop' in command
    assert '--cycles' in command and '40' in command


def test_free_set_input_keeps_the_engine_early_stop(tmp_path, monkeypatch):
    command = _capture_handdm_command(tmp_path, monkeypatch, free_rows=300)
    assert '--no-auto-stop' not in command
