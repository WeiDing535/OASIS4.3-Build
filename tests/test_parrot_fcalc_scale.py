"""Regression tests for the trace-feedback Fcalc column choice and amplitude scaling.

The default ``fc`` mode undoes REFMAC's overall B factor so the partial-model
amplitude enters Parrot at the conventional Fcalc scale. Without this, the
amplitude is ~1.66x too large, the trace collapses, and R-free degrades by ~17 pp
on 7E1D. These tests pin both the log parser and the amplitude formula.
"""
from pathlib import Path

import gemmi
import numpy as np
import pytest

from phasebuild import parrot_dm


REFMAC_LOG = """Overall scale should be multiplied by       1.0000
Overall               : scale =     1.5629, B  =  -13.5775
Babinet's bulk solvent: scale =     0.5019, B  =  160.5187
"""


def test_overall_b_parses_refmac_log(tmp_path):
    log = tmp_path / 'refmac.log'
    log.write_text(REFMAC_LOG)
    assert parrot_dm._overall_b(log) == pytest.approx(-13.5775)


def test_overall_b_missing_or_malformed_returns_none(tmp_path):
    log = tmp_path / 'refmac.log'
    log.write_text('no scaling data here\n')
    assert parrot_dm._overall_b(log) is None
    assert parrot_dm._overall_b(tmp_path / 'absent.log') is None


def test_b_amplitude_factor_matches_fmodel_scale():
    """exp(B*s^2/4) must reproduce the fmodel/REFMAC ratio per shell.

    Measured on 7E1D cycle_02 (fmodel_fcalc_matched.mtz vs refmac fcalc.mtz):
    shell 2.0-2.3 A ratio 0.482 (factor 0.484), 4.0-6.0 A ratio 0.857
    (factor 0.854), median after/f 0.996, Pearson 0.99995.
    """
    B = -13.5775
    for d, expected in [(2.15, 0.482), (2.5, 0.577), (3.0, 0.675),
                        (3.6, 0.767), (5.0, 0.857)]:
        factor = np.exp(B / (4.0 * d * d))
        assert factor == pytest.approx(expected, abs=0.02), (d, factor)
    # Monotonic and always <= 1 for negative B (high resolution shrinks most).
    ds = np.array([2.0, 3.0, 5.0, 10.0])
    factors = np.exp(B / (4.0 * ds * ds))
    assert np.all(factors <= 1.0)
    assert np.all(np.diff(factors) > 0)


@pytest.mark.parametrize('mode', ['fc', 'fc_all', 'refined'])
def test_attach_fc_selects_requested_columns(tmp_path, monkeypatch, mode):
    """The fcalc option must actually choose REFMAC columns (it used to be dead)."""
    from phasebuild import parrot_dm

    staged = gemmi.Mtz(with_base=True)
    staged.spacegroup = gemmi.SpaceGroup('P 21 21 21')
    staged.cell = gemmi.UnitCell(35.8, 52.8, 97.4, 90, 90, 90)
    for col in ['FP', 'SIGFP']:
        staged.add_column(col, 'F' if col == 'FP' else 'Q')
    staged.set_data(np.array([[1, 0, 0, 20.0, 2.0], [2, 0, 0, 30.0, 2.0]], np.float32))
    staged.write_to_file(str(tmp_path / 'staged.mtz'))

    _fixture_model(tmp_path)

    def fake_run(command, **kwargs):
        out = Path(command[command.index('HKLOUT') + 1])
        m = gemmi.Mtz(with_base=True)
        m.spacegroup = gemmi.SpaceGroup('P 21 21 21')
        m.cell = gemmi.UnitCell(35.8, 52.8, 97.4, 90, 90, 90)
        for label, typ in [('FP', 'F'), ('SIGFP', 'Q'),
                           ('FC', 'F'), ('PHIC', 'P'),
                           ('FC_ALL', 'F'), ('PHIC_ALL', 'P'),
                           ('FC_ALL_LS', 'F'), ('PHIC_ALL_LS', 'P')]:
            m.add_column(label, typ)
        vals = {
            'FC': [50.0, 60.0], 'PHIC': [10.0, 20.0],
            'FC_ALL': [40.0, 48.0], 'PHIC_ALL': [11.0, 21.0],
            'FC_ALL_LS': [30.0, 36.0], 'PHIC_ALL_LS': [12.0, 22.0],
        }
        m.set_data(np.array([[1, 0, 0, 20, 2] + [vals[c][0] for c in
                             ['FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', 'FC_ALL_LS', 'PHIC_ALL_LS']],
                            [2, 0, 0, 30, 2] + [vals[c][1] for c in
                             ['FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', 'FC_ALL_LS', 'PHIC_ALL_LS']]],
                           np.float32))
        m.write_to_file(str(out))
        (Path(command[command.index('HKLOUT') + 1]).parent / 'refmac.log').write_text(REFMAC_LOG)
        return type('R', (), {'returncode': 0})()

    monkeypatch.setattr(parrot_dm.subprocess, 'run', fake_run)
    monkeypatch.setattr(parrot_dm, '_refmac_binary', lambda runtime=None: '/usr/bin/true')
    monkeypatch.setattr('phasebuild.runtime.runtime_environment', lambda *a, **k: {})

    out = tmp_path / 'with_model.mtz'
    info = parrot_dm._attach_fc(tmp_path / 'staged.mtz',
                                tmp_path / 'model.pdb', out, fcalc=mode)
    result = gemmi.read_mtz_file(str(out))
    fc = np.asarray(result.column_with_label('FC'))
    ph = np.asarray(result.column_with_label('PHIC'))
    expected = {'fc': ([50.0, 60.0], [10.0, 20.0]),
                'fc_all': ([40.0, 48.0], [11.0, 21.0]),
                'refined': ([30.0, 36.0], [12.0, 22.0])}[mode]
    assert info['fcalc'] == mode
    assert info['columns'] == {'fc': ['FC', 'PHIC'], 'fc_all': ['FC_ALL', 'PHIC_ALL'],
                               'refined': ['FC_ALL_LS', 'PHIC_ALL_LS']}[mode]
    assert ph == pytest.approx(expected[1], abs=1e-3)  # phase never changes
    assert np.all(fc <= expected[0])  # fc mode rescales down; others unchanged
    if mode == 'fc':
        assert info['amplitude_normalised'] is True
        # B<0: high-res (smaller d for h=2 here) shrinks more than h=1
        assert fc[1] / expected[0][1] < fc[0] / expected[0][0]


def _fixture_model(tmp_path):
    path = tmp_path / 'model.pdb'
    path.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n')
    return path


def test_invalid_fcalc_rejected(tmp_path):
    with pytest.raises(ValueError, match='fcalc must be'):
        parrot_dm._attach_fc(tmp_path / 'staged.mtz', tmp_path / 'model.pdb',
                             tmp_path / 'out.mtz', fcalc='bogus')


_7E1D_BASE = Path('/Volumes/Fast_4T/Workbuddy_Oasis/OASIS4.3-20261117/work/reproduce/7E1D')
@pytest.mark.skipif(not (_7E1D_BASE / 'fmodel_compare/cycle_02/fmodel_fcalc_matched.mtz').is_file(),
                    reason='7E1D fmodel reference MTZ not available (work/ is gitignored)')
def test_real_7e1d_b_factor_reproduces_fmodel_amplitude():
    """End-to-end sanity check against the measured 7E1D reference data."""
    base = Path('/Volumes/Fast_4T/Workbuddy_Oasis/OASIS4.3-20261117/work/reproduce/7E1D')
    fmodel = gemmi.read_mtz_file(str(base / 'fmodel_compare/cycle_02/fmodel_fcalc_matched.mtz'))
    refmac = gemmi.read_mtz_file(str(base / 'endpoint_fb3_s0.50/03_dm/cycle_02/refmac_fc/fcalc.mtz'))
    f = np.asarray(fmodel.column_with_label('FC'))
    r = np.asarray(refmac.column_with_label('FC'))
    hkl = np.array([[h, k, l] for h, k, l, *_ in np.array(fmodel)], np.int32)
    d = fmodel.cell.calculate_d_array(hkl)
    factor = np.exp(parrot_dm._overall_b(
        base / 'endpoint_fb3_s0.50/03_dm/cycle_02/refmac_fc/refmac.log') / (4.0 * d * d))
    mask = (f > 1) & (r > 1)
    assert np.corrcoef(r[mask] * factor[mask], f[mask])[0, 1] > 0.999
    assert np.median((r[mask] * factor[mask]) / f[mask]) == pytest.approx(1.0, abs=0.02)
