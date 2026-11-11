"""Observation-only scale recovery and fixed-likelihood seed contracts."""
from dataclasses import replace
import numpy as np
import gemmi
import pytest
from oasis_py.data import Reflections, Atoms, structure_factors
from oasis_py.phasing import iterate_pplus
from phasebuild.pplus_experimental import prepare_pplus_calibration
from phasebuild.pplus_observed_scale import phase_pplus_observed_scale


OPTIONS = dict(max_neighbors=30, max_per_target=8)
COMPOSITION = {'C': 200., 'N': 50., 'O': 60., 'Se': 1.}
FPP = {'Se': 6.}


def synthetic(alpha=.23, with_centric=False, placeholder=1.):
    # E[DANO^2] = sigmaD^2 + 2 alpha^2 |Q|^2 exactly per row, so the
    # marginal-variance optimum is known without protein phase information.
    hkl = np.array([[h, k, 1] for h in range(1, 9) for k in range(1, 5)])
    if with_centric:
        hkl = np.vstack([hkl, [1, 0, 1], [2, 0, 1]])
    n = len(hkl)
    data = Reflections(hkl, 20. + np.arange(n) / 2, np.ones(n),
        gemmi.UnitCell(30, 32, 34, 90, 90, 90),
        gemmi.SpaceGroup('P 21' if with_centric else 'P 1'),
        dano=np.ones(n), sigdano=np.linspace(.4, 1.2, n))
    atoms = Atoms([[.13, .17, .19]], ['Se'], [1.], [12.])
    raw = structure_factors(data, atoms, FPP)
    signs = np.where(np.arange(n) % 2, -1, 1)
    data.dano[:] = signs * np.sqrt(data.sigdano**2 + 2 * alpha**2 * abs(raw)**2)
    single = data.centric.copy()
    data.dano[single] = 0.
    data.sigdano[single] = placeholder
    context = prepare_pplus_calibration(data, COMPOSITION, atoms, FPP,
        centric_singletons=single, **OPTIONS)
    return data, atoms, single, context, raw


def run(data, atoms, single, context):
    return phase_pplus_observed_scale(data, COMPOSITION, atoms, FPP,
        calibration_context=context, centric_singletons=single, **OPTIONS)


def test_known_anomalous_amplitude_scale_is_recovered_without_true_phases():
    data, atoms, single, context, raw = synthetic()
    result = run(data, atoms, single, context)
    alpha = result.metadata['geometry_anomalous_scale']
    assert alpha == pytest.approx(.23, rel=1e-5)
    difference = np.arccos(np.clip(data.dano / (2 * abs(alpha * raw)), -1, 1))
    expected = iterate_pplus(np.angle(raw), difference, context.error,
                             context.sim, context.graph, cycles=5)
    np.testing.assert_allclose(result.phase, expected[0], rtol=0, atol=1e-14)
    np.testing.assert_allclose(result.fom, expected[1], rtol=0, atol=1e-14)
    assert result.metadata['seed_only_geometry_ablation']
    assert not result.metadata['fully_recalibrated_likelihood']
    assert result.metadata['cosine_rank'] is False


def test_null_anomalous_signal_is_rejected_instead_of_inventing_scale():
    data, atoms, single, _, _ = synthetic()
    data.dano[:] = 0.
    context = prepare_pplus_calibration(data, COMPOSITION, atoms, FPP,
        centric_singletons=single, **OPTIONS)
    with pytest.raises(ValueError, match='No positive observed anomalous scale'):
        run(data, atoms, single, context)


def test_frozen_likelihood_fields_graph_and_old_observations_are_unchanged():
    data, atoms, single, context, _ = synthetic()
    arrays = [context.reference, context.difference, context.error,
              context.observations, context.sim, context.graph.kappa]
    before = [a.tobytes() for a in arrays]
    original_data = [a.tobytes() for a in [data.f, data.sigf, data.dano, data.sigdano]]
    result = run(data, atoms, single, context)
    assert [a.tobytes() for a in arrays] == before
    assert [a.tobytes() for a in [data.f, data.sigf, data.dano, data.sigdano]] == original_data
    assert result.hl.tobytes() == context.observations.tobytes()
    assert not np.shares_memory(result.hl, context.observations)
    result.hl[0, 0] += 1
    assert context.observations.tobytes() == before[3]
    assert result.metadata['observation_HL_unchanged']


def test_centric_singleton_sigma_placeholder_cannot_affect_scale_or_seed():
    a = synthetic(with_centric=True, placeholder=1.)
    b = synthetic(with_centric=True, placeholder=1e12)
    assert a[2].sum() == 2
    first, second = run(*a[:4]), run(*b[:4])
    for key,value in first.metadata['anomalous_scale_fit'].items():
        other=second.metadata['anomalous_scale_fit'][key]
        if isinstance(value,float):
            assert value == pytest.approx(other,rel=1e-12,abs=1e-11)
        elif isinstance(value,list):
            np.testing.assert_allclose(value,other,rtol=1e-12,atol=1e-11)
        else:
            assert value == other
    # Independent direct summations and Wilson fits may differ by a few ULPs.
    # The separate frozen-context test checks unchanged HL byte-for-byte;
    # this repeated-fit comparison uses a numerical tolerance instead.
    np.testing.assert_allclose(first.fom * np.exp(1j * first.phase),
                               second.fom * np.exp(1j * second.phase),
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(first.hl, second.hl, rtol=1e-12, atol=1e-11)
    np.testing.assert_array_equal(first.hl[a[2]], 0.)
    assert first.metadata['anomalous_scale_fit']['acentric_observations'] == 32


def test_changed_inputs_and_mismatched_calibration_are_rejected_before_fit(monkeypatch):
    import phasebuild.pplus_observed_scale as module
    data, atoms, single, context, _ = synthetic()
    def forbidden(*args, **kwargs):
        raise AssertionError('Mismatched inputs reached scale fitting')
    monkeypatch.setattr(module, 'fit_anomalous_scale', forbidden)
    with pytest.raises(ValueError, match='does not match'):
        run(data, atoms, single, replace(context, identity='incorrect'))
    data.dano[0] *= 1.01
    with pytest.raises(ValueError, match='does not match'):
        run(data, atoms, single, context)
    with pytest.raises(ValueError, match='frozen PplusCalibration'):
        run(data, atoms, single, None)
