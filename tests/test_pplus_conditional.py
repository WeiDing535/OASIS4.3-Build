import gemmi
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import expit

from oasis_py.data import Reflections
from oasis_py.phasing import iterate_pplus
from oasis_py.clipper_likelihood import reflection_classes
from phasebuild.pplus_conditional import (
    ConditionalBranches, prepare_conditional_branches, iterate_conditional_pplus)


def data(n=2, sg='P 1', hkl=None):
    return Reflections(
        hkl if hkl is not None else [[i+1, 0, 1] for i in range(n)],
        np.ones(n), np.ones(n), gemmi.UnitCell(30, 30, 30, 90, 90, 90),
        gemmi.SpaceGroup(sg))


class ConstantField:
    def __init__(self, values):
        self.values = np.asarray(values, dtype=complex)
        self.nreflections = len(self.values)
        self.calls = 0

    def field(self, phase, fom):
        self.calls += 1
        assert np.isfinite(phase).all() and np.isfinite(fom).all()
        return self.values


def run(d, branches, field=None, sim=None, **kwargs):
    n = len(d.f)
    return iterate_conditional_pplus(
        d, branches, np.zeros(n, complex) if sim is None else sim,
        ConstantField(np.zeros(n) if field is None else field), **kwargs)


def test_delta_limit_recovers_original_pplus_tanh_iteration():
    d = data(3)
    ref = np.array([.2, -.4, .8])
    delta = np.array([.3, 1.2, 2.1])
    first = np.exp(1j*(ref[:, None]+delta[:, None]*np.array([1., -1.])))
    branches = ConditionalBranches(np.full((3, 2), -np.log(2)), first, first**2)
    sim = np.array([.2+.3j, -.7+.1j, .4-.5j])
    field = np.array([.4-.2j, .1+.2j, -.3+.8j])
    expected = iterate_pplus(ref, delta, np.ones(3), sim,
                             ConstantField(field), cycles=5)
    actual = run(d, branches, field, sim, cycles=5)
    np.testing.assert_allclose(actual.probability, expected[2], atol=3e-15)
    np.testing.assert_allclose(actual.first,
                              expected[1]*np.exp(1j*expected[0]), atol=3e-15)
    np.testing.assert_allclose(actual.second,
                              actual.probability*first[:, 0]**2
                              +(1-actual.probability)*first[:, 1]**2, atol=2e-15)


def test_uniform_no_field_has_zero_moments_and_balanced_probability():
    d = data()
    branches = prepare_conditional_branches(d, np.zeros((2, 4)), [.37, .91])
    result = run(d, branches, cycles=7)
    np.testing.assert_allclose(result.probability, .5, atol=2e-15)
    np.testing.assert_allclose(result.first, 0, atol=2e-15)
    np.testing.assert_allclose(result.second, 0, atol=2e-15)
    assert result.metadata['observation_likelihood_accumulated'] is False
    assert 'not exact continuous posterior' in result.metadata['approximation']


def test_first_cycle_is_observation_only_even_with_nonzero_fields():
    d = data()
    hl = np.array([[.7, -.4, .3, .2], [-.2, .8, -.4, .1]])
    branches = prepare_conditional_branches(d, hl, [.2, .6])
    graph = ConstantField([100j, -100j])
    result = iterate_conditional_pplus(d, branches, [20, 30], graph, cycles=1)
    assert graph.calls == 0
    mass = np.exp(branches.log_mass)
    np.testing.assert_allclose(result.first, np.sum(mass*branches.first, axis=1))
    np.testing.assert_allclose(result.second, np.sum(mass*branches.second, axis=1))


def test_finite_shape_update_minimizes_two_state_variational_objective():
    d = data(1)
    b = prepare_conditional_branches(d, [[.6, -.4, .7, -.2]], [.3])
    field = np.array([.9+.8j])
    result = run(d, b, sim=field, cycles=2)
    def objective(p):
        q = np.array([p, 1-p])
        return np.sum(q*(np.log(q)-b.log_mass[0]))-np.real(
            field[0]*np.conj(np.sum(q*b.first[0])))
    p = result.probability[0]
    assert objective(p) < objective(p-.01)
    assert objective(p) < objective(p+.01)
    np.testing.assert_allclose(result.second,
        p*b.second[:, 0]+(1-p)*b.second[:, 1], atol=1e-15)
    assert abs(result.second[0]-result.first[0]**2) > .01


def test_probability_damping_and_triplet_weight_have_explicit_semantics():
    d = data(1)
    b = prepare_conditional_branches(d, [[.5, -.1, .2, .3]], [.2])
    field, sim = np.array([.7j]), np.array([.3])
    p0 = np.exp(b.log_mass[:, 0])
    target = expit(b.log_mass[:, 0]-b.log_mass[:, 1]
                  +np.real((sim+.25*field)*np.conj(b.first[:, 0]-b.first[:, 1])))
    result = run(d, b, field, sim, cycles=2, triplet_weight=.25, damping=.3)
    np.testing.assert_allclose(result.probability, .7*p0+.3*target)


def test_grid_aligned_rotation_and_branch_exchange_preserve_moments():
    d = data()
    hl = np.array([[.8, -.3, .6, -.1], [-.2, .7, -.4, .3]])
    ref = np.array([.31, .63])
    sim = np.array([.4+.7j, -.2+.3j])
    graph = np.array([-.1+.2j, .3-.2j])
    b = prepare_conditional_branches(d, hl, ref)
    original = run(d, b, graph, sim, cycles=4)
    angle = 37*2*np.pi/720
    rotated = hl.copy()
    for harmonic, columns in [(1, (0, 1)), (2, (2, 3))]:
        z = (hl[:, columns[0]]+1j*hl[:, columns[1]])*np.exp(1j*harmonic*angle)
        rotated[:, columns[0]], rotated[:, columns[1]] = z.real, z.imag
    br = prepare_conditional_branches(d, rotated, ref+angle)
    rr = run(d, br, graph*np.exp(1j*angle), sim*np.exp(1j*angle), cycles=4)
    np.testing.assert_allclose(rr.probability, original.probability, atol=2e-14)
    np.testing.assert_allclose(rr.first, original.first*np.exp(1j*angle), atol=2e-14)
    np.testing.assert_allclose(rr.second, original.second*np.exp(2j*angle), atol=2e-14)
    swapped = prepare_conditional_branches(d, hl, ref+np.pi)
    sr = run(d, swapped, graph, sim, cycles=4)
    np.testing.assert_allclose(sr.probability, 1-original.probability, atol=2e-14)
    np.testing.assert_allclose(sr.first, original.first, atol=2e-14)
    np.testing.assert_allclose(sr.second, original.second, atol=2e-14)


def test_conditional_integrals_match_independent_high_precision_quadrature():
    d = data()
    hl = np.array([[.8, -.4, .5, -.3], [-.3, .6, -.2, .4]])
    ref = np.array([0., np.pi/2])
    b = prepare_conditional_branches(d, hl, ref, samples=8192)
    for i in range(2):
        def density(theta):
            return np.exp(hl[i] @ np.array(
                [np.cos(theta), np.sin(theta), np.cos(2*theta), np.sin(2*theta)]))
        total = quad(density, 0, 2*np.pi, epsabs=1e-12, epsrel=1e-12)[0]
        for branch in range(2):
            start = ref[i]+branch*np.pi
            mass = quad(density, start, start+np.pi, epsabs=1e-12)[0]
            np.testing.assert_allclose(np.exp(b.log_mass[i, branch]), mass/total,
                                       atol=5e-8)
            for harmonic, actual in [(1, b.first), (2, b.second)]:
                real = quad(lambda x: density(x)*np.cos(harmonic*x),
                            start, start+np.pi, epsabs=1e-12)[0]/mass
                imag = quad(lambda x: density(x)*np.sin(harmonic*x),
                            start, start+np.pi, epsabs=1e-12)[0]/mass
                np.testing.assert_allclose(actual[i, branch], real+1j*imag, atol=1e-7)


# Independent continuous quadrature replaces the old absolute-grid identity.


def test_centric_states_and_boundary_tie_break_use_allowed_phases():
    d = data(3, 'P 21 21 21', [[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    _, centric, allowed = reflection_classes(d)
    assert centric.all()
    hl = np.array([[.2, .3, .5, .7], [-.4, .6, .2, .1], [.3, -.2, .4, .1]])
    b = prepare_conditional_branches(d, hl, allowed)
    result = run(d, b, sim=[.4j, .5, -.2j], cycles=3)
    np.testing.assert_allclose(np.sin(result.phase-allowed), 0, atol=2e-15)
    np.testing.assert_allclose(result.second, np.exp(2j*allowed), atol=2e-15)
    assert b.metadata['centric_boundary_rows'] == 3
    swapped = prepare_conditional_branches(d, hl, allowed+np.pi)
    sr = run(d, swapped, sim=[.4j, .5, -.2j], cycles=3)
    np.testing.assert_allclose(sr.probability, 1-result.probability, atol=2e-15)
    np.testing.assert_allclose(sr.first, result.first, atol=2e-15)


def test_observations_are_frozen_and_never_accumulate_between_cycles():
    d = data()
    hl = np.array([[.8, -.4, .5, -.3], [-.3, .6, -.2, .4]])
    before = hl.copy()
    ref = np.array([.1, .7])
    b = prepare_conditional_branches(d, hl, ref)
    arrays = [x.copy() for x in [b.log_mass, b.first, b.second]]
    a = run(d, b, sim=[.7j, -.2], cycles=2)
    z = run(d, b, sim=[.7j, -.2], cycles=9)
    np.testing.assert_array_equal(a.probability, z.probability)
    for actual, expected in zip([b.log_mass, b.first, b.second], arrays):
        np.testing.assert_array_equal(actual, expected)
        assert not actual.flags.writeable
    np.testing.assert_array_equal(hl, before)
    np.testing.assert_array_equal(ref, [.1, .7])
    assert np.isfinite(z.probability).all()
    assert np.all((z.probability >= 0) & (z.probability <= 1))
    assert np.all(z.fom <= 1)
    with pytest.raises(ValueError):
        b.first[0, 0] = 0


@pytest.mark.parametrize('samples', [True, 0, 71, 73, 720.])
def test_invalid_angular_grid_is_rejected(samples):
    with pytest.raises(ValueError):
        prepare_conditional_branches(data(1), [[0., 0., 0., 0.]], [0.], samples=samples)


def test_invalid_inputs_and_reflection_identity_are_rejected():
    d = data()
    with pytest.raises(ValueError):
        prepare_conditional_branches(d, np.zeros((2, 3)), [0, 0])
    with pytest.raises(ValueError):
        prepare_conditional_branches(d, np.zeros((2, 4)), [0, np.nan])
    with pytest.raises(ValueError):
        ConditionalBranches([[0, 0]], [[2, 0]], [[1, 0]])
    with pytest.raises(ValueError):
        ConditionalBranches([[-np.inf, -np.inf]], [[0, 0]], [[0, 0]])
    with pytest.raises(ValueError):
        ConditionalBranches([[0, 0]], [[1, 0]], [[0, 0]])
    b = prepare_conditional_branches(d, np.zeros((2, 4)), [0, 0])
    mismatch = data(2, hkl=[[2, 0, 1], [3, 0, 1]])
    with pytest.raises(ValueError, match='identity'):
        run(mismatch, b)
    for kwargs in [{'cycles': 0}, {'damping': 0}, {'damping': 1.1},
                   {'triplet_weight': -1}, {'triplet_weight': np.nan}]:
        with pytest.raises(ValueError):
            run(d, b, **kwargs)
    with pytest.raises(ValueError):
        run(d, b, sim=[np.nan, 0])
    with pytest.raises(ValueError):
        run(d, b, field=[np.inf, 0], cycles=2)
    with pytest.raises(ValueError):
        run(d, b, field=[0], cycles=2)

