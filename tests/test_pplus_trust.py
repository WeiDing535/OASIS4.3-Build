import gemmi
import numpy as np
import pytest
from oasis_py.data import Reflections
from oasis_py.phasing import iterate_pplus
from phasebuild.pplus_conditional import ConditionalBranches, prepare_conditional_branches
from phasebuild.pplus_trust import bernoulli_kl, bounded_branch_mass, iterate_trusted_pplus


def data(hkl=None):
    return Reflections([[1, 0, 1], [2, 0, 1]] if hkl is None else hkl,
                       np.ones(2), np.ones(2), gemmi.UnitCell(30, 30, 30, 90, 90, 90),
                       gemmi.SpaceGroup('P 1'))


class Field:
    nreflections = 2

    def field(self, phase, fom):
        return np.full(2, 1j)


@pytest.mark.parametrize('budget', [0., .0005, .0025, .01, .04])
def test_kl_bound_and_update_direction(budget):
    p = np.array([.01, .2, .5, .8, .99])
    lm = np.log(np.column_stack([p, 1-p]))
    score = np.array([100, -100, 100, -1e6, 1e6])
    q = bounded_branch_mass(lm, score, budget)
    assert np.max(bernoulli_kl(q, lm)) <= budget + 1e-12
    assert np.all((q-p)*score >= -1e-14)
    if budget == 0:
        np.testing.assert_allclose(q, p, atol=1e-15)


def test_branch_exchange_and_exact_zero_mass():
    p = np.array([.2, .5, .9])
    lm = np.log(np.column_stack([p, 1-p]))
    s = np.array([2, -3, 10])
    np.testing.assert_allclose(bounded_branch_mass(lm, s, .01),
                               1-bounded_branch_mass(lm[:, ::-1], -s, .01), atol=1e-13)
    np.testing.assert_array_equal(bounded_branch_mass([[0, -np.inf], [-np.inf, 0]],
                                                       [-1000, 1000], .01), [1, 0])


def test_real_branch_update_converges_with_valid_conditional_moments():
    branches = ConditionalBranches(np.log(np.full((2, 2), .5)),
                                   np.array([[1j, -1j], [1j, -1j]]),
                                   np.full((2, 2), -1+0j))
    result = iterate_trusted_pplus(data(), branches, Field(), kl_budget=.01)
    assert result.metadata['converged']
    assert np.all(result.probability > .55)
    assert np.max(bernoulli_kl(result.probability, branches.log_mass)) <= .01 + 1e-12
    np.testing.assert_allclose(result.first, 1j*(2*result.probability-1))
    np.testing.assert_allclose(result.fom, abs(result.first))
    assert np.all(abs(result.second-result.first**2) <= 1-abs(result.first)**2+1e-12)
    assert result.metadata['observation_likelihood_accumulated'] is False
    assert result.metadata['posterior_HL_exported'] is False


def test_loose_trust_bound_recovers_legacy_Pplus_in_doublet_limit():
    reference = np.array([.3, -.5])
    difference = np.array([.6, 1.2])
    first = np.exp(1j*(reference[:, None] + difference[:, None]*np.array([1, -1])))
    branches = ConditionalBranches(np.full((2, 2), -np.log(2)), first, first**2)
    original = iterate_pplus(reference, difference, np.ones(2), np.zeros(2),
                             Field(), cycles=3)
    actual = iterate_trusted_pplus(data(), branches, Field(), kl_budget=100,
                                   damping=1)
    np.testing.assert_allclose(actual.probability, original[2], atol=2e-15)
    np.testing.assert_allclose(actual.first, original[1]*np.exp(1j*original[0]), atol=2e-15)


def test_crystal_identity_guard():
    d = data()
    branches = prepare_conditional_branches(d, np.zeros((2, 4)), [.3, .7])
    with pytest.raises(ValueError, match='identity'):
        iterate_trusted_pplus(data([[3, 0, 1], [4, 0, 1]]), branches, Field())


@pytest.mark.parametrize('budget', [0, -.1, np.nan, np.inf, None, True])
def test_iteration_rejects_disabled_or_invalid_budget(budget):
    d = data()
    b = prepare_conditional_branches(d, np.zeros((2, 4)), [.3, .7])
    with pytest.raises(ValueError, match='controls'):
        iterate_trusted_pplus(d, b, Field(), kl_budget=budget)


@pytest.mark.parametrize('mass', [[[np.nan, 0]], [[0, 0]], [[-np.inf, -np.inf]], []])
def test_invalid_mass_rejected(mass):
    with pytest.raises(ValueError, match='log masses'):
        bounded_branch_mass(mass, [1], .01)


def test_invalid_field_and_probability_rejected():
    with pytest.raises(ValueError, match='field'):
        bounded_branch_mass(np.log([[.5, .5]]), [np.nan], .01)
    with pytest.raises(ValueError, match='probabilities'):
        bernoulli_kl([1.1], np.log([[.5, .5]]))
