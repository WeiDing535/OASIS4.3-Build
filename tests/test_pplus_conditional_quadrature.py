"""Branch-aligned integration checks using independent adaptive quadrature."""
import numpy as np
import pytest
from scipy.integrate import quad

from oasis_py.phasing import sad_hl
from phasebuild.pplus_conditional import prepare_conditional_branches
from test_pplus_conditional import data


def independent_integrals(hl, reference):
    """Adaptive continuous integrals; branch-local rescaling avoids tiny densities.

    A coarse set only picks an arbitrary safe exponent offset. It is not
    used to approximate any integral, which is evaluated independently by
    scipy's adaptive Gauss-Kronrod quadrature.
    """
    hl = np.asarray(hl)
    def log_density(theta):
        return hl @ [np.cos(theta), np.sin(theta),
                     np.cos(2*theta), np.sin(2*theta)]
    def integrate(function, low, high):
        return quad(function, low, high, epsabs=2e-13, epsrel=2e-13,
                    limit=300)[0]
    log_mass, first, second = [], [], []
    for sign in [1, -1]:
        bound = max(log_density(reference+sign*u)
                    for u in np.linspace(0, np.pi, 65))
        def weight(u):
            return np.exp(log_density(reference+sign*u)-bound)
        value = integrate(weight, 0, np.pi)
        log_mass.append(np.log(value)+bound)
        for harmonic, target in [(1, first), (2, second)]:
            re = integrate(lambda u: weight(u)*np.cos(harmonic*(reference+sign*u)), 0, np.pi)
            im = integrate(lambda u: weight(u)*np.sin(harmonic*(reference+sign*u)), 0, np.pi)
            target.append((re+1j*im)/value)
    mass = np.exp(np.array(log_mass)-max(log_mass))
    return mass/mass.sum(), np.array(first), np.array(second)


@pytest.mark.parametrize('hl,reference', [
    ([.8, -.4, .5, -.3], .217),
    ([-.3, .6, -.2, .4], 1.119),
    ([3.2, -2.1, .9, -.7], -.719),
    ([210., -70., -20., 8.], .217),
])
def test_high_precision_asymmetric_branches_match_adaptive_quad(hl, reference):
    b = prepare_conditional_branches(data(1), [hl], [reference], samples=16384)
    mass, first, second = independent_integrals(hl, reference)
    np.testing.assert_allclose(np.exp(b.log_mass[0]), mass, atol=3e-6, rtol=0)
    np.testing.assert_allclose(b.first[0], first, atol=3e-6, rtol=0)
    np.testing.assert_allclose(b.second[0], second, atol=5e-6, rtol=0)


def test_sad_mirror_branch_masses_are_half_even_with_arbitrary_references():
    rng = np.random.default_rng(517209)
    n = 240
    reference = rng.uniform(-4*np.pi, 4*np.pi, n)
    amplitude = np.geomspace(.001, 80., n)
    difference = rng.uniform(-2, 2, n)*amplitude
    hl = sad_hl(reference, amplitude, difference, np.ones(n))
    b = prepare_conditional_branches(data(n), hl, reference, samples=1440)
    np.testing.assert_allclose(np.exp(b.log_mass), .5, atol=5e-11, rtol=0)


def test_arbitrary_rotation_preserves_mass_and_rotates_both_moments():
    d = data(2)
    hl = np.array([[.8, -.3, .6, -.1], [-.2, .7, -.4, .3]])
    reference = np.array([.217, -.631])
    b = prepare_conditional_branches(d, hl, reference, samples=1440)
    angle = .123456789
    rotated = hl.copy()
    for harmonic, columns in [(1, (0, 1)), (2, (2, 3))]:
        value = (hl[:, columns[0]]+1j*hl[:, columns[1]])*np.exp(1j*harmonic*angle)
        rotated[:, columns[0]], rotated[:, columns[1]] = value.real, value.imag
    r = prepare_conditional_branches(d, rotated, reference+angle, samples=1440)
    np.testing.assert_allclose(r.log_mass, b.log_mass, atol=3e-14, rtol=0)
    np.testing.assert_allclose(r.first, b.first*np.exp(1j*angle), atol=3e-14, rtol=0)
    np.testing.assert_allclose(r.second, b.second*np.exp(2j*angle), atol=3e-14, rtol=0)


def test_combined_first_and_second_moments_match_continuous_circle():
    hl = np.array([1.4, .3, .2, -.8])
    reference = .317
    b = prepare_conditional_branches(data(1), [hl], [reference], samples=1440)
    mass, first, second = independent_integrals(hl, reference)
    np.testing.assert_allclose(np.exp(b.log_mass[0])@b.first[0], mass@first, atol=2e-13, rtol=0)
    np.testing.assert_allclose(np.exp(b.log_mass[0])@b.second[0], mass@second, atol=2e-13, rtol=0)


def test_half_circle_integrals_show_second_order_midpoint_convergence():
    hl, reference = np.array([1.4, .3, .2, -.8]), .317
    mass, first, second = independent_integrals(hl, reference)
    errors = []
    for samples in [72, 144, 288]:
        b = prepare_conditional_branches(data(1), [hl], [reference], samples=samples)
        errors.append(max(abs(np.exp(b.log_mass[0])-mass).max(),
                          abs(b.first[0]-first).max(), abs(b.second[0]-second).max()))
    np.testing.assert_array_less(3.5, np.array(errors[:-1])/errors[1:])
    np.testing.assert_array_less(np.array(errors[:-1])/errors[1:], 4.5)




def test_total_sample_semantics_are_preserved():
    b = prepare_conditional_branches(data(1), [[0., 0., 0., 0.]], [.217], samples=1440)
    assert b.metadata['angular_samples'] == 1440
    assert 'samples total nodes' in b.metadata['angular_grid']
    assert 'no boundary nodes' in b.metadata['acentric_boundary_rule']
