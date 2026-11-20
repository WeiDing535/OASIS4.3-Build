import gemmi
import numpy as np
import pytest
from phasebuild.geometry import template, prepend, flexible_template


@pytest.mark.parametrize('phi,psi', [(-63., -42.), (-120., 130.), (60., 40.)])
def test_ideal_backbone_geometry(phi, psi):
    xyz = template(phi, psi)
    assert np.linalg.norm(xyz[1, 1]) == pytest.approx(0.)
    for residue in xyz:
        assert np.linalg.norm(residue[0] - residue[1]) == pytest.approx(1.458)
        assert np.linalg.norm(residue[1] - residue[2]) == pytest.approx(1.525)
        assert np.linalg.norm(residue[2] - residue[3]) == pytest.approx(1.231)
    for i in range(2):
        assert np.linalg.norm(xyz[i, 2] - xyz[i + 1, 0]) == pytest.approx(1.329)
        for points, expected in [([xyz[i, 0], xyz[i, 1], xyz[i, 2], xyz[i+1, 0]], psi),
                                 ([xyz[i, 1], xyz[i, 2], xyz[i+1, 0], xyz[i+1, 1]], 180.),
                                 ([xyz[i, 2], xyz[i+1, 0], xyz[i+1, 1], xyz[i+1, 2]], phi)]:
            actual = np.rad2deg(gemmi.calculate_dihedral(*[gemmi.Position(*p) for p in points]))
            delta = (actual - expected + 180.) % 360. - 180.
            assert delta == pytest.approx(0., abs=1e-6)


def test_backward_extension_recovers_ideal_chain():
    xyz = template(-63., -42.)
    np.testing.assert_allclose(prepend(xyz[1], -63., -42.), xyz[0], atol=1e-6)


def test_flexible_template_preserves_ideal_geometry():
    xyz = flexible_template([-42., -63., -42., -63., -42., 0.])
    np.testing.assert_allclose(xyz, template(-63., -42.), atol=1e-6)
    perturbed = flexible_template([-30., -80., 130., -120., 140., 65.])
    for r in perturbed:
        np.testing.assert_allclose(np.linalg.norm(r[0]-r[1]), 1.458, atol=1e-6)
        np.testing.assert_allclose(np.linalg.norm(r[1]-r[2]), 1.525, atol=1e-6)
    for i in range(2):
        np.testing.assert_allclose(np.linalg.norm(perturbed[i, 2]-perturbed[i+1, 0]), 1.329, atol=1e-6)
