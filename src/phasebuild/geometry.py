"""Ideal trans-peptide geometry. Coordinates are always orthogonal Angstroms."""
import numpy as np


def place(a, b, c, length, angle, dihedral):
    """Append an atom using a bond length, bond angle and torsion (degrees)."""
    bc = np.asarray(c) - b
    bc = bc / np.linalg.norm(bc)
    normal = np.cross(np.asarray(b) - a, bc)
    normal /= np.linalg.norm(normal)
    side = np.cross(normal, bc)
    theta, tau = np.deg2rad([angle, dihedral])
    return c + length * (-np.cos(theta) * bc + np.sin(theta) *
                         (np.cos(tau) * side + np.sin(tau) * normal))


def oxygen(n, ca, c, psi):
    return place(n, ca, c, 1.231, 120.8, psi + 180.0)


def extend(last, psi, phi):
    n, ca, c = last[:3]
    nn = place(n, ca, c, 1.329, 116.2, psi)
    aa = place(ca, c, nn, 1.458, 121.7, 180.0)
    cc = place(c, nn, aa, 1.525, 111.2, phi)
    return np.array([nn, aa, cc, oxygen(nn, aa, cc, -40.0)])


def template(phi, psi, length=3):
    n = np.array([-1.458, 0., 0.])
    ca = np.zeros(3)
    c = 1.525 * np.array([np.cos(np.deg2rad(68.8)), np.sin(np.deg2rad(68.8)), 0.])
    residues = [np.array([n, ca, c, oxygen(n, ca, c, psi)])]
    for _ in range(length - 1):
        residues.append(extend(residues[-1], psi, phi))
        residues[-1][3] = oxygen(*residues[-1][:3], psi)
    xyz = np.array(residues)
    return xyz - xyz[length // 2, 1]


# General conformations, not learned from the target/reference structure.
CONFORMATIONS = [(-63., -42.), (-120., 130.), (-80., 150.), (-70., 0.), (60., 40.)]


def background_probes(xyz):
    """Off-atom probes penalize fits into indiscriminately dense regions."""
    atoms = xyz.reshape(-1, 3)
    offsets = 1.4 * np.concatenate([np.eye(3), -np.eye(3)])
    probes = (atoms[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
    distances = np.linalg.norm(probes[:, None, :] - atoms[None, :, :], axis=-1)
    return probes[distances.min(axis=1) > 1.25]


def flexible_template(angles):
    """Tripeptide with independent psi0/phi1/psi1/phi2/psi2 and terminal N."""
    psi0, phi1, psi1, phi2, psi2, n_twist = angles
    first = template(phi1, psi0, length=1)[0]
    second = extend(first, psi0, phi1)
    second[3] = oxygen(*second[:3], psi1)
    third = extend(second, psi1, phi2)
    third[3] = oxygen(*third[:3], psi2)
    xyz = np.array([first, second, third])
    # The missing preceding residue leaves the first N's torsion independent.
    axis = xyz[0, 2] - xyz[0, 1]
    axis /= np.linalg.norm(axis)
    vec = xyz[0, 0] - xyz[0, 1]
    theta = np.deg2rad(n_twist)
    xyz[0, 0] = (xyz[0, 1] + vec * np.cos(theta) + np.cross(axis, vec) * np.sin(theta)
                 + axis * np.dot(axis, vec) * (1. - np.cos(theta)))
    return xyz - xyz[1, 1]


def prepend(first, phi, psi):
    """Build a preceding residue from the current residue's phi and previous psi."""
    n, ca, c = first[:3]
    pc = place(c, ca, n, 1.329, 121.7, phi)
    pca = place(ca, n, pc, 1.525, 116.2, 180.)
    pn = place(n, pc, pca, 1.458, 111.2, psi)
    return np.array([pn, pca, pc, oxygen(pn, pca, pc, psi)])
