"""NumPy restrained reciprocal-space refinement; no external program is called.

This is the internal replacement for the REFMAC stage of the OASIS pipeline.
Everything here is NumPy/SciPy/Gemmi: structure factors are summed directly
(no FFT), the bulk solvent is a flat mask model, and the restraints are bond,
angle and non-bonded repulsion terms built from the local monomer library.

Documented limits of this prototype, all of which are ordinary engineering
work rather than new science:

* direct summation, so the cost grows as n_reflections x n_atoms x n_symmetry
  (a production version needs an FFT/NUFFT structure-factor engine);
* no planarity, chirality, torsion or NCS restraints yet;
* maximum-likelihood target not implemented; the target is least squares with
  sigma weights and an analytic overall scale;
* no bulk-solvent gradient with respect to atomic coordinates;
* centric reflections are not special-cased (all ten calibration space groups
  are primitive or body-centred acentric).
"""
import json
import math
from pathlib import Path

import gemmi
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from .chemistry import BACKBONE, monomer
from .data import fingerprint
from .refinement import reflection_partition

# Cromer-Mann four-Gaussian coefficients (International Tables for
# Crystallography Vol. C, table 6.1.1.4) for the elements that occur in the
# calibration set; the fallback is a single Gaussian with the correct electron
# count, which is adequate only at low resolution.
CROMER_MANN = {
    'H': (0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305),
    'C': (2.31000, 20.8439, 1.02000, 10.2075, 1.58860, 0.568700, 0.865000, 51.6512, 0.215600),
    'N': (12.2126, 0.005700, 3.13220, 9.89330, 2.01250, 28.9975, 1.16630, 0.582600, -11.5290),
    'O': (3.04850, 13.2771, 2.28680, 5.70110, 1.54630, 0.323900, 0.867000, 32.9089, 0.250800),
    'F': (3.53920, 10.2825, 2.64120, 4.29440, 1.51700, 0.261500, 1.02430, 26.1476, 0.277600),
    'P': (6.43450, 1.90670, 4.17910, 27.1570, 1.78000, 0.526000, 1.64420, 68.6700, 1.90670),
    'S': (6.90530, 1.46790, 5.20340, 22.2151, 1.43790, 0.253600, 1.58630, 56.1720, 0.866900),
    'CL': (11.4604, 0.010400, 7.19640, 1.16620, 6.25560, 18.5194, 1.64550, 47.7784, -9.5574),
    'SE': (17.0006, 2.40980, 5.81960, 0.272600, 3.97310, 15.2372, 4.35430, 43.8163, 2.84090),
    'BR': (17.4504, 2.66400, 7.15610, 0.082400, 5.44540, 16.4441, 3.78460, 42.9615, 2.83210),
    'ZN': (14.0743, 3.26550, 7.03180, 0.233300, 5.16520, 10.3163, 2.41000, 58.7097, 1.30410),
    'MG': (5.42040, 2.82750, 2.17350, 79.2611, 1.22690, 0.380800, 2.30730, 7.19370, 0.858400),
    'CA': (8.62660, 10.4421, 7.38730, 0.659900, 1.58990, 85.7484, 1.02110, 178.437, 1.37510),
    'MN': (11.2819, 5.34090, 7.35730, 0.343200, 3.01930, 17.8674, 2.24410, 83.7543, 1.08960),
    'FE': (11.7695, 4.76110, 7.35730, 0.307200, 3.52220, 15.3535, 2.30450, 76.8805, 1.03690),
    'CU': (13.3380, 3.58280, 7.16760, 0.247000, 5.61580, 11.3966, 1.67350, 64.8126, 1.19100),
    'CD': (19.4238, 1.90140, 6.74170, 9.83910, 4.27670, 0.310200, 2.67860, 53.4941, 1.83130),
    'TB': (29.9551, 1.89130, 15.9692, 0.344000, 11.1969, 13.4809, 5.94060, 59.2863, 7.79550),
    'HG': (35.5078, 1.57800, 22.5787, 8.15670, 12.8107, 0.295400, 5.70440, 45.8290, 10.5148),
    'AU': (16.8819, 0.461100, 18.5913, 8.62160, 25.5582, 1.48260, 5.86000, 36.3956, 12.0658),
}


def _cromer_mann(element, s):
    """Scattering factor f(s), s = sin(theta)/lambda, in electrons."""
    element = element.upper()
    if element in CROMER_MANN:
        a1, b1, a2, b2, a3, b3, a4, b4, c = CROMER_MANN[element]
        s2 = s * s
        return (a1 * np.exp(-b1 * s2) + a2 * np.exp(-b2 * s2)
                + a3 * np.exp(-b3 * s2) + a4 * np.exp(-b4 * s2) + c)
    electrons = float(gemmi.Element(element).atomic_number)
    return electrons * np.exp(-6.0 * s * s)


def _phase_factors(arguments):
    """exp(2*pi*i*x) evaluated as cos/sin in single precision.

    Single precision is deliberate: the phase argument is a large multiple of
    2*pi, so the absolute error is dominated by the argument itself, and the
    float32 trigonometric path is roughly four times faster than the complex
    exponential for the array sizes used here.
    """
    angle = (2.0 * np.pi) * arguments.astype(np.float32, copy=False)
    phase = np.empty(angle.shape, dtype=np.complex64)
    np.cos(angle, out=phase.real)
    np.sin(angle, out=phase.imag)
    return phase


class Atoms:
    """Flat atom table with per-atom element, occupancy and B factor."""

    def __init__(self, structure, chain_filter=None):
        self.names, self.residues, self.elements = [], [], []
        self.xyz, self.occ, self.b = [], [], []
        self.residue_start, self.residue_names, self.chain_of_residue = [], [], []
        for chain in structure[0]:
            if chain_filter is not None and chain.name not in chain_filter:
                continue
            for residue in chain:
                start = len(self.names)
                for atom in residue:
                    if atom.altloc not in ('', 'A', '\0'):
                        continue
                    if atom.element.name in ('H', 'D'):
                        continue
                    self.names.append(atom.name)
                    self.residues.append(len(self.residue_names))
                    self.elements.append(atom.element.name.upper())
                    self.xyz.append([atom.pos.x, atom.pos.y, atom.pos.z])
                    self.occ.append(float(atom.occ))
                    self.b.append(float(atom.b_iso))
                if len(self.names) > start:
                    self.residue_start.append(start)
                    self.residue_names.append(residue.name.upper())
                    self.chain_of_residue.append(chain.name)
        self.names = np.array(self.names)
        self.residues = np.array(self.residues, dtype=int)
        self.elements = np.array(self.elements)
        self.xyz = np.array(self.xyz, dtype=float)
        self.occ = np.array(self.occ, dtype=float)
        self.b = np.array(self.b, dtype=float)
        self.residue_start = np.array(self.residue_start, dtype=int)
        self.unique_elements = sorted(set(self.elements.tolist()))

    def __len__(self):
        return len(self.names)

    def residue_span(self, index):
        start = self.residue_start[index]
        stop = self.residue_start[index + 1] if index + 1 < len(self.residue_start) else len(self.names)
        return start, stop


def _residue_atom_index(atoms, residue_index):
    start, stop = atoms.residue_span(residue_index)
    return {atoms.names[i]: i for i in range(start, stop)}


def build_restraints(atoms, monomer_root, bond_sigma=0.02, angle_sigma=2.0,
                     peptide_sigma=0.02, repulsion_radius=3.0, repulsion_weight=1.0,
                     skip_unknown=True):
    """Bond, angle, peptide-bond and non-bonded repulsion restraints."""
    bonds, angles, unknown = [], [], []
    cache = {}
    for residue_index, name in enumerate(atoms.residue_names):
        index = _residue_atom_index(atoms, residue_index)
        template = None
        if name in cache:
            template = cache[name]
        elif name not in ('HOH', 'WAT', 'DOD'):
            try:
                template = monomer(name, monomer_root)
            except Exception:
                template = None
            cache[name] = template
        if template is None:
            if name not in ('HOH', 'WAT', 'DOD'):
                unknown.append(name)
            continue
        for first, second, _order, distance in template['bonds']:
            if first in index and second in index:
                bonds.append((index[first], index[second], float(distance), bond_sigma))
        for first, second, third, angle in template['angles']:
            if first in index and second in index and third in index:
                angles.append((index[first], index[second], index[third], float(angle), angle_sigma))
    # Peptide bonds and the C-N-CA angle are inter-residue and absent from the
    # monomer dictionary.
    chain_residues = {}
    for residue_index, chain in enumerate(atoms.chain_of_residue):
        chain_residues.setdefault(chain, []).append(residue_index)
    for chain, order in chain_residues.items():
        for left, right in zip(order, order[1:]):
            first = _residue_atom_index(atoms, left)
            second = _residue_atom_index(atoms, right)
            if 'C' in first and 'N' in second:
                distance = np.linalg.norm(atoms.xyz[first['C']] - atoms.xyz[second['N']])
                if distance < 2.0:
                    bonds.append((first['C'], second['N'], 1.329, peptide_sigma))
                    if 'CA' in first:
                        angles.append((first['CA'], first['C'], second['N'], 116.2, angle_sigma))
                    if 'CA' in second:
                        angles.append((first['C'], second['N'], second['CA'], 121.7, angle_sigma))
    bonds = np.array(bonds, dtype=float).reshape(-1, 4)
    angles = np.array(angles, dtype=float).reshape(-1, 5)
    # Non-bonded pairs: close contacts that are neither bonded nor 1-3.
    tree = cKDTree(atoms.xyz)
    pairs = np.array(sorted(tree.query_pairs(repulsion_radius)), dtype=int).reshape(-1, 2)
    excluded = set()
    for i, j, _d, _s in bonds:
        excluded.add((min(int(i), int(j)), max(int(i), int(j))))
    neighbours = {}
    for i, j, k, _t, _s in angles:
        for a, b in ((i, j), (j, i), (k, j), (j, k), (i, k), (k, i)):
            neighbours.setdefault((int(a), int(b)), set()).add(int(j))
    keep = []
    for i, j in pairs:
        key = (min(int(i), int(j)), max(int(i), int(j)))
        if key in excluded:
            continue
        if int(j) in neighbours.get((int(i), int(j)), ()):
            continue
        keep.append(key)
    repulsion = np.array(sorted(set(keep)), dtype=int).reshape(-1, 2)
    return {'bonds': bonds, 'angles': angles, 'repulsion': repulsion,
            'repulsion_radius': float(repulsion_radius),
            'repulsion_weight': float(repulsion_weight),
            'unknown_residues': sorted(set(unknown))}


class Geometry:
    """Vectorised restraint energy and gradient over Cartesian coordinates."""

    def __init__(self, restraints):
        self.bonds = restraints['bonds']
        self.angles = restraints['angles']
        self.repulsion = restraints['repulsion']
        self.radius = restraints['repulsion_radius']
        self.repulsion_weight = restraints['repulsion_weight']
        self.bond_weight = 1.0 / np.maximum(self.bonds[:, 3], 1e-6) ** 2 if len(self.bonds) else np.zeros(0)
        self.angle_weight = (1.0 / np.deg2rad(np.maximum(self.angles[:, 4], 1e-6)) ** 2
                             if len(self.angles) else np.zeros(0))

    def energy_gradient(self, xyz):
        energy = 0.0
        gradient = np.zeros_like(xyz)
        if len(self.bonds):
            i = self.bonds[:, 0].astype(int)
            j = self.bonds[:, 1].astype(int)
            delta = xyz[i] - xyz[j]
            distance = np.linalg.norm(delta, axis=1)
            deviation = distance - self.bonds[:, 2]
            weight = self.bond_weight
            energy += float(np.sum(weight * deviation ** 2))
            factor = (2.0 * weight * deviation / np.maximum(distance, 1e-9))[:, None] * delta
            np.add.at(gradient, i, factor)
            np.add.at(gradient, j, -factor)
        if len(self.angles):
            i = self.angles[:, 0].astype(int)
            j = self.angles[:, 1].astype(int)
            k = self.angles[:, 2].astype(int)
            first = xyz[i] - xyz[j]
            second = xyz[k] - xyz[j]
            first_length = np.linalg.norm(first, axis=1)
            second_length = np.linalg.norm(second, axis=1)
            cosine = np.clip(np.sum(first * second, axis=1) / (first_length * second_length), -1.0, 1.0)
            angle = np.rad2deg(np.arccos(cosine))
            deviation = angle - self.angles[:, 3]
            weight = self.angle_weight
            energy += float(np.sum(weight * deviation ** 2))
            # The energy uses degrees, so the radian derivative needs the
            # 180/pi factor: dE/dx = 2 w dtheta_deg (180/pi) dtheta_rad/dx.
            sine = np.maximum(np.sin(np.deg2rad(angle)), 1e-6)
            scale = (-2.0 * weight * deviation * (180.0 / np.pi) / sine)[:, None]
            first_unit = first / np.maximum(first_length, 1e-9)[:, None]
            second_unit = second / np.maximum(second_length, 1e-9)[:, None]
            d_first = (second_unit - cosine[:, None] * first_unit) / np.maximum(first_length, 1e-9)[:, None]
            d_second = (first_unit - cosine[:, None] * second_unit) / np.maximum(second_length, 1e-9)[:, None]
            np.add.at(gradient, i, scale * d_first)
            np.add.at(gradient, k, scale * d_second)
            np.add.at(gradient, j, -scale * (d_first + d_second))
        if len(self.repulsion):
            i = self.repulsion[:, 0].astype(int)
            j = self.repulsion[:, 1].astype(int)
            delta = xyz[i] - xyz[j]
            distance = np.linalg.norm(delta, axis=1)
            violation = np.maximum(self.radius - distance, 0.0)
            active = violation > 0
            if np.any(active):
                weight = self.repulsion_weight
                energy += float(np.sum(weight * violation[active] ** 2))
                factor = (-2.0 * weight * violation[active]
                          / np.maximum(distance[active], 1e-9))[:, None] * delta[active]
                np.add.at(gradient, i[active], factor)
                np.add.at(gradient, j[active], -factor)
        return energy, gradient


class Scattering:
    """Per-reflection, per-element scattering factor table."""

    def __init__(self, atoms, reflections):
        self.elements = atoms.unique_elements
        index = {element: position for position, element in enumerate(self.elements)}
        self.atom_element = np.array([index[element] for element in atoms.elements], dtype=int)
        s = reflections['s']
        self.table = np.stack([_cromer_mann(element, s) for element in self.elements], axis=1)
        self.s2 = s * s


class Reflections:
    """Observed data with a pre-existing free set."""

    def __init__(self, path, free_value=0, resolution=None, fp='F', sigfp='SIGF', free='FreeR_flag'):
        mtz, partition = reflection_partition(path, fp=fp, sigfp=sigfp, free=free, free_value=free_value)
        values = np.array(mtz, copy=True)
        labels = mtz.column_labels()
        f = values[:, labels.index(fp)]
        sigma = values[:, labels.index(sigfp)]
        flag = values[:, labels.index(free)]
        usable = np.isfinite(f) & np.isfinite(sigma) & (f >= 0) & (sigma > 0)
        spacing = mtz.make_d_array()
        if resolution is not None:
            usable &= spacing >= resolution
        hkl = np.rint(values[:, :3]).astype(int)
        self.hkl = hkl[usable]
        self.f = f[usable]
        self.sigma = sigma[usable]
        self.free = flag[usable] == free_value
        self.spacing = spacing[usable]
        self.partition = partition
        self.s = 0.5 / self.spacing
        self.mtz = mtz

    def __len__(self):
        return len(self.f)


class BulkSolvent:
    """Flat bulk-solvent model on a periodic FFT grid.

    The grid holds the solvent region (1 outside the molecular envelope, 0
    inside), so its Fourier transform is the complement of the envelope
    transform and interferes destructively with the atoms at low resolution,
    which is what a bulk-solvent correction has to do.  The envelope is built
    from atoms wrapped into the cell *and* their 26 periodic images: without
    the images a molecule that crosses a cell face loses part of its envelope
    and the transform picks up a spurious boundary term.

    ``probe`` and ``shrinkage`` follow REFMAC: a grid point is inside the
    envelope when it is within ``probe * shrinkage`` of an atom.
    """

    def __init__(self, cell, spacegroup, atoms, spacing=1.0, probe=1.2,
                 shrinkage=0.8, smoothing=0.8):
        self.cell = cell
        self.spacing = spacing
        lengths = [float(length) for length in cell.parameters[:3]]
        self.grid_shape = tuple(max(8, int(math.ceil(length / spacing)))
                                for length in lengths)
        # The mesh has to span exactly one cell: a grid that overshoots the
        # cell edge shears the phases and decorrelates the solvent transform.
        self.grid_step = tuple(length / count
                               for length, count in zip(lengths, self.grid_shape))
        n1, n2, n3 = self.grid_shape
        orthogonal = np.array(cell.orth.mat.tolist(), dtype=float)
        fractional = atoms.xyz @ np.linalg.inv(orthogonal).T
        fractional -= np.floor(fractional)
        shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1)
                           for k in (-1, 0, 1)], dtype=float)
        images = ((fractional[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
                  @ orthogonal.T)
        radii = np.tile(self._radii(atoms, probe, shrinkage), len(shifts))
        axis = [np.arange(n) * step for n, step in zip(self.grid_shape, self.grid_step)]
        mesh = np.stack(np.meshgrid(*axis, indexing='ij'), axis=-1).reshape(-1, 3)
        tree = cKDTree(images)
        # A single radius per atom keeps the envelope from developing the holes
        # that a flat skin leaves between bonded atoms; a holey envelope has a
        # rapidly decaying transform and cannot correct the low-resolution
        # amplitudes the way a real bulk-solvent term must.
        candidates = tree.query_ball_point(mesh, float(np.max(radii)))
        counts = np.fromiter((len(entry) for entry in candidates), dtype=int,
                             count=len(candidates))
        point_index = np.repeat(np.arange(len(mesh)), counts)
        atom_index = (np.concatenate([np.asarray(entry, dtype=int) for entry in candidates])
                      if counts.sum() else np.zeros(0, dtype=int))
        distance = np.sum((mesh[point_index] - images[atom_index]) ** 2, axis=1)
        inside = distance < radii[atom_index] ** 2
        mask = np.ones(len(mesh), dtype=float)
        mask[point_index[inside]] = 0.0
        mask = mask.reshape(self.grid_shape)
        if smoothing:
            mask = gaussian_filter(mask, float(smoothing), mode='wrap')
        self.grid = mask
        self.solvent_fraction = float(mask.mean())
        self.volume_element = float(np.prod(cell.parameters[:3])) / float(n1 * n2 * n3)
        self.spacegroup = spacegroup

    @staticmethod
    def _radii(atoms, probe, shrinkage):
        """Per-atom envelope radius: shrunk van der Waals radius plus probe."""
        table = {}
        radii = np.empty(len(atoms), dtype=float)
        for index, element in enumerate(atoms.elements):
            key = element.upper()
            if key not in table:
                try:
                    table[key] = float(gemmi.Element(key).vdw_r)
                except Exception:
                    table[key] = 1.7
            radii[index] = shrinkage * table[key] + probe
        return radii

    def structure_factors(self, hkl):
        transform = np.fft.ifftn(self.grid)
        n1, n2, n3 = self.grid_shape
        index = (hkl[:, 0] % n1, hkl[:, 1] % n2, hkl[:, 2] % n3)
        return self.volume_element * n1 * n2 * n3 * transform[index]


class Refiner:
    """Least-squares restrained refinement of coordinates and isotropic B.

    Structure factors are summed directly in chunks of reflections, so the
    peak memory stays proportional to chunk_size x n_atoms.
    """

    def __init__(self, atoms, reflections, spacegroup, cell, restraints,
                 refine_b=True, refine_solvent=True, bulk_solvent=True,
                 geometry_weight=1.0, b_restraint_sigma=5.0, b_prior=None,
                 chunk_size=1500):
        self.atoms = atoms
        self.data = reflections
        self.spacegroup = spacegroup
        self.cell = cell
        self.refine_b = refine_b
        self.refine_solvent = refine_solvent
        self.geometry_weight = geometry_weight
        self.b_restraint_sigma = b_restraint_sigma
        self.b_prior = atoms.b.copy() if b_prior is None else np.array(b_prior, dtype=float)
        self.chunk_size = int(chunk_size)
        self.geometry = Geometry(restraints)
        self.restraint_summary = {key: int(len(restraints[key])) for key in ('bonds', 'angles', 'repulsion')}
        self.scattering = Scattering(atoms, {'s': reflections.s})
        self.orthogonal = np.array(cell.orth.mat.tolist(), dtype=float)
        self.ops = self._symmetry_operators(spacegroup)
        # REFMAC-style overall B applied to the atomic part only; it is a
        # scaling parameter, never folded into the per-atom B values.
        self.b_add = 0.0
        self.bulk = BulkSolvent(cell, spacegroup, atoms) if bulk_solvent else None
        self.f_bulk = (self.bulk.structure_factors(reflections.hkl) if self.bulk is not None
                       else np.zeros(len(reflections), dtype=complex))
        self.hkl_rotated = [reflections.hkl @ rotation for rotation, _t in self.ops]
        self.phase_offsets = [np.exp(2j * np.pi * (reflections.hkl @ translation)) for _r, translation in self.ops]

    def _to_fractional(self, xyz):
        return xyz @ np.linalg.inv(self.orthogonal).T

    def _to_cartesian(self, fractional):
        return fractional @ self.orthogonal.T

    def _symmetry_operators(self, spacegroup):
        operators = []
        for op in spacegroup.operations():
            rotation = np.array([[op.rot[i][j] for j in range(3)] for i in range(3)], dtype=float) / op.DEN
            translation = np.array([op.tran[i] / op.DEN for i in range(3)], dtype=float)
            operators.append((rotation, translation))
        return operators

    def _pack(self, xyz, b, scale, k_sol, b_sol):
        parts = [self._to_fractional(xyz).reshape(-1)]
        if self.refine_b:
            parts.append(b)
        parts.append(np.array([math.log(max(scale, 1e-6))]))
        if self.refine_solvent:
            parts.append(np.array([k_sol, b_sol]))
        return np.concatenate(parts)

    def _unpack(self, parameters):
        count = len(self.atoms)
        xyz_fractional = parameters[:3 * count].reshape(count, 3)
        offset = 3 * count
        if self.refine_b:
            b = parameters[offset:offset + count]
            offset += count
        else:
            b = self.atoms.b
        scale = math.exp(float(parameters[offset]))
        offset += 1
        if self.refine_solvent:
            k_sol = float(parameters[offset])
            b_sol = float(parameters[offset + 1])
        else:
            k_sol, b_sol = 0.35, 50.0
        return xyz_fractional, b, scale, k_sol, b_sol

    def _amplitude(self, b):
        attenuation = np.exp(-0.25 * self.scattering.s2[:, None] * b[None, :])
        return self.scattering.table[:, self.scattering.atom_element] * attenuation * self.atoms.occ[None, :]

    def _structure_factors(self, xyz_fractional, b, scale, k_sol, b_sol):
        amplitude = self._amplitude(b)
        bulk = k_sol * np.exp(-0.25 * self.scattering.s2 * b_sol) * self.f_bulk
        total = np.zeros(len(self.data), dtype=complex)
        for rotated, offset in zip(self.hkl_rotated, self.phase_offsets):
            phase = _phase_factors(rotated @ xyz_fractional.T)
            total += offset * np.sum(amplitude * phase, axis=1)
        total = self.atomic_attenuation * total
        return scale * (total + bulk), amplitude, bulk, total

    @property
    def atomic_attenuation(self):
        """Overall isotropic B applied to the atomic contribution.

        Always an array: ``_target`` indexes it per reflection chunk, so the
        zero-B case has to return ones rather than the scalar 1.0.
        """
        return np.exp(-0.25 * self.scattering.s2 * self.b_add)

    def _target(self, parameters, gradient=True):
        xyz_fractional, b, scale, k_sol, b_sol = self._unpack(parameters)
        xyz = self._to_cartesian(xyz_fractional)
        count = len(self.atoms)
        amplitude = self._amplitude(b)
        solvent_shape = np.exp(-0.25 * self.scattering.s2 * b_sol)
        bulk = k_sol * solvent_shape * self.f_bulk
        atomic_attenuation = self.atomic_attenuation
        weight = 1.0 / np.maximum(self.data.sigma, 1e-6) ** 2
        observed = self.data.f
        n = len(observed)
        fcalc = np.zeros(n, dtype=complex)
        energy = 0.0
        gradient_xyz = np.zeros((count, 3))
        gradient_b = np.zeros(count) if self.refine_b else None
        gradient_scale = 0.0
        gradient_ksol = 0.0
        gradient_bsol = 0.0
        for start in range(0, n, self.chunk_size):
            stop = min(start + self.chunk_size, n)
            sl = slice(start, stop)
            chunk_amplitude = amplitude[sl]
            chunk_phases = []
            total = np.zeros(stop - start, dtype=complex)
            for rotated, offset in zip(self.hkl_rotated, self.phase_offsets):
                phase = _phase_factors(rotated[sl] @ xyz_fractional.T) * offset[sl][:, None].astype(np.complex64)
                chunk_phases.append(phase)
                total += np.sum(chunk_amplitude * phase, axis=1)
            chunk_fcalc = scale * (atomic_attenuation[sl] * total + bulk[sl])
            fcalc[sl] = chunk_fcalc
            amplitude_fcalc = np.maximum(np.abs(chunk_fcalc), 1e-9)
            residual = observed[sl] - amplitude_fcalc
            energy += float(np.sum(weight[sl] * residual ** 2))
            if not gradient:
                continue
            d_fcalc = -2.0 * weight[sl] * residual * chunk_fcalc / amplitude_fcalc
            d_total = scale * atomic_attenuation[sl] * d_fcalc
            conjugate = np.conj(d_total)[:, None]
            for (rotated, offset), phase in zip(zip(self.hkl_rotated, self.phase_offsets), chunk_phases):
                product = conjugate * chunk_amplitude * phase
                gradient_xyz += -2.0 * np.pi * (np.imag(product).T @ rotated[sl])
                if self.refine_b:
                    gradient_b += -0.25 * np.sum(
                        self.scattering.s2[sl][:, None] * np.real(product), axis=0)
            gradient_scale += float(np.real(np.sum(np.conj(d_fcalc) * chunk_fcalc)))
            if self.refine_solvent:
                solvent = solvent_shape[sl] * self.f_bulk[sl]
                gradient_ksol += float(np.real(np.sum(np.conj(d_total) * solvent)))
                gradient_bsol += float(np.real(np.sum(
                    np.conj(d_total) * k_sol * (-0.25 * self.scattering.s2[sl]) * solvent)))
        energy_geometry, gradient_geometry = self.geometry.energy_gradient(xyz)
        energy += self.geometry_weight * energy_geometry
        # The X-ray term is already differentiated with respect to fractional
        # coordinates; the restraint term is Cartesian and must be converted.
        gradient_fractional = gradient_xyz + self.geometry_weight * (gradient_geometry @ self.orthogonal.T)
        parameters_gradient = np.zeros_like(parameters)
        parameters_gradient[:3 * count] = gradient_fractional.reshape(-1)
        offset = 3 * count
        if self.refine_b:
            parameters_gradient[offset:offset + count] = (
                gradient_b + 2.0 * (b - self.b_prior) / self.b_restraint_sigma ** 2)
            offset += count
        parameters_gradient[offset] = gradient_scale
        offset += 1
        if self.refine_solvent:
            parameters_gradient[offset] = gradient_ksol
            parameters_gradient[offset + 1] = gradient_bsol
        return energy, parameters_gradient

    def metrics(self, fcalc):
        amplitude = np.abs(fcalc)
        observed = self.data.f
        work = ~self.data.free
        free = self.data.free

        def r_factor(selection):
            denominator = float(np.sum(np.abs(observed[selection])))
            if denominator <= 0:
                return None
            return float(np.sum(np.abs(observed[selection] - amplitude[selection])) / denominator)

        bond_rmsd = angle_rmsd = None
        xyz = self._to_cartesian(self._unpack(self.parameters)[0])
        if len(self.geometry.bonds):
            i = self.geometry.bonds[:, 0].astype(int)
            j = self.geometry.bonds[:, 1].astype(int)
            distance = np.linalg.norm(xyz[i] - xyz[j], axis=1)
            bond_rmsd = float(np.sqrt(np.mean((distance - self.geometry.bonds[:, 2]) ** 2)))
        if len(self.geometry.angles):
            i = self.geometry.angles[:, 0].astype(int)
            j = self.geometry.angles[:, 1].astype(int)
            k = self.geometry.angles[:, 2].astype(int)
            first = xyz[i] - xyz[j]
            second = xyz[k] - xyz[j]
            cosine = np.clip(np.sum(first * second, axis=1)
                             / (np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)), -1.0, 1.0)
            angle = np.rad2deg(np.arccos(cosine))
            angle_rmsd = float(np.sqrt(np.mean((angle - self.geometry.angles[:, 3]) ** 2)))
        return {'r_work': r_factor(work), 'r_free': r_factor(free),
                'bond_rmsd': bond_rmsd, 'angle_rmsd': angle_rmsd,
                'reflections_work': int(work.sum()), 'reflections_free': int(free.sum())}

    def _atom_structure_factors(self, xyz_fractional, b):
        """Unscaled contribution of the atoms alone (bulk solvent excluded)."""
        amplitude = self._amplitude(b)
        total = np.zeros(len(self.data), dtype=complex)
        for rotated, offset in zip(self.hkl_rotated, self.phase_offsets):
            phase = _phase_factors(rotated @ xyz_fractional.T)
            total += offset * np.sum(amplitude * phase, axis=1)
        return total

    def fit_scaling(self, iterations=300):
        """Fit overall scale, overall B and the flat bulk-solvent parameters.

        REFMAC treats the overall scale and an overall B as least-squares
        scaling parameters that are *not* written into the per-atom B values;
        without that second number the atomic contribution cannot be damped at
        high resolution and the residual stalls roughly 6 points of R above
        what the same model can reach.  These four numbers are badly
        conditioned next to thousands of coordinates, so they are solved in
        their own small problem.  The atom part is computed once, which makes
        every trial evaluation cost one pass over the reflections.
        """
        xyz_fractional, b, scale, k_sol, b_sol = self._unpack(self.parameters)
        total = self._atom_structure_factors(xyz_fractional, b)
        observed = self.data.f
        weight = 1.0 / np.maximum(self.data.sigma, 1e-6) ** 2
        s2 = self.scattering.s2
        attenuation = np.exp(-0.25 * s2)
        bulk_shape = self.f_bulk

        def energy(q):
            factor, b_add, solvent_scale, solvent_b = math.exp(q[0]), q[1], math.exp(q[2]), math.exp(q[3])
            atomic = np.exp(-0.25 * s2 * b_add)
            solvent = solvent_scale * attenuation ** solvent_b * bulk_shape
            z = atomic * total + solvent
            fcalc = factor * z
            magnitude = np.abs(fcalc)
            residual = observed - magnitude
            value = float(np.sum(weight * residual ** 2))
            unit = np.conj(z) / np.maximum(magnitude / factor, 1e-12)
            slope = np.empty((4, len(observed)))
            slope[0] = magnitude
            slope[1] = factor * np.real(unit * atomic * total * (-0.25 * s2))
            slope[2] = factor * np.real(unit * solvent)
            slope[3] = factor * np.real(unit * solvent * (-0.25 * s2))
            gradient = -2.0 * (weight * residual) @ slope.T
            return value, gradient

        start = np.array([math.log(max(scale, 1e-6)), self.b_add,
                          math.log(max(k_sol, 1e-6)), math.log(max(b_sol, 1e-3))])
        bounds = [(math.log(0.2), math.log(5.0)),
                  (-30.0, 200.0),
                  (math.log(0.05), math.log(1.0)),
                  (math.log(5.0), math.log(400.0))]
        result = minimize(energy, start, jac=True, method='L-BFGS-B', bounds=bounds,
                          options={'maxiter': iterations, 'ftol': 1e-14})
        factor, b_add = math.exp(result.x[0]), float(result.x[1])
        solvent_scale, solvent_b = math.exp(result.x[2]), math.exp(result.x[3])
        self.b_add = b_add
        count = len(self.atoms)
        offset = 3 * count + (count if self.refine_b else 0)
        self.parameters[offset] = math.log(max(factor, 1e-6))
        self.parameters[offset + 1] = solvent_scale
        self.parameters[offset + 2] = solvent_b
        return factor, b_add, solvent_scale, solvent_b

    def run(self, cycles=10, max_iterations=40, tolerance=1e-4, verbose=True,
            target_bond_rmsd=0.02, scale_iterations=40):
        self.parameters = self._pack(self.atoms.xyz, self.atoms.b, 1.0, 0.35, 50.0)
        count = len(self.atoms)
        free = 3 * count + (count if self.refine_b else 0)
        history = []

        def objective(subset):
            parameters = self.parameters.copy()
            parameters[:free] = subset
            energy, gradient = self._target(parameters)
            return energy, gradient[:free]

        if target_bond_rmsd and cycles > 1:
            # One short run shows how far the diffraction term alone pulls the
            # geometry; the harmonic balance then fixes the relative weight in
            # a single step instead of creeping towards it cycle by cycle.
            probe = minimize(objective, self.parameters[:free].copy(), jac=True,
                             method='L-BFGS-B',
                             options={'maxiter': 5, 'maxfun': 10, 'maxls': 20})
            self.parameters[:free] = probe.x
            self.fit_scaling(iterations=40)
            xyz_fractional, b, scale, k_sol, b_sol = self._unpack(self.parameters)
            fcalc, _a, _b, _t = self._structure_factors(xyz_fractional, b, scale, k_sol, b_sol)
            probe_metrics = self.metrics(fcalc)
            if probe_metrics['bond_rmsd']:
                ratio = probe_metrics['bond_rmsd'] / target_bond_rmsd
                self.geometry_weight *= float(np.clip(ratio ** 2, 1.0, 400.0))
            if verbose:
                print(json.dumps({'calibration': True, 'bond_rmsd': probe_metrics['bond_rmsd'],
                                  'r_work': probe_metrics['r_work'],
                                  'geometry_weight': self.geometry_weight}), flush=True)

        for cycle in range(1, cycles + 1):
            def objective(subset):
                parameters = self.parameters.copy()
                parameters[:free] = subset
                energy, gradient = self._target(parameters)
                return energy, gradient[:free]
            _energy, gradient = objective(self.parameters[:free])
            gradient_norm = max(float(np.linalg.norm(gradient)), 1.0)
            result = minimize(objective, self.parameters[:free].copy(), jac=True,
                              method='L-BFGS-B',
                              options={'maxiter': max_iterations, 'maxfun': 2 * max_iterations,
                                       'ftol': 1e-10, 'gtol': tolerance * gradient_norm,
                                       'maxls': 20})
            self.parameters[:free] = result.x
            self.fit_scaling(iterations=scale_iterations)
            xyz_fractional, b, scale, k_sol, b_sol = self._unpack(self.parameters)
            fcalc, _a, _b, _t = self._structure_factors(xyz_fractional, b, scale, k_sol, b_sol)
            metrics = self.metrics(fcalc)
            metrics.update({'cycle': cycle, 'energy': float(result.fun),
                            'iterations': int(result.nit), 'function_evaluations': int(result.nfev),
                            'converged': bool(result.success), 'message': str(result.message),
                            'scale': scale, 'b_add': self.b_add,
                            'k_sol': k_sol, 'b_sol': b_sol,
                            'mean_b': float(np.mean(b)),
                            'geometry_weight': self.geometry_weight})
            history.append(metrics)
            if verbose:
                print(json.dumps({key: metrics[key] for key in
                                  ('cycle', 'r_work', 'r_free', 'bond_rmsd', 'angle_rmsd',
                                   'iterations', 'scale', 'b_add', 'k_sol', 'b_sol',
                                   'geometry_weight')}), flush=True)
            # The relative weight of the diffraction term is not known a
            # priori; adjust it so the model keeps the geometry it started
            # with, which is what a restrained refinement is asked to do.
            if target_bond_rmsd and metrics['bond_rmsd']:
                ratio = metrics['bond_rmsd'] / target_bond_rmsd
                self.geometry_weight *= min(max(ratio ** 2, 0.2), 5.0)
        xyz_fractional, b, scale, k_sol, b_sol = self._unpack(self.parameters)
        self.atoms.xyz = self._to_cartesian(xyz_fractional)
        self.atoms.b = np.array(b)
        return {'history': history, 'scale': scale, 'b_add': self.b_add,
                'k_sol': k_sol, 'b_sol': b_sol}


def refine_structure(model_path, native_path, output_dir, monomer_root,
                     cycles=5, max_iterations=40, free_value=0, resolution=None,
                     refine_b=True, refine_solvent=True, geometry_weight=1.0,
                     b_restraint_sigma=5.0, chunk_size=1500, verbose=True,
                     target_bond_rmsd=0.02):
    """Refine ``model_path`` against ``native_path`` and write the result."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    input_fingerprint = fingerprint(model_path)
    structure = gemmi.read_structure(str(model_path))
    structure.remove_alternative_conformations()
    structure.remove_hydrogens()
    structure.remove_empty_chains()
    spacegroup = structure.find_spacegroup()
    if spacegroup is None or structure.spacegroup_hm in ('', 'P 1'):
        raise ValueError('Model has no usable space group; refusing to guess one')
    atoms = Atoms(structure)
    if not len(atoms):
        raise ValueError('Model contains no heavy atoms')
    reflections = Reflections(native_path, free_value=free_value, resolution=resolution)
    restraints = build_restraints(atoms, monomer_root)
    refiner = Refiner(atoms, reflections, spacegroup, structure.cell, restraints,
                      refine_b=refine_b, refine_solvent=refine_solvent,
                      geometry_weight=geometry_weight, b_restraint_sigma=b_restraint_sigma,
                      chunk_size=chunk_size)
    report = refiner.run(cycles=cycles, max_iterations=max_iterations, verbose=verbose,
                         target_bond_rmsd=target_bond_rmsd)
    # Write the refined coordinates back onto the original hierarchy.
    positions = iter(atoms.xyz)
    b_values = iter(atoms.b)
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.altloc not in ('', 'A', '\0') or atom.element.name in ('H', 'D'):
                    continue
                xyz = next(positions)
                atom.pos = gemmi.Position(float(xyz[0]), float(xyz[1]), float(xyz[2]))
                atom.b_iso = float(next(b_values))
    model_path = output / 'refined.cif'
    structure.write_pdb(str(output / 'refined.pdb'))
    structure.make_mmcif_document().write_file(str(model_path))
    report.update({'input_model': input_fingerprint,
                   'native': reflections.partition['input'],
                   'atoms': len(atoms), 'restraints': refiner.restraint_summary,
                   'unknown_residues': restraints['unknown_residues'],
                   'free_value': int(free_value), 'resolution_cutoff': resolution,
                   'cycles': cycles, 'refine_b': bool(refine_b),
                   'geometry_weight': geometry_weight,
                   'initial': {k: v for k, v in report['history'][0].items() if k != 'message'},
                   'final': {k: v for k, v in report['history'][-1].items() if k != 'message'}})
    (output / 'refinement_report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
