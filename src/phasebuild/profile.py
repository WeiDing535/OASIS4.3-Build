"""Independent local density profile trained on an unrelated known structure.

Uses an orthonormal N/CA/C frame and inverse-variance weighted correlation.
This is not Buccaneer's LLK implementation or a calibrated probability model.
"""
import json
import gemmi
import numpy as np
from .data import fingerprint


def frames(xyz):
    ca = xyz[..., 1, :]
    x = xyz[..., 2, :] - ca
    x /= np.maximum(np.linalg.norm(x, axis=-1, keepdims=True), 1e-8)
    y = xyz[..., 0, :] - ca
    y -= (y * x).sum(axis=-1, keepdims=True) * x
    y /= np.maximum(np.linalg.norm(y, axis=-1, keepdims=True), 1e-8)
    return ca, np.stack([x, y, np.cross(x, y)], axis=-1)


def patches(density, xyz, offsets):
    ca, axes = frames(xyz)
    points = np.einsum('...ij,pj->...pi', axes, offsets) + ca[..., None, :]
    return density.sample(points)


class DensityProfile:
    def __init__(self, offsets, mean, variance, metadata):
        self.offsets, self.mean = np.asarray(offsets), np.asarray(mean)
        self.variance = np.asarray(variance)
        if (self.offsets.ndim != 2 or self.offsets.shape[1] != 3 or
                self.mean.shape != (len(self.offsets),) or self.variance.shape != self.mean.shape or
                not all(np.isfinite(a).all() for a in [self.offsets, self.mean, self.variance]) or
                np.any(self.variance < 0)):
            raise ValueError('Invalid density-profile arrays')
        self.weight = 1. / (self.variance + .25)
        self.weight /= self.weight.sum()
        self.centered = self.mean - np.sum(self.mean * self.weight)
        self.norm = np.sqrt(np.sum(self.centered**2 * self.weight))
        self.metadata = metadata
        if self.norm < 1e-6:
            raise ValueError('Density profile has no spatial contrast')

    def correlation(self, values):
        values = values - np.sum(values * self.weight, axis=-1, keepdims=True)
        covariance = np.sum(values * self.centered * self.weight, axis=-1)
        norm = np.sqrt(np.sum(values**2 * self.weight, axis=-1)) * self.norm
        return covariance / np.maximum(norm, 1e-8)

    def score(self, density, xyz):
        atom_density = np.clip(density.sample(xyz), -3., 6.).mean(axis=-1)
        return 4. * self.correlation(patches(density, xyz, self.offsets)) + .15 * atom_density

    def save(self, path):
        np.savez_compressed(path, offsets=self.offsets, mean=self.mean,
                            variance=self.variance,
                            metadata=json.dumps(self.metadata))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            model = cls(data['offsets'], data['mean'], data['variance'], json.loads(str(data['metadata'])))
            return model


def train(density, structure_path, map_path, radius=3., spacing=.8):
    extent = int(np.floor(radius / spacing))
    grid = np.arange(-extent, extent + 1) * spacing
    offsets = np.array(np.meshgrid(grid, grid, grid, indexing='ij')).reshape(3, -1).T
    offsets = offsets[np.linalg.norm(offsets, axis=1) <= radius]
    structure = gemmi.read_structure(str(structure_path))
    xyz = []
    for chain in structure[0]:
        for residue in chain:
            atoms = [residue.find_atom(name, '*') for name in ['N', 'CA', 'C', 'O']]
            if any(atom is None for atom in atoms):
                continue
            coords = np.array([list(atom.pos) for atom in atoms])
            if not (1.2 < np.linalg.norm(coords[0]-coords[1]) < 1.7 and
                    1.2 < np.linalg.norm(coords[1]-coords[2]) < 1.8):
                continue
            xyz.append(coords)
    if len(xyz) < 30:
        raise ValueError('At least 30 valid backbone groups required for training')
    values = patches(density, np.array(xyz), offsets).astype(float)
    values -= values.mean(axis=-1, keepdims=True)
    values /= np.maximum(values.std(axis=-1, keepdims=True), 1e-6)
    variance = values.var(axis=0)
    metadata = {'method': 'inverse-variance weighted local density correlation',
                'training_structure': fingerprint(structure_path), 'training_map': fingerprint(map_path),
                'residues': len(xyz), 'radius': radius, 'spacing': spacing,
                'resolution_high': density.mtz.resolution_high(),
                'warning': 'Not a probability; training reference must be unrelated to the test protein.'}
    model = DensityProfile(offsets, values.mean(axis=0), variance, metadata)
    scores = model.correlation(values)
    metadata['training_positive_correlation_quantiles'] = np.quantile(scores, [.05,.5,.95]).tolist()
    return model
