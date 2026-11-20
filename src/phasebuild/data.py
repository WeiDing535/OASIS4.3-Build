import hashlib
from pathlib import Path
import gemmi
import numpy as np


def fingerprint(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path.resolve()), 'sha256': digest.hexdigest()}


def audit(path):
    mtz = gemmi.read_mtz_file(str(path))
    cols = []
    for col in mtz.columns:
        values = np.asarray(col)
        finite = values[np.isfinite(values)]
        dataset = mtz.dataset(col.dataset_id)
        item = {'label': col.label, 'type': col.type, 'dataset_id':col.dataset_id,
                'dataset':dataset.dataset_name, 'valid': int(len(finite)),
                'missing': int(len(values) - len(finite))}
        if len(finite):
            item.update(min=float(finite.min()), max=float(finite.max()))
        if col.label == 'FreeR_flag':
            keys, counts = np.unique(finite, return_counts=True)
            item['counts'] = {str(int(k)): int(v) for k, v in zip(keys, counts)}
        cols.append(item)
    return dict(input=fingerprint(path), title=mtz.title,
                spacegroup=mtz.spacegroup.hm, cell=list(mtz.cell.parameters),
                resolution_high=mtz.resolution_high(), reflections=mtz.nreflections,
                columns=cols)


class Density:
    def __init__(self, path, amplitude='FWT', phase='PHWT', sample_rate=3., resolution=None):
        self.source = fingerprint(path)
        self.mtz = gemmi.read_mtz_file(str(path))
        self.profile = None
        if resolution is not None:
            data = np.array(self.mtz, copy=True)
            self.mtz.set_data(data[self.mtz.make_d_array() >= resolution])
            self.mtz.update_reso()
        labels = self.mtz.column_labels()
        for label, kind in [(amplitude, 'F'), (phase, 'P')]:
            if labels.count(label)>1:
                raise ValueError('Ambiguous duplicate MTZ label: '+label+'; use unique column labels')
            if label not in labels:
                raise ValueError('Missing MTZ column: ' + label)
            if self.mtz.column_with_label(label).type != kind:
                raise ValueError('Unexpected MTZ column type for ' + label)
        arr = np.array(self.mtz, copy=False)
        valid = np.isfinite(arr[:, labels.index(amplitude)]) & np.isfinite(arr[:, labels.index(phase)])
        if not np.any(valid):
            raise ValueError('No finite amplitude/phase pairs')
        # FWT is already a map amplitude: do not multiply by FOM a second time.
        self.grid = self.mtz.transform_f_phi_to_map(amplitude, phase, sample_rate=sample_rate)
        self.cell = self.grid.unit_cell
        self.sg = self.grid.spacegroup
        self.shape = tuple(self.grid.shape)
        values = np.array(self.grid, copy=False)
        self.mean, self.sigma = float(values.mean()), float(values.std())
        if self.sigma <= 0 or not np.isfinite(self.sigma):
            raise ValueError('Map has zero or nonfinite variance')
        values[:] = (values - self.mean) / self.sigma

    def sample(self, xyz):
        shape = xyz.shape[:-1]
        points = np.ascontiguousarray(xyz.reshape(-1, 3), dtype=np.float64)
        return self.grid.interpolate_position_array(points, order=1).reshape(shape)

    def score(self, xyz):
        """Density score per residue; a heuristic, not a calibrated probability."""
        if self.profile is not None:
            return self.profile.score(self, xyz)
        values = np.clip(self.sample(xyz), -3., 6.)
        return np.average(values, axis=-1, weights=[1., 1., 1., .8])

    def canonical(self, xyz):
        frac = self.cell.fractionalize(gemmi.Position(*map(float, xyz)))
        positions = [np.mod(op.apply_to_xyz(list(frac)), 1.) for op in self.sg.operations()]
        return min(tuple(np.round(p, 6)) for p in positions)

    def distance(self, a, b):
        return self.cell.find_nearest_image(gemmi.Position(*map(float, a)),
                                           gemmi.Position(*map(float, b))).dist()

    def write_map(self, path):
        out = gemmi.Ccp4Map()
        out.grid = self.grid
        out.update_ccp4_header()
        out.write_ccp4_map(str(path))
