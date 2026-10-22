"""Train an opt-in patch discriminator on one unrelated reference only.

Defaults use CCP4's 1TQW; 1AJR is evaluated only after fitting and calibration.
Neither these defaults nor the API read the four SAD development references.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import gemmi
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild.data import Density, fingerprint
from phasebuild.backbone_discriminant import (features, patch_offsets, fit_discriminant,
                                              discrimination_summary)


def backbone_groups(path):
    structure = gemmi.read_structure(str(path))
    coords = []
    for chain in structure[0]:
        for residue in chain:
            atoms = [residue.find_atom(name, '*') for name in ['N', 'CA', 'C', 'O']]
            if any(atom is None for atom in atoms):
                continue
            xyz = np.array([list(atom.pos) for atom in atoms])
            if (1.2 < np.linalg.norm(xyz[0] - xyz[1]) < 1.7 and
                    1.2 < np.linalg.norm(xyz[2] - xyz[1]) < 1.8):
                coords.append(xyz)
    return np.array(coords)


def perturb_groups(xyz, rng):
    groups = [xyz + rng.normal(0., .08, (len(xyz), 1, 3))]
    labels = [np.ones(len(xyz), dtype=int)]
    for kind in ['rotation', 'translation', 'near_translation', 'random_location']:
        ca = xyz[:, 1:2]
        if kind == 'rotation':
            axis = rng.normal(size=(len(xyz), 3))
            axis /= np.linalg.norm(axis, axis=-1, keepdims=True)
            angle = rng.uniform(np.deg2rad(60.), np.pi, size=(len(xyz), 1))
            matrices = Rotation.from_rotvec(axis * angle).as_matrix()
            candidate = np.einsum('nij,naj->nai', matrices, xyz - ca) + ca
        else:
            direction = rng.normal(size=(len(xyz), 3))
            direction /= np.linalg.norm(direction, axis=-1, keepdims=True)
            bounds = {'translation': (2., 5.), 'near_translation': (1.2, 2.),
                      'random_location': (5., 15.)}[kind]
            displacement = direction * rng.uniform(*bounds, size=(len(xyz), 1))
            candidate = xyz + displacement[:, None]
        # Same-CA rotations about the N/C bisector can be almost unchanged.
        # Force those rare ambiguous cases to be genuine coordinate decoys.
        wrong = np.max(np.linalg.norm(candidate[:, :3] - xyz[:, :3], axis=-1), axis=-1) >= 1.1
        candidate[~wrong] += np.array([2., 0., 0.])
        groups.append(candidate)
        labels.append(np.zeros(len(xyz), dtype=int))
    return np.concatenate(groups), np.concatenate(labels), np.tile(np.arange(len(xyz)), len(groups))


def dataset(mtz_path, structure, output, resolution, phase_errors, limit, seed, offsets):
    xyz = backbone_groups(structure)
    if len(xyz) < 50:
        raise ValueError('Unrelated structure requires at least 50 valid backbone groups')
    rng = np.random.RandomState(seed)
    if len(xyz) > limit:
        xyz = xyz[np.sort(rng.choice(len(xyz), limit, replace=False))]
    rows, labels, owners, conditions = [], [], [], []
    for error in phase_errors:
        mtz = gemmi.read_mtz_file(str(mtz_path))
        data = np.array(mtz, copy=True)
        phase = mtz.column_labels().index('FC_BEST.F_phi.phi')
        data[:, phase] = (data[:, phase] + rng.normal(0., error, len(data))) % 360.
        mtz.set_data(data)
        path = output / ('augmented_%g_%g.mtz' % (resolution, error))
        mtz.write_to_file(str(path))
        density = Density(path, 'FC_BEST.F_phi.F', 'FC_BEST.F_phi.phi', resolution=resolution)
        candidates, truth, ids = perturb_groups(xyz, rng)
        rows.append(features(density, candidates, offsets))
        labels.append(truth)
        owners.append(ids)
        conditions.append({'phase_rms_degrees': float(error), 'rows': len(truth),
                           'map': fingerprint(path)})
    return np.concatenate(rows), np.concatenate(labels), np.concatenate(owners), conditions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    folder = '/Applications/ccp4-9/lib/data/reference_structures/'
    parser.add_argument('--training-structure', default=folder + 'reference-1tqw.pdb')
    parser.add_argument('--training-map', default=folder + 'reference-1tqw.mtz')
    parser.add_argument('--validation-structure', default=folder + 'reference-1ajr.pdb')
    parser.add_argument('--validation-map', default=folder + 'reference-1ajr.mtz')
    parser.add_argument('--resolution', type=float, default=2.8)
    parser.add_argument('--output', required=True)
    parser.add_argument('--limit', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=731)
    parser.add_argument('--regularization', type=float, default=20.)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    offsets = patch_offsets()
    x, y, owners, conditions = dataset(args.training_map, args.training_structure,
        out, args.resolution, [0., 35., 60., 85.], args.limit, args.seed, offsets)
    # All perturbed copies of a residue stay in the same split, across maps.
    holdout = owners % 5 == 0
    metadata = {'method': 'L2 logistic discrimination of local backbone versus synthetic coordinate decoys',
        'training_structure': fingerprint(args.training_structure),
        'training_map': fingerprint(args.training_map), 'resolution': args.resolution,
        'seed': args.seed, 'conditions': conditions,
        'training_rows': int((~holdout).sum()), 'calibration_rows': int(holdout.sum()),
        'calibration_rule': '95th percentile of held-out unrelated-reference decoy logits',
        'warning': 'Synthetic-task logits, not calibrated experimental probabilities; opt-in research only.'}
    model = fit_discriminant(x[~holdout], y[~holdout], offsets, metadata,
                             regularization=args.regularization)
    calibration = discrimination_summary(y[holdout], model.score_features(x[holdout]))
    thresholds = {key: value['logit'] for key, value in calibration['thresholds'].items()}
    model.metadata['calibration'] = calibration
    model.save(out / 'discriminant.npz')
    print(json.dumps({'calibration': calibration, 'fitting': model.metadata['optimization']}), flush=True)
    validation, truth, _, validation_conditions = dataset(args.validation_map, args.validation_structure,
        out, args.resolution, [20., 50., 80.], args.limit, args.seed + 1, offsets)
    report = dict(model.metadata, validation_structure=fingerprint(args.validation_structure),
                  validation_map=fingerprint(args.validation_map), validation_conditions=validation_conditions,
                  independent_validation=discrimination_summary(truth, model.score_features(validation), thresholds),
                  seconds=time.time() - started)
    (out / 'training_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'independent_validation': report['independent_validation'], 'seconds': report['seconds']}), flush=True)


if __name__ == '__main__':
    main()
