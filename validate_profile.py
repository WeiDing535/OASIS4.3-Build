"""Synthetic orientation discrimination on an unrelated validation structure."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
import gemmi
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import rankdata
from phasebuild.data import Density, fingerprint
from phasebuild.profile import DensityProfile


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--profile', required=True)
    p.add_argument('--mtz', required=True)
    p.add_argument('--structure', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--resolution', type=float, default=1.9)
    p.add_argument('--seed', type=int, default=1)
    args = p.parse_args()
    density = Density(args.mtz, 'FC_BEST.F_phi.F', 'FC_BEST.F_phi.phi', resolution=args.resolution)
    density.profile = DensityProfile.load(args.profile)
    groups = []
    for chain in gemmi.read_structure(args.structure)[0]:
        for residue in chain:
            atoms = [residue.find_atom(name,'*') for name in ['N','CA','C','O']]
            if all(atom is not None for atom in atoms):
                groups.append([list(atom.pos) for atom in atoms])
    if not groups:
        raise ValueError('No complete backbone groups for validation')
    xyz = np.array(groups)
    ca = xyz[:,1:2]
    negative = np.einsum('rij,raj->rai', Rotation.random(len(xyz), random_state=args.seed).as_matrix(), xyz-ca)+ca
    positive_scores, negative_scores = density.score(xyz), density.score(negative)
    ranks = rankdata(np.r_[positive_scores,negative_scores])
    n = len(xyz)
    auc = (ranks[:n].sum()-n*(n+1)/2)/(n*n)
    result = {'validation_structure': fingerprint(args.structure), 'map': fingerprint(args.mtz),
              'profile': fingerprint(args.profile), 'seed': args.seed, 'residues': n,
              'quantiles': [.05,.1,.5,.9,.95],
              'positive_scores': np.quantile(positive_scores,[.05,.1,.5,.9,.95]).tolist(),
              'random_rotation_negative_scores': np.quantile(negative_scores,[.05,.1,.5,.9,.95]).tolist(),
              'auc_vs_random_rotation': float(auc),
              'scope': 'Synthetic orientation test at known CA sites. Not a model-building benchmark.'}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
