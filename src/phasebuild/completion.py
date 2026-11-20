"""Recover short internal gaps supported by an archived parent trace.

Both anchors must match one parent under one symmetry operation and lattice
translation. This does not invent loops or extend unanchored termini.
"""
import gemmi
import numpy as np
from .validation import peptide_metrics, peptide_valid


def register(fragment, parent, density, tolerance=.6):
    matches = []
    for start in range(len(parent)-len(fragment)+1):
        image = density.cell.find_nearest_image(gemmi.Position(*parent[start,1]),
                                                gemmi.Position(*fragment[0,1]))
        if image.dist() > tolerance:
            continue
        transformed = np.array([list(density.cell.orthogonalize(density.cell.fract_image(
            image,density.cell.fractionalize(gemmi.Position(*p))))) for p in fragment.reshape(-1,3)]).reshape(fragment.shape)
        delta = transformed[:,:3]-parent[start:start+len(fragment),:3]
        distances = np.linalg.norm(delta,axis=-1)
        rms = float(np.sqrt(np.mean(distances**2)))
        if rms <= tolerance and distances.max() <= 1.:
            matches.append((start,transformed,rms))
    return matches


def complete_internal_gaps(chains, parents, density, max_gap=2, minimum_score=2.):
    if max_gap < 0 or max_gap > 5:
        raise ValueError('max_gap must be between 0 and 5')
    chains = [c.copy() for c in chains]
    accepted = []
    while True:
        proposals = []
        for parent_index,parent in enumerate(parents):
            registrations = [register(c,parent,density) for c in chains]
            for i,left_options in enumerate(registrations):
                for j,right_options in enumerate(registrations):
                    if i == j:
                        continue
                    for left_start,left,left_rms in left_options:
                        for right_start,right,right_rms in right_options:
                            end = left_start+len(left)
                            gap = right_start-end
                            if not 0 <= gap <= max_gap:
                                continue
                            missing = parent[end:right_start].copy()
                            scores = density.score(missing) if gap else np.array([])
                            if gap and (scores.min() < minimum_score or density.sample(missing[:,:3]).min() < .5):
                                continue
                            merged = np.concatenate([left,missing,right])
                            new_links = [peptide_metrics(merged[k],merged[k+1])
                                         for k in range(len(left)-1,len(left)+gap)]
                            if not all(peptide_valid(m) for m in new_links):
                                continue
                            other_ca = [r[1] for k,c in enumerate(chains) if k not in [i,j] for r in c]
                            if any(density.distance(r[1],p) < 2.8 for r in missing for p in other_ca):
                                continue
                            if any(density.distance(merged[a,1],merged[b,1]) < 2.8
                                   for a in range(len(merged)) for b in range(a-1)):
                                continue
                            score = -(left_rms+right_rms) + (.01*float(scores.mean()) if gap else 0.)
                            proposals.append((score,i,j,merged,{'parent':parent_index+1,'gap_residues':gap,
                                  'left_anchor_rmsd':left_rms,'right_anchor_rmsd':right_rms,
                                  'missing_residue_scores':scores.tolist(),'new_links':new_links}))
        if not proposals:
            break
        proposals.sort(key=lambda item:item[0],reverse=True)
        _,i,j,merged,record = proposals[0]
        accepted.append(record)
        chains = [c for k,c in enumerate(chains) if k not in [i,j]] + [merged]
    return chains, {'method':'two-anchor recovery from archived parent traces',
                    'max_gap':max_gap,'minimum_score':minimum_score,'joins':accepted,
                    'added_residues':sum(r['gap_residues'] for r in accepted),
                    'scope':'No target reference or sequence used. Unsupported gaps remain unbuilt.'}
