"""Three-residue contextual seed search with bounded-memory profile sampling.

Uses the complete candidate peptide context during coarse orientation ranking;
no reference model or sequence hypothesis is an input.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from .geometry import CONFORMATIONS,template,background_probes,flexible_template
from .profile import frames
from .trace import peaks


def score_oriented_profiles(density,oriented,centers,block=128):
    if density.profile is None:raise ValueError('Context search requires a trained profile')
    p=density.profile;ca,axes=frames(oriented)
    patches=np.einsum('...ij,pj->...pi',axes,p.offsets)+ca[...,None,:]
    output=np.empty((len(centers),len(oriented)))
    for start in range(0,len(centers),block):
        center=centers[start:start+block]
        values=density.sample(center[:,None,None,None,:]+patches[None])
        correlation=p.correlation(values)
        atoms=np.clip(density.sample(center[:,None,None,None,:]+oriented[None]),-3.,6.).mean(axis=-1)
        scores=4.*correlation+.15*atoms
        output[start:start+len(center)]=scores.mean(axis=-1)-.15*scores.std(axis=-1)
    return output


def find_context_seeds(density, rotations=512, refine_count=96, seed=0, progress=print):
    centers = peaks(density)
    if len(centers) == 0:
        return [], {'peak_count': 0, 'coarse_candidates': 0}
    matrices = Rotation.random(rotations, random_state=seed).as_matrix()
    candidates = []
    background_weight = .55 if density.profile is None else 0.
    for phi, psi in CONFORMATIONS:
        base = template(phi, psi)
        probes = background_probes(base)
        best = np.full((len(centers), 2), -np.inf)
        best_rot = np.zeros((len(centers), 2), dtype=int)
        for start in range(0, rotations, 16):
            rot = matrices[start:start + 16]
            oriented = np.einsum('rij,akj->raki', rot, base)
            xyz = centers[:, None, None, None, :] + oriented[None, ...]
            scores = score_oriented_profiles(density,oriented,centers)
            if background_weight:
                background = centers[:, None, None, :] + np.einsum('rij,aj->rai', rot, probes)[None, ...]
                scores -= background_weight * np.clip(density.sample(background), -3., 6.).mean(axis=-1)
            merged = np.concatenate([best, scores], axis=1)
            ids = np.broadcast_to(np.arange(start, start + len(rot)), scores.shape)
            merged_rot = np.concatenate([best_rot, ids], axis=1)
            top = np.argsort(merged, axis=1)[:, -2:]
            best = np.take_along_axis(merged, top, axis=1)
            best_rot = np.take_along_axis(merged_rot, top, axis=1)
        for i, center in enumerate(centers):
            for j in range(2):
                candidates.append((float(best[i, j]), center, matrices[best_rot[i, j]], base, probes, phi, psi))
    candidates.sort(key=lambda x: x[0], reverse=True)
    # Limit domination by repeated orientations at one density maximum.
    chosen, counts = [], {}
    chosen_centers=[]
    for candidate in candidates:
        key = tuple(candidate[1])
        if counts.get(key,0)==0 and any(density.distance(candidate[1],c)<2.2 for c in chosen_centers):
            continue
        if counts.get(key, 0) < 2:
            if counts.get(key,0)==0:chosen_centers.append(candidate[1])
            chosen.append(candidate)
            counts[key] = counts.get(key, 0) + 1
        if len(chosen) >= refine_count:
            break
    refined = []
    for _, center, rotation, base, probes, phi, psi in chosen:
        def coordinates(params):
            matrix = Rotation.from_rotvec(params[3:]).as_matrix() @ rotation
            return base @ matrix.T + center + params[:3]
        def loss(params):
            scores = density.score(coordinates(params))
            matrix = Rotation.from_rotvec(params[3:]).as_matrix() @ rotation
            background = probes @ matrix.T + center + params[:3]
            return -float(scores.mean() - background_weight * np.clip(density.sample(background), -3., 6.).mean())
        fit = minimize(loss, np.zeros(6), method='L-BFGS-B',
                       bounds=[(-1.2, 1.2)] * 3 + [(-.6, .6)] * 3,
                       options={'maxiter': 45, 'eps': 1e-3, 'ftol': 1e-7})
        rigid_matrix = Rotation.from_rotvec(fit.x[3:]).as_matrix() @ rotation
        rigid_center = center + fit.x[:3]
        start_angles = np.array([psi, phi, psi, phi, psi, 0.])
        def flex_coordinates(params):
            matrix = Rotation.from_rotvec(params[3:6]).as_matrix() @ rigid_matrix
            return flexible_template(start_angles + np.rad2deg(params[6:])) @ matrix.T + rigid_center + params[:3]
        def flex_loss(params):
            xyz = flex_coordinates(params)
            matrix = Rotation.from_rotvec(params[3:6]).as_matrix() @ rigid_matrix
            background = probes @ matrix.T + rigid_center + params[:3]
            return -float(density.score(xyz).mean() - background_weight * np.clip(density.sample(background), -3., 6.).mean())
        bounds = ([(-.6, .6)] * 3 + [(-.25, .25)] * 3
                  + [(-np.deg2rad(35.), np.deg2rad(35.))] * 4
                  + [(-np.pi, np.pi)] * 2)
        fitted = minimize(flex_loss, np.zeros(12), method='L-BFGS-B', bounds=bounds,
                          options={'maxiter': 40, 'eps': 1e-3, 'ftol': 1e-6})
        xyz = flex_coordinates(fitted.x)
        refined.append((-float(flex_loss(fitted.x)), xyz, float(start_angles[4] + np.rad2deg(fitted.x[10]))))
        if len(refined) % 32 == 0:
            progress('Fitted %d/%d seeds' % (len(refined), len(chosen)))
    refined.sort(key=lambda x: x[0], reverse=True)
    progress('Refined %d seed candidates from %d symmetry-unique peaks' % (len(refined), len(centers)))
    return refined, {'peak_count': len(centers),
                     'coarse_candidates': len(centers) * rotations * len(CONFORMATIONS),
                     'refined_candidates': len(refined)}

