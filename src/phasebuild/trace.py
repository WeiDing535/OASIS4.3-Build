"""Independent map-only rigid-tripeptide search and forward beam growth.

This intentionally small research prototype does not assign sequence or refine.
"""
import time
import json
import gemmi
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from .geometry import CONFORMATIONS, template, extend, oxygen, background_probes, flexible_template, prepend


def peaks(density, threshold=1.):
    values = np.array(density.grid, copy=False)
    mask = (values == maximum_filter(values, size=3, mode='wrap')) & (values > threshold)
    indices = np.argwhere(mask)
    indices = indices[np.argsort(values[tuple(indices.T)])[::-1]]
    unique, points = set(), []
    for idx in indices:
        frac = gemmi.Fractional(*(idx / np.array(values.shape)))
        pos = density.cell.orthogonalize(frac)
        xyz = np.array(list(pos))
        key = density.canonical(xyz)
        if key not in unique:
            unique.add(key)
            points.append(xyz)
    return np.array(points)


def find_seeds(density, rotations=512, refine_count=96, seed=0, progress=print):
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
            scores = density.score(xyz if density.profile is None else xyz[..., 1:2, :, :]).mean(axis=-1)
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
    for candidate in candidates:
        key = tuple(candidate[1])
        if counts.get(key, 0) < 2:
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


def grow(density, xyz, last_psi, max_length=30, beam_width=4, min_score=.9):
    # Each state carries its best prefix: a bad extension cannot force acceptance.
    beams = [(0., xyz.copy(), last_psi)]
    best = (0., xyz.copy(), last_psi)
    # A compact, general Ramachandran proposal set. Cis peptides are unsupported.
    proposals = [(psi, phi) for phi, psi in CONFORMATIONS]
    proposals += [(psi, phi) for phi in [-150., -120., -90., -60.]
                  for psi in [-60., -30., 0., 120., 150., 180.]]
    proposals = list(dict.fromkeys(proposals))
    for _ in range(max(0, max_length - len(xyz))):
        branches = []
        for cumulative, chain, _ in beams:
            variants = []
            for psi, phi in proposals:
                residue = extend(chain[-1], psi, phi)
                previous = chain[-1].copy()
                previous[3] = oxygen(*previous[:3], psi)
                # Test the new carbonyl in the next cycle; score its N/CA/C now.
                samples = density.sample(residue[:3])
                quality = (float(density.score(residue)) if density.profile is not None
                           else float(np.clip(samples, -3., 6.).mean()))
                peptide_score = float(density.score(previous))
                if (quality < min_score or peptide_score < min_score or
                        samples.min() < .5 * min_score or
                        density.sample(previous[3]) < .5 * min_score):
                    continue
                if len(chain) > 1:
                    # Local nonbonded clash check (same unwrapped chain).
                    delta = residue[:3, None, :] - chain[:-1, :3].reshape(1, -1, 3)
                    if np.min(np.linalg.norm(delta, axis=-1)) < 2.1:
                        continue
                    # Detect closure onto own symmetry copies, excluding neighbour.
                    if any(density.distance(residue[1], ca) < 2.8 for ca in chain[:-1, 1]):
                        continue
                variants.append((quality + .25 * peptide_score, residue, previous, psi))
            variants.sort(key=lambda x: x[0], reverse=True)
            for quality, residue, previous, psi in variants[:beam_width]:
                # Refine the outgoing psi, incoming phi and new terminal psi.
                # This relaxes torsion discretization while retaining ideal bonds.
                # Recover the starting phi by testing the proposal coordinates.
                phi = min(proposals, key=lambda pair: np.linalg.norm(
                    extend(chain[-1], pair[0], pair[1])[:3] - residue[:3]))[1]
                def proposal(params):
                    new_residue = extend(chain[-1], params[0], params[1])
                    new_residue[3] = oxygen(*new_residue[:3], params[2])
                    prev = chain[-1].copy()
                    prev[3] = oxygen(*prev[:3], params[0])
                    return new_residue, prev
                def loss(params):
                    next_residue, prev = proposal(params)
                    if density.profile is not None:
                        return -float(.8 * density.score(next_residue) + .2 * density.score(prev))
                    points = np.concatenate([next_residue, prev[3:4]])
                    return -float(np.clip(density.sample(points), -3., 6.).mean())
                fit = minimize(loss, [psi, phi, -40.], method='L-BFGS-B',
                               bounds=[(psi - 20., psi + 20.), (phi - 20., phi + 20.), (-180., 180.)],
                               options={'maxiter': 18, 'eps': .15, 'ftol': 1e-5})
                residue, previous = proposal(fit.x)
                # Re-check constraints after torsion optimization.
                if np.min(density.sample(residue[:3])) < .5 * min_score:
                    continue
                if len(chain) > 1:
                    delta = residue[:3, None, :] - chain[:-1, :3].reshape(1, -1, 3)
                    if np.min(np.linalg.norm(delta, axis=-1)) < 2.1:
                        continue
                    if any(density.distance(residue[1], ca) < 2.8 for ca in chain[:-1, 1]):
                        continue
                quality = -1.25 * float(loss(fit.x))
                if quality < 1.25 * min_score:
                    continue
                new = np.concatenate([chain[:-1], previous[None], residue[None]])
                gain = quality - 1.25 * min_score
                branches.append((cumulative + gain, new, float(fit.x[2])))
        if not branches:
            break
        branches.sort(key=lambda x: x[0], reverse=True)
        beams = branches[:beam_width]
        if beams[0][0] > best[0]:
            best = beams[0]
    return best[1]


def grow_backward(density, xyz, max_length, beam_width, min_score, used):
    beams = [(0., xyz.copy())]
    best = beams[0]
    for _ in range(max(0, max_length - len(xyz))):
        branches = []
        for cumulative, chain in beams:
            variants = []
            for phi in [-150., -120., -90., -63., 60.]:
                for psi in [-60., -42., 0., 40., 120., 150., 180.]:
                    residue = prepend(chain[0], phi, psi)
                    score = float(density.score(residue))
                    variants.append((score, phi, psi))
            variants.sort(reverse=True)
            for _, phi, psi in variants[:beam_width]:
                def loss(params):
                    return -float(density.score(prepend(chain[0], *params)))
                fit = minimize(loss, [phi, psi], method='L-BFGS-B',
                               bounds=[(phi - 25., phi + 25.), (psi - 25., psi + 25.)],
                               options={'maxiter': 20, 'eps': .15, 'ftol': 1e-5})
                residue = prepend(chain[0], *fit.x)
                quality = -float(fit.fun)
                if quality < min_score or np.min(density.sample(residue)) < .5 * min_score:
                    continue
                delta = residue[:3, None, :] - chain[1:, :3].reshape(1, -1, 3)
                if np.min(np.linalg.norm(delta, axis=-1)) < 2.1:
                    continue
                if any(density.distance(residue[1], ca) < 2.8 for ca in list(chain[1:, 1]) + used):
                    continue
                new = np.concatenate([residue[None], chain])
                branches.append((cumulative + quality - min_score, new))
        if not branches:
            break
        branches.sort(key=lambda x: x[0], reverse=True)
        beams = branches[:beam_width]
        if beams[0][0] > best[0]:
            best = beams[0]
    return best[1]


def write_backbone(path, chains, density):
    structure = gemmi.Structure()
    structure.name = 'phasebuild_unassigned_backbone'
    structure.cell = density.cell
    structure.spacegroup_hm = density.sg.hm
    model = gemmi.Model('1')
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    for i, xyz in enumerate(chains):
        chain = gemmi.Chain(alphabet[i] if i < len(alphabet) else str(i + 1))
        for j, points in enumerate(xyz):
            residue = gemmi.Residue()
            residue.name = 'UNK'
            residue.seqid = gemmi.SeqId(j + 1, ' ')
            for name, element, point in zip(['N', 'CA', 'C', 'O'], ['N', 'C', 'C', 'O'], points):
                atom = gemmi.Atom()
                atom.name, atom.element = name, gemmi.Element(element)
                atom.pos = gemmi.Position(*map(float, point))
                atom.occ, atom.b_iso = 1., 30.
                residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.make_mmcif_document().write_file(str(path))


def trace(density, rotations=512, refine_count=96, max_length=30,
          beam_width=4, max_residues=65, min_score=.9, seed=0, progress=print,
          candidate_output=None, selection_cutoff=None):
    started = time.perf_counter()
    seeds, stats = find_seeds(density, rotations, refine_count, seed, progress)
    chains = []
    seed_scores = []
    for score, xyz, psi in seeds:
        if np.min(density.score(xyz)) < min_score:
            continue
        used = [ca for chain in chains for ca in chain[:, 1]]
        if used and any(density.distance(ca, old) < 2.5 for ca in xyz[:, 1] for old in used):
            continue
        remaining = max_residues - sum(len(c) for c in chains)
        if remaining < 3:
            break
        chain = grow(density, xyz, psi, min(max_length, remaining), beam_width, min_score)
        # Stop the proposed chain at its first overlap with accepted fragments.
        keep = len(chain)
        for j in range(3, len(chain)):
            if any(density.distance(chain[j, 1], old) < 2.8 for old in used):
                keep = j
                break
        chain = chain[:keep]
        chain = grow_backward(density, chain, min(max_length, remaining), beam_width, min_score, used)
        chains.append(chain)
        seed_scores.append(score)
        progress('Accepted fragment %d: %d residues, seed score %.3f' % (len(chains), len(chain), score))
    if candidate_output is not None:
        np.savez_compressed(candidate_output,
                            seeds=np.array([s[1] for s in seeds]).reshape(-1,3,4,3),
                            seed_scores=np.array([s[0] for s in seeds]),
                            seed_psi=np.array([s[2] for s in seeds]),
                            grown=np.concatenate(chains) if chains else np.empty((0,4,3)),
                            grown_lengths=np.array([len(c) for c in chains], dtype=int),
                            metadata=json.dumps({'input': density.source,
                                                 'profile_training': density.profile.metadata if density.profile else None}))
    if selection_cutoff is not None:
        from .selection import select_fragments
        stats['before_selection_residues'] = sum(map(len, chains))
        chains, stats['selection'] = select_fragments(density, chains + [s[1] for s in seeds],
                                                      max_residues, selection_cutoff)
        progress('Global selection retained %d residues in %d fragments' % (sum(map(len,chains)),len(chains)))
    stats.update(elapsed_seconds=time.perf_counter() - started,
                 fragments=len(chains), residues=sum(map(len, chains)),
                 fragment_lengths=[len(c) for c in chains], seed_scores=seed_scores,
                 limitations=['unassigned UNK residues', 'growth proposals are greedy; no global joining',
                              'no sidechains, refinement or sequence assignment',
                              'density score is heuristic; false fragments are possible'])
    return chains, stats
