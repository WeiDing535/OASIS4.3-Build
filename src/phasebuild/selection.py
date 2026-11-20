"""Bounded beam search over mutually conflicting backbone fragments.

This optimizes a map-derived surrogate, not reference correctness. It is not an
exact maximum-weight independent-set solver and does not join peptide chains.
"""
import numpy as np


def independent_set(weights, lengths, conflicts, capacity, width=256):
    # States: (score, residue count, selected bitset, blocked future bitset).
    states = [(0., 0, 0, 0)]
    for i, (weight, length) in enumerate(zip(weights, lengths)):
        expanded = list(states)
        for score, count, selected, blocked in states:
            if not (blocked >> i) & 1 and count + length <= capacity:
                expanded.append((score + weight, count + length, selected | (1 << i), blocked | conflicts[i]))
        # Future-equivalent states can safely retain only the best accumulated score.
        unique = {}
        future_mask = ~((1 << (i + 1)) - 1)
        for state in expanded:
            key = (state[1], state[3] & future_mask)
            if key not in unique or state[0] > unique[key][0]:
                unique[key] = state
        states = sorted(unique.values(), key=lambda s: (-s[0], s[1], s[2]))[:width]
    best = max(states, key=lambda s: s[0])
    return [i for i in range(len(weights)) if (best[2] >> i) & 1]


def select_fragments(density, candidates, capacity=65, cutoff=2.5, width=256):
    if density.profile is None:
        raise ValueError('Global profile selection requires a trained density profile')
    pieces = []
    for xyz in candidates:
        scores = density.score(xyz)
        # Split at weak residues instead of allowing a strong average to hide them.
        good = scores >= cutoff
        starts = np.flatnonzero(good & ~np.r_[False, good[:-1]])
        ends = np.flatnonzero(good & ~np.r_[good[1:], False]) + 1
        for start, end in zip(starts, ends):
            if end - start >= 3:
                pieces.append((float(np.sum(scores[start:end] - cutoff)), xyz[start:end]))
    pieces.sort(key=lambda item: item[0], reverse=True)
    # Limit search size explicitly; all input candidates remain in the saved archive.
    pieces = pieces[:160]
    conflicts = [0] * len(pieces)
    for i in range(len(pieces)):
        for j in range(i):
            if any(density.distance(a, b) < 2.5 for a in pieces[i][1][:,1] for b in pieces[j][1][:,1]):
                conflicts[i] |= 1 << j
                conflicts[j] |= 1 << i
    indices = independent_set([p[0] for p in pieces], [len(p[1]) for p in pieces], conflicts, capacity, width)
    return [pieces[i][1] for i in indices], {
        'method': 'bounded beam search on fragment conflict graph', 'beam_width': width,
        'score_cutoff': cutoff, 'input_candidates': len(candidates), 'eligible_pieces': len(pieces),
        'selected_pieces': len(indices), 'objective': float(sum(pieces[i][0] for i in indices)),
        'optimality': 'not guaranteed; optimizes a map-score surrogate',
        'scope': 'No sequence information, target reference, or peptide joining'}
