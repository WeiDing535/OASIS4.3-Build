"""Opt-in local sequence anchors with the searched maximum in the null.

A noisy trace need not have the same sequence register over its full length.
This module finds ungapped local intervals, tests them against composition
shuffles searched over *all* the same intervals, and splits only accepted
anchors. It never uses a target structure or a previously solved model.
"""
import numpy as np
from .residue_profile import ALPHABET
from .selection import independent_set
from .sequence_null import required_shuffles


def _codes(sequence):
    return np.array([ALPHABET.index(c) if c != 'X' else 20 for c in sequence], dtype=np.int32)


def local_maxima(scores, catalogs, minimum_length=8, residue_penalty=.5, details=False):
    """Maximum local sum over every sequence offset and eligible subinterval.

    ``catalogs`` is a list of arrays with shape (batch, sequence_length).
    The optional candidates have fields (score, record, delta, begin, end),
    where sequence_position = model_position + delta and ends are exclusive.
    Prefix minima give the exact variable-length maximum, including searches
    which begin or end within either the model or the sequence.
    """
    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 2 or scores.shape[1] != 20 or not np.isfinite(scores).all():
        raise ValueError('Expected finite residue scores with shape (length,20)')
    if minimum_length < 1 or not np.isfinite(residue_penalty) or residue_penalty < 0:
        raise ValueError('Invalid local sequence search parameters')
    if not catalogs:
        return (np.empty(0), []) if details else np.empty(0)
    batch = len(catalogs[0]); length = len(scores)
    if any(np.asarray(c).ndim != 2 or len(c) != batch for c in catalogs):
        raise ValueError('Sequence catalogs must have the same batch size')
    if details and batch != 1:
        raise ValueError('Candidate details require one unshuffled catalog')
    maxima = np.full(batch, -np.inf); candidates = []
    if length < minimum_length:
        return (maxima, candidates) if details else maxima
    padded = np.column_stack([scores, np.zeros(length)])
    model_positions = np.arange(length)
    for ri, catalog in enumerate(catalogs):
        sequence_length = catalog.shape[1]
        if sequence_length < minimum_length:
            continue
        deltas = np.arange(minimum_length-length, sequence_length-minimum_length+1)
        positions = deltas[:, None]+model_positions[None, :]
        valid = (positions >= 0) & (positions < sequence_length)
        codes = catalog[:, np.clip(positions, 0, sequence_length-1)]
        values = padded[model_positions[None, None, :], codes]-residue_penalty
        values = np.where(valid[None, :, :], values, 0.)
        prefixes = np.concatenate([np.zeros((batch, len(deltas), 1)), values.cumsum(axis=-1)], axis=-1)
        begin_valid = np.maximum(0, -deltas)
        end_valid = np.minimum(length, sequence_length-deltas)
        prefix_positions = np.arange(length+1)
        valid_starts = ((prefix_positions[None, :] >= begin_valid[:, None]) &
                        (prefix_positions[None, :] <= (end_valid-minimum_length)[:, None]))
        eligible_prefixes = np.where(valid_starts[None, :, :], prefixes, np.inf)
        minima = np.minimum.accumulate(eligible_prefixes, axis=-1)
        ends = np.arange(minimum_length, length+1)
        totals = prefixes[..., ends]-minima[..., ends-minimum_length]
        valid_ends = ((ends[None, :] >= (begin_valid+minimum_length)[:, None]) &
                      (ends[None, :] <= end_valid[:, None]))
        totals = np.where(valid_ends[None, :, :], totals, -np.inf)
        maxima = np.maximum(maxima, totals.max(axis=(1, 2)))
        if details:
            # Keep the best begin for every possible end and offset. Searching
            # additional begins cannot exceed this score, so the same null
            # maximum controls these candidate anchor scores.
            for di, delta in enumerate(deltas):
                for ei in np.flatnonzero(valid_ends[di]):
                    end = int(ends[ei]); stop = end-minimum_length+1
                    begin = int(np.argmin(eligible_prefixes[0, di, :stop]))
                    candidates.append((float(totals[0, di, ei]), ri, int(delta), begin, end))
    return (maxima, candidates) if details else maxima


def local_permutation_maxima(scores, sequences, trials, random, minimum_length=8,
                             residue_penalty=.5, batch_size=16):
    """Repeat the complete local search on each composition-preserving shuffle."""
    codes = [_codes(sequence) for sequence in sequences]
    length = len(scores)
    largest = max((len(c)+length for c in codes), default=1)
    block = max(1, min(batch_size, 1_000_000//max(1, largest*length)))
    result = []
    for start in range(0, trials, block):
        count = min(block, trials-start)
        # Preserve deterministic trial-major, then record-major random draws.
        decoys = [[random.permutation(code) for code in codes] for _ in range(count)]
        catalogs = [np.stack([trial[ri] for trial in decoys]) for ri in range(len(codes))]
        result.extend(local_maxima(scores, catalogs, minimum_length, residue_penalty).tolist())
    return np.array(result)


def assign_local_sequence(chains, density, target, records, minimum_length=8,
                          shuffles=499, alpha=.05, seed=0, residue_penalty=.5,
                          max_anchors_per_chain=3, minimum_margin=1.):
    """Return split backbone chains and a standard annotation-compatible report.

    The family threshold is alpha / eligible input chains. For each chain,
    *every* anchor is compared to maxima from the full local search, rather
    than a null for its subsequently selected length. Identical declared
    sequence copies provide capacity but are one sequence hypothesis.
    """
    if not records or minimum_length < 1 or max_anchors_per_chain < 1:
        raise ValueError('Invalid local sequence assignment configuration')
    copies = {}; representatives = []
    for ri, (_, sequence) in enumerate(records):
        if not sequence or set(sequence)-set(ALPHABET+'X'):
            raise ValueError('Invalid protein sequence')
        if sequence not in copies:
            representatives.append(ri)
        copies.setdefault(sequence, []).append(ri)
    unique_sequences = [records[ri][1] for ri in representatives]
    maximum_sequence = max(map(len, unique_sequences))
    eligible = sum(len(chain) >= minimum_length and maximum_sequence >= minimum_length for chain in chains)
    threshold = alpha/max(1, eligible)
    trials = required_shuffles(shuffles, eligible, alpha)
    random = np.random.RandomState(seed)
    chain_reports = []; proposals = []
    for ci, chain in enumerate(chains):
        audit = {'source_fragment_index': ci, 'length': len(chain), 'eligible': len(chain) >= minimum_length}
        chain_reports.append(audit)
        if len(chain) < minimum_length or maximum_sequence < minimum_length:
            audit.update(reason='too short or no sequence window', candidates_tested=0)
            continue
        scores = target.scores(density, chain)
        observed, candidates = local_maxima(scores, [_codes(s)[None, :] for s in unique_sequences],
                                            minimum_length, residue_penalty, True)
        null = local_permutation_maxima(scores, unique_sequences, trials, random,
                                        minimum_length, residue_penalty)
        null95 = float(np.quantile(null, .95))
        audit.update(best_local_score=float(observed[0]), null_95=null95,
                     best_tail_fraction=float((1+np.count_nonzero(null >= observed[0]))/(trials+1)),
                     corrected_tail_threshold=threshold, candidates_tested=len(candidates))
        candidates.sort(reverse=True)
        if not candidates or (1+np.count_nonzero(null >= candidates[0][0]))/(trials+1) > threshold:
            audit['reason'] = 'insufficient local sequence evidence'
            continue
        array = np.array(candidates)
        kept = []
        for score, uri, delta, begin, end in candidates:
            tail = float((1+np.count_nonzero(null >= score))/(trials+1))
            if tail > threshold:
                break
            # A competitor in another register over the same model interval
            # is ambiguous. End-trim variants of one register are synonymous.
            overlaps = np.maximum(0, np.minimum(array[:, 4], end)-np.maximum(array[:, 3], begin))
            competing = ((array[:, 1] != uri) | (array[:, 2] != delta)) & (overlaps >= .5*(end-begin))
            runner = float(array[competing, 0].max()) if competing.any() else -np.inf
            margin = score-runner
            if margin < minimum_margin or any(max(begin, a) < min(end, b) for a, b in kept):
                continue
            kept.append((begin, end))
            ri = representatives[uri]
            for copy_ri in copies[records[ri][1]]:
                proposals.append({'source_fragment_index': ci, 'source_begin': begin, 'source_end': end,
                                  'sequence_record': copy_ri, 'start': begin+delta+1, 'end': end+delta,
                                  'length': end-begin, 'best_score': score, 'runner_up_margin': margin,
                                  'shuffle_tail_fraction': tail, 'corrected_tail_threshold': threshold,
                                  'weight': max(.01, score-null95)})
            if len(kept) >= max_anchors_per_chain:
                break
        audit['reason'] = 'local candidates available' if kept else 'local register is ambiguous'
    conflicts = [0]*len(proposals)
    for i, a in enumerate(proposals):
        for j, b in enumerate(proposals[:i]):
            model_overlap = (a['source_fragment_index'] == b['source_fragment_index'] and
                             max(a['source_begin'], b['source_begin']) < min(a['source_end'], b['source_end']))
            sequence_overlap = (a['sequence_record'] == b['sequence_record'] and
                                max(a['start'], b['start']) <= min(a['end'], b['end']))
            if model_overlap or sequence_overlap:
                conflicts[i] |= 1 << j; conflicts[j] |= 1 << i
    chosen = independent_set([p['weight'] for p in proposals], [p['length'] for p in proposals],
                             conflicts, sum(len(sequence) for _, sequence in records))
    selected = [proposals[i] for i in chosen]
    output = []; decisions = []
    for ci, chain in enumerate(chains):
        anchors = sorted((p for p in selected if p['source_fragment_index'] == ci), key=lambda p: p['source_begin'])
        cursor = 0
        for anchor in anchors:
            begin, end = anchor['source_begin'], anchor['source_end']
            if begin > cursor:
                output.append(np.array(chain[cursor:begin], copy=True))
                decisions.append({'fragment_index': len(output)-1, 'source_fragment_index': ci,
                                  'source_begin': cursor, 'source_end': begin, 'length': begin-cursor,
                                  'accepted': False, 'reason': 'outside accepted local anchor'})
            output.append(np.array(chain[begin:end], copy=True))
            decision = {k: v for k, v in anchor.items() if k != 'weight'}
            decision.update(fragment_index=len(output)-1, accepted=True,
                            record_id=records[anchor['sequence_record']][0],
                            reason='local permutation, register and nonoverlap checks passed')
            decisions.append(decision); cursor = end
        if cursor < len(chain):
            output.append(np.array(chain[cursor:], copy=True))
            decisions.append({'fragment_index': len(output)-1, 'source_fragment_index': ci,
                              'source_begin': cursor, 'source_end': len(chain), 'length': len(chain)-cursor,
                              'accepted': False, 'reason': 'no accepted local anchor'})
    known = sum(sum(c != 'X' for c in records[d['sequence_record']][1][d['start']-1:d['end']])
                for d in decisions if d['accepted'])
    report = {'decisions': decisions, 'chain_searches': chain_reports, 'assigned_residues': known,
              'assigned_catalog_positions': sum(d['length'] for d in decisions if d['accepted']),
              'identical_sequence_copy_groups': list(copies.values()), 'requested_shuffles': shuffles,
              'shuffles': trials, 'minimum_attainable_tail_fraction': 1/(trials+1),
              'familywise_tail_threshold': threshold, 'alpha': alpha, 'seed': seed,
              'minimum_length': minimum_length, 'residue_penalty': residue_penalty,
              'maximum_anchors_per_chain': max_anchors_per_chain, 'minimum_margin': minimum_margin,
              'sequence_records': [{'id': name, 'length': len(seq)} for name, seq in records],
              'null_search': 'Maximum over all eligible contiguous model subsegments and all sequence offsets; identical search for each composition shuffle.',
              'scope': 'Composition-shuffle tails are heuristic evidence, not calibrated correctness probabilities. No target reference used.'}
    return output, report
