"""Conservative, opt-in cleanup of extra protein fragments in one-copy models.

This is a coordinate-subset operation, not a phasing method. Never infer the
ASU copy count from sequence coverage or silently apply it during building.
"""
import gemmi


def _protein(residue):
    return gemmi.find_tabulated_residue(residue.name).is_amino_acid()


def cleanup_single_copy(structure, sequence, *, asu_copies):
    """Keep a unique exact full-sequence chain and preserve nonprotein residues.

    Require an explicitly declared single-copy ASU and exactly one complete
    sequence match. No map, reference coordinates or R factors are accepted.
    Return a clone and an audit; never mutate the supplied structure.
    """
    if isinstance(asu_copies, bool) or not isinstance(asu_copies, int) or asu_copies != 1:
        raise ValueError('Cleanup requires an explicitly declared single-copy ASU')
    if not isinstance(sequence, str):
        raise ValueError('Sequence must be a protein sequence string')
    sequence = ''.join(sequence.split()).upper()
    if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('A nonempty, unambiguous standard protein sequence is required')
    if len(structure) != 1 or not len(structure[0]):
        raise ValueError('Exactly one populated coordinate model is required')

    matches = []
    for index, chain in enumerate(structure[0]):
        residues = [r for r in chain if _protein(r)]
        letters = ''.join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                          for r in residues)
        complete = all(r.find_atom('CA', '*') is not None for r in residues)
        if letters == sequence and complete:
            matches.append(index)
    if len(matches) != 1:
        raise ValueError('Cleanup requires exactly one complete, exact full-sequence chain')

    keep = matches[0]
    result = structure.clone()
    removed = []
    for index in range(len(result[0]) - 1, -1, -1):
        if index == keep:
            continue
        chain = result[0][index]
        for ri in range(len(chain) - 1, -1, -1):
            residue = chain[ri]
            if _protein(residue):
                removed.append({'chain': chain.name, 'seqid': str(residue.seqid),
                                'resname': residue.name, 'atoms': len(residue)})
                del chain[ri]
        if not len(chain):
            del result[0][index]

    before = sum(_protein(r) for c in structure[0] for r in c)
    after = sum(_protein(r) for c in result[0] for r in c)
    if after != len(sequence):
        raise RuntimeError('Internal single-copy residue-count invariant failed')
    report = {'operation': 'explicit-single-copy-full-sequence-cleanup',
              'asu_copies': asu_copies, 'sequence_length': len(sequence),
              'matched_source_chain': structure[0][keep].name,
              'matched_source_chain_index': keep,
              'input_protein_residues': before, 'output_protein_residues': after,
              'removed_protein_residues': len(removed), 'removed': removed,
              'selection': 'Unique exact complete sequence chain; no reference/map/Rfree selection',
              'preservation': 'Input not mutated; retained atoms and nonprotein residues cloned unchanged',
              'limitation': 'Exact sequence coverage does not establish coordinate correctness; refine and validate separately'}
    return result, report
