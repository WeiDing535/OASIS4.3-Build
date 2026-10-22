"""Reference access is confined to this independent evaluation module."""
import gemmi
import numpy as np
from scipy.optimize import linear_sum_assignment
from .data import fingerprint


def ca_atoms(structure):
    result = []
    if len(structure) == 0:
        return result
    for chain in structure[0]:
        for residue in chain:
            atom = residue.find_atom('CA', '*')
            if atom is not None:
                result.append((chain.name, str(residue.seqid), residue.name, atom.pos, residue))
    return result


def evaluate(model_path, reference_path, cutoff=1.):
    if cutoff <= 0:
        raise ValueError('Matching cutoff must be positive')
    model = gemmi.read_structure(str(model_path))
    reference = gemmi.read_structure(str(reference_path))
    if model.spacegroup_hm != reference.spacegroup_hm:
        a, b = gemmi.find_spacegroup_by_name(model.spacegroup_hm), gemmi.find_spacegroup_by_name(reference.spacegroup_hm)
        if a is None or b is None or a.number != b.number:
            raise ValueError('Model and reference space groups differ')
    if not np.allclose(model.cell.parameters, reference.cell.parameters, rtol=1e-4, atol=1e-3):
        raise ValueError('Model and reference unit cells differ')
    built, target = ca_atoms(model), ca_atoms(reference)
    if not target:
        raise ValueError('Reference contains no C-alpha atoms')
    reference.setup_cell_images()
    cell = reference.cell
    distances = np.empty((len(built), len(target)))
    for i, item in enumerate(built):
        for j, other in enumerate(target):
            distances[i, j] = cell.find_nearest_image(item[3], other[3]).dist()
    assignments = []
    if len(built):
        # Invalid pairs cost 1. Valid costs sum to <1: maximize count first.
        cost = np.where(distances < cutoff, distances / (cutoff * (len(built) + len(target) + 1)), 1.)
        rows, cols = linear_sum_assignment(cost)
        assignments = [(int(i), int(j), float(distances[i, j]))
                       for i, j in zip(rows, cols) if distances[i, j] < cutoff]
    matches = []
    backbone_correct = 0
    sequence_correct = 0
    sequence_and_backbone_correct = 0
    assigned_matched = 0
    image_keys = {}
    for i, j, distance in assignments:
        image = cell.find_nearest_image(built[i][3], target[j][3])
        image_keys[i] = (image.sym_idx,tuple(image.pbc_shift))
        atom_distances = {}
        for name in ['N', 'CA', 'C', 'O']:
            a, b = built[i][4].find_atom(name, '*'), target[j][4].find_atom(name, '*')
            if a is not None and b is not None:
                transformed = cell.orthogonalize(cell.fract_image(image, cell.fractionalize(b.pos)))
                atom_distances[name] = a.pos.dist(transformed)
        correct = all(name in atom_distances and atom_distances[name] < cutoff for name in ['N', 'CA', 'C'])
        backbone_correct += int(correct)
        letter = gemmi.find_tabulated_residue(built[i][2]).one_letter_code.upper()
        target_letter = gemmi.find_tabulated_residue(target[j][2]).one_letter_code.upper()
        assigned = letter in 'ACDEFGHIKLMNPQRSTVWY'
        type_correct = assigned and letter == target_letter
        assigned_matched += int(assigned)
        sequence_correct += int(type_correct)
        sequence_and_backbone_correct += int(type_correct and correct)
        matches.append({'built_chain': built[i][0], 'built_residue': built[i][1],
                        'reference_chain': target[j][0], 'reference_residue': target[j][1],
                        'ca_distance': distance, 'atom_distances': atom_distances,
                        'built_residue_name':built[i][2],'reference_residue_name':target[j][2],
                        'sequence_assigned':assigned,'sequence_type_correct':type_correct,
                        'n_ca_c_all_within_cutoff': correct})
    indexed = {i: j for i, j, _ in assignments}
    links = {'both_ca_matched': 0, 'reference_forward': 0, 'reference_reverse': 0, 'other': 0}
    for i in range(len(built) - 1):
        if built[i][0] != built[i + 1][0] or i not in indexed or i + 1 not in indexed:
            continue
        a, b = indexed[i], indexed[i + 1]
        links['both_ca_matched'] += 1
        if target[a][0] == target[b][0] and b == a + 1:
            links['reference_forward'] += 1
        elif target[a][0] == target[b][0] and b == a - 1:
            links['reference_reverse'] += 1
        else:
            links['other'] += 1
    chain_registration = []
    for chain_name in dict.fromkeys(item[0] for item in built):
        indices = [i for i in indexed if built[i][0]==chain_name]
        reference_chains = sorted({target[indexed[i]][0] for i in indices})
        transforms = {image_keys[i] for i in indices}
        positions = [indexed[i] for i in sorted(indices)]
        chain_registration.append({'built_chain':chain_name,'matched_ca':len(indices),
            'reference_chains':reference_chains,'distinct_symmetry_and_translation_images':len(transforms),
            'one_reference_chain_and_image':bool(indices) and len(reference_chains)==1 and len(transforms)==1,
            'matched_reference_order_strictly_increasing':bool(positions) and all(b>a for a,b in zip(positions,positions[1:]))})
    assigned_total = sum(gemmi.find_tabulated_residue(r[2]).one_letter_code.upper() in
                         'ACDEFGHIKLMNPQRSTVWY' for r in built)
    return {'model': fingerprint(model_path), 'reference': fingerprint(reference_path),
            'cutoff_angstrom': cutoff, 'built_residues': len(built),
            'reference_residues': len(target), 'matched_residues': len(matches),
            'ca_completeness': len(matches) / len(target),
            'ca_precision': len(matches) / len(built) if built else 0.,
            'n_ca_c_correct_residues': backbone_correct,
            'n_ca_c_completeness': backbone_correct / len(target),
            'n_ca_c_precision': backbone_correct / len(built) if built else 0.,
            'matched_links': links,
            'chain_registration':chain_registration,
            'sequence_evaluation':{'assigned_residues':assigned_total,'assigned_with_ca_match':assigned_matched,
                'correct_type_with_ca_match':sequence_correct,
                'correct_type_and_backbone':sequence_and_backbone_correct,
                'identity_accuracy_among_assigned_ca_matches':sequence_correct/assigned_matched if assigned_matched else None,
                'joint_sequence_backbone_accuracy_among_all_assigned':sequence_and_backbone_correct/assigned_total if assigned_total else None,
                'scope':'Unknown residues are not guesses. Unmatched assigned residues cannot be verified by this CA-cutoff metric.'},
            'matched_ca_rmsd': float(np.sqrt(np.mean([d*d for _, _, d in assignments]))) if matches else None,
            'matches': matches,
            'metric_scope': 'One-to-one CA assignment allowing crystallographic symmetry and lattice translations. '
                            'N/CA/C correctness uses the same symmetry image selected for CA. '
                            'Link counts test reference residue order, not full peptide geometry. '
                            'Sequence identities are checked only on spatially matched residues. '
                            'Chain registration reports consistency of matched images; no arbitrary-origin search. '
                            'Matched-only RMSD is not an all-residue RMSD.'}
