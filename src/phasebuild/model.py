"""Backbone model IO; no target-reference access."""
import gemmi
import numpy as np


def read_backbone(path, density):
    structure = gemmi.read_structure(str(path))
    if not np.allclose(structure.cell.parameters, density.cell.parameters, rtol=1e-4, atol=1e-3):
        raise ValueError('Model and map unit cells differ')
    sg = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
    if sg is None or sg.number != density.sg.number:
        raise ValueError('Model and map space groups differ')
    chains = []
    for chain in structure[0]:
        residues = []
        for residue in chain:
            atoms = [residue.find_atom(name, '*') for name in ['N','CA','C','O']]
            if not all(a is not None for a in atoms):
                raise ValueError('This autonomous improvement stage requires N/CA/C/O in every residue')
            residues.append([list(a.pos) for a in atoms])
        if residues:
            chains.append(np.array(residues, dtype=float))
    return chains
