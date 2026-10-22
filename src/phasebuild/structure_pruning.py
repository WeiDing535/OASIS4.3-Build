"""Structure-preserving experimental classifier pruning; no reference input."""
import string
import gemmi
import numpy as np
from .backbone_discriminant import supported_segments


def prune_structure(structure,density,classifier,minimum_length=4,bridge=1):
    """Retain complete original residues and split at classifier-removed runs.

    Residue names, seqids, MSE chemistry, sidechains and all atom attributes are
    copied. The first retained segment keeps its chain name. Further segments
    receive unused one-character IDs so the standard PDB/REFMAC path remains
    unambiguous. The input object is never mutated.
    """
    if len(structure)!=1 or not len(structure[0]):
        raise ValueError('Exactly one populated structure model is required')
    if minimum_length!=4 or bridge!=1:
        raise ValueError('This focused pruning protocol fixes minimum_length4 and bridge1')
    if not np.allclose(structure.cell.parameters,density.cell.parameters,rtol=1e-4,atol=1e-3):
        raise ValueError('Structure and scoring map unit cells differ')
    sg=gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
    if sg is None or sg.number!=density.sg.number:
        raise ValueError('Structure and scoring map space groups differ')
    threshold=float(classifier.metadata['calibration']['thresholds']['0.05']['logit'])
    result=structure.clone()
    for i in range(len(result[0])-1,-1,-1):
        del result[0][i]
    reserved={chain.name for chain in structure[0]}
    free_names=iter([name for name in string.ascii_uppercase+string.ascii_lowercase+string.digits
                     if name not in reserved])
    decisions=[];output_residues=0
    for chain in structure[0]:
        xyz=[]
        for residue in chain:
            atoms=[residue.find_atom(name,'*') for name in ['N','CA','C','O']]
            if not all(atom is not None for atom in atoms):
                raise ValueError('Complete N/CA/C/O required for every scored residue; input preserved')
            xyz.append([list(atom.pos) for atom in atoms])
        if not xyz:
            continue
        logits=classifier.score(density,np.asarray(xyz))
        intervals=supported_segments(logits,threshold,minimum_length=4,bridge=1)
        segments=[]
        for number,(a,b) in enumerate(intervals):
            if number==0:
                name=chain.name
                if len(name)>1:
                    raise ValueError('Original chain IDs must fit the one-character PDB refinement path')
            else:
                try:name=next(free_names)
                except StopIteration:raise ValueError('Too many split fragments for unique one-character PDB chain IDs')
            kept=gemmi.Chain(name)
            for residue in list(chain)[a:b]:
                kept.add_residue(residue.clone())
            result[0].add_chain(kept)
            output_residues+=b-a
            segments.append({'source_start_index':int(a),'source_stop_index':int(b),
                             'source_first_seqid':str(chain[a].seqid),
                             'source_last_seqid':str(chain[b-1].seqid),'output_chain':name})
        decisions.append({'source_chain':chain.name,'input_residues':len(chain),
                          'logits':logits.tolist(),'retained_segments':segments})
    if not output_residues:
        raise ValueError('No supported segments remain; original structure is preserved')
    report={'classifier_logit_threshold':threshold,'calibration_entry':'0.05',
            'minimum_length':4,'internal_weak_bridge':1,
            'input_residues':sum(len(c) for c in structure[0]),'output_residues':output_residues,
            'output_fragments':len(result[0]),'decisions':decisions,
            'preserved':'Whole retained residues including names, seqids, MSE, sidechains, B, occupancy, '
                        'altlocs and atom attributes. First retained chain segment preserves its name; '
                        'additional segments use new unused one-character chain IDs.',
            'warning':'Unrelated synthetic-decoy calibration is not an experimental false-positive rate'}
    return result,report
