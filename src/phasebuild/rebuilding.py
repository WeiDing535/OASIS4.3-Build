"""Map/B-factor supported fragment pruning and property transfer after joins."""
import gemmi
import numpy as np
from .chemistry import BACKBONE


def prune_weak_fragments(structure,density,max_length=3,maximum_profile=3.,
                         minimum_b=45.,minimum_relative_b=1.5):
    known=[a.b_iso for c in structure[0] for r in c if r.name!='UNK'
           for a in r if a.name in BACKBONE]
    baseline=float(np.median(known)) if len(known)>=48 else None
    removed=[];decisions=[]
    for ci,chain in enumerate(structure[0]):
        if not chain or len(chain)>max_length or any(r.name!='UNK' for r in chain):
            continue
        xyz=np.array([[list(r.find_atom(n,'*').pos) for n in ['N','CA','C','O']] for r in chain])
        profile=float(np.mean(density.score(xyz)))
        b=float(np.median([a.b_iso for r in chain for a in r if a.name in BACKBONE]))
        weak=(baseline is not None and profile<maximum_profile and b>=minimum_b and b>=minimum_relative_b*baseline)
        decisions.append({'chain':chain.name,'length':len(chain),'mean_profile':profile,'median_b':b,
                          'known_backbone_median_b':baseline,'removed':bool(weak)})
        if weak:
            removed.append(ci)
    for ci in reversed(removed):
        del structure[0][ci]
    return {'method':'conjunctive short UNK, weak map and elevated refined B-factor gate',
            'thresholds':{'maximum_length':max_length,'profile_below':maximum_profile,
                          'median_b_at_least':minimum_b,'relative_b_at_least':minimum_relative_b},
            'decisions':decisions,'removed_fragments':len(removed),
            'removed_residues':sum(d['length'] for d in decisions if d['removed']),
            'scope':'Development heuristic, not a calibrated error probability; no Rfree or target reference used.'}


def transfer_properties(previous,current,maximum_ca_distance=.6):
    """Transfer B factors and compatible sidechain seeds using one crystal image."""
    previous.setup_cell_images();cell=previous.cell
    residues=[r for c in previous[0] for r in c]
    used=set();matched=0
    for chain in current[0]:
        for residue in chain:
            ca=residue.find_atom('CA','*')
            options=[]
            for index,old in enumerate(residues):
                if index in used:
                    continue
                image=cell.find_nearest_image(ca.pos,old.find_atom('CA','*').pos)
                if image.dist()<=maximum_ca_distance:
                    options.append((image.dist(),index,old,image))
            if not options:
                continue
            _,index,old,image=min(options,key=lambda v:v[0]);used.add(index);matched+=1
            lookup={a.name:a for a in old}
            for atom in residue:
                if atom.name in lookup:
                    atom.b_iso=lookup[atom.name].b_iso
            if residue.name!=old.name:
                continue
            # One image for the whole residue; never choose each atom's closest image.
            for atom in old:
                if atom.name in BACKBONE or atom.element.name in ['H','D']:
                    continue
                copied=atom.clone()
                copied.pos=cell.orthogonalize(cell.fract_image(image,cell.fractionalize(atom.pos)))
                residue.add_atom(copied)
    return {'matched_previous_residues':matched,'maximum_ca_distance':maximum_ca_distance,
            'scope':'B factors transferred to matched backbone atoms; same-type sidechains are fitting seeds.'}
