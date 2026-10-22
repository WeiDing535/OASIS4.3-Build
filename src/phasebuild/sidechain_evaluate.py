"""Independent atom-level reference evaluation; never imported by builders."""
import gemmi
import numpy as np
from .evaluate import evaluate
from .chemistry import BACKBONE
from .data import fingerprint

SYMMETRIES={'ASP':{'OD1':'OD2','OD2':'OD1'},'GLU':{'OE1':'OE2','OE2':'OE1'},
            'ARG':{'NH1':'NH2','NH2':'NH1'},'VAL':{'CG1':'CG2','CG2':'CG1'},
            'LEU':{'CD1':'CD2','CD2':'CD1'},
            'PHE':{'CD1':'CD2','CD2':'CD1','CE1':'CE2','CE2':'CE1'},
            'TYR':{'CD1':'CD2','CD2':'CD1','CE1':'CE2','CE2':'CE1'}}


def side_atoms(residue):
    return {a.name:a for a in residue if a.name not in BACKBONE and a.element.name not in ['H','D']}


def evaluate_sidechains(model_path,reference_path,cutoff=1.):
    backbone=evaluate(model_path,reference_path,cutoff)
    model=gemmi.read_structure(str(model_path));reference=gemmi.read_structure(str(reference_path))
    reference.setup_cell_images();cell=reference.cell
    models={(c.name,str(r.seqid)):r for c in model[0] for r in c}
    targets={(c.name,str(r.seqid)):r for c in reference[0] for r in c}
    mapping={(m['built_chain'],m['built_residue']):m for m in backbone['matches']}
    built=0;correct=0;distances=[];records=[]
    for key,residue in models.items():
        atoms=side_atoms(residue)
        built+=len(atoms)
        match=mapping.get(key)
        if match is None or not atoms:
            if atoms:
                records.append({'chain':key[0],'residue':key[1],'built_atoms':len(atoms),'within_cutoff':0,
                                'reason':'no CA match'})
            continue
        target=targets[(match['reference_chain'],match['reference_residue'])]
        if residue.name!=target.name:
            records.append({'chain':key[0],'residue':key[1],'built_atoms':len(atoms),'within_cutoff':0,
                            'reason':'residue type differs'})
            continue
        target_atoms=side_atoms(target)
        image=cell.find_nearest_image(residue.find_atom('CA','*').pos,target.find_atom('CA','*').pos)
        options=[{}]
        if residue.name in SYMMETRIES:
            options.append(SYMMETRIES[residue.name])
        outcomes=[]
        for permute in options:
            errors={}
            for name,atom in atoms.items():
                target_atom=target_atoms.get(permute.get(name,name))
                if target_atom is not None:
                    point=cell.orthogonalize(cell.fract_image(image,cell.fractionalize(target_atom.pos)))
                    errors[name]=atom.pos.dist(point)
            outcomes.append((sum(d<cutoff for d in errors.values()),sum(d*d for d in errors.values()),errors,permute))
        outcomes.sort(key=lambda x:(-x[0],x[1]))
        count,_,errors,permutation=outcomes[0]
        correct+=count;distances+=list(errors.values())
        records.append({'chain':key[0],'residue':key[1],'type':residue.name,'built_atoms':len(atoms),
                        'reference_atoms':len(target_atoms),'within_cutoff':count,'atom_distances':errors,
                        'equivalent_name_permutation':permutation,
                        'all_expected_atoms_present_and_within_cutoff':len(atoms)==len(target_atoms) and count==len(target_atoms)})
    expected_all=sum(len(side_atoms(r)) for c in reference[0] for r in c
                     if r.find_atom('CA','*') is not None)
    return {'model':fingerprint(model_path),'reference':fingerprint(reference_path),'cutoff_angstrom':cutoff,
            'sidechain_atoms_built':built,'sidechain_atoms_within_cutoff':correct,
            'sidechain_atom_precision':correct/built if built else None,
            'reference_sidechain_atoms_whole_protein':expected_all,
            'whole_protein_sidechain_atom_recovery':correct/expected_all if expected_all else None,
            'atom_rmsd_on_name_matched_atoms':float(np.sqrt(np.mean(np.square(distances)))) if distances else None,
            'fully_correct_complete_sidechains':sum(r.get('all_expected_atoms_present_and_within_cutoff',False) for r in records),
            'residues':records,'backbone_summary':{k:v for k,v in backbone.items() if k not in ['matches']},
            'scope':'One CA-selected crystal image per residue; chemically equivalent heavy-atom names allowed. '
                    'All built sidechain atoms are in the precision denominator, including spatially unmatched atoms. '
                    'RMSD uses all atom-name matches, not only those below the cutoff. Hydrogens excluded.'}
