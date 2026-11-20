"""Independent geometry/contact audit without target-reference coordinates."""
import gemmi
import numpy as np
from .chemistry import BACKBONE,close_pairs
from .sidechains import restraint_issues


def atom_geometry_report(structure,library):
    structure.setup_cell_images()
    model=structure[0];cell=structure.cell
    search=gemmi.NeighborSearch(model,cell,3.).populate(include_h=False)
    issues=[];clashes=[];sulfur=[];seen=set();atom_count=0
    for ci,chain in enumerate(model):
        for ri,residue in enumerate(chain):
            entry=library['entries'].get(residue.name)
            names=[a.name for a in residue if a.element.name not in ['H','D']]
            xyz=np.array([list(a.pos) for a in residue if a.element.name not in ['H','D']])
            if entry:
                excluded=close_pairs(entry['component'],names)
            else:
                # UNK still has the known N-CA-C-O backbone connectivity.
                pairs=[('N','CA'),('CA','C'),('C','O'),('N','C'),('CA','O')]
                excluded=set(pairs)|{(b,a) for a,b in pairs}
            if entry:
                for issue in restraint_issues(xyz,names,entry['component']):
                    issues.append(dict(chain=chain.name,residue=str(residue.seqid),**issue))
            for ai,atom in enumerate(residue):
                if atom.element.name in ['H','D']:
                    continue
                atom_count+=1
                owner=(ci,ri,ai)
                for mark in search.find_neighbors(atom,min_dist=0.,max_dist=2.3):
                    other=model[mark.chain_idx][mark.residue_idx][mark.atom_idx]
                    other_residue=model[mark.chain_idx][mark.residue_idx]
                    partner=(mark.chain_idx,mark.residue_idx,mark.atom_idx)
                    point=cell.find_nearest_pbc_position(atom.pos,mark.pos,0)
                    distance=atom.pos.dist(point)
                    identity=mark.image_idx==0 and point.dist(other.pos)<1e-5
                    if identity and owner==partner:
                        continue
                    if identity and owner[:2]==partner[:2] and (atom.name,other.name) in excluded:
                        continue
                    if identity and ci==mark.chain_idx and abs(ri-mark.residue_idx)==1:
                        left,right=(atom.name,other.name) if ri<mark.residue_idx else (other.name,atom.name)
                        if (left,right) in [('C','N'),('CA','N'),('O','N'),('C','CA')]:
                            continue
                    # Opposite views of the same pair/image distance are counted once.
                    key=(tuple(sorted((owner,partner))),round(distance,5))
                    if key in seen:
                        continue
                    seen.add(key)
                    contact={'atom_1':[chain.name,str(residue.seqid),atom.name],
                             'atom_2':[model[mark.chain_idx].name,str(other_residue.seqid),other.name],
                             'distance':distance,'symmetry_image':int(mark.image_idx),
                             'involves_sidechain':atom.name not in BACKBONE or other.name not in BACKBONE}
                    if distance<1.8:
                        clashes.append(contact)
                    if residue.name==other_residue.name=='CYS' and atom.name==other.name=='SG' and 1.8<=distance<=2.3:
                        sulfur.append(contact)
    return {'heavy_atoms':atom_count,'sidechain_restraint_issues':issues,
            'severe_nonbonded_contacts_below_1p8A':clashes,
            'severe_contacts_involving_sidechains':sum(c['involves_sidechain'] for c in clashes),
            'candidate_disulfide_contacts':sulfur,
            'scope':'Symmetry/PBC-aware heavy-atom contact audit; covalent 1-2/1-3 pairs excluded. '
                    'Sidechain bond and angle tolerances are 0.15 A and 20 degrees. '
                    'SG-SG proximity does not establish disulfide geometry or chemistry. '
                    'This is not complete stereochemical validation or crystallographic refinement.'}
