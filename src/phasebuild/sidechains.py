"""Empirical sidechain placement, torsion fitting and conservative atom pruning."""
import gemmi
import numpy as np
from scipy.optimize import minimize
from .chemistry import BACKBONE, graph, close_pairs, rotatable_bonds, ring_groups
from .profile import frames
from .validation import angle


def neighbor_atoms(structure,owner,density,center,radius=12.):
    positions=[];labels=[]
    for ci,chain in enumerate(structure[0]):
        for ri,residue in enumerate(chain):
            same=(ci,ri)==owner
            for atom in residue:
                if atom.element.name in ['H','D']:
                    continue
                fractional=density.cell.fractionalize(atom.pos)
                for op in density.sg.operations():
                    transformed=density.cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(list(fractional))))
                    pos=density.cell.find_nearest_pbc_position(gemmi.Position(*center),transformed,0)
                    if same and pos.dist(atom.pos)<1e-5:
                        continue
                    xyz=np.array(list(pos))
                    if np.linalg.norm(xyz-center)<radius:
                        positions.append(xyz);labels.append((residue.name,atom.name))
    return np.array(positions).reshape(-1,3),labels


def rotate_chis(base,rotors,names,angles):
    xyz=base.copy();lookup={name:i for i,name in enumerate(names)}
    for (a,b,moving),theta in zip(rotors,angles):
        origin=xyz[lookup[b]]
        axis=xyz[lookup[b]]-xyz[lookup[a]];axis/=np.linalg.norm(axis)
        indices=[lookup[n] for n in moving]
        vector=xyz[indices]-origin
        xyz[indices]=(origin+vector*np.cos(theta)+np.cross(axis,vector)*np.sin(theta)
                      +np.outer(vector@axis,axis)*(1.-np.cos(theta)))
    return xyz


def restraint_issues(xyz,names,component):
    lookup={n:i for i,n in enumerate(names)};issues=[]
    for a,b,_,ideal in component['bonds']:
        if a not in lookup or b not in lookup or a in BACKBONE and b in BACKBONE:
            continue
        value=float(np.linalg.norm(xyz[lookup[a]]-xyz[lookup[b]]))
        if abs(value-ideal)>.15:
            issues.append({'kind':'bond','atoms':[a,b],'value':value,'ideal':ideal})
    for a,b,c,ideal in component['angles']:
        if not all(n in lookup for n in [a,b,c]) or all(n in BACKBONE for n in [a,b,c]):
            continue
        value=angle(xyz[lookup[a]],xyz[lookup[b]],xyz[lookup[c]])
        if abs(value-ideal)>20.:
            issues.append({'kind':'angle','atoms':[a,b,c],'value':value,'ideal':ideal})
    return issues


def collision_data(xyz,names,component,neighbors,labels):
    sc=xyz[4:];n=len(sc)
    internal=np.linalg.norm(sc[:,None,:]-xyz[None,:,:],axis=-1)
    excluded=close_pairs(component,names)
    for i,a in enumerate(names[4:]):
        for j,b in enumerate(names):
            if (a,b) in excluded:
                internal[i,j]=100.
    external=np.linalg.norm(sc[:,None,:]-neighbors[None,:,:],axis=-1)
    limits=np.full_like(external,2.3)
    if component['name']=='CYS' and 'SG' in names:
        i=names[4:].index('SG')
        for j,label in enumerate(labels):
            if label==('CYS','SG'):
                limits[i,j]=1.85 # Candidate disulfides are reported separately.
    penalty=np.maximum(2.3-internal,0.)**2
    external_penalty=np.maximum(limits-external,0.)**2
    per_atom=np.minimum(internal.min(axis=1),external.min(axis=1) if len(neighbors) else np.full(n,100.))
    return float((penalty.sum()+external_penalty.sum())/max(n,1)),per_atom


def atom_map_score(density,sc):
    raw=density.sample(sc)
    shell=sc[:,None,:]+1.2*np.concatenate([np.eye(3),-np.eye(3)])[None,:,:]
    contrast=np.clip(raw,-3.,6.)-.3*np.clip(density.sample(shell),-3.,6.).mean(axis=-1)
    return float(contrast.mean()),raw


def retained_atoms(raw,clearance,names,component,minimum_density=.7):
    good={n for n,r,d in zip(names[4:],raw,clearance) if r>=minimum_density and d>=1.8}
    for ring in ring_groups(component):
        if not set(ring)<=good:
            good-=set(ring)
    adjacency=graph(component);connected={'CA'};queue=['CA']
    while queue:
        node=queue.pop()
        for child in adjacency[node]:
            if child in good and child not in connected:
                connected.add(child);queue.append(child)
    return [i for i,name in enumerate(names[4:]) if name in connected]


def build_sidechains(structure,density,library,progress=print,maximum_seeds=8):
    records=[]
    for ci,chain in enumerate(structure[0]):
        for ri,residue in enumerate(chain):
            name=residue.name
            if name not in library['entries'] or name=='GLY':
                continue
            entry=library['entries'][name];component=entry['component']
            atoms=[residue.find_atom(n,'*') for n in ['N','CA','C','O']]
            if any(a is None for a in atoms):
                continue
            backbone=np.array([list(a.pos) for a in atoms])
            ca,axes=frames(backbone)
            names=['N','CA','C','O']+component['sidechain_names']
            neighbors,labels=neighbor_atoms(structure,(ci,ri),density,ca)
            backbone_b=float(np.mean([a.b_iso for a in atoms]))
            rotors=rotatable_bonds(component)
            def score(xyz):
                contrast,raw=atom_map_score(density,xyz[4:])
                penalty,clearance=collision_data(xyz,names,component,neighbors,labels)
                return contrast-4.*penalty,raw,clearance
            candidates=[]
            for local in entry['conformers']:
                xyz=np.array(local)@axes.T+ca
                xyz[:4]=backbone
                candidates.append((score(xyz)[0],xyz))
            candidates.sort(key=lambda r:r[0],reverse=True)
            seeds=[xyz for _,xyz in candidates[:maximum_seeds]]
            current=[residue.find_atom(n,'*') for n in names]
            has_current=all(a is not None for a in current)
            if has_current:
                seeds.append(np.array([list(a.pos) for a in current]))
            fits=[]
            for base in seeds:
                if rotors:
                    fit=minimize(lambda x:-score(rotate_chis(base,rotors,names,x))[0],
                                 np.zeros(len(rotors)),method='L-BFGS-B',
                                 bounds=[(-np.deg2rad(35.),np.deg2rad(35.))]*len(rotors),
                                 options={'maxiter':40,'eps':1e-3,'ftol':1e-6})
                    xyz=rotate_chis(base,rotors,names,fit.x)
                else:
                    xyz=base
                quality,raw,clearance=score(xyz)
                if restraint_issues(xyz,names,component):
                    continue
                keep=retained_atoms(raw,clearance,names,component)
                if not keep or np.mean(raw[keep])<1.2:
                    continue
                utility=float(np.maximum(raw[keep]-.7,0.).sum())
                fits.append((utility,quality,xyz,keep,raw))
            record={'chain':chain.name,'residue':str(residue.seqid),'type':name,
                    'expected_sidechain_atoms':len(names)-4,'library_conformers':len(candidates),
                    'current_complete_conformer_seed':has_current,'seeds_fitted':len(seeds)}
            if fits:
                fits.sort(key=lambda r:(r[0],r[1]),reverse=True)
                _,quality,xyz,keep,raw=fits[0]
                for index in range(len(residue)-1,-1,-1):
                    if residue[index].name not in BACKBONE:
                        del residue[index]
                for index in keep:
                    atom=gemmi.Atom();atom.name=names[4+index]
                    atom.element=gemmi.Element(component['elements'][atom.name])
                    atom.pos=gemmi.Position(*xyz[4+index]);atom.occ=1.;atom.b_iso=backbone_b
                    residue.add_atom(atom)
                record.update(built_atoms=[names[4+i] for i in keep],score=quality,
                              atom_density={names[4+i]:float(raw[i]) for i in keep},
                              complete=len(keep)==len(names)-4)
            else:
                existing=[a for a in residue if a.name not in BACKBONE and a.element.name not in ['H','D']]
                record.update(built_atoms=[a.name for a in existing],complete=len(existing)==len(names)-4,
                              reason='no supported new conformer; existing atoms retained')
            records.append(record)
        progress('Sidechain fitting finished fragment %s'%chain.name)
    return {'method':'independent empirical conformers, graph-derived chi refinement, connected atom pruning',
            'maximum_library_seeds_per_residue':maximum_seeds,
            'residues_attempted':len(records),'residues_with_atoms':sum(bool(r['built_atoms']) for r in records),
            'complete_sidechains':sum(r['complete'] for r in records),
            'sidechain_atoms_built':sum(len(r['built_atoms']) for r in records),'residues':records,
            'scope':'Main-chain coordinates fixed. Partial sidechains retained; no crystallographic refinement. '
                    'SG-SG contacts are allowed provisionally and are not imposed as known disulfide bonds.'}
