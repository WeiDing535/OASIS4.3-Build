"""Local CCP4 monomer restraints and graph-derived single-bond rotations."""
from pathlib import Path
from collections import deque
import gemmi

BACKBONE={'N','CA','C','O','OXT'}


def monomer(name,root):
    path=Path(root)/name[0].lower()/(name+'.cif')
    document=gemmi.cif.read_file(str(path))
    block=document.find_block('comp_'+name)
    if block is None:
        raise ValueError('Missing monomer block: '+name)
    atoms=block.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol'])
    elements={r[0]:r[1].upper() for r in atoms if r[1].upper() not in ['H','D']}
    bonds=[]
    order='_chem_comp_bond.type' if block.find_values('_chem_comp_bond.type') else '_chem_comp_bond.value_order'
    for row in block.find(['_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2',
                           order,'_chem_comp_bond.value_dist']):
        if row[0] in elements and row[1] in elements:
            bonds.append([row[0],row[1],row[2].upper(),float(row[3])])
    if len(elements)>1 and not bonds:
        raise ValueError('No heavy-atom bonds found in monomer '+name)
    angles=[]
    for row in block.find(['_chem_comp_angle.atom_id_1','_chem_comp_angle.atom_id_2',
                           '_chem_comp_angle.atom_id_3','_chem_comp_angle.value_angle']):
        if all(row[i] in elements for i in range(3)):
            angles.append([row[0],row[1],row[2],float(row[3])])
    return {'name':name,'elements':elements,'sidechain_names':[a for a in elements if a not in BACKBONE],
            'bonds':bonds,'angles':angles,'path':str(path.resolve())}


def graph(component):
    result={a:set() for a in component['elements']}
    for a,b,_,_ in component['bonds']:
        result[a].add(b);result[b].add(a)
    return result


def component_after_cut(adjacency,start,a,b):
    seen={start};queue=[start]
    while queue:
        node=queue.pop()
        for neighbor in adjacency[node]:
            if {node,neighbor}=={a,b} or neighbor in seen:
                continue
            seen.add(neighbor);queue.append(neighbor)
    return seen


def rotatable_bonds(component):
    adjacency=graph(component);names=set(component['sidechain_names'])
    depth={'CA':0};queue=deque(['CA'])
    while queue:
        node=queue.popleft()
        for neighbor in sorted(adjacency[node]):
            if neighbor not in depth:
                depth[neighbor]=depth[node]+1;queue.append(neighbor)
    result=[]
    for a,b,kind,_ in component['bonds']:
        if kind!='SINGLE' or not ({a,b}&names):
            continue
        if depth.get(a,999)>depth.get(b,999):
            a,b=b,a
        moving=component_after_cut(adjacency,b,a,b)
        if a in moving or moving & (BACKBONE-{'OXT'}):
            continue # Ring or backbone bond, including the proline closure.
        moving=sorted((moving & names)-{a,b})
        if moving:
            result.append((a,b,moving))
    result.sort(key=lambda r:depth.get(r[0],999))
    return result


def close_pairs(component,names):
    adjacency=graph(component);result=set()
    for a in names:
        seen={a};frontier={a}
        for _ in range(2):
            frontier={v for u in frontier for v in adjacency[u]}-seen
            seen|=frontier
        for b in seen:
            if b in names:
                result.add((a,b))
    return result


def ring_groups(component):
    adjacency=graph(component);edges=[]
    for a,b,_,_ in component['bonds']:
        if a in component_after_cut(adjacency,b,a,b):
            edges.append((a,b))
    groups=[]
    for a,b in edges:
        intersect=[g for g in groups if a in g or b in g]
        merged={a,b}
        for g in intersect:
            merged|=g;groups.remove(g)
        groups.append(merged)
    return [sorted(g&set(component['sidechain_names'])) for g in groups if g&set(component['sidechain_names'])]
