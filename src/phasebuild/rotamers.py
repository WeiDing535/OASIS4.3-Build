"""Empirical conformers from an unrelated structure; no target coordinates."""
import json
from pathlib import Path
import gemmi
import numpy as np
from .chemistry import monomer
from .data import fingerprint
from .profile import frames
from .residue_profile import THREE


def train_rotamers(structure_path,monomer_root,maximum=32):
    components={name:monomer(name,monomer_root) for name in THREE.values()}
    groups={name:[] for name in components}
    for chain in gemmi.read_structure(str(structure_path))[0]:
        for residue in chain:
            name=residue.name
            if name not in components:
                continue
            names=['N','CA','C','O']+components[name]['sidechain_names']
            atoms=[residue.find_atom(n,'*') for n in names]
            if any(a is None or a.occ<.5 for a in atoms):
                continue
            xyz=np.array([list(a.pos) for a in atoms])
            ca,axes=frames(xyz[:4])
            local=(xyz-ca)@axes
            groups[name].append((float(np.mean([a.b_iso for a in atoms])),local))
    entries={}
    for name,values in groups.items():
        values.sort(key=lambda r:r[0])
        chosen=[]
        for _,local in values:
            if not chosen or (len(local)>4 and min(np.sqrt(np.mean(np.sum((local[4:]-old[4:])**2,axis=-1)))
                                                       for old in chosen)>.6):
                chosen.append(local)
            if len(chosen)>=maximum:
                break
        if not chosen:
            raise ValueError('No complete independent conformer for '+name)
        entries[name]={'component':components[name],'conformers':[c.tolist() for c in chosen],
                       'training_residues':len(values)}
    return {'metadata':{'source_structure':fingerprint(structure_path),
                         'monomer_restraints':[fingerprint(c['path']) for c in components.values()],
                         'maximum_conformers_per_type':maximum,'clustering_rmsd_angstrom':.6,
                         'scope':'Empirical conformer library from unrelated 1TQW; no INS reference.'},
            'entries':entries}


def save_library(library,path):
    Path(path).write_text(json.dumps(library,indent=2,allow_nan=False)+'\n')


def load_library(path):
    library=json.loads(Path(path).read_text())
    for name,entry in library['entries'].items():
        array=np.asarray(entry['conformers'])
        if array.ndim!=3 or array.shape[1:]!=(4+len(entry['component']['sidechain_names']),3) or not np.isfinite(array).all():
            raise ValueError('Invalid rotamer library entry '+name)
    return library
