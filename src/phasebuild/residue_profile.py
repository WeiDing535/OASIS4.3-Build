"""Independent residue-type density targets trained on unrelated protein maps."""
import json
import gemmi
import numpy as np
from .profile import frames, patches
from .data import fingerprint

ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
THREE = dict(zip(ALPHABET,['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                          'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']))


def training_groups(structure_path):
    groups, labels, cb = [], [], []
    for chain in gemmi.read_structure(str(structure_path))[0]:
        for residue in chain:
            letter = gemmi.find_tabulated_residue(residue.name).one_letter_code
            atoms = [residue.find_atom(n,'*') for n in ['N','CA','C','O']]
            if letter not in ALPHABET or not all(a is not None for a in atoms):
                continue
            xyz = np.array([list(a.pos) for a in atoms])
            groups.append(xyz)
            labels.append(ALPHABET.index(letter))
            atom = residue.find_atom('CB','*')
            if atom is not None:
                ca,axes = frames(xyz)
                cb.append((np.array(list(atom.pos))-ca)@axes)
    return np.array(groups),np.array(labels),np.array(cb)


def features(density, xyz, offsets):
    scale = np.clip(density.sample(xyz[...,:3,:]).mean(axis=-1),1.,6.)
    return np.clip(patches(density,xyz,offsets),-2.,8.) / scale[...,None]


class ResidueProfile:
    def __init__(self, offsets, means, variance, metadata):
        self.offsets,self.means,self.variance = map(np.asarray,[offsets,means,variance])
        if (self.offsets.ndim!=2 or self.offsets.shape[1]!=3 or
            self.means.shape!=(20,len(self.offsets)) or self.variance.shape!=(len(self.offsets),) or
            any(not np.isfinite(v).all() for v in [self.offsets,self.means,self.variance]) or
            np.any(self.variance<0)):
            raise ValueError('Invalid residue density profile')
        self.metadata = metadata

    def scores(self,density,xyz):
        values = features(density,xyz,self.offsets)
        scores = -.5*np.mean((values[...,None,:]-self.means)**2/(self.variance+.02),axis=-1)
        # Equal class prior. Center each site so sequence composition is not a prior.
        return (scores-scores.mean(axis=-1,keepdims=True))/np.maximum(scores.std(axis=-1,keepdims=True),.03)

    def save(self,path):
        np.savez_compressed(path,offsets=self.offsets,means=self.means,variance=self.variance,
                            metadata=json.dumps(self.metadata))

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:
            return cls(z['offsets'],z['means'],z['variance'],json.loads(str(z['metadata'])))


def train_residue_profile(density,structure_path,map_path,spacing=.8,pseudocount=5.):
    xyz,labels,cb = training_groups(structure_path)
    if len(cb)<30:
        raise ValueError('Too few sidechain anchors for training')
    cb_center = np.median(cb,axis=0)
    direction = cb_center/np.linalg.norm(cb_center)
    grid = np.arange(-6,7)*spacing
    offsets = np.array(np.meshgrid(grid,grid,grid,indexing='ij')).reshape(3,-1).T
    offsets = offsets[(np.linalg.norm(offsets-cb_center,axis=1)<=3.6)&(offsets@direction>.8)]
    # Remove the common backbone stem; retain sidechain-facing density.
    core = np.array([[0.,0.,0.],[1.525,0.,0.],[-.53,1.36,0.]])
    offsets = offsets[np.linalg.norm(offsets[:,None,:]-core[None,:,:],axis=-1).min(axis=1)>1.2]
    values = features(density,xyz,offsets)
    global_mean = values.mean(axis=0)
    counts = np.bincount(labels,minlength=20)
    if np.min(counts)<5:
        raise ValueError('At least 5 training residues of each amino-acid type are required')
    empirical = np.array([values[labels==i].mean(axis=0) for i in range(20)])
    means = np.array([(values[labels==i].sum(axis=0)+pseudocount*global_mean)/(counts[i]+pseudocount)
                      for i in range(20)])
    variance = np.mean((values-empirical[labels])**2,axis=0)
    metadata = {'alphabet':ALPHABET,'training_structure':fingerprint(structure_path),
                'training_map':fingerprint(map_path),'residues':len(xyz),
                'class_counts':{letter:int(counts[i]) for i,letter in enumerate(ALPHABET)},
                'resolution_high':density.mtz.resolution_high(),'spacing':spacing,
                'pseudocount':pseudocount,'local_cb_center':cb_center.tolist(),
                'method':'sidechain-facing local density, regularized class means, shared diagonal variance',
                'scope':'Heuristic residue-type scores, not posterior probabilities; no target protein used.'}
    return ResidueProfile(offsets,means,variance,metadata)
