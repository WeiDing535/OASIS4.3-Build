"""Map-only fitting of free terminal torsions; not crystallographic refinement."""
import numpy as np


def rotate(points,origin,axis,angles):
    axis=axis/np.linalg.norm(axis)
    vector=points-origin
    theta=np.deg2rad(angles).reshape(-1,1,1)
    return (origin+vector[None]*np.cos(theta)+np.cross(axis,vector)[None]*np.sin(theta)
            +axis*np.sum(vector*axis,axis=-1)[None,:,None]*(1.-np.cos(theta)))


def fit_terminals(chains,density,minimum_gain=.15):
    chains=[c.copy() for c in chains]
    changes=[]
    angles=np.arange(-180.,180.,5.)
    for ci,chain in enumerate(chains):
        for which in ['first_N','last_C_O']:
            ri=0 if which=='first_N' else len(chain)-1
            original=chain[ri].copy()
            old_score=float(density.score(original))
            if which=='first_N':
                variants=np.repeat(original[None],len(angles),axis=0)
                variants[:,0]=rotate(original[0:1],original[1],original[2]-original[1],angles)[:,0]
            else:
                phi=np.repeat(original[None],len(angles),axis=0)
                phi[:,2:]=rotate(original[2:],original[1],original[1]-original[0],angles)
                variants=[]
                for residue in phi:
                    candidates=np.repeat(residue[None],len(angles),axis=0)
                    candidates[:,3]=rotate(residue[3:4],residue[2],residue[2]-residue[1],angles)[:,0]
                    variants.append(candidates)
                variants=np.concatenate(variants)
            scores=density.score(variants)
            order=np.argsort(scores)[::-1]
            other_atoms=[point for cj,c in enumerate(chains) for rj,r in enumerate(c)
                         if (cj,rj)!=(ci,ri) for point in r]
            chosen=None
            for index in order:
                if scores[index]<old_score+minimum_gain:
                    break
                candidate=variants[index]
                changed=candidate[0:1] if which=='first_N' else candidate[2:]
                # Same residue is excluded; its ideal geometry is preserved analytically.
                if any(density.distance(point,other)<1.8 for point in changed for other in other_atoms):
                    continue
                chosen=int(index)
                break
            if chosen is not None:
                chain[ri]=variants[chosen]
                changes.append({'fragment_index':ci,'residue_index':ri,'kind':which,
                    'old_score':old_score,'new_score':float(scores[chosen]),
                    'maximum_atom_displacement':float(np.linalg.norm(chain[ri]-original,axis=-1).max())})
    return chains,{'method':'grid search of free terminal N and C/O torsions',
                   'minimum_gain':minimum_gain,'grid_degrees':5.,'changes':changes,
                   'scope':'CA positions and internal peptide links remain fixed. '
                           'Map score improvements do not guarantee coordinate correctness.'}
