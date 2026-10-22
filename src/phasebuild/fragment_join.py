"""Direction-aware overlap assembly of map-derived peptide fragments.

This experimental module accepts coordinates and a map, never a reference or
residue identity. It requires at least two N/CA/C overlap residues. Appended
residues are rebuilt from their torsions, retaining ideal trans-peptide bonds
instead of splicing independently fitted coordinates across a broken bond.
"""
import itertools
import gemmi
import numpy as np
from scipy.spatial import cKDTree
from .geometry import extend, oxygen
from .selection import independent_set


def _dihedral(points):
    return float(np.rad2deg(gemmi.calculate_dihedral(*[
        gemmi.Position(*map(float, p)) for p in points])))


def append_torsions(left, right, overlap):
    """Continue left with right's internal torsions, without a rigid-fit seam."""
    if overlap < 2 or overlap > min(len(left), len(right)):
        raise ValueError('An overlap of at least two residues is required')
    merged = [r.copy() for r in left]
    for i in range(overlap, len(right)):
        previous, current = right[i-1], right[i]
        psi = _dihedral([previous[0], previous[1], previous[2], current[0]])
        phi = _dihedral([previous[2], current[0], current[1], current[2]])
        merged[-1][3] = oxygen(*merged[-1][:3], psi)
        added = extend(merged[-1], psi, phi)
        if i+1 < len(right):
            next_psi = _dihedral([current[0], current[1], current[2], right[i+1,0]])
        else:
            next_psi = _dihedral([current[0], current[1], current[2], current[3]])-180.
        added[3] = oxygen(*added[:3], next_psi)
        merged.append(added)
    return np.asarray(merged)


class ImageIndex:
    """Spatial index of first residues, with symmetry and lattice images."""
    def __init__(self, fragments, density):
        self.fragments = fragments
        self.orth = np.asarray(density.cell.orth.mat.tolist())
        self.frac = np.asarray(density.cell.frac.mat.tolist())
        self.operations = [np.asarray(op.float_seitz()) for op in density.sg.operations()]
        self.frac_fragments = [c @ self.frac.T for c in fragments]
        points, labels = [], []
        for ci, xyz in enumerate(self.frac_fragments):
            for oi, operation in enumerate(self.operations):
                anchor = xyz[0,1] @ operation[:3,:3].T + operation[:3,3]
                lattice = -np.floor(anchor)
                for extra in itertools.product((-1,0,1), repeat=3):
                    shift = lattice + extra
                    points.append((anchor+shift) @ self.orth.T)
                    labels.append((ci, oi, shift))
        self.labels = labels
        self.tree = cKDTree(np.asarray(points))

    def query(self, anchor, radius):
        fractional = np.asarray(anchor) @ self.frac.T
        cell_shift = np.floor(fractional)
        wrapped = (fractional-cell_shift) @ self.orth.T
        for hit in self.tree.query_ball_point(wrapped, radius):
            ci, oi, shift = self.labels[hit]
            operation = self.operations[oi]
            xyz = self.frac_fragments[ci] @ operation[:3,:3].T + operation[:3,3]
            yield ci, (xyz+shift+cell_shift) @ self.orth.T


def overlap_graph(fragments, density, minimum_overlap=2, maximum_overlap=40,
                  rms_limit=.8, atom_limit=1.5, image_index=None):
    """Only matching atom names and residue order can establish an overlap."""
    if minimum_overlap < 2 or maximum_overlap < minimum_overlap:
        raise ValueError('Overlap bounds must satisfy 2 <= minimum <= maximum')
    index = image_index if image_index is not None else ImageIndex(fragments, density)
    edges = [[] for _ in fragments]
    for li, left in enumerate(fragments):
        matches = {}
        for overlap in range(minimum_overlap, min(maximum_overlap,len(left))+1):
            for ri, right in index.query(left[-overlap,1], atom_limit):
                if ri == li or len(right) <= overlap:
                    continue
                delta = np.linalg.norm(left[-overlap:,:3]-right[:overlap,:3], axis=-1)
                rms = float(np.sqrt(np.mean(delta**2)))
                if rms > rms_limit or float(delta.max()) > atom_limit:
                    continue
                key = (ri, overlap)
                if key not in matches or rms < matches[key][0]:
                    matches[key] = (rms, right)
        edges[li] = [(ri, overlap, value[0], value[1])
                     for (ri, overlap), value in matches.items()]
        edges[li].sort(key=lambda e: (e[2], -e[1], e[0]))
    return edges


def _self_clear(chain, previous_length):
    ca = chain[:,1]
    for i in range(previous_length, len(chain)):
        if i > 1 and np.min(np.linalg.norm(ca[:i-1]-ca[i],axis=-1)) < 2.8:
            return False
        if i > 1 and np.min(np.linalg.norm(chain[:i-1,:3].reshape(-1,3)[:,None]-
                                          chain[i,:3][None],axis=-1)) < 2.1:
            return False
    return True


def assemble_fragments(density, fragments, capacity, maximum_length=40,
                       beam_width=3, cutoff=2.5, progress=None):
    """Generate overlap paths, then select compatible complete peptide chains.

    A bounded search retains alternatives from every starting fragment. Map
    evidence counts each residue once, so duplicate seeds cannot inflate it.
    """
    if capacity < 3 or maximum_length < 3 or beam_width < 1:
        raise ValueError('Capacity/length must be >=3 and beam width positive')
    fragments = [np.asarray(c,dtype=float) for c in fragments if len(c)>=3]
    if any(c.ndim!=3 or c.shape[1:]!=(4,3) or not np.isfinite(c).all() for c in fragments):
        raise ValueError('Fragments must be finite arrays of N/CA/C/O coordinates')
    if not fragments:
        return [], {'input_fragments':0,'directed_overlap_edges':0,'selected_fragments':0}
    image_index=ImageIndex(fragments,density)
    edges = overlap_graph(fragments,density,image_index=image_index)
    candidates, accepted_joins, rejected = [], 0, {'displacement':0,'clash':0,'map':0}
    for start, initial in enumerate(fragments):
        initial_score = float(np.sum(density.score(initial)-cutoff))
        beam = [(initial_score, initial, start, frozenset([start]))]
        best = beam[0]
        for depth in range(maximum_length-len(initial)):
            branches=[]
            for total, chain, last, visited in beam:
                for following, overlap, rms, right in edges[last]:
                    if following in visited or len(chain)+len(right)-overlap > maximum_length:
                        continue
                    # The edge was measured in the original previous fragment's
                    # image. Put the new fragment near the assembled chain end.
                    options = [(ci,c) for ci,c in image_index.query(chain[-overlap,1],1.5)
                               if ci==following]
                    if not options:
                        continue
                    current_right = min((c for _,c in options), key=lambda c:
                                        np.mean((chain[-overlap:,:3]-c[:overlap,:3])**2))
                    mismatch=np.linalg.norm(chain[-overlap:,:3]-current_right[:overlap,:3],axis=-1)
                    if np.sqrt(np.mean(mismatch**2))>.8 or mismatch.max()>1.5:
                        continue
                    merged=append_torsions(chain,current_right,overlap)
                    new=merged[len(chain):]
                    movement=np.linalg.norm(new[:,:3]-current_right[overlap:,:3],axis=-1)
                    if movement.max()>1.0:
                        rejected['displacement']+=1;continue
                    if not _self_clear(merged,len(chain)):
                        rejected['clash']+=1;continue
                    scores=density.score(new)
                    old_scores=density.score(current_right[overlap:])
                    # Appending changes the preceding terminal O as well. Score
                    # the complete resulting chain so that atom cannot silently
                    # lose map evidence while only new residues are counted.
                    proposal_score=float(np.sum(density.score(merged)-cutoff))
                    gain=proposal_score-total
                    if scores.min()<1.5 or scores.mean()<old_scores.mean()-.25 or gain<=0:
                        rejected['map']+=1;continue
                    accepted_joins+=1
                    branches.append((proposal_score,merged,following,visited|{following}))
            if not branches:
                break
            branches.sort(key=lambda state:(-state[0],-len(state[1]),state[2]))
            beam=branches[:beam_width]
            if beam[0][0]>best[0]:
                best=beam[0]
        candidates.append((best[0],best[1]))
        if progress is not None and (start+1)%50==0:
            progress('Overlap paths %d/%d; best length %d'%(start+1,len(fragments),len(best[1])))
    candidates.sort(key=lambda item:item[0],reverse=True)
    # Keep longer paths and original high-scoring traces as competing candidates.
    candidates=candidates[:min(512,max(160,capacity))]
    conflicts=[0]*len(candidates)
    for i,(_,a) in enumerate(candidates):
        for j,(_,b) in enumerate(candidates[:i]):
            if any(density.distance(x,y)<2.5 for x in a[:,1] for y in b[:,1]):
                conflicts[i]|=1<<j;conflicts[j]|=1<<i
    selected=independent_set([max(0.,x[0]) for x in candidates],
                             [len(x[1]) for x in candidates],conflicts,capacity)
    chains=[candidates[i][1] for i in selected]
    return chains, {'input_fragments':len(fragments),
                    'directed_overlap_edges':sum(map(len,edges)),
                    'accepted_join_extensions':accepted_joins,'rejected_extensions':rejected,
                    'path_candidates':len(candidates),'selected_fragments':len(chains),
                    'fragment_lengths':list(map(len,chains)),
                    'map_objective':float(sum(candidates[i][0] for i in selected)),
                    'scope':'Map-only direction-consistent overlap assembly; no sequence or reference input. '
                            'Bounded path and conflict searches, with no optimality guarantee.'}
