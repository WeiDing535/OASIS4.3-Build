"""Opt-in overlap assembly with unrelated-model evidence before conflicts.

The classifier is frozen elsewhere from unrelated training data. Its logits are
synthetic-task scores, not probabilities of experimental residue correctness.
Map support and peptide geometry still constrain every continuation.
"""
import numpy as np
from .backbone_discriminant import supported_segments
from .fragment_join import ImageIndex,overlap_graph,append_torsions,_self_clear
from .selection import independent_set


def positive_evidence(logits,threshold):
    """Nonnegative margin above a threshold fixed by unrelated calibration."""
    return float(np.maximum(np.asarray(logits)-threshold,0.).sum())


def assemble_discriminant_fragments(density,fragments,classifier,capacity,
                                    maximum_length=40,beam_width=3,
                                    calibration_decoy_rate='0.05',progress=None):
    """Prune/rank alternatives before selecting geometrically conflicting paths.

    Fixed recipe: held-out-training decoy95th-percentile threshold, at least3
    residues, at most1 internal weak residue, and no weak terminal tails.
    Tripeptides remain eligible because they establish directional overlaps.
    """
    if calibration_decoy_rate!='0.05':
        raise ValueError('This focused assembly protocol uses the fixed 0.05 calibration entry')
    threshold=float(classifier.metadata['calibration']['thresholds'][calibration_decoy_rate]['logit'])
    if capacity<3 or maximum_length<3 or beam_width<1 or not np.isfinite(threshold):
        raise ValueError('Invalid assembly capacity, length, width, or calibration threshold')
    raw=[np.asarray(c,dtype=float) for c in fragments if len(c)>=3]
    if any(c.ndim!=3 or c.shape[1:]!=(4,3) or not np.isfinite(c).all() for c in raw):
        raise ValueError('Fragments must contain finite N/CA/C/O coordinates')
    pieces=[]
    for chain in raw:
        logits=classifier.score(density,chain)
        for a,b in supported_segments(logits,threshold,minimum_length=3,bridge=1):
            pieces.append(chain[a:b])
    report={'input_fragments':len(raw),'input_residues':sum(map(len,raw)),
            'classifier_supported_pieces':len(pieces),
            'classifier_supported_residues':sum(map(len,pieces)),
            'calibration_decoy_rate':calibration_decoy_rate,'classifier_logit_threshold':threshold,
            'minimum_supported_length':3,'maximum_internal_weak_island':1,
            'objective':'Sum of max(classifier_logit - fixed calibration threshold, 0); each residue once',
            'warning':'Synthetic-task logits/decoy quantile do not imply experimental probabilities or false-positive rates'}
    if not pieces:
        return [],dict(report,directed_overlap_edges=0,selected_fragments=0,fragment_lengths=[])
    image_index=ImageIndex(pieces,density)
    edges=overlap_graph(pieces,density,image_index=image_index)
    candidates=[];accepted=0;rejected={'displacement':0,'clash':0,'map':0,'classifier_continuity':0}
    for start,initial in enumerate(pieces):
        logits=classifier.score(density,initial)
        # Classifier evidence ranks paths; the original map objective guards
        # every continuation independently, retaining the prior map protocol.
        beam=[(positive_evidence(logits,threshold),float(np.sum(density.score(initial)-2.5)),
               initial,start,frozenset([start]))]
        best=beam[0]
        for _ in range(maximum_length-len(initial)):
            branches=[]
            for evidence,map_total,chain,last,visited in beam:
                for following,overlap,rms,right in edges[last]:
                    if following in visited or len(chain)+len(right)-overlap>maximum_length:
                        continue
                    options=[c for ci,c in image_index.query(chain[-overlap,1],1.5) if ci==following]
                    if not options:
                        continue
                    right=min(options,key=lambda c:np.mean((chain[-overlap:,:3]-c[:overlap,:3])**2))
                    mismatch=np.linalg.norm(chain[-overlap:,:3]-right[:overlap,:3],axis=-1)
                    if np.sqrt(np.mean(mismatch**2))>.8 or mismatch.max()>1.5:
                        continue
                    merged=append_torsions(chain,right,overlap);new=merged[len(chain):]
                    if np.linalg.norm(new[:,:3]-right[overlap:,:3],axis=-1).max()>1.:
                        rejected['displacement']+=1;continue
                    if not _self_clear(merged,len(chain)):
                        rejected['clash']+=1;continue
                    scores=density.score(new);old_scores=density.score(right[overlap:])
                    proposal_map=float(np.sum(density.score(merged)-2.5))
                    if scores.min()<1.5 or scores.mean()<old_scores.mean()-.25 or proposal_map<=map_total:
                        rejected['map']+=1;continue
                    logits=classifier.score(density,merged)
                    supported=supported_segments(logits,threshold,minimum_length=3,bridge=1)
                    if supported!=[(0,len(merged))]:
                        rejected['classifier_continuity']+=1;continue
                    proposal_evidence=positive_evidence(logits,threshold)
                    if proposal_evidence<=evidence:
                        rejected['classifier_continuity']+=1;continue
                    accepted+=1
                    branches.append((proposal_evidence,proposal_map,merged,following,visited|{following}))
            if not branches:
                break
            branches.sort(key=lambda state:(-state[0],-state[1],-len(state[2]),state[3]))
            beam=branches[:beam_width]
            if beam[0][:2]>best[:2]:
                best=beam[0]
        candidates.append((best[0],best[2]))
        if progress is not None and (start+1)%50==0:
            progress('Classifier overlap paths %d/%d; best length %d'%(start+1,len(pieces),len(best[2])))
    candidates.sort(key=lambda item:item[0],reverse=True)
    candidates=candidates[:min(512,max(160,capacity))]
    conflicts=[0]*len(candidates)
    for i,(_,a) in enumerate(candidates):
        for j,(_,b) in enumerate(candidates[:i]):
            if any(density.distance(x,y)<2.5 for x in a[:,1] for y in b[:,1]):
                conflicts[i]|=1<<j;conflicts[j]|=1<<i
    chosen=independent_set([item[0] for item in candidates],
                           [len(item[1]) for item in candidates],conflicts,capacity)
    chains=[candidates[i][1] for i in chosen]
    report.update(directed_overlap_edges=sum(map(len,edges)),accepted_join_extensions=accepted,
                  rejected_extensions=rejected,path_candidates=len(candidates),
                  selected_fragments=len(chains),fragment_lengths=list(map(len,chains)),
                  classifier_objective=float(sum(candidates[i][0] for i in chosen)),
                  scope='Map, cached native candidates and frozen unrelated-model classifier only; no target labels. '
                        'Bounded path/conflict search with no optimality guarantee.')
    return chains,report
