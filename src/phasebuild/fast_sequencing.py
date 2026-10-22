"""Opt-in scalable sequencing with an attainable familywise permutation threshold.

Candidate hypotheses, copy capacity, evidence and acceptance criteria are shared
with legacy sequencing; null scoring is batched, not statistically relaxed.
"""
import numpy as np
from .sequencing import windows,placement_scores
from .selection import independent_set
from .sequence_null import permutation_maxima,required_shuffles

def assign_sequence(chains,density,target,records,minimum_length=6,shuffles=499,alpha=.05,seed=0):
    if shuffles<19:
        raise ValueError('At least 19 composition shuffles are required')
    requested_shuffles=shuffles
    random = np.random.RandomState(seed)
    # Identical full sequence records are declared ASU copies. They must not
    # compete as distinct sequence hypotheses, but each has its own capacity.
    copies = {}
    representatives = []
    for ri, record in enumerate(records):
        if record[1] not in copies:
            representatives.append(ri)
        copies.setdefault(record[1], []).append(ri)
    unique_records = [records[ri] for ri in representatives]
    eligible = sum(minimum_length<=len(c)<=max(len(s) for _,s in records) for c in chains)
    threshold = alpha/max(1,eligible)
    shuffles=required_shuffles(shuffles,eligible,alpha)
    decisions,proposals = [],[]
    for ci,chain in enumerate(chains):
        record = {'fragment_index':ci,'length':len(chain)}
        choices = windows(unique_records,len(chain))
        if len(chain)<minimum_length or not choices:
            record.update(accepted=False,reason='too short or no sequence window')
            decisions.append(record)
            continue
        scores = target.scores(density,chain)
        values = placement_scores(scores,choices)
        order = np.argsort(values)[::-1]
        best = int(order[0]); best_value = float(values[best])
        margin = best_value-float(values[order[1]]) if len(order)>1 else best_value
        null = permutation_maxima(scores,[seq for _,seq in unique_records],shuffles,random)
        tail = (1+sum(v>=best_value for v in null))/(shuffles+1)
        unique_ri,start,_ = choices[best]
        ri=representatives[unique_ri]
        record.update(sequence_record=ri,record_id=records[ri][0],start=start+1,end=start+len(chain),
                      best_score=best_value,runner_up_margin=margin,shuffle_tail_fraction=tail,
                      corrected_tail_threshold=threshold,null_95=float(np.quantile(null,.95)))
        good = tail<=threshold and margin>=1.
        record.update(accepted=False,reason='candidate' if good else 'insufficient sequence evidence')
        decisions.append(record)
        if good:
            for copy_ri in copies[records[ri][1]]:
                proposals.append((ci,copy_ri,start,start+len(chain),max(.01,best_value-float(np.quantile(null,.95)))))
    conflicts = [0]*len(proposals)
    for i,a in enumerate(proposals):
        for j,b in enumerate(proposals[:i]):
            if a[0]==b[0] or (a[1]==b[1] and max(a[2],b[2])<min(a[3],b[3])):
                conflicts[i] |= 1<<j
                conflicts[j] |= 1<<i
    selected = independent_set([p[4] for p in proposals],[p[3]-p[2] for p in proposals],conflicts,
                                sum(len(seq) for _,seq in records))
    for i in selected:
        ci,ri,start,end,_ = proposals[i]
        decisions[ci].update(accepted=True,sequence_record=ri,record_id=records[ri][0],
                             start=start+1,end=end,reason='sequence and nonoverlap checks passed')
    for decision in decisions:
        if decision['reason']=='candidate':
            decision['reason']='sequence interval conflicts with a better-supported placement'
    assigned_positions=sum(d['length'] for d in decisions if d['accepted'])
    assigned_known=sum(sum(c!='X' for c in records[d['sequence_record']][1][d['start']-1:d['end']])
                       for d in decisions if d['accepted'])
    return {'decisions':decisions,'assigned_residues':assigned_known,'assigned_catalog_positions':assigned_positions,
            'identical_sequence_copy_groups':list(copies.values()),
            'requested_shuffles':requested_shuffles,'minimum_attainable_tail_fraction':1/(shuffles+1),
            'minimum_length':minimum_length,'shuffles':shuffles,'seed':seed,
            'familywise_tail_threshold':threshold,'alpha':alpha,
            'sequence_records':[{'id':name,'length':len(seq)} for name,seq in records],
            'assumptions':'One copy per sequence record; no sequence interval reused. A record is an indexing '
                          'catalog, not evidence that every adjacent entry forms a peptide bond.',
            'scope':'Composition-shuffle tail fractions are heuristic evidence, not calibrated correctness probabilities.'}

