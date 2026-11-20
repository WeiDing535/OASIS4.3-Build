"""Batched composition-preserving sequence nulls with bounded working memory."""
import numpy as np
from .residue_profile import ALPHABET

def permutation_maxima(scores,sequences,trials,random,batch_size=32):
    scores=np.column_stack([scores,np.zeros(len(scores))])
    length=len(scores)
    codes=[np.array([ALPHABET.index(x) if x!='X' else 20 for x in seq],dtype=np.int32) for seq in sequences]
    indices=[np.arange(len(seq)-length+1)[:,None]+np.arange(length)[None,:] if len(seq)>=length else None for seq in codes]
    results=[]
    # Cap the indexed score tensor to about 2 million values per batch.
    windows=max((len(x) for x in indices if x is not None),default=1)
    block=max(1,min(batch_size,2_000_000//max(1,windows*length)))
    for start in range(0,trials,block):
        n=min(block,trials-start)
        permutations=[[random.permutation(code) for code in codes] for _ in range(n)]
        maxima=np.full(n,-np.inf)
        for ri,index in enumerate(indices):
            if index is None:continue
            catalog=np.stack([p[ri] for p in permutations])[:,index]
            totals=scores[np.arange(length)[None,None,:],catalog].sum(axis=-1)
            maxima=np.maximum(maxima,totals.max(axis=-1))
        results.extend(maxima.tolist())
    return np.array(results)


def required_shuffles(requested,eligible,alpha=.05):
    """Make the minimum attainable tail probability reach the stated threshold."""
    if not 0<alpha<=1 or requested<19:raise ValueError('Invalid permutation test configuration')
    return max(int(requested),int(np.ceil(max(1,eligible)/alpha))-1)
