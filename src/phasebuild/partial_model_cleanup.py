"""Explicit one-copy, sequence-anchored reduction of partial protein models.

Coordinates are only selected, never fitted or generated. Sequence and an
independent work-only DM map decide overlapping fragments; Rfree is not input.
"""
import math
import gemmi
import numpy as np


def _letter(residue):
    info=gemmi.find_tabulated_residue(residue.name)
    return info.one_letter_code.upper() if info.is_amino_acid() else ''


def work_density(map_file,native_file,free_value=0):
    from .sad_bridge import keys,crystal,columns
    m=gemmi.read_mtz_file(str(map_file));native=gemmi.read_mtz_file(str(native_file));crystal(m,native)
    columns(m,{'FWT':'F','PHWT':'P'});columns(native,{'FreeR_flag':'I'})
    lookup=dict(zip(keys(native),native.column_with_label('FreeR_flag').array));a=np.array(m,copy=True);labels=m.column_labels();retained=0;free_removed=0;keep=[]
    for i,h in enumerate(keys(m)):
        flag=lookup.get(h,float('nan'))
        if not np.isfinite(flag) or flag==free_value:
            free_removed+=bool(flag==free_value);a[i,[labels.index('FWT'),labels.index('PHWT')]]=np.nan
        elif np.isfinite(a[i,[labels.index('FWT'),labels.index('PHWT')]]).all():retained+=1;keep.append(i)
    if not retained:raise ValueError('No measured working map coefficients')
    m.set_data(a[keep]);grid=m.transform_f_phi_to_map('FWT','PHWT',sample_rate=3)
    if not np.isfinite(np.array(grid)).all() or np.std(np.array(grid))<=1e-12:raise ValueError('Working map has no finite density variation')
    grid.normalize()
    return grid,{'measured_work_coefficients':retained,'free_coefficients_removed':free_removed,'unknown_HKL_excluded':True,'normalization':'Map mean0 variance1; trilinear backbone sampling'}


def _connected(left,right):
    c=left.find_atom('C','*');n=right.find_atom('N','*');a=left.find_atom('CA','*');b=right.find_atom('CA','*')
    return bool(c and n and a and b and 1.1<=c.pos.dist(n.pos)<=1.7 and 2.7<=a.pos.dist(b.pos)<=4.5)


def cleanup_partial_single_copy(structure,sequence,grid,*,asu_copies,min_anchor=12,min_coverage=.8,min_long_anchor=40):
    if isinstance(asu_copies,bool) or not isinstance(asu_copies,int) or asu_copies!=1:raise ValueError('Explicit single-copy ASU required')
    sequence=''.join(sequence.split()).upper()
    if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('Unambiguous standard sequence required')
    if min_long_anchor is None:min_long_anchor=min(40,max(min_anchor,math.ceil(len(sequence)/2)))
    if isinstance(min_long_anchor,bool) or not isinstance(min_long_anchor,int) or min_long_anchor<min_anchor:raise ValueError('Long anchor must be an integer at least min_anchor')
    if len(structure)!=1 or grid is None:raise ValueError('One model and a work-only density map are required')
    kmers={}
    for j in range(len(sequence)-min_anchor+1):kmers.setdefault(sequence[j:j+min_anchor],[]).append(j)
    candidates={};anchors=[];residue_lookup={};protein_count=0
    for ci,chain in enumerate(structure[0]):
        runs=[];current=[]
        for ri,residue in enumerate(chain):
            if _letter(residue):protein_count+=1
            if _letter(residue) in 'ACDEFGHIKLMNPQRSTVWY' and _letter(residue) and residue.find_atom('CA','*'):
                current.append((ri,residue))
            else:
                if current:runs.append(current);current=[]
        if current:runs.append(current)
        found=set()
        for run in runs:
            letters=''.join(_letter(r) for _,r in run)
            for k in range(len(letters)-min_anchor+1):
                matches=kmers.get(letters[k:k+min_anchor],[])
                if len(matches)!=1:continue
                j=matches[0];left=0;right=min_anchor
                while k-left>0 and j-left>0 and letters[k-left-1]==sequence[j-left-1]:left+=1
                while k+right<len(letters) and j+right<len(sequence) and letters[k+right]==sequence[j+right]:right+=1
                token=(k-left,j-left,left+right);identity=(run[0][0],*token)
                if identity in found:continue
                found.add(identity);start,target,length=token
                text=letters[start:start+length]
                if sequence.find(text)!=sequence.rfind(text):continue
                anchors.append({'source_chain':chain.name,'source_chain_index':ci,'sequence_start':target+1,'length':length})
                for offset in range(length):
                    ri,res=run[start+offset];key=(ci,ri);residue_lookup[key]=res
                    values=[grid.interpolate_value(res.find_atom(n,'*').pos) for n in ['N','CA','C','O'] if res.find_atom(n,'*')]
                    if len(values)<3 or not np.isfinite(values).all():continue
                    score=float(np.mean(values))+.02*math.log1p(length)
                    entry={'source':key,'position':target+offset,'score':score,'anchor_length':length}
                    per=candidates.setdefault(target+offset,{})
                    if key not in per or per[key]['score']<score:per[key]=entry
    if not anchors or max(a['length'] for a in anchors)<min_long_anchor:raise ValueError('No sufficiently long unique sequence anchor')
    coverage=len(candidates)/len(sequence)
    if coverage<min_coverage:raise ValueError('Insufficient uniquely anchored sequence coverage')
    positions=sorted(candidates);layers=[];previous=[];previous_scores=[]
    for index,pos in enumerate(positions):
        options=[candidates[pos][k] for k in sorted(candidates[pos])];scores=[];backs=[]
        for option in options:
            if not previous:scores.append(option['score']);backs.append(-1);continue
            values=[]
            for j,old in enumerate(previous):
                penalty=0.
                if positions[index-1]+1==pos:
                    l=residue_lookup[old['source']];r=residue_lookup[option['source']]
                    same=old['source'][0]==option['source'][0] and old['source'][1]+1==option['source'][1]
                    if not (same and _connected(l,r)):
                        c=l.find_atom('C','*');n=r.find_atom('N','*');distance=c.pos.dist(n.pos) if c and n else 10.
                        penalty=1.+.25*min(abs(distance-1.33),4.)
                values.append(previous_scores[j]+option['score']-penalty)
            best=int(np.argmax(values));scores.append(values[best]);backs.append(best)
        layers.append((options,backs));previous=options;previous_scores=scores
    chosen=[];last=int(np.argmax(previous_scores))
    for options,backs in reversed(layers):chosen.append(options[last]);last=backs[last]
    chosen.reverse();result=structure.clone();result[0]=gemmi.Model(1)
    chain=None;previous=None;segments=[]
    for item in chosen:
        residue=residue_lookup[item['source']]
        joined=bool(previous and previous['position']+1==item['position'] and _connected(residue_lookup[previous['source']],residue))
        if not joined:
            name=chr(65+len(segments)) if len(segments)<26 else 'F'+str(len(segments)+1);chain=gemmi.Chain(name);result[0].add_chain(chain);chain=result[0][-1];segments.append({'chain':name,'sequence_start':item['position']+1,'residues':0})
        copied=residue.clone();copied.seqid=gemmi.SeqId(item['position']+1,' ');chain.add_residue(copied);segments[-1]['residues']+=1;previous=item
    # Nonprotein ligands/ions are retained without using any reference model.
    for original in structure[0]:
        extra=gemmi.Chain('L'+original.name)
        for residue in original:
            if not _letter(residue):extra.add_residue(residue.clone())
        if len(extra):result[0].add_chain(extra)
    report={'operation':'explicit-single-copy-sequence-anchored-partial-model-reduction','asu_copies':1,'sequence_length':len(sequence),'anchored_positions':len(chosen),'coverage':coverage,'input_protein_residues':protein_count,'output_protein_residues':len(chosen),'removed_protein_residues':protein_count-len(chosen),'overlapping_sequence_positions':sum(len(c)>1 for c in candidates.values()),'anchors':anchors,'segments':segments,'selection':[{'sequence_position':c['position']+1,'source_chain_index':c['source'][0],'source_residue_index':c['source'][1],'work_density_score':c['score']} for c in chosen],'coordinate_generation':False,'forced_nonphysical_joins':False,'selection_uses_Rfree':False,'thresholds':{'min_anchor':min_anchor,'min_coverage':min_coverage,'min_long_anchor':min_long_anchor},'limitation':'Exact sequence and work-map support do not establish coordinate correctness. Unanchored fragments excluded only under declared single-copy coverage gates.'}
    return result,report
