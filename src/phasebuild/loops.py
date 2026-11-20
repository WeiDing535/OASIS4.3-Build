"""De novo short-loop torsion search with fixed backbone anchors.

No parent trace or target-reference coordinates are used. This bounded search
is conservative and can return no accepted loop.
"""
import gemmi
import numpy as np
from scipy.optimize import least_squares,minimize
from .geometry import extend,oxygen,CONFORMATIONS
from .validation import peptide_metrics,peptide_valid
from .residue_profile import ALPHABET


def sequence_bindings(structure,report,records):
    bindings=[None]*len(structure[0])
    for decision in report['decisions']:
        if not decision['accepted']:
            continue
        ci=decision['fragment_index'];chain=structure[0][ci]
        start=decision['start']-1;end=decision['end'];record=decision['sequence_record']
        seq=records[record][1]
        letters=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in chain)
        if len(chain)!=end-start or letters!=seq[start:end] or any(r.seqid.num!=start+i+1 for i,r in enumerate(chain)):
            raise ValueError('Sequence report is inconsistent with the supplied model')
        bindings[ci]={'record':record,'start':start,'end':end}
    return bindings


def sequence_hypothesis(left_length,right_length,left,right,gap,records):
    if left is None and right is None:
        return None
    if left is not None and right is not None:
        if left['record']!=right['record'] or right['start']-left['end']!=gap:
            return None
        ri=left['record'];start=left['start'];end=right['end']
    elif left is not None:
        ri=left['record'];start=left['start'];end=left['end']+gap+right_length
    else:
        ri=right['record'];end=right['end'];start=right['start']-gap-left_length
    if start<0 or end>len(records[ri][1]):
        return None
    return {'record':ri,'start':start,'end':end}


def transformed_images(fragment,reference,density,maximum_distance):
    output=[]
    for op in density.sg.operations():
        xyz=np.array([list(density.cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(
            list(density.cell.fractionalize(gemmi.Position(*p))))))) for p in fragment.reshape(-1,3)]).reshape(fragment.shape)
        nearest=density.cell.find_nearest_pbc_position(gemmi.Position(*reference),gemmi.Position(*xyz[0,1]),0)
        xyz+=np.array(list(nearest))-xyz[0,1]
        distance=float(np.linalg.norm(xyz[0,1]-reference))
        if 2.8<distance<=maximum_distance and not any(np.linalg.norm(xyz[0,1]-old[1][0,1])<.1 for old in output):
            output.append((distance,xyz,op.triplet()))
    output.sort(key=lambda item:item[0])
    return output[:2]


def kinematic_bridge(left_last,angles):
    residues=[left_last.copy()]
    for psi,phi in np.rad2deg(np.asarray(angles).reshape(-1,2)):
        residues[-1][3]=oxygen(*residues[-1][:3],psi)
        residues.append(extend(residues[-1],psi,phi))
    return np.array(residues)


def loop_candidates(left,right,gap,density,target,missing_sequence,anchor_sequence_gain=0.,starts=48,seed=0,
                    right_letter='X',flexible_right=0,flexible_sequence=''):
    random=np.random.RandomState(seed)
    flex_target=right[:flexible_right]
    right=right[flexible_right:]
    search_sequence=missing_sequence+flexible_sequence
    search_gap=gap+flexible_right
    sigma=np.full((flexible_right,3,1),.3)
    if flexible_right:
        sigma[0,0]=.7 # Terminal N lacks an upstream peptide restraint in a fragment.
    bounds=[]
    for i in range(search_gap+1):
        letter=search_sequence[i] if i<search_gap else right_letter
        maximum=180. if letter=='G' else (80. if letter in 'DNSTX?' else 0.)
        bounds.extend([(-np.pi,np.pi),(-np.pi,np.deg2rad(maximum))])
    lower=np.array([b[0] for b in bounds]);upper=np.array([b[1] for b in bounds])
    initial=[]
    for mode in [0,1,2]:
        phi,psi=CONFORMATIONS[mode]
        initial.append(np.tile(np.deg2rad([psi,phi]),search_gap+1))
    while len(initial)<starts:
        values=[]
        for _ in range(search_gap+1):
            phi,psi=CONFORMATIONS[random.choice(len(CONFORMATIONS),p=[.4,.25,.15,.15,.05])]
            values.extend(np.deg2rad([psi+random.normal(0.,20.),phi+random.normal(0.,20.)]))
        initial.append(np.array(values))
    feasible=[]
    def closure(x):
        bridge=kinematic_bridge(left[-1],x)
        result=(bridge[-1,:3]-right[0,:3]).reshape(-1)
        if flexible_right:
            prior=.1*(bridge[1+gap:-1,:3]-flex_target[:,:3])/sigma
            result=np.r_[result,prior.reshape(-1)]
        return result
    for x in initial[:starts]:
        x=np.clip(x,lower+1e-6,upper-1e-6)
        fitted=least_squares(closure,x,bounds=(lower,upper),max_nfev=55,ftol=.01,xtol=1e-4,gtol=1e-4)
        bridge=kinematic_bridge(left[-1],fitted.x)
        rms=float(np.sqrt(np.mean(np.sum((bridge[-1,:3]-right[0,:3])**2,axis=-1))))
        if rms>.4:
            continue
        missing=bridge[1:-1]
        if feasible and any(np.sqrt(np.mean(np.sum((missing[:,1]-old[2][:,1])**2,axis=-1)))<.5 for old in feasible):
            continue
        quality=float(density.score(missing).mean())
        feasible.append((quality,fitted.x,missing))
    feasible.sort(key=lambda item:item[0],reverse=True)
    output=[]
    rejected={key:0 for key in ['closure','missing_residue_density','anchor_displacement',
                              'anchor_density','left_oxygen_density','peptide_geometry']}
    codes=np.array([ALPHABET.index(c) if c in ALPHABET else -1 for c in search_sequence])
    def identity_score(missing):
        if target is None:
            return 0.
        values=target.scores(density,missing)
        return float(np.mean([values[i,c] if c>=0 else 0. for i,c in enumerate(codes)]))
    for _,parameters,_ in feasible[:8]:
        def loss(x):
            bridge=kinematic_bridge(left[-1],x);missing=bridge[1:-1]
            penalty=100.*np.mean(np.sum((bridge[-1,:3]-right[0,:3])**2,axis=-1))
            if flexible_right:
                penalty+=np.mean(((missing[gap:,:3]-flex_target[:,:3])/sigma)**2)
            return penalty-float(density.score(missing).mean())-.15*identity_score(missing)
        fitted=minimize(loss,parameters,method='L-BFGS-B',bounds=bounds,
                        options={'maxiter':45,'eps':1e-3,'ftol':1e-6})
        bridge=kinematic_bridge(left[-1],fitted.x);missing=bridge[1:-1]
        rms=float(np.sqrt(np.mean(np.sum((bridge[-1,:3]-right[0,:3])**2,axis=-1))))
        profile=density.score(missing)
        if rms>.25:
            rejected['closure']+=1
            continue
        if profile[:gap].min()<2.0 or profile[:gap].mean()<2.6 or density.sample(missing[:gap,:3]).min()<.5:
            rejected['missing_residue_density']+=1
            continue
        displacement=None
        if flexible_right:
            delta=np.linalg.norm(missing[gap:,:3]-flex_target[:,:3],axis=-1)
            allowed=np.full_like(delta,.55);allowed[0,0]=1.4
            if np.any(delta>allowed):
                rejected['anchor_displacement']+=1
                continue
            if np.mean(profile[gap:])<np.mean(density.score(flex_target))-.2:
                rejected['anchor_density']+=1
                continue
            displacement=float(delta.max())
        modified_left=left.copy();modified_left[-1,3]=bridge[0,3]
        if density.sample(modified_left[-1,3])<.5:
            rejected['left_oxygen_density']+=1
            continue
        merged=np.concatenate([modified_left,missing,right])
        links=[peptide_metrics(merged[k],merged[k+1]) for k in range(len(left)-1,len(left)+search_gap)]
        if not all(peptide_valid(m) for m in links):
            rejected['peptide_geometry']+=1
            continue
        quality=float(profile.mean())+.15*identity_score(missing)+.15*anchor_sequence_gain
        output.append({'xyz':merged,'missing':missing[:gap],'quality':quality,'gap':gap,
                       'flexible_right_residues':flexible_right,'maximum_anchor_atom_displacement':displacement,
                       'closure_rmsd':rms,'profile_scores':profile.tolist(),'identity_score':identity_score(missing),
                       'anchor_sequence_gain':anchor_sequence_gain,'links':links})
    return output,{'closure_starts':len(initial[:starts]),'distinct_feasible_closures':len(feasible),
                   'density_and_geometry_passed':len(output),'flexible_right_residues':flexible_right,
                   'refined_candidates':min(len(feasible),8),'first_rejection_reasons':rejected}


def build_short_loops(chains,density,records,bindings,target,starts=48,seed=0,progress=print):
    if starts<3:
        raise ValueError('At least three closure starts required')
    chains=[c.copy() for c in chains];bindings=list(bindings)
    accepted=[];attempts=[]
    while True:
        proposals=[]
        for i,left in enumerate(chains):
            if len(left)<4 or (bindings[i] is None and np.mean(density.score(left))<3.):
                continue
            for j,right_original in enumerate(chains):
                if i==j or len(right_original)<4 or (bindings[j] is None and np.mean(density.score(right_original))<3.):
                    continue
                for gap in [1,2,3]:
                    binding=sequence_hypothesis(len(left),len(right_original),bindings[i],bindings[j],gap,records)
                    if binding is None:
                        continue
                    if any(k not in [i,j] and b is not None and b['record']==binding['record'] and
                           max(b['start'],binding['start'])<min(b['end'],binding['end']) for k,b in enumerate(bindings)):
                        continue
                    seq=records[binding['record']][1]
                    missing_seq=seq[binding['start']+len(left):binding['start']+len(left)+gap]
                    anchor_gains=[]
                    for fragment,b,offset in [(left,bindings[i],binding['start']),
                                               (right_original,bindings[j],binding['end']-len(right_original))]:
                        if b is None:
                            values=target.scores(density,fragment)
                            expected=seq[offset:offset+len(fragment)]
                            gain=np.mean([values[k,ALPHABET.index(c)] if c in ALPHABET else 0. for k,c in enumerate(expected)])
                            anchor_gains.append(float(gain))
                    if anchor_gains and min(anchor_gains)<.25:
                        continue
                    for distance,right,operator in transformed_images(right_original,left[-1,1],density,min(12.5,3.9*(gap+1))):
                        progress('Loop hypothesis %d -> %d: %d missing residues, CA %.2f A'%(i+1,j+1,gap,distance))
                        candidates,stats=loop_candidates(left,right,gap,density,target,missing_seq,
                            float(np.mean(anchor_gains)) if anchor_gains else 0.,starts,seed+len(attempts),
                            right_letter=seq[binding['end']-len(right)])
                        attempts.append({'left_fragment':i+1,'right_fragment':j+1,'gap':gap,'operator':operator,**stats})
                        if not candidates and bindings[j] is not None and len(right)>=5:
                            flex=2
                            right_start=binding['end']-len(right)
                            progress('Trying bounded flexibility in first two right-anchor residues')
                            candidates,stats=loop_candidates(left,right,gap,density,target,missing_seq,
                                float(np.mean(anchor_gains)) if anchor_gains else 0.,starts,seed+len(attempts),
                                right_letter=seq[right_start+flex],flexible_right=flex,
                                flexible_sequence=seq[right_start:right_start+flex])
                            attempts.append({'left_fragment':i+1,'right_fragment':j+1,'gap':gap,'operator':operator,**stats})
                        for candidate in candidates:
                            merged=candidate['xyz'];other=[c for k,c in enumerate(chains) if k not in [i,j]]
                            if any(density.distance(a[1],b[1])<2.8 for a in candidate['missing'] for c in other for b in c):
                                continue
                            if any(density.distance(merged[a,1],merged[b,1])<2.8 for a in range(len(merged)) for b in range(a-1)):
                                continue
                            # Severe nonbonded backbone clashes; adjacent residues use peptide checks above.
                            bad=False
                            for a in range(len(left),len(left)+gap):
                                neighbors=[r for b,r in enumerate(merged) if abs(b-a)>1]+[r for c in other for r in c]
                                if any(density.distance(p,q)<1.8 for p in merged[a] for r in neighbors for q in r):
                                    bad=True;break
                            if not bad:
                                proposals.append((candidate['quality'],i,j,candidate,binding,operator))
        if not proposals:
            break
        proposals.sort(key=lambda item:item[0],reverse=True)
        quality,i,j,candidate,binding,operator=proposals[0]
        competitors=[p[0] for p in proposals[1:] if p[1:3]==(i,j) and
                     (p[3]['gap']!=candidate['gap'] or np.sqrt(np.mean(np.sum(
                         (p[3]['missing'][:,1]-candidate['missing'][:,1])**2,axis=-1)))>.7)]
        margin=quality-max(competitors) if competitors else None
        if margin is not None and margin<.15:
            attempts.append({'reason':'best loop is ambiguous','score_margin':margin})
            break
        accepted.append({k:v for k,v in candidate.items() if k not in ['xyz','missing']}|
                        {'sequence_binding':binding,'operator':operator,'score_margin':margin})
        chains=[c for k,c in enumerate(chains) if k not in [i,j]]+[candidate['xyz']]
        bindings=[b for k,b in enumerate(bindings) if k not in [i,j]]+[binding]
    return chains,{'method':'de novo torsion closure with fixed outer anchors and bounded proximal flexibility',
                   'max_gap':3,
                   'starts_per_hypothesis':starts,'attempts':attempts,'accepted':accepted,
                   'added_residues':sum(r['gap'] for r in accepted),
                   'scope':'No parent trace or target-reference coordinates. Weak or ambiguous loops remain absent.'}
