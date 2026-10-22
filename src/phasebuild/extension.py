"""Sequence-anchored terminal beam search; no target-reference coordinates."""
import gemmi
import numpy as np
from scipy.optimize import minimize
from .geometry import extend,prepend,oxygen
from .validation import peptide_metrics,peptide_valid
from .residue_profile import ALPHABET


def free_positions(binding,bindings,records,direction,maximum=12):
    word=records[binding['record']][1]
    occupied={i for b in bindings if b is not None and b is not binding and b['record']==binding['record']
              for i in range(b['start'],b['end'])}
    start=binding['end'] if direction=='C' else binding['start']-1
    step=1 if direction=='C' else -1
    result=[]
    while 0<=start<len(word) and start not in occupied and len(result)<maximum:
        result.append(start);start+=step
    return result


def phi_grid(letter):
    if letter=='P':return [-85.,-65.,-45.]
    result=[-150.,-120.,-90.,-63.]
    if letter in 'GDNSTX':result += [60.,90.,120.]
    return result


def phi_bounds(letter):
    return (-100.,-30.) if letter=='P' else (-179.9,179.9 if letter=='G' else 80. if letter in 'DNSTX' else 0.)


def image_atoms(points,center,density):
    """Expand complete crystal images around an anchor, preserving each atom's image."""
    points=np.asarray(points).reshape(-1,3)
    if not len(points):return points
    output=[]
    for p in points:
        f=density.cell.fractionalize(gemmi.Position(*p))
        for op in density.sg.operations():
            image=density.cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(list(f))))
            nearest=density.cell.find_nearest_pbc_position(gemmi.Position(*center),image,0)
            xyz=np.array(list(nearest))
            if np.linalg.norm(xyz-center)<55.:output.append(xyz)
    return np.array(output).reshape(-1,3)


def terminal_candidates(chain,positions,word,direction,density,target,external=None,
                        beam_width=6,minimum_profile=2.5,minimum_atom=.2,maximum_refinements=8):
    """Keep alternate prefixes; a weak residue must be supported by later growth."""
    if not positions:return [],{'steps':[],'reason':'no free sequence positions'}
    external=np.empty((0,3)) if external is None else np.asarray(external).reshape(-1,3)
    beams=[(0.,chain.copy(),[],[])];accepted=[];steps=[]
    psi_values=[-180.,-150.,-120.,-90.,-60.,-42.,0.,40.,80.,120.,150.,180.]
    def identity(xyz,letter):
        if target is None or letter not in ALPHABET:return np.zeros(xyz.shape[:-2])
        return target.scores(density,xyz)[...,ALPHABET.index(letter)]
    for depth,position in enumerate(positions):
        branches=[];letter=word[position]
        for cumulative,current,profiles,identities in beams:
            anchor=current[-1] if direction=='C' else current[0]
            anchor_letter=word[position-1] if direction=='C' else word[position+1]
            def coordinates(x):
                if direction=='C':
                    new=extend(anchor,x[0],x[1]);new[3]=oxygen(*new[:3],x[2])
                    changed=anchor.copy();changed[3]=oxygen(*changed[:3],x[0])
                    return new,changed
                return prepend(anchor,x[0],x[1]),anchor
            if direction=='C':
                seeds=[[psi,phi,terminal] for psi in psi_values for phi in phi_grid(letter) for terminal in [-42.,40.,130.]]
                bounds=[(-180.,180.),phi_bounds(letter),(-180.,180.)]
            else:
                seeds=[[phi,psi] for phi in phi_grid(anchor_letter) for psi in psi_values]
                bounds=[phi_bounds(anchor_letter),(-180.,180.)]
            seeds=np.asarray(seeds)
            variants=np.array([coordinates(x)[0] for x in seeds])
            quality=density.score(variants)+.25*identity(variants,letter)
            order=np.argsort(quality)[::-1][:maximum_refinements*(2 if depth==0 else 1)]
            for index in order:
                def loss(x):
                    new,changed=coordinates(x)
                    return -float(density.score(new)+.25*identity(new,letter)+
                                  (.15*density.score(changed) if direction=='C' else 0.))
                seed=np.array([np.clip(v,lo+1e-4,hi-1e-4) for v,(lo,hi) in zip(seeds[index],bounds)])
                fitted=minimize(loss,seed,method='L-BFGS-B',bounds=bounds,
                                options={'maxiter':22,'eps':.12,'ftol':1e-5})
                new,changed=coordinates(fitted.x)
                profile=float(density.score(new));seq=float(identity(new,letter))
                if profile<1.5 or np.min(density.sample(new[:3]))<minimum_atom:
                    continue
                if direction=='C' and float(density.sample(changed[3]))<minimum_atom:
                    continue
                link=peptide_metrics(changed,new) if direction=='C' else peptide_metrics(new,anchor)
                if not peptide_valid(link):continue
                nonadjacent=current[:-1] if direction=='C' else current[1:]
                environment=np.concatenate([external,nonadjacent.reshape(-1,3)])
                if len(environment) and np.min(np.linalg.norm(new[:,None,:]-environment[None,:,:],axis=-1))<1.8:
                    continue
                # Explicit symmetry-aware CA exclusion, including growing chain images.
                if any(density.distance(new[1],r[1])<2.8 for r in nonadjacent):continue
                if direction=='C':
                    merged=np.concatenate([current[:-1],changed[None],new[None]])
                else:
                    merged=np.concatenate([new[None],current])
                values=profiles+[profile];types=identities+[seq]
                if len(values)>=2 and values[-1]<minimum_profile and values[-2]<minimum_profile:continue
                gain=cumulative+profile-minimum_profile+.25*seq
                # Preserve CA-distinct alternatives; duplicated torsion solutions waste beam slots.
                if any(np.sqrt(np.mean(np.sum((merged[:,1]-b[1][:,1])**2,axis=-1)))<.35 for b in branches):continue
                branches.append((gain,merged,values,types))
        if not branches:
            steps.append({'position':position+1,'candidates':0});break
        branches.sort(key=lambda x:x[0],reverse=True)
        width=beam_width*2 if depth==0 else beam_width
        beams=branches[:width]
        steps.append({'position':position+1,'candidates':len(branches),'kept':len(beams),
                      'best_profile':beams[0][2][-1],'best_gain':beams[0][0]})
        for gain,xyz,profiles,types in beams:
            n=len(profiles)
            if n>=min(2,len(positions)) and profiles[-1]>=minimum_profile and np.mean(profiles)>=minimum_profile and gain>.8:
                accepted.append({'xyz':xyz,'added':n,'gain':gain,'profiles':profiles,'identities':types,
                                 'positions':positions[:n]})
    accepted.sort(key=lambda c:c['gain'],reverse=True)
    return accepted,{'steps':steps,'direction':direction,'beam_width':beam_width,
                     'minimum_profile':minimum_profile,'minimum_atom':minimum_atom,
                     'accepted_prefixes':len(accepted)}


def extend_anchors(chains,bindings,records,density,target,maximum=12,beam_width=6,progress=print):
    chains=[c.copy() for c in chains];bindings=[None if b is None else dict(b) for b in bindings]
    events=[];attempts=[]
    # Extend ends with the largest free intervals first; occupancy is recomputed after every acceptance.
    jobs=sorted([(len(free_positions(b,bindings,records,d,maximum)),i,d)
                 for i,b in enumerate(bindings) if b is not None for d in ['C','N']],reverse=True)
    for _,index,direction in jobs:
        binding=bindings[index]
        positions=free_positions(binding,bindings,records,direction,maximum)
        if not positions:continue
        word=records[binding['record']][1]
        progress('Terminal growth fragment %d %s: at most %d free positions'%(index+1,direction,len(positions)))
        original=chains[index]
        # Try the current terminal and a one-residue repair with the next interior residue fixed.
        variants=[(original,positions,0)]
        if len(original)>=5:
            repair_position=binding['end']-1 if direction=='C' else binding['start']
            variants.append((original[:-1] if direction=='C' else original[1:],
                             [repair_position]+positions,1))
        proposals=[]
        for anchor,free,repair in variants:
            external=np.concatenate([c.reshape(-1,3) for j,c in enumerate(chains) if j!=index]) if len(chains)>1 else np.empty((0,3))
            center=anchor[-1,1] if direction=='C' else anchor[0,1]
            external=image_atoms(external,center,density)
            candidates,stats=terminal_candidates(anchor,free,word,direction,density,target,external,
                                                 beam_width=beam_width)
            attempts.append({'fragment':index+1,'direction':direction,'repair_residues':repair,**stats})
            for candidate in candidates:
                added=candidate['added']-repair
                if added<1:continue
                if repair:
                    old=original[-1] if direction=='C' else original[0]
                    replacement=candidate['xyz'][len(original)-1] if direction=='C' else candidate['xyz'][added]
                    if float(density.score(replacement))<float(density.score(old))+.1:continue
                proposals.append((candidate['gain']-repair*max(0.,float(density.score(original[-1 if direction=='C' else 0]))-2.5),
                                  candidate,added,repair))
        if not proposals:continue
        proposals.sort(key=lambda x:x[0],reverse=True)
        gain,candidate,added,repair=proposals[0]
        chains[index]=candidate['xyz']
        binding['end' if direction=='C' else 'start']+=added*(1 if direction=='C' else -1)
        events.append({'fragment':index+1,'direction':direction,'added_residues':added,
                       'repaired_residues':repair,'gain':gain,'profiles':candidate['profiles'],
                       'sequence_binding':dict(binding)})
    return chains,bindings,{'method':'sequence-anchored bidirectional terminal beam search with one-residue repair',
                            'events':events,'attempts':attempts,'added_residues':sum(e['added_residues'] for e in events),
                            'scope':'No target reference, Buccaneer coordinates or Rfree. Fixed interior anchors; trans peptides only.'}
