# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei and contributors.
"""Protein sequence composition and Matthews solvent estimate for a specified ASU.

Vm = Vcell / (Nsym * Masu); solvent fraction = 1 - 1.23/Vm.
Protein residue formulas are polymerized (free amino acid minus water);
one H2O is restored per FASTA/PIR chain, including repeated records.
"""
from collections import Counter
from pathlib import Path
import hashlib
import math
import re
import gemmi

# C, H, N, O, S in a peptide residue. Charges, disulfides, ligands and
# non-sequence modifications are not inferred from a primary sequence.
RESIDUE_FORMULAS={
    'A':(3,5,1,1,0),'R':(6,12,4,1,0),'N':(4,6,2,2,0),'D':(4,5,1,3,0),
    'C':(3,5,1,1,1),'E':(5,7,1,3,0),'Q':(5,8,2,2,0),'G':(2,3,1,1,0),
    'H':(6,7,3,1,0),'I':(6,11,1,1,0),'L':(6,11,1,1,0),'K':(6,12,2,1,0),
    'M':(5,9,1,1,1),'F':(9,9,1,1,0),'P':(5,7,1,1,0),'S':(3,5,1,2,0),
    'T':(4,7,1,2,0),'W':(11,10,2,1,0),'Y':(9,9,1,2,0),'V':(5,9,1,1,0),
    'O':(12,19,3,2,0),
}
MODIFIED_FORMULAS={'MSE':{'C':5,'H':9,'N':1,'O':1,'Se':1},
                   'SEC':{'C':3,'H':5,'N':1,'O':1,'Se':1}}


def protein_sequence(path):
    path=Path(path).expanduser().resolve()
    if not path.is_file():raise ValueError('Please select an existing sequence file (FASTA, PIR, or plain text)')
    if path.stat().st_size>5_000_000:raise ValueError('Sequence file exceeds 5 MB; check that the correct file was selected')
    raw=path.read_bytes()
    try:text=raw.decode('utf-8-sig').strip()
    except UnicodeDecodeError as e:raise ValueError('Sequence file must be UTF-8 text; Word/PDF files are not accepted') from e
    if not text:raise ValueError('Sequence file is empty')
    try:records=gemmi.read_pir_or_fasta(text if text.startswith('>') else '>sequence\n'+text)
    except (RuntimeError,ValueError) as e:raise ValueError('Sequence format error: '+str(e)) from e
    composition=Counter();chains=[];total=0
    for record in records:
        header=record.header.splitlines()[0] if record.header else f'chain_{len(chains)+1}'
        seq=record.seq.upper();tokens=re.findall(r'\([A-Z0-9]+\)|[A-Z]',seq)
        if ''.join(tokens)!=seq or not tokens:raise ValueError('Sequence contains empty entries, alignment gaps, or unsupported characters; provide an ungapped protein sequence')
        for index,token in enumerate(tokens,1):
            if token=='U':formula=MODIFIED_FORMULAS['SEC']
            elif token.startswith('('):
                name=token[1:-1]
                if name in MODIFIED_FORMULAS:formula=MODIFIED_FORMULAS[name]
                else:
                    residue=gemmi.find_tabulated_residue(name)
                    letter=residue.one_letter_code
                    if not residue.is_amino_acid() or letter not in RESIDUE_FORMULAS:
                        raise ValueError(f'unsupported modified residue {token}; standard amino acids and (MSE)/(SEC) are supported')
                    # Only accept canonical names, never approximate other modifications.
                    if gemmi.expand_one_letter(letter,gemmi.ResidueKind.AA)!=name:raise ValueError(f'cannot infer accurate elemental composition from {token}')
                    formula=dict(zip(['C','H','N','O','S'],RESIDUE_FORMULAS[letter]))
            elif token in RESIDUE_FORMULAS:formula=dict(zip(['C','H','N','O','S'],RESIDUE_FORMULAS[token]))
            else:raise ValueError(f'sequence {header} position {index} is an unknown residue {token}; specify the residue type')
            composition.update(formula)
        composition.update({'H':2,'O':1})
        chains.append(dict(name=header,residues=len(tokens)))
        total+=len(tokens)
    if not chains:raise ValueError('No protein sequence found in the sequence file')
    composition={k:v for k,v in sorted(composition.items()) if v>0}
    weight=sum(gemmi.Element(k).weight*v for k,v in composition.items())
    return dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest(),chains=chains,
        chains_per_copy=len(chains),residues_per_copy=total,composition_per_copy=composition,mass_per_copy_da=weight)


MATTHEWS_SOLVENT_RANGE=(0.20,0.65)


CHAIN_DECLARATION=re.compile(r'Chains?\s+([A-Za-z0-9]+(?:[,\s]+[A-Za-z0-9]+)*)')


def declared_copies(path):
    """Number of chains declared in the RCSB FASTA header, if present.

    A header such as ``>6HF7_1|Chains A, B, C|...`` states that the deposited
    asymmetric unit holds three copies of this sequence. This is an input
    annotation, not a reference structure: it is read from the sequence file the
    user supplied.
    """
    path=Path(path).expanduser().resolve()
    text=path.read_text(errors='ignore')
    total=0;found=False
    for record in gemmi.read_pir_or_fasta(text if text.startswith('>') else '>sequence\n'+text):
        header=(record.header or '').replace('\r','')
        match=CHAIN_DECLARATION.search(header)
        if not match:continue
        letters=[token for token in re.split(r'[,\s]+',match.group(1)) if token]
        if letters:
            total+=len(letters);found=True
    return total if found else None


def infer_copies(path,cell,spacegroup,maximum=8):
    """Copies plus the evidence used: FASTA chain annotation or Matthews range."""
    declared=declared_copies(path)
    candidates=plausible_copies(path,cell,spacegroup,maximum)
    if declared is not None:
        row=next((item for item in candidates if item['copies']==declared),None)
        return dict(copies=declared,source='fasta_chain_annotation',
                    solvent=row['solvent'] if row else None,
                    matthews_plausible=bool(row and row['plausible']),
                    candidates=candidates)
    plausible=[item for item in candidates if item['plausible']]
    chosen=plausible[0]['copies'] if plausible else 1
    return dict(copies=chosen,source='matthews_plausible_range',
                solvent=next(item['solvent'] for item in candidates if item['copies']==chosen),
                matthews_plausible=bool(plausible),candidates=candidates)


def plausible_copies(path,cell,spacegroup,maximum=8):
    """Copy numbers giving a Matthews-plausible solvent fraction.

    Uses only the sequence, unit cell and space-group multiplicity: no
    experimental phases, reference coordinates or free-set information.
    """
    info=protein_sequence(path);nsym=len(spacegroup.operations())
    results=[]
    for copies in range(1,int(maximum)+1):
        mass=info['mass_per_copy_da']*copies
        vm=cell.volume/(nsym*mass)
        solvent=1-1.23/vm
        results.append(dict(copies=copies,vm=vm,solvent=solvent,
            plausible=bool(MATTHEWS_SOLVENT_RANGE[0]<=solvent<=MATTHEWS_SOLVENT_RANGE[1])))
    return results


def resolve_copies(path,cell,spacegroup,value):
    """Resolve 'auto' using the FASTA chain annotation, else the Matthews range.

    An explicit integer is returned unchanged; ``auto`` falls back to 1 when
    neither source is available so a warning, not a silent guess, drives the
    user's decision.
    """
    if value!='auto':return value
    return infer_copies(path,cell,spacegroup)['copies']


def sequence_asu(path,copies,cell,spacegroup):
    if copies=='auto':
        copies=resolve_copies(path,cell,spacegroup,'auto')
    if isinstance(copies,bool) or not isinstance(copies,int) or copies<1:
        raise ValueError('ASU copies must be a positive integer giving how many times the sequence file repeats in one ASU')
    if spacegroup is None or not cell.is_crystal() or not math.isfinite(cell.volume) or cell.volume<=0:
        raise ValueError('Experimental MTZ lacks a valid unit cell or space group; solvent content cannot be estimated')
    info=protein_sequence(path);nsym=len(spacegroup.operations())
    mass=info['mass_per_copy_da']*copies;vm=cell.volume/(nsym*mass);solvent=1-1.23/vm
    if not math.isfinite(solvent) or not 0<solvent<1:
        raise ValueError(f'sequence/copies incompatible with the unit cell: estimated solvent content {solvent:.1%}; check copies, sequence, and space group')
    warnings=[]
    if solvent<.2 or solvent>.85:warnings.append('Estimated solvent content is outside the typical protein-crystal range; check sequence and copies (copies were not changed automatically)')
    elif not MATTHEWS_SOLVENT_RANGE[0]<=solvent<=MATTHEWS_SOLVENT_RANGE[1]:
        plausible=[row for row in plausible_copies(path,cell,spacegroup) if row['plausible']]
        hint=', '.join('%d copy -> solvent %.0f%%'%(row['copies'],100*row['solvent']) for row in plausible[:3])
        warnings.append('current copies=%d gives a high estimated solvent fraction %.0f%%; ASU copies may be missing (%s). Use --copies auto or specify explicitly'
                       %(copies,100*solvent,hint or 'no suitable copy number'))
    if any(k in info['composition_per_copy'] for k in ['Se']):warnings.append('selenomethionine residues use explicit sequence tags for composition')
    return {**info,'copies':copies,'asu_residues':info['residues_per_copy']*copies,
        'asu_chains':info['chains_per_copy']*copies,'asu_mass_da':mass,
        'asu_composition':{k:v*copies for k,v in info['composition_per_copy'].items()},
        'symmetry_multiplicity':nsym,'cell':list(cell.parameters),'spacegroup':spacegroup.hm,
        'matthews_vm':vm,'solvent_fraction':solvent,'warnings':warnings,
        'method':'Vm=Vcell/(Nsym*ASU_mass); solvent=1-1.23/Vm; protein partial specific volume approximation',
        'assumptions':['Each file record is a distinct chain; copies repeats the entire file',
            'Unmodified linear protein chains; native methionine remains M, SeMet requires (MSE)',
            'Ligands, glycosylation, bound solvent and unspecified modifications are not included',
            'Solvent fraction is an estimate, not a measured value; copies is supplied, not inferred']}
