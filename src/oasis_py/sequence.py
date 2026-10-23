# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
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
    if not path.is_file():raise ValueError('请选择存在的蛋白质序列文件（FASTA、PIR 或单序列文本）')
    if path.stat().st_size>5_000_000:raise ValueError('序列文件超过 5 MB，请检查是否误选了其他文件')
    raw=path.read_bytes()
    try:text=raw.decode('utf-8-sig').strip()
    except UnicodeDecodeError as e:raise ValueError('序列文件须为 UTF-8 文本，不能使用 Word/PDF 文件') from e
    if not text:raise ValueError('序列文件为空')
    try:records=gemmi.read_pir_or_fasta(text if text.startswith('>') else '>sequence\n'+text)
    except (RuntimeError,ValueError) as e:raise ValueError('序列格式错误：'+str(e)) from e
    composition=Counter();chains=[];total=0
    for record in records:
        header=record.header.splitlines()[0] if record.header else f'chain_{len(chains)+1}'
        seq=record.seq.upper();tokens=re.findall(r'\([A-Z0-9]+\)|[A-Z]',seq)
        if ''.join(tokens)!=seq or not tokens:raise ValueError('序列含空条目、比对缺口或不支持的字符；请提供未比对的真实蛋白序列')
        for index,token in enumerate(tokens,1):
            if token=='U':formula=MODIFIED_FORMULAS['SEC']
            elif token.startswith('('):
                name=token[1:-1]
                if name in MODIFIED_FORMULAS:formula=MODIFIED_FORMULAS[name]
                else:
                    residue=gemmi.find_tabulated_residue(name)
                    letter=residue.one_letter_code
                    if not residue.is_amino_acid() or letter not in RESIDUE_FORMULAS:
                        raise ValueError(f'不支持的修饰残基 {token}；当前支持标准氨基酸及 (MSE)/(SEC)')
                    # Only accept canonical names, never approximate other modifications.
                    if gemmi.expand_one_letter(letter,gemmi.ResidueKind.AA)!=name:raise ValueError(f'不能从 {token} 推断准确元素组成')
                    formula=dict(zip(['C','H','N','O','S'],RESIDUE_FORMULAS[letter]))
            elif token in RESIDUE_FORMULAS:formula=dict(zip(['C','H','N','O','S'],RESIDUE_FORMULAS[token]))
            else:raise ValueError(f'序列 {header} 第 {index} 位为未知残基 {token}，请先明确残基种类')
            composition.update(formula)
        composition.update({'H':2,'O':1})
        chains.append(dict(name=header,residues=len(tokens)))
        total+=len(tokens)
    if not chains:raise ValueError('序列文件中没有蛋白质序列')
    composition={k:v for k,v in sorted(composition.items()) if v>0}
    weight=sum(gemmi.Element(k).weight*v for k,v in composition.items())
    return dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest(),chains=chains,
        chains_per_copy=len(chains),residues_per_copy=total,composition_per_copy=composition,mass_per_copy_da=weight)


def sequence_asu(path,copies,cell,spacegroup):
    if isinstance(copies,bool) or not isinstance(copies,int) or copies<1:
        raise ValueError('ASU copies 必须是正整数，表示整个序列文件在一个 ASU 中重复的份数')
    if spacegroup is None or not cell.is_crystal() or not math.isfinite(cell.volume) or cell.volume<=0:
        raise ValueError('实验 MTZ 中缺少有效晶胞或空间群，不能估算溶剂含量')
    info=protein_sequence(path);nsym=len(spacegroup.operations())
    mass=info['mass_per_copy_da']*copies;vm=cell.volume/(nsym*mass);solvent=1-1.23/vm
    if not math.isfinite(solvent) or not 0<solvent<1:
        raise ValueError(f'序列/copies 与晶胞不相容：估算溶剂含量 {solvent:.1%}，请检查 copies、序列和空间群')
    warnings=[]
    if solvent<.2 or solvent>.85:warnings.append('估算溶剂含量偏离常见蛋白晶体范围，请核对序列与 copies；程序未自动改变 copies')
    if any(k in info['composition_per_copy'] for k in ['Se']):warnings.append('硒代残基使用显式序列标记计算组成')
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
