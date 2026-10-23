# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Validate desktop SAD jobs and create shell-free CLI arguments."""
from pathlib import Path
import math
import sys
import gemmi
from .data import parse_composition,read_mtz,read_atoms
from .reference_data import bundled_reference,read_reference
from .sequence import sequence_asu


def positive(value, name):
    number=float(value)
    if not math.isfinite(number) or number<=0:raise ValueError(f'{name} 必须是正数')
    return number


def build_sad_command(spec):
    source=Path(spec['input']).expanduser().resolve()
    heavy=Path(spec['heavy']).expanduser().resolve()
    out=Path(spec['output']).expanduser().resolve()
    if not source.is_file() or not heavy.is_file():raise ValueError('请选择存在的实验 MTZ 和重原子坐标文件')
    if out.suffix.lower()!='.mtz':raise ValueError('输出文件名必须以 .mtz 结尾')
    for path in [out,out.with_suffix('.json'),out.with_suffix('.run.log')]:
        if path.exists():raise ValueError(f'文件已存在，请换一个输出名称：{path}')
    if out in [source,heavy]:raise ValueError('输出不能覆盖输入文件')
    labels=spec.get('labels',{})
    dmin=positive(spec['dmin'],'分辨率') if spec.get('dmin') else None
    data,_=read_mtz(source,'sad',labels,dmin=dmin)
    if len(data.f)<100:raise ValueError('有效 SAD 反射不足 100 条，请检查数据列或分辨率')
    atoms=read_atoms(heavy,data)
    sequence_meta=None
    if 'sequence' in spec:
        sequence_meta=sequence_asu(spec['sequence'],spec.get('copies',1),data.cell,data.spacegroup)
        composition=None;residues=0
    else:
        composition=spec.get('composition','').strip()
        residues=int(spec.get('residues') or 0)
        parse_composition(composition or None,residues)
    fpp={}
    for item in spec['fpp'].split(','):
        el,value=item.split('=',1);element=gemmi.Element(el.strip())
        if element.atomic_number==0:raise ValueError(f'未知元素：{el}')
        fpp[element.name]=positive(value,"f''")
    missing=set(atoms.elements)-set(fpp)
    if missing:raise ValueError("重原子文件的以下元素没有 f''："+', '.join(sorted(missing)))
    method=spec.get('method','likelihood')
    if method not in {'likelihood','pplus'}:raise ValueError('未知 SAD 方法')
    cycles=int(spec.get('cycles',5));dmcycles=int(spec.get('dm_cycles',5))
    if not 1<=cycles<=100 or not 1<=dmcycles<=100:raise ValueError('迭代轮数必须在 1–100 之间')
    command=[sys.executable,'-u','-m','oasis_py','phase',str(source),'--mode','sad',
        '--sad-method',method,'--heavy',str(heavy),'--fpp',','.join(f'{k}={v}' for k,v in fpp.items()),
        '--cycles',str(cycles),'--seed',str(int(spec.get('seed',20261016))),
        '--max-neighbors','600','--max-per-target',str(max(1,min(96,1950000//len(data.f)))),
        '--out',str(out)]
    if sequence_meta:command+=['--sequence',sequence_meta['path'],'--copies',str(sequence_meta['copies'])]
    elif composition:command+=['--composition',composition]
    else:command+=['--residues',str(residues)]
    if labels:command+=['--labels',','.join(f'{k}={v}' for k,v in labels.items())]
    if dmin:command+=['--dmin',str(dmin)]
    if spec.get('density_modification',True):
        solvent=sequence_meta['solvent_fraction'] if sequence_meta else positive(spec.get('solvent',''),'溶剂比例')
        if solvent>=1:raise ValueError('溶剂比例必须在 0 和 1 之间（例如 0.5）')
        reference=Path(spec['reference']).resolve() if spec.get('reference') else bundled_reference()
        read_reference(reference)
        command+=['--native-dm','--dm-engine','mlhl-best' if method=='likelihood' else 'legacy',
            '--dm-cycles',str(dmcycles),'--reference-mtz',str(reference)]
        if sequence_meta is None:command+=['--solvent-fraction',str(solvent)]
    return command
