# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Optional mature algorithms accessed via isolated external Python runtimes."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import gemmi
import numpy as np


def prepare_phaser_pairs(hklin, labels, output):
    """Remove unusable Bijvoet pairs before EP_DAT; retain the full MTZ schema."""
    if len(labels)!=4:raise ValueError('Phaser requires four Bijvoet amplitude/sigma labels')
    mtz=gemmi.read_mtz_file(str(hklin))
    columns=[]
    for label,types in zip(labels,['FG','QL','FG','QL']):
        column=mtz.column_with_label(label)
        if column is None or column.type not in types:
            raise ValueError(f'Invalid amplitude/sigma column for Phaser: {label}')
        columns.append(np.asarray(column.array))
    pairs=np.column_stack(columns)
    valid=np.all(np.isfinite(pairs)&(pairs>0),axis=1)
    if not valid.any():raise ValueError('No complete positive Bijvoet pairs for Phaser')
    stats={'input_rows':len(valid),'retained_rows':int(valid.sum()),
           'removed_rows':int((~valid).sum()),'labels':list(labels)}
    if not valid.all():
        mtz.set_data(np.array(mtz.array[valid],copy=True))
        mtz.write_to_file(str(output))
        return str(Path(output).resolve()),stats
    return str(Path(hklin).resolve()),stats


def phaser_sad(phenix_python,hklin,sequence,out_dir,wavelength,heavy=None,
               sites=None,element='S',copies=1,labels=None,jobs=1,dmin=None,timeout=1800):
    if wavelength<=0 or copies<1 or jobs<1:raise ValueError('Positive wavelength, copies and jobs required')
    if (heavy is None)==(sites is None):raise ValueError('Supply either heavy atoms or a site count for EP_SSD search')
    if sites is not None and sites<1:raise ValueError('sites must be positive')
    if timeout<=0:raise ValueError('timeout must be positive')
    out=Path(out_dir).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Phaser output directory must be empty')
    for path in [phenix_python,hklin,sequence]+([heavy] if heavy else []):
        if not Path(path).is_file():raise ValueError(f'File not found: {path}')
    out.mkdir(parents=True,exist_ok=True)
    pair_labels=labels or ['F(+)','SIGF(+)','F(-)','SIGF(-)']
    prepared,input_filter=prepare_phaser_pairs(hklin,pair_labels,out/'input_filtered.mtz')
    (out/'input_filter.json').write_text(json.dumps(input_filter,indent=2)+'\n')
    job={'hklin':prepared,'sequence':str(Path(sequence).resolve()),
         'wavelength':wavelength,'copies':copies,'jobs':jobs,'dmin':dmin,
         'heavy':str(Path(heavy).resolve()) if heavy else None,'search':sites is not None,
         'sites':sites,'element':element,'elements':[element],
         'labels':pair_labels,
         'root':str(out/'phaser'),'summary':str(out/'summary.json')}
    job_path=out/'job.json';job_path.write_text(json.dumps(job,indent=2)+'\n')
    driver=Path(__file__).with_name('phaser_driver.py')
    command=[str(Path(phenix_python).resolve()),str(driver),str(job_path)]
    with (out/'phaser.log').open('w') as log:
        try:
            run=subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=False,timeout=timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f'Phaser exceeded {timeout}s and was terminated; see {out}/phaser.log')
    if run.returncode:raise RuntimeError(f'Phaser failed with exit {run.returncode}; see {out}/phaser.log')
    if not (out/'summary.json').is_file():raise RuntimeError('Phaser did not write its completion summary')
    summary=json.loads((out/'summary.json').read_text())
    summary['input_filter']=input_filter
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
