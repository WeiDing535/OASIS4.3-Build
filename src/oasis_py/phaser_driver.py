# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Run with phenix.python, using Phaser's documented experimental-phasing API.

This is an external maximum-likelihood backend adapter. The likelihood engine
is Phaser/CCTBX; it is not represented as a Python reimplementation.
"""
from __future__ import print_function
import json
import sys
import phaser

def native_strings(value):
    # Some installed Phenix releases still embed Python 2 and their Boost
    # bindings accept byte strings rather than JSON's unicode objects.
    if sys.version_info[0] < 3:
        if isinstance(value, unicode):
            return value.encode('utf-8')
        if isinstance(value, dict):
            return {native_strings(k):native_strings(v) for k,v in value.items()}
        if isinstance(value, list):
            return [native_strings(v) for v in value]
    return value

def require_success(result,stage):
    if not result.Success():
        raise RuntimeError('%s: %s: %s' % (stage,result.ErrorName(),result.ErrorMessage()))


def main():
    with open(sys.argv[1]) as handle:job=native_strings(json.load(handle))
    xtal,wave='crystal','sad'
    inp=phaser.InputEP_DAT();inp.setHKLI(job['hklin'])
    inp.addCRYS_ANOM_LABI(xtal,wave,*job['labels'])
    inp.setMUTE(False)
    data=phaser.runEP_DAT(inp);require_success(data,'EP_DAT')
    if job.get('search'):
        search=phaser.InputEP_SSD()
        search.setCRYS_DATA(data.getCrysData())
        search.setSPAC_HALL(data.getSpaceGroupHall());search.setCELL6(data.getUnitCell())
        search.setWAVE(job['wavelength']);search.setFIND_SCAT(job['element']);search.setFIND_NUMB(job['sites'])
        search.setJOBS(job['jobs']);search.setROOT(job['root']+'_search')
        if job.get('dmin'):search.setRESO_HIGH(job['dmin'])
        found=phaser.runEP_SSD(search);require_success(found,'EP_SSD')
        if not found.number_of_solutions():raise RuntimeError('EP_SSD found no substructure solutions')
        substructure=found.getDotSol()
        # EP_SSD in some installed releases mutates reflection selections
        # through shared data objects. Re-read observations before EP_AUTO.
        fresh=phaser.InputEP_DAT();fresh.setHKLI(job['hklin'])
        fresh.addCRYS_ANOM_LABI(xtal,wave,*job['labels']);fresh.setMUTE(True)
        data=phaser.runEP_DAT(fresh);require_success(data,'EP_DAT reload')
    inp=phaser.InputEP_AUTO()
    inp.setSPAC_HALL(data.getSpaceGroupHall());inp.setCELL6(data.getUnitCell())
    # Preserve Phaser's complete reflection dataset. EP_SSD may update its
    # internal reflection selection; reconstructing individual arrays can
    # otherwise mix selections with different lengths.
    inp.setCRYS_DATA(data.getCrysData())
    if job.get('search'):inp.setATOM(substructure)
    else:inp.setATOM_PDB(xtal,job['heavy'])
    inp.addCOMP_PROT_SEQ_NUM(job['sequence'],float(job['copies']))
    inp.setWAVE(job['wavelength']);inp.setHAND('both');inp.setJOBS(job['jobs'])
    inp.setROOT(job['root']);inp.setHKLO(True);inp.setXYZO(True);inp.setMUTE(False)
    inp.setLLGC_COMP(True)
    for element in job['elements']:inp.addLLGC_SCAT(element)
    if job.get('dmin'):inp.setRESO_HIGH(job['dmin'])
    result=phaser.runEP_AUTO(inp);require_success(result,'EP_AUTO')
    summary={'backend':'Phaser EP_AUTO maximum likelihood (external)',
             'substructure_backend':'Phaser EP_SSD (external)' if job.get('search') else 'supplied model',
             'log_likelihood':float(result.getLogLikelihood()),'mtz':result.getMtzFile(),
             'hands':[]}
    for hand in [False,True]:
        summary['hands'].append({'inverted':hand,'fom':float(result.stats_fom(hand)),
                                 'reflections':int(result.stats_num(hand))})
    with open(job['summary'],'w') as handle:json.dump(summary,handle,indent=2)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print('PHASER BACKEND ERROR:',str(exc),file=sys.stderr);sys.exit(1)
