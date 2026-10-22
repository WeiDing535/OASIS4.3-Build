import json
import xml.etree.ElementTree as ET
import gemmi
import numpy as np
import pytest
from phasebuild.refinement import reflection_partition,verify_free_flags,parse_refmac,refmac_refine


def mtz(rows):
    m=gemmi.Mtz(with_base=True);m.spacegroup=gemmi.find_spacegroup_by_name('P 1')
    m.set_cell_for_all(gemmi.UnitCell(40,40,40,90,90,90));m.add_dataset('test')
    for name,kind in [('F','F'),('SIGF','Q'),('FreeR_flag','I')]:
        m.add_column(name,kind)
    m.set_data(np.array(rows,dtype=np.float32));m.update_reso()
    return m


def test_partition_keeps_existing_flag_zero_and_rejects_missing_observed_flags(tmp_path):
    path=tmp_path/'data.mtz'
    m=mtz([[1,1,0,10,1,0],[2,1,0,20,1,5],[3,1,0,np.nan,np.nan,0]])
    m.write_to_file(str(path))
    _,p=reflection_partition(path)
    assert p['usable_free']==1 and p['usable_work']==1 and not p['flags_regenerated']
    m=mtz([[1,1,0,10,1,np.nan],[2,1,0,20,1,5]])
    m.write_to_file(str(path))
    with pytest.raises(ValueError,match='integer free-set flag'):
        reflection_partition(path)


def test_audit_allows_only_unobserved_omissions_and_no_flag_changes():
    a=mtz([[1,1,0,10,1,0],[2,1,0,20,1,5],[3,1,0,np.nan,np.nan,3]])
    b=mtz([[1,1,0,10,1,0],[2,1,0,20,1,5]])
    report=verify_free_flags(a,b,'FreeR_flag','F','SIGF')
    assert report['removed_unobserved_hkls']==[[3,1,0]]
    json.dumps(report,allow_nan=False)
    b=mtz([[1,1,0,10,1,0]])
    with pytest.raises(RuntimeError,match='observed reflection'):
        verify_free_flags(a,b,'FreeR_flag','F','SIGF')
    b=mtz([[1,1,0,10,1,1],[2,1,0,20,1,5]])
    with pytest.raises(RuntimeError,match='changed an original'):
        verify_free_flags(a,b,'FreeR_flag','F','SIGF')


def test_r_factors_must_be_consistent_between_log_and_xml(tmp_path):
    path=tmp_path/'refmac.xml'
    root=ET.Element('REFMAC');ET.SubElement(root,'err_level').text='0'
    overall=ET.SubElement(root,'Overall_stats')
    for key,value in [('n_reflections_all',100),('n_reflections_free',5),('resolution_low',20),('resolution_high',2)]:
        ET.SubElement(overall,key).text=str(value)
    vs=ET.SubElement(overall,'stats_vs_cycle');cycle=ET.SubElement(vs,'new_cycle')
    for key,value in [('cycle',0),('r_factor',.437),('r_free',.433)]:
        ET.SubElement(cycle,key).text=str(value)
    ET.ElementTree(root).write(str(path))
    log='Ncyc    Rfact    Rfree\n0 0.4375 0.4333 0.58 3000 150 0.0148 1.27 1.56 0.89 0.1\n$TEXT:'
    report=parse_refmac(log,path)
    assert report['final']['r_work']==.4375 and report['reflections_used_work']==95
    with pytest.raises(RuntimeError,match='values disagree'):
        parse_refmac(log.replace('0.4333','0.6000'),path)
    with pytest.raises(RuntimeError,match='cycle table missing'):
        parse_refmac('normal exit without statistics',path)


def test_refinement_adapter_exports_auditable_zero_cycle_result(tmp_path,monkeypatch):
    from pathlib import Path
    from types import SimpleNamespace
    from phasebuild.geometry import template
    from phasebuild.trace import write_backbone
    import phasebuild.refinement as refinement
    data=mtz([[1,1,0,10,1,0],[2,1,0,20,1,5],[3,1,0,np.nan,np.nan,3]])
    native=tmp_path/'native.mtz';data.write_to_file(str(native))
    model=tmp_path/'input.cif'
    write_backbone(model,[template(-63.,-42.,length=2)+[20.,20.,20.]],
                   SimpleNamespace(cell=data.cell,sg=data.spacegroup))
    ccp4=tmp_path/'ccp4';(ccp4/'bin').mkdir(parents=True)
    (ccp4/'bin/refmac5').write_bytes(b'fake executable, not run')
    (ccp4/'bin/refmac5').chmod(0o755)
    def fake_run(command,**kwargs):
        assert 'FREE 0\n' in kwargs['input'] and 'SCAL MLSC WORK' in kwargs['input']
        assert 'NCYC 0' in kwargs['input'] and 'MAKE HYDR N HOUT N BUIL N' in kwargs['input']
        out=Path(kwargs['cwd'])
        gemmi.read_structure(str(out/'input.pdb')).write_pdb(str(out/'refined.pdb'))
        mtz([[1,1,0,10,1,0],[2,1,0,20,1,5]]).write_to_file(str(out/'refined.mtz'))
        (out/'refmac.xml').write_text('''<REFMAC><err_level>0</err_level><Overall_stats>
<n_reflections_all>2</n_reflections_all><n_reflections_free>1</n_reflections_free>
<resolution_low>28.3</resolution_low><resolution_high>17.9</resolution_high>
<stats_vs_cycle><new_cycle><cycle>0</cycle><r_factor>0.437</r_factor><r_free>0.433</r_free>
</new_cycle></stats_vs_cycle></Overall_stats></REFMAC>''')
        kwargs['stdout'].write('Ncyc    Rfact    Rfree\n0 0.4375 0.4333 0.58 3000 150 0.0148 1.27 1.56 0.89 0.1\n$TEXT:')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(refinement.subprocess,'run',fake_run)
    out=tmp_path/'result'
    report=refmac_refine(native,model,out,cycles=0,ccp4=ccp4)
    saved=json.loads((out/'refinement_report.json').read_text())
    assert saved['status']=='success' and report['statistics']['final']['r_free']==.4333
    assert saved['free_flag_audit']['removed_unobserved_hkls']==[[3,1,0]]
    assert not saved['atom_geometry']['severe_nonbonded_contacts_below_1p8A']
    assert (out/'model.cif').is_file()
