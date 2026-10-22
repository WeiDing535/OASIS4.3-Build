import json
import subprocess
import gemmi
import numpy as np
import pytest
from phasebuild.rfactors import atomic_factors,fit_amplitudes,calculate_rfactors


def synthetic(tmp_path):
    st=gemmi.Structure();st.cell=gemmi.UnitCell(18,20,22,90,90,90);st.spacegroup_hm='P 21 21 21'
    model=gemmi.Model('1');chain=gemmi.Chain('A');res=gemmi.Residue();res.name='ALA';res.seqid=gemmi.SeqId(1,' ')
    for name,element,pos in [('N','N',(3.1,4.3,5.2)),('CA','C',(4.4,4.6,5.7)),('C','C',(5.1,5.7,6.3)),('O','O',(5.3,6.5,5.6))]:
        atom=gemmi.Atom();atom.name=name;atom.element=gemmi.Element(element);atom.pos=gemmi.Position(*pos);atom.b_iso=20;atom.occ=1;res.add_atom(atom)
    chain.add_residue(res);model.add_chain(chain);st.add_model(model);st.setup_cell_images()
    hkl=gemmi.make_miller_array(st.cell,st.find_spacegroup(),2.,7.)[:100]
    fc=atomic_factors(st,hkl,2.,method='direct');x=np.array([st.cell.calculate_1_d2(h) for h in hkl])/4
    fo=2.4*np.exp(-11*x)*np.abs(fc)
    flags=np.ones(len(hkl));flags[::5]=0
    mtz=gemmi.Mtz(with_base=True);mtz.spacegroup=st.find_spacegroup();mtz.set_cell_for_all(st.cell);mtz.add_dataset('test')
    for label,kind in [('FP','F'),('SIGFP','Q'),('TEST','I')]:mtz.add_column(label,kind)
    mtz.set_data(np.column_stack([hkl,fo,np.ones(len(hkl)),flags]).astype(np.float32))
    mp=tmp_path/'data.mtz';sp=tmp_path/'model.cif';mtz.write_to_file(str(mp));st.make_mmcif_document().write_file(str(sp))
    return st,hkl,fc,x,fo,flags,mp,sp


def test_fft_matches_direct_with_crystal_symmetry(tmp_path):
    st,hkl,fc,*_=synthetic(tmp_path)
    fft=atomic_factors(st,hkl,2.,'fft')
    assert np.linalg.norm(fft-fc)/np.linalg.norm(fc)<.002


def test_work_only_scale_recovers_known_model_and_ignores_free_amplitudes(tmp_path):
    _,_,fc,x,fo,flags,*_=synthetic(tmp_path);work=flags!=0
    predicted,p=fit_amplitudes(fo,fc,x,work)
    assert p['k_overall']==pytest.approx(2.4,rel=1e-5)
    assert p['b_overall']==pytest.approx(11,abs=1e-4)
    np.testing.assert_allclose(predicted,fo,rtol=1e-6)
    changed=fo.copy();changed[~work]*=3
    other,q=fit_amplitudes(changed,fc,x,work)
    assert p==q
    np.testing.assert_array_equal(predicted,other)


def test_bulk_solvent_scale_recovers_synthetic_parameters_without_free_fitting(tmp_path):
    _,_,fc,x,_,flags,*_=synthetic(tmp_path)
    mask=np.exp(1j*np.arange(len(fc))*.7)*70
    fo=1.7*np.exp(-7*x)*np.abs(fc+.32*np.exp(-40*x)*mask)
    pred,params=fit_amplitudes(fo,fc,x,flags!=0,mask)
    np.testing.assert_allclose(pred,fo,rtol=1e-5)
    assert params['k_sol']==pytest.approx(.32,abs=1e-4)
    assert params['b_sol']==pytest.approx(40,abs=.01)


def test_standalone_report_needs_no_external_process_and_preserves_model(tmp_path,monkeypatch):
    st,_,_,_,_,_,mp,sp=synthetic(tmp_path)
    monkeypatch.setattr(subprocess,'Popen',lambda *a,**kw:pytest.fail('No external process permitted'))
    report=calculate_rfactors(mp,sp,tmp_path/'out',fp='FP',sigfp='SIGFP',free='TEST',
                             bulk_solvent=False,method='direct',progress=lambda s:None)
    assert report['statistics']['final']['r_work']<1e-6
    assert report['statistics']['final']['r_free']<1e-6
    assert report['external_refmac_used'] is False and report['coordinates_refined'] is False
    original=[list(cra.atom.pos) for cra in st[0].all()]
    output=[list(cra.atom.pos) for cra in gemmi.read_structure(str(tmp_path/'out/model.cif'))[0].all()]
    np.testing.assert_allclose(output,original,atol=1e-6)
    assert json.loads((tmp_path/'out/r_factor_report.json').read_text())['status']=='success'
