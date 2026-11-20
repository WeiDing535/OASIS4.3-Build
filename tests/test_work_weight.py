from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,gemmi,pytest
from phasebuild.work_weight import mask_test_observations,select_geometry_weight,parse_work_geometry
import phasebuild.work_weight as ww
from test_sad_pipeline import mtz

def test_pilot_mask_preserves_original_flags_and_working_observations(tmp_path):
 path=tmp_path/'native.mtz';rows=np.array([[1,0,0,10,1,0],[0,1,0,20,2,3]],float);mtz(path,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],rows);before=path.read_bytes();report=mask_test_observations(path,tmp_path/'work.mtz');out=gemmi.read_mtz_file(str(tmp_path/'work.mtz'));a=np.array(out)
 assert np.isnan(a[0,3:5]).all();np.testing.assert_array_equal(a[:,5],rows[:,5]);np.testing.assert_array_equal(a[1],rows[1]);assert path.read_bytes()==before and report['original_test_rows_masked']==1

def test_geometry_policy_cannot_choose_by_Rfree():
 reports=[{'status':'completed','cycle':30,'weight':w,'bond_rmsd_angstrom':.003,'angle_rmsd_degrees':a,'r_free':rf} for w,a,rf in [(.005,.9,.2),(.01,1.,.5),(.02,1.2,.6),(.04,1.5,.1)]]
 assert select_geometry_weight(reports)==.02
 for r in reports:r['r_free']=0 if r['weight']==.04 else 1
 assert select_geometry_weight(reports)==.02
 assert select_geometry_weight([])==.005

def xml(tmp_path,nfree):
 r=ET.Element('REFMAC');ET.SubElement(r,'err_level').text='0';overall=ET.SubElement(r,'Overall_stats');ET.SubElement(overall,'n_reflections_free').text=str(nfree);cycles=ET.SubElement(overall,'stats_vs_cycle')
 for i in range(31):
  c=ET.SubElement(cycles,'new_cycle')
  for k,v in [('cycle',i),('r_factor',.25),('rmsBOND',.004),('rmsANGLE',1.1)]:ET.SubElement(c,k).text=str(v)
 p=tmp_path/'refmac.xml';ET.ElementTree(r).write(p);log=tmp_path/'refmac.log';log.write_text('Ncyc    Rfact    Rfree\n30 0.25 0 0 0 0 0.004 0 1.1 0 0\n$TEXT:');return log,p

def test_work_parser_allows_absent_free_sentinel_and_never_returns_Rfree(tmp_path):
 log,p=xml(tmp_path,-999);r=parse_work_geometry(log,p);assert r['cycle']==30 and 'r_free' not in r

def test_work_parser_rejects_unmasked_test_data(tmp_path):
 log,p=xml(tmp_path,642)
 with pytest.raises(RuntimeError,match='observed test'):parse_work_geometry(log,p)


def test_failed_pilots_are_excluded_not_fatal():
 # A numerically divergent pilot is recorded as failed; healthy pilots still
 # drive selection, and all-failed falls back to the predeclared default.
 reports=[{'status':'failed','weight':.02,'reason':'timeout','error':'CGMAT diverged'},
          {'status':'completed','cycle':30,'weight':.005,'bond_rmsd_angstrom':.003,'angle_rmsd_degrees':.9},
          {'status':'completed','cycle':30,'weight':.01,'bond_rmsd_angstrom':.004,'angle_rmsd_degrees':1.0}]
 assert select_geometry_weight(reports)==.01
 assert select_geometry_weight([r for r in reports if r['status']=='failed'])==.005


def test_pilot_timeout_is_isolated_to_that_weight(tmp_path,monkeypatch):
 # The pilot in work_weight is a nested function; drive it through the real
 # calibrate_work_weight with a fake Popen that times out only for w0p02.
 import subprocess as sp
 path=tmp_path/'native.mtz'
 rows=np.array([[1,0,0,10,1,7],[0,1,0,20,2,3]],float)
 mtz(path,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],rows)
 model=tmp_path/'model.pdb'
 model.write_text('CRYST1   16.000   16.000   16.000  90.00  90.00  90.00 P 1           1\nEND\n')
 runtime=tmp_path/'runtime';(runtime/'bin').mkdir(parents=True)
 r=tmp_path/'calib'

 class FakeProc:
  def __init__(self,cmd,*a,**k):
   self.cmd=cmd;self.pid=99999;self.returncode=0
   self.cwd=Path(k['cwd'])
   self.weight='.02' if 'w0p02' in str(self.cwd) else None
  def communicate(self,input=None,timeout=None):
   if self.weight:raise sp.TimeoutExpired(self.cmd,timeout)
   xml(self.cwd,0)  # healthy pilot writes a valid 30-cycle REFMAC report
   return ('','')
  def wait(self,timeout=None):return 0
 def fake_popen(cmd,*a,**k):
  return FakeProc(cmd,*a,**k)
 monkeypatch.setattr(ww.subprocess,'Popen',fake_popen)
 monkeypatch.setattr(ww,'_terminate_group',lambda p:None)
 selected,report=ww.calibrate_work_weight(str(path),str(model),str(runtime),str(r),free_value=7)
 fail=[x for x in report['pilots'] if x['status']!='completed']
 assert len(fail)==1 and fail[0]['weight']==.02 and fail[0]['reason']=='timeout'
 assert len([x for x in report['pilots'] if x['status']=='completed'])==3
 assert (r/'FAILED_PILOTS.json').exists()
 assert selected==.04  # failed .02 excluded; largest healthy weight wins
