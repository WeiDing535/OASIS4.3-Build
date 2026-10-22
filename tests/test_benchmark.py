import json
from types import SimpleNamespace
from pathlib import Path
import gemmi
import numpy as np
from phasebuild.benchmark import buccaneer_full
from phasebuild.geometry import template
from phasebuild.trace import write_backbone


def test_control_enables_sequence_and_sidechains_without_using_native_builder(tmp_path,monkeypatch):
    import phasebuild.benchmark as benchmark
    mtz=gemmi.Mtz(with_base=True);mtz.spacegroup=gemmi.find_spacegroup_by_name('P 1')
    mtz.set_cell_for_all(gemmi.UnitCell(40,40,40,90,90,90));mtz.add_dataset('test')
    for name,kind in [('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('FOM','W'),('FreeR_flag','I')]:
        mtz.add_column(name,kind)
    mtz.set_data(np.array([[1,1,0,10,1,10,20,.8,0]],dtype=np.float32));mtz.update_reso()
    data=tmp_path/'phase.mtz';mtz.write_to_file(str(data))
    seq=tmp_path/'sequence.fasta';seq.write_text('>protein\nACDEFGHIK\n')
    root=tmp_path/'ccp4';(root/'bin').mkdir(parents=True)
    (root/'bin/cbuccaneer').write_bytes(b'not run')
    (root/'bin/cbuccaneer').chmod(0o755)
    library=root/'lib/data/reference_structures';library.mkdir(parents=True)
    for name in ['reference-1tqw.pdb','reference-1tqw.mtz']:(library/name).write_bytes(b'training placeholder')
    def fake_run(command,**kwargs):
        assert '-sequence' in command and '-rebuild' in command and '-correlation-mode' in command
        assert command[command.index('-cycles')+1]=='3'
        if '-pdbin' in command:
            assert Path(command[command.index('-pdbin')+1])==tmp_path/'control/model.cif'
        sequence=Path(command[command.index('-seqin')+1]).read_text()
        assert sequence=='>protein\nACDEFGHIK\n'
        path=Path(kwargs['cwd'])/'temporary.cif'
        write_backbone(path,[template(-63.,-42.,length=2)+[20.,20.,20.]],
                       SimpleNamespace(cell=mtz.cell,sg=mtz.spacegroup))
        gemmi.read_structure(str(path)).write_pdb(str(Path(kwargs['cwd'])/'model.pdb'))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(benchmark.subprocess,'run',fake_run)
    report=buccaneer_full(data,seq,tmp_path/'control',ccp4=root,correlation=True)
    assert report['residues']==2 and report['status']=='success'
    json.dumps(report,allow_nan=False)
    seeded=buccaneer_full(data,seq,tmp_path/'seeded_control',ccp4=root,correlation=True,
                         starting_model=tmp_path/'control/model.cif')
    assert seeded['starting_model']['sha256']==report['output']['sha256']
