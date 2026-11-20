import gemmi
import pytest
import importlib.util
import json
from pathlib import Path
from phasebuild.model_cleanup import cleanup_single_copy


def model():
    s = gemmi.Structure()
    s.cell = gemmi.UnitCell(20, 20, 20, 90, 90, 90)
    s.spacegroup_hm = 'P 1'
    m = gemmi.Model('1')
    for name, residues in [('A', ['ALA', 'MSE']), ('B', ['UNK', 'ALA', 'HOH', 'CA'])]:
        c = gemmi.Chain(name)
        for i, resname in enumerate(residues, 1):
            r = gemmi.Residue()
            r.name = resname
            r.seqid = gemmi.SeqId(i, ' ')
            a = gemmi.Atom()
            a.name = 'O' if resname == 'HOH' else 'CA'
            a.element = gemmi.Element('O' if resname == 'HOH' else 'Ca' if resname == 'CA' else 'C')
            a.pos = gemmi.Position(i, 2, 3)
            a.occ = .7
            a.b_iso = 31
            r.add_atom(a)
            c.add_residue(r)
        m.add_chain(c)
    s.add_model(m)
    return s


def test_cleanup_preserves_source_atoms_chemistry_and_nonprotein():
    s = model()
    before = s.make_mmcif_document().as_string()
    out, report = cleanup_single_copy(s, 'AM', asu_copies=1)
    assert s.make_mmcif_document().as_string() == before
    assert report['removed_protein_residues'] == 2
    assert [r.name for r in out[0]['A']] == ['ALA', 'MSE']
    assert [r.name for r in out[0]['B']] == ['HOH', 'CA']
    a = out[0]['A'][1][0]
    assert (a.pos.x, a.pos.y, a.pos.z) == (2, 2, 3)
    assert a.occ == pytest.approx(.7)
    assert a.b_iso == 31
    assert out.cell.parameters == s.cell.parameters
    assert out.spacegroup_hm == s.spacegroup_hm
    assert str(out[0]['A'][1].seqid) == str(s[0]['A'][1].seqid)


@pytest.mark.parametrize('copies', [0, 2, True, 1.0, None])
def test_requires_explicit_single_copy(copies):
    with pytest.raises(ValueError, match='single-copy'):
        cleanup_single_copy(model(), 'AM', asu_copies=copies)


@pytest.mark.parametrize('sequence', ['', 'AX', 'AMZ', '>fasta\nAM', None])
def test_ambiguous_sequence_rejected(sequence):
    with pytest.raises(ValueError):
        cleanup_single_copy(model(), sequence, asu_copies=1)


def test_two_complete_copies_are_ambiguous():
    s = model()
    c = s[0]['A'].clone()
    c.name = 'C'
    s[0].add_chain(c)
    with pytest.raises(ValueError, match='exactly one complete'):
        cleanup_single_copy(s, 'AM', asu_copies=1)


def test_partial_chain_and_missing_ca_rejected():
    with pytest.raises(ValueError, match='exactly one complete'):
        cleanup_single_copy(model(), 'AMA', asu_copies=1)
    s = model()
    del s[0]['A'][0][0]
    with pytest.raises(ValueError, match='exactly one complete'):
        cleanup_single_copy(s, 'AM', asu_copies=1)


def test_idempotent_and_no_extra_fragments():
    out, _ = cleanup_single_copy(model(), 'AM', asu_copies=1)
    again, report = cleanup_single_copy(out, 'AM', asu_copies=1)
    assert report['removed_protein_residues'] == 0
    assert again.make_mmcif_document().as_string() == out.make_mmcif_document().as_string()


def cli_main():
    path = Path(__file__).resolve().parents[1] / 'run_single_copy_cleanup.py'
    spec = importlib.util.spec_from_file_location('cleanup_cli_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.main


def test_cli_writes_audit_and_refuses_existing_output(tmp_path):
    source = tmp_path / 'source.cif'
    model().make_mmcif_document().write_file(str(source))
    before = source.read_bytes()
    sequence = tmp_path / 'seq.fasta'
    sequence.write_text('>single\nAM\n')
    output = tmp_path / 'new'
    args = ['--model', str(source), '--sequence', str(sequence), '--asu-copies', '1',
            '--output', str(output)]
    assert cli_main()(args) == 0
    audit = (output / 'cleanup_report.json').read_bytes()
    assert json.loads(audit)['removed_protein_residues'] == 2
    assert source.read_bytes() == before
    with pytest.raises(FileExistsError):
        cli_main()(args)
    assert (output / 'cleanup_report.json').read_bytes() == audit


def test_cli_rejects_multiple_sequences_before_writing(tmp_path):
    source = tmp_path / 'source.cif'
    model().make_mmcif_document().write_file(str(source))
    sequence = tmp_path / 'seq.fasta'
    sequence.write_text('>first\nAM\n>second\nAM\n')
    output = tmp_path / 'new'
    with pytest.raises(ValueError, match='one sequence record'):
        cli_main()(['--model', str(source), '--sequence', str(sequence),
                    '--asu-copies', '1', '--output', str(output)])
    assert not output.exists()
