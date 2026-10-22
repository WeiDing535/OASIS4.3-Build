import gemmi
import pytest
from phasebuild.evaluate import evaluate


def write_structure(path, positions, spacegroup='I 21 3'):
    s = gemmi.Structure()
    s.cell = gemmi.UnitCell(40, 40, 40, 90, 90, 90)
    s.spacegroup_hm = spacegroup
    m, c = gemmi.Model('1'), gemmi.Chain('A')
    for i, pos in enumerate(positions):
        r, a = gemmi.Residue(), gemmi.Atom()
        r.name, r.seqid = 'ALA', gemmi.SeqId(i + 1, ' ')
        a.name, a.element, a.pos = 'CA', gemmi.Element('C'), gemmi.Position(*pos)
        r.add_atom(a)
        c.add_residue(r)
    m.add_chain(c)
    s.add_model(m)
    s.make_mmcif_document().write_file(str(path))


def test_symmetry_translation_and_one_to_one_matching(tmp_path):
    ref, built = tmp_path / 'ref.cif', tmp_path / 'built.cif'
    points = [(3., 6., 9.), (5., 8., 11.)]
    write_structure(ref, points)
    cell = gemmi.UnitCell(40, 40, 40, 90, 90, 90)
    op = list(gemmi.find_spacegroup_by_name('I 21 3').operations())[3]
    transformed = []
    for point in points:
        f = cell.fractionalize(gemmi.Position(*point))
        p = cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(list(f))))
        transformed.append((p.x + 40, p.y, p.z))
    # A duplicate must not match the same reference residue twice.
    write_structure(built, transformed + [transformed[0]])
    result = evaluate(built, ref)
    assert result['matched_residues'] == 2
    assert result['ca_completeness'] == 1.
    assert result['ca_precision'] == pytest.approx(2 / 3)
    assert result['matched_ca_rmsd'] < 1e-5


def test_reject_incompatible_cell(tmp_path):
    ref, built = tmp_path / 'ref.cif', tmp_path / 'built.cif'
    write_structure(ref, [(3., 6., 9.)])
    write_structure(built, [(3., 6., 9.)])
    s = gemmi.read_structure(str(built))
    s.cell = gemmi.UnitCell(41, 40, 40, 90, 90, 90)
    s.make_mmcif_document().write_file(str(built))
    with pytest.raises(ValueError, match='unit cells differ'):
        evaluate(built, ref)


def test_ca_match_does_not_imply_correct_backbone_direction(tmp_path):
    ref, built = tmp_path / 'ref.cif', tmp_path / 'built.cif'
    write_structure(ref, [(3., 6., 9.)])
    s = gemmi.read_structure(str(ref))
    residue = s[0][0][0]
    for name, offset in [('N', -1.), ('C', 1.)]:
        a = gemmi.Atom()
        a.name, a.element = name, gemmi.Element('N' if name == 'N' else 'C')
        a.pos = gemmi.Position(3. + offset, 6., 9.)
        residue.add_atom(a)
    s.make_mmcif_document().write_file(str(ref))
    n, c = residue.find_atom('N', '*'), residue.find_atom('C', '*')
    n.pos, c.pos = gemmi.Position(4., 6., 9.), gemmi.Position(2., 6., 9.)
    s.make_mmcif_document().write_file(str(built))
    result = evaluate(built, ref)
    assert result['matched_residues'] == 1
    assert result['n_ca_c_correct_residues'] == 0
    assert result['sequence_evaluation']['correct_type_with_ca_match'] == 1
    assert result['sequence_evaluation']['correct_type_and_backbone'] == 0
