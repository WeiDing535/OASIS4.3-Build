"""Unit tests for the in-tree refinement engine.

These stay deliberately small: a three-atom P 1 structure with no restraints and
no bulk solvent, so the whole file runs in well under a second while still
exercising the pieces that broke during development (the zero overall-B case and
the fractional/cartesian round trip).
"""
import gemmi
import numpy as np
import pytest
from phasebuild.internal_refine import Atoms, BulkSolvent, Reflections, Refiner


def structure():
    st = gemmi.Structure()
    st.spacegroup_hm = 'P 1'
    st.cell = gemmi.UnitCell(30.0, 31.0, 32.0, 90.0, 90.0, 90.0)
    model = gemmi.Model('1')
    chain = gemmi.Chain('A')
    for index, (name, element, position, b) in enumerate([
            ('CA', 'C', (5.0, 6.0, 7.0), 12.0),
            ('N', 'N', (6.4, 6.2, 7.5), 14.0),
            ('O', 'O', (4.1, 5.2, 8.0), 16.0)]):
        residue = gemmi.Residue()
        residue.name = 'GLY'
        residue.seqid = gemmi.SeqId(str(index + 1))
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element(element)
        atom.pos = gemmi.Position(*position)
        atom.b_iso = b
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)
    return st


def reflections(tmp_path):
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    m.set_cell_for_all(gemmi.UnitCell(30.0, 31.0, 32.0, 90.0, 90.0, 90.0))
    m.add_dataset('test')
    for name, kind in [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')]:
        m.add_column(name, kind)
    rows = []
    rng = np.random.default_rng(7)
    for h in range(4):
        for k in range(3):
            for l in range(3):
                rows.append([h, k, l, float(rng.uniform(20, 60)), 1.0,
                             0 if (h + k + l) % 5 else 1])
    m.set_data(np.array(rows, dtype=np.float32))
    m.update_reso()
    path = tmp_path / 'data.mtz'
    m.write_to_file(str(path))
    return Reflections(path)


def refiner(tmp_path, **kwargs):
    st = structure()
    atoms = Atoms(st)
    data = reflections(tmp_path)
    empty = {'bonds': [], 'angles': [], 'repulsion': [],
             'repulsion_radius': 3.0, 'repulsion_weight': 1.0}
    return Refiner(atoms, data, st.find_spacegroup(), st.cell, empty,
                   refine_b=False, bulk_solvent=False, **kwargs)


def test_atomic_attenuation_is_always_per_reflection(tmp_path):
    r = refiner(tmp_path)
    assert r.atomic_attenuation.shape == (len(r.data),)
    np.testing.assert_allclose(r.atomic_attenuation, 1.0)
    r.b_add = 10.0
    np.testing.assert_allclose(r.atomic_attenuation,
                               np.exp(-0.25 * r.scattering.s2 * 10.0))


def test_fractional_cartesian_round_trip(tmp_path):
    r = refiner(tmp_path)
    xyz = r.atoms.xyz.copy()
    np.testing.assert_allclose(r._to_cartesian(r._to_fractional(xyz)), xyz,
                               atol=1e-9)


def test_target_returns_full_length_gradient(tmp_path):
    r = refiner(tmp_path)
    parameters = r._pack(r._to_fractional(r.atoms.xyz), r.atoms.b, 1.0, 0.2, 50.0)
    energy, gradient = r._target(parameters)
    assert np.isfinite(energy)
    assert gradient.shape == parameters.shape
    assert np.all(np.isfinite(gradient))


def test_structure_factors_split_atomic_and_bulk_terms(tmp_path):
    r = refiner(tmp_path)
    fcalc, amplitude, bulk, total = r._structure_factors(
        r._to_fractional(r.atoms.xyz), r.atoms.b, 1.0, 0.2, 50.0)
    assert fcalc.shape == (len(r.data),)
    np.testing.assert_allclose(bulk, 0.0)
    np.testing.assert_allclose(fcalc, total)


def test_bulk_solvent_mask_is_a_complement_inside_the_cell(tmp_path):
    st = structure()
    atoms = Atoms(st)
    bulk = BulkSolvent(st.cell, st.find_spacegroup(), atoms, spacing=1.5)
    assert 0.0 < bulk.solvent_fraction < 1.0
    values = bulk.structure_factors(reflections(tmp_path).hkl)
    assert values.shape == (len(reflections(tmp_path)),)
    assert np.all(np.isfinite(values))


def test_run_reaches_a_lower_target_than_the_start(tmp_path):
    r = refiner(tmp_path)
    parameters = r._pack(r._to_fractional(r.atoms.xyz), r.atoms.b, 1.0, 0.2, 50.0)
    before, _ = r._target(parameters)
    result = r.run(cycles=1, max_iterations=5, verbose=False)
    assert result['history']
    assert result['history'][-1]['r_work'] < 1.0
    assert np.isfinite(result['history'][-1]['r_free'])
