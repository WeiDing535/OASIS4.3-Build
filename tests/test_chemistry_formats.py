import pytest
from phasebuild.chemistry import monomer,close_pairs,ring_groups


@pytest.mark.parametrize('order',['type','value_order'])
def test_monomer_accepts_both_ccp4_bond_order_column_names(tmp_path,order):
    folder=tmp_path/'p';folder.mkdir()
    (folder/'PRO.cif').write_text('''data_comp_PRO
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
N N
CA C
C C
O O
CB C
CG C
CD C
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.'''+order+'''
_chem_comp_bond.value_dist
N CA SINGLE 1.496
CA C SINGLE 1.525
C O DOUBLE 1.252
CA CB SINGLE 1.526
CB CG SINGLE 1.516
CG CD SINGLE 1.497
CD N SINGLE 1.499
''')
    comp=monomer('PRO',tmp_path)
    assert len(comp['bonds'])==7 and ('N','CA') in close_pairs(comp,list(comp['elements']))
    assert set(ring_groups(comp)[0])=={'CB','CG','CD'}
