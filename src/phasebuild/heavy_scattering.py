"""Heavy-only anomalous factors with Cartesian anisotropic displacement tensors.

The isotropic route delegates unchanged to the existing direct summation.
No model-derived protein phases or isotropization are introduced.
"""
from dataclasses import dataclass
import gemmi
import numpy as np
from oasis_py.data import read_atoms,read_structure_checked,structure_factors

@dataclass
class AnisotropicHeavy:
    structure: object
    elements: list
    anisotropic_atoms: int


def read_heavy(path,data):
    st,_=read_structure_checked(path)
    atoms=[a for model in st for chain in model for residue in chain for a in residue if a.occ>0]
    if not any(a.aniso.nonzero() for a in atoms):return read_atoms(path,data)
    if len(st)!=1:raise ValueError('Use a single-model coordinate file')
    if st.cell.is_crystal() and not np.allclose(st.cell.parameters,data.cell.parameters,rtol=1e-4,atol=1e-3):raise ValueError('Model and reflection unit cells differ')
    if st.spacegroup_hm and gemmi.SpaceGroup(st.spacegroup_hm).operations()!=data.spacegroup.operations():raise ValueError('Model and reflection space-group settings differ')
    for a in atoms:
        if not a.element.atomic_number:raise ValueError('Unknown heavy-atom element')
        if not np.isfinite([a.pos.x,a.pos.y,a.pos.z,a.occ,a.b_iso]).all():raise ValueError('Nonfinite heavy-atom parameters')
        if a.aniso.nonzero():
            u=np.asarray(a.aniso.as_mat33())
            if not np.isfinite(u).all() or np.linalg.eigvalsh(u).min()<-1e-6:raise ValueError('Invalid anisotropic heavy-atom displacement tensor')
    if not atoms:raise ValueError('No occupied heavy atoms')
    return AnisotropicHeavy(st,[a.element.name for a in atoms],sum(a.aniso.nonzero() for a in atoms))


def anomalous_factors(data,heavy,fpp):
    if not isinstance(heavy,AnisotropicHeavy):return structure_factors(data,heavy,fpp)
    # q_cart is the reciprocal vector (cycles / Angstrom); U is Cartesian A^2.
    frac=np.asarray(data.cell.frac.mat);orth=np.asarray(data.cell.orth.mat)
    q=np.asarray(data.hkl,dtype=float)@frac
    result=np.zeros(len(data.f),complex)
    for chain in heavy.structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.occ<=0:continue
                fdouble=float(fpp[atom.element.name])
                if not np.isfinite(fdouble) or fdouble<=0:raise ValueError('Invalid anomalous scattering factor')
                p=data.cell.fractionalize(atom.pos);p=np.array([p.x,p.y,p.z]);seen=[]
                for op in data.spacegroup.operations():
                    position=np.mod(op.apply_to_xyz(p.tolist()),1.)
                    if any(np.max(abs((position-v+.5)%1-.5))<1e-7 for v in seen):continue
                    seen.append(position)
                    if atom.aniso.nonzero():
                        rotation=orth@(np.asarray(op.rot,dtype=float)/op.DEN)@frac
                        u=rotation@np.asarray(atom.aniso.as_mat33())@rotation.T
                        thermal=np.exp(-2*np.pi**2*np.einsum('ni,ij,nj->n',q,u,q))
                    else:thermal=np.exp(-data.stol2*atom.b_iso)
                    result+=1j*fdouble*atom.occ*thermal*np.exp(2j*np.pi*(data.hkl@position))
    return result


def scattering_metadata(heavy):
    return {'ADP_model':'Cartesian anisotropic U rotated with each symmetry operator; isotropic B retained for other sites' if isinstance(heavy,AnisotropicHeavy) else 'Original isotropic direct summation (unchanged)',
            'anisotropic_atoms':heavy.anisotropic_atoms if isinstance(heavy,AnisotropicHeavy) else 0,'anisotropization_or_isotropization_performed':False}
