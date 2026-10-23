# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Crystallographic inputs, with explicit conventions and validation."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import warnings
import gemmi
import numpy as np
from scipy.optimize import least_squares


@dataclass
class Reflections:
    hkl: np.ndarray
    f: np.ndarray
    sigf: np.ndarray
    cell: gemmi.UnitCell
    spacegroup: gemmi.SpaceGroup
    dano: np.ndarray | None = None
    sigdano: np.ndarray | None = None
    fderivative: np.ndarray | None = None
    sigderivative: np.ndarray | None = None
    source: str = ''

    def __post_init__(self):
        self.hkl = np.asarray(self.hkl, dtype=np.int32)
        self.f = np.asarray(self.f, dtype=float)
        self.sigf = np.asarray(self.sigf, dtype=float)
        n = len(self.f)
        if self.hkl.shape != (n, 3) or self.sigf.shape != (n,):
            raise ValueError('Inconsistent reflection array shapes')
        if n == 0 or not np.all(np.isfinite(self.f)) or np.any(self.f <= 0):
            raise ValueError('Native phasing requires finite positive amplitudes')
        if not np.all(np.isfinite(self.sigf)) or np.any(self.sigf <= 0):
            raise ValueError('Finite positive amplitude uncertainties are required')
        if self.cell.volume <= 0 or not self.cell.is_crystal():
            raise ValueError('A valid crystal unit cell is required')
        if len(np.unique(self.hkl, axis=0)) != n:
            raise ValueError('Duplicate HKL rows; merge observations before phasing')
        if np.any(np.all(self.hkl == 0, axis=1)):
            raise ValueError('F(000) must not be included in measured reflections')
        for name in ['dano', 'sigdano', 'fderivative', 'sigderivative']:
            value = getattr(self, name)
            if value is not None:
                value = np.asarray(value, dtype=float)
                if value.shape != (n,) or not np.all(np.isfinite(value)):
                    raise ValueError(f'Invalid {name} array')
                if name.startswith('sig') and np.any(value <= 0):
                    raise ValueError(f'{name} must be positive')
                setattr(self, name, value)

    @property
    def stol2(self):
        return self.cell.calculate_1_d2_array(self.hkl) / 4

    @property
    def epsilon(self):
        return self.spacegroup.operations().epsilon_factor_array(self.hkl).astype(float)

    @property
    def centric(self):
        return self.spacegroup.operations().centric_flag_array(self.hkl)

    def subset(self, mask):
        values = {name: (None if getattr(self, name) is None else getattr(self, name)[mask])
                  for name in ['dano', 'sigdano', 'fderivative', 'sigderivative']}
        return Reflections(self.hkl[mask], self.f[mask], self.sigf[mask],
                           self.cell, self.spacegroup, source=self.source, **values)


def read_mtz(path, mode='dmr', labels=None, dmin=None, dmax=None,
             free_label=None, free_value=None):
    """Read amplitudes (not intensities). Missing required observations are removed.

    SAD can use FP/SIGFP/DANO/SIGDANO or a complete Friedel pair. SIR uses
    native FP and derivative FPH. Free-set exclusion is opt-in and explicit.
    """
    mode = mode.lower()
    if mode not in {'sad', 'sir', 'dmr'}:
        raise ValueError('mode must be sad, sir or dmr')
    mtz = gemmi.read_mtz_file(str(path))
    mtz.ensure_asu()
    labels = dict(labels or {})
    present = {c.label for c in mtz.columns}
    def column(key, candidates, types):
        name = labels.get(key)
        if name is None:
            name = next((c for c in candidates if c in present), None)
        if name is None or name not in present:
            raise ValueError(f'Missing {key}; specify --labels. Available: {sorted(present)}')
        c = mtz.column_with_label(name)
        if c.type not in types:
            raise ValueError(f'{name} has MTZ type {c.type}; expected {types} (amplitudes, not intensities)')
        return np.array(c.array, dtype=float)
    hkl = mtz.make_miller_array()
    extra = {}
    pair = mode == 'sad' and ('fplus' in labels or
            (not any(c in present for c in ['DANO','DANOM']) and 'F(+)' in present))
    if pair:
        fp = column('fplus',['F(+)'], 'FG')
        fm = column('fminus',['F(-)'], 'FG')
        sp = column('sigfplus',['SIGF(+)'], 'QL')
        sm = column('sigfminus',['SIGF(-)'], 'QL')
        f, sigf = (fp+fm)/2, np.hypot(sp, sm)/2
        extra = {'dano': fp-fm, 'sigdano': np.hypot(sp,sm)}
    else:
        f = column('f', ['FP','F','Fobs','FOBS'], 'FG')
        sigf = column('sigf', ['SIGFP','SIGF','SIGFobs','SIGFOBS'], 'QL')
        if mode == 'sad':
            extra = {'dano': column('dano',['DANO','DANOM'],'D'),
                     'sigdano': column('sigdano',['SIGDANO','SIGDANOM'],'Q')}
        elif mode == 'sir':
            extra = {'fderivative': column('fderivative',['FPH'],'FG'),
                     'sigderivative': column('sigderivative',['SIGFPH'],'QL')}
    keep = np.isfinite(f) & np.isfinite(sigf) & (f > 0) & (sigf > 0)
    for name,a in extra.items():
        keep &= np.isfinite(a)
        if name.startswith('sig') or name == 'fderivative':
            keep &= a > 0
    ops = mtz.spacegroup.operations()
    keep &= ~ops.systematic_absences(hkl)
    keep &= np.any(hkl != 0, axis=1)
    d = 1/np.sqrt(mtz.cell.calculate_1_d2_array(hkl))
    if dmin is not None: keep &= d >= dmin
    if dmax is not None: keep &= d <= dmax
    if free_label is not None:
        if free_value is None:
            raise ValueError('Specify free_value explicitly; conventions differ across datasets')
        if free_label not in present:
            raise ValueError(f'Missing free-set column {free_label}')
        flag = np.array(mtz.column_with_label(free_label).array)
        keep &= np.isfinite(flag) & (flag != free_value)
    result = Reflections(hkl[keep],f[keep],sigf[keep],mtz.cell,mtz.spacegroup,
                         source=str(Path(path).resolve()),
                         **{name:a[keep] for name,a in extra.items()})
    return result, {'input_rows':len(hkl), 'retained_rows':int(keep.sum()),
                    'removed_rows':int((~keep).sum()), 'mode':mode,
                    'free_set_excluded':free_label is not None}


def parse_composition(text=None, residues=None):
    if text:
        result = {}
        for field in text.split(','):
            element,count = field.split(':')
            el = gemmi.Element(element.strip())
            if el.atomic_number == 0 or el.it92 is None:
                raise ValueError(f'Unsupported element: {element}')
            number=float(count)
            if not np.isfinite(number) or number <= 0:
                raise ValueError('Composition counts must be finite and positive')
            result[el.name] = result.get(el.name,0)+number
        return result
    if residues is None or residues <= 0:
        raise ValueError('Provide ASU composition or a positive residue count')
    return {'C':5*residues,'N':1.2*residues,'O':1.5*residues,'H':8*residues}


def form_factor(element, stol2):
    c = gemmi.Element(element).it92
    if c is None:
        raise ValueError(f'No X-ray scattering table for {element}')
    return np.exp(-np.asarray(stol2)[...,None]*np.asarray(c.b)) @ np.asarray(c.a)+c.c


def normalize(data, composition, bins=20):
    """Robust isotropic Wilson fit; composition counts are per ASU.

    E=F/sqrt(epsilon*sum_cell(n*f^2)*scale*exp(-2B*s^2)).
    This intentionally does not reproduce legacy index-group rescaling.
    """
    s2=data.stol2
    order=len(data.spacegroup.operations())
    scatter=sum(count*form_factor(el,s2)**2 for el,count in composition.items())*order
    reduced=data.f**2/(data.epsilon*scatter)
    groups=np.array_split(np.argsort(s2),min(bins,max(2,len(s2)//30)))
    x=np.array([s2[g].mean() for g in groups if len(g)])
    y=np.log([reduced[g].mean() for g in groups if len(g)])
    if len(x)<2 or np.ptp(x)<1e-12:
        scale=float(np.mean(reduced)); b=0.
    else:
        fit=least_squares(lambda p:p[0]-2*p[1]*x-y,[float(np.mean(y)),10.],
                          loss='soft_l1', bounds=([-100,-100],[100,500]))
        scale,b=float(np.exp(fit.x[0])),float(fit.x[1])
    expected=data.epsilon*scatter*scale*np.exp(-2*b*s2)
    if np.any(expected<=0) or not np.all(np.isfinite(expected)):
        raise ValueError('Wilson model produced invalid expected intensities')
    e=data.f/np.sqrt(expected)
    # Number of atoms in the entire unit cell, same convention as Cochran.
    sigma2=sum(n*gemmi.Element(el).atomic_number**2 for el,n in composition.items())*order
    sigma3=sum(n*gemmi.Element(el).atomic_number**3 for el,n in composition.items())*order
    return e,expected,2*sigma3/sigma2**1.5, {'scale':scale,'wilson_b':b,
                                          'mean_e2':float(np.mean(e**2))}


@dataclass
class Atoms:
    fractional: np.ndarray
    elements: list[str]
    occupancy: np.ndarray
    b_iso: np.ndarray

    def __post_init__(self):
        self.fractional=np.asarray(self.fractional,dtype=float)
        self.occupancy=np.asarray(self.occupancy,dtype=float)
        self.b_iso=np.asarray(self.b_iso,dtype=float)
        n=len(self.elements)
        if self.fractional.shape!=(n,3) or self.occupancy.shape!=(n,) or self.b_iso.shape!=(n,):
            raise ValueError('Inconsistent atom array shapes')
        if n==0 or not all(np.all(np.isfinite(a)) for a in [self.fractional,self.occupancy,self.b_iso]):
            raise ValueError('An empty or non-finite atom model is not supported')
        if np.any(self.occupancy<0) or np.any(self.b_iso<0):
            raise ValueError('Negative occupancy or isotropic B factor')
        self.fractional%=1

    def subset(self,index):
        index=np.atleast_1d(index)
        return Atoms(self.fractional[index],[self.elements[i] for i in index],
                     self.occupancy[index],self.b_iso[index])


def read_structure_checked(path):
    """Read coordinates, repairing a narrowly identifiable PDB column error.

    Some exported PDB files put a two-letter element in columns 78–79
    instead of 77–78 (Se becomes S to a standard reader). Repair only when
    the complete valid symbol also agrees with the atom name, the preceding
    column is blank, and no charge/trailing text makes the shift ambiguous.
    The original file is never modified; every repair is reported.
    """
    path=Path(path)
    repairs=[]
    if path.suffix.lower() not in {'.pdb','.ent'}:
        return gemmi.read_structure(str(path)),repairs
    lines=path.read_text().splitlines(keepends=True)
    for i,line in enumerate(lines):
        if not line.startswith(('ATOM  ','HETATM')):continue
        old=line[76:78].strip();shifted=line[77:79].strip()
        name=line[12:16].strip()
        if (len(shifted)==2 and line[76:77]==' ' and not line[79:].strip()
                and gemmi.Element(shifted).atomic_number>0
                and gemmi.Element(shifted).name.upper()==name.upper()
                and gemmi.Element(shifted).name!=gemmi.Element(old).name):
            symbol=gemmi.Element(shifted).name
            lines[i]=line[:76]+symbol.rjust(2)+line[79:]
            repairs.append({'line':i+1,'old_element_field':old,'element':symbol,
                            'operation':'repair one-column-right PDB element field',
                            'evidence':'complete element symbol and matching atom name'})
    if repairs:
        warnings.warn(f'{path}: repaired {len(repairs)} shifted PDB element fields; '
                      f'first repair at line {repairs[0]["line"]} '
                      f'{repairs[0]["old_element_field"]} -> {repairs[0]["element"]}',
                      UserWarning,stacklevel=2)
        return gemmi.read_pdb_string(''.join(lines)),repairs
    return gemmi.read_structure(str(path)),repairs


def read_atoms(path,data):
    st,_=read_structure_checked(path)
    if len(st)!=1: raise ValueError('Use a single-model coordinate file')
    if st.cell.is_crystal() and not np.allclose(st.cell.parameters,data.cell.parameters,rtol=1e-4,atol=1e-3):
        raise ValueError('Model and reflection unit cells differ; reindex/transform first')
    if st.spacegroup_hm and gemmi.SpaceGroup(st.spacegroup_hm).operations()!=data.spacegroup.operations():
        raise ValueError('Model and reflection space-group settings differ')
    xyz=[];elements=[];occ=[];b=[]
    for chain in st[0]:
        for residue in chain:
            for atom in residue:
                if atom.occ<=0: continue
                if atom.element.atomic_number==0:
                    raise ValueError(f'Unknown element at {chain.name}/{residue.seqid}/{atom.name}')
                if atom.aniso.nonzero():
                    raise ValueError('Native prototype supports isotropic ADPs; use reference backend for anisotropic models')
                p=data.cell.fractionalize(atom.pos)
                xyz.append([p.x,p.y,p.z]);elements.append(atom.element.name)
                occ.append(atom.occ);b.append(atom.b_iso)
    return Atoms(xyz,elements,occ,b)


def structure_factors(data,atoms,anomalous=None,chunk=128):
    """Vectorized direct summation, F=sum f exp(+2pi i h.r), isotropic ADPs.

    anomalous=None: normal X-ray form factors. Otherwise the values are f'',
    and the result includes the factor i. Symmetry special positions are
    deduplicated per input atom; atoms must be an ASU model.
    """
    xyz=[];el=[];oc=[];bi=[]
    for p,e,o,b in zip(atoms.fractional,atoms.elements,atoms.occupancy,atoms.b_iso):
        orbit=[]
        for op in data.spacegroup.operations():
            q=np.mod(op.apply_to_xyz(p.tolist()),1.)
            if not any(np.max(np.abs((q-r+.5)%1-.5))<1e-7 for r in orbit):
                orbit.append(q)
        xyz.extend(orbit);el.extend([e]*len(orbit));oc.extend([o]*len(orbit));bi.extend([b]*len(orbit))
    xyz=np.asarray(xyz);oc=np.asarray(oc);bi=np.asarray(bi)
    elements=sorted(set(el));ids={e:np.flatnonzero(np.array(el)==e) for e in elements}
    out=np.zeros(len(data.f),complex)
    all_stol2=data.stol2
    for start in range(0,len(data.f),chunk):
        sl=slice(start,start+chunk);s2=all_stol2[sl]
        # Bound temporaries in both reflection and atom dimensions.
        for element,index in ids.items():
            scattering=(form_factor(element,s2) if anomalous is None
                        else np.full(len(s2),1j*float(anomalous.get(element,0)),complex))
            for offset in range(0,len(index),1024):
                ix=index[offset:offset+1024]
                phases=np.exp(2j*np.pi*(data.hkl[sl]@xyz[ix].T))
                thermal=np.exp(-s2[:,None]*bi[ix])
                out[sl]+=scattering*(np.sum(phases*thermal*oc[ix],axis=1))
    return out


def write_mtz(path,data,phase,fom,hl=None,probability=None):
    phase=np.asarray(phase);fom=np.asarray(fom)
    if phase.shape!=data.f.shape or fom.shape!=data.f.shape:
        raise ValueError('Output shapes differ from reflection data')
    if not np.all(np.isfinite(phase)) or not np.all(np.isfinite(fom)) or np.any((fom<0)|(fom>1)):
        raise ValueError('Invalid output phases or FOM')
    mtz=gemmi.Mtz(with_base=True)
    mtz.spacegroup=data.spacegroup;mtz.set_cell_for_all(data.cell);mtz.add_dataset('OASIS_PY_RESEARCH')
    cols=[('FP','F'),('SIGFP','Q'),('PHIB','P'),('FOM','W')]
    arrays=[data.hkl,data.f[:,None],data.sigf[:,None],np.rad2deg(phase[:,None])%360,fom[:,None]]
    if probability is not None:
        p=np.asarray(probability)
        if p.shape!=data.f.shape or not np.all(np.isfinite(p)) or np.any((p<0)|(p>1)):
            raise ValueError('Invalid PPLUS')
        cols.append(('PPLUS','R'));arrays.append(p[:,None])
    if hl is not None:
        hl=np.asarray(hl)
        if hl.shape!=(len(data.f),4) or not np.all(np.isfinite(hl)):
            raise ValueError('Invalid HL coefficients')
        cols.extend((name,'A') for name in ['HLA','HLB','HLC','HLD']);arrays.append(hl)
    for name,typ in cols:mtz.add_column(name,typ)
    mtz.set_data(np.column_stack(arrays).astype(np.float32))
    mtz.sort();mtz.update_reso()
    mtz.history=['OASIS-PY research prototype; not legacy-validated phasing',
                 'Original input flags/columns are not copied; retain source MTZ']
    mtz.write_to_file(str(path))
