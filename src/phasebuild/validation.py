"""Map-independent geometric checks, with explicit conservative thresholds."""
import gemmi
import numpy as np


def angle(a, b, c):
    x, y = a-b, c-b
    return float(np.rad2deg(np.arccos(np.clip(np.dot(x,y) / (np.linalg.norm(x)*np.linalg.norm(y)), -1.,1.))))


def peptide_metrics(left, right):
    torsion = np.rad2deg(gemmi.calculate_dihedral(*[gemmi.Position(*p)
                           for p in [left[1],left[2],right[0],right[1]]]))
    return {'c_n': float(np.linalg.norm(left[2]-right[0])),
            'ca_ca': float(np.linalg.norm(left[1]-right[1])),
            'ca_c_n': angle(left[1],left[2],right[0]),
            'c_n_ca': angle(left[2],right[0],right[1]),
            'omega_trans_deviation': abs(180.-abs(float(torsion)))}


def peptide_valid(metrics):
    return (1.15 <= metrics['c_n'] <= 1.50 and 3.3 <= metrics['ca_ca'] <= 4.2 and
            abs(metrics['ca_c_n']-116.2) <= 15. and abs(metrics['c_n_ca']-121.7) <= 15. and
            metrics['omega_trans_deviation'] <= 25.)


def geometry_report(chains, density=None):
    bonds, links, issues = [], [], []
    ideals = [(0,1,1.458,'N-CA'),(1,2,1.525,'CA-C'),(2,3,1.231,'C-O')]
    for ci, chain in enumerate(chains):
        for ri, residue in enumerate(chain):
            for a,b,ideal,name in ideals:
                length = float(np.linalg.norm(residue[a]-residue[b]))
                bonds.append(length-ideal)
                if abs(length-ideal) > .12:
                    issues.append({'fragment':ci+1,'residue':ri+1,'kind':name,'length':length})
        for ri in range(len(chain)-1):
            metrics = peptide_metrics(chain[ri],chain[ri+1])
            links.append(metrics)
            if not peptide_valid(metrics):
                issues.append({'fragment':ci+1,'residue':ri+1,'kind':'peptide_link',**metrics})
    clashes = []
    if density is not None:
        ca = [(ci,ri,r[1]) for ci,c in enumerate(chains) for ri,r in enumerate(c)]
        for i,(ci,ri,a) in enumerate(ca):
            for cj,rj,b in ca[:i]:
                if ci == cj and abs(ri-rj) <= 1:
                    continue
                distance = density.distance(a,b)
                if distance < 2.8:
                    clashes.append({'fragment_1':ci+1,'residue_1':ri+1,'fragment_2':cj+1,
                                    'residue_2':rj+1,'ca_distance':distance})
    return {'residues':sum(map(len,chains)),'fragments':len(chains),'peptide_links':len(links),
            'bond_length_rms_deviation':float(np.sqrt(np.mean(np.square(bonds)))) if bonds else None,
            'max_omega_trans_deviation':max((m['omega_trans_deviation'] for m in links),default=None),
            'severe_geometry_issues':issues,'nonadjacent_ca_clashes':clashes,
            'scope':'Ideal bond/peptide checks and symmetry-aware CA clashes. '
                    'Not a complete all-atom validation; cis peptides are not supported.'}
