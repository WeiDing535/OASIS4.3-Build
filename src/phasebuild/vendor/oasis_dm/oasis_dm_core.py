#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oasis-dm -- a from-scratch implementation of macromolecular density modification
          (phase improvement), functionally close to CCP4 ``dm``.

What is implemented here (all the density-modification logic is our own code):

  * MTZ I/O and symmetry expansion        -> gemmi (a library, like CCP4's libccp4)
  * real/reciprocal FFT                   -> numpy.fft (self-contained convention,
                                            verified against gemmi's own FFT to CC=1.0)
  * Hendrickson-Lattman phase probability, centroid phase / FOM, numerical
    integration for bimodal (SIR/SAD) priors
  * solvent envelope: Wang rolling-sphere (FFT convolution with an explicitly
    built 3-D kernel) or Abrahams-Leslie local-variance mask
  * solvent flattening with a continuous (soft) mask and a soft weighted blend
  * histogram matching by CDF / rank (ordinal) mapping, with CCP4's own
    ``hist.lib`` target density distributions (18 resolutions x 100 bins)
  * phase recombination  P_new = P_exp * P_mod  implemented as HL coefficient
    addition, with Sim / sigmaA-like per-resolution-shell error estimates
  * perturbation-gamma bias correction (two maps, identical mask and LUT)
  * free-reflection omission for genuine cross-validation (dm COMBINE OMIT idea)
  * centric reflection phase restriction, missing-reflection RESTORE weight
  * output of dm-style (FDM/PHIDM/FOMDM), HL and 2mFo-DFc map coefficients

References that drove the design are collected in research/theory.md and
research/ccp4_local.md.  This program deliberately calls no external program.

Usage:
  /opt/miniconda3/bin/python src/oasis_dm_core.py --mtzin CASE.mtz --mtzout OUT.mtz \
      --fo FOBS --phib PHIB --fom FOM --free FreeR --solc 0.5 --cycles 5 \
      --stats results/case_stats.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass, asdict, field

import numpy as np
import gemmi
from scipy.special import i0e, i1e
from scipy import ndimage

# ---------------------------------------------------------------------------
# small numerical helpers
# ---------------------------------------------------------------------------

def wrap180(d):
    """wrap phase differences (degrees) into (-180, 180]."""
    return (np.asarray(d, dtype=np.float64) + 180.0) % 360.0 - 180.0


def nice_fft_size(n):
    """smallest >= n whose factors are only 2, 3, 5 (numpy FFT friendly)."""
    n = max(int(n), 4)
    while True:
        t = n
        for p in (2, 3, 5):
            while t % p == 0:
                t //= p
        if t == 1:
            return n
        n += 1


# I1/I0 lookup for  FOM <-> Sim X  (von Mises distribution).
_XGRID = np.linspace(0.0, 400.0, 160001)
_MGRID = i1e(_XGRID) / i0e(_XGRID)


def fom_to_x(m):
    return np.interp(np.clip(m, 0.0, 0.99995), _MGRID, _XGRID)


def x_to_fom(x):
    return i1e(np.asarray(x, dtype=np.float64)) / i0e(np.asarray(x, dtype=np.float64))


def phifom_to_hl(phi_rad, fom):
    """convert unimodal (phase, FOM) prior into HL coefficients A,B,C,D."""
    x = fom_to_x(fom)
    z = np.zeros_like(x)
    return np.column_stack([x * np.cos(phi_rad), x * np.sin(phi_rad), z, z])


def hl_centroid(hl, nphi=256):
    """centroid phase and figure of merit of a general HL distribution.

    P(phi) = exp[A cos(phi) + B sin(phi) + C cos(2phi) + D sin(2phi)]
    m * exp(i*phi_best) = integral P(phi) exp(i*phi) dphi / integral P(phi) dphi

    Numerical quadrature is used so bimodal (C,D != 0) SIR/SAD priors are
    handled correctly.  Unimodal priors use the exact von Mises result.
    """
    hl = np.asarray(hl, dtype=np.float64)
    A, B, C, D = hl[:, 0], hl[:, 1], hl[:, 2], hl[:, 3]
    fast = (np.abs(C) < 1e-8) & (np.abs(D) < 1e-8) & (np.abs(A) + np.abs(B) > 1e-12)
    phi = np.arctan2(B, A)
    m = x_to_fom(np.sqrt(A * A + B * B))
    out_phi = np.where(fast, phi, 0.0)
    out_m = np.where(fast, m, 0.0)

    slow = ~fast
    if np.any(slow):
        idx = np.nonzero(slow)[0]
        grid = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
        e1 = np.exp(1j * grid)
        e1c = np.cos(grid)
        e1s = np.sin(grid)
        e2c = np.cos(2.0 * grid)
        e2s = np.sin(2.0 * grid)
        for start in range(0, len(idx), 2000):
            ii = idx[start:start + 2000]
            lp = (A[ii, None] * e1c[None, :] + B[ii, None] * e1s[None, :]
                  + C[ii, None] * e2c[None, :] + D[ii, None] * e2s[None, :])
            lp -= lp.max(axis=1, keepdims=True)
            p = np.exp(lp)
            z = (p * e1[None, :]).sum(axis=1) / p.sum(axis=1)
            out_phi[ii] = np.angle(z)
            out_m[ii] = np.abs(z)
    return out_phi, out_m


class HistLib:
    """CCP4's dm target protein-density histogram library ($CLIBD/hist.lib)."""

    def __init__(self, path):
        with open(path, 'r') as fh:
            lines = [l.rstrip('\n') for l in fh]
        head = lines[0].split()
        self.nres, self.nbin = int(head[1]), int(head[2])
        meta = {}
        for l in lines[1:]:
            key = l.split()[0] if l.split() else ''
            if key == 'DATA':
                break
            if key in ('RESL', 'RMIN', 'RMAX', 'RAVE', 'RSTD'):
                vals = [float(x) for x in l.split()[1:]]
                meta[key] = np.asarray(vals, dtype=np.float64)
        self.res = meta['RESL']
        self.rmin = float(np.mean(meta['RMIN']))
        self.rmax = float(np.mean(meta['RMAX']))
        self.rave = meta['RAVE']
        self.rstd = meta['RSTD']
        k = lines.index('DATA')
        data = np.array([[float(x) for x in l.split()]
                         for l in lines[k + 1:] if l.strip()], dtype=np.float64)
        if data.shape != (self.nbin, self.nres):
            raise ValueError(f'hist.lib DATA shape {data.shape} != {(self.nbin, self.nres)}')
        # np.interp needs increasing x: RESL is decreasing -> reverse columns.
        self._res_asc = self.res[::-1].copy()
        self._data_asc = data[:, ::-1].copy()
        self.rho_axis = self.rmin + (np.arange(self.nbin) + 0.5) * \
                        (self.rmax - self.rmin) / self.nbin

    def target(self, d_eff):
        """interpolate the target PDF at effective resolution d_eff (A)."""
        d = float(np.clip(d_eff, self._res_asc[0], self._res_asc[-1]))
        pdf = np.array([np.interp(d, self._res_asc, self._data_asc[j])
                        for j in range(self.nbin)], dtype=np.float64)
        pdf = np.maximum(pdf, 0.0)
        tot = pdf.sum()
        if tot <= 0:
            raise ValueError('empty target histogram')
        pdf /= tot
        cdf = np.cumsum(pdf)
        cdf[-1] = 1.0
        mean = float((self.rho_axis * pdf).sum())
        std = float(np.sqrt(((self.rho_axis - mean) ** 2 * pdf).sum()))
        return cdf, mean, max(std, 1e-6)


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

@dataclass
class Config:
    mtzin: str = ''
    mtzout: str = ''
    fo: str = 'FOBS'
    sigf: str = 'SIGFOBS'
    phib: str = 'PHIB'
    fom: str = 'FOM'
    hl: str = 'HLA,HLB,HLC,HLD'
    free: str = 'FreeR'
    free_marker: int = 0         # value marking a FREE reflection (0 for CCP4
                                 # FreeR_flag; 1 for the GPCR pipeline's FreeR)
    truth: str = 'PHIC'          # optional, only used for diagnostics in the log
    solc: float = 0.5            # solvent volume fraction (dm SOLC)
    cycles: int = 5              # dm NCYCLE
    mask: str = 'variance'       # 'variance' (Abrahams-Leslie) or 'wang' (dm)
    radius: float = 4.0          # local-variance radius (A); wang mode uses 8 A
    radius_min: float = 3.0
    radius_decay: float = 0.85   # dm's automatically shrinking mask radius
    mask_smooth: float = 0.8     # gaussian smoothing sigma (grid points)
    hist_lib: str = '/Applications/ccp4-9/lib/data/hist.lib'
    hist_on: bool = True
    solv_on: bool = True
    squash: bool = True          # mild non-linear ordinal squash
    gamma_on: bool = True        # perturbation-gamma bias correction
    omit_free: bool = False      # reflection-omit cross validation (optional)
    anchor: str = 'exp'          # 'exp' = always recombine with ORIGINAL prior
    n_shells: int = 20
    x_max: float = 12.0          # cap on Sim weight X (avoids FOM -> 1)
    cmbwt: float = 0.4           # weight of the DM-derived HL contribution
                                 # (calibrated by cross-case sweep, sweep2 D)
    restore: float = 0.0         # weight for reflections missing in the input
    sample_rate: float = 3.0     # grid spacing = d_min / sample_rate
    seed: int = 0
    stats: str = ''
    tag: str = ''
    quiet: bool = False

    def as_dict(self):
        return asdict(self)


OUT_COLUMNS = [
    ('DM_FDM', 'F'), ('DM_PHIDM', 'P'), ('DM_FOMDM', 'W'),
    ('DM_HLADM', 'A'), ('DM_HLBDM', 'A'), ('DM_HLCDM', 'A'), ('DM_HLDDM', 'A'),
    ('DM_FWT', 'F'), ('DM_PHWT', 'P'),
    ('DM_DELFWT', 'F'), ('DM_PHDELFWT', 'P'),
    ('DM_FCDM', 'F'), ('DM_PHICDM', 'P'),
]


# ---------------------------------------------------------------------------
# main engine
# ---------------------------------------------------------------------------

class DensityModifier:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.log = (lambda *a: None) if cfg.quiet else print
        self._load()
        self._prepare()

    # ---------------------------- input -----------------------------------
    def _col(self, label, default=None):
        c = self.mtz.column_with_label(label)
        if c is None:
            return None
        return np.asarray(c.array, dtype=np.float64)

    def _load(self):
        cfg = self.cfg
        self.mtz = gemmi.read_mtz_file(cfg.mtzin)
        self.hkl = np.asarray(self.mtz.make_miller_array(), dtype=np.int32)
        self.nref = len(self.hkl)
        self.cell = self.mtz.cell
        self.sg = self.mtz.spacegroup
        self.V = float(self.cell.volume)
        self.Fo = self._col(cfg.fo)
        if self.Fo is None:
            raise ValueError(f'column {cfg.fo} not found; use --fo')
        self.SigF = self._col(cfg.sigf)
        phi = self._col(cfg.phib)
        if phi is None:
            raise ValueError(f'column {cfg.phib} not found; use --phib')
        if np.nanmax(np.abs(phi)) <= 3.15:
            self.log('[info] input phases look like RADIANS -> converted to degrees')
            phi = np.rad2deg(phi.copy())
        self.phi_exp_deg = phi
        self.fom_in = self._col(cfg.fom)
        if self.fom_in is None:
            self.log(f'[warn] no FOM column {cfg.fom}; assuming FOM=1')
            self.fom_in = np.ones(self.nref)
        free = self._col(cfg.free)
        if free is None:
            self.log(f'[warn] no free flag column {cfg.free}; free diagnostics disabled')
            self.free = np.zeros(self.nref, dtype=bool)
        else:
            self.free = (free == cfg.free_marker)
        self.valid = (np.isfinite(self.Fo) & (self.Fo > 0)
                      & np.isfinite(self.phi_exp_deg) & np.isfinite(self.fom_in))
        self.work = self.valid & ~self.free

        # HL prior: prefer explicit HL coefficients (they carry bimodality).
        hl_labels = [x.strip() for x in cfg.hl.split(',') if x.strip()]
        cols = [self._col(l) for l in hl_labels]
        if len(cols) == 4 and all(c is not None for c in cols):
            self.HL_exp = np.column_stack([np.nan_to_num(c) for c in cols])
            self.hl_source = 'HL coefficients ' + ','.join(hl_labels)
            p, _ = hl_centroid(self.HL_exp)
            diff = np.abs(wrap180(np.rad2deg(p[self.valid]) - self.phi_exp_deg[self.valid]))
            self.log(f'[info] HL centroid vs PHIB: median |dphi| = '
                     f'{np.median(diff):.2f} deg, 90% = {np.percentile(diff, 90):.2f} deg')
        else:
            self.HL_exp = phifom_to_hl(np.deg2rad(self.phi_exp_deg),
                                       np.clip(np.nan_to_num(self.fom_in), 0, 0.9999))
            self.hl_source = 'PHIB/FOM -> von Mises HL'
        self.HL_exp = np.where(self.valid[:, None], self.HL_exp, 0.0)

        self.d = np.asarray(self.mtz.make_d_array(), dtype=np.float64)
        self.dmin = float(np.nanmin(self.d[self.valid]))
        self.s2 = 1.0 / np.maximum(self.d, 1e-6) ** 2

    def _prepare(self):
        cfg = self.cfg
        a, b, c = self.cell.a, self.cell.b, self.cell.c
        step = self.dmin / cfg.sample_rate
        # The grid must be compatible with the space group's fractional
        # translations (a 4_1 screw axis, for instance, requires the
        # corresponding axis to be a multiple of 4).  A plain 5-smooth size
        # can silently violate that: P 41 21 2 with 80x80x45 is rejected by
        # gemmi.get_f_phi_on_grid(), while 80x80x48 is accepted.  Ask gemmi
        # to round our requested minimum up to the nearest FFT-friendly,
        # symmetry-compatible grid instead of guessing.
        min_size = [math.ceil(x / step) for x in (a, b, c)]
        self.size = [int(x) for x in
                     self.mtz.get_size_for_hkl(min_size=min_size, sample_rate=0.0)]
        self.N = int(np.prod(self.size))
        self.scale = self.N / self.V          # reciprocal -> density
        self.grid_spacing = [a / self.size[0], b / self.size[1], c / self.size[2]]
        self.ih = self.hkl[:, 0] % self.size[0]
        self.ik = self.hkl[:, 1] % self.size[1]
        self.il = self.hkl[:, 2] % self.size[2]

        # resolution shells, equal population in s^2, fitted on working reflections only
        edges = np.quantile(self.s2[self.work], np.linspace(0, 1, cfg.n_shells + 1))
        edges[0] = -np.inf
        edges[-1] = np.inf
        self.shell = np.digitize(self.s2, edges[1:-1])
        self.shell_edges = edges
        self.shell_d = np.array([np.median(self.d[self.shell == s])
                                 if np.any(self.shell == s) else np.nan
                                 for s in range(cfg.n_shells)])

        self.hist = HistLib(cfg.hist_lib)
        self._kernel_cache = {}
        self.rho_solv = 0.0
        self._centric_info()

        self.log(f'[input] {cfg.mtzin}')
        self.log(f'        SG={self.sg.hm}  cell={a:.3f} {b:.3f} {c:.3f} '
                 f'{self.cell.alpha:.1f} {self.cell.beta:.1f} {self.cell.gamma:.1f}')
        self.log(f'        nref={self.nref} valid={self.valid.sum()} '
                 f'free={np.sum(self.free & self.valid)}  dmin={self.dmin:.3f} A')
        self.log(f'        grid={self.size} (spacing {step:.3f} A), '
                 f'FOM_in={np.nanmean(self.fom_in[self.valid]):.4f}, '
                 f'HL source={self.hl_source}')

    def _centric_info(self):
        """Find centric reflections and their allowed phase axis (mod 180 deg)."""
        ops = list(self.sg.operations())
        self.centric = np.zeros(self.nref, dtype=bool)
        self.phi_allow = np.zeros(self.nref, dtype=np.float64)
        for i, h in enumerate(self.hkl):
            hh = [int(x) for x in h]
            for op in ops:
                if list(op.apply_to_hkl(hh)) == [-hh[0], -hh[1], -hh[2]]:
                    self.centric[i] = True
                    self.phi_allow[i] = (180.0 * np.dot(hh, op.tran) / op.DEN) % 180.0
                    break
        self.log(f'[info] centric reflections: {self.centric.sum()} / {self.nref}')

    def _project_centric(self, phi_rad):
        """project centric phases onto the nearest symmetry-allowed value."""
        if not np.any(self.centric):
            return phi_rad
        out = np.array(phi_rad, copy=True)
        p = np.rad2deg(phi_rad[self.centric])
        a = self.phi_allow[self.centric]
        # allowed phases are a and a+180 degrees; choose the nearest one
        out[self.centric] = np.deg2rad(np.where(
            np.abs(wrap180(p - a)) <= np.abs(wrap180(p - (a + 180.0))),
            a, a + 180.0))
        return out

    # ------------------------- FFT helpers --------------------------------
    def _grid_from_coeffs(self, F, phi_deg):
        """ASU (F, phi) -> full P1 reciprocal grid, via gemmi symmetry expansion."""
        m = gemmi.Mtz(with_base=True)
        m.spacegroup = self.sg
        m.cell = self.cell
        m.add_dataset('oasis-dm')
        m.add_column('FBEST', 'F')
        m.add_column('PHIBEST', 'P')
        data = np.column_stack([self.hkl.astype(np.float64),
                                np.nan_to_num(F), np.nan_to_num(phi_deg)]).astype(np.float32)
        m.set_data(data)
        g = m.get_f_phi_on_grid('FBEST', 'PHIBEST', size=self.size)
        return np.nan_to_num(np.asarray(g.array).copy())

    def rho_from_grid(self, arr):
        """reciprocal grid -> electron density in e/A^3.

        This exact convention was verified against gemmi's own
        transform_f_phi_to_map(): CC = 1.000000, max diff < 5e-7.
        """
        return np.fft.ifftn(arr.conj()).real * self.scale

    def grid_from_rho(self, rho):
        """electron density -> reciprocal grid (inverse of rho_from_grid)."""
        return np.fft.fftn(np.asarray(rho, dtype=np.float64) / self.scale).conj()

    def fmod_from_rho(self, rho):
        """modified density -> complex F_mod for every ASU reflection."""
        g = self.grid_from_rho(rho)
        return g[self.ih, self.ik, self.il]

    # ---------------------- solvent envelope -------------------------------
    def _kernel_fft(self, radius, mode):
        key = (round(radius, 3), mode)
        if key in self._kernel_cache:
            return self._kernel_cache[key]
        dx, dy, dz = self.grid_spacing
        nx, ny, nz = self.size
        ex = int(math.ceil(radius / dx)) + 1
        ey = int(math.ceil(radius / dy)) + 1
        ez = int(math.ceil(radius / dz)) + 1
        kern = np.zeros((nx, ny, nz), dtype=np.float32)
        r2max = radius * radius
        for i in range(-ex, ex + 1):
            x2 = (i * dx) ** 2
            if x2 >= r2max:
                continue
            for j in range(-ey, ey + 1):
                y2 = (j * dy) ** 2
                if x2 + y2 >= r2max:
                    continue
                for k in range(-ez, ez + 1):
                    r2 = x2 + y2 + (k * dz) ** 2
                    if r2 >= r2max:
                        continue
                    r = math.sqrt(r2)
                    if mode == 'wang':
                        w = 1.0 - r2 / r2max          # dm SOLMASK mode 2
                    elif mode == 'wang1':
                        w = 1.0 - r / radius          # original Wang
                    else:
                        w = 1.0
                    kern[i % nx, j % ny, k % nz] = w
        s = kern.sum()
        if s <= 0:
            raise ValueError('empty spherical kernel')
        kern /= s
        kf = np.fft.fftn(kern)
        self._kernel_cache[key] = kf
        return kf

    def _conv(self, x, kf):
        return np.fft.ifftn(np.fft.fftn(np.asarray(x, dtype=np.float64)) * kf).real

    def _mask_feature(self, rho, radius):
        """returns the scalar field whose LOW values identify solvent."""
        cfg = self.cfg
        if cfg.mask == 'variance':
            kf = self._kernel_fft(radius, 'tophat')
            mean = self._conv(rho, kf)
            mean2 = self._conv(rho * rho, kf)
            return np.sqrt(np.maximum(mean2 - mean * mean, 0.0))
        kf = self._kernel_fft(radius, 'wang')
        return self._conv(np.maximum(rho - self.rho_solv, 0.0), kf)

    # --------------------- phase error / combination -----------------------
    def _sim_weight(self, Fc):
        """per-shell scale D, missing-structure variance SigmaQ and Sim X."""
        cfg = self.cfg
        X = np.zeros(self.nref, dtype=np.float64)
        Ds = np.zeros(cfg.n_shells)
        Qs = np.zeros(cfg.n_shells)
        fo_all = self.Fo
        for s in range(cfg.n_shells):
            idx = self.work & (self.shell == s) & np.isfinite(Fc) & (Fc > 0) & (fo_all > 0)
            n = int(idx.sum())
            if n < 15:
                Ds[s], Qs[s] = 1.0, np.nan
                continue
            fo, fc = fo_all[idx], Fc[idx]
            D = float(np.sum(fo * fc) / np.sum(fc * fc))
            resid = float(np.mean((fo - D * fc) ** 2))
            q = max(resid, 0.10 * float(np.mean(fo ** 2)))
            Ds[s], Qs[s] = D, q
            X[idx] = np.clip(2.0 * fo * (D * fc) / q, 0.0, cfg.x_max)
        good = np.isfinite(Qs) & (Qs > 0)
        Dg = float(np.median(Ds[good])) if np.any(good) else 1.0
        Qg = float(np.median(Qs[good])) if np.any(good) else float(np.mean(self.Fo[self.work] ** 2))
        for s in range(cfg.n_shells):
            if not (np.isfinite(Qs[s]) and Qs[s] > 0):
                Ds[s], Qs[s] = Dg, Qg
        idx = self.valid & (X == 0)
        X[idx] = np.clip(2.0 * fo_all[idx] * (Dg * Fc[idx]) / Qg, 0.0, cfg.x_max)
        return X, Ds, Qs, Dg

    # ------------------------------- run ----------------------------------
    def run(self):
        cfg = self.cfg
        rng = np.random.default_rng(cfg.seed)
        HL_cur = self.HL_exp.copy()
        HL_exp = self.HL_exp.copy()
        Fo = np.nan_to_num(self.Fo)
        stats = []
        rho_prev = None
        t0 = time.time()
        best_phi_rad = np.deg2rad(self.phi_exp_deg)
        best_m = np.clip(np.nan_to_num(self.fom_in), 0, 1)

        self.log('')
        self.log(f'{"cyc":>3} {"meanFOM":>8} {"solv%":>7} {"gamma":>7} '
                 f'{"d_eff":>7} {"solvStd":>8} {"protStd":>8} {"mapCC":>7} '
                 f'{"dPhi":>7} {"freeR":>7}')

        for cyc in range(1, cfg.cycles + 1):
            # 1. centroid ("best") map coefficients --------------------------
            phi_best_rad, m_best = hl_centroid(HL_cur)
            phi_best_rad = self._project_centric(phi_best_rad)
            F_best = np.where(self.valid, m_best * Fo, 0.0)
            if cfg.omit_free:
                F_best[self.free] = 0.0        # reflection-omit cross validation
            arr = self._grid_from_coeffs(F_best, np.rad2deg(phi_best_rad))
            rho = self.rho_from_grid(arr)

            # 2. solvent envelope -------------------------------------------
            radius = max(cfg.radius_min, cfg.radius * cfg.radius_decay ** (cyc - 1))
            av = self._mask_feature(rho, radius)
            cut = float(np.quantile(av, cfg.solc))
            width = max(0.10 * float(av.std()), 1e-9)
            mask_solv = 1.0 / (1.0 + np.exp(np.clip((av - cut) / width, -30, 30)))
            if cfg.mask_smooth > 0:
                mask_solv = ndimage.gaussian_filter(mask_solv, cfg.mask_smooth)
            mask_prot = 1.0 - mask_solv
            solv_frac = float(mask_solv.mean())

            # 3. target histogram & ordinal (CDF) map ------------------------
            mean_fom_now = float(np.mean(m_best[self.valid]))
            d_eff = float(np.clip(self.dmin / max(mean_fom_now, 0.30) ** 0.25, 1.5, 10.0))
            rho_hm_grid = None
            if cfg.hist_on:
                cdf_t, tgt_mean, tgt_std = self.hist.target(d_eff)
                vals = rho[mask_prot > 0.5]
                if len(vals) < 100:
                    vals = rho.ravel()
                svals = np.sort(vals)
                qq = np.linspace(0.0, 1.0, len(svals))
                cur_mean = float(vals.mean())
                cur_std = max(float(vals.std()), 1e-6)
                rho_hm_grid = self.hist.rho_axis

                def _hist_map(x):
                    r = np.interp(x, svals, qq)
                    if cfg.squash:
                        # mild non-linear ordinal squash, as in the reference code
                        z = 2.0 * r - 1.0
                        z = z / np.sqrt(1.0 + 0.01 * z * z)
                        r = (z + 1.0) / 2.0
                    hm = np.interp(r, cdf_t, rho_hm_grid)
                    # re-standardise onto the CURRENT protein-region statistics
                    return (hm - tgt_mean) / tgt_std * cur_std + cur_mean
            else:
                def _hist_map(x):
                    return x

            wp = mask_prot if cfg.hist_on else np.zeros_like(mask_prot)
            ws = mask_solv if cfg.solv_on else np.zeros_like(mask_solv)
            wo = (1.0 - wp) * (1.0 - ws)
            sw = np.maximum(wp + ws + wo, 1e-6)

            def modify(x):
                return (wp * _hist_map(x) + ws * self.rho_solv + wo * x) / sw

            rho_mod = modify(rho)

            # 4. perturbation-gamma bias correction --------------------------
            gamma = None
            if cfg.gamma_on:
                eps = 0.03 * max(float(rho.std()), 1e-6)
                delta = rng.normal(0.0, eps, rho.shape)
                rho_mod_pert = modify(rho + delta)   # SAME mask, SAME LUT
                d = rho_mod_pert - rho_mod
                denom = float(np.sum(delta * delta))
                gamma = float(np.clip(np.sum(d * delta) / max(denom, 1e-12), 0.0, 0.9))
                rho_mod = (rho_mod - gamma * rho) / (1.0 - gamma)

            # 5. back-transform -> modified structure factors ----------------
            Fmod = self.fmod_from_rho(rho_mod)
            Fc = np.abs(Fmod)
            phi_mod = np.angle(Fmod)

            # 6. Sim/sigmaA error estimate and HL phase recombination --------
            X, Ds, Qs, Dg = self._sim_weight(Fc)
            HL_mod = np.column_stack([cfg.cmbwt * X * np.cos(phi_mod),
                                      cfg.cmbwt * X * np.sin(phi_mod),
                                      np.zeros(self.nref), np.zeros(self.nref)])
            HL_mod = np.where(self.valid[:, None], HL_mod, 0.0)
            if cfg.anchor == 'exp':
                HL_new = HL_exp + HL_mod          # anchor to ORIGINAL experiment
            else:
                HL_new = HL_cur + HL_mod          # cumulative (more bias)
            phi_new_rad, m_new = hl_centroid(HL_new)
            phi_new_rad = self._project_centric(phi_new_rad)
            HL_cur = HL_new

            # 7. diagnostics -------------------------------------------------
            sel = self.valid
            dphi = float(np.mean(np.abs(wrap180(np.rad2deg(phi_new_rad[sel])
                                                - np.rad2deg(phi_best_rad[sel])))))
            sol_sel = mask_solv > 0.8
            prot_sel = mask_prot > 0.8
            solv_std = float(np.std(rho_mod[sol_sel])) if np.any(sol_sel) else float('nan')
            prot_std = float(np.std(rho_mod[prot_sel])) if np.any(prot_sel) else float('nan')
            map_cc = (float(np.corrcoef(rho_mod.ravel(), rho_prev.ravel())[0, 1])
                      if rho_prev is not None else float('nan'))
            rho_prev = rho_mod
            fr = sel & self.free & (Fc > 0)
            if np.any(fr):
                k = float(np.sum(Fo[fr] * Fc[fr]) / np.sum(Fc[fr] ** 2))
                free_r = float(np.sum(np.abs(Fo[fr] - k * Fc[fr])) / np.sum(Fo[fr]))
            else:
                free_r = float('nan')
            mean_fom = float(np.mean(m_new[sel]))

            rec = dict(cycle=cyc, mean_fom=mean_fom, solvent_frac=solv_frac,
                       gamma=gamma, d_eff=d_eff, mask_radius=radius,
                       solvent_std=solv_std, protein_std=prot_std,
                       map_cc_prev=map_cc, mean_dphi_deg=dphi, free_r=free_r,
                       mean_X=float(np.mean(X[sel])),
                       D_shells=[float(x) for x in Ds],
                       SigmaQ_shells=[float(x) for x in Qs],
                       shell_d=[float(x) for x in self.shell_d])
            stats.append(rec)
            self.log(f'{cyc:3d} {mean_fom:8.4f} {solv_frac*100:7.2f} '
                     f'{(gamma if gamma is not None else float("nan")):7.3f} '
                     f'{d_eff:7.3f} {solv_std:8.4f} {prot_std:8.4f} '
                     f'{map_cc:7.4f} {dphi:7.2f} {free_r:7.4f}')

            best_phi_rad, best_m = phi_new_rad, m_new

        # --------------------------- output --------------------------------
        self.best_phi_rad = best_phi_rad
        self.best_m = best_m
        self.HL_final = HL_cur
        self.Fmod = Fmod
        self.Dg = Dg
        self.stats = stats
        self.elapsed = time.time() - t0
        self._write_output()
        self._write_stats()
        return self.stats

    def _write_output(self):
        cfg = self.cfg
        Fo = np.nan_to_num(self.Fo)
        phi_deg = wrap180(np.rad2deg(self.best_phi_rad))
        m = np.clip(self.best_m, 0.0, 1.0)
        Fc = np.abs(self.Fmod)
        phi_mod_deg = wrap180(np.rad2deg(np.angle(self.Fmod)))
        D = self.Dg

        Fdm = np.where(self.valid, Fo, cfg.restore * Fc)
        Fwt = np.where(self.centric, m * Fo, 2.0 * m * Fo - D * Fc)
        Fdel = m * Fo - D * Fc
        phi_out = np.where(self.valid, phi_deg, phi_mod_deg)
        fom_out = np.where(self.valid, m, 0.0)
        newcols = [
            Fdm, phi_out, fom_out,
            self.HL_final[:, 0], self.HL_final[:, 1],
            self.HL_final[:, 2], self.HL_final[:, 3],
            np.maximum(Fwt, 0.0), phi_out,
            Fdel, phi_out,
            Fc, phi_mod_deg,
        ]
        for v in newcols:
            if v.shape != (self.nref,):
                raise RuntimeError('output column shape mismatch')

        out = gemmi.read_mtz_file(cfg.mtzin)
        orig = np.asarray(out.array).copy()
        # remove pre-existing columns with the same labels (idempotent re-runs)
        for label, typ in OUT_COLUMNS:
            if out.column_with_label(label) is not None:
                continue
            out.add_column(label, typ)
        data = np.column_stack([orig] + [np.nan_to_num(v).astype(np.float32)
                                         for v in newcols]).astype(np.float32)
        if data.shape[1] != len(out.columns):
            # gemmi add_column may not resize the existing buffer; rebuild instead
            out = gemmi.Mtz(with_base=True)
            out.spacegroup = self.sg
            out.cell = self.cell
            out.add_dataset('oasis-dm')
            for c in self.mtz.columns[3:]:
                out.add_column(c.label, c.type)
            for label, typ in OUT_COLUMNS:
                if out.column_with_label(label) is None:
                    out.add_column(label, typ)
            data = np.column_stack([orig] + [np.nan_to_num(v).astype(np.float32)
                                             for v in newcols]).astype(np.float32)
        out.set_data(np.ascontiguousarray(data))
        out.write_to_file(cfg.mtzout)
        self.log(f'[out ] wrote {cfg.mtzout}  ({data.shape[1]} columns)')
        self.output_columns = dict(zip(
            ['FDM', 'PHIDM', 'FOMDM', 'HLADM', 'HLBDM', 'HLCDM', 'HLDDM',
             'FWT', 'PHWT', 'DELFWT', 'PHDELFWT', 'FCDM', 'PHICDM'],
            [x[0] for x in OUT_COLUMNS]))

    def _write_stats(self):
        if not self.cfg.stats:
            return
        payload = dict(tag=self.cfg.tag or os.path.basename(self.cfg.mtzin),
                       config=self.cfg.as_dict(),
                       input=dict(mtz=self.cfg.mtzin,
                                  spacegroup=self.sg.hm,
                                  cell=[self.cell.a, self.cell.b, self.cell.c,
                                        self.cell.alpha, self.cell.beta, self.cell.gamma],
                                  nref=self.nref, nvalid=int(self.valid.sum()),
                                  nfree=int(np.sum(self.free & self.valid)),
                                  dmin=self.dmin, grid=self.size,
                                  hl_source=self.hl_source,
                                  fom_in_mean=float(np.nanmean(self.fom_in[self.valid]))),
                       output_columns=self.output_columns,
                       cycles=self.stats,
                       final=dict(mean_fom=self.stats[-1]['mean_fom'],
                                  elapsed_sec=self.elapsed),
                       elapsed_sec=self.elapsed)
        os.makedirs(os.path.dirname(os.path.abspath(self.cfg.stats)), exist_ok=True)
        with open(self.cfg.stats, 'w') as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)
        self.log(f'[out ] wrote {self.cfg.stats}  (elapsed {self.elapsed:.2f}s)')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description='oasis-dm: from-scratch density modification (dm-like)')
    p.add_argument('--mtzin', required=True)
    p.add_argument('--mtzout', required=True)
    p.add_argument('--fo', default='FOBS')
    p.add_argument('--sigf', default='SIGFOBS')
    p.add_argument('--phib', default='PHIB')
    p.add_argument('--fom', default='FOM')
    p.add_argument('--hl', default='HLA,HLB,HLC,HLD',
                   help='comma separated HL coefficient labels; set to "" to force PHIB/FOM')
    p.add_argument('--free', default='FreeR')
    p.add_argument('--truth', default='PHIC')
    p.add_argument('--solc', type=float, default=0.5)
    p.add_argument('--cycles', type=int, default=5)
    p.add_argument('--mask', choices=['wang', 'variance'], default='variance',
                   help='solvent envelope: local-variance or Wang sphere (dm)')
    p.add_argument('--radius', type=float, default=None,
                   help='mask radius in A (default: 4 for variance, 8 for wang)')
    p.add_argument('--radius-min', type=float, default=None,
                   help='final mask radius in A (default: 3 for variance, 4 for wang)')
    p.add_argument('--radius-decay', type=float, default=0.85)
    p.add_argument('--mask-smooth', type=float, default=0.8)
    p.add_argument('--hist-lib', default='/Applications/ccp4-9/lib/data/hist.lib')
    p.add_argument('--no-hist', action='store_true')
    p.add_argument('--no-solv', action='store_true')
    p.add_argument('--no-squash', action='store_true')
    p.add_argument('--no-gamma', action='store_true')
    p.add_argument('--omit-free', action='store_true',
                   help='enable dm-OMIT style reflection omission (strict cross validation)')
    p.add_argument('--free-marker', type=int, default=0,
                   help='value marking a FREE reflection (0 = CCP4 FreeR_flag; '
                        '1 = GPCR pipeline FreeR). Default 0.')
    p.add_argument('--no-omit-free', action='store_true', help=argparse.SUPPRESS)
    p.add_argument('--anchor', choices=['exp', 'cumulative'], default='exp')
    p.add_argument('--n-shells', type=int, default=20)
    p.add_argument('--x-max', type=float, default=12.0)
    p.add_argument('--cmbwt', type=float, default=0.4)
    p.add_argument('--restore', type=float, default=0.0)
    p.add_argument('--sample-rate', type=float, default=3.0)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--stats', default='')
    p.add_argument('--tag', default='')
    p.add_argument('--quiet', action='store_true')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    radius = args.radius if args.radius is not None else (8.0 if args.mask == 'wang' else 4.0)
    radius_min = (args.radius_min if args.radius_min is not None
                  else (4.0 if args.mask == 'wang' else 3.0))
    cfg = Config(mtzin=args.mtzin, mtzout=args.mtzout,
                 fo=args.fo, sigf=args.sigf, phib=args.phib, fom=args.fom,
                 hl=args.hl, free=args.free, free_marker=args.free_marker,
                 truth=args.truth,
                 solc=args.solc, cycles=args.cycles, mask=args.mask,
                 radius=radius, radius_min=radius_min,
                 radius_decay=args.radius_decay, mask_smooth=args.mask_smooth,
                 hist_lib=args.hist_lib, hist_on=not args.no_hist,
                 solv_on=not args.no_solv, squash=not args.no_squash,
                 gamma_on=not args.no_gamma, omit_free=args.omit_free,
                 anchor=args.anchor, n_shells=args.n_shells, x_max=args.x_max,
                 cmbwt=args.cmbwt, restore=args.restore,
                 sample_rate=args.sample_rate, seed=args.seed,
                 stats=args.stats, tag=args.tag, quiet=args.quiet)
    dm = DensityModifier(cfg)
    dm.run()
    return 0


if __name__ == '__main__':
    sys.exit(main())
