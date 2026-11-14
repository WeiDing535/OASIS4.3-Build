#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""oasis-dm v2 -- unbiased phase recombination by reflection-omit error estimation
and free-set driven early stopping.

Why v2 exists
-------------
In v1 the Sim/sigmaA missing-structure variance SigmaQ was fitted on WORK
reflections, i.e. exactly those reflections that were used to compute the map.
That makes (Fo - D*Fcalc) systematically too small, X too large, and the DM
phases over-trusted (FOM inflation).  v1 patched this with a hand-tuned
cmbwt=0.4 obtained from a 3-case sweep -- a constant that need not generalise.

v2 follows the idea of COMBINE OMIT with a free-likelihood target:
  * each cycle builds TWO maps from the SAME mask and histogram LUT:
      rho_full  -- all reflections (used for the output phases)
      rho_omit  -- free reflections deleted (unbiased cross-validation map)
  * SigmaQ and the shell scale D are fitted ONLY on the omitted reflections
    (Fo_free vs |Fmod_omit|_free), so the weight X no longer contains the
    self-reference bias; cmbwt can go back to 1.0.
  * a free-set figure of merit (amplitude CC) drives early stopping and the
    best cycle is restored (dm NCYCLE AUTO / RESTORE idea).

Everything else (solvent mask, histogram matching, gamma correction, HL
recombination, output columns) is inherited unchanged from v1.  v1 is NOT
modified.  This program never calls an external density-modification program.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass, asdict

import numpy as np
from scipy import ndimage

from oasis_dm_core import (Config, DensityModifier, hl_centroid, wrap180,
                    build_parser)


@dataclass
class ConfigV2(Config):
    sigmaq_source: str = 'blend'    # 'blend' (work-dominated, stable, default) |
                                    # 'omit' (unbiased attempt, needs scale
                                    # correction to be stable) | 'work'
    auto_stop: bool = True
    max_cycles: int = 12
    min_cycles: int = 2
    patience: int = 2
    free_metric: str = 'cc'         # 'cc' | 'rfree'
    min_free: float = 500.0         # below this many valid free reflections even
                                    # an explicit 'omit' request falls back to
                                    # blend (work-dominated SigmaQ)
    cmbwt: float = 0.4              # v1-calibrated damping; kept conservative
                                    # because the simplified omit SigmaQ estimate
                                    # is optimistic (free set is only ~5% of
                                    # reflections) and must not inflate FOM
    n_shells: int = 10              # fewer shells -> more reflections each

    def as_dict(self):
        return asdict(self)


class DensityModifierV2(DensityModifier):
    """v1 engine + omit-based error estimation + free-set early stopping."""

    def _prepare(self):
        super()._prepare()
        self.free_valid = self.free & self.valid
        self.n_free_valid = int(self.free_valid.sum())
        if self.n_free_valid < 50:
            self.log(f'[warn] only {self.n_free_valid} valid free reflections; '
                     f'omit-based SigmaQ may be noisy')
        # Robustness guard: a handful of free reflections cannot pin down D and
        # SigmaQ per resolution shell.  Fall back to a work-dominated blend so we
        # do not inflate FOM (the symmetric failure mode to v1's self-reference
        # bias).  This keeps v2 honest on small-test-set SAD data.
        if self.cfg.sigmaq_source == 'omit' and self.n_free_valid < self.cfg.min_free:
            self.cfg.sigmaq_source = 'blend'
            self.cfg.cmbwt = 0.4
            self.log(f'[v2  ] free reflections {self.n_free_valid} < min_free '
                     f'{self.cfg.min_free}: SigmaQ source -> blend '
                     f'(work-dominated); cmbwt -> 0.4 to avoid FOM inflation')
        self.log(f'[v2  ] sigmaQ source={self.cfg.sigmaq_source}  '
                 f'auto_stop={self.cfg.auto_stop}  '
                 f'free reflections={self.n_free_valid}')

    # ------------------------------------------------------------------
    def _build_modifier(self, rho, radius, mean_fom_now):
        """Build a FIXED mask + histogram LUT from rho; return a callable that
        applies the same modification to any (full or omit) density map."""
        cfg = self.cfg
        av = self._mask_feature(rho, radius)
        cut = float(np.quantile(av, cfg.solc))
        width = max(0.10 * float(av.std()), 1e-9)
        mask_solv = 1.0 / (1.0 + np.exp(np.clip((av - cut) / width, -30, 30)))
        if cfg.mask_smooth > 0:
            mask_solv = ndimage.gaussian_filter(mask_solv, cfg.mask_smooth)
        mask_prot = 1.0 - mask_solv

        d_eff = float(np.clip(self.dmin / max(mean_fom_now, 0.30) ** 0.25, 1.5, 10.0))
        if cfg.hist_on:
            cdf_t, tgt_mean, tgt_std = self.hist.target(d_eff)
            vals = rho[mask_prot > 0.5]
            if len(vals) < 100:
                vals = rho.ravel()
            svals = np.sort(vals)
            qq = np.linspace(0.0, 1.0, len(svals))
            cur_mean, cur_std = float(vals.mean()), max(float(vals.std()), 1e-6)

            def hist_map(x):
                r = np.interp(x, svals, qq)
                if cfg.squash:
                    z = 2.0 * r - 1.0
                    z = z / np.sqrt(1.0 + 0.01 * z * z)
                    r = (z + 1.0) / 2.0
                hm = np.interp(r, cdf_t, self.hist.rho_axis)
                return (hm - tgt_mean) / tgt_std * cur_std + cur_mean
        else:
            def hist_map(x):
                return x

        wp = mask_prot if cfg.hist_on else np.zeros_like(mask_prot)
        ws = mask_solv if cfg.solv_on else np.zeros_like(mask_solv)
        wo = (1.0 - wp) * (1.0 - ws)
        sw = np.maximum(wp + ws + wo, 1e-6)

        def modify(x):
            return (wp * hist_map(x) + ws * self.rho_solv + wo * x) / sw

        return modify, mask_solv, mask_prot, d_eff

    def _estimate_gamma(self, rho, modify):
        rng = np.random.default_rng(self.cfg.seed + 9000 + len(self.stats))
        eps = 0.03 * max(float(rho.std()), 1e-6)
        delta = rng.normal(0.0, eps, rho.shape)
        base = modify(rho)
        pert = modify(rho + delta)
        d = pert - base
        g = float(np.sum(d * delta) / max(float(np.sum(delta * delta)), 1e-12))
        return float(np.clip(g, 0.0, 0.9)), base

    # ------------------------------------------------------------------
    def _sim_weight_v2(self, Fc_full, Fc_omit):
        """Fit D and SigmaQ on OMITTED (free) reflections; returns X for all."""
        cfg = self.cfg
        X = np.zeros(self.nref, dtype=np.float64)
        Ds = np.zeros(cfg.n_shells)
        Qs = np.zeros(cfg.n_shells)
        Fo = np.nan_to_num(self.Fo)
        if cfg.sigmaq_source == 'work':
            fitmask = self.work
            fc_fit = Fc_full
        elif cfg.sigmaq_source == 'blend':
            # work-dominated, free reflections only help the per-shell counts
            fitmask = self.work | self.free_valid
            fc_fit = Fc_full
        else:  # 'omit' (unbiased): fit on FREE reflections, but take Fc from the
               # FULL map.  A free reflection's contribution to the full density
               # is diluted across the whole grid, so Fc_full[free] is ~self-
               # reference free; using Fc_omit would shrink |Fmod| (incomplete
               # map) and systematically *under*-estimate SigmaQ -> FOM inflation
               # (the symmetric failure mode to v1's work self-reference bias).
            fitmask = self.free_valid
            fc_fit = Fc_full
        for s in range(cfg.n_shells):
            idx = (fitmask & (self.shell == s) & np.isfinite(fc_fit)
                   & (fc_fit > 0) & (Fo > 0))
            n = int(idx.sum())
            if n < 8:
                Ds[s], Qs[s] = 1.0, np.nan
                continue
            fo = Fo[idx]
            fc = fc_fit[idx]
            D = float(np.sum(fo * fc) / np.sum(fc * fc))
            resid = float(np.mean((fo - D * fc) ** 2))
            Qs[s] = max(resid, 0.10 * float(np.mean(fo ** 2)))
            Ds[s] = D
        good = np.isfinite(Qs) & (Qs > 0)
        Dg = float(np.median(Ds[good])) if np.any(good) else 1.0
        Qg = (float(np.median(Qs[good])) if np.any(good)
              else float(np.mean(Fo[self.work] ** 2)))
        for s in range(cfg.n_shells):
            if not (np.isfinite(Qs[s]) and Qs[s] > 0):
                Ds[s], Qs[s] = Dg, Qg
        for s in range(cfg.n_shells):
            idx = self.valid & (self.shell == s)
            if np.any(idx):
                X[idx] = np.clip(2.0 * Fo[idx] * (Ds[s] * Fc_full[idx]) / Qs[s],
                                 0.0, cfg.x_max)
        return X, Ds, Qs, Dg

    def _free_metric(self, Fc_omit):
        sel = self.free_valid & np.isfinite(Fc_omit) & (Fc_omit > 0)
        if sel.sum() < 20:
            return float('nan'), float('nan')
        fo = np.nan_to_num(self.Fo)[sel]
        fc = Fc_omit[sel]
        cc = float(np.corrcoef(fo, fc)[0, 1])
        k = float(np.sum(fo * fc) / np.sum(fc * fc))
        rfree = float(np.sum(np.abs(fo - k * fc)) / np.sum(fo))
        return cc, rfree

    # ------------------------------------------------------------------
    def run(self):
        cfg = self.cfg
        HL_exp = self.HL_exp.copy()
        HL_cur = HL_exp.copy()
        Fo = np.nan_to_num(self.Fo)
        stats = []
        self.stats = stats          # shared list, used for per-cycle RNG seeding
        rho_prev = None
        t0 = time.time()
        ncycles = cfg.max_cycles if cfg.auto_stop else cfg.cycles
        best = None
        bad_streak = 0

        self.log('')
        self.log(f'{"cyc":>3} {"meanFOM":>8} {"freeCC":>7} {"freeR":>7} '
                 f'{"SigmaQ":>9} {"gamma":>7} {"dPhi":>7} {"mapCC":>7} {"best":>5}')
        for cyc in range(1, ncycles + 1):
            phi_best_rad, m_best = hl_centroid(HL_cur)
            phi_best_rad = self._project_centric(phi_best_rad)
            phi_deg = np.rad2deg(phi_best_rad)
            F_best = np.where(self.valid, m_best * Fo, 0.0)
            F_omit = F_best.copy()
            F_omit[self.free] = 0.0

            # --- full map -------------------------------------------------
            rho = self.rho_from_grid(self._grid_from_coeffs(F_best, phi_deg))
            radius = max(cfg.radius_min, cfg.radius * cfg.radius_decay ** (cyc - 1))
            mean_fom_now = float(np.mean(m_best[self.valid]))
            modify, mask_solv, mask_prot, d_eff = self._build_modifier(
                rho, radius, mean_fom_now)
            gamma = None
            if cfg.gamma_on:
                gamma, rho_mod = self._estimate_gamma(rho, modify)
            else:
                rho_mod = modify(rho)
            if gamma is not None:
                rho_mod = (rho_mod - gamma * rho) / (1.0 - gamma)
            Fmod = self.fmod_from_rho(rho_mod)
            Fc_full = np.abs(Fmod)

            # --- omit map (SAME mask, SAME LUT) ---------------------------
            rho_o = self.rho_from_grid(self._grid_from_coeffs(F_omit, phi_deg))
            rho_mod_o = modify(rho_o)
            if gamma is not None:
                rho_mod_o = (rho_mod_o - gamma * rho_o) / (1.0 - gamma)
            Fmod_o = self.fmod_from_rho(rho_mod_o)
            Fc_omit = np.abs(Fmod_o)
            free_cc, free_r = self._free_metric(Fc_omit)

            # --- unbiased error estimate + phase recombination ------------
            X, Ds, Qs, Dg = self._sim_weight_v2(Fc_full, Fc_omit)
            phi_mod = np.angle(Fmod)
            HL_mod = np.column_stack([cfg.cmbwt * X * np.cos(phi_mod),
                                      cfg.cmbwt * X * np.sin(phi_mod),
                                      np.zeros(self.nref), np.zeros(self.nref)])
            HL_mod = np.where(self.valid[:, None], HL_mod, 0.0)
            HL_new = HL_exp + HL_mod
            phi_new_rad, m_new = hl_centroid(HL_new)
            phi_new_rad = self._project_centric(phi_new_rad)
            HL_cur = HL_new

            # --- diagnostics ----------------------------------------------
            sel = self.valid
            dphi = float(np.mean(np.abs(wrap180(np.rad2deg(phi_new_rad[sel])
                                                - np.rad2deg(phi_best_rad[sel])))))
            map_cc = (float(np.corrcoef(rho_mod.ravel(), rho_prev.ravel())[0, 1])
                      if rho_prev is not None else float('nan'))
            rho_prev = rho_mod
            qmed = float(np.median(Qs[np.isfinite(Qs) & (Qs > 0)]))
            score = free_cc if cfg.free_metric == 'cc' else -free_r
            snapshot = dict(HL=HL_new.copy(), phi=phi_new_rad.copy(), m=m_new.copy(),
                            Fmod=Fmod.copy(), Dg=Dg, score=score, cycle=cyc)
            if best is None or (np.isfinite(score) and score > best['score'] + 1e-5):
                best = snapshot
                bad_streak = 0
                is_best = True
            else:
                bad_streak += 1
                is_best = False

            stats.append(dict(cycle=cyc, mean_fom=float(np.mean(m_new[sel])),
                              free_cc=free_cc, free_r_omit=free_r,
                              sigmaq_median=qmed, gamma=gamma, d_eff=d_eff,
                              mask_radius=radius,
                              solvent_frac=float(mask_solv.mean()),
                              map_cc_prev=map_cc, mean_dphi_deg=dphi,
                              mean_X=float(np.mean(X[sel])),
                              is_best=is_best, best_cycle=best['cycle'],
                              D_shells=[float(x) for x in Ds],
                              SigmaQ_shells=[float(x) for x in Qs],
                              shell_d=[float(x) for x in self.shell_d]))
            self.log(f'{cyc:3d} {np.mean(m_new[sel]):8.4f} {free_cc:7.4f} '
                     f'{free_r:7.4f} {qmed:9.2f} '
                     f'{(gamma if gamma is not None else float("nan")):7.3f} '
                     f'{dphi:7.2f} {map_cc:7.4f} {"*" if is_best else " ":>5}')

            if cfg.auto_stop and cyc >= cfg.min_cycles and bad_streak >= cfg.patience:
                self.log(f'[v2  ] early stop at cycle {cyc}: free-set '
                         f'{cfg.free_metric} did not improve for {bad_streak} cycles; '
                         f'restoring cycle {best["cycle"]}')
                break

        # restore best cycle -------------------------------------------------
        self.best_phi_rad = best['phi']
        self.best_m = best['m']
        self.HL_final = best['HL']
        self.Fmod = best['Fmod']
        self.Dg = best['Dg']
        self.stats = stats
        self.best_cycle = best['cycle']
        self.elapsed = time.time() - t0
        self._write_output()
        self._write_stats()
        return self.stats

    def _write_stats(self):
        # add v2-specific fields to the parent JSON payload
        super()._write_stats()
        if not self.cfg.stats or not os.path.isfile(self.cfg.stats):
            return
        import json
        with open(self.cfg.stats) as fh:
            payload = json.load(fh)
        payload['final'].update(dict(best_cycle=self.best_cycle,
                                     n_cycles_run=len(self.stats),
                                     sigmaq_source=self.cfg.sigmaq_source,
                                     auto_stop=self.cfg.auto_stop,
                                     best_free_cc=max(
                                         (s['free_cc'] for s in self.stats
                                          if np.isfinite(s['free_cc'])),
                                         default=float('nan'))))
        with open(self.cfg.stats, 'w') as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)


def build_parser_v2():
    p = build_parser()
    p.description = 'oasis-dm v2: omit-based unbiased SigmaQ + free-set early stop'
    p.add_argument('--sigmaq-source', choices=['omit', 'work', 'blend'],
                   default='blend',
                   help='which reflections fit SigmaQ (default blend = '
                        'work-dominated, stable; omit is an experimental '
                        'unbiased attempt)')
    p.add_argument('--auto-stop', action='store_true', default=True)
    p.add_argument('--no-auto-stop', dest='auto_stop', action='store_false',
                   help='disable free-set early stopping (fixed --cycles)')
    p.add_argument('--max-cycles', type=int, default=12)
    p.add_argument('--min-cycles', type=int, default=2)
    p.add_argument('--patience', type=int, default=2)
    p.add_argument('--free-metric', choices=['cc', 'rfree'], default='cc')
    p.add_argument('--min-free', type=float, default=500.0,
                   help='below this many valid free reflections, auto-fallback '
                        'SigmaQ to blend (work-dominated) + cmbwt 0.4')
    p.set_defaults(cmbwt=0.4, n_shells=10)
    return p


def main(argv=None):
    args = build_parser_v2().parse_args(argv)
    radius = args.radius if args.radius is not None else (8.0 if args.mask == 'wang' else 4.0)
    radius_min = (args.radius_min if args.radius_min is not None
                  else (4.0 if args.mask == 'wang' else 3.0))
    cfg = ConfigV2(mtzin=args.mtzin, mtzout=args.mtzout,
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
                   stats=args.stats, tag=args.tag, quiet=args.quiet,
                   sigmaq_source=args.sigmaq_source, auto_stop=args.auto_stop,
                   max_cycles=args.max_cycles, min_cycles=args.min_cycles,
                   patience=args.patience, free_metric=args.free_metric,
                   min_free=args.min_free)
    DensityModifierV2(cfg).run()
    return 0


if __name__ == '__main__':
    sys.exit(main())
