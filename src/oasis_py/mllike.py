"""Multivariate likelihood for substructure refinement and SAD phasing.

For every reflection of every dataset the two Bijvoet mates are treated as two
independent observations, split into the channels that actually carry
information:

* the *native* amplitude ``fn = (|F+| + |F-|)/2``, whose calculated value is
  ``(|Fc+| + |Fc-|)/2``.  The native structure factor ``F_P`` is unknown; it is
  integrated out with a Wilson amplitude prior and a phase that is uniform for
  acentric reflections and two-valued for centric ones.  Its error is the
  Wilson variance of the resolution bin, large enough to absorb the model error
  of the native part and to keep the phase integral smooth.
* the *anomalous difference* ``df = |F+| - |F-|``, whose calculated value is
  ``|Fc+| - |Fc-|`` and whose error is the experimental variance of the two
  amplitudes.  Only this channel is sensitive to the anomalous scattering of
  the substructure, and it must not be diluted by the native error.

The substructure itself enters through ``F_H`` and its Friedel mate: writing
``F_H = S_r + i S_i`` the mate is ``S_r - i S_i`` (the anomalous term changes
sign, it is not the complex conjugate), so

    Fc+ = A exp(i alpha) + F_H,      Fc- = A exp(i alpha) + conj-mate.

A Gaussian prior on ``F_H`` with a variance ``sigh`` per resolution bin fixes
the otherwise arbitrary scale of the substructure, because the free native
phase can absorb any ``|F_H|``.  Refined parameters are therefore the atomic
coordinates, occupancies and B factors plus one ``sigh`` per bin.

Because amplitudes -- not complex values -- are observed, the phase integral
does not collapse to a Bessel function, and the amplitude and phase are
integrated on a two-dimensional quadrature that adapts to each reflection.
This is the same statistical model a BP3-class program uses; the implementation
is independent and the external program serves only as a numerical oracle.
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

__all__ = ['Dataset', 'Substructure', 'MultivariateLikelihood']

# Gauss-Legendre nodes for the amplitude integral.
_GL, _GLW = np.polynomial.legendre.leggauss(20)


class Dataset:
    """One native dataset: Friedel pair amplitudes with their standard deviations."""

    def __init__(self, fp, sigfp, fm, sigfm, name='native'):
        self.fp = np.asarray(fp, dtype=float)
        self.sigfp = np.asarray(sigfp, dtype=float)
        self.fm = np.asarray(fm, dtype=float)
        self.sigfm = np.asarray(sigfm, dtype=float)
        self.name = name


class Substructure:
    """Substructure atoms with their scattering factors and symmetry."""

    def __init__(self, elements, xyz, occ, b, hkl, cell, spacegroup, s2):
        self.elements = list(elements)
        self.xyz = np.asarray(xyz, dtype=float)
        self.occ = np.asarray(occ, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self.hkl = np.asarray(hkl, dtype=int)
        self.cell = cell
        self.spacegroup = spacegroup
        self.s2 = np.asarray(s2, dtype=float)
        self.ops = [(np.asarray(op.rot, dtype=float) / op.DEN,
                     np.asarray(op.tran, dtype=float) / op.DEN)
                    for op in spacegroup.operations()]
        self.f0 = np.stack([cromer_mann(e, 0.5 * np.sqrt(self.s2)) for e in self.elements],
                           axis=1)
        self.fp = np.zeros(len(self.elements))
        self.fpp = np.zeros(len(self.elements))

    def set_anomalous(self, fp, fpp):
        self.fp = np.broadcast_to(np.asarray(fp, dtype=float), self.occ.shape).copy()
        self.fpp = np.broadcast_to(np.asarray(fpp, dtype=float), self.occ.shape).copy()

    def factors(self, xyz=None, occ=None, b=None, gradient=False, coef=None):
        """F_H, and optionally dF_H/dp for every refined parameter.

        Parameter order: x, y, z of every atom, then occupancies, then B
        factors, i.e. 5*nat columns.  ``coef`` overrides the per-atom
        scattering factor, which is how the real and imaginary parts of the
        heavy-atom contribution are separated (see :meth:`factors_split`).
        """
        xyz = self.xyz if xyz is None else np.asarray(xyz, dtype=float)
        occ = self.occ if occ is None else np.asarray(occ, dtype=float)
        b = self.b if b is None else np.asarray(b, dtype=float)
        nat = len(occ)
        attenuation = np.exp(-0.25 * self.s2[:, None] * b[None, :])
        if coef is None:
            form = self.f0 + self.fp[None, :] + 1j * self.fpp[None, :]
        else:
            form = np.asarray(coef, dtype=complex)
        amp = occ[None, :] * attenuation * form
        total = np.zeros(len(self.hkl), dtype=complex)
        deriv = np.zeros((len(self.hkl), 5 * nat), dtype=complex) if gradient else None
        for rot, tran in self.ops:
            phase = np.exp(2j * np.pi * (self.hkl @ (xyz @ rot.T + tran).T))
            total += np.sum(amp * phase, axis=1)
            if not gradient:
                continue
            rotated = self.hkl @ rot
            for axis in range(3):
                deriv[:, 3 * np.arange(nat) + axis] += (
                    2j * np.pi * rotated[:, axis, None] * amp * phase)
            deriv[:, 3 * nat:4 * nat] += form * attenuation * phase
            deriv[:, 4 * nat:] += -0.25 * self.s2[:, None] * amp * phase
        return (total, deriv) if gradient else total

    def factors_split(self, xyz=None, occ=None, b=None, gradient=False):
        """Separate the real and imaginary parts of the heavy-atom contribution.

        Writing ``F_H = S_r + i S_i`` with

            S_r = sum_j (f0 + f') exp(2 pi i h x_j)
            S_i = sum_j f''      exp(2 pi i h x_j)

        the Friedel mate is ``F_H(-h) = S_r - i S_i``: the anomalous term enters
        with the opposite sign, it is *not* the complex conjugate of ``F_H(h)``.
        Both ``S_r`` and ``S_i`` are ordinary structure factors with real
        coefficients, so their parameter derivatives have the usual form.
        """
        sr = self.factors(xyz, occ, b, gradient,
                          coef=self.f0 + self.fp[None, :])
        si = self.factors(xyz, occ, b, gradient,
                          coef=self.fpp[None, :].astype(complex))
        if gradient:
            return sr[0], si[0], sr[1], si[1]
        return sr, si


def cromer_mann(element, s):
    """Elastic scattering factor in electrons."""
    from phasebuild.internal_refine import _cromer_mann
    return _cromer_mann(element, s)


class MultivariateLikelihood:
    """Joint likelihood of all datasets for a fixed substructure."""

    def __init__(self, datasets, substructure, bins=9, n_alpha=16):
        self.sub = substructure
        self.datasets = datasets
        self.bins = bins
        self.n_alpha = n_alpha
        d = 1.0 / np.sqrt(substructure.s2)
        edges = np.quantile(d, np.linspace(0.0, 1.0, bins + 1))
        edges[0] -= 1e-6
        self.bin = np.clip(np.digitize(d, edges[1:-1]), 0, bins - 1)
        self.sigma_p = self._wilson()
        self.centric, self.allowed = _centric_phases(substructure)

    def _wilson(self):
        """Wilson amplitude scale (rms |F|) of every resolution bin."""
        amp = np.mean([0.5 * (d.fp + d.fm) for d in self.datasets], axis=0)
        out = np.ones(self.bins)
        for i in range(self.bins):
            mask = self.bin == i
            if mask.any():
                out[i] = np.sqrt(np.mean(amp[mask] ** 2))
        return out

    # ---------------------------------------------------------------- core
    def core(self, sr, si, fp, fm, vp, vm, sigh, index=None, marginal=False):
        """log P(|F+|, |F-| | F_H) for every reflection, one dataset.

        The two Bijvoet mates are split into the two channels that carry
        independent information:

        * the *native* amplitude ``fn = (|F+| + |F-|)/2`` with the Wilson
          variance of its resolution bin.  Its model counterpart
          ``(|Fc+| + |Fc-|)/2`` is dominated by the integrated-out native
          amplitude, so a large error is appropriate and the phase integral
          stays smooth;
        * the *anomalous difference* ``df = |F+| - |F-|`` with the experimental
          variance ``vp + vm``.  Only this channel is sensitive to the anomalous
          part of ``F_H``; diluting it with the native error would make the
          substructure invisible.

        ``sr``/``si`` are the real and imaginary parts of the heavy-atom
        contribution, ``F_H = sr + i si``, so that the Friedel mate is
        ``sr - i si``.  ``sigh`` is the variance of the Gaussian prior on
        ``F_H`` in the resolution bin of each reflection; that prior is what
        fixes the scale of the substructure, because the native phase can
        otherwise absorb an arbitrary ``|F_H|``.

        With ``marginal=True`` the per-reflection native phase distribution is
        returned as well, as ``(logp, logmarginal)``: the joint log-probability
        over the amplitude grid and the phase grid, and its marginal over the
        amplitude.  :meth:`phase_moments` turns that into a phase and a FOM.
        """
        n = len(sr)
        fh = sr + 1j * si
        fhn = sr - 1j * si
        b = self.bin if index is None else self.bin[index]
        S = self.sigma_p[b]
        cen = self.centric if index is None else self.centric[index]

        fn = 0.5 * (fp + fm)
        df = fp - fm
        vn = S ** 2
        vd = vp + vm

        # --- amplitude grid ------------------------------------------------
        # Rayleigh prior (rms S) has mode S/sqrt(2) and a spread of about 0.6 S;
        # the native observation constrains the amplitude to about sqrt(vn).
        # Their product is a Gaussian whose mode and width set the grid.
        m0 = S / np.sqrt(2.0)
        s0 = 0.6 * S
        prec = 1.0 / s0 ** 2 + 1.0 / vn
        centre = (m0 / s0 ** 2 + fn / vn) / prec
        width = 1.0 / np.sqrt(prec)
        lo = np.maximum(centre - 6.0 * width, 0.0)
        hi = centre + 6.0 * width
        half = np.maximum(0.5 * (hi - lo), 1e-6)
        mid = 0.5 * (hi + lo)
        A = mid[:, None] + half[:, None] * _GL[None, :]
        logw = np.log(half[:, None] * _GLW[None, :])
        log_prior = np.log(2.0 * A / S[:, None] ** 2) - (A / S[:, None]) ** 2
        log_prior_h = -np.abs(fh) ** 2 / sigh - np.log(np.pi * sigh)
        base = -0.5 * np.log(2.0 * np.pi * vn) - 0.5 * np.log(2.0 * np.pi * vd)
        fh2 = np.abs(fh) ** 2
        fhn2 = np.abs(fhn) ** 2

        out = np.empty(n)
        if marginal:
            log_marg = np.full((n, max(self.n_alpha, 2)), -np.inf)
        if not cen.all():
            # --- acentric: average over the native phase -------------------
            al = (np.arange(self.n_alpha) + 0.5) * (2.0 * np.pi / self.n_alpha)
            rot = np.exp(-1j * al)
            # |A exp(i a) + F_H|^2 = A^2 + |F_H|^2 + 2 A Re(F_H exp(-i a))
            cross_p = np.real(fh[:, None] * rot[None, :])          # (n, na)
            cross_m = np.real(fhn[:, None] * rot[None, :])
            A2 = A ** 2
            cp = np.sqrt(np.maximum(A2[:, :, None]
                                    + fh2[:, None, None]
                                    + 2.0 * A[:, :, None] * cross_p[:, None, :], 0.0))
            cm = np.sqrt(np.maximum(A2[:, :, None]
                                    + fhn2[:, None, None]
                                    + 2.0 * A[:, :, None] * cross_m[:, None, :], 0.0))
            cn = 0.5 * (cp + cm)
            dfc = cp - cm
            lp = (-(fn[:, None, None] - cn) ** 2 / (2.0 * vn[:, None, None])
                  - (df[:, None, None] - dfc) ** 2 / (2.0 * vd[:, None, None]))
            lp += (base[:, None, None] + log_prior_h[:, None, None]
                   + log_prior[:, :, None] + logw[:, :, None]
                   - np.log(self.n_alpha))
            flat = lp.reshape(n, -1)
            out[:] = logsumexp(flat, axis=1)
            if marginal:
                log_marg[:] = logsumexp(lp, axis=1)
        if cen.any():
            # --- centric: two discrete native phases, no anomalous difference
            sel = np.where(cen)[0]
            phi = self.allowed if index is None else self.allowed[index]
            hrot = fh[sel] * np.exp(-1j * phi[sel])
            As = A[sel]
            # the two signs are kept apart so that the marginal can report the
            # probability of each allowed phase
            common = (-0.5 * np.log(2.0 * np.pi * vn[sel, None])
                      + log_prior_h[sel, None] + log_prior[sel] + logw[sel]
                      - np.log(2.0))
            terms = []
            for sgn in (1.0, -1.0):
                mag = np.abs(sgn * As + hrot[:, None])
                terms.append(-(fn[sel, None] - mag) ** 2 / (2.0 * vn[sel, None])
                             + common)
            qmin = np.minimum(terms[0], terms[1])
            lp = -qmin + np.log(np.exp(-(terms[0] - qmin))
                                + np.exp(-(terms[1] - qmin)))
            out[sel] = logsumexp(lp, axis=1)
            if marginal:
                sign = np.stack([logsumexp(t, axis=1) for t in terms], axis=1)
                log_marg[sel, :2] = sign - logsumexp(sign, axis=1)[:, None]
        return (out, log_marg) if marginal else out

    def phase_moments(self, data, sr, si, sigh, index=None):
        """Native phase and figure of merit implied by the full likelihood.

        The joint probability over the native amplitude and phase is
        marginalised over the amplitude, and the first circular moment of the
        remaining phase distribution gives the phase and the FOM.  A centric
        reflection has only two allowed phases, so its FOM is the difference of
        the two probabilities.

        Returns ``(phase, fom)`` in radians and ``[0, 1]``.
        """
        vp, vm = data.sigfp ** 2, data.sigfm ** 2
        if index is not None:
            vp, vm = vp[index], vm[index]
        _, log_marg = self.core(sr, si, data.fp, data.fm, vp, vm,
                                sigh[self.bin if index is None else index],
                                index=index, marginal=True)
        n = len(sr)
        al = (np.arange(self.n_alpha) + 0.5) * (2.0 * np.pi / self.n_alpha)
        phase = np.zeros(n)
        fom = np.zeros(n)
        acen = ~(self.centric if index is None else self.centric[index])
        if acen.any():
            w = np.exp(log_marg[acen])
            w = w / np.maximum(w.sum(axis=1)[:, None], 1e-300)
            z = np.sum(w * np.exp(1j * al)[None, :], axis=1)
            phase[acen] = np.angle(z) % (2.0 * np.pi)
            fom[acen] = np.abs(z)
        if not acen.all():
            sel = np.where(~acen)[0]
            phi = self.allowed if index is None else self.allowed[index]
            p = np.exp(log_marg[sel, :2])
            p = p / np.maximum(p.sum(axis=1)[:, None], 1e-300)
            phase[sel] = phi[sel] % (2.0 * np.pi)
            fom[sel] = np.abs(p[:, 0] - p[:, 1])
        return phase, np.clip(fom, 0.0, 1.0)

    def default_sigh(self, fraction=0.1):
        """Starting value of the substructure prior variance, per bin.

        The expected intensity of the heavy-atom contribution is a fraction of
        the Wilson intensity; ``fraction`` is that ratio, so the prior standard
        deviation is ``sqrt(fraction) * rms|F|``.
        """
        return fraction * self.sigma_p ** 2

    def sigh_bounds(self, low=0.02, high=0.5):
        """Bounds for the substructure prior variance when it is refined.

        ``Sigma_H`` may not be refined freely.  The prior
        ``-|F_H|^2/Sigma_H - log(pi Sigma_H)`` is minimised at
        ``Sigma_H = <|F_H|^2>``, so a free ``Sigma_H`` makes the prior penalty
        scale-free: growing the substructure by a factor k grows the penalty
        by k^2 but grows its own optimum by k^2 as well.  The likelihood then
        acquires a flat direction and the occupancies run away (observed on
        insulin, whose anomalous signal is at the noise level).  Keeping
        ``Sigma_H`` inside a fraction of the Wilson intensity removes the
        degeneracy while still letting it adapt to the data.

        Returns ``(low, high)`` arrays of ``log(Sigma_H)`` per bin.
        """
        return (np.log(low * self.sigma_p ** 2),
                np.log(high * self.sigma_p ** 2))

    def dataset_logp(self, data, sr, si, sigh, index=None):
        fp, fm = data.fp, data.fm
        vp, vm = data.sigfp ** 2, data.sigfm ** 2
        sg = sigh[self.bin]
        if index is not None:
            fp, fm = fp[index], fm[index]
            vp, vm = vp[index], vm[index]
            sr, si = sr[index], si[index]
            sg = sg[index]
        return self.core(sr, si, fp, fm, vp, vm, sg, index)

    def logp(self, sr, si, sigh):
        return float(np.sum(self.logp_given(sr, si, sigh)))

    def logp_given(self, sr, si, sigh, index=None):
        """Log likelihood with the structure factors supplied (no recomputation)."""
        return float(np.sum([np.sum(self.dataset_logp(d, sr, si, sigh, index))
                             for d in self.datasets]))

    # ----------------------------------------------------------- gradient
    def gradient(self, sr, si, dsr=None, dsi=None, sigh=None, atoms=True, bins=True,
                 step=1e-4):
        """Gradient of :meth:`nll` with respect to every parameter.

        ``core`` is a smooth function of the two complex structure factors
        ``S_r`` and ``S_i`` and each reflection is independent, so the
        derivative with respect to the structure factors is obtained by a
        central difference on the real and on the imaginary part separately;
        the chain rule to the atomic parameters is then exact.  ``S_r`` and
        ``S_i`` carry four real degrees of freedom (``F_H`` and its Friedel
        mate vary independently in that representation), so all four
        directions are needed.
        """
        nat = len(self.sub.occ)
        g = np.zeros(5 * nat + self.bins)
        if atoms:
            gr = np.zeros(len(sr), dtype=complex)
            gi = np.zeros(len(si), dtype=complex)
            for data in self.datasets:
                vp, vm = data.sigfp ** 2, data.sigfm ** 2
                sg = sigh[self.bin]

                def f(a, b):
                    return self.core(a, b, data.fp, data.fm, vp, vm, sg)

                gr += ((f(sr + step, si) - f(sr - step, si))
                       + 1j * (f(sr + 1j * step, si) - f(sr - 1j * step, si))) \
                    / (2 * step)
                gi += ((f(sr, si + step) - f(sr, si - step))
                       + 1j * (f(sr, si + 1j * step) - f(sr, si - 1j * step))) \
                    / (2 * step)
            g[:5 * nat] = -(np.real(gr[:, None] * np.conj(dsr)).sum(axis=0)
                            + np.real(gi[:, None] * np.conj(dsi)).sum(axis=0))
        if bins:
            for i in range(self.bins):
                mask = self.bin == i
                if not mask.any():
                    continue
                h = 1e-5 * max(1.0, abs(sigh[i]))

                def total(v):
                    vv = sigh.copy()
                    vv[i] = v
                    return self.logp_given(sr, si, vv, index=mask)

                g[5 * nat + i] = -sigh[i] * (total(sigh[i] + h)
                                             - total(sigh[i] - h)) / (2 * h)
        return g

    # ------------------------------------------------------------ parameters
    def pack(self, xyz, occ, b, sigh):
        return np.concatenate([np.ravel(xyz), occ, b, np.log(sigh)])

    def unpack(self, p):
        nat = (len(p) - self.bins) // 5
        xyz = p[:3 * nat].reshape(nat, 3)
        occ = p[3 * nat:4 * nat]
        b = p[4 * nat:5 * nat]
        sigh = np.exp(p[5 * nat:])
        return xyz, occ, b, sigh

    def nll(self, p):
        xyz, occ, b, sigh = self.unpack(p)
        sr, si = self.sub.factors_split(xyz, occ, b)
        return -np.sum(self.logp(sr, si, sigh))


def _centric_phases(sub):
    """Which reflections are centric, and their allowed phase."""
    centric = np.zeros(len(sub.hkl), bool)
    allowed = np.zeros(len(sub.hkl))
    for i, h in enumerate(sub.hkl.tolist()):
        for rot, tran in sub.ops:
            if list(np.dot(rot.T, h).astype(int)) == [-v for v in h]:
                centric[i] = True
                allowed[i] = (-np.pi * np.dot(h, tran)) % np.pi
                break
    return centric, allowed
