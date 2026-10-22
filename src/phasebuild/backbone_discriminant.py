"""Experimental backbone/decoy patch discrimination, independent of Buccaneer.

Training uses an unrelated model and controlled coordinate perturbations. Scores
are classifier logits for that synthetic task, not calibrated probabilities of
experimental structure correctness. This module is opt-in and does not alter
DensityProfile, tracing, sequence placement, or the default pipeline.
"""
import json
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from .profile import patches


def patch_offsets(radius=3., spacing=.8):
    grid = np.arange(-int(radius / spacing), int(radius / spacing) + 1) * spacing
    offsets = np.array(np.meshgrid(grid, grid, grid, indexing='ij')).reshape(3, -1).T
    return offsets[np.linalg.norm(offsets, axis=1) <= radius]


def features(density, xyz, offsets):
    values = np.clip(patches(density, xyz, offsets), -6., 12.)
    mean = values.mean(axis=-1, keepdims=True)
    std = np.maximum(values.std(axis=-1, keepdims=True), .15)
    contrast = np.clip((values - mean) / std, -5., 5.)
    atoms = np.clip(density.sample(xyz), -6., 12.)
    # Raw density retains the signal lost by local correlation normalization;
    # normalized first and second moments distinguish local spatial patterns.
    summary = np.concatenate([mean, std, atoms,
                              atoms.min(axis=-1, keepdims=True),
                              atoms.mean(axis=-1, keepdims=True)], axis=-1)
    return np.concatenate([contrast, contrast ** 2, .25 * values, summary], axis=-1)


class BackboneDiscriminant:
    def __init__(self, offsets, center, scale, coefficient, intercept, metadata):
        self.offsets = np.asarray(offsets, dtype=float)
        self.center = np.asarray(center, dtype=float)
        self.scale = np.asarray(scale, dtype=float)
        self.coefficient = np.asarray(coefficient, dtype=float)
        self.intercept = float(intercept)
        self.metadata = metadata
        expected = 3 * len(self.offsets) + 8
        if (self.offsets.ndim != 2 or self.offsets.shape[1:] != (3,) or
                any(a.shape != (expected,) for a in
                    [self.center, self.scale, self.coefficient]) or
                any(not np.isfinite(a).all() for a in
                    [self.offsets, self.center, self.scale, self.coefficient]) or
                np.any(self.scale <= 0) or not np.isfinite(self.intercept)):
            raise ValueError('Invalid discriminant arrays')

    def score_features(self, values):
        return ((values - self.center) / self.scale) @ self.coefficient + self.intercept

    def score(self, density, xyz, block_size=256):
        xyz = np.asarray(xyz)
        if xyz.shape[-2:] != (4, 3):
            raise ValueError('Expected N/CA/C/O groups')
        shape = xyz.shape[:-2]
        flat = xyz.reshape(-1, 4, 3)
        scores = [self.score_features(features(density, flat[i:i + block_size], self.offsets))
                  for i in range(0, len(flat), block_size)]
        return np.concatenate(scores).reshape(shape) if scores else np.empty(shape)

    def save(self, path):
        np.savez_compressed(path, offsets=self.offsets, center=self.center,
                            scale=self.scale, coefficient=self.coefficient,
                            intercept=self.intercept, metadata=json.dumps(self.metadata))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            return cls(data['offsets'], data['center'], data['scale'],
                       data['coefficient'], float(data['intercept']),
                       json.loads(str(data['metadata'])))


def fit_discriminant(values, labels, offsets, metadata, regularization=1., maxiter=180):
    values, labels = np.asarray(values, dtype=float), np.asarray(labels, dtype=float)
    if values.ndim != 2 or labels.shape != (len(values),) or set(np.unique(labels)) != {0., 1.}:
        raise ValueError('Training requires finite feature rows and both binary classes')
    if not np.isfinite(values).all() or regularization <= 0:
        raise ValueError('Training requires finite features and positive regularization')
    center, scale = values.mean(axis=0), np.maximum(values.std(axis=0), .05)
    standardized = (values - center) / scale
    # Balance the two synthetic classes; experimental priors are unknown.
    weights = np.where(labels > .5, .5 / labels.mean(), .5 / (1. - labels.mean()))
    normalizer = weights.sum()
    def loss(params):
        logits = standardized @ params[:-1] + params[-1]
        penalty = .5 * regularization * np.sum(params[:-1] ** 2) / normalizer
        value = np.sum(weights * (np.logaddexp(0., logits) - labels * logits)) / normalizer + penalty
        residual = weights * (expit(logits) - labels) / normalizer
        gradient = np.r_[standardized.T @ residual + regularization * params[:-1] / normalizer,
                         residual.sum()]
        return float(value), gradient
    optimized = minimize(loss, np.zeros(values.shape[1] + 1), jac=True, method='L-BFGS-B',
                         options={'maxiter': maxiter, 'ftol': 1e-10, 'gtol': 1e-6})
    info = dict(metadata, optimization={'success': bool(optimized.success),
                'message': str(optimized.message), 'iterations': int(optimized.nit),
                'loss': float(optimized.fun), 'regularization': float(regularization)})
    return BackboneDiscriminant(offsets, center, scale, optimized.x[:-1], optimized.x[-1], info)


def discrimination_summary(labels, scores, thresholds=None):
    from scipy.stats import rankdata
    labels, scores = np.asarray(labels, dtype=bool), np.asarray(scores, dtype=float)
    positive, negative = scores[labels], scores[~labels]
    if not len(positive) or not len(negative):
        raise ValueError('Both classes required for discrimination summary')
    ranks = rankdata(scores)
    auc = (ranks[labels].sum() - len(positive) * (len(positive) + 1) / 2.) / (len(positive) * len(negative))
    if thresholds is None:
        thresholds = {str(rate): float(np.quantile(negative, 1. - rate)) for rate in [.01, .05, .10, .20]}
    return {'synthetic_auc': float(auc), 'positives': len(positive), 'negatives': len(negative),
            'positive_logit_quantiles': np.quantile(positive, [.05, .5, .95]).tolist(),
            'negative_logit_quantiles': np.quantile(negative, [.05, .5, .95]).tolist(),
            'thresholds': {key: {'logit': float(threshold),
                'positive_retention': float(np.mean(positive >= threshold)),
                'negative_retention': float(np.mean(negative >= threshold))}
                for key, threshold in thresholds.items()}}


def supported_segments(scores, threshold, minimum_length=4, bridge=1):
    """Retain supported runs and bounded internal gaps; no terminal weak tails."""
    keep = np.asarray(scores) >= threshold
    start = None
    for i, value in enumerate(keep):
        if not value and start is None:
            start = i
        elif value and start is not None:
            if start > 0 and i - start <= bridge:
                keep[start:i] = True
            start = None
    edges = np.diff(np.r_[False, keep, False].astype(int))
    return [(int(a), int(b)) for a, b in zip(np.where(edges == 1)[0], np.where(edges == -1)[0])
            if b - a >= minimum_length]
