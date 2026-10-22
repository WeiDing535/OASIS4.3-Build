"""Optional residue target fitted to unrelated density patches only.

Class balancing implements an equal amino-acid prior. The output retains the
learned logit magnitude instead of forcing every site to unit class variance.
"""
import json
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.spatial.transform import Rotation
from .residue_profile import ALPHABET, features


class ResidueLogisticProfile:
    def __init__(self, offsets, center, scale, coefficients, intercept, metadata):
        self.offsets, self.center, self.scale = map(np.asarray, [offsets, center, scale])
        self.coefficients, self.intercept = map(np.asarray, [coefficients, intercept])
        self.metadata = metadata
        size = len(self.offsets)
        if (self.offsets.shape != (size, 3) or self.center.shape != (size,) or
                self.scale.shape != (size,) or self.coefficients.shape != (size, 20) or
                self.intercept.shape != (20,) or np.any(self.scale <= 0) or
                any(not np.isfinite(v).all() for v in [self.offsets, self.center, self.scale,
                                                       self.coefficients, self.intercept])):
            raise ValueError('Invalid residue logistic profile')

    def scores_from_features(self, values):
        values = np.asarray(values)
        logits = ((values-self.center)/self.scale)@self.coefficients+self.intercept
        return logits-logits.mean(axis=-1, keepdims=True)

    def scores(self, density, xyz):
        return self.scores_from_features(features(density, xyz, self.offsets))

    def save(self, path):
        np.savez_compressed(path, offsets=self.offsets, center=self.center, scale=self.scale,
                            coefficients=self.coefficients, intercept=self.intercept,
                            metadata=json.dumps(self.metadata))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as z:
            return cls(z['offsets'], z['center'], z['scale'], z['coefficients'], z['intercept'],
                       json.loads(str(z['metadata'])))


def perturb_backbones(xyz, random, translation_sigma=.25, rotation_sigma_degrees=5.):
    """Rigid local-frame errors for augmentation; no target-derived parameters."""
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 3 or xyz.shape[1:] != (4, 3):
        raise ValueError('Expected N/CA/C/O local backbone coordinates')
    translations = random.normal(scale=translation_sigma, size=(len(xyz), 1, 3))
    vectors = random.normal(scale=np.deg2rad(rotation_sigma_degrees), size=(len(xyz), 3))
    matrices = Rotation.from_rotvec(vectors).as_matrix()
    ca = xyz[:, 1:2, :]
    return np.einsum('nij,nkj->nki', matrices, xyz-ca)+ca+translations


def fit_residue_logistic(values, labels, offsets, metadata=None, regularization=.1, max_iterations=300):
    """Fit class-balanced softmax; returned metadata audits all fit parameters."""
    values = np.asarray(values, dtype=float); labels = np.asarray(labels, dtype=int)
    if (values.ndim != 2 or labels.shape != (len(values),) or not np.isfinite(values).all() or
            not np.isfinite(regularization) or regularization <= 0 or
            labels.min(initial=0) < 0 or labels.max(initial=0) >= 20):
        raise ValueError('Invalid residue training data or regularization')
    counts = np.bincount(labels, minlength=20)
    if np.any(counts == 0):
        raise ValueError('All twenty residue classes are required in training')
    center = values.mean(axis=0)
    scale = np.maximum(values.std(axis=0), .05)
    x = (values-center)/scale
    examples, dimensions = x.shape
    weights = 1./(20*counts[labels])
    def objective(vector):
        coefficients = vector[:dimensions*20].reshape(dimensions, 20)
        intercept = vector[dimensions*20:]
        logits = x@coefficients+intercept
        log_norm = logsumexp(logits, axis=1)
        loss = np.dot(weights, log_norm-logits[np.arange(examples), labels])
        loss += .5*regularization*np.square(coefficients).sum()
        residual = np.exp(logits-log_norm[:, None])
        residual[np.arange(examples), labels] -= 1.
        residual *= weights[:, None]
        gradient = np.concatenate([(x.T@residual+regularization*coefficients).ravel(), residual.sum(axis=0)])
        return float(loss), gradient
    result = minimize(objective, np.zeros(dimensions*20+20), jac=True, method='L-BFGS-B',
                      options={'maxiter': max_iterations, 'ftol': 1e-10, 'gtol': 1e-6})
    fitted = dict(metadata or {})
    fitted.update(alphabet=ALPHABET, method='class-balanced ridge multinomial logistic density patches',
                  regularization=regularization, training_examples=examples,
                  class_counts={letter: int(counts[i]) for i, letter in enumerate(ALPHABET)},
                  optimizer_success=bool(result.success), optimizer_message=str(result.message),
                  iterations=int(result.nit), objective=float(result.fun),
                  score_scale='centered raw logits; no per-site variance normalization',
                  scope='Trained on supplied unrelated examples; no target-reference input in this API.')
    return ResidueLogisticProfile(offsets, center, scale,
                                   result.x[:dimensions*20].reshape(dimensions, 20), result.x[dimensions*20:], fitted)
