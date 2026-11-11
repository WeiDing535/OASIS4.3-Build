"""Fixed-conditional-shape, variational two-state P+ research approximation.

Observation-conditional branch shapes stay frozen; only their masses change.
This is not an exact continuous posterior or a joint Bijvoet likelihood.
No scale fitting, external phases or protein models are accepted.
"""
from dataclasses import dataclass, field
import hashlib
import json
import numpy as np
from scipy.special import expit, logsumexp
from oasis_py.clipper_likelihood import reflection_classes


def _identity(data):
    digest = hashlib.sha256(np.ascontiguousarray(data.hkl, dtype=np.int32).tobytes())
    digest.update(json.dumps({'cell': list(data.cell.parameters),
        'symmetry': [op.triplet() for op in data.spacegroup.operations()]},
        sort_keys=True).encode())
    return digest.hexdigest()


def _frozen(value, dtype):
    out = np.array(value, dtype=dtype, copy=True)
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class ConditionalBranches:
    """Normalized log masses and conditional moments, each shape (n, 2).

    Column 0 is the positive-sine semicircle relative to reference. Factory
    summaries carry crystal/HKL identity. Manually constructed summaries may
    omit identity; their caller is responsible for matching reflection order.
    """
    log_mass: np.ndarray
    first: np.ndarray
    second: np.ndarray
    reference: np.ndarray | None = None
    centric: np.ndarray | None = None
    data_identity: str | None = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        lm = np.array(self.log_mass, dtype=float, copy=True)
        first = np.asarray(self.first, dtype=complex)
        second = np.asarray(self.second, dtype=complex)
        if (lm.ndim != 2 or lm.shape[1] != 2 or not len(lm)
                or first.shape != lm.shape or second.shape != lm.shape):
            raise ValueError('Branch arrays must have nonempty shape (n, 2)')
        if (np.isnan(lm).any() or np.isposinf(lm).any()
                or np.any(np.all(np.isneginf(lm), axis=1))
                or not np.isfinite(first).all() or not np.isfinite(second).all()):
            raise ValueError('Invalid branch masses or moments')
        if np.any(abs(first) > 1 + 1e-10) or np.any(abs(second) > 1 + 1e-10):
            raise ValueError('Conditional moments exceed the unit disk')
        if np.any(abs(second-first**2) > 1-abs(first)**2+1e-9):
            raise ValueError('Inconsistent first and second conditional moments')
        lm -= logsumexp(lm, axis=1)[:, None]
        object.__setattr__(self, 'log_mass', _frozen(lm, float))
        object.__setattr__(self, 'first', _frozen(first, complex))
        object.__setattr__(self, 'second', _frozen(second, complex))
        if self.reference is not None:
            ref = np.asarray(self.reference, dtype=float)
            if ref.shape != (len(lm),) or not np.isfinite(ref).all():
                raise ValueError('Invalid branch references')
            object.__setattr__(self, 'reference', _frozen(ref % (2*np.pi), float))
        if self.centric is not None:
            cen = np.asarray(self.centric)
            if cen.shape != (len(lm),) or cen.dtype.kind != 'b':
                raise ValueError('Invalid branch centric mask')
            object.__setattr__(self, 'centric', _frozen(cen, bool))
        if self.data_identity is not None and not isinstance(self.data_identity, str):
            raise ValueError('Invalid branch data identity')
        object.__setattr__(self, 'metadata', dict(self.metadata))


@dataclass(frozen=True)
class ConditionalPplusResult:
    phase: np.ndarray
    fom: np.ndarray
    probability: np.ndarray
    first: np.ndarray
    second: np.ndarray
    history: list
    metadata: dict


def prepare_conditional_branches(data, hl, reference, samples=720):
    """Integrate frozen HL on paired branch-aligned midpoint grids.

    ``samples`` remains the total number of acentric nodes. Each half circle
    uses samples//2 midpoint nodes at reference +/- (j+.5)*pi/(samples//2).
    Equal-weight half-circle integrals are normalized together. This avoids
    a discontinuous mask on an absolute grid and preserves mirror-symmetric
    observation masses up to roundoff. Centric observations still use their
    two allowed states, with the existing cosine boundary tie break.
    """
    n = len(data.f)
    coefficients = np.asarray(hl, dtype=float)
    ref = np.asarray(reference, dtype=float)
    if (coefficients.shape != (n, 4) or ref.shape != (n,)
            or not np.isfinite(coefficients).all() or not np.isfinite(ref).all()):
        raise ValueError('Invalid observation HL or branch reference')
    if type(samples) is not int or samples < 72 or samples % 2:
        raise ValueError('Use an even angular grid of at least 72 samples')
    _, centric, allowed = reflection_classes(data)
    lm = np.empty((n, 2))
    first = np.empty((n, 2), dtype=complex)
    second = np.empty((n, 2), dtype=complex)
    ref = ref % (2*np.pi)
    half_samples = samples//2
    offset = (np.arange(half_samples)+.5)*(np.pi/half_samples)
    positive_basis = np.array([np.cos(offset), np.sin(offset),
                               np.cos(2*offset), np.sin(2*offset)])
    negative_basis = positive_basis * np.array([1., -1., 1., -1.])[:, None]
    positive_e1 = np.exp(1j*offset)
    positive_e2 = np.exp(2j*offset)
    indices = np.flatnonzero(~centric)
    for start in range(0, len(indices), 128):
        ix = indices[start:start+128]
        angle = ref[ix]
        rotation1, rotation2 = np.exp(1j*angle), np.exp(2j*angle)
        # Express the same absolute HL in coordinates relative to reference.
        relative = np.empty_like(coefficients[ix])
        for harmonic, columns in [(1, (0, 1)), (2, (2, 3))]:
            value = (coefficients[ix, columns[0]]
                     +1j*coefficients[ix, columns[1]])*np.exp(-1j*harmonic*angle)
            relative[:, columns[0]], relative[:, columns[1]] = value.real, value.imag
        for branch, (basis, e1, e2) in enumerate([
                (positive_basis, positive_e1, positive_e2),
                (negative_basis, positive_e1.conj(), positive_e2.conj())]):
            lp = relative @ basis
            if not np.isfinite(lp).all():
                raise ValueError('Observation HL overflows the angular grid')
            total = logsumexp(lp, axis=1)
            weights = np.exp(lp-total[:, None])
            # Both semicircles have the same node weights; the common
            # factor pi/half_samples cancels when normalizing their masses.
            lm[ix, branch] = total
            first[ix, branch] = rotation1*(weights @ e1)
            second[ix, branch] = rotation2*(weights @ e2)
        lm[ix] -= logsumexp(lm[ix], axis=1)[:, None]
    boundary_rows = 0
    for i in np.flatnonzero(centric):
        states = allowed[i]+np.array([0., np.pi])
        cbasis = np.array([np.cos(states), np.sin(states),
                           np.cos(2*states), np.sin(2*states)])
        lp = coefficients[i] @ cbasis
        if not np.isfinite(lp).all():
            raise ValueError('Observation HL overflows the centric states')
        sine = np.sin(states-ref[i])
        boundary = abs(sine) <= 1e-12
        boundary_rows += int(boundary.any())
        positive = np.where(boundary, np.cos(states-ref[i]) > 0, sine > 0)
        if int(positive.sum()) != 1:
            raise ValueError('Centric branch assignment is degenerate')
        order = [int(np.flatnonzero(positive)[0]), int(np.flatnonzero(~positive)[0])]
        lm[i] = lp[order]-logsumexp(lp)
        state_first = np.exp(1j*allowed[i])*np.array([1., -1.])
        first[i] = state_first[order]
        second[i] = np.exp(2j*allowed[i])
    return ConditionalBranches(lm, first, second, ref, centric, _identity(data), {
        'method': 'fixed_observation_conditional_branch_summaries',
        'angular_samples': samples,
        'angular_grid': 'paired branch-aligned semicircle midpoint grids; samples total nodes',
        'acentric_boundary_rule': 'no boundary nodes; equal-weight paired half circles',
        'centric_boundary_rule': 'cosine sign tie break',
        'centric_boundary_rows': boundary_rows,
        'observation_HL_modified': False,
        'quadrature_scope': 'Fixed grid; convergence must be checked separately'})


def iterate_conditional_pplus(data, branches, sim_field, relations, cycles=5,
                              triplet_weight=1., damping=1.):
    """Update P+ masses with frozen observation-conditional phase shapes.

    The first cycle uses observations only. Later cycles use
    logit(P+) = log(m+/m-) + Re(field*conjugate(z+ - z-)).
    This is a two-state variational update, not the exact continuous posterior.
    Damping mixes branch probabilities; observations never accumulate.
    """
    if not isinstance(branches, ConditionalBranches):
        raise ValueError('ConditionalBranches are required')
    n = len(data.f)
    if branches.log_mass.shape != (n, 2):
        raise ValueError('Branch/reflection counts differ')
    if branches.data_identity is not None and branches.data_identity != _identity(data):
        raise ValueError('Branch summaries do not match reflection/crystal identity')
    if branches.centric is not None and not np.array_equal(branches.centric, data.centric):
        raise ValueError('Branch centric classes differ from reflections')
    if type(cycles) is not int or cycles < 1:
        raise ValueError('P+ cycles must be a positive integer')
    if (not np.isfinite(triplet_weight) or triplet_weight < 0
            or not np.isfinite(damping) or not 0 < damping <= 1):
        raise ValueError('Invalid triplet weight or probability damping')
    sim = np.asarray(sim_field, dtype=complex)
    if sim.shape != (n,) or not np.isfinite(sim).all():
        raise ValueError('Invalid partial-structure field')
    if not callable(getattr(relations, 'field', None)):
        raise ValueError('A relation-field provider is required')
    if getattr(relations, 'nreflections', n) != n:
        raise ValueError('Relation/reflection counts differ')
    _, centric, allowed = reflection_classes(data)
    axis = np.exp(1j*allowed[centric])
    if (np.any(abs(np.imag(branches.first[centric]/axis[:, None])) > 1e-9)
            or np.any(abs(branches.second[centric]-axis[:, None]**2) > 1e-9)):
        raise ValueError('Conditional moments violate allowed centric states')
    base_logit = branches.log_mass[:, 0]-branches.log_mass[:, 1]
    probability = expit(base_logit)
    history = []
    for cycle in range(cycles):
        previous = probability.copy()
        if cycle:
            graph_field = np.asarray(relations.field(phase, fom), dtype=complex)
            if graph_field.shape != (n,) or not np.isfinite(graph_field).all():
                raise ValueError('Invalid relation field')
            current_field = sim+triplet_weight*graph_field
            if not np.isfinite(current_field).all():
                raise ValueError('Combined phase field is not finite')
            score = np.real(current_field*np.conj(
                branches.first[:, 0]-branches.first[:, 1]))
            proposed = expit(base_logit+score)
            probability = (1-damping)*probability+damping*proposed
        first = probability*branches.first[:, 0]+(1-probability)*branches.first[:, 1]
        second = probability*branches.second[:, 0]+(1-probability)*branches.second[:, 1]
        fom = np.clip(abs(first), 0, 1)
        phase = np.angle(first) % (2*np.pi)
        signed = np.real(first[centric]/axis)
        phase[centric] = (allowed[centric]+np.where(signed < 0, np.pi, 0)) % (2*np.pi)
        history.append({'cycle': cycle+1, 'mean_fom': float(fom.mean()),
                        'mean_pplus': float(probability.mean()),
                        'maximum_probability_change': float(abs(probability-previous).max())})
    return ConditionalPplusResult(
        _frozen(phase, float), _frozen(fom, float), _frozen(probability, float),
        _frozen(first, complex), _frozen(second, complex), history, {
        'method': 'conditional_shape_variational_two_state_Pplus_research',
        'approximation': 'Fixed observation-conditional shapes; not exact continuous posterior',
        'cycles': cycles, 'triplet_weight': float(triplet_weight),
        'probability_damping': float(damping),
        'observation_likelihood_accumulated': False,
        'first_cycle': 'Observation branch masses and moments only',
        'posterior_HL_exported': False})

