"""Proper ZINB likelihood, deterministic histogram fitting and scalar scaling.

Unlike the inspected TRANSIT objective (DeJesus et al. 2015), the zero mass includes (1-pi).
The empirical catalogue defines the population of zeros; it is not ground truth.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import nbinom

LOG_BOUNDS = ((-18., 25.), (-14., 18.))  # NB mean and size, numerical safeguards
START_SIZES = (.2, 2., 20.)


def profile(values, frequencies, zeros, log_parameters):
    """Mean negative log likelihood with pi maximized analytically."""
    mu, size = np.exp(log_parameters)
    logp = -np.log1p(mu / size)
    logq = -np.log1p(size / mu)
    logzero = size * logp
    nbzero = np.exp(logzero)
    count = zeros + frequencies.sum()
    fraction = zeros / count
    pi = float(np.clip((fraction - nbzero) / (-np.expm1(logzero)), 0., 1.-1e-12))
    logmixzero = np.logaddexp(np.log(pi) if pi else -np.inf, np.log1p(-pi) + logzero)
    logs = (gammaln(values + size) - gammaln(size) - gammaln(values + 1)
            + size * logp + values * logq)
    loss = -(zeros * logmixzero + np.dot(frequencies, np.log1p(-pi) + logs)) / count
    return float(loss), pi


def fit(positive, universe):
    positive = np.asarray(positive, dtype=float)
    if (not len(positive) or universe < len(positive) or np.any(~np.isfinite(positive))
            or np.any(positive <= 0) or np.any(positive != np.floor(positive))):
        raise ValueError("ZINB needs positive integer support and a valid catalogue size")
    values, frequencies = np.unique(positive, return_counts=True)
    zeros = int(universe - len(positive))
    starts = [np.log([positive.mean(), size]) for size in START_SIZES]
    candidates = [minimize(lambda v: profile(values, frequencies, zeros, v)[0], start,
                          method="L-BFGS-B", bounds=LOG_BOUNDS,
                          options=dict(maxiter=1000, ftol=1e-12, gtol=1e-7)) for start in starts]
    valid = [r for r in candidates if r.success and np.isfinite(r.fun)]
    if not valid:
        raise ValueError("ZINB optimizer failed for all deterministic starts; no fallback scaling")
    best = min(valid, key=lambda r: r.fun)
    mu, size = np.exp(best.x)
    loss, pi = profile(values, frequencies, zeros, best.x)
    p = size / (size + mu)
    fitted_zero = pi + (1-pi)*np.exp(size*np.log(p))
    # Discrete CDF discrepancy checked on both sides of each occupied count.
    empirical = (zeros + np.cumsum(frequencies)) / universe
    before = empirical - frequencies / universe
    discrepancy = max(abs(fitted_zero - zeros/universe),
        np.max(abs(pi + (1-pi)*nbinom.cdf(values, size, p) - empirical)),
        np.max(abs(pi + (1-pi)*nbinom.cdf(values-1, size, p) - before)))
    boundary = any(min(abs(v-lo), abs(v-hi)) < 1e-3 for v, (lo, hi) in zip(best.x, LOG_BOUNDS))
    return dict(factor=float(1/mu), zinb_pi=pi, zinb_nb_mean=float(mu),
        zinb_size=float(size), zinb_probability=float(p), zinb_observed_zero_fraction=zeros/universe,
        zinb_fitted_zero_fraction=float(fitted_zero), zinb_cdf_max_discrepancy=float(discrepancy),
        zinb_mean_negative_log_likelihood=loss, zinb_converged=True,
        zinb_successful_starts=len(valid), zinb_total_starts=len(starts),
        zinb_iterations=int(best.nit), zinb_boundary=bool(boundary or pi == 0),
        zinb_warning="Parameter boundary; interpret fit cautiously" if boundary or pi == 0 else "",
        zinb_optimizer_message=str(best.message))
