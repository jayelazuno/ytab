"""Beta-geometric count correction (BGC), adapted from TRANSIT's BetaGeomNorm with audited parameter estimators.

Primary evaluation uses the exact beta-geometric CDF; a seeded simulation engine
retains the upstream approximation for sensitivity checks. Zeros are preserved.
"""
import hashlib

import numpy as np
from scipy import sparse
from scipy.special import betaln
from scipy.stats import geom, trim_mean

from .ttr import ttr_parameters


def fit_parameters(positive, universe):
    positive = np.asarray(positive, dtype=float)
    if (positive.ndim != 1 or not len(positive) or universe < len(positive) or
            np.any(~np.isfinite(positive)) or np.any(positive <= 0)):
        raise ValueError("BGC requires nonempty finite positive counts and a common universe")
    ordered = np.r_[np.zeros(universe-len(positive)), np.sort(positive)]
    shifted_mean = float(trim_mean(ordered+1, .001))
    rho = max(1 / shifted_mean, .0001)
    logs = np.log1p(-np.arange(universe, dtype=float)/universe)
    denominator = float(np.dot(ordered, logs))
    if denominator == 0 or not 0 < rho < 1:
        raise ValueError("BGC fit is undefined for these counts; no silent distribution fallback")
    A = float(np.dot(logs, logs)/denominator)
    exp_A = np.exp(A)
    divisor = exp_A+rho-1
    if divisor == 0:
        raise ValueError("BGC concentration estimator is singular; inspect the input distribution")
    kappa = max(float((2*exp_A-1)/divisor), 10.)
    if not np.isfinite(kappa):
        raise ValueError("BGC concentration is nonfinite")
    return dict(rho=rho, kappa=kappa, alpha=kappa*rho, beta=kappa*(1-rho), A=A,
        shifted_trimmed_mean=shifted_mean, fit_trim_each_tail=.001, rho_lower_bound=.0001,
        kappa_lower_bound=10., universe_sites=int(universe), occupied_sites=len(positive))


def beta_geometric_cdf(x, alpha, beta):
    """P(X<=x)=1-B(alpha,beta+floor(x))/B(alpha,beta), support 1,2,..."""
    if not np.isfinite(alpha) or not np.isfinite(beta) or alpha <= 0 or beta <= 0:
        raise ValueError("Beta-geometric shape parameters must be finite and positive")
    x = np.asarray(x, dtype=float)
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("CDF counts must be finite and nonnegative")
    return np.clip(-np.expm1(np.minimum(betaln(alpha,beta+np.floor(x))-betaln(alpha,beta),0)),0,1)


def correct(matrix, groups, cdf="analytic", seed=0, draws=200000, post_ttr=True, trim=.05, target=100):
    matrix = sparse.csc_matrix(matrix, dtype=float, copy=True)
    matrix.eliminate_zeros()
    if (cdf not in {"analytic","simulation"} or not isinstance(seed,(int,np.integer)) or seed < 0 or
            not isinstance(draws,(int,np.integer)) or draws < 1 or matrix.shape[1] != len(groups) or
            not matrix.shape[0] or np.any(~np.isfinite(matrix.data)) or np.any(matrix.data < 0)):
        raise ValueError("Invalid BGC engine, nonnegative integer seed, positive draws, matrix or groups")
    output = matrix.copy()
    occupied = np.diff(matrix.indptr)
    if np.any(occupied == 0):
        raise ValueError("Empty samples cannot define BGC distributions")
    means = np.asarray(matrix.sum(axis=0)).ravel()/occupied
    details = []
    for j, group in enumerate(groups):
        positive = matrix.data[matrix.indptr[j]:matrix.indptr[j+1]]
        fit = fit_parameters(positive, matrix.shape[0])
        target_mean = float(np.mean(means[np.array([g==group for g in groups])]))
        p_target = 1/target_mean
        if not 0 < p_target <= 1:
            raise ValueError("BGC target geometric probability must be in (0,1]")
        unique, inverse, frequency = np.unique(positive, return_inverse=True, return_counts=True)
        # Derive an order-independent sample RNG from its fitted count distribution.
        digest = hashlib.sha256(np.sort(positive).astype("<f8").tobytes()).hexdigest()
        child_seed = int.from_bytes(hashlib.sha256(f"{seed}:{matrix.shape[0]}:{digest}".encode()).digest()[:8],"little")
        if cdf == "analytic":
            probability = beta_geometric_cdf(unique, fit["alpha"], fit["beta"])
        else:
            rng = np.random.default_rng(child_seed)
            p = rng.beta(fit["alpha"], fit["beta"], size=draws)
            if np.any(p <= 0) or np.any(p > 1):
                raise ValueError("Simulated BGC beta probabilities are invalid; use analytic CDF or inspect fit")
            simulated = np.sort(rng.geometric(p))
            if np.any(simulated < 1):
                raise ValueError("Simulated BGC counts overflowed; use analytic CDF or inspect fit")
            probability = np.searchsorted(simulated, unique, side="right")/draws
        upper = np.nextafter(1.,0.)
        lower_count = int(frequency[probability <= 0].sum())
        upper_count = int(frequency[probability >= 1].sum())
        # Preserve positive observations even when finite simulation yields ECDF=0.
        bounded = np.clip(probability, np.nextafter(0.,1.), upper)
        transformed = np.ones(len(unique)) if p_target == 1 else geom.ppf(bounded,p_target)
        transformed = np.maximum(transformed,1.)
        if np.any(~np.isfinite(transformed)):
            raise ValueError("BGC produced nonfinite transformed counts")
        values = transformed[inverse]
        scale = ttr_parameters(values, matrix.shape[0], trim, target)["factor"] if post_ttr else 1.
        output.data[output.indptr[j]:output.indptr[j+1]] = values*scale
        details.append(dict(bgc_rho=fit["rho"], bgc_kappa=fit["kappa"], bgc_alpha=fit["alpha"], bgc_beta=fit["beta"],
            bgc_A=fit["A"], bgc_shifted_trimmed_mean=fit["shifted_trimmed_mean"], bgc_fit_trim_each_tail=.001,
            bgc_target_geometric_mean=target_mean, bgc_cdf=cdf, bgc_seed=seed,
            bgc_sample_seed=str(child_seed) if cdf=="simulation" else "not_used",
            bgc_draws=draws if cdf=="simulation" else 0, bgc_post_ttr=post_ttr, bgc_post_ttr_factor=scale,
            bgc_post_ttr_trim_each_tail=trim if post_ttr else "", bgc_post_ttr_target=target if post_ttr else "",
            bgc_lower_probability_safeguard_sites=lower_count, bgc_upper_probability_safeguard_sites=upper_count,
            bgc_raw_top_1pct_support_fraction=float(np.sort(positive)[-max(1,int(np.ceil(len(positive)*.01))):].sum()/positive.sum()),
            bgc_corrected_top_1pct_support_fraction=float(np.sort(values)[-max(1,int(np.ceil(len(values)*.01))):].sum()/values.sum())))
    return output.tocsr(), details
