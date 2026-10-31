"""Experimental control-conditioned Hermes Gumbel (count-aware longest-run extension), adapted from TRANSIT's Gumbel method (DeJesus et al. 2013).

Input is a YTAB scalar site normalization (ttr, totreads or nzmean); no normalization
is computed here. The model uses the catalogue's raw integer counts with each
observation's exposure taken from its normalization factor (1/factor).
Opportunities are control-supported coordinates. A site is encoded as depleted
when its treatment count is improbably low under a shared-dispersion negative
binomial model conditioned on the matched controls. The longest run of depleted
sites in each gene is scored against a discrete Gumbel fitted to simulated
no-effect runs and checked on independent simulations and the exact conditional
run distribution. Experimental run-based depletion evidence, not a validated
essentiality or differential-fitness classifier.
"""
from dataclasses import asdict, dataclass
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp
from scipy.stats import nbinom, norm

from ytab.normalization.ttr import sha256, write_table

MODEL = "experimental control-conditioned Hermes Gumbel"
REVISION = "ccg_v5_normalized_input"
NORMALIZATIONS = ("ttr", "totreads", "nzmean")   # scalar methods whose factors are valid count-model exposures
DEFAULT_NORMALIZATION = "ttr"
DEFAULT_SIMULATIONS = {"exact":1000, "gumbel-mc":20000}   # exact: Gumbel diagnostic only; gumbel-mc: tested setting
LOG_TAU_BOUNDS = (np.log(1e-3), np.log(20.))
LOG_PHI_BOUNDS = (np.log(1e-6), np.log(1e2))
LOG_LINEAR_BOUNDS = (np.log(1e-6), np.log(1e4))
BETA_BOUNDS = (1e-3, 1e3)
_GRID = np.linspace(-8, 8, 1601)
_GRID_WEIGHTS = norm.pdf(_GRID)/norm.pdf(_GRID).sum()


@dataclass(frozen=True)
class Settings:
    alpha_site: float = .05
    min_control_reads: int = 1
    eligibility: str = "control-support"
    min_sites: int = 3
    calibration_sims: int | None = None   # None: DEFAULT_SIMULATIONS for the inference mode
    validation_sims: int | None = None
    seed: int = 2026
    nodes: int = 32
    fdr: float = .05
    tail_tolerance: float = 2.
    ks_confidence: float = .99
    nterm: float = 0.
    cterm: float = 0.
    ignore_stop: bool = True
    treatment_dispersion: str = "replicates"
    dispersion_model: str = "linear-quadratic"
    inference: str = "exact"

    def __post_init__(self):
        default = DEFAULT_SIMULATIONS.get(self.inference, 1000)
        for name in ("calibration_sims", "validation_sims"):
            if getattr(self, name) is None:object.__setattr__(self, name, default)

    def validate(self):
        integers = (self.min_control_reads, self.min_sites, self.calibration_sims,
                    self.validation_sims, self.seed, self.nodes)
        if (any(not isinstance(v, int) or isinstance(v, bool) for v in integers) or
                self.min_control_reads < 1 or self.min_sites < 1 or self.calibration_sims < 50 or
                self.validation_sims < 50 or self.seed < 0 or not 4 <= self.nodes <= 200):
            raise ValueError("Require min_control_reads>=1, min_sites>=1, >=50 calibration/validation simulations, seed>=0 and 4..200 nodes")
        if (not 0 < self.alpha_site < .5 or not 0 < self.fdr < 1 or not self.tail_tolerance > 1 or
                not .5 <= self.ks_confidence < 1 or self.eligibility not in {"control-support", "detectable"} or
                self.treatment_dispersion not in {"replicates", "control"} or
                self.dispersion_model not in {"linear-quadratic", "quadratic"} or
                self.inference not in DEFAULT_SIMULATIONS or
                not 0 <= self.nterm < 100 or not 0 <= self.cterm < 100 or self.nterm+self.cterm >= 100):
            raise ValueError("Invalid site alpha, FDR, tail tolerance, KS confidence, eligibility, dispersion, inference settings or terminal percentages")


def annotations(path):
    """Resolve CDS parents through transcript ancestry, using each gene's own CDS."""
    records = {}
    cds = []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):break
            if line.startswith("#") or not line.strip():continue
            fields=line.rstrip().split("\t")
            if len(fields)!=9:raise ValueError("Expected nine-column GFF3")
            attrs={unquote(k):unquote(v) for k,v in
                   (entry.split("=",1) for entry in fields[8].split(";") if "=" in entry)}
            rec=dict(chrom=fields[0],kind=fields[2],start=int(fields[3]),end=int(fields[4]),
                     strand=fields[6],attrs=attrs)
            if rec["start"]<1 or rec["end"]<rec["start"]:raise ValueError("Invalid inclusive GFF coordinates")
            if fields[2]=="CDS":cds.append(rec)
            elif attrs.get("ID"):
                if attrs["ID"] in records:raise ValueError("Duplicate non-CDS GFF ID")
                records[attrs["ID"]]=rec
    genes={}
    def ancestor(identifier,seen):
        if identifier in seen:raise ValueError("Cyclic GFF ancestry")
        if identifier not in records:raise ValueError(f"Unresolved GFF CDS parent: {identifier}")
        record=records[identifier]
        if record["kind"] in {"gene","pseudogene"} or not record["attrs"].get("Parent"):
            return [identifier]
        return [g for p in record["attrs"]["Parent"].split(",") for g in ancestor(p,seen|{identifier})]
    for record in cds:
        parents=record["attrs"].get("Parent","").split(",")
        if not all(parents):raise ValueError("CDS must have explicit GFF Parent")
        for gid in {g for p in parents for g in ancestor(p,set())}:
            parent=records[gid]
            if parent["chrom"]!=record["chrom"] or parent["strand"]!=record["strand"]:
                raise ValueError("CDS and gene disagree on chromosome/strand")
            identifier=parent["attrs"].get("locus_tag",gid)
            gene=genes.setdefault(gid,dict(id=identifier,gff_id=gid,name=parent["attrs"].get("gene") or parent["attrs"].get("Gene") or parent["attrs"].get("Name",identifier),
                description=parent["attrs"].get("description",parent["attrs"].get("Note","")),
                chrom=parent["chrom"],strand=parent["strand"],blocks=[]))
            gene["blocks"].append((record["start"],record["end"]))
    if not genes:raise ValueError("No protein-coding CDS groups found in GFF3")
    if len({g["id"] for g in genes.values()})!=len(genes):raise ValueError("Duplicate gene identifiers")
    for gene in genes.values():
        merged=[]
        for a,b in sorted(set(gene["blocks"])):
            if merged and a<=merged[-1][1]+1:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
            else:merged.append((a,b))
        gene["blocks"]=merged
    return sorted(genes.values(),key=lambda g:g["id"])


def _integer_counts(values, name):
    values = np.asarray(values)
    if (not np.issubdtype(values.dtype, np.number) or np.any(~np.isfinite(values)) or
            np.any(values < 0) or np.any(values != np.floor(values))):
        raise ValueError(f"{name} must be raw nonnegative integer counts")
    return values.astype(np.int64)


def nb_logpmf(counts, mean, phi):
    """NB with E=mean and Var=mean+phi*mean^2 (phi may vary with mean)."""
    size = 1/phi
    return (gammaln(counts+size)-gammaln(size)-gammaln(counts+1)
            - counts*np.log1p(size/mean)-size*np.log1p(mean/size))


def effective_phi(mean, dispersion):
    """Dispersion (phi, linear) -> NB phi at this mean: Var = mean*(1+linear) + phi*mean^2.

    The linear term is extra variance proportional to the mean (bottleneck-like sampling of
    clones); a scalar dispersion means no linear term.
    """
    phi, linear = (dispersion, 0.) if np.ndim(dispersion) == 0 else dispersion
    return phi+linear/mean if linear else phi


def selection_log_probability(m, tau, dispersion, exposures, min_reads):
    """log P(sum of control counts >= min_reads) under the lognormal-NB background."""
    mean = np.exp(m+tau*_GRID)[:, None]*exposures
    if min_reads == 1:
        below = np.exp(nb_logpmf(0, mean, effective_phi(mean, dispersion)).sum(1))
    else:
        pmf = np.exp(nb_logpmf(np.arange(min_reads)[None, None, :], mean[..., None], effective_phi(mean[..., None], dispersion)))
        below = pmf[:, 0, :]
        for j in range(1, len(exposures)):
            below = np.stack([(below[:, :k+1]*pmf[:, j, k::-1]).sum(1) for k in range(min_reads)], axis=1)
        below = below.sum(1)
    return float(np.log1p(-min(float(np.dot(_GRID_WEIGHTS, below)), 1.)))


def _log_joint(u, counts, exposures, m, tau, dispersion):
    mean = np.exp(u)[..., None]*exposures
    return nb_logpmf(counts, mean, effective_phi(mean, dispersion)).sum(-1)+norm.logpdf(u, m, tau)


def site_posterior(counts, exposures, m, tau, dispersion, nodes, start=None):
    """Adaptive Gauss-Hermite nodes for p(log lambda | counts).

    Returns node locations (n,Q), normalized log posterior weights (n,Q) and log marginal likelihoods.
    `start` optionally warm-starts the mode search (e.g. from a nearby parameter value); results do not depend on it.
    """
    counts = np.asarray(counts, dtype=float)
    phi, linear = (dispersion, 0.) if np.ndim(dispersion) == 0 else dispersion
    u = np.log((counts.sum(1)+.5)/exposures.sum()) if start is None else np.array(start, dtype=float)
    def derivatives(u):
        if not linear:
            size = 1/phi; frac = np.exp(u)[:, None]*exposures; frac = frac/(frac+size)
            return (counts-(counts+size)*frac).sum(1)-(u-m)/tau**2, -((counts+size)*frac*(1-frac)).sum(1)-1/tau**2
        h = 1e-3   # mean-dependent dispersion: central differences of the exact log joint
        f0, fp, fm = (_log_joint(u+d, counts, exposures, m, tau, dispersion) for d in (0, h, -h))
        return (fp-fm)/(2*h), np.minimum((fp-2*f0+fm)/h**2, -1/tau**2)
    for _ in range(200):
        gradient, curvature = derivatives(u)
        step = np.clip(-gradient/curvature, -2, 2); u = u+step
        # Finite-difference steps are noisy near 1e-8; 1e-6 is ample for centring adaptive quadrature.
        if np.max(np.abs(step), initial=0) < (1e-6 if linear else 1e-8):break
    else:raise ValueError("Site posterior mode did not converge")
    sd = 1/np.sqrt(-derivatives(u)[1])
    x, w = hermgauss(nodes)
    locations = u[:, None]+np.sqrt(2)*sd[:, None]*x
    terms = np.log(w)+x**2+_log_joint(locations, counts[:, None, :], exposures, m, tau, dispersion)
    total = logsumexp(terms, axis=1)
    return locations, terms-total[:, None], total+np.log(np.sqrt(2)*sd)


def _unpack(theta, settings):
    m, log_tau, log_phi = theta[:3]
    linear = float(np.exp(theta[3])) if settings.dispersion_model == "linear-quadratic" else 0.
    return float(m), float(np.exp(log_tau)), (float(np.exp(log_phi)), linear)


def _background_objective(theta, rows, weights, exposures, settings, modes):
    m, tau, dispersion = _unpack(theta, settings)
    logsel = selection_log_probability(m, tau, dispersion, exposures, settings.min_control_reads)
    if not np.isfinite(logsel):return 1e10
    locations, _, marginal = site_posterior(rows, exposures, m, tau, dispersion, settings.nodes, modes.get("u"))
    if np.all(np.isfinite(locations)):modes["u"] = locations.mean(1)   # node centre = posterior mode
    value = -(weights*marginal).sum()/weights.sum()+logsel
    return float(value) if np.isfinite(value) else 1e10


def _replicate_summary(controls, exposures):
    controls = np.asarray(controls); total = controls.sum(1)
    both = (controls[:, 0] > 0) & (controls[:, 1] > 0)
    ratio = np.log2((controls[both, 0]/exposures[0])/(controls[both, 1]/exposures[1]))
    mean = total/2; levels = {}
    for lo, hi in ((5, 20), (20, 100), (100, 1000)):
        keep = both & (mean >= lo) & (mean < hi)
        r = np.log2((controls[keep, 0]/exposures[0])/(controls[keep, 1]/exposures[1]))
        levels[f"{lo}-{hi}"] = float(np.median(np.abs(r))) if len(r) else None
    return dict(any_control_zero=float(np.mean((controls == 0).any(1))),
        control_sum_quantiles={str(q):float(np.percentile(total, q)) for q in (10, 50, 90, 99)},
        median_abs_log2_replicate_ratio=float(np.median(np.abs(ratio))) if len(ratio) else None,
        median_abs_log2_replicate_ratio_by_mean_count=levels)


def fit_background(controls, exposures, settings):
    """Selection-aware empirical-Bayes fit of log lambda~N(m,tau^2) and shared dispersion from replicate tracks."""
    rows, weights = np.unique(np.asarray(controls, dtype=np.int64), axis=0, return_counts=True)
    if weights.sum() < 20:raise ValueError("Need at least 20 supported sites to fit the background")
    naive = np.log(rows.sum(1)/exposures.sum())
    m0 = float(np.average(naive, weights=weights))
    tau0 = float(np.clip(np.sqrt(np.average((naive-m0)**2, weights=weights)), .1, 10))
    bounds = [(None, None), LOG_TAU_BOUNDS, LOG_PHI_BOUNDS]
    starts = [(m0, np.log(tau0), np.log(.1)), (m0-.5, np.log(2*tau0), np.log(1.))]
    if settings.dispersion_model == "linear-quadratic":
        bounds.append(LOG_LINEAR_BOUNDS); starts = [starts[0]+(np.log(.1),), starts[1]+(np.log(2.),)]
    fits = [minimize(_background_objective, start, args=(rows, weights, exposures, settings, {}), method="L-BFGS-B", bounds=bounds)
            for start in starts]
    converged = [f for f in fits if f.success and np.all(np.isfinite(f.x))]
    if not converged:raise ValueError("Background fit failed to converge from either start; no fallback")
    best = min(converged, key=lambda f:f.fun)
    m, tau, dispersion = _unpack(best.x, settings)
    coarse = site_posterior(rows, exposures, m, tau, dispersion, settings.nodes)[2]
    fine = site_posterior(rows, exposures, m, tau, dispersion, 2*settings.nodes)[2]
    named = [("tau", best.x[1], LOG_TAU_BOUNDS), ("phi", best.x[2], LOG_PHI_BOUNDS)]
    if settings.dispersion_model == "linear-quadratic":named.append(("linear", best.x[3], LOG_LINEAR_BOUNDS))
    boundary = [name for name, value, (lo, hi) in named if value-lo < 1e-3 or hi-value < 1e-3]
    # Posterior predictive check from the fitted, selection-conditioned generative model.
    rng = np.random.default_rng(np.random.SeedSequence(settings.seed, spawn_key=(0,)))
    draws = min(4*int(weights.sum()), 400000)
    mean = np.exp(m+tau*rng.standard_normal(draws))[:, None]*exposures
    size = 1/effective_phi(mean, dispersion)
    simulated = rng.negative_binomial(size, size/(size+mean))
    simulated = simulated[simulated.sum(1) >= settings.min_control_reads]
    observed = np.repeat(rows, weights, axis=0)
    return dict(m=m, tau=tau, phi=dispersion[0], linear=dispersion[1], dispersion_model=settings.dispersion_model, converged=True,
        starts_converged=len(converged), starts=len(fits),
        objective_per_site=float(best.fun), start_parameter_spread=float(np.ptp([f.x for f in converged], axis=0).max()) if len(converged) > 1 else None,
        start_objective_spread=float(np.ptp([f.fun for f in converged])) if len(converged) > 1 else None,
        boundary=boundary, sites=int(weights.sum()), unique_profiles=len(rows),
        selection_probability=float(np.exp(selection_log_probability(m, tau, dispersion, exposures, settings.min_control_reads))),
        integration_check=dict(nodes=settings.nodes, comparison_nodes=2*settings.nodes,
            max_abs_log_marginal_difference=float(np.max(np.abs(coarse-fine))),
            total_log_likelihood_difference=float(np.abs((weights*(coarse-fine)).sum()))),
        posterior_predictive=dict(observed=_replicate_summary(observed, exposures),
            simulated=_replicate_summary(simulated, exposures), simulated_selected_sites=len(simulated)))


def predictive_cdf(t, locations, log_weights, treatment_exposure, dispersion):
    """Control-conditioned predictive P(T<=t) under no treatment effect."""
    mean = treatment_exposure*np.exp(locations); size = 1/effective_phi(mean, dispersion)
    cdf = nbinom.cdf(np.asarray(t)[:, None], size, size/(size+mean))
    return np.clip((np.exp(log_weights)*cdf).sum(1), 0, 1)


def depletion_thresholds(locations, log_weights, treatment_exposure, dispersion, alpha):
    """Largest t with F(t)<=alpha (-1 when a zero count would not be flagged) and null flag probability F(t*)."""
    n = len(locations); threshold = np.full(n, -1, dtype=np.int64); probability = np.zeros(n)
    args = (treatment_exposure, dispersion)
    index = np.flatnonzero(predictive_cdf(np.zeros(n), locations, log_weights, *args) <= alpha)
    if len(index):
        u, w = locations[index], log_weights[index]
        lo = np.zeros(len(index), dtype=np.int64)
        hi = np.maximum(1, np.ceil((np.exp(w)*treatment_exposure*np.exp(u)).sum(1))).astype(np.int64)
        while True:
            low = predictive_cdf(hi, u, w, *args) <= alpha
            if not low.any():break
            lo[low] = hi[low]; hi[low] *= 2
        while np.any(hi-lo > 1):
            active = np.flatnonzero(hi-lo > 1); middle = (lo[active]+hi[active])//2
            low = predictive_cdf(middle, u[active], w[active], *args) <= alpha
            lo[active[low]] = middle[low]; hi[active[~low]] = middle[~low]
        threshold[index] = lo
        probability[index] = predictive_cdf(lo, u, w, *args)
    return threshold, probability


def longest_run(flags):
    """Longest consecutive True run; the last tied maximum wins, as in the source convention."""
    best = run = 0; start = end = None
    for i, flag in enumerate(flags):
        run = run+1 if flag else 0
        if flag and run >= best:best, start, end = run, i-run+1, i
    return best, start, end


def longest_runs(flags):
    current = np.zeros(len(flags), dtype=np.int64); best = current.copy()
    for column in np.asarray(flags).T:
        current = (current+1)*column; np.maximum(best, current, out=best)
    return best


def exact_run_tails(probabilities, kmax):
    """P(R>=k), k=1..kmax, for independent heterogeneous Bernoulli flags (no subtraction)."""
    size = max(int(kmax), 1); state = np.zeros((size, size)); state[:, 0] = 1.
    allowed = np.tril(np.ones((size, size), dtype=bool)); diagonal = np.arange(size); tails = np.zeros(size)
    for p in probabilities:
        tails += p*state[diagonal, diagonal]
        updated = np.empty_like(state); updated[:, 0] = (1-p)*state.sum(1); updated[:, 1:] = p*state[:, :-1]
        state = updated*allowed
    return tails


def _log_gumbel_cdf(x, mu, beta):
    with np.errstate(over="ignore"):
        return -np.exp(-(np.asarray(x, dtype=float)-mu)/beta)


def _log_difference(upper, lower):
    """log(exp(upper)-exp(lower)) for upper>=lower; -inf when both CDF values underflow."""
    with np.errstate(divide="ignore", invalid="ignore"):
        value = upper+np.log(-np.expm1(np.minimum(lower-upper, 0)))
    return np.where(np.isneginf(upper), -np.inf, value)


def discrete_gumbel_logpmf(r, n, mu, beta):
    """log P(R=r) for a Gumbel rounded to integers and truncated to 0..n."""
    r = np.asarray(r, dtype=float)
    normalizer = _log_difference(_log_gumbel_cdf(n+.5, mu, beta), _log_gumbel_cdf(-.5, mu, beta))
    with np.errstate(invalid="ignore"):
        value = _log_difference(_log_gumbel_cdf(r+.5, mu, beta), _log_gumbel_cdf(r-.5, mu, beta))-normalizer
    return np.where((r >= 0) & (r <= n) & np.isfinite(value), value, -np.inf)


def discrete_gumbel_cdf(k, n, mu, beta):
    k = np.minimum(np.asarray(k, dtype=float), n)
    normalizer = _log_difference(_log_gumbel_cdf(n+.5, mu, beta), _log_gumbel_cdf(-.5, mu, beta))
    return np.where(k < 0, 0., np.exp(_log_difference(_log_gumbel_cdf(k+.5, mu, beta), _log_gumbel_cdf(-.5, mu, beta))-normalizer))


def discrete_gumbel_logsf(r, n, mu, beta):
    """log P(R>=r) for 1<=r<=n, kept in log space so tiny tails stay differentiable."""
    normalizer = _log_difference(_log_gumbel_cdf(n+.5, mu, beta), _log_gumbel_cdf(-.5, mu, beta))
    return float(_log_difference(_log_gumbel_cdf(n+.5, mu, beta), _log_gumbel_cdf(r-.5, mu, beta))-normalizer)


def _hessian(f, x, h=1e-3):
    x = np.asarray(x, dtype=float); steps = np.eye(len(x))*h; H = np.empty((len(x), len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            H[i, j] = (f(x+steps[i]+steps[j])-f(x+steps[i]-steps[j])-f(x-steps[i]+steps[j])+f(x-steps[i]-steps[j]))/(4*h*h)
    return (H+H.T)/2


def gumbel_log_tail_se(fit, r, n, h=1e-5):
    """Delta-method Monte Carlo standard error of log P(R>=r) from the ML covariance of (mu, log beta).

    None when the observed information is not positive definite (the caller then allows no Monte Carlo slack).
    """
    if fit.get("covariance") is None or not 1 <= r <= n:return None
    theta = np.array([fit["mu"], np.log(fit["beta"])]); gradient = np.empty(2)
    for i in range(2):
        step = np.eye(2)[i]*h
        gradient[i] = (discrete_gumbel_logsf(r, n, theta[0]+step[0], np.exp(theta[1]+step[1]))-
                       discrete_gumbel_logsf(r, n, theta[0]-step[0], np.exp(theta[1]-step[1])))/(2*h)
    variance = float(gradient@np.asarray(fit["covariance"])@gradient)
    return float(np.sqrt(variance)) if np.isfinite(variance) and variance >= 0 else None


def discrete_gumbel_sf(r, n, mu, beta):
    """Inclusive upper tail P(R>=r); 1 for r<=0 and 0 for r>n."""
    r = float(r)
    if r <= 0:return 1.
    if r > n:return 0.
    normalizer = _log_difference(_log_gumbel_cdf(n+.5, mu, beta), _log_gumbel_cdf(-.5, mu, beta))
    return float(np.exp(_log_difference(_log_gumbel_cdf(n+.5, mu, beta), _log_gumbel_cdf(r-.5, mu, beta))-normalizer))


def fit_discrete_gumbel(runs, n):
    """Maximum-likelihood rounded/truncated Gumbel on calibration runs; failures are reported, not replaced."""
    values, counts = np.unique(np.asarray(runs, dtype=np.int64), return_counts=True)
    if len(values) < 2:
        return dict(converged=False, reason=f"degenerate calibration null (all runs={int(values[0])})", mu=None, beta=None)
    mean = float(np.average(values, weights=counts)); sd = float(np.sqrt(np.average((values-mean)**2, weights=counts)))
    def objective(theta):
        # Floor log-probabilities so finite-difference gradients stay finite at extreme trial parameters.
        return float(-(counts*np.maximum(discrete_gumbel_logpmf(values, n, theta[0], np.exp(theta[1])), -700.)).sum())
    beta0 = max(sd*np.sqrt(6)/np.pi, .1)
    bounds = [(None, None), tuple(np.log(BETA_BOUNDS))]
    fits = [minimize(objective, start, method="L-BFGS-B", bounds=bounds)
            for start in ((mean-np.euler_gamma*beta0, np.log(beta0)), (float(np.median(np.repeat(values, counts))), 0.))]
    good = [f for f in fits if f.success and np.isfinite(f.fun)]
    if not good:return dict(converged=False, reason="optimizer did not converge", mu=None, beta=None)
    best = min(good, key=lambda f:f.fun); mu, beta = float(best.x[0]), float(np.exp(best.x[1]))
    if not BETA_BOUNDS[0]*1.001 < beta < BETA_BOUNDS[1]/1.001:
        return dict(converged=False, reason=f"scale at boundary ({beta:.3g})", mu=mu, beta=beta)
    spread = float(np.ptp([f.x for f in good], axis=0).max()) if len(good) > 1 else None
    # Observed information in (mu, log beta); its inverse is the Monte Carlo covariance of the fit.
    H = _hessian(objective, best.x)
    covariance = np.linalg.inv(H).tolist() if np.all(np.isfinite(H)) and np.all(np.linalg.eigvalsh(H) > 0) else None
    return dict(converged=True, reason="", mu=mu, beta=beta, negative_log_likelihood=float(best.fun), start_spread=spread,
                covariance=covariance)


def _gene_rng(seed, gene_id, stream):
    key = int.from_bytes(hashlib.sha256(str(gene_id).encode()).digest()[:8], "big")
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(key, stream)))


def simulate_runs(rng, locations, log_weights, threshold, treatment_exposure, dispersion, sims):
    """No-effect treatment counts drawn from the control-conditioned predictive, encoded exactly as observed."""
    n, nodes = locations.shape; flags = np.zeros((sims, n), dtype=bool)
    active = np.flatnonzero(threshold >= 0)
    if len(active):
        cumulative = np.cumsum(np.exp(log_weights[active]), axis=1); cumulative[:, -1] = 1.
        offsets = np.arange(len(active))
        flat = (cumulative+offsets[:, None]).ravel()
        chunk = max(1, int(2e6//len(active)))
        for begin in range(0, sims, chunk):
            rows = min(chunk, sims-begin)
            node = np.searchsorted(flat, rng.random((rows, len(active)))+offsets, side="right")-offsets*nodes
            mean = treatment_exposure*np.exp(locations[active][offsets, np.minimum(node, nodes-1)])
            size = 1/effective_phi(mean, dispersion)
            counts = rng.negative_binomial(size, size/(size+mean))
            flags[begin:begin+rows, active] = counts <= threshold[active]
    return longest_runs(flags)


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float); order = np.argsort(p); ranked = p[order]*len(p)/np.arange(1, len(p)+1)
    adjusted = np.empty_like(p); adjusted[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return adjusted


def _log2_ratio(treatment, control, treatment_exposure, control_exposure):
    if control <= 0 or treatment <= 0:return None
    return float(np.log2((treatment/treatment_exposure)/(control/control_exposure)))


def gene_sites(genes, sites, settings):
    """Catalogue indices in each gene's own CDS union (N/C-terminal trimming and last 3 bases excluded, as in TRANSIT's Gumbel method)."""
    by_chrom = {}
    for i, (chrom, p) in enumerate(sites):by_chrom.setdefault(chrom, []).append((p, i))
    chroms = {}
    for chrom, entries in by_chrom.items():
        entries.sort(); positions, indices = np.asarray(entries, dtype=np.int64).T; chroms[chrom] = (positions, indices)
    result = []
    for gene in genes:
        positions, indices = chroms.get(gene["chrom"], (np.array([], dtype=np.int64), np.array([], dtype=np.int64)))
        selected = np.asarray([k for a, b in gene["blocks"]
                               for k in range(np.searchsorted(positions, a), np.searchsorted(positions, b, side="right"))], dtype=int)
        p = positions[selected]; idx = indices[selected]
        lengths = np.array([b-a+1 for a, b in gene["blocks"]]); length = int(lengths.sum())
        offset = np.zeros(len(p), dtype=int); cumulative = 0
        for (a, b), width in zip(gene["blocks"], lengths):
            mask = (p >= a) & (p <= b); offset[mask] = cumulative+p[mask]-a; cumulative += int(width)
        oriented = offset if gene["strand"] != "-" else length-1-offset
        keep = (oriented >= length*settings.nterm/100) & (oriented < length*(1-settings.cterm/100))
        if settings.ignore_stop:keep &= oriented < length-3
        result.append(idx[keep])
    return result


def combined_dispersion(dispersion, exposures):
    """Moment-matched dispersion of a sum of independent replicate tracks.

    Var_k = e_k*lambda*(1+linear) + phi*(e_k*lambda)^2 sums to mu*(1+linear) + phi*sum(e^2)/sum(e)^2*mu^2
    with mu = sum(e)*lambda: the linear term is unchanged and phi is rescaled.
    """
    phi, linear = (dispersion, 0.) if np.ndim(dispersion) == 0 else dispersion
    exposures = np.asarray(exposures, dtype=float)
    scaled = float(phi*np.sum(exposures**2)/exposures.sum()**2)
    return scaled if np.ndim(dispersion) == 0 else (scaled, float(linear))


def _gene_result(gene, idx, control, treatment, positions, threshold, flag_probability,
                 locations, log_weights, control_exposure, treatment_exposure, treatment_dispersion, settings, tail_floor):
    n = len(idx)
    row = dict(ORF=gene["id"], Name=gene["name"], Description=gene["description"], chromosome=gene["chrom"],
               strand=gene["strand"], n_eligible=n)
    if not n:
        return row | dict(status="no control-supported sites", in_family=False)
    flags = treatment[idx] <= threshold[idx]
    r, a, b = longest_run(flags)
    control_reads = int(control[idx].sum()); treatment_reads = int(treatment[idx].sum())
    row |= dict(n_detectable=int((threshold[idx] >= 0).sum()), n_depleted=int(flags.sum()),
        expected_depleted_null=float(flag_probability[idx].sum()), r=r,
        run_start=int(positions[idx[a]]) if r else None, run_end=int(positions[idx[b]]) if r else None,
        run_span_bp=int(positions[idx[b]]-positions[idx[a]]+1) if r else 0,
        control_reads=control_reads, treatment_reads=treatment_reads, treatment_reads_zero=treatment_reads == 0,
        log2_scaled_ratio=_log2_ratio(treatment_reads, control_reads, treatment_exposure, control_exposure))
    if r:
        run = idx[a:b+1]
        row |= dict(run_control_reads=int(control[run].sum()), run_treatment_reads=int(treatment[run].sum()),
            run_log2_scaled_ratio=_log2_ratio(int(treatment[run].sum()), int(control[run].sum()), treatment_exposure, control_exposure))
    if n < settings.min_sites:
        return row | dict(status="too few control-supported sites", in_family=False)
    args = (locations[idx], log_weights[idx], threshold[idx], treatment_exposure, treatment_dispersion)
    calibration = simulate_runs(_gene_rng(settings.seed, gene["id"], 1), *args, settings.calibration_sims)
    validation = simulate_runs(_gene_rng(settings.seed, gene["id"], 2), *args, settings.validation_sims)
    fit = fit_discrete_gumbel(calibration, n)
    kmax = int(max(r, calibration.max(), validation.max()))+1
    tails = exact_run_tails(flag_probability[idx], kmax)
    exact_cdf = 1-tails; k = np.arange(kmax)
    p_exact = float(tails[r-1]) if r else 1.
    p_sim = float((1+(validation >= r).sum())/(settings.validation_sims+1))
    mc = settings.inference == "gumbel-mc"
    # DKW covers validation sampling only; gumbel-mc also allows for the calibration sample the fit came from.
    bound = float(np.sqrt(np.log(2/(1-settings.ks_confidence))/2*(1/settings.validation_sims+(1/settings.calibration_sims if mc else 0))))
    null_max = int(max(calibration.max(), validation.max()))
    row |= dict(null_mean_run=float(calibration.mean()), null_max_run=null_max,
        gumbel_mu=fit["mu"], gumbel_beta=fit["beta"], fit_converged=fit["converged"], ks_bound=bound,
        exact_tail_p=p_exact, simulation_tail_p=p_sim, tail_extrapolated=bool(r > null_max))
    reasons = []
    if fit["converged"]:
        cdf = discrete_gumbel_cdf(k, n, fit["mu"], fit["beta"])
        ecdf = (validation[None, :] <= k[:, None]).mean(1)
        p_fit = discrete_gumbel_sf(r, n, fit["mu"], fit["beta"])
        row |= dict(validation_ks=float(np.abs(cdf-ecdf).max()), exact_cdf_max_deviation=float(np.abs(cdf-exact_cdf).max()),
                    gumbel_p_fitted=p_fit, gumbel_to_exact_tail_ratio=p_fit/p_exact if r and p_exact > 0 else None)
        if row["validation_ks"] > bound:reasons.append("validation CDF outside DKW bound")
        if r:
            # Relative tail agreement is required wherever it can change inference; below the
            # Bonferroni-scale floor both tails imply rejection and the magnitude is reported as unvalidated.
            # gumbel-mc widens the tolerance by z Monte Carlo standard errors of the fitted log tail.
            se = gumbel_log_tail_se(fit, r, n) if mc else None
            slack = float(norm.ppf(.5+settings.ks_confidence/2))*se if se else 0.
            row["gumbel_log_tail_se"] = se
            agrees = p_exact > 0 and p_fit > 0 and abs(np.log(p_fit/p_exact)) <= np.log(settings.tail_tolerance)+slack
            beyond = p_fit <= tail_floor and p_exact <= tail_floor
            row["tail_check"] = "agrees" if agrees else "both below inference floor; magnitude unvalidated" if beyond else "disagrees"
            if not (agrees or beyond):reasons.append("fitted tail disagrees with exact conditional tail")
    else:reasons.append(fit["reason"])
    row["null_fit_adequate"] = not reasons
    if not row["n_detectable"]:
        return row | dict(status="no depletion-detectable sites", in_family=True, p_family=1.)
    if not r:
        return row | dict(status="no depleted site", in_family=True, p_family=1.)
    row["gumbel_check"] = "; ".join(reasons) or "adequate"
    if settings.inference == "exact":
        # Deterministic p from the exact conditional run distribution; the Gumbel fit is a reported diagnostic only.
        return row | dict(status="exact conditional tail", in_family=True, p_value=p_exact, p_family=p_exact)
    if reasons:
        return row | dict(status="Gumbel null inadequate: "+"; ".join(reasons), in_family=True, p_family=1.)
    return row | dict(status="calibrated", in_family=True, gumbel_p=row["gumbel_p_fitted"], p_value=row["gumbel_p_fitted"], p_family=row["gumbel_p_fitted"])


def analyze(controls, treatments, exposures, genes, gene_indices, positions, settings=None):
    """Run the extension on in-memory arrays.

    controls: sites x J (J>=2) and treatments: sites x K replicate tracks; exposures: J control then K treatment values.
    Replicates are combined (summed) within condition: the combined control defines the null universe and
    baseline, the combined treatment is compared. Each condition's dispersion comes from its own replicate
    disagreement before combining; shared drug effects do not create replicate disagreement.
    """
    settings = settings or Settings(); settings.validate()
    controls = _integer_counts(controls, "Control counts"); treatments = _integer_counts(treatments, "Treatment counts")
    if treatments.ndim == 1:treatments = treatments[:, None]
    exposures = np.asarray(exposures, dtype=float); positions = np.asarray(positions, dtype=np.int64)
    if (controls.ndim != 2 or controls.shape[1] < 2 or treatments.ndim != 2 or len(treatments) != len(controls)
            or len(positions) != len(controls)):
        raise ValueError("Need sites x >=2 control tracks, matching treatment tracks and positions")
    J, K = controls.shape[1], treatments.shape[1]
    if settings.treatment_dispersion == "replicates" and K < 2:
        raise ValueError("Treatment dispersion from replicates needs >=2 treatment tracks; use treatment_dispersion='control' to assume control dispersion explicitly")
    if exposures.shape != (J+K,) or np.any(~np.isfinite(exposures)) or np.any(exposures <= 0):
        raise ValueError("Need one finite positive exposure per control track then per treatment track")
    if len(genes) != len(gene_indices) or len({g["id"] for g in genes}) != len(genes):
        raise ValueError("Gene identities and site lists must align and be unique")
    control_exposures, treatment_exposures = exposures[:J], exposures[J:]
    control, treatment = controls.sum(1), treatments.sum(1)
    control_exposure, treatment_exposure = float(control_exposures.sum()), float(treatment_exposures.sum())
    eligible = control >= settings.min_control_reads
    background = fit_background(controls[eligible], control_exposures, settings)
    if settings.treatment_dispersion == "replicates":
        supported = treatment >= settings.min_control_reads
        treatment_fit = fit_background(treatments[supported], treatment_exposures, settings)
        treatment_replicate = (treatment_fit["phi"], treatment_fit["linear"])
    else:
        treatment_fit = None; treatment_replicate = (background["phi"], background["linear"])
    control_replicate = (background["phi"], background["linear"])
    combined_control_dispersion = combined_dispersion(control_replicate, control_exposures)
    combined_treatment_dispersion = combined_dispersion(treatment_replicate, treatment_exposures)
    needed = np.unique(np.concatenate([np.asarray(i, dtype=np.int64) for i in gene_indices] or [np.array([], dtype=np.int64)]))
    needed = needed[eligible[needed]]
    profiles, inverse = np.unique(control[needed], return_inverse=True); inverse = np.ravel(inverse)
    threshold = np.full(len(control), -1, dtype=np.int64); flag_probability = np.zeros(len(control))
    locations = np.zeros((len(control), settings.nodes)); log_weights = np.full((len(control), settings.nodes), -np.inf)
    observed_cdf = np.full(len(control), np.nan)
    if len(needed):
        u, w, _ = site_posterior(profiles[:, None], np.array([control_exposure]), background["m"], background["tau"], combined_control_dispersion, settings.nodes)
        t, p = depletion_thresholds(u, w, treatment_exposure, combined_treatment_dispersion, settings.alpha_site)
        threshold[needed], flag_probability[needed] = t[inverse], p[inverse]
        locations[needed], log_weights[needed] = u[inverse], w[inverse]
        observed_cdf[needed] = predictive_cdf(treatment[needed], locations[needed], log_weights[needed], treatment_exposure, combined_treatment_dispersion)
    if settings.eligibility == "detectable":eligible &= threshold >= 0
    gene_indices = [np.asarray(idx, dtype=np.int64) for idx in gene_indices]
    # Family fixed from control-defined eligibility before any treatment-dependent result is computed.
    family_size = sum(int(eligible[idx].sum()) >= settings.min_sites for idx in gene_indices)
    tail_floor = settings.fdr/max(family_size, 1)
    rows = []
    for gene, idx in zip(genes, gene_indices):
        rows.append(_gene_result(gene, idx[eligible[idx]], control, treatment, positions, threshold, flag_probability,
                                 locations, log_weights, control_exposure, treatment_exposure, combined_treatment_dispersion, settings, tail_floor))
    family = [r for r in rows if r["in_family"]]
    if family:
        for row, q in zip(family, benjamini_hochberg([r["p_family"] for r in family])):row["q_bh"] = float(q)
    for row in rows:
        row["call"] = ("not tested" if not row["in_family"] else
                       "depleted run" if row["status"] in {"calibrated", "exact conditional tail"} and row["q_bh"] <= settings.fdr else "not significant")
    owners = {}
    for gene, idx in zip(genes, gene_indices):
        for i in idx[eligible[idx]]:owners.setdefault(int(i), []).append(gene["id"])
    site_rows = [dict(site_index=i, genes=";".join(owners[i]), position=int(positions[i]),
                      **{f"control_{j+1}":int(controls[i, j]) for j in range(J)},
                      **{f"treatment_{k+1}":int(treatments[i, k]) for k in range(K)},
                      control_combined=int(control[i]), treatment_combined=int(treatment[i]),
                      predictive_cdf=float(observed_cdf[i]), depletion_threshold=int(threshold[i]),
                      null_flag_probability=float(flag_probability[i]), depleted=bool(treatment[i] <= threshold[i]))
                 for i in sorted(owners)]
    gene_site_index = np.array(sorted(owners), dtype=np.int64)
    statuses = {}
    for row in rows:
        key = row["status"].split(":")[0]; statuses[key] = statuses.get(key, 0)+1
    tested = [r for r in family if r.get("r")]
    diagnostics = dict(background=background, treatment_replicate_fit=treatment_fit,
        dispersion=dict(source=settings.treatment_dispersion, model=settings.dispersion_model,
            variance="Var = mu*(1+linear) + phi*mu^2 per condition",
            control_replicate=dict(phi=control_replicate[0], linear=control_replicate[1]),
            treatment_replicate=dict(phi=treatment_replicate[0], linear=treatment_replicate[1]),
            combined_control=dict(phi=combined_control_dispersion[0], linear=combined_control_dispersion[1]),
            combined_treatment=dict(phi=combined_treatment_dispersion[0], linear=combined_treatment_dispersion[1])),
        combined_exposures=dict(control=control_exposure, treatment=treatment_exposure),
        eligible_sites=int(eligible.sum()), gene_assigned_eligible_sites=len(gene_site_index),
        observed_depleted_sites=int((treatment[gene_site_index] <= threshold[gene_site_index]).sum()) if len(gene_site_index) else 0,
        expected_depleted_sites_null=float(flag_probability[gene_site_index].sum()),
        detectable_sites=int((threshold[gene_site_index] >= 0).sum()) if len(gene_site_index) else 0,
        status_counts=statuses, family_size=len(family), tail_inference_floor=tail_floor,
        family_tail_magnitude_unvalidated=sum(r.get("tail_check", "").startswith("both below") for r in family),
        family_with_depleted_run=len(tested), family_null_fit_inadequate=sum(r["status"].startswith("Gumbel null inadequate") for r in family),
        inference=settings.inference, family_gumbel_check_failed=sum(r.get("gumbel_check", "adequate") != "adequate" for r in tested),
        family_tail_extrapolated=sum(bool(r.get("tail_extrapolated")) for r in tested),
        calls={c:sum(r["call"] == c for r in rows) for c in ("depleted run", "not significant", "not tested")})
    return rows, site_rows, diagnostics


def synthetic_experiment(rng, sites_per_gene, exposures, m=np.log(30), tau=1.5, phi=.2, effects=None, spacing=40,
                         treatments=1, treatment_phi=None, linear=0., treatment_linear=None):
    """Counts generated from the assumed lognormal-NB model, for validation only.

    exposures: control tracks then the last `treatments` treatment tracks; effects: optional per-site
    multiplicative treatment-mean factors (1 = no effect); treatment_phi: replicate dispersion of treatment
    tracks (default phi); linear / treatment_linear: mean-proportional extra variance (bottleneck-like).
    Returns controls, treatments (sites x treatments, squeezed when one), genes, indices, positions.
    """
    sites_per_gene = np.asarray(sites_per_gene, dtype=int); total = int(sites_per_gene.sum())
    exposures = np.asarray(exposures, dtype=float); effects = np.ones(total) if effects is None else np.asarray(effects, dtype=float)
    lam = np.exp(m+tau*rng.standard_normal(total))
    def draw(mean, dispersion):
        size = 1/effective_phi(np.maximum(mean, 1e-12), dispersion); return rng.negative_binomial(size, size/(size+mean))
    treatment_dispersion = (phi if treatment_phi is None else treatment_phi, linear if treatment_linear is None else treatment_linear)
    controls = np.column_stack([draw(e*lam, (phi, linear)) for e in exposures[:-treatments]])
    treated = np.column_stack([draw(e*lam*effects, treatment_dispersion) for e in exposures[-treatments:]])
    genes, indices, positions, start, cursor = [], [], np.zeros(total, dtype=np.int64), 1, 0
    for g, n in enumerate(sites_per_gene):
        positions[cursor:cursor+n] = start+spacing*np.arange(n)
        genes.append(dict(id=f"g{g:04d}", name=f"g{g:04d}", description="", chrom="chr1", strand="+",
                          blocks=[(start, start+spacing*max(n, 1)+3)]))
        indices.append(np.arange(cursor, cursor+n)); start += spacing*max(n, 1)+1000; cursor += n
    return controls, treated[:, 0] if treatments == 1 else treated, genes, indices, positions


def _json_ready(value):
    if isinstance(value, dict):return {str(k):_json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):return [_json_ready(v) for v in value]
    if isinstance(value, (np.bool_,)):return bool(value)
    if isinstance(value, np.integer):return int(value)
    if isinstance(value, (float, np.floating)):return float(value) if np.isfinite(value) else None
    return value


NONLINEAR_NORMALIZATIONS = ("betageom", "aBGC")   # resampling only: count-reshaping, no scalar factor


def load_normalization(normalized, normalization=None, nonlinear=False):
    """Locate one scalar method inside a YTAB `normalize` output (run directory or its method folder).

    nonlinear=True (resampling only) also accepts the beta-geometric corrections, which reshape counts.
    """
    normalized = Path(normalized).resolve()
    if (normalized/"method.json").is_file():
        folder, root = normalized, normalized.parent
    else:
        root = normalized
        folder = root/(normalization or DEFAULT_NORMALIZATION)
    if not (root/"run_manifest.json").is_file() or not (folder/"method.json").is_file():
        raise ValueError(f"Expected a YTAB normalization output with {folder.name}/method.json and run_manifest.json: {normalized}")
    method = json.loads((folder/"method.json").read_text())["id"]
    if normalization and method != normalization:
        raise ValueError(f"Normalized input is {method!r}, not the requested {normalization!r}")
    if method not in NORMALIZATIONS+(NONLINEAR_NORMALIZATIONS if nonlinear else ()):
        raise ValueError(f"{method!r} is not supported by the control-conditioned Gumbel; use one of {', '.join(NORMALIZATIONS)} "
                         "(scalar factors only: distribution-reshaping methods cannot be expressed as count-model exposures)")
    manifest = json.loads((root/"run_manifest.json").read_text())
    if manifest.get("status") != "success" or method not in manifest.get("methods", []):
        raise ValueError(f"Normalization run is incomplete or does not list {method!r}")
    with (root/"observations.tsv").open() as handle:observations = [r["observation_id"] for r in csv.DictReader(handle, delimiter="\t")]
    with (folder/"factors.tsv").open() as handle:factors = list(csv.DictReader(handle, delimiter="\t"))
    if [r["observation_id"] for r in factors] != observations or len(set(observations)) != len(observations):
        raise ValueError("Normalization factors and observations disagree")
    return dict(method=method, root=root, folder=folder, manifest=manifest, observations=observations,
                factors={r["observation_id"]:float(r["factor"]) if r["factor"] != "" else None for r in factors},
                groups={r["observation_id"]:r["normalization_group"] for r in factors})


def load_design(normalized, normalization, reference, output, controls, treatments, min_controls=2, min_treatments=1, stratified=False,
                nonlinear=False, check_roles=True, independent=False):
    """Validated controls/treatments from a YTAB scalar normalization output and its source catalogue.

    stratified=True (resampling only) allows several independent libraries, each with >=1 control and
    >=1 treatment sharing one normalization group; comparisons are then made within library.

    Shared by the control-conditioned Gumbel and resampling so both read normalized input identically.
    Returns raw integer counts and normalized counts (selected columns), factors, sites and provenance.
    nonlinear=True (resampling only) accepts betageom/aBGC output: factors are None, and the normalized
    counts are checked to keep every raw zero and the within-sample order of raw counts.
    check_roles=False (multi-condition tests only): conditions are user-assigned, so catalogue roles are not compared.
    independent=True (unstratified resampling only): the samples may come from different, independently made libraries
    (e.g. one library per strain, as in SATAY mutant screens); labels are permuted among samples within each site, as in
    TRANSIT's resampling. Differences in library composition then add to the noise, so the test is for library-level designs.
    """
    norm = load_normalization(normalized, normalization, nonlinear)
    catalogue = Path(norm["manifest"].get("catalogue", "")).resolve()
    reference, output = Path(reference).resolve(), Path(output).resolve()
    if any(output == x or x in output.parents or output in x.parents for x in (catalogue, norm["root"])) or (output.exists() and any(output.iterdir())):
        raise ValueError("Choose a new output directory, separate from the normalization and its catalogue")
    controls = list(controls); treatments = [treatments] if isinstance(treatments, str) else list(treatments)
    selected = controls+treatments
    if len(controls) < min_controls or len(treatments) < min_treatments or len(set(selected)) != len(selected):
        raise ValueError(f"Select >={min_controls} unique control observations and distinct treatment observations explicitly")
    if not (catalogue/"run_manifest.json").is_file():
        raise ValueError(f"Normalization source catalogue not found: {catalogue}")
    stale = [name for name, digest in norm["manifest"].get("source_hashes", {}).items() if sha256(catalogue/name) != digest]
    if stale or not norm["manifest"].get("source_hashes"):
        raise ValueError(f"Catalogue changed since normalization (checksum mismatch: {stale}); renormalize")
    source = json.loads((catalogue/"run_manifest.json").read_text())
    if source.get("status") != "success" or sha256(reference) != source.get("reference_sha256") or norm["manifest"].get("reference_sha256") != source.get("reference_sha256"):
        raise ValueError("Successful catalogue and checksum-identical reference required")
    with (catalogue/"observations.tsv").open() as handle:metadata = list(csv.DictReader(handle, delimiter="\t"))
    ids = [r["observation_id"] for r in metadata]
    if len(set(ids)) != len(ids):raise ValueError("Duplicate observation IDs in catalogue")
    missing = [x for x in selected if x not in ids or x not in norm["factors"]]
    if missing:raise ValueError(f"Observations missing from the catalogue or normalization: {missing}")
    chosen = [metadata[ids.index(x)] for x in selected]
    for field in ("project",) if stratified or independent else ("project", "library_id", "pool"):
        if len({r.get(field, "") for r in chosen}) != 1 or not chosen[0].get(field):
            raise ValueError(f"Controls and treatments must share one explicit {field}")
    if stratified:
        libraries = {r.get("library_id", "") for r in chosen}
        for library in libraries:
            members = [x for x, r in zip(selected, chosen) if r.get("library_id", "") == library]
            if not library or not set(members) & set(controls) or not set(members) & set(treatments):
                raise ValueError(f"Each library needs an explicit library_id with >=1 control and >=1 treatment: {library!r}")
            if len({norm["groups"][x] for x in members}) != 1:
                raise ValueError(f"Samples of library {library!r} must come from one normalization group")
    elif len({norm["groups"][x] for x in selected}) != 1:
        raise ValueError("Controls and treatments must come from one normalization group")
    control_roles = {r.get("role", "") for r in chosen[:len(controls)]}
    treatment_roles = {r.get("role", "") for r in chosen[len(controls):]}
    if check_roles and (len(control_roles) != 1 or len(treatment_roles) != 1 or control_roles == treatment_roles):
        raise ValueError("Controls must share one role, treatments another, and the roles must differ")
    raw = sparse.load_npz(catalogue/"raw_counts.npz").tocsc(); raw.sum_duplicates()
    if (raw.shape[1] != len(ids) or raw.shape[0] != source.get("universe_sites") or np.any(~np.isfinite(raw.data))
            or np.any(raw.data < 0) or np.any(raw.data != np.floor(raw.data))):
        raise ValueError("Raw integer matrix dimensions/counts disagree with catalogue")
    with gzip.open(catalogue/"sites.tsv.gz", "rt") as handle:
        sites = [(r["chromosome"], int(r["position_1based"])) for r in csv.DictReader(handle, delimiter="\t")]
    if len(sites) != raw.shape[0] or len(set(sites)) != len(sites) or any(p < 1 for _, p in sites):
        raise ValueError("Coordinate catalogue duplicates or dimensions invalid")
    counts = raw[:, [ids.index(x) for x in selected]].toarray().astype(np.int64)
    scaled = sparse.load_npz(norm["folder"]/"normalized_counts.npz").tocsc()
    if scaled.shape[0] != raw.shape[0]:raise ValueError("Normalized matrix and catalogue universe disagree")
    scaled = scaled[:, [norm["observations"].index(x) for x in selected]].toarray()
    if norm["method"] in NONLINEAR_NORMALIZATIONS:
        factors = None
        for j in range(counts.shape[1]):
            order = np.lexsort((scaled[:, j], counts[:, j]))
            if (np.any((scaled[:, j] > 0) != (counts[:, j] > 0)) or np.any(~np.isfinite(scaled[:, j]))
                    or np.any(np.diff(scaled[order, j]) < -1e-9*max(1., scaled[:, j].max()))):
                raise ValueError("Normalized counts do not preserve the raw zeros and count order of this catalogue")
    else:
        factors = np.array([norm["factors"][x] for x in selected])
        if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
            raise ValueError("Normalization factors must be finite and positive")
        if not np.allclose(scaled, counts*factors, rtol=1e-9, atol=1e-9):
            raise ValueError("Normalized counts are not raw catalogue counts times the reported factors; input is not a valid scalar normalization")
    return dict(norm=norm, catalogue=catalogue, output=output, controls=controls, treatments=treatments, selected=selected,
                libraries=[r.get("library_id", "") for r in chosen],
                chosen=chosen, counts=counts, scaled=scaled, factors=factors, sites=sites)


def run(normalized, gff, reference, output, controls, treatments, settings=None, normalization=None):
    """File-based runner on a YTAB scalar normalization output (ttr default, totreads or nzmean).

    The normalization's source catalogue supplies raw integer counts and coordinates; each observation's
    exposure is 1/factor relative to the control mean. The normalized matrix must equal raw*factor.
    """
    settings = settings or Settings(); settings.validate()
    design = load_design(normalized, normalization, reference, output, controls, treatments)
    norm, catalogue, output, controls, treatments, selected, chosen, counts, factors, sites = (design[k] for k in
        ("norm", "catalogue", "output", "controls", "treatments", "selected", "chosen", "counts", "factors", "sites"))
    gff = Path(gff).resolve()
    J = len(controls); inverse = 1/factors; exposures = inverse/inverse[:J].mean()
    genes = annotations(gff)
    positions = np.array([p for _, p in sites], dtype=np.int64)
    rows, site_rows, diagnostics = analyze(counts[:, :J], counts[:, J:], exposures, genes, gene_sites(genes, sites, settings), positions, settings)
    for row in site_rows:row["chromosome"] = sites[row["site_index"]][0]
    output.mkdir(parents=True, exist_ok=True)
    leading = ["ORF", "Name", "Description", "n_eligible", "n_depleted", "r", "run_span_bp", "p_value", "q_bh", "call", "status"]
    write_table(output/"genes.tsv", leading+list(dict.fromkeys(k for r in rows for k in r if k not in leading)), rows)
    first = ["site_index", "chromosome", "position", "genes"]
    write_table(output/"sites.tsv.gz", first+[k for k in site_rows[0] if k not in first] if site_rows else ["site_index"], site_rows)
    manifest = dict(status="success", analysis=MODEL, algorithm_revision=REVISION,
        controls=chosen[:J], treatments=chosen[J:], parameters=asdict(settings),
        normalization=dict(method=norm["method"], input=str(norm["folder"]), catalogue=str(catalogue),
            factors={x:float(f) for x, f in zip(selected, factors)}, group=norm["groups"][selected[0]],
            hashes={name:sha256(norm["folder"]/name) for name in ("normalized_counts.npz", "factors.tsv", "method.json")}
                   | dict(run_manifest=sha256(norm["root"]/"run_manifest.json")),
            check="normalized_counts == raw catalogue counts * factor (verified)"),
        exposures=dict(rule="exposure_j = (1/factor_j) / mean over controls of (1/factor); counts are the catalogue's raw integers",
                       values={x:float(e) for x, e in zip(selected, exposures)},
                       total_site_support={x:float(t) for x, t in zip(selected, counts.sum(0))}),
        model=dict(combination="replicate tracks summed within condition; combined control defines the null universe and baseline, combined treatment is compared",
            eligibility=f"combined control count >= {settings.min_control_reads}" + ("; and a zero treatment count would be flagged" if settings.eligibility == "detectable" else ""),
            background="C_ij ~ NB(e_j*lambda_i), Var = mu*(1+linear_C) + phi_C*mu^2" + (" (linear term fixed at 0)" if settings.dispersion_model == "quadratic" else "") + "; log lambda_i ~ N(m, tau^2); m, tau, phi_C, linear_C fitted on control replicate tracks with zero-truncated (selection-conditioned) marginal likelihood",
            treatment_dispersion=("phi_T, linear_T fitted the same way on the treatment replicate tracks (replicate disagreement only; shared effects are not dispersion)"
                                  if settings.treatment_dispersion == "replicates" else "treatment dispersion assumed equal to control dispersion (explicit setting)"),
            combined_tracks="sum of independent replicates approximated by NB with the same mean, linear term unchanged and phi*sum(e^2)/sum(e)^2",
            encoding="D_i = 1[F_i(T_i) <= alpha_site], F_i = predictive NB CDF of the combined treatment given the combined control (adaptive Gauss-Hermite posterior)",
            statistic="R_g = longest run of consecutive D_i=1 across eligible sites in genomic order; last tied maximum; span = end-start+1",
            null="per-gene seeded count-level simulation from the same predictive, identical encoding; discrete rounded/truncated Gumbel fitted by ML on calibration runs",
            adequacy=("validation-run CDF within a two-sample KS bound covering calibration and validation sampling; |log fitted/exact tail| <= log(tail_tolerance) + z*SE (delta-method Monte Carlo SE of the fitted log tail), or both below fdr/family_size"
                      if settings.inference == "gumbel-mc" else
                      "not gating: p is the exact conditional tail P(R>=r) from the per-site null flag probabilities; Gumbel fit and checks reported in gumbel_check"),
            inference=settings.inference,
            multiple_testing="BH over genes with >=min_sites eligible sites; uncalibrated or run-free genes enter as p=1"),
        genes=len(genes), diagnostics=diagnostics,
        reference_sha256=sha256(reference),
        source_hashes={str(p):sha256(p) for p in [catalogue/"raw_counts.npz", catalogue/"observations.tsv",
                       catalogue/"sites.tsv.gz", catalogue/"run_manifest.json", gff, reference]},
        implementation_sha256=sha256(__file__),
        limitations=["Experimental run-based depletion evidence; not validated essentiality, fitness or a published Hermes method",
            "Dispersion is shared across sites within each condition as a two-term mean-variance curve; site-specific dispersion is not modelled",
            "Treatment dispersion from two replicates can absorb replicate-specific (non-shared) effects; it cannot separate them from noise",
            "Summed replicate tracks use a moment-matched NB approximation",
            "Sites are conditionally independent given lambda; shared clone, positional or depth-estimation dependence is not modelled",
            "Plug-in background hyperparameters; their estimation uncertainty is not propagated into site flags",
            "Lognormal abundance prior is an assumption; see posterior_predictive diagnostics",
            "Scalar normalization exposures (ttr/totreads/nzmean) assume most sites are unaffected by treatment",
            "Discrete Gumbel is an approximation; inadequate fits give no significance, not evidence of no depletion",
            "Overlapping genes share sites; BH assumes valid p-values and tolerable dependence",
            "Fitness, classifier and UI are unchanged"])
    (output/"run_manifest.json").write_text(json.dumps(_json_ready(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
