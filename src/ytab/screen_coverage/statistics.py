"""Screen planning mathematics. All coverage values are fractions, not percents.

Multinomial draws model colonies; hypergeometric draws model a finite site
catalogue. Neither sequencing reads nor distinct sites are colony counts.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson, binom, hypergeom


def _integer(value, name, minimum=0):
    if not np.isfinite(value) or value != int(value) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _probabilities(p):
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or not len(p) or not np.all(np.isfinite(p)) or np.any(p < 0) or np.any(p > 1) or p.sum() > 1 + 1e-10:
        raise ValueError("Gene probabilities must be finite, nonnegative and sum to at most one")
    return p


def poisson_probability(x, lam):
    if not np.isfinite(lam) or lam < 0:
        raise ValueError("lambda must be finite and nonnegative")
    return poisson.pmf(x, lam)


def redundancy_probability(lam, r=1):
    _integer(r, "r", 1)
    lam = np.asarray(lam, dtype=float)
    if np.any(~np.isfinite(lam)) or np.any(lam < 0):
        raise ValueError("lambda must be finite and nonnegative")
    return poisson.sf(r - 1, lam)


def poisson_coverage(k, n, r=1):
    _integer(n, "N", 1)
    return float(redundancy_probability(float(k) / n, r))


def colonies_for_target_coverage(n, target):
    return colonies_for_redundancy_target(n, target, 1)


def colonies_for_redundancy_target(n, target, r=1):
    _integer(n, "N", 1)
    _integer(r, "r", 1)
    if not 0 <= target < 1:
        raise ValueError("Target must be >=0 and <1")
    if target == 0:
        return 0
    if r == 1:
        return math.ceil(-n * math.log1p(-target))
    high = float(r)
    while redundancy_probability(high, r) < target:
        high *= 2
    return math.ceil(n * brentq(lambda lam: redundancy_probability(lam, r) - target, 0, high, xtol=1e-12))


def empirical_gene_probabilities(abundance, total=None):
    a = np.asarray(abundance, dtype=float)
    if a.ndim != 1 or not len(a) or np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError("Invalid abundance")
    total = a.sum() if total is None else float(total)
    if not np.isfinite(total) or total < a.sum() - 1e-8 or total < 0:
        raise ValueError("Denominator cannot be smaller than screenable abundance")
    return _probabilities(a / total if total else np.zeros_like(a))


def coverage_ceiling(p):
    p = _probabilities(p)
    return float(np.count_nonzero(p) / len(p))


def weighted_poisson_coverage(k, p, r=1):
    p = _probabilities(p)
    return float(np.mean(redundancy_probability(float(k) * p, r)))


def exact_empirical_coverage(k, p, r=1):
    k = _integer(k, "k")
    _integer(r, "r", 1)
    return float(binom.sf(r - 1, k, _probabilities(p)).mean())


def empirical_target(p, target, r=1, exact=False, maximum=10**12):
    p = _probabilities(p)
    _integer(r, "r", 1)
    if not 0 <= target < 1:
        raise ValueError("Target must be >=0 and <1")
    if target == 0:
        return {"colonies": 0, "status": "achievable"}
    ceiling = coverage_ceiling(p)
    if target > ceiling or (target == ceiling and not (exact and np.all(p[p > 0] == 1))):
        return {"colonies": None, "status": "Not achievable under current library representation"}
    fn = exact_empirical_coverage if exact else weighted_poisson_coverage
    lo, hi = 0, max(1, len(p))
    while fn(hi, p, r) < target and hi < maximum:
        hi = min(hi * 2, maximum)
    if fn(hi, p, r) < target:
        return {"colonies": None, "status": f"Exceeds numerical search limit ({maximum:,} colonies)"}
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if fn(mid, p, r) >= target:
            hi = mid
        else:
            lo = mid
    return {"colonies": hi, "status": "achievable"}


def _summary(values, k, r, n):
    lo, hi = np.quantile(values, [0.025, 0.975])
    return dict(k=int(k), r=int(r), mean=float(values.mean()), sd=float(values.std(ddof=1)),
                lower=float(lo), upper=float(hi), represented=float(n * values.mean()),
                missed=float(n * (1 - values.mean())))


def monte_carlo_hermes(p, sizes, r=(1, 2, 3), iterations=1000, seed=2026):
    """Exact joint multinomial occupancy; aggregated insertion-to-gene draws.

    The final category is unscreenable mass. Incremental trajectories preserve
    monotonicity and avoid allocating one array element per picked colony.
    """
    p = _probabilities(p)
    _integer(iterations, "iterations", 2)
    rng = np.random.default_rng(_integer(seed, "seed"))
    thresholds = sorted(set(_integer(x, "r", 1) for x in np.atleast_1d(r)))
    sizes = sorted(set(_integer(x, "k") for x in sizes))
    full = np.r_[p, max(0., 1 - p.sum())]
    full /= full.sum()
    summaries = {rr: np.empty((iterations,len(sizes))) for rr in thresholds}
    # Bound working memory even for the advanced 10,000-iteration setting.
    for start in range(0,iterations,128):
        stop = min(start+128,iterations)
        counts = np.zeros((stop-start,len(full)),dtype=np.int64)
        previous = 0
        for j,k in enumerate(sizes):
            counts += rng.multinomial(k-previous,full,size=stop-start)
            for rr in thresholds:
                summaries[rr][start:stop,j] = (counts[:,:-1]>=rr).mean(axis=1)
            previous = k
    return [_summary(summaries[rr][:,j],k,rr,len(p)) for j,k in enumerate(sizes) for rr in thresholds]


def monte_carlo_uniform(n, sizes, r=(1, 2, 3), iterations=1000, seed=2026):
    _integer(n, "N", 1)
    return monte_carlo_hermes(np.full(n, 1 / n), sizes, r, iterations, seed)


def monte_carlo_alleles(allele_genes, allele_probabilities, n, sizes, iterations=1000, seed=2026):
    """Exact with-replacement colony draws, counting distinct observed alleles.

    Binomial thinning removes non-target draws without changing the sampling
    distribution. A coordinate's repeat draws never increment its gene count.
    Memory is bounded by the catalogue plus at most 100,000 draws at a time.
    """
    _integer(n,"N",1);_integer(iterations,"iterations",2)
    gids=np.asarray(allele_genes,dtype=int)
    p=np.asarray(allele_probabilities,dtype=float)
    if len(p):
        _probabilities(p)
    if gids.shape!=p.shape or np.any(gids<0) or np.any(gids>=n):
        raise ValueError("Invalid allele-to-gene mapping")
    rng=np.random.default_rng(_integer(seed,"seed"))
    sizes=sorted(set(_integer(k,"k") for k in sizes))
    if sizes and max(sizes)>1_000_000:
        raise ValueError("Distinct-allele simulations support at most 1,000,000 colonies; analytical models remain available above this limit")
    mass=min(1.,float(p.sum()))
    positive=p>0;gids=gids[positive];p=p[positive]
    cdf=np.cumsum(p/mass) if mass else np.array([])
    if len(cdf):cdf[-1]=1.
    results={rr:np.empty((iterations,len(sizes))) for rr in (1,2,3)}
    for it in range(iterations):
        seen=np.zeros(len(p),dtype=bool);counts=np.zeros(n,dtype=int);previous=0
        for j,k in enumerate(sizes):
            target_draws=rng.binomial(k-previous,mass)
            while target_draws:
                chunk=min(target_draws,100000)
                ids=np.unique(np.searchsorted(cdf,rng.random(chunk),side="right"))
                new=ids[~seen[ids]];seen[new]=True
                counts+=np.bincount(gids[new],minlength=n)
                target_draws-=chunk
            for rr in results:results[rr][it,j]=np.mean(counts>=rr)
            previous=k
    return [_summary(results[rr][:,j],k,rr,n) for j,k in enumerate(sizes) for rr in results]


def empirical_rarefaction(site_counts, sizes, r=(1, 2, 3), iterations=1000, seed=2026, other_sites=0):
    """Without-replacement catalogue draws; unit is distinct insertion sites.

    Counts represent disjoint genes; an extra category retains sites outside
    the chosen universe. This does NOT infer viable-colony sampling.
    """
    a = np.asarray([_integer(x, "site count") for x in site_counts], dtype=np.int64)
    if not len(a):
        raise ValueError("Empty universe")
    _integer(iterations, "iterations", 2)
    full = np.r_[a, _integer(other_sites, "other sites")]
    total = int(full.sum())
    rng = np.random.default_rng(_integer(seed, "seed"))
    sizes = sorted(set(_integer(x, "sites") for x in sizes))
    if sizes and max(sizes) > total:
        raise ValueError("Cannot sample more sites than the catalogue contains")
    thresholds = sorted(set(_integer(x, "r", 1) for x in np.atleast_1d(r)))
    values = {rr: np.empty((iterations, len(sizes))) for rr in thresholds}
    # One random permutation per replicate; catalogue sizes are bounded by
    # observed sites, never expanded to sequencing-read or cell counts.
    catalogue = np.repeat(np.arange(len(full)), full)
    for it in range(iterations):
        shuffled = rng.permutation(catalogue)
        counts = np.zeros(len(full), dtype=int)
        previous = 0
        for j, k in enumerate(sizes):
            counts += np.bincount(shuffled[previous:k], minlength=len(full))
            for rr in thresholds:
                values[rr][it, j] = np.mean(counts[:-1] >= rr)
            previous = k
    result = []
    for j, k in enumerate(sizes):
        for rr in thresholds:
            row = _summary(values[rr][:, j], k, rr, len(a))
            row["exact_mean"] = float(hypergeom.sf(rr - 1, total, a, k).mean()) if total else 0.
            result.append(row)
    return result
