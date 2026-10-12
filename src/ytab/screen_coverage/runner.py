"""JSON worker for the Shiny Screen Coverage tab; never changes project data."""
import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import binom, poisson

from .audited_library import load_audited_library, BASES, MODEL_LABELS
from .statistics import (poisson_coverage, colonies_for_redundancy_target,
    weighted_poisson_coverage, exact_empirical_coverage, empirical_target,
    coverage_ceiling, monte_carlo_uniform, monte_carlo_hermes, empirical_rarefaction, monte_carlo_alleles)


def calculate_technical(request, cache_dir=None):
    audit_dir=request.get("audit_dir",str(Path(__file__).resolve().parents[3]/"docs/codex/parent_union_audit"))
    library = load_audited_library(audit_dir,request.get("selected_sources"),request.get("universe","informative"),
                                  request.get("basis","sites"),request.get("conditioning","physical"))
    allele_genes=library.pop("_allele_genes");allele_p=library.pop("_allele_p")
    n = len(library["genes"])
    p = np.array([g["p"] for g in library["genes"]])
    for name, default in (("k",10000),("r",1),("iterations",1000),("seed",2026),("max_k",100000)):
        value = request.get(name,default)
        if not isinstance(value,(int,float)) or not math.isfinite(value) or value != int(value):
            raise ValueError(f"{name} must be a finite integer")
    k, r = int(request.get("k", 10000)), int(request.get("r", 1))
    target = float(request.get("target", .95))
    if k < 0 or k > 10**12 or not 1 <= r <= 10000 or not 0 < target < 1:
        raise ValueError("Use k between 0 and 10^12, r between 1 and 10000 and a target between 0 and 100%")
    iterations, seed = int(request.get("iterations", 1000)), int(request.get("seed", 2026))
    if not 2 <= iterations <= 10000:
        raise ValueError("Iterations must be between 2 and 10000")
    thresholds = sorted({1, 2, 3, 5, r})
    targets = sorted({.90, .95, .99, .999, target})
    inverse, curves = [], []
    ps = {basis: np.array([g[f"p_{basis}"] for g in library["genes"]]) for basis in BASES}
    models = {"Equal-probability Poisson": (None, False)}
    for b in BASES:
        models[MODEL_LABELS[b]+" — Poisson"]=(ps[b],False)
        models[MODEL_LABELS[b]+" — exact"]=(ps[b],True)
    for method, (probs, exact) in models.items():
        for rr in thresholds:
            for t in targets:
                res = (dict(colonies=colonies_for_redundancy_target(n, t, rr), status="achievable") if probs is None
                    else empirical_target(probs, t, rr, exact))
                kk = res["colonies"]
                coverage = None if kk is None else (poisson_coverage(kk, n, rr) if probs is None else
                    (exact_empirical_coverage if exact else weighted_poisson_coverage)(kk, probs, rr))
                inverse.append(dict(method=method, target=t, r=rr, depth=None if kk is None else kk/n,
                    represented=None if coverage is None else n*coverage,
                    missed=None if coverage is None else n*(1-coverage), **res))
    max_k = int(request.get("max_k", max(k * 2, colonies_for_redundancy_target(n, .999, r))))
    if not 1 <= max_k <= 10**12:
        raise ValueError("Maximum curve size must be between 1 and 10^12")
    practical_sizes=[1000,2500,5000,7500,10000,12500,15000,20000]
    grid = np.unique(np.r_[np.linspace(0, max_k, 100).astype(np.int64), k, practical_sizes,
        [x["colonies"] for x in inverse if x["colonies"] is not None and x["colonies"] <= max_k and x["r"] == r]])
    for method, (probs, exact) in models.items():
        for rr in thresholds:
            for kk in grid:
                cov = poisson_coverage(kk, n, rr) if probs is None else (exact_empirical_coverage if exact else weighted_poisson_coverage)(int(kk), probs, rr)
                curves.append(dict(method=method, k=int(kk), r=rr, mean=cov, represented=n*cov, missed=n*(1-cov)))
    planned = []
    for method, (probs, exact) in models.items():
        values = {rr: poisson_coverage(k, n, rr) if probs is None else (exact_empirical_coverage if exact else weighted_poisson_coverage)(k, probs, rr) for rr in {1,2,3,r}}
        entry = dict(method=method, k=k, depth=k/n, represented=n*values[1], missed=n*(1-values[1]),
                     coverage=values[1], at_least_2=values[2], at_least_3=values[3], selected_r=values[r])
        for t in (.95,.99):
            found = next(x for x in inverse if x["method"] == method and x["r"] == r and x["target"] == t)
            entry[f"additional_{round(t*100)}"] = max(0, found["colonies"] - k) if found["colonies"] is not None else found["status"]
        planned.append(entry)
    for gene in library["genes"]:
        prob = gene["p"]
        gene.update(expected_first=None if prob == 0 else 1/prob,
            recovery=float(binom.sf(0,k,prob)), recovery_2=float(binom.sf(1,k,prob)))
    rare_counts = [g["informative_sites"] for g in library["genes"]]
    catalogue_size = library["unique_sites"]
    simulations, rarefaction, alleles = [], [], []
    seed_map={"uniform":seed,"sites":seed+1,"h2o2":seed+2,"zn":seed+3,"rarefaction":seed+4,
              "alleles_sites":seed+5,"alleles_h2o2":seed+6,"alleles_zn":seed+7}
    allele_note="Distinct-allele simulations use practical colony counts and the planned k up to 1,000,000; analytical and colony-redundancy calculations are not limited by this cap."
    if request.get("simulate", False):
        engine_hash = hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name("statistics.py").read_bytes()).hexdigest()
        sim_request = dict(engine=engine_hash,ps={b:v.tolist() for b,v in ps.items()}, n=n, max_k=max_k, k=k, r=thresholds, iterations=iterations, seed=seed,
                           catalogue=rare_counts, other=library["rarefaction_other_sites"],sources=library["sources"],
                           allele_hash=hashlib.sha256(allele_genes.tobytes()+b"".join(v.tobytes() for v in allele_p.values())).hexdigest())
        key = hashlib.sha256(json.dumps(sim_request, sort_keys=True).encode()).hexdigest()
        cache = Path(cache_dir)/f"screen-coverage-{key}.json" if cache_dir else None
        if cache and cache.exists():
            cached = json.loads(cache.read_text())
            simulations, rarefaction, alleles = cached["simulations"], cached["rarefaction"], cached["alleles"]
        else:
            sizes = sorted(set(np.linspace(0,max_k,31).astype(int).tolist()+[k]+practical_sizes))
            simulations.extend(dict(method="Uniform Monte Carlo",**row) for row in monte_carlo_uniform(n,sizes,thresholds,iterations,seed))
            for b in BASES:
                simulations.extend(dict(method=MODEL_LABELS[b]+" — Monte Carlo",**row)
                    for row in monte_carlo_hermes(ps[b],sizes,thresholds,iterations,seed_map[b]))
                allele_sizes=sorted(set(practical_sizes+([k] if k<=1_000_000 else [])))
                alleles.extend(dict(method=MODEL_LABELS[b],**row) for row in monte_carlo_alleles(
                    allele_genes,allele_p[b],n,allele_sizes,iterations,seed_map["alleles_"+b]))
            rarefaction = empirical_rarefaction(rare_counts,np.unique(np.linspace(0,catalogue_size,31).astype(int)),
                thresholds,iterations,seed_map["rarefaction"],library["rarefaction_other_sites"])
            if cache:
                cache.parent.mkdir(parents=True,exist_ok=True)
                cache.write_text(json.dumps(dict(simulations=simulations,rarefaction=rarefaction,alleles=alleles)))
        for method in ("Uniform Monte Carlo",*[MODEL_LABELS[b]+" — Monte Carlo" for b in BASES],"Site catalogue rarefaction"):
            rows = rarefaction if method == "Site catalogue rarefaction" else [x for x in simulations if x["method"] == method]
            for rr in thresholds:
                ceiling = (float(np.mean(np.array(rare_counts)>=rr)) if method == "Site catalogue rarefaction" else
                    1. if method == "Uniform Monte Carlo" else coverage_ceiling(ps[next(b for b in BASES if method==MODEL_LABELS[b]+" — Monte Carlo")]))
                for t in targets:
                    reached = [x for x in rows if x["r"] == rr and x["mean"] >= t]
                    kk = min(x["k"] for x in reached) if reached else None
                    status = ("First simulated grid point reaching target; not a precise inverse" if reached else
                        "Not achievable under current library representation" if t > ceiling else "Not reached in simulated range")
                    inverse.append(dict(method=method,target=t,r=rr,colonies=kk,status=status,
                        unit="distinct sites" if method == "Site catalogue rarefaction" else "colonies"))
    lam = k/n
    upper = int(poisson.ppf(.999999,lam))
    # Only the visible PMF range is limited, never the calculation of coverage.
    start = max(0,int(poisson.ppf(.000001,lam)))
    xs = np.arange(start,upper+1) if upper-start < 500 else np.unique(np.linspace(start,upper,500).astype(int))
    pmf = [dict(x=int(x), probability=float(poisson.pmf(x,lam))) for x in xs]
    return dict(library=library, settings=dict(request, k=k, r=r, target=target, max_k=max_k, iterations=iterations,seed=seed,
        universe=request.get("universe","informative"),basis=request.get("basis","sites"),conditioning=request.get("conditioning","physical")),
        software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        n=n, ceiling=coverage_ceiling(p), useful_sampling_mass=float(p.sum()), lambda_value=lam,
        inverse=inverse,curves=curves,planned=planned,simulations=simulations,rarefaction=rarefaction,pmf=pmf,
        alleles=alleles,allele_simulation_note=allele_note,seed_map=seed_map,
        model_ceilings={b:coverage_ceiling(ps[b]) for b in BASES},
        pmf_display_note="PMF omits extreme tails outside 0.000001–0.999999 quantiles; very wide ranges use at most 500 discrete x values.")


def calculate(request, cache_dir=None):
    # Preserve earlier analytical/rarefaction tools in Statistical Details.
    # They never supply the pool-stratified headline results.
    from .pool_runner import calculate_pool_screen
    technical = dict(request, simulate=bool(request.get("technical_simulations", False)))
    result = calculate_technical(technical, cache_dir)
    result["pool_screen"] = calculate_pool_screen(request, cache_dir)
    result["settings"].update(simulate=bool(request.get("simulate", True)))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request",required=True)
    parser.add_argument("--output",required=True)
    parser.add_argument("--cache-dir")
    args = parser.parse_args()
    try:
        result = calculate(json.loads(Path(args.request).read_text()), args.cache_dir)
    except Exception as exc:
        Path(args.output).write_text(json.dumps({"error":str(exc)}))
        raise
    Path(args.output).write_text(json.dumps(result,allow_nan=False))


if __name__ == "__main__":
    main()
