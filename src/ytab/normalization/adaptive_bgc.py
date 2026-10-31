"""Adaptive beta-geometric correction with explicit search and tail conventions.

Retains TRANSIT's positive-count estimators and maximum discrepancy score
(AdaptiveBGCNorm; the beta-geometric correction is DeJesus & Ioerger 2016).
Analytic CDF default; upstream-style seeded simulation is an optional engine.
"""
import hashlib

import numpy as np
from scipy import sparse
from scipy.stats import geom, trim_mean

from .betageom import beta_geometric_cdf
from .ttr import ttr_parameters

DEFAULT_FRACTIONS = tuple(i/100 for i in range(0,31,2))


def fit_candidate(positive, fraction):
    z = np.sort(np.asarray(positive,dtype=float))
    if z.ndim != 1 or len(z)<2 or np.any(~np.isfinite(z)) or np.any(z<=0) or not 0<=fraction<.5:
        raise ValueError("Adaptive BGC requires at least two finite positive counts and 0<=fraction<0.5")
    rho = float(1/trim_mean(z,fraction))
    logs = np.log1p(-np.arange(len(z))/len(z))
    A = float(np.dot(logs,logs)/np.dot(z,logs))
    divisor = np.exp(A)+rho-1
    if divisor == 0:raise ValueError("Singular concentration estimator")
    kappa = float((2*np.exp(A)-1)/divisor)
    alpha,beta = kappa*rho,kappa*(1-rho)
    if not (0<rho<1 and np.isfinite(kappa) and kappa>0 and alpha>0 and beta>0):
        raise ValueError("Invalid beta-geometric shape parameters; no fallback distribution")
    return dict(rho=rho,kappa=kappa,alpha=alpha,beta=beta,A=A)


def mapped_counts(unique, fit, engine, seed, draws, fingerprint, fraction, phase):
    if engine=="analytic":
        probability=beta_geometric_cdf(unique,fit["alpha"],fit["beta"])
    else:
        token=f"{seed}:{fingerprint}:{fraction:.17g}:{phase}"
        child=int.from_bytes(hashlib.sha256(token.encode()).digest()[:8],"little")
        rng=np.random.default_rng(child)
        p=rng.beta(fit["alpha"],fit["beta"],size=draws)
        if np.any(p<=0) or np.any(p>1):raise ValueError("Invalid simulated beta probabilities")
        simulated=np.sort(rng.geometric(p))
        if np.any(simulated<1):raise ValueError("Simulated geometric counts overflowed")
        probability=np.searchsorted(simulated,unique,side="right")/draws
    bounded=np.clip(probability,np.nextafter(0.,1.),np.nextafter(1.,0.))
    if engine=="analytic":
        # F_BG(1)=E[p]=rho exactly. Log-beta subtraction can place that
        # percentile just above a discrete geometric step and incorrectly
        # turn every one-read observation into two reads.
        probability[unique==1]=fit["rho"]
        bounded=np.clip(probability,np.nextafter(0.,1.),np.nextafter(1.,0.))
        inverse=np.log1p(-bounded)/np.log1p(-fit["rho"])
        rounded=np.rint(inverse)
        inverse=np.where(np.isclose(inverse,rounded,rtol=1e-12,atol=1e-10),rounded,inverse)
        values=np.maximum(np.ceil(inverse),1.)
        # Jensen's inequality implies F_BG(x)<=F_Geom,rho(x) for integer x,
        # so the exact mapping cannot exceed the original positive count.
        values=np.minimum(values,unique)
    else:
        values=np.maximum(geom.ppf(bounded,fit["rho"]),1.)
    if np.any(~np.isfinite(values)):raise ValueError("Nonfinite corrected counts")
    return values,probability


# "endpoint": scoring quantiles r/n, the convention of the source implementation; "transit" (its earlier name) is accepted.
TAIL_ALIASES={"transit":"endpoint"}


def discrepancy(ordered, tail="endpoint"):
    tail=TAIL_ALIASES.get(tail,tail)
    y=np.asarray(ordered,dtype=float)
    if len(y)<2 or np.any(~np.isfinite(y)) or np.any(y<=0) or tail not in {"endpoint","midpoint"}:
        raise ValueError("Invalid corrected counts or tail convention")
    n=len(y)
    q=np.arange(1,n+1)/n if tail=="endpoint" else (np.arange(n)+.5)/n
    theoretical=geom.ppf(q,1/y.mean())
    infinite=np.isinf(theoretical)
    if np.any(infinite):
        finite=theoretical[~infinite]
        if not len(finite):raise ValueError("No finite target quantiles")
        theoretical[infinite]=finite.max()+100
    if np.any(theoretical<=0) or np.any(~np.isfinite(theoretical)):
        raise ValueError("Invalid theoretical quantiles")
    return float(np.max((y-theoretical)**2/theoretical))


def correct(matrix,groups,cdf="analytic",seed=0,draws=200000,post_ttr=True,tail="endpoint",fractions=None,trim=.05,target=100):
    tail=TAIL_ALIASES.get(tail,tail)
    matrix=sparse.csc_matrix(matrix,dtype=float,copy=True)
    matrix.eliminate_zeros()
    grid=DEFAULT_FRACTIONS if fractions is None else tuple(fractions)
    if (cdf not in {"analytic","simulation"} or tail not in {"endpoint","midpoint"} or
            not isinstance(seed,(int,np.integer)) or seed<0 or not isinstance(draws,(int,np.integer)) or draws<1 or
            not grid or len(set(grid))!=len(grid) or any(not np.isfinite(f) or not 0<=f<.5 for f in grid) or
            matrix.shape[1]!=len(groups) or not matrix.shape[0] or np.any(~np.isfinite(matrix.data)) or
            np.any(matrix.data<0) or np.any(matrix.data!=np.floor(matrix.data))):
        raise ValueError("Invalid adaptive BGC settings, fraction grid or raw integer matrix")
    output=matrix.copy()
    details=[]
    for j in range(matrix.shape[1]):
        positive=matrix.data[matrix.indptr[j]:matrix.indptr[j+1]]
        if len(positive)<2:raise ValueError(f"Sample {j}: adaptive BGC requires at least two occupied coordinates")
        unique,inverse,frequency=np.unique(positive,return_inverse=True,return_counts=True)
        fingerprint=hashlib.sha256(np.sort(positive).astype("<f8").tobytes()).hexdigest()
        candidates=[]
        for fraction in sorted(grid):
            row=dict(trim_fraction=float(fraction),status="invalid",reason="",rho="",kappa="",alpha="",beta="",A="",score="",selected=False)
            try:
                fit=fit_candidate(positive,fraction)
                row.update(fit)
                values,_=mapped_counts(unique,fit,cdf,seed,draws,fingerprint,fraction,"candidate")
                ordered=np.repeat(values,frequency)
                score=discrepancy(ordered,tail)
                if not np.isfinite(score):raise ValueError("Nonfinite candidate discrepancy")
                row.update(status="valid",score=score)
            except (ValueError,FloatingPointError,OverflowError) as error:
                row["reason"]=str(error)
            candidates.append(row)
        valid=[r for r in candidates if r["status"]=="valid"]
        if not valid:
            raise ValueError(f"Sample {j}: no valid adaptive BGC candidate; "+"; ".join(sorted({r['reason'] for r in candidates})))
        chosen=min(valid,key=lambda r:(r["score"],r["trim_fraction"]))
        chosen["selected"]=True
        values,probability=mapped_counts(unique,chosen,cdf,seed,draws,fingerprint,chosen["trim_fraction"],"final")
        mapped=values[inverse]
        scale=ttr_parameters(mapped,matrix.shape[0],trim,target)["factor"] if post_ttr else 1.
        output.data[output.indptr[j]:output.indptr[j+1]]=mapped*scale
        top=max(1,int(np.ceil(len(positive)*.01)))
        details.append(dict(abgc_selected_trim=chosen["trim_fraction"],abgc_score=chosen["score"],abgc_rho=chosen["rho"],
            abgc_kappa=chosen["kappa"],abgc_alpha=chosen["alpha"],abgc_beta=chosen["beta"],abgc_A=chosen["A"],
            abgc_valid_candidates=len(valid),abgc_invalid_candidates=len(candidates)-len(valid),abgc_tail=tail,abgc_cdf=cdf,
            abgc_seed=seed,abgc_draws=draws if cdf=="simulation" else 0,abgc_post_ttr=post_ttr,abgc_post_ttr_factor=scale,
            abgc_post_ttr_trim_each_tail=trim if post_ttr else "",abgc_post_ttr_target=target if post_ttr else "",
            abgc_lower_probability_safeguard_sites=int(frequency[probability<=0].sum()),
            abgc_upper_probability_safeguard_sites=int(frequency[probability>=1].sum()),
            abgc_raw_top_1pct_support_fraction=float(np.sort(positive)[-top:].sum()/positive.sum()),
            abgc_corrected_top_1pct_support_fraction=float(np.sort(mapped)[-top:].sum()/mapped.sum()),
            _abgc_candidates=candidates))
    return output.tocsr(),details
