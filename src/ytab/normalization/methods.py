"""Shared method descriptions and explicit, opt-in site-normalization choices.

This registry is also the future normalization UI contract. Planned methods are
never silently dispatched or advertised as ready. Existing workflows are separate.
"""
from dataclasses import asdict, dataclass
import csv
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from .ttr import sha256, ttr_parameters, write_table
from .quantile import normalize_positive_quantiles
from .betageom import correct as correct_betageom
from .adaptive_bgc import correct as correct_adaptive_bgc, DEFAULT_FRACTIONS
from .zinfnb import fit as fit_zinfnb, LOG_BOUNDS, START_SIZES
from .midlc import thin as thin_midlc, SEED as MIDLC_SEED

# The methods adapted here are published in DeJesus et al. 2015 (TRANSIT, PLoS Comput Biol 11:e1004401) and, for the
# beta-geometric corrections, DeJesus & Ioerger 2016 (J Bioinform Comput Biol 14:1642004). See docs/references.md.
TRANSIT_SOURCE = "https://github.com/ioerger/transit/blob/50e23d6bec2605a24835fbd67f33691e330f3361/src/pytransit/norm_tools.py"


@dataclass(frozen=True)
class Method:
    id: str
    name: str
    status: str
    input_level: str
    description: str
    formula: str
    may_suit: str
    caution: str
    diagnostics: tuple
    source: str = TRANSIT_SOURCE


METHODS = {m.id:m for m in (
    Method("nonorm","Unnormalized baseline","available","raw insertion-site counts",
      "Preserves raw support and zeros without scaling.","factor_j = 1",
      "Baseline for evaluating what normalization changes.",
      "Unequal sequencing support remains; not an automatic recommendation for comparisons.",("support depth","occupied sites")),
    Method("ttr","Trimmed total reads (TTR)","available","raw insertion-site counts",
      "Scales by site density times the trimmed positive-site mean; defaults trim 5% per tail and target 100.",
      "factor_j = target / (density_j * trimmed_positive_mean_j)",
      "Evaluation of robust scalar scaling when a few high-count sites influence read totals.",
      "Does not correct missing mutants or all distribution skew; density can include genuine library-complexity differences.",
      ("density","trimmed positive-site mean","outlier influence","matched-control agreement")),
    Method("totreads","Total-read scaling","available","raw insertion-site counts",
      "Scales each library's total reads to the mean total within the selected normalization group.",
      "factor_j = mean_group(total_support) / total_support_j",
      "A simple baseline when support-depth differences dominate and a few sites do not dominate totals.",
      "Sensitive to extreme counts and global composition changes; group membership sets the scale.",
      ("total support","high-count concentration","matched-control agreement")),
    Method("nzmean","Nonzero-mean scaling","available","raw insertion-site counts",
      "Scales each library's mean reads per occupied site to the group average.",
      "factor_j = mean_group(positive_site_mean) / positive_site_mean_j",
      "A comparison method when the intended quantity is support per occupied site rather than density-adjusted support.",
      "Does not adjust insertion density; ignores zeros when estimating the scale and remains sensitive to high counts.",
      ("positive-site mean","density differences","high-count concentration")),
    Method("cpm","Counts per million (total reads, gene level; z-score input)","existing_workflow","gene-level feature read counts",
      "The treated-vs-control z-score scales gene-level counts to counts per million (ytab.fitness.zscore); same scaling as total reads.",
      "CPM_gj = feature_reads_gj / sum_genes(feature_reads_j) * 1000000",
      "Reproducing the existing fitness method.",
      "Not interchangeable with site-level counts; do not feed pre-normalized counts into a workflow that normalizes again.",
      ("feature-level library size","control-control behavior","composition"),"src/ytab/fitness/zscore.py"),
    Method("midlc","MidLC common saturation","available","raw insertion-site counts",
      "Thins every library by binomial sampling to the same depth relative to its own complexity: the group's smallest depth-to-MidLC ratio. MidLC is the depth at which a library is expected to show half of its occupied sites. Downsampling only; counts stay whole numbers.",
      "MidLC_j: sum_s[1-(1-p)^c_sj] = occupied_j/2, MidLC_j = p*total_j; p_j = min_group(total/MidLC) * MidLC_j / total_j; counts ~ Binomial(c, p_j)",
      "Comparing libraries sequenced to different saturation, and essentiality inputs that should not depend on sequencing depth.",
      "Discards reads (random, seeded); the most saturated libraries lose the most. Does not correct insertion-site preference or skew. Units are thinned reads, not a common scale.",
      ("MidLC","saturation ratio","thinning fraction","sites retained"),"src/ytab/normalization/midlc.py"),
    Method("quantile","Hermes zero-preserving quantile (experimental)","available","raw insertion-site counts",
      "Maps positive counts toward a within-library mean empirical quantile function, averaging tied rank targets and preserving every zero. Ties take the average of their rank targets rather than a dense rank.",
      "Q_group(p) = mean_j Q_positive,j(p); map positive ranks to Q_group; average targets within tied counts; retain zeros",
      "Experimental sensitivity comparison when positive-count distribution matching is scientifically justified.",
      "Nonlinear; may erase real abundance/composition differences. Unequal occupied-site counts and ties prevent exact full-distribution matching. Not a physical abundance estimate or a default method.",
      ("positive-count distributions","ties and zeros","control behavior","gene-effect preservation")),
    Method("betageom","Beta-geometric correction (experimental)","available","raw insertion-site counts",
      "Fits beta-geometric parameters including catalogue zeros, maps counts to a group geometric target and optionally applies explicit post-TTR scaling.",
      "p~Beta(kappa*rho,kappa*(1-rho)); X|p~Geom(p); corrected=Geom_target^-1(F_BG(X)); zeros retained",
      "Sensitivity analysis for heavy-tailed insertion support dominated by a few sites when a geometric target is scientifically defensible.",
      "May suppress genuine enrichment; catalogue zeros affect fitted parameters; not a viable-mutant abundance model. Analytic CDF default; seeded simulation available. No automatic preference.",
      ("fit parameters","catalogue sensitivity","upper-tail concentration","effect preservation","CDF safeguards")),
    Method("aBGC","Adaptive beta-geometric correction (experimental)","available","raw insertion-site counts",
      "Searches positive-count trimming fractions for the smallest maximum geometric-quantile discrepancy; preserves zeros; optional explicit final TTR.",
      "choose f minimizing max_r((corrected_r-geometric_quantile_r)^2/geometric_quantile_r); fit positive counts only",
      "Sensitivity option for skewed counts when a geometric target and adaptive distribution-based fitting are appropriate.",
      "The score is not a p-value or biological validation. Tail convention and grid affect selection; can compress real effects; invalid candidates are reported, not silently replaced.",
      ("candidate table","tail sensitivity","selected parameters","effect preservation","seed reproducibility")),
    Method("zinfnb","Zero-inflated negative-binomial scaling (experimental)","available","raw insertion-site counts",
      "Fits a proper ZINB mixture per sample over the full empirical catalogue; the zero term includes the (1-pi) factor.",
      "P(0)=pi+(1-pi)*NB(0); P(x>0)=(1-pi)*NB(x); factor=1/fitted_NB_mean",
      "Sensitivity analysis when a ZINB distribution and fitted NB-component mean are defensible scale references.",
      "Catalogue-dependent zeros are not biological structural zeros. Fit convergence does not validate the model; no automatic preference or subsequent CPM.",
      ("zero interpretation","convergence","parameter boundaries","CDF discrepancy","catalogue sensitivity")),
    Method("emphist","Empirical-histogram normalization","planned","insertion-site counts",
      "Empirical count-distribution adjustment.","Awaiting adaptation and validation.",
      "An exploratory option only after its reference and distribution assumptions are specified.",
      "Not selectable until its exact implementation and Hermes behavior are validated.",("reference distribution","zero preservation","control behavior")),
    Method("QMgeom","Geometric quantile matching","planned","insertion-site counts",
      "Matches counts to a geometric target distribution and applies scaling.","Awaiting adaptation and validation.",
      "Potential diagnostic when a geometric target is scientifically defensible.",
      "Can change distribution shape; not a correction for all Hermes insertion preferences.",("distribution fit","zeros","control behavior")),
)}


def method_catalogue():
    return [asdict(method) for method in METHODS.values()]


def normalize(matrix, method, groups, trim=.05, target=100, bgc_options=None, abgc_options=None):
    if method not in METHODS:
        raise ValueError(f"Unknown method: {method}")
    if METHODS[method].status != "available":
        raise ValueError(f"{method} is {METHODS[method].status}; it is not selectable in the shared site-level runner")
    matrix=sparse.csc_matrix(matrix,copy=True)
    matrix.eliminate_zeros()
    if (not matrix.shape[0] or matrix.shape[1]!=len(groups) or any(not str(g).strip() for g in groups) or
        np.any(~np.isfinite(matrix.data)) or np.any(matrix.data<0) or np.any(matrix.data!=np.floor(matrix.data))):
        raise ValueError("Supply a nonempty site-by-sample raw integer-count matrix and one explicit group per sample")
    total=np.asarray(matrix.sum(axis=0)).ravel()
    occupied=np.diff(matrix.indptr)
    if np.any(total<=0):
        raise ValueError("Every selected sample must have positive raw support; empty samples need review, not a fallback factor")
    positive_mean=total/occupied
    factors=np.ones(matrix.shape[1])
    bgc_output,bgc_details = correct_betageom(matrix,groups,trim=trim,target=target,**(bgc_options or {})) if method=="betageom" else (None,None)
    abgc_output,abgc_details = correct_adaptive_bgc(matrix,groups,trim=trim,target=target,**(abgc_options or {})) if method=="aBGC" else (None,None)
    midlc_output,midlc_details = thin_midlc(matrix,groups) if method=="midlc" else (None,None)
    rows=[]
    for i,group in enumerate(groups):
        members=np.array([g==group for g in groups])
        if method=="ttr":
            positive=matrix.data[matrix.indptr[i]:matrix.indptr[i+1]]
            params=ttr_parameters(positive,matrix.shape[0],trim,target)
            factors[i]=params["factor"]
        else:
            params={}
            if method=="totreads":factors[i]=total[members].mean()/total[i]
            if method=="nzmean":factors[i]=positive_mean[members].mean()/positive_mean[i]
            if method=="zinfnb":
                params=fit_zinfnb(matrix.data[matrix.indptr[i]:matrix.indptr[i+1]],matrix.shape[0])
                factors[i]=params["factor"]
        rows.append(dict(sample_index=i,method=method,normalization_group=group,group_samples=int(members.sum()),
          factor="" if method in {"quantile","betageom","aBGC","midlc"} else float(factors[i]),raw_total_support=int(total[i]),occupied_sites=int(occupied[i]),
          density=float(occupied[i]/matrix.shape[0]),positive_site_mean=float(positive_mean[i]),
          trimmed_nonzero_mean=params.get("trimmed_nonzero_mean",""),
          trim_fraction_each_tail=trim if method=="ttr" else "",target=target if method=="ttr" else ""))
        if bgc_details:rows[-1].update(bgc_details[i])
        if abgc_details:rows[-1].update(abgc_details[i])
        if midlc_details:rows[-1].update(midlc_details[i])
        if method=="zinfnb":rows[-1].update(params)
    output=(midlc_output if method=="midlc" else abgc_output if method=="aBGC" else bgc_output if method=="betageom" else normalize_positive_quantiles(matrix,groups) if method=="quantile" else
        matrix.copy() if method=="nonorm" else matrix.astype(float).multiply(factors).tocsc())
    return output.tocsr(),rows


def run_choices(catalogue, output, methods, group_by, observations=None, trim=.05, target=100, bgc_options=None, abgc_options=None):
    """Apply explicit choices to an existing raw matrix; never normalize twice."""
    catalogue,output=Path(catalogue).resolve(),Path(output).resolve()
    if output==catalogue or output in catalogue.parents:
        raise ValueError("Output must not replace the input catalogue")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is nonempty; choose a new version")
    if not methods or len(set(methods))!=len(methods):
        raise ValueError("Choose one or more distinct methods explicitly")
    with (catalogue/"observations.tsv").open() as handle:
        metadata=list(csv.DictReader(handle,delimiter="\t"))
    ids=[row["observation_id"] for row in metadata]
    if len(set(ids))!=len(ids):raise ValueError("Duplicate observation IDs in input catalogue")
    selected=ids if observations is None else observations
    if not selected or len(set(selected))!=len(selected) or not set(selected)<=set(ids):
        raise ValueError("Choose unique observation IDs present in the input catalogue")
    if group_by not in {"library_id","project","dataset","cohort"}:
        raise ValueError("Choose an explicit grouping: library_id, project, dataset or cohort")
    indices=[ids.index(value) for value in selected]
    raw=sparse.load_npz(catalogue/"raw_counts.npz")
    if raw.shape[1]!=len(ids):raise ValueError("Raw matrix and observation metadata disagree")
    raw=raw[:,indices]
    metadata=[metadata[i] for i in indices]
    groups=[row[group_by] if group_by!="cohort" else "explicit_selected_cohort" for row in metadata]
    source_manifest=json.loads((catalogue/"run_manifest.json").read_text())
    if source_manifest.get("status")!="success" or source_manifest.get("universe_sites")!=raw.shape[0]:
        raise ValueError("Input catalogue is incomplete or its universe size disagrees with the matrix")
    # Validate all choices before creating any output. Each method always starts
    # from the same raw matrix, never another method's normalized output.
    for method in methods:
        normalize(raw,method,groups,trim,target,bgc_options,abgc_options)
    output.mkdir(parents=True,exist_ok=True)
    (output/"implementation_snapshot.py").write_bytes(Path(__file__).read_bytes())
    all_rows=[]
    for method in methods:
        normalized,rows=normalize(raw,method,groups,trim,target,bgc_options,abgc_options)
        folder=output/method
        folder.mkdir()
        sparse.save_npz(folder/"normalized_counts.npz",normalized)
        candidate_rows=[]
        for row,meta in zip(rows,metadata):
            candidate_rows.extend(dict(observation_id=meta["observation_id"],**candidate) for candidate in row.pop("_abgc_candidates",[]))
            row["observation_id"]=meta["observation_id"]
            row["display_name"]=meta["display_name"]
        write_table(folder/"factors.tsv",list(rows[0]),rows)
        if candidate_rows:write_table(folder/"candidates.tsv",list(candidate_rows[0]),candidate_rows)
        (folder/"method.json").write_text(json.dumps(asdict(METHODS[method]),indent=2)+"\n")
        all_rows.extend(rows)
    write_table(output/"method_comparison.tsv",list(dict.fromkeys(k for row in all_rows for k in row)),all_rows)
    write_table(output/"observations.tsv",list(metadata[0]),metadata)
    manifest=dict(status="success",methods=methods,group_by=group_by,observations=selected,
      universe_sites=raw.shape[0],catalogue=str(catalogue),reference_sha256=source_manifest["reference_sha256"],
      source_hashes={name:sha256(catalogue/name) for name in ("raw_counts.npz","observations.tsv","sites.tsv.gz","run_manifest.json")},
      implementation_sha256=sha256(__file__),ttr_implementation_sha256=sha256(Path(__file__).with_name("ttr.py")),
      quantile_implementation_sha256=sha256(Path(__file__).with_name("quantile.py")) if "quantile" in methods else None,
      zinfnb_implementation_sha256=sha256(Path(__file__).with_name("zinfnb.py")) if "zinfnb" in methods else None,
      zinfnb_parameters=dict(likelihood="proper ZINB; zero term includes (1-pi)",fit_universe="full empirical catalogue",optimizer="L-BFGS-B",maxiter=1000,ftol=1e-12,gtol=1e-7,log_mean_size_bounds=LOG_BOUNDS,start_sizes=START_SIZES,pi="profiled analytically; [0,1-1e-12]",target="NB component mean 1; no group target") if "zinfnb" in methods else None,
      betageom_implementation_sha256=sha256(Path(__file__).with_name("betageom.py")) if "betageom" in methods else None,
      bgc_parameters=dict(cdf="analytic",seed=0,draws=200000,post_ttr=True) | (bgc_options or {}) | dict(trim=trim,target=target,fit_universe="full catalogue including zeros",fit_trim_each_tail=.001,rho_lower_bound=.0001,kappa_lower_bound=10) if "betageom" in methods else None,
      adaptive_bgc_implementation_sha256=sha256(Path(__file__).with_name("adaptive_bgc.py")) if "aBGC" in methods else None,
      abgc_parameters=dict(cdf="analytic",seed=0,draws=200000,post_ttr=True,tail="endpoint",fractions=list(DEFAULT_FRACTIONS)) | (abgc_options or {}) | dict(trim=trim,target=target,fit_universe="positive counts only",score="maximum normalized squared quantile discrepancy; not a p-value") if "aBGC" in methods else None,
      midlc_implementation_sha256=sha256(Path(__file__).with_name("midlc.py")) if "midlc" in methods else None,
      midlc_parameters=dict(target="common saturation: smallest total/MidLC ratio in the normalization group",thinning="binomial, downsample only",seed=MIDLC_SEED,midlc="depth with half the occupied sites expected under binomial thinning (exact)") if "midlc" in methods else None,
      software_versions={name:version(name) for name in ("numpy","scipy")},
      ttr_parameters=dict(trim_each_tail=trim,target=target) if "ttr" in methods else None,
      scale_warning="Units differ by method. Compare within-method relative results; larger absolute normalized counts are not evidence of better normalization.",
      workflow_policy="Opt-in evaluation only; the fitness z-score default (CPM) is unchanged; no downstream double normalization")
    (output/"run_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest
