"""Hermes site-restricted resampling, adapted from TRANSIT's resampling test (DeJesus et al. 2015) to empirical Hermes sites.

Per gene, T = mean(treatment) - mean(control) over eligible site x sample normalized counts.
Under the null, sample labels are exchangeable within each site, independently across sites
(site-restricted resampling). Then T* is an affine function of U = sum_i S_i, where S_i is
the sum of K values drawn without replacement from site i's J+K values, uniformly over the
C(J+K, K) label subsets. The p-value is either exact (convolution of the per-site subset-sum
distributions on a grid, with a conservative rounding bound) or Monte Carlo (seeded draws).
Input is a YTAB scalar normalization (ttr, totreads or nzmean). Options (v3):
- alternative: two-sided (default) or one-sided depleted / enriched;
- recenter: one scalar rescaling of the treatment, estimated from genes with pooled mean >= the median,
  against the shift a scalar normalization picks up when many genes change. "median-gene": the median
  gene log ratio becomes 0; "trimmed-sum" (TMM-like): genes in either recenter_trim log-ratio tail are
  dropped and the remaining genes' summed means are matched. The factor is reported;
- winsorize: label-invariant winsorization (adapted from TRANSIT's -winz option): in each gene, the largest site x sample
  value (over both conditions) is replaced by the second largest, so the permutation stays exact.
"""
from dataclasses import asdict, dataclass
from itertools import combinations
import hashlib
import json
from pathlib import Path

import numpy as np

from ytab.essentiality.control_conditioned_gumbel import annotations, benjamini_hochberg, gene_sites, load_design, _json_ready
from ytab.normalization.ttr import sha256, write_table

MODEL = "Hermes site-restricted resampling"
REVISION = "hermes_resampling_v3"


@dataclass(frozen=True)
class Settings:
    inference: str = "exact"            # exact | monte-carlo
    opportunities: str = "any-support"  # any-support (label-invariant) | control-support
    samples: int = 10000                # Monte Carlo permutations
    grid_bins: int = 65536              # exact: grid resolution for the convolution
    min_sites: int = 1
    pseudocount: float = 1.
    fdr: float = .05
    seed: int = 2026
    nterm: float = 0.
    cterm: float = 0.
    ignore_stop: bool = True
    alternative: str = "two-sided"      # two-sided | depleted | enriched
    recenter: str = "none"              # none | median-gene | trimmed-sum
    recenter_trim: float = .2           # trimmed-sum: fraction of genes trimmed from each log-ratio tail
    winsorize: bool = False
    asymmetry_warning: float = 2.        # warn when replicate noise differs more than this between conditions

    def validate(self):
        if (self.alternative not in {"two-sided", "depleted", "enriched"} or self.recenter not in {"none", "median-gene", "trimmed-sum"}
                or not 0 <= self.recenter_trim < .5 or not isinstance(self.winsorize, bool) or not self.asymmetry_warning >= 1):
            raise ValueError("Invalid alternative, recenter or winsorize")
        if (self.inference not in {"exact", "monte-carlo"} or self.opportunities not in {"any-support", "control-support"} or
                not isinstance(self.samples, int) or self.samples < 100 or not isinstance(self.grid_bins, int) or self.grid_bins < 256 or
                not isinstance(self.min_sites, int) or self.min_sites < 1 or not self.pseudocount >= 0 or not 0 < self.fdr < 1 or
                not isinstance(self.seed, int) or self.seed < 0 or not 0 <= self.nterm < 100 or not 0 <= self.cterm < 100 or
                self.nterm+self.cterm >= 100):
            raise ValueError("Invalid inference, opportunities, samples (>=100), grid_bins (>=256), min_sites, pseudocount, FDR, seed or terminal percentages")


def subset_sums(control, treatment):
    """Per-site sums of every K-subset of the J+K values (treatment label sets) and the observed treatment sum.

    Columns are ordered controls then treatments, so the observed labelling is the last subset.
    """
    values = np.concatenate([control, treatment], axis=1); K = treatment.shape[1]
    subsets = np.array(list(combinations(range(values.shape[1]), K)))
    return values[:, subsets].sum(2), treatment.sum(1)


def _tail(t, t_obs, slack, alternative):
    """Mask of permutation statistics at least as extreme as t_obs (slack widens it conservatively)."""
    if alternative == "depleted":return t <= t_obs+slack
    if alternative == "enriched":return t >= t_obs-slack
    return np.abs(t) >= abs(t_obs)-slack


def exact_units(blocks, n_control, n_treatment, bins, alternative="two-sided"):
    """Exact p over independent permutation units (two-sided or one-sided).

    blocks: list of (sums, totals) per library, sums units x m (equally likely treatment-label subset
    sums, observed labelling last) and totals the units' all-sample sums. T = mean(treatment) -
    mean(control) = a*U - b*sum(totals) with U = sum of chosen subset sums, a = 1/NB + 1/NA, b = 1/NA.
    Subset sums are rounded to a grid of step delta = range/bins (total error <= units*delta/2).
    Returns the conservative p (within that margin of |T_obs|, ties included), the exact p of the
    rounded problem (ties included), and T_obs.
    """
    a, b = 1/n_treatment+1/n_control, 1/n_control
    total = sum(float(t.sum()) for _, t in blocks); observed = sum(float(s[:, -1].sum()) for s, _ in blocks)
    t_obs = a*observed-b*total
    low = [s.min(1) for s, _ in blocks]; spread = sum(float((s.max(1)-l).sum()) for (s, _), l in zip(blocks, low))
    if spread <= 0:return 1., 1., t_obs
    delta = spread/bins; units = sum(len(s) for s, _ in blocks)
    indices = [np.rint((s-l[:, None])/delta).astype(np.int64) for (s, _), l in zip(blocks, low)]
    distribution = np.zeros(int(sum(i.max(1).sum() for i in indices))+1); distribution[0] = 1.; top = 0
    for index in indices:
        for row in index:
            values, counts = np.unique(row, return_counts=True); updated = np.zeros_like(distribution)
            for v, c in zip(values, counts):updated[v:v+top+1] += distribution[:top+1]*(c/len(row))
            distribution = updated; top += int(row.max())
    base = sum(float(l.sum()) for l in low)
    t = a*(base+delta*np.arange(len(distribution)))-b*total
    tolerance = 1e-9*max(1., abs(t_obs), a*spread)
    upper = float(distribution[_tail(t, t_obs, a*units*delta/2+tolerance, alternative)].sum())
    t_rounded = a*(base+delta*sum(int(i[:, -1].sum()) for i in indices))-b*total
    rounded = float(distribution[_tail(t, t_rounded, tolerance, alternative)].sum())
    return min(1., upper), min(1., rounded), t_obs


def monte_carlo_units(rng, blocks, n_control, n_treatment, samples, alternative="two-sided"):
    """Seeded permutations of independent units; p = (b+1)/(S+1)."""
    a, b = 1/n_treatment+1/n_control, 1/n_control
    total = sum(float(t.sum()) for _, t in blocks); t_obs = a*sum(float(s[:, -1].sum()) for s, _ in blocks)-b*total
    units = sum(len(s) for s, _ in blocks); exceed = 0; chunk = max(1, int(4e6//max(units, 1)))
    tolerance = 1e-12*max(1., abs(t_obs))
    for begin in range(0, samples, chunk):
        rows = min(chunk, samples-begin); draws = np.zeros(rows)
        for s, _ in blocks:draws += s[np.arange(len(s)), rng.integers(0, s.shape[1], size=(rows, len(s)))].sum(1)
        exceed += int(_tail(a*draws-b*total, t_obs, tolerance, alternative).sum())
    return (exceed+1)/(samples+1), t_obs


def exact_p(sums, observed, J, K, bins):
    """Single-library form: sites x C(J+K, K) subset sums with the observed labelling last."""
    n = len(sums); totals = sums.mean(1)*(J+K)/K
    return exact_units([(sums, totals)], n*J, n*K, bins)


def monte_carlo_p(rng, sums, observed, J, K, samples):
    n = len(sums)
    return monte_carlo_units(rng, [(sums, sums.mean(1)*(J+K)/K)], n*J, n*K, samples)


def replicate_spread(values):
    """Robust SD of log2 ratios between replicate columns (mean over pairs) at mid-abundance sites; None if < 2 replicates."""
    values = np.asarray(values, dtype=float)
    if values.shape[1] < 2:return None
    spreads = []
    for i in range(values.shape[1]):
        for j in range(i+1, values.shape[1]):
            pair = values[:, [i, j]]; both = (pair > 0).all(1)
            if both.sum() < 50:continue
            mean = pair[both].mean(1); lo, hi = np.quantile(mean, [.25, .9]); band = (mean >= lo) & (mean <= hi)
            r = np.log2(pair[both][band, 0]/pair[both][band, 1]); spreads.append(1.4826*np.median(np.abs(r-np.median(r))))
    return float(np.mean(spreads)) if spreads else None


def noise_asymmetry(control, treatment, threshold=2.):
    """Replicate-noise asymmetry between conditions: max(spread_T/spread_C, spread_C/spread_T) (>= 1), with both spreads.

    Pre-registered v3 evaluation (2 vs 2, fresh seeds): asymmetry ~1-2 calibrated (null p<=0.05 rate
    4.1-5.5%); ~3.6-3.8 (noisy control or treatment) 6.1-6.2% with BGC + recentering, 8.4-8.5% with TTR.
    Gale FLU: 3.5-3.6. Above `threshold` a warning is returned: p-values are mildly anti-conservative.
    """
    c, t = replicate_spread(control), replicate_spread(treatment)
    if not c or not t:
        return dict(control_spread=c, treatment_spread=t, asymmetry=None, warning="fewer than 2 replicates in a condition; asymmetry not measurable")
    a = float(max(t/c, c/t)); noisier = "treatment" if t > c else "control"
    warning = (f"{noisier} replicates are {a:.1f}x noisier than the other condition's (threshold {threshold:g}); permutation p-values are "
               "mildly anti-conservative in this regime (about 6% instead of 5% at p<=0.05 with BGC + recentering, 8-9% with TTR); "
               "treat calls near the q cutoff with caution and prefer BGC normalization") if a > threshold else None
    return dict(control_spread=c, treatment_spread=t, asymmetry=a, threshold=threshold, warning=warning)


def _winsorized(units):
    """Label-invariant winsorization of one gene: its largest site x sample value (both conditions) becomes the second largest."""
    values = np.concatenate([np.concatenate([st["control"][i], st["treatment"][i]], axis=1).ravel() for st, i in units])
    if len(values) < 2:return units
    top, second = np.sort(values)[-1], np.sort(values)[-2]
    if top == second:return units
    out = []
    for st, i in units:
        if not len(i):out.append((st, i)); continue
        c, t = np.minimum(st["control"][i], second), np.minimum(st["treatment"][i], second)
        sums, _ = subset_sums(c, t)   # gene-local unit arrays, indexed 0..len(i)-1
        out.append((dict(library=st["library"], control=c, treatment=t, sums=sums, totals=c.sum(1)+t.sum(1)), np.arange(len(i))))
    return out


def _gene_rng(seed, gene_id):
    key = int.from_bytes(hashlib.sha256(str(gene_id).encode()).digest()[:8], "big")
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(key, 7)))


def analyze(control, treatment, genes, gene_indices, positions, settings=None, control_libraries=None, treatment_libraries=None):
    """control: sites x J and treatment: sites x K normalized counts (same scalar normalization per library).

    Optional per-column library labels make permutation units (site, library): labels are permuted only
    among the samples of one library at one site. One library (default) gives plain site-restricted resampling.
    """
    settings = settings or Settings(); settings.validate()
    control = np.asarray(control, dtype=float); treatment = np.asarray(treatment, dtype=float)
    if control.ndim != 2 or treatment.ndim != 2 or len(control) != len(treatment) or control.shape[1] < 1 or treatment.shape[1] < 1:
        raise ValueError("Need sites x >=1 control and sites x >=1 treatment normalized counts")
    if np.any(~np.isfinite(control)) or np.any(~np.isfinite(treatment)) or np.any(control < 0) or np.any(treatment < 0):
        raise ValueError("Normalized counts must be finite and nonnegative")
    if len(genes) != len(gene_indices) or len({g["id"] for g in genes}) != len(genes):
        raise ValueError("Gene identities and site lists must align and be unique")
    control_libraries = ["library"]*control.shape[1] if control_libraries is None else list(control_libraries)
    treatment_libraries = ["library"]*treatment.shape[1] if treatment_libraries is None else list(treatment_libraries)
    libraries = list(dict.fromkeys(control_libraries+treatment_libraries))
    if (len(control_libraries) != control.shape[1] or len(treatment_libraries) != treatment.shape[1] or
            set(control_libraries) != set(treatment_libraries)):
        raise ValueError("Every library needs >=1 control and >=1 treatment column")
    rescale = 1.
    if settings.recenter != "none":
        a_all, b_all = control.mean(1), treatment.mean(1); ratios, pooled = [], []
        for idx in gene_indices:
            idx = np.asarray(idx, dtype=np.int64)
            if len(idx):ratios.append((a_all[idx].mean(), b_all[idx].mean())); pooled.append(a_all[idx].mean()+b_all[idx].mean())
        ratios, pooled = np.array(ratios), np.array(pooled); keep = (pooled >= np.median(pooled)) & (ratios[:, 0] > 0) & (ratios[:, 1] > 0)
        if keep.sum() >= 10:
            log_ratio = np.log2(ratios[keep, 1]/ratios[keep, 0])
            if settings.recenter == "median-gene":rescale = float(2**-np.median(log_ratio))
            else:   # TMM-like: drop genes in either log-ratio tail, then match the remaining genes' total means (mean-consistent)
                lo, hi = np.quantile(log_ratio, [settings.recenter_trim, 1-settings.recenter_trim]); core = (log_ratio >= lo) & (log_ratio <= hi)
                rescale = float(ratios[keep][core, 0].sum()/ratios[keep][core, 1].sum())
        treatment = treatment*rescale
    strata = []
    for library in libraries:
        a = control[:, [j for j, x in enumerate(control_libraries) if x == library]]
        b = treatment[:, [k for k, x in enumerate(treatment_libraries) if x == library]]
        sums, _ = subset_sums(a, b)
        # any-support depends only on the pooled values, so it is invariant to label permutation; control-support is not.
        eligible = (a.sum(1)+b.sum(1) > 0) if settings.opportunities == "any-support" else a.sum(1) > 0
        strata.append(dict(library=library, control=a, treatment=b, sums=sums, totals=a.sum(1)+b.sum(1), eligible=eligible))
    rows = []
    for gene, idx in zip(genes, gene_indices):
        idx = np.asarray(idx, dtype=np.int64)
        units = [(st, idx[st["eligible"][idx]]) for st in strata]; n = sum(len(i) for _, i in units)
        row = dict(ORF=gene["id"], Name=gene["name"], Description=gene["description"], chromosome=gene["chrom"], strand=gene["strand"], n_sites=n)
        if n < settings.min_sites:
            rows.append(row | dict(status="too few eligible sites" if n else "no eligible sites", in_family=False)); continue
        if settings.winsorize:
            units = _winsorized(units)
        a = np.concatenate([st["control"][i].ravel() for st, i in units]); b = np.concatenate([st["treatment"][i].ravel() for st, i in units])
        mean_a, mean_b = float(a.mean()), float(b.mean())
        log2 = lambda x, y:float(np.log2((y+settings.pseudocount)/(x+settings.pseudocount))) if settings.pseudocount > 0 or (x > 0 and y > 0) else None
        row |= dict(mean_control=mean_a, mean_treatment=mean_b, sum_control=float(a.sum()), sum_treatment=float(b.sum()), log2FC=log2(mean_a, mean_b))
        if len(strata) > 1:
            row["library_log2FC"] = ";".join(f"{st['library']}={log2(st['control'][i].mean(), st['treatment'][i].mean()):.3f}" if len(i) else f"{st['library']}=NA"
                                             for st, i in units)
            row["libraries_with_sites"] = sum(len(i) > 0 for _, i in units)
        blocks = [(st["sums"][i], st["totals"][i]) for st, i in units if len(i)]
        if settings.inference == "exact":
            p, p_rounded, t_obs = exact_units(blocks, len(a), len(b), settings.grid_bins, settings.alternative)
            row |= dict(delta_mean=t_obs, p_value=p, p_value_rounded_problem=p_rounded)
        else:
            p, t_obs = monte_carlo_units(_gene_rng(settings.seed, gene["id"]), blocks, len(a), len(b), settings.samples, settings.alternative)
            row |= dict(delta_mean=t_obs, p_value=p)
        rows.append(row | dict(status="tested", in_family=True))
    family = [r for r in rows if r["in_family"]]
    if family:
        for row, q in zip(family, benjamini_hochberg([r["p_value"] for r in family])):row["q_bh"] = float(q)
    for row in rows:
        row["call"] = ("not tested" if not row["in_family"] else "not significant" if row["q_bh"] > settings.fdr else
                       "depleted" if row["delta_mean"] < 0 and settings.alternative != "enriched" else
                       "enriched" if row["delta_mean"] > 0 and settings.alternative != "depleted" else "not significant")
    diagnostics = dict(inference=settings.inference, opportunities=settings.opportunities, libraries=libraries, alternative=settings.alternative,
                       recenter=settings.recenter, treatment_rescale=rescale, winsorize=settings.winsorize,
                       replicate_noise=noise_asymmetry(control, treatment, settings.asymmetry_warning) if len(strata) == 1 else None,
                       eligible_units=int(sum(st["eligible"].sum() for st in strata)),
                       family_size=len(family), calls={c:sum(r["call"] == c for r in rows) for c in ("depleted", "enriched", "not significant", "not tested")})
    if settings.inference == "exact" and family:
        diagnostics["max_discretization_gap"] = float(max(r["p_value"]-r["p_value_rounded_problem"] for r in family))
    return rows, diagnostics


def run(normalized, gff, reference, output, controls, treatments, settings=None, normalization=None, stratified=False):
    """File-based runner on a YTAB scalar normalization output (ttr default, totreads or nzmean).

    stratified=True: several independent libraries (e.g. paired pools), permuted within library.
    """
    settings = settings or Settings(); settings.validate()
    design = load_design(normalized, normalization, reference, output, controls, treatments, min_controls=1, stratified=stratified, nonlinear=True)
    output, sites, scaled = design["output"], design["sites"], design["scaled"]; J = len(design["controls"])
    genes = annotations(Path(gff).resolve())
    positions = np.array([p for _, p in sites], dtype=np.int64)
    libraries = design["libraries"] if stratified else None
    indices = gene_sites(genes, sites, settings)
    rows, diagnostics = analyze(scaled[:, :J], scaled[:, J:], genes, indices, positions, settings,
                                libraries[:J] if libraries else None, libraries[J:] if libraries else None)
    output.mkdir(parents=True, exist_ok=True)
    leading = ["ORF", "Name", "Description", "n_sites", "mean_control", "mean_treatment", "log2FC", "delta_mean", "p_value", "q_bh", "call", "status"]
    write_table(output/"genes.tsv", leading+list(dict.fromkeys(k for r in rows for k in r if k not in leading)), rows)
    # Normalized counts at each gene's sites, per library (for site-level plots); before recentering.
    K = scaled.shape[1]-J
    names = [f"control_{i+1}" for i in range(J)]+[f"treatment_{i+1}" for i in range(K)]
    site_rows = [dict(site_index=int(i), chromosome=sites[i][0], position=sites[i][1], genes=g["id"], **{n:float(scaled[i, j]) for j, n in enumerate(names)})
                 for g, idx in zip(genes, indices) for i in sorted(idx)]
    write_table(output/"sites.tsv.gz", ["site_index", "chromosome", "position", "genes"]+names, site_rows)
    norm = design["norm"]
    manifest = dict(status="success", analysis=MODEL, algorithm_revision=REVISION, parameters=asdict(settings),
        controls=design["chosen"][:J], treatments=design["chosen"][J:],
        normalization=dict(method=norm["method"], input=str(norm["folder"]), catalogue=str(design["catalogue"]),
            factors=None if design["factors"] is None else {x:float(f) for x, f in zip(design["selected"], design["factors"])},
            group=norm["groups"][design["selected"][0]],
            hashes={name:sha256(norm["folder"]/name) for name in ("normalized_counts.npz", "factors.tsv", "method.json")}
                   | dict(run_manifest=sha256(norm["root"]/"run_manifest.json")),
            check=("normalized_counts preserve raw zeros and within-sample count order (verified; nonlinear method)" if design["factors"] is None
                   else "normalized_counts == raw catalogue counts * factor (verified)")),
        model=dict(statistic="T = mean(treatment) - mean(control) over eligible site x sample normalized counts",
            null=("site-restricted within library: labels exchangeable among one library's samples at each site, independently across (site, library) units"
                  if stratified else "site-restricted: labels exchangeable among samples within each site, independently across sites"),
            p_value=("exact: convolution of per-site uniform K-subset-sum distributions on a grid; conservative bound after rounding"
                     if settings.inference == "exact" else "Monte Carlo: seeded per-gene site-restricted permutations, p=(b+1)/(S+1)"),
            opportunities=("sites with >0 normalized support in any compared sample (label-invariant)" if settings.opportunities == "any-support"
                           else "sites with >0 combined control support (not label-invariant)"),
            multiple_testing="BH over genes with >= min_sites eligible sites; call direction from sign(T)"),
        diagnostics=diagnostics, implementation_sha256=sha256(__file__),
        limitations=["Unequal replicate noise between conditions (see diagnostics.replicate_noise) makes permutation p-values mildly anti-conservative",
            "Sites are treated as independent within a gene; sample-level effects shared across sites are not modelled",
            "Exchangeability across conditions is assumed; unequal replicate noise between conditions can violate it",
            "Scalar normalization assumes most sites are unaffected by treatment",
            "Overlapping genes share sites; BH assumes valid p-values"])
    (output/"run_manifest.json").write_text(json.dumps(_json_ready(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
