"""Treated-vs-control local-abundance z-score: call genes depleted or enriched in a treated library against its matched control.

Port of the established YTAB fitness analysis (formerly scripts/ytab_treated_vs_parent_screen.R); the calculation is
unchanged:

1. Gene values per library: by default CPM from each sample's raw SummaryTable feature table (reads / total * 1e6).
   Any site normalization can be used instead (`normalized=`): the gene value is then the sum of the normalized site
   counts in the gene, from a Site Normalization run, expressed per million of the mean library total (one constant
   for all libraries, so the pseudocount has the same meaning as with CPM).
2. Control noise model per background, from two control (parent) libraries: A = mean of the two, M = log2 ratio with a
   pseudocount (0.5). Genes sorted by A (descending); running mean of A and running SD of M over 41 genes; the SD is
   fitted as SD(A) = m1 + m2 * A^m3 (bounded least squares; start from a log-log line).
3. Each treated/control pair: A = mean, log2FC = log2((treated + pc) / (control + pc)), z = log2FC / SD(A). The call
   threshold is the 0.99 quantile of |z| between the two controls.
4. Calls. combine="replicates" (default, as in Gale et al. 2020, G3 10:3859): the replicate libraries of each background x
   treated condition are averaged (control mean, treated mean) and each gene gets one log2FC and one z against the same
   threshold; final_call is depleted, enriched or none. combine="pairs" (the earlier YTAB rule): per gene across pairs,
   consistently depleted/enriched (2+ pairs, none opposite), single-pool, mixed, or none. The per-pair table is always
   written.
5. Optional recentering (recenter="trimmed-sum", off by default): one scalar per treated library or treated mean, from genes
   with pooled mean >= the median after dropping the recenter_trim log-ratio tails (TMM-like, as in resampling).

Writes the same tables as before into OUTPUT/tables, so the app and the runner read them unchanged.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, asdict
import json
import math
from pathlib import Path
import platform
import re

import numpy as np
import pandas as pd
from scipy import optimize

FEATURE_CANDIDATES = ["standard_name", "feature_name", "feature", "gene", "gene_name", "locus_tag", "orf", "qng_id"]
READS_CANDIDATES = ["reads", "read_count", "feature_reads", "reads_per_feature", "total_reads", "sum_reads"]


CONTRACT_VERSION = 1  # tables and manifest fields listed in ytab.fitness.contract
UNIT_NAMES = dict(betageom="BGC", aBGC="aBGC", ttr="TTR", totreads="Total", nzmean="Nonzero-mean", quantile="Quantile",
                  zinfnb="ZINB", midlc="MidLC", nonorm="Unnormalized")


@dataclass
class Settings:
    pseudocount: float = .5
    z_quantile: float = .99
    window: int = 40           # running window: window/2 genes on each side (41 genes)
    min_abundance: float = .1  # running A below this is not used in the noise fit
    combine: str = "replicates"  # replicates (average replicates, one z per gene) | pairs (per-pair calls, consistency rule)
    recenter: str = "none"       # none | trimmed-sum
    recenter_trim: float = .2
    hit_min_log2FC: float = 1.     # top hits: called genes with |log2FC| at least this (2-fold, Gale et al. 2020)
    hit_min_abundance: float = 0.  # optional abundance floor for top hits (per million); off by default

    def __post_init__(self):
        if self.combine not in {"replicates", "pairs"} or self.recenter not in {"none", "trimmed-sum"} or not 0 <= self.recenter_trim < .5:
            raise ValueError("Invalid combine, recenter or recenter_trim")


def clean_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(x).lower())


def pick_col(columns, candidates, label):
    clean = [clean_name(c) for c in columns]
    for candidate in candidates:
        if clean_name(candidate) in clean:return columns[clean.index(clean_name(candidate))]
    raise ValueError(f"Could not detect {label} column among: {', '.join(map(str, columns))}")


def read_feature_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[""])
    cols = [clean_name(c) for c in df.columns]
    if "reads" not in cols and "standard_name" not in cols:
        alt = pd.read_csv(path, skiprows=1, dtype=str, keep_default_na=False, na_values=[""])
        if alt.shape[1] > df.shape[1]:df = alt
    return df


def sample_gene_reads(sample: str, path: Path) -> pd.DataFrame:
    df = read_feature_table(path)
    feature, reads = pick_col(list(df.columns), FEATURE_CANDIDATES, "feature/gene id"), pick_col(list(df.columns), READS_CANDIDATES, "reads")
    out = df.copy(); out.insert(0, "sample", sample)
    out["feature_id"] = df[feature].astype("string"); out["reads"] = pd.to_numeric(df[reads], errors="coerce")
    return out[out.feature_id.notna() & (out.feature_id != "") & out.reads.notna()]


def first_nonempty(values: pd.Series):
    v = values.dropna().astype(str); v = v[v != ""]
    return v.iloc[0] if len(v) else None


def running(values: np.ndarray, half: int, fun) -> np.ndarray:
    """Centered running statistic over 2*half+1 values; NaN where the window is incomplete."""
    n = len(values); out = np.full(n, np.nan)
    for i in range(half, n-half):out[i] = fun(values[i-half:i+half+1])
    return out


def quantile7(x: np.ndarray, p: float) -> float:
    x = np.sort(np.asarray(x, dtype=float)[~np.isnan(x)])
    if not len(x):return math.nan
    h = (len(x)-1)*p; lo = math.floor(h); hi = min(lo+1, len(x)-1)
    return float(x[lo] + (h-lo)*(x[hi]-x[lo])) if x[hi] != x[lo] else float(x[lo])


def sd_from_fit(coef: dict, a: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return np.maximum(coef["m1"]+coef["m2"]*np.power(a, coef["m3"]), 1e-8)


def control_model(values: pd.DataFrame, c1: str, c2: str, s: Settings) -> dict:
    a = (values[c1].to_numpy()+values[c2].to_numpy())/2
    m = np.log2((values[c1].to_numpy()+s.pseudocount)/(values[c2].to_numpy()+s.pseudocount))
    order = np.argsort(-a, kind="stable")  # descending A; ties keep gene order
    tab = pd.DataFrame(dict(feature_id=values.index.to_numpy()[order], A_ctrl=a[order], M_ctrl=m[order]))
    half = s.window//2
    tab["run_A"] = running(tab.A_ctrl.to_numpy(), half, np.mean)
    tab["run_SD"] = running(tab.M_ctrl.to_numpy(), half, lambda v:np.std(v, ddof=1))
    fit = tab[tab.run_A.notna() & tab.run_SD.notna() & (tab.run_A >= s.min_abundance) & (tab.run_SD > 0)]
    if len(fit) < 20:raise ValueError(f"Too few points to fit the control-control noise model for {c1} vs {c2}")
    x, y = fit.run_A.to_numpy(), fit.run_SD.to_numpy()
    slope = np.polyfit(np.log(x), np.log(y), 1)[0]
    m3 = max(min(slope, .3), -2.5); m1 = quantile7(y, .02)
    m2 = max((float(np.median(y))-m1)/(float(np.median(x))**m3), 1e-6)
    result = optimize.least_squares(lambda p:p[0]+p[1]*x**p[2]-y, [m1, m2, m3], bounds=([0, 0, -5], [np.inf, np.inf, .5]),
                                    method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
    coef = dict(m1=float(result.x[0]), m2=float(result.x[1]), m3=float(result.x[2]))
    return dict(table=tab, coef=coef, start=dict(m1=m1, m2=m2, m3=m3), converged=bool(result.success))


def recenter_factor(control: np.ndarray, treated: np.ndarray, trim: float) -> float:
    """Scalar for the treated values: genes with pooled mean >= the median and both values > 0, the trim tails of their
    log2 ratios dropped, remaining sums matched (as resampling's trimmed-sum recentering)."""
    pooled = control+treated; keep = (pooled >= np.median(pooled)) & (control > 0) & (treated > 0)
    if keep.sum() < 10:return 1.
    ratio = np.log2(treated[keep]/control[keep]); lo, hi = np.quantile(ratio, [trim, 1-trim]); core = (ratio >= lo) & (ratio <= hi)
    return float(control[keep][core].sum()/treated[keep][core].sum())


def call_of(z: np.ndarray, thr: float) -> np.ndarray:
    return np.where(z >= thr, "enriched", np.where(z <= -thr, "depleted", "none"))


def final_call(enriched: int, depleted: int) -> str:
    if enriched >= 2 and depleted == 0:return "consistently_enriched"
    if depleted >= 2 and enriched == 0:return "consistently_depleted"
    if enriched >= 1 and depleted >= 1:return "mixed"
    if enriched == 1:return "single_pool_enriched"
    if depleted == 1:return "single_pool_depleted"
    return "none"


def load_design(path: Path, comparisons: list[str] | None = None):
    design = pd.read_csv(path, dtype=str, keep_default_na=False)
    design = design[design["include"].str.strip().str.lower().isin(["1", "true", "yes", "y"])]
    ids = comparisons or list(design.comparison_id)
    contrasts = design[design.comparison_id.isin(ids)]
    if contrasts.empty:raise ValueError("No selected comparisons were found in the design.")
    samples = pd.concat([design.rename(columns=dict(parent_sample="sample", parent_condition="condition"))[["sample", "condition", "background", "pool"]],
                         design.rename(columns=dict(treated_sample="sample", treated_condition="condition"))[["sample", "condition", "background", "pool"]]])
    return contrasts, samples.drop_duplicates("sample", keep="first").reset_index(drop=True)


def gene_values_from_feature_tables(sample_map: pd.DataFrame, sample_files: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(long raw feature table, gene x sample reads) from raw SummaryTable files."""
    missing = [s for s in sample_map["sample"] if s not in sample_files or not Path(sample_files[s]).is_file()]
    if missing:raise ValueError("Every design sample requires an existing raw SummaryTable input: " + ", ".join(missing))
    long = pd.concat([sample_gene_reads(s, Path(sample_files[s])) for s in sample_map["sample"]], ignore_index=True)
    long = long.merge(sample_map, on="sample", how="left")
    reads = long.groupby(["feature_id", "sample", "condition", "background", "pool"], sort=True, as_index=False)["reads"].sum()
    return long, reads


def _run_gene_sums(run_dir: Path, counts_path: Path, observations_path: Path, gff: Path, samples: list[str]) -> pd.DataFrame:
    """Gene x sample sums of a site-count matrix of a Site Normalization run over each gene's sites (the CDS, with the
    same N/C-terminal trimming as the other gene-level analyses)."""
    import gzip
    from scipy import sparse
    from ytab.essentiality.control_conditioned_gumbel import annotations, gene_sites
    from ytab.essentiality.library_correlation import SITE_RULES, gene_key
    observations = pd.read_csv(observations_path, sep="\t", dtype=str)
    missing = [s for s in samples if s not in set(observations["sample"])]
    if missing:raise ValueError(f"Normalization run {Path(run_dir).name} does not include: {', '.join(missing)}")
    with gzip.open(Path(run_dir)/"experiment/matrix/sites.tsv.gz", "rt") as handle:
        sites = pd.read_csv(handle, sep="\t", dtype={"chromosome":str})
    counts = sparse.load_npz(counts_path).tocsc()
    genes = annotations(Path(gff)); idx = gene_sites(genes, list(zip(sites.chromosome, sites.position_1based.astype(int))), SITE_RULES)
    columns = {s:int(np.flatnonzero(observations["sample"].to_numpy() == s)[0]) for s in samples}
    values = {s:np.asarray(counts[:, j].todense()).ravel() for s, j in columns.items()}
    rows = {gene_key(g["id"]):{s:float(values[s][i].sum()) if len(i) else 0. for s in samples} for g, i in zip(genes, idx)}
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("feature_id")


def gene_values_from_normalization(run_dir: Path, method: str, gff: Path, samples: list[str]) -> pd.DataFrame:
    """Gene x sample values from a Site Normalization run: the method's normalized site counts summed per gene."""
    methods = Path(run_dir)/"experiment"/"methods"; path = methods/method/"normalized_counts.npz"
    if not path.is_file():raise ValueError(f"Normalization run {Path(run_dir).name} has no {method} counts")
    return _run_gene_sums(run_dir, path, methods/"observations.tsv", gff, samples)


def raw_gene_reads_from_run(run_dir: Path, gff: Path, samples: list[str]) -> pd.DataFrame:
    """Gene x sample raw reads (before any normalization) from the run's site matrix: library sizes."""
    matrix = Path(run_dir)/"experiment"/"matrix"
    return _run_gene_sums(run_dir, matrix/"raw_counts.npz", matrix/"observations.tsv", gff, samples)


def gene_annotation_from_gff(gff: Path) -> pd.DataFrame:
    """feature_id, standard_name, common_name, chrom, start, stop, strand from the GFF (any species)."""
    from ytab.essentiality.control_conditioned_gumbel import annotations
    from ytab.essentiality.library_correlation import gene_key
    rows = [dict(feature_id=gene_key(g["id"]), standard_name=g["id"], common_name=g.get("name") or g["id"], chrom=g["chrom"],
                 start=min(b[0] for b in g["blocks"]), stop=max(b[1] for b in g["blocks"]), strand=g["strand"]) for g in annotations(Path(gff))]
    return pd.DataFrame(rows).drop_duplicates("feature_id")


def pool_label(pool: pd.Series) -> pd.Series:
    """Display label of a pool: "pool3" for a bare number, the pool's own name otherwise (e.g. "Cg1")."""
    text = pool.astype(str)
    return text.where(~text.str.fullmatch(r"\d+"), "pool"+text)


def _rank_hits(frame: pd.DataFrame, s: Settings) -> pd.DataFrame:
    """Top hits per direction, as in published Tn-seq hit lists. Gale et al. 2020 kept genes called by the z-score with
    at least a 2-fold change and reported them by z; resampling/DESeq2-style lists order significant genes by their
    statistic with a fold-change filter.
    - Candidates: called in the direction (call ends with depleted/enriched) and |log2FC| >= hit_min_log2FC (default 1,
      i.e. 2-fold), with mean abundance >= hit_min_abundance (default 0: the noise model already widens the noise of
      low-count genes, which lowers their z).
    - Order: larger |z| in the direction, then larger |log2FC|, then feature ID (stable).
    Ranks are 1..n among candidates, empty otherwise."""
    lfc, z, a = frame.log2FC.to_numpy(float), frame.z.to_numpy(float), frame.mean_abundance.to_numpy(float)
    call = frame["call"].astype(str).str.lower()
    frame["passes_abundance"] = np.isfinite(a) & (a >= s.hit_min_abundance)
    ok = np.isfinite(lfc) & frame.passes_abundance.to_numpy()
    called_dep, called_enr = call.str.endswith("depleted").to_numpy(), call.str.endswith("enriched").to_numpy()
    for direction, eligible, zdir in (("depleted", ok & called_dep & (lfc <= -s.hit_min_log2FC), np.where(z < 0, -z, 0.)),
                                      ("enriched", ok & called_enr & (lfc >= s.hit_min_log2FC), np.where(z > 0, z, 0.)),
                                      ("both", ok & ((called_dep & (lfc <= -s.hit_min_log2FC)) | (called_enr & (lfc >= s.hit_min_log2FC))), np.abs(z))):
        zdir = np.where(np.isfinite(zdir), zdir, 0.)
        idx = np.flatnonzero(eligible)
        order = idx[np.lexsort((frame.feature_id.to_numpy(str)[idx], -np.abs(lfc[idx]), -zdir[idx]))]
        rank = np.full(len(frame), np.nan); rank[order] = np.arange(1, len(order)+1)
        frame[f"z_strength_{direction}"] = zdir; frame[f"rank_{direction}"] = rank
    return frame


def hit_evidence(summary: pd.DataFrame, by_pair: pd.DataFrame, s: Settings) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-gene evidence for the app's top hits: the effect and z behind each call, the abundance, and how many
    control-treated pairs support or call each direction (combined view), plus the same per pair."""
    pairs = by_pair[np.isfinite(by_pair.log2FC) & np.isfinite(by_pair.z) & (by_pair.parent_cpm+by_pair.treated_cpm > 0)].copy()
    pairs["pool_label"] = pool_label(pairs["pool"])
    dep = (pairs.z < 0) & (pairs.log2FC < 0); enr = (pairs.z > 0) & (pairs.log2FC > 0)
    g = pairs.assign(dep=dep, enr=enr, dep_call=pairs.call == "depleted", enr_call=pairs.call == "enriched").groupby("feature_id")
    counts = g.agg(n_pairs_valid=("z", "size"), n_pairs_depleted_support=("dep", "sum"), n_pairs_enriched_support=("enr", "sum"),
                   n_pairs_depleted_called=("dep_call", "sum"), n_pairs_enriched_called=("enr_call", "sum"))
    join = lambda mask:pairs[mask].groupby("feature_id").pool_label.agg(lambda x:";".join(dict.fromkeys(x)))
    combined = pd.DataFrame(dict(feature_id=summary.feature_id, standard_name=summary.standard_name, common_name=summary.common_name,
                                 call=summary.final_call, log2FC=summary.call_log2FC, z=summary.call_z, mean_abundance=summary.mean_abundance))
    combined = combined.merge(counts, left_on="feature_id", right_index=True, how="left")
    combined["pools_depleted_support"] = join(dep).reindex(combined.feature_id).fillna("").to_numpy()
    combined["pools_enriched_support"] = join(enr).reindex(combined.feature_id).fillna("").to_numpy()
    for c in counts.columns:combined[c] = combined[c].fillna(0).astype(int)
    combined = _rank_hits(combined.reset_index(drop=True), s)
    per_pair = []
    for contrast, x in by_pair.groupby("contrast", sort=False):
        valid = np.isfinite(x.log2FC) & np.isfinite(x.z) & (x.parent_cpm+x.treated_cpm > 0)
        f = pd.DataFrame(dict(contrast=contrast, feature_id=x.feature_id, standard_name=x.standard_name, common_name=x.common_name, call=x.call,
                              log2FC=x.log2FC, z=x.z, mean_abundance=x.A_exp, n_pairs_valid=valid.astype(int),
                              n_pairs_depleted_support=(valid & (x.z < 0) & (x.log2FC < 0)).astype(int),
                              n_pairs_enriched_support=(valid & (x.z > 0) & (x.log2FC > 0)).astype(int),
                              n_pairs_depleted_called=(x.call == "depleted").astype(int), n_pairs_enriched_called=(x.call == "enriched").astype(int)))
        f["pools_depleted_support"] = np.where(f.n_pairs_depleted_support > 0, pool_label(x["pool"]), "")
        f["pools_enriched_support"] = np.where(f.n_pairs_enriched_support > 0, pool_label(x["pool"]), "")
        f.loc[~valid.to_numpy(), "mean_abundance"] = np.nan
        per_pair.append(_rank_hits(f.reset_index(drop=True), s))
    return combined, pd.concat(per_pair, ignore_index=True)


def run(comparison_design: Path, output: Path, sample_files: dict | None = None, normalized: dict | None = None,
        comparisons: list[str] | None = None, settings: Settings | None = None, project_id: str = "", analysis_id: str = "") -> dict:
    """Run the z-score analysis. `sample_files`: sample -> raw SummaryTable (CPM input, the default).
    `normalized`: {"genes": DataFrame (feature_id x sample normalized gene values), "label": method} for any other
    normalization."""
    s = settings or Settings(); out = Path(output); tables = out/"tables"; tables.mkdir(parents=True, exist_ok=True)
    contrasts, sample_map = load_design(Path(comparison_design), comparisons)
    # Library sizes and gene annotation: from the raw SummaryTables when given (always for CPM); otherwise from the
    # Site Normalization run's raw site matrix and the GFF (normalized["raw"], normalized["annotation"]).
    have_tables = bool(sample_files) and all(x in sample_files and Path(sample_files[x]).is_file() for x in sample_map["sample"])
    reads, annot, sizes, sizes_source = None, None, None, ""
    if normalized is None or have_tables:
        sizes_source = "SummaryTable gene reads"
        long, reads = gene_values_from_feature_tables(sample_map, sample_files or {})
        long.to_csv(tables/"feature_table.long.raw_reads.csv", index=False)
        annot_cols = [c for c in long.columns if c not in {"sample", "condition", "background", "pool", "reads"}]
        annot = long[annot_cols].groupby("feature_id", sort=True).agg(first_nonempty).reset_index()
        annot = annot[["feature_id"]+[c for c in annot_cols if c != "feature_id"]]
    elif normalized.get("raw") is not None:
        sizes_source = "Site Normalization run raw site counts, summed per gene"
        reads = normalized["raw"].rename_axis("feature_id").reset_index().melt(id_vars="feature_id", var_name="sample", value_name="reads")
        reads = reads.merge(sample_map, on="sample", how="inner")
    if reads is not None:
        sizes = reads.groupby(["sample", "condition", "background", "pool"], as_index=False).agg(
            total_feature_reads=("reads", "sum"), features_with_reads=("reads", lambda r:int((r > 0).sum())))
        sizes = sizes.sort_values(["background", "pool", "condition"], kind="stable")
    if normalized is None:
        cpm_long = reads.merge(sizes[["sample", "total_feature_reads"]], on="sample", how="left")
        cpm_long["cpm"] = np.where(cpm_long.total_feature_reads > 0, cpm_long.reads/cpm_long.total_feature_reads*1e6, 0.)
        cpm_long.to_csv(tables/"feature_reads_cpm.long.csv", index=False)
        values = cpm_long.pivot_table(index="feature_id", columns="sample", values="cpm", aggfunc="sum", fill_value=0.)
        values = values.loc[sorted(values.index)]
        value_label, value_scale = "CPM", 1.
    else:
        values = normalized["genes"].astype(float).fillna(0.); values = values.loc[sorted(values.index)]
        # Per million of the mean library total: one constant for all libraries, so the method's between-library
        # correction is kept and the pseudocount means the same as with CPM (amendment 1 of the pre-registration).
        totals = values.sum(axis=0); value_scale = 1e6/float(totals.mean()) if totals.mean() > 0 else 1.
        if abs(value_scale-1) < 1e-9:value_scale = 1.
        values = values*value_scale
        if annot is None and normalized.get("annotation") is not None:annot = normalized["annotation"]
        value_label = normalized.get("label", "normalized")
    if annot is None:annot = pd.DataFrame(dict(feature_id=values.index))
    for column in ("standard_name", "common_name"):   # every result names its genes; the ID when no name is known
        if column not in annot:annot[column] = annot.feature_id
        annot[column] = annot[column].where(annot[column].notna() & (annot[column].astype(str).str.strip() != ""), annot.feature_id)
    value_unit = "CPM" if value_label == "CPM" else f"{UNIT_NAMES.get(value_label, value_label)} reads per million"
    # The gene values the z-score used, for every normalization (and the raw reads when known).
    gene_values = values.rename_axis("feature_id").reset_index().melt(id_vars="feature_id", var_name="sample", value_name="value")
    gene_values = gene_values.merge(sample_map, on="sample", how="inner")
    if reads is not None:gene_values = gene_values.merge(reads[["feature_id", "sample", "reads"]], on=["feature_id", "sample"], how="left")
    gene_values.to_csv(tables/"gene_values.long.csv", index=False)
    if sizes is not None:
        sizes = sizes.merge(values.sum(axis=0).rename("value_total").rename_axis("sample").reset_index(), on="sample", how="left")
        sizes.to_csv(tables/"library_sizes.feature_reads.csv", index=False)
    missing = [x for x in sample_map["sample"] if x not in values.columns]
    if missing:raise ValueError("No gene values for: " + ", ".join(missing))

    models, thresholds, coefs, control_rows = {}, [], [], []
    controls_by_bg = sample_map[sample_map.condition == "parent"].groupby("background", sort=True)["sample"].apply(lambda x:list(dict.fromkeys(x)))
    for bg, controls in controls_by_bg.items():
        if len(controls) < 2:raise ValueError("Each background requires at least two control (parent) samples for the noise model.")
        c1, c2 = controls[:2]; model = control_model(values, c1, c2, s)
        z = model["table"].M_ctrl.to_numpy()/sd_from_fit(model["coef"], model["table"].A_ctrl.to_numpy())
        model["threshold"] = quantile7(np.abs(z), s.z_quantile); model["controls"] = [c1, c2]; models[bg] = model
        thresholds.append(dict(background=bg, z_quantile=s.z_quantile, z_threshold=model["threshold"], parent_comparison=f"{c1}_vs_{c2}"))
        coefs += [dict(background=bg, parameter=k, value=v) for k, v in model["coef"].items()]
        control_rows.append(pd.DataFrame(dict(background=bg, parent_comparison=f"{c1}_vs_{c2}", feature_id=model["table"].feature_id,
                                              A_ctrl=model["table"].A_ctrl, M_ctrl=model["table"].M_ctrl, z=z)))
    pd.DataFrame(thresholds).to_csv(tables/"background_specific_z_thresholds.csv", index=False)
    pd.DataFrame(coefs).to_csv(tables/"background_noise_model_coefficients.csv", index=False)
    pd.concat(control_rows).to_csv(tables/"control_control_z_scores.csv", index=False)

    pairs = []
    for row in contrasts.itertuples():
        if row.background not in models:raise ValueError(f"No control-control noise model for background: {row.background}")
        p, t = values[row.parent_sample].to_numpy(), values[row.treated_sample].to_numpy()
        if s.recenter != "none":t = t*recenter_factor(p, t, s.recenter_trim)
        a = (t+p)/2; m = np.log2((t+s.pseudocount)/(p+s.pseudocount)); model = models[row.background]
        sd = sd_from_fit(model["coef"], a)
        with np.errstate(divide="ignore", invalid="ignore"):z = m/sd
        thr = model["threshold"]
        pairs.append(pd.DataFrame(dict(feature_id=values.index, contrast=row.comparison_id, background=row.background, pool=row.pool,
            parent_sample=row.parent_sample, treated_sample=row.treated_sample, parent_cpm=p, treated_cpm=t, A_exp=a, log2FC=m,
            sd_local=sd, z=z, Zthr=thr, call=call_of(z, thr))))
    by_pair = pd.concat(pairs, ignore_index=True).merge(annot, on="feature_id", how="left")
    by_pair.to_csv(tables/"treated_vs_parent.by_pool.log2fc_z.csv", index=False)

    combined = []   # replicates averaged per background x treated condition
    for (bg, condition), group in contrasts.groupby(["background", "treated_condition"], sort=True):
        controls, treated = list(dict.fromkeys(group.parent_sample)), list(dict.fromkeys(group.treated_sample))
        p, t = values[controls].to_numpy().mean(1), values[treated].to_numpy().mean(1)
        factor = recenter_factor(p, t, s.recenter_trim) if s.recenter != "none" else 1.; t = t*factor
        a = (t+p)/2; m = np.log2((t+s.pseudocount)/(p+s.pseudocount)); model = models[bg]; sd = sd_from_fit(model["coef"], a)
        with np.errstate(divide="ignore", invalid="ignore"):z = m/sd
        combined.append(pd.DataFrame(dict(feature_id=values.index, background=bg, treated_condition=condition, controls=";".join(controls),
            treated=";".join(treated), parent_mean=p, treated_mean=t, recenter_factor=factor, A_exp=a, log2FC=m, sd_local=sd, z=z,
            Zthr=model["threshold"], call=call_of(z, model["threshold"]))))
    combined = pd.concat(combined, ignore_index=True)
    combined.merge(annot, on="feature_id", how="left").to_csv(tables/"treated_vs_parent.combined.log2fc_z.csv", index=False)

    g = by_pair.assign(enr=by_pair.call == "enriched", dep=by_pair.call == "depleted", abs_z=by_pair.z.abs()).groupby("feature_id", sort=True)
    summary = g.agg(mean_log2FC=("log2FC", "mean"), median_log2FC=("log2FC", "median"), sd_log2FC=("log2FC", "std"), mean_z=("z", "mean"),
                    max_abs_z=("abs_z", "max"), n_pairs=("log2FC", "count"), n_enriched_pairs=("enr", "sum"), n_depleted_pairs=("dep", "sum"))
    for column, flag in (("pools_enriched", "enr"), ("pools_depleted", "dep")):
        hits = by_pair[by_pair.call == ("enriched" if flag == "enr" else "depleted")]
        joined = hits.groupby("feature_id").pool.agg(lambda x:";".join(sorted(set(x.astype(str)))))
        summary[column] = joined.reindex(summary.index).fillna("")
    summary["final_call"] = [final_call(int(e), int(d)) for e, d in zip(summary.n_enriched_pairs, summary.n_depleted_pairs)]
    summary = summary.astype(dict(n_pairs=int, n_enriched_pairs=int, n_depleted_pairs=int)).reset_index()
    c = combined.assign(enr=combined.call == "enriched", dep=combined.call == "depleted").groupby("feature_id", sort=True)
    pooled = c.agg(combined_log2FC=("log2FC", "mean"), combined_z=("z", "mean"), combined_enriched=("enr", "sum"), combined_depleted=("dep", "sum"))
    summary = summary.merge(pooled.reset_index(), on="feature_id", how="left")
    if s.combine == "replicates":
        summary["final_call"] = [("mixed" if e and d else "enriched" if e else "depleted" if d else "none")
                                 for e, d in zip(summary.combined_enriched, summary.combined_depleted)]
    # Plot-ready abundances (in value_unit): mean over all design libraries, and the control / treated replicate means.
    design_samples = [x for x in sample_map["sample"] if x in values.columns]
    summary["mean_abundance"] = values[design_samples].mean(axis=1).reindex(summary.feature_id).to_numpy()
    means = combined.groupby("feature_id", sort=True)[["parent_mean", "treated_mean"]].mean()
    summary["control_mean"] = means.parent_mean.reindex(summary.feature_id).to_numpy()
    summary["treated_mean"] = means.treated_mean.reindex(summary.feature_id).to_numpy()
    # The effect and z that final_call is based on (combined for the replicate rule, the pair mean for the pair rule).
    replicates = s.combine == "replicates"
    summary["call_log2FC"] = summary.combined_log2FC if replicates else summary.mean_log2FC
    summary["call_z"] = summary.combined_z if replicates else summary.mean_z
    summary = summary.merge(annot, on="feature_id", how="left").sort_values("mean_log2FC", kind="stable")
    summary.to_csv(tables/"treated_vs_parent.summary_by_feature.csv", index=False)
    hits, hits_by_pair = hit_evidence(summary, by_pair, s)
    hits.to_csv(tables/"top_hit_evidence.csv", index=False); hits_by_pair.to_csv(tables/"top_hit_evidence.by_pair.csv", index=False)
    summary.head(100).to_csv(tables/"top100_depleted_in_treated.csv", index=False)
    summary.sort_values("mean_log2FC", ascending=False, kind="stable").head(100).to_csv(tables/"top100_enriched_in_treated.csv", index=False)
    counts = by_pair.groupby(["background", "pool", "contrast", "call"], sort=False).size().rename("n").reset_index()
    counts["call"] = pd.Categorical(counts.call, ["depleted", "none", "enriched"])
    counts.sort_values(["background", "pool", "call"], kind="stable").to_csv(tables/"call_counts.by_pool.csv", index=False)
    pd.DataFrame([dict(implementation="ytab.fitness.zscore", python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                       scipy=__import__("scipy").__version__)]).to_csv(tables/"analysis_runtime_environment.csv", index=False)
    written = [*(["feature_table.long.raw_reads.csv"] if normalized is None or have_tables else []),
               *(["library_sizes.feature_reads.csv"] if sizes is not None else []), *(["feature_reads_cpm.long.csv"] if normalized is None else []),
               "gene_values.long.csv", "background_specific_z_thresholds.csv", "background_noise_model_coefficients.csv",
               "control_control_z_scores.csv", "treated_vs_parent.by_pool.log2fc_z.csv", "treated_vs_parent.combined.log2fc_z.csv",
               "treated_vs_parent.summary_by_feature.csv", "top100_depleted_in_treated.csv", "top100_enriched_in_treated.csv",
               "top_hit_evidence.csv", "top_hit_evidence.by_pair.csv", "call_counts.by_pool.csv", "analysis_runtime_environment.csv"]
    manifest = dict(status="success", contract=CONTRACT_VERSION, tables=[f"tables/{x}" for x in written], method="local-abundance z-score", project_id=project_id, analysis_id=analysis_id,
                    normalization=value_label, normalization_run=str((normalized or {}).get("run") or ""), value_unit=value_unit,
                    value_scale=value_scale, n_features=int(len(summary)), library_sizes=sizes is not None, library_sizes_source=sizes_source,
                    settings=asdict(s), call_rule=s.combine, comparisons=list(contrasts.comparison_id),
                    backgrounds={bg:dict(controls=m["controls"], threshold=m["threshold"], coefficients=m["coef"], converged=m["converged"])
                                 for bg, m in models.items()},
                    calls={k:int(v) for k, v in summary.final_call.value_counts().items()})
    (out/"zscore_manifest.json").write_text(json.dumps(manifest, indent=2, default=str)+"\n")
    return manifest
