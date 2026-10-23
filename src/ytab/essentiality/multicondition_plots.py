"""Response correlation and response heatmap for Hermes multi-condition runs, adapted from TRANSIT's corrplot and heatmap (DeJesus et al. 2015).

Reads a `ytab multicondition` output (gene_sample_sums.tsv, anova_genes.tsv and/or pairs_genes.tsv).
Values are the run's normalized, recentered site x sample values; a gene's mean in any set of samples is
sum(S)/sum(N) over its eligible sites.

Adaptations to Hermes split designs:
- TRANSIT's corrplot correlates mean counts. Samples split from one pool share founder abundances, so mean-count
  correlations measure library structure (about 0.93 with no biological effect). Here the correlation is
  between condition *responses*, the log2 fold changes of non-reference conditions against the reference.
- Fold changes against a shared reference correlate through the reference's own noise (simulation, no
  effects: Spearman +0.47). With >= 2 reference replicates in every library, the split-reference correlation
  uses a different reference replicate for each condition. Replicate numbers are also matched (one
  replicate of the condition against one of the reference): a mean of 2 replicates over a single replicate
  has a gene-dependent log bias (Jensen) shared by both conditions, which alone gave +0.11 with no effects.
  The estimate averages over all replicate combinations (at most 64, chosen deterministically). With one
  reference replicate, the shared-reference value is reported and flagged.
- Heatmap: log2 fold change against the pooled reference on a zero-centred scale, for responsive genes
  (ANOVA "varies", or any all-pairs call when ANOVA's noise warning fired or ANOVA was not run), genes with
  >= min_sites eligible sites only (in place of TRANSIT's scale-dependent mean >= 5 filter). Genes are
  clustered on their fold-change profiles (average linkage, Euclidean); conditions are ordered by the
  split-reference correlation when available.
"""
import csv
import textwrap
from itertools import combinations, permutations
import json
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import spearmanr

REVISION = "multicondition_plots_v1"


def load_run(run):
    """Arrays from a multicondition output: S, N (genes x samples), genes, samples, reference, anova/pairs rows."""
    run = Path(run)
    manifest = json.loads((run/"run_manifest.json").read_text())
    if not (run/"gene_sample_sums.tsv").is_file():
        raise ValueError(f"{run} has no gene_sample_sums.tsv; rerun `ytab multicondition` (outputs before this file was added cannot be plotted)")
    samples = [s["observation_id"] for s in manifest["samples"]]
    genes, index = [], {}
    S = N = None
    with (run/"gene_sample_sums.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for r in rows:
        if r["ORF"] not in index:index[r["ORF"]] = len(genes); genes.append(dict(ORF=r["ORF"], Name=r["Name"]))
    S = np.zeros((len(genes), len(samples))); N = np.zeros_like(S); column = {x:j for j, x in enumerate(samples)}
    for r in rows:
        i, j = index[r["ORF"]], column[r["observation_id"]]; S[i, j] = float(r["sum"]); N[i, j] = float(r["eligible_sites"])
    read = lambda name:list(csv.DictReader((run/name).open(), delimiter="\t")) if (run/name).is_file() else None
    return dict(manifest=manifest, S=S, N=N, genes=genes, samples=samples,
                conditions=[s["condition"] for s in manifest["samples"]], libraries=[s["library"] for s in manifest["samples"]],
                reference=manifest["reference_condition"], anova=read("anova_genes.tsv"), pairs=read("pairs_genes.tsv"))


def group_mean(S, N, cols):
    cols = np.asarray(cols, dtype=int); n = N[:, cols].sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, S[:, cols].sum(1)/np.where(n > 0, n, 1), np.nan)


def replicate_sets(conditions, libraries, condition, minimum=1):
    """Replicate r of a condition = its r-th sample in every library; [] when some library has < minimum."""
    per_library = {l:[] for l in dict.fromkeys(libraries)}
    for j, (c, l) in enumerate(zip(conditions, libraries)):
        if c == condition:per_library[l].append(j)
    R = min(len(v) for v in per_library.values())
    return [[v[r] for v in per_library.values()] for r in range(R)] if R >= minimum else []


def reference_replicates(conditions, libraries, reference):
    return replicate_sets(conditions, libraries, reference, minimum=2)


def lfc(a, b, pseudocount):
    return np.log2((a+pseudocount)/(b+pseudocount))


def _rho(x, y):
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3 or np.ptp(x[keep]) == 0 or np.ptp(y[keep]) == 0:return np.nan
    return float(spearmanr(x[keep], y[keep])[0])


def response_correlation(S, N, conditions, libraries, reference, mask=None, pseudocount=1., boot=2000, seed=2026):
    """Spearman correlation between the responses (log2FC vs reference) of each pair of non-reference conditions.
    Returns rows with the split-reference estimate (when available), a gene-bootstrap 95% CI, and the
    shared-reference value for comparison."""
    names = [c for c in dict.fromkeys(conditions) if c != reference]
    cond_cols = {c:[j for j, x in enumerate(conditions) if x == c] for c in dict.fromkeys(conditions)}
    mask = np.ones(len(S), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    splits = reference_replicates(conditions, libraries, reference)
    means = {c:group_mean(S, N, cols) for c, cols in cond_cols.items()}
    ref_split = [group_mean(S, N, cols) for cols in splits]
    reps = {c:[group_mean(S, N, cols) for cols in replicate_sets(conditions, libraries, c)] for c in names}
    rng = np.random.default_rng(seed); rows = []
    for a, b in combinations(names, 2):
        shared = (lfc(means[a], means[reference], pseudocount), lfc(means[b], means[reference], pseudocount))
        if splits and reps[a] and reps[b]:
            combos = [(r1, r2, i, j) for r1, r2 in permutations(range(len(splits)), 2) for i in range(len(reps[a])) for j in range(len(reps[b]))]
            if len(combos) > 64:combos = [combos[k] for k in np.linspace(0, len(combos)-1, 64).round().astype(int)]
            pairs = [(lfc(reps[a][i], ref_split[r1], pseudocount), lfc(reps[b][j], ref_split[r2], pseudocount)) for r1, r2, i, j in combos]
        else:pairs = [shared]
        valid = mask & np.all([np.isfinite(x) & np.isfinite(y) for x, y in pairs+[shared]], axis=0)
        idx = np.flatnonzero(valid)
        estimate = lambda sel:float(np.nanmean([_rho(x[sel], y[sel]) for x, y in pairs]))
        rho = estimate(idx) if len(idx) >= 3 else np.nan
        draws = [estimate(rng.choice(idx, len(idx))) for _ in range(boot)] if len(idx) >= 3 and boot else []
        rows.append(dict(condition_a=a, condition_b=b, genes=int(len(idx)), method="split reference" if len(pairs) > 1 or splits else "shared reference (biased upward)",
                         rho=rho, ci_low=float(np.nanquantile(draws, .025)) if draws else np.nan,
                         ci_high=float(np.nanquantile(draws, .975)) if draws else np.nan,
                         rho_shared_reference=_rho(shared[0][idx], shared[1][idx]) if len(idx) >= 3 else np.nan,
                         reference_splits=len(splits)))
    return rows


def responsive_mask(data):
    """ANOVA 'varies', unless ANOVA's noise warning fired or ANOVA was not run: then any all-pairs call."""
    orf = [g["ORF"] for g in data["genes"]]
    warning = ((data["manifest"].get("diagnostics") or {}).get("anova") or {}).get("warning")
    if data["anova"] is not None and not warning:
        hit = {r["ORF"] for r in data["anova"] if r["call"] == "varies"}; source = "ANOVA (q <= FDR)"
    elif data["pairs"] is not None:
        hit = {r["ORF"] for r in data["pairs"] if r["call_all_pairs"] in ("depleted", "enriched")}
        source = "all-pairs (any pair called)" + (" because ANOVA's replicate-noise warning fired" if warning else "")
    else:raise ValueError("The run has neither ANOVA nor all-pairs results")
    return np.array([o in hit for o in orf]), source


def tested_mask(data):
    rows = data["anova"] or data["pairs"]
    tested = {r["ORF"] for r in rows if r.get("in_family") in ("True", "true", True)}
    return np.array([g["ORF"] in tested for g in data["genes"]])


def heatmap_matrix(S, N, conditions, reference, mask, min_sites=3, pseudocount=1., column_order=None):
    """log2FC vs the pooled reference for masked genes with >= min_sites eligible sites; rows clustered."""
    names = [c for c in dict.fromkeys(conditions) if c != reference]
    if column_order:names = [c for c in column_order if c in names]
    cols = {c:[j for j, x in enumerate(conditions) if x == c] for c in dict.fromkeys(conditions)}
    ref = group_mean(S, N, cols[reference])
    M = np.column_stack([lfc(group_mean(S, N, cols[c]), ref, pseudocount) for c in names])
    keep = np.asarray(mask, dtype=bool) & (N.max(1) >= min_sites) & np.isfinite(M).all(1)
    rows = np.flatnonzero(keep); M = M[rows]
    order = leaves_list(linkage(M, "average", metric="euclidean", optimal_ordering=True)) if len(rows) > 2 else np.arange(len(rows))
    return M[order], rows[order], names


def condition_order(corr_rows, names):
    """Order conditions by average-linkage clustering on 1 - rho (split reference), when >= 3 conditions."""
    if len(names) < 3:return names
    D = np.zeros((len(names), len(names)))
    for r in corr_rows:
        i, j = names.index(r["condition_a"]), names.index(r["condition_b"]); D[i, j] = D[j, i] = 1-(r["rho"] if np.isfinite(r["rho"]) else 0.)
    condensed = D[np.triu_indices(len(names), 1)]
    return [names[i] for i in leaves_list(linkage(condensed, "average", optimal_ordering=True))]


def plot_correlation(rows, names, path, title):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(names); R = np.eye(n); text = [["1" if i == j else "" for j in range(n)] for i in range(n)]
    for r in rows:
        i, j = names.index(r["condition_a"]), names.index(r["condition_b"]); R[i, j] = R[j, i] = r["rho"]
        text[i][j] = text[j][i] = f"{r['rho']:.2f}\n[{r['ci_low']:.2f}, {r['ci_high']:.2f}]" if np.isfinite(r["ci_low"]) else f"{r['rho']:.2f}"
    fig, ax = plt.subplots(figsize=(1.6+1.2*n, 1.2+1.1*n))
    im = ax.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1)
    for i in range(n):
        for j in range(n):ax.text(j, i, text[i][j], ha="center", va="center", fontsize=8)
    ax.set_xticks(range(n), names, rotation=45, ha="right"); ax.set_yticks(range(n), names)
    fig.colorbar(im, ax=ax, label="Spearman rho of responses"); ax.set_title(textwrap.fill(title, 45), fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_heatmap(M, labels, names, path, title, max_labels=80):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    limit = max(1., float(np.quantile(np.abs(M), .99))) if M.size else 1.
    height = min(30., 2.5+.14*max(len(M), 10)); fig, ax = plt.subplots(figsize=(2.5+.9*len(names), height))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit, interpolation="nearest")
    ax.set_xticks(range(len(names)), names, rotation=45, ha="right")
    if len(M) <= max_labels:ax.set_yticks(range(len(M)), labels, fontsize=6)
    else:ax.set_yticks([]); ax.set_ylabel(f"{len(M)} genes (clustered)")
    fig.colorbar(im, ax=ax, label="log2FC vs reference\n(limit: 99th percentile of |log2FC|)", shrink=.6)
    ax.set_title(textwrap.fill(title, 18+12*len(names)), fontsize=8); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def _clean(x):
    if isinstance(x, dict):return {k:_clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):return [_clean(v) for v in x]
    if isinstance(x, (float, np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):return int(x)
    return x


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)


def run(run_dir, output, min_sites=3, pseudocount=1., boot=2000, seed=2026, top=None):
    """Writes response_correlation.tsv/.png (>= 2 non-reference conditions), heatmap_matrix.tsv/heatmap.png and plots_manifest.json."""
    output = Path(output)
    if output.exists() and any(output.iterdir()):raise ValueError("Choose a new output directory")
    data = load_run(run_dir); output.mkdir(parents=True, exist_ok=True)
    S, N, conditions, libraries, reference = data["S"], data["N"], data["conditions"], data["libraries"], data["reference"]
    names = [c for c in dict.fromkeys(conditions) if c != reference]
    tested = tested_mask(data); responsive, source = responsive_mask(data)
    files, notes, corr_all = {}, [], []
    if len(names) >= 2:
        corr_all = response_correlation(S, N, conditions, libraries, reference, tested, pseudocount, boot, seed)
        corr_resp = response_correlation(S, N, conditions, libraries, reference, responsive, pseudocount, boot, seed)
        rows = [r | dict(gene_set="all tested genes") for r in corr_all]+[r | dict(gene_set="responsive genes") for r in corr_resp]
        write_tsv(output/"response_correlation.tsv", ["gene_set", "condition_a", "condition_b", "genes", "method", "rho", "ci_low", "ci_high",
                                                      "rho_shared_reference", "reference_splits"], rows)
        method = corr_all[0]["method"]
        plot_correlation(corr_all, names, output/"response_correlation.png",
                         f"Response correlation, all tested genes ({method}; reference {reference})")
        files |= dict(correlation_table="response_correlation.tsv", correlation_plot="response_correlation.png")
        if method != "split reference":notes.append("Fewer than 2 reference replicates in some library: correlations use a shared reference and are biased upward")
    else:notes.append("Response correlation needs >= 2 non-reference conditions")
    if top:
        q = {r["ORF"]:float(r["q_bh"]) for r in (data["anova"] or []) if r.get("q_bh")}
        ranked = sorted((i for i in np.flatnonzero(responsive)), key=lambda i:q.get(data["genes"][i]["ORF"], 1.))[:top]
        responsive = np.isin(np.arange(len(S)), ranked)
    M, rows, cols = heatmap_matrix(S, N, conditions, reference, responsive, min_sites, pseudocount, condition_order(corr_all, names) if corr_all else names)
    labels = [data["genes"][i]["Name"] or data["genes"][i]["ORF"] for i in rows]
    write_tsv(output/"heatmap_matrix.tsv", ["ORF", "Name"]+[f"log2FC_{c}" for c in cols],
              [dict(ORF=data["genes"][i]["ORF"], Name=data["genes"][i]["Name"]) | {f"log2FC_{c}":float(v) for c, v in zip(cols, m)} for i, m in zip(rows, M)])
    plot_heatmap(M, labels, cols, output/"heatmap.png", f"{len(M)} responsive genes ({source}; >= {min_sites} sites), log2FC vs {reference}")
    files |= dict(heatmap_matrix="heatmap_matrix.tsv", heatmap="heatmap.png")
    manifest = dict(status="success", revision=REVISION, run=str(Path(run_dir).resolve()), reference_condition=reference, conditions=names,
                    responsive_genes=dict(source=source, selected=int(responsive.sum()), plotted=int(len(M))),
                    settings=dict(min_sites=min_sites, pseudocount=pseudocount, bootstrap=boot, seed=seed, top=top),
                    correlation=corr_all, files=files, notes=notes,
                    caveat="Correlations on responsive genes are descriptive: selecting genes on significance can itself induce correlation.")
    (output/"plots_manifest.json").write_text(json.dumps(_clean(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
