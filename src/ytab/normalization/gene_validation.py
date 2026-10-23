"""Separate gene-level evaluation of raw, whole-library TTR and feature CPM.

Reuses SummaryTable annotation/counting, verifies raw aggregation, and does not
produce significance calls or change the existing R fitness workflow.
"""
import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import pearsonr, spearmanr

from ytab.core.RangeSet import RangeSet
from ytab.pipeline.project_config import load_project_config
from ytab.summary.SummaryTable import load_feature_db
from .ttr import sha256, write_table


def rows(path, delimiter="\t"):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def summary_rows(path):
    with Path(path).open() as handle:
        first = handle.readline()
        if not first.startswith("RDF,"):
            handle.seek(0)
        result = list(csv.DictReader(handle))
    if not result or not {"standard_name", "reads"} <= result[0].keys():
        raise ValueError(f"Not a raw SummaryTable feature table: {path}")
    if len({r["standard_name"] for r in result}) != len(result):
        raise ValueError("Duplicate feature IDs in raw summary")
    return result


def projection(features, sites):
    """Match SummaryTable's global exon mask within each feature's inclusive span.

    Preserve existing overlapping-feature behavior rather than introducing a
    different CDS assignment rule into this normalization comparison.
    """
    by_chrom = {}
    for i, (chrom, position) in enumerate(sites):
        by_chrom.setdefault(chrom, []).append((position, i))
    union = {}
    for feature in features:
        union.setdefault(feature.chromosome, []).extend(feature.exons)
    union = {chrom: list(RangeSet(blocks)) for chrom, blocks in union.items()}
    chromosome_data = {}
    for chrom, entries in by_chrom.items():
        entries.sort()
        positions, indices = np.asarray(entries, dtype=np.int64).T
        coding = np.zeros(len(entries), dtype=bool)
        for start, end in union.get(chrom, []):
            coding[np.searchsorted(positions, start):np.searchsorted(positions, end, side="right")] = True
        chromosome_data[chrom] = positions, indices, coding
    gene_indices, site_indices = [], []
    for i, feature in enumerate(features):
        if feature.chromosome not in chromosome_data:
            continue
        positions, indices, coding = chromosome_data[feature.chromosome]
        a, b = np.searchsorted(positions, feature.start), np.searchsorted(positions, feature.stop, side="right")
        selected = indices[a:b][coding[a:b]]
        gene_indices.extend([i] * len(selected)); site_indices.extend(selected)
    return sparse.csr_matrix((np.ones(len(site_indices), dtype=np.int64), (gene_indices, site_indices)), shape=(len(features), len(sites)))


def method_matrices(raw, ttr_factors):
    raw = np.asarray(raw)
    factors = np.asarray(ttr_factors, dtype=float)
    if raw.ndim != 2 or factors.shape != (raw.shape[1],) or np.any(raw < 0) or np.any(~np.isfinite(raw)) or np.any(raw != np.floor(raw)):
        raise ValueError("Expected nonnegative integer gene counts and one TTR factor per sample")
    totals = raw.sum(axis=0)
    if np.any(totals <= 0) or np.any(~np.isfinite(factors)) or np.any(factors <= 0):
        raise ValueError("Positive gene support and finite positive factors required")
    return {"nonorm": raw.copy(), "ttr": raw * factors, "cpm": raw / totals * 1e6}


def metric(values, mask):
    selected = values[mask]
    return dict(genes=int(mask.sum()), median_log2_ratio=float(np.median(selected)) if len(selected) else None,
        median_absolute_log2_ratio=float(np.median(np.abs(selected))) if len(selected) else None,
        log2_ratio_iqr=float(np.percentile(selected, 75) - np.percentile(selected, 25)) if len(selected) else None)


def correlation(x, y, kind):
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float((spearmanr if kind == "spearman" else pearsonr)(x, y).statistic)


def agreement_scope(left, right, library_lookup, compare_independent_pools=False):
    if library_lookup[left["parent_sample"]] == library_lookup[right["parent_sample"]]:
        return "within_library"
    if (compare_independent_pools and left.get("background") and left.get("comparison_group") and
            left["background"] == right.get("background") and left["comparison_group"] == right.get("comparison_group")):
        return "between_independent_pools"
    return None


def evaluate(validation, output, min_count=20, fold_threshold=2, existing_fitness=None, plots=True, compare_independent_pools=False, site_normalization=None):
    validation, output = Path(validation).resolve(), Path(output).resolve()
    if output == validation or output in validation.parents:
        raise ValueError("Output cannot replace validation inputs")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is nonempty; choose a new version")
    if min_count < 1 or not np.isfinite(fold_threshold) or fold_threshold <= 1:
        raise ValueError("Require min_count >= 1 and fold_threshold > 1")
    source = json.loads((validation / "run_manifest.json").read_text())
    if source.get("status") != "success":
        raise ValueError("Site validation must be complete")
    cfg = load_project_config(Path(source["project_config"]))
    root = Path(cfg["_repo_root"])
    resolve = lambda p: Path(p) if Path(p).is_absolute() else root / p
    project = resolve(cfg["output_project_dir"])
    gff, fasta = resolve(cfg["reference"]["gff"]), resolve(cfg["reference"]["fasta"])
    if sha256(fasta) != source["reference_sha256"] or sha256(source["project_config"]) != source["project_config_sha256"]:
        raise ValueError("Reference/config changed since site validation; rebuild validation")
    matrix_dir = validation / "matrix"
    metadata = rows(matrix_dir / "observations.tsv")
    samples = [r["sample"] for r in metadata]
    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate sample metadata")
    for row in metadata:
        if sha256(row["input"]) != row["input_sha256"]:
            raise ValueError("Raw hit inputs changed since site validation")
    raw_sites = sparse.load_npz(matrix_dir / "raw_counts.npz")
    with gzip.open(matrix_dir / "sites.tsv.gz", "rt") as handle:
        sites = [(r["chromosome"], int(r["position_1based"])) for r in csv.DictReader(handle, delimiter="\t")]
    if raw_sites.shape != (len(sites), len(samples)):
        raise ValueError("Matrix dimensions disagree with sites/observations")
    db = load_feature_db(cfg["species"], str(gff), str(fasta))
    features = sorted(db.get_all_features(), key=lambda f: f.standard_name)
    ids = [f.standard_name for f in features]
    mapping = projection(features, sites)
    raw = (mapping @ raw_sites).toarray()
    hashes, parity = {}, []
    annotations = {}
    # Fail closed if existing summaries used different filters/annotation rules.
    for j, sample in enumerate(samples):
        folder = project / "summary" / sample
        candidates = list(folder.glob(sample + ".feature_table.RDF_1.csv"))
        if len(candidates) != 1:
            raise ValueError(f"Need exactly one raw RDF_1 summary for {sample}")
        path = candidates[0]
        records = summary_rows(path)
        lookup = {r["standard_name"]: r for r in records}
        if set(lookup) != set(ids):
            raise ValueError("Existing raw summary feature universe differs from current annotation")
        expected = np.array([float(lookup[gid]["reads"]) for gid in ids])
        if not np.array_equal(raw[:, j], expected):
            raise ValueError(f"Raw aggregation differs from accepted summary for {sample}; inspect filters/annotation before comparing normalization")
        hashes[str(path)] = sha256(path)
        parity.append(dict(sample=sample, features=len(ids), mismatched_features=0, raw_summary_exact_match=True))
        annotations.update(lookup)
    factor_rows = rows(validation / "methods/ttr/factors.tsv")
    factor_lookup = {r["observation_id"]: float(r["factor"]) for r in factor_rows}
    factors = np.array([factor_lookup[r["observation_id"]] for r in metadata])
    methods = method_matrices(raw, factors)
    scales = {"nonorm": np.ones(len(samples)), "ttr": factors, "cpm": 1e6 / raw.sum(axis=0)}
    # Independently verify scaling sites then projecting == scaling raw genes.
    normalized_sites = sparse.load_npz(validation / "methods/ttr/normalized_counts.npz")
    np.testing.assert_allclose((mapping @ normalized_sites).toarray(), methods["ttr"], rtol=1e-12)
    if site_normalization:
        extra = Path(site_normalization).resolve()
        extra_manifest = json.loads((extra / "run_manifest.json").read_text())
        if extra_manifest.get("status") != "success" or extra_manifest.get("group_by") != "library_id":
            raise ValueError("Additional site normalization must be complete and grouped by library_id")
        if extra_manifest.get("observations") != [r["observation_id"] for r in metadata]:
            raise ValueError("Additional site normalization sample selection/order differs")
        for name in ("raw_counts.npz", "observations.tsv", "sites.tsv.gz"):
            if extra_manifest.get("source_hashes", {}).get(name) != sha256(matrix_dir / name):
                raise ValueError("Additional site normalization does not use the same raw matrix/metadata/sites")
        for method in extra_manifest["methods"]:
            if method in methods:
                continue
            if method not in {"quantile","betageom","aBGC","nzmean","totreads","zinfnb"}:
                raise ValueError("Only explicitly supported site normalization methods can be added here")
            path = extra / method / "normalized_counts.npz"
            transformed = sparse.load_npz(path).tocsr()
            transformed.eliminate_zeros()
            canonical = raw_sites.tocsr().copy()
            canonical.eliminate_zeros()
            if (transformed.shape != canonical.shape or np.any(~np.isfinite(transformed.data)) or np.any(transformed.data < 0) or
                    not np.array_equal(transformed.indptr, canonical.indptr) or not np.array_equal(transformed.indices, canonical.indices)):
                raise ValueError("Site normalization must preserve the exact positive site pattern")
            # Nonlinear methods MUST precede aggregation; scalar methods are
            # independently checked against their factor formula and gene sums.
            methods[method] = (mapping @ transformed).toarray()
            if method in {"nzmean","totreads","zinfnb"}:
                from .methods import normalize
                expected,factor_rows=normalize(raw_sites,method,[r["library_id"] for r in metadata])
                expected=expected.tocsr()
                np.testing.assert_allclose(transformed.data,expected.data,rtol=1e-12)
                scales[method]=np.array([r["factor"] for r in factor_rows])
                np.testing.assert_allclose(methods[method],raw*scales[method],rtol=1e-12)
            hashes[str(path)] = sha256(path)
        hashes[str(extra / "run_manifest.json")] = sha256(extra / "run_manifest.json")
    design = rows(validation / "comparisons.tsv")
    library_lookup = {r["sample"]: r["library_id"] for r in metadata}
    comparisons = [dict(comparison_id=r["comparison_id"], kind="treated_control", left=r["parent_sample"], right=r["treated_sample"]) for r in design]
    controls = [r for r in metadata if r["role"] == "control"]
    comparisons.extend(dict(comparison_id=a["sample"] + "_vs_" + b["sample"], kind="control_control", left=a["sample"], right=b["sample"])
        for a, b in itertools.combinations(controls, 2) if a["library_id"] == b["library_id"])
    if compare_independent_pools:
        contexts = {r["sample"]: {(d.get("background"), d.get("comparison_group")) for d in design
            if d["parent_sample"] == r["sample"] and d.get("background") and d.get("comparison_group")} for r in controls}
        comparisons.extend(dict(comparison_id=a["sample"] + "_vs_" + b["sample"], kind="between_pool_controls", left=a["sample"], right=b["sample"])
            for a, b in itertools.combinations(controls, 2) if a["library_id"] != b["library_id"] and
            contexts[a["sample"]] & contexts[b["sample"]])
    effects, diagnostics, agreements = {}, [], []
    threshold = float(np.log2(fold_threshold))
    thresholds = sorted({1, 10, min_count, 50})
    for method, values in methods.items():
        for comp in comparisons:
            a, b = samples.index(comp["left"]), samples.index(comp["right"])
            shared = (raw[:, a] > 0) & (raw[:, b] > 0)
            change = np.full(len(ids), np.nan)
            # A common scalar shift preserves exact raw-ratio ties. Dividing
            # independently scaled floats can otherwise split tied ranks.
            change[shared] = (np.log2(raw[shared, b] / raw[shared, a]) + np.log2(scales[method][b] / scales[method][a])
                if method in scales else np.log2(values[shared, b] / values[shared, a]))
            effects[method, comp["comparison_id"]] = (change, shared, a, b)
            for cutoff in thresholds:
                mask = (raw[:, a] >= cutoff) & (raw[:, b] >= cutoff)
                diagnostics.append(dict(method=method, **comp, minimum_raw_count_each=cutoff, **metric(change, mask),
                    descriptive_large_changes=int((np.abs(change[mask]) >= threshold).sum()),
                    fold_threshold=fold_threshold, left_zero_right_positive=int(((raw[:, a] == 0) & (raw[:, b] > 0)).sum()),
                    right_zero_left_positive=int(((raw[:, b] == 0) & (raw[:, a] > 0)).sum()),
                    both_zero=int(((raw[:, a] == 0) & (raw[:, b] == 0)).sum())))
        # Agreement is descriptive; comparisons from one induced pool are not independent inductions.
        for left, right in itertools.combinations(design, 2):
            scope = agreement_scope(left, right, library_lookup, compare_independent_pools)
            if scope is None:
                continue
            x, mx, a, b = effects[method, left["comparison_id"]]
            y, my, c, d = effects[method, right["comparison_id"]]
            for cutoff in thresholds:
                mask = np.all(raw[:, [a,b,c,d]] >= cutoff, axis=1)
                agreements.append(dict(method=method, comparison1=left["comparison_id"], comparison2=right["comparison_id"],
                    scope=scope,
                    minimum_raw_count_each=cutoff, genes=int(mask.sum()), pearson=correlation(x[mask], y[mask], "pearson"),
                    spearman=correlation(x[mask], y[mask], "spearman"),
                    median_absolute_effect_difference=float(np.median(np.abs(x[mask] - y[mask]))) if mask.any() else None,
                    same_direction_fraction=float(np.mean(np.sign(x[mask]) == np.sign(y[mask]))) if mask.any() else None,
                    descriptive_large_same_direction=int(((np.abs(x[mask]) >= threshold) & (np.abs(y[mask]) >= threshold) & (np.sign(x[mask]) == np.sign(y[mask]))).sum())))
    cpm_parity = None
    if existing_fitness:
        fitness = Path(existing_fitness).resolve()
        exported = fitness / "tables/feature_reads_cpm.long.csv"
        count, seen = 0, set()
        id_index = {gid:i for i,gid in enumerate(ids)}
        for row in rows(exported, ","):
            if row["sample"] not in samples or row["feature_id"] not in annotations:
                raise ValueError("Existing fitness has a different sample/feature scope")
            key = row["feature_id"], row["sample"]
            if key in seen:
                raise ValueError("Duplicate rows in existing fitness CPM export")
            seen.add(key)
            i, j = id_index[row["feature_id"]], samples.index(row["sample"])
            np.testing.assert_allclose(methods["cpm"][i,j], float(row["cpm"]), rtol=1e-12, atol=1e-10)
            count += 1
        if count != raw.size:
            raise ValueError("Existing fitness CPM export is incomplete")
        contrast_path = fitness / "tables/treated_vs_parent.by_pool.log2fc_z.csv"
        count_effects, seen_effects = 0, set()
        for row in rows(contrast_path, ","):
            key = row["feature_id"], row["contrast"]
            if key in seen_effects or not any(r["comparison_id"] == row["contrast"] and r["parent_sample"] == row["parent_sample"] and r["treated_sample"] == row["treated_sample"] for r in design):
                raise ValueError("Existing fitness effects have duplicate rows or different comparison scope")
            seen_effects.add(key)
            i, a, b = id_index[row["feature_id"]], samples.index(row["parent_sample"]), samples.index(row["treated_sample"])
            # Reproduce the current default 0.5 CPM pseudocount as a parity check only.
            expected = np.log2((methods["cpm"][i,b] + .5) / (methods["cpm"][i,a] + .5))
            np.testing.assert_allclose(expected, float(row["log2FC"]), rtol=1e-12, atol=1e-10)
            count_effects += 1
        if count_effects != len(ids) * len(design):
            raise ValueError("Existing fitness effects export is incomplete")
        cpm_parity = dict(exact_cpm_rows=count, exact_default_pseudocount_effect_rows=count_effects, default_cpm_pseudocount=.5)
        hashes.update({str(p):sha256(p) for p in (exported, contrast_path)})
    output.mkdir(parents=True, exist_ok=True)
    for method, values in methods.items():
        write_table(output / (method + "_gene_counts.tsv.gz"), ["feature_id", *samples],
            (dict(feature_id=gid, **dict(zip(samples, values[i]))) for i,gid in enumerate(ids)))
    write_table(output / "gene_annotations.tsv", ["feature_id", "common_name", "chromosome", "start", "stop", "coding_length"],
        (dict(feature_id=f.standard_name, common_name=annotations[f.standard_name].get("common_name", ""),
            chromosome=f.chromosome, start=f.start, stop=f.stop, coding_length=f.coding_length) for f in features))
    write_table(output / "aggregation_parity.tsv", list(parity[0]), parity)
    write_table(output / "comparison_diagnostics.tsv", list(diagnostics[0]), diagnostics)
    if agreements:
        write_table(output / "treatment_agreement.tsv", list(agreements[0]), agreements)
    for nonlinear in [name for name in methods if name not in scales]:
        retention = []
        for left, right in itertools.combinations(design, 2):
            if agreement_scope(left, right, library_lookup, compare_independent_pools) is None:
                continue
            x, _, a, b = effects["cpm", left["comparison_id"]]
            y, _, c, d = effects["cpm", right["comparison_id"]]
            qx = effects[nonlinear, left["comparison_id"]][0]
            qy = effects[nonlinear, right["comparison_id"]][0]
            for cutoff in thresholds:
                mask = np.all(raw[:,[a,b,c,d]] >= cutoff, axis=1) & (np.abs(x) >= threshold) & (np.abs(y) >= threshold) & (np.sign(x) == np.sign(y))
                retained = mask & (np.abs(qx) >= threshold) & (np.abs(qy) >= threshold) & (np.sign(qx) == np.sign(x)) & (np.sign(qy) == np.sign(y))
                changed_sign = mask & ((np.sign(qx) != np.sign(x)) | (np.sign(qy) != np.sign(y)))
                retention.append(dict(comparison1=left["comparison_id"], comparison2=right["comparison_id"], minimum_raw_count_each=cutoff,
                    reference="CPM descriptive consistent large effects; not ground truth", fold_threshold=fold_threshold,
                    cpm_consistent_large_genes=int(mask.sum()), **{nonlinear+"_retained_large_same_direction":int(retained.sum()),
                    nonlinear+"_changed_sign_at_least_one":int(changed_sign.sum()),
                    "median_absolute_effect_ratio_"+nonlinear+"_to_cpm":float(np.median(np.r_[np.abs(qx[mask]/x[mask]),np.abs(qy[mask]/y[mask])])) if mask.any() else None}))
        if retention:
            write_table(output / ("effect_retention.tsv" if nonlinear=="quantile" else nonlinear+"_effect_retention.tsv"), list(retention[0]), retention)
    effect_rows = []
    for (method, comp), (change, shared, a, b) in effects.items():
        for i, gid in enumerate(ids):
            effect_rows.append(dict(method=method, comparison_id=comp, feature_id=gid,
                left_raw=int(raw[i,a]), right_raw=int(raw[i,b]), log2_ratio=float(change[i]) if shared[i] else "",
                evidence="both_zero" if not raw[i,a] and not raw[i,b] else "dropout_or_gain" if not shared[i] else "low_count" if min(raw[i,a],raw[i,b]) < min_count else "shared_supported"))
    write_table(output / "gene_effects.tsv.gz", list(effect_rows[0]), effect_rows)
    write_table(output / "observations.tsv", list(metadata[0]), metadata)
    if plots:
        plot(output, raw, methods, comparisons, effects, min_count, metadata)
    hashes.update({str(p):sha256(p) for p in (gff, fasta, validation / "run_manifest.json", matrix_dir / "raw_counts.npz",
        matrix_dir / "sites.tsv.gz", matrix_dir / "observations.tsv", validation / "comparisons.tsv",
        validation / "methods/ttr/factors.tsv", validation / "methods/ttr/normalized_counts.npz")})
    manifest = dict(status="success", validation=str(validation), project=cfg["project_id"], methods=list(methods), features=len(ids),
        samples=samples, minimum_raw_count_each=min_count, count_sensitivity_cutoffs=thresholds,
        descriptive_fold_threshold=fold_threshold, source_hashes=hashes, implementation_sha256=sha256(__file__),
        compare_independent_pools=compare_independent_pools,
        additional_site_normalization=str(Path(site_normalization).resolve()) if site_normalization else None,
        raw_aggregation_exact=True, scaling_before_after_aggregation_agrees=True, existing_fitness_cpm_parity=cpm_parity,
        annotation_policy="Exact SummaryTable global-exon-mask within inclusive feature span; overlapping features retain existing behavior",
        zero_policy="No pseudocount in cross-method diagnostics; report dropout/gain separately; shared-zero genes excluded",
        limitations=["No p-values, FDR estimates, z-score refits or significant fitness calls", "Large changes in controls are descriptive background excursions, not a measured false-positive rate",
            "Culture contrasts sharing a library do not constitute independent induction pools", "Scalar factors change centering, not shared-positive ranks or spread",
            "Between-pool controls, if requested, measure library differences and are not technical control-noise replicates",
            "Method selection cannot be based on treatment centering alone", "Quantile uses a jointly estimated reference and can suppress real distribution shifts; agreement is not independent validation",
            "CPM effect retention is a sensitivity check, not biological ground truth", "No app defaults or existing fitness results changed"])
    if "betageom" in methods:
        manifest["limitations"].append("BGC changes distribution shape and depends on catalogue zeros; it may suppress genuine enrichment; post-TTR is an explicit method step, not CPM re-normalization")
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "implementation_snapshot.py").write_bytes(Path(__file__).read_bytes())
    return manifest


def plot(output, raw, methods, comparisons, effects, cutoff, metadata):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = {r["sample"]:r["display_name"] for r in metadata}
    fig, axes = plt.subplots(len(comparisons), len(methods), figsize=(5*len(methods), 3.6*len(comparisons)), squeeze=False)
    for col, method in enumerate(methods):
        for row, comp in enumerate(comparisons):
            change, shared, a, b = effects[method, comp["comparison_id"]]
            mask = (raw[:,a] >= cutoff) & (raw[:,b] >= cutoff)
            ax = axes[row,col]
            if mask.any():
                abundance = (np.log10(raw[mask,a]) + np.log10(raw[mask,b])) / 2
                ax.hexbin(abundance, change[mask], gridsize=50, bins="log", mincnt=1)
            ax.axhline(0, color="red", linestyle="--", linewidth=.7)
            ax.set(title=f"{method}: {names[comp['right']]}\nvs {names[comp['left']]}", xlabel="log10 raw geometric-mean gene support", ylabel="log2 right / left")
    for axes_row in axes:
        limits = [ax.get_ylim() for ax in axes_row]
        for ax in axes_row:
            ax.set_ylim(min(v[0] for v in limits), max(v[1] for v in limits))
    fig.suptitle(f"Gene-level normalization diagnostics: raw count ≥{cutoff} in both samples; no significance calls")
    fig.tight_layout(rect=(0,0,1,.96))
    fig.savefig(output / "gene_comparisons.png", dpi=140)
    plt.close(fig)
