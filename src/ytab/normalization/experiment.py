"""Opt-in, design-driven site normalization evaluation for an experiment.

Does not modify fitness results. The parent catalogue remains immutable; sites
observed only in the experiment are retained in a separately versioned matrix.
"""
import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from ytab.pipeline.project_config import load_project_config
from .methods import normalize, run_choices
from .ttr import read_counts, sha256, write_table


def included(row):
    return str(row.get("include", "true")).strip().lower() in {"1", "true", "yes"}


def shared_log2_change(raw, normalized, a, b):
    """log2(right/left) of normalized values at sites positive in both raw samples.

    Equals the raw ratio times the factor ratio for scalar methods; nonlinear methods (quantile, BGC, aBGC)
    have no single factor, so their normalized values are compared directly. Zeros are preserved by every method.
    """
    x, y = raw[:, a].toarray().ravel(), raw[:, b].toarray().ravel()
    nx, ny = normalized[:, a].toarray().ravel(), normalized[:, b].toarray().ravel()
    shared = (x > 0) & (y > 0) & (nx > 0) & (ny > 0)
    return shared, np.log2(ny[shared] / nx[shared])


def samples_by_id(records):
    return {row["sample"]: row for row in records}


def evaluate(project_config, catalogue, output, methods, trim=.05, target=100, plots=True, sample_list=None, together=False):
    """Select assay samples from comparison_design.csv, never from names; or, with `sample_list`, exactly the samples listed
    (for comparisons between strains or conditions that are not pool-matched pairs: multi-condition tests). With a
    sample list, the sample sheet's role says control (parent/control) or treated, and no pairing is required. together=True (sample list only) normalizes all listed samples as one group
    (the project) instead of within each library, so independently made libraries are on one scale (TRANSIT's
    convention); library_id still records each sample's library."""
    if together and not sample_list:raise ValueError("Normalizing together needs an explicit sample list")
    catalogue, output = Path(catalogue).resolve(), Path(output).resolve()
    if output == catalogue or output in catalogue.parents:
        raise ValueError("Output cannot replace the parent catalogue")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is nonempty; choose a new version")
    cfg = load_project_config(Path(project_config).resolve())
    root = Path(cfg["_repo_root"])
    resolve = lambda path: Path(path) if Path(path).is_absolute() else root / path
    project = resolve(cfg["output_project_dir"])
    sheet = resolve(cfg["sample_sheet"])
    design_path = project / "config/comparison_design.csv"
    with sheet.open() as handle:
        records = list(csv.DictReader(handle))
    samples = {row["sample"]: row for row in records}
    if len(samples) != len(records):
        raise ValueError("Duplicate sample IDs in sample sheet")
    if sample_list:
        missing = [s for s in sample_list if s not in samples_by_id(records)]
        if missing:raise ValueError("Samples not in the sample sheet: " + ", ".join(missing))
        design = []
        role = lambda s: (samples_by_id(records)[s].get("library_role") or samples_by_id(records)[s].get("condition") or "").strip().lower()
        controls = {s for s in sample_list if role(s) in ("parent", "control")}
        treatments = set(sample_list) - controls
    else:
        with design_path.open() as handle:
            design = [row for row in csv.DictReader(handle) if included(row)]
        if not design or len({r["comparison_id"] for r in design}) != len(design):
            raise ValueError("Supply nonempty comparisons with unique comparison IDs")
        controls = {r["parent_sample"] for r in design}
        treatments = {r["treated_sample"] for r in design}
    if controls & treatments:
        raise ValueError("A sample cannot be both a control and treatment in this evaluation")
    selected = sorted(controls | treatments)
    for sample in selected:
        if sample not in samples or not included(samples[sample]):
            raise ValueError(f"Comparison references missing/excluded sample: {sample}")
        if not samples[sample].get("pool", "").strip():
            raise ValueError(f"Explicit pool metadata required: {sample}")
    group = lambda sample: (samples[sample].get("background", ""), samples[sample]["pool"])
    for row in design:
        if group(row["parent_sample"]) != group(row["treated_sample"]):
            raise ValueError("Matched comparisons must share background and pool metadata")
    reference = resolve(cfg["reference"]["fasta"])
    reference_hash = sha256(reference)
    source = json.loads((catalogue / "run_manifest.json").read_text())
    if source.get("status") != "success" or source["reference_sha256"] != reference_hash:
        raise ValueError("Parent catalogue must be complete and use the same reference FASTA")
    with gzip.open(catalogue / "sites.tsv.gz", "rt") as handle:
        sites = [(row["chromosome"], int(row["position_1based"])) for row in csv.DictReader(handle, delimiter="\t")]
    if len(sites) != source["universe_sites"] or len(set(sites)) != len(sites):
        raise ValueError("Parent catalogue sites disagree with its manifest")
    counts, metadata = [], []
    for sample in selected:
        row = samples[sample]
        hit = resolve(row["hit_file"]) if row.get("hit_file") else project / "create_hit_file" / sample / (sample + "_hits.txt")
        values, _ = read_counts(hit)
        counts.append(values)
        metadata.append(dict(observation_id=cfg["project_id"] + "::" + sample, sample=sample,
            display_name=row.get("display_name") or sample, project=cfg["project_id"], dataset=cfg["project_id"],
            library_id="::".join((cfg["project_id"], *group(sample))), pool=row["pool"],
            role="control" if sample in controls else "treated", input=str(hit.resolve()), input_sha256=sha256(hit)))
    parent_set = set(sites)
    additions = sorted(set().union(*(set(values) for values in counts)) - parent_set)
    sites.extend(additions)
    index = {key: i for i, key in enumerate(sites)}
    rows, cols, data = [], [], []
    for col, values in enumerate(counts):
        for key, value in values.items():
            rows.append(index[key]); cols.append(col); data.append(value)
    raw = sparse.csr_matrix((np.asarray(data, dtype=np.int64), (rows, cols)), shape=(len(sites), len(selected)))
    # Fail on invalid method/parameters before producing any files.
    if not methods or len(set(methods)) != len(methods):
        raise ValueError("Choose distinct methods explicitly")
    for method in methods:
        normalize(raw, method, [m["library_id"] for m in metadata], trim, target)
    output.mkdir(parents=True, exist_ok=True)
    matrix_dir = output / "matrix"
    matrix_dir.mkdir()
    sparse.save_npz(matrix_dir / "raw_counts.npz", raw)
    write_table(matrix_dir / "observations.tsv", list(metadata[0]), metadata)
    write_table(matrix_dir / "sites.tsv.gz", ["site_index", "chromosome", "position_1based", "in_parent_catalogue"],
        (dict(site_index=i, chromosome=key[0], position_1based=key[1], in_parent_catalogue=i < len(parent_set)) for i, key in enumerate(sites)))
    manifest = dict(status="success", project=cfg["project_id"], reference_sha256=reference_hash,
        universe_sites=len(sites), parent_catalogue_sites=len(parent_set), experiment_additional_sites=len(additions),
        project_config=str(Path(project_config).resolve()), project_config_sha256=sha256(project_config),
        sample_sheet_sha256=sha256(sheet), design_sha256=sha256(design_path) if design else None, sample_list=sorted(sample_list) if sample_list else None, normalization_group="project" if together else "library_id",
        parent_catalogue=str(catalogue), parent_manifest_sha256=sha256(catalogue / "run_manifest.json"),
        parent_sites_sha256=sha256(catalogue / "sites.tsv.gz"), implementation_sha256=sha256(__file__),
        site_policy="Parent empirical catalogue plus experiment-observed additions; zeros mean undetected, not impossible sites",
        scope="Exploratory site-level validation; no fitness calls, no downstream CPM, no independent-induction claim")
    (matrix_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if design:write_table(output / "comparisons.tsv", list(design[0]), design)
    run_choices(matrix_dir, output / "methods", methods, "project" if together else "library_id", trim=trim, target=target)
    comparison_rows = [dict(comparison_id=r["comparison_id"], kind="treated_control", left=r["parent_sample"], right=r["treated_sample"]) for r in design]
    comparison_rows.extend(dict(comparison_id=a + "_vs_" + b, kind="control_control", left=a, right=b)
        for a, b in itertools.combinations(sorted(controls), 2) if group(a) == group(b))
    raw = raw.tocsc()
    metrics = []
    for method in methods:
        normalized = normalize(raw, method, [m["library_id"] for m in metadata], trim, target)[0].tocsc()
        for row in comparison_rows:
            a, b = selected.index(row["left"]), selected.index(row["right"])
            x, y = raw[:, a].toarray().ravel(), raw[:, b].toarray().ravel()
            observed = (x > 0) | (y > 0)
            # Positive shared sites: no pseudocount or arbitrary normalization-unit artifact.
            shared, change = shared_log2_change(raw, normalized, a, b)
            rank = float(spearmanr(x[observed], y[observed]).statistic) if observed.sum() > 1 and np.ptp(x[observed]) and np.ptp(y[observed]) else None
            metrics.append(dict(method=method, **row, shared_positive_sites=int(shared.sum()), observed_union_sites=int(observed.sum()),
                left_only_sites=int(((x > 0) & (y == 0)).sum()), right_only_sites=int(((x == 0) & (y > 0)).sum()),
                spearman_observed_union=rank,
                median_log2_ratio_shared=float(np.median(change)) if len(change) else None,
                median_absolute_log2_ratio_shared=float(np.median(np.abs(change))) if len(change) else None,
                log2_ratio_iqr_shared=float(np.percentile(change, 75) - np.percentile(change, 25)) if len(change) else None))
    if metrics:write_table(output / "comparison_diagnostics.tsv", list(metrics[0]), metrics)  # none when no pairs share a pool
    if plots and comparison_rows:
        plot_comparisons(output, raw, selected, comparison_rows, methods, metadata, trim, target)
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "implementation_snapshot.py").write_bytes(Path(__file__).read_bytes())
    return manifest


def plot_comparisons(output, raw, selected, comparisons, methods, metadata, trim, target):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(comparisons), len(methods), figsize=(5 * len(methods), 3.5 * len(comparisons)), squeeze=False)
    for col, method in enumerate(methods):
        normalized = normalize(raw, method, [m["library_id"] for m in metadata], trim, target)[0].tocsc()
        for row, comparison in enumerate(comparisons):
            a, b = selected.index(comparison["left"]), selected.index(comparison["right"])
            x, y = raw[:, a].toarray().ravel(), raw[:, b].toarray().ravel()
            shared, change = shared_log2_change(raw, normalized, a, b)
            # Raw geometric-mean abundance is deliberately identical across methods.
            abundance = (np.log10(x[shared]) + np.log10(y[shared])) / 2
            ax = axes[row, col]
            if shared.any():
                ax.hexbin(abundance, change, gridsize=65, bins="log", mincnt=1, cmap="viridis")
            ax.axhline(0, color="red", linestyle="--", linewidth=.7)
            ax.set(title=f"{method}: {comparison['comparison_id']}", xlabel="log10 raw geometric-mean site support", ylabel="log2 right / left (shared positive sites)")
    for row_axes in axes:
        limits = [ax.get_ylim() for ax in row_axes]
        for ax in row_axes:
            ax.set_ylim(min(v[0] for v in limits), max(v[1] for v in limits))
    fig.suptitle("Site-level normalization diagnostics — not gene fitness calls", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, .97))
    fig.savefig(output / "comparison_diagnostics.png", dpi=140)
    plt.close(fig)
