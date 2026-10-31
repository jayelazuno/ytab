"""Per-library Tn-seq count statistics (definitions adapted from TRANSIT's tnseq_stats; DeJesus et al. 2015) extended
with existing, explicitly linked Hermes QC."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy import sparse, stats

from ytab.normalization.ttr import sha256, write_table

TNSEQ_STATS_COLUMNS = ["dataset", "density", "mean_ct", "NZmean", "NZmedian", "max_ct",
                   "total_cts", "skewness", "kurtosis", "pickands_tail_index"]


def statistics(positive, universe):
    x = np.asarray(positive, dtype=float)
    if (universe < 1 or universe < len(x) or x.ndim != 1 or
            np.any(~np.isfinite(x)) or np.any(x <= 0)):
        raise ValueError("Supply finite positive counts and a valid nonempty catalogue")
    tail = np.sort(x)[::-1]
    estimates = []
    # Exact upstream zero-based order-statistic indices; implicit trailing zeros.
    for m in range(10, 100):
        if 4*m >= universe:
            continue
        a, b, c = [tail[k] if k < len(tail) else 0. for k in (m, 2*m, 4*m)]
        if a > b > c:
            estimates.append(float(np.log2((a-b)/(b-c))))
    variable = len(x) > 1 and np.ptp(x) > 0
    return dict(density=len(x)/universe, mean_ct=float(x.sum()/universe),
        NZmean=float(x.mean()) if len(x) else None,
        NZmedian=float(np.median(x)) if len(x) else None,
        max_ct=float(x.max()) if len(x) else 0., total_cts=float(x.sum()),
        skewness=float(stats.skew(x, bias=True)) if variable else None,
        kurtosis=float(stats.kurtosis(x, fisher=True, bias=True)) if variable else None,
        pickands_tail_index=float(np.median(estimates)) if estimates else None,
        pti_valid_terms=len(estimates), catalogue_sites=universe, occupied_sites=len(x),
        stats_notes="; ".join(note for note in (
            "No occupied sites" if not len(x) else "",
            "Positive-count moments undefined" if not variable else "",
            "PTI unavailable: too few sites or tied order statistics" if not estimates else "" ) if note))


def run(catalogue, output, diagnostics=None):
    catalogue, output = Path(catalogue).resolve(), Path(output).resolve()
    diagnostics = {key: Path(value).resolve() for key,value in (diagnostics or {}).items()}
    if output == catalogue or output in catalogue.parents:
        raise ValueError("Output must not replace the input catalogue")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose a new, empty output directory")
    inputs = {str(catalogue/name):sha256(catalogue/name) for name in
              ("raw_counts.npz", "observations.tsv", "sites.tsv.gz", "run_manifest.json")}
    manifest = json.loads((catalogue/"run_manifest.json").read_text())
    matrix = sparse.load_npz(catalogue/"raw_counts.npz").tocsc()
    matrix.sum_duplicates(); matrix.eliminate_zeros()
    with (catalogue/"observations.tsv").open() as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    if (manifest.get("status") != "success" or manifest.get("universe_sites") != matrix.shape[0]
            or len(metadata) != matrix.shape[1] or not matrix.shape[1]
            or np.any(~np.isfinite(matrix.data)) or np.any(matrix.data < 0)
            or np.any(matrix.data != np.floor(matrix.data))):
        raise ValueError("Expected a successful raw integer-count catalogue matching its metadata")
    ids = [r["observation_id"] for r in metadata]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate observation IDs")
    if not set(diagnostics) <= {r["project"] for r in metadata}:
        raise ValueError("Diagnostic project keys must occur in catalogue metadata")
    summaries = {}
    for project, root in diagnostics.items():
        path = root/"library_diagnostics.summary.csv"
        inputs[str(path)] = sha256(path)
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        if len({r["sample"] for r in rows}) != len(rows):
            raise ValueError("Duplicate samples in existing QC summary")
        summaries[project] = {r["sample"]:r for r in rows}
    output_rows = []
    for j, meta in enumerate(metadata):
        calculated = statistics(matrix.data[matrix.indptr[j]:matrix.indptr[j+1]], matrix.shape[0])
        row = dict(dataset=meta.get("display_name") or meta["observation_id"], **calculated)
        row.update({key:meta.get(key, "") for key in
                    ("observation_id", "sample", "project", "library_id", "pool", "background", "role")})
        project, sample = meta["project"], meta["sample"]
        existing = summaries.get(project, {}).get(sample)
        row["hermes_qc_status"] = "not supplied" if project not in diagnostics else "sample unavailable"
        if existing:
            if (float(existing["total_reads"]) != calculated["total_cts"] or
                    int(existing["unique_sites"]) != calculated["occupied_sites"]):
                row["hermes_qc_status"] = "count mismatch; saved metrics withheld"
                row["hermes_saved_total_reads"] = existing["total_reads"]
                row["hermes_saved_unique_sites"] = existing["unique_sites"]
                existing = None
        if existing:
            row["hermes_qc_status"] = "reused; totals agree; historical source identity not revalidated"
            for key, value in existing.items():
                if key not in {"sample", "input_path", "total_reads", "unique_sites"}:
                    row["hermes_"+key] = value
            root = diagnostics[project]/sample
            for suffix in ("midlc.csv", "seqbias_2_7.tsv", "centromere_bins.tsv", "centromere_bias.png",
                           "tss_metaplot.tsv", "tts_metaplot.tsv", "trna_metaplot.tsv", "metaplots.png"):
                path = root/(sample+"."+suffix)
                row["hermes_"+suffix.replace(".", "_")+"_file"] = ""
                if path.is_file():
                    inputs[str(path)] = sha256(path)
                    row["hermes_"+suffix.replace(".", "_")+"_file"] = str(path)
                    if suffix == "seqbias_2_7.tsv":
                        with path.open() as handle:
                            for pair in csv.DictReader(handle, delimiter="\t"):
                                # Copy the existing convention, do not reinterpret as adjacent TA.
                                if pair["pair"] not in {a+b for a in "ACGT" for b in "ACGT"}:
                                    raise ValueError("Unexpected sequence-bias pair")
                                for field in ("obs_frac", "bg_frac", "enrichment"):
                                    row["hermes_2_7_"+pair["pair"]+"_"+field] = pair[field]
        output_rows.append(row)
    columns = TNSEQ_STATS_COLUMNS + list(dict.fromkeys(k for r in output_rows for k in r if k not in TNSEQ_STATS_COLUMNS))
    output.mkdir(parents=True, exist_ok=True)
    write_table(output/"tnseq_stats.tsv", columns, output_rows)
    result = dict(status="success", samples=len(metadata), catalogue=str(catalogue),
        reference_sha256=manifest.get("reference_sha256"), universe_sites=matrix.shape[0],
        source_hashes=inputs, implementation_sha256=sha256(__file__), count_state="raw; no normalization",
        statistics_conventions="tnseq_stats columns (adapted from TRANSIT, DeJesus et al. 2015); positive-only biased skew and excess kurtosis; guarded PTI m=10..99, zero-based indices",
        limitations=["Density is empirical catalogue occupancy, not genome saturation",
            "No absolute QC thresholds, automatic exclusions or biological pass/fail calls",
            "Existing Hermes outputs are reused, not recomputed; historical source identity is not revalidated",
            "MidLC remains a diagnostic; this table does not validate classifier targets",
            "Existing +2/+7 sequence-bias convention is copied, not a new strand-aware motif analysis"])
    result["qc_status_counts"] = {status:sum(r["hermes_qc_status"]==status for r in output_rows)
                                  for status in sorted({r["hermes_qc_status"] for r in output_rows})}
    (output/"run_manifest.json").write_text(json.dumps(result, indent=2)+"\n")
    return result
