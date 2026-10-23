#!/usr/bin/env python3
"""Treated-vs-control local-abundance z-score (ytab.fitness.zscore): depleted and enriched genes, replicates averaged.

Default input: raw per-sample SummaryTable feature tables, as CPM (--sample-table SAMPLE=FILE for every design sample).
Any site normalization instead: --normalized RUN --normalization METHOD --gff GFF (gene totals of the method's
normalized site counts from a Site Normalization run).

  python scripts/local/ytab.py fitness-zscore --comparison-design output/projects/P/config/comparison_design.csv \\
      --output-dir OUT --sample-table s1=output/projects/P/summary/s1/s1.feature_table.RDF_1.csv ...
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.fitness import contract, zscore  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-id", default=""); p.add_argument("--analysis-id", default="")
    p.add_argument("--comparison-design", required=True, type=Path); p.add_argument("--comparisons", help="Comma-separated comparison IDs")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--sample-table", action="append", default=[], help="SAMPLE=raw SummaryTable feature table (CPM input)")
    p.add_argument("--normalized", type=Path, help="Site Normalization run directory (instead of CPM)")
    p.add_argument("--normalization", help="Method in that run, e.g. betageom or ttr")
    p.add_argument("--gff", type=Path, help="GFF of the reference (with --normalized)")
    p.add_argument("--pseudocount", type=float, default=.5); p.add_argument("--z-quantile", type=float, default=.99)
    p.add_argument("--combine", choices=["replicates", "pairs"], default="replicates",
                   help="replicates: average each condition's replicates, one z per gene (default); pairs: per-pair calls")
    p.add_argument("--recenter", choices=["none", "trimmed-sum"], default="none", help="Scalar recentering of the treated values")
    a = p.parse_args()
    settings = zscore.Settings(pseudocount=a.pseudocount, z_quantile=a.z_quantile, combine=a.combine, recenter=a.recenter)
    comparisons = [x for x in (a.comparisons or "").split(",") if x] or None
    try:
        normalized = None
        if a.normalized:
            if not (a.normalization and a.gff):raise ValueError("--normalized needs --normalization and --gff")
            _, samples = zscore.load_design(a.comparison_design, comparisons)
            names = list(samples["sample"])
            normalized = dict(genes=zscore.gene_values_from_normalization(a.normalized, a.normalization, a.gff, names),
                              raw=zscore.raw_gene_reads_from_run(a.normalized, a.gff, names),
                              annotation=zscore.gene_annotation_from_gff(a.gff), label=a.normalization, run=str(Path(a.normalized).resolve()))
        files = {x.split("=", 1)[0]:x.split("=", 1)[1] for x in a.sample_table}
        manifest = zscore.run(a.comparison_design, a.output_dir, sample_files=files, normalized=normalized, comparisons=comparisons,
                              settings=settings, project_id=a.project_id, analysis_id=a.analysis_id)
    except (ValueError, KeyError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr); return 1
    problems = contract.validate(a.output_dir)
    if problems:
        print("ERROR: the result does not meet the output contract: " + "; ".join(problems), file=sys.stderr); return 1
    print(json.dumps({k:manifest[k] for k in ("status", "normalization", "value_unit", "comparisons", "calls")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
