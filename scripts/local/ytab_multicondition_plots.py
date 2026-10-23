#!/usr/bin/env python3
"""Response correlation and response heatmap for a `ytab multicondition` run, adapted from TRANSIT's corrplot and heatmap
(DeJesus et al. 2015).

  response_correlation.png/.tsv  Spearman correlation between conditions' responses (log2FC vs the reference),
                                 split-reference (unbiased) when every library has >= 2 reference replicates,
                                 with gene-bootstrap 95% CIs; on all tested genes and (table only) responsive genes
  heatmap.png / heatmap_matrix.tsv  log2FC vs reference for responsive genes (ANOVA, or all-pairs when ANOVA's
                                 noise warning fired), zero-centred colours, genes clustered on their profiles

Example:
  ytab multicondition-plots --run output/multicondition/RUN --output output/multicondition/RUN_plots
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality.multicondition_plots import run


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, help="Output directory of `ytab multicondition`")
    parser.add_argument("--output", required=True, help="New output directory")
    parser.add_argument("--min-sites", type=int, default=3, help="Heatmap: minimum eligible sites per gene (default 3)")
    parser.add_argument("--top", type=int, help="Heatmap: only the top N responsive genes by ANOVA q")
    parser.add_argument("--pseudocount", type=float, default=1.)
    parser.add_argument("--bootstrap", type=int, default=2000, help="Gene-bootstrap resamples for correlation CIs (0 = none)")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    try:
        manifest = run(args.run, args.output, args.min_sites, args.pseudocount, args.bootstrap, args.seed, args.top)
    except (ValueError, KeyError, OSError) as error:parser.exit(1, f"ERROR: {error}\n")
    print(json.dumps({k:manifest[k] for k in ("status", "reference_condition", "conditions", "responsive_genes", "files", "notes")}, indent=2))
    for note in manifest["notes"]:print(f"NOTE: {note}", file=sys.stderr)


if __name__ == "__main__":main()
