#!/usr/bin/env python3
"""Compare gene-level raw, TTR and CPM results without changing the fitness workflow."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ytab.normalization.gene_validation import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", required=True, help="Completed normalize-experiment output directory")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-count", type=int, default=20)
    parser.add_argument("--fold-threshold", type=float, default=2)
    parser.add_argument("--existing-fitness", help="Optional existing fitness directory for exact default CPM parity checks")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--site-normalization", help="Optional normalize output with zero-preserving quantile results from the identical raw matrix, ordered samples and library_id groups")
    parser.add_argument("--compare-independent-pools", action="store_true", help="Also compare independent pools sharing an explicit background and comparison_group; between-pool controls are not technical noise")
    args = parser.parse_args()
    try:
        result = evaluate(args.validation, args.output, args.min_count, args.fold_threshold, args.existing_fitness, not args.no_plots, args.compare_independent_pools, args.site_normalization)
    except (ValueError, KeyError, OSError, AssertionError) as error:
        parser.exit(1, f"ERROR: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
