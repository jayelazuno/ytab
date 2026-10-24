#!/usr/bin/env python3
"""Evaluate explicit site-normalization choices on design-defined treatments and controls."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ytab.normalization.experiment import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--catalogue", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--method", action="append", required=True)
    parser.add_argument("--trim", type=float, default=.05)
    parser.add_argument("--target", type=float, default=100)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--normalize-together", action="store_true", help="With --sample: one normalization group for all samples (independent libraries)")
    parser.add_argument("--sample", action="append", help="Normalize exactly these samples (repeat; roles from the sample sheet) instead of the comparison design")
    args = parser.parse_args()
    try:
        result = evaluate(args.project_config, args.catalogue, args.output, args.method, args.trim, args.target, not args.no_plots, sample_list=args.sample, together=args.normalize_together)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"ERROR: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
