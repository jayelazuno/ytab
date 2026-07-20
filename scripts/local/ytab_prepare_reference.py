#!/usr/bin/env python3
"""Prepare FASTA and Bowtie2 indexes for one local YTAB reference."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.project_config import default_threads
from ytab.pipeline.reference_prepare import prepare_reference
from ytab.pipeline.reference_registry import resolve_reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True)
    parser.add_argument("--threads", type=int, default=default_threads())
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        reference = resolve_reference(args.species, args.repo_root)
        annotation = reference.gff or reference.gtf or reference.feature_table
        print(f"Selected species: {args.species}")
        print(f"FASTA path: {reference.fasta or 'missing'}")
        print(f"Annotation path: {annotation or 'missing'}")
        print(
            "Bowtie2 index complete: "
            f"{'yes' if reference.bowtie2_index_complete else 'no'}"
        )
        result = prepare_reference(
            args.species, args.repo_root, threads=args.threads, force=args.force
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if result["bowtie2_index_complete_before"] and not args.force:
        print("Bowtie2 index found; indexing skipped.")
    else:
        print("Bowtie2 index preparation complete.")
    print(f"FASTA index created: {'yes' if result['fasta_index_created'] else 'no'}")
    print(f"Manifest: {result['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
