#!/usr/bin/env python3
"""Initialize (but do not run) a local YTAB project."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.project_config import create_project_config, default_threads, detected_cpu_count, validate_project_config
from ytab.pipeline.reference_registry import list_available_species


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--force", action="store_true", help="Overwrite existing config files; never delete outputs.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.repo_root.expanduser().resolve()
    config_path = root / "output" / "projects" / args.project_id / "config" / "project.yaml"
    if config_path.exists() and not args.force:
        print(f"ERROR: config already exists: {config_path}\nUse --force to overwrite config files without deleting outputs.", file=sys.stderr)
        return 2
    detected = detected_cpu_count()
    selected = args.threads if args.threads is not None else (default_threads() if detected >= 2 else 0)
    print("Available species:", ", ".join(list_available_species(root / "resources" / "species")) or "none")
    try:
        config = create_project_config(args.project_id, args.fastq_dir, args.species, root, selected)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    reference = config["reference"]
    samples = config["samples"]
    fastq_count = sum(bool(row.get("fastq_1")) + bool(row.get("fastq_2")) for row in samples)
    lines = {
        "Selected species": args.species, "Reference FASTA": reference.get("fasta") or "not found",
        "Feature table": reference.get("feature_table") or "not found", "GFF": reference.get("gff") or "not found",
        "gffutils DB": reference.get("gffutils_db") or "not found",
        "Bowtie2 index prefix": reference.get("bowtie2_index_prefix") or "not found",
        "Bowtie2 index type": reference.get("bowtie2_index_type") or "not found",
        "Bowtie2 index complete": "yes" if reference.get("bowtie2_index_complete") else "no",
        "Indexing skipped": "yes" if reference.get("bowtie2_index_complete") else "no; required later",
        "Detected CPU count": detected, "Selected thread count": selected,
        "Thread selection": "default" if args.threads is None else "user-selected", "FASTQs found": fastq_count,
        "Samples detected": len(samples), "Included samples": sum(bool(row.get("include")) for row in samples),
        "Output project dir": config["output_project_dir"], "Output export dir": config["output_export_dir"],
        "project.yaml": str(config_path), "sample_sheet.csv": str(config_path.with_name("sample_sheet.csv")),
        "reference_resolved.json": str(config_path.with_name("reference_resolved.json")),
    }
    for label, value in lines.items():
        print(f"{label}: {value}")
    messages = validate_project_config(config)
    print("Validation:")
    print("  PASS" if not messages else "\n".join(f"  {message}" for message in messages))
    print("Alignment run: no (Step 1 initialization only)")
    return 1 if any(message.startswith("ERROR:") for message in messages) else 0


if __name__ == "__main__":
    raise SystemExit(main())
