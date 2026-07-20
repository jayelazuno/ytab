#!/usr/bin/env python3
"""Lightweight smoke test for local reference and CPU discovery."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.project_config import default_threads, detected_cpu_count
from ytab.pipeline.reference_registry import list_available_species, resolve_reference


def summarize(species: str):
    reference = resolve_reference(species, REPO_ROOT)
    print(f"{species}: FASTA={reference.fasta or 'not found'}; index_complete={reference.bowtie2_index_complete}; index_type={reference.bowtie2_index_type or 'none'}")
    for warning in reference.warnings:
        print(f"  WARNING: {warning}")
    return reference


def main() -> int:
    try:
        species = list_available_species(REPO_ROOT / "resources" / "species")
        print("Available species:", ", ".join(species))
        if "glabrata" in species:
            reference = summarize("glabrata")
            assert reference.fasta is not None and reference.fasta.is_file(), "glabrata FASTA not found"
            if list(reference.reference_dir.glob("*.bt2")):
                assert reference.bowtie2_index_complete, "glabrata .bt2 files exist but index is incomplete"
        if "albicans" in species:
            reference = summarize("albicans")
            if list(reference.reference_dir.glob("*.bt2")):
                assert reference.bowtie2_index_complete, "albicans .bt2 files exist but index is incomplete"
        detected = detected_cpu_count()
        default = default_threads()
        print(f"Detected CPU count: {detected}")
        print(f"Default thread count: {default}")
        assert default >= 2
        print("PASS: reference discovery smoke test")
        return 0
    except Exception as exc:
        print(f"FAIL: reference discovery smoke test: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
