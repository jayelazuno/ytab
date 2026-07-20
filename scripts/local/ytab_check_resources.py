#!/usr/bin/env python3
"""Report whether local species references are ready for the YTAB pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.reference_registry import (
    list_available_species,
    reference_can_be_prepared,
    reference_is_runnable,
    reference_is_selectable,
    resolve_reference,
)


def yes_no(value: object) -> str:
    return "yes" if value else "no"


def main() -> int:
    resources_dir = REPO_ROOT / "resources" / "species"
    species_names = list_available_species(resources_dir)
    print("Available species:", ", ".join(species_names) or "none")
    if not species_names:
        print("FAIL: no species directories found under resources/species/.")
        return 1

    runnable: list[str] = []
    preparable: list[str] = []
    for species in species_names:
        reference = resolve_reference(species, REPO_ROOT)
        selectable = reference_is_selectable(reference)
        can_prepare = reference_can_be_prepared(reference)
        can_run = reference_is_runnable(reference)
        if can_run:
            runnable.append(species)
        if can_prepare:
            preparable.append(species)
        print(f"\n{species}:")
        fields = (
            ("FASTA found", reference.fasta),
            ("GFF found", reference.gff),
            ("GTF found", reference.gtf),
            ("feature table found", reference.feature_table),
            ("gffutils DB found", reference.gffutils_db),
            ("Bowtie2 index complete", reference.bowtie2_index_complete),
            ("centromere BED found", reference.centromere_bed),
            ("tRNA BED found", reference.trna_bed),
            ("orthology file found", reference.orthology_file),
        )
        for label, value in fields:
            print(f"  {label}: {yes_no(value)}")
        print(f"  selectable: {yes_no(selectable)}")
        print(f"  can_prepare: {yes_no(can_prepare)}")
        print(f"  runnable: {yes_no(can_run)}")
        if reference.fasta and not reference.bowtie2_index_complete:
            print("  Reference can be prepared; run ytab_prepare_reference.py.")
        elif not reference.fasta:
            print("  Reference files must be restored/downloaded/provided.")
        for warning in reference.warnings:
            print(f"  WARNING: {warning}")

    if runnable or preparable:
        if runnable:
            print(f"\nRunnable species: {', '.join(runnable)}")
        if preparable:
            print(f"References that can be prepared: {', '.join(preparable)}")
        print("PASS: at least one reference can run or be prepared.")
        return 0
    print("\nFAIL: no species can run and no species can be prepared.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
