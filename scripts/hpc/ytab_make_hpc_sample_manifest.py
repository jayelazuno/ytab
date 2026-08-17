#!/usr/bin/env python3
"""Create HPC-ready Zn_toxicity_screen YTAB manifests and config files.

This script performs only project setup and FASTQ discovery. It does not copy
FASTQs, run mapping, or run any downstream scientific calculation.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Iterable

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.reference_registry import resolve_reference  # noqa: E402

PROJECT_ID = "Zn_toxicity_screen"
PROJECT_TITLE = "Zinc Toxicity Screen"
SEQUENCING_RUN_DEFAULT = "He_26110"
COMPARISON_GROUP = "Zn_1_5mM_vs_mock"
COMPARISON_LABEL = "1.5 mM Zn-treated vs mock"

MANIFEST_COLUMNS = [
    "sample", "original_sample_id", "barcode", "lane", "sequencing_run",
    "species", "read_layout", "r1_fastq", "r2_fastq",
    "scratch_mapfastq_dir", "scratch_bam_dir", "scratch_tmp_dir",
    "project_hit_dir", "condition", "condition_label", "treatment",
    "treatment_label", "treatment_concentration", "treatment_concentration_mM",
    "treatment_unit", "control_or_treated", "library_role", "pool",
    "pool_id", "replicate_id", "original_background", "comparison_group",
    "comparison_label", "classifier_role", "fitness_role", "included",
]

SAMPLE_SHEET_COLUMNS = [
    "sample", "fastq_1", "fastq_2", "layout", "condition", "condition_label",
    "background", "pool", "pool_id", "replicate_id", "treatment",
    "treatment_label", "control_or_treated", "library_role", "classifier_role",
    "fitness_role", "species", "hit_file", "guessed_condition",
    "guessed_background", "guessed_pool", "include", "warnings",
]

COMPARISON_COLUMNS = [
    "comparison_id", "display_comparison", "comparison_group",
    "comparison_label", "control_sample", "treated_sample", "parent_sample",
    "treatment", "treatment_label", "treatment_concentration_mM",
    "control_condition", "treated_condition", "pool", "pool_id",
    "replicate_id", "original_background", "included", "status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", default=PROJECT_ID)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--scratch-work-dir", required=True, type=Path)
    parser.add_argument("--species", default="glabrata")
    parser.add_argument("--read-layout", choices=("paired_end",), default="paired_end")
    parser.add_argument("--repo-root", default=REPO_ROOT, type=Path)
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_submission_sheet(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"sample sheet is missing: {path}")
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for index, row in enumerate(reader, 1):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) < 5:
                fail(f"sample sheet row {index} has fewer than 5 columns: {row}")
            rows.append({
                "sample": row[0].strip(),
                "barcode": row[1].strip(),
                "lane": row[3].strip(),
                "sequencing_run": row[4].strip() or SEQUENCING_RUN_DEFAULT,
            })
    if not rows:
        fail(f"no samples found in {path}")
    return rows


def discover_fastqs(directory: Path) -> list[Path]:
    if not directory.is_dir():
        fail(f"FASTQ directory is missing: {directory}")
    patterns = ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq")
    fastqs: list[Path] = []
    for pattern in patterns:
        fastqs.extend(directory.rglob(pattern))
    return sorted(path.resolve() for path in fastqs if path.is_file())


def _mate(path: Path) -> str | None:
    name = path.name
    if re.search(r"(^|[_.])R1($|[_.])", name, re.IGNORECASE) or "_R1_" in name:
        return "R1"
    if re.search(r"(^|[_.])R2($|[_.])", name, re.IGNORECASE) or "_R2_" in name:
        return "R2"
    return None


def _sample_matches(sample: str, path: Path) -> bool:
    name = path.name
    stem = name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return sample in name or stem.startswith(sample)


def resolve_fastq_pair(sample: str, fastqs: Iterable[Path]) -> tuple[Path, Path]:
    matches = [path for path in fastqs if _sample_matches(sample, path)]
    r1 = [path for path in matches if _mate(path) == "R1"]
    r2 = [path for path in matches if _mate(path) == "R2"]
    if not r1 or not r2:
        fail(
            f"missing paired FASTQs for sample {sample}; "
            f"R1 matches={len(r1)} R2 matches={len(r2)} under FASTQ directory"
        )
    if len(r1) > 1 or len(r2) > 1:
        details = "\n  R1: " + "\n      ".join(str(path) for path in r1)
        details += "\n  R2: " + "\n      ".join(str(path) for path in r2)
        fail(f"ambiguous paired FASTQs for sample {sample}:{details}")
    return r1[0], r2[0]


def infer_metadata(sample: str) -> dict[str, str]:
    pool_match = re.search(r"pool([1-4])", sample)
    pool_id = pool_match.group(1) if pool_match else ""
    background_match = re.match(r"(yH\d+)-", sample)
    background = background_match.group(1) if background_match else ""
    if "1_5mM-Zn-treated" in sample:
        return {
            "condition": "Zn-treated",
            "condition_label": "1.5 mM Zn-treated",
            "treatment": "1_5mM-Zn",
            "treatment_label": "1.5 mM Zn",
            "treatment_concentration": "1.5",
            "treatment_concentration_mM": "1.5",
            "treatment_unit": "mM",
            "control_or_treated": "treated",
            "library_role": "treated",
            "classifier_role": "exclude",
            "fitness_role": "treated",
            "pool": f"pool{pool_id}" if pool_id else "",
            "pool_id": pool_id,
            "replicate_id": pool_id,
            "original_background": background,
        }
    if sample.endswith("-mock") or "mock" in sample:
        return {
            "condition": "mock",
            "condition_label": "mock",
            "treatment": "mock",
            "treatment_label": "mock",
            "treatment_concentration": "0",
            "treatment_concentration_mM": "0",
            "treatment_unit": "mM",
            "control_or_treated": "control",
            "library_role": "control",
            "classifier_role": "classifier_control",
            "fitness_role": "control",
            "pool": f"pool{pool_id}" if pool_id else "",
            "pool_id": pool_id,
            "replicate_id": pool_id,
            "original_background": background,
        }
    fail(f"could not infer Zn/mock role from sample name: {sample}")


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def portable_reference(reference: dict, repo_root: Path) -> dict:
    portable = dict(reference)
    for key, value in list(portable.items()):
        if not isinstance(value, str) or not value:
            continue
        path = Path(value)
        if not path.is_absolute():
            continue
        try:
            portable[key] = str(path.relative_to(repo_root))
        except ValueError:
            portable[key] = value
    return portable


def make_comparisons(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_pool: dict[str, dict[str, dict[str, str]]] = {}
    for row in rows:
        by_pool.setdefault(row["pool_id"], {})[row["fitness_role"]] = row
    comparisons: list[dict[str, str]] = []
    for pool_id in ("1", "2", "3", "4"):
        pair = by_pool.get(pool_id, {})
        control = pair.get("control")
        treated = pair.get("treated")
        if not control or not treated:
            fail(f"missing matched mock/Zn-treated comparison for pool {pool_id}")
        comparisons.append({
            "comparison_id": f"Zn_1_5mM_pool{pool_id}_vs_mock",
            "display_comparison": f"1.5 mM Zn-treated pool{pool_id} vs mock pool{pool_id}",
            "comparison_group": COMPARISON_GROUP,
            "comparison_label": COMPARISON_LABEL,
            "control_sample": control["sample"],
            "treated_sample": treated["sample"],
            "parent_sample": control["sample"],
            "treatment": "1_5mM-Zn",
            "treatment_label": "1.5 mM Zn",
            "treatment_concentration_mM": "1.5",
            "control_condition": "mock",
            "treated_condition": "Zn-treated",
            "pool": f"pool{pool_id}",
            "pool_id": pool_id,
            "replicate_id": pool_id,
            "original_background": treated["original_background"],
            "included": "true",
            "status": "active",
        })
    return comparisons


def main() -> int:
    args = parse_args()
    if args.project_id != PROJECT_ID:
        fail(f"project id must be {PROJECT_ID}, not {args.project_id}")
    repo_root = args.repo_root.expanduser().resolve()
    fastq_dir = args.fastq_dir.expanduser().resolve()
    scratch = args.scratch_work_dir.expanduser()
    project_dir = repo_root / "output" / "projects" / args.project_id
    export_dir = repo_root / "output" / "exports" / args.project_id
    config_dir = project_dir / "config"

    submission_rows = read_submission_sheet(args.sample_sheet)
    fastqs = discover_fastqs(fastq_dir)
    if len(submission_rows) != 8:
        fail(f"expected 8 samples, found {len(submission_rows)}")

    manifest_rows: list[dict[str, str]] = []
    for row in submission_rows:
        sample = row["sample"]
        r1, r2 = resolve_fastq_pair(sample, fastqs)
        meta = infer_metadata(sample)
        manifest_rows.append({
            "sample": sample,
            "original_sample_id": sample,
            "barcode": row["barcode"],
            "lane": row["lane"],
            "sequencing_run": row["sequencing_run"],
            "species": args.species,
            "read_layout": args.read_layout,
            "r1_fastq": str(r1),
            "r2_fastq": str(r2),
            "scratch_mapfastq_dir": str(scratch / "mapfastq" / sample),
            "scratch_bam_dir": str(scratch / "bam" / sample),
            "scratch_tmp_dir": str(scratch / "tmp" / sample),
            "project_hit_dir": str(project_dir / "create_hit_file" / sample),
            "comparison_group": COMPARISON_GROUP,
            "comparison_label": COMPARISON_LABEL,
            "included": "true",
            **meta,
        })

    comparisons = make_comparisons(manifest_rows)
    sample_rows: list[dict[str, str]] = []
    for row in manifest_rows:
        sample_rows.append({
            "sample": row["sample"],
            "fastq_1": row["r1_fastq"],
            "fastq_2": row["r2_fastq"],
            "layout": "paired",
            "condition": row["condition"],
            "condition_label": row["condition_label"],
            "background": row["original_background"],
            "pool": row["pool_id"],
            "pool_id": row["pool_id"],
            "replicate_id": row["replicate_id"],
            "treatment": row["treatment"],
            "treatment_label": row["treatment_label"],
            "control_or_treated": row["control_or_treated"],
            "library_role": row["library_role"],
            "classifier_role": row["classifier_role"],
            "fitness_role": row["fitness_role"],
            "species": row["species"],
            "hit_file": f"output/projects/{PROJECT_ID}/create_hit_file/{row['sample']}/{row['sample']}_hits.txt",
            "guessed_condition": "parent" if row["classifier_role"] == "classifier_control" else "treated",
            "guessed_background": row["original_background"],
            "guessed_pool": row["pool_id"],
            "include": "true",
            "warnings": "",
        })

    reference = portable_reference(resolve_reference(args.species, repo_root).to_dict(), repo_root)
    project_config = {
        "project_id": PROJECT_ID,
        "display_name": PROJECT_ID,
        "project_title": PROJECT_TITLE,
        "species": args.species,
        "read_layout": args.read_layout,
        "sequencing_run": SEQUENCING_RUN_DEFAULT,
        "start_stage": "fastq",
        "entry_mode": "hpc_fastq",
        "hpc_project": True,
        "release_project": False,
        "public_pilot_project": False,
        "software_validation_project": False,
        "smoke_test_project": False,
        "analysis_profile": "all",
        "repo_root": "../../..",
        "fastq_dir": str(fastq_dir),
        "scratch_fastq_dir": str(fastq_dir),
        "scratch_work_dir": str(scratch),
        "output_project_dir": f"output/projects/{PROJECT_ID}",
        "output_export_dir": f"output/exports/{PROJECT_ID}",
        "sample_sheet": f"output/projects/{PROJECT_ID}/config/sample_sheet.csv",
        "hpc_sample_manifest": f"output/projects/{PROJECT_ID}/config/hpc_sample_manifest.csv",
        "comparison_design": f"output/projects/{PROJECT_ID}/config/comparison_design.csv",
        "classifier_input_role": "mock_controls_only",
        "fitness_design": "Zn_treated_vs_mock",
        "treatment_label": "1.5 mM Zn",
        "comparison_label": COMPARISON_LABEL,
        "threads": 8,
        "reference": reference,
        "samples": sample_rows,
    }

    write_csv(config_dir / "hpc_sample_manifest.csv", manifest_rows, MANIFEST_COLUMNS)
    write_csv(config_dir / "sample_sheet.csv", sample_rows, SAMPLE_SHEET_COLUMNS)
    write_csv(config_dir / "comparison_design.csv", comparisons, COMPARISON_COLUMNS)
    config_dir.mkdir(parents=True, exist_ok=True)
    (config_dir / "reference_resolved.json").write_text(json.dumps(reference, indent=2) + "\n", encoding="utf-8")
    (config_dir / "project.yaml").write_text(yaml.safe_dump(project_config, sort_keys=False), encoding="utf-8")

    for directory in (
        project_dir / "logs", project_dir / "manifests", project_dir / "create_hit_file",
        project_dir / "summary", project_dir / "library_diagnostics",
        project_dir / "sample_normalization", project_dir / "summary_normalized",
        project_dir / "combined_hits", project_dir / "summary_combined",
        project_dir / "classifier", project_dir / "treated_vs_parent", export_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    print(f"Project ID: {PROJECT_ID}")
    print(f"Samples: {len(manifest_rows)}")
    print(f"Mock controls: {sum(row['classifier_role'] == 'classifier_control' for row in manifest_rows)}")
    print(f"Zn-treated: {sum(row['fitness_role'] == 'treated' for row in manifest_rows)}")
    print("Paired-end handling: R1 and R2 are recorded separately and passed to MapFastq as --r1/--r2; MapFastq uses Bowtie2 -1/-2 and junction-mate R1 by default.")
    print(f"Project config: {config_dir / 'project.yaml'}")
    print(f"HPC sample manifest: {config_dir / 'hpc_sample_manifest.csv'}")
    print(f"Comparison design: {config_dir / 'comparison_design.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
