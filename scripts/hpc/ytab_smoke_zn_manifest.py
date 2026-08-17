#!/usr/bin/env python3
"""Smoke-test Zn_toxicity_screen HPC manifest/config generation."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ID = "Zn_toxicity_screen"
SAMPLES = [
    "yH298-parent-pool1-1_5mM-Zn-treated",
    "yH298-parent-pool1-mock",
    "yH298-parent-pool2-1_5mM-Zn-treated",
    "yH298-parent-pool2-mock",
    "yH299-parent-pool3-1_5mM-Zn-treated",
    "yH299-parent-pool3-mock",
    "yH299-parent-pool4-1_5mM-Zn-treated",
    "yH299-parent-pool4-mock",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", default=PROJECT_ID)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--sample-sheet", type=Path)
    parser.add_argument("--scratch-work-dir", type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--fixture-mode", action="store_true",
                        help="Create temporary fake paired FASTQs/reference files for local smoke validation.")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def local_sample_sheet(repo_root: Path, fastq_dir: Path, explicit: Path | None) -> Path:
    candidates = []
    if explicit:
        candidates.append(explicit)
    candidates.extend([
        fastq_dir / "sampleSheetUsed.csv",
        repo_root / "codex" / "sampleSheetUsed.csv",
        repo_root / "docs" / "codex" / "sampleSheetUsed.csv",
    ])
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit("ERROR: sampleSheetUsed.csv was not found in the supplied paths.")


def make_fixture_repo(source_sheet: Path, fastq_dir: Path) -> tuple[tempfile.TemporaryDirectory, Path, Path, Path]:
    tmp = tempfile.TemporaryDirectory(prefix="ytab_zn_manifest_")
    root = Path(tmp.name) / "repo"
    ref = root / "resources" / "species" / "glabrata" / "reference_genome"
    ref.mkdir(parents=True)
    (ref / "fixture.fna").write_text(">chrI\nACGTACGTACGT\n", encoding="utf-8")
    (ref / "fixture.gff").write_text("##gff-version 3\n", encoding="utf-8")
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        (ref / f"fixture{suffix}").write_text("", encoding="utf-8")
    sheet = root / "codex" / "sampleSheetUsed.csv"
    sheet.parent.mkdir(parents=True)
    sheet.write_text(source_sheet.read_text(encoding="utf-8"), encoding="utf-8")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    for sample in SAMPLES:
        (fastq_dir / f"{sample}_R1.fastq.gz").write_text("", encoding="utf-8")
        (fastq_dir / f"{sample}_R2.fastq.gz").write_text("", encoding="utf-8")
    scratch = Path(tmp.name) / "scratch" / PROJECT_ID
    return tmp, root, sheet, scratch


def validate(repo_root: Path, fastq_dir: Path) -> list[str]:
    project = repo_root / "output" / "projects" / PROJECT_ID
    manifest = project / "config" / "hpc_sample_manifest.csv"
    sample_sheet = project / "config" / "sample_sheet.csv"
    comparisons = project / "config" / "comparison_design.csv"
    config = project / "config" / "project.yaml"
    errors: list[str] = []
    for path in (manifest, sample_sheet, comparisons, config):
        if not path.is_file():
            errors.append(f"missing generated file: {path}")
    if errors:
        return errors
    rows = read_csv(manifest)
    comp = read_csv(comparisons)
    sample_rows = read_csv(sample_sheet)
    text = config.read_text(encoding="utf-8") + manifest.read_text(encoding="utf-8") + comparisons.read_text(encoding="utf-8")

    if PROJECT_ID not in text:
        errors.append("project ID Zn_toxicity_screen missing from generated files")
    forbidden_prefix = "EO46" "_Zn" "_tox"
    if f"{forbidden_prefix}_v1" in text or forbidden_prefix in text:
        errors.append("forbidden legacy Zn project ID/path text found")
    if "H2O2" "_screen_v1" in text:
        errors.append("H2O2 project path/name found in Zn generated files")
    if len(rows) != 8:
        errors.append(f"expected 8 samples, found {len(rows)}")
    if sum(row["classifier_role"] == "classifier_control" for row in rows) != 4:
        errors.append("expected 4 mock classifier controls")
    if sum(row["fitness_role"] == "treated" for row in rows) != 4:
        errors.append("expected 4 Zn-treated fitness samples")
    if {row["pool_id"] for row in rows} != {"1", "2", "3", "4"}:
        errors.append("expected pools 1-4")
    for row in rows:
        sample = row["sample"]
        if row["original_sample_id"] != sample:
            errors.append(f"sample ID changed for {sample}")
        if not Path(row["r1_fastq"]).is_file() or not Path(row["r2_fastq"]).is_file():
            errors.append(f"paired FASTQs missing in manifest for {sample}")
        if "/ytab_work/Zn_toxicity_screen/" not in row["scratch_bam_dir"] and "scratch" not in row["scratch_bam_dir"]:
            errors.append(f"BAM path does not point to scratch for {sample}")
        if not row["project_hit_dir"].endswith(f"output/projects/{PROJECT_ID}/create_hit_file/{sample}"):
            errors.append(f"processed hit dir does not point to Zn project for {sample}")
        if "1_5mM-Zn-treated" in sample:
            if row["treatment_label"] != "1.5 mM Zn" or row["treatment_concentration_mM"] != "1.5":
                errors.append(f"Zn label/concentration incorrect for {sample}")
            if row["classifier_role"] != "exclude" or row["fitness_role"] != "treated":
                errors.append(f"Zn roles incorrect for {sample}")
        if sample.endswith("-mock"):
            if row["treatment_concentration_mM"] != "0" or row["classifier_role"] != "classifier_control":
                errors.append(f"mock role/concentration incorrect for {sample}")
    expected_comparisons = {f"Zn_1_5mM_pool{i}_vs_mock" for i in range(1, 5)}
    if {row["comparison_id"] for row in comp} != expected_comparisons:
        errors.append("comparison IDs do not match expected file-safe Zn_1_5mM IDs")
    if any("1.5 mM Zn" not in row["display_comparison"] for row in comp):
        errors.append("comparison labels do not use human-readable 1.5 mM Zn")
    if any("1_5mM" not in row["comparison_id"] for row in comp):
        errors.append("comparison IDs do not preserve file-safe 1_5mM")
    if any(row["hit_file"].startswith("/") for row in sample_rows):
        errors.append("sample_sheet hit_file should be repo-relative for local app portability")
    return errors


def main() -> int:
    args = parse_args()
    if args.project_id != PROJECT_ID:
        print(f"FAIL: project ID must be {PROJECT_ID}, not {args.project_id}")
        return 1
    repo_root = args.repo_root.resolve()
    source_sheet = local_sample_sheet(repo_root, args.fastq_dir, args.sample_sheet)
    cleanup = None
    fastq_dir = args.fastq_dir.resolve()
    scratch = args.scratch_work_dir or Path("/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/ytab_work/Zn_toxicity_screen")
    if args.fixture_mode:
        cleanup, repo_root, source_sheet, scratch = make_fixture_repo(source_sheet, fastq_dir)
    try:
        project_manifest = repo_root / "output" / "projects" / PROJECT_ID / "config" / "hpc_sample_manifest.csv"
        if args.fixture_mode or not project_manifest.is_file():
            cmd = [
                sys.executable, str(Path(__file__).with_name("ytab_make_hpc_sample_manifest.py")),
                "--project-id", PROJECT_ID,
                "--sample-sheet", str(source_sheet),
                "--fastq-dir", str(fastq_dir),
                "--scratch-work-dir", str(scratch),
                "--species", "glabrata",
                "--read-layout", "paired_end",
                "--repo-root", str(repo_root),
            ]
            completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if completed.returncode != 0:
                print(completed.stdout)
                print("FAIL: manifest generator failed")
                return 1
        errors = validate(repo_root, fastq_dir)
    finally:
        if cleanup is not None:
            cleanup.cleanup()
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
