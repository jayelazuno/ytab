#!/usr/bin/env python3
"""Initialize a YTAB project from existing CreateHitFile outputs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.project_config import _reference_dict, _version  # existing config helpers
from ytab.pipeline.reference_registry import resolve_reference
from ytab.pipeline.project_status import build_project_status, write_project_status


MAPFASTQ_FIELDS = ["sample", "status", "layout", "fastq_1", "fastq_2", "output_dir", "log_file", "elapsed_seconds", "message"]
CREATE_HIT_FIELDS = ["sample", "status", "input_bam", "hits_file", "output_dir", "log_file", "elapsed_seconds", "message"]
DESIGN_FIELDS = ["comparison_id", "pair_id", "parent_sample", "treated_sample", "background", "pool", "parent_condition", "treated_condition", "include", "notes"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--hit-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--sample-metadata", type=Path)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--force", action="store_true", help="Overwrite config/manifests only; never delete hit outputs.")
    return parser.parse_args()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def portable(path: Path, root: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def resolve_path(path: Path, root: Path) -> Path:
    path = Path(path).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def validate_project_id(project_id: str) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_-]+", project_id or ""):
        raise ValueError("Project ID may contain only letters, numbers, hyphens, and underscores.")


def discover_hit_files(hit_dir: Path) -> list[dict]:
    paths = sorted(path for path in hit_dir.rglob("*_hits.txt") if path.is_file() and path.stat().st_size > 0)
    rows = []
    for path in paths:
        sample = re.sub(r"_hits\.txt$", "", path.name)
        if path.parent.name and path.parent.name == sample:
            sample = path.parent.name
        if Path(sample).name != sample:
            raise ValueError(f"Sample name contains path separators: {sample}")
        rows.append({"sample": sample, "hit_file": path})
    dedup = {}
    for row in rows:
        dedup.setdefault(row["sample"], row)
    if len(dedup) != len(rows):
        repeated = sorted(sample for sample in dedup if sum(r["sample"] == sample for r in rows) > 1)
        raise ValueError("Multiple hit files resolve to the same sample: " + ", ".join(repeated))
    if not rows:
        raise ValueError(f"No non-empty *_hits.txt files found under {hit_dir}")
    return rows


def infer_metadata(sample: str) -> dict:
    lower = sample.lower()
    condition = "parent" if "parent" in lower else "treated" if any(x in lower for x in ("treated", "h2o2", "facs", "selected", "selection")) else ""
    background = ""
    hit = re.search(r"(?i)(?:^|[-_])(y?H\d+)(?:[-_]|$)", sample)
    if hit:
        background = hit.group(1)
    pool = ""
    pool_hit = re.search(r"(?i)(?:^|[-_])pool[-_]?([A-Za-z0-9.]+)", sample)
    if pool_hit:
        pool = pool_hit.group(1)
    treatment = "H2O2" if "h2o2" in lower else ""
    return {"condition": condition, "background": background, "pool": pool, "treatment": treatment}


def read_metadata(path: Path | None, root: Path) -> dict[str, dict]:
    if path is None:
        return {}
    resolved = resolve_path(path, root)
    if not resolved.is_file():
        raise FileNotFoundError(f"Sample metadata CSV not found: {resolved}")
    with resolved.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "sample" not in rows[0]:
        raise ValueError("Sample metadata CSV must include a sample column.")
    return {row["sample"]: row for row in rows}


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in rows)


def build_sample_rows(hit_rows: list[dict], metadata: dict[str, dict], root: Path, project_dir: Path, species: str) -> list[dict]:
    rows = []
    for row in hit_rows:
        sample = row["sample"]
        inferred = infer_metadata(sample)
        supplied = metadata.get(sample, {})
        condition = supplied.get("condition") or supplied.get("guessed_condition") or inferred["condition"]
        background = supplied.get("background") or supplied.get("guessed_background") or inferred["background"]
        pool = supplied.get("pool") or supplied.get("guessed_pool") or inferred["pool"]
        hit_file = row["hit_file"]
        rows.append({
            "sample": sample,
            "fastq_1": "",
            "fastq_2": "",
            "layout": supplied.get("layout") or "imported_hit_file",
            "condition": condition,
            "background": background,
            "pool": pool,
            "treatment": supplied.get("treatment") or inferred["treatment"],
            "species": supplied.get("species") or species,
            "hit_file": portable(hit_file, root),
            "guessed_condition": condition,
            "guessed_background": background,
            "guessed_pool": pool,
            "include": supplied.get("include") if supplied.get("include") not in (None, "") else "true",
            "import_status": "imported_success",
            "warnings": supplied.get("warnings") or "",
        })
    return rows


def build_comparison_design(samples: list[dict]) -> list[dict]:
    parents = [row for row in samples if str(row.get("condition", "")).lower() == "parent" and str(row.get("include", "")).lower() in {"true", "1", "yes", "y"}]
    treated = [row for row in samples if str(row.get("condition", "")).lower() == "treated" and str(row.get("include", "")).lower() in {"true", "1", "yes", "y"}]
    rows = []
    for child in treated:
        background = str(child.get("background") or "")
        pool = str(child.get("pool") or "")
        matches = [parent for parent in parents if str(parent.get("background") or "") == background and str(parent.get("pool") or "") == pool]
        parent = matches[0] if len(matches) == 1 else {}
        pair = f"{background}_pool{pool}" if background and pool else f"pair_{len(rows)+1}"
        rows.append({
            "comparison_id": f"{pair}_treated_vs_parent",
            "pair_id": pair,
            "parent_sample": parent.get("sample", ""),
            "treated_sample": child.get("sample", ""),
            "background": background,
            "pool": pool,
            "parent_condition": "parent",
            "treated_condition": "treated",
            "include": "true",
            "notes": "" if len(matches) == 1 else f"Expected one matching parent; found {len(matches)}.",
        })
    return rows


def main() -> int:
    args = parse_args()
    try:
        root = args.repo_root.expanduser().resolve()
        validate_project_id(args.project_id)
        if args.threads < 1:
            raise ValueError("Threads must be at least 1.")
        project_dir = root / "output" / "projects" / args.project_id
        export_dir = root / "output" / "exports" / args.project_id
        config_dir = project_dir / "config"
        config_path = config_dir / "project.yaml"
        if config_path.exists() and not args.force:
            raise FileExistsError(f"Project config already exists: {config_path}. Use --force to refresh import metadata without deleting outputs.")
        hit_dir = resolve_path(args.hit_dir, root)
        hit_rows = discover_hit_files(hit_dir)
        metadata = read_metadata(args.sample_metadata, root)
        reference = _reference_dict(resolve_reference(args.species, root), root)
        for directory in (project_dir, export_dir, project_dir / "logs", config_dir, project_dir / "manifests" / "mapfastq", project_dir / "manifests" / "create_hit_file"):
            directory.mkdir(parents=True, exist_ok=True)

        samples = build_sample_rows(hit_rows, metadata, root, project_dir, args.species)
        sample_path = config_dir / "sample_sheet.csv"
        sample_fields = ["sample", "fastq_1", "fastq_2", "layout", "condition", "background", "pool", "treatment", "species", "hit_file", "guessed_condition", "guessed_background", "guessed_pool", "include", "import_status", "warnings"]
        write_csv(sample_path, samples, sample_fields)
        (config_dir / "reference_resolved.json").write_text(json.dumps(reference, indent=2) + "\n", encoding="utf-8")
        design = build_comparison_design(samples)
        write_csv(config_dir / "comparison_design.csv", design, DESIGN_FIELDS)
        write_csv(config_dir / "comparison_design_issues.csv", [row for row in design if row.get("notes")], ["comparison_id", "pair_id", "issue", "message"])

        config = {
            "project_id": args.project_id,
            "display_name": args.project_id,
            "project_type": "release" if args.project_id == "H2O2_screen_v1" else "user",
            "show_in_launcher": True,
            "repo_root": str(root),
            "start_stage": "create_hit_file",
            "entry_mode": "import_existing_hit_files",
            "fastq_dir": "",
            "hit_file_import_dir": portable(hit_dir, root),
            "output_project_dir": portable(project_dir, root),
            "output_export_dir": portable(export_dir, root),
            "species": args.species,
            "threads": args.threads,
            "reference": reference,
            "samples": samples,
            "sample_sheet": portable(sample_path, root),
            "created_at": now(),
            "ytab_version_or_git_commit_if_available": _version(root),
        }
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

        map_rows = [{
            "sample": row["sample"], "status": "imported_or_not_required", "layout": row["layout"],
            "fastq_1": "", "fastq_2": "", "output_dir": "", "log_file": "",
            "elapsed_seconds": "0", "message": "Project starts from existing CreateHitFile outputs; MapFastq is not required.",
        } for row in samples]
        hit_status_rows = [{
            "sample": row["sample"], "status": "imported_success", "input_bam": "",
            "hits_file": row["hit_file"], "output_dir": str(Path(row["hit_file"]).parent if Path(row["hit_file"]).is_absolute() else (root / row["hit_file"]).parent),
            "log_file": "", "elapsed_seconds": "0", "message": "Existing CreateHitFile output imported.",
        } for row in samples]
        write_csv(project_dir / "manifests" / "mapfastq" / "mapfastq_status.csv", map_rows, MAPFASTQ_FIELDS)
        write_csv(project_dir / "manifests" / "create_hit_file" / "create_hit_file_status.csv", hit_status_rows, CREATE_HIT_FIELDS)
        import_manifest = {
            "project_id": args.project_id, "status": "success", "start_stage": "create_hit_file",
            "hit_file_count": len(samples), "hit_file_import_dir": portable(hit_dir, root),
            "imported_at": now(), "samples": samples,
            "mapfastq": "imported_or_not_required", "create_hit_file": "imported_success",
        }
        (project_dir / "manifests" / "create_hit_file" / "import_existing_hit_project_manifest.json").write_text(json.dumps(import_manifest, indent=2) + "\n", encoding="utf-8")
        status = build_project_status(config_path)
        write_project_status(config_path, status)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Project ID: {args.project_id}")
    print(f"Start stage: create_hit_file")
    print(f"Imported hit files: {len(samples)}")
    print(f"project.yaml: {config_path}")
    print(f"sample_sheet.csv: {sample_path}")
    print(f"reference_resolved.json: {config_dir / 'reference_resolved.json'}")
    print(f"comparison_design.csv: {config_dir / 'comparison_design.csv'}")
    print("MapFastq: imported_or_not_required")
    print("CreateHitFile: imported_success")
    print("Next stage: SummaryTable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
