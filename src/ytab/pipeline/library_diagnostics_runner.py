"""Restartable local orchestration for cross-library raw QC diagnostics."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .mapfastq_runner import get_included_samples as _get_included_samples
from .mapfastq_runner import load_project_for_mapping


STATUS_FIELDS = [
    "project_id", "status", "samples", "input_hit_files_count", "output_dir",
    "export_dir", "log_file", "elapsed_seconds", "message",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def load_project_for_library_diagnostics(project_config: Path) -> dict:
    """Load project, sample sheet, and resolved reference configuration."""
    return load_project_for_mapping(project_config)


def get_included_samples(config: dict) -> list[dict]:
    return _get_included_samples(config)


def _paths(config: dict) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project_dir = _resolve(config.get("output_project_dir"), root)
    export_root = _resolve(config.get("output_export_dir"), root)
    if project_dir is None or export_root is None:
        raise ValueError("Project and export output directories must be configured.")
    manifest_dir = project_dir / "manifests" / "library_diagnostics"
    return {
        "project": project_dir,
        "output": project_dir / "library_diagnostics",
        "export": export_root / "qc" / "diagnostics",
        "log": project_dir / "logs" / "library_diagnostics" / "library_diagnostics.log",
        "manifest": manifest_dir / "library_diagnostics_manifest.json",
        "status": manifest_dir / "library_diagnostics_status.csv",
    }


def find_sample_hits_file(sample: dict, config: dict) -> Path:
    name = str(sample["sample"])
    return _paths(config)["project"] / "create_hit_file" / name / f"{name}_hits.txt"


def _select_samples(config: dict, samples: list[str] | None) -> list[dict]:
    included = get_included_samples(config)
    by_name = {str(sample.get("sample") or ""): sample for sample in included}
    names = list(by_name) if samples is None else samples
    if not names:
        raise ValueError("No included samples were selected for library diagnostics.")
    if len(names) != len(set(names)):
        raise ValueError("Selected sample names must be unique.")
    unknown = [name for name in names if name not in by_name]
    if unknown:
        raise ValueError(
            "Selected samples are not included or do not exist: " + ", ".join(unknown)
        )
    return [by_name[name] for name in names]


def _reference_files(config: dict) -> dict[str, Path | None]:
    root = Path(config["_repo_root"])
    reference = config.get("reference") or {}
    return {
        "fasta": _resolve(reference.get("fasta"), root),
        "gff": _resolve(reference.get("gff"), root),
        "gtf": _resolve(reference.get("gtf"), root),
        "feature_table": _resolve(reference.get("feature_table"), root),
        "centromere_bed": _resolve(
            reference.get("centromere_bed") or reference.get("centromeres_bed"), root
        ),
        "trna_bed": _resolve(reference.get("trna_bed"), root),
    }


def validate_library_diagnostics_inputs(
    config: dict, samples: list[str] | None = None
) -> list[str]:
    errors: list[str] = []
    try:
        selected = _select_samples(config, samples)
    except ValueError as exc:
        return [str(exc)]
    for sample in selected:
        name = str(sample.get("sample") or "")
        if not name:
            errors.append("An included sample has no sample name.")
        elif Path(name).name != name:
            errors.append(f"Sample name contains path separators: {name}")
        hits = find_sample_hits_file(sample, config)
        if not hits.is_file():
            errors.append(f"Hit file missing for sample {name}: {hits}")
    return errors


def _resource_warnings(config: dict) -> list[str]:
    reference = _reference_files(config)
    warnings: list[str] = []
    if not reference["fasta"] or not reference["fasta"].is_file():
        warnings.append("Reference FASTA is missing; sequence-bias diagnostics will be skipped.")
    if not reference["centromere_bed"] or not reference["centromere_bed"].is_file():
        warnings.append("Centromere BED is missing; centromere-bias diagnostics will be skipped.")
    if not reference["trna_bed"] or not reference["trna_bed"].is_file():
        warnings.append("tRNA BED is missing; no separate tRNA BED diagnostics will be available.")
    if not reference["gff"] or not reference["gff"].is_file():
        warnings.append("Reference GFF is missing; TSS/TTS/tRNA metaplots will be skipped.")
    return warnings


def build_library_diagnostics_command(
    config: dict, samples: list[str] | None = None, threads: int | None = None
) -> list[str]:
    if threads is not None and int(threads) < 1:
        raise ValueError("Threads must be at least 1.")
    selected = _select_samples(config, samples)
    root = Path(config["_repo_root"])
    reference = _reference_files(config)
    command = [
        sys.executable,
        str(root / "src" / "ytab" / "qc" / "LibraryDiagnostics.py"),
        "--hits-txt",
        *[str(find_sample_hits_file(sample, config)) for sample in selected],
        "--outdir",
        str(_paths(config)["output"]),
    ]
    if reference["fasta"] and reference["fasta"].is_file():
        command.extend(["--fasta", str(reference["fasta"])])
    if reference["centromere_bed"] and reference["centromere_bed"].is_file():
        command.extend(["--centromeres-bed", str(reference["centromere_bed"])])
    if reference["gff"] and reference["gff"].is_file():
        command.extend(["--metaplot-gff", str(reference["gff"])])
    return command


def _detected_outputs(output_dir: Path) -> list[str]:
    if not output_dir.is_dir():
        return []
    return sorted(str(path) for path in output_dir.rglob("*") if path.is_file())


def _write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _write_status(path: Path, manifest: dict) -> None:
    row = {
        "project_id": manifest["project_id"],
        "status": manifest["status"],
        "samples": ";".join(manifest["selected_samples"]),
        "input_hit_files_count": len(manifest["input_hit_files"]),
        "output_dir": manifest["output_dir"],
        "export_dir": manifest["export_dir"],
        "log_file": manifest["log_file"],
        "elapsed_seconds": manifest["elapsed_seconds"],
        "message": manifest.get("error_message") or "; ".join(manifest.get("warnings") or []),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerow(row)


def _collect_exports(output_dir: Path, export_dir: Path) -> None:
    for source in output_dir.rglob("*"):
        if not source.is_file():
            continue
        target = export_dir / source.relative_to(output_dir)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def run_library_diagnostics_project(
    project_config: Path, samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
) -> dict:
    config = load_project_for_library_diagnostics(project_config)
    errors = validate_library_diagnostics_inputs(config, samples)
    if errors:
        raise ValueError("\n".join(errors))
    selected = _select_samples(config, samples)
    selected_names = [str(sample["sample"]) for sample in selected]
    hit_files = [find_sample_hits_file(sample, config) for sample in selected]
    paths = _paths(config)
    command = build_library_diagnostics_command(config, selected_names, threads)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"),
        "species": config.get("species"),
        "selected_samples": selected_names,
        "input_hit_files": [str(path) for path in hit_files],
        "output_dir": str(paths["output"]),
        "export_dir": str(paths["export"]),
        "log_file": str(paths["log"]),
        "status": "failed",
        "start_time": _now(),
        "end_time": None,
        "elapsed_seconds": 0.0,
        "command_run": command,
        "detected_outputs": [],
        "warnings": _resource_warnings(config),
        "error_message": None,
    }

    expected = paths["output"] / "library_diagnostics.summary.csv"
    if not force and paths["manifest"].is_file() and expected.is_file():
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        same_inputs = previous.get("selected_samples") == selected_names
        if previous.get("status") in {"success", "skipped"} and same_inputs:
            manifest["status"] = "skipped"
            manifest["warnings"].append("Successful cached diagnostics found; LibraryDiagnostics skipped.")
    elif dry_run:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Dry run; command was not executed.")
    else:
        paths["output"].mkdir(parents=True, exist_ok=True)
        paths["log"].parent.mkdir(parents=True, exist_ok=True)
        try:
            with paths["log"].open("w", encoding="utf-8") as log:
                log.write("Command: " + " ".join(command) + "\n\n")
                log.flush()
                completed = subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, text=True
                )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"LibraryDiagnostics exited with status {completed.returncode}; see {paths['log']}"
                )
            if not expected.is_file():
                raise RuntimeError("LibraryDiagnostics completed without its combined summary CSV.")
            _collect_exports(paths["output"], paths["export"])
            manifest["status"] = "success"
        except Exception as exc:  # keep a durable failed manifest for app/status readers
            manifest["error_message"] = str(exc)

    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    manifest["detected_outputs"] = _detected_outputs(paths["output"])
    _write_manifest(paths["manifest"], manifest)
    _write_status(paths["status"], manifest)
    return {
        **manifest,
        "included_samples": len(get_included_samples(config)),
        "manifest_path": str(paths["manifest"]),
        "status_path": str(paths["status"]),
        "dry_run": dry_run,
    }
