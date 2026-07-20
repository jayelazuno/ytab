"""YAML-backed local YTAB project initialization and validation."""

from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .reference_registry import ReferenceInfo, resolve_reference
from .sample_discovery import discover_fastqs, write_sample_sheet


def detected_cpu_count() -> int:
    return os.cpu_count() or 2


def default_threads(max_threads: int | None = None) -> int:
    detected = detected_cpu_count()
    if detected < 2:
        raise RuntimeError("YTAB local pipeline requires at least 2 CPU threads.")
    default = max(2, min(4, detected))
    if max_threads is not None:
        default = min(default, max_threads)
        if default < 2:
            raise RuntimeError("YTAB local pipeline requires at least 2 CPU threads.")
    return default


def _root(repo_root: Path | None) -> Path:
    return (Path(repo_root) if repo_root else Path(__file__).resolve().parents[3]).expanduser().resolve()


def _portable(path: Path | None, root: Path) -> str | None:
    if path is None:
        return None
    path = Path(path).resolve()
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _reference_dict(reference: ReferenceInfo, root: Path) -> dict[str, Any]:
    data = reference.to_dict()
    for key, value in data.items():
        if key != "warnings" and value is not None and (key.endswith("dir") or key in {
            "fasta", "feature_table", "gff", "gtf", "gffutils_db", "bowtie2_index_prefix",
            "centromere_bed", "trna_bed", "orthology_file",
        }):
            data[key] = _portable(Path(value), root)
    return data


def _version(root: Path) -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root,
                              text=True, capture_output=True, check=True).stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def create_project_config(project_id, fastq_dir, species, repo_root=None, threads=None) -> dict:
    root = _root(repo_root)
    fastq_path = Path(fastq_dir).expanduser().resolve()
    selected_threads = default_threads() if threads is None else int(threads)
    reference = resolve_reference(str(species), root)
    samples = discover_fastqs(fastq_path)
    project_dir = root / "output" / "projects" / str(project_id)
    export_dir = root / "output" / "exports" / str(project_id)
    config_dir = project_dir / "config"
    for directory in (project_dir, export_dir, project_dir / "logs", config_dir):
        directory.mkdir(parents=True, exist_ok=True)

    sample_path = config_dir / "sample_sheet.csv"
    reference_path = config_dir / "reference_resolved.json"
    config_path = config_dir / "project.yaml"
    write_sample_sheet(samples, sample_path)
    reference_data = _reference_dict(reference, root)
    reference_path.write_text(json.dumps(reference_data, indent=2) + "\n", encoding="utf-8")
    sample_records = samples.where(pd.notnull(samples), None).to_dict(orient="records")
    config = {
        "project_id": str(project_id), "repo_root": str(root), "fastq_dir": _portable(fastq_path, root),
        "output_project_dir": _portable(project_dir, root), "output_export_dir": _portable(export_dir, root),
        "species": str(species), "threads": selected_threads, "reference": reference_data,
        "samples": sample_records, "sample_sheet": _portable(sample_path, root),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "ytab_version_or_git_commit_if_available": _version(root),
    }
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config


def load_project_config(config_path: Path) -> dict:
    with Path(config_path).open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"Project config is not a YAML mapping: {config_path}")
    return loaded


def _resolve_config_path(value: str | None, root: Path) -> Path | None:
    if not value:
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else root / path


def validate_project_config(config: dict) -> list[str]:
    messages: list[str] = []
    root = Path(config.get("repo_root") or _root(None)).expanduser().resolve()
    detected = detected_cpu_count()
    if detected < 2:
        messages.append("ERROR: YTAB local pipeline requires at least 2 CPU threads.")
    try:
        threads = int(config.get("threads", 0))
    except (TypeError, ValueError):
        threads = 0
    if threads < 2:
        messages.append("ERROR: threads must be >= 2 for local YTAB runs")
    if threads > detected:
        messages.append("ERROR: requested threads exceed detected CPU count")
    if threads > 4:
        messages.append("WARNING: Using more than 4 threads may increase memory use on low-memory machines.")

    fastq_dir = _resolve_config_path(config.get("fastq_dir"), root)
    if fastq_dir is None or not fastq_dir.is_dir():
        messages.append("ERROR: FASTQ directory does not exist.")
    samples = config.get("samples") or []
    if not samples:
        messages.append("ERROR: no samples were detected.")
    if not any(bool(sample.get("fastq_1") or sample.get("fastq_2")) for sample in samples):
        messages.append("ERROR: no FASTQs were detected.")
    if not any(bool(sample.get("include")) for sample in samples):
        messages.append("ERROR: at least one sample must have include=True.")

    reference = config.get("reference") or {}
    fasta = _resolve_config_path(reference.get("fasta"), root)
    feature = _resolve_config_path(reference.get("feature_table"), root)
    gff = _resolve_config_path(reference.get("gff"), root)
    if fasta is None or not fasta.is_file():
        messages.append("ERROR: reference FASTA is missing.")
    if not ((feature and feature.is_file()) or (gff and gff.is_file())):
        messages.append("ERROR: reference feature table or GFF is required.")
    if not reference.get("bowtie2_index_complete"):
        messages.append("WARNING: No complete Bowtie2 index found. Indexing will be required before alignment.")
    for warning in reference.get("warnings") or []:
        rendered = f"WARNING: {warning}"
        if rendered not in messages:
            messages.append(rendered)

    for key in ("output_project_dir", "output_export_dir"):
        output = _resolve_config_path(config.get(key), root)
        if output is None:
            messages.append(f"ERROR: {key} is not configured.")
            continue
        try:
            output.mkdir(parents=True, exist_ok=True)
            if not os.access(output, os.W_OK):
                messages.append(f"ERROR: {key} is not writable: {output}")
        except OSError as exc:
            messages.append(f"ERROR: {key} is not writable: {exc}")
    return messages
