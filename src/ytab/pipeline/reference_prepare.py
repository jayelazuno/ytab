"""Explicit preparation of local reference indexes (never alignment)."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .reference_registry import find_bowtie2_index_prefix, resolve_reference


def _run(command: list[str]) -> None:
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required command not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Command failed with exit status {exc.returncode}: {' '.join(command)}"
        ) from exc


def _index_prefix_complete(index_prefix: Path) -> bool:
    return any(
        all(Path(f"{index_prefix}{suffix}.{extension}").is_file() for suffix in (
            ".1", ".2", ".3", ".4", ".rev.1", ".rev.2"
        ))
        for extension in ("bt2", "bt2l")
    )


def build_bowtie2_index(
    fasta: Path, index_prefix: Path, threads: int = 2, force: bool = False
) -> dict:
    fasta = Path(fasta).expanduser().resolve()
    index_prefix = Path(index_prefix).expanduser().resolve()
    if not fasta.is_file():
        raise FileNotFoundError("Reference FASTA missing. Cannot build Bowtie2 index.")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if _index_prefix_complete(index_prefix) and not force:
        return {"created": False, "skipped": True, "command": None}
    command = ["bowtie2-build", "--threads", str(threads), str(fasta), str(index_prefix)]
    _run(command)
    return {"created": True, "skipped": False, "command": command}


def ensure_fasta_index(fasta: Path, force: bool = False) -> dict:
    fasta = Path(fasta).expanduser().resolve()
    if not fasta.is_file():
        raise FileNotFoundError("Reference FASTA missing. Cannot build Bowtie2 index.")
    fai = Path(f"{fasta}.fai")
    if fai.is_file() and not force:
        return {"created": False, "skipped": True, "path": fai, "command": None}
    command = ["samtools", "faidx", str(fasta)]
    _run(command)
    return {"created": True, "skipped": False, "path": fai, "command": command}


def prepare_reference(
    species: str,
    repo_root: Path | None = None,
    threads: int = 2,
    force: bool = False,
) -> dict:
    info = resolve_reference(species, repo_root)
    if info.fasta is None:
        raise FileNotFoundError("Reference FASTA missing. Cannot build Bowtie2 index.")

    fasta = info.fasta.resolve()
    index_prefix = fasta.with_suffix("")
    commands_run: list[list[str]] = []
    warnings = list(info.warnings)
    complete_before = info.bowtie2_index_complete

    fasta_result = ensure_fasta_index(fasta, force=force)
    if fasta_result["command"]:
        commands_run.append(fasta_result["command"])

    bowtie_result = build_bowtie2_index(
        fasta, index_prefix, threads=threads, force=force
    )
    if bowtie_result["command"]:
        commands_run.append(bowtie_result["command"])

    _, complete_after, _ = find_bowtie2_index_prefix(info.reference_dir)
    if info.gff and not info.gffutils_db:
        warnings.append(
            "gffutils DB missing; DB generation will be added in a later reference-prep step."
        )

    manifest = {
        "species": species,
        "fasta": str(fasta),
        "index_prefix": str(index_prefix),
        "bowtie2_index_complete_before": complete_before,
        "bowtie2_index_complete_after": complete_after,
        "fasta_index_created": fasta_result["created"],
        "commands_run": commands_run,
        "warnings": warnings,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = info.reference_dir / "reference_prepare_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    manifest["manifest"] = str(manifest_path)
    return manifest
