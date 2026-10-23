"""Output contract of an Essentiality run (scripts/local/ytab_essentiality.py): the files and fields every run has,
whatever the species or libraries. The app only reads these; anything it shows is written here by the pipeline.

`validate(run_dir)` lists what is missing; the essentiality CLI runs it after every run and fails on any problem.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

VERSION = 1
STATES = ("ES", "GD", "NE", "GA")

MANIFEST_FIELDS = ("status", "contract", "project", "species", "libraries", "hmm_run", "models", "hybrid", "hmm", "calls",
                   "agreement", "chromosomes", "files")
GENES = ("ORF", "Name", "chromosome", "hmm_score", "hmm_call")
HYBRID_GENES = ("hybrid_score", "hybrid_essential", "agreement")
# one row per window: the HMM posteriors, the most probable state ("unmappable" for masked windows) and the insertions
WINDOWS = ("chromosome", "start", *(f"P_{s}" for s in STATES), "state", "insertions", "mappable")


def _header(path: Path, sep: str = "\t") -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.readline().rstrip("\n").split(sep)


def validate(run_dir: Path) -> list[str]:
    """Problems with an Essentiality run folder (empty when it meets the contract)."""
    out = Path(run_dir); path = out/"manifest.json"
    if not path.is_file():return ["missing manifest.json"]
    manifest = json.loads(path.read_text())
    problems = [f"manifest.json: missing field {k}" for k in MANIFEST_FIELDS if k not in manifest]
    if manifest.get("contract") != VERSION:problems.append(f"manifest.json: contract {manifest.get('contract')} is not {VERSION}")
    files = manifest.get("files") or {}
    for key, columns in (("genes", GENES+(HYBRID_GENES if manifest.get("hybrid") else ())), ("windows", WINDOWS)):
        name = files.get(key)
        if not name or not (out/name).is_file():problems.append(f"missing {key} file ({name})"); continue
        header = _header(out/name)
        problems += [f"{name}: missing column {c}" for c in columns if c not in header]
    return problems
