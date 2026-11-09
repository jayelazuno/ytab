"""Filename-only FASTQ sample discovery."""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

import pandas as pd


COLUMNS = ["sample", "fastq_1", "fastq_2", "layout", "guessed_condition",
           "guessed_background", "guessed_pool", "include", "warnings"]
FASTQ_RE = re.compile(r"(?i)(\.fastq|\.fq)(\.gz)?$")
READ_RE = re.compile(r"(?i)(?P<sample>.+?)(?P<separator>_|\.)(?P<read>R?[12])$")


def _guesses(sample: str) -> tuple[str, str, str]:
    lower = sample.lower()
    condition = ""
    if "parent" in lower:
        condition = "parent"
    elif any(token in lower for token in ("treated", "h2o2", "facs", "selected", "selection")):
        condition = "treated"
    background_match = re.search(r"(?i)(?:y)?H(298|299)", sample)
    background = background_match.group(0) if background_match else ""
    pool_match = re.search(r"(?i)pool[\s_.-]*(\d+)", sample)
    return condition, background, pool_match.group(1) if pool_match else ""


def discover_fastqs(fastq_dir: Path, include_files: Iterable[str] | None = None) -> pd.DataFrame:
    """Group the FASTQs in one directory into samples.

    ``include_files`` names the FASTQ files to analyse, as filenames relative to
    ``fastq_dir``. Samples built from any other file are still recorded, with
    ``include`` set to False, so the discovered set stays complete and a sample
    can be taken up later without rebuilding the project.
    """
    fastq_dir = Path(fastq_dir).expanduser()
    if not fastq_dir.is_dir():
        return pd.DataFrame(columns=COLUMNS)
    wanted = None if include_files is None else {Path(name).name for name in include_files}
    files = sorted(path.resolve() for path in fastq_dir.iterdir() if path.is_file() and FASTQ_RE.search(path.name))
    groups: dict[str, dict[str, list[Path]]] = {}
    for path in files:
        stem = FASTQ_RE.sub("", path.name)
        match = READ_RE.fullmatch(stem)
        if match:
            sample = match.group("sample")
            read = match.group("read")[-1]
            groups.setdefault(sample, {"1": [], "2": [], "single": []})[read].append(path)
        else:
            groups.setdefault(stem, {"1": [], "2": [], "single": []})["single"].append(path)

    rows = []
    for sample, reads in sorted(groups.items()):
        warnings = []
        paired = bool(reads["1"] or reads["2"])
        if paired and (not reads["1"] or not reads["2"]):
            warnings.append("Paired-end filename is missing R1 or R2.")
        if any(len(reads[key]) > 1 for key in reads):
            warnings.append("Multiple FASTQs map to the same sample/read; first file selected.")
        fastq_1 = (reads["1"] or reads["single"] or [None])[0]
        fastq_2 = (reads["2"] or [None])[0]
        condition, background, pool = _guesses(sample)
        if wanted is None:
            include = True
        else:
            names = {path.name for path in (fastq_1, fastq_2) if path}
            include = bool(names & wanted)
        rows.append({"sample": sample, "fastq_1": str(fastq_1) if fastq_1 else "",
                     "fastq_2": str(fastq_2) if fastq_2 else "", "layout": "paired" if paired else "single",
                     "guessed_condition": condition, "guessed_background": background,
                     "guessed_pool": pool, "include": include, "warnings": " ".join(warnings)})
    return pd.DataFrame(rows, columns=COLUMNS)


def write_sample_sheet(df: pd.DataFrame, out_csv: Path) -> None:
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
