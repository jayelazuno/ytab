#!/usr/bin/env python3
"""Fetch a read subsample for a project's libraries from ENA, at the paths its sample sheet expects.

A project whose reads are public can ship a small, honest sample of them, so the mapping stage can be run without
the full libraries: the first N read pairs of each run, taken straight from ENA with HTTP range requests. Only the
head of each file is transferred, so fetching 200,000 pairs costs a few megabytes per library rather than the whole
run, and nothing has to be cleaned up afterwards.

The subsample is for trying the workflow. It is a small, non-random slice of the start of each file, so its insertion
counts are far below the full library's and its results do not reproduce the project's stored results.

    python scripts/local/ytab_fetch_demo_reads.py --project-config output/projects/<ID>/config/project.yaml

Run accessions come from the sample sheet's `sample` column, and each file is written to the `fastq_1`/`fastq_2`
path recorded there, under the project's `fastq_dir`.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
import zlib
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.project_config import load_project_config  # noqa: E402

ENA_REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
# Measured on Gale et al. 2020: ~43 compressed bytes per read. Ask for more than the estimate so one request
# usually suffices, and grow it if the file turns out to be denser than that.
BYTES_PER_READ = 60
MAX_ATTEMPTS = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--reads", type=int, default=200_000, help="read pairs per library (default 200000)")
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--force", action="store_true", help="overwrite files that are already there")
    parser.add_argument("--dry-run", action="store_true", help="say what would be fetched and stop")
    return parser.parse_args()


def ena_urls(run: str) -> list[str]:
    """The run's FASTQ URLs at ENA, in mate order."""
    query = urllib.parse.urlencode(dict(accession=run, result="read_run", fields="fastq_ftp", format="tsv"))
    with urllib.request.urlopen(f"{ENA_REPORT}?{query}", timeout=120) as response:
        rows = response.read().decode().splitlines()
    if len(rows) < 2 or not rows[1].strip():
        raise RuntimeError(f"ENA has no FASTQ for {run}")
    return [f"https://{path}" for path in rows[1].split("\t")[-1].split(";") if path]


def head_records(url: str, reads: int) -> bytes:
    """The first `reads` complete FASTQ records of a gzipped file, as uncompressed bytes.

    Only the head of the remote file is requested. A gzip stream cut mid-way still decompresses up to the cut, so
    the decompressor is fed what arrived and the trailing partial record is dropped.
    """
    wanted = reads * 4
    for attempt in range(MAX_ATTEMPTS):
        span = reads * BYTES_PER_READ * (2 ** attempt)
        request = urllib.request.Request(url, headers={"Range": f"bytes=0-{span}"})
        with urllib.request.urlopen(request, timeout=300) as response:
            raw = response.read()
        text = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw)
        lines = text.split(b"\n")
        if len(lines) > wanted:
            return b"\n".join(lines[:wanted]) + b"\n"
        if len(raw) < span:  # the whole file arrived and it is simply shorter than asked for
            return b"\n".join(lines[:len(lines) // 4 * 4]) + b"\n"
    raise RuntimeError(f"{url}: {reads} reads not reached after {MAX_ATTEMPTS} attempts")


def main() -> int:
    args = parse_args()
    root = args.repo_root.expanduser().resolve()
    config = load_project_config(args.project_config.expanduser().resolve())
    fastq_dir = Path(config["fastq_dir"])
    if not fastq_dir.is_absolute():
        fastq_dir = root / fastq_dir
    samples = [s for s in config.get("samples", []) if s.get("fastq_1")]
    if not samples:
        print("ERROR: the project has no samples with FASTQ paths", file=sys.stderr)
        return 2

    print(f"Project {config['project_id']}: {len(samples)} libraries, {args.reads} read pairs each")
    print(f"Writing under {fastq_dir}")
    written = []
    for sample in samples:
        run = str(sample["sample"])
        targets = [sample.get("fastq_1"), sample.get("fastq_2")]
        targets = [fastq_dir / t for t in targets if t]
        if all(t.exists() for t in targets) and not args.force:
            print(f"  {run}: already present, skipping")
            continue
        if args.dry_run:
            print(f"  {run}: would write {', '.join(str(t.relative_to(fastq_dir)) for t in targets)}")
            continue
        urls = ena_urls(run)
        if len(urls) < len(targets):
            print(f"  WARNING {run}: ENA lists {len(urls)} file(s) for {len(targets)} expected; skipping",
                  file=sys.stderr)
            continue
        for url, target in zip(urls, targets):
            target.parent.mkdir(parents=True, exist_ok=True)
            records = head_records(url, args.reads)
            target.write_bytes(_gzip(records))
            written.append(dict(run=run, file=str(target.relative_to(root)), source=url,
                                reads=records.count(b"\n") // 4, bytes=target.stat().st_size))
            print(f"  {run}: {target.name} ({written[-1]['reads']} reads, {written[-1]['bytes'] / 1048576:.1f} MB)")

    if written:
        note = fastq_dir / "SUBSAMPLE.json"
        note.write_text(json.dumps(dict(
            project=config["project_id"], fetched_at=datetime.now(timezone.utc).isoformat(),
            reads_per_library=args.reads, source="ENA (https://www.ebi.ac.uk/ena)",
            note="The first N read pairs of each run, not a random sample. For trying the mapping stage; "
                 "insertion counts are far below the full libraries' and do not reproduce the stored results.",
            files=written), indent=2) + "\n", encoding="utf-8")
        print(f"Provenance: {note}")
    return 0


def _gzip(data: bytes) -> bytes:
    import gzip as _g
    import io
    buffer = io.BytesIO()
    with _g.GzipFile(fileobj=buffer, mode="wb", compresslevel=6, mtime=0) as handle:
        handle.write(data)
    return buffer.getvalue()


if __name__ == "__main__":
    raise SystemExit(main())
