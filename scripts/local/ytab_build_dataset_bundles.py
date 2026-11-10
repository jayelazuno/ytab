#!/usr/bin/env python3
"""Package each processed dataset as one archive, ready to publish and to fetch back into a YTAB installation.

A dataset's archive holds its project directory and its exports, minus what a machine should not carry to another:
logs, caches, job state and the generated export archives. Legacy normalization targets are left out as well, so a
published dataset contains only runs of the current pipeline.

The archives are gzipped tar files; the tables inside are repetitive text and compress about six-fold, so a dataset
downloads at roughly a sixth of its size on disk. Each archive is written with a sha256 and a row in
`bundles_index.tsv`, which carries the fields `resources/datasets/catalog.tsv` needs to point the app at them:
`archive_name`, `archive_bytes`, `archive_sha256`.

    python scripts/local/ytab_build_dataset_bundles.py --out-dir output/zenodo
    python scripts/local/ytab_build_dataset_bundles.py --out-dir output/zenodo --dataset Michel_2017_SATAY

Datasets listed in `resources/datasets/catalog.tsv` are published; anything else in `output/projects` is unpublished
work and is skipped unless named explicitly with --dataset.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Machine-specific or regenerable, so never published. The genome-browser tracks are not in this list: the browser
# reads the bedGraph, WIG and bigWig files directly, so a dataset that arrived without them would open without its
# tracks.
SKIP_DIRS = {"logs", "cache", "tmp", "bundles", "jobs", "argon", "__pycache__"}
SKIP_SUFFIXES = {".pyc", ".lock"}
# The superseded MidLC-value targets: a directory named T<number>, as against the current <method>__<run>__<hash>.
LEGACY_TARGET = re.compile(r"^T\d[\dA-Za-z]*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=Path("output/zenodo"))
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--dataset", action="append", help="only this dataset (repeatable)")
    parser.add_argument("--exclude", action="append", default=[], help="skip this dataset (repeatable)")
    parser.add_argument("--dry-run", action="store_true", help="report what each archive would hold and stop")
    parser.add_argument("--index-only", action="store_true",
                        help="rewrite the index from the archives already in --out-dir, building nothing")
    return parser.parse_args()


def published_projects(root: Path) -> list[str]:
    """The project folders the dataset catalogue lists, in catalogue order."""
    path = root / "resources" / "datasets" / "catalog.tsv"
    names: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            names += [p for p in (row.get("project") or "").split() if p and p not in names]
    return names


def wanted(path: Path) -> bool:
    parts = set(path.parts)
    if parts & SKIP_DIRS or path.suffix in SKIP_SUFFIXES:
        return False
    return not any(LEGACY_TARGET.match(part) for part in path.parts)


def members(root: Path, project: str) -> list[tuple[Path, str]]:
    """(file on disk, path inside the archive) for everything the dataset should carry."""
    found: list[tuple[Path, str]] = []
    for base in (Path("output/projects") / project, Path("output/exports") / project):
        directory = root / base
        if not directory.is_dir():
            continue
        for item in sorted(directory.rglob("*")):
            if item.is_file() and wanted(item.relative_to(directory)):
                found.append((item, str(base / item.relative_to(directory))))
    return found


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    root = args.repo_root.expanduser().resolve()
    out_dir = (root / args.out_dir) if not args.out_dir.is_absolute() else args.out_dir
    projects = args.dataset or [p for p in published_projects(root) if p not in set(args.exclude)]
    if not projects:
        print("ERROR: no datasets to package", file=sys.stderr)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for project in ([] if args.index_only else projects):
        files = members(root, project)
        if not files:
            print(f"  {project}: nothing to package, skipping", file=sys.stderr)
            continue
        plain = sum(f.stat().st_size for f, _ in files)
        if args.dry_run:
            print(f"  {project:34} {len(files):5} files  {plain / 1048576:8.0f} MB on disk")
            continue
        archive = out_dir / f"ytab_{project}.tar.gz"
        started = time.time()
        with tarfile.open(archive, "w:gz", compresslevel=6) as tar:
            for source, name in files:
                tar.add(source, arcname=name, recursive=False)
        size = archive.stat().st_size
        row = dict(project=project, archive_name=archive.name, archive_bytes=size,
                   archive_sha256=sha256(archive), files=len(files), bytes_on_disk=plain)
        rows.append(row)
        (archive.parent / f"{archive.name}.sha256").write_text(f"{row['archive_sha256']}  {archive.name}\n",
                                                               encoding="utf-8")
        print(f"  {project:34} {len(files):5} files  {plain / 1048576:8.0f} MB -> "
              f"{size / 1048576:7.0f} MB  ({plain / size:.1f}x, {time.time() - started:.0f}s)")

    # An archive that this run did not rebuild is still published, so the index describes everything in the
    # directory; rebuilding one dataset must not drop the other seventeen rows.
    built = {r["project"] for r in rows}
    for archive in sorted(out_dir.glob("ytab_*.tar.gz")):
        project = archive.name[len("ytab_"):-len(".tar.gz")]
        if project in built:
            continue
        files = members(root, project)
        rows.append(dict(project=project, archive_name=archive.name, archive_bytes=archive.stat().st_size,
                         archive_sha256=sha256(archive), files=len(files),
                         bytes_on_disk=sum(f.stat().st_size for f, _ in files)))
    rows.sort(key=lambda r: r["project"])

    if rows:
        index = out_dir / "bundles_index.tsv"
        with index.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        (out_dir / "bundles_index.json").write_text(json.dumps(dict(
            built_at=datetime.now(timezone.utc).isoformat(), datasets=rows), indent=2) + "\n", encoding="utf-8")
        total = sum(r["archive_bytes"] for r in rows)
        print(f"\n{len(rows)} archives, {total / 1073741824:.2f} GB total")
        print(f"Index: {index}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
