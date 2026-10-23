#!/usr/bin/env python3
"""Build the hybrid essentiality model's gene features from the parents, normalized with any method.

The parents are normalized with the chosen method, written as hit files, combined and summarized (SummaryTable); the
result (summary_combined/<tag>/combined_feature_table.<tag>.txt) is what `ytab essentiality --features` takes.
--with-classifier also scores it with the legacy classifier (classifier/<tag>/essentiality_predictions.<tag>.csv),
for comparison only. Every method, MidLC common saturation included, is computed from the newest Site
Normalization run whose parent catalogue holds the parents (see ytab.pipeline.hybrid_inputs for the whole-number
rules). Results already built for the same method, run and parents are reused without being touched.

  python scripts/local/ytab.py hybrid-inputs --project-config output/projects/P/config/project.yaml \\
      --normalization betageom --samples s1,s2,s3
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))

from ytab.pipeline.hybrid_inputs import catalogue_runs, normalize_parents  # noqa: E402
from ytab.pipeline.mapfastq_runner import get_included_samples, load_project_for_mapping  # noqa: E402
from ytab.pipeline.summary_table_runner import find_sample_hits_file  # noqa: E402

SCRIPTS = ROOT/"scripts/local"


def run(script: str, args: list[str], label: str) -> None:
    print(f"=== {label}", file=sys.stderr, flush=True)
    result = subprocess.run([sys.executable, str(SCRIPTS/script), *args], cwd=ROOT, capture_output=True, text=True)
    sys.stderr.write(result.stdout+result.stderr)
    if result.returncode:
        lines = [x for x in (result.stderr+result.stdout).strip().splitlines() if x.strip()]
        raise ValueError(f"{label} failed: {(lines or ['no output'])[-1].removeprefix('ERROR: ')}")


def build(project_config: Path, normalization: str, samples: list[str], normalization_run: Path | None = None,
          threads: int | None = None, force: bool = False, with_classifier: bool = False) -> dict:
    """Gene features of the parents normalized with the method: combined, summarized (SummaryTable). With
    with_classifier, the legacy classifier also scores them (its prediction table, for comparison only)."""
    config = load_project_for_mapping(project_config)
    root = Path(config["_repo_root"]); project = Path(config["output_project_dir"])
    project = project if project.is_absolute() else root/project
    included = {str(s["sample"]):s for s in get_included_samples(config)}
    unknown = [s for s in samples if s not in included]
    if unknown:raise ValueError(f"Not included samples: {', '.join(unknown)}")
    project_id = str(config.get("project_id"))
    runs = [Path(normalization_run)] if normalization_run else catalogue_runs(project, [f"{project_id}::{s}" for s in samples])
    if not runs:raise ValueError("No Site Normalization run holds these parents; run Normalization > Site Normalization first")
    hits = {s:find_sample_hits_file(included[s], config) for s in samples}
    source = normalize_parents(project, runs[0], normalization, project_id, samples, hits, force=force)
    tag = source["tag"]
    features = project/"summary_combined"/tag/f"combined_feature_table.{tag}.txt"
    predictions = project/"classifier"/tag/f"essentiality_predictions.{tag}.csv"
    common = ["--project-config", str(project_config), "--target", tag, *(["--force"] if force else [])]
    status = "reused"
    if force or not features.is_file():
        run("ytab_run_combine_hits.py", [*common, "--samples", ",".join(samples)], f"Combining {len(samples)} parents ({tag})")
        run("ytab_run_summary_combined.py", [*common, *(["--threads", str(threads)] if threads else [])], "Combined SummaryTable")
        status = "success"
    if not features.is_file():raise ValueError(f"The combined SummaryTable wrote no {features.name}")
    if with_classifier and (force or not predictions.is_file()):
        run("ytab_run_classifier.py", common, "Legacy classifier (comparison only)")
        if not predictions.is_file():raise ValueError(f"The classifier wrote no {predictions.name}")
    return dict(status=status, normalization=normalization, tag=tag, features=str(features),
                predictions=str(predictions) if with_classifier else None, source=source)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--normalization", default="betageom", help="Normalization method id (see ytab normalize --list-methods), e.g. betageom, ttr or midlc")
    p.add_argument("--normalization-run", type=Path, help="Site Normalization run to take the parent catalogue from (default: newest that has the parents)")
    p.add_argument("--samples", required=True, help="Comma-separated parent libraries (one per independent pool)")
    p.add_argument("--threads", type=int); p.add_argument("--force", action="store_true")
    p.add_argument("--with-classifier", action="store_true", help="Also score the features with the legacy classifier (comparison only)")
    a = p.parse_args()
    try:
        result = build(a.project_config, a.normalization, [s.strip() for s in a.samples.split(",") if s.strip()],
                       a.normalization_run, a.threads, a.force, a.with_classifier)
    except (ValueError, FileNotFoundError, OSError) as error:
        p.exit(1, f"ERROR: {error}\n")
    print(json.dumps(result, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
