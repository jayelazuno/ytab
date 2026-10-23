#!/usr/bin/env python3
"""Map reads, create hit files and build SummaryTables for selected samples in one run (the app's Run > Mapping).

Each stage uses its own cache-aware runner (ytab_run_mapfastq.py, ytab_run_create_hit_file.py,
ytab_run_summary_table.py), so finished samples are skipped. A sample whose later outputs already exist is not sent
to the earlier stages: a hit file means it was mapped (e.g. on the HPC, with the BAM left there), and a SummaryTable
means its hit file was made; a sample with a SummaryTable is not run at all. A sample that fails a stage does not
continue to the next one.

  python scripts/local/ytab.py preprocess --project-config output/projects/P/config/project.yaml --samples s1,s2
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))

from ytab.pipeline.create_hit_file_runner import validate_create_hit_file_inputs  # noqa: E402
from ytab.pipeline.mapfastq_runner import get_included_samples, load_project_for_mapping, validate_mapping_inputs  # noqa: E402
from ytab.pipeline.summary_table_runner import validate_summary_table_inputs  # noqa: E402
from ytab.pipeline.progress_tracker import ProgressTracker, load_progress_state, make_job_id  # noqa: E402

STAGES = (("mapfastq", "Mapping", "ytab_run_mapfastq.py"), ("create_hit_file", "Hit files", "ytab_run_create_hit_file.py"),
          ("summary", "SummaryTable", "ytab_run_summary_table.py"))
DONE = {"success", "skipped", "cached", "dry_run_success", "dry-run", "dry_run"}
VALIDATE = {"mapfastq":validate_mapping_inputs, "create_hit_file":validate_create_hit_file_inputs, "summary":validate_summary_table_inputs}


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def resolve(value, root: Path) -> Path | None:
    if not value:return None
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else root/path


def has_summary(project: Path, name: str) -> bool:
    folder = project/"summary"/name
    return folder.is_dir() and any(nonempty(p) for pattern in ("*.feature_table*.csv", "*.feature_table*.tsv", "*.feature_table*.txt")
                                   for p in folder.glob(pattern))


def has_hits(project: Path, sample: dict, root: Path) -> bool:
    name = str(sample["sample"])
    imported = resolve(sample.get("hit_file"), root)
    return nonempty(project/"create_hit_file"/name/f"{name}_hits.txt") or bool(imported and nonempty(imported))


def has_mapping(project: Path, name: str) -> bool:
    folder = project/"mapfastq"/name; bam = folder/f"{name}.sorted.bam"
    stats = nonempty(folder/f"{name}.mapping_stats.csv")
    indexed = nonempty(bam) and (nonempty(Path(f"{bam}.bai")) or nonempty(bam.with_suffix(".bai")))
    return stats and (indexed or nonempty(folder/"scratch_bam_path.txt"))


def plan(config: dict, selected: list[str], force: bool) -> dict[str, list[str]]:
    """Samples each stage receives before any stage runs; later stages also drop samples that fail earlier ones.
    Finished samples are never sent to a runner (without --force): a runner's cache check rewrites its manifest,
    which would replace the record of the original run."""
    root = Path(config["_repo_root"]); project = resolve(config.get("output_project_dir"), root)
    rows = {str(s["sample"]):s for s in get_included_samples(config)}
    summary = {n:has_summary(project, n) for n in selected}; hits = {n:summary[n] or has_hits(project, rows[n], root) for n in selected}
    mapped = {n:hits[n] or has_mapping(project, n) for n in selected}
    imported = str(config.get("start_stage") or "fastq").strip().lower() == "create_hit_file"
    return {"mapfastq":[] if imported else [n for n in selected if force or not mapped[n]],
            "create_hit_file":[n for n in selected if force or not hits[n]],
            "summary":[n for n in selected if force or not summary[n]]}


def input_errors(config: dict, stage: str, samples: list[str]) -> dict[str, str]:
    """Per-sample input problems for a stage (a runner checks all its samples together and stops on the first)."""
    shared = set(VALIDATE[stage](config, []))-{"No included samples were found in the sample sheet."}
    if shared:return {n:sorted(shared)[0] for n in samples}
    return {n:errors[0] for n in samples if (errors := VALIDATE[stage](config, [n]))}


def run_stage(script: str, args: list[str], child_progress: Path, tracker: ProgressTracker, label: str, poll: float = 2.0) -> dict:
    """Run one stage script, mirroring its current sample and phase into this job's progress record."""
    command = [sys.executable, str(ROOT/"scripts/local"/script), *args, "--progress-file", str(child_progress)]
    print(f"\n=== {label}: {' '.join(command[1:])}", flush=True)
    process = subprocess.Popen(command, cwd=ROOT)
    while process.poll() is None:
        state = load_progress_state(child_progress) or {}
        item = state.get("current_item")
        tracker.state.update(current_item=item, current_phase=f"{label}: {state.get('current_phase') or 'starting'}",
                             message=f"{label}: {state.get('processed_items', 0)}/{state.get('total_items', 0)} samples"+(f" (now {item})" if item else ""))
        tracker.heartbeat(); time.sleep(poll)
    state = load_progress_state(child_progress) or {}
    return {"returncode":process.returncode, "items":{x["item"]:x for x in state.get("items", [])}}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--samples", help="Comma-separated included sample names (default: all included)")
    p.add_argument("--threads", type=int)
    p.add_argument("--force", action="store_true", help="Rerun every stage for the selected samples")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--keep-going", action="store_true")
    p.add_argument("--job-id"); p.add_argument("--progress-file", type=Path)
    p.add_argument("--junction-mate", choices=("R1", "R2"))
    p.add_argument("--require-mate-alignment", action=argparse.BooleanOptionalAction, default=None)
    a = p.parse_args()
    config = load_project_for_mapping(a.project_config)
    if a.threads:config["threads"] = a.threads
    for sample in config.get("samples", []):  # the mapping options apply to the input check as they do to the run
        for key, value in (("junction_mate", a.junction_mate), ("require_mate_alignment", a.require_mate_alignment)):
            if value is not None:sample[key] = value
    included = [str(s["sample"]) for s in get_included_samples(config)]
    selected = [x.strip() for x in a.samples.split(",") if x.strip()] if a.samples else included
    unknown = sorted(set(selected)-set(included))
    if unknown:p.exit(2, f"ERROR: not included samples: {', '.join(unknown)}\n")
    project = resolve(config.get("output_project_dir"), Path(config["_repo_root"]))
    job_id = a.job_id or make_job_id("preprocessing")
    progress = a.progress_file or project/"manifests"/"jobs"/f"{job_id}.progress.json"
    tracker = ProgressTracker.create(progress, job_id, str(config.get("project_id")), "preprocessing", selected)
    tracker.start("Mapping, hit files and SummaryTable")
    stages = plan(config, selected, a.force)
    failed: dict[str, str] = {}
    common = ["--project-config", str(a.project_config), *(["--threads", str(a.threads)] if a.threads else []),
              *(["--force"] if a.force else []), *(["--dry-run"] if a.dry_run else []), *(["--keep-going"] if a.keep_going else [])]
    mapping = [*(["--junction-mate", a.junction_mate] if a.junction_mate else []),
               *([] if a.require_mate_alignment is None else ["--require-mate-alignment" if a.require_mate_alignment else "--no-require-mate-alignment"])]
    dry_ran: set[str] = set()  # in a dry run, later stages have no inputs for samples an earlier stage only previewed
    for stage, label, script in STAGES:
        samples = [n for n in stages[stage] if n not in failed and n not in dry_ran]
        skipped = [n for n in selected if n not in stages[stage] and n not in failed]
        if skipped:print(f"\n=== {label}: already done for {', '.join(skipped)}", flush=True)
        for n, error in input_errors(config, stage, samples).items():
            failed[n] = f"{label}: {error}"; print(f"ERROR: {error}", file=sys.stderr, flush=True)
        samples = [n for n in samples if n not in failed]
        if not samples or (failed and not a.keep_going):
            if failed and not a.keep_going:break
            continue
        result = run_stage(script, [*common, "--samples", ",".join(samples), *(mapping if stage == "mapfastq" else [])],
                           progress.with_name(f"{job_id}.{stage}.progress.json"), tracker, label)
        for n in samples:
            status = (result["items"].get(n) or {}).get("status", "failed" if result["returncode"] else "success")
            if status not in DONE:failed[n] = f"{label} {status}"
        if a.dry_run:dry_ran.update(samples)
        if failed and not a.keep_going:break
    for index, n in enumerate(selected):
        ran = any(n in stages[s] for s, _, _ in STAGES)
        tracker.start_item(n, index, "recording")
        tracker.finish_item(n, "failed" if n in failed else "success" if ran else "skipped",
                            message=failed.get(n, ("Previewed (dry run)" if a.dry_run else "Mapping, hit file and SummaryTable done") if ran else "Already done"), error_message=failed.get(n, ""))
    status = "dry_run_success" if a.dry_run and not failed else ("success" if not failed else "partial" if len(failed) < len(selected) else "failed")
    tracker.finalize(status, "Preprocessing complete" if not failed else f"{len(failed)} sample(s) failed: "+"; ".join(f"{k} ({v})" for k, v in failed.items()))
    print(json.dumps({"status":status, "samples":len(selected), "failed":failed, "stages":stages}, indent=1))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
