"""Restartable, one-sample-at-a-time local MapFastq orchestration."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml
from .progress_tracker import ProgressTracker, make_job_id


STATUS_FIELDS = [
    "sample", "status", "layout", "fastq_1", "fastq_2", "output_dir",
    "log_file", "elapsed_seconds", "message",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _index_complete(prefix: Path) -> bool:
    return any(
        all(Path(f"{prefix}{suffix}.{extension}").is_file() for suffix in (
            ".1", ".2", ".3", ".4", ".rev.1", ".rev.2"
        ))
        for extension in ("bt2", "bt2l")
    )


def load_project_for_mapping(project_config: Path) -> dict:
    config_path = Path(project_config).expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Project config is not a YAML mapping: {config_path}")

    repo_root = _resolve(config.get("repo_root"), config_path.parent) or Path(__file__).resolve().parents[3]
    config["_project_config"] = str(config_path)
    config["_repo_root"] = str(repo_root)

    sample_sheet = _resolve(config.get("sample_sheet"), repo_root)
    if sample_sheet:
        if not sample_sheet.is_file():
            raise FileNotFoundError(f"Sample sheet is missing: {sample_sheet}")
        with sample_sheet.open(newline="", encoding="utf-8") as handle:
            config["samples"] = list(csv.DictReader(handle))

    reference_path = config_path.with_name("reference_resolved.json")
    if reference_path.is_file():
        with reference_path.open(encoding="utf-8") as handle:
            reference = json.load(handle)
        if not isinstance(reference, dict):
            raise ValueError(f"Reference file is not a JSON object: {reference_path}")
        config["reference"] = reference
    config["_reference_resolved"] = str(reference_path)
    return config


def get_included_samples(config: dict) -> list[dict]:
    return [dict(sample) for sample in config.get("samples", []) if _as_bool(sample.get("include", True))]


def validate_mapping_inputs(config: dict) -> list[str]:
    errors: list[str] = []
    root = Path(config["_repo_root"])
    try:
        threads = int(config.get("threads", 0))
    except (TypeError, ValueError):
        threads = 0
    if threads < 2:
        errors.append("Threads must be at least 2.")

    reference = config.get("reference") or {}
    prefix = _resolve(reference.get("bowtie2_index_prefix"), root)
    if prefix is None or not _index_complete(prefix):
        errors.append("No complete Bowtie2 index found. Run reference preparation before alignment.")

    samples = get_included_samples(config)
    if not samples:
        errors.append("No included samples were found in the sample sheet.")
    for sample in samples:
        settings = {**(config.get("mapping") or {}), **{k:v for k,v in sample.items() if v not in (None, "")}}
        name = str(sample.get("sample") or "")
        if not name:
            errors.append("An included sample has no sample name.")
        elif Path(name).name != name:
            errors.append(f"Sample name contains path separators: {name}")
        layout = str(sample.get("layout") or "single").lower()
        fastq_1 = _resolve(sample.get("fastq_1"), root)
        fastq_2 = _resolve(sample.get("fastq_2"), root)
        mate=str(settings.get("junction_mate") or settings.get("genome_mapping_mate") or "R1").upper()
        paired=layout=="paired" or str(settings.get("sequencing_layout") or settings.get("read_layout") or "").lower() in {"paired","paired_end"}
        agreement=settings.get("require_mate_alignment")
        agreement=str(settings.get("mapping_strategy") or "").lower()=="paired_genome" if agreement in (None, "") else _as_bool(agreement)
        if (not paired or agreement or mate=="R1") and (fastq_1 is None or not fastq_1.is_file()):
            errors.append(f"FASTQ file missing for sample {name}: {fastq_1 or 'fastq_1 not set'}")
        if paired and (agreement or mate=="R2") and (fastq_2 is None or not fastq_2.is_file()):
            errors.append(f"Required FASTQ R2 missing for sample {name}: {fastq_2 or 'fastq_2 not set'}")
        if layout not in {"single", "paired"}:
            errors.append(f"Unsupported layout for sample {name}: {layout}")
    return errors


def _paths(config: dict, sample_name: str) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project_dir = _resolve(config.get("output_project_dir"), root)
    export_dir = _resolve(config.get("output_export_dir"), root)
    if project_dir is None or export_dir is None:
        raise ValueError("Project and export output directories must be configured.")
    return {
        "output": project_dir / "mapfastq" / sample_name,
        "log": project_dir / "logs" / "mapfastq" / f"{sample_name}.log",
        "manifest": project_dir / "manifests" / "mapfastq" / f"{sample_name}.mapfastq_manifest.json",
        "status": project_dir / "manifests" / "mapfastq" / "mapfastq_status.csv",
        "qc": export_dir / "qc" / "mapping",
    }


def build_mapfastq_command(sample: dict, config: dict, threads: int | None = None) -> list[str]:
    sample = {**(config.get("mapping") or {}), **{k:v for k,v in sample.items() if v not in (None, "")}}
    root = Path(config["_repo_root"])
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 2:
        raise ValueError("Threads must be at least 2.")

    name = str(sample["sample"])
    layout = str(sample.get("layout") or "single").lower()

    # Separate sequencing layout from genome-mapping strategy.
    # Paired-end input means both FASTQs exist and are tracked.
    # It does not automatically mean Bowtie2 paired-end genome alignment.
    sequencing_layout = str(
        sample.get("sequencing_layout")
        or sample.get("read_layout")
        or ("paired_end" if layout in {"paired", "paired_end"} else "single_end")
    ).lower()

    raw_mapping_strategy = sample.get("mapping_strategy") or "single_mate"
    raw_genome_mapping_mate = (
        sample.get("junction_mate") or sample.get("genome_mapping_mate")
        or sample.get("resolved_genome_mapping_mate")
        or "R1"
    )

    if sequencing_layout in {"paired", "paired_end"}:
        mapping_strategy = str(raw_mapping_strategy).lower()
        if mapping_strategy == "auto":
            mapping_strategy = "single_mate"

        if mapping_strategy not in {"single_mate", "paired_genome"}:
            raise ValueError(
                f"Unsupported mapping_strategy for {name}: {mapping_strategy}. "
                "Use single_mate, paired_genome, or auto."
            )

        genome_mapping_mate = str(raw_genome_mapping_mate).upper()

    else:
        mapping_strategy = str(raw_mapping_strategy or "single_mate").lower()
        genome_mapping_mate = str(raw_genome_mapping_mate or "R1").upper()

    if genome_mapping_mate == "AUTO":
        raise ValueError("Automatic mate guessing is not used: declare junction_mate R1 or R2 (default R1)")
    if genome_mapping_mate not in {"R1","R2"}:
        raise ValueError("Select one junction_mate: R1 or R2, never both insertion endpoints")
    require = sample.get("require_mate_alignment")
    require = mapping_strategy == "paired_genome" if require in (None, "") else _as_bool(require)
    if require and sequencing_layout not in {"paired","paired_end"}:
        raise ValueError("Mate agreement requires paired-end inputs")

    prefix = _resolve((config.get("reference") or {}).get("bowtie2_index_prefix"), root)
    output = _paths(config, name)["output"]
    script = root / "src" / "ytab" / "mapping" / "MapFastq.py"

    command = [
        sys.executable, str(script), "--out-dir", str(output), "--bt2-index",
        str(prefix), "--threads", str(selected_threads), "--sample-name", name,
    ]
    if _as_bool(sample.get("hermes_cleanup",True)):
        command.append("--hermes-preprocess")
        for field, flag in (("primer_catalogue","--primer-catalogue"),("hermes_reference","--hermes-reference")):
            if sample.get(field):
                command.extend([flag,str(_resolve(sample[field],root))])
        for field, flag in (("hermes_overlap","--hermes-overlap"),("adapter_error_rate","--adapter-error-rate"),("min_read_length","--min-read-length"),("junction_read_end","--junction-read-end")):
            if sample.get(field) not in (None, ""):
                command.extend([flag,str(sample[field])])
    for field,flag in (("mapq","--mapq"),("drop_flags","--drop-flags")):
        if sample.get(field) not in (None, ""):
            command.extend([flag,str(sample[field])])

    fastq_1 = _resolve(sample.get("fastq_1") or sample.get("r1_fastq"), root)
    fastq_2 = _resolve(sample.get("fastq_2") or sample.get("r2_fastq"), root)

    if require:
        if fastq_2 is None:
            raise ValueError(f"paired_genome mapping requested for {name}, but fastq_2 is missing.")
        command.extend(["--r1", str(fastq_1), "--r2", str(fastq_2),"--junction-mate",genome_mapping_mate,"--require-mate-alignment"])
        return command

    if sequencing_layout in {"paired", "paired_end"}:
        if genome_mapping_mate == "R2":
            if fastq_2 is None:
                raise ValueError(f"R2 genome mapping requested for {name}, but fastq_2 is missing.")
            map_fastq = fastq_2
        else:
            map_fastq = fastq_1
    else:
        map_fastq = fastq_1

    command.extend(["--input-file-name", str(map_fastq),"--source-junction-mate",genome_mapping_mate])
    return command


def mapping_fingerprint(command: list[str]) -> str:
    """Invalidate old outputs when reads, settings, index, code or primer assets change."""
    root=Path(__file__).resolve().parents[3]
    payload={"schema":2,"command":command,"files":{}}
    payload["tools"]={}
    for tool in ("bowtie2","samtools"):
        binary=shutil.which(tool)
        payload["tools"][tool]=[binary,Path(binary).stat().st_size,Path(binary).stat().st_mtime_ns] if binary else None
    for flag in ("--input-file-name","--r1","--r2"):
        if flag in command:
            path=Path(command[command.index(flag)+1])
            payload["files"][str(path)]=[path.stat().st_size,path.stat().st_mtime_ns] if path.is_file() else None
    prefix=command[command.index("--bt2-index")+1]
    for path in Path(prefix).parent.glob(Path(prefix).name+"*.bt2*"):
        payload["files"][str(path)]=[path.stat().st_size,path.stat().st_mtime_ns]
    files=[Path(__file__),root/"src/ytab/mapping/MapFastq.py"]
    if "--hermes-preprocess" in command:
        import cutadapt
        payload["tools"]["cutadapt_version"]=cutadapt.__version__
        from ytab.mapping.hermes_preprocess import DEFAULT_CATALOGUE
        catalogue=Path(command[command.index("--primer-catalogue")+1]) if "--primer-catalogue" in command else DEFAULT_CATALOGUE
        data=json.loads(catalogue.read_text())
        reference=Path(command[command.index("--hermes-reference")+1]) if "--hermes-reference" in command else catalogue.parent/data["reference"]
        files += [catalogue,reference,root/"src/ytab/mapping/hermes_preprocess.py"]
    for path in files:
        payload["files"][str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()

def _detected_outputs(output_dir: Path) -> list[str]:
    return sorted(str(path) for path in output_dir.iterdir() if path.is_file()) if output_dir.is_dir() else []


def _expected_outputs(output_dir: Path, sample: str) -> bool:
    bam = output_dir / f"{sample}.sorted.bam"
    return bam.is_file() and (Path(f"{bam}.bai").is_file() or bam.with_suffix(".bai").is_file())


def _write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_mapfastq_sample(
    sample: dict, config: dict, threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    root = Path(config["_repo_root"])
    name = str(sample["sample"])
    paths = _paths(config, name)
    selected_threads = int(config.get("threads") if threads is None else threads)
    command = build_mapfastq_command(sample, config, selected_threads)
    fingerprint = mapping_fingerprint(command)
    fastq_1 = _resolve(sample.get("fastq_1"), root)
    fastq_2 = _resolve(sample.get("fastq_2"), root)
    prefix = _resolve((config.get("reference") or {}).get("bowtie2_index_prefix"), root)
    start_time = _now()
    started = time.monotonic()
    warnings = [str(sample.get("warnings"))] if sample.get("warnings") else []
    manifest = {
        "project_id": config.get("project_id"), "sample": name,
        "fastq_1": str(fastq_1) if fastq_1 else None,
        "fastq_2": str(fastq_2) if fastq_2 else None,
        "layout": str(sample.get("layout") or "single").lower(),
        "sequencing_layout": str(
            sample.get("sequencing_layout")
            or sample.get("read_layout")
            or ("paired_end" if str(sample.get("layout") or "single").lower() in {"paired", "paired_end"} else "single_end")
        ).lower(),
        "mapping_strategy": str(sample.get("mapping_strategy") or "").lower(),
        "genome_mapping_mate": str(sample.get("genome_mapping_mate") or sample.get("resolved_genome_mapping_mate") or "").upper(),
        "junction_mate": str(sample.get("junction_mate") or "").upper(),
        "require_mate_alignment": str(sample.get("require_mate_alignment") or "false").lower(),
        "species": config.get("species"), "bowtie2_index_prefix": str(prefix),
        "threads": selected_threads, "output_dir": str(paths["output"]),
        "log_file": str(paths["log"]), "status": "failed",
        "start_time": start_time, "end_time": None, "elapsed_seconds": 0.0,
        "command_run": command, "detected_outputs": [], "warnings": warnings,
        "request_fingerprint": fingerprint,
        "error_message": None,
    }
    mate_flag="--junction-mate" if "--junction-mate" in command else "--source-junction-mate"
    manifest.update(mapping_strategy="paired_genome" if "--require-mate-alignment" in command else "single_mate",
        junction_mate=command[command.index(mate_flag)+1],genome_mapping_mate=command[command.index(mate_flag)+1],
        require_mate_alignment="--require-mate-alignment" in command,hermes_cleanup="--hermes-preprocess" in command)

    if not force and paths["manifest"].is_file() and _expected_outputs(paths["output"], name):
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        if previous.get("status") in {"success", "skipped"} and previous.get("request_fingerprint")==fingerprint:
            manifest.update(previous)
            manifest.update(status="skipped", end_time=_now())
            manifest["warnings"].append("Successful cached outputs found; mapping skipped.")
            manifest["detected_outputs"] = _detected_outputs(paths["output"])
            _write_manifest(paths["manifest"], manifest)
            return manifest

    if dry_run:
        manifest.update(status="dry-run", end_time=_now())
        manifest["warnings"].append("Dry run; command was not executed.")
        _write_manifest(paths["manifest"], manifest)
        return manifest

    paths["output"].mkdir(parents=True, exist_ok=True)
    paths["log"].parent.mkdir(parents=True, exist_ok=True)
    try:
        with paths["log"].open("w", encoding="utf-8") as log:
            log.write("Command: " + " ".join(command) + "\n\n")
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
        if completed.returncode != 0:
            raise RuntimeError(f"MapFastq exited with status {completed.returncode}; see {paths['log']}")
        if not _expected_outputs(paths["output"], name):
            raise RuntimeError("MapFastq completed without the expected sorted BAM and index.")
        stats = paths["output"] / f"{name}.mapping_stats.csv"
        if stats.is_file():
            paths["qc"].mkdir(parents=True, exist_ok=True)
            shutil.copy2(stats, paths["qc"] / stats.name)
        manifest["status"] = "success"
        preprocessing=paths["output"]/f"{name}.hermes_preprocessing.json"
        if "--hermes-preprocess" in command:
            if not preprocessing.is_file():
                raise RuntimeError("Hermes preprocessing report is missing")
            manifest["hermes_preprocessing"]=json.loads(preprocessing.read_text())
            manifest["genome_mapping_mate"]=manifest["hermes_preprocessing"]["junction_mate"]
            manifest["junction_mate"]=manifest["genome_mapping_mate"]
            manifest["require_mate_alignment"]=manifest["hermes_preprocessing"]["require_mate_alignment"]
        manifest["status"] = "success"
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error_message"] = str(exc)
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    manifest["detected_outputs"] = _detected_outputs(paths["output"])
    _write_manifest(paths["manifest"], manifest)
    return manifest


def _write_status(path: Path, results: list[dict]) -> None:
    existing: dict[str, dict] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = {row["sample"]: row for row in csv.DictReader(handle)}
    for result in results:
        existing[result["sample"]] = {
            "sample": result["sample"], "status": result["status"],
            "layout": result["layout"], "fastq_1": result["fastq_1"] or "",
            "fastq_2": result["fastq_2"] or "", "output_dir": result["output_dir"],
            "log_file": result["log_file"], "elapsed_seconds": result["elapsed_seconds"],
            "message": result.get("error_message") or "; ".join(result.get("warnings") or []),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[name] for name in sorted(existing))


def run_mapfastq_project(
    project_config: Path, samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
    keep_going: bool = False, job_id: str | None = None,
    progress_file: Path | None = None,
    junction_mate: str | None = None, require_mate_alignment: bool | None = None,
    hermes_cleanup: bool | None = None,
) -> dict:
    config = load_project_for_mapping(project_config)
    if threads is not None:
        config["threads"] = threads
    for sample in config.get("samples",[]):
        for key,value in (("junction_mate",junction_mate),("require_mate_alignment",require_mate_alignment),("hermes_cleanup",hermes_cleanup)):
            if value is not None:
                sample[key]=value
    errors = validate_mapping_inputs(config)
    if errors:
        raise ValueError("\n".join(errors))
    included = get_included_samples(config)
    by_name = {str(sample["sample"]): sample for sample in included}
    selected_names = list(by_name) if samples is None else samples
    unknown = [name for name in selected_names if name not in by_name]
    if unknown:
        raise ValueError(f"Selected samples are not included or do not exist: {', '.join(unknown)}")
    project_dir = _resolve(config.get("output_project_dir"), Path(config["_repo_root"]))
    job_id = job_id or make_job_id("mapfastq")
    progress_file = Path(progress_file) if progress_file else project_dir / "manifests" / "jobs" / f"{job_id}.progress.json"
    tracker = ProgressTracker.create(progress_file, job_id, str(config.get("project_id")), "mapfastq", selected_names,
        input_files={name:str(_resolve(by_name[name].get("fastq_1"),Path(config["_repo_root"]))) for name in selected_names})
    tracker.state["dry_run"]=dry_run
    tracker.start("Mapping started")
    results: list[dict] = []
    for index, name in enumerate(selected_names, 1):
        tracker.start_item(name,index,"cache check")
        paths=_paths(config,name)
        running={"sample":name,"status":"running","layout":str(by_name[name].get("layout") or "single"),"fastq_1":str(_resolve(by_name[name].get("fastq_1"),Path(config["_repo_root"]))) if by_name[name].get("fastq_1") else None,"fastq_2":str(_resolve(by_name[name].get("fastq_2"),Path(config["_repo_root"]))) if by_name[name].get("fastq_2") else None,"output_dir":str(paths["output"]),"log_file":str(paths["log"]),"elapsed_seconds":0.0,"warnings":[],"error_message":"Mapping command starting."}
        _write_status(paths["status"],results+[running]);tracker.phase("mapping / BAM processing","MapFastq command running")
        result = run_mapfastq_sample(by_name[name], config, threads, force, dry_run)
        results.append(result)
        _write_status(_paths(config, name)["status"], results)
        message=result.get("error_message") or "; ".join(result.get("warnings") or [])
        tracker.finish_item(name,result["status"],result.get("detected_outputs"),message,result.get("error_message") or "")
        if result["status"] == "failed" and not keep_going:
            break
    final_status="success" if not any(x["status"]=="failed" for x in results) else ("partial" if any(x["status"] in {"success","skipped"} for x in results) else "failed")
    tracker.finalize("dry_run_success" if dry_run else final_status,"Commands were validated but alignment was not executed." if dry_run else "Mapping complete")
    return {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "bowtie2_index_prefix": str(_resolve((config.get("reference") or {}).get("bowtie2_index_prefix"), Path(config["_repo_root"]))),
        "threads": int(config["threads"]), "included_samples": len(included),
        "selected_samples": selected_names, "results": results,
        "success": sum(row["status"] == "success" for row in results),
        "failed": sum(row["status"] == "failed" for row in results),
        "skipped": sum(row["status"] == "skipped" for row in results),
        "dry_run": dry_run,
        "job_id":job_id,"progress_file":str(progress_file),
    }
