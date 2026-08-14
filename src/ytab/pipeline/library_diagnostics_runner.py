"""Restartable orchestration and sample-set-aware caching for LibraryDiagnostics."""
from __future__ import annotations

import csv, hashlib, json, shutil, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

from .mapfastq_runner import get_included_samples as _get_included_samples
from .mapfastq_runner import load_project_for_mapping
from .progress_tracker import ProgressTracker, make_job_id

STATUS_FIELDS = ["project_id","run_id","sample_set_label","cache_signature","status","samples","input_hit_files_count","output_dir","export_dir","log_file","elapsed_seconds","message"]

def _now() -> str: return datetime.now(timezone.utc).isoformat()
def _resolve(value, root):
    if value in (None, ""): return None
    path=Path(value).expanduser(); return path.resolve() if path.is_absolute() else (root/path).resolve()
def sha256_file(path: Path) -> str:
    digest=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""): digest.update(block)
    return digest.hexdigest()
def canonical_sha256(value) -> str:
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()

def load_project_for_library_diagnostics(project_config: Path) -> dict: return load_project_for_mapping(project_config)
def get_included_samples(config: dict) -> list[dict]: return _get_included_samples(config)
def _base_paths(config):
    root=Path(config["_repo_root"]); project=_resolve(config.get("output_project_dir"),root); export=_resolve(config.get("output_export_dir"),root)
    if project is None or export is None: raise ValueError("Project and export output directories must be configured.")
    return {"root":root,"project":project,"stable_output":project/"library_diagnostics","export_root":export/"qc"/"diagnostics"}
def find_sample_hits_file(sample,config):
    name=str(sample["sample"]); return _base_paths(config)["project"]/"create_hit_file"/name/f"{name}_hits.txt"
def _select_samples(config,samples):
    included=get_included_samples(config); by_name={str(x.get("sample") or ""):x for x in included}; names=list(by_name) if samples is None else list(samples)
    if not names: raise ValueError("No included samples were selected for library diagnostics.")
    if len(names)!=len(set(names)): raise ValueError("Selected sample names must be unique.")
    unknown=[x for x in names if x not in by_name]
    if unknown: raise ValueError("Selected samples are not included or do not exist: "+", ".join(unknown))
    return [by_name[x] for x in names]
def _reference_files(config):
    root=Path(config["_repo_root"]); ref=config.get("reference") or {}
    return {"fasta":_resolve(ref.get("fasta"),root),"gff":_resolve(ref.get("gff"),root),"gtf":_resolve(ref.get("gtf"),root),"feature_table":_resolve(ref.get("feature_table"),root),"centromere_bed":_resolve(ref.get("centromere_bed") or ref.get("centromeres_bed"),root),"trna_bed":_resolve(ref.get("trna_bed"),root)}
def validate_library_diagnostics_inputs(config,samples=None):
    try: selected=_select_samples(config,samples)
    except ValueError as exc: return [str(exc)]
    errors=[]
    for sample in selected:
        name=str(sample.get("sample") or ""); hits=find_sample_hits_file(sample,config)
        if not name: errors.append("An included sample has no sample name.")
        elif Path(name).name!=name: errors.append(f"Sample name contains path separators: {name}")
        if not hits.is_file() or hits.stat().st_size<=0: errors.append(f"Valid non-empty hit file missing for sample {name}: {hits}")
    return errors
def _resource_warnings(config):
    ref=_reference_files(config); pairs=(("fasta","Reference FASTA is missing; sequence-bias diagnostics will be skipped."),("centromere_bed","Centromere BED is missing; centromere-bias diagnostics will be skipped."),("trna_bed","tRNA BED is missing; no separate tRNA BED diagnostics will be available."),("gff","Reference GFF is missing; TSS/TTS/tRNA metaplots will be skipped."))
    return [message for key,message in pairs if not ref[key] or not ref[key].is_file()]
def sample_set_label(config,names):
    selected=set(names); included=get_included_samples(config); all_names={str(x["sample"]) for x in included}; conditions={str(x["sample"]):str(x.get("guessed_condition") or "").lower() for x in included}
    if selected==all_names: return "All eligible samples"
    if selected and all(conditions.get(x)=="parent" for x in selected): return "Parents"
    if selected and all(conditions.get(x)=="treated" for x in selected): return "Treated"
    return "Custom selection"
def build_cache_context(config,samples=None,threads=None,scientific_parameters=None):
    selected=_select_samples(config,samples); names=sorted(str(x["sample"]) for x in selected); hits=[find_sample_hits_file(next(x for x in selected if str(x["sample"])==name),config) for name in names]
    refs={k:v for k,v in _reference_files(config).items() if v and v.is_file()}; script=Path(config["_repo_root"])/"src"/"ytab"/"qc"/"LibraryDiagnostics.py"
    parameters={"threads":threads,"normalization_flags":[],"scientific_cli_options":scientific_parameters or {}}
    payload={"project_id":str(config.get("project_id")),"species":str(config.get("species")),"selected_samples_sorted":names,"selected_samples_sha256":canonical_sha256(names),"input_hit_files":[str(x) for x in hits],"input_hit_file_hashes":{str(x):sha256_file(x) for x in hits},"reference_hashes":{k:sha256_file(v) for k,v in sorted(refs.items())},"library_diagnostics_script_hash":sha256_file(script),"scientific_parameters":parameters,"run_mode":"library_diagnostics"}
    payload["input_files_sha256"]=canonical_sha256({"hits":payload["input_hit_file_hashes"],"references":payload["reference_hashes"]}); payload["cache_signature"]=canonical_sha256(payload)
    label=sample_set_label(config,names); run_id={"All eligible samples":"all_samples","Parents":"parents","Treated":"treated"}.get(label,f"custom_{payload['selected_samples_sha256'][:10]}")
    payload.update(run_id=run_id,sample_set_label=label); return payload
def cache_reusable(manifest,context,force=False): return not force and manifest.get("status") in {"success","cached"} and manifest.get("cache_signature")==context.get("cache_signature") and Path(manifest.get("output_dir","")).is_dir()
def find_matching_cache(config,context,force=False):
    base=_base_paths(config); path=base["project"]/"manifests"/"library_diagnostics"/"runs"/context["run_id"]/"manifest.json"
    try: manifest=json.loads(path.read_text())
    except (OSError,json.JSONDecodeError): return None
    return manifest if cache_reusable(manifest,context,force) else None
def _run_paths(config,context):
    base=_base_paths(config); run=context["run_id"]; manifests=base["project"]/"manifests"/"library_diagnostics"
    return {**base,"output":base["stable_output"]/"runs"/run,"export":base["export_root"]/run,"log":base["project"]/"logs"/"library_diagnostics"/f"{run}.log","manifest":manifests/"runs"/run/"manifest.json","status":manifests/"library_diagnostics_status.csv"}
def build_library_diagnostics_command(config,samples=None,threads=None,output_dir=None):
    if threads is not None and int(threads)<1: raise ValueError("Threads must be at least 1.")
    selected=_select_samples(config,samples); ref=_reference_files(config); out=Path(output_dir) if output_dir else _base_paths(config)["stable_output"]
    command=[sys.executable,str(Path(config["_repo_root"])/"src"/"ytab"/"qc"/"LibraryDiagnostics.py"),"--hits-txt",*[str(find_sample_hits_file(x,config)) for x in selected],"--outdir",str(out)]
    if ref["fasta"] and ref["fasta"].is_file(): command += ["--fasta",str(ref["fasta"])]
    if ref["centromere_bed"] and ref["centromere_bed"].is_file(): command += ["--centromeres-bed",str(ref["centromere_bed"])]
    if ref["gff"] and ref["gff"].is_file(): command += ["--metaplot-gff",str(ref["gff"])]
    return command
def _detected_outputs(path): return sorted(str(x) for x in Path(path).rglob("*") if x.is_file()) if Path(path).is_dir() else []
def _write_manifest(path,value): path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(value,indent=2)+"\n")
def _write_status(path,m):
    row={"project_id":m["project_id"],"run_id":m["run_id"],"sample_set_label":m["sample_set_label"],"cache_signature":m["cache_signature"],"status":m["status"],"samples":";".join(m["selected_samples"]),"input_hit_files_count":len(m["input_files"]),"output_dir":m["output_dir"],"export_dir":m["export_dir"],"log_file":m["log_file"],"elapsed_seconds":m["elapsed_seconds"],"message":m.get("error_message") or "; ".join(m.get("warnings") or [])}
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as handle: writer=csv.DictWriter(handle,fieldnames=STATUS_FIELDS);writer.writeheader();writer.writerow(row)
def _copy_tree(source,target):
    target.mkdir(parents=True,exist_ok=True)
    for item in source.rglob("*"):
        if item.is_file(): dest=target/item.relative_to(source);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(item,dest)
def _publish_stable(run_dir,stable_dir):
    stable_dir.mkdir(parents=True,exist_ok=True)
    for item in run_dir.iterdir():
        dest=stable_dir/item.name
        if item.is_dir():
            if dest.exists(): shutil.rmtree(dest)
            shutil.copytree(item,dest)
        else: shutil.copy2(item,dest)

def run_library_diagnostics_project(project_config,samples=None,threads=None,force=False,dry_run=False,job_id=None,progress_file=None,scientific_parameters=None):
    config=load_project_for_library_diagnostics(project_config); errors=validate_library_diagnostics_inputs(config,samples)
    if errors: raise ValueError("\n".join(errors))
    context=build_cache_context(config,samples,threads,scientific_parameters); paths=_run_paths(config,context); command=build_library_diagnostics_command(config,context["selected_samples_sorted"],threads,paths["output"])
    job_id=job_id or make_job_id("library_diagnostics"); progress_file=Path(progress_file) if progress_file else paths["project"]/"manifests"/"jobs"/f"{job_id}.progress.json"
    tracker=ProgressTracker.create(progress_file,job_id,str(config.get("project_id")),"Library Diagnostics",[context["sample_set_label"]],command=command,input_files={context["sample_set_label"]:context["input_hit_files"][0]});tracker.state.update(dry_run=dry_run,execution_mode="preview" if dry_run else "run",sample_set_label=context["sample_set_label"],selected_sample_count=len(context["selected_samples_sorted"]));tracker.start("Validating inputs");tracker.start_item(context["sample_set_label"],1,"validating inputs")
    started=time.monotonic(); cached=find_matching_cache(config,context,force)
    manifest={**context,"project_id":config.get("project_id"),"species":config.get("species"),"selected_samples":context["selected_samples_sorted"],"input_files":context["input_hit_files"],"input_hashes":context["input_hit_file_hashes"],"output_dir":str(paths["output"]),"export_dir":str(paths["export"]),"log_file":str(paths["log"]),"manifest_path":str(paths["manifest"]),"progress_file":str(progress_file),"status":"failed","cached_from":None,"start_time":_now(),"end_time":None,"elapsed_seconds":0.0,"command_run":command,"detected_outputs":[],"warnings":_resource_warnings(config),"error_message":None,"scientific_parameters":context["scientific_parameters"]}
    if dry_run:
        manifest["status"]="preview"; tracker.phase("collecting outputs","Preview complete; diagnostics were not executed")
    elif cached:
        manifest={**cached,"status":"cached","cached_from":str(paths["manifest"]),"command_run":command,"progress_file":str(progress_file)}; tracker.phase("collecting outputs","Matching cached diagnostics found")
    else:
        paths["output"].mkdir(parents=True,exist_ok=True);paths["log"].parent.mkdir(parents=True,exist_ok=True)
        try:
            tracker.phase("loading hit files","Loading selected hit files");tracker.phase("running diagnostics","LibraryDiagnostics running")
            with paths["log"].open("w") as log:
                log.write("Command: "+" ".join(command)+"\n\n");log.flush();completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,text=True)
            if completed.returncode: raise RuntimeError(f"LibraryDiagnostics exited with status {completed.returncode}; see {paths['log']}")
            expected=paths["output"]/"library_diagnostics.summary.csv"
            if not expected.is_file(): raise RuntimeError("LibraryDiagnostics completed without its combined summary CSV.")
            tracker.phase("collecting outputs","Collecting diagnostic outputs");_copy_tree(paths["output"],paths["export"]);tracker.phase("exporting results","Publishing stable current result");_publish_stable(paths["output"],paths["stable_output"]);manifest["status"]="success"
        except Exception as exc: manifest["error_message"]=str(exc)
    manifest["end_time"]=_now();manifest["elapsed_seconds"]=round(time.monotonic()-started,3);manifest["detected_outputs"]=_detected_outputs(paths["output"])
    if not dry_run: _write_manifest(paths["manifest"],manifest);_write_status(paths["status"],manifest)
    final="failed" if manifest["status"]=="failed" else "dry_run_success" if dry_run else "success"; message="Preview complete; diagnostics were not executed." if dry_run else "Cached result reused" if manifest["status"]=="cached" else "Library Diagnostics complete";tracker.finish_item(context["sample_set_label"],"failed" if final=="failed" else "skipped" if manifest["status"]=="cached" else "success",manifest["detected_outputs"],message,manifest.get("error_message") or "");tracker.finalize(final,message)
    return {**manifest,"included_samples":len(get_included_samples(config)),"status_path":str(paths["status"]),"dry_run":dry_run,"job_id":job_id}
