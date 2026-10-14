"""Sequential, shell-free master workflows over registered YTAB CLI modules."""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from .registry import ROOT, command, module


def now():
    return datetime.now(timezone.utc).isoformat()


def has_option(arguments, option):
    return any(arg == option or arg.startswith(option+"=") for arg in arguments)


def build_plan(config, base, project=None, threads=None, force=False):
    if not isinstance(config, dict) or config.get("version", 1) != 1:
        raise ValueError("Workflow must be a mapping with version: 1")
    steps = config.get("steps")
    if not isinstance(steps, list) or not steps:
        raise ValueError("Workflow requires a nonempty steps list")
    project = project or config.get("project_config")
    if project:
        project = Path(project)
        project = (base/project).resolve() if not project.is_absolute() else project.resolve()
        if not project.is_file():
            raise ValueError(f"Project/reference configuration not found: {project}")
    result = []
    known = set()
    for i, step in enumerate(steps):
        if not isinstance(step, dict):
            raise ValueError(f"Step {i+1} must be a mapping")
        selected = module(step.get("module", ""))
        if not selected.workflow:
            raise ValueError(f"Interactive module {selected.name} cannot run inside a master workflow")
        identifier = step.get("id", f"{i+1}-{selected.name}")
        if not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", identifier) or identifier in known:
            raise ValueError(f"Invalid or duplicate step id: {identifier!r}")
        arguments = step.get("args", [])
        if not isinstance(arguments, list) or not all(isinstance(v, str) for v in arguments):
            raise ValueError(f"{identifier}: args must be a list of strings, not a shell command")
        arguments = list(arguments)
        if project and selected.project_option and not has_option(arguments, selected.project_option) and not has_option(arguments, "--project-config"):
            arguments += [selected.project_option, str(project)]
        if selected.project_option and not has_option(arguments, selected.project_option) and not has_option(arguments, "--project-config") and not any(a in ("--help", "-h") for a in arguments):
            raise ValueError(f"{identifier}: supply {selected.project_option} or a workflow project_config")
        if threads is not None and selected.threads and not has_option(arguments, "--threads"):
            arguments += ["--threads", str(threads)]
        if force and selected.force_option and not has_option(arguments, selected.force_option):
            arguments.append(selected.force_option)
        depends = step.get("depends_on", [result[-1]["id"]] if result else [])
        if not isinstance(depends, list) or not all(isinstance(d, str) and d in known for d in depends):
            raise ValueError(f"{identifier}: depends_on must name earlier steps; forward/cyclic dependencies are not allowed")
        result.append(dict(id=identifier, module=selected.name, depends_on=depends, command=command(selected.name, arguments), cwd=str(base)))
        known.add(identifier)
    return result


def write_state(path, state):
    temporary = path.with_suffix(".json.pending")
    temporary.write_text(json.dumps(state, indent=2)+"\n")
    temporary.replace(path)


def execute_plan(plan, output, keep_going=False):
    output = Path(output).resolve()
    run = output / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8])
    run.mkdir(parents=True)
    state = dict(schema_version=1, status="running", started=now(), plan=plan, results=[])
    manifest = run/"workflow.json"
    write_state(manifest, state)
    outcomes = {}
    env = dict(os.environ)
    env.setdefault("MPLBACKEND", "Agg")
    started = time.monotonic()
    try:
        for step in plan:
            if any(outcomes.get(dep) != "success" for dep in step["depends_on"]):
                entry = dict(id=step["id"], module=step["module"], status="blocked", reason="A required earlier step did not succeed")
            else:
                log = run/(step["id"]+".log")
                print(f"[{step['id']}] {shlex.join(step['command'])}", flush=True)
                try:
                    with log.open("w") as handle:
                        completed = subprocess.run(step["command"], cwd=step["cwd"], env=env, stdout=handle, stderr=subprocess.STDOUT)
                    code = completed.returncode
                    entry = dict(id=step["id"], module=step["module"], status="success" if code == 0 else "failed", returncode=code, log=str(log))
                except OSError as e:
                    entry = dict(id=step["id"], module=step["module"], status="failed", returncode=1, error=str(e), log=str(log))
                print(f"[{step['id']}] {entry['status']}; log: {log}", flush=True)
            outcomes[step["id"]] = entry["status"]
            state["results"].append(entry)
            write_state(manifest, state)
            if entry["status"] == "failed" and not keep_going:
                break
        finished = {e["id"] for e in state["results"]}
        for step in plan:
            if step["id"] not in finished:
                state["results"].append(dict(id=step["id"], module=step["module"], status="not-run", reason="Master run stopped after failure"))
        state["status"] = "failed" if any(e["status"] != "success" for e in state["results"]) else "success"
    except KeyboardInterrupt:
        state["status"] = "cancelled"
    state.update(finished=now(), elapsed_seconds=round(time.monotonic()-started, 3))
    write_state(manifest, state)
    state["manifest"] = str(manifest)
    return state


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run a YTAB master pipeline or ordered, extensible YAML workflow. --dry-run never executes modules or writes outputs.")
    parser.add_argument("--workflow", type=Path, help="YAML workflow; argument paths are relative to this file's directory")
    parser.add_argument("--project-config", type=Path)
    parser.add_argument("--profile", choices=("core", "classifier", "fitness", "all"), default="all")
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--analysis-id", default="treated_vs_parent_raw_cpm")
    parser.add_argument("--stages", help="Comma-separated existing pipeline stages")
    parser.add_argument("--from-stage")
    parser.add_argument("--to-stage")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT/"output/cli_workflows", help="Master manifests/logs, not an override of module scientific outputs")
    args = parser.parse_args(argv)
    try:
        if args.threads is not None and args.threads < 1:
            raise ValueError("--threads must be positive")
        project = args.project_config.resolve() if args.project_config else None
        if args.workflow:
            import yaml
            path = args.workflow.resolve()
            try:
                config = yaml.safe_load(path.read_text())
            except yaml.YAMLError as error:
                raise ValueError(f"Invalid workflow YAML: {error}") from error
            base = path.parent
            if any((args.stages, args.from_stage, args.to_stage)):
                raise ValueError("Pipeline stage selectors belong in a workflow pipeline step's args")
        else:
            if not project:
                raise ValueError("Supply --project-config for a standard master run, or --workflow")
            arguments = ["--profile", args.profile, "--target", args.target, "--analysis-id", args.analysis_id]
            for flag, value in (("--stages",args.stages),("--from-stage",args.from_stage),("--to-stage",args.to_stage)):
                if value:arguments += [flag, value]
            if args.keep_going:arguments.append("--keep-going")
            config = dict(version=1, steps=[dict(id="pipeline", module="pipeline", args=arguments)])
            base = Path.cwd()
        plan = build_plan(config, base, project, args.threads, args.force)
        if args.dry_run:
            print(json.dumps(dict(status="dry-run", plan=plan), indent=2))
            return 0
        result = execute_plan(plan, args.output, args.keep_going)
        print(f"Master status: {result['status']}\nManifest: {result['manifest']}")
        return 130 if result["status"] == "cancelled" else 1 if result["status"] == "failed" else 0
    except (ValueError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
