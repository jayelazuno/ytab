"""CLI dispatch, scientific parity and extensible master workflow tests."""
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from ytab.cli import main as dispatcher
from ytab.cli.coverage import main as coverage
from ytab.cli.registry import ROOT, MODULES, command
from ytab.cli.workflow import build_plan, execute_plan, main as workflow
from ytab.screen_coverage.statistics import colonies_for_target_coverage, poisson_coverage


def test_registry_paths_and_aliases():
    for name, selected in MODULES.items():
        assert Path(command(name, ["--help"])[1]).is_file(), selected.script
    assert command("screen-coverage", []) == command("coverage", [])
    assert command("mutant-identification", []) == command("mutants", [])


def test_top_level_help_and_json(capsys):
    assert dispatcher.main([]) == 0
    assert "modular analysis" in capsys.readouterr().out
    assert dispatcher.main(["modules", "--json"]) == 0
    items = json.loads(capsys.readouterr().out)
    assert {m["name"] for m in items} == set(MODULES)
    assert dispatcher.main(["does-not-exist"]) == 2


def test_dispatch_forwards_backend_exit_and_arguments():
    with patch.object(dispatcher.subprocess, "run", return_value=SimpleNamespace(returncode=7)) as runner:
        assert dispatcher.main(["fitness", "--help"]) == 7
    args, kwargs = runner.call_args
    assert args[0] == command("fitness", ["--help"])
    assert not kwargs.get("shell", False)
    assert "MPLBACKEND" in kwargs["env"]


@pytest.mark.parametrize("n,q", [(1,95), (4379,95), (5220,95), (5220,99.9)])
def test_poisson_exports_use_unchanged_math(tmp_path, capsys, n, q):
    out = tmp_path/"coverage"
    args = ["--genes",str(n),"--target",str(q),"--output",str(out),"--no-plot","--json"]
    assert coverage(args) == 0
    result = json.loads(capsys.readouterr().out)
    expected = colonies_for_target_coverage(n,q/100)
    assert result["required_picks"] == expected
    assert result["expected_recovery_percent"] == 100*poisson_coverage(expected,n)
    assert not result["reused"]
    assert coverage(args) == 0
    assert json.loads(capsys.readouterr().out)["reused"]
    assert (out/"curve.csv").exists() and (out/"summary.csv").exists()
    if n == 5220 and q == 95:
        assert expected == 15638


def test_coverage_plot_and_overwrite_guard(tmp_path):
    args = ["--genes","10","--output",str(tmp_path)]
    assert coverage(args) == 0
    assert (tmp_path/"coverage.png").read_bytes().startswith(b"\x89PNG")
    before = (tmp_path/"coverage.json").read_bytes()
    assert coverage(args+["--target","90"]) == 1
    assert (tmp_path/"coverage.json").read_bytes() == before
    assert coverage(args+["--target","90","--force"]) == 0


@pytest.mark.parametrize("args", [["--genes","0"], ["--genes","10","--target","100"], ["--genes","10","--target","nan"]])
def test_invalid_coverage_no_writes(tmp_path,args):
    with pytest.raises(SystemExit) as error:
        coverage(args+["--output",str(tmp_path/"absent")])
    assert error.value.code == 2
    assert not (tmp_path/"absent").exists()


def test_workflow_configuration_injection(tmp_path):
    project=tmp_path/"project.yaml"; project.touch()
    config=dict(version=1,project_config="project.yaml",steps=[
        dict(id="qc",module="library-qc"),
        dict(id="mutant",module="mutants",args=["--input","reads.ab1"]),
        dict(id="design",module="comparison-design"),
        dict(id="benchmark",module="coverage",depends_on=[],args=["--genes","10"])])
    plan=build_plan(config,tmp_path,threads=3,force=True)
    assert ["--project-config",str(project)] == plan[0]["command"][-5:-3]
    assert "--threads" in plan[0]["command"] and "--force" in plan[0]["command"]
    assert "--reference-config" in plan[1]["command"]
    assert "--link-project-config" not in plan[1]["command"]
    assert "--force" not in plan[1]["command"]
    assert "--overwrite" in plan[2]["command"]
    assert "--threads" not in plan[3]["command"]
    assert plan[1]["depends_on"] == ["qc"]
    assert plan[3]["depends_on"] == []


@pytest.mark.parametrize("steps", [
    [dict(module="no-module")], [dict(module="app")],
    [dict(module="coverage",args="echo hello")],
    [dict(module="coverage",depends_on=["future"])],
    [dict(id="../unsafe",module="coverage")],
    [dict(id="x",module="coverage"),dict(id="x",module="coverage")],
    [dict(module="fitness")],
])
def test_workflow_rejects_invalid_config(tmp_path,steps):
    with pytest.raises(ValueError):
        build_plan(dict(steps=steps),tmp_path)


def test_keep_going_blocks_dependents_but_runs_independent(tmp_path):
    plan=build_plan(dict(steps=[dict(id="fail",module="coverage",args=["--help"]),
        dict(id="blocked",module="coverage",args=["--help"]),
        dict(id="independent",module="coverage",depends_on=[],args=["--help"])]),tmp_path)
    with patch("ytab.cli.workflow.subprocess.run",side_effect=[SimpleNamespace(returncode=2),SimpleNamespace(returncode=0)]) as runner:
        state=execute_plan(plan,tmp_path/"logs",keep_going=True)
    assert runner.call_count == 2
    assert [r["status"] for r in state["results"]] == ["failed","blocked","success"]
    assert state["status"] == "failed"
    assert json.loads(Path(state["manifest"]).read_text())["status"] == "failed"


def test_stop_and_cancellation(tmp_path):
    plan=build_plan(dict(steps=[dict(module="coverage",args=["--help"]),dict(module="coverage",args=["--help"])]),tmp_path)
    with patch("ytab.cli.workflow.subprocess.run",return_value=SimpleNamespace(returncode=1)) as runner:
        state=execute_plan(plan,tmp_path/"logs")
    assert runner.call_count == 1
    assert [r["status"] for r in state["results"]] == ["failed","not-run"]
    with patch("ytab.cli.workflow.subprocess.run",side_effect=KeyboardInterrupt):
        state=execute_plan(plan,tmp_path/"logs")
    assert state["status"] == "cancelled"


def test_dry_run_no_output_and_bad_yaml(tmp_path,capsys):
    config=tmp_path/"workflow.yaml"
    config.write_text('version: 1\nsteps:\n  - module: coverage\n    args: ["--genes", "10"]\n')
    out=tmp_path/"no-output"
    assert workflow(["--workflow",str(config),"--dry-run","--output",str(out)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "dry-run"
    assert not out.exists()
    config.write_text('steps: [bad YAML')
    assert workflow(["--workflow",str(config)]) == 2


def test_real_master_workflow(tmp_path):
    config=tmp_path/"workflow.yaml"
    config.write_text('version: 1\nsteps:\n  - module: coverage\n    args: ["--genes", "10", "--no-plot", "--output", "benchmark"]\n  - module: mutants\n    args: ["--help"]\n')
    assert workflow(["--workflow",str(config),"--output",str(tmp_path/"logs")]) == 0
    assert (tmp_path/"benchmark/coverage.json").is_file()
    manifest=next((tmp_path/"logs").glob("*/workflow.json"))
    assert json.loads(manifest.read_text())["status"] == "success"


def test_module_entrypoint_headless():
    result=subprocess.run([sys.executable,"-m","ytab","modules","--json"],capture_output=True,text=True)
    assert result.returncode == 0, result.stderr
    assert any(m["name"]=="mutants" for m in json.loads(result.stdout))
