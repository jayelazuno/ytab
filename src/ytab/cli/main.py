"""One CLI for existing YTAB modules, ordered workflows, and the Shiny app."""
import json
import os
import subprocess
import sys
from .registry import MODULES, command, inventory


def help_text():
    rows = "\n".join(f"  {m.name:20} {m.description}" for m in MODULES.values())
    return f"""YTAB — modular analysis CLI

Usage:
  ytab modules [--json]
  ytab <module> --help
  ytab <module> <module options>
  ytab run --project-config PROJECT.yaml --profile core|classifier|fitness|all
  ytab run --workflow WORKFLOW.yaml [--dry-run]

Use python scripts/local/ytab.py or PYTHONPATH=src python -m ytab.
The app and CLI use the same analysis backends and result formats.
Run --help on any module for its scientific settings and output options.

  run                  Master run: existing pipeline or ordered YAML workflow
{rows}
"""


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        print(help_text())
        return 0
    if argv[0] == "modules":
        if argv[1:] not in ([], ["--json"]):
            print("Usage: ytab modules [--json]", file=sys.stderr)
            return 2
        print(json.dumps(inventory(), indent=2) if argv[1:] else help_text())
        return 0
    if argv[0] == "run":
        from .workflow import main as run_main
        return run_main(argv[1:])
    try:
        args = command(argv[0], argv[1:])
        env = dict(os.environ)
        env.setdefault("MPLBACKEND", "Agg")
        result = subprocess.run(args, env=env)
        return result.returncode if result.returncode >= 0 else 128-result.returncode
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except OSError as e:
        print(f"ERROR: could not start module: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130
