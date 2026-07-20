#!/usr/bin/env python3
"""Check the local YTAB Python, command-line, and R dependencies."""

from __future__ import annotations

import importlib
import shutil
import subprocess
import sys


PYTHON_MODULES = {
    "pandas": "pandas",
    "PyYAML (yaml)": "yaml",
    "Biopython (Bio)": "Bio",
    "pysam": "pysam",
    "gffutils": "gffutils",
    "numpy": "numpy",
    "scipy": "scipy",
    "scikit-learn (sklearn)": "sklearn",
    "matplotlib": "matplotlib",
    "openpyxl": "openpyxl",
}
COMMANDS = ("bowtie2", "samtools", "Rscript")
R_PACKAGES = ("shiny", "processx", "jsonlite", "yaml", "DT", "bslib", "tidyverse")


def report(ok: bool, label: str, detail: str = "") -> bool:
    suffix = f" ({detail})" if detail else ""
    print(f"{'PASS' if ok else 'FAIL'}: {label}{suffix}")
    return ok


def check_r_package(rscript: str, package: str) -> bool:
    expression = f"quit(status=if (requireNamespace('{package}', quietly=TRUE)) 0 else 1)"
    result = subprocess.run(
        [rscript, "--vanilla", "-e", expression],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def main() -> int:
    results: list[bool] = []
    version_ok = sys.version_info >= (3, 10)
    results.append(report(version_ok, "Python >= 3.10", sys.version.split()[0]))

    for label, module_name in PYTHON_MODULES.items():
        try:
            importlib.import_module(module_name)
        except Exception as exc:  # imports may fail because a shared library is absent
            results.append(report(False, f"Python import {label}", str(exc)))
        else:
            results.append(report(True, f"Python import {label}"))

    found_commands: dict[str, str | None] = {}
    for command in COMMANDS:
        found_commands[command] = shutil.which(command)
        results.append(report(found_commands[command] is not None, f"command {command}", found_commands[command] or "not on PATH"))

    rscript = found_commands["Rscript"]
    for package in R_PACKAGES:
        ok = bool(rscript) and check_r_package(rscript, package)
        results.append(report(ok, f"R package {package}", "installed" if ok else "missing"))

    if all(results):
        print("PASS: local YTAB environment is ready.")
        return 0
    print("FAIL: local YTAB environment is missing critical dependencies.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
