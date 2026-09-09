#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ytab.go_analysis.runner import (  # noqa: E402
    LOOKUP_RELATIVE,
    build_universe,
    load_lookup,
    load_project,
    map_genes_to_scer,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Smoke-test YTAB GO analysis mapping.")
    parser.add_argument(
        "--project-config",
        type=Path,
        default=ROOT / "output/projects/Zn_toxicity_screen/config/project.yaml",
    )
    args = parser.parse_args()
    config, repo_root, project_root = load_project(args.project_config)
    lookup = load_lookup(repo_root)
    assert (repo_root / LOOKUP_RELATIVE).is_file()
    assert "cgla_gene_name_from_deseq" not in {
        "input_gene",
        "matched",
        "gwk60_id_clean",
        "cagl_id",
        "gene_name",
        "cg_to_sc_relationship",
        "scer_gene_id",
    }
    usable = lookup[lookup["scer_gene_id"].astype(str).str.len() > 0].copy()
    assert len(usable) > 0
    row = usable.iloc[0]
    probes = [
        row["gwk60_id_clean"],
        row["cagl_id"],
        row["scer_gene_id"].split(";")[0].split(",")[0],
    ]
    if str(row.get("scer_gene_name", "")).strip():
        probes.append(row["scer_gene_name"])
    mapped = map_genes_to_scer(probes, lookup)
    assert len(mapped) == len(probes)
    assert mapped["scer_gene_id"].astype(str).str.len().all()
    universe = build_universe(project_root, lookup)
    assert len(universe) > 0
    assert "scer_gene_id" in universe.columns
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
