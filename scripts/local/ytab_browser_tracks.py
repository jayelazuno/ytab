#!/usr/bin/env python3
"""Genome browser files for a project (Gene Explorer > Genome browser reads only these): reference with display
chromosome names, aliases, gene annotation BED and one insertion track per sample with a CreateHitFile hit table.

  ytab browser-tracks --project-config output/projects/P/config/project.yaml [--force]

Writes <project>/genome_browser/ (see ytab.domain_explorer.browser_tracks); rebuilds only when an input is newer.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.domain_explorer import browser_tracks  # noqa: E402
from ytab.pipeline.project_config import _resolve_config_path, load_project_config  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path); p.add_argument("--force", action="store_true")
    a = p.parse_args()
    try:
        cfg = load_project_config(a.project_config); root = Path(cfg["_repo_root"])
        project = _resolve_config_path(cfg.get("output_project_dir") or f"output/projects/{cfg['project_id']}", root)
        ref = cfg.get("reference") or {}
        fasta, gff = _resolve_config_path(ref.get("fasta"), root), _resolve_config_path(ref.get("gff"), root)
        if not (fasta and fasta.is_file() and gff and gff.is_file()):raise ValueError("The project reference needs a FASTA and a GFF")
        species = str(cfg.get("species") or ref.get("species") or "")
        manifest = browser_tracks.build(project, root, species, fasta, gff, force=a.force)
        problems = browser_tracks.validate(project)
        if problems:raise ValueError("the browser files do not meet the contract: " + "; ".join(problems))
    except (ValueError, KeyError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr); return 1
    print(json.dumps(dict(status="success", output=str(project/"genome_browser"), genes=manifest["gene_count"],
                          tracks=[t["sample"] for t in manifest["tracks"]]), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
