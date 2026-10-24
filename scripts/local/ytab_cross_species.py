#!/usr/bin/env python3
"""Essentiality across species: label-free HMM calls of one-to-one orthologs in two or more projects (used by the app's
Gene Explorer > Across Species > Essentiality across species).

Genes are matched through one hub, S. cerevisiae systematic names (YAL001C ...): each project's reference config names
an ortholog table (`reference.orthology_file`) that maps its genes to S. cerevisiae; an S. cerevisiae project is the
hub itself. Only one-to-one orthologs are kept (a gene mapped to exactly one S. cerevisiae gene and back, per species),
and only hub genes with an ortholog scored by the HMM in every species, so each species is compared on the same genes.

Ortholog tables: tab-separated, '#' comments allowed. With a header that has a `scer_id` column, the first column is the
species' gene and `scer_id` its ortholog; without one, each row's first field is the species' gene and the first later
field that is a single S. cerevisiae systematic name is its ortholog (rows naming several genes, e.g. "YBR265W+YDR302W",
are ambiguous and skipped). This covers CGD, PomBase and YTAB-curated tables; a new species needs one such table.

Calls come from each project's essentiality run (newest unless --essentiality-run is given): P_ES_run and the HMM call
(ES or not); genes the HMM did not score (unmappable) are left out.

Outputs in OUTPUT/:
  orthologs.tsv      hub gene and, per species, the gene, its name, P_ES_run and ES (1/0)
  intersections.tsv  genes by the set of species in which they are essential (shared and species-specific sets)
  pairwise.tsv       per species pair: genes, both essential, either essential, Cohen's kappa of the calls (Cohen
                     1960) and the Spearman correlation of P_ES_run
  manifest.json

  ytab cross-species --project-config output/projects/Gale_2020/config/project.yaml \\
      --project-config output/projects/Pombe_Hermes/config/project.yaml --output output/projects/Gale_2020/cross_species/run1
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from itertools import combinations
import json
from pathlib import Path
import re

import numpy as np
from scipy.stats import spearmanr
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCER = re.compile(r"^(Y[A-P][LR]\d{3}[WC](-[A-Z])?|Q\d{4})$")


def species_label(key):
    """'glabrata' -> 'C. glabrata', from the comparative species manifest; the key itself when it is not listed."""
    path = ROOT/"resources/comparative/manifest/comparative_species_manifest.yaml"
    labels = {k:(v or {}).get("label", "") for k, v in (yaml.safe_load(path.read_text()) or {}).items()} if path.exists() else {}
    words = (labels.get(key) or "").split()
    return f"{words[0][0]}. {' '.join(words[1:])}" if len(words) >= 2 else key


def read_orthologs(path):
    """{species gene: S. cerevisiae gene}, one-to-one pairs only."""
    pairs = set()
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle if line.strip() and not line.startswith("#")]
    if lines and "scer_id" in lines[0]:
        col = lines[0].index("scer_id"); lines = lines[1:]
        pairs = {(f[0].strip(), f[col].strip()) for f in lines if len(f) > col and SCER.match(f[col].strip())}
    else:
        for f in lines:
            hub = next((x.strip() for x in f[1:] if SCER.match(x.strip())), None)
            if hub and f[0].strip():pairs.add((f[0].strip(), hub))
    na, nb = Counter(a for a, _ in pairs), Counter(b for _, b in pairs)
    return {a:b for a, b in pairs if na[a] == 1 and nb[b] == 1}


def newest_run(project_dir):
    runs = sorted((p for p in (project_dir/"essentiality").glob("*/manifest.json")), key=lambda p:p.stat().st_mtime, reverse=True)
    if not runs:raise ValueError(f"{project_dir.name}: no essentiality run; run essentiality first")
    return runs[0].parent


def hmm_calls(run):
    out = {}
    with open(run/"models/genes.tsv") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r.get("P_ES_run") not in ("", "nan", None):
                out[r["ORF"].removeprefix("gene:")] = (float(r["P_ES_run"]), r.get("hmm_call") == "ES", r.get("Name") or r["ORF"])
    return out


def kappa(a, b):
    po = (a == b).mean(); pe = a.mean()*b.mean()+(1-a.mean())*(1-b.mean())
    return float((po-pe)/(1-pe)) if pe < 1 else float("nan")


def write(path, rows, fields=None):
    with open(path, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def build(configs, runs, output):
    species, maps, calls, sources = [], {}, {}, []
    for i, config in enumerate(configs):
        cfg = yaml.safe_load(Path(config).read_text()); key = cfg.get("species") or cfg["reference"]["species"]
        label = species_label(key)
        if label in species:raise ValueError(f"{label} is given twice; choose one project per species")
        project_dir = ROOT/(cfg.get("output_project_dir") or f"output/projects/{cfg['project_id']}")
        run = Path(runs[i]).resolve() if runs and runs[i] else newest_run(project_dir)
        orthology = (cfg.get("reference") or {}).get("orthology_file")
        if key == "cerevisiae":maps[label] = None
        elif orthology:maps[label] = read_orthologs(ROOT/orthology)
        else:raise ValueError(f"{cfg['project_id']}: no ortholog table to S. cerevisiae (reference.orthology_file in its project.yaml)")
        species.append(label); calls[label] = hmm_calls(run)
        sources.append(dict(project=cfg["project_id"], species=label, essentiality_run=str(run), orthology_file=orthology,
                            one_to_one_orthologs=len(maps[label]) if maps[label] is not None else None, scored_genes=len(calls[label])))
    if len(species) < 2:raise ValueError("Give at least two projects of different species")
    to_species = {s:({b:a for a, b in m.items()} if m is not None else None) for s, m in maps.items()}
    hubs = set.intersection(*(set(v) if v is not None else set(calls[s]) for s, v in to_species.items()))
    rows = []
    for hub in sorted(hubs):
        genes = {s:(hub if to_species[s] is None else to_species[s][hub]) for s in species}
        if not all(genes[s] in calls[s] for s in species):continue
        row = dict(scer_id=hub)
        for s in species:
            p, es, name = calls[s][genes[s]]
            row.update({f"{s} gene":genes[s], f"{s} name":name, f"{s} P_ES_run":round(p, 4), f"{s} ES":int(es)})
        rows.append(row)
    if not rows:raise ValueError("No one-to-one orthologs scored in every species")
    output.mkdir(parents=True)
    write(output/"orthologs.tsv", rows)
    combos = Counter(tuple(s for s in species if r[f"{s} ES"]) for r in rows)
    write(output/"intersections.tsv", [dict(essential_in=" + ".join(c) or "none", n_species=len(c), genes=n, **{s:int(s in c) for s in species})
                                       for c, n in combos.most_common()])
    pairwise = []
    for a, b in combinations(species, 2):
        ca, cb = (np.array([r[f"{s} ES"] for r in rows]) for s in (a, b))
        pa, pb = (np.array([r[f"{s} P_ES_run"] for r in rows]) for s in (a, b))
        pairwise.append(dict(species_a=a, species_b=b, genes=len(rows), both_ES=int((ca & cb).sum()), either_ES=int((ca | cb).sum()),
                             kappa_ES=round(kappa(ca, cb), 4), spearman_P_ES_run=round(float(spearmanr(pa, pb).statistic), 4)))
    write(output/"pairwise.tsv", pairwise)
    manifest = dict(status="success", species=species, projects=sources, orthologs=len(rows),
                    essential_in_all=combos.get(tuple(species), 0), essential_in_none=combos.get((), 0),
                    hub="S. cerevisiae systematic names; one-to-one orthologs scored in every species",
                    created=datetime.now(timezone.utc).isoformat(),
                    files=dict(orthologs="orthologs.tsv", intersections="intersections.tsv", pairwise="pairwise.tsv"))
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", action="append", required=True, help="Repeat: one project per species (at least two)")
    p.add_argument("--essentiality-run", action="append", help="Repeat in the same order as --project-config (default: each project's newest run)")
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    if a.essentiality_run and len(a.essentiality_run) != len(a.project_config):p.error("give one --essentiality-run per --project-config")
    out = a.output.resolve()
    try:
        if out.exists() and any(out.iterdir()):raise ValueError("Output directory is nonempty; choose a new one")
        if out.exists():out.rmdir()
        m = build(a.project_config, a.essentiality_run, out)
    except (ValueError, KeyError, OSError) as error:
        p.exit(1, f"ERROR: {error}\n")
    print(json.dumps({k:m[k] for k in ("status", "species", "orthologs", "essential_in_all", "essential_in_none")}, indent=1))


if __name__ == "__main__":main()
