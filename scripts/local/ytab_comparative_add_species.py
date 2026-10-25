#!/usr/bin/env python3
"""Add a species to the app's comparative orthology tables (resources/comparative/orthology/) through S. cerevisiae.

The comparative tables (gene_lookup.csv, orthogroups_long.csv) were imported from an OrthoFinder analysis that did not
include every species YTAB supports. A species with a gene annotation (resources/species/<species>/annotation/
gene_annotation.tsv, whose scer_gene_id column holds its S. cerevisiae ortholog) joins each orthogroup that contains its
S. cerevisiae ortholog; its rows are marked source "via S. cerevisiae ortholog" so they can be told from the
OrthoFinder groups. Rerunning replaces the species' earlier rows. The species is then enabled in
resources/comparative/manifest/comparative_species_manifest.yaml.

  python scripts/local/ytab_comparative_add_species.py --species pombe --label "Schizosaccharomyces pombe"
"""
import argparse
import csv
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
ORTHO = ROOT/"resources/comparative/orthology"
MANIFEST = ROOT/"resources/comparative/manifest/comparative_species_manifest.yaml"
SOURCE = "via S. cerevisiae ortholog"


def read(path):
    """Header, parsed rows and the raw lines (existing lines are written back unchanged)."""
    lines = Path(path).read_text().splitlines()
    reader = csv.DictReader(lines); rows = list(reader)
    return reader.fieldnames, rows, lines[1:]


def line(fields, row):
    """One row in the tables' style: values quoted, empty values bare."""
    return ",".join('"%s"' % str(row.get(f, "")).replace('"', '""') if row.get(f, "") != "" else "" for f in fields)


def write(path, fields, header, kept, new):
    Path(path).write_text("\n".join([header]+kept+[line(fields, r) for r in new])+"\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True); p.add_argument("--label", required=True, help="Full species name, e.g. Schizosaccharomyces pombe")
    a = p.parse_args()
    with open(ROOT/"resources/species"/a.species/"annotation/gene_annotation.tsv") as handle:
        genes = [r for r in csv.DictReader(handle, delimiter="\t") if r.get("scer_gene_id")]
    og_fields, groups, og_lines = read(ORTHO/"orthogroups_long.csv")
    lk_fields, lookup, lk_lines = read(ORTHO/"gene_lookup.csv")
    og_header = (ORTHO/"orthogroups_long.csv").read_text().splitlines()[0]; lk_header = (ORTHO/"gene_lookup.csv").read_text().splitlines()[0]
    og_kept = [l for r, l in zip(groups, og_lines) if r["species"] != a.species]; lk_kept = [l for r, l in zip(lookup, lk_lines) if r["species"] != a.species]
    new_groups, new_lookup = [], []
    og_of = {}
    for r in groups:
        if r["species"] == "cerevisiae":og_of.setdefault(r["gene_id"], []).append(r["orthogroup_id"])
    hog_of = {r["gene_id"]:r.get("hog_id", "") for r in lookup if r["species"] == "cerevisiae"}
    added = 0
    for g in genes:
        scer = [s for s in g["scer_gene_id"].split(";") if s]
        name = g.get("gene_name") or g["systematic_id"]
        ogs = sorted({og for s in scer for og in og_of.get(s, [])})
        for og in ogs:new_groups.append(dict(orthogroup_id=og, species=a.species, gene_id=g["systematic_id"], gene_name=name, source_column=SOURCE))
        new_lookup.append(dict(gene_id=g["systematic_id"], species=a.species, expression_id=g["gene_id"], id_type="systematic", source_info=SOURCE,
                           gene_name=g.get("gene_name", ""), hog_id=hog_of.get(scer[0], "") if scer else "", og_id=";".join(ogs)))
        added += bool(ogs)
    write(ORTHO/"orthogroups_long.csv", og_fields, og_header, og_kept, new_groups)
    write(ORTHO/"gene_lookup.csv", lk_fields, lk_header, lk_kept, new_lookup)
    # the species' manifest entry is replaced; the other entries are kept as written
    text = MANIFEST.read_text(); blocks = [b for b in text.strip().split("\n\n") if not b.startswith(f"{a.species}:")]
    genus = a.label.split()[0]
    blocks.append(f"{a.species}:\n  label: {a.label}\n  enabled: true\n  ytab_species_key: {a.species}\n"
                  f"  rnacross_aliases: [{genus[0].lower()}{a.species[0]}, {genus[0].lower()}{a.species[:3]}, {genus[0]}_{a.species}, {a.label}]\n"
                  f"  note: Orthogroups through S. cerevisiae orthologs.")
    MANIFEST.write_text("\n\n".join(blocks)+"\n")
    assert yaml.safe_load(MANIFEST.read_text())[a.species]["enabled"] is True
    print(f"{a.species}: {len(genes)} genes with an S. cerevisiae ortholog; {added} placed in orthogroups")


if __name__ == "__main__":main()
