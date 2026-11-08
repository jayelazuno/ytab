#!/usr/bin/env python3
"""Curate K. phaffii essentiality labels that are independent of transposon data.

Why not use the published calls. Zhu et al. 2018 (Sci Rep 8:10223) classify every gene from a logistic regression on
the gene insertion index of their own TcBuster/Sleeping beauty libraries, trained on orthologs of S. cerevisiae and
S. pombe essential genes. Those calls are therefore unusable as truth twice over: they come from transposon insertion
data, as the libraries YTAB analyses do, and their essential training set is itself ortholog-derived, so it overlaps
the independent labels built here. They are kept as `zhu2018_comparator` and are never labels - the rule the
C. glabrata labels already follow for Gale et al. 2020.

Evidence used instead, in order of precedence:
  gold   nonessential: the 171 K. phaffii genes Zhu et al. report as successfully knocked out (Dataset S3, the
         "non-essential genes" column). A viable knockout is direct evidence, obtained without transposons.
  silver both labels: the viability of the 1:1 S. cerevisiae ortholog in SGD, from this repository's own curated
         cerevisiae labels (gold tier only). Orthologs come from Zhu's Dataset S6, which pairs each P. pastoris
         locus tag with an S. cerevisiae systematic name; a pair is 1:1 when the locus tag has exactly one partner
         and that partner is listed for no other locus tag.

Zhu's own ortholog-derived essential training set (Dataset S3, "essential genes" column) is recorded as
`zhu2018_training_essential` for transparency and is not a label either: it is built the same way as the silver tier
and would double-count the same evidence.

  python scripts/local/ytab_curate_phaffii_labels.py

Outputs go to resources/species/phaffii/essentiality_labels/; the raw supplements must already be in raw/.
"""
import csv
import collections
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LABELS_DIR = ROOT/"resources/species/phaffii/essentiality_labels"
RAW = LABELS_DIR/"raw/zhu_2018_scirep"
S3 = RAW/"Dataset S3 training gene sets.xlsx"
S4 = RAW/"Dataset S4 gene insertion index and essentiality.xlsx"
S6 = RAW/"Dataset S6 S cerevisiae ortholog comparison.xlsx"
SCER_LABELS = ROOT/"resources/species/cerevisiae/essentiality_labels/labels.tsv"
GFF = ROOT/"resources/species/phaffii/reference_genome/GCF_000027005.1_ASM2700v1_genomic.gff"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows_of(path, sheet, skip=0):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True)[sheet]
    return list(ws.iter_rows(values_only=True))[skip:]


def cerevisiae_labels():
    out = {}
    with SCER_LABELS.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["tier"] == "gold" and r["label"] in ("essential", "nonessential"):out[r["gene_id"]] = r["label"]
    return out


def annotated_genes():
    """{locus tag: gene name} for the genes of the reference the libraries are mapped to."""
    names = {}
    for line in GFF.open():
        if line.startswith("#"):continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        tag = attrs.get("locus_tag") or attrs.get("ID", "").removeprefix("gene-")
        if tag:names[tag] = attrs.get("gene", "")
    return names


def main():
    scer, annotated = cerevisiae_labels(), annotated_genes()

    # silver: viability of the 1:1 S. cerevisiae ortholog (Dataset S6 gives systematic names directly)
    pairs = [(r[0], r[4]) for r in rows_of(S6, "Sc ortholog comparison", skip=2) if r[0] and r[4]]
    per_gene, per_ortholog = collections.Counter(p for p, _ in pairs), collections.Counter(s for _, s in pairs)
    one_to_one = {p:s for p, s in pairs if per_gene[p] == 1 and per_ortholog[s] == 1}
    silver = {p:(scer[s], s) for p, s in one_to_one.items() if s in scer}

    # gold: a viable knockout was made; and their ortholog-derived essential training set, for transparency only
    training = rows_of(S3, "Dataset S3", skip=1)
    knocked_out = {r[1] for r in training if r[1]}
    training_essential = {r[0] for r in training if r[0]}

    comparator = {r[0]:dict(category=str(r[5] or ""), nep=r[4]) for r in rows_of(S4, "Gene essentiality", skip=1) if r[0]}

    rows, counts = [], collections.Counter()
    for gene in sorted(set(annotated) | set(one_to_one) | knocked_out):
        gold = "nonessential" if gene in knocked_out else ""
        silver_label, ortholog = silver.get(gene, ("", ""))
        label, tier = (gold, "gold") if gold else (silver_label, "silver" if silver_label else "")
        conflict = silver_label if gold and silver_label and silver_label != gold else ""
        reason = "" if label else ("no 1:1 S. cerevisiae ortholog in Zhu Dataset S6" if gene not in one_to_one
                                   else "1:1 ortholog has no curated SGD viability")
        rows.append(dict(gene_id=gene, name=annotated.get(gene, ""), label=label, tier=tier,
                         unlabelled_reason=reason, scer_ortholog=ortholog, conflict=conflict,
                         viable_knockout=str(gene in knocked_out).lower(),
                         zhu2018_comparator=comparator.get(gene, {}).get("category", ""),
                         zhu2018_nep=comparator.get(gene, {}).get("nep", ""),
                         zhu2018_training_essential=str(gene in training_essential).lower(),
                         in_reference_annotation=str(gene in annotated).lower()))
        counts[(tier or "unlabelled", label or "none")] += 1

    LABELS_DIR.mkdir(parents=True, exist_ok=True)
    with (LABELS_DIR/"labels.tsv").open("w") as handle:
        w = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    their = {"Essential":"essential", "Non-essential":"nonessential"}
    agreement = collections.Counter((r["label"], their.get(r["zhu2018_comparator"], "")) for r in rows
                                    if r["label"] and their.get(r["zhu2018_comparator"], ""))
    agreed = sum(n for (a, b), n in agreement.items() if a == b)
    summary = dict(
        genes=len(rows), reference_annotation=len(annotated),
        by_tier={f"{t}:{l}":n for (t, l), n in sorted(counts.items())},
        gold_source="viable knockout reported by Zhu et al. 2018 Dataset S3 (non-essential column)",
        silver_source="1:1 S. cerevisiae ortholog viability (SGD, via this repository's curated cerevisiae labels)",
        one_to_one_orthologs=len(one_to_one), silver_labels=len(silver),
        conflicts_gold_over_silver=sum(1 for r in rows if r["conflict"]),
        comparator_agreement=dict(labelled_and_called=sum(agreement.values()), agree=agreed,
                                  fraction=agreed/sum(agreement.values()) if agreement else None,
                                  table={f"{a} vs {b}":n for (a, b), n in sorted(agreement.items())}))
    (LABELS_DIR/"summary.json").write_text(json.dumps(summary, indent=1)+"\n")
    (LABELS_DIR/"sources.json").write_text(json.dumps(dict(
        zhu_2018=dict(citation="Zhu J, Gong R, Zhu Q et al. (2018) Genome-wide determination of gene essentiality by "
                               "transposon insertion sequencing in yeast Pichia pastoris. Sci Rep 8:10223",
                      doi="10.1038/s41598-018-28217-z",
                      files={p.name:sha(p) for p in sorted(RAW.glob("*.xlsx"))}),
        cerevisiae_labels=dict(path=str(SCER_LABELS.relative_to(ROOT)), sha256=sha(SCER_LABELS)),
        reference_gff=dict(path=str(GFF.relative_to(ROOT)), sha256=sha(GFF)),
        script=dict(path=str(Path(__file__).relative_to(ROOT)), sha256=sha(__file__))), indent=1)+"\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
