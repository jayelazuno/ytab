#!/usr/bin/env python3
"""Species resources for C. albicans and S. cerevisiae, built from SGD/CGD files with fixed rules.

For each species writes, under resources/species/<species>/:
- reference_genome/: (cerevisiae) the S288C FASTA and GFF from an SGD release archive (--sgd-release, the
  S288C_reference_genome_*.tgz); centromeres.<source>.bed and tRNAs.bed from the GFF; chromosome_names.tsv
  (contig, display name, order) for plots.
- annotation/gene_annotation.tsv: every protein-coding gene of the mapped reference with its name, description,
  S. cerevisiae ortholog (1:1 or not) and that ortholog's SGD null-mutant viability.
- essentiality_labels/: labels.tsv, summary.json, sources.json, README.md and the raw evidence files.

Essentiality label rules (they never look at any YTAB model output):
- S. cerevisiae: SGD null-mutant viable/inviable; systematic deletion sets (haploid or homozygous diploid) decide
  when they agree, otherwise classical genetics when it agrees; mixed -> conflict. Tier gold (direct evidence).
- C. albicans direct (gold), from CGD: inviable null in a homozygous diploid or haploid, or inviable repressible
  allele under repression -> essential; viable null in a homozygous diploid or haploid -> nonessential; inviable only
  under another condition -> conditional (not labelled); both -> conflict (not labelled). Heterozygous-diploid
  viability is not evidence (one working copy remains). Transposon/machine-learning survey calls (PMID 30377286) are
  kept as a comparator column, never as labels: they come from the same kind of data YTAB models are run on.
- C. albicans silver: the S. cerevisiae call through 1:1 orthologs (CGD ortholog table) when there is no direct call.

  python scripts/local/ytab_build_species_resources.py --species cerevisiae --sgd-release S288C_reference_genome_Current_Release.tgz \\
      --sgd-phenotypes phenotype_data.tab.gz
  python scripts/local/ytab_build_species_resources.py --species albicans --cgd-phenotypes C_albicans_SC5314_phenotype_data.tab \\
      --sgd-phenotypes phenotype_data.tab.gz --sgd-gff saccharomyces_cerevisiae_R64-5-1_20240529.gff
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import date
import gzip
import json
from pathlib import Path
import re
import shutil
import sys
import tarfile
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality.control_conditioned_gumbel import annotations  # noqa: E402
from ytab.normalization.ttr import sha256  # noqa: E402

SPECIES = ROOT/"resources/species"
TRANSPOSON_PMIDS = {"30377286"}
URLS = dict(sgd_release="http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz",
            sgd_phenotypes="https://downloads.yeastgenome.org/curation/literature/phenotype_data.tab",
            cgd_phenotypes="http://www.candidagenome.org/download/phenotype/C_albicans_SC5314_phenotype_data.tab",
            cgd_orthologs="CGD C_albicans_SC5314_S_cerevisiae_orthologs.txt (repository copy; header records its date)")
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]


def opener(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def gff_rows(path: Path):
    for line in opener(path):
        if line.startswith("##FASTA"):break
        if line.startswith("#") or not line.strip():continue
        f = line.rstrip("\n").split("\t")
        if len(f) == 9:yield f, {unquote(k):unquote(v) for k, v in (x.split("=", 1) for x in f[8].split(";") if "=" in x)}


def fasta_names(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in open(path) if line.startswith(">")]


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as handle:
        w = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def install_sgd_release(archive: Path, reference: Path) -> tuple[Path, Path, str]:
    """Extract the release GFF; its ##FASTA section (chrI..chrmt, the GFF's own sequence names) becomes the FASTA."""
    with tarfile.open(archive) as tar:
        member = next(m for m in tar.getmembers() if m.isfile() and not Path(m.name).name.startswith("._")
                      and re.search(r"saccharomyces_cerevisiae_R[\d-]+_\d+\.gff\.gz$", m.name))
        release = re.search(r"(R[\d-]+_\d+)", member.name).group(1)
        text = gzip.decompress(tar.extractfile(member).read()).decode()
    annotation, sequence = text.split("##FASTA\n", 1)
    reference.mkdir(parents=True, exist_ok=True)
    gff, fasta = reference/f"saccharomyces_cerevisiae_{release}.gff", reference/f"S288C_reference_{release}.fasta"
    gff.write_text(annotation); fasta.write_text(sequence)
    return fasta, gff, release


def chromosome_display(contig: str) -> tuple[str, int]:
    if re.search(r"(chrmt|chrM)", contig):return "Mito", 99
    m = re.match(r"chr([IVX]+)$", contig)
    if m and m.group(1) in ROMAN:return f"Chr {m.group(1)}", ROMAN.index(m.group(1))+1
    m = re.match(r"Ca\d+chr([0-9R])[AB]?_", contig)
    if m:return f"Chr {m.group(1)}", 8 if m.group(1) == "R" else int(m.group(1))
    return contig, 98


def tracks(gff: Path, contigs: set[str], reference: Path, source: str):
    centromeres, trnas = [], []
    for f, a in gff_rows(gff):
        if f[0] not in contigs:continue
        if f[2] == "centromere":centromeres.append((f[0], int(f[3])-1, int(f[4]), a.get("Name", a.get("ID", ""))))
        if f[2] in ("tRNA", "tRNA_gene") and (f[2] == "tRNA" or not any(r[3] == a.get("ID") for r in trnas)):
            trnas.append((f[0], int(f[3])-1, int(f[4]), a.get("ID", ""), ".", f[6]))
    # SGD lists each tRNA twice (tRNA_gene and its tRNA); keep one per interval.
    trnas = sorted({(c, s, e): (c, s, e, n, d, st) for c, s, e, n, d, st in trnas}.values())
    # SGD has a centromere feature and its centromere child at the same interval; keep one per interval.
    centromeres = sorted({(c, s, e): (c, s, e, n) for c, s, e, n in sorted(centromeres, key=lambda r:len(r[3]), reverse=True)}.values())
    (reference/f"centromeres.{source}.bed").write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in centromeres))
    (reference/"tRNAs.bed").write_text("".join("\t".join(map(str, r))+"\n" for r in trnas))
    names = [dict(contig=c, display=chromosome_display(c)[0], order=chromosome_display(c)[1]) for c in sorted(contigs, key=lambda c:chromosome_display(c)[1])]
    write_tsv(reference/"chromosome_names.tsv", names)
    return len(centromeres), len(trnas)


def sgd_calls(path: Path) -> dict[str, tuple[str, str]]:
    calls = defaultdict(lambda:defaultdict(set))
    for r in csv.reader(opener(path), delimiter="\t"):
        if len(r) < 10 or r[6] != "null" or r[9] not in ("viable", "inviable"):continue
        source = "systematic" if "systematic mutation set" in r[5] else "classical" if "classical genetics" in r[5] else None
        if source:calls[r[0]][source].add(r[9])
    out = {}
    for gene, c in calls.items():
        for source in ("systematic", "classical"):
            if c[source]:
                out[gene] = ("conflict" if len(c[source]) > 1 else "essential" if c[source] == {"inviable"} else "nonessential", source); break
    return out


def sgd_genes(gff: Path) -> dict[str, dict]:
    return {a["ID"]:dict(name=a.get("gene", a["ID"]), description=a.get("display", a.get("Note", "")))
            for f, a in gff_rows(gff) if f[2] == "gene" and a.get("ID")}


def cgd_direct(path: Path):
    direct, comparator = defaultdict(list), defaultdict(set)
    for r in csv.DictReader(opener(path), delimiter="\t"):
        if r["Phenotype"] not in ("viable", "inviable"):continue
        pmid = (re.search(r"PMID: (\d+)", r["Reference"]) or [None, ""])[1]
        if pmid in TRANSPOSON_PMIDS:comparator[r["Feature Name"]].add(r["Phenotype"]); continue
        ploidy, mutant, cond = r["Experiment Type"].lower(), r["Mutant Type"], r["Condition"].strip()
        complete = ploidy.startswith("homozygous diploid") or ploidy.startswith("haploid")
        repressed = "repress" in cond.lower()
        if r["Phenotype"] == "inviable" and mutant == "repressible":e = "essential" if repressed or not cond else "conditional"
        elif r["Phenotype"] == "inviable" and mutant == "null" and complete:e = "conditional" if cond else "essential"
        elif r["Phenotype"] == "viable" and mutant == "null" and complete:e = "nonessential"
        else:continue
        direct[r["Feature Name"]].append(dict(evidence=e, pmid=pmid, condition=cond))
    resolved = {}
    for gene, ev in direct.items():
        kinds = {x["evidence"] for x in ev}
        call = ("conflict" if {"essential", "nonessential"} <= kinds else "essential" if "essential" in kinds else
                "nonessential" if "nonessential" in kinds else "conditional")
        resolved[gene] = (call, ";".join(sorted({x["pmid"] for x in ev if x["pmid"]})), ";".join(sorted({x["condition"] for x in ev if x["condition"]})))
    comparator = {g:("conflict" if len(p) > 1 else "essential" if p == {"inviable"} else "nonessential") for g, p in comparator.items()}
    return resolved, comparator


def cgd_orthologs(path: Path) -> dict[str, tuple[str, list[str], list[str]]]:
    """C. albicans feature -> (relationship, S. cerevisiae systematic IDs, names); 1:1 when each side has one partner."""
    pairs = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#") and l.strip()]
    ca_to_sc, sc_to_ca = defaultdict(set), defaultdict(set)
    names = {}
    for ca, _, _, sc, sc_name, _ in pairs:
        ca_to_sc[ca].add(sc); sc_to_ca[sc].add(ca); names[sc] = sc_name
    out = {}
    for ca, scs in ca_to_sc.items():
        n_ca = max(len(sc_to_ca[s]) for s in scs)
        out[ca] = (f"{n_ca}:{len(scs)}", sorted(scs), [names[s] for s in sorted(scs)])
    return out


def readme(species: str, summary: dict) -> str:
    gold, silver = summary["labels"]["gold"], summary["labels"].get("silver", {"essential":0, "nonessential":0})
    if species == "cerevisiae":
        rows = f"| gold | S. cerevisiae null-mutant viability (SGD) | {gold['essential']:,} | {gold['nonessential']:,} |"
    else:
        rows = (f"| gold | C. albicans null (homozygous diploid or haploid) or repressible allele (CGD), excluding transposon surveys | {gold['essential']:,} | {gold['nonessential']:,} |\n"
                f"| silver | S. cerevisiae null-mutant viability, 1:1 ortholog | {silver['essential']:,} | {silver['nonessential']:,} |")
    note = ("\n- **`transposon_comparator`:** transposon mutagenesis and machine-learning calls in a haploid strain as recorded in "
            "CGD (PMID 30377286). Never labels: they come from the same kind of data YTAB models are run on.") if species == "albicans" else ""
    return (f"# Curated {'S. cerevisiae' if species == 'cerevisiae' else 'C. albicans'} essentiality labels\n\n"
            "Built by `scripts/local/ytab_build_species_resources.py`; rules are in its docstring, sources and SHA-256 hashes in "
            "`sources.json`, counts in `summary.json`.\n\n| tier | evidence | essential | nonessential |\n|---|---|---|---|\n"
            f"{rows}\n\n- **Unlabelled ({sum(summary['unlabelled'].values()):,} genes):** reasons in `summary.json`.{note}\n")


def build(species: str, args) -> dict:
    base = SPECIES/species; reference = base/"reference_genome"; raw = base/"essentiality_labels/raw"
    raw.mkdir(parents=True, exist_ok=True)
    sources = dict(built=str(date.today()), script_sha256=sha256(Path(__file__)), urls={}, files={})
    if species == "cerevisiae":
        fasta, gff, release = install_sgd_release(args.sgd_release, reference)
        sources["urls"]["reference"] = URLS["sgd_release"]; sources["release"] = release
        sources["files"][args.sgd_release.name] = sha256(args.sgd_release)
        sgd_gff = gff
    else:
        from ytab.pipeline.reference_registry import resolve_reference
        info = resolve_reference(species, ROOT); fasta, gff = info.fasta, info.gff; sgd_gff = args.sgd_gff
    contigs = set(fasta_names(fasta))
    n_cen, n_trna = tracks(gff, contigs, reference, "SGD" if species == "cerevisiae" else "CGD")
    shutil.copyfile(args.sgd_phenotypes, raw/"SGD_phenotype_data.tab.gz") if str(args.sgd_phenotypes).endswith(".gz") else \
        (raw/"SGD_phenotype_data.tab.gz").write_bytes(gzip.compress(Path(args.sgd_phenotypes).read_bytes()))
    sources["urls"]["sgd_phenotypes"] = URLS["sgd_phenotypes"]; sources["files"]["SGD_phenotype_data.tab.gz"] = sha256(raw/"SGD_phenotype_data.tab.gz")
    sc = sgd_calls(raw/"SGD_phenotype_data.tab.gz"); scg = sgd_genes(sgd_gff)
    genes = [g for g in annotations(gff) if g["chrom"] in contigs]
    rows, annotation = [], []
    if species == "cerevisiae":
        for g in genes:
            call, src = sc.get(g["id"], ("", ""))
            label, tier, reason = (call, "gold", "") if call in ("essential", "nonessential") else ("", "", "S. cerevisiae conflict" if call == "conflict" else "no viability data")
            name = scg.get(g["id"], {}).get("name", g["name"])
            rows.append(dict(gene_id=g["id"], name=name, chromosome=g["chrom"], start=g["blocks"][0][0], end=g["blocks"][-1][1],
                             sc_call=call, sc_source=src, label=label, tier=tier, unlabelled_reason=reason))
            annotation.append(dict(gene_id=g["id"], systematic_id=g["id"], gene_name=name, scer_relationship="self", scer_gene_id=g["id"],
                                   scer_gene_name=name, scer_essentiality=call, description=scg.get(g["id"], {}).get("description", g["description"])))
    else:
        shutil.copyfile(args.cgd_phenotypes, raw/"C_albicans_SC5314_phenotype_data.tab.gz") if str(args.cgd_phenotypes).endswith(".gz") else \
            (raw/"C_albicans_SC5314_phenotype_data.tab.gz").write_bytes(gzip.compress(Path(args.cgd_phenotypes).read_bytes()))
        ortho_path = base/"C_albicans_SC5314_S_cerevisiae_orthologs.txt"
        sources["urls"] |= dict(cgd_phenotypes=URLS["cgd_phenotypes"], cgd_orthologs=URLS["cgd_orthologs"])
        sources["files"] |= {"C_albicans_SC5314_phenotype_data.tab.gz":sha256(raw/"C_albicans_SC5314_phenotype_data.tab.gz"), ortho_path.name:sha256(ortho_path),
                             gff.name:sha256(gff), Path(sgd_gff).name:sha256(Path(sgd_gff))}
        direct, comparator = cgd_direct(args.cgd_phenotypes); ortho = cgd_orthologs(ortho_path)
        for g in genes:
            rel, scs, sc_names = ortho.get(g["id"], ("no ortholog", [], []))
            one = rel == "1:1"
            sc_call, sc_src = sc.get(scs[0], ("", "")) if one else ("", "")
            d_call, pmids, conds = direct.get(g["id"], ("", "", ""))
            if d_call in ("essential", "nonessential"):label, tier, reason = d_call, "gold", ""
            elif d_call in ("conflict", "conditional"):label, tier, reason = "", "", f"C. albicans direct evidence {d_call}"
            elif not one:label, tier, reason = "", "", f"orthology {rel}"
            elif sc_call in ("essential", "nonessential"):label, tier, reason = sc_call, "silver", ""
            else:label, tier, reason = "", "", "S. cerevisiae conflict" if sc_call == "conflict" else "S. cerevisiae ortholog without viability data"
            rows.append(dict(gene_id=g["id"], name=g["name"], chromosome=g["chrom"], start=g["blocks"][0][0], end=g["blocks"][-1][1],
                             orthology=rel, scer_id=";".join(scs), scer_name=";".join(sc_names), sc_call=sc_call, sc_source=sc_src,
                             direct_call=d_call, direct_pmids=pmids, direct_conditions=conds, label=label, tier=tier, unlabelled_reason=reason,
                             transposon_comparator=comparator.get(g["id"], "")))
            annotation.append(dict(gene_id=g["id"], systematic_id=re.sub(r"_[AB]$", "", g["id"]), gene_name=g["name"] if g["name"] != g["id"] else "",
                                   scer_relationship=rel if scs else "", scer_gene_id=";".join(scs), scer_gene_name=";".join(sc_names),
                                   scer_essentiality=";".join(sc.get(s, ("", ""))[0] for s in scs) if scs else "",
                                   description=g["description"] or ";".join(scg.get(s, {}).get("description", "") for s in scs)))
    write_tsv(base/"essentiality_labels/labels.tsv", rows)
    write_tsv(base/"annotation/gene_annotation.tsv", annotation)
    count = lambda pred:sum(1 for r in rows if pred(r))
    summary = dict(genes=len(rows), centromeres=n_cen, trnas=n_trna, chromosomes=len(contigs),
                   labels={t:{l:count(lambda r, t=t, l=l:r["tier"] == t and r["label"] == l) for l in ("essential", "nonessential")} for t in ("gold", "silver")},
                   unlabelled=dict(sorted({x:count(lambda r, x=x:not r["label"] and r["unlabelled_reason"] == x)
                                           for x in {r["unlabelled_reason"] for r in rows if not r["label"]}}.items(), key=lambda kv:-kv[1])))
    if species == "albicans":
        both = [r for r in rows if r["direct_call"] in ("essential", "nonessential") and r["sc_call"] in ("essential", "nonessential")]
        labelled = [r for r in rows if r["label"] and r["transposon_comparator"] in ("essential", "nonessential")]
        summary |= dict(gold_vs_sc_ortholog=dict(genes=len(both), agree=sum(r["direct_call"] == r["sc_call"] for r in both)),
                        comparator_only_transposon=dict(labelled_genes_with_call=len(labelled), agree=sum(r["label"] == r["transposon_comparator"] for r in labelled)))
    (base/"essentiality_labels/summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    (base/"essentiality_labels/sources.json").write_text(json.dumps(sources, indent=2)+"\n")
    (base/"essentiality_labels/README.md").write_text(readme(species, summary))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, choices=["cerevisiae", "albicans"])
    p.add_argument("--sgd-release", type=Path); p.add_argument("--sgd-phenotypes", type=Path, required=True)
    p.add_argument("--cgd-phenotypes", type=Path); p.add_argument("--sgd-gff", type=Path)
    a = p.parse_args()
    if a.species == "cerevisiae" and not a.sgd_release:p.error("--sgd-release is required for cerevisiae")
    if a.species == "albicans" and not (a.cgd_phenotypes and a.sgd_gff):p.error("--cgd-phenotypes and --sgd-gff are required for albicans")
    print(json.dumps(build(a.species, a), indent=1))


if __name__ == "__main__":
    main()
