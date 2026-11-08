#!/usr/bin/env python3
"""Build a transposon donor reference from its parts (resources/transposons/donors/<donor>/recipe.yaml -> donor.fasta).

Most donor plasmids are not deposited; labs work from files pieced together from known elements. The reads that matter
start at a transposon end and run 20-150 bp into whatever flanks it in the donor, so the reference needs every part
that can flank an end (vector backbone and MCS, markers, promoters, terminators, host-genome fragments), not the exact
junctions. Each part becomes one contig of donor.fasta; the order in the recipe follows the published construct.

Recipe items (each one contig):
  - {name: pRS416 backbone, genbank: U03450}                      a GenBank record (components/<accession>.gb, fetched once)
  - {name: Hermes, fasta: resources/transposons/hermes/L34807.1.fasta}
  - {name: ADE2, species: cerevisiae, gene: YOR128C, flank: 600}   a host gene with flanks, from the species' reference
  - {name: GAL1-10 promoter, species: cerevisiae, region: chrII:278352-279021}
Validation (`--compare FILE.gbk|fasta`): the real plasmid is cut into 100-bp windows, each aligned (bowtie2, local) to
the built reference; reported: windows covered, and windows within 150 bp outside each transposon end.

  python scripts/local/ytab_build_donor.py resources/transposons/donors/<donor> [--compare real.gbk]
"""
from __future__ import annotations
import argparse, json, re, subprocess, sys, tempfile
from pathlib import Path
import yaml
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]; COMPONENTS = ROOT/"resources/transposons/donors/components"
GENOMES = {"cerevisiae":"resources/species/cerevisiae/reference_genome", "glabrata":"resources/species/glabrata/reference_genome",
           "glabrata_CBS138":"resources/species/glabrata/strains/CBS138/reference_genome", "pombe":"resources/species/pombe/reference_genome",
           "albicans":"resources/species/albicans/reference_genome", "auris":"resources/species/auris/reference_genome",
           "neoformans":"resources/species/neoformans/reference_genome", "phaffii":"resources/species/phaffii/reference_genome",
           "lipolytica":"resources/species/lipolytica/reference_genome"}


def genbank(acc):
    path = COMPONENTS/f"{acc}.gb"
    if not path.exists():
        subprocess.run(["curl", "-s", "-G", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", "--data-urlencode", "db=nuccore",
                        "--data-urlencode", f"id={acc}", "--data-urlencode", "rettype=gb", "--data-urlencode", "retmode=text", "-o", str(path)], check=True)
    return str(SeqIO.read(path, "genbank").seq).upper()


def genome_files(species):
    d = ROOT/GENOMES[species]
    fasta = next(p for p in sorted(d.iterdir()) if p.suffix in (".fa", ".fna", ".fasta") and not p.name.endswith(".fai"))
    gff = next((p for p in sorted(d.iterdir()) if p.suffix in (".gff", ".gff3")), None)
    return fasta, gff


_seqs = {}
def chrom(species, name):
    if species not in _seqs:_seqs[species] = {r.id:str(r.seq).upper() for r in SeqIO.parse(genome_files(species)[0], "fasta")}
    return _seqs[species][name]


def gene_region(species, gene, flank):
    _, gff = genome_files(species)
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":continue
        attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        ids = {attrs.get(k, "") for k in ("ID", "Name", "gene", "locus_tag", "Gene", "Alias")} | {attrs.get("ID", "").removeprefix("gene:").removeprefix("gene-")}
        if gene in ids or any(gene in v.split(",") for v in ids):
            seq = chrom(species, f[0]); start, end = max(0, int(f[3])-1-flank), int(f[4])+flank
            return seq[start:end], f"{f[0]}:{start+1}-{end}"
    raise SystemExit(f"{gene} not found in {gff}")


def build(folder: Path):
    recipe = yaml.safe_load(open(folder/"recipe.yaml")); contigs, log = [], []
    for part in recipe["parts"]:
        if "genbank" in part:seq, where = genbank(part["genbank"]), f"GenBank {part['genbank']}"
        elif "fasta" in part:seq, where = str(next(SeqIO.parse(ROOT/part["fasta"], "fasta")).seq).upper(), part["fasta"]
        elif "gene" in part:seq, where = gene_region(part["species"], part["gene"], int(part.get("flank", 0))); where = f"{part['species']} {part['gene']} {where}"
        elif "region" in part:
            c, span = part["region"].split(":"); a, b = map(int, span.split("-")); seq, where = chrom(part["species"], c)[a-1:b], f"{part['species']} {part['region']}"
        else:raise SystemExit(f"unknown part {part}")
        contigs.append((re.sub(r"\W+", "_", part["name"]), seq)); log.append(dict(name=part["name"], source=where, length=len(seq), role=part.get("role", "")))
    (folder/"donor.fasta").write_text("".join(f">{n}\n{s}\n" for n, s in contigs))
    (folder/"parts.json").write_text(json.dumps(dict(donor=recipe.get("donor", folder.name), parts=log), indent=1))
    print(f"{folder/'donor.fasta'}: {len(contigs)} parts, {sum(len(s) for _, s in contigs):,} bp")
    return folder/"donor.fasta"


def compare(reference: Path, real: Path, ends: list[str]):
    rec = SeqIO.read(real, "genbank" if real.suffix in (".gb", ".gbk") else "fasta"); s = str(rec.seq).upper(); L = len(s); ss = s+s
    rc = lambda x:x[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    end_pos = sorted(p for e in ends for q in (e, rc(e)) for p in [ss.find(q)] if 0 <= p < L)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["bowtie2-build", "-q", str(reference), f"{tmp}/ref"], check=True, stdout=subprocess.DEVNULL)
        fq = Path(tmp)/"w.fq"; fq.write_text("".join(f"@w{i}\n{ss[i:i+100]}\n+\n{'I'*100}\n" for i in range(0, L, 10)))
        out = subprocess.run(f"bowtie2 --local -x {tmp}/ref -U {fq} --no-unal --no-hd -p 4", shell=True, capture_output=True, text=True).stdout
    hit = {int(l.split("\t")[0][1:]) for l in out.splitlines()}
    windows = list(range(0, L, 10))
    near = [w for w in windows if any(0 < (w-e) % L < 150 or 0 < (e-w-100) % L < 150 for e in end_pos)]
    result = dict(real=real.name, length=L, windows=len(windows), covered=round(len(hit)/len(windows), 4),
                  end_positions=[e+1 for e in end_pos], flank_windows=len(near), flank_covered=round(sum(w in hit for w in near)/max(len(near), 1), 4))
    print(json.dumps(result, indent=1)); return result


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("folder", type=Path); a.add_argument("--compare", type=Path); a = a.parse_args()
    ref = build(a.folder)
    if a.compare:
        ends = yaml.safe_load(open(a.folder/"recipe.yaml")).get("transposon_ends", [])
        (a.folder/"validation.json").write_text(json.dumps(compare(ref, a.compare, ends), indent=1))
