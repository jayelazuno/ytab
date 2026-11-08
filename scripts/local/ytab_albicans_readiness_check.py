#!/usr/bin/env python3
"""C. albicans readiness check: run the species-agnostic Hermes tools on the real C. albicans SC5314 reference.

No C. albicans Hermes data exist in the repository, so two Hermes libraries are simulated on the real
haploid reference (A22 HapA FASTA) with the real A22 annotation and a known set of essential genes:
- genes: haplotype-A genes on HapA chromosomes; 15% drawn as essential (seed 31);
- insertions: 250,000 distinct sites per library, uniform over the genome, 2x enriched outside CDS, none in
  essential CDS; reads lognormal. Libraries are independent draws (independent pools).
Then, with the real commands:
1. ytab hermes-hmm on HapA + A22 GFF + the two hit files -> recovery of the simulated essential genes (AUC);
2. the same with the full diploid FASTA (haplotypes A and B) -> must stop with the unmappable-reference error;
3. ytab library-correlation of the two C. albicans libraries and one C. glabrata Gale library through the
   four-species 1:1 orthologs.
Writes output/essentiality/species_readiness/albicans/. This checks that the tools run correctly on
C. albicans inputs; it is not a validation of performance on real C. albicans data.
"""
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import hermes_hmm as hh

REF = ROOT/"resources/species/albicans/reference_genome"
HAPA = REF/"C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA.fasta"
DIPLOID = REF/"C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes.fasta"
GFF = REF/"C_albicans_SC5314_version_A22-s07-m01-r08_features.gff"
OUT = ROOT/"output/essentiality/species_readiness/albicans"
SEED, SITES, ESSENTIAL = 31, 250_000, .15


def auc(score, label):
    score, label = np.asarray(score, float), np.asarray(label, bool)
    order = np.argsort(score, kind="mergesort"); ranks = np.empty(len(score)); ranks[order] = np.arange(1, len(score)+1)
    _, inv, counts = np.unique(score, return_inverse=True, return_counts=True); ranks = (np.bincount(inv, ranks)/counts)[inv]
    n1, n0 = label.sum(), (~label).sum(); return float((ranks[label].sum()-n1*(n1+1)/2)/(n1*n0))


def simulate():
    rng = np.random.default_rng(SEED)
    seqs = hh.read_fasta(HAPA); lengths = {c:len(s) for c, s in seqs.items()}
    genes = [g for g in ccg.annotations(GFF) if g["chrom"] in lengths]
    essential = {g["id"] for g in genes if rng.random() < ESSENTIAL}
    files = []
    for lib in ("Ca_sim_pool1", "Ca_sim_pool2"):
        rows = []
        for c, n in lengths.items():
            w = np.ones(n); cds = np.zeros(n, dtype=bool)
            for g in genes:
                if g["chrom"] != c:continue
                for a, b in g["blocks"]:cds[a-1:b] = True
            w[~cds] = 2.
            for g in genes:
                if g["chrom"] == c and g["id"] in essential:
                    for a, b in g["blocks"]:w[a-1:b] = 0.
            k = int(round(SITES*n/sum(lengths.values())))
            pos = np.sort(rng.choice(n, k, replace=False, p=w/w.sum()))+1
            reads = np.maximum(1, np.round(np.exp(rng.normal(2.5, 1.2, k)))).astype(int)
            rows += [(c, p, r) for p, r in zip(pos, reads)]
        path = OUT/"hits"/f"{lib}_hits.txt"; path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as h:
            h.write("Chromosome\tPosition\tReads\n"); h.writelines(f"{c}\t{p}\t{r}\n" for c, p, r in rows)
        files.append(path)
    (OUT/"truth.tsv").write_text("ORF\tessential\n"+"".join(f"{g['id']}\t{g['id'] in essential}\n" for g in genes))
    return files, essential


def ytab(*args):
    return subprocess.run([sys.executable, str(ROOT/"scripts/local/ytab.py"), *map(str, args)], capture_output=True, text=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True); result = {}
    files, essential = simulate()
    hmm = ytab("hermes-hmm", OUT/"hmm_hapA", "--fasta", HAPA, "--gff", GFF, *[x for f in files for x in ("--hits", f)])
    result["hmm_hapA_returncode"] = hmm.returncode; result["hmm_hapA_notes"] = [l for l in hmm.stderr.splitlines() if l.startswith("NOTE")]
    if hmm.returncode == 0:
        fit = json.loads((OUT/"hmm_hapA/joint/fit.json").read_text())
        rows = [r for r in csv.DictReader((OUT/"hmm_hapA/joint/genes.tsv").open(), delimiter="\t") if r["status"] == "called"]
        y = np.array([r["ORF"] in essential for r in rows]); s = np.array([float(r["P_ES_run"]) for r in rows])
        called = np.array([r["call"] == "ES" for r in rows])
        result["hmm_hapA"] = dict(masked_fraction=fit["masked_fraction"], converged=fit["converged"], iterations=fit["iterations"], r=fit["r"],
                                  genes_called=len(rows), simulated_essential_called=int(y.sum()), auc=auc(s, y),
                                  es_precision=float(y[called].mean()) if called.any() else None, es_recall=float(called[y].mean()), notes=fit["notes"])
    else:result["hmm_hapA_error"] = hmm.stderr[-2000:]
    diploid = ytab("hermes-hmm", OUT/"hmm_diploid", "--fasta", DIPLOID, "--gff", GFF, *[x for f in files for x in ("--hits", f)])
    result["diploid_reference_guard"] = dict(returncode=diploid.returncode, message=diploid.stderr.strip()[-600:],
                                             stopped_as_intended=diploid.returncode != 0 and "unmappable" in diploid.stderr)
    table = OUT/"libraries.tsv"
    with table.open("w") as h:
        h.write("name\thits\tspecies\tgff\tstrain\tpool\tcondition\n")
        for f in files:h.write(f"albicans::{f.stem.removesuffix('_hits')}\t{f}\talbicans\t{GFF}\tsimulated\t{f.stem.split('_')[2]}\tparent\n")
        h.write(f"glabrata::Gale_Cg1_parent\t{ROOT}/output/projects/Gale_2020/create_hit_file/SRR11561508/SRR11561508_hits.txt\tglabrata\t"
                f"{ROOT}/resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff\tBG14\tCg1\tparent\n")
    corr = ytab("library-correlation", "--samples", table, "--orthologs", ROOT/"resources/species/orthologs/glabrata_cerevisiae_pombe_albicans_1to1.tsv",
                "--bootstrap", "100", "--output", OUT/"library_correlation")
    result["library_correlation_returncode"] = corr.returncode
    if corr.returncode == 0:
        result["library_correlation_pairs"] = [{k:r[k] for k in ("sample_a", "sample_b", "pair_type", "genes", "rho", "shared_site_read_rho", "clone_evidence")}
                                               for r in csv.DictReader((OUT/"library_correlation/pairs.tsv").open(), delimiter="\t")]
    else:result["library_correlation_error"] = corr.stderr[-2000:]
    (OUT/"result.json").write_text(json.dumps(result, indent=2)+"\n"); print(json.dumps(result, indent=1))


if __name__ == "__main__":main()
