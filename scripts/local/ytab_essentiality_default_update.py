#!/usr/bin/env python3
"""Rescore the paper's essentiality fits with the current HMM gene-scoring default and re-evaluate them.

After R20/R21 (output/essentiality/hmm_proportional_run, hmm_run_density) the default essential run became at least
20% of a gene's mappable windows (floor 3). The fits do not change with the run rule, so this copies each fit used in
the paper (fit.json, windows.npz, target_model.tsv) to output/essentiality/default_update/hmm/<name>/joint, rewrites
its genes.tsv with `ytab hermes-hmm --rescore`, and repeats the four-species evaluation (same labels, masks, gene
filters and Segal et al.'s RF scores, seed 81723) there. The pre-registered folders are not touched, so their
registered results can still be regenerated.

  python scripts/local/ytab_essentiality_default_update.py rescore
  python scripts/local/ytab_essentiality_default_update.py evaluate   # writes evaluation.tsv, result.json
"""
import csv
import json
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
import ytab_four_species_benchmark as fsb
import ytab_hermes_hmm as runner
from ytab.essentiality import hermes_hmm as hh

OUT = ROOT/"output/essentiality/default_update"
SOURCES = {name:fsb.BENCH/"hmm"/name/"joint" for name in fsb.FITS} | {
    "glabrata_gale":ROOT/"output/essentiality/hermes_hmm/benchmark_v3/gale_blinded/joint_run1"}


def rescore():
    log = []
    for name, src in SOURCES.items():
        dst = OUT/"hmm"/name/"joint"; dst.mkdir(parents=True, exist_ok=True)
        for f in ("fit.json", "windows.npz", "target_model.tsv"):
            if (src/f).exists() and not (dst/f).exists():shutil.copy2(src/f, dst/f)
        info = runner.rescore(dst)
        log.append(dict(name=name, source=str(src.relative_to(ROOT)), converged=info["converged"], **info["rescored"], gene_calls=info["gene_calls"]))
        print(name, info["gene_calls"], f"max posterior difference {info['rescored']['max_posterior_difference']:.1e}", flush=True)
    (OUT/"rescore.json").write_text(json.dumps(dict(settings=hh._json_ready(vars(hh.Settings())), fits=log), indent=1))


def evaluate():
    """fsb.evaluate() on the rescored gene tables, writing here; plus C. glabrata on S. cerevisiae ortholog labels (silver)."""
    fits = json.loads((fsb.BENCH/"fits.json").read_text())
    (OUT/"fits.json").write_text(json.dumps(fits, indent=1))
    fsb.BENCH = OUT
    fsb.check_preregistration = lambda: dict(note="default update, not the registered analysis; registered result: output/essentiality/four_species/result.json")
    fsb.evaluate()
    import ytab_hmm_sensitivity as sens
    lab = sens.labels("glabrata"); genes = sens.genes(OUT/"hmm/glabrata_gale/joint/genes.tsv")
    common = sorted(set(lab) & set(genes)); y = np.array([lab[g] for g in common]); s = np.array([genes[g] for g in common])
    result = json.loads((OUT/"result.json").read_text())
    calls = {r["ORF"].removeprefix("gene:"):r["call"] for r in csv.DictReader(open(OUT/"hmm/glabrata_gale/joint/genes.tsv"), delimiter="\t")}
    result["glabrata_silver"] = dict(genes=len(common), essential=int(y.sum()), auc_hmm=fsb.auc(s, y),
                                     calls_hmm=fsb.calls([calls[g] == "ES" for g in common], y))
    (OUT/"result.json").write_text(json.dumps(result, indent=1)); print(json.dumps(result["glabrata_silver"], indent=1))


if __name__ == "__main__":
    {"rescore":rescore, "evaluate":evaluate}[sys.argv[1]]()
