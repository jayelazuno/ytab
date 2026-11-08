#!/usr/bin/env python3
"""aBGC arm of the R17 condition comparisons (output/conditions/ADDENDUM_abgc.md; descriptive, not pre-registered).

Same tests, libraries and endpoints as scripts/local/ytab_conditions_benchmark.py and the S. pombe between-strain
addendum, with aBGC normalization runs of their own.
  run:       normalization runs and the four tests -> output/conditions/<test>/aBGC/
  evaluate:  P1 recovery, P2 SGD rapamycin AUC (aBGC - BGC), P3 heterochromatin, between-strain recovery
             -> output/conditions/abgc_addendum/
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("r17", ROOT/"scripts/local/ytab_conditions_benchmark.py")
r17 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(r17)

ARM = "aBGC"
OUT = r17.BENCH/"abgc_addendum"
SATAY_NORM = ROOT/"output/projects"/r17.SATAY/"normalization/conditions_abgc"
POMBE_NORM = ROOT/"output/projects"/r17.POMBE/"normalization/conditions_abgc"
BETWEEN_NORM = ROOT/"output/projects"/r17.POMBE/"normalization/between_strain_abgc"
BETWEEN = {"SRR8502436":"WT", "SRR8502435":"WT", "SRR8502434":"otr1R", "SRR8502433":"otr1R"}


def run():
    S, step, config = r17.SCRIPTS, r17.step, r17.config
    normalize = lambda name, project, out, samples, together:step(name, [S/"ytab_normalize_project.py", "--project-config", config(project), "--output", out,
        "--method", ARM, *[x for s in samples for x in ("--sample", s)], *(["--normalize-together"] if together else [])])
    if not (SATAY_NORM/"manifest.json").exists():normalize("normalize_satay_abgc", r17.SATAY, SATAY_NORM, [*r17.SATAY_CONDITIONS, "ERR1533149"], True)
    if not (POMBE_NORM/"manifest.json").exists():normalize("normalize_pombe_abgc", r17.POMBE, POMBE_NORM, [s for s, _ in r17.POMBE_GI], False)
    if not (BETWEEN_NORM/"manifest.json").exists():normalize("normalize_pombe_between_abgc", r17.POMBE, BETWEEN_NORM, list(BETWEEN), True)
    out = r17.BENCH/"rapamycin"/ARM
    if not (out/"genes.tsv").exists():
        step(f"rapamycin_{ARM}", [S/"ytab_treated_vs_control.py", "--project-config", config(r17.SATAY), "--method", "resampling",
                                  "--normalization-run", SATAY_NORM, "--normalization", ARM, "--comparisons", r17.RAPAMYCIN, "--output", out])
    out = r17.BENCH/"satay_mutants"/ARM
    if not out.exists():
        step(f"satay_mutants_{ARM}", [S/"ytab_multicondition_project.py", "--project-config", config(r17.SATAY), "--normalized", SATAY_NORM/"experiment/methods",
                                      "--normalization", ARM, "--reference-condition", "WT", "--independent-libraries",
                                      *[x for s, c in r17.SATAY_CONDITIONS.items() for x in ("--sample", f"{r17.SATAY}::{s}={c}")], "--output", out])
    out = r17.BENCH/"pombe_gi"/ARM
    if not out.exists():
        fasta, gff = r17.reference(r17.POMBE)
        step(f"pombe_gi_{ARM}", [S/"ytab_genetic_interactions.py", "--normalized", POMBE_NORM/"experiment/methods", "--normalization", ARM,
                                 "--gff", gff, "--reference", fasta, "--output", out, "--reference-strain", "WT", "--query-strain", "otr1R",
                                 *[x for s, role in r17.POMBE_GI for x in ("--sample", f"{r17.POMBE}::{s}={role}")]])
    out = r17.BENCH/"pombe_between_strain"/ARM
    if not out.exists():
        step(f"pombe_between_{ARM}", [S/"ytab_multicondition_project.py", "--project-config", config(r17.POMBE), "--normalized", BETWEEN_NORM/"experiment/methods",
                                      "--normalization", ARM, "--reference-condition", "WT", "--independent-libraries",
                                      *[x for s, c in BETWEEN.items() for x in ("--sample", f"{r17.POMBE}::{s}={c}")], "--output", out])
    print("done")


def between(norm):
    """Between-strain endpoints: named-hit recovery and canonical heterochromatin genes among calls depleted in otr1R."""
    import pandas as pd
    from scipy.stats import fisher_exact
    d = pd.read_csv(r17.BENCH/"pombe_between_strain"/norm/"tests/pairs_genes.tsv", sep="\t")
    d = d[(d.status == "tested") & ~d.Name.isin(r17.MASKED)]
    dep = d.call_all_pairs == "depleted"; canon = d.Name.str.lower().isin(r17.CANONICAL)
    named = r17.POSITIVES[-1][3].split()
    table = [[int((canon & dep).sum()), int((canon & ~dep).sum())], [int((~canon & dep).sum()), int((~canon & ~dep).sum())]]
    odds, p = fisher_exact(table, alternative="greater")
    return dict(normalization=norm, tested=len(d), depleted=int(dep.sum()), named_recovered=int(d[dep].Name.str.lower().isin([g.lower() for g in named]).sum()),
                named=len(named), canonical_depleted=table[0][0], canonical_tested=int(canon.sum()), odds_ratio=float(odds), p=float(p))


def evaluate():
    import pandas as pd
    sys.path.append(str(r17.SCRIPTS)); OUT.mkdir(parents=True, exist_ok=True)
    r17.NORMS = ["betageom", ARM, "totreads", "ttr"]   # BGC stays the bootstrap reference of rapamycin_auc(): delta_bgc_vs_aBGC = BGC - aBGC
    rec = pd.DataFrame([r for n in r17.NORMS for r in r17.recovery(n)]); rec.to_csv(OUT/"recovery.tsv", sep="\t", index=False)
    by_test = rec.groupby(["normalization", "test"], sort=False).agg(positives=("gene", "size"), called=("called", "sum")).reset_index()
    by_test.to_csv(OUT/"recovery_by_test.tsv", sep="\t", index=False)
    sets = r17.sgd_rapamycin(); auc = pd.DataFrame(r17.rapamycin_auc(sets)); auc.to_csv(OUT/"rapamycin_sgd_auc.tsv", sep="\t", index=False)
    het = [r17.heterochromatin(n) for n in r17.NORMS]
    pd.DataFrame(het[1]["candidates"]).to_csv(OUT/f"pombe_candidates_{ARM}.tsv", sep="\t", index=False)
    bs = pd.DataFrame([between(n) for n in ("betageom", ARM)]); bs.to_csv(OUT/"between_strain.tsv", sep="\t", index=False)
    counts = {n:dict(rapamycin=r17.tables(n)[0].call.value_counts().to_dict(),
                     satay_mutants={p:g.call.value_counts().to_dict() for p, g in r17.tables(n)[1].groupby("pair")},
                     pombe_gi=r17.tables(n)[2].call.value_counts().to_dict()) for n in r17.NORMS}
    result = dict(addendum="output/conditions/ADDENDUM_abgc.md", recovery=by_test.to_dict("records"), rapamycin_sgd_auc=auc.to_dict("records"),
                  heterochromatin=[{k:v for k, v in h.items() if k != "candidates"} | dict(candidates=len(h["candidates"])) for h in het],
                  between_strain=bs.to_dict("records"), call_counts=counts)
    (OUT/"result.json").write_text(json.dumps(result, indent=1, default=str))
    print(by_test.to_string()); print(auc.to_string()); print(bs.to_string())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("phase", choices=["run", "evaluate"])
    a = p.parse_args(); run() if a.phase == "run" else evaluate()
