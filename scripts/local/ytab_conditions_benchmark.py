#!/usr/bin/env python3
"""Condition comparisons on public Ac/Ds and Hermes data (rapamycin, SATAY mutants, S. pombe heterochromatin).

Implements output/conditions/PREREGISTRATION.md (frozen; both phases refuse to run if it changed).
  run:       normalization runs and the tests, with ordinary YTAB commands (reads no published positives)
  evaluate:  recovery of the published positives, the SGD rapamycin AUCs, the S. pombe heterochromatin enrichment
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.ttr import sha256

BENCH = ROOT/"output/conditions"
SCRIPTS = ROOT/"scripts/local"
NORMS = ["betageom", "totreads", "ttr"]
SATAY, POMBE = "Michel_2017_SATAY", "Pombe_Hermes"
SATAY_CONDITIONS = {"ERR1533147":"WT", "ERR1533148":"WT", "ERR1533143":"dpl1del", "ERR1533144":"dpl1del_psd2del",
                    "ERR1533146":"VPS13_D716H", "ERR1533145":"mmm1del_VPS13_D716H"}
RAPAMYCIN = "BY4741_ade2_poolWT2_treated_vs_parent"
POMBE_GI = [("SRR8502440", "WT:control"), ("SRR8502436", "WT:treated"), ("SRR8502439", "WT:control"), ("SRR8502435", "WT:treated"),
            ("SRR8502438", "otr1R:control"), ("SRR8502434", "otr1R:treated"), ("SRR8502437", "otr1R:control"), ("SRR8502433", "otr1R:treated")]


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def config(project):
    return ROOT/"output/projects"/project/"config/project.yaml"


def step(name, cmd):
    log = BENCH/"logs"/f"{name}.log"; log.parent.mkdir(parents=True, exist_ok=True)
    print(name, flush=True)
    with open(log, "w") as handle:
        r = subprocess.run([sys.executable]+[str(c) for c in cmd], stdout=handle, stderr=subprocess.STDOUT, cwd=ROOT)
    if r.returncode:raise SystemExit(f"{name} failed; see {log}")


def reference(project):
    import yaml
    cfg = yaml.safe_load(config(project).read_text())["reference"]
    return ROOT/cfg["fasta"], ROOT/cfg["gff"]


def run():
    prereg = check_preregistration(); done = {}
    satay_norm = ROOT/"output/projects"/SATAY/"normalization/conditions_v1"
    pombe_norm = ROOT/"output/projects"/POMBE/"normalization/conditions_v1"
    samples = [*SATAY_CONDITIONS, "ERR1533149"]
    if not (satay_norm/"manifest.json").exists():
        step("normalize_satay", [SCRIPTS/"ytab_normalize_project.py", "--project-config", config(SATAY), "--output", satay_norm,
                                 *[x for m in NORMS for x in ("--method", m)], *[x for s in samples for x in ("--sample", s)], "--normalize-together"])
    if not (pombe_norm/"manifest.json").exists():
        step("normalize_pombe", [SCRIPTS/"ytab_normalize_project.py", "--project-config", config(POMBE), "--output", pombe_norm,
                                 *[x for m in NORMS for x in ("--method", m)], *[x for s, _ in POMBE_GI for x in ("--sample", s)]])
    for norm in NORMS:
        out = BENCH/"rapamycin"/norm
        if not (out/"genes.tsv").exists():
            step(f"rapamycin_{norm}", [SCRIPTS/"ytab_treated_vs_control.py", "--project-config", config(SATAY), "--method", "resampling",
                                       "--normalization-run", satay_norm, "--normalization", norm, "--comparisons", RAPAMYCIN, "--output", out])
        out = BENCH/"satay_mutants"/norm
        if not out.exists():
            step(f"satay_mutants_{norm}", [SCRIPTS/"ytab_multicondition_project.py", "--project-config", config(SATAY),
                                           "--normalized", satay_norm/"experiment/methods", "--normalization", norm, "--reference-condition", "WT", "--independent-libraries",
                                           *[x for s, c in SATAY_CONDITIONS.items() for x in ("--sample", f"{SATAY}::{s}={c}")], "--output", out])
        out = BENCH/"pombe_gi"/norm
        if not out.exists():
            fasta, gff = reference(POMBE)
            step(f"pombe_gi_{norm}", [SCRIPTS/"ytab_genetic_interactions.py", "--normalized", pombe_norm/"experiment/methods", "--normalization", norm,
                                      "--gff", gff, "--reference", fasta, "--output", out, "--reference-strain", "WT", "--query-strain", "otr1R",
                                      *[x for s, role in POMBE_GI for x in ("--sample", f"{POMBE}::{s}={role}")]])
    (BENCH/"run.json").write_text(json.dumps(dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), normalizations=NORMS,
                                                 satay_normalization=str(satay_norm), pombe_normalization=str(pombe_norm)), indent=1))
    print("done")


# ---------------------------------------------------------------- evaluation (reads the published positives)

POSITIVES = [  # (test, pair in the output, expected call, genes) — PREREGISTRATION.md table
    ("rapamycin", None, "enriched", "FPR1 RRD1 TIP41 GLN3 SAP190 PSP2 CCS1 ESL2 PPH21 PPM1"),
    ("dpl1del vs WT", "dpl1del vs WT", "depleted", "LCB3"),
    ("dpl1del vs WT", "dpl1del vs WT", "enriched", "HIP1"),
    ("dpl1del_psd2del vs WT", "dpl1del_psd2del vs WT", "depleted", "LCB3 VPS13 PSD1"),
    ("dpl1del_psd2del vs WT", "dpl1del_psd2del vs WT", "enriched", "CDC10"),
    ("mmm1del_VPS13_D716H vs WT", "mmm1del_VPS13_D716H vs WT", "enriched", "MDM10 MDM12 MDM34"),
    ("mmm1del_VPS13_D716H vs WT", "mmm1del_VPS13_D716H vs WT", "depleted", "VPS13 MCP1 VAM6 TOM70 VPS41 YPT7 VMS1 YME1"),
    ("S. pombe GI", None, "Aggravating", "mit1 sir2 msi2 rna14 iss1 ysh1 ssu72"),
]
CANONICAL = "clr4 rik1 raf1 raf2 cul4 chp1 tas3 ago1 dcr1 rdp1 hrr1 cid12 swi6 clr3 mit1 sir2 stc1".split()
MASKED = {"ADE2", "URA3", "CYC1", "GAL1", "GAL10", "HIS3", "LEU2", "MET15", "MET17", "LYS2", "ura4"}
BOOT_SEED, BOOTSTRAPS, FDR = 31415, 2000, .05


def tables(norm):
    import pandas as pd
    rap = pd.read_csv(BENCH/"rapamycin"/norm/"genes.tsv", sep="\t")
    pairs = pd.read_csv(BENCH/"satay_mutants"/norm/"tests/pairs_genes.tsv", sep="\t")
    gi = pd.read_csv(BENCH/"pombe_gi"/norm/"genes.tsv", sep="\t"); gi["ORF"] = gi.ORF.str.removeprefix("gene:")
    return rap, pairs, gi


def recovery(norm):
    import numpy as np
    rap, pairs, gi = tables(norm); rows = []
    for test, pair, expected, genes in POSITIVES:
        if test == "rapamycin":d, stat = rap, "log2FC"
        elif test == "S. pombe GI":d, stat = gi, "delta_log2FC"
        else:d, stat = pairs[pairs.pair == pair], "log2FC"
        d = d.copy(); sign = 1 if expected == "enriched" else -1
        d["rank"] = (-sign*d[stat]).rank(method="min", na_option="bottom")
        index = {str(n).lower():r for n, r in zip(d.Name, d.to_dict("records"))}
        for gene in genes.split():
            r = index.get(gene.lower())
            called = r is not None and r["call"] == expected
            rows.append(dict(normalization=norm, test=test, expected=expected, gene=gene, present=r is not None, called=bool(called),
                             call=None if r is None else r["call"], statistic=None if r is None else r[stat],
                             q=None if r is None else r["q_bh"], rank=None if r is None else int(r["rank"]), tested_genes=int(d[stat].notna().sum())))
    return rows


def sgd_rapamycin():
    import gzip
    import pandas as pd
    import ytab_four_species_benchmark as fs
    cols = list(range(15))
    with gzip.open(ROOT/"resources/species/cerevisiae/essentiality_labels/raw/SGD_phenotype_data.tab.gz", "rt") as h:
        ph = pd.read_csv(h, sep="\t", header=None, names=cols, dtype=str, quoting=3, keep_default_na=False)
    ph = ph[ph[10].str.contains("rapamycin|sirolimus", case=False, na=False) & (ph[6] == "null")
            & ph[5].isin(["classical genetics", "systematic mutation set", "homozygous diploid, systematic mutation set"])]
    sensitive = set(ph[ph[9] == "resistance to chemicals: decreased"][0]); resistant = set(ph[ph[9] == "resistance to chemicals: increased"][0])
    annotated = set(ph[0]); conflict = sensitive & resistant
    gold = fs.labels("cerevisiae"); masked = fs.masked("cerevisiae")
    negatives = {g for g, v in gold.items() if v == 0 and g not in annotated}
    return sensitive-conflict-masked, resistant-conflict-masked, negatives-masked, len(conflict)


def rapamycin_auc(sets):
    import numpy as np
    from sklearn.metrics import roc_auc_score
    sensitive, resistant, negatives, _ = sets; scores = {}
    for norm in NORMS:
        rap = tables(norm)[0].dropna(subset=["log2FC"]); scores[norm] = dict(zip(rap.ORF, rap.log2FC))
    common = set.intersection(*[set(v) for v in scores.values()])
    rows, rng = [], np.random.default_rng(BOOT_SEED)
    for name, positives, sign in (("sensitive (depleted)", sensitive, -1), ("resistant (enriched)", resistant, 1)):
        genes = sorted((positives | negatives) & common); y = np.array([g in positives for g in genes], int)
        s = {n:sign*np.array([scores[n][g] for g in genes]) for n in NORMS}
        aucs = {n:float(roc_auc_score(y, s[n])) for n in NORMS}
        boots = {n:[] for n in NORMS[1:]}
        for _ in range(BOOTSTRAPS):
            i = rng.integers(0, len(y), len(y))
            if 0 < y[i].sum() < len(y):
                base = roc_auc_score(y[i], s["betageom"][i])
                for n in NORMS[1:]:boots[n].append(base-roc_auc_score(y[i], s[n][i]))
        row = dict(phenotype=name, positives=int(y.sum()), negatives=int(len(y)-y.sum()), **{f"auc_{n}":aucs[n] for n in NORMS})
        for n in NORMS[1:]:
            b = np.array(boots[n]); row |= {f"delta_bgc_vs_{n}":aucs["betageom"]-aucs[n], f"ci_low_{n}":float(np.percentile(b, 2.5)), f"ci_high_{n}":float(np.percentile(b, 97.5))}
        rows.append(row)
    return rows


def heterochromatin(norm):
    import numpy as np
    from scipy.stats import fisher_exact
    from sklearn.metrics import roc_auc_score
    gi = tables(norm)[2]; gi = gi[~gi.Name.isin(MASKED)]
    tested = gi[gi.status == "tested"].dropna(subset=["delta_log2FC"])
    canon = tested.Name.str.lower().isin(CANONICAL); agg = tested.call == "Aggravating"
    table = [[int((canon & agg).sum()), int((canon & ~agg).sum())], [int((~canon & agg).sum()), int((~canon & ~agg).sum())]]
    odds, p = fisher_exact(table, alternative="greater")
    named = {g.lower() for _, _, _, genes in POSITIVES[-1:] for g in genes.split()} | set(CANONICAL)
    candidates = tested[agg & ~tested.Name.str.lower().isin(named)].sort_values("q_bh")
    counts = tested.call.value_counts().to_dict()
    # exploratory (not pre-registered): do canonical factors rank toward one end of delta_log2FC without a call?
    auc_down = float(roc_auc_score(canon, -tested.delta_log2FC)) if canon.any() else None
    return dict(normalization=norm, tested=len(tested), calls=counts, fisher_table=table, odds_ratio=float(odds), p=float(p),
                canonical_tested=int(canon.sum()), exploratory_canonical_auc_aggravating_direction=auc_down,
                canonical_directions=tested[canon].interaction.value_counts().to_dict(),
                candidates=candidates[["ORF", "Name", "delta_log2FC", "q_bh"]].to_dict("records"))


def evaluate():
    import pandas as pd
    sys.path.append(str(SCRIPTS))
    prereg = check_preregistration()
    rec = pd.DataFrame([r for n in NORMS for r in recovery(n)]); rec.to_csv(BENCH/"recovery.tsv", sep="\t", index=False)
    summary = (rec.groupby(["normalization", "test", "expected"], sort=False)
                  .agg(positives=("gene", "size"), called=("called", "sum")).reset_index())
    by_test = summary.groupby(["normalization", "test"], sort=False)[["positives", "called"]].sum().reset_index()
    by_test["recovers"] = by_test.called >= by_test.positives/2; by_test.to_csv(BENCH/"recovery_by_test.tsv", sep="\t", index=False)
    sets = sgd_rapamycin(); p2 = rapamycin_auc(sets); pd.DataFrame(p2).to_csv(BENCH/"rapamycin_sgd_auc.tsv", sep="\t", index=False)
    p3 = [heterochromatin(n) for n in NORMS]
    for h in p3:pd.DataFrame(h["candidates"]).to_csv(BENCH/f"pombe_candidates_{h['normalization']}.tsv", sep="\t", index=False)
    counts = {}
    for n in NORMS:
        rap, pairs, gi = tables(n)
        counts[n] = dict(rapamycin=rap.call.value_counts().to_dict(),
                         satay_mutants={p:g.call.value_counts().to_dict() for p, g in pairs.groupby("pair")},
                         pombe_gi=gi.call.value_counts().to_dict())
    result = dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), recovery=by_test.to_dict("records"),
                  sgd_sets=dict(sensitive=len(sets[0]), resistant=len(sets[1]), negatives=len(sets[2]), conflicting_excluded=sets[3]),
                  rapamycin_sgd_auc=p2, heterochromatin=[{k:v for k, v in h.items() if k != "candidates"} | dict(candidates=len(h["candidates"])) for h in p3],
                  call_counts=counts)
    (BENCH/"result.json").write_text(json.dumps(result, indent=1, default=str))
    print(by_test.to_string()); print(json.dumps({k:result[k] for k in ("sgd_sets", "rapamycin_sgd_auc", "heterochromatin")}, indent=1, default=str))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("phase", choices=["run", "evaluate"])
    a = p.parse_args()
    if a.phase == "run":run()
    else:evaluate()


if __name__ == "__main__":main()
