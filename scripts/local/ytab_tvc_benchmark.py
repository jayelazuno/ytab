#!/usr/bin/env python3
"""Treated-vs-control benchmark (docs/manuscript/preregistration_treated_vs_control.md).

  ytab tvc-benchmark calibrate --output DIR                calibration from the base project (default: Gale 2020 FLU)
  ytab tvc-benchmark dev --output DIR [--experiments 2]     development run on seeds 7001-7100 (never the benchmark seed)
  ytab tvc-benchmark blind --output DIR                     the pre-registered run; refuses unless the pre-registration is frozen
  ytab tvc-benchmark amend --output DIR                     amendment 1: rerun the amended z-score arms on the same experiments
  ytab tvc-benchmark extend --output DIR                    R27: the arms R13 never measured, on the same experiments
  ytab tvc-benchmark extension-report --output DIR          R27's non-inferiority tests and gates -> extension_report.md
  ytab tvc-benchmark report --output DIR                    metrics, paired tests, calibration gates -> report.md

The base project needs one Site Normalization run with two control and two treated libraries in its experiment matrix.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import json
from pathlib import Path
import pickle
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.benchmarks import treated_vs_control as tvc  # noqa: E402
from ytab.essentiality import control_conditioned_gumbel as ccg, resampling  # noqa: E402
from ytab.normalization.ttr import sha256  # noqa: E402

PREREG = ROOT/"docs/manuscript/preregistration_treated_vs_control.md"
FROZEN = ROOT/"docs/manuscript/preregistration_treated_vs_control.sha256"
AMENDMENT = ROOT/"docs/manuscript/preregistration_treated_vs_control.amendment1.md"
AMENDMENT_FROZEN = ROOT/"docs/manuscript/preregistration_treated_vs_control.amendment1.sha256"
# Amendment 1: z-score arms whose input is not already per million (all other arms are unchanged by it).
AMENDED_ARMS = ("zscore_ttr", "zscore_bgc", "zscore_ttr_pairs", "zscore_bgc_pairs")
BASE = dict(run=ROOT/"output/projects/Gale_2020/normalization/norm_20261020_231352",
            controls=["SRR11561505", "SRR11561504"], treated=["SRR11561503", "SRR11561502"],
            gff=ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff")
BENCHMARK_SEED, DEV_SEEDS = 27182, range(7001, 7101)
# R27 (preregistration_normalization_arms.md): the arms each statistic accepts that R13 never measured, each tested
# for non-inferior ranking against that statistic's recommended arm, and against the same three calibration gates.
EXTENSION = ROOT/"docs/manuscript/preregistration_normalization_arms.md"
EXTENSION_FROZEN = ROOT/"docs/manuscript/preregistration_normalization_arms.sha256"
EXTENSION_MARGIN = -0.02
EXTENSION_REFERENCE = {"zscore":"zscore_bgc", "resampling":"resampling_bgc", "gumbel":"gumbel_ttr"}
EXTENSION_FAMILIES = {   # family -> {hypothesis: (arm, reference arm, metric)}
    family:{f"H1 {arm} vs {EXTENSION_REFERENCE[family]} (AUC)":(arm, EXTENSION_REFERENCE[family], "auc")
            for arm, (method, _, _) in tvc.EXTENSION_ARMS.items() if method == family}
    for family in ("zscore", "resampling", "gumbel")
}
FAMILIES = {
    "normalization": {"H1a z-score: BGC vs CPM": ("zscore_bgc", "zscore_cpm", "auc"),
                      "H1b resampling: BGC vs total reads": ("resampling_bgc", "resampling_totreads", "auc"),
                      "H1c Gumbel: TTR vs total reads": ("gumbel_ttr", "gumbel_totreads", "auc"),
                      "H1d z-score: TTR vs CPM": ("zscore_ttr", "zscore_cpm", "auc"),
                      "H1e resampling: TTR vs total reads": ("resampling_ttr", "resampling_totreads", "auc"),
                      "H1f z-score: BGC vs TTR": ("zscore_bgc", "zscore_ttr", "auc"),
                      "H1g resampling: BGC vs TTR": ("resampling_bgc", "resampling_ttr", "auc")},
    "method": {"H2a resampling BGC vs z-score CPM": ("resampling_bgc", "zscore_cpm", "precision_at_k"),
               "H2b resampling BGC vs z-score BGC": ("resampling_bgc", "zscore_bgc", "precision_at_k"),
               "H2c Gumbel TTR vs z-score BGC (depletion)": ("gumbel_ttr", "zscore_bgc", "precision_at_k_depletion")},
    "calls": {"H3a resampling BGC vs z-score CPM": ("resampling_bgc", "zscore_cpm", "fdp"),
              "H3b resampling BGC vs z-score BGC": ("resampling_bgc", "zscore_bgc", "fdp"),
              "H3c Gumbel TTR vs z-score BGC (depletion)": ("gumbel_ttr", "zscore_bgc", "fdp_depletion")},
}
GATED = ("resampling_bgc", "gumbel_ttr")


def calibrate(out: Path):
    matrix = BASE["run"]/"experiment/matrix"
    obs = [r["sample"] for r in csv.DictReader(open(matrix/"observations.tsv"), delimiter="\t")]
    with gzip.open(matrix/"sites.tsv.gz", "rt") as fh:
        sites = [(r["chromosome"], int(r["position_1based"])) for r in csv.DictReader(fh, delimiter="\t")]
    raw = sparse.load_npz(matrix/"raw_counts.npz")
    cal = tvc.calibrate(raw, sites, [obs.index(s) for s in BASE["controls"]], [obs.index(s) for s in BASE["treated"]])
    genes = ccg.annotations(BASE["gff"])
    gene_r = resampling_sites(genes, cal.sites); gene_g = ccg.gene_sites(genes, cal.sites, ccg.Settings())
    out.mkdir(parents=True, exist_ok=True)
    pickle.dump(dict(cal=cal, genes=genes, gene_r=gene_r, gene_g=gene_g), open(out/"calibration.pkl", "wb"))
    info = dict(base={k:str(v) for k, v in BASE.items()}, sites=len(cal.sites), exposures=cal.exposures.round(4).tolist(),
                control_phi=dict(log_mu=cal.control_phi[0].round(3).tolist(), log_phi=cal.control_phi[1].round(3).tolist()),
                treated_phi=dict(log_mu=cal.treated_phi[0].round(3).tolist(), log_phi=cal.treated_phi[1].round(3).tolist()),
                sources={p.name:sha256(p) for p in (matrix/"raw_counts.npz", matrix/"sites.tsv.gz", BASE["gff"], Path(tvc.__file__))})
    (out/"calibration.json").write_text(json.dumps(info, indent=1)+"\n")
    return info


def resampling_sites(genes, sites):
    return ccg.gene_sites(genes, sites, resampling.Settings())


def task(args):
    out, scenario, experiment, seed, *arms = args
    data = pickle.load(open(Path(out)/"calibration.pkl", "rb"))
    rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(list(tvc.SCENARIOS).index(scenario), experiment)))
    raw, truth = tvc.simulate(data["cal"], data["gene_r"], scenario, rng)
    names = [g["id"] for g in data["genes"]]; truth = pd.Series(truth, index=names)
    positions = np.array([p for _, p in data["cal"].sites], dtype=np.int64)
    results = tvc.run_arms(raw, data["cal"], data["genes"], data["gene_r"], data["gene_g"], positions,
                           {k:tvc.ALL_ARMS[k] for k in arms[0]} if arms else tvc.ARMS)
    rows = []
    for arm, r in results.items():
        row = dict(scenario=scenario, experiment=experiment, arm=arm) | tvc.score(r, truth, depletion_only=arm.startswith("gumbel"))
        if not arm.startswith("gumbel"):
            dep = tvc.score(r.assign(call=r.call.where(r.call < 0, 0)), truth.where(truth <= 0, 0), depletion_only=True)
            row["precision_at_k_depletion"], row["fdp_depletion"] = dep.get("precision_at_k"), dep["fdp"]
        else:
            row["precision_at_k_depletion"], row["fdp_depletion"] = row.get("precision_at_k"), row["fdp"]
        rows.append(row)
    return rows


def run(out: Path, seed: int, experiments: int, null_experiments: int, workers: int):
    if not (out/"calibration.pkl").exists():calibrate(out)
    tasks = [(str(out), s, e, seed) for s in tvc.SCENARIOS for e in range(null_experiments if s == "null" else experiments)]
    rows = []
    with ProcessPoolExecutor(workers) as pool:
        for i, result in enumerate(pool.map(task, tasks), 1):
            rows += result; print(f"{i}/{len(tasks)} {result[0]['scenario']} #{result[0]['experiment']}", file=sys.stderr, flush=True)
            pd.DataFrame(rows).to_csv(out/"metrics.tsv", sep="\t", index=False)
    return pd.DataFrame(rows)


def amend(out: Path, seed: int, experiments: int, null_experiments: int, workers: int):
    """Rerun AMENDED_ARMS on the same experiments and replace their rows; the original metrics are kept."""
    original = out/"metrics.original.tsv"
    if not original.exists():(out/"metrics.tsv").rename(original)
    m = pd.read_csv(original, sep="\t", keep_default_na=False, na_values=[""])
    tasks = [(str(out), s, e, seed, AMENDED_ARMS) for s in tvc.SCENARIOS for e in range(null_experiments if s == "null" else experiments)]
    rows = []
    with ProcessPoolExecutor(workers) as pool:
        for i, result in enumerate(pool.map(task, tasks), 1):
            rows += result; print(f"amend {i}/{len(tasks)} {result[0]['scenario']} #{result[0]['experiment']}", file=sys.stderr, flush=True)
    new = pd.DataFrame(rows)
    pd.concat([m[~m.arm.isin(AMENDED_ARMS)], new], ignore_index=True).to_csv(out/"metrics.tsv", sep="\t", index=False)


def extend(out: Path, seed: int, experiments: int, null_experiments: int, workers: int):
    """R27: run the extension arms on the same experiments as the blind run and add their rows to metrics.tsv."""
    metrics = out/"metrics.tsv"
    m = pd.read_csv(metrics, sep="\t", keep_default_na=False, na_values=[""])
    kept = m[~m.arm.isin(tvc.EXTENSION_ARMS)]
    names = tuple(tvc.EXTENSION_ARMS)
    tasks = [(str(out), s, e, seed, names) for s in tvc.SCENARIOS for e in range(null_experiments if s == "null" else experiments)]
    rows = []
    with ProcessPoolExecutor(workers) as pool:
        for i, result in enumerate(pool.map(task, tasks), 1):
            rows += result; print(f"extend {i}/{len(tasks)} {result[0]['scenario']} #{result[0]['experiment']}", file=sys.stderr, flush=True)
            pd.concat([kept, pd.DataFrame(rows)], ignore_index=True).to_csv(metrics, sep="\t", index=False)


def extension_report(out: Path):
    """R27's own report: non-inferiority of each extension arm's ranking, and the three gates per arm."""
    m = pd.read_csv(out/"metrics.tsv", sep="\t", keep_default_na=False, na_values=[""])
    signal = m[m.scenario != "null"]; null = m[m.scenario == "null"]
    missing = [a for a in tvc.EXTENSION_ARMS if a not in set(m.arm)]
    if missing:raise SystemExit("extension arms not in metrics.tsv: " + ", ".join(missing))
    lines = ["# R27 — every normalization each statistic accepts", "",
             f"Margin for non-inferior ranking: {EXTENSION_MARGIN} AUC against the statistic's recommended arm.", ""]
    tests = []
    for family, hyps in EXTENSION_FAMILIES.items():
        ps = {}
        for name, (a, b, metric) in hyps.items():
            wide = signal.pivot_table(index=["scenario", "experiment"], columns="arm", values=metric)
            diff = (wide[a]-wide[b]).to_numpy(); ps[name] = tvc.sign_flip(diff)
            lo, hi = tvc.bootstrap_ci(diff)
            tests.append(dict(family=family, arm=a, reference=b, metric=metric, mean_difference=np.nanmean(diff),
                              ci_low=lo, ci_high=hi, p=ps[name], experiments=int(np.isfinite(diff).sum())))
        adjusted = tvc.holm(ps)
        for t, name in zip([t for t in tests if t["family"] == family], hyps):t["p_holm"] = adjusted[name]
    tests = pd.DataFrame(tests)
    tests["ranking"] = np.where(tests.ci_low > EXTENSION_MARGIN, np.where((tests.p_holm < .05) & (tests.ci_low > 0), "better", "not worse"),
                                np.where(tests.ci_high < EXTENSION_MARGIN, "worse", "no claim"))
    tests.round(5).to_csv(out/"extension_tests.tsv", sep="\t", index=False)
    lines += ["## H1 ranking, non-inferiority within family (Holm)", "", markdown(tests.round(4)), ""]
    gates = []
    for arm in tvc.EXTENSION_ARMS:
        n = null[null.arm == arm]; tested, p05 = n.unaffected_tested.sum(), n.unaffected_p05.sum()
        zscore_arm = tvc.EXTENSION_ARMS[arm][0] == "zscore"     # no p-value: gate 1 does not apply, as in R13
        lower = beta.ppf(.025, p05, tested-p05+1) if p05 else 0.
        any_call = (n.calls > 0).mean(); s = signal[signal.arm == arm]
        fdr = s.false_calls.sum()/max(s.calls.sum(), 1)
        gates.append(dict(arm=arm, null_p05_rate="" if zscore_arm else p05/max(tested, 1),
                          null_p05_lower95="" if zscore_arm else lower, gate1="n/a" if zscore_arm else bool(lower <= .05),
                          null_experiments_with_calls=any_call, gate2=bool(any_call <= .10), pooled_fdr=fdr, gate3=bool(fdr <= .10)))
    gates = pd.DataFrame(gates)
    gates["pass"] = gates.gate2 & gates.gate3 & gates.gate1.isin([True, "n/a"])
    gates.round(4).to_csv(out/"extension_gates.tsv", sep="\t", index=False)
    lines += ["## H2 calibration gates per arm", "", markdown(gates.round(4)), "",
              "## Cleared for", ""]
    verdict = tests.set_index("arm").ranking.to_dict(); passed = gates.set_index("arm")["pass"].to_dict()
    rows = [dict(arm=a, ranking=verdict[a], gates="pass" if passed[a] else "fail",
                 cleared_for="ranking and calls" if verdict[a] in ("not worse", "better") and passed[a]
                 else "ranking only" if verdict[a] in ("not worse", "better") else "not cleared")
            for a in tvc.EXTENSION_ARMS]
    lines += [markdown(pd.DataFrame(rows)), ""]
    (out/"extension_report.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines))


def frozen_digest(document: Path, record: Path, p, verify_sources=False) -> tuple[str, str]:
    """The first line of the record is the document's own digest and the UTC freeze time; every later line is
    `digest path`, a source file the design depends on. R27 checks those too, as its text promises."""
    if not record.is_file():p.error(f"{document.name} is not frozen ({record.name} missing)")
    lines = [l.split() for l in record.read_text().splitlines() if l.strip()]
    digest, when = lines[0][:2]
    if sha256(document) != digest:p.error(f"{document.name} changed after it was frozen; refusing to run")
    if verify_sources:
        for entry in lines[1:]:
            if len(entry) < 2:continue
            recorded, path = entry[0], ROOT/entry[1]
            if not path.is_file():p.error(f"{entry[1]} is recorded in {record.name} but missing; refusing to run")
            if sha256(path) != recorded:p.error(f"{entry[1]} changed after {document.name} was frozen; refusing to run")
    return digest, when


def markdown(df: pd.DataFrame) -> str:
    fmt = lambda v: f"{v:.4g}" if isinstance(v, float) else str(v)
    rows = [list(map(str, df.columns))] + [[fmt(v) for v in r] for r in df.itertuples(index=False)]
    return "\n".join(["| " + " | ".join(rows[0]) + " |", "|" + "---|"*len(rows[0])] + ["| " + " | ".join(r) + " |" for r in rows[1:]])


def report(out: Path):
    m = pd.read_csv(out/"metrics.tsv", sep="\t", keep_default_na=False, na_values=[""]); signal = m[m.scenario != "null"]  # "null" is a scenario name, not NA
    lines = ["# Treated-vs-control benchmark report", "", f"Experiments: {m.groupby('scenario').experiment.nunique().to_dict()}", ""]
    summary = m.groupby(["scenario", "arm"]).agg(auc=("auc", "mean"), precision_at_k=("precision_at_k", "mean"), power=("power", "mean"),
                                                 calls=("calls", "sum"), false_calls=("false_calls", "sum")).reset_index()
    summary["empirical_fdr"] = summary.false_calls/summary.calls.where(summary.calls > 0)
    summary.round(4).to_csv(out/"summary.tsv", sep="\t", index=False)
    lines += ["## Mean metrics by scenario and arm", "", markdown(summary.round(3)), ""]
    tests = []
    for family, hyps in FAMILIES.items():
        ps = {}
        for name, (a, b, metric) in hyps.items():
            wide = signal.pivot_table(index=["scenario", "experiment"], columns="arm", values=metric)
            diff = (wide[a]-wide[b]).to_numpy(); ps[name] = tvc.sign_flip(diff)
            lo, hi = tvc.bootstrap_ci(diff)
            tests.append(dict(family=family, hypothesis=name, metric=metric, mean_difference=np.nanmean(diff), ci_low=lo, ci_high=hi, p=ps[name], experiments=int(np.isfinite(diff).sum())))
        adjusted = tvc.holm(ps)
        for t in tests:
            if t["family"] == family:t["p_holm"] = adjusted[t["hypothesis"]]
    tests = pd.DataFrame(tests); tests["significant"] = (tests.p_holm < .05) & ((tests.ci_low > 0) | (tests.ci_high < 0))
    tests.round(5).to_csv(out/"tests.tsv", sep="\t", index=False)
    lines += ["## Pre-registered tests (Holm within family)", "", markdown(tests.round(4)), ""]
    gates = []
    null = m[m.scenario == "null"]
    for arm in GATED:
        n = null[null.arm == arm]; tested, p05 = n.unaffected_tested.sum(), n.unaffected_p05.sum()
        lower = beta.ppf(.025, p05, tested-p05+1) if p05 else 0.
        any_call = (n.calls > 0).mean()
        s = signal[signal.arm == arm]; fdr = s.false_calls.sum()/max(s.calls.sum(), 1)
        gates.append(dict(arm=arm, null_p05_rate=p05/max(tested, 1), null_p05_lower95=lower, gate1=lower <= .05,
                          null_experiments_with_calls=any_call, gate2=any_call <= .10, pooled_fdr=fdr, gate3=fdr <= .10))
    gates = pd.DataFrame(gates); gates["pass"] = gates.gate1 & gates.gate2 & gates.gate3
    gates.round(4).to_csv(out/"gates.tsv", sep="\t", index=False)
    lines += ["## Calibration gates", "", markdown(gates.round(4)), ""]
    (out/"report.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("phase", choices=["calibrate", "dev", "blind", "amend", "extend", "extension-report", "report"])
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--experiments", type=int, default=20, help="per signal scenario")
    p.add_argument("--null-experiments", type=int, default=40); p.add_argument("--workers", type=int, default=6)
    p.add_argument("--dev-seed", type=int, default=DEV_SEEDS[0])
    a = p.parse_args(); out = a.output.resolve()
    if a.phase == "calibrate":print(json.dumps(calibrate(out), indent=1)[:2000]); return
    if a.phase == "report":report(out); return
    if a.phase == "extension-report":extension_report(out); return
    if a.phase == "dev":
        if a.dev_seed not in DEV_SEEDS:p.error("development seeds are 7001-7100")
        run(out, a.dev_seed, a.experiments, a.null_experiments, a.workers); report(out); return
    if a.phase == "extend":
        digest, when = frozen_digest(EXTENSION, EXTENSION_FROZEN, p, verify_sources=True)
        if not (out/"frozen.json").is_file():p.error("extend applies to a completed blind run")
        info = json.loads((out/"frozen.json").read_text()) | dict(r27_sha256=digest, r27_frozen_utc=when)
        (out/"frozen.json").write_text(json.dumps(info, indent=1))
        extend(out, BENCHMARK_SEED, 20, 40, a.workers); extension_report(out); return
    if a.phase == "amend":
        digest, when = frozen_digest(AMENDMENT, AMENDMENT_FROZEN, p)
        if not (out/"frozen.json").is_file():p.error("amend applies to a completed blind run")
        info = json.loads((out/"frozen.json").read_text()) | dict(amendment1_sha256=digest, amendment1_frozen_utc=when)
        (out/"frozen.json").write_text(json.dumps(info, indent=1))
        amend(out, BENCHMARK_SEED, 20, 40, a.workers); report(out); return
    digest, when = frozen_digest(PREREG, FROZEN, p)
    (out/"frozen.json").parent.mkdir(parents=True, exist_ok=True)
    (out/"frozen.json").write_text(json.dumps(dict(preregistration_sha256=digest, frozen_utc=when, seed=BENCHMARK_SEED), indent=1))
    if (a.experiments, a.null_experiments) != (20, 40):p.error("The pre-registered run uses 20 experiments per signal scenario and 40 null experiments")
    run(out, BENCHMARK_SEED, a.experiments, a.null_experiments, a.workers); report(out)


if __name__ == "__main__":
    main()
