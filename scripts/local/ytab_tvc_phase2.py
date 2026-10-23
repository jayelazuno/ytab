#!/usr/bin/env python3
"""Treated-vs-control benchmark, phase 2 (pre-registration, "Real data, after the simulation"; descriptive).

  ytab tvc-phase2 --output DIR [--workers 6]

Runs every pre-registered arm (ytab.benchmarks.treated_vs_control.ARMS, the code of the amended blind run) on real
libraries, each design with 2 control and 2 treated libraries:
  gale_flu        Gale 2020 Cg1: no-drug controls 1-2 vs fluconazole 1-2
  gale_swap_null  control 1 + FLU 1 vs control 2 + FLU 2 (each side has one control and one FLU culture: no calls expected)
  h2o2_<bg>, zn_<bg>   H2O2 and Zn screens per background (illustration only)
and writes:
  calls.tsv            per design, arm and gene: score, call (+1/-1/0), p
  call_summary.tsv     calls per design and arm
  overlap.tsv          Jaccard overlap of depleted and of enriched calls between arms (gale_flu)
  publication.tsv      gale_flu arms vs Gale's published fluconazole scores
  control_control.tsv  the z-score's control-control comparison (fraction beyond its own call threshold)
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.benchmarks import treated_vs_control as tvc  # noqa: E402
from ytab.essentiality import control_conditioned_gumbel as ccg, resampling  # noqa: E402
from ytab.fitness import zscore  # noqa: E402
from ytab.normalization.ttr import sha256  # noqa: E402

GFF = ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff"
PUBLICATION = ROOT/"output/essentiality/control_conditioned_gumbel/Gale_2020/publication_comparison.tsv"
PROJECTS = ROOT/"output/projects"
# design -> (normalization run, [control 1, control 2, treated 1, treated 2]); samples paired by position
DESIGNS = {
    "gale_flu": (PROJECTS/"Gale_2020/normalization/norm_20261020_231352", ["SRR11561505", "SRR11561504", "SRR11561503", "SRR11561502"]),
    "gale_swap_null": (PROJECTS/"Gale_2020/normalization/norm_20261020_231352", ["SRR11561505", "SRR11561503", "SRR11561504", "SRR11561502"]),
    "h2o2_yH298": (PROJECTS/"H2O2_screen_v1/normalization/norm_20261022_all_methods",
                   ["yH298-parent-pool1", "yH298-parent-pool2", "yH298-H2O2-treated-facs-pool1", "yH298-H2O2-treated-facs-pool2"]),
    "h2o2_yH299": (PROJECTS/"H2O2_screen_v1/normalization/norm_20261022_all_methods",
                   ["yH299-parent-pool3", "yH299-parent-pool4", "yH299-H2O2-treated-facs-pool3", "yH299-H2O2-treated-facs-pool4"]),
    "zn_yH298": (PROJECTS/"Zn_toxicity_screen/normalization/norm_20261022_all_methods",
                 ["yH298-parent-pool1-mock", "yH298-parent-pool2-mock", "yH298-parent-pool1-1_5mM-Zn-treated", "yH298-parent-pool2-1_5mM-Zn-treated"]),
    "zn_yH299": (PROJECTS/"Zn_toxicity_screen/normalization/norm_20261022_all_methods",
                 ["yH299-parent-pool3-mock", "yH299-parent-pool4-mock", "yH299-parent-pool3-1_5mM-Zn-treated", "yH299-parent-pool4-1_5mM-Zn-treated"]),
}
def _rows_tolerant(rows, signed):
    """tvc._rows (frozen with the benchmark) for real data, where a method can leave every gene untested: its table then
    lacks p_value/log2FC/call columns, which are taken as missing (no call). Identical output whenever they exist."""
    d = pd.DataFrame(rows)
    for column, default in (("p_value", np.nan), ("q_bh", np.nan), ("log2FC", np.nan), ("call", "not tested"), ("in_family", False)):
        if column not in d:d[column] = default
    return _frozen_rows(d.to_dict("records"), signed)


_frozen_rows = tvc._rows
tvc._rows = _rows_tolerant

PRIMARY = ("zscore_cpm", "zscore_ttr", "zscore_bgc", "resampling_totreads", "resampling_ttr", "resampling_bgc", "gumbel_totreads", "gumbel_ttr")


def load(run: Path, samples: list[str]):
    matrix = run/"experiment/matrix"
    obs = [r["sample"] for r in csv.DictReader(open(matrix/"observations.tsv"), delimiter="\t")]
    with gzip.open(matrix/"sites.tsv.gz", "rt") as fh:
        sites = [(r["chromosome"], int(r["position_1based"])) for r in csv.DictReader(fh, delimiter="\t")]
    raw = sparse.load_npz(matrix/"raw_counts.npz").tocsc()[:, [obs.index(s) for s in samples]].toarray()
    keep = raw.sum(1) > 0   # sites observed in any of the four libraries
    return raw[keep], [s for s, k in zip(sites, keep) if k]


def run_design(name: str, out: Path | None = None) -> pd.DataFrame:
    saved = out/f"calls_{name}.tsv" if out else None
    if saved and saved.is_file():return pd.read_csv(saved, sep="\t")
    run, samples = DESIGNS[name]
    raw, sites = load(run, samples)
    genes = ccg.annotations(GFF)
    gene_r = ccg.gene_sites(genes, sites, resampling.Settings()); gene_g = ccg.gene_sites(genes, sites, ccg.Settings())
    positions = np.array([p for _, p in sites], dtype=np.int64)
    results = tvc.run_arms(raw, None, genes, gene_r, gene_g, positions)
    rows = [r.assign(design=name, arm=arm).rename_axis("gene").reset_index() for arm, r in results.items()]
    calls = pd.concat(rows, ignore_index=True)
    if saved:calls.to_csv(saved, sep="\t", index=False)
    print(f"done {name}", file=sys.stderr, flush=True)
    return calls


def control_control(name: str) -> list[dict]:
    """The z-score's own control-control comparison: share of genes beyond its call threshold (1% by construction)."""
    run, samples = DESIGNS[name]
    raw, sites = load(run, samples)
    genes = ccg.annotations(GFF); idx = ccg.gene_sites(genes, sites, resampling.Settings()); names = [g["id"] for g in genes]
    out = []
    for norm in ("cpm", "ttr", "betageom"):
        if norm == "cpm":
            totals = tvc._gene_totals(raw, idx, names); values = totals/totals.sum()*1e6
        else:
            values = tvc._gene_totals(tvc.normalize(sparse.csc_matrix(raw), norm, ["sim"]*4)[0].toarray(), idx, names)
        with tempfile.TemporaryDirectory() as tmp:
            design = Path(tmp)/"design.csv"; design.write_text(tvc.DESIGN)
            zscore.run(design, Path(tmp)/"out", normalized=dict(genes=values.rename_axis("feature_id"), label="phase2"))
            cc = pd.read_csv(Path(tmp)/"out/tables/control_control_z_scores.csv")
            thr = pd.read_csv(Path(tmp)/"out/tables/background_specific_z_thresholds.csv").z_threshold.iloc[0]
        z = cc.z.replace([np.inf, -np.inf], np.nan).dropna()
        out.append(dict(design=name, normalization=norm, genes=len(z), threshold=thr, beyond=int((z.abs() > thr).sum()),
                        fraction_beyond=float((z.abs() > thr).mean())))
    return out


def jaccard(a: set, b: set) -> float:
    return len(a & b)/len(a | b) if a | b else float("nan")


def publication(calls: pd.DataFrame) -> pd.DataFrame:
    pub = pd.read_csv(PUBLICATION, sep="\t").set_index("ORF")
    pub = pub[pub.published_z.notna()]
    strict = set(pub.index[pub.publication_class.astype(str).str.startswith("A")])   # z <= -3 and log2FC <= -1
    z3 = set(pub.index[pub.published_z <= -3]); z3e = set(pub.index[pub.published_z >= 3])
    rows = []
    for arm, d in calls[calls.design == "gale_flu"].groupby("arm"):
        d = d.set_index("gene").reindex(pub.index)
        score = d.score.fillna(0.)
        ranked = score.sort_values(kind="stable").index   # most depleted first (all arms: lower = more depleted)
        dep = set(d.index[d.call == -1]); enr = set(d.index[d.call == 1])
        rows.append(dict(arm=arm, genes=len(pub), spearman_vs_published_z=float(score.corr(pub.published_z, method="spearman")),
                         published_strict=len(strict), strict_in_top_k=len(strict & set(ranked[:len(strict)])),
                         strict_called_depleted=len(strict & dep), published_z_le_3=len(z3), z_le_3_in_top_k=len(z3 & set(ranked[:len(z3)])),
                         z_le_3_called_depleted=len(z3 & dep), depleted_calls=len(dep),
                         published_z_ge_3=len(z3e), z_ge_3_called_enriched=len(z3e & enr), enriched_calls=len(enr),
                         upc2a_rank=int(list(ranked).index("GWK60_C00759"))+1 if "GWK60_C00759" in ranked else None))
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", required=True, type=Path); p.add_argument("--workers", type=int, default=6)
    a = p.parse_args(); out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    sources = {str(f.relative_to(ROOT)):sha256(f) for f in (ROOT/"src/ytab/benchmarks/treated_vs_control.py", ROOT/"src/ytab/fitness/zscore.py",
                                                            ROOT/"src/ytab/essentiality/resampling.py", ROOT/"src/ytab/essentiality/control_conditioned_gumbel.py",
                                                            Path(__file__), PUBLICATION)}
    (out/"sources.json").write_text(json.dumps(dict(sources=sources, designs={k:dict(run=str(v[0]), samples=v[1]) for k, v in DESIGNS.items()}), indent=1))
    with ProcessPoolExecutor(a.workers) as pool:
        calls = pd.concat(pool.map(run_design, DESIGNS, [out]*len(DESIGNS)), ignore_index=True)
        cc = [r for rows in pool.map(control_control, ["gale_flu"]) for r in rows]
    calls.to_csv(out/"calls.tsv", sep="\t", index=False)
    summary = calls.groupby(["design", "arm"]).agg(depleted=("call", lambda c:int((c == -1).sum())), enriched=("call", lambda c:int((c == 1).sum())),
                                                   tested=("p", lambda x:int(x.notna().sum()))).reset_index()
    summary.to_csv(out/"call_summary.tsv", sep="\t", index=False)
    flu = calls[calls.design == "gale_flu"]
    sets = {arm:(set(d.gene[d.call == -1]), set(d.gene[d.call == 1])) for arm, d in flu.groupby("arm") if arm in PRIMARY}
    overlap = pd.DataFrame([dict(arm_a=x, arm_b=y, depleted_jaccard=jaccard(sets[x][0], sets[y][0]), enriched_jaccard=jaccard(sets[x][1], sets[y][1]),
                                 depleted_shared=len(sets[x][0] & sets[y][0]), enriched_shared=len(sets[x][1] & sets[y][1]))
                            for i, x in enumerate(PRIMARY) for y in PRIMARY[i+1:] if x in sets and y in sets])
    overlap.to_csv(out/"overlap.tsv", sep="\t", index=False)
    publication(calls).to_csv(out/"publication.tsv", sep="\t", index=False)
    pd.DataFrame(cc).to_csv(out/"control_control.tsv", sep="\t", index=False)
    print(summary.pivot(index="arm", columns="design", values="depleted").to_string())


if __name__ == "__main__":
    main()
