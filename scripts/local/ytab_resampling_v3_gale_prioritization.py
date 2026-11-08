#!/usr/bin/env python3
"""POST-HOC exploratory comparison (not acceptance): resampling v3 Gale calls (BGC + trimmed-sum recentering,
output/essentiality/resampling/Gale_2020_bgc_v3/genes.tsv) vs Gale's published S2 scores and the control-conditioned
Gumbel, and how a candidate ranking behaves against the published scores.

The prioritization rule is fixed here before reading the published scores and is not tuned to them:
  depleted   tier 1 = resampling depleted AND Gumbel depleted run (two independent tests agree)
             tier 2 = resampling depleted only
             tier 3 = Gumbel depleted run only (resampling not significant)
  enriched   tier 1 = resampling enriched with >= 5 sites; tier 2 = resampling enriched with < 5 sites
  within a tier: q_bh ascending, then |log2FC| descending, then n_sites descending.
Published scores are used only to evaluate. Writes output/essentiality/resampling/v3_gale_prioritization/.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"output/essentiality/resampling/v3_gale_prioritization"
V3 = ROOT/"output/essentiality/resampling/Gale_2020_bgc_v3/genes.tsv"
PUB = ROOT/"output/essentiality/control_conditioned_gumbel/Gale_2020/publication_comparison.tsv"
MIN_SITES_ENRICHED = 5


def ranked(d, direction):
    if direction == "depleted":
        res, gum = d.call == "depleted", d.gumbel_call == "depleted run"
        tier = np.select([res & gum, res, gum & ~res], [1, 2, 3], 0)
    else:
        res = d.call == "enriched"; tier = np.select([res & (d.n_sites >= MIN_SITES_ENRICHED), res], [1, 2], 0)
    c = d.assign(tier=tier, abs_lfc=d.log2FC.abs())[tier > 0]
    c = c.sort_values(["tier", "q_bh", "abs_lfc", "n_sites"], ascending=[True, True, False, False], na_position="last")
    return c.assign(rank=np.arange(1, len(c)+1)).drop(columns="abs_lfc")


def hit(c, direction):
    return c.published_z <= -2 if direction == "depleted" else c.published_z >= 2


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = pd.read_csv(V3, sep="\t")
    pub = pd.read_csv(PUB, sep="\t")[["ORF", "publication_class", "published_z", "published_log2fc", "call", "q_bh"]] \
        .rename(columns={"call":"gumbel_call", "q_bh":"gumbel_q"})
    d = res.merge(pub, on="ORF", how="left")
    scored = d[d.published_z.notna()]
    summary = dict(genes=len(d), tested=int((d.status == "tested").sum()), with_published_score=len(scored),
                   base_rate_published_z_le_minus2=float((scored.published_z <= -2).mean()),
                   base_rate_published_z_ge_2=float((scored.published_z >= 2).mean()))
    t = scored[scored.status == "tested"]
    rho = spearmanr(t.log2FC, t.published_log2fc); summary["spearman_log2fc_vs_published_tested"] = dict(rho=float(rho[0]), n=len(t))
    tiers, tables = [], []
    for direction in ("depleted", "enriched"):
        c = ranked(d, direction); h = hit(c, direction); opposite = c.published_z >= 0 if direction == "depleted" else c.published_z <= 0
        c = c.assign(published_agrees=h, published_opposite=opposite); tables.append(c.assign(direction=direction))
        for tier, g in c.groupby("tier"):
            s = g.published_z.notna()
            tiers.append(dict(direction=direction, tier=int(tier), genes=len(g), with_published_score=int(s.sum()),
                              published_agrees=int(g.published_agrees.sum()), published_opposite=int(g.published_opposite.sum()),
                              agree_fraction=round(g.published_agrees.sum()/max(s.sum(), 1), 3),
                              median_published_z=None if not s.any() else round(float(g.published_z.median()), 2),
                              median_log2FC=round(float(g.log2FC.median()), 2), median_sites=float(g.n_sites.median())))
        # precision at k along the ranked list (genes without a published score count as not agreeing)
        summary[f"{direction}_precision_at_k"] = {k:round(float(c.published_agrees.head(k).mean()), 3) for k in (10, 25, 50, 100, 200) if k <= len(c)}
        # how many strong published genes appear in the list and where
        strong = scored[scored.published_z <= -3] if direction == "depleted" else scored[scored.published_z >= 3]
        found = c[c.ORF.isin(strong.ORF)]
        summary[f"{direction}_published_strong"] = dict(threshold="z<=-3" if direction == "depleted" else "z>=3", genes=len(strong),
            in_ranked_list=len(found), median_rank=None if found.empty else float(found["rank"].median()),
            untested=int(d[d.ORF.isin(strong.ORF)].status.ne("tested").sum()))
    strict = d[d.publication_class.fillna("").str.startswith("A")]
    summary["strict_published_19"] = dict(tested=int((strict.status == "tested").sum()), resampling_depleted=int((strict.call == "depleted").sum()),
                                          gumbel_depleted=int((strict.gumbel_call == "depleted run").sum()))
    allc = pd.concat(tables)
    cols = ["direction", "rank", "tier", "ORF", "Name", "n_sites", "log2FC", "p_value", "q_bh", "call", "gumbel_call", "gumbel_q",
            "published_z", "published_log2fc", "publication_class", "published_agrees", "published_opposite", "Description"]
    allc[cols].to_csv(OUT/"ranked_candidates.tsv", sep="\t", index=False)
    pd.DataFrame(tiers).to_csv(OUT/"tiers_vs_published.tsv", sep="\t", index=False)
    strict[["Name", "ORF", "n_sites", "status", "log2FC", "q_bh", "call", "gumbel_call", "published_z", "published_log2fc"]] \
        .to_csv(OUT/"strict_published_genes.tsv", sep="\t", index=False)
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=1)); print(pd.DataFrame(tiers).to_string(index=False))
    print(strict[["Name", "n_sites", "status", "log2FC", "q_bh", "call", "gumbel_call", "published_z"]].to_string(index=False))
    for direction in ("depleted", "enriched"):
        print(f"\nTop 15 {direction}"); print(allc[allc.direction == direction][["rank", "tier", "Name", "ORF", "n_sites", "log2FC", "q_bh", "gumbel_call", "published_z"]].head(15).to_string(index=False))


if __name__ == "__main__":main()
