"""Agreement between treated-vs-control methods on the same comparisons: the local-abundance z-score, resampling and
the control-conditioned Gumbel (and, optionally, published calls). Writes every table the app plots; the app reads
them and derives nothing.

Per background x treated condition ("group"), genes tested by both methods of a pair are compared on
- effect: Spearman of the effects (z-score and resampling log2FC; Gumbel log2 scaled read ratio). Rank correlation,
  because the scales and tails differ between methods;
- score: Spearman of signed significance (z for the z-score; sign(effect) x -log10 p for resampling; -(-log10 p) for
  the depletion-only Gumbel, so depletion is negative for every method);
- calls, per direction: how many each calls, how many both call, and recovery both ways, |A n B| / |A| (of A's calls,
  the share B also calls) and |A n B| / |B|. With published calls as A this is the recall of the published hits.
  No Jaccard: it is capped by unequal list sizes (37 published vs 250 called caps it at 0.15 with every published gene
  recovered), so it reads as disagreement where there is none. The Gumbel calls depletion only, so it is left out of
  enriched-call comparisons.

Outputs (output directory): genes.tsv (one row per group x gene: each method's effect, score, call), pairs.tsv (one row
per group x method pair), sets.tsv (exclusive call intersections per group and direction, for Venn/UpSet), manifest.json.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ytab.fitness.contract import VERSION

LABELS = {"zscore": "Local-abundance z-score", "resampling": "Resampling", "gumbel": "Control-conditioned Gumbel", "published": "Published"}
DEPLETION_ONLY = {"gumbel"}


def _zscore(path: Path) -> pd.DataFrame:
    table = path/"tables/treated_vs_parent.combined.log2fc_z.csv"
    x = pd.read_csv(table if table.exists() or not table.with_suffix(".csv.gz").exists() else table.with_suffix(".csv.gz"))   # archived runs keep it gzipped
    return pd.DataFrame(dict(background=x.background.astype(str), condition=x.treated_condition.astype(str), ORF=x.feature_id.astype(str),
                             Name=x.get("common_name", x.feature_id).astype(str), effect=x.log2FC, score=x.z,
                             call=np.where(x.call.isin(["depleted", "enriched"]), x.call, "none"), tested=np.isfinite(x.z)))


def _tvc(path: Path, method: str) -> pd.DataFrame:
    x = pd.read_csv(path/"genes.tsv", sep="\t", low_memory=False)
    tested = x.status.eq("tested") | (x["in_family"].astype(str).eq("True") if "in_family" in x else False)
    p = pd.to_numeric(x.p_value, errors="coerce").clip(lower=1e-300); lp = -np.log10(p)
    if method == "gumbel":
        effect = pd.to_numeric(x.log2_scaled_ratio, errors="coerce"); score = -lp
        call = np.where(x.call.eq("depleted run"), "depleted", "none")
    else:
        effect = pd.to_numeric(x.log2FC, errors="coerce"); score = np.sign(effect)*lp
        call = np.where(x.call.isin(["depleted", "enriched"]), x.call, "none")
    return pd.DataFrame(dict(background=x.background.astype(str), condition=x.condition.astype(str), ORF=x.ORF.astype(str),
                             Name=x.Name.astype(str), effect=effect, score=score, call=call, tested=tested & np.isfinite(score)))


def read_method(path: Path, method: str) -> pd.DataFrame:
    """One method's genes in the common shape: background, condition, ORF, Name, effect, score, call, tested. Gene IDs
    are matched without a GFF type prefix ("gene:SPAC1002.01" in resampling output from PomBase GFFs)."""
    x = _zscore(Path(path)) if method == "zscore" else _tvc(Path(path), method)
    x["ORF"] = x.ORF.str.replace(r"^gene:", "", regex=True); return x


def read_published(path: Path) -> pd.DataFrame:
    """Published calls: a TSV with ORF, call (depleted / enriched / none) and optional effect, score, background,
    condition. Without background/condition the calls apply to every group."""
    x = pd.read_csv(path, sep="\t")
    out = pd.DataFrame(dict(ORF=x.ORF.astype(str), Name=x.get("Name", x.ORF).astype(str), effect=pd.to_numeric(x.get("effect"), errors="coerce"),
                            score=pd.to_numeric(x.get("score"), errors="coerce"), call=x.call.fillna("none").astype(str), tested=True))
    for k in ("background", "condition"):
        if k in x:out[k] = x[k].astype(str)
    return out


def _spearman(a: pd.Series, b: pd.Series) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    return float(a[ok].rank().corr(b[ok].rank())) if ok.sum() >= 3 else float("nan")


def _set_stats(a: set, b: set) -> dict:
    both = len(a & b)
    return dict(n_a=len(a), n_b=len(b), n_both=both, a_recovered_by_b=both/len(a) if a else float("nan"),
                b_recovered_by_a=both/len(b) if b else float("nan"))


def agreement(methods: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """methods: name -> read_method()/read_published() frame. Returns (genes, pairs, sets)."""
    names = list(methods)
    groups = sorted({(b, c) for n in names if "background" in methods[n] for b, c in zip(methods[n].background, methods[n].condition)})
    genes_rows, pair_rows, set_rows = [], [], []
    for b, c in groups:
        per = {}
        for n in names:
            x = methods[n]
            per[n] = (x[(x.background == b) & (x.condition == c)] if "background" in x else x).drop_duplicates("ORF").set_index("ORF")
        orfs = sorted(set().union(*[set(v.index) for v in per.values()]))
        g = pd.DataFrame(index=pd.Index(orfs, name="ORF"))
        g["Name"] = pd.concat([v.Name for v in per.values()]).groupby(level=0).first().reindex(orfs)
        for n, v in per.items():
            g[f"{n}_effect"] = v.effect.reindex(orfs); g[f"{n}_score"] = v.score.reindex(orfs)
            g[f"{n}_call"] = v.call.reindex(orfs).fillna("not tested"); g[f"{n}_tested"] = v.tested.reindex(orfs).fillna(False).astype(bool)
        g.insert(0, "condition", c); g.insert(0, "background", b); g.insert(0, "group", f"{b}__{c}")
        genes_rows.append(g.reset_index())
        for m1, m2 in itertools.combinations(names, 2):
            both = g[g[f"{m1}_tested"] & g[f"{m2}_tested"]]
            row = dict(group=f"{b}__{c}", background=b, condition=c, method_a=m1, method_b=m2, genes_tested_by_both=len(both),
                       spearman_effect=_spearman(both[f"{m1}_effect"], both[f"{m2}_effect"]),
                       spearman_score=_spearman(both[f"{m1}_score"], both[f"{m2}_score"]))
            for d in ("depleted", "enriched"):
                if d == "enriched" and (m1 in DEPLETION_ONLY or m2 in DEPLETION_ONLY):
                    row |= {f"{d}_{k}": float("nan") for k in ("n_a", "n_b", "n_both", "a_recovered_by_b", "b_recovered_by_a")}; continue
                st = _set_stats(set(both.index[both[f"{m1}_call"] == d]), set(both.index[both[f"{m2}_call"] == d]))
                row |= {f"{d}_{k}": v for k, v in st.items()}
            pair_rows.append(row)
        for d in ("depleted", "enriched"):
            use = [n for n in names if not (d == "enriched" and n in DEPLETION_ONLY)]
            if len(use) < 2:continue
            tested = g[np.logical_and.reduce([g[f"{n}_tested"] for n in use])]
            member = pd.DataFrame({n: tested[f"{n}_call"].eq(d) for n in use})
            key = member.apply(lambda r: "&".join(n for n in use if r[n]), axis=1)
            for k, count in key[key != ""].value_counts().items():
                set_rows.append(dict(group=f"{b}__{c}", background=b, condition=c, direction=d, methods=k, n_methods=k.count("&")+1, genes=int(count)))
    genes = pd.concat(genes_rows, ignore_index=True) if genes_rows else pd.DataFrame()
    return genes, pd.DataFrame(pair_rows), pd.DataFrame(set_rows)


def run(inputs: dict[str, Path], output: Path, published: Path | None = None, sources: dict | None = None) -> dict:
    """inputs: method -> result folder (z-score treated_vs_parent folder, resampling or Gumbel run folder)."""
    methods = {m: read_method(p, m) for m, p in inputs.items()}
    if published is not None:methods["published"] = read_published(published)
    if len(methods) < 2:raise SystemExit("Agreement needs results of at least two methods.")
    genes, pairs, sets = agreement(methods)
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    genes.to_csv(out/"genes.tsv", sep="\t", index=False); pairs.to_csv(out/"pairs.tsv", sep="\t", index=False); sets.to_csv(out/"sets.tsv", sep="\t", index=False)
    manifest = dict(status="success", contract=VERSION, kind="tvc_agreement", methods=list(methods), labels={m: LABELS.get(m, m) for m in methods},
                    inputs={m: str(p) for m, p in inputs.items()} | ({"published": str(published)} if published else {}), sources=sources or {},
                    groups=sorted(genes.group.unique().tolist()) if len(genes) else [], depletion_only=sorted(DEPLETION_ONLY & set(methods)),
                    tables=["genes.tsv", "pairs.tsv", "sets.tsv"])
    (out/"manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest
