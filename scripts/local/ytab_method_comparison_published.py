#!/usr/bin/env python3
"""Method-comparison study, step 2: published calls of the Cunningham-lab C. glabrata screens in the shape
ytab.fitness.agreement reads (ORF, Name, call, score, effect, background, condition), one TSV per dataset in
output/comparison/<dataset>/published.tsv.

Sources
- Gale et al. 2020 (G3): per-gene published z and log2FC (output/essentiality/control_conditioned_gumbel/Gale_2020/
  publication_comparison.tsv, from the paper's supplement).
- Nickels & Cunningham 2025 (mSphere) Table S1 (resources/species/glabrata/tnseq_data/nickels_2025_echinocandins/
  msphere.00270-25-s0003.xlsx): z-scores of every gene in every strain x drug x dose, which also carries the
  Nickels et al. 2024 CBS138 micafungin screens (its "ref. 2" columns).

Published call (one rule for every study): z <= -3 depleted, z >= 3 enriched (the |z| threshold of Gale's publication
classes). Genes are CBS138 ORFs (CAGL...) in Table S1; BG2 datasets are matched through the BG2 ortholog table.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT/"output/comparison"
S1 = ROOT/"resources/species/glabrata/tnseq_data/nickels_2025_echinocandins/msphere.00270-25-s0003.xlsx"
BG2_MAP = ROOT/"resources/species/glabrata/orthologs/C_glabrata_BG2_S_cerevisiae_orthologs.txt"
GALE = ROOT/"output/essentiality/control_conditioned_gumbel/Gale_2020/publication_comparison.tsv"
THRESHOLD = 3.0
STRAINS = {"BG2u (BG14)": "BG2", "pdr1∆ (CGM1094)": "pdr1_KO", "PDR1-GOF (TJN88)": "PDR1_GOF", "2001u (TJN32)": "CBS138_2001", "CBS138u (TJN69)": "CBS138"}
DRUGS = {"M": "micafungin", "I": "ibrexafungerp", "R": "rezafungin", "A": "anidulafungin", "C": "caspofungin"}
# which dataset each (background, drug) of Table S1 belongs to
DATASET = {("BG2", d): "Nickels_2025_glabrata" for d in DRUGS.values()} | {("pdr1_KO", d): "Nickels_2025_glabrata" for d in DRUGS.values()} | \
          {("PDR1_GOF", d): "Nickels_2025_glabrata" for d in DRUGS.values()} | {("CBS138_2001", "caspofungin"): "Nickels_2025_glabrata_CBS138",
          ("CBS138_2001", "micafungin"): "Nickels_2024_glabrata", ("CBS138", "micafungin"): "Nickels_2024_glabrata"}


def call(z): return np.where(z <= -THRESHOLD, "depleted", np.where(z >= THRESHOLD, "enriched", "none"))


def table_s1() -> pd.DataFrame:
    x = pd.read_excel(S1, header=None)
    strain = x.iloc[15].ffill(); dose = x.iloc[16]; body = x.iloc[18:]
    rows = []
    for j in range(9, x.shape[1]):
        code = str(dose.iloc[j]); s = STRAINS.get(str(strain.iloc[j]))
        if s is None or len(code) < 3 or code[-1] != "F" or code[-2] not in DRUGS:continue   # fluconazole (Gale) columns are skipped here
        drug = DRUGS[code[-2]]; ng = code[:-2]
        z = pd.to_numeric(body.iloc[:, j], errors="coerce")
        rows.append(pd.DataFrame(dict(CAGL=body.iloc[:, 1].astype(str), Name=body.iloc[:, 2].fillna(body.iloc[:, 4]).astype(str), score=z,
                                      background=s, condition=f"{drug}_{ng}ng", dataset=DATASET.get((s, drug), ""))))
    d = pd.concat(rows, ignore_index=True); d = d[np.isfinite(d.score) & d.dataset.ne("")]
    d["call"] = call(d.score.to_numpy()); d["effect"] = np.nan
    return d


def main():
    d = table_s1(); bg2 = pd.read_csv(BG2_MAP, sep="\t", dtype=str)
    to_bg2 = bg2.dropna(subset=["cglab_cagl_id"]).drop_duplicates("cglab_cagl_id").set_index("cglab_cagl_id").cglab_id
    for dataset, g in d.groupby("dataset"):
        g = g.copy()
        g["ORF"] = g.CAGL.map(to_bg2) if dataset == "Nickels_2025_glabrata" else g.CAGL
        lost = g.ORF.isna().mean(); g = g.dropna(subset=["ORF"])
        out = OUT/dataset; out.mkdir(parents=True, exist_ok=True)
        g[["ORF", "Name", "call", "score", "effect", "background", "condition"]].to_csv(out/"published.tsv", sep="\t", index=False)
        print(f"{dataset}: {g.groupby(['background', 'condition']).ngroups} published conditions, {len(g):,} gene rows"
              f"{f', {lost:.1%} of CBS138 genes without a BG2 ortholog' if lost else ''}")
    p = pd.read_csv(GALE, sep="\t"); p = p[np.isfinite(p.published_z)]
    out = OUT/"Gale_2020"; out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dict(ORF=p.ORF, Name=p.Name, call=call(p.published_z.to_numpy()), score=p.published_z, effect=p.published_log2fc,
                      background="BG14", condition="treated")).to_csv(out/"published.tsv", sep="\t", index=False)
    print(f"Gale_2020: 1 published condition, {len(p):,} genes")


if __name__ == "__main__":main()
