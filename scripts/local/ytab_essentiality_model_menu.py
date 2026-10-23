#!/usr/bin/env python3
"""Side-by-side essentiality model menu (HMM, hybrid, current RF) for one set of libraries.

HMM only (no labels, e.g. a non-model organism):
    ytab_essentiality_model_menu.py --hmm genes.tsv --out DIR
With gene features and curated labels (adds the hybrid):
    ytab_essentiality_model_menu.py --hmm genes.tsv --features combined_feature_table.TAG.txt --labels labels.tsv --out DIR
--rf essentiality_predictions.csv adds the legacy classifier's score for comparison (and its gene features when
--features is not given).
--reproduce-gale: the recorded Gale 2020 Cg1-Cg3 inputs -> output/essentiality/model_menu/Gale_2020/.
Writes genes.tsv, summary.json and MODELS.md.
"""
import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import model_menu as mm
from ytab.normalization.ttr import sha256

DEFAULTS = dict(hmm=ROOT/"output/essentiality/hermes_hmm/benchmark_v2/gale_blinded/joint_run1/genes.tsv",
                # recorded input of the pre-registered reproduction (moved with the other non-thinning MidLC outputs)
                rf=ROOT/"archive/ui_reorganization/backup/old_midlc/projects/Gale_2020/classifier/T409600/essentiality_predictions.T409600.csv",
                labels=ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv",
                out=ROOT/"output/essentiality/model_menu/Gale_2020")

GUIDE = """# Essentiality model menu

`genes.tsv` has one row per gene with every available model side by side. Choose the model that fits the
organism and the evidence available; the methods, validation and limitations are in
`output/essentiality/hermes_hmm/hermes_hmm_methods.md`.

| model | columns | needs labels? | use when | main caveat |
|---|---|---|---|---|
| **HMM** (primary) | `hmm_score` (P_ES_run), `hmm_call` | no | always; the only option without essentiality labels | lower precision (accessibility cold spots, severe growth defects) |
| **hybrid** | `hybrid_score`, `hybrid_essential` | yes | labels exist for the species or a close relative | inherits the bias of its labels |
| **legacy RF** (comparison only) | `rf_existing_score`, `rf_existing_essential` | trained on *S. cerevisiae* | showing the limits of transferring a model-organism classifier | trained on another species' dense libraries; over-calls at low insertion density |

**Columns:**
- `hmm_call`: ES, GD, NE, GA (insertion-density classes; only ES means essential), or `uncertain`
  (< 50% of the trimmed CDS mappable).
- `hybrid_source`: `out-of-fold (training gene)` means the gene's label was used for training and its
  score comes from a model that did not see it; `trained model` means the gene was not in training.
- `agreement` (HMM ES call vs hybrid score ≥ 0.5):
  - `consensus essential`: both models; the most reliable calls.
  - `HMM only`: candidates for species-specific essentiality, or accessibility artefacts.
  - `hybrid only`: essential-like by gene features but without a confident essential run.
  - `not essential`: neither.
- `reference_label`, `reference_tier`: curated labels, shown for inspection only.

Thresholds are 0.5 for every model, as benchmarked. No setting was tuned to these labels.
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hmm", type=Path); ap.add_argument("--rf", type=Path); ap.add_argument("--labels", type=Path)
    ap.add_argument("--features", type=Path, help="SummaryTable feature table (gene features for the hybrid)")
    ap.add_argument("--out", type=Path); ap.add_argument("--train-tiers", default="silver")
    ap.add_argument("--reproduce-gale", action="store_true", help="Use the recorded Gale 2020 inputs (validation reproduction only)")
    a = ap.parse_args()
    if a.reproduce_gale:a.hmm, a.rf, a.labels, a.out = a.hmm or DEFAULTS["hmm"], a.rf or DEFAULTS["rf"], a.labels or DEFAULTS["labels"], a.out or DEFAULTS["out"]
    if a.hmm is None or a.out is None:raise SystemExit("give --hmm and --out (or --reproduce-gale)")
    hmm = pd.read_csv(a.hmm, sep="\t"); rf = pd.read_csv(a.rf) if a.rf else None
    features = pd.read_csv(a.features, skiprows=1) if a.features and open(a.features).readline().startswith("RDF") else (pd.read_csv(a.features) if a.features else None)
    labels = pd.read_csv(a.labels, sep="\t", dtype=str) if a.labels else None
    menu, summary = mm.build(hmm, rf, labels, tuple(a.train_tiers.split(",")), features=features)
    a.out.mkdir(parents=True, exist_ok=True)
    menu.to_csv(a.out/"genes.tsv", sep="\t", index=False)
    summary["inputs"] = {k:dict(path=str(p), sha256=sha256(p)) for k, p in (("hmm", a.hmm), ("features", a.features), ("rf", a.rf), ("labels", a.labels)) if p}
    summary["script_sha256"] = sha256(Path(__file__)); summary["module_sha256"] = sha256(Path(mm.__file__))
    (a.out/"summary.json").write_text(json.dumps(summary, indent=2, default=str)+"\n"); (a.out/"MODELS.md").write_text(GUIDE)
    print(json.dumps({k:v for k, v in summary.items() if k != "inputs"}, indent=1, default=str))


if __name__ == "__main__":main()
