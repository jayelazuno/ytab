# Essentiality model menu

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
