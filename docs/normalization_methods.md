# Selectable normalization methods

YTAB is developing a shared normalization interface, not replacing every workflow
with TTR. Users select a method for a stated analysis purpose. The future single
normalization tab will use the same machine-readable method descriptions as the
CLI: status, input level, calculation, potential use, cautions and diagnostics.
The app's Normalization tab now uses this interface (see "App Normalization tab" below). The treated-vs-control
z-score defaults to total reads (CPM) and accepts any method; the hybrid essentiality model's gene features can use
any method, MidLC being one of them (see "Hybrid gene features" below).

## Availability and input scope

| Method | State | What it changes / important caveat |
|---|---|---|
| Unnormalized baseline (`nonorm`) | Available in shared site runner | Preserves raw counts; useful for evaluating changes, not an automatic recommendation for unequal-depth comparisons. |
| TTR (`ttr`) | Available in shared site runner | Density × trimmed nonzero-mean scaling; robust scalar estimate, not a correction for all skew or absent mutants. |
| Total-read scaling (`totreads`) | Available in shared site runner | Equalizes totals to the group's average total; sensitive to dominant counts and composition changes. |
| Nonzero-mean scaling (`nzmean`) | Available in shared site runner | Equalizes means among occupied sites; does not adjust density differences. |
| CPM (`cpm`) | Z-score default (`ytab fitness`) | Gene-level total-read scaling (the same scaling as `totreads`, at gene level); listed for reference, not a separate site method. |
| MidLC common saturation (`midlc`) | Available in shared site runner | Binomial thinning of every library to the group's smallest depth-to-MidLC ratio (MidLC: depth with half the occupied sites expected, solved exactly); downsampling only, whole-number counts. A strongly selected library looks highly saturated and is thinned most, so it suits unselected libraries (parents) better than treated–control comparisons. |
| Zero-preserving positive-count quantile (`quantile`) | Available, experimental | Hermes-specific adaptation: preserves zeros, averages tied rank targets and transforms positive counts within explicit groups; not a literal TRANSIT dense-rank port or an established improvement. |
| Beta-geometric correction (`betageom`) | Available, experimental | Fits and transforms heavy-tailed site counts toward a geometric target; optional explicit final TTR. Catalogue-dependent fit; can compress genuine enrichment. |
| Adaptive BGC (`aBGC`) | Available, experimental | Searches positive-count trimming fractions using a maximum geometric-quantile discrepancy; exposes candidate validity, tail convention and simulation sensitivity. Not biological optimization. |
| Zero-inflated NB (`zinfnb`) | Available, experimental | Proper ZINB likelihood, correcting the inspected TRANSIT zero term; scales by inverse fitted NB mean. Catalogue- and model-dependent; convergence is not biological validation. |
| Empirical histogram, geometric quantile matching | Planned | Require separate adaptation, parameter review and validation before becoming selectable. |

Descriptions are guidance on assumptions, **not automatic method recommendations**.
Dataset evaluations characterize trade-offs and suitable data characteristics;
they do not establish a preferred method. Users choose. Any later recommendations
require broader, method-appropriate validation after the adaptations are complete.
Methods marked `planned` or `existing_workflow` cannot be executed by the shared
site-level runner. The `normalization` command (parent-only MidLC thinning written as hit files) now uses the same
rule as `midlc`: automatic target = common saturation of the selected parents, MidLC computed exactly. Its earlier
automatic targets were read counts passed where a multiple of MidLC is expected (e.g. `T409600` meant
409,600 × MidLC), so those runs did not thin at all; such targets are now rejected with an explanation.

## Inspect and choose

```sh
python scripts/local/ytab.py normalize --list-methods
python scripts/local/ytab.py normalize --describe ttr
python scripts/local/ytab.py normalize --describe totreads

python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_our_parents_ttr_v1 \
  --output output/normalization/gale_our_method_choices_v1 \
  --group-by library_id \
  --method nonorm --method ttr --method totreads --method nzmean
```

There is no default method. The user must supply `--method` and `--group-by`.
Repeat `--method` to compare methods; repeat `--observation PROJECT::SAMPLE` to
select observations explicitly. Available grouping keys are `library_id`,
`project`, `dataset` and `cohort`. The last intentionally treats the explicit
selection as one scaling group and must not be mistaken for independent-pool
stratification. `library_id` retains the independent-pool identities from metadata.

This example evaluates existing **parent observations**, not fitness treatment
effects. It does not add treatment samples to the parent catalogue or feed these
outputs into an analysis that already performs CPM. Experimental matrices now
have a separate opt-in runner below. Downstream method-aware fitness analysis
and the single normalization UI remain subsequent integration steps.

## Evaluate an experiment against its controls

```sh
python scripts/local/ytab.py normalize-experiment \
  --project-config output/projects/Gale_2020/config/project.yaml \
  --catalogue output/normalization/gale_our_parents_ttr_v1 \
  --output output/normalization/gale_flu_ttr_validation_v1 \
  --method nonorm --method ttr
```

This generic runner reads included comparisons from the project's
`config/comparison_design.csv` and resolves sample identities, display names,
backgrounds and pools from its sample sheet. `parent_sample` in the comparison
design denotes the comparison control, not necessarily an original induced
parent library. Only design-selected control/treatment samples enter the count
matrix. Matched samples must share background and pool; control–control
diagnostics are computed only within that same background/pool stratum.

The raw matrix starts with the supplied empirical parent catalogue, retaining
its row ordering, and appends any experiment-observed coordinates absent from
it. The original parent catalogue is never modified. This preserves all observed
counts instead of silently dropping treatment/control-only sites; it is an
explicit **experiment-expanded empirical catalogue**, not a saturated target-site
universe. All four assay samples share the same rows and have zeros where
undetected. Adding rows that are zero in every selected sample rescales absolute
TTR factors equally, without changing relative factors or log ratios.

Outputs include the raw matrix, coordinate/observation metadata, independent
per-method results, factor tables, comparison metrics and plots. Site log ratios
are evaluated at shared positive sites without pseudocounts. Left-only and
right-only site counts are reported separately: shared-site metrics alone cannot
describe dropout or establish improved fitness calls. The plots use identical
raw geometric-mean abundance on the horizontal axis across methods. Scalar
normalization cannot change ranks, zero patterns or log-ratio spread (IQR); it
can shift overall log-ratio centering. Control/control centering is an offset
diagnostic, not proof that real treatment effects are correctly estimated.

No gene-level significance testing, fitness calls, independent-induction claims,
or downstream CPM re-normalization occur in this validation runner. The three
original Gale induced pools remain the parent-library catalogue/classifier
inputs; Gale's two controls and two treated Cg1-derived cultures are the assay
validation samples. See the generated
[Gale validation report](../output/normalization/gale_flu_ttr_validation_v1/INTERPRETATION.md).

## Calculation and grouping

### Gene-level validation after the experimental site run

```sh
python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm
```

This separate diagnostic command compares raw gene counts, whole-library TTR
and the current feature-total CPM method. It reuses the existing feature database
and SummaryTable's annotation/counting policy, including overlapping-feature
behavior. Raw counts reconstructed from the experimental site matrix must match
every feature in each existing RDF_1 summary exactly; otherwise the command
fails before producing results. Normalized-site aggregation must also agree
with scalar scaling after raw-gene aggregation. TTR and CPM each start from raw
counts independently: CPM is never applied after TTR.

Optional `--existing-fitness` verifies the complete existing CPM export and its
default log2 fold changes with a 0.5 CPM pseudocount. This is a parity check, not
a new fit of the existing control-noise/z-score model. Other existing pseudocount
settings fail this explicit default-setting check rather than being silently
treated as equivalent.

Cross-method diagnostics deliberately use **no pseudocount** and compare shared
positive genes, retaining zero/gain/dropout flags separately. A fixed 0.5 added
to raw, TTR and CPM counts would have different meanings because those counts
have different units. Count filtering is always on raw support, with sensitivity
checks at 1, 10, 20 and 50 counts (plus the selected `--min-count`, default 20).
The optional `--fold-threshold` defaults to twofold and is descriptive only:
control excursions are not p-values, significant hits, or an estimated FDR.

By default, treatment agreement and control–control diagnostics stay within the
same library. Experiments with independent pools can explicitly request
`--compare-independent-pools`. This adds effect-agreement comparisons only when
the design specifies the same nonempty `background` and `comparison_group`.
Mock–mock contrasts from different pools are labeled `between_pool_controls`,
not technical control noise; agreement rows carry `between_independent_pools`.
No pools or backgrounds are merged to estimate a normalization factor.

For the zinc validation:

```sh
python scripts/local/ytab.py normalize-experiment \
  --project-config output/projects/Zn_toxicity_screen/config/project.yaml \
  --catalogue output/normalization/gale_our_parents_ttr_v1 \
  --output output/normalization/zn_ttr_validation_v1 \
  --method nonorm --method ttr

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/zn_ttr_validation_v1 \
  --output output/normalization/zn_gene_validation_v1 \
  --existing-fitness output/projects/Zn_toxicity_screen/treated_vs_parent/Zn_1_5mM_vs_mock \
  --compare-independent-pools
```

The [zinc report](../output/normalization/zn_gene_validation_v1/INTERPRETATION.md)
records the different control-replicate structure and mixed method results.

Outputs include per-method gene matrices, annotations, count-stratified
control/treatment metrics, treatment-effect agreement within library strata,
gene effects with evidence flags, plots and hash-based provenance. There is no
UI integration or change to existing fitness/classifier defaults. Results of the
first four-sample Gale evaluation are in the
[gene-level report](../output/normalization/gale_flu_gene_validation_v1/INTERPRETATION.md).

### Shared site-method calculations

Every method starts independently from `raw_counts.npz`, never from another
method's output. Counts must be finite, nonnegative integers; empty samples fail
explicitly rather than receiving an invented factor. All seven available methods
preserve the zero pattern. Duplicate selections and incomplete catalogues fail.

The four scalar methods return sample factors. Experimental `quantile` and
`betageom`/`aBGC` are nonlinear and leave `factor` blank rather than inventing a scalar. See
[quantile mathematics and evaluation](quantile_normalization.md). Its normalized
sites must be summed into genes **after** transformation, never approximated by
scaling raw gene totals. `normalize-genes --site-normalization DIRECTORY` accepts
the matching site output and checks identical raw-matrix/coordinate hashes,
ordered observations, library grouping and zero preservation before aggregation.
BGC's own post-TTR step is explicitly recorded, not a subsequent CPM workflow.
Its numerical engine, model assumptions and catalogue sensitivity are described in
[BGC mathematics and guidance](betageom_normalization.md).
Adaptive BGC has a separate [principle, parameter-search and sensitivity guide](adaptive_bgc_normalization.md).
The [ZINB guide](zinfnb_normalization.md) documents its proper mixture likelihood,
profiled zero component, numerical bounds, per-sample scalar fitting and explicit
departure from TRANSIT. It does not equalize group totals or use a group target.

For total-read scaling, `factor_j = mean_group(total_support) / total_support_j`.
For nonzero-mean scaling, `factor_j = mean_group(nonzero_mean) / nonzero_mean_j`.
These follow the inspected TRANSIT implementations. Both have factor 1 in a
single-sample group: there is no other sample from which to estimate a relative
group scale. TTR instead uses the explicit target, default 100, and a common site
universe; its default trim is 5% per tail. See [TTR details](ttr_normalization.md).
The dedicated [nzmean guide](nzmean_normalization.md) documents positive-site
scaling, grouping, outlier sensitivity and the Gale/zinc evaluations. Gene
validation accepts `nzmean` site results, checks the factor formula independently,
and verifies that scaling sites before aggregation equals scaling raw gene sums.
It does not recompute a nonzero mean across genes or apply CPM afterward.
The [totreads guide](totreads_normalization.md) documents whole-site total scaling
and its separate Gale/zinc evaluations. Gene validation also accepts `totreads`,
verifies the site-derived factor and aggregation identity, and distinguishes its
whole-site denominator from the existing gene-level CPM denominator.

Output count units differ across methods. Larger normalized totals do not mean
better normalization. Compare relative treatment/control results, matched-control
behavior, and agreement within appropriate pool strata—not absolute output counts
between different normalization methods.

## Outputs and provenance

Each selected method has its own folder with `normalized_counts.npz`, `factors.tsv`
and `method.json`. Matrix rows retain the original catalogue's `sites.tsv.gz`
ordering; columns follow the exported `observations.tsv` order. A comparison
table and run manifest record selection, grouping, method parameters, input hashes,
reference identity, implementation hashes and NumPy/SciPy versions. New output
directories are required; input catalogues and existing results are not overwritten.

This interface does not alter the fitness/essentiality run caches, accepted
biological results, empirical site catalogue, or scientific defaults. It does not
infer viable mutant abundance, introduce motif weights, or claim saturation.

## App Normalization tab and `ytab normalize-project`

The app's **Normalization → Site Normalization** tab runs one command for the open project:

```sh
ytab normalize-project --project-config output/projects/PROJECT/config/project.yaml \
    --output output/projects/PROJECT/normalization/RUN          # every available method; --method M to restrict
```

It runs the existing stages in order:
1. a parent site catalogue from every included parent/control sample;
2. the experiment matrix from `comparison_design.csv`, with each method applied to the same raw matrix within each
   library; unnormalized counts are always added as the baseline;
3. plot summaries in `plots/`: positive-count percentiles, matched-pair diagnostics, and each method's native
   before/after view (`native_*.tsv`): TRANSIT-style per-sample statistics (density, mean, nonzero mean and median,
   max, total, skew, kurtosis); for scaling methods the statistic they equalize (TTR: density × trimmed nonzero mean;
   total reads; nonzero mean; ZINB: fitted NB mean), before and after; for BGC, aBGC and unnormalized counts the
   quantiles of the nonzero counts against a geometric distribution fitted to them, and the ranked counts; for MidLC the
   exact saturation curve of each library with its MidLC point and the point after thinning.

Normalization is per sample: each library gets its own factor (scaling methods) or its own fitted transform (BGC,
aBGC, quantile), within its background and pool. `--refresh-summaries RUN` adds the summaries to a run made before
they existed (the app does this on first view).

In the app, a run computes every method. The **Show methods** boxes (full names, Unnormalized included) choose which
methods are shown. View **Each normalization, before and after** gives one panel per ticked method with its own plot
and axes (the method's native view, above); **Replicate agreement** is the one view comparing methods side by side
(matched pairs: median log2 ratio and IQR at shared sites, short method names at 40°). The per-sample statistics table
follows TRANSIT's quality-control table.

`RUN/experiment/methods` is the input for `ytab multicondition`, `ytab genetic-interactions` and `ytab resampling`.

**Suggested method in the app: BGC.** The CLI still has no default. The app shows BGC (and TTR for comparison) first
because, in the pre-registered Hermes benchmarks, BGC then TTR with recentering kept false discoveries lowest
(resampling v3: FDR 7.0% vs 12.4% with TTR), and the multi-condition and GI tests were validated with BGC. Unnormalized
counts can be shown alongside for reference.

**Fix (2026-10-20).** `normalize-experiment` computed matched-pair diagnostics from per-sample scale factors, which
nonlinear methods (quantile, BGC, aBGC) do not have, so it failed for them (`TypeError`). Diagnostics now compare
normalized values directly. For scalar methods this equals the previous formula (tested).

## Hybrid gene features from any method (`ytab hybrid-inputs`)

MidLC is no longer special to the classifier branch. The hybrid essentiality model's gene features (hits, reads,
freedom and neighbourhood indices, longest free interval, and the current random forest's scores) are built from the
parents normalized with the chosen method, default BGC:

```sh
ytab hybrid-inputs --project-config output/projects/PROJECT/config/project.yaml --normalization betageom --samples P1,P2,P3
```

`ytab essentiality --labels ...` runs it automatically with the HMM's own parent libraries (`--normalization`, default
`betageom`; `--rf` still takes an existing prediction table). Steps: the parents are normalized with the method from
the newest Site Normalization run whose parent catalogue holds them; written as hit files under
`sample_normalization/<method>__<run>__<parents>/`; then combined, summarized (SummaryTable) and scored by the
classifier exactly as for MidLC. Results for the same method, run and parents are reused untouched.

- **Whole numbers.** SummaryTable reads whole-number counts. A normalized count is rounded to the nearest whole number
  and a site with insertions keeps at least 1, so site-based features match the data and read-based features carry
  the method's scaling. A site's W and C strand rows share its normalized count in proportion to their raw reads.
  Each written file is checked against the catalogue's raw counts; a mismatch stops the step.
- **MidLC.** Common saturation over all chosen parents (one group). Its counts are already whole numbers; a site's
  thinned count is shared between its strand rows by largest remainder and thinned-out sites are dropped (not raised
  to 1). On Gale's five parents the least saturated library (Cg3 parent, 13× MidLC) is unchanged and the most
  saturated (Cg1 parent, 73×) keeps 18% of its reads.
- **Gale check (2026-10-21, descriptive, not pre-registered).** Same HMM fit, parents and labels; silver/gold hybrid
  AUC: real MidLC (common saturation) 0.947/0.972, BGC 0.947/0.972, TTR 0.947/0.972, the earlier non-thinning
  "MidLC" target 0.947/0.974 (all within 0.002 of real MidLC, bootstrap intervals include 0). The current RF, trained
  on S. cerevisiae and not retrained, loses 0.002–0.007 with BGC/TTR. Details:
  `output/essentiality/hybrid/normalization_comparison_real_midlc_2026-10-21/RESULT.md`.

## Validation

```sh
PYTHONPATH=src python -m pytest -q \
  tests/test_normalization_genes.py tests/test_normalization_experiment.py tests/test_normalization_choices.py \
  tests/test_zinfnb_normalization.py tests/test_totreads_normalization.py tests/test_nzmean_normalization.py tests/test_adaptive_bgc.py tests/test_betageom_normalization.py tests/test_quantile_normalization.py tests/test_ttr.py tests/test_cli.py
```

The tests check TRANSIT numerical parity when the optional source checkout is
available, pool-stratified scaling, raw/zero preservation, independent method
execution, selected column order, provenance and rejection of unavailable methods.
