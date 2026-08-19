# Fitness candidate-ranking inventory

## Executive summary

The scientific calculation is in `scripts/ytab_treated_vs_parent_screen.R`. The Python runner (`src/ytab/pipeline/treated_vs_parent_runner.py`) validates inputs, runs that R script, and copies its summary output; it does not add a ranking algorithm.

Each feature/pool row is calculated from CPM, treated/control log2FC, a local parent-parent noise-model z-score, and a stored z-threshold call. The combined summary averages log2FC and z across valid pool rows and records directional pool counts. The legacy `top100` files then sort the combined summary by mean log2FC and take 100 rows without filtering by call.

The current Shiny MA plot is separate from the legacy `top100` files. Combined mode uses `treated_vs_parent.summary_by_feature.csv`; individual mode filters `treated_vs_parent.by_pool.log2fc_z.csv` to the selected contrast. Both use stored call categories for eligibility and log2FC for direction-specific ordering.

Strong single-pool effects are retained, but combined averaging can dilute them. The legacy top100 files can also contain `final_call = none` because they are not call-filtered.

## Files inspected

- `scripts/ytab_treated_vs_parent_screen.R`
  - `read_sample_feature_table()`
  - `make_control_running()`
  - `compute_contrast()`
  - `summary_results` aggregation
  - `top_depleted` / `top_enriched` creation
  - legacy MA plot generation (`p_ma_all`, `base_ma`, `sig_labels`)
- `src/ytab/pipeline/treated_vs_parent_runner.py`
  - `build_treated_vs_parent_command()`
  - `choose_primary_treated_vs_parent_table()`
  - `create_stable_treated_vs_parent_outputs()`
  - schema/cache validation
- `app/shiny/R/fitness_plot_utils.R`
  - `fitness_ma_data()`
  - `fitness_ma_pair_choices()`
  - `fitness_ma_highlight_rows()`
  - `fitness_ma_annotation_rows()`
- `app/shiny/app.R` and `app/shiny/R/ui_fitness.R`
  - Fitness selectors, MA renderer, and download handlers
- `app/shiny/R/fitness_result_state.R`
- `scripts/local/ytab_run_treated_vs_parent.py`
- `scripts/local/ytab_smoke_fitness_ui.R`

No alternate ranking implementation was found in `pipeline/hpc/` or the Zn local wrappers.

## Output tables inspected

Inspected both:

- `output/projects/Zn_toxicity_screen/treated_vs_parent/Zn_1_5mM_vs_mock/tables/`
- `output/projects/H2O2_screen_v1/treated_vs_parent/H2O2_vs_parent/tables/`

| Table | Rows | Relevant columns |
|---|---:|---|
| `treated_vs_parent.by_pool.log2fc_z.csv` | 20,936 | `feature_id`, `contrast`, `pool`, `parent_cpm`, `treated_cpm`, `A_exp`, `log2FC`, `sd_local`, `z`, `Zthr`, `call` |
| `treated_vs_parent.summary_by_feature.csv` | 5,234 | `mean_log2FC`, `median_log2FC`, `sd_log2FC`, `mean_z`, `max_abs_z`, `n_pairs`, `n_enriched_pairs`, `n_depleted_pairs`, pool lists, `final_call` |
| `top100_depleted_in_treated.csv` | 100 | Summary columns; mean-log2FC ascending |
| `top100_enriched_in_treated.csv` | 100 | Summary columns; mean-log2FC descending |
| `feature_reads_cpm.long.csv` | 41,872 | `feature_id`, `sample`, `reads`, `total_feature_reads`, `cpm` |
| `background_specific_z_thresholds.csv` | 2 | background-specific quantile and z threshold |
| `background_noise_model_coefficients.csv` | 6 | fitted local-noise parameters |

Both projects have four valid pool comparisons; every summary feature has `n_pairs = 4`.

## By-pool calculation

`compute_contrast()` creates one row per feature and treated/control pool comparison:

```text
A_exp = (treated_cpm + parent_cpm) / 2
log2FC = log2((treated_cpm + pseudocount) / (parent_cpm + pseudocount))
z = log2FC / sd_from_fit(parent-parent local-noise model, A_exp)
```

Calls are assigned from the stored local z threshold:

- `z >= z_thr` → `enriched`
- `z <= -z_thr` → `depleted`
- otherwise → `none`

The table has no candidate rank column and no explicit low-support flag. It retains CPM/read-support fields, local z fields, pool/contrast metadata, and annotations. `read_sample_feature_table()` only requires a nonempty feature ID and numeric reads value; no later minimum CPM/read-count filter is applied before ranking.

## Summary-by-feature construction

The `summary_results` block groups by `feature_id` and calculates:

- `mean_log2FC`, `median_log2FC`, `sd_log2FC`
- `mean_z`, `max_abs_z`
- `n_pairs`
- `n_enriched_pairs`, `n_depleted_pairs`
- `pools_enriched`, `pools_depleted`

`final_call` is:

1. ≥2 enriched and no depleted: `consistently_enriched`
2. ≥2 depleted and no enriched: `consistently_depleted`
3. both directions present: `mixed`
4. exactly one enriched: `single_pool_enriched`
5. exactly one depleted: `single_pool_depleted`
6. otherwise: `none`

Single-pool calls are retained. Means/medians ignore missing values rather than penalizing them. The current output has four valid rows per feature.

## Top100 generation

The exact source is `scripts/ytab_treated_vs_parent_screen.R`:

```r
top_depleted <- summary_results %>% arrange(mean_log2FC) %>% slice_head(n = 100)
top_enriched <- summary_results %>% arrange(desc(mean_log2FC)) %>% slice_head(n = 100)
```

Therefore these files:

- use the combined summary, not the by-pool table;
- rank only by `mean_log2FC` (ascending/descending);
- do not use CPM, reads, z-score magnitude, pool count, or consistency as ranking keys;
- do not filter by `final_call`.

Observed call composition:

| Project/file | Single-pool | Consistent | Mixed | None |
|---|---:|---:|---:|---:|
| Zn depleted | 7 | 0 | 0 | 93 |
| Zn enriched | 20 | 1 | 0 | 79 |
| H2O2 depleted | 40 | 13 | 0 | 47 |
| H2O2 enriched | 30 | 1 | 0 | 69 |

## Current Shiny MA ranking

### Combined mode

`fitness_ma_data(result, "combined")` reads `treated_vs_parent.summary_by_feature.csv` and uses existing CPM rows to derive display abundance. `fitness_ma_highlight_rows()` selects rows whose stored `final_call` contains `depleted` or `enriched`, then sorts by `mean_log2FC` (depleted ascending, enriched descending). `single_pool_depleted` and `single_pool_enriched` are eligible; `mixed` is not matched by the current string test.

### Individual-pair mode

`fitness_ma_data(result, "individual", pair)` filters `treated_vs_parent.by_pool.log2fc_z.csv` to the selected `contrast`. It uses that row’s `parent_cpm`, `treated_cpm`, pair-specific `log2FC`, and pair-specific stored `call`. The same helper ranks the selected pair by its own log2FC. It does not use combined rankings or top100 files.

### Annotation

`fitness_ma_annotation_rows()` matches custom feature IDs/names against the current plotted table. Custom annotation bypasses ranking for matched entries; top-hit annotation uses only the current highlight result.

## Single-pool answer

| Path | Answer | Evidence |
|---|---|---|
| `top100_depleted_in_treated.csv` | **YES, but diluted and not call-filtered** | Single-pool entries occur; sorting is combined `mean_log2FC` only. |
| `top100_enriched_in_treated.csv` | **YES, but diluted and not call-filtered** | Single-pool entries occur; sorting is combined `mean_log2FC` only. |
| Combined MA | **YES, partially** | Single-pool final calls are eligible, but four-pool averaging can dilute the effect; mixed calls are excluded by string matching. |
| Individual pair MA | **YES** | The selected pair is filtered first and ranked directly by that pair’s log2FC/call. |

## Pool/replicate consistency gate

There is a consistency classification, not a universal ranking gate. `n_enriched_pairs`, `n_depleted_pairs`, pool lists, and `final_call` record directional support. The top100 files ignore `final_call`; combined MA accepts single-pool calls but not `mixed`; individual MA has no cross-pool requirement. There is no pool correlation gate, all-pools-agree gate, minimum support-count gate, or minimum number of valid comparisons in the ranking path. Parent-parent controls calibrate local z thresholds, but are not a candidate-consistency gate.

## Plain-language summary

YTAB calculates CPM-based treated/control log2FC per pool, calibrates a local z-score from parent-parent controls, and stores a depleted/enriched/none call. Legacy top100 files average log2FC across pools and take the 100 most negative or positive features, regardless of call. Combined MA uses the summary mean log2FC plus stored call eligibility. Individual MA uses the selected pair’s own log2FC plus stored pair call. CPM/read support affects abundance and z calibration but does not rank candidates. z-score determines calls but is not a secondary rank key. Pool support is recorded/classified but not required for top100 ranking; strong single-pool effects can remain visible, although combined means may dilute them.

## Future ranking modes supported by existing columns

| Mode | Current status | Existing columns |
|---|---|---|
| Combined evidence | Partially implemented | `mean_log2FC`, `median_log2FC`, `mean_z`, `max_abs_z`, `n_pairs`, directional counts, pool lists, `final_call` |
| Strongest single-pool effect | Partially implemented | By-pool `log2FC`, `z`, CPM, `A_exp`, `call`, `pool`, `contrast`; summary has `max_abs_z` but not max-absolute-log2FC or pool-of-maximum-effect |
| Multi-pool consensus | Partially implemented | Directional counts, pool lists, `final_call`; consensus labels exist but are not required by top100 ranking |

## Gaps and design questions

- The `top100` names imply called candidates, but the files contain many `final_call = none` rows.
- Combined MA uses a string test, so `mixed` features are excluded from either direction.
- Combined mean log2FC does not expose strongest-pool effect or dilution.
- `max_abs_z` has no associated pool column, and there is no `max_abs_log2FC` summary column.
- No explicit abundance/read-support floor is applied at ranking time.
- `treated_vs_parent_results.csv` is a stable copy of the summary table, not an independent rank table.

Questions for a future design decision: should top100 be call-filtered; should MA offer mean-effect/strongest-pool/consensus modes; should mixed calls be selectable; should abundance support be a display-only or ranking filter; and should single-pool candidates receive a separate view?

## Audit scope

This report was produced by reading existing source and Zn/H2O2 output tables. No code, UI, scientific table, project output, or pipeline result was modified or recomputed.
