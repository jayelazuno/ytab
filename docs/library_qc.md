# Library QC table (`ytab qc`)

One table with every per-library QC metric of a project, each flagged green / yellow / red (or grey: not flagged) and
commented, so users can judge their libraries and choose a normalization. It runs on its own and the app's
**Quality Control** tab only displays its output.

```sh
python scripts/local/ytab.py qc --project-config output/projects/PROJECT/config/project.yaml \
    --output output/projects/PROJECT/qc/RUN
```

Code: `src/ytab/qc/library_qc.py` (package module) and `scripts/local/ytab_qc.py` (CLI). Takes seconds.

## Outputs

| file | contents |
|---|---|
| `qc_table.tsv` | one row per library, one column per metric |
| `qc_flags.tsv` | one row per library × metric: value, flag, reason, source of the rule, metric comment, normalization hint |
| `qc_table.xlsx` | the table coloured by flag; each cell's note holds its reason and hint; header notes hold the metric comments; a Read me and a Flags sheet |
| `metrics.tsv` | the metric definitions (theme, label, bad direction, rule, comment) |
| `plots/histogram.tsv.gz`, `qq.tsv.gz`, `ranked.tsv.gz` | count histogram, QQ against a geometric fit and ranked counts for every library (adapted from TRANSIT's QC plots) |
| `manifest.json` | counts, catalogue size, genome size, the rules and their parameters |

## Metrics

- **Reads & Mapping:** reads sequenced, retained after Hermes cleanup, with 3' adaptor (Hermes cleanup record);
  mapped %, mapped with MAPQ above the threshold %, mean MAPQ (MapFastq).
- **Insertions** (from the hit file; strands combined; definitions adapted from TRANSIT's `tnseq_stats`, see `docs/tnseq_stats.md`):
  reads at insertion sites (total_cts), insertion sites, insertions per kb of genome, catalogue density, mean reads
  per catalogue site (mean_ct), NZmean, NZmedian, max reads at a site (and where), top site's share of reads, top 1%
  of sites' share, skewness, kurtosis, Pickands tail index, strand correlation (TPP FR_corr).
- **Genes** (SummaryTable): genes hit %, insertions per gene, insertions in genes %, insertion–neighbourhood
  Spearman correlation.
- **Saturation** (saved library diagnostics): MidLC estimate, depth / MidLC, reads in jackpot sites, jackpot sites.

The catalogue for density and mean_ct is the project's empirical site catalogue: every position with an insertion in
any included library. Missing inputs (e.g. no cleanup record for libraries mapped on the HPC) give an empty, grey
cell, never a zero.

## Flags

No cut-off is invented:

- **Rule of thumb** — only where published Tn-seq QC guidance is scale-free and carries over to Hermes: skewness of
  the positive counts, above 30 yellow ("often signals a problem"), above 50 red ("< 50 desired"); skew can be higher
  under strong selection. Source: the quality-control guidance of TRANSIT (DeJesus et al. 2015), written for Himar1.
- **Project-relative** — every other flagged metric compares a library with the project's libraries **of the same
  role** (parents with parents, treated with treated, from the sample sheet's `library_role`), so selection in treated
  libraries is not mistaken for a defect. Robust z = (value − median) / (1.4826 × MAD), in the metric's bad direction
  only; yellow beyond 2.5 and red beyond 3 (Leys et al. 2013, *J Exp Soc Psychol* 49:764: "moderately" and "very"
  conservative). At least 4 libraries of the role are needed, otherwise grey with that reason.
  The spread is at least 10% of the role median, so near-identical libraries are not flagged for trivial
  differences. **This 10% floor is a YTAB choice, not a published value** (`MIN_RELATIVE_SPREAD`).
- **Not flagged (descriptive):** mean_ct, max reads at a site (depends on depth; the top-site share is flagged
  instead), kurtosis (not very meaningful on its own), strand correlation (for Himar1 both orientations share a TA site;
  for Hermes they are separate insertions), MidLC estimate and depth / MidLC.

**Hermes.** Hermes has no fixed set of possible insertion sites, so the absolute density (> 30%), NZmean (30–200) and
total-read (> 1 M) thresholds of the Himar1 guidance are not applied; insertions per kb and the project-relative comparison take their place.

**Normalization hints** (from the same QC guidance), shown on yellow and red cells: skewed or jackpot-dominated
libraries suit BGC (betageom); TTR is robust to isolated outlier sites. In the QQ plot against a geometric fit,
strong curvature away from the diagonal also suggests BGC.

## In the app

The app uses neutral wording throughout; the sources of the rules are given here, not in the interface.

- **Overview:** build or choose a QC table; the coloured table (hover a cell for its reason, the metric's meaning and
  the hint; click a library for its count histogram, QQ against a geometric fit and ranked counts); downloads (Excel,
  TSV, flags).
- **Reads & Mapping, Insertions & Genes, Saturation & Bias:** one plot at a time from a menu: every QC-table metric of
  the theme across libraries (coloured by flag), plus the theme's curves and summaries (total and mapped reads; count
  distributions, genes by insertion count, genes hit by role, insertions per gene by role, insertion signal along
  chromosomes; saturation curves, jackpots and depth, centromere bias, TSS/TTS/tRNA metaplots, +2/+7 sequence bias).
  Every plot takes the shared display options (style, height, width, text size, label mode and angle, grid, point
  size and opacity, bar orientation, value labels, heatmap palette and labels), a role filter and colour-by where they
  apply, and PNG/PDF/plotted-data downloads. Saturation & Bias also holds the library-diagnostics run controls.
- **Library Agreement:** run library correlation **within this project** (e.g. Cg1 against Cg3) or **with other
  projects** (e.g. this project's libraries against Gale's); then the correlation heatmap, a gene insertion heatmap
  (most variable genes) or two libraries gene by gene, plus the pair table. `ytab library-correlation` now also writes
  `gene_values.tsv.gz` (per library and gene: sites, reads, density, reads per kb) for the gene-level plots.

The former Mapping QC, Summary QC and Library Diagnostics plot panels were removed (backup:
`output/ui_reorganization/backup/qc_legacy_plots/`); their plots are in the menus above. Old links open the matching
subtab.

## First runs (2026-10-21)

- **H2O2:** one treated library (yH299 pool 4) is red on six metrics that point the same way: NZmean, top site's
  share, top 1% share, skewness (116), Pickands index, reads in jackpots — a jackpot-dominated library (BGC hint).
- **Zn:** one control library (yH299 pool 4) has NZmedian 170 against 14 for the other controls.
- **Gale:** fewer than 4 libraries per role, so only the skewness rule of thumb applies: five libraries are red (skewness 113–168) and the
  two parent libraries from pools Cg2 and Cg3 yellow (36 and 43), consistent with BGC being suggested for Gale.
