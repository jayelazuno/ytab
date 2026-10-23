# Empirical Hermes catalogue and TTR evaluation

This standalone, opt-in analysis builds a common observed insertion-coordinate
universe and evaluates trimmed total reads (TTR) normalization, adapted from TRANSIT's TTRNorm
(DeJesus et al. 2015; see [references](references.md)).
It does **not** replace classifier MidLC, fitness CPM, or the existing
`normalization` CLI module, and it does not add a Shiny tab.

## Run

```sh
python scripts/local/ytab.py ttr \
  --manifest resources/normalization/gale_our_parent_catalogue_v1.yaml \
  --output output/normalization/gale_our_parents_ttr_v1
```

The direct runner is `scripts/local/ytab_ttr.py`. Run `ytab.py ttr --help` for
options. `--trim` defaults to 0.05 **per tail**; `--target` defaults to 100.
Choose an empty/new output directory for each version; existing results are not
overwritten. The same module can be an explicit step in a CLI master workflow.

## Inputs and identity

The YAML manifest explicitly lists source project configurations, included parent
or untreated/control observations, dataset groups, and independent `library_id`
values. Source locations and dataset identities live in metadata, not in the
algorithm. An observation ID is `project_id::sample`; sample accessions remain
stable. Repeated sequencing or outgrowth observations of one induced pool must
share the same independent-library ID. Treatment samples are rejected.

All observations must reference an identical FASTA checksum. The analysis uses
the existing unnormalized YTAB `_hits.txt` files and their reported chromosome
and 1-based `Hit position`. Raw nonnegative `Hit count` values at identical
coordinates are added, retaining W/C strand provenance. There is no additional
coordinate shift, TSD correction, motif filter, or remapping.

YTAB browser WIG files encode strand with the sign of the count. They are **not**
summed as signed values here: opposite strands must not cancel. The canonical
nonnegative hit files avoid that ambiguity.

The universe is the **observed union**, not every genomic base or every possible
Hermes target. Within that union, an absent sample/site value is zero. This means
not observed at the available depth, not biological impossibility. The output
retains sample, pool, independent-library and dataset provenance for each site.
No gene-length, motif, bottleneck or loss-of-function correction is fitted.

## TTR calculation

For sample `j`, over one common catalogue of `U` sites:

```
theta_j = number of sites with positive raw counts / U
mu_j    = mean of the positive counts after trimming floor(t*n_positive)
          observations from each tail
factor_j = target / (theta_j * mu_j)
normalized_count_ij = raw_count_ij * factor_j
```

Trimming estimates the factor; it does not discard sites or alter their raw
counts. Zero counts remain zero. Counts stay floating-point after scaling, with
no stochastic read thinning or rounding to integers. Normalized total support
need not be equal across samples because the estimator trims outliers.

For a common universe, adding sites that are zero in every sample multiplies
all factors by the same amount; it does not change factor ratios. A larger
catalogue alone therefore does not demonstrate better relative normalization.
Differences in true pool complexity also affect density: TTR scaling should not
be interpreted as recovering missing mutants or as proven to preserve every
biological abundance difference across independently induced pools.

The formula and defaults were checked directly against
[TRANSIT's TTRNorm implementation](https://github.com/ioerger/transit/blob/50e23d6bec2605a24835fbd67f33691e330f3361/src/pytransit/norm_tools.py#L148).
The local optional upstream checkout is used for numerical parity tests, not
as a runtime dependency or copied normalization class.

## Outputs

- `observations.tsv`: explicit sources, library identities and input/reference hashes.
- `sites.tsv.gz` and `site_provenance.tsv.gz`: coordinate universe and per-observation support.
- `raw_counts.npz` and `raw_counts.tsv.gz`: identical site-by-observation matrices;
  sparse zeros in NPZ are explicitly written in the tabular export.
- `ttr_normalized_counts.npz` and `.tsv.gz`: scaled matrices, retaining zeros.
- `ttr_parameters.tsv`: occupied sites, density, trimmed mean, factor and target check.
- `dataset_overlap.tsv`, `library_overlap.tsv`, `observation_overlap.tsv`:
  intersections, unique contributions and Jaccard overlap at each provenance level.
- `dataset_contributions.tsv`: sites added beyond the other dataset groups.
- `library_discovery.tsv`: exact expected union over subsets of the observed
  independent libraries, not an extrapolated saturation estimate.
- `depth_rarefaction.tsv`: expected observed sites at shared count depths under
  hypergeometric subsampling without replacement; no reads are actually discarded.
- `ttr_diagnostics.png`, `site_discovery.png`: compact diagnostic figures.
- `run_manifest.json`, `source_manifest.yaml`: parameters, identities and limitations.

Depth in this QC is **usable called-insertion support**, not total FASTQ reads.
PCR/sequencing support is not a direct measurement of viable clone abundance.
Repeated observations improve discovery but are not extra independent inductions.
The empirical catalogue, overlap, and discovery curves do not by themselves prove
that all recoverable nonessential-genome insertion sites have been saturated.

## Tests

```sh
PYTHONPATH=src python -m pytest -q tests/test_ttr.py tests/test_cli.py
```

These cover numerical parity with TRANSIT when its optional local checkout is
present, trimmed-mean rules, common-universe scaling, zero preservation, matrix
exports, strand merging, library provenance, treatment rejection, and exact small
rarefaction examples. No production analysis consumes these TTR outputs yet.
