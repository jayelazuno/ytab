# Gale 2020 local downstream analysis

## Validation hold: automatic normalization units

The raw summaries and library QC completed successfully for all seven libraries.
The subsequent workflow executed through the classifier/report/export, but its
normalization-dependent results are **provisional, not validated**. The existing
automatic selector used raw MidLC read-count estimates as candidate multipliers,
selected T409600, and produced no thinning for any of the three parent libraries.
For Cg1, MidLC=409600 and multiplier=409600 produced reads_target=167772160000,
although the input contained only 23743484 reads. The normalization comparison
table confirms unchanged read totals in all three parent samples. The selector's
units and fallback/retention behavior require correction/review before a validated
normalized-parent/classifier analysis can be completed. No final classifier
target has been accepted. Do not rerun/promote this automatic workflow until that
issue has been addressed.

The affected site's and feature-level automatic recommendations have been
invalidated (their former target is preserved for audit). Raw summaries, QC,
original hit files and the diagnostic T409600 outputs remain intact. The first
report/archive predates this validation hold and must not be treated as validated
analysis. See `VALIDATION_STATUS.json` for the machine-readable status.

Source: Gale et al. (2020), https://doi.org/10.1534/g3.120.401595,
SRA BioProject PRJNA625944 and the accession/sample mapping supplied for this project.

| Run | Library | Independent induction pool | Use |
| --- | --- | --- | --- |
| SRR11561508 | Cg1_parent | Cg1 | Parent essentiality branch |
| SRR11561507 | Cg2_parent | Cg2 | Parent essentiality branch |
| SRR11561506 | Cg3_parent | Cg3 | Parent essentiality branch |
| SRR11561505 | Cg1_control1 | Cg1 | Regrown no-drug screen control |
| SRR11561504 | Cg1_control2 | Cg1 | Regrown no-drug screen control |
| SRR11561503 | Cg1_FLU1 | Cg1 | FLU screen culture 1 |
| SRR11561502 | Cg1_FLU2 | Cg1 | FLU screen culture 2 |

All libraries use strain BG14, a ura3-deficient BG2 derivative. Mapping was already
completed on Argon against the project's BG2 reference, not the CBS138 reference
used for the published analysis. This is a YTAB reanalysis, not an exact
reproduction of the paper's mapping/counting and statistical settings.

The original HPC mapping and hit-file manifests are retained as provenance.
The local run starts at raw SummaryTable generation from seven completed hit files;
no FASTQs or BAMs are required and neither stage is rerun.

`downstream_workflow.yaml` explicitly selects only the three original independent
parent pools for MidLC normalization, normalized summaries and parent combination.
The two regrown control cultures are not additional independently induced pools.
Do not replace this workflow with a generic parents-only selection without
checking it: the existing backend also recognizes screen controls as parents.

The existing automatic MidLC/site-retention selection and feature-retention check
are used without changing scientific methods. The existing essentiality classifier
uses seed 0; no random-forest improvements are introduced here, and no final target
is marked as manually reviewed. Read normalization warnings before interpreting
the selected target. All QC/summary denominators and metrics come from the inputs.

The FLU experiment uses exactly four cultures derived from Cg1: control1
(SRR11561505), control2 (SRR11561504), FLU1 (SRR11561503), and FLU2 (SRR11561502).
The user reconfirmed this separation from the three original induction libraries.
`comparison_design.csv` defines numbered FLU1/control1 and FLU2/control2 contrasts
for the existing contrast-based runner. These labels do not assert matched
biological induction pairs: both contrasts retain the same source pool, Cg1.
The two no-drug controls jointly define the control/control noise model.
Original Cg1/Cg2/Cg3 parent sequencing observations are not fitness controls and
do not enter the fitness analysis. Replicate agreement is culture-level agreement,
not evidence of replication across independent insertion-pool inductions.

`fitness_workflow.yaml` runs this branch independently from raw per-culture
SummaryTables using CPM in the established fitness R script. It does not consume
MidLC-normalized tables and does not annotate with the provisional classifier.
Thus the classifier normalization hold does not block the FLU fitness branch.

The four-culture fitness run completed successfully under analysis ID
`Cg1_FLU_vs_no_drug_raw_cpm`. Its recorded raw-summary inputs contain exactly the
four screen cultures above, and its control/control noise model records
SRR11561505 versus SRR11561504. Classifier annotation is disabled. Results are in
`treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm/`; plots and provenance use the same
analysis directory. The separate parent-classifier validation hold still applies.

Reproduce/cache-resume from the repository root:

```bash
python scripts/local/ytab.py run \
  --workflow output/projects/Gale_2020/config/downstream_workflow.yaml --threads 4
```

Scientific results use the canonical `output/projects/Gale_2020/` and
`output/exports/Gale_2020/` routes and can be opened through the desktop app.
