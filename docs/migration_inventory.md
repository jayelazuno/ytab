

# YTAB migration history

## 2026-03-12 — initial core migration from HermesTnseq

### Goal
Start a clean migration of the core Hermes transposon-seq pipeline into the new YTAB repo structure, while preserving behavior in the first pass.

### New repo identity
- Project name: YTAB
- Expanded name: Yeast Transposon Analysis Browser
- Repo: `ytab`

### High-level decisions
- Use a new repo rather than reorganizing the legacy repo in place.
- Migrate only stable core pipeline code first.
- Leave species/reference asset reorganization for a later phase.
- Leave Docker, Apptainer, and Nextflow setup until the core code structure is stable.
- Run the app locally first; HPC will be used for the pipeline/runtime side.

### Legacy repo assessment
Legacy repo contains:
- numbered stage/output directories (`02.results`, `03.mapping`, `04.hits`, `05.summary`, `05.summary_normalized`, `06.classifier`, `07.QC`)
- HPC job wrapper scripts
- main source code in `transposon-pipeline/`

Conclusion:
- `transposon-pipeline/` is the real migration source
- numbered directories are mainly workflow/output references, not the new code structure

### Core files migrated
Copied from `HermesTnseq/transposon-pipeline/` into `ytab/src/ytab/`:

#### Pipeline-facing modules
- `MapFastq.py` -> `src/ytab/mapping/MapFastq.py`
- `CreateHitFile.py` -> `src/ytab/insertions/CreateHitFile.py`
- `SummaryTable.py` -> `src/ytab/summary/SummaryTable.py`
- `Classifier.py` -> `src/ytab/essentiality/Classifier.py`
- `LibraryDiagnostics.py` -> `src/ytab/qc/LibraryDiagnostics.py`

#### Core support modules
- `GenomicFeatures.py` -> `src/ytab/core/GenomicFeatures.py`
- `Organisms.py` -> `src/ytab/core/Organisms.py`
- `Shared.py` -> `src/ytab/core/Shared.py`
- `RangeSet.py` -> `src/ytab/core/RangeSet.py`
- `SortedCollection.py` -> `src/ytab/core/SortedCollection.py`

### Package structure created
Created package directories and `__init__.py` files under:
- `src/ytab/`
- `src/ytab/core/`
- `src/ytab/mapping/`
- `src/ytab/insertions/`
- `src/ytab/summary/`
- `src/ytab/essentiality/`
- `src/ytab/qc/`

### Import refactor completed
Updated old flat-directory imports to package-style imports, e.g.:
- `import Shared` -> `from ytab.core import Shared`
- `import GenomicFeatures` -> `from ytab.core import GenomicFeatures`
- `from RangeSet import RangeSet` -> `from ytab.core.RangeSet import RangeSet`
- `import SummaryTable` -> `from ytab.summary import SummaryTable`

### Current status
Internal package migration is working.
Import failures have moved from internal module resolution to external dependency/environment issues.

### Import test results
Main remaining issues on local Mac environment:
- `ModuleNotFoundError: No module named 'gffutils'`
- `ModuleNotFoundError: No module named 'pysam'`
- `matplotlib` / `PIL` import failure due to missing `libtiff` in local base conda environment

Interpretation:
- YTAB package structure is valid
- Remaining issues are environment/dependency related, not code-layout related

### Environment notes
- Local Mac conda environment has solving/runtime issues
- Existing `tnseq` conda environment already exists on HPC and is the preferred place to continue
- Next step is to clone `ytab` on HPC and continue testing there

### Deferred work
Not yet done:
- species/reference asset reorganization
- manifest-based organism/resource handling
- Nextflow scaffold
- Docker image
- Apptainer/HPC runtime profile
- CLI wrappers
- lowercase renaming / code style cleanup

### Important notes
- There is a space in the current local parent directory path: `Repositories /`
- This may be annoying later for shell commands and mounts, but does not block migration
- Keep original capitalized filenames for now; rename later only after the migrated code is stable

### Immediate next step
On HPC:
1. clone `ytab`
2. activate existing `tnseq` environment
3. run import tests from repo root with `PYTHONPATH=src`
4. identify any remaining missing packages or path assumptions
5. continue cleanup there

## 2026-11-14 — packages by question asked, not by where the code grew

### Goal
`src/ytab/essentiality/` had become the home of most analysis code, including every treated-vs-control method and
the shared helpers they all import. A reader looking for the comparative methods found them under "essentiality",
and `qc`, `normalization` and `screen_coverage` imported genome annotation from a Gumbel module. The packages now
say what each module answers.

### Shared helpers extracted to `src/ytab/core/`
Both `control_conditioned_gumbel.py` and `library_correlation.py` had become de facto utility modules. What the
rest of the repo imported from them now lives on its own:
- `annotations`, `gene_sites`, `SITE_RULES`, `gene_key` -> `src/ytab/core/annotation.py`
- `load_design`, `load_normalization`, `NORMALIZATIONS`, `NONLINEAR_NORMALIZATIONS`, `DEFAULT_NORMALIZATION`
  -> `src/ytab/core/design.py`
- `benjamini_hochberg`, `_json_ready` (now `json_ready`) -> `src/ytab/core/stats.py`
- `read_hits` -> `src/ytab/core/hits.py`

### Comparative methods moved to `src/ytab/fitness/`
Each needs at least two conditions, so none of them is an essentiality classifier:
- `src/ytab/essentiality/control_conditioned_gumbel.py` -> `src/ytab/fitness/control_conditioned_gumbel.py`
- `src/ytab/essentiality/pooled_gumbel.py` -> `src/ytab/fitness/pooled_gumbel.py`
- `src/ytab/essentiality/resampling.py` -> `src/ytab/fitness/resampling.py`
- `src/ytab/essentiality/pooled_gene_ql.py` -> `src/ytab/fitness/pooled_gene_ql.py`
- `src/ytab/essentiality/multicondition.py` -> `src/ytab/fitness/multicondition.py`
- `src/ytab/essentiality/multicondition_plots.py` -> `src/ytab/fitness/multicondition_plots.py`
- `src/ytab/essentiality/genetic_interactions.py` -> `src/ytab/fitness/genetic_interactions.py`

### Library agreement moved to `src/ytab/qc/`
It measures how alike libraries are, not what a treatment did, and it powers Quality Control > Library Agreement:
- `src/ytab/essentiality/library_correlation.py` -> `src/ytab/qc/library_correlation.py`

`src/ytab/essentiality/` now holds only single-condition essentiality: `hermes_hmm.py`, `Classifier.py`,
`model_menu.py` and `contract.py`. Both package `__init__.py` files state the boundary.

### Pre-registration checksums
`docs/manuscript/preregistration_*.sha256` record source paths that no longer exist. They are left exactly as
written: they are the provenance of a registered analysis, not a live index. Their recorded hashes had already
diverged from the code before this migration (the files were developed further after registration), so verifying
one means checking out the registered commit, where both the old paths and the recorded contents are intact.
