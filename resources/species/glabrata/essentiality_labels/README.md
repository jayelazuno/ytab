# Curated C. glabrata essentiality labels

Built by `scripts/local/ytab_curate_essentiality_labels.py` from CGD (direct C. glabrata evidence) and
S. cerevisiae (SGD) via 1:1 orthologs only. Rules are in the script docstring; sources, retrieval dates
and SHA-256 hashes are in `sources.json`; counts in `summary.json`.

| tier | evidence | essential | nonessential |
|---|---|---|---|
| gold | C. glabrata deletion / repressible allele (CGD, excluding Gale 2020) | 20 | 256 |
| silver | S. cerevisiae null-mutant viability, 1:1 orthogroup | 924 | 2,367 |

- **Unlabelled (1,653 genes):** not 1:1 (paralog groups, mainly WGD ohnolog pairs), no orthogroup, S.
  cerevisiae conflict or no data, or C. glabrata evidence that is condition-dependent (RIB1–5, RIB7
  inviable only without riboflavin).
- **Gold vs S. cerevisiae ortholog:** 141 of 151 agree. The 10 discordant genes (e.g. ERG11, PGS1,
  CBF1, OPI1, GSH1) are real species differences and are labelled by the C. glabrata evidence.
- **`gale2020_comparator`:** Gale et al. 2020 Hermes calls as recorded in CGD (PMID 32819971). These
  are never labels: they come from the same data that YTAB models are run on. Agreement with the
  curated labels is 91% (gold 92%, silver 91%).

Use: train on silver, evaluate on gold plus held-out silver; report Gale agreement separately.
