# Curated K. phaffii essentiality labels

Built by `scripts/local/ytab_curate_phaffii_labels.py` from direct knockout evidence and from *S. cerevisiae* (SGD)
viability through 1:1 orthologs only. Rules are in the script docstring; sources and SHA-256 hashes are in
`sources.json`; counts in `summary.json`.

| tier | evidence | essential | nonessential |
|---|---|---|---|
| gold | a viable knockout was made (Zhu et al. Dataset S3) | — | 171 |
| silver | *S. cerevisiae* null-mutant viability (SGD), 1:1 orthogroup | 908 | 2,351 |

5,040 annotated genes; 3,430 labelled, 1,610 unlabelled.

- **Unlabelled (1,610 genes):** no 1:1 *S. cerevisiae* ortholog in Zhu's Dataset S6 (which pairs 4,132 genes, 3,528
  of them 1:1), or a 1:1 ortholog whose SGD viability is not curated.
- **Asymmetry, stated plainly:** as for *C. neoformans*, there is no gold *essential* source. Every essential label
  is silver, inferred from the *S. cerevisiae* ortholog; the gold labels are all nonessential. An evaluation on these
  labels is more trustworthy about what is **not** essential than about what is.
- **Conflicts (2 genes):** a viable knockout exists although the *S. cerevisiae* ortholog is essential; the gold
  label wins and the silver one is kept in `conflict`.
- **`zhu2018_comparator` and `zhu2018_nep`:** Zhu et al.'s own category and non-essential probability. **These are
  never labels**, for two reasons. They come from a logistic regression on the gene insertion index of their own
  TcBuster/Sleeping beauty transposon libraries — the same kind of data YTAB models are run on — and the model's
  essential training set was itself built from *S. cerevisiae* and *S. pombe* orthologs, which is how the silver tier
  here is built, so using them would double-count one line of evidence. Agreement with the curated labels is 83.7%
  (566 of 591 essential, 1,388 of 1,743 nonessential).
- **`zhu2018_training_essential`:** the 161 genes Zhu et al. used as essential training examples, recorded for
  transparency. Not a label, for the same reason.

## Source

Zhu J, Gong R, Zhu Q et al. (2018) Genome-wide determination of gene essentiality by transposon insertion sequencing
in yeast *Pichia pastoris*. Sci Rep 8:10223. doi:10.1038/s41598-018-28217-z. Datasets S3, S4 and S6 are in
`raw/zhu_2018_scirep/`.
