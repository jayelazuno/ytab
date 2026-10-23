# Curated C. albicans essentiality labels

Built by `scripts/local/ytab_build_species_resources.py`; rules are in its docstring, sources and SHA-256 hashes in `sources.json`, counts in `summary.json`.

| tier | evidence | essential | nonessential |
|---|---|---|---|
| gold | C. albicans null (homozygous diploid or haploid) or repressible allele (CGD), excluding transposon surveys | 271 | 1,695 |
| silver | S. cerevisiae null-mutant viability, 1:1 ortholog | 729 | 1,815 |

- **Unlabelled (1,691 genes):** reasons in `summary.json`.
- **`transposon_comparator`:** transposon mutagenesis and machine-learning calls in a haploid strain as recorded in CGD (PMID 30377286). Never labels: they come from the same kind of data YTAB models are run on.
