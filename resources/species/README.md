
# YTAB species resources

This directory contains the organism-specific reference resources used by YTAB.

The goal of this layout is to keep all species assets in one stable location outside of `src/`, so pipeline code, QC scripts, app code, and future workflow managers all resolve the same reference files from a single place.

## Purpose

Each species subdirectory stores the files needed for one or more stages of the YTAB workflow, including:

- genome sequence
- genome annotation
- feature tables
- ortholog maps
- aligner indexes
- any additional species-specific support files

These resources are treated as **inputs** to the pipeline and app, not generated outputs.

## Directory layout

Every species uses the same subdirectories; code resolves files only through them.

```text
resources/species/
├── README.md
├── orthologs/                 cross-species 1:1 ortholog groups (glabrata, cerevisiae, pombe, albicans)
└── <species>/                 albicans, cerevisiae, glabrata, pombe, or any new species
    ├── reference_genome/      FASTA, GFF/GTF, feature tables, Bowtie2 index, chromosome_names.tsv,
    │                          centromere and tRNA BEDs (resolved by ytab.pipeline.reference_registry)
    ├── annotation/            gene_annotation.tsv (names, S. cerevisiae orthologs, descriptions), GO files,
    │                          paralogs, repeated/deleted regions, genetic-interaction and fitness tables
    ├── orthologs/             the species' ortholog tables (S. cerevisiae first; S. pombe where available)
    ├── essentiality_labels/   labels.tsv, sources.json, summary.json, README.md; raw/ holds the source files
    │                          (phenotype downloads, viability annotations, published essential-gene lists)
    └── tnseq_data/            published transposon insertion data, one folder per study
                               (e.g. cerevisiae/tnseq_data/michel_2017_satay, pombe/tnseq_data/guo_2013_hermes)
```

A strain whose genome differs enough to need its own reference (for example *C. glabrata* CBS138, about 1% divergent
from the BG2 reference) has the same layout under `<species>/strains/<strain>/` and is selected as species
`<species>_<strain>` (e.g. `glabrata_CBS138`); its reads are mapped to, and its genes annotated from, that folder.

Only the folders a species needs exist. A species becomes selectable when `reference_genome/` contains at least one
recognized reference file. Files are inputs to the pipelines and app, never outputs.

