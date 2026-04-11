
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

```text
resources/species/
├── README.md
├── glabrata/
│   ├── C_glabrata_BG2_S_cerevisiae_orthologs.txt
│   └── reference_genome/
│       ├── GCA_014217725.1_ASM1421772v1_feature_table.txt
│       ├── GCA_014217725.1_ASM1421772v1_genomic.fna
│       ├── GCA_014217725.1_ASM1421772v1_genomic.gff
│       ├── GCA_014217725.1_ASM1421772v1_genomic.gtf
│       ├── GCA_014217725.1_ASM1421772v1_protein.faa
│       ├── GCA_014217725.1_ASM1421772v1_genomic.1.bt2
│       ├── GCA_014217725.1_ASM1421772v1_genomic.2.bt2
│       ├── GCA_014217725.1_ASM1421772v1_genomic.3.bt2
│       ├── GCA_014217725.1_ASM1421772v1_genomic.4.bt2
│       ├── GCA_014217725.1_ASM1421772v1_genomic.rev.1.bt2
│       └── GCA_014217725.1_ASM1421772v1_genomic.rev.2.bt2
├── cerevisiae/
├── pombe/
├── albicans/
└── Kornmann/
