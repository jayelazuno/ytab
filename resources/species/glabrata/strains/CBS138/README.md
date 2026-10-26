# *C. glabrata* CBS138 (species `glabrata_CBS138`)

The reference strain CBS138 (ATCC 2001) and its derivatives (e.g. CBS138u, 2001u of Nickels et al. 2024) are mapped
here instead of to the BG2 reference in `resources/species/glabrata/` (about 1% sequence divergence; on 2,000 reads of
SRR29743993, 87% of reads mapped at MAPQ >= 20 and 78% without mismatch here, against 81% and 41% on BG2).

Sources (CGD, downloaded 2026-10-26; SHA-256, first 16 characters):
- `reference_genome/C_glabrata_CBS138_version_s05-m03-r48_chromosomes.fasta`: http://www.candidagenome.org/download/sequence/C_glabrata_CBS138/current/ (gunzipped) — `66b388dee96dc2d4`
- `reference_genome/C_glabrata_CBS138_version_s05-m03-r48_features.gff`: http://www.candidagenome.org/download/gff/C_glabrata_CBS138/ — `2ea9405aa4fd0497`
- `orthologs/C_glabrata_CBS138_S_cerevisiae_orthologs.txt`: http://www.candidagenome.org/download/homology/orthologs/C_glabrata_CBS138_S_cerevisiae_by_CGOB/ — `60fbd8041618b298`

Built from them:
- `reference_genome/centromeres.CGD.bed`, `tRNAs.bed`, `chromosome_names.tsv` and `annotation/gene_annotation.tsv`:
  `python scripts/local/ytab_build_species_resources.py --species glabrata_CBS138 --annotation-only --tracks`
- Bowtie2 index and FASTA index: `python scripts/local/ytab.py reference --species glabrata_CBS138`

No essentiality labels are built for this strain yet.
