# Hermes read selection and automatic cleanup

The project mapping runner uses **R1 as the designated Hermes/junction read** by
default. Users can declare R2 instead; it does not guess mate identity from mapping
rate. Both sequencing mates may map well, but only one starts at the insertion.
The default does not require the other mate to map, or even to be present.

The existing CLI, master pipeline and Shiny mapping controls all use this runner.
No changes to normalization, fitness, Poisson planning or reference annotation are
introduced here.

## Running

```sh
# Default: clean and align R1 only.
python scripts/local/ytab.py mapfastq --project-config PROJECT.yaml --threads 8

# Select R2; do not require the other mate to align.
python scripts/local/ytab.py mapfastq --project-config PROJECT.yaml \
  --junction-mate R2 --no-require-mate-alignment --threads 8

# Explicitly require concordant paired-end placement.
python scripts/local/ytab.py mapfastq --project-config PROJECT.yaml \
  --junction-mate R1 --require-mate-alignment --threads 8
```

The Mapping tab provides equivalent junction-read and mate-agreement selectors.
Their default respects project/sample settings. CLI overrides apply to the selected
mapping run without rewriting the user's sample sheet.

Reusable project defaults go under `mapping` in project YAML:

```yaml
mapping:
  junction_mate: R1
  require_mate_alignment: false
  hermes_cleanup: true
  junction_read_end: start
  hermes_overlap: 16
  adapter_error_rate: 0.10
  min_read_length: 25
```

Sample-sheet columns with the same names override project defaults; blank cells
inherit defaults. `genome_mapping_mate` remains a supported alias for the mate.
An explicitly selected legacy `mapping_strategy: paired_genome` implies requiring
concordance unless `require_mate_alignment` explicitly overrides it. Strategy
`auto` means the normal single-mate default, **not evidence-based mate inference**;
`genome_mapping_mate: auto` is rejected with instructions to choose R1 or R2.

## Automatic sequence cleanup

The bundled `resources/transposons/hermes/primers.json` contains the library PCR
primers, sequencing primers, splinkerette sequences, inline indices and common
Illumina adaptor sequences. Source: Gale et al. 2020, Supplemental Table S3,
doi:10.25387/g3.12445373; paper doi:10.1534/g3.120.401595. The accessioned natural
Hermes sequence is `L34807.1.fasta`; it is not the complete engineered launchpad.

The streaming cleaner uses Cutadapt's alignment library to:

- remove detected outward-facing left/right Hermes terminal sequence from the
  junction read's 5' region, including known retained primer-plus-terminal tails;
- remove complete retained known primer/adaptor sequence at the 5' start;
- remove detected 3' adaptor/Hermes read-through, including reverse complements;
- discard reads shorter than the configured minimum after cleanup.

Default minimum overlap is 16 nt, maximum error rate 0.10, with no indels in the
adapter match. Hermes matches must reach the transposon boundary, with a start
within the first 80 read bases. Six-base barcodes alone are **not** clipping
anchors. Low-complexity 3' matches containing fewer than three distinct bases are
ignored. The percentage/overlap thresholds and reference/catalogue are recorded.
Cleanup currently streams in one process; the thread setting applies to alignment
and BAM processing.

Crucially, **a Hermes match is not mandatory**. Gale's Hermes sequencing primer's
56-nt terminal suffix ends at reference position 2749, the transposon boundary;
the first sequenced base is therefore expected to be genomic. HR2 PCR priming
instead leaves 33 Hermes bases, but there is no unconditional 33-bp cut. Only
detected sequence is removed. Low-quality genomic bases at the junction-facing
start are not silently trimmed, since that would shift the insertion coordinate.

For a reverse-complemented single-read export, `junction_read_end: end` reverses
both sequence and quality before cleanup. This is not required for ordinary
reverse-strand genomic alignments: Bowtie2 already searches both genomic strands.
Nonstandard end-facing exports are not accepted in concordant-pair mode.

Custom `primer_catalogue` and `hermes_reference` paths may be provided in the same
mapping settings. The catalogue must follow the bundled schema. For a non-Hermes
dataset or an externally precleaned input, use `hermes_cleanup: false` or the CLI
flag `--no-hermes-cleanup`. Legacy direct `MapFastq.py` calls remain unchanged unless
`--hermes-preprocess` is passed.

## Alignment and reporting

Default alignment passes only the designated cleaned read to Bowtie2. Thus an
unmapped opposite read cannot remove a valid junction alignment. End-to-end
mapping preserves the cleaned read start; normal MAPQ filtering still applies.

With explicit agreement enabled, Bowtie2 uses `--no-mixed --no-discordant` and the
final BAM requires the proper-pair flag. **Only the designated R1 or R2 record is
emitted**, never both endpoints. This preserves downstream insertion counting:
forward reference start +1, or reverse reference end, in 1-based coordinates.
Paired FASTQs must have matching names/counts and both cleaned reads must pass
minimum length. No UMI/PCR deduplication or target-site duplication correction is
introduced by this cleanup.

Each sample writes `<sample>.hermes_preprocessing.json` beside its sorted BAM,
including raw input paths, declared mate, orientation setting, retained counts,
detected signatures, clipped-base counts, Cutadapt version and resource hashes.
The report is included in the mapping manifest; job logs report cleanup progress.
Raw FASTQs are never modified or deleted by this path.

Mapping cache keys include selected input file size/mtime, commands, index
size/mtime, mapping/preprocessing code, resource contents, Cutadapt version and
aligner-tool paths/size/mtime. Older unversioned mapping caches are not silently
reused. A preview is not a successful cache entry.
The master runner delegates mapping cache checks to the sample runner.

When intentionally replacing existing alignments, regenerate downstream hit files,
summaries and analyses as well (for example with a deliberately scoped master
`--force` run). This change does not automatically overwrite existing downstream
scientific results. This implementation has been tested on synthetic junctions;
the downloaded Gale FASTQs on the HPC still require their first production run.
