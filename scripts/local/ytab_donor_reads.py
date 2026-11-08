#!/usr/bin/env python3
"""Count the reads that come from the transposon donor, without remapping the library to the genome.

The reads are put against the donor's own index (a few kilobases, so one quick pass); only the reads that hit it, or
that carry a recurring junction flank, are then checked against the genome. A library already mapped is not mapped
again: `ytab mapfastq` writes the same table while it runs.

  ytab donor-reads --project-config output/projects/P/config/project.yaml           # every included library
  ytab donor-reads --project-config ... --samples A,B --threads 8
  ytab donor-reads --fastq reads.fastq.gz --sample NAME --donor-fasta D.fasta --output DIR

Writes <sample>.donor_reads.csv into each library's mapfastq folder (or --output), which `ytab qc` reads for the
"Reads from the transposon donor (%)" metric.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.mapping import donor_reads  # noqa: E402
from ytab.pipeline.mapfastq_runner import get_included_samples, load_project_for_mapping  # noqa: E402
from ytab.mapping.transposons import resolve_design  # noqa: E402


def resolve(value, root):
    if value in (None, ""):return None
    p = Path(value)
    return p if p.is_absolute() else root/p


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", type=Path); p.add_argument("--samples", help="Comma-separated sample names")
    p.add_argument("--fastq", type=Path); p.add_argument("--sample"); p.add_argument("--output", type=Path)
    p.add_argument("--donor-fasta", type=Path); p.add_argument("--no-genome", action="store_true",
                   help="Skip the genome check of donor-hit reads (no split into donor-only and shared)")
    p.add_argument("--genome-index", help="Bowtie2 index prefix, to split donor reads into donor-only and shared (with --fastq)")
    p.add_argument("--transposon-fasta", type=Path, help="Transposon reference, masked inside the donor so only the flanks align (default: the project's)")
    p.add_argument("--mapq", type=int, default=20); p.add_argument("--threads", type=int, default=4)
    a = p.parse_args()

    if a.fastq:
        if not a.sample or not a.output:p.error("--fastq needs --sample and --output")
        a.output.mkdir(parents=True, exist_ok=True)
        row = donor_reads.from_fastq(str(a.fastq), a.sample, str(a.output), str(a.donor_fasta) if a.donor_fasta else None,
                                     a.genome_index, a.mapq, a.threads, str(a.transposon_fasta) if a.transposon_fasta else None)
        value = row.get("percent_donor") if row.get("percent_donor") is not None else row.get("percent_recurrent_flank")
        print(f"{a.sample}: {value}% donor reads, donor-only {row.get('percent_donor_only')}%, "
              f"placed on the genome {row.get('percent_donor_shared')}% ({row['source']})")
        return

    if not a.project_config:p.error("--project-config or --fastq is required")
    config = load_project_for_mapping(a.project_config); root = Path(config["_repo_root"])
    wanted = {s.strip() for s in a.samples.split(",")} if a.samples else None
    donor = resolve((config.get("mapping") or {}).get("donor_fasta"), root)
    index = None if a.no_genome else (config.get("reference") or {}).get("bowtie2_index_prefix")
    index = str(resolve(index, root)) if index else None
    project = resolve(config.get("output_project_dir") or f"output/projects/{config['project_id']}", root)
    for sample in get_included_samples(config):
        name = str(sample["sample"])
        if wanted and name not in wanted:continue
        fastq = resolve(sample.get("fastq_1") or sample.get("r1_fastq"), root)
        if fastq is None or not fastq.is_file():print(f"{name}: no FASTQ ({fastq})"); continue
        out = project/"mapfastq"/name; out.mkdir(parents=True, exist_ok=True)
        design = resolve_design({**(config.get("mapping") or {}), **{k: v for k, v in sample.items() if v not in (None, "")}}, root)
        tn = str(a.transposon_fasta) if a.transposon_fasta else next((str(x) for x in sorted((ROOT/"resources/transposons"/design["transposon"]).glob("*.fasta"))), None)
        row = donor_reads.from_fastq(str(fastq), name, str(out), str(donor) if donor else None, index, a.mapq, a.threads, tn)
        value = row.get("percent_donor") if row.get("percent_donor") is not None else row.get("percent_recurrent_flank")
        print(f"{name}: {value}% donor reads, {row.get('percent_donor_shared')}% placed on the genome ({row['source']})", flush=True)


if __name__ == "__main__":main()
