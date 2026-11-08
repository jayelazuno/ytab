#!/usr/bin/env python3
"""Unique-molecule (PCR-duplicate-free) complexity of a paired-end Hermes library, from sonication break points.

Each original DNA molecule carrying a Hermes junction ends at a random sonication break point, recorded
by the mate of the junction read. Pairs with the same junction site and strand and the same fragment
length are PCR copies of one molecule. Distinct (site, strand, length) keys estimate unique molecules
(a slight underestimate when two molecules from one clone share a length by chance).

Prepare a paired BAM on the HPC (a subsample is enough: zinc libraries saturate near 1M reads), e.g.
    seqtk sample -s2026 SAMPLE_R1.fastq.gz 2000000 > r1.fq; seqtk sample -s2026 SAMPLE_R2.fastq.gz 2000000 > r2.fq
    bowtie2 --local -X 1500 -p 8 -x BT2_INDEX -1 r1.fq -2 r2.fq | samtools view -b -o SAMPLE.pairs.bam -
then
    python ytab_fragment_complexity.py --bam SAMPLE.pairs.bam --out SAMPLE.fragments [--junction-mate R1]
Outputs SAMPLE.fragments.summary.json and SAMPLE.fragments.sites.tsv.gz (site, strand, pairs, molecules).
"""
import argparse
from collections import Counter
import gzip
import json

import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True); parser.add_argument("--out", required=True)
    parser.add_argument("--junction-mate", choices=("R1", "R2"), default="R1"); parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--max-fragment", type=int, default=1500)
    args = parser.parse_args()
    pairs, molecules = Counter(), set(); seen = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            seen["records"] += 1
            if read.is_secondary or read.is_supplementary or not read.is_paired:continue
            if read.is_read1 != (args.junction_mate == "R1"):continue
            seen["junction_reads"] += 1
            if read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair or read.mapping_quality < args.mapq:continue
            length = abs(read.template_length)
            if not 0 < length <= args.max_fragment:continue
            strand = "-" if read.is_reverse else "+"
            site = (read.reference_name, read.reference_end if read.is_reverse else read.reference_start+1, strand)
            pairs[site] += 1; molecules.add(site+(length,))
    per_site = Counter(m[:3] for m in molecules)
    total_pairs, total_molecules = sum(pairs.values()), len(molecules)
    with gzip.open(args.out+".sites.tsv.gz", "wt") as handle:
        handle.write("chrom\tposition\tstrand\tpairs\tmolecules\n")
        for site, n in sorted(pairs.items()):handle.write(f"{site[0]}\t{site[1]}\t{site[2]}\t{n}\t{per_site[site]}\n")
    multi = [per_site[s] for s in pairs]
    summary = dict(bam=args.bam, junction_mate=args.junction_mate, mapq=args.mapq, records=seen["records"], junction_reads=seen["junction_reads"],
                   usable_pairs=total_pairs, sites=len(pairs), unique_molecules=total_molecules,
                   duplication_rate=1-total_molecules/total_pairs if total_pairs else None,
                   pairs_per_molecule=total_pairs/total_molecules if total_molecules else None,
                   molecules_per_site_mean=total_molecules/len(pairs) if pairs else None,
                   sites_with_one_molecule=sum(m == 1 for m in multi)/len(multi) if multi else None)
    with open(args.out+".summary.json", "w") as handle:json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":main()
