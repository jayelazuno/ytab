"""Reads that come from the transposon donor rather than from an insertion in the genome.

In induced libraries the donor (a plasmid, or a chromosomal donor locus) is not always lost: reads that start at the
transposon end and run into the donor's flanking sequence are sequenced along with the genomic junctions. Mapping to
the genome drops those whose flank is donor-only sequence (vector, foreign markers, terminators) and places those
whose flank the donor shares with the host genome (a promoter, marker or CEN/ARS fragment) at that locus, where they
look like insertions.

Counting them does not need the genome mapping redone. The work is one pass of the reads against the donor, whose
index is a few kilobases; only the small subsets that follow are ever put against the genome:
  donor sequence known   every read is aligned to the donor; the reads that hit it (and only those) are aligned to
                         the genome, which splits them into donor-only and shared (the donor's genomic fragments).
  no donor sequence      a donor junction is one fixed sequence, so donor reads share their first PREFIX bases: the
                         prefixes recurring in at least MIN_SHARE of reads are candidates, and the reads carrying them
                         (again a small subset) are checked against the genome. An estimate for any transposon.

Entry points:
  from_fastq(...)   standalone, from the FASTQ (ytab donor-reads / scripts/local/ytab_donor_reads.py);
  count(...)        during mapping, reusing the SAM that already exists, so the genome pass is free.
Both write <sample>.donor_reads.csv. Insertion calls are never changed here.
"""
from __future__ import annotations

import csv
import gzip
import os
import shlex
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

PREFIX = 20
MIN_SHARE = 1e-4   # a recurring flank must hold at least 0.01% of the reads to be a donor-junction candidate
MIN_ALIGN = 30     # bases of a read that must align to donor flank sequence before it counts as a donor read
FIELDS = ["sample", "donor", "reads", "donor_reads", "donor_only_reads", "donor_shared_reads", "percent_donor", "percent_donor_only",
          "percent_donor_shared", "recurrent_flank_reads", "percent_recurrent_flank", "top_flanks", "source"]


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _sam_records(sam: str):
    with open(sam) as h:
        for line in h:
            if not line.startswith("@"):yield line.split("\t", 11)


def _run(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def _stream(cmd: str):
    """Yield a command's stdout line by line.

    An alignment against the donor returns one SAM record per aligned read, and a deeply sequenced library has tens
    of millions of them: holding that output as a string, and again as a list of lines, needs tens of gigabytes and
    is what the job gets killed for. Only one line is in memory at a time here.
    """
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, bufsize=1 << 20)
    try:
        for line in process.stdout:
            yield line
    finally:
        process.stdout.close()
        error = process.stderr.read()
        process.stderr.close()
        if process.wait():
            raise subprocess.CalledProcessError(process.returncode, cmd, stderr=error)


def _index(fasta: str, where: str) -> str:
    prefix = os.path.join(where, "ref")
    subprocess.run(["bowtie2-build", "-q", str(fasta), prefix], check=True, stdout=subprocess.DEVNULL)
    return prefix


def _sequences(fasta: str) -> list[tuple[str, str]]:
    out, name, seq = [], None, []
    with _open(fasta) as h:
        for line in h:
            if line.startswith(">"):
                if name:out.append((name, "".join(seq)))
                name, seq = line[1:].strip(), []
            else:seq.append(line.strip().upper())
    if name:out.append((name, "".join(seq)))
    return out


def mask_transposon(donor_fasta: str, transposon_fasta: str | None, where: str) -> str:
    """The donor carries the transposon, so a read that keeps the transposon end would match the donor whatever its
    flank. Masking the transposon's own sequence inside the donor leaves only the flanks to align to, which is what
    tells a donor read from a genomic one. Returns the donor to index (unchanged when no transposon is given)."""
    if not transposon_fasta or not os.path.exists(transposon_fasta):return donor_fasta
    tn = _sequences(transposon_fasta)[0][1]
    rc = tn[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    out, masked = [], 0
    for name, seq in _sequences(donor_fasta):
        for probe in (tn, rc):
            for size in (len(probe), 400, 150):   # the whole element, else its end: donors carry a mini-transposon
                for end in (probe[:size], probe[-size:]):
                    start = seq.find(end)
                    while start >= 0:
                        seq = seq[:start] + "N"*len(end) + seq[start+len(end):]; masked += len(end)
                        start = seq.find(end)
        out.append((name, seq))
    if not masked:return donor_fasta
    out = [(name, seq) for name, seq in out if set(seq) != {"N"}]   # a contig that is only the transposon has no flank
    path = os.path.join(where, "donor_masked.fasta")
    with open(path, "w") as h:
        for name, seq in out:h.write(f">{name}\n{seq}\n")
    return path


def genome_classes(sam: str, mapq: int, mate: str = "R1") -> tuple[set, Counter, int]:
    """From an existing genome SAM: (names placed at MAPQ >= mapq, first-PREFIX-base counts of unmapped reads,
    primary reads). Paired data: only the junction mate (`mate`) is counted."""
    placed, flanks, reads = set(), Counter(), 0
    for f in _sam_records(sam):
        flag = int(f[1])
        if flag & 0x900:continue            # secondary / supplementary
        if flag & 0x1 and not flag & (0x40 if mate == "R1" else 0x80):continue   # paired: the junction mate only
        reads += 1
        if flag & 0x4:flanks[f[9][:PREFIX]] += 1
        elif int(f[4]) >= mapq:placed.add(f[0])
    return placed, flanks, reads


def _aligned_bases(cigar: str) -> int:
    total = number = 0
    for ch in cigar:
        if ch.isdigit():number = number*10+int(ch)
        else:
            if ch in "M=X":total += number
            number = 0
    return total


def donor_hits(fastq: str, donor_fasta: str, threads: int, keep_reads: str | None = None,
               transposon_fasta: str | None = None) -> set:
    """Read names whose flank aligns to the donor (local alignment over masked transposon sequence, at least
    MIN_ALIGN bases). `keep_reads` also writes those reads, for the genome check."""
    with tempfile.TemporaryDirectory() as tmp:
        reference = mask_transposon(donor_fasta, transposon_fasta, tmp)
        index = _index(reference, tmp)
        al = f"--al {shlex.quote(keep_reads)} " if keep_reads else ""
        # Consumed inside the temporary directory's lifetime, because the index lives there.
        names = set()
        for line in _stream(f"bowtie2 --local -p {threads} --no-unal --no-hd {al}"
                            f"-x {shlex.quote(index)} -U {shlex.quote(fastq)}"):
            f = line.split("\t")
            if len(f) > 5 and _aligned_bases(f[5]) >= MIN_ALIGN:names.add(f[0])
    return names


def placed_on_genome(fastq: str, genome_index: str, mapq: int, threads: int) -> set:
    """Read names in `fastq` that the genome places at MAPQ >= mapq (used on small subsets only)."""
    if not genome_index or not os.path.getsize(fastq):return set()
    stream = _stream(f"bowtie2 -p {threads} --no-unal --no-hd -x {shlex.quote(genome_index)} -U {shlex.quote(fastq)}")
    return {f[0] for f in (l.split("\t") for l in stream if l) if len(f) > 4 and int(f[4]) >= mapq}


# A flank only matters at MIN_SHARE of the reads, so the table is thinned as it grows: every PRUNE_EVERY reads the
# keys seen fewer than PRUNE_BELOW times are dropped. Without this the table holds one key per distinct read start,
# which for a library merged from sixteen runs is hundreds of millions of keys and tens of gigabytes. A flank that
# survives to matter is seen far more often than the threshold, so the recurring set is unchanged.
PRUNE_EVERY = 20_000_000
PRUNE_BELOW = 2


def flank_counts(fastq: str) -> tuple[Counter, int]:
    """First-PREFIX-base counts over a FASTQ, and the number of reads."""
    counts, reads = Counter(), 0
    with _open(fastq) as h:
        for i, line in enumerate(h):
            if i % 4 == 1:
                counts[line[:PREFIX]] += 1
                reads += 1
                if reads % PRUNE_EVERY == 0:
                    for key in [k for k, n in counts.items() if n < PRUNE_BELOW]:del counts[key]
    return counts, reads


def _recurring(flanks: Counter, reads: int) -> list[tuple[str, int]]:
    return [(k, n) for k, n in flanks.most_common(50) if reads and n/reads >= MIN_SHARE and "N" not in k and len(k) == PREFIX]


def _write(row: dict, sample: str, out_dir: str) -> Path:
    for key in ("donor", "donor_only", "donor_shared", "recurrent_flank"):
        n = row.get(f"{key}_reads"); reads = row.get("reads")
        row[f"percent_{key}"] = round(100*n/reads, 4) if n is not None and reads else None
    path = Path(out_dir)/f"{sample}.donor_reads.csv"
    with open(path, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=FIELDS); w.writeheader(); w.writerow({k:row.get(k) for k in FIELDS})
    return path


def from_fastqs(fastqs: list[str], sample: str, out_dir: str, donor_fasta: str | None = None,
                genome_index: str | None = None, mapq: int = 20, threads: int = 4,
                transposon_fasta: str | None = None) -> dict:
    """Donor reads for one library whose reads are spread over several runs.

    A library merged from runs has no FASTQ of its own; its share of donor reads is the share across all the reads
    that went into it, so each run is counted and the counts are added. The flank tables are merged before the
    recurring ones are chosen, so a flank is judged against the library's reads, not one run's.
    """
    totals: dict[str, int | None] = {}
    merged: Counter = Counter()
    row: dict = {}
    for fastq in fastqs:
        part = _count(fastq, sample, donor_fasta, genome_index, mapq, threads, transposon_fasta)
        merged.update(part.pop("_flanks"))
        row = row or dict(sample=sample, donor=part["donor"], source=part["source"])
        for key in ("reads", "donor_reads", "donor_only_reads", "donor_shared_reads", "recurrent_flank_reads"):
            value = part.get(key)
            if value is None:
                totals.setdefault(key, None)
            else:
                totals[key] = (totals.get(key) or 0) + value
    row |= totals
    recurring = _recurring(merged, row.get("reads") or 0)
    row["top_flanks"] = ";".join(f"{k}:{n}" for k, n in recurring[:5])
    return row | dict(path=str(_write(row, sample, out_dir)))


def from_fastq(fastq: str, sample: str, out_dir: str, donor_fasta: str | None = None, genome_index: str | None = None,
               mapq: int = 20, threads: int = 4, transposon_fasta: str | None = None) -> dict:
    """Donor reads straight from a FASTQ, without remapping the library to the genome."""
    row = _count(fastq, sample, donor_fasta, genome_index, mapq, threads, transposon_fasta)
    row.pop("_flanks", None)
    return row | dict(path=str(_write(row, sample, out_dir)))


def _count(fastq: str, sample: str, donor_fasta: str | None, genome_index: str | None,
           mapq: int, threads: int, transposon_fasta: str | None) -> dict:
    """The counts for one FASTQ, with its flank table under `_flanks`; nothing is written."""
    flanks, reads = flank_counts(fastq)
    recurring = _recurring(flanks, reads)
    row = dict(sample=sample, donor=Path(donor_fasta).parent.name if donor_fasta else "", reads=reads,
               top_flanks=";".join(f"{k}:{n}" for k, n in recurring[:5]),
               source="donor alignment" if donor_fasta else "recurring flanks (no donor sequence)")
    with tempfile.TemporaryDirectory() as tmp:
        if donor_fasta:
            hit_fastq = os.path.join(tmp, "donor_hits.fq")
            names = donor_hits(fastq, donor_fasta, threads, hit_fastq, transposon_fasta)
            placed = placed_on_genome(hit_fastq, genome_index, mapq, threads) if genome_index else set()
            row |= dict(donor_reads=len(names), donor_shared_reads=len(names & placed) if genome_index else None,
                        donor_only_reads=len(names - placed) if genome_index else None)
        if recurring:
            # only the reads carrying a recurring flank go against the genome; those it cannot place are donor reads
            wanted = {k for k, _ in recurring}; sub = os.path.join(tmp, "flanks.fq"); kept = 0
            with _open(fastq) as h, open(sub, "w") as out:
                block = []
                for i, line in enumerate(h):
                    block.append(line)
                    if i % 4 == 3:
                        if block[1][:PREFIX] in wanted:out.write("".join(block)); kept += 1
                        block = []
            unplaced = kept-len(placed_on_genome(sub, genome_index, mapq, threads)) if genome_index else kept
            row["recurrent_flank_reads"] = unplaced
    return row | {"_flanks": flanks}


def count(sam: str, clean_fastq: str | None, sample: str, out_dir: str, mapq: int, threads: int, donor_fasta: str | None = None,
          mate: str = "R1", transposon_fasta: str | None = None) -> dict:
    """During mapping: the genome SAM already exists, so both numbers come out of it plus one donor pass."""
    placed, flanks, reads = genome_classes(sam, mapq, mate)
    recurring = _recurring(flanks, reads)
    row = dict(sample=sample, donor=Path(donor_fasta).parent.name if donor_fasta else "", reads=reads,
               recurrent_flank_reads=sum(n for _, n in recurring), top_flanks=";".join(f"{k}:{n}" for k, n in recurring[:5]),
               source="donor alignment" if donor_fasta else "recurring flanks (no donor sequence)")
    if donor_fasta and clean_fastq and os.path.exists(clean_fastq):
        names = donor_hits(clean_fastq, donor_fasta, threads, transposon_fasta=transposon_fasta)
        row |= dict(donor_reads=len(names), donor_shared_reads=len(names & placed), donor_only_reads=len(names - placed))
    _write(row, sample, out_dir)
    return row
