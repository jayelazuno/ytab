"""Evidence-based cleanup of a user-designated Hermes junction read.

Mate identity is declared (default R1), not guessed from alignment rate. Cutadapt
aligners recognize terminal Hermes/primer overhang and adapter read-through.
"""
from collections import Counter
from contextlib import ExitStack
import hashlib
import json
from itertools import zip_longest
from pathlib import Path

import dnaio
from Bio import SeqIO
from Bio.Seq import Seq
from cutadapt.adapters import BackAdapter, FrontAdapter

DEFAULT_CATALOGUE = Path(__file__).resolve().parents[3]/"resources/transposons/hermes/primers.json"


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


class HermesCleaner:
    def __init__(self, catalogue=DEFAULT_CATALOGUE, reference=None, overlap=16, error_rate=.1, min_length=25):
        if overlap < 12 or not 0 <= error_rate <= .2 or min_length < 1:
            raise ValueError("Require overlap >=12, error rate in [0,.2], and min length >=1")
        self.path = Path(catalogue).resolve()
        self.catalogue = json.loads(self.path.read_text())
        self.reference = Path(reference).resolve() if reference else self.path.parent/self.catalogue["reference"]
        self.hermes = str(SeqIO.read(self.reference, "fasta").seq).upper()
        self.min_length = min_length
        self.parameters = dict(overlap=overlap, error_rate=error_rate, min_length=min_length)
        options = dict(max_errors=error_rate, min_overlap=overlap, indels=False)
        # Each sequence ends at a genuine transposon boundary in outward-read orientation.
        terminals = dict(Hermes_right=self.hermes[-120:], Hermes_left=reverse_complement(self.hermes[:120]))
        self.front = [(name, FrontAdapter(seq, **options), 0) for name,seq in terminals.items()]
        # Retained amplification/sequencing oligos may include a non-Hermes 5' tail.
        # Extend their matched Hermes suffix to the physical boundary, rather than
        # removing only the primer and accidentally leaving the 33-bp HR2 tail.
        for name in ("HR-1","HR-2-PCR","Hermes-sequencing"):
            primer=self.catalogue["primers"].get(name,"")
            for length in range(len(primer),19,-1):
                position=self.hermes.rfind(primer[-length:])
                if position >= len(self.hermes)-200:
                    sequence=primer+self.hermes[position+length:]
                    self.front.append((name,FrontAdapter(sequence,**options),len(primer)-length))
                    break
        # Complete supplied primers are useful for retained PCR/adapter sequence, but
        # only the boundary-reaching Hermes suffix determines junction-side clipping.
        # Back adapters use no short homopolymer tails or six-base barcode-only matches.
        sequences = dict(self.catalogue["adaptors"])
        sequences.update(self.catalogue["primers"])
        sequences.update(terminals)
        self.adapter_front=[(name,FrontAdapter(seq,**options)) for name,seq in sequences.items() if name not in terminals and len(seq)>=overlap]
        self.back = []
        seen = set()
        for name, seq in sequences.items():
            for orientation, text in (("forward",seq),("rc",reverse_complement(seq))):
                if text in seen or len(text)<overlap:
                    continue
                seen.add(text)
                self.back.append((f"{name}:{orientation}",BackAdapter(text, **options)))

    def provenance(self):
        return dict(parameters=self.parameters, catalogue=str(self.path), reference=str(self.reference),
                    catalogue_sha256=hashlib.sha256(self.path.read_bytes()).hexdigest(),
                    reference_sha256=hashlib.sha256(self.reference.read_bytes()).hexdigest())

    def clean(self, record, junction=True, read_end="start"):
        sequence, quality = record.sequence.upper(), record.qualities
        stats = Counter()
        if junction and read_end == "end":
            sequence, quality = reverse_complement(sequence), quality[::-1]
            stats["reverse_complemented"] += 1
        if junction:
            matches = []
            for name, adapter, hermes_start in self.front:
                match = adapter.match_to(sequence)
                # Do not mistake an internal genomic motif far into the read for its start.
                if match and match.rstart <= 80 and match.astop == len(adapter.sequence):
                    matches.append((match.score, name, match, hermes_start))
            if matches:
                _, name, match, hermes_start = max(matches,key=lambda x:x[0])
                stats["hermes_5p_reads"] += 1
                stats["hermes_overlap_bases"] += match.astop-max(match.astart,hermes_start)
                stats["trimmed_5p_bases"] += match.rstop
                stats[f"signature:{name}"] += 1
                sequence,quality = sequence[match.rstop:],quality[match.rstop:]
            else:
                # A complete retained primer/adapter at the start is removable;
                # a short prefix resemblance is not enough to delete genomic DNA.
                candidates=[]
                for name,adapter in self.adapter_front:
                    match=adapter.match_to(sequence)
                    if match and match.astart==0 and match.astop==len(adapter.sequence) and match.rstart==0:
                        candidates.append((name,match))
                if candidates:
                    name,match=max(candidates,key=lambda x:x[1].score)
                    stats["adapter_5p_reads"] += 1
                    stats["trimmed_5p_bases"] += match.rstop
                    stats[f"signature:{name}"] += 1
                    sequence,quality=sequence[match.rstop:],quality[match.rstop:]
        # Clip only demonstrated adapter read-through, not arbitrary low-quality
        # genomic bases at the junction-facing start (which would shift insertion sites).
        matches = [(name,a.match_to(sequence)) for name,a in self.back]
        matches = [(name,m) for name,m in matches if m is not None and len(set(sequence[m.rstart:m.rstop]))>=3]
        if matches:
            name, match = max(matches,key=lambda x:x[1].score)
            stats["adapter_3p_reads"] += 1
            stats["trimmed_3p_bases"] += len(sequence)-match.rstart
            stats[f"signature:{name}"] += 1
            sequence,quality=sequence[:match.rstart],quality[:match.rstart]
        if len(sequence) < self.min_length:
            stats["too_short"] += 1
            return None,stats
        return dnaio.SequenceRecord(record.name,sequence,quality),stats


def prepare_fastqs(inputs, outputs, junction_index=0, read_end="start", **options):
    if len(inputs) != len(outputs) or len(inputs) not in (1,2) or junction_index not in range(len(inputs)):
        raise ValueError("Invalid input/output mates")
    if read_end not in {"start","end"}:
        raise ValueError("Junction read end must be start or end")
    if len(inputs)==2 and read_end=="end":
        raise ValueError("Concordant mapping assumes the normal start-facing paired-read layout; use single-mate mapping for reverse-complemented exports")
    cleaner = HermesCleaner(**options)
    counts = [Counter() for _ in inputs]
    kept = 0
    with ExitStack() as stack:
        readers = [stack.enter_context(dnaio.open(str(p))) for p in inputs]
        writers = [stack.enter_context(dnaio.open(str(p),mode="w")) for p in outputs]
        for records in zip_longest(*readers):
            if any(r is None for r in records):
                raise ValueError("Paired FASTQs have unequal record counts")
            if len(records)==2:
                names=[r.name.split()[0].removesuffix("/1").removesuffix("/2") for r in records]
                if names[0] != names[1]:
                    raise ValueError(f"Paired FASTQ names disagree: {names}")
            cleaned=[]
            for index,record in enumerate(records):
                counts[index]["input_reads"] += 1
                if index==0 and counts[index]["input_reads"]%1000000==0:
                    print(f"Hermes cleanup: {counts[index]['input_reads']:,} reads/pairs examined; {kept:,} retained",flush=True)
                value,stats=cleaner.clean(record,junction=index==junction_index,read_end=read_end)
                counts[index].update(stats)
                cleaned.append(value)
            # Paired mode deliberately retains only pairs with both reads usable.
            if all(r is not None for r in cleaned):
                kept += 1
                for writer,record in zip(writers,cleaned):
                    writer.write(record)
    if not kept:
        raise ValueError("No reads remain after Hermes/adapter cleanup")
    import cutadapt
    return dict(schema_version=1, mode="designated-junction-mate", junction_index=junction_index,
                junction_read_end=read_end, retained_reads_or_pairs=kept,
                input_files=[str(p) for p in inputs], mates=[dict(c) for c in counts],
                cutadapt_version=cutadapt.__version__, **cleaner.provenance(),
                note="Absent Hermes sequence is allowed: the sequencing primer may end at the junction. No genomic strand is imposed; Bowtie2 searches both strands.")
