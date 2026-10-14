"""Adapt MapFastq trimming/alignment to Sanger FASTA and individual mutant calls."""
from __future__ import annotations
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import fcntl
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ytab.domain_explorer.domain_insertion_explorer import load_project_gene_context, REPO_ROOT, _parse_attributes
from ytab.pipeline.mapfastq_runner import _index_complete

SUPPORTED={".fa",".fasta",".fna",".fas",".seq",".txt",".fq",".fastq",".ab1",".abi"}
MAX_BYTES=100_000_000
MAX_MEMBER=10_000_000
PRIMERS={"HR1":"AATACTTGCACTCAAAAGGCTTGACACC","HR2":"ACACTTGTGCTTATCTATGTGGC","ARB2":"GGCACGCGTCGACTAGTAC"}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def input_members(paths):
    """Read archive members in memory, with no extraction or link following."""
    total=0
    def checked(name,data):
        nonlocal total
        p=PurePosixPath(name.replace("\\","/"))
        if p.is_absolute() or ".." in p.parts or len(data)>MAX_MEMBER:
            raise ValueError(f"Unsafe or oversized input member: {name}")
        total+=len(data)
        if total>MAX_BYTES:raise ValueError("Uncompressed input exceeds 100 MB")
        return name,data
    for source in paths:
        path=Path(source).expanduser().resolve()
        if path.is_dir():
            children=sorted(p for p in path.rglob("*") if p.is_file() and not p.is_symlink()
                            and (p.suffix.lower() in SUPPORTED or p.suffix.lower() in {".gz",".tar",".tgz",".zip"}))
            if not children:raise ValueError(f"No supported sequence files in {path}")
            for name,data in input_members(children):yield checked(name,data)
        elif not path.is_file():raise FileNotFoundError(f"Input does not exist: {path}")
        elif zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as arc:
                for m in sorted(arc.infolist(),key=lambda m:m.filename):
                    if m.is_dir():continue
                    if ((m.external_attr>>16)&0o170000)==0o120000:raise ValueError("Archive links are not supported")
                    if m.file_size>MAX_MEMBER:raise ValueError("Archive member exceeds size limit")
                    if Path(m.filename).suffix.lower() in SUPPORTED:yield checked(m.filename,arc.read(m))
        elif tarfile.is_tarfile(path):
            with tarfile.open(path,"r:*") as arc:
                for m in sorted(arc.getmembers(),key=lambda m:m.name):
                    if m.issym() or m.islnk():raise ValueError("Archive links are not supported")
                    if not m.isfile():continue
                    if m.size>MAX_MEMBER:raise ValueError("Archive member exceeds size limit")
                    if Path(m.name).suffix.lower() in SUPPORTED:yield checked(m.name,arc.extractfile(m).read())
        elif path.suffix.lower()==".gz":
            with gzip.open(path,"rb") as f:yield checked(path.stem,f.read(MAX_MEMBER+1))
        elif path.suffix.lower() in SUPPORTED:
            with path.open("rb") as f:yield checked(path.name,f.read(MAX_MEMBER+1))
        else:raise ValueError(f"Unsupported input format: {path.name}")

def load_reads(paths,quality=20):
    reads=[]
    for name,data in input_members(paths):
        if Path(name).suffix.lower() in {".abi",".ab1"}:records=[SeqIO.read(io.BytesIO(data),"abi")]
        else:
            text=data.decode("utf-8-sig").strip()
            if text.startswith(">"):records=list(SeqIO.parse(io.StringIO(text),"fasta"))
            elif text.startswith("@"):records=list(SeqIO.parse(io.StringIO(text),"fastq"))
            else:
                seq="".join(text.split()).upper()
                if not re.fullmatch("[ACGTRYSWKMBDHVN]+",seq):raise ValueError(f"Not a DNA sequence file: {name}")
                records=[SeqRecord(Seq(seq),id=Path(name).stem)]
        for record in records:
            seq=str(record.seq).upper()
            if not re.fullmatch("[ACGTRYSWKMBDHVN]+",seq):raise ValueError(f"Invalid DNA characters in {name}")
            q=record.letter_annotations.get("phred_quality");a,b=0,len(seq)
            if q:
                while a<b and q[a]<quality:a+=1
                while b>a and q[b-1]<quality:b-=1
            while a<b and seq[a]=="N":a+=1
            while b>a and seq[b-1]=="N":b-=1
            reads.append(dict(read_id=f"r{len(reads)+1:05d}",mutant=record.id,source=name,
                source_member_sha256=hashlib.sha256(data).hexdigest(),raw_length=len(seq),
                sequence=seq[a:b],quality_available=bool(q),end_trimmed_bases=a+len(seq)-b,
                input_format=Path(name).suffix.lower(),pair_key=Path(name).stem.casefold() if len(records)==1 else "",
                raw_sequence=seq))
            if len(reads)>5000:raise ValueError("At most 5,000 reads per run")
    if not reads:raise ValueError("No sequence records found")
    return merge_reads(reads)

def merge_reads(reads):
    """Pair same-stem trace/text exports, never merge unrelated clones by sequence."""
    groups=defaultdict(list)
    for r in reads:
        groups[r["pair_key"] or (r["mutant"],r["source"])].append(r)
    merged=[]
    for group in groups.values():
        # Exact duplicate sources may appear via both folder and explicit input.
        unique={(r["source_member_sha256"],r["mutant"]):r for r in group}
        group=list(unique.values())
        traces=[r for r in group if r["input_format"] in {".ab1",".abi"}]
        texts=[r for r in group if r["input_format"]==".seq"]
        if len(traces)==1 and len(texts)==1 and len(group)==2:
            trace,text=traces[0],texts[0];r=dict(trace)
            same=trace["raw_sequence"]==text["raw_sequence"]
            r.update(mutant=text["mutant"],source=trace["source"]+"; "+text["source"],
                     input_merge="Paired ABI + sequence export; identical base calls" if same else
                     "Paired ABI + sequence export disagree; ABI base calls/qualities used, review required",
                     sequence_export_agrees=same)
            r["input_sources"]=[{k:v[k] for k in ("source","mutant","source_member_sha256","input_format","raw_length")} for v in (trace,text)]
            merged.append(r)
        else:
            for r in group:
                r=dict(r);r["input_merge"]="Single input; no paired ABI/sequence export"
                if len(traces)>1 or len(texts)>1:r["input_merge"]="Multiple same-stem inputs; pairing ambiguous, retained separately"
                r["input_sources"]=[{k:r[k] for k in ("source","mutant","source_member_sha256","input_format","raw_length")}]
                merged.append(r)
    merged.sort(key=lambda r:(re.sub(r"\d+",lambda m:f"{int(m[0]):020d}",r["mutant"].casefold()),r["source"]))
    for i,r in enumerate(merged,1):
        r["read_id"]=f"r{i:05d}"
        r.pop("pair_key",None);r["pre_qc_sequence"]=r.pop("raw_sequence")
    return merged

def request_fingerprint(reads,project_id,reference_hash,hermes_hash,genes,parameters,links,versions):
    payload=dict(schema="mutant-identification-explicit-cds-v5",project_id=project_id,reference=reference_hash,
                 hermes=hermes_hash,genes=genes,parameters=parameters,versions=versions,
                 inputs=sorted((Path(s["source"]).name,s["mutant"],s["source_member_sha256"])
                               for r in reads for s in r["input_sources"]),
                 linked_results=sorted((str(p),digest(p)) for p in links))
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()).hexdigest()

def cached_result(output,fingerprint):
    files=sorted(set([output/"results.json",*output.parent.glob("*/results.json")]))
    for file in files:
        try:
            data=json.loads(file.read_text())
            if not isinstance(data,dict):continue
            if data.get("request_fingerprint")!=fingerprint or data.get("status")!="success":continue
            if Path(data["output_dir"]).resolve()!=file.parent.resolve():continue
            if not all((file.parent/name).is_file() for name in ("mutant_identifications.csv","genomic_flanks.fasta")):continue
            return dict(data,reused=True)
        except (OSError,ValueError,KeyError,TypeError):continue
    return None

def command(args,output,name,commands):
    commands.append([str(a) for a in args]);r=subprocess.run(args,capture_output=True,text=True)
    (output/f"{name}.log").write_text(r.stdout+"\n"+r.stderr)
    if r.returncode:raise RuntimeError(f"{name} failed: {r.stderr[-2500:]}")
    return r.stdout

def adapters_from_reference(path,window=200):
    with Path(path).open() as f:records=list(SeqIO.parse(f,"fasta"))
    if len(records)!=1 or len(records[0])<window:raise ValueError("Construct FASTA needs one complete element with both terminal ends")
    seq=str(records[0].seq).upper()
    return {"hermes_front_right":seq[-window:],"hermes_front_left":str(Seq(seq[:window]).reverse_complement())}

def clip_reads(reads,hermes,output,commands,overlap=12,error_rate=.15,window=200):
    raw=output/"end_trimmed.fasta";clean=output/"genomic_flanks.fasta"
    with raw.open("w") as f:
        for r in reads:
            if r["sequence"]:f.write(f">{r['read_id']}\n{r['sequence']}\n")
    adapters=adapters_from_reference(hermes,window)
    if not any(r["sequence"] for r in reads):
        clean.write_text("")
        for r in reads:r.update(clean_sequence="",clipped_bases=0,clipping=[])
        return adapters
    args=[sys.executable,"-m","cutadapt","--overlap",str(overlap),"--error-rate",str(error_rate),"--times","2",
          "--info-file",str(output/"clipping.tsv"),"--json",str(output/"cutadapt.json")]
    for name,seq in adapters.items():args.extend(["-g",f"{name}={seq}","-a",f"{name.replace('front','back')}={Seq(seq).reverse_complement()}"])
    # Vendor .seq reads often lose the earliest Hermes bases. A short exact
    # terminal signature can remove the remaining tail, but does not by itself
    # authorize an exact junction call (20 matching bases are required below).
    for name,seq in adapters.items():
        args.extend(["-g",f"{name}_short={seq[-12:]};min_overlap=12;max_error_rate=0;noindels",
                     "-a",f"{name.replace('front','back')}_short={Seq(seq[-12:]).reverse_complement()};min_overlap=12;max_error_rate=0;noindels"])
    args.extend(["-a",f"ARB2_reverse={Seq(PRIMERS['ARB2']).reverse_complement()};min_overlap=17","-o",str(clean),str(raw)])
    command(args,output,"cutadapt",commands)
    with clean.open() as f:seqs={r.id:str(r.seq) for r in SeqIO.parse(f,"fasta")}
    evidence=defaultdict(list)
    for line in (output/"clipping.tsv").read_text().splitlines():
        c=line.split("\t")
        if len(c)>=8 and c[1]!="-1":evidence[c[0]].append(dict(adapter=c[7],errors=int(c[1]),start=int(c[2]),end=int(c[3])))
    for r in reads:
        r["clean_sequence"]=seqs.get(r["read_id"],"");r["clipped_bases"]=len(r["sequence"])-len(r["clean_sequence"])
        r["clipping"]=evidence[r["read_id"]]
    return adapters

def junction_position(aln,evidence,min_overlap=20):
    sides={"left" if "front" in e["adapter"] else "right" for e in evidence
           if e["adapter"].startswith("hermes_") and e["end"]-e["start"]>=min_overlap}
    if len(sides)!=1:return None
    ref_left=(next(iter(sides))=="left")!=aln.is_reverse
    if aln.cigartuples[0 if ref_left else -1][0] not in {0,7,8}:return None
    return aln.reference_start+1 if ref_left else aln.reference_end

def junction_review(aln,evidence,min_overlap=20):
    detected=[e for e in evidence if e["adapter"].startswith("hermes_")]
    overlap=max((e["end"]-e["start"] for e in detected),default=0)
    if not detected:return overlap,"Hermes end not detected; a noisy read start or construct difference is possible."
    if overlap<min_overlap:return overlap,f"Hermes overlap {overlap} bp (<{min_overlap} bp); genomic alignment retained for review."
    if aln is None:return overlap,"Hermes end detected; no unambiguous genomic alignment passes the selected thresholds."
    if junction_position(aln,evidence,min_overlap) is None:
        return overlap,"Hermes overlap meets the review threshold; local clipping or competing terminal evidence prevents an exact junction-position call."
    return overlap,f"Hermes overlap meets the {min_overlap} bp review threshold and the adjacent genomic base aligns."

def junction_status(position,evidence,min_overlap=20):
    """Separate supported Hermes sequence from resolution of its exact boundary."""
    if position is not None:return "Hermes junction supported"
    if any(e["adapter"].startswith("hermes_") and e["end"]-e["start"]>=min_overlap for e in evidence):
        return "Hermes overlap supported; exact junction position unresolved"
    return "Genomic flank mapped; junction unconfirmed"

def annotate_gene_structure(genes,path):
    """Use the selected genome's feature hierarchy, never gene span as CDS."""
    genes=[dict(g) for g in genes];by_id={g["gene_id"]:g for g in genes}
    metadata=dict(coordinates="1-based, inclusive",annotation=str(path))
    features=[];parents=defaultdict(set);roots=defaultdict(set)
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#!genome-build-accession"):
                metadata["assembly"]=line.split()[-1].split(":")[-1]
            if line.startswith("#"):continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<9:continue
            attrs=_parse_attributes(parts[8]);identifier=attrs.get("ID")
            parent_ids=attrs.get("Parent","").split(",") if attrs.get("Parent") else []
            if attrs.get("gene_id"):parent_ids.append(attrs["gene_id"])
            if identifier:parents[identifier].update(parent_ids)
            gene_id=attrs.get("locus_tag") or attrs.get("gene_id") or identifier
            if parts[2] in {"gene","pseudogene"} and gene_id in by_id:
                roots[identifier or gene_id].add(gene_id);roots[gene_id].add(gene_id)
            if parts[2]=="region" and attrs.get("strain"):metadata["strain"]=attrs["strain"]
            features.append((parts,attrs,parent_ids))
    def ancestors(identifier,seen=None):
        seen=set() if seen is None else seen
        if identifier in seen:return set()
        seen=seen|{identifier}
        return set(roots[identifier])|set().union(*(ancestors(p,seen) for p in parents[identifier]))
    for g in genes:
        g.update(cds_blocks=[],transcript_exons=[],explicit_introns=[],utr_blocks=[],structure_source=str(path))
    kinds={"CDS":"cds_blocks","exon":"transcript_exons","intron":"explicit_introns",
           "five_prime_UTR":"utr_blocks","three_prime_UTR":"utr_blocks","UTR":"utr_blocks"}
    for parts,attrs,parent_ids in features:
        if parts[2] not in kinds:continue
        direct=attrs.get("locus_tag") or attrs.get("gene_id")
        targets={direct} if direct in by_id else set()
        for parent in parent_ids:targets.update(ancestors(parent))
        for gid in targets:
            if by_id[gid]["chromosome"]==parts[0]:by_id[gid][kinds[parts[2]]].append([int(parts[3]),int(parts[4])])
    for g in genes:
        for field in kinds.values():g[field]=sorted(set(tuple(block) for block in g[field]))
    return genes,metadata

def coding_blocks(gene):
    # Explicitly empty CDS annotation means noncoding, not a missing fallback.
    return gene["cds_blocks"] if "cds_blocks" in gene else gene.get("exons") or [(gene["start"],gene["end"])]

def gene_region(gene,position):
    contains=lambda blocks:any(a<=position<=b for a,b in blocks)
    if contains(coding_blocks(gene)):return "CDS"
    if contains(gene.get("explicit_introns",[])):return "Intron"
    exons=gene.get("transcript_exons",[])
    if contains(gene.get("utr_blocks",[])):return "UTR"
    if contains(exons):return "Noncoding exon / UTR"
    if len(exons)>1 and min(a for a,b in exons)<position<max(b for a,b in exons):return "Intron"
    return "Noncoding gene span" if gene["start"]<=position<=gene["end"] else "Intergenic"

def annotated_candidates(genes,chrom,start,end,position,include_neighbors=True):
    found=[]
    for g in genes:
        if g["chromosome"]!=chrom:continue
        hit=g["start"]<=position<=g["end"] if position is not None else max(start,g["start"])<=min(end,g["end"])
        if hit:
            row={k:g.get(k,"") for k in ("gene_id","display_name","common_name","systematic_name","product","strand")}
            region=gene_region(g,position) if position is not None else None
            row.update(gene_start=g["start"],gene_end=g["end"],cds_blocks=list(coding_blocks(g)),region_type=region,
                context="Mapped flank overlaps gene; junction unconfirmed" if position is None else
                ("Coding exon" if region=="CDS" else region))
            found.append(row)
    if not found and include_neighbors:
        distance=lambda g:min(abs(position-g["start"]),abs(position-g["end"])) if position is not None else max(g["start"]-end,start-g["end"],0)
        for g in sorted((g for g in genes if g["chromosome"]==chrom),key=distance)[:2]:
            found.append(dict(gene_id=g["gene_id"],display_name=g["display_name"],product=g.get("product",""),
                context="Intergenic: nearby gene, not an insertion within this gene" if position is not None else "Mapped intergenic flank: nearby gene; junction unconfirmed",distance_bp=distance(g)))
    return found

def junction_proximal_candidates(aln,genes,evidence,position,window=100,min_bases=50,default_end="start"):
    """Assign at the junction-facing edge, not all genes crossed by a flank."""
    sides={"start" if "front" in e["adapter"] else "end" for e in evidence if e["adapter"].startswith("hermes_")}
    info=dict(window_bases=window,minimum_aligned_bases=min_bases,confirmed_junction=position is not None)
    if len(sides)>1:
        return [],dict(info,reason="Competing Hermes ends; junction-facing end is ambiguous")
    if sides:
        end=next(iter(sides));basis="Detected Hermes terminal orientation"
    elif default_end=="require":
        return [],dict(info,reason="Hermes orientation undetected; no assumed junction end")
    else:
        end=default_end;basis=f"Assumed read {end} from the selected sequencing-primer direction"
    left=(end=="start")!=aln.is_reverse
    length=aln.query_length
    pairs=[(q,r) for q,r in aln.get_aligned_pairs() if q is not None and r is not None and
           (q<window if left else q>=length-window)]
    info.update(read_end=end,orientation_basis=basis,aligned_bases=len(pairs))
    if len(pairs)<min_bases:
        return [],dict(info,reason="Insufficient aligned sequence at the junction-facing end; no gene assigned")
    reference_positions=[r+1 for q,r in pairs]
    anchor=position if position is not None else (min(reference_positions) if left else max(reference_positions))
    info.update(proximal_reference_start=min(reference_positions),proximal_reference_end=max(reference_positions),
                putative_adjacent_base=anchor,reason="Confirmed junction" if position is not None else "Provisional junction-proximal assignment; exact junction unconfirmed")
    candidates=annotated_candidates(genes,aln.reference_name,anchor,anchor,anchor,include_neighbors=False)
    info["boundary_gene_regions"]=[dict(gene_id=g["gene_id"],region=g["region_type"]) for g in candidates]
    candidates=[g for g in candidates if g["region_type"]=="CDS"]
    for g in candidates:
        support=sum(any(a<=r<=b for a,b in g["cds_blocks"]) for r in reference_positions)
        g["junction_proximal_aligned_bases"]=support
        if position is None:g["context"]="Provisional junction-proximal gene; exact junction unconfirmed"
    # An unresolved edge needs substantial same-gene support; a confirmed
    # insertion near a gene boundary must not be discarded merely for being near its end.
    if position is None:candidates=[g for g in candidates if g["junction_proximal_aligned_bases"]>=min_bases]
    if not candidates:info["reason"]="No sufficiently supported gene at the junction-facing edge; downstream/nearby genes are not assigned"
    return candidates,info

def assign_flank_candidates(aln,genes,evidence,position,window=100,min_bases=50,default_end="start"):
    """Keep single-CDS flank matches; resolve only multiple-CDS overlaps."""
    spans=annotated_candidates(genes,aln.reference_name,aln.reference_start+1,aln.reference_end,None,include_neighbors=False)
    aligned_positions=[r+1 for q,r in aln.get_aligned_pairs() if q is not None and r is not None]
    overlaps=[g for g in spans if any(any(a<=r<=b for a,b in g["cds_blocks"]) for r in aligned_positions)]
    info=dict(overlapping_cds_count=len(overlaps),confirmed_junction=position is not None,
              full_flank_gene_ids=[g["gene_id"] for g in overlaps],gene_span_ids=[g["gene_id"] for g in spans])
    if position is not None:
        # Actual junction coordinates take precedence over a flank association.
        candidates=annotated_candidates(genes,aln.reference_name,position,position,position,include_neighbors=False)
        regions=[dict(gene_id=g["gene_id"],region=g["region_type"]) for g in candidates]
        return [g for g in candidates if g["region_type"]=="CDS"],dict(info,rule="confirmed_junction",boundary_gene_regions=regions,reason="Classification at the confirmed junction base")
    if len(overlaps)<=1:
        return overlaps,dict(info,rule="single_cds_flank" if overlaps else "noncoding_gene_overlap" if spans else "no_cds_overlap",
            reason="Single-CDS flank match retained; exact junction unconfirmed" if overlaps else
            "Mapped flank does not overlap a CDS; nearby genes are not insertion targets")
    candidates,resolution=junction_proximal_candidates(aln,genes,evidence,position,window,min_bases,default_end)
    return candidates,dict(info,**resolution,rule="multiple_cds_junction_proximal")

def link_results(root,genes,selected_fitness=None):
    root=Path(root);wanted={g["gene_id"] for g in genes};out=[]
    if not wanted:return out
    if selected_fitness:
        p=Path(selected_fitness).resolve()
        if not p.is_relative_to(root.resolve()) or not p.is_file():raise ValueError("Fitness result must be an existing file in the active project")
        files=[("Fitness screen",p)]
    else:files=[("Fitness screen",p) for p in sorted((root/"treated_vs_parent").glob("*/treated_vs_parent_results.csv"))]
    files.extend(("Raw sample summary",p) for p in sorted((root/"summary").glob("*/*.feature_table.RDF_1.csv")))
    for kind,path in files:
        with path.open() as f:
            if not f.readline().startswith("RDF,"):f.seek(0)
            for row in csv.DictReader(f):
                gid=row.get("feature_id") or row.get("standard_name") or row.get("gene_id")
                if gid in wanted:
                    keys=("final_call","fitness_call","mean_log2FC","mean_log2fc","mean_z","max_abs_z","n_pairs",
                          "n_enriched_pairs","n_depleted_pairs","pools_enriched","pools_depleted","classifier_label",
                          "hits","reads","insertion_index","neighborhood_index")
                    out.append(dict(gene_id=gid,kind=kind,analysis_or_sample=path.parent.name,source=str(path),**{k:row[k] for k in keys if k in row}))
    return out

def identify(project_config,inputs,output,hermes=None,threads=2,min_length=80,min_identity=.9,min_mapq=20,
             quality=20,overlap=12,error_rate=.15,window=200,junction_min_overlap=20,fitness_result=None,condition="",pool="",
             link_project_config=None,replace_current=False,junction_gene_window=100,junction_gene_min_bases=50,junction_read_end="start"):
    output=Path(output).resolve()
    if threads<1 or min_length<20 or not 0<min_identity<=1 or min_mapq<0 or not 0<=quality<=60:raise ValueError("Invalid alignment/quality settings")
    if overlap<8 or not 0<=error_rate<=.3 or window<overlap:raise ValueError("Invalid clipping settings")
    if not 12<=junction_min_overlap<=window:raise ValueError("Junction overlap must be between 12 bp and the terminal window")
    if not 20<=junction_gene_min_bases<=junction_gene_window<=500:raise ValueError("Invalid junction-proximal window/minimum aligned bases")
    if junction_read_end not in {"start","end","require"}:raise ValueError("Invalid junction-facing read end")
    context=load_project_gene_context(project_config);ref=context["reference"]
    resolve=lambda p:Path(p) if Path(p).is_absolute() else REPO_ROOT/p
    annotation_path=resolve(ref.get("gff") or ref.get("gtf") or "")
    if not annotation_path.is_file():raise ValueError("A matching GFF/GTF with explicit CDS features is required for mutant gene assignment")
    genes,annotation_metadata=annotate_gene_structure(context["genes"],annotation_path)
    if not any(g["cds_blocks"] for g in genes):raise ValueError("No CDS features could be associated with the selected reference genes; check annotation/assembly compatibility")
    fasta=resolve(ref["fasta"]);index=resolve(ref["bowtie2_index_prefix"])
    if not fasta.is_file() or not _index_complete(index):raise ValueError("Prepare the selected reference FASTA and Bowtie2 index first")
    bowtie=shutil.which("bowtie2") or str(Path(sys.executable).parent/"bowtie2")
    if not Path(bowtie).is_file():raise ValueError("bowtie2 is missing from the environment")
    hermes=Path(hermes or REPO_ROOT/"resources/transposons/hermes/L34807.1.fasta")
    if not hermes.is_file():raise ValueError("Hermes reference FASTA is missing")
    reads=load_reads(inputs,quality)
    parameters=dict(threads=threads,min_length=min_length,min_identity=min_identity,min_mapq=min_mapq,quality=quality,
                    overlap=overlap,error_rate=error_rate,terminal_window=window,junction_min_overlap=junction_min_overlap,condition=condition,pool=pool,
                    junction_gene_window=junction_gene_window,junction_gene_min_bases=junction_gene_min_bases,junction_read_end=junction_read_end)
    versions={"bowtie2":subprocess.run([bowtie,"--version"],capture_output=True,text=True).stdout.splitlines()[0],
              "cutadapt":subprocess.run([sys.executable,"-m","cutadapt","--version"],capture_output=True,text=True).stdout.strip()}
    link_context=load_project_gene_context(link_project_config) if link_project_config else None
    root=Path(link_context["project_root"]) if link_context else None
    if link_context:
        linked_fasta=resolve(link_context["reference"]["fasta"])
        if digest(linked_fasta)!=digest(fasta):raise ValueError("Linked project must use the selected alignment reference genome")
    if fitness_result:
        if root is None:raise ValueError("Select a project explicitly before linking a Fitness result")
        fitness_result=Path(fitness_result).resolve()
        if not fitness_result.is_relative_to(root.resolve()) or not fitness_result.is_file():
            raise ValueError("Fitness result must be an existing file in the active project")
    link_files=[] if root is None else ([Path(fitness_result)] if fitness_result else sorted((root/"treated_vs_parent").glob("*/treated_vs_parent_results.csv")))
    if root is not None:link_files+=sorted((root/"summary").glob("*/*.feature_table.RDF_1.csv"))
    reference_hash,hermes_hash=digest(fasta),digest(hermes)
    fingerprint=request_fingerprint(reads,link_context["project_id"] if link_context else "universal",reference_hash,hermes_hash,genes,parameters,link_files,versions)
    cached=cached_result(output,fingerprint)
    if cached:return cached
    final_output=output
    if replace_current:
        output.parent.mkdir(parents=True,exist_ok=True)
        output=Path(tempfile.mkdtemp(prefix=".pending-",dir=output.parent))
    elif output.exists() and any(output.iterdir()):raise ValueError("Output contains a different or incomplete run; choose an empty directory")
    output.mkdir(parents=True,exist_ok=True);commands=[]
    # Measure Hermes in the original base calls without changing QC/alignment.
    # This evidence is diagnostic, not authorization for an exact junction call.
    pre_qc=output/"pre_qc_hermes";pre_qc.mkdir()
    original=[dict(read_id=r["read_id"],sequence=r["pre_qc_sequence"]) for r in reads]
    clip_reads(original,hermes,pre_qc,commands,overlap,error_rate,window)
    for r,raw in zip(reads,original):
        r["hermes_clipping_before_qc"]=[e for e in raw["clipping"] if e["adapter"].startswith("hermes_")]
        r["hermes_overlap_before_qc_bases"]=max((e["end"]-e["start"] for e in r["hermes_clipping_before_qc"]),default=0)
    adapters=clip_reads(reads,hermes,output,commands,overlap,error_rate,window)
    eligible={r["read_id"] for r in reads if len(r["clean_sequence"])>=min_length};query=output/"alignment_input.fasta"
    with query.open("w") as f:
        for r in reads:
            if r["read_id"] in eligible:f.write(f">{r['read_id']}\n{r['clean_sequence']}\n")
    alignments=defaultdict(list)
    if eligible:
        sam=output/"alignments.sam"
        # MapFastq's Bowtie2 index/read/output/thread conventions; local FASTA
        # alignment tolerates Sanger ends and -k checks alternate placements.
        command([bowtie,"--very-sensitive-local","-f","-k","5","--seed","2026","-p",str(threads),
                 "-x",str(index),"-U",str(query),"-S",str(sam)],output,"bowtie2",commands)
        with pysam.AlignmentFile(sam,"r") as f:
            for a in f:
                if not a.is_unmapped:alignments[a.query_name].append(a)
    results=[];flat=[]
    for r in reads:
        row={k:v for k,v in r.items() if k not in {"sequence","clean_sequence","pre_qc_sequence"}}
        row.update(condition=condition,pool=pool,clean_length=len(r["clean_sequence"]),candidates=[],alignments=[])
        row["hermes_overlap_bases"],row["junction_review_reason"]=junction_review(None,r["clipping"],junction_min_overlap)
        passing=[]
        for a in sorted(alignments[r["read_id"]],key=lambda a:a.get_tag("AS"),reverse=True):
            columns=sum(n for op,n in a.cigartuples if op in {0,1,2,7,8});identity=max(0,1-a.get_tag("NM")/max(columns,1))
            e=dict(chromosome=a.reference_name,start=a.reference_start+1,end=a.reference_end,strand="-" if a.is_reverse else "+",
                   mapq=None if a.mapping_quality==255 else a.mapping_quality,mapq_available=a.mapping_quality!=255,
                   identity=round(identity,5),aligned_bases=a.query_alignment_length,score=a.get_tag("AS"),cigar=a.cigarstring)
            row["alignments"].append(e)
            if identity>=min_identity and a.query_alignment_length>=min_length:passing.append((a,e))
        if r["read_id"] not in eligible:row["status"]="Insufficient genomic sequence"
        elif not passing:row["status"]="No alignment passes thresholds"
        else:
            a,e=passing[0];row.update(e)
            rivals=[b for b in alignments[r["read_id"]] if b.get_tag("AS")>=e["score"]-10 and (b.reference_name,b.reference_start)!=(a.reference_name,a.reference_start)]
            if rivals or (a.mapping_quality!=255 and a.mapping_quality<min_mapq):row["status"]="Ambiguous mapping; review alternate alignments"
            else:
                pos=junction_position(a,r["clipping"],junction_min_overlap);row["junction_adjacent_base"]=pos
                row["hermes_overlap_bases"],row["junction_review_reason"]=junction_review(a,r["clipping"],junction_min_overlap)
                row["status"]=junction_status(pos,r["clipping"],junction_min_overlap)
                evidence=r["clipping"]+r["hermes_clipping_before_qc"]
                row["candidates"],row["gene_assignment"]=assign_flank_candidates(a,genes,evidence,pos,
                    junction_gene_window,junction_gene_min_bases,junction_read_end)
        results.append(row)
        flat.append(dict(mutant=row["mutant"],input_file=row["source"],status=row["status"],condition=condition,pool=pool,
            genes="; ".join(g["gene_id"] for g in row["candidates"]),gene_names="; ".join(g["display_name"] for g in row["candidates"]),
            chromosome=row.get("chromosome",""),junction_adjacent_base=row.get("junction_adjacent_base"),strand=row.get("strand",""),
            identity=row.get("identity"),mapq=row.get("mapq"),raw_length=row["raw_length"],clean_length=row["clean_length"],aligned_bases=row.get("aligned_bases"),
            clipped_bases=row["clipped_bases"],quality_available=row["quality_available"],
            hermes_overlap_bases=row["hermes_overlap_bases"],hermes_overlap_before_qc_bases=row["hermes_overlap_before_qc_bases"],junction_review_reason=row["junction_review_reason"]))
    linked=link_results(root,[g for r in results for g in r["candidates"]],fitness_result) if root is not None else []
    for row in results:
        gids={g["gene_id"] for g in row["candidates"]}
        row["screen_results"]=[r for r in linked if r["gene_id"] in gids]
    with (output/"mutant_identifications.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    result=dict(status="success",created=datetime.now(timezone.utc).isoformat(),project_id="universal",project_config=str(Path(project_config).resolve()),
        annotation_metadata=annotation_metadata,linked_project_id=link_context["project_id"] if link_context else None,
        linked_project_config=str(Path(link_project_config).resolve()) if link_project_config else None,
        request_fingerprint=fingerprint,reused=False,
        reference=str(fasta),reference_species=context["species"],reference_sha256=reference_hash,hermes_reference=str(hermes),hermes_sha256=hermes_hash,
        inputs=[str(Path(p).resolve()) for p in inputs],input_hashes={str(Path(p).resolve()):digest(p) for p in inputs if Path(p).is_file()},
        parameters=parameters,
        adapters=adapters,versions=versions,commands=commands,results=results,output_dir=str(final_output),
        limitations=["Plain sequence files have no per-base qualities; ABI/FASTQ qualities are used only when available.",
          "A mapped flank alone does not establish an exact insertion junction or loss of function.",
          "Without Hermes end evidence, junction-facing orientation is a user-selected assumption; proximal gene assignments remain provisional.",
          "Junction position is the adjacent reference base; no unverified target-site-duplication offset is applied.",
          "Wild-type Hermes termini may differ from an engineered cassette; supply its complete sequence when necessary.",
          "Condition and pool are user-supplied; filename and gel lane ordering do not infer clone identity.",
          "Fitness links are gene-level screen associations, not measured phenotypes of this individual clone."])
    pending=output/"results.json.tmp"
    pending.write_text(json.dumps(result,indent=2));pending.replace(output/"results.json")
    if replace_current:
        previous=None
        if final_output.exists():
            archive=final_output.parent/"archive";archive.mkdir(exist_ok=True)
            previous=archive/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")+"_previous")
            final_output.rename(previous)
        try:output.rename(final_output)
        except OSError:
            if previous is not None:previous.rename(final_output)
            raise
    return result

def identify_current(**kwargs):
    """Serialize replacement; only publish a new current result after success."""
    output=Path(kwargs["output"]).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    with (output.parent/".current.lock").open("a") as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError("Another mutant identification is running; wait for it to finish")
        return identify(**kwargs,replace_current=True)
