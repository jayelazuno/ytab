"""Genome browser files for one project (the app's Gene Explorer > Genome browser reads only these).

Writes <project>/genome_browser/:
  reference/reference.fna (+ .fai)     the reference FASTA with display chromosome names, 80 bases per line
  reference/chromosome_aliases.tsv     display name, then every other name of each chromosome (for loaded tracks)
  reference/genes.bed.gz               genes from the GFF; label: common name (see gene_labels), then gene ID, then GFF name/locus/ID
  tracks/<sample>.insertions.bed.gz    one interval per insertion site, from each sample's CreateHitFile hit table
  tracks/<sample>.insertions.bedgraph.gz  the sample's insertion bedGraph, renamed to display chromosome names (if present)
  igv_session.xml                      an IGV desktop session of all of the above (paths relative to this folder)
  manifest.json                        the files above and, per sample, its track files (labels, roles and pools are
                                       sample-sheet metadata the app adds)

Chromosome names come from resources/species/<species>/reference_genome: chr_to_cp.tsv (name, contig) or
chromosome_names.tsv (contig, display, order). Gene labels come from the comparative annotation lookup or
resources/species/<species>/annotation/gene_annotation.tsv; without either, labels are the GFF's own. Files are rebuilt
only when an input is newer than the manifest.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path
from urllib.parse import unquote

VERSION = 2  # 2: tracks and genes gzip-compressed (igv.js and IGV desktop read .bed.gz directly)
LOOKUP = Path("resources/comparative/orthology/20260610_Cgla_CAGL_to_Scer_annotation_lookup.csv")
EMPTY = {"", "0", "na", "n/a", "none", "null", ".", "-"}


def _clean(value) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in EMPTY else text


def _first(*values) -> str:
    return next((v for v in map(_clean, values) if v), "")


def chromosome_map(repo_root: Path, species: str) -> list[dict]:
    """[{contig, display, order, aliases}] from the species resources (empty when there are none)."""
    ref = Path(repo_root)/"resources/species"/species/"reference_genome"
    rows = []
    if (ref/"chr_to_cp.tsv").is_file():
        for name, contig in csv.reader(open(ref/"chr_to_cp.tsv"), delimiter="\t"):
            letter = re.match(r"^Chr([A-Z])", name)
            if letter:rows.append(dict(contig=contig, display=f"Chr {letter.group(1)}", order=ord(letter.group(1))-64, aliases=[name]))
    elif (ref/"chromosome_names.tsv").is_file():
        for r in csv.DictReader(open(ref/"chromosome_names.tsv"), delimiter="\t"):
            rows.append(dict(contig=r["contig"], display=r["display"], order=int(float(r["order"])), aliases=[]))
    return sorted(rows, key=lambda r:r["order"])


def gene_labels(repo_root: Path, species: str) -> dict[str, str]:
    """feature ID (lower case) -> display label: the S. cerevisiae name for C. glabrata (ortholog lookup), the species' own
    common name for species with their own annotation; the gene ID when there is no name."""
    labels = {}
    lookup = Path(repo_root)/LOOKUP
    if lookup.is_file():
        for r in csv.DictReader(open(lookup, newline="")):
            key = _clean(r.get("gwk60_id_clean")).lower()
            if key and key not in labels:labels[key] = _first(r.get("scer_gene_name"), r.get("cagl_id"))
    other = Path(repo_root)/"resources/species"/species/"annotation/gene_annotation.tsv"
    if other.is_file():
        for r in csv.DictReader(open(other, newline=""), delimiter="\t"):
            key = _clean(r.get("gene_id")).lower()
            if key and key not in labels:labels[key] = _first(r.get("gene_name"), r.get("systematic_id"))
    return {k:v for k, v in labels.items() if v}


def _attribute(attributes: str, key: str) -> str:
    hit = re.search(rf"(^|;){key}=([^;]+)", attributes)
    return unquote(hit.group(2)) if hit else ""


def write_reference(fasta: Path, out: Path, display: dict[str, str]) -> None:
    """FASTA with display names (original name kept in the header) and its .fai index."""
    index, offset = [], 0
    with open(fasta) as source, open(out, "wb") as target:
        def flush(name, original, chunks):
            nonlocal offset
            if name is None:return
            header = f">{name} {original}\n".encode(); target.write(header); offset += len(header)
            seq = "".join(chunks); start = offset
            for i in range(0, len(seq), 80):
                line = (seq[i:i+80]+"\n").encode(); target.write(line); offset += len(line)
            index.append(f"{name}\t{len(seq)}\t{start}\t80\t81")
        name = original = None; chunks = []
        for line in source:
            if line.startswith(">"):
                flush(name, original, chunks)
                original = line[1:].split()[0]; name = display.get(original, original); chunks = []
            else:chunks.append("".join(line.split()))
        flush(name, original, chunks)
    Path(str(out)+".fai").write_text("\n".join(index)+"\n")


def write_genes(gff: Path, out: Path, display: dict[str, str], labels: dict[str, str]) -> int:
    rows = []
    for line in open(gff):
        if line.startswith("#"):continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":continue
        locus, name, gid = _attribute(f[8], "locus_tag"), _attribute(f[8], "Name"), _attribute(f[8], "ID")
        feature = _first(locus, gid, name)
        label = _first(labels.get(feature.lower()), name, locus, gid)
        if label:rows.append((display.get(f[0], f[0]), max(0, int(f[3])-1), int(f[4]), label, 0, f[6]))
    with gzip.open(out, "wt") as handle:
        for r in rows:handle.write("\t".join(map(str, r))+"\n")
    return len(rows)


def write_insertions(hits: Path, out: Path, display: dict[str, str]) -> int:
    """One BED interval per insertion site of a CreateHitFile hit table (columns Chromosome, Hit position)."""
    sites = set()
    with open(hits, newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            try:pos = int(r["Hit position"])
            except (KeyError, ValueError):continue
            if pos > 0:sites.add((display.get(r["Chromosome"], r["Chromosome"]), pos))
    with gzip.open(out, "wt") as handle:
        for chrom, pos in sorted(sites):handle.write(f"{chrom}\t{pos-1}\t{pos}\n")
    return len(sites)


def write_bedgraph(source: Path, out: Path, display: dict[str, str]) -> None:
    with open(source) as src, gzip.open(out, "wt") as dst:
        for line in src:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:dst.write("\t".join([display.get(f[0], f[0])]+f[1:4])+"\n")


def write_igv_session(out: Path, manifest: dict) -> None:
    """IGV desktop session: the display reference, the gene annotation and every insertion track (relative paths)."""
    from xml.sax.saxutils import quoteattr
    first = manifest["chromosome_order"][0] if manifest["chromosome_order"] else ""
    tracks = [f'    <Resource path={quoteattr(manifest["genes"])} name="Genes"/>']
    tracks += [f'    <Resource path={quoteattr(t["file"])} name={quoteattr(t.get("label") or t["sample"])}/>' for t in manifest["tracks"]]
    (out/"igv_session.xml").write_text('<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<Session genome={quoteattr(manifest["reference"])} locus={quoteattr(first)} version="8">\n  <Resources>\n'
        + "\n".join(tracks) + "\n  </Resources>\n</Session>\n")


def build(project_root: Path, repo_root: Path, species: str, fasta: Path, gff: Path, force: bool = False) -> dict:
    """Every sample with a CreateHitFile hit table gets a track. Returns the manifest."""
    out = Path(project_root)/"genome_browser"; manifest_path = out/"manifest.json"
    hit_root = Path(project_root)/"create_hit_file"
    samples = sorted(p.parent.name for p in hit_root.glob("*/*_hits.txt") if p.name == f"{p.parent.name}_hits.txt")
    inputs = [fasta, gff, Path(repo_root)/LOOKUP] + sorted(hit_root.glob("*/*_hits.txt")) + sorted(hit_root.glob("*/*.insertions.bedgraph"))
    newest = max((p.stat().st_mtime for p in inputs if p.is_file()), default=0)
    if not force and manifest_path.is_file() and manifest_path.stat().st_mtime >= newest:
        old = json.loads(manifest_path.read_text())
        if old.get("contract") == VERSION and [t["sample"] for t in old["tracks"]] == samples:
            return old
    (out/"reference").mkdir(parents=True, exist_ok=True); (out/"tracks").mkdir(exist_ok=True)
    chromosomes = chromosome_map(repo_root, species); display = {c["contig"]:c["display"].replace(" ", "") for c in chromosomes}
    write_reference(fasta, out/"reference/reference.fna", display)
    with open(out/"reference/chromosome_aliases.tsv", "w") as handle:
        for c in chromosomes:handle.write("\t".join(dict.fromkeys([c["display"].replace(" ", ""), c["display"], c["contig"], *c["aliases"]]))+"\n")
    genes = write_genes(gff, out/"reference/genes.bed.gz", display, gene_labels(repo_root, species))
    tracks = []
    sheet = Path(project_root)/"config/sample_sheet.csv"   # display names, when the sample sheet has them
    names = {r["sample"]:_clean(r.get("display_name")) for r in csv.DictReader(open(sheet, newline=""))} if sheet.is_file() else {}
    for sample in samples:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample)
        n = write_insertions(hit_root/sample/f"{sample}_hits.txt", out/"tracks"/f"{safe}.insertions.bed.gz", display)
        track = dict(sample=sample, label=names.get(sample) or sample, file=f"tracks/{safe}.insertions.bed.gz", format="bed", sites=n)
        bedgraph = hit_root/sample/f"{sample}.insertions.bedgraph"
        if bedgraph.is_file():
            write_bedgraph(bedgraph, out/"tracks"/f"{safe}.insertions.bedgraph.gz", display); track["bedgraph"] = f"tracks/{safe}.insertions.bedgraph.gz"
        tracks.append(track)
    manifest = dict(contract=VERSION, species=species, reference="reference/reference.fna", index="reference/reference.fna.fai",
                    aliases="reference/chromosome_aliases.tsv", genes="reference/genes.bed.gz", gene_count=genes,
                    chromosome_order=[c["display"].replace(" ", "") for c in chromosomes], tracks=tracks, igv_session="igv_session.xml")
    write_igv_session(out, manifest)
    manifest_path.write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def validate(project_root: Path) -> list[str]:
    out = Path(project_root)/"genome_browser"; path = out/"manifest.json"
    if not path.is_file():return ["missing genome_browser/manifest.json"]
    m = json.loads(path.read_text()); problems = []
    if m.get("contract") != VERSION:problems.append(f"contract {m.get('contract')} is not {VERSION}")
    for key in ("reference", "index", "aliases", "genes", "igv_session"):
        if not (out/m.get(key, "")).is_file():problems.append(f"missing {key} file")
    problems += [f"missing track {t['file']}" for t in m.get("tracks", []) if not (out/t["file"]).is_file()]
    return problems
