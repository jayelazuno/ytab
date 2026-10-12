"""Build a NEW parent catalogue from any number of explicit projects.

Never overwrites an accepted audit. Only included untreated parents contribute.
"""
import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import yaml

from .library import is_parent_sample, read_coding_genes
from .pool_coverage import pool_identity


def audit(repo, output, project_names):
    repo, output = Path(repo).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use a new empty output directory; accepted audits are immutable")
    if not project_names or len(set(project_names))!=len(project_names):
        raise ValueError("Provide distinct project IDs to include")
    output.mkdir(parents=True, exist_ok=True)
    source_paths = []
    def resolve(path):
        p = Path(path)
        return p if p.is_absolute() else repo / p
    projects = {}
    for name in project_names:
        config = repo / "output/projects" / name / "config/project.yaml"
        cfg = yaml.safe_load(config.read_text())
        root = resolve(cfg["output_project_dir"])
        ref_file = root / "config/reference_resolved.json"
        ref = json.loads(ref_file.read_text())
        sheet, gff = resolve(cfg["sample_sheet"]), resolve(ref["gff"])
        source_paths.extend([config, ref_file, sheet, gff])
        with sheet.open() as f:
            parents = [r for r in csv.DictReader(f) if is_parent_sample(r)
                       and r.get("include", "true").lower() in {"true", "yes", "1"}]
        projects[name] = dict(root=root, gff=gff, parents=parents)
    gffs = [p["gff"] for p in projects.values()]
    assert len({hashlib.sha256(p.read_bytes()).hexdigest() for p in gffs}) == 1
    genes = read_coding_genes(gffs[0])  # Includes CDS-bearing pseudogenes; classify explicitly below.
    biotypes, record_types, transcripts, cds_records = {}, Counter(), defaultdict(set), []
    for line in gffs[0].read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) != 9:
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        record_types[f[2]] += 1
        gid = a.get("locus_tag")
        if f[2] in {"gene", "pseudogene"}:
            biotypes[gid] = a.get("gene_biotype", f[2])
        if f[2] in {"mRNA", "transcript"}:
            transcripts[gid].add(a["ID"])
        if f[2] == "CDS":
            cds_records.append(line)
    protein_genes = {g for g in genes if biotypes.get(g) == "protein_coding"}
    pseudo_genes = {g for g in genes if biotypes.get(g) == "pseudogene"}
    assert protein_genes | pseudo_genes == set(genes)
    sites = defaultdict(dict)
    sample_details = []
    for name, project in projects.items():
        for row in project["parents"]:
            sample = row["sample"]
            path = resolve(row.get("hit_file", ""))
            if not path.is_file():
                path = project["root"] / "create_hit_file" / sample / f"{sample}_hits.txt"
            source_paths.append(path)
            key = f"{name}::{sample}"
            sample_sites = defaultdict(float)
            with path.open() as f:
                for hit in csv.DictReader(f, delimiter="\t"):
                    coord = (hit["Chromosome"], int(hit["Hit position"]))
                    reads = float(hit["Hit count"])
                    assert coord[1] > 0 and np.isfinite(reads) and reads >= 0
                    if reads > 0:
                        sample_sites[coord] += reads
            for coord, reads in sample_sites.items():
                sites[coord][key] = reads
            sample_details.append(dict(project=name, sample=sample, pool=row.get("pool_id") or row.get("pool"),
                                       physical_pool_id=pool_identity(row,hashlib.sha256(gffs[0].read_bytes()).hexdigest()[:12]),
                                       parent_replicate_provenance=row.get("background"),
                                       unique_sites=len(sample_sites), reads=sum(sample_sites.values())))
    assignment, observed = defaultdict(set), defaultdict(set)
    chrom_genes = defaultdict(list)
    for gid, g in genes.items():
        chrom_genes[g["chromosome"]].append((gid, g))
    for chrom in {c for c, _ in sites}:
        positions = np.array(sorted(p for c, p in sites if c == chrom))
        for gid, g in chrom_genes[chrom]:
            for pos in positions[np.searchsorted(positions,g["start"]):np.searchsorted(positions,g["end"],side="right")]:
                observed[gid].add((chrom,int(pos)))
            for a,b in g["blocks"]:
                for pos in positions[np.searchsorted(positions,a):np.searchsorted(positions,b,side="right")]:
                    assignment[(chrom,int(pos))].add(gid)
    unique_assignment = {coord: next(iter(gids)) for coord,gids in assignment.items() if len(gids)==1}
    all_by_gene = defaultdict(set)
    for coord,gid in unique_assignment.items():
        all_by_gene[gid].add(coord)
    project_sites = {name: {coord for coord,sources in sites.items() if any(key.startswith(name+"::") for key in sources)} for name in projects}
    def inventory(coords):
        coding = {coord for coord in coords if unique_assignment.get(coord) in protein_genes}
        coding_loci = {coord for coord in coords if coord in unique_assignment}
        gene_counts = Counter(unique_assignment[coord] for coord in coding)
        return dict(unique_sites=len(coords), unambiguous_protein_coding_sites=len(coding),
                    unambiguous_all_CDS_locus_sites=len(coding_loci),
                    protein_genes_with_coding_insertion=len(gene_counts),
                    all_CDS_loci_with_insertion=len({unique_assignment[c] for c in coding_loci}),
                    protein_genes_with_ge2_sites=sum(n>=2 for n in gene_counts.values()),
                    protein_genes_with_ge3_sites=sum(n>=3 for n in gene_counts.values()),
                    protein_genes_without_coding_insertion=len(protein_genes)-len(gene_counts),
                    protein_genes_with_observed_ORF_insertion=sum(bool(observed[g]&coords) for g in protein_genes),
                    ambiguous_CDS_sites=sum(len(assignment[c])>1 for c in coords))
    # Retain historical diagnostic columns for compatibility only; neither
    # inclusion nor the primary model depends on these two project names.
    h=project_sites.get("H2O2_screen_v1",set());z=project_sites.get("Zn_toxicity_screen",set())
    hgenes = {unique_assignment[c] for c in h if unique_assignment.get(c) in protein_genes}
    zgenes = {unique_assignment[c] for c in z if unique_assignment.get(c) in protein_genes}
    summary = dict(reference=str(gffs[0]), gene_annotation=dict(CDS_containing_loci=len(genes),
        protein_coding_genes=len(protein_genes), CDS_containing_pseudogenes=len(pseudo_genes),
        CDS_rows=len(cds_records), duplicate_CDS_rows=len(cds_records)-len(set(cds_records)),
        loci_with_multiple_transcripts=sum(len(v)>1 for g,v in transcripts.items() if g in genes),
        chromosomes=sorted(chrom_genes), feature_counts=dict(record_types)),
        samples=sample_details, by_project={name:inventory(coords) for name,coords in project_sites.items()},
        union=inventory(set(sites)), overlap=dict(shared_sites=len(h&z),H2O2_only_sites=len(h-z),Zn_only_sites=len(z-h),
            site_jaccard=len(h&z)/len(h|z) if h|z else None, shared_protein_genes=len(hgenes&zgenes),
            H2O2_only_protein_genes=len(hgenes-zgenes),Zn_only_protein_genes=len(zgenes-hgenes),
            protein_gene_jaccard=len(hgenes&zgenes)/len(hgenes|zgenes) if hgenes|zgenes else None),
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in dict.fromkeys(source_paths)])
    with (output/"parent_union_site_provenance.csv").open("w") as f:
        w=csv.writer(f);w.writerow(["chromosome","position","unambiguous_CDS_locus","biotype","ambiguous_CDS_loci","parent_sources","source_read_support_json"])
        for coord, sources in sorted(sites.items()):
            gid=unique_assignment.get(coord,"")
            w.writerow([*coord,gid,biotypes.get(gid,""),";".join(sorted(assignment[coord])) if len(assignment[coord])>1 else "",
                        ";".join(sorted(sources)),json.dumps(sources,sort_keys=True)])
    with (output/"parent_union_gene_inventory.csv").open("w") as f:
        w=csv.writer(f);w.writerow(["gene_id","biotype","chromosome","start","end","coding_length","union_distinct_CDS_sites","H2O2_distinct_CDS_sites","Zn_distinct_CDS_sites"])
        for gid,g in sorted(genes.items()):
            w.writerow([gid,biotypes[gid],g["chromosome"],g["start"],g["end"],g["length"],len(all_by_gene[gid]),len(all_by_gene[gid]&h),len(all_by_gene[gid]&z)])
    (output/"parent_union_audit.json").write_text(json.dumps(summary,indent=2))
    return summary


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo",required=True);parser.add_argument("--output",required=True)
    parser.add_argument("--projects",nargs="+",required=True)
    args=parser.parse_args()
    result=audit(args.repo,args.output,args.projects)
    print(json.dumps(result["union"],indent=2))
