#!/usr/bin/env python3
"""Exercise real clipping/alignment with synthetic known genomic junctions."""
import argparse
import tempfile
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from Bio import SeqIO
from Bio.Seq import Seq
from ytab.domain_explorer.domain_insertion_explorer import load_project_gene_context, REPO_ROOT
from ytab.mutant_screen.identify import identify, identify_current, PRIMERS

p=argparse.ArgumentParser();p.add_argument("--project-config",required=True,type=Path);a=p.parse_args()
context=load_project_gene_context(a.project_config)
fasta=Path(context["reference"]["fasta"])
if not fasta.is_absolute():fasta=REPO_ROOT/fasta
with fasta.open() as handle:genome={r.id:str(r.seq).upper() for r in SeqIO.parse(handle,"fasta")}
with (REPO_ROOT/"resources/transposons/hermes/L34807.1.fasta").open() as handle:
    hermes=str(SeqIO.read(handle,"fasta").seq)
start=hermes.find(PRIMERS["HR2"])+len(PRIMERS["HR2"])
assert start>=len(PRIMERS["HR2"])
tail=hermes[start:]
chosen=None
for g in context["genes"]:
    for left,right in g.get("exons",[]):
        if right-left<900 or g["chromosome"] not in genome:continue
        pos=left+100;flank=genome[g["chromosome"]][pos-1:pos-1+500]
        if len(flank)!=500 or "N" in flank:continue
        if tail[-12:] in flank or str(Seq(tail[-12:]).reverse_complement()) in flank:continue
        if sum(s.count(flank)+s.count(str(Seq(flank).reverse_complement())) for s in genome.values())!=1:continue
        chosen=(g,pos,flank);break
    if chosen:break
assert chosen,"No unique reference segment available for the test"
g,pos,flank=chosen
with tempfile.TemporaryDirectory(prefix="ytab-junction-test-") as tmp:
    root=Path(tmp);reads=root/"known_junctions.fasta"
    forward=tail+flank
    reads.write_text(f">forward\n{forward}\n>reverse\n{Seq(forward).reverse_complement()}\n>no_Hermes\n{flank}\n>empty\nNNNN\n")
    d=identify(a.project_config,[reads],root/"run")
    repeated=identify(a.project_config,[reads],root/"duplicate-request")
    assert repeated["reused"] and repeated["output_dir"]==d["output_dir"]
    assert not (root/"duplicate-request").exists()
    results={r["mutant"]:r for r in d["results"]}
    for name in ("forward","reverse"):
        r=results[name]
        assert r["status"]=="Hermes junction supported",(name,r["status"],r["clipping"])
        assert r["junction_adjacent_base"]==pos,(name,pos,r.get("junction_adjacent_base"))
        assert g["gene_id"] in {x["gene_id"] for x in r["candidates"]}
        assert r["clipped_bases"]==len(tail),(name,r["clipped_bases"],len(tail))
    assert results["no_Hermes"]["status"]=="Genomic flank mapped; junction unconfirmed"
    assert results["no_Hermes"]["junction_adjacent_base"] is None
    assert results["empty"]["status"]=="Insufficient genomic sequence"
    assert d["project_id"]=="universal" and d["linked_project_id"] is None
    assert all(not r["screen_results"] for r in d["results"])
    linked=identify(a.project_config,[reads],root/"linked",link_project_config=a.project_config)
    assert linked["linked_project_id"]==context["project_id"]
    assert any(x["kind"]=="Raw sample summary" for r in linked["results"] for x in r["screen_results"])
    current=root/"shared/current"
    first=identify_current(project_config=a.project_config,inputs=[reads],output=current)
    assert current.is_dir() and first["project_id"]=="universal"
    assert identify_current(project_config=a.project_config,inputs=[reads],output=current)["reused"]
    updated=identify_current(project_config=a.project_config,inputs=[reads],output=current,quality=21)
    assert updated["parameters"]["quality"]==21
    assert len(list((current.parent/"archive").iterdir()))==1
    try:identify_current(project_config=a.project_config,inputs=[reads],output=current,min_identity=2)
    except ValueError:pass
    else:raise AssertionError("Invalid request unexpectedly succeeded")
    import json
    assert json.loads((current/"results.json").read_text())["parameters"]["quality"]==21
print("PASS: real Bowtie2/Cutadapt, both orientations, correct junction/gene, no fabricated junction and empty-read outcome")
