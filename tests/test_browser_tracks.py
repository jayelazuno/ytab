import gzip
from pathlib import Path

from ytab.domain_explorer import browser_tracks as bt


def _project(tmp_path):
    repo, project = tmp_path/"repo", tmp_path/"project"
    ref = repo/"resources/species/testspecies/reference_genome"; ref.mkdir(parents=True)
    (ref/"chromosome_names.tsv").write_text("contig\tdisplay\torder\nctg1\tChr I\t1\n")
    fasta = tmp_path/"ref.fna"; fasta.write_text(">ctg1 test\n" + "ACGT"*50 + "\n")
    gff = tmp_path/"ref.gff"
    gff.write_text("ctg1\tx\tgene\t11\t60\t.\t+\t.\tID=gene-G1;Name=ABC1;locus_tag=G1\n"
                   "ctg1\tx\tgene\t81\t120\t.\t-\t.\tID=gene-G2;locus_tag=G2\n")
    hits = project/"create_hit_file/s1"; hits.mkdir(parents=True)
    (hits/"s1_hits.txt").write_text("Chromosome\tSource\tHit position\tHit count\nctg1\tW\t15\t3\nctg1\tC\t15\t2\nctg1\tW\t90\t1\n")
    (project/"config").mkdir(); (project/"config/sample_sheet.csv").write_text("sample,display_name\ns1,Sample one\n")
    return repo, project, fasta, gff


def test_browser_files_for_any_species(tmp_path):
    repo, project, fasta, gff = _project(tmp_path)
    m = bt.build(project, repo, "testspecies", fasta, gff)
    out = project/"genome_browser"
    assert bt.validate(project) == []
    assert (out/m["reference"]).read_text().startswith(">ChrI ctg1")          # display names from the species resource
    assert (out/m["index"]).read_text().split("\t")[:2] == ["ChrI", "200"]
    genes = [line.split("\t") for line in gzip.open(out/m["genes"], "rt").read().splitlines()]
    assert [g[3] for g in genes] == ["ABC1", "G2"] and genes[0][:3] == ["ChrI", "10", "60"]   # GFF label fallback
    assert gzip.open(out/m["tracks"][0]["file"], "rt").read().splitlines() == ["ChrI\t14\t15", "ChrI\t89\t90"]  # one row per site
    session = (out/m["igv_session"]).read_text()
    assert 'genome="reference/reference.fna"' in session and 'name="Sample one"' in session
    assert bt.build(project, repo, "testspecies", fasta, gff) == m             # up to date: not rebuilt


def test_saved_gene_queries_are_capped(tmp_path):
    import os, time
    from ytab.domain_explorer.domain_insertion_explorer import _limit_saved_queries
    base = tmp_path/"gene_domain_explorer"
    for d in ("manifests", "figures", "tables"):(base/d).mkdir(parents=True)
    for i in range(5):
        stem = f"G{i}.raw.abc{i}"
        (base/"manifests"/f"{stem}.json").write_text("{}"); (base/"figures"/f"{stem}.png").write_text("x"); (base/"tables"/f"{stem}.insertions.csv").write_text("x")
        os.utime(base/"manifests"/f"{stem}.json", (time.time()-100+i, time.time()-100+i))
    os.utime(base/"manifests"/"G0.raw.abc0.json")           # G0 was just reused
    assert _limit_saved_queries(base, keep=3) == 2
    kept = sorted(p.name.split(".")[0] for p in (base/"manifests").glob("*.json"))
    assert kept == ["G0", "G3", "G4"] and not (base/"figures"/"G1.raw.abc1.png").exists() and not (base/"tables"/"G2.raw.abc2.insertions.csv").exists()
