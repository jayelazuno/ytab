import numpy as np

from ytab.essentiality import library_correlation as lc


def write_gff(path, chrom, genes):
    with path.open("w") as h:
        h.write("##gff-version 3\n")
        for gid, a, b in genes:
            h.write(f"{chrom}\t.\tgene\t{a}\t{b}\t.\t+\t.\tID=gene:{gid};Name={gid}\n")
            h.write(f"{chrom}\t.\tCDS\t{a}\t{b}\t.\t+\t0\tID=cds-{gid};Parent=gene:{gid}\n")


def write_hits(path, chrom, positions, reads, ytab=True):
    with path.open("w") as h:
        if ytab:
            h.write("Chromosome\tSource\tUp feature type\tUp feature name\tUp gene name\tUp feature dist\tDown feature type\tDown feature name\t"
                    "Down gene name\tDown feature dist\tIntergenicType\tHit position\tHit count\n")
            for p, r in zip(positions, reads):h.write(f"{chrom}\tW\t\t\t\t\t\t\t\t\t\t{p}\t{r}\n")
        else:
            h.write("Chromosome,Strand,Position,Reads\n")
            for p, r in zip(positions, reads):h.write(f"{chrom},W,{p},{r}\n")


def library(rng, essential, clone, sites, keep=1.):
    """Insertions at sites outside essential genes; reads from shared (clone) or fresh abundances."""
    ok = ~np.isin(sites//1000, essential) & (rng.random(len(sites)) < keep)
    return sites[ok], np.maximum(1, rng.poisson(clone[ok])).astype(int)


def test_pairs_typed_matched_across_species_and_clone_evidence(tmp_path):
    rng = np.random.default_rng(1)
    genes = [(f"G{g:03d}", g*1000+1, g*1000+900) for g in range(300)]
    write_gff(tmp_path/"a.gff", "chrA", genes); write_gff(tmp_path/"b.gff", "chrB", [(f"P{g:03d}", a, b) for g, (_, a, b) in enumerate(genes)])
    essential = rng.choice(300, 60, replace=False); sites = np.arange(1, 300_000, 37)
    clone = np.exp(rng.normal(3, 1.5, len(sites)))
    samples = []
    for name, ab, keep in (("split1", clone, .9), ("split2", clone, .9), ("independent", np.exp(rng.normal(3, 1.5, len(sites))), .9)):
        pos, reads = library(rng, essential, ab, sites, keep); write_hits(tmp_path/f"{name}.txt", "chrA", pos, reads)
        samples.append(lc.Sample(name, str(tmp_path/f"{name}.txt"), "spA", str(tmp_path/"a.gff"), "S1", "p1" if name != "independent" else "p2", "parent"))
    pos, reads = library(rng, essential, clone, sites); write_hits(tmp_path/"other.csv", "chrB", pos, reads, ytab=False)
    samples.append(lc.Sample("other", str(tmp_path/"other.csv"), "spB", str(tmp_path/"b.gff"), "S9", "1", "parent"))
    with (tmp_path/"orth.tsv").open("w") as h:
        h.write("group\tspecies\tgene\n")
        for g in range(250):h.write(f"OG{g}\tspA\tG{g:03d}\nOG{g}\tspB\tP{g:03d}\n")
        h.write("OGdup\tspA\tG299\nOGdup\tspA\tG298\nOGdup\tspB\tP299\n")      # not 1:1: excluded
    manifest = lc.run(samples, tmp_path/"out", orthologs=tmp_path/"orth.tsv", boot=50)
    import csv
    rows = {(r["sample_a"], r["sample_b"]):r for r in csv.DictReader((tmp_path/"out"/"pairs.tsv").open(), delimiter="\t")}
    split = rows[("split1", "split2")]; indep = rows[("split1", "independent")]; cross = rows[("split1", "other")]
    assert split["pair_type"] == "same pool (shared clones)" and indep["pair_type"] == "same strain, different pool"
    assert cross["pair_type"] == "different species" and int(cross["genes"]) == 250
    assert split["clone_evidence"] == "shared clones" and indep["clone_evidence"] == "no shared clone abundance"
    assert min(float(split["rho"]), float(indep["rho"]), float(cross["rho"])) > .4    # shared essential genes
    assert not split["flag"] and not indep["flag"] and manifest["flagged_pairs"] == []


def test_conflicting_declaration_is_flagged():
    assert lc.conflict("same pool (shared clones)", "no shared clone abundance").startswith("declared same pool")
    assert lc.conflict("different strain", "shared clones").startswith("declared independent")
    assert lc.conflict("same pool (shared clones)", "ambiguous") == ""


def test_gene_key_strips_annotation_prefix():
    assert lc.gene_key("gene:SPAC1002.01") == "SPAC1002.01" and lc.gene_key("GWK60_A00011") == "GWK60_A00011"
