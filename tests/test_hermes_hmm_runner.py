import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]/"scripts/local"))
import ytab_hermes_hmm as runner
from ytab.essentiality import hermes_hmm as hh


def genome(tmp_path, names, length=40_000, duplicate=False, seed=0):
    rng = np.random.default_rng(seed); seq = "".join(rng.choice(list("ACGT"), length))
    with (tmp_path/"g.fna").open("w") as h:
        for i, n in enumerate(names):h.write(f">{n} description\n{seq if duplicate or i == 0 else ''.join(rng.choice(list('ACGT'), length))}\n")
    return tmp_path/"g.fna"


def gff(tmp_path, entries):
    with (tmp_path/"g.gff").open("w") as h:
        h.write("##gff-version 3\n")
        for gid, chrom, a, b in entries:
            h.write(f"{chrom}\t.\tgene\t{a}\t{b}\t.\t+\t.\tID={gid};Name={gid}\n{chrom}\t.\tCDS\t{a}\t{b}\t.\t+\t0\tID=cds-{gid};Parent={gid}\n")
    return tmp_path/"g.gff"


def hits(tmp_path, name, chrom_positions):
    path = tmp_path/f"{name}_hits.txt"
    with path.open("w") as h:
        h.write("Chromosome\tPosition\tReads\n")
        for c, ps in chrom_positions.items():h.writelines(f"{c}\t{p}\t5\n" for p in ps)
    return path


@pytest.fixture(autouse=True)
def private_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "CACHE", tmp_path/"cache")


def test_feature_cache_is_keyed_by_genome(tmp_path):
    settings = hh.Settings()
    a = genome(tmp_path, ["c1"], seed=1); fa = runner.features(settings, [], a)
    (tmp_path/"b").mkdir(); b = genome(tmp_path/"b", ["c1"], seed=2); fb = runner.features(settings, [], b)
    assert not np.array_equal(fa["c1"]["at"], fb["c1"]["at"]) and len(list((tmp_path/"cache").glob("features_*"))) == 2


def test_duplicated_haplotypes_are_refused(tmp_path):
    fasta = genome(tmp_path, ["chr1A", "chr1B"], duplicate=True); annotation = gff(tmp_path, [("g1", "chr1A", 1001, 3000)])
    h = hits(tmp_path, "lib", {"chr1A": range(10, 40_000, 37)})
    with pytest.raises(ValueError, match="unmappable"):
        runner.run(tmp_path/"out", ["lib"], hh.Settings(), "t", fasta=fasta, gff=annotation, hits=[h])


def test_chromosome_name_mismatch_is_refused(tmp_path):
    fasta = genome(tmp_path, ["chrI"]); annotation = gff(tmp_path, [("g1", "chrI", 1001, 3000)])
    h = hits(tmp_path, "lib", {"I": range(10, 40_000, 37)})
    with pytest.raises(ValueError, match="not in the FASTA"):
        runner.run(tmp_path/"out", ["lib"], hh.Settings(), "t", fasta=fasta, gff=annotation, hits=[h])


def test_genes_off_reference_are_excluded_and_reported(tmp_path):
    fasta = genome(tmp_path, ["chr1A"], length=60_000)
    genes = [(f"A{i}", "chr1A", 2000+i*3000, 3500+i*3000) for i in range(15)]+[("B0", "chr1B", 1001, 3000)]
    annotation = gff(tmp_path, genes)
    rng = np.random.default_rng(3); positions = np.sort(rng.choice(np.arange(1, 60_000), 3000, replace=False))
    essential = (positions >= 2000) & (positions <= 3500)                    # gene A0 has no insertions
    h1 = hits(tmp_path, "lib1", {"chr1A": positions[~essential]}); h2 = hits(tmp_path, "lib2", {"chr1A": positions[~essential][::2]})
    summary, rows = runner.run(tmp_path/"out", ["lib1", "lib2"], hh.Settings(max_iterations=30), "t", fasta=fasta, gff=annotation, hits=[h1, h2])
    assert {r["ORF"] for r in rows} == {f"A{i}" for i in range(15)}
    assert any("1 genes on 1 chromosomes absent from the FASTA" in n for n in summary["notes"])
    assert summary["inputs"]["fasta_sha256"] and (tmp_path/"out"/"genes.tsv").is_file()
