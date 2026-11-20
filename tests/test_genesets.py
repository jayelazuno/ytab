"""Gene sets from GO: ontology parsing, upward propagation, and what must not propagate."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "local"))

import importlib.util
spec = importlib.util.spec_from_file_location("gs", ROOT / "scripts/local/ytab_build_genesets.py")
gs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gs)


@pytest.fixture
def obo(tmp_path):
    # child --is_a--> parent --part_of--> grandparent; sibling regulates parent (must NOT be swept in);
    # old is an obsolete term; GO:0000099 is an alias of child.
    p = tmp_path / "t.obo"
    p.write_text(
        "[Term]\nid: GO:0000001\nname: grandparent\n\n"
        "[Term]\nid: GO:0000002\nname: parent\nrelationship: part_of GO:0000001\n\n"
        "[Term]\nid: GO:0000003\nname: child\nis_a: GO:0000002\nalt_id: GO:0000099\n\n"
        "[Term]\nid: GO:0000004\nname: regulator\nrelationship: regulates GO:0000002\n\n"
        "[Term]\nid: GO:0000005\nname: old\nis_obsolete: true\nis_a: GO:0000002\n")
    return p


def gaf_line(symbol, term, qualifier=""):
    f = ["SGD", f"S{symbol}", symbol, qualifier, term, "ref", "IDA", "", "P", "", "", "gene", "", "taxon:4932", "20260101", "SGD"]
    return "\t".join(f)


@pytest.fixture
def gaf(tmp_path):
    p = tmp_path / "t.gaf"
    p.write_text("!comment\n" + "\n".join([
        gaf_line("GRAND", "GO:0000001"),
        gaf_line("PAR", "GO:0000002"),
        gaf_line("CHILD", "GO:0000003"),
        gaf_line("REG", "GO:0000004"),
        gaf_line("NOPE", "GO:0000003", "NOT"),
    ]) + "\n")
    return p


def test_part_of_and_is_a_both_carry_annotations_upward(obo, gaf):
    parents, _ = gs.read_obo(obo)
    by_term = gs.read_gaf(gaf)
    members = set().union(*(by_term.get(t, set()) for t in gs.descendants(parents, {"GO:0000001"})))
    # the grandparent collects itself, the part_of parent and the is_a child
    assert members == {"GRAND", "PAR", "CHILD"}


def test_regulates_does_not_carry_annotations(obo, gaf):
    parents, _ = gs.read_obo(obo)
    by_term = gs.read_gaf(gaf)
    members = set().union(*(by_term.get(t, set()) for t in gs.descendants(parents, {"GO:0000002"})))
    # REG regulates the parent; a gene that regulates a process is not in it.
    assert "REG" not in members
    assert members == {"PAR", "CHILD"}


def test_not_qualified_annotations_are_skipped(gaf):
    # "NOT" asserts the gene is explicitly not in the term; it must never become a member.
    assert "NOPE" not in gs.read_gaf(gaf).get("GO:0000003", set())


def test_obsolete_terms_are_dropped(obo):
    parents, _ = gs.read_obo(obo)
    assert "GO:0000005" not in parents


def test_alt_ids_resolve_to_the_current_term(obo):
    _, alias = gs.read_obo(obo)
    assert alias["GO:0000099"] == "GO:0000003"


def test_descendants_includes_the_term_itself(obo):
    parents, _ = gs.read_obo(obo)
    assert "GO:0000003" in gs.descendants(parents, {"GO:0000003"})


def test_gmt_is_written_with_name_description_members(tmp_path, obo, gaf):
    out = tmp_path / "out.gmt"
    sets = tmp_path / "sets.tsv"
    sets.write_text("# a comment\nTEST\tGO:0000001\tA description\n")
    gs.main([str(gaf), str(obo), str(out), "--sets", str(sets), "--min-genes", "1"])
    body = [l for l in out.read_text().splitlines() if not l.startswith("#")]
    assert len(body) == 1
    name, description, *members = body[0].split("\t")
    assert name == "TEST"
    assert "GO:0000001" in description
    assert sorted(members) == ["CHILD", "GRAND", "PAR"]


def test_a_set_below_min_genes_is_skipped(tmp_path, obo, gaf):
    out = tmp_path / "out.gmt"
    sets = tmp_path / "sets.tsv"
    sets.write_text("TINY\tGO:0000004\n")        # only REG, one gene
    gs.main([str(gaf), str(obo), str(out), "--sets", str(sets), "--min-genes", "5"])
    assert [l for l in out.read_text().splitlines() if not l.startswith("#")] == []


def test_sets_file_needs_a_go_term(tmp_path, obo, gaf):
    sets = tmp_path / "bad.tsv"
    sets.write_text("JUSTANAME\n")
    with pytest.raises(SystemExit, match="name and GO term"):
        gs.main([str(gaf), str(obo), str(tmp_path / "o.gmt"), "--sets", str(sets)])
