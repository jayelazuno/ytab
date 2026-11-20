"""Comet table: set assignment, ortholog resolution, label choice and the plotted geometry."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ytab.fitness.comet import Settings, UNASSIGNED, assign_sets, comet_table, read_gmt


@pytest.fixture
def gmt(tmp_path):
    p = tmp_path / "sets.gmt"
    p.write_text("# a comment line is ignored\n"
                 "ERG\tErgosterol\tERG11\tERG3\n"
                 "\n"
                 "EFFLUX\tEfflux pumps\tCDR1\tPDR16\n")
    return p


@pytest.fixture
def results():
    return pd.DataFrame(dict(
        feature_id=["g1", "g2", "g3", "g4", "g5"],
        common_name=["ERG11", "", "g3", "LOCAL", "g5"],
        control_mean=[100.0, 50.0, 10.0, 4.0, 1.0],
        treated_mean=[1.0, 200.0, 10.0, 4.0, 1.0],
        final_call=["depleted", "enriched", "none", "none", "none"]))


def test_read_gmt_skips_comments_and_blanks(gmt):
    sets = read_gmt(gmt)
    assert set(sets) == {"ERG", "EFFLUX"}
    assert sets["ERG"]["members"] == {"ERG11", "ERG3"}


def test_read_gmt_rejects_a_duplicate_set(tmp_path):
    p = tmp_path / "dup.gmt"
    p.write_text("A\td\tG1\nA\td\tG2\n")
    with pytest.raises(ValueError, match="twice"):
        read_gmt(p)


def test_members_match_identifier_own_name_and_ortholog_name(gmt):
    frame = pd.DataFrame(dict(feature_id=["g1", "CDR1", "g3"],
                              display_name=["ERG11", "x", "y"],
                              ortholog_name=["", "", "PDR16"]))
    got = assign_sets(frame, read_gmt(gmt), ("ERG", "EFFLUX"))
    # g1 by its own name, CDR1 by its identifier, g3 only through its ortholog.
    assert list(got) == ["ERG", "EFFLUX", "EFFLUX"]


def test_set_priority_is_the_declared_order(tmp_path):
    p = tmp_path / "overlap.gmt"
    p.write_text("FIRST\td\tSHARED\nSECOND\td\tSHARED\n")
    sets = read_gmt(p)
    frame = pd.DataFrame(dict(feature_id=["SHARED"], display_name=[""], ortholog_name=[""]))
    assert assign_sets(frame, sets, ("FIRST", "SECOND")).iloc[0] == "FIRST"
    assert assign_sets(frame, sets, ("SECOND", "FIRST")).iloc[0] == "SECOND"


def test_unknown_set_in_order_is_an_error(gmt):
    frame = pd.DataFrame(dict(feature_id=["g"], display_name=[""], ortholog_name=[""]))
    with pytest.raises(ValueError, match="not in the GMT"):
        assign_sets(frame, read_gmt(gmt), ("ERG", "NOPE"))


def test_off_diagonal_is_the_plotted_distance(results, gmt):
    t = comet_table(results, Settings(gene_sets=gmt, label_count=0))
    row = t[t.feature_id == "g3"].iloc[0]
    assert row.off_diagonal == pytest.approx(0.0)          # equal reads sit on the diagonal
    down = t[t.feature_id == "g1"].iloc[0]
    assert down.off_diagonal < 0                           # depleted falls below it


def test_pseudocount_keeps_zero_read_genes_on_a_log_axis(gmt):
    r = pd.DataFrame(dict(feature_id=["z"], common_name=["Z"],
                          control_mean=[0.0], treated_mean=[0.0]))
    t = comet_table(r, Settings(gene_sets=gmt, label_count=0))
    assert t.x.iloc[0] > 0 and t.y.iloc[0] > 0


def test_genes_in_no_set_are_unassigned(results, gmt):
    t = comet_table(results, Settings(gene_sets=gmt, label_count=0))
    assert t[t.feature_id == "g5"].gene_set.iloc[0] == UNASSIGNED


def test_labels_require_support_in_a_condition(results, gmt):
    orth = {"g4": "ERG3"}
    # g4 is in ERG and far off the diagonal only in ratio terms, but has 4 reads: below the floor.
    loud = comet_table(results, Settings(gene_sets=gmt, label_count=5, label_min_reads=10.0),
                       orthologs=orth)
    assert not loud[loud.feature_id == "g4"].label.iloc[0]
    quiet = comet_table(results, Settings(gene_sets=gmt, label_count=5, label_min_reads=0.0),
                        orthologs=orth)
    assert quiet[quiet.feature_id == "g4"].label.iloc[0]


def test_explicit_label_genes_are_always_labelled(results, gmt):
    t = comet_table(results, Settings(gene_sets=gmt, label_count=0, label_genes=("g5",)))
    assert t[t.feature_id == "g5"].label.iloc[0]


def test_label_name_prefers_the_ortholog(results, gmt):
    t = comet_table(results, Settings(gene_sets=gmt, label_count=0, label_genes=("g3",)),
                    orthologs={"g3": "SCER1"})
    assert t[t.feature_id == "g3"].label_name.iloc[0] == "SCER1"


def test_missing_columns_are_rejected():
    with pytest.raises(ValueError, match="Missing columns|missing columns"):
        comet_table(pd.DataFrame(dict(feature_id=["g"])), Settings())


def test_labels_are_capped_per_set_not_across_all_sets(gmt):
    # Two sets of unequal size. With one shared pool the big set takes every label; per set, each gets its own.
    r = pd.DataFrame(dict(
        feature_id=["ERG11", "ERG3", "CDR1", "PDR16"],
        common_name=["", "", "", ""],
        control_mean=[100.0, 100.0, 100.0, 100.0],
        treated_mean=[1.0, 2.0, 50.0, 60.0]))      # ERG moves far, EFFLUX moves little
    t = comet_table(r, Settings(gene_sets=gmt, label_count=1))
    labelled = t[t.label]
    assert set(labelled.gene_set) == {"ERG", "EFFLUX"}     # both sets represented
    assert len(labelled) == 2                              # one each, not two from ERG


def test_direction_restricts_which_hits_are_labelled(gmt):
    r = pd.DataFrame(dict(
        feature_id=["ERG11", "ERG3"], common_name=["", ""],
        control_mean=[100.0, 1.0], treated_mean=[1.0, 100.0]))   # one down, one up
    down = comet_table(r, Settings(gene_sets=gmt, label_count=5, direction="depleted"))
    assert list(down[down.label].feature_id) == ["ERG11"]
    up = comet_table(r, Settings(gene_sets=gmt, label_count=5, direction="enriched"))
    assert list(up[up.label].feature_id) == ["ERG3"]
    both = comet_table(r, Settings(gene_sets=gmt, label_count=5, direction="both"))
    assert set(both[both.label].feature_id) == {"ERG11", "ERG3"}


def test_direction_does_not_change_which_genes_are_coloured(gmt):
    # The overlay shows where a whole family sits; hiding half of it would misrepresent the family.
    r = pd.DataFrame(dict(
        feature_id=["ERG11", "ERG3"], common_name=["", ""],
        control_mean=[100.0, 1.0], treated_mean=[1.0, 100.0]))
    for direction in ("both", "depleted", "enriched"):
        t = comet_table(r, Settings(gene_sets=gmt, label_count=5, direction=direction))
        assert (t.gene_set == "ERG").sum() == 2
