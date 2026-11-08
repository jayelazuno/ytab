import pytest

from ytab.normalization import site_filter

GENES = [dict(id="gene:A", name="A", chrom="c1", strand="+", blocks=[(1, 100)]),
         dict(id="gene:B", name="B", chrom="c1", strand="+", blocks=[(201, 300)])]
COUNTS = {("c1", 10):900, ("c1", 20):100, ("c1", 30):40, ("c1", 210):50, ("c1", 500):7}


def test_none_keeps_counts():
    out, s = site_filter.apply(COUNTS, GENES, "none")
    assert out == COUNTS and s["genes_changed"] == 0


def test_winsorize_caps_top_site_at_second():
    out, s = site_filter.apply(COUNTS, GENES, "winsorize")
    assert out[("c1", 10)] == 100 and out[("c1", 210)] == 50 and out[("c1", 500)] == 7   # one-site gene and intergenic unchanged
    assert s == dict(site_filter="winsorize", genes_changed=1, reads_removed=800)


def test_drop_top_removes_top_site():
    out, s = site_filter.apply(COUNTS, GENES, "drop-top")
    assert ("c1", 10) not in out and ("c1", 210) not in out and out[("c1", 20)] == 100 and out[("c1", 500)] == 7
    assert s["genes_changed"] == 2 and s["reads_removed"] == 950


def test_unknown_filter():
    with pytest.raises(ValueError):site_filter.apply(COUNTS, GENES, "trim")
