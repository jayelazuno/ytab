"""Parent hit files for the hybrid's gene features from any normalization method (ytab.pipeline.hybrid_inputs)."""
import csv
import gzip
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest
from scipy import sparse

from ytab.pipeline import hybrid_inputs as hi

HEADER = ["Chromosome", "Source", "Up feature type", "Up feature name", "Up gene name", "Up feature dist", "Down feature type",
          "Down feature name", "Down gene name", "Down feature dist", "IntergenicType", "Hit position", "Hit count"]


def _hits(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(HEADER)
        for chrom, strand, pos, count in rows:writer.writerow([chrom, strand, "gene", "g1", "G1", "5", "gene", "g2", "G2", "9", "W-W", pos, count])


def _project(tmp_path: Path, normalized: np.ndarray, raw_rows):
    """A catalogue of 3 sites x 1 observation, and a fake engine that writes the given normalized counts."""
    run = tmp_path/"project/normalization/norm_1"; cat = run/"catalogue"; cat.mkdir(parents=True)
    with gzip.open(cat/"sites.tsv.gz", "wt") as handle:
        handle.write("site_index\tchromosome\tposition_1based\n0\tc1\t10\n1\tc1\t20\n2\tc1\t30\n")
    (cat/"observations.tsv").write_text("observation_id\tsample\nP::s1\ts1\n")
    site_raw = {}
    for chrom, _, pos, count in raw_rows:site_raw[pos] = site_raw.get(pos, 0)+count
    sparse.save_npz(cat/"raw_counts.npz", sparse.csr_matrix(np.array([[site_raw.get(p, 0)] for p in (10, 20, 30)])))
    hits = tmp_path/"hits/s1_hits.txt"; _hits(hits, raw_rows)

    def engine(command, **kwargs):
        out = Path(command[command.index("--output")+1]); method = command[command.index("--method")+1]
        (out/method).mkdir(parents=True); (out/"observations.tsv").write_text("observation_id\tsample\nP::s1\ts1\n")
        sparse.save_npz(out/method/"normalized_counts.npz", sparse.csr_matrix(normalized.reshape(-1, 1)))
        return subprocess.CompletedProcess(command, 0, "", "")
    return run, hits, engine


def test_round_keep_observed():
    assert hi.round_keep_observed(np.array([.2, .5, 1.49, 1.5, 7.6])).tolist() == [1, 1, 1, 2, 8]


def test_strands_share_site_count_and_rows_keep_annotation(tmp_path, monkeypatch):
    # site 10 has W=3 and C=1 reads; its normalized count 8 is shared 6/2. Site 20 normalized to 0.2 keeps 1.
    rows = [("c1", "W", 10, 3), ("c1", "C", 10, 1), ("c1", "W", 20, 2)]
    run, hits, engine = _project(tmp_path, np.array([8., .2, 0.]), rows)
    monkeypatch.setattr(hi.subprocess, "run", engine)
    source = hi.normalize_parents(tmp_path/"project", run, "ttr", "P", ["s1"], {"s1":hits})
    assert source["tag"].startswith("ttr__norm_1__")
    written = list(csv.reader(open(tmp_path/"project/sample_normalization"/source["tag"]/"s1/s1_normalized_hits.txt"), delimiter="\t"))
    assert written[0] == HEADER and [r[-1] for r in written[1:]] == ["6", "2", "1"]
    assert [r[:-1] for r in written[1:]] == [r[:-1] for r in list(csv.reader(open(hits), delimiter="\t"))[1:]]
    assert source["totals"]["s1"]["raised_to_one"] == 1 and source["rounding"] == hi.ROUNDING
    # A second call with the same parents reuses the result without running the engine.
    monkeypatch.setattr(hi.subprocess, "run", lambda *a, **k:pytest.fail("engine ran again"))
    assert hi.normalize_parents(tmp_path/"project", run, "ttr", "P", ["s1"], {"s1":hits})["tag"] == source["tag"]


def test_hit_file_that_does_not_match_the_catalogue_is_refused(tmp_path, monkeypatch):
    run, hits, engine = _project(tmp_path, np.array([8., .2, 0.]), [("c1", "W", 10, 4), ("c1", "W", 20, 2)])
    _hits(hits, [("c1", "W", 10, 5), ("c1", "W", 20, 2)])  # the hit file changed after the catalogue was built
    monkeypatch.setattr(hi.subprocess, "run", engine)
    with pytest.raises(ValueError, match="does not match the catalogue"):
        hi.normalize_parents(tmp_path/"project", run, "ttr", "P", ["s1"], {"s1":hits})
    assert not list((tmp_path/"project/sample_normalization").glob("*/normalization_source.json"))


def test_tag_depends_on_method_run_and_parent_set():
    a = hi.method_tag("betageom", Path("x/norm_1"), ["s2", "s1"])
    assert a == hi.method_tag("betageom", Path("y/norm_1"), ["s1", "s2"])
    assert a != hi.method_tag("betageom", Path("x/norm_1"), ["s1"]) and a != hi.method_tag("ttr", Path("x/norm_1"), ["s1", "s2"])


def test_combine_resolves_a_method_tag(tmp_path):
    from ytab.pipeline.combine_hits_runner import resolve_combine_target
    tag_dir = tmp_path/"sample_normalization/ttr__norm_1__abcdef"; tag_dir.mkdir(parents=True)
    (tag_dir/hi.SOURCE_FILE).write_text(json.dumps({"method":"ttr"}))
    config = {"_repo_root":str(tmp_path), "output_project_dir":str(tmp_path), "output_export_dir":str(tmp_path/"exports")}
    info = resolve_combine_target(config, "ttr__norm_1__abcdef")
    assert info["target_tag"] == "ttr__norm_1__abcdef" and info["target"] == "ttr"
