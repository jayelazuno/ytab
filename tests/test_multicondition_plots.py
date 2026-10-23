import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ytab.essentiality import multicondition as mc
from ytab.essentiality import multicondition_plots as mp


def split_design(rng, genes=600, sites=12, a_effect=None, b_effect=None, ref_noise=.6):
    """One pool split into reference R (2 replicates), A and B (2 each); per-gene effects multiply A and B."""
    clone = np.exp(rng.normal(3, 1.2, genes*sites))
    noise = lambda:np.exp(rng.normal(0, ref_noise, genes*sites))
    gene_of = np.repeat(np.arange(genes), sites)
    fa = np.ones(genes) if a_effect is None else a_effect; fb = np.ones(genes) if b_effect is None else b_effect
    cols = [clone*noise(), clone*noise(), clone*fa[gene_of]*noise(), clone*fa[gene_of]*noise(), clone*fb[gene_of]*noise(), clone*fb[gene_of]*noise()]
    idx = [np.arange(g*sites, (g+1)*sites) for g in range(genes)]
    S, N = mc.gene_sample_sums(np.column_stack(cols), idx)
    return S, N, ["R", "R", "A", "A", "B", "B"], ["L"]*6


def test_gene_sample_sums_respects_library_eligibility():
    V = np.array([[1., 2., 0., 0.], [0., 0., 5., 1.], [3., 0., 0., 0.]])
    S, N = mc.gene_sample_sums(V, [np.array([0, 1, 2])], ["L1", "L1", "L2", "L2"])
    assert S.tolist() == [[4., 2., 5., 1.]] and N.tolist() == [[2., 2., 1., 1.]]


def test_split_reference_removes_shared_reference_correlation():
    rows = [mp.response_correlation(*split_design(np.random.default_rng(1000+s)), "R", boot=100, seed=s)[0] for s in range(8)]
    assert all(r["method"] == "split reference" and r["reference_splits"] == 2 for r in rows)
    assert np.mean([r["rho_shared_reference"] for r in rows]) > .3     # shared reference noise alone
    assert abs(np.mean([r["rho"] for r in rows])) < .03                  # split + replicate-matched: unbiased
    assert sum(r["ci_low"] <= 0 <= r["ci_high"] for r in rows) >= 7


def test_split_reference_keeps_true_sign():
    rng = np.random.default_rng(4)
    effect = np.exp(rng.normal(0, 1, 600))
    S, N, cond, libs = split_design(rng, a_effect=effect, b_effect=1/effect)
    row, = mp.response_correlation(S, N, cond, libs, "R", boot=0)
    assert row["rho"] < -.5
    S, N, cond, libs = split_design(rng, a_effect=effect, b_effect=effect)
    row, = mp.response_correlation(S, N, cond, libs, "R", boot=0)
    assert row["rho"] > .5


def test_single_reference_replicate_is_flagged():
    rng = np.random.default_rng(5)
    S, N, cond, libs = split_design(rng)
    keep = [1, 2, 3, 4, 5]
    row, = mp.response_correlation(S[:, keep], N[:, keep], [cond[j] for j in keep], [libs[j] for j in keep], "R", boot=0)
    assert row["method"].startswith("shared reference") and row["reference_splits"] == 0


def test_heatmap_rows_filtered_clustered_and_zero_centred_input():
    rng = np.random.default_rng(6)
    effect = np.ones(600); effect[:30] = 8.; effect[30:60] = 1/8
    S, N, cond, _ = split_design(rng, a_effect=effect, ref_noise=.1)
    mask = np.zeros(600, dtype=bool); mask[:60] = True
    M, rows, names = mp.heatmap_matrix(S, N, cond, "R", mask, min_sites=3)
    assert names == ["A", "B"] and sorted(rows.tolist()) == list(range(60))
    up, down = set(range(30)), set(range(30, 60))
    order = rows.tolist(); first = order[:30]
    assert set(first) in (up, down)                          # the two response groups are contiguous after clustering
    assert np.all(M[np.isin(rows, list(up)), 0] > 1) and np.all(M[np.isin(rows, list(down)), 0] < -1)


def test_cli_end_to_end(tmp_path):
    from test_multicondition import normalization_output
    lost = lambda i:np.where((i >= 60) & (i < 90), 0., 1.)
    one = lambda i:np.ones(len(i))
    names = [("c1", one), ("c2", one), ("f1", lost), ("f2", lost), ("z1", one), ("z2", one)]
    meta = [dict(observation_id=f"P::{x}", display_name=x, project="P", library_id="P::L1", pool="L1", role="any") for x, _ in names]
    gff, reference = normalization_output(tmp_path, [f for _, f in names], meta)
    samples = [("P::c1", "control"), ("P::c2", "control"), ("P::f1", "FLU"), ("P::f2", "FLU"), ("P::z1", "Zn"), ("P::z2", "Zn")]
    mc.run(tmp_path/"norm", gff, reference, tmp_path/"run", samples)
    script = Path(__file__).resolve().parents[1]/"scripts/local/ytab_multicondition_plots.py"
    done = subprocess.run([sys.executable, str(script), "--run", str(tmp_path/"run"), "--output", str(tmp_path/"plots"), "--bootstrap", "100"],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    manifest = json.loads((tmp_path/"plots"/"plots_manifest.json").read_text())
    assert manifest["correlation"][0]["method"] == "split reference"
    assert all((tmp_path/"plots"/f).is_file() for f in ("heatmap.png", "heatmap_matrix.tsv", "response_correlation.png", "response_correlation.tsv"))
    assert "g02" in (tmp_path/"plots"/"heatmap_matrix.tsv").read_text()
