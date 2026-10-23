from itertools import product

import numpy as np

from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import multicondition as mc
from ytab.essentiality import resampling


def gene(i):
    return dict(id=f"g{i}", name=f"g{i}", description="", chrom="chr1", strand="+", blocks=[(1, 2)])


def mc_tolerance(p, B):
    return 5*np.sqrt(p*(1-p)/B)+3/B


def test_assignments_are_distinct_multiset_arrangements():
    A = mc.assignments([0, 0, 1, 1, 2, 2], 1000)
    assert A.shape == (90, 6) and len({tuple(r) for r in A}) == 90 and all(sorted(r) == [0, 0, 1, 1, 2, 2] for r in A.tolist())
    assert mc.assignments([0, 1, 2, 3, 4, 5, 6, 7], 1000) is None


def test_two_conditions_match_resampling_v3_exact():
    rng = np.random.default_rng(5)
    sites = rng.integers(3, 25, 40)
    C, T, genes, idx, pos = ccg.synthetic_experiment(rng, sites, (1., 1.2, .9, 1.1), treatments=2,
                                                     effects=np.where(rng.random(sites.sum()) < .15, .3, 1.))
    v3, _ = resampling.analyze(C.astype(float), T.astype(float), genes, idx, pos, resampling.Settings(recenter="trimmed-sum"))
    B = 99999
    rows, _ = mc.anova(np.column_stack([C, T]).astype(float), ["c", "c", "t", "t"], genes, idx,
                       settings=mc.Settings(stages=(B,), thresholds=(), recenter="trimmed-sum"))
    for a, b in zip(v3, rows):
        assert a["n_sites"] == b["n_sites"]
        assert abs(a["p_value"]-b["p_value"]) <= mc_tolerance(a["p_value"], B)+1e-6
        assert np.isclose(a["mean_treatment"], b["mean_t"]) and np.isclose(a["mean_control"], b["mean_c"])


def exact_three_condition_p(V, labels):
    A = mc.assignments(labels, 10**6); K = int(labels.max())+1
    per_site = [np.array([[V[s, a == k].sum() for k in range(K)] for a in A]) for s in range(len(V))]
    N = np.array([len(V)*(labels == k).sum() for k in range(K)], dtype=float)
    obs = mc._statistic(np.array([V[:, labels == k].sum() for k in range(K)]), N)
    stats = [mc._statistic(sum(per_site[s][c] for s, c in enumerate(choice)), N) for choice in product(range(len(A)), repeat=len(V))]
    return float(np.mean(np.array(stats) >= obs-1e-9*max(1., obs)))


def test_three_conditions_match_full_enumeration_both_paths():
    V = np.array([[5., 9., 2., 1., 14., 11.], [3., 4., 0., 2., 8., 10.]])
    labels = np.array([0, 0, 1, 1, 2, 2]); exact = exact_three_condition_p(V, labels)
    B = 99999
    for limit in (5000, 1):   # enumerated and sampled permutation paths
        rows, d = mc.anova(V, ["a", "a", "b", "b", "c", "c"], [gene(0)], [np.arange(2)],
                           settings=mc.Settings(stages=(B,), thresholds=(), recenter="none", max_enumerate=limit))
        assert abs(rows[0]["p_value"]-exact) <= mc_tolerance(exact, B)


def test_labels_permuted_only_within_library():
    # each library has identical values in all its samples, so every within-library arrangement gives the observed T
    V = np.array([[10., 10., 10., 0., 0., 0.], [0., 0., 0., 50., 50., 50.]])
    rows, _ = mc.anova(V, ["a", "b", "c", "a", "b", "c"], [gene(0)], [np.arange(2)], libraries=["L1"]*3+["L2"]*3,
                       settings=mc.Settings(recenter="none"))
    assert rows[0]["p_value"] == 1.


def test_staged_monte_carlo_is_deterministic_and_escalates():
    rng = np.random.default_rng(1)
    V = rng.poisson(20, (30, 6)).astype(float); V[:, 4:] *= .1
    run = lambda:mc.anova(V, ["a", "a", "b", "b", "c", "c"], [gene(0)], [np.arange(30)], settings=mc.Settings(recenter="none"))[0][0]
    first, second = run(), run()
    assert first["p_value"] == second["p_value"] and first["permutations"] == 99999 and first["top_condition"] == "c"


def test_all_pairs_one_fdr_family():
    rng = np.random.default_rng(2)
    sites = rng.integers(3, 15, 30)
    C, T, genes, idx, pos = ccg.synthetic_experiment(rng, sites, (1., 1., 1., 1.), treatments=2)
    V = np.column_stack([C, T, T[:, ::-1]]).astype(float)
    rows, diag = mc.all_pairs(V, ["c", "c", "x", "x", "y", "y"], genes, idx, pos)
    assert set(diag) == {"recenter_factors", "x vs c", "y vs c", "y vs x"} and len(rows) == 3*len(genes)
    family = [r for r in rows if r["in_family"]]
    assert np.allclose([r["q_all_pairs"] for r in family], ccg.benjamini_hochberg([r["p_value"] for r in family]))


def test_iterative_median_recovers_shift_despite_one_sided_contamination():
    rng = np.random.default_rng(4)
    genes, n = 400, 20
    clone = np.repeat(np.exp(rng.normal(3, 1, genes)), n)[:, None]*np.exp(rng.normal(0, 1, (genes*n, 1)))
    ref = clone*np.exp(rng.normal(0, .25, clone.shape))
    other = clone*np.exp(rng.normal(0, .25, clone.shape))*2**-.3   # composition shift: other is 0.3 log2 too low
    idx = [np.arange(g*n, (g+1)*n) for g in range(genes)]
    for g in range(40):other[idx[g]] *= 8.                   # 10% of genes strongly enriched
    f = mc.iterative_median_factor(ref, other, idx)
    assert abs(np.log2(f)-.3) < .03
    assert abs(np.log2(mc.trimmed_sum_factor(ref, other, idx, .2))-.3) > abs(np.log2(f)-.3)


def test_all_pairs_v2_recenters_once_to_the_reference():
    rng = np.random.default_rng(6)
    sites = rng.integers(3, 15, 40)
    C, T, genes, idx, pos = ccg.synthetic_experiment(rng, sites, (1., 1., 1., 1.), treatments=2)
    V = np.column_stack([C, T, 2*T]).astype(float)
    rows, diag = mc.all_pairs(V, ["c", "c", "x", "x", "y", "y"], genes, idx, pos)
    assert diag["recenter_factors"]["c"] == 1. and np.isclose(diag["recenter_factors"]["y"]*2, diag["recenter_factors"]["x"], rtol=1e-6)
    assert all(diag[p]["treatment_rescale"] == 1. for p in ("x vs c", "y vs c", "y vs x"))


def normalization_output(tmp_path, columns, meta):
    """Synthetic YTAB catalogue + TTR normalization run (as in tests/test_resampling.py)."""
    import csv, gzip, json
    from scipy import sparse
    from ytab.normalization.methods import run_choices
    from ytab.normalization.ttr import sha256
    rng = np.random.default_rng(7)
    sites, n_genes = 30, 40
    clone = np.exp(np.log(30)+1.5*rng.standard_normal(sites*n_genes))
    pos = 1+40*np.arange(sites*n_genes)+np.repeat(np.arange(n_genes)*1000, sites)
    raw = np.column_stack([rng.negative_binomial(5, 5/(5+clone*f(np.arange(sites*n_genes)))) for f in columns])
    gff = tmp_path/"genes.gff"; reference = tmp_path/"ref.fna"; reference.write_text(">chr1\nACGT\n")
    with gff.open("w") as handle:
        handle.write("##gff-version 3\n")
        for g in range(n_genes):
            a, b = int(pos[g*sites]), int(pos[g*sites+sites-1])+3
            handle.write(f"chr1\t.\tgene\t{a}\t{b}\t.\t+\t.\tID=g{g:02d};locus_tag=g{g:02d};Name=g{g:02d}\n")
            handle.write(f"chr1\t.\tCDS\t{a}\t{b}\t.\t+\t0\tID=cds-g{g:02d};Parent=g{g:02d}\n")
    catalogue = tmp_path/"cat"; catalogue.mkdir()
    sparse.save_npz(catalogue/"raw_counts.npz", sparse.csc_matrix(raw))
    with gzip.open(catalogue/"sites.tsv.gz", "wt") as handle:
        handle.write("site_index\tchromosome\tposition_1based\n")
        for i, p in enumerate(pos):handle.write(f"{i}\tchr1\t{p}\n")
    with (catalogue/"observations.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, list(meta[0]), delimiter="\t"); writer.writeheader(); writer.writerows(meta)
    (catalogue/"run_manifest.json").write_text(json.dumps(dict(status="success", reference_sha256=sha256(reference), universe_sites=len(pos))))
    run_choices(catalogue, tmp_path/"norm", ["ttr"], "library_id")
    return gff, reference


def test_cli_runs_three_conditions_from_one_library(tmp_path):
    import csv, json, subprocess, sys
    from pathlib import Path
    lost = lambda i:np.where((i >= 60) & (i < 90), 0., 1.)            # gene g02 loses its insertions under FLU
    one = lambda i:np.ones(len(i))
    names = [("p1", "parent", one), ("c1", "control", one), ("c2", "control", one), ("f1", "treated", lost), ("f2", "treated", lost)]
    meta = [dict(observation_id=f"P::{x}", display_name=x, project="P", library_id="P::L1", pool="L1", role=r) for x, r, _ in names]
    gff, reference = normalization_output(tmp_path, [f for *_, f in names], meta)
    script = Path(__file__).resolve().parents[1]/"scripts/local/ytab_multicondition.py"
    args = [sys.executable, str(script), "--normalized", str(tmp_path/"norm"), "--gff", str(gff), "--reference", str(reference),
            "--output", str(tmp_path/"out"), "--reference-condition", "control",
            "--sample", "P::p1=parent", "--sample", "P::c1=control", "--sample", "P::c2=control", "--sample", "P::f1=FLU", "--sample", "P::f2=FLU"]
    done = subprocess.run(args, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    summary = json.loads(done.stdout)
    assert summary["conditions"] == ["parent", "control", "FLU"] and summary["reference"] == "control"
    with (tmp_path/"out"/"anova_genes.tsv").open() as handle:anova = {r["ORF"]:r for r in csv.DictReader(handle, delimiter="\t")}
    assert anova["g02"]["call"] == "varies" and anova["g02"]["top_condition"] == "FLU"
    with (tmp_path/"out"/"pairs_genes.tsv").open() as handle:pairs = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["pair"] for r in pairs} == {"parent vs control", "FLU vs control", "FLU vs parent"}   # reference condition first
    assert next(r for r in pairs if r["ORF"] == "g02" and r["pair"] == "FLU vs control")["call_all_pairs"] == "depleted"
    manifest = json.loads((tmp_path/"out"/"run_manifest.json").read_text())
    assert manifest["algorithm_revision"] == "hermes_anova_v2" and manifest["normalization"]["method"] == "ttr"
    again = subprocess.run(args, capture_output=True, text=True)
    assert again.returncode != 0 and "new output directory" in again.stderr     # never overwrites
    bad = subprocess.run(args[:-2]+["--sample", "P::f2"], capture_output=True, text=True)
    assert bad.returncode != 0 and "OBSERVATION_ID=CONDITION" in bad.stderr


def test_runner_stratified_libraries(tmp_path):
    one = lambda i:np.ones(len(i))
    names = [("a_c", "L1", one), ("a_x", "L1", one), ("a_y", "L1", one), ("b_c", "L2", one), ("b_x", "L2", one), ("b_y", "L2", one)]
    meta = [dict(observation_id=f"P::{x}", display_name=x, project="P", library_id=f"P::{l}", pool=l, role="any") for x, l, _ in names]
    gff, reference = normalization_output(tmp_path, [f for *_, f in names], meta)
    samples = [(f"P::{x}", x.split("_")[1]) for x, *_ in names]
    manifest = mc.run(tmp_path/"norm", gff, reference, tmp_path/"out", samples, stratified=True)
    assert manifest["diagnostics"]["anova"]["libraries"] == ["P::L1", "P::L2"] and manifest["conditions"] == ["c", "x", "y"]
    import pytest
    with pytest.raises(ValueError, match="one explicit library_id|library_id"):
        mc.run(tmp_path/"norm", gff, reference, tmp_path/"out2", samples)                  # unstratified across two pools is refused
