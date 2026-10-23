import csv
import gzip
import itertools
import json

import numpy as np
import pytest
from scipy import sparse

from ytab.essentiality.control_conditioned_gumbel import synthetic_experiment
from ytab.essentiality.resampling import Settings, analyze, exact_p, monte_carlo_p, run, subset_sums
from ytab.normalization.methods import normalize, run_choices
from ytab.normalization.ttr import sha256

EXPOSURES = np.array([1., 1.3, .8, 1.1])


def brute_force_p(control, treatment):
    """Enumerate every site-restricted label assignment."""
    J, K = control.shape[1], treatment.shape[1]; values = np.concatenate([control, treatment], 1)
    t_obs = treatment.mean()-control.mean(); per_site = list(itertools.combinations(range(J+K), K)); hits = total = 0
    for choice in itertools.product(per_site, repeat=len(values)):
        mask = np.zeros(values.shape, dtype=bool)
        for i, subset in enumerate(choice):mask[i, list(subset)] = True
        t = values[mask].mean()-values[~mask].mean(); total += 1; hits += abs(t) >= abs(t_obs)-1e-9
    return hits/total, t_obs


def test_settings_guard():
    for bad in (dict(inference="pooled"), dict(opportunities="all"), dict(samples=10), dict(grid_bins=10), dict(fdr=0)):
        with pytest.raises(ValueError):Settings(**bad).validate()
    assert Settings().inference == "exact" and Settings().samples == 10000


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_exact_matches_enumeration_and_bounds_it(seed):
    rng = np.random.default_rng(seed)
    control = rng.poisson(20, (5, 2))*1.37; treatment = rng.poisson(12, (5, 2))*.91   # non-lattice scaled counts
    sums, observed = subset_sums(control, treatment)
    upper, rounded, t_obs = exact_p(sums, observed, 2, 2, 65536)
    truth, t_true = brute_force_p(control, treatment)
    assert t_obs == pytest.approx(t_true) and truth <= upper+1e-12 and rounded == pytest.approx(truth) and upper-truth < .01


def test_unequal_group_sizes_and_monte_carlo_agree_with_exact():
    rng = np.random.default_rng(4)
    control = rng.poisson(30, (40, 3))*1.1; treatment = rng.poisson(31, (40, 2))*1.02   # modest effect: p within Monte Carlo resolution
    sums, observed = subset_sums(control, treatment)
    assert sums.shape == (40, 10)
    upper, rounded, t_obs = exact_p(sums, observed, 3, 2, 65536)
    mc, t_mc = monte_carlo_p(np.random.default_rng(0), sums, observed, 3, 2, 200000)
    assert t_mc == pytest.approx(t_obs) and t_obs == pytest.approx(treatment.mean()-control.mean())
    assert abs(mc-rounded) < 4*np.sqrt(rounded*(1-rounded)/200000) and rounded <= upper < rounded+.005


def test_exact_ties_count_toward_p():
    # Single site: all six label assignments give |T| >= |T_obs|, so the exact p is 1 (ties included).
    upper, rounded, t_obs = exact_p(*subset_sums(np.array([[10., 0.]]), np.array([[0., 3.]])), 2, 2, 65536)
    assert t_obs == -3.5 and upper == pytest.approx(1) and rounded == pytest.approx(1)


def test_identical_values_give_p_one():
    values = np.full((6, 2), 5.)
    assert exact_p(*subset_sums(values, values), 2, 2, 1024)[0] == 1


def experiment(seed, effects=None, genes=120, sites=30, **model):
    rng = np.random.default_rng(seed)
    C, T, genes, idx, pos = synthetic_experiment(rng, np.full(genes, sites), EXPOSURES, effects=effects, treatments=2, **model)
    scaled, _ = normalize(sparse.csc_matrix(np.column_stack([C, T])), "ttr", ["g"]*4)
    scaled = scaled.toarray()
    return scaled[:, :2], scaled[:, 2:], genes, idx, pos


def test_null_calibration_and_label_invariant_opportunities():
    C, T, genes, idx, pos = experiment(11, genes=200, m=np.log(3))
    for inference in ("exact", "monte-carlo"):
        rows, diagnostics = analyze(C, T, genes, idx, pos, Settings(inference=inference, samples=2000))
        p = np.array([r["p_value"] for r in rows if r["in_family"]])
        assert np.mean(p <= .05) < .1 and diagnostics["calls"]["depleted"]+diagnostics["calls"]["enriched"] <= 2
    # Swapping the condition labels mirrors any-support results exactly (label invariance).
    forward = analyze(C, T, genes, idx, pos, Settings())[0]; reverse = analyze(T, C, genes, idx, pos, Settings())[0]
    assert [r.get("p_value") for r in forward] == pytest.approx([r.get("p_value") for r in reverse])


def test_detects_depletion_and_enrichment_and_exact_is_seed_free():
    effects = np.ones(120*30); effects[:30*10] = .2; effects[30*60:30*70] = 5.
    C, T, genes, idx, pos = experiment(21, effects=effects)
    rows = analyze(C, T, genes, idx, pos, Settings())[0]
    depleted = {r["ORF"] for r in rows if r["call"] == "depleted"}; enriched = {r["ORF"] for r in rows if r["call"] == "enriched"}
    assert len(depleted & {genes[g]["id"] for g in range(10)}) >= 8 and len(enriched & {genes[g]["id"] for g in range(60, 70)}) >= 8
    other = analyze(C, T, genes, idx, pos, Settings(seed=99))[0]
    assert [r.get("p_value") for r in rows] == [r.get("p_value") for r in other]
    mc = [analyze(C, T, genes, idx, pos, Settings(inference="monte-carlo", samples=1000, seed=s))[0] for s in (1, 1, 2)]
    assert [r.get("p_value") for r in mc[0]] == [r.get("p_value") for r in mc[1]]
    assert [r.get("p_value") for r in mc[0]] != [r.get("p_value") for r in mc[2]]


def test_file_runner_on_normalization_output(tmp_path):
    effects = np.ones(40*30); effects[30*2:30*3] = 0.
    rng = np.random.default_rng(61)
    C, T, genes, idx, pos = synthetic_experiment(rng, np.full(40, 30), EXPOSURES, effects=effects, treatments=2)
    sites = [("chr1", int(p)) for p in pos]
    gff = tmp_path/"genes.gff"; reference = tmp_path/"ref.fna"; reference.write_text(">chr1\nACGT\n")
    with gff.open("w") as handle:
        handle.write("##gff-version 3\n")
        for g in genes:
            a, b = g["blocks"][0]
            handle.write(f"chr1\t.\tgene\t{a}\t{b}\t.\t+\t.\tID={g['id']};locus_tag={g['id']};Name={g['id']}\n")
            handle.write(f"chr1\t.\tCDS\t{a}\t{b}\t.\t+\t0\tID=cds-{g['id']};Parent={g['id']}\n")
    catalogue = tmp_path/"cat"; catalogue.mkdir()
    sparse.save_npz(catalogue/"raw_counts.npz", sparse.csc_matrix(np.column_stack([C, T])))
    with gzip.open(catalogue/"sites.tsv.gz", "wt") as handle:
        handle.write("site_index\tchromosome\tposition_1based\n")
        for i, (chrom, p) in enumerate(sites):handle.write(f"{i}\t{chrom}\t{p}\n")
    meta = [dict(observation_id=f"P::{x}", display_name=x, project="P", library_id="P::L1", pool="L1", role=r)
            for x, r in (("c1", "control"), ("c2", "control"), ("t1", "treated"), ("t2", "treated"))]
    with (catalogue/"observations.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, list(meta[0]), delimiter="\t"); writer.writeheader(); writer.writerows(meta)
    (catalogue/"run_manifest.json").write_text(json.dumps(dict(status="success", reference_sha256=sha256(reference), universe_sites=len(sites))))
    run_choices(catalogue, tmp_path/"norm", ["ttr", "quantile"], "library_id")
    manifest = run(tmp_path/"norm", gff, reference, tmp_path/"out", ["P::c1", "P::c2"], ["P::t1", "P::t2"])
    assert manifest["normalization"]["method"] == "ttr" and manifest["diagnostics"]["inference"] == "exact"
    with (tmp_path/"out"/"genes.tsv").open() as handle:rows = {r["ORF"]:r for r in csv.DictReader(handle, delimiter="\t")}
    assert rows["g0002"]["call"] == "depleted"
    with pytest.raises(ValueError, match="not supported"):run(tmp_path/"norm"/"quantile", gff, reference, tmp_path/"q", ["P::c1"], ["P::t1"])
    with pytest.raises(ValueError, match="role"):run(tmp_path/"norm", gff, reference, tmp_path/"r", ["P::c1", "P::t1"], ["P::c2", "P::t2"])


def test_stratified_exact_matches_enumeration_within_libraries():
    rng = np.random.default_rng(8)
    control = rng.poisson(15, (4, 2))*np.array([1.1, .9]); treatment = rng.poisson(9, (4, 2))*np.array([.8, 1.3])
    rows = analyze(control, treatment, [dict(id="g", name="g", description="", chrom="c", strand="+")], [np.arange(4)], np.arange(4),
                   Settings(), control_libraries=["L1", "L2"], treatment_libraries=["L1", "L2"])[0]
    # Brute force: each (site, library) unit swaps its single control/treatment pair or not.
    pairs = [(control[i, L], treatment[i, L]) for i in range(4) for L in range(2)]
    t_obs = np.mean([t for _, t in pairs])-np.mean([c for c, _ in pairs]); hits = 0
    for flips in itertools.product((0, 1), repeat=len(pairs)):
        diff = np.mean([(c if f else t)-(t if f else c) for (c, t), f in zip(pairs, flips)]); hits += abs(diff) >= abs(t_obs)-1e-9
    assert rows[0]["delta_mean"] == pytest.approx(t_obs) and rows[0]["p_value"] == pytest.approx(hits/2**len(pairs), abs=2e-3)
    assert rows[0]["p_value_rounded_problem"] == pytest.approx(hits/2**len(pairs)) and rows[0]["libraries_with_sites"] == 2
    with pytest.raises(ValueError, match="library"):
        analyze(control, treatment, [dict(id="g", name="g", description="", chrom="c", strand="+")], [np.arange(4)], np.arange(4),
                Settings(), control_libraries=["L1", "L1"], treatment_libraries=["L1", "L2"])


def brute_force_one_sided(control, treatment, alternative):
    J, K = control.shape[1], treatment.shape[1]; values = np.concatenate([control, treatment], 1)
    t_obs = treatment.mean()-control.mean(); per_site = list(itertools.combinations(range(J+K), K)); hits = total = 0
    for choice in itertools.product(per_site, repeat=len(values)):
        mask = np.zeros(values.shape, dtype=bool)
        for i, subset in enumerate(choice):mask[i, list(subset)] = True
        t = values[mask].mean()-values[~mask].mean(); total += 1
        hits += (t <= t_obs+1e-9) if alternative == "depleted" else (t >= t_obs-1e-9)
    return hits/total


@pytest.mark.parametrize("alternative", ["depleted", "enriched"])
def test_one_sided_exact_matches_enumeration(alternative):
    from ytab.essentiality.resampling import exact_units
    rng = np.random.default_rng(5)
    control = rng.poisson(20, (5, 2))*1.37; treatment = rng.poisson(12, (5, 2))*.91
    sums, _ = subset_sums(control, treatment)
    upper, rounded, _ = exact_units([(sums, control.sum(1)+treatment.sum(1))], 10, 10, 65536, alternative)
    truth = brute_force_one_sided(control, treatment, alternative)
    assert rounded == pytest.approx(truth, abs=1e-12) and upper >= truth-1e-12


def test_winsorize_is_label_invariant_and_exact():
    rng = np.random.default_rng(6)
    control = rng.poisson(10, (4, 2))*1.37; treatment = rng.poisson(10, (4, 2))*.91; treatment[1, 0] = 500.   # non-lattice values
    genes = [dict(id="g", name="g", description="", chrom="c", strand="+")]
    rows, _ = analyze(control, treatment, genes, [np.arange(4)], np.arange(4), Settings(winsorize=True))
    values = np.sort(np.concatenate([control, treatment], 1).ravel()); second = values[-2]
    expected, _ = brute_force_p(np.minimum(control, second), np.minimum(treatment, second))
    assert rows[0]["p_value_rounded_problem"] == pytest.approx(expected, abs=1e-12) and rows[0]["p_value"] >= expected-1e-12
    swapped, _ = analyze(treatment, control, genes, [np.arange(4)], np.arange(4), Settings(winsorize=True))
    # same winsorized data either labelling (mirror-image problem; equal up to grid rounding, still bounded)
    assert swapped[0]["p_value_rounded_problem"] == pytest.approx(expected, abs=5e-3) and swapped[0]["p_value"] >= expected-1e-12


def test_recenter_removes_a_global_scale_shift():
    C, T, genes, idx, pos = experiment(21)
    rows, diag = analyze(C, T*1.8, genes, idx, pos, Settings(recenter="median-gene"))
    assert diag["treatment_rescale"] == pytest.approx(1/1.8, rel=.08)
    assert np.mean([r["p_value"] <= .05 for r in rows if r["in_family"]]) < .08


def test_noise_asymmetry_warning():
    from ytab.essentiality.resampling import noise_asymmetry
    rng = np.random.default_rng(8); mu = np.exp(rng.normal(3, 1, 5000))
    quiet = np.column_stack([rng.poisson(mu), rng.poisson(mu)]).astype(float)
    noisy = np.column_stack([rng.negative_binomial(2, 2/(2+mu)), rng.negative_binomial(2, 2/(2+mu))]).astype(float)
    assert noise_asymmetry(quiet, quiet[:, ::-1])["warning"] is None
    result = noise_asymmetry(quiet, noisy)
    assert result["asymmetry"] > 2 and "treatment replicates" in result["warning"]
    assert noise_asymmetry(quiet[:, :1], noisy)["asymmetry"] is None
