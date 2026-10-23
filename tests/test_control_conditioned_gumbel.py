import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.stats import gumbel_r

from ytab.essentiality.control_conditioned_gumbel import (
    Settings, analyze, benjamini_hochberg, combined_dispersion, depletion_thresholds, discrete_gumbel_cdf, discrete_gumbel_logpmf,
    discrete_gumbel_sf, exact_run_tails, fit_background, fit_discrete_gumbel, gene_sites, gumbel_log_tail_se, longest_run,
    longest_runs, predictive_cdf, run, selection_log_probability, simulate_runs, site_posterior,
    synthetic_experiment)
from ytab.normalization.methods import run_choices
from ytab.normalization.ttr import sha256

FAST = Settings(calibration_sims=300, validation_sims=300)
EXPOSURES = np.array([1., 1.3, .8, 1.1])   # two controls then two treatment replicates
CONTROL_EXPOSURES = EXPOSURES[:2]


def experiment(seed, sites=40, genes=120, effects=None, exposures=EXPOSURES, **model):
    rng = np.random.default_rng(seed)
    return synthetic_experiment(rng, np.full(genes, sites), exposures, effects=effects, treatments=2, **model)


def block_effects(genes, sites, affected, fold, start=5, length=25, scattered=False):
    effects = np.ones(genes*sites)
    for g in affected:
        offsets = np.arange(start, start+2*length, 2)[:length] if scattered else np.arange(start, start+length)
        effects[g*sites+offsets] = fold
    return effects


def test_settings_and_input_guards():
    for bad in (dict(alpha_site=.6), dict(min_control_reads=0), dict(calibration_sims=10),
                dict(eligibility="treated"), dict(nterm=60, cterm=50), dict(tail_tolerance=1.), dict(inference="gumbel")):
        with pytest.raises(ValueError):Settings(**bad).validate()
    C, T, genes, idx, pos = experiment(1, genes=10)
    with pytest.raises(ValueError, match="integer"):analyze(C-1, T, EXPOSURES, genes, idx, pos, FAST)
    with pytest.raises(ValueError, match="integer"):analyze(C+.5, T, EXPOSURES, genes, idx, pos, FAST)
    with pytest.raises(ValueError, match=">=2 control"):analyze(C[:, :1], T, EXPOSURES[1:], genes, idx, pos, FAST)
    with pytest.raises(ValueError, match="exposure"):analyze(C, T, [1, 0, 1, 1], genes, idx, pos, FAST)
    with pytest.raises(ValueError, match=">=2 treatment"):analyze(C, T[:, 0], EXPOSURES[:3], genes, idx, pos, FAST)
    assert analyze(C, T[:, 0], EXPOSURES[:3], genes, idx, pos, Settings(calibration_sims=300, validation_sims=300, treatment_dispersion="control"))[2]["dispersion"]["source"] == "control"
    with pytest.raises(ValueError, match="unique"):analyze(C, T, EXPOSURES, genes+genes[:1], idx+idx[:1], pos, FAST)
    with pytest.raises(ValueError, match="20.*supported"):
        analyze(np.zeros_like(C), T, EXPOSURES, genes, idx, pos, FAST)
    assert Settings().inference == "exact" and Settings().calibration_sims == Settings().validation_sims == 1000
    assert Settings(inference="gumbel-mc").calibration_sims == 20000 and Settings(inference="gumbel-mc", validation_sims=500).validation_sims == 500


def test_background_fit_recovers_parameters_and_selection_matters():
    rng = np.random.default_rng(3)
    C, *_ = synthetic_experiment(rng, [6000], EXPOSURES[:3], m=np.log(4), tau=1.6, phi=.3)
    fit = fit_background(C[C.sum(1) >= 1], CONTROL_EXPOSURES, FAST)
    assert fit["m"] == pytest.approx(np.log(4), abs=.15)
    assert fit["tau"] == pytest.approx(1.6, abs=.15) and fit["phi"] == pytest.approx(.3, abs=.08)
    assert fit["integration_check"]["max_abs_log_marginal_difference"] < 1e-6
    assert set(fit["boundary"]) <= {"linear"}
    # The extra linear component need not estimate exactly zero in a finite
    # sample. Check the recovered variance curve, not an arbitrary coefficient
    # cutoff; the separate linear-noise test checks identifiability away from zero.
    for mean in (1., 4., 30.):
        recovered_variance = mean*(1+fit["linear"])+fit["phi"]*mean**2
        assert recovered_variance == pytest.approx(mean+.3*mean**2, rel=.1)
    # Ignoring the zero-truncation (treating observed-only sites as unselected) biases baseline abundance upward.
    naive = np.mean(np.log(C[C.sum(1) >= 1].sum(1)/CONTROL_EXPOSURES.sum()))
    assert naive > fit["m"]+.1


def test_selection_probability_matches_monte_carlo():
    rng = np.random.default_rng(9); m, tau, phi = np.log(2), 1.2, .5
    lam = np.exp(m+tau*rng.standard_normal(400000))
    counts = np.column_stack([rng.negative_binomial(1/phi, (1/phi)/(1/phi+e*lam)) for e in CONTROL_EXPOSURES])
    for minimum in (1, 3):
        expected = np.mean(counts.sum(1) >= minimum)
        assert np.exp(selection_log_probability(m, tau, phi, CONTROL_EXPOSURES, minimum)) == pytest.approx(expected, abs=3e-3)


def test_predictive_flags_low_control_sites_and_partial_depletion():
    profiles = np.array([[1, 0], [0, 1], [200, 260], [5, 7]])
    u, w, _ = site_posterior(profiles, CONTROL_EXPOSURES, np.log(30), 1.5, .2, 32)
    threshold, probability = depletion_thresholds(u, w, EXPOSURES[2], .2, .05)
    assert threshold[0] == threshold[1] == -1 and probability[0] == 0   # weak controls give uncertainty, not depletion
    assert threshold[2] > 0 and 0 < probability[2] <= .05                 # reads retained can still be depleted
    t = np.arange(threshold[2]+3)
    cdf = predictive_cdf(t, np.repeat(u[2:3], len(t), 0), np.repeat(w[2:3], len(t), 0), EXPOSURES[2], .2)
    assert np.all(np.diff(cdf) >= 0) and cdf[threshold[2]] <= .05 < cdf[threshold[2]+1]


def test_exact_run_tails_match_enumeration():
    p = np.array([.3, .6, .1, .8, .5, .2, .7])
    tails = exact_run_tails(p, 7)
    expected = np.zeros(8)
    for bits in itertools.product((0, 1), repeat=len(p)):
        weight = np.prod(np.where(bits, p, 1-p)); expected[longest_run(bits)[0]] += weight
    assert np.allclose(tails, expected[::-1].cumsum()[::-1][1:])


def test_longest_run_tie_convention_and_vector_form():
    assert longest_run([1, 1, 0, 1, 1, 0, 1]) == (2, 3, 4)
    assert longest_run([0, 0]) == (0, None, None)
    flags = np.array([[1, 1, 0, 1, 1, 1], [0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 1, 1]], dtype=bool)
    assert longest_runs(flags).tolist() == [3, 0, 2]


def test_count_simulation_matches_exact_bernoulli_null():
    profiles = np.array([[40, 50], [3, 2], [120, 90], [12, 20], [1, 0], [60, 60]]*5)
    u, w, _ = site_posterior(profiles, CONTROL_EXPOSURES, np.log(30), 1.5, .2, 32)
    threshold, probability = depletion_thresholds(u, w, EXPOSURES[2], .2, .2)
    runs = simulate_runs(np.random.default_rng(5), u, w, threshold, EXPOSURES[2], .2, 40000)
    tails = exact_run_tails(probability, 6)
    assert np.allclose([(runs >= k).mean() for k in range(1, 7)], tails, atol=4e-3)


def test_discrete_gumbel_fit_tail_and_failures():
    draws = np.clip(np.round(gumbel_r.rvs(3, 1.2, size=8000, random_state=3)), 0, 40)
    fit = fit_discrete_gumbel(draws, 40)
    assert fit["converged"] and fit["mu"] == pytest.approx(3, abs=.1) and fit["beta"] == pytest.approx(1.2, abs=.08)
    pmf = np.exp(discrete_gumbel_logpmf(np.arange(41), 40, 3, 1.2))
    assert pmf.sum() == pytest.approx(1) and discrete_gumbel_cdf(40, 40, 3, 1.2) == pytest.approx(1)
    assert discrete_gumbel_sf(5, 40, 3, 1.2) == pytest.approx(pmf[5:].sum())
    assert discrete_gumbel_sf(0, 40, 3, 1.2) == 1 and discrete_gumbel_sf(41, 40, 3, 1.2) == 0
    degenerate = fit_discrete_gumbel(np.zeros(500), 10)
    assert not degenerate["converged"] and "degenerate" in degenerate["reason"]


def test_no_effect_false_positives_and_unequal_depths():
    for exposures in (EXPOSURES, np.array([1., 1.2, 3., 2.5]), np.array([1., .9, .3, .4])):
        C, T, genes, idx, pos = experiment(11, exposures=exposures)
        rows, sites, diagnostics = analyze(C, T, exposures, genes, idx, pos, FAST)
        family = [r for r in rows if r["in_family"]]
        assert diagnostics["calls"]["depleted run"] == 0
        assert np.mean([r["p_family"] <= .05 for r in family]) <= .1
        expected, observed = diagnostics["expected_depleted_sites_null"], diagnostics["observed_depleted_sites"]
        assert abs(observed-expected) < 5*np.sqrt(expected)


@pytest.mark.parametrize("fold", [0., .2])
def test_complete_and_partial_depletion_detected(fold):
    affected = range(0, 120, 12)
    C, T, genes, idx, pos = experiment(21, effects=block_effects(120, 40, affected, fold))
    rows, sites, _ = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)
    called = {r["ORF"] for r in rows if r["call"] == "depleted run"}
    truth = {genes[g]["id"] for g in affected}
    assert len(called & truth) >= 7 and len(called-truth) <= 1
    if fold:
        flagged = [s for s in sites if s["depleted"] and s["genes"] in truth]
        assert sum(s["treatment_combined"] > 0 for s in flagged) > len(flagged)/4   # depletion recognized with reads remaining


def test_adjacent_depletion_beats_isolated_and_zero_counts_reported():
    affected = range(0, 120, 12)
    C, T, genes, idx, pos = experiment(31, effects=block_effects(120, 40, affected, 0., length=14))
    adjacent = {r["ORF"]:r for r in analyze(C, T, EXPOSURES, genes, idx, pos, FAST)[0]}
    C2, T2, *_ = experiment(31, effects=block_effects(120, 40, affected, 0., length=14, scattered=True))
    isolated = {r["ORF"]:r for r in analyze(C2, T2, EXPOSURES, genes, idx, pos, FAST)[0]}
    truth = [genes[g]["id"] for g in affected]
    assert np.mean([adjacent[g]["r"] for g in truth]) > 2*np.mean([isolated[g]["r"] for g in truth])
    assert sum(adjacent[g]["call"] == "depleted run" for g in truth) > sum(isolated[g]["call"] == "depleted run" for g in truth)
    effects = np.ones(120*40); effects[:40] = 0
    C, T, genes, idx, pos = experiment(32, effects=effects)
    row = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)[0][0]
    assert row["treatment_reads_zero"] and row["log2_scaled_ratio"] is None and row["call"] == "depleted run"


def test_few_and_many_eligible_coordinates():
    rng = np.random.default_rng(41)
    C, T, genes, idx, pos = synthetic_experiment(rng, [0, 2, 3, 400]+[40]*60, EXPOSURES, treatments=2)
    rows = {r["ORF"]:r for r in analyze(C, T, EXPOSURES, genes, idx, pos, FAST)[0]}
    assert rows["g0000"]["status"] == "no control-supported sites" and rows["g0000"]["call"] == "not tested"
    assert rows["g0001"]["status"] == "too few control-supported sites" and not rows["g0001"]["in_family"]
    assert rows["g0002"]["in_family"] and rows["g0003"]["n_eligible"] > 350 and "gumbel_mu" in rows["g0003"]


def test_detectable_eligibility_drops_uninformative_sites():
    C, T, genes, idx, pos = experiment(45, m=np.log(3))
    base = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)
    strict = analyze(C, T, EXPOSURES, genes, idx, pos, Settings(calibration_sims=300, validation_sims=300, eligibility="detectable"))
    assert strict[2]["gene_assigned_eligible_sites"] < base[2]["gene_assigned_eligible_sites"]
    assert all(s["depletion_threshold"] >= 0 for s in strict[1])


def test_reproducible_order_independent_and_no_mutation():
    C, T, genes, idx, pos = experiment(51, genes=40, effects=block_effects(40, 40, [3], .1))
    copies = C.copy(), T.copy(); state = np.random.get_state()[1].copy()
    first = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)
    assert np.array_equal(C, copies[0]) and np.array_equal(T, copies[1]) and np.array_equal(state, np.random.get_state()[1])
    second = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)
    reversed_rows = analyze(C, T, EXPOSURES, genes[::-1], idx[::-1], pos, FAST)[0]
    assert first[0] == second[0]
    by_gene = {r["ORF"]:r for r in reversed_rows}
    for row in first[0]:
        assert by_gene[row["ORF"]]["gumbel_mu"] == row["gumbel_mu"] and by_gene[row["ORF"]]["q_bh"] == row["q_bh"]
    other = analyze(C, T, EXPOSURES, genes, idx, pos, Settings(calibration_sims=300, validation_sims=300, seed=7))[0]
    assert any(a.get("gumbel_mu") != b.get("gumbel_mu") for a, b in zip(first[0], other))


def test_benjamini_hochberg():
    assert np.allclose(benjamini_hochberg([.01, .04, .03, .5]), [.04, .16/3, .16/3, .5])


def test_gene_sites_follow_transit_trimming_conventions():
    blocks = [(1, 300), (401, 800)]; length = 700   # two-exon CDS; spliced length 700
    positions = np.arange(1, 801, 7)
    sites = [("chr1", int(p)) for p in positions]
    def spliced(p):return p-1 if p <= 300 else 300+p-401 if 401 <= p <= 800 else None
    for strand in "+-":
        gene = [dict(id="g", name="g", description="", chrom="chr1", strand=strand, blocks=blocks)]
        for nterm, cterm, stop in ((0., 0., True), (10., 5., False), (0., 0., False)):
            chosen = gene_sites(gene, sites, Settings(nterm=nterm, cterm=cterm, ignore_stop=stop))[0]
            expected = []
            for k, p in enumerate(positions):
                offset = spliced(p)
                if offset is None:continue
                oriented = offset if strand == "+" else length-1-offset
                if length*nterm/100 <= oriented < length*(1-cterm/100) and (not stop or oriented < length-3):expected.append(k)
            assert sorted(chosen.tolist()) == expected and np.all(np.diff([sites[k][1] for k in chosen]) > 0)


def write_catalogue(root, counts, sites, metadata, reference):
    root.mkdir()
    sparse.save_npz(root/"raw_counts.npz", sparse.csc_matrix(counts))
    with gzip.open(root/"sites.tsv.gz", "wt") as handle:
        handle.write("site_index\tchromosome\tposition_1based\n")
        for i, (chrom, p) in enumerate(sites):handle.write(f"{i}\t{chrom}\t{p}\n")
    fields = ["observation_id", "display_name", "project", "library_id", "pool", "role"]
    with (root/"observations.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(metadata)
    (root/"run_manifest.json").write_text(json.dumps(dict(status="success", reference_sha256=sha256(reference), universe_sites=len(sites))))


def test_file_runner_on_normalization_outputs_and_guards(tmp_path):
    C, T, genes, idx, pos = experiment(61, genes=40, effects=block_effects(40, 40, [2], 0.))
    sites = [("chr1", int(p)) for p in pos]
    gff = tmp_path/"genes.gff"; reference = tmp_path/"ref.fna"; reference.write_text(">chr1\nACGT\n")
    with gff.open("w") as handle:
        handle.write("##gff-version 3\n")
        for g in genes:
            a, b = g["blocks"][0]
            handle.write(f"chr1\t.\tgene\t{a}\t{b}\t.\t+\t.\tID={g['id']};locus_tag={g['id']};Name={g['id']}\n")
            handle.write(f"chr1\t.\tCDS\t{a}\t{b}\t.\t+\t0\tID=cds-{g['id']};Parent={g['id']}\n")
    meta = [dict(observation_id=f"P::{x}", display_name=x, project="P", library_id="P::L1", pool="L1", role=r)
            for x, r in (("c1", "control"), ("c2", "control"), ("t1", "treated"), ("t2", "treated"), ("o1", "treated"))]
    counts = np.column_stack([C, T, T[:, 0]])
    catalogue = tmp_path/"cat"; write_catalogue(catalogue, counts, sites, meta, reference)
    normalized = tmp_path/"norm"; run_choices(catalogue, normalized, ["ttr", "totreads", "nzmean", "quantile"], "library_id")
    controls, treated = ["P::c1", "P::c2"], ["P::t1", "P::t2"]
    out = tmp_path/"out"
    manifest = run(normalized, gff, reference, out, controls, treated, FAST)   # run directory: ttr by default
    assert manifest["status"] == "success" and manifest["normalization"]["method"] == "ttr"
    with (out/"genes.tsv").open() as handle:rows = {r["ORF"]:r for r in csv.DictReader(handle, delimiter="\t")}
    assert rows["g0002"]["call"] == "depleted run" and (out/"sites.tsv.gz").exists()
    with (normalized/"ttr"/"factors.tsv").open() as handle:factors = {r["observation_id"]:float(r["factor"]) for r in csv.DictReader(handle, delimiter="\t")}
    inverse = {k:1/v for k, v in factors.items()}
    assert manifest["exposures"]["values"]["P::t1"] == pytest.approx(inverse["P::t1"]/np.mean([inverse[c] for c in controls]))
    assert manifest["diagnostics"]["dispersion"]["source"] == "replicates" and manifest["diagnostics"]["inference"] == "exact"
    # totreads exposures are each observation's total support relative to the control mean.
    totreads = run(normalized/"totreads", gff, reference, tmp_path/"tot", controls, treated, FAST)
    assert totreads["exposures"]["values"]["P::t2"] == pytest.approx(T[:, 1].sum()/C.sum(0).mean())
    assert run(normalized, gff, reference, tmp_path/"nz", controls, treated, FAST, "nzmean")["normalization"]["method"] == "nzmean"
    with pytest.raises(ValueError, match="not supported"):run(normalized/"quantile", gff, reference, tmp_path/"q", controls, treated, FAST)
    with pytest.raises(ValueError, match="not the requested"):run(normalized/"ttr", gff, reference, tmp_path/"m", controls, treated, FAST, "nzmean")
    with pytest.raises(ValueError, match="new output"):run(normalized, gff, reference, out, controls, treated, FAST)
    with pytest.raises(ValueError, match="distinct treatment"):run(normalized, gff, reference, tmp_path/"o2", controls, ["P::c1"], FAST)
    with pytest.raises(ValueError, match=">=2"):run(normalized, gff, reference, tmp_path/"o3", ["P::c1"], treated, FAST)
    with pytest.raises(ValueError, match="role"):run(normalized, gff, reference, tmp_path/"o4", ["P::c1", "P::t1"], ["P::t2", "P::o1"], FAST)
    with pytest.raises(ValueError, match="checksum"):run(normalized, gff, gff, tmp_path/"o5", controls, treated, FAST)
    with pytest.raises(ValueError, match="normalization output"):run(catalogue, gff, reference, tmp_path/"o6", controls, treated, FAST)
    # A normalized matrix that is not raw*factor (or a catalogue changed after normalization) is refused.
    tampered = tmp_path/"tampered"; run_choices(catalogue, tampered, ["ttr"], "library_id")
    matrix = sparse.load_npz(tampered/"ttr"/"normalized_counts.npz").tolil(); matrix[0, 0] += 1; sparse.save_npz(tampered/"ttr"/"normalized_counts.npz", matrix.tocsr())
    with pytest.raises(ValueError, match="not raw catalogue counts"):run(tampered, gff, reference, tmp_path/"o7", controls, treated, FAST)
    mixed_meta = [dict(m) for m in meta]; mixed_meta[4]["library_id"] = "P::L2"
    mixed = tmp_path/"mixed"; write_catalogue(mixed, counts, sites, mixed_meta, reference)
    run_choices(mixed, tmp_path/"mixed_norm", ["ttr"], "library_id")
    with pytest.raises(ValueError, match="library_id"):run(tmp_path/"mixed_norm", gff, reference, tmp_path/"o8", controls, ["P::t1", "P::o1"], FAST)
    sparse.save_npz(catalogue/"raw_counts.npz", sparse.csc_matrix(counts+1))
    with pytest.raises(ValueError, match="Catalogue changed"):run(normalized, gff, reference, tmp_path/"o9", controls, treated, FAST)


def test_combined_dispersion_matches_sum_of_replicates():
    rng = np.random.default_rng(71); e = np.array([1., 1.6]); phi = .3; lam = 20.
    total = sum(rng.negative_binomial(1/phi, (1/phi)/(1/phi+x*lam), 400000) for x in e)
    combined = combined_dispersion(phi, e)
    assert total.mean() == pytest.approx(e.sum()*lam, rel=.01)
    assert total.var() == pytest.approx(e.sum()*lam+combined*(e.sum()*lam)**2, rel=.02)


def test_noisy_treatment_replicates_need_treatment_dispersion():
    """Gale-like null: treatment replicates far noisier than controls, no true effect."""
    C, T, genes, idx, pos = experiment(81, phi=.02, treatment_phi=.4)
    calibrated = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)[2]
    assumed = analyze(C, T, EXPOSURES, genes, idx, pos, Settings(calibration_sims=300, validation_sims=300, treatment_dispersion="control"))[2]
    assert calibrated["dispersion"]["treatment_replicate"]["phi"] > 8*calibrated["dispersion"]["control_replicate"]["phi"]
    assert calibrated["calls"]["depleted run"] == 0
    assert abs(calibrated["observed_depleted_sites"]-calibrated["expected_depleted_sites_null"]) < 5*np.sqrt(calibrated["expected_depleted_sites_null"])
    assert assumed["observed_depleted_sites"] > 3*assumed["expected_depleted_sites_null"] and assumed["calls"]["depleted run"] > 0


def test_linear_dispersion_recovered_and_needed_for_bottleneck_noise():
    rng = np.random.default_rng(91)
    C, *_ = synthetic_experiment(rng, [6000], EXPOSURES[:3], m=np.log(20), tau=1.4, phi=.05, linear=2.)
    fit = fit_background(C[C.sum(1) >= 1], CONTROL_EXPOSURES, FAST)
    assert fit["linear"] == pytest.approx(2., rel=.2) and fit["phi"] == pytest.approx(.05, abs=.03)
    assert fit["integration_check"]["max_abs_log_marginal_difference"] < 1e-5
    # Bottleneck-like treatment noise (variance proportional to the mean) with no true effect.
    C, T, genes, idx, pos = experiment(92, phi=.02, treatment_phi=.05, treatment_linear=4.)
    full = analyze(C, T, EXPOSURES, genes, idx, pos, FAST)[2]
    quadratic = analyze(C, T, EXPOSURES, genes, idx, pos, Settings(calibration_sims=300, validation_sims=300, dispersion_model="quadratic"))[2]
    assert full["calls"]["depleted run"] == 0
    assert abs(full["observed_depleted_sites"]-full["expected_depleted_sites_null"]) < 5*np.sqrt(full["expected_depleted_sites_null"])
    assert quadratic["observed_depleted_sites"] > full["observed_depleted_sites"] and quadratic["dispersion"]["treatment_replicate"]["linear"] == 0


def test_combined_dispersion_keeps_linear_term():
    rng = np.random.default_rng(93); e = np.array([1., 1.6]); dispersion = (.2, 3.); lam = 15.
    def draw(mean):
        size = 1/(dispersion[0]+dispersion[1]/mean); return rng.negative_binomial(size, size/(size+mean), 400000)
    total = draw(e[0]*lam)+draw(e[1]*lam); phi, linear = combined_dispersion(dispersion, e); mu = e.sum()*lam
    assert linear == 3. and total.var() == pytest.approx(mu*(1+linear)+phi*mu**2, rel=.02)


def test_exact_inference_is_seed_independent_and_detects_depletion():
    C, T, genes, idx, pos = experiment(61, genes=60, effects=block_effects(60, 40, [2, 9, 30], .1))
    runs = [analyze(C, T, EXPOSURES, genes, idx, pos, Settings(calibration_sims=300, validation_sims=300, seed=s, inference="exact"))
            for s in (1, 2)]
    (a, _, _), (b, _, _) = runs
    assert [r.get("p_family") for r in a] == [r.get("p_family") for r in b]
    assert [r["call"] for r in a] == [r["call"] for r in b]
    rows, _, diagnostics = runs[0]
    assert diagnostics["inference"] == "exact" and diagnostics["family_null_fit_inadequate"] == 0
    called = {r["ORF"] for r in rows if r["call"] == "depleted run"}
    assert {genes[g]["id"] for g in (2, 9, 30)} <= called
    for row in rows:
        if row["status"] == "exact conditional tail":
            assert row["p_value"] == row["exact_tail_p"] and row["gumbel_check"]


def test_monte_carlo_tail_se_shrinks_and_matches_refits():
    rng = np.random.default_rng(3); n = 30
    probabilities = np.full(n, .1)
    def fits(sims, k):
        return [fit_discrete_gumbel(longest_runs(rng.random((sims, n)) < probabilities), n) for _ in range(k)]
    small, large = fits(1000, 40), fits(16000, 3)
    se_small = np.median([gumbel_log_tail_se(f, 5, n) for f in small])
    se_large = np.median([gumbel_log_tail_se(f, 5, n) for f in large])
    assert 3 < se_small/se_large < 5.5   # ~sqrt(16)
    spread = np.std([np.log(discrete_gumbel_sf(5, n, f["mu"], f["beta"])) for f in small])
    assert .6 < se_small/spread < 1.6    # delta-method SE agrees with refit-to-refit variability


def test_monte_carlo_aware_gate_is_reproducible():
    C, T, genes, idx, pos = experiment(71, genes=60, effects=block_effects(60, 40, [4, 20], .1))
    calls = []
    for seed in (1, 2):
        rows, _, diagnostics = analyze(C, T, EXPOSURES, genes, idx, pos,
                                       Settings(calibration_sims=4000, validation_sims=4000, seed=seed, inference="gumbel-mc"))
        assert diagnostics["inference"] == "gumbel-mc"
        assert all(r.get("gumbel_log_tail_se") is None or r["gumbel_log_tail_se"] > 0 for r in rows)
        calls.append({r["ORF"] for r in rows if r["call"] == "depleted run"})
    assert calls[0] == calls[1] and {genes[4]["id"], genes[20]["id"]} <= calls[0]
    with pytest.raises(ValueError, match="inference"):
        Settings(inference="bootstrap").validate()
