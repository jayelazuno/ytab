import numpy as np
import pandas as pd

from ytab.fitness import zscore


def _project(tmp_path, scale=1.0, seed=3):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(800)]
    mu = rng.lognormal(4, 1.2, len(genes))
    effect = np.ones(len(genes)); effect[:40] = 1 / 16; effect[40:80] = 16
    rows, files = [], {}
    for pool in ("1", "2"):
        for role, fold in (("parent", np.ones(len(genes))), ("treated", effect)):
            sample = f"bg-{role}-{pool}"
            reads = rng.negative_binomial(20, 20 / (20 + mu * fold))
            if role == "treated" and pool == "1": reads = np.round(reads * scale)
            path = tmp_path / f"{sample}.csv"
            pd.DataFrame({"standard_name": genes, "reads": reads}).to_csv(path, index=False)
            files[sample] = path
        rows.append(dict(comparison_id=f"p{pool}", pair_id=f"p{pool}", parent_sample=f"bg-parent-{pool}",
                         treated_sample=f"bg-treated-{pool}", background="bg", pool=pool, parent_condition="parent",
                         treated_condition="treated", include="true", notes=""))
    design = tmp_path / "design.csv"; pd.DataFrame(rows).to_csv(design, index=False)
    return design, files, genes


def test_planted_effects_are_called(tmp_path):
    design, files, genes = _project(tmp_path)
    zscore.run(design, tmp_path / "out", sample_files=files)
    calls = pd.read_csv(tmp_path / "out" / "tables" / "treated_vs_parent.summary_by_feature.csv").set_index("feature_id").final_call
    # Default: replicates averaged, one z per gene (Gale et al. 2020).
    assert (calls[genes[:40]] == "depleted").mean() > .8
    assert (calls[genes[40:80]] == "enriched").mean() > .8
    assert (calls[genes[80:]] != "none").mean() < .02


def test_pairs_rule_and_recentering(tmp_path):
    design, files, genes = _project(tmp_path)
    zscore.run(design, tmp_path / "out", sample_files=files, settings=zscore.Settings(combine="pairs", recenter="trimmed-sum"))
    calls = pd.read_csv(tmp_path / "out" / "tables" / "treated_vs_parent.summary_by_feature.csv").set_index("feature_id").final_call
    assert (calls[genes[:40]] == "consistently_depleted").mean() > .8
    assert (calls[genes[40:80]] == "consistently_enriched").mean() > .8
    # Null genes: no consistent calls (single-pool calls can follow from CPM composition shifts caused by the planted effects).
    assert calls[genes[80:]].str.startswith("consistently").mean() < .01
    combined = pd.read_csv(tmp_path / "out" / "tables" / "treated_vs_parent.combined.log2fc_z.csv")
    # The planted 16x enrichment shrinks every other gene's CPM; recentering puts unaffected genes back at log2FC ~ 0.
    unaffected = combined.set_index("feature_id").loc[genes[80:]]
    assert combined.recenter_factor.nunique() == 1 and combined.recenter_factor.iloc[0] > 1.2
    assert abs(unaffected[unaffected.A_exp > 50].log2FC.median()) < .1


def test_cpm_makes_calls_depth_invariant(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"; a.mkdir(); b.mkdir()
    design, files, _ = _project(a)
    zscore.run(design, a / "out", sample_files=files)
    design2, files2, _ = _project(b, scale=3.0)
    zscore.run(design2, b / "out", sample_files=files2)
    one = pd.read_csv(a / "out" / "tables" / "treated_vs_parent.summary_by_feature.csv").set_index("feature_id").final_call
    two = pd.read_csv(b / "out" / "tables" / "treated_vs_parent.summary_by_feature.csv").set_index("feature_id").final_call
    assert (one == two.loc[one.index]).mean() > .97


def test_running_and_quantile_match_r_conventions():
    x = np.arange(1, 11, dtype=float)
    assert zscore.quantile7(x, .99) == np.quantile(x, .99)
    mean, sd = zscore.running(x, 1, np.mean), zscore.running(x, 1, lambda v: np.std(v, ddof=1))
    assert np.isnan(mean[0]) and mean[1] == 2 and np.isclose(sd[1], 1) and np.isnan(sd[-1])


def _normalized(files, genes, scale):
    """Gene x sample values standing in for a site normalization, on an arbitrary scale."""
    values = pd.DataFrame({s: pd.read_csv(p).reads.to_numpy()*scale for s, p in files.items()}, index=pd.Index(genes, name="feature_id"))
    return dict(genes=values, raw=values/scale, label="betageom", run="run",
                annotation=pd.DataFrame(dict(feature_id=genes, standard_name=genes, common_name=[f"N{g}" for g in genes])))


def test_normalized_values_are_per_million_so_scale_does_not_matter(tmp_path):
    from ytab.fitness import contract
    design, files, genes = _project(tmp_path)
    for name, scale in (("one", 1.), ("big", 37.)):
        zscore.run(design, tmp_path/name, normalized=_normalized(files, genes, scale))
        assert contract.validate(tmp_path/name) == []
    one, big = (pd.read_csv(tmp_path/n/"tables/treated_vs_parent.summary_by_feature.csv").set_index("feature_id") for n in ("one", "big"))
    assert (one.final_call == big.final_call.loc[one.index]).all()
    assert np.allclose(one.combined_z, big.combined_z.loc[one.index], rtol=0, atol=.05)   # noise fit converges within ~0.4%
    manifest = __import__("json").loads((tmp_path/"big/zscore_manifest.json").read_text())
    assert manifest["value_unit"] == "BGC reads per million" and manifest["library_sizes"]
    values = pd.read_csv(tmp_path/"big/tables/gene_values.long.csv")
    assert np.isclose(values.groupby("sample").value.sum().mean(), 1e6)
    assert one.common_name.loc[genes[0]] == f"N{genes[0]}"


def test_cpm_result_meets_contract(tmp_path):
    from ytab.fitness import contract
    design, files, _ = _project(tmp_path)
    zscore.run(design, tmp_path/"out", sample_files=files)
    assert contract.validate(tmp_path/"out") == []
    summary = pd.read_csv(tmp_path/"out/tables/treated_vs_parent.summary_by_feature.csv")
    assert (summary.call_log2FC == summary.combined_log2FC).all()   # the replicate rule's own effect


def test_top_hit_ranking_is_written_by_the_pipeline(tmp_path):
    design, files, genes = _project(tmp_path)
    zscore.run(design, tmp_path/"out", sample_files=files)
    hits = pd.read_csv(tmp_path/"out/tables/top_hit_evidence.csv").set_index("feature_id")
    top = hits.rank_depleted.dropna().sort_values()
    assert set(top.index[:30]) <= set(genes[:40])                 # planted 16x depletions lead the depleted ranks
    assert hits.loc[genes[40:80]].rank_enriched.notna().mean() > .8
    # candidates (Gale et al. 2020): called genes with at least a 2-fold change; ranks 1..n ordered by |z|
    cand = hits[hits.rank_both.notna()].sort_values("rank_both")
    assert (cand.log2FC.abs() >= 1).all() and cand["call"].isin(["depleted", "enriched"]).all()
    assert list(cand.rank_both) == list(range(1, len(cand)+1)) and cand.z_strength_both.is_monotonic_decreasing
    assert hits[hits["call"] == "none"].rank_both.isna().all()
    assert (hits.n_pairs_depleted_support <= hits.n_pairs_valid).all() and hits.n_pairs_valid.max() == 2
    pairs = pd.read_csv(tmp_path/"out/tables/top_hit_evidence.by_pair.csv")
    assert set(pairs.contrast) == {"p1", "p2"} and (pairs.n_pairs_valid <= 1).all()
