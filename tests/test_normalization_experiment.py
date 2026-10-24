import csv
import gzip
import json

import numpy as np
import pytest
from scipy import sparse
import yaml

from ytab.normalization.experiment import evaluate
from ytab.normalization.ttr import sha256


def fixture(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 50 + "\n")
    project = tmp_path / "study"
    config = project / "config"
    config.mkdir(parents=True)
    sheet = config / "samples.csv"
    sheet.write_text("sample,include,pool,background,display_name\nc1,true,p1,bg,Control 1\nc2,true,p1,bg,Control 2\nt1,true,p1,bg,Treated 1\nt2,true,p1,bg,Treated 2\noriginal_parent,true,p1,bg,Original parent\n")
    (config / "comparison_design.csv").write_text("comparison_id,parent_sample,treated_sample,include\ncontrast1,c1,t1,true\ncontrast2,c2,t2,true\n")
    cfg = config / "project.yaml"
    cfg.write_text(yaml.safe_dump(dict(repo_root=str(tmp_path), project_id="study", output_project_dir=str(project), sample_sheet=str(sheet), reference=dict(fasta=str(reference)))))
    values = {"c1": [(1, 10), (2, 20)], "c2": [(1, 20), (2, 40)], "t1": [(1, 5), (3, 5)], "t2": [(1, 10), (3, 10)]}
    for sample, counts in values.items():
        hit = project / "create_hit_file" / sample / (sample + "_hits.txt")
        hit.parent.mkdir(parents=True)
        hit.write_text("Chromosome\tHit position\tHit count\n" + "".join(f"chr1\t{p}\t{n}\n" for p,n in counts))
    catalogue = tmp_path / "catalogue"
    catalogue.mkdir()
    with gzip.open(catalogue / "sites.tsv.gz", "wt") as handle:
        handle.write("chromosome\tposition_1based\nchr1\t1\nchr1\t2\nchr1\t4\n")
    (catalogue / "run_manifest.json").write_text(json.dumps(dict(status="success", reference_sha256=sha256(reference), universe_sites=3)))
    return cfg, catalogue, sheet


def test_experimental_counts_design_controls_and_no_double_scaling(tmp_path):
    cfg, catalogue, _ = fixture(tmp_path)
    original = sha256(catalogue / "sites.tsv.gz")
    output = tmp_path / "results"
    result = evaluate(cfg, catalogue, output, ["nonorm", "ttr"], plots=False)
    assert result["parent_catalogue_sites"] == 3 and result["experiment_additional_sites"] == 1
    raw = sparse.load_npz(output / "matrix/raw_counts.npz").toarray()
    np.testing.assert_array_equal(raw, [[10,20,5,10],[20,40,0,0],[0,0,0,0],[0,0,5,10]])
    np.testing.assert_array_equal(sparse.load_npz(output / "methods/nonorm/normalized_counts.npz").toarray(), raw)
    normalized = sparse.load_npz(output / "methods/ttr/normalized_counts.npz").toarray()
    assert np.array_equal(raw == 0, normalized == 0)
    np.testing.assert_allclose(normalized[:,0], normalized[:,1])
    with (output / "comparison_diagnostics.tsv").open() as handle:
        metrics = list(csv.DictReader(handle, delimiter="\t"))
    controls = [r for r in metrics if r["kind"] == "control_control"]
    assert len(metrics) == 6 and len(controls) == 2
    assert float(controls[0]["median_log2_ratio_shared"]) == 1
    assert float(controls[1]["median_log2_ratio_shared"]) == pytest.approx(0)
    assert controls[0]["spearman_observed_union"] == controls[1]["spearman_observed_union"]
    assert sha256(catalogue / "sites.tsv.gz") == original
    with pytest.raises(ValueError, match="nonempty"):
        evaluate(cfg, catalogue, output, ["ttr"], plots=False)


def test_reference_mismatch_and_cross_pool_pairs_fail(tmp_path):
    cfg, catalogue, sheet = fixture(tmp_path)
    sheet.write_text(sheet.read_text().replace("t1,true,p1", "t1,true,p2"))
    with pytest.raises(ValueError, match="share background and pool"):
        evaluate(cfg, catalogue, tmp_path / "bad", ["ttr"], plots=False)
    assert not (tmp_path / "bad").exists()
    sheet.write_text(sheet.read_text().replace("t1,true,p2", "t1,true,p1"))
    manifest = catalogue / "run_manifest.json"
    data = json.loads(manifest.read_text()); data["reference_sha256"] = "wrong"
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="same reference"):
        evaluate(cfg, catalogue, tmp_path / "bad", ["ttr"], plots=False)


def test_nonlinear_methods_get_matched_diagnostics_and_scalar_results_are_unchanged(tmp_path):
    """Quantile, BGC and aBGC have no single factor; their diagnostics compare normalized values directly.
    For scalar methods this equals the raw ratio times the factor ratio (the previous formula)."""
    cfg, catalogue, _ = fixture(tmp_path)
    output = tmp_path / "results"
    evaluate(cfg, catalogue, output, ["nonorm", "ttr", "quantile", "betageom", "aBGC"], plots=True)
    with (output / "comparison_diagnostics.tsv").open() as handle:
        metrics = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["method"] for r in metrics} == {"nonorm", "ttr", "quantile", "betageom", "aBGC"}
    assert all(r["median_log2_ratio_shared"] != "" for r in metrics)
    ttr = {r["comparison_id"]: float(r["median_log2_ratio_shared"]) for r in metrics if r["method"] == "ttr"}
    raw = sparse.load_npz(output / "matrix/raw_counts.npz").toarray()
    with (output / "methods/ttr/factors.tsv").open() as handle:
        factors = [float(r["factor"]) for r in csv.DictReader(handle, delimiter="\t")]
    shared = (raw[:, 0] > 0) & (raw[:, 1] > 0)
    previous = np.median(np.log2(raw[shared, 1] / raw[shared, 0]) + np.log2(factors[1] / factors[0]))
    assert ttr["c1_vs_c2"] == pytest.approx(previous)
    assert (output / "comparison_diagnostics.png").stat().st_size > 0


def test_explicit_sample_list_needs_no_pairs(tmp_path):
    """Strain comparisons (multi-condition): the listed samples, roles from the sample sheet, pools need not match."""
    cfg, catalogue, sheet = fixture(tmp_path)
    sheet.write_text("sample,include,pool,background,library_role,display_name\nc1,true,p1,bg,parent,Control 1\nc2,true,p2,bg,parent,Control 2\n"
                     "t1,true,m1,mutA,treated,Mutant A\nt2,true,m2,mutB,treated,Mutant B\noriginal_parent,true,p1,bg,parent,Original parent\n")
    output = tmp_path / "strains"
    result = evaluate(cfg, catalogue, output, ["nonorm", "ttr"], plots=False, sample_list=["c1", "c2", "t1", "t2"])
    assert result["sample_list"] == ["c1", "c2", "t1", "t2"] and result["design_sha256"] is None
    with (output / "matrix/observations.tsv").open() as handle:
        roles = {r["sample"]:r["role"] for r in csv.DictReader(handle, delimiter="\t")}
    assert roles == {"c1":"control", "c2":"control", "t1":"treated", "t2":"treated"}
    assert not (output / "comparisons.tsv").exists()
    with pytest.raises(ValueError, match="not in the sample sheet"):
        evaluate(cfg, catalogue, tmp_path / "bad", ["ttr"], plots=False, sample_list=["c1", "nope"])


def test_normalize_together_uses_one_group(tmp_path):
    cfg, catalogue, sheet = fixture(tmp_path)
    sheet.write_text("sample,include,pool,background,library_role\nc1,true,p1,bg,parent\nc2,true,p2,bg,parent\nt1,true,m1,mutA,treated\n"
                     "t2,true,m2,mutB,treated\noriginal_parent,true,p1,bg,parent\n")
    output = tmp_path / "together"
    result = evaluate(cfg, catalogue, output, ["ttr"], plots=False, sample_list=["c1", "c2", "t1", "t2"], together=True)
    assert result["normalization_group"] == "project"
    with (output / "methods/ttr/factors.tsv").open() as handle:
        groups = {r["normalization_group"] for r in csv.DictReader(handle, delimiter="\t")}
    assert len(groups) == 1
    with pytest.raises(ValueError, match="explicit sample list"):
        evaluate(cfg, catalogue, tmp_path / "bad", ["ttr"], plots=False, together=True)
