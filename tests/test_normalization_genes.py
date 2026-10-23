import csv
import gzip
import json
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse
import yaml

from ytab.normalization import gene_validation as gv
from ytab.normalization.experiment import evaluate as evaluate_sites
from ytab.normalization.ttr import sha256


def feature(gid, start, stop, exons):
    return SimpleNamespace(standard_name=gid, chromosome="chr1", start=start, stop=stop,
        exons=exons, coding_length=sum(b-a+1 for a,b in exons))


def test_projection_preserves_inclusive_boundaries_introns_and_existing_overlap():
    features = [feature("a", 2, 8, [(2,3),(7,8)]), feature("b", 5, 5, [(5,5)])]
    sites = [("chr1", p) for p in range(1,10)] + [("chr2", 5)]
    actual = gv.projection(features, sites).toarray()
    # SummaryTable counts the chromosome-wide exon mask inside each gene span.
    # Preserve that accepted behavior at an overlapping CDS inside an intron.
    np.testing.assert_array_equal(np.flatnonzero(actual[0]), [1,2,4,6,7])
    np.testing.assert_array_equal(np.flatnonzero(actual[1]), [4])


def test_independent_scaling_and_cpm_no_double_normalization():
    raw = np.array([[10,20],[30,40],[0,0]])
    result = gv.method_matrices(raw, [2,3])
    np.testing.assert_array_equal(result["nonorm"], raw)
    np.testing.assert_array_equal(result["ttr"], raw * [2,3])
    np.testing.assert_allclose(result["cpm"].sum(axis=0), [1e6,1e6])
    assert np.array_equal(result["ttr"] == 0, raw == 0)
    assert np.array_equal(result["cpm"] == 0, raw == 0)
    np.testing.assert_array_equal(raw, [[10,20],[30,40],[0,0]])
    with pytest.raises(ValueError): gv.method_matrices([[0,1]],[2,3])
    with pytest.raises(ValueError): gv.method_matrices([[1.5,2]],[2,3])


def test_independent_pool_comparisons_are_explicit_and_context_matched():
    left = dict(parent_sample="c1", background="a", comparison_group="zinc")
    right = dict(parent_sample="c2", background="a", comparison_group="zinc")
    libraries = {"c1":"pool1", "c2":"pool2"}
    assert gv.agreement_scope(left, right, libraries) is None
    assert gv.agreement_scope(left, right, libraries, True) == "between_independent_pools"
    assert gv.agreement_scope(left, dict(right, background="b"), libraries, True) is None
    assert gv.agreement_scope(left, dict(right, comparison_group="other"), libraries, True) is None
    assert gv.agreement_scope(left, dict(right, comparison_group=""), libraries, True) is None
    assert gv.agreement_scope(left, right, {"c1":"pool", "c2":"pool"}) == "within_library"


def fixture(tmp_path, monkeypatch, independent_pools=False):
    project = tmp_path / "study"
    config = project / "config"
    config.mkdir(parents=True)
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A"*100 + "\n")
    gff = tmp_path / "ref.gff"
    gff.write_text("##gff-version 3\n")
    sheet = config / "samples.csv"
    sheet.write_text("sample,include,pool,background,display_name\nc1,true,p1,bg,Control1\nc2,true,p1,bg,Control2\nt1,true,p1,bg,Treated1\nt2,true,p1,bg,Treated2\n")
    (config / "comparison_design.csv").write_text("comparison_id,parent_sample,treated_sample,include\ncontrast1,c1,t1,true\ncontrast2,c2,t2,true\n")
    if independent_pools:
        sheet.write_text(sheet.read_text().replace("c2,true,p1", "c2,true,p2").replace("t2,true,p1", "t2,true,p2"))
        (config / "comparison_design.csv").write_text("comparison_id,parent_sample,treated_sample,include,background,comparison_group\ncontrast1,c1,t1,true,bg,assay\ncontrast2,c2,t2,true,bg,assay\n")
    cfg = config / "project.yaml"
    cfg.write_text(yaml.safe_dump(dict(repo_root=str(tmp_path), project_id="study", species="glabrata", output_project_dir=str(project), sample_sheet=str(sheet), reference=dict(fasta=str(reference),gff=str(gff)))))
    values = {"c1":[10,20,0], "c2":[20,40,0], "t1":[5,5,0], "t2":[10,10,0]}
    for sample, counts in values.items():
        hit = project / "create_hit_file" / sample / (sample + "_hits.txt")
        hit.parent.mkdir(parents=True)
        hit.write_text("Chromosome\tHit position\tHit count\n" + "".join(f"chr1\t{i+1}\t{n}\n" for i,n in enumerate(counts)))
        summary = project / "summary" / sample / (sample + ".feature_table.RDF_1.csv")
        summary.parent.mkdir(parents=True)
        summary.write_text("RDF,1\nstandard_name,reads,common_name\n" + "".join(f"gene{i+1},{n},Name{i+1}\n" for i,n in enumerate(counts)))
    catalogue = tmp_path / "catalogue"
    catalogue.mkdir()
    with gzip.open(catalogue / "sites.tsv.gz", "wt") as handle:
        handle.write("chromosome\tposition_1based\nchr1\t1\nchr1\t2\nchr1\t3\n")
    (catalogue / "run_manifest.json").write_text(json.dumps(dict(status="success", reference_sha256=sha256(reference), universe_sites=3)))
    validation = tmp_path / "validation"
    evaluate_sites(cfg, catalogue, validation, ["nonorm","ttr"], plots=False)
    features = [feature(f"gene{i}",i,i,[(i,i)]) for i in range(1,4)]
    monkeypatch.setattr(gv, "load_feature_db", lambda *args: SimpleNamespace(get_all_features=lambda:features))
    return project, validation


def test_end_to_end_summary_parity_dropout_and_control_diagnostics(tmp_path, monkeypatch):
    project, validation = fixture(tmp_path, monkeypatch)
    output = tmp_path / "genes"
    result = gv.evaluate(validation, output, min_count=10, plots=False)
    assert result["raw_aggregation_exact"] and result["scaling_before_after_aggregation_agrees"]
    assert result["methods"] == ["nonorm","ttr","cpm"]
    metrics = gv.rows(output / "comparison_diagnostics.tsv")
    control = [r for r in metrics if r["kind"] == "control_control" and r["minimum_raw_count_each"] == "10"]
    assert float(control[0]["median_log2_ratio"]) == 1
    assert float(control[1]["median_log2_ratio"]) == pytest.approx(0)
    assert all(r["both_zero"] == "1" for r in metrics)
    with gzip.open(output / "gene_effects.tsv.gz", "rt") as handle:
        effects = list(csv.DictReader(handle, delimiter="\t"))
    assert all(r["log2_ratio"] == "" and r["evidence"] == "both_zero" for r in effects if r["feature_id"] == "gene3")
    with pytest.raises(ValueError, match="nonempty"):gv.evaluate(validation, output, plots=False)
    path = project / "summary/c1/c1.feature_table.RDF_1.csv"
    path.write_text(path.read_text().replace("gene1,10", "gene1,11"))
    with pytest.raises(ValueError, match="differs from accepted summary"):
        gv.evaluate(validation, tmp_path / "bad", plots=False)
    assert not (tmp_path / "bad").exists()


def test_modified_raw_input_rejected(tmp_path, monkeypatch):
    project, validation = fixture(tmp_path, monkeypatch)
    hit = project / "create_hit_file/c1/c1_hits.txt"
    hit.write_text(hit.read_text().replace("\t10", "\t11"))
    with pytest.raises(ValueError, match="Raw hit inputs changed"):
        gv.evaluate(validation, tmp_path / "bad", plots=False)


def test_between_pool_controls_are_not_technical_noise(tmp_path, monkeypatch):
    _, validation = fixture(tmp_path, monkeypatch, independent_pools=True)
    output = tmp_path / "genes"
    gv.evaluate(validation, output, plots=False, compare_independent_pools=True)
    metrics = gv.rows(output / "comparison_diagnostics.tsv")
    assert {r["kind"] for r in metrics} == {"treated_control", "between_pool_controls"}
    agreement = gv.rows(output / "treatment_agreement.tsv")
    assert {r["scope"] for r in agreement} == {"between_independent_pools"}
    gv.evaluate(validation, tmp_path / "strict", plots=False)
    assert {r["kind"] for r in gv.rows(tmp_path / "strict/comparison_diagnostics.tsv")} == {"treated_control"}
    assert not (tmp_path / "strict/treatment_agreement.tsv").exists()


def test_existing_cpm_and_default_effects_parity(tmp_path, monkeypatch):
    project, validation = fixture(tmp_path, monkeypatch)
    fitness = tmp_path / "fitness/tables"
    fitness.mkdir(parents=True)
    samples = [r["sample"] for r in gv.rows(validation / "matrix/observations.tsv")]
    raw = np.array([[10,20,5,10],[20,40,5,10],[0,0,0,0]])
    cpm = gv.method_matrices(raw, np.ones(4))["cpm"]
    with (fitness / "feature_reads_cpm.long.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["feature_id","sample","cpm"])
        writer.writeheader()
        for i in range(3):
            for j, sample in enumerate(samples):
                writer.writerow(dict(feature_id=f"gene{i+1}", sample=sample, cpm=cpm[i,j]))
    with (fitness / "treated_vs_parent.by_pool.log2fc_z.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["feature_id","contrast","parent_sample","treated_sample","log2FC"])
        writer.writeheader()
        for row in gv.rows(validation / "comparisons.tsv"):
            a, b = samples.index(row["parent_sample"]), samples.index(row["treated_sample"])
            for i in range(3):
                writer.writerow(dict(feature_id=f"gene{i+1}", contrast=row["comparison_id"], parent_sample=row["parent_sample"],
                    treated_sample=row["treated_sample"], log2FC=np.log2((cpm[i,b]+.5)/(cpm[i,a]+.5))))
    result = gv.evaluate(validation, tmp_path / "genes", existing_fitness=fitness.parent, plots=False)
    assert result["existing_fitness_cpm_parity"]["exact_cpm_rows"] == 12
    assert result["existing_fitness_cpm_parity"]["exact_default_pseudocount_effect_rows"] == 6
    with (fitness / "feature_reads_cpm.long.csv").open("a") as handle:
        handle.write("gene1,c1,333333.3333333333\n")
    with pytest.raises(ValueError, match="Duplicate rows"):
        gv.evaluate(validation, tmp_path / "bad", existing_fitness=fitness.parent, plots=False)


def test_quantile_site_results_are_projected_and_provenance_checked(tmp_path, monkeypatch):
    from ytab.normalization.methods import run_choices
    _, validation = fixture(tmp_path, monkeypatch)
    choices = tmp_path / "choices"
    run_choices(validation / "matrix", choices, ["quantile"], "library_id")
    output = tmp_path / "genes"
    result = gv.evaluate(validation, output, site_normalization=choices, plots=False)
    assert result["methods"] == ["nonorm","ttr","cpm","quantile"]
    assert (output / "effect_retention.tsv").exists()
    actual = gv.rows(output / "comparison_diagnostics.tsv")
    assert {r["method"] for r in actual} == {"nonorm","ttr","cpm","quantile"}
    with gzip.open(output / "quantile_gene_counts.tsv.gz", "rt") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    expected = sparse.load_npz(choices / "quantile/normalized_counts.npz").toarray()
    np.testing.assert_allclose([[float(r[k]) for k in ["c1","c2","t1","t2"]] for r in records], expected)
    manifest = choices / "run_manifest.json"
    value = json.loads(manifest.read_text()); value["source_hashes"]["raw_counts.npz"] = "wrong"
    manifest.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="same raw matrix"):
        gv.evaluate(validation, tmp_path / "bad", site_normalization=choices, plots=False)


def test_bgc_site_results_are_aggregated_without_cpm_rescaling(tmp_path, monkeypatch):
    from ytab.normalization.methods import run_choices
    _, validation = fixture(tmp_path, monkeypatch)
    choices = tmp_path / "choices"
    run_choices(validation / "matrix", choices, ["betageom"], "library_id")
    output = tmp_path / "genes"
    result = gv.evaluate(validation, output, site_normalization=choices, plots=False)
    assert result["methods"] == ["nonorm","ttr","cpm","betageom"]
    assert (output / "betageom_effect_retention.tsv").exists()
    with gzip.open(output / "betageom_gene_counts.tsv.gz", "rt") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    expected = sparse.load_npz(choices / "betageom/normalized_counts.npz").toarray()
    np.testing.assert_allclose([[float(r[k]) for k in ["c1","c2","t1","t2"]] for r in records], expected)


def test_adaptive_bgc_gene_aggregation_and_retention(tmp_path,monkeypatch):
    from ytab.normalization.methods import run_choices
    _,validation=fixture(tmp_path,monkeypatch)
    choices=tmp_path/"adaptive"
    run_choices(validation/"matrix",choices,["aBGC"],"library_id")
    output=tmp_path/"genes"
    result=gv.evaluate(validation,output,site_normalization=choices,plots=False)
    assert result["methods"]==["nonorm","ttr","cpm","aBGC"]
    assert (output/"aBGC_effect_retention.tsv").exists()
    with gzip.open(output/"aBGC_gene_counts.tsv.gz","rt") as handle:
        records=list(csv.DictReader(handle,delimiter="\t"))
    expected=sparse.load_npz(choices/"aBGC/normalized_counts.npz").toarray()
    np.testing.assert_allclose([[float(r[k]) for k in ["c1","c2","t1","t2"]] for r in records],expected)
