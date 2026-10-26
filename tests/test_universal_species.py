"""Species-agnostic classifier inputs, SummaryTable feature database and MidLC target units."""
from pathlib import Path

import pytest

from ytab.core.GenomicFeatures import GffFeatureDB
from ytab.essentiality import Classifier
from ytab.pipeline import classifier_runner, normalization_runner, summary_table_runner
from ytab.normalization.midlc import midlc


def test_orthology_formats(tmp_path):
    annotation = tmp_path/"gene_annotation.tsv"
    annotation.write_text("gene_id\tsystematic_id\tgene_name\tscer_relationship\tscer_gene_id\tscer_gene_name\tscer_essentiality\tdescription\n"
                          "C1_00060W_A\tC1_00060W\tTUP1\t1:1\tYCR084C\tTUP1\tnonessential\tx\n"
                          "C1_00070C_A\tC1_00070C\t\t\t\t\t\t\n")
    cgd = tmp_path/"cgd.txt"
    cgd.write_text("# CGD header\nC1_00060W_A\tTUP1\tCAL1\tYCR084C\tTUP1\tS000000680\nC2_1_A\t\tCAL2\tYAL001C\tTFC3\tS1\nC2_1_A\t\tCAL2\tYBR001C\tNTH2\tS2\n")
    bg2 = tmp_path/"bg2.txt"
    bg2.write_text("cglab_id\tscer_id\tscer_gene_name\tOrthogroup\tscer_sgdid\nGWK60_A1\tYNL318C\tHXT14\tOG1\tS5\n")
    a = Classifier.load_scer_orthology_map(annotation)
    assert a["C1_00060W_A"]["scer_id_list"] == ["YCR084C"] and a["C1_00070C_A"]["scer_id_list"] == []
    c = Classifier.load_scer_orthology_map(cgd)
    assert c["C2_1_A"]["scer_id_list"] == ["YAL001C", "YBR001C"] and c["C1_00060W_A"]["scer_gene_name"] == "TUP1"
    assert Classifier.load_scer_orthology_map(bg2)["GWK60_A1"]["scer_id"] == "YNL318C"


def test_gff_feature_db_any_style(tmp_path):
    fasta = tmp_path/"ref.fasta"; fasta.write_text(">chrI\n" + "A"*2000 + "\n")
    gff = tmp_path/"ref.gff"
    gff.write_text("##gff-version 3\n"
                   "chrI\tSGD\tgene\t100\t400\t.\t+\t.\tID=YAL001C;Name=YAL001C;gene=TFC3\n"
                   "chrI\tSGD\tmRNA\t100\t400\t.\t+\t.\tID=YAL001C_mRNA;Parent=YAL001C\n"
                   "chrI\tSGD\tCDS\t100\t200\t.\t+\t0\tParent=YAL001C_mRNA\n"
                   "chrI\tSGD\tCDS\t300\t400\t.\t+\t0\tParent=YAL001C_mRNA\n"
                   "chrI\tNCBI\tgene\t600\t900\t.\t-\t.\tID=gene-X;locus_tag=LT_1\n"
                   "chrI\tNCBI\tCDS\t600\t900\t.\t-\t0\tID=cds-X;Parent=gene-X;protein_id=P1\n"
                   "chrII\tSGD\tgene\t1\t90\t.\t+\t.\tID=OTHER\n"
                   "chrII\tSGD\tCDS\t1\t90\t.\t+\t0\tParent=OTHER\n")
    db = GffFeatureDB(str(gff), str(fasta))
    genes = {f.standard_name: f for f in db.get_all_features()}
    assert set(genes) == {"YAL001C", "LT_1"}          # chrII is not in the FASTA
    assert genes["YAL001C"].common_name == "TFC3" and (genes["YAL001C"].start, genes["YAL001C"].stop) == (100, 400)
    assert genes["LT_1"].strand == "C" and genes["LT_1"].protein_id == "P1"


def test_summary_table_uses_reference_gff_for_any_species(tmp_path):
    gff = tmp_path/"a.gff"; gff.write_text(""); fasta = tmp_path/"a.fasta"; fasta.write_text(">c\nA\n")
    reference = {"gff": gff, "fasta": fasta, "orthology_file": None}
    assert summary_table_runner.feature_db_args({"species": "newyeast"}, reference)[:2] == ["--feature-db", "gff"]
    assert summary_table_runner.feature_db_args({"species": "glabrata"}, reference)[:2] == ["--feature-db", "glabrata"]
    assert summary_table_runner.feature_db_errors({"species": "newyeast"}, reference) == []
    assert summary_table_runner.feature_db_errors({"species": "newyeast"}, {"gff": None, "fasta": None})


def test_classifier_runner_accepts_any_species(tmp_path):
    table = tmp_path/"combined_feature_table.T.txt"
    table.write_text("RDF,1\n" + ",".join(sorted(classifier_runner.REQUIRED_COLUMNS)) + "\n")
    errors = classifier_runner.validate_classifier_inputs({"species": "albicans"}, {"target_tag": "T"}, table, {})
    assert not any("species" in e and "only" in e for e in errors)


def test_midlc_target_units(tmp_path):
    hits = tmp_path/"s_hits.txt"
    hits.write_text("Chromosome\tSource\tHit position\tHit count\nc\tW\t10\t5\nc\tC\t10\t3\nc\tW\t20\t1\n")
    assert sorted(normalization_runner.site_counts(hits)) == [1, 8]
    saturation = {"a": {"ratio": 13.2}, "b": {"ratio": 73.0}}
    normalization_runner.check_target_units(20, saturation)              # thins b
    with pytest.raises(ValueError, match="multiple of MidLC"):
        normalization_runner.check_target_units(409600, saturation)       # the old read-count targets
    est, p = midlc([8, 1])
    assert 0 < p < 1 and est == pytest.approx(p*9)


def test_strain_reference_folders(tmp_path):
    """A strain with its own genome lives in <species>/strains/<strain> and is named <species>_<strain>."""
    from ytab.pipeline.reference_registry import list_available_species, species_resource_dir
    for folder in ("glabrata", "glabrata/strains/CBS138", "pombe"):
        (tmp_path/folder/"reference_genome").mkdir(parents=True)
        (tmp_path/folder/"reference_genome"/"ref.fasta").write_text(">c\nA\n")
    assert species_resource_dir(tmp_path, "glabrata") == tmp_path/"glabrata"
    assert species_resource_dir(tmp_path, "glabrata_CBS138") == tmp_path/"glabrata/strains/CBS138"
    assert species_resource_dir(tmp_path, "glabrata_cbs138") == tmp_path/"glabrata/strains/CBS138"
    assert species_resource_dir(tmp_path, "glabrata_BG2") == tmp_path/"glabrata_BG2"   # unknown strain: not found, as before
    assert list_available_species(tmp_path) == ["glabrata", "glabrata_CBS138", "pombe"]
