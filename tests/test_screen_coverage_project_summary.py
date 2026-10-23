import csv
import json
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from ytab.screen_coverage.project_summary import parent_rows, project_summary


def test_parent_selection_is_metadata_driven():
    rows = [dict(sample="parent", include="True", library_role="parent"),
            dict(sample="control", include="True", library_role="control"),
            dict(sample="parent-treated", include="True", library_role="treated"),
            dict(sample="excluded", include="False", library_role="parent")]
    assert [r["sample"] for r in parent_rows(rows)] == ["parent"]
    assert [r["sample"] for r in parent_rows(rows[1:])] == ["control"]
    assert parent_rows([rows[2]]) == []


def make_project(root, name):
    project = root / name
    config = project / "config"
    config.mkdir(parents=True)
    gff = root / "reference.gff"
    gff.write_text("chr\ttest\tgene\t1\t5\t.\t+\t.\tlocus_tag=pseudo;gene_biotype=pseudogene\n")
    sheet = config / "sample_sheet.csv"
    with sheet.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "include", "library_role"])
        writer.writeheader()
        writer.writerow(dict(sample=name+"-parent", include="True", library_role="parent"))
    hits = project / "create_hit_file" / (name+"-parent") / (name+"-parent_hits.txt")
    hits.parent.mkdir(parents=True)
    hits.write_text("placeholder")
    path = config / "project.yaml"
    path.write_text(yaml.safe_dump(dict(project_id=name, repo_root=str(root),
        output_project_dir=name, sample_sheet=str(sheet), reference=dict(gff=str(gff)))))
    return path, hits


def test_inventory_cache_is_project_scoped_and_never_uses_other_audits(tmp_path):
    first, hit = make_project(tmp_path, "study-A")
    second, _ = make_project(tmp_path, "study-B")
    library = dict(genes=[dict(gene_id="coding",informative_sites=2),dict(gene_id="pseudo",informative_sites=3)],
                   unique_sites=7, sources=[],informative_definition="coding blocks")
    with patch("ytab.screen_coverage.project_summary.load_library", return_value=library) as adapter:
        a = project_summary(first)
        b = project_summary(second)
        assert a["project_id"] == "study-A" and b["project_id"] == "study-B"
        assert a["samples"] == ["study-A-parent"] and b["samples"] == ["study-B-parent"]
        assert a["annotated"] == a["represented"] == a["ge2"] == 1
        assert a["ge3"] == 0 and a["coding_sites"] == 2
        assert a["signature"] != b["signature"]
        assert project_summary(first) == a and adapter.call_count == 2
        hit.write_text("changed insertion catalogue")
        assert project_summary(first)["signature"] != a["signature"]
        assert adapter.call_count == 3
    hit.unlink()
    with pytest.raises(ValueError, match="missing"):
        project_summary(first)
