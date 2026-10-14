"""Discoverable CLI modules; all entries delegate to existing YTAB backends."""
from dataclasses import asdict, dataclass
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class Module:
    name: str
    script: str
    description: str
    project_option: str | None = "--project-config"
    threads: bool = False
    force_option: str | None = "--force"
    workflow: bool = True
    interpreter: str = "python"


MODULES = {m.name: m for m in (
    Module("pipeline", "ytab_run_pipeline.py", "Existing core/classifier/fitness/all project pipeline", threads=True),
    Module("reference", "ytab_prepare_reference.py", "Prepare species reference resources", project_option=None, threads=True),
    Module("mapfastq", "ytab_run_mapfastq.py", "FASTQ clipping and reference alignment", threads=True),
    Module("hit-files", "ytab_run_create_hit_file.py", "Create unnormalized insertion hit files", threads=True),
    Module("summary", "ytab_run_summary_table.py", "Raw per-sample gene summaries", threads=True),
    Module("library-qc", "ytab_run_library_diagnostics.py", "Library diagnostics and insertion architecture", threads=True),
    Module("normalization", "ytab_run_sample_normalization.py", "Existing parent-only MidLC normalization", threads=True),
    Module("summary-normalized", "ytab_run_summary_normalized.py", "Summaries of normalized parent libraries", threads=True),
    Module("combine-hits", "ytab_run_combine_hits.py", "Combine independent parent insertion tracks"),
    Module("summary-combined", "ytab_run_summary_combined.py", "Combined-parent gene summaries", threads=True),
    Module("classifier", "ytab_run_classifier.py", "Existing essentiality classifier"),
    Module("comparison-design", "ytab_init_comparison_design.py", "Create explicit fitness comparisons", force_option="--overwrite"),
    Module("fitness", "ytab_run_treated_vs_parent.py", "Raw-summary CPM treated-versus-control fitness"),
    Module("go", "ytab_run_go_analysis.py", "GO enrichment and optional semantic reduction"),
    Module("gene-insertions", "ytab_draw_gene_domain_insertions.py", "Gene/domain insertion plots and queries"),
    Module("coverage", "ytab_poisson_coverage.py", "Analytical Poisson / equal-gene screen benchmark", project_option=None),
    Module("mutants", "ytab_identify_mutants.py", "Universal paired Sanger/ABI mutant identification", project_option="--reference-config", threads=True, force_option=None),
    Module("status", "ytab_project_status.py", "Inspect resumable project stage status", force_option=None),
    Module("report", "ytab_build_project_report.py", "Local HTML report of existing project outputs"),
    Module("export", "ytab_export_project.py", "Portable project-result export"),
    Module("init", "ytab_init_local_project.py", "Initialize a project from local FASTQs", project_option=None, force_option=None),
    Module("import-hits", "ytab_import_existing_hit_project.py", "Initialize a project from existing hit files", project_option=None, force_option=None),
    Module("app", "ytab_launch_app.sh", "Open the Shiny app on existing results", workflow=False, force_option=None, interpreter="bash"),
)}

ALIASES = {name.replace("-", "_"): name for name in MODULES if "-" in name}
ALIASES.update({"treated-vs-parent": "fitness", "screen-coverage": "coverage", "mutant-identification": "mutants"})


def module(name):
    name = ALIASES.get(name, name)
    if name not in MODULES:
        raise ValueError(f"Unknown module {name!r}. Use 'ytab modules' to list available modules.")
    return MODULES[name]


def command(name, arguments):
    selected = module(name)
    path = ROOT / "scripts" / "local" / selected.script
    if not path.is_file():
        raise ValueError(f"Missing runner for {name}: {path}")
    interpreter = sys.executable if selected.interpreter == "python" else selected.interpreter
    return [interpreter, str(path), *arguments]


def inventory():
    return [asdict(m) for m in MODULES.values()]
