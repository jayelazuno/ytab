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
    Module("preprocess", "ytab_run_preprocessing.py", "Mapping, hit files and SummaryTable for selected samples in one run (skips stages already done)", threads=True),
    Module("mapfastq", "ytab_run_mapfastq.py", "FASTQ clipping and reference alignment", threads=True),
    Module("hit-files", "ytab_run_create_hit_file.py", "Create unnormalized insertion hit files", threads=True),
    Module("summary", "ytab_run_summary_table.py", "Raw per-sample gene summaries", threads=True),
    Module("library-qc", "ytab_run_library_diagnostics.py", "Library diagnostics and insertion architecture", threads=True),
    Module("qc", "ytab_qc.py", "Library QC table for a project: all QC metrics, flagged (rule of thumb / project-relative) and commented, plus plot data", force_option=None),
    Module("tnseq-stats", "ytab_tnseq_stats.py", "Per-library Tn-seq count statistics plus saved Hermes QC", project_option=None, force_option=None),
    Module("normalization", "ytab_run_sample_normalization.py", "Legacy parent-only MidLC thinning to a fixed target (superseded by the midlc method of normalize-project)", threads=True),
    Module("ttr", "ytab_ttr.py", "Empirical site catalogue and TTR normalization evaluation", project_option=None, force_option=None),
    Module("normalize", "ytab_normalize.py", "Explicit site-normalization choices and method guidance", project_option=None, force_option=None),
    Module("hybrid-inputs", "ytab_hybrid_inputs.py", "Hybrid essentiality gene features from parents normalized with any method (default BGC; reuses earlier results)", threads=True, force_option=None),
    Module("normalize-project", "ytab_normalize_project.py", "Normalize a project's sites with chosen methods in one step (catalogue, experiment matrix, plot summaries)", project_option=None, force_option=None),
    Module("normalize-experiment", "ytab_normalize_experiment.py", "Evaluate normalization on design-defined controls and treatments", project_option=None, force_option=None),
    Module("normalize-genes", "ytab_normalize_genes.py", "Gene-level raw/TTR/CPM validation without changing fitness defaults", project_option=None, force_option=None),
    Module("bgc-audit", "ytab_bgc_audit.py", "BGC numerical, seed and catalogue-sensitivity diagnostics", project_option=None, force_option=None),
    Module("abgc-audit", "ytab_abgc_audit.py", "Adaptive BGC search, tail and simulation sensitivity diagnostics", project_option=None, force_option=None),
    Module("summary-normalized", "ytab_run_summary_normalized.py", "Summaries of normalized parent libraries", threads=True),
    Module("combine-hits", "ytab_run_combine_hits.py", "Combine independent parent insertion tracks"),
    Module("summary-combined", "ytab_run_summary_combined.py", "Combined-parent gene summaries", threads=True),
    Module("classifier", "ytab_run_classifier.py", "Existing essentiality classifier"),
    Module("resampling", "ytab_resampling.py", "Hermes site-restricted resampling on a ttr/totreads/nzmean normalization output", project_option=None, force_option=None),
    Module("control-gumbel", "ytab_control_conditioned_gumbel.py", "Experimental control-conditioned count-aware Hermes Gumbel on a ttr/totreads/nzmean normalization output", project_option=None, force_option=None),
    Module("multicondition", "ytab_multicondition.py", "Hermes multi-condition tests (permutation ANOVA and all-pairs resampling) on a normalization output", project_option=None, force_option=None),
    Module("multicondition-project", "ytab_multicondition_project.py", "Multi-condition tests and their plots in one step for a project (reference from the project config)", project_option=None, force_option=None),
    Module("multicondition-plots", "ytab_multicondition_plots.py", "Response correlation (corrplot) and response heatmap for a multicondition run", project_option=None, force_option=None),
    Module("genetic-interactions", "ytab_genetic_interactions.py", "Hermes genetic interactions: strain x treatment response differences", project_option=None, force_option=None),
    Module("essentiality", "ytab_essentiality.py", "Essentiality for a project in one step: Hermes HMM (label-free), current RF scores and the hybrid when labels are given", project_option=None, force_option=None),
    Module("cross-species", "ytab_cross_species.py", "Essentiality across species: HMM calls of one-to-one orthologs (S. cerevisiae hub) in two or more projects", project_option=None, force_option=None),
    Module("hermes-hmm", "ytab_hermes_hmm.py", "Label-free Hermes window HMM essentiality for any yeast (haploid reference FASTA + GFF + hit files)", project_option=None, force_option=None),
    Module("resource-bundle", "ytab_resource_bundle.py", "Bundle reprocessed projects (results, no FASTQ/BAM) with a checksum manifest, to open them in the app elsewhere", project_option=None, force_option=None),
    Module("datasets", "ytab_datasets.py", "Catalogue of the reprocessed public datasets (study, strain, transposon, accession, QC, results in the app) for the Datasets page", project_option=None, force_option=None),
    Module("library-correlation", "ytab_library_correlation.py", "Correlate Hermes libraries across pools, strains, conditions and species (gene insertion density)", project_option=None, force_option=None),
    Module("comparison-design", "ytab_init_comparison_design.py", "Create explicit fitness comparisons", force_option="--overwrite"),
    Module("fitness", "ytab_run_treated_vs_parent.py", "Treated vs control, local-abundance z-score for a project (CPM by default, or any site normalization)"),
    Module("fitness-zscore", "ytab_fitness_zscore.py", "Local-abundance z-score on explicit inputs (feature tables or a normalization run)", project_option=None, force_option=None),
    Module("tvc-benchmark", "ytab_tvc_benchmark.py", "Pre-registered treated-vs-control benchmark: Gale-calibrated simulation, method x normalization arms, paired tests", project_option=None, force_option=None),
    Module("browser-tracks", "ytab_browser_tracks.py", "Genome browser files for a project: reference with display chromosome names, gene BED, insertion tracks and an IGV session", project_option=None, force_option=None),
    Module("tvc-phase2", "ytab_tvc_phase2.py", "Treated-vs-control benchmark phase 2: every arm on real Gale FLU, a swap null, H2O2/Zn, and Gale's published scores", project_option=None, force_option=None),
    Module("treated-vs-control", "ytab_treated_vs_control.py", "Treated vs control for a project by resampling or the control-conditioned Gumbel (normalization run + comparison design)", project_option=None, force_option=None),
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
ALIASES.update({"treated-vs-parent": "fitness", "screen-coverage": "coverage", "mutant-identification": "mutants",
                "anova": "multicondition", "all-pairs": "multicondition", "corrplot": "multicondition-plots", "heatmap": "multicondition-plots", "gi": "genetic-interactions"})


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
