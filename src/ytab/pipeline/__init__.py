"""Local YTAB pipeline project setup utilities."""

from .project_config import (
    create_project_config,
    default_threads,
    detected_cpu_count,
    load_project_config,
    validate_project_config,
)
from .mapfastq_runner import (
    build_mapfastq_command,
    get_included_samples,
    load_project_for_mapping,
    run_mapfastq_project,
    run_mapfastq_sample,
    validate_mapping_inputs,
)
from .create_hit_file_runner import (
    build_create_hit_file_command,
    find_sample_bam,
    load_project_for_create_hit_file,
    run_create_hit_file_project,
    run_create_hit_file_sample,
    validate_create_hit_file_inputs,
)
from .summary_table_runner import (
    build_summary_table_command,
    collect_summary_stats,
    find_sample_hits_file,
    load_project_for_summary_table,
    run_summary_table_project,
    run_summary_table_sample,
    validate_summary_table_inputs,
)
from .reference_prepare import build_bowtie2_index, ensure_fasta_index, prepare_reference
from .reference_registry import (
    ReferenceInfo,
    list_available_species,
    reference_can_be_prepared,
    reference_is_runnable,
    reference_is_selectable,
    resolve_reference,
)
from .sample_discovery import discover_fastqs, write_sample_sheet

__all__ = [
    "ReferenceInfo",
    "build_bowtie2_index",
    "build_create_hit_file_command",
    "build_mapfastq_command",
    "build_summary_table_command",
    "collect_summary_stats",
    "create_project_config",
    "default_threads",
    "detected_cpu_count",
    "discover_fastqs",
    "ensure_fasta_index",
    "find_sample_bam",
    "find_sample_hits_file",
    "get_included_samples",
    "list_available_species",
    "load_project_config",
    "load_project_for_create_hit_file",
    "load_project_for_mapping",
    "load_project_for_summary_table",
    "prepare_reference",
    "reference_can_be_prepared",
    "reference_is_runnable",
    "reference_is_selectable",
    "resolve_reference",
    "run_mapfastq_project",
    "run_mapfastq_sample",
    "run_create_hit_file_project",
    "run_create_hit_file_sample",
    "run_summary_table_project",
    "run_summary_table_sample",
    "validate_create_hit_file_inputs",
    "validate_mapping_inputs",
    "validate_summary_table_inputs",
    "validate_project_config",
    "write_sample_sheet",
]
