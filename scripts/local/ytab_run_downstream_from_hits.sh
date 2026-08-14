#!/usr/bin/env bash
set -o pipefail

usage() {
  echo "Usage: $0 --project-config output/projects/<id>/config/project.yaml [--threads N] [--force] [--keep-going]"
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON_BIN="${YTAB_PYTHON:-python}"
PROJECT_CONFIG=""
THREADS="2"
FORCE=0
KEEP_GOING=0
FAILURES=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-config) PROJECT_CONFIG="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force|--force-downstream) FORCE=1; shift ;;
    --keep-going) KEEP_GOING=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$PROJECT_CONFIG" ]]; then
  usage >&2
  exit 2
fi

run_stage() {
  local label="$1"
  shift
  echo
  echo "== $label =="
  "$@"
  local status=$?
  if [[ "$status" -ne 0 ]]; then
    FAILURES=$((FAILURES + 1))
    echo "Stage failed ($status): $label" >&2
    if [[ "$KEEP_GOING" -ne 1 ]]; then
      exit "$status"
    fi
  fi
}

force_arg=()
if [[ "$FORCE" -eq 1 ]]; then force_arg=(--force); fi
keep_arg=()
if [[ "$KEEP_GOING" -eq 1 ]]; then keep_arg=(--keep-going); fi

cd "$ROOT"

run_stage "Raw SummaryTable" "$PYTHON_BIN" scripts/local/ytab_run_summary_table.py --project-config "$PROJECT_CONFIG" --threads "$THREADS" "${force_arg[@]}" "${keep_arg[@]}"
run_stage "Raw LibraryDiagnostics" "$PYTHON_BIN" scripts/local/ytab_run_library_diagnostics.py --project-config "$PROJECT_CONFIG" --threads "$THREADS" "${force_arg[@]}"
run_stage "Parent-only MidLC normalization" "$PYTHON_BIN" scripts/local/ytab_run_sample_normalization.py --project-config "$PROJECT_CONFIG" --targets auto --sample-mode parents --threads "$THREADS" "${force_arg[@]}" "${keep_arg[@]}"
run_stage "Normalized SummaryTable target evaluation" "$PYTHON_BIN" scripts/local/ytab_run_summary_normalized.py --project-config "$PROJECT_CONFIG" --targets recommended --sample-mode parents --threads "$THREADS" "${force_arg[@]}" "${keep_arg[@]}"
run_stage "Combine normalized parent libraries" "$PYTHON_BIN" scripts/local/ytab_run_combine_hits.py --project-config "$PROJECT_CONFIG" --target recommended "${force_arg[@]}"
run_stage "Combined parent SummaryTable" "$PYTHON_BIN" scripts/local/ytab_run_summary_combined.py --project-config "$PROJECT_CONFIG" --target recommended --threads "$THREADS" "${force_arg[@]}"
run_stage "Essentiality classifier" "$PYTHON_BIN" scripts/local/ytab_run_classifier.py --project-config "$PROJECT_CONFIG" --target recommended "${force_arg[@]}"
run_stage "Working classifier target record" "$PYTHON_BIN" scripts/local/ytab_record_working_classifier_target.py --project-config "$PROJECT_CONFIG" --target recommended
run_stage "Fitness comparison design" "$PYTHON_BIN" scripts/local/ytab_init_comparison_design.py --project-config "$PROJECT_CONFIG" --overwrite
run_stage "Treated-versus-parent fitness" "$PYTHON_BIN" scripts/local/ytab_run_treated_vs_parent.py --project-config "$PROJECT_CONFIG" --analysis-id H2O2_vs_parent --input-mode raw-summary "${force_arg[@]}" "${keep_arg[@]}"
run_stage "Project status refresh" "$PYTHON_BIN" scripts/local/ytab_project_status.py --project-config "$PROJECT_CONFIG" --show-next
run_stage "Project report" "$PYTHON_BIN" scripts/local/ytab_build_project_report.py --project-config "$PROJECT_CONFIG" --force
run_stage "Project export" "$PYTHON_BIN" scripts/local/ytab_export_project.py --project-config "$PROJECT_CONFIG" --force --include-hit-files

if [[ "$FAILURES" -ne 0 ]]; then
  echo
  echo "Completed with $FAILURES downstream failure(s)." >&2
  exit 1
fi

echo
echo "Downstream analyses completed from existing hit files."
