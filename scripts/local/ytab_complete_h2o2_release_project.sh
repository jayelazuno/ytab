#!/usr/bin/env bash
set -o pipefail

usage() {
  echo "Usage: $0 [--threads N] [--keep-going] [--force-downstream]"
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON_BIN="${YTAB_PYTHON:-python}"
THREADS="2"
KEEP_GOING=0
FORCE_DOWNSTREAM=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2 ;;
    --keep-going) KEEP_GOING=1; shift ;;
    --force-downstream) FORCE_DOWNSTREAM=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

cd "$ROOT"
PROJECT_CONFIG="output/projects/H2O2_screen_v1/config/project.yaml"
IMPORT_ARGS=(--project-id H2O2_screen_v1 --species glabrata --hit-dir output/projects/H2O2_screen_v1/create_hit_file --threads "$THREADS" --force)
DOWNSTREAM_ARGS=(--project-config "$PROJECT_CONFIG" --threads "$THREADS")
if [[ "$KEEP_GOING" -eq 1 ]]; then DOWNSTREAM_ARGS+=(--keep-going); fi
if [[ "$FORCE_DOWNSTREAM" -eq 1 ]]; then DOWNSTREAM_ARGS+=(--force); fi

YTAB_PYTHON="$PYTHON_BIN" "$PYTHON_BIN" scripts/local/ytab_import_existing_hit_project.py "${IMPORT_ARGS[@]}" || exit $?
YTAB_PYTHON="$PYTHON_BIN" bash scripts/local/ytab_run_downstream_from_hits.sh "${DOWNSTREAM_ARGS[@]}" || exit $?
"$PYTHON_BIN" scripts/local/ytab_project_status.py --project-config "$PROJECT_CONFIG" --show-next
