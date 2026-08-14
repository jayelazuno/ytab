#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --project-config output/projects/<id>/config/project.yaml [--threads N] [--force] [--keep-going]"
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON_BIN="${YTAB_PYTHON:-python}"
PROJECT_CONFIG=""
THREADS=""
FORCE=0
KEEP_GOING=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-config) PROJECT_CONFIG="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --keep-going) KEEP_GOING=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$PROJECT_CONFIG" ]]; then
  usage >&2
  exit 2
fi

common=(--project-config "$PROJECT_CONFIG")
if [[ -n "$THREADS" ]]; then common+=(--threads "$THREADS"); fi
if [[ "$FORCE" -eq 1 ]]; then common+=(--force); fi
if [[ "$KEEP_GOING" -eq 1 ]]; then common+=(--keep-going); fi

cd "$ROOT"
"$PYTHON_BIN" scripts/local/ytab_run_mapfastq.py "${common[@]}"
"$PYTHON_BIN" scripts/local/ytab_run_create_hit_file.py "${common[@]}"
"$PYTHON_BIN" scripts/local/ytab_project_status.py --project-config "$PROJECT_CONFIG" --show-next
