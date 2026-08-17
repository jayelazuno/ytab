#!/bin/bash
# Submit the full Zn_toxicity_screen YTAB pipeline on SGE/qsub.

set -euo pipefail

PROJECT_ID="Zn_toxicity_screen"
SAMPLE_SHEET=""
FASTQ_DIR="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110"
SCRATCH_WORK_DIR="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/ytab_work/Zn_toxicity_screen"
REPO_ROOT="/Users/jayelazuno/workspace/ytab"
THREADS="8"
QUEUE="UI"
KEEP_GOING="0"
DRY_RUN="0"
FORCE_DOWNSTREAM="0"
YTAB_CONDA_ENV="${YTAB_CONDA_ENV:-ytab-local}"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/hpc/ytab_submit_zn_toxicity_screen_hpc.sh [options]

Options:
  --project-id Zn_toxicity_screen
  --sample-sheet PATH
  --fastq-dir PATH
  --scratch-work-dir PATH
  --repo-root PATH
  --threads 8
  --queue UI
  --conda-env ytab-local
  --keep-going
  --dry-run
  --force-downstream
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-id) PROJECT_ID="$2"; shift 2 ;;
    --sample-sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --scratch-work-dir) SCRATCH_WORK_DIR="$2"; shift 2 ;;
    --repo-root) REPO_ROOT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --queue) QUEUE="$2"; shift 2 ;;
    --conda-env) YTAB_CONDA_ENV="$2"; shift 2 ;;
    --keep-going) KEEP_GOING="1"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    --force-downstream) FORCE_DOWNSTREAM="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$PROJECT_ID" != "Zn_toxicity_screen" ]]; then
  echo "ERROR: project ID must be Zn_toxicity_screen, not $PROJECT_ID" >&2
  exit 2
fi

if [[ -z "$SAMPLE_SHEET" ]]; then
  if [[ -s "$FASTQ_DIR/sampleSheetUsed.csv" ]]; then
    SAMPLE_SHEET="$FASTQ_DIR/sampleSheetUsed.csv"
  else
    SAMPLE_SHEET="$REPO_ROOT/codex/sampleSheetUsed.csv"
  fi
fi

PROJECT_DIR="$REPO_ROOT/output/projects/$PROJECT_ID"
EXPORT_DIR="$REPO_ROOT/output/exports/$PROJECT_ID"
PROJECT_CONFIG="$PROJECT_DIR/config/project.yaml"
MANIFEST="$PROJECT_DIR/config/hpc_sample_manifest.csv"
HPC_SCRIPT_DIR="$REPO_ROOT/pipeline/hpc/zn_toxicity_screen"
PYTHON="${PYTHON:-python}"

if [[ "$DRY_RUN" != "1" && "$(pwd -P)" != "$REPO_ROOT" ]]; then
  echo "ERROR: run from $REPO_ROOT or pass --repo-root for this clone." >&2
  echo "Current directory: $(pwd -P)" >&2
  exit 2
fi

mkdir -p "$SCRATCH_WORK_DIR/logs" "$SCRATCH_WORK_DIR/mapfastq" "$SCRATCH_WORK_DIR/bam" "$SCRATCH_WORK_DIR/tmp"
mkdir -p "$PROJECT_DIR/logs" "$EXPORT_DIR"

echo "Project ID:       $PROJECT_ID"
echo "Project title:    Zinc Toxicity Screen"
echo "Repo root:        $REPO_ROOT"
echo "Sample sheet:     $SAMPLE_SHEET"
echo "FASTQ dir:        $FASTQ_DIR"
echo "Scratch work dir: $SCRATCH_WORK_DIR"
echo "Threads:          $THREADS"
echo "Queue:            $QUEUE"
echo "Conda env:        $YTAB_CONDA_ENV"
echo "Dry run:          $DRY_RUN"
echo

if [[ "$DRY_RUN" == "1" && ! -s "$SAMPLE_SHEET" ]]; then
  echo "[DRY-RUN WARNING] sample sheet is not visible here: $SAMPLE_SHEET"
  echo "[DRY-RUN WARNING] qsub dependency chain will still be printed with N=8."
else
  "$PYTHON" "$REPO_ROOT/scripts/hpc/ytab_make_hpc_sample_manifest.py" \
    --project-id "$PROJECT_ID" \
    --sample-sheet "$SAMPLE_SHEET" \
    --fastq-dir "$FASTQ_DIR" \
    --scratch-work-dir "$SCRATCH_WORK_DIR" \
    --species glabrata \
    --read-layout paired_end \
    --repo-root "$REPO_ROOT"
fi

if [[ -s "$MANIFEST" ]]; then
  SAMPLE_COUNT="$("$PYTHON" "$REPO_ROOT/scripts/hpc/ytab_hpc_manifest_row.py" --manifest "$MANIFEST" --task-id 1 --count)"
else
  SAMPLE_COUNT="8"
fi

if [[ "$SAMPLE_COUNT" != "8" ]]; then
  echo "ERROR: expected 8 Zn toxicity samples, found $SAMPLE_COUNT" >&2
  exit 2
fi

COMMON_VARS="PROJECT_ID=$PROJECT_ID,REPO_ROOT=$REPO_ROOT,PROJECT_CONFIG=$PROJECT_CONFIG,MANIFEST=$MANIFEST,THREADS=$THREADS,FORCE=$FORCE_DOWNSTREAM,YTAB_CONDA_ENV=$YTAB_CONDA_ENV"

submit_job() {
  local label="$1"
  local script="$2"
  local hold="${3:-}"
  local qsub_cmd=(qsub -q "$QUEUE" -pe smp "$THREADS" -v "$COMMON_VARS")
  if [[ "$script" == *array.sge ]]; then
    qsub_cmd+=(-t "1-$SAMPLE_COUNT")
  fi
  if [[ -n "$hold" ]]; then
    qsub_cmd+=(-hold_jid "$hold")
  fi
  qsub_cmd+=("$script")
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "DRY-RUN $label:" >&2
    printf '  %q' "${qsub_cmd[@]}" >&2
    echo >&2
    echo "$label"
  else
    echo "Submitting $label..." >&2
    local output
    output="$("${qsub_cmd[@]}")"
    echo "$output" >&2
    awk '{print $3}' <<<"$output"
  fi
}

map_job="$(submit_job YTAB_Zn_MapFastq "$HPC_SCRIPT_DIR/01_mapfastq_array.sge")"
hit_job="$(submit_job YTAB_Zn_CreateHitFile "$HPC_SCRIPT_DIR/02_create_hit_file_array.sge" "$map_job")"
summary_job="$(submit_job YTAB_Zn_SummaryRaw "$HPC_SCRIPT_DIR/03_summary_raw.sge" "$hit_job")"
qc_job="$(submit_job YTAB_Zn_LibraryDiagnostics "$HPC_SCRIPT_DIR/04_library_diagnostics.sge" "$summary_job")"
classifier_job="$(submit_job YTAB_Zn_ClassifierBranch "$HPC_SCRIPT_DIR/05_classifier_branch.sge" "$summary_job")"
fitness_job="$(submit_job YTAB_Zn_Fitness "$HPC_SCRIPT_DIR/06_fitness_zn_vs_mock.sge" "$summary_job")"
final_hold="$qc_job,$classifier_job,$fitness_job"
final_job="$(submit_job YTAB_Zn_Finalize "$HPC_SCRIPT_DIR/07_report_export_finalize.sge" "$final_hold")"

echo
echo "Submitted dependency chain:"
echo "  manifest/config"
echo "    -> $map_job"
echo "    -> $hit_job"
echo "    -> $summary_job"
echo "    -> $qc_job"
echo "    -> $classifier_job"
echo "    -> $fitness_job"
echo "    -> $final_job"
echo
echo "Follow with:"
echo "  qstat -u jayelazuno"
echo "  qstat -j $final_job"
echo
echo "After completion, pull app-facing outputs only. Do not pull scratch BAMs by default:"
echo '  rsync -avP \'
echo "    jayelazuno@<HPC_HOST>:$REPO_ROOT/output/projects/$PROJECT_ID \\"
echo '    "/Volumes/rdss_bhe2/User/Joshua-Ayelazuno/Repositories /ytab/output/projects/"'
echo
echo '  rsync -avP \'
echo "    jayelazuno@<HPC_HOST>:$REPO_ROOT/output/exports/$PROJECT_ID \\"
echo '    "/Volumes/rdss_bhe2/User/Joshua-Ayelazuno/Repositories /ytab/output/exports/"'
