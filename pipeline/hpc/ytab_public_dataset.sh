#!/bin/bash
# Download or process one public dataset described in pipeline/hpc/datasets/<dataset>/ (samples.csv, dataset.env).
# Run on the Argon login node from the YTAB checkout:
#   bash pipeline/hpc/ytab_public_dataset.sh <dataset> download          # qsub the download job
#   bash pipeline/hpc/ytab_public_dataset.sh <dataset> submit [options]  # create the project, submit its stage chain
#   bash pipeline/hpc/ytab_public_dataset.sh <dataset> status            # files downloaded and verified
# submit options go to argon_pipeline.py submit (default --profile core: reference, mapping, hit files, SummaryTable,
# library diagnostics, report and export). Add --dry-run to preview.
set -euo pipefail
YTAB="${YTAB:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}"
DATASET="${1:?Usage: ytab_public_dataset.sh <dataset> download|submit|status [options]}"
ACTION="${2:?Usage: ytab_public_dataset.sh <dataset> download|submit|status [options]}"
shift 2
DIR="$YTAB/pipeline/hpc/datasets/$DATASET"
[[ -f "$DIR/samples.csv" && -f "$DIR/dataset.env" ]] || { echo "ERROR: Unknown dataset: $DATASET ($(ls "$YTAB/pipeline/hpc/datasets" | tr '\n' ' '))" >&2; exit 1; }
source "$DIR/dataset.env"
DATA_ROOT="${DATA_ROOT:-/nfsscratch/jayelazuno/workspace/public_data}"
FASTQ_DIR="$DATA_ROOT/$DATASET/fastq"

case "$ACTION" in
  download)
    mkdir -p "$DATA_ROOT/$DATASET/logs"
    cd "$YTAB"
    qsub -v "DATASET=$DATASET,YTAB=$YTAB,DATA_ROOT=$DATA_ROOT" -N "YTAB_dl_${DATASET:0:20}" \
         -o "$DATA_ROOT/$DATASET/logs/qsub.log" "$YTAB/pipeline/hpc/ytab_download_public.txt"
    ;;
  status)
    files="$DATA_ROOT/$DATASET/metadata/files.tsv"
    [[ -f "$files" ]] || { echo "Not resolved yet: $files"; exit 0; }
    expected=$(( $(wc -l < "$files") - 1 ))
    done_count=$(find "$FASTQ_DIR" -name '*.fastq.gz.complete' 2>/dev/null | wc -l)
    echo "$DATASET: $done_count of $expected files downloaded and verified ($FASTQ_DIR)"
    ;;
  submit)
    files="$DATA_ROOT/$DATASET/metadata/files.tsv"
    expected=$(( $(wc -l < "$files") - 1 ))
    done_count=$(find "$FASTQ_DIR" -name '*.fastq.gz.complete' | wc -l)
    (( done_count == expected )) || { echo "ERROR: Only $done_count of $expected files are verified; finish the download first." >&2; exit 1; }
    exec bash "$YTAB/pipeline/hpc/submit_pipeline.sh" --project-id "$PROJECT_ID" --fastq-dir "$FASTQ_DIR" \
         --species "$SPECIES" --sample-sheet "$DIR/samples.csv" "$@"
    ;;
  *) echo "ERROR: Action must be download, submit or status." >&2; exit 1 ;;
esac
