#!/bin/bash
# Submit the Zn toxicity R2-vs-Hermes BLAST diagnostic on SGE/qsub.

set -euo pipefail

REPO_ROOT="/old_Users/jayelazuno/workspace/ytab"
R2_DIR="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110"
SAMPLE_SHEET="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110/sampleSheetUsed.csv"
BLAST_ROOT="/old_Users/jayelazuno/workspace/ytab/blast"
ALIGN_ROOT="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/ytab_work/Zn_toxicity_screen/blast_alignments"
HERMES_CDS="/old_Users/jayelazuno/workspace/ytab/blast/Hermes-cds.txt"
YTAB_CONDA_ENV="${YTAB_CONDA_ENV:-ytab-hermes-blast}"
THREADS="8"
QUEUE="UI"
MIN_PIDENT="90"
MIN_ALIGN_LEN="30"
MAX_EVALUE="1e-10"
MAX_READS="0"
DRY_RUN="0"
FORCE="0"
SMOKE="0"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/hpc/ytab_submit_zn_r2_hermes_blast.sh [options]

Options:
  --repo-root PATH
  --r2-dir PATH
  --sample-sheet PATH
  --blast-root PATH
  --align-root PATH
  --hermes-cds PATH
  --conda-env ytab-hermes-blast
  --threads 8
  --queue UI
  --min-pident 90
  --min-align-len 30
  --max-evalue 1e-10
  --max-reads N
  --smoke
  --force
  --dry-run
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-root) REPO_ROOT="$2"; shift 2 ;;
    --r2-dir) R2_DIR="$2"; shift 2 ;;
    --sample-sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --blast-root) BLAST_ROOT="$2"; shift 2 ;;
    --align-root) ALIGN_ROOT="$2"; shift 2 ;;
    --hermes-cds) HERMES_CDS="$2"; shift 2 ;;
    --conda-env) YTAB_CONDA_ENV="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --queue) QUEUE="$2"; shift 2 ;;
    --min-pident) MIN_PIDENT="$2"; shift 2 ;;
    --min-align-len) MIN_ALIGN_LEN="$2"; shift 2 ;;
    --max-evalue) MAX_EVALUE="$2"; shift 2 ;;
    --max-reads) MAX_READS="$2"; shift 2 ;;
    --smoke) SMOKE="1"; shift ;;
    --force) FORCE="1"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$SMOKE" == "1" && "$MAX_READS" == "0" ]]; then
  MAX_READS="10000"
fi

activate_env() {
  if [[ "${YTAB_SKIP_CONDA:-0}" == "1" ]]; then return 0; fi
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
  else
    for conda_sh in "${CONDA_SH:-}" \
      /old_Users/jayelazuno/miniforge3/etc/profile.d/conda.sh \
      /Users/jayelazuno/miniforge3/etc/profile.d/conda.sh; do
      if [[ -n "$conda_sh" && -s "$conda_sh" ]]; then
        # shellcheck disable=SC1090
        source "$conda_sh"
        break
      fi
    done
  fi
  conda activate "$YTAB_CONDA_ENV"
}

find_hermes() {
  for candidate in "$HERMES_CDS" \
    "/old_Users/jayelazuno/workspace/ytab/blast/Hermes-cds.txt" \
    "$REPO_ROOT/docs/codex/Hermes-cds.txt"; do
    if [[ -s "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

validate_and_copy_hermes() {
  local source="$1"
  local dest="$2"
  python - "$source" "$dest" <<'PY'
from pathlib import Path
import re
import sys
source = Path(sys.argv[1])
dest = Path(sys.argv[2])
if not source.is_file() or source.stat().st_size == 0:
    raise SystemExit(f"ERROR: Hermes CDS file missing or empty: {source}")
text = source.read_text(encoding="utf-8", errors="replace").strip()
if not text.startswith(">"):
    raise SystemExit(f"ERROR: Hermes CDS does not appear to be FASTA: {source}")
seq = "".join(line.strip() for line in text.splitlines() if line and not line.startswith(">"))
if not seq:
    raise SystemExit("ERROR: Hermes CDS FASTA has no sequence")
letters = re.sub(r"[^A-Za-z]", "", seq).upper()
nucleotide = sum(1 for c in letters if c in "ACGTUN")
if not letters or nucleotide / len(letters) < 0.95:
    raise SystemExit("ERROR: Hermes CDS does not appear nucleotide-like enough for blastn")
dest.parent.mkdir(parents=True, exist_ok=True)
dest.write_text(text + "\n", encoding="utf-8")
print(f"Validated nucleotide FASTA: {source}")
PY
}

if [[ "$DRY_RUN" != "1" && "$(pwd -P)" != "$REPO_ROOT" ]]; then
  echo "ERROR: run from $REPO_ROOT or pass --repo-root for this clone." >&2
  echo "Current directory: $(pwd -P)" >&2
  exit 2
fi

if [[ "$DRY_RUN" != "1" ]]; then
  activate_env
  for tool in python makeblastdb blastn seqkit Rscript; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing required tool: $tool" >&2; exit 2; }
  done
fi

mkdir -p \
  "$BLAST_ROOT/reference" "$BLAST_ROOT/db" "$BLAST_ROOT/manifest" "$BLAST_ROOT/summary/per_sample" \
  "$BLAST_ROOT/figures" "$BLAST_ROOT/logs" "$BLAST_ROOT/scripts_metadata" \
  "$ALIGN_ROOT/fasta_queries" "$ALIGN_ROOT/blast_outfmt6" "$ALIGN_ROOT/unique_hit_reads" "$ALIGN_ROOT/seqkit_stats"

HERMES_SOURCE="$(find_hermes || true)"
if [[ -z "$HERMES_SOURCE" ]]; then
  echo "ERROR: Hermes CDS file not found. Checked --hermes-cds, blast/Hermes-cds.txt, and docs/codex/Hermes-cds.txt." >&2
  exit 2
fi
HERMES_FASTA="$BLAST_ROOT/reference/Hermes-cds.fa"
DB_PREFIX="$BLAST_ROOT/db/hermes_cds"
MANIFEST="$BLAST_ROOT/manifest/r2_fastq_manifest.tsv"

if [[ "$DRY_RUN" != "1" ]]; then
  validate_and_copy_hermes "$HERMES_SOURCE" "$HERMES_FASTA"
  if [[ "$FORCE" == "1" ]] || ! compgen -G "${DB_PREFIX}.n*" >/dev/null; then
    makeblastdb -in "$HERMES_FASTA" -dbtype nucl -parse_seqids -out "$DB_PREFIX"
  else
    echo "[SKIP] Existing BLAST database found at $DB_PREFIX"
  fi
  python "$REPO_ROOT/scripts/hpc/ytab_make_r2_hermes_blast_manifest.py" \
    --r2-dir "$R2_DIR" \
    --sample-sheet "$SAMPLE_SHEET" \
    --output "$MANIFEST"
fi

if [[ -s "$MANIFEST" ]]; then
  SAMPLE_COUNT="$(python - "$MANIFEST" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    print(sum(1 for _ in csv.DictReader(handle, delimiter="\t")))
PY
)"
else
  SAMPLE_COUNT="0"
fi
if [[ "$DRY_RUN" == "1" && "$SAMPLE_COUNT" -lt 1 ]]; then
  SAMPLE_COUNT="8"
fi
if [[ "$DRY_RUN" != "1" && "$SAMPLE_COUNT" -lt 1 ]]; then
  echo "ERROR: R2 manifest has no rows: $MANIFEST" >&2
  exit 2
fi

cat > "$BLAST_ROOT/scripts_metadata/run_config.tsv" <<EOF
repo_root	$REPO_ROOT
r2_dir	$R2_DIR
sample_sheet	$SAMPLE_SHEET
blast_root	$BLAST_ROOT
align_root	$ALIGN_ROOT
hermes_cds_source	$HERMES_SOURCE
hermes_fasta	$HERMES_FASTA
db_prefix	$DB_PREFIX
manifest	$MANIFEST
min_pident	$MIN_PIDENT
min_align_len	$MIN_ALIGN_LEN
max_evalue	$MAX_EVALUE
max_reads	$MAX_READS
threads	$THREADS
conda_env	$YTAB_CONDA_ENV
EOF

python - "$BLAST_ROOT/scripts_metadata/run_config.tsv" "$BLAST_ROOT/scripts_metadata/run_config.json" <<'PY'
import json
import sys
items = {}
with open(sys.argv[1], encoding="utf-8") as handle:
    for line in handle:
        key, value = line.rstrip("\n").split("\t", 1)
        items[key] = value
with open(sys.argv[2], "w", encoding="utf-8") as handle:
    json.dump(items, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

COMMON_VARS="REPO_ROOT=$REPO_ROOT,BLAST_ROOT=$BLAST_ROOT,ALIGN_ROOT=$ALIGN_ROOT,MANIFEST=$MANIFEST,DB_PREFIX=$DB_PREFIX,THREADS=$THREADS,MIN_PIDENT=$MIN_PIDENT,MIN_ALIGN_LEN=$MIN_ALIGN_LEN,MAX_EVALUE=$MAX_EVALUE,MAX_READS=$MAX_READS,YTAB_CONDA_ENV=$YTAB_CONDA_ENV"
ARRAY_SCRIPT="$REPO_ROOT/pipeline/hpc/zn_toxicity_screen/08_r2_hermes_blast_array.sge"
SUMMARY_SCRIPT="$REPO_ROOT/pipeline/hpc/zn_toxicity_screen/09_r2_hermes_blast_summary.sge"

array_cmd=(qsub -q "$QUEUE" -pe smp "$THREADS" -t "1-$SAMPLE_COUNT" -v "$COMMON_VARS" "$ARRAY_SCRIPT")
summary_cmd=(qsub -q "$QUEUE" -pe smp 1 -hold_jid "YTAB_Zn_R2_Hermes_BLAST" -v "$COMMON_VARS" "$SUMMARY_SCRIPT")

echo "Repo root:      $REPO_ROOT"
echo "R2 directory:   $R2_DIR"
echo "Sample sheet:   $SAMPLE_SHEET"
echo "Hermes CDS:     $HERMES_SOURCE"
echo "Manifest:       $MANIFEST"
echo "R2 FASTQs:      $SAMPLE_COUNT"
echo "Alignment root: $ALIGN_ROOT"
echo "Summary root:   $BLAST_ROOT/summary"
echo "Figure root:    $BLAST_ROOT/figures"
echo

if [[ "$DRY_RUN" == "1" ]]; then
  echo "DRY-RUN array qsub:"
  printf '  %q' "${array_cmd[@]}"; echo
  echo "DRY-RUN summary qsub:"
  printf '  %q' "${summary_cmd[@]}"; echo
  exit 0
fi

echo "Submitting R2 Hermes BLAST array..."
array_output="$("${array_cmd[@]}")"
echo "$array_output"
array_job_id="$(awk '{print $3}' <<<"$array_output")"

echo "Submitting held summary/plot job..."
summary_cmd=(qsub -q "$QUEUE" -pe smp 1 -hold_jid "$array_job_id" -v "$COMMON_VARS" "$SUMMARY_SCRIPT")
summary_output="$("${summary_cmd[@]}")"
echo "$summary_output"
summary_job_id="$(awk '{print $3}' <<<"$summary_output")"

echo
echo "Submitted jobs:"
echo "  BLAST array: $array_job_id"
echo "  Summary/plot: $summary_job_id"
echo
echo "After jobs finish:"
echo "  ls -lh $BLAST_ROOT/summary/"
echo "  ls -lh $BLAST_ROOT/figures/"
