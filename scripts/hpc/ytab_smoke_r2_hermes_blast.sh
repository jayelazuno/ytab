#!/bin/bash
# Smoke-test the Zn R2 Hermes BLAST diagnostic without submitting the full SGE run.

set -euo pipefail

REPO_ROOT="/old_Users/jayelazuno/workspace/ytab"
R2_DIR="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110"
SAMPLE_SHEET="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110/sampleSheetUsed.csv"
BLAST_ROOT="/old_Users/jayelazuno/workspace/ytab/blast"
ALIGN_ROOT="/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/ytab_work/Zn_toxicity_screen/blast_alignments"
HERMES_CDS="/old_Users/jayelazuno/workspace/ytab/blast/Hermes-cds.txt"
YTAB_CONDA_ENV="${YTAB_CONDA_ENV:-ytab-hermes-blast}"
MAX_READS="10000"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/hpc/ytab_smoke_r2_hermes_blast.sh [options]

Options:
  --repo-root PATH
  --r2-dir PATH
  --sample-sheet PATH
  --blast-root PATH
  --align-root PATH
  --hermes-cds PATH
  --conda-env ytab-hermes-blast
  --max-reads 10000
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
    --max-reads) MAX_READS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

cd "$REPO_ROOT"

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

for tool in blastn makeblastdb seqkit python Rscript; do
  command -v "$tool" >/dev/null 2>&1 || { echo "FAIL: missing required tool: $tool" >&2; exit 1; }
done

stamp="$(mktemp)"
mkdir -p \
  "$BLAST_ROOT/reference" "$BLAST_ROOT/db" "$BLAST_ROOT/manifest" "$BLAST_ROOT/summary/per_sample" \
  "$BLAST_ROOT/figures" "$BLAST_ROOT/logs" "$BLAST_ROOT/scripts_metadata" \
  "$ALIGN_ROOT/fasta_queries" "$ALIGN_ROOT/blast_outfmt6" "$ALIGN_ROOT/unique_hit_reads" "$ALIGN_ROOT/seqkit_stats"

bash "$REPO_ROOT/scripts/hpc/ytab_submit_zn_r2_hermes_blast.sh" \
  --repo-root "$REPO_ROOT" \
  --r2-dir "$R2_DIR" \
  --sample-sheet "$SAMPLE_SHEET" \
  --blast-root "$BLAST_ROOT" \
  --align-root "$ALIGN_ROOT" \
  --hermes-cds "$HERMES_CDS" \
  --conda-env "$YTAB_CONDA_ENV" \
  --max-reads "$MAX_READS" \
  --smoke \
  --dry-run >/dev/null

hermes_source=""
for candidate in "$HERMES_CDS" "$BLAST_ROOT/Hermes-cds.txt" "$REPO_ROOT/docs/codex/Hermes-cds.txt"; do
  if [[ -s "$candidate" ]]; then hermes_source="$candidate"; break; fi
done
if [[ -z "$hermes_source" ]]; then
  echo "FAIL: Hermes CDS file could not be found" >&2
  exit 1
fi

python - "$hermes_source" "$BLAST_ROOT/reference/Hermes-cds.fa" <<'PY'
from pathlib import Path
import re
import sys
source = Path(sys.argv[1])
dest = Path(sys.argv[2])
text = source.read_text(encoding="utf-8", errors="replace").strip()
if not text.startswith(">"):
    raise SystemExit("FAIL: Hermes CDS is not FASTA")
seq = "".join(line.strip() for line in text.splitlines() if line and not line.startswith(">"))
letters = re.sub(r"[^A-Za-z]", "", seq).upper()
if not letters or sum(c in "ACGTUN" for c in letters) / len(letters) < 0.95:
    raise SystemExit("FAIL: Hermes CDS is not nucleotide-like")
dest.write_text(text + "\n", encoding="utf-8")
PY

makeblastdb -in "$BLAST_ROOT/reference/Hermes-cds.fa" -dbtype nucl -parse_seqids -out "$BLAST_ROOT/db/hermes_cds" >/dev/null

manifest="$BLAST_ROOT/manifest/r2_fastq_manifest.tsv"
python "$REPO_ROOT/scripts/hpc/ytab_make_r2_hermes_blast_manifest.py" \
  --r2-dir "$R2_DIR" \
  --sample-sheet "$SAMPLE_SHEET" \
  --output "$manifest" >/dev/null
row_count="$(python - "$manifest" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    print(sum(1 for _ in csv.DictReader(handle, delimiter="\t")))
PY
)"
if [[ "$row_count" -lt 1 ]]; then
  echo "FAIL: manifest has no R2 FASTQ rows" >&2
  exit 1
fi

SGE_TASK_ID=1 \
REPO_ROOT="$REPO_ROOT" \
BLAST_ROOT="$BLAST_ROOT" \
ALIGN_ROOT="$ALIGN_ROOT" \
MANIFEST="$manifest" \
DB_PREFIX="$BLAST_ROOT/db/hermes_cds" \
MAX_READS="$MAX_READS" \
YTAB_SKIP_CONDA=1 \
bash "$REPO_ROOT/pipeline/hpc/zn_toxicity_screen/08_r2_hermes_blast_array.sge" >/dev/null

sample_summary_count="$(find "$BLAST_ROOT/summary/per_sample" -name '*.r2_hermes_blast_summary.tsv' | wc -l | tr -d ' ')"
if [[ "$sample_summary_count" -lt 1 ]]; then
  echo "FAIL: per-sample BLAST summary was not created" >&2
  exit 1
fi

python "$REPO_ROOT/scripts/hpc/ytab_summarize_r2_hermes_blast.py" \
  --per-sample-dir "$BLAST_ROOT/summary/per_sample" \
  --summary-dir "$BLAST_ROOT/summary" >/dev/null

Rscript "$REPO_ROOT/scripts/hpc/ytab_plot_r2_hermes_blast_summary.R" \
  --summary "$BLAST_ROOT/summary/r2_hermes_blast_summary.csv" \
  --figures-dir "$BLAST_ROOT/figures" >/dev/null

for path in \
  "$BLAST_ROOT/summary/r2_hermes_blast_summary.csv" \
  "$BLAST_ROOT/summary/r2_hermes_blast_summary.tsv" \
  "$BLAST_ROOT/summary/r2_hermes_blast_compact.csv" \
  "$BLAST_ROOT/summary/r2_hermes_blast_summary.txt" \
  "$BLAST_ROOT/figures/r2_hermes_percent_by_sample.png" \
  "$BLAST_ROOT/figures/r2_hermes_percent_by_sample.pdf" \
  "$BLAST_ROOT/figures/r2_hermes_percent_by_role.png" \
  "$BLAST_ROOT/figures/r2_hermes_percent_by_role.pdf"; do
  if [[ ! -s "$path" ]]; then
    echo "FAIL: expected smoke output missing: $path" >&2
    exit 1
  fi
done

if [[ -d "$REPO_ROOT/output/projects" ]] && find "$REPO_ROOT/output/projects" "$REPO_ROOT/output/exports" -type f -newer "$stamp" 2>/dev/null | grep -q .; then
  echo "FAIL: smoke touched output/projects or output/exports" >&2
  exit 1
fi

rm -f "$stamp"
echo "PASS"
