#!/bin/bash
# Smoke-test Zn_toxicity_screen HPC wrapper and qsub scripts.

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd -P)}"
SCRIPT_DIR="$REPO_ROOT/scripts/hpc"
HPC_DIR="$REPO_ROOT/pipeline/hpc/zn_toxicity_screen"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ytab_zn_hpc_smoke.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

required=(
  "$SCRIPT_DIR/ytab_make_hpc_sample_manifest.py"
  "$SCRIPT_DIR/ytab_submit_zn_toxicity_screen_hpc.sh"
  "$SCRIPT_DIR/ytab_hpc_manifest_row.py"
  "$SCRIPT_DIR/ytab_hpc_config_value.py"
  "$SCRIPT_DIR/ytab_hpc_write_sample_status.py"
  "$HPC_DIR/01_mapfastq_array.sge"
  "$HPC_DIR/02_create_hit_file_array.sge"
  "$HPC_DIR/03_summary_raw.sge"
  "$HPC_DIR/04_library_diagnostics.sge"
  "$HPC_DIR/05_classifier_branch.sge"
  "$HPC_DIR/06_fitness_zn_vs_mock.sge"
  "$HPC_DIR/07_report_export_finalize.sge"
)

for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "FAIL: missing required script: $path"; exit 1; }
done

for path in "$SCRIPT_DIR"/*.sh "$HPC_DIR"/*.sge; do
  bash -n "$path"
done
python -B -c '
import pathlib
files = list(pathlib.Path("'"$SCRIPT_DIR"'").glob("ytab_*zn*.py")) + list(pathlib.Path("'"$SCRIPT_DIR"'").glob("ytab_hpc_*.py"))
for path in files:
    compile(path.read_text(encoding="utf-8"), str(path), "exec")
'

bad_project="EO46""_Zn""_tox"
bad_h2o2_project="H2O2""_screen_v1"
if grep -R "$bad_h2o2_project\|yH298-H2O2\|yH299-H2O2\|${bad_project}_v1\|$bad_project" \
  "$SCRIPT_DIR/ytab_make_hpc_sample_manifest.py" \
  "$SCRIPT_DIR/ytab_submit_zn_toxicity_screen_hpc.sh" \
  "$SCRIPT_DIR/ytab_hpc_manifest_row.py" \
  "$SCRIPT_DIR/ytab_hpc_config_value.py" \
  "$SCRIPT_DIR/ytab_hpc_write_sample_status.py" \
  "$HPC_DIR" >/dev/null; then
  echo "FAIL: forbidden H2O2 or legacy Zn project text found in HPC scripts"
  exit 1
fi

grep -q '#\$ -t 1-8' "$HPC_DIR/01_mapfastq_array.sge" || { echo "FAIL: MapFastq is not an array job"; exit 1; }
grep -q '#\$ -t 1-8' "$HPC_DIR/02_create_hit_file_array.sge" || { echo "FAIL: CreateHitFile is not an array job"; exit 1; }
grep -q 'SGE_TASK_ID' "$HPC_DIR/01_mapfastq_array.sge" || { echo "FAIL: MapFastq does not use SGE_TASK_ID"; exit 1; }
grep -q 'SGE_TASK_ID' "$HPC_DIR/02_create_hit_file_array.sge" || { echo "FAIL: CreateHitFile does not use SGE_TASK_ID"; exit 1; }
grep -q 'hpc_sample_manifest.csv' "$SCRIPT_DIR/ytab_submit_zn_toxicity_screen_hpc.sh" || { echo "FAIL: submit wrapper does not reference hpc manifest"; exit 1; }
grep -q -- '--r1' "$HPC_DIR/01_mapfastq_array.sge" || { echo "FAIL: MapFastq R1 not passed"; exit 1; }
grep -q -- '--r2' "$HPC_DIR/01_mapfastq_array.sge" || { echo "FAIL: MapFastq R2 not passed"; exit 1; }
grep -q 'scratch BAM' "$HPC_DIR/02_create_hit_file_array.sge" || { echo "FAIL: CreateHitFile does not document scratch BAM input"; exit 1; }
grep -q 'ytab_run_sample_normalization.py' "$HPC_DIR/05_classifier_branch.sge" || { echo "FAIL: classifier branch missing normalization"; exit 1; }
grep -q 'sample-mode parents' "$HPC_DIR/05_classifier_branch.sge" || { echo "FAIL: classifier branch does not use mock/control-only parent mode"; exit 1; }
grep -q 'Zn_1_5mM_vs_mock' "$HPC_DIR/06_fitness_zn_vs_mock.sge" || { echo "FAIL: fitness branch analysis ID missing"; exit 1; }
grep -q 'ytab_export_project.py' "$HPC_DIR/07_report_export_finalize.sge" || { echo "FAIL: final export job missing"; exit 1; }
grep -q -- '-hold_jid' "$SCRIPT_DIR/ytab_submit_zn_toxicity_screen_hpc.sh" || { echo "FAIL: hold_jid dependency support missing"; exit 1; }

bash "$SCRIPT_DIR/ytab_submit_zn_toxicity_screen_hpc.sh" \
  --dry-run \
  --repo-root "$TMP_DIR/repo" \
  --fastq-dir "$TMP_DIR/fastqs" \
  --scratch-work-dir "$TMP_DIR/scratch" \
  --sample-sheet "$TMP_DIR/fastqs/sampleSheetUsed.csv" \
  > "$TMP_DIR/dry_run.out" 2>&1

grep -q 'DRY-RUN YTAB_Zn_MapFastq' "$TMP_DIR/dry_run.out" || { echo "FAIL: dry-run MapFastq qsub not printed"; exit 1; }
grep -q 'DRY-RUN YTAB_Zn_CreateHitFile' "$TMP_DIR/dry_run.out" || { echo "FAIL: dry-run CreateHitFile qsub not printed"; exit 1; }
grep -q 'YTAB_Zn_ClassifierBranch' "$TMP_DIR/dry_run.out" || { echo "FAIL: dry-run classifier dependency not printed"; exit 1; }
grep -q 'YTAB_Zn_Fitness' "$TMP_DIR/dry_run.out" || { echo "FAIL: dry-run fitness dependency not printed"; exit 1; }
grep -q 'rsync -avP' "$TMP_DIR/dry_run.out" || { echo "FAIL: rsync pull instructions missing"; exit 1; }

echo "PASS"
