#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-resume400-to1k-hotspot
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=0:25:00
#SBATCH --output=sbatch-output/output/bogner_resume400_to1k_hotspot_%j.out
#SBATCH --error=sbatch-output/error/bogner_resume400_to1k_hotspot_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Deep-dive resume: loads the step_400 checkpoint from the all-snapshots diagnostic
# run (job 4261093, staged at ~/mt-dennys-mirror2/apps/showcases/Antidunes/resume_staging_step400/,
# md5-verified after transfer) and reruns from t=401 with a widened VTK region covering the two
# hot-spot columns (X=8,Z=26) and (X=0,Z=25), Y=0-8, at writeFrequency=1 -- a full per-step field dump
# leading up to the (expected, deterministic) crash at step 444. See this run's own prm header for the
# full trace/rationale. Only ~44-50 steps expected before the crash; 20 min is generous margin.
# Queue: independent, confirmed clean/free before submission.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
STAGE_DIR="$REPO/apps/showcases/Antidunes/resume_staging_step400"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "## Staging step_400 checkpoint into $BIN_DIR/snapshots"
mkdir -p "$BIN_DIR/snapshots"
cp "$STAGE_DIR"/* "$BIN_DIR/snapshots/"
ls -la "$BIN_DIR/snapshots"

echo "############################################################"
echo "## Bogner resume from step 400, hot-spot region, writeFrequency=1"
echo "############################################################"
set +e
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_rho0767_tilt26p565_resume400_to1k_hotspot_gas.prm"
RC=$?
set -e
echo "## run exited with status $RC (a crash around step 444 is the EXPECTED outcome here)"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_resume400_to1k_hotspot"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
exit 0
