#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-resume400-full50k
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=1:00:00
#SBATCH --output=sbatch-output/output/bogner_resume400_full50k_%j.out
#SBATCH --error=sbatch-output/error/bogner_resume400_full50k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-24] Continuation of the 26.565deg-seed rho=0.7665 box from the step_400 checkpoint
# (same one staged for job 4261250's hot-spot diagnostic, reused here), now running the FULL remaining
# ~49,600 steps with NORMAL logging (no per-step keepAllSnapshots, no widened VTK region, no
# writeFrequency=1) -- see this run's own prm header and THESIS_PLAN.md for why: the diagnostic run
# did NOT reproduce the deterministic step-444 crash, so this continuation checks where the box
# actually settles once that crash is avoided.
# Timing: ~49,600 remaining steps at this box's own observed rate (100k steps ~62min => ~31min for 50k)
# -> ~30 min expected + archive overhead; 1h limit for margin.
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
echo "## Bogner resume from step 400, FULL 50k steps, normal logging"
echo "############################################################"
"$SH" --run --np 65 --prm "$PRM_DIR/bogner_exact_3box_s16_rho0767_tilt26p565_resume400_full50k_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_resume400_full50k"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
