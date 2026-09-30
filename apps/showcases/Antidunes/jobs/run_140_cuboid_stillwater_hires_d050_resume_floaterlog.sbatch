#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d050-resume-floaterlog
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=140
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d050_resume_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d050_resume_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-18] Thesis item 1a (THESIS_PLAN.md) -- rotation-angle logging rerun of the
# already-completed d050 resume (job 4244926), per explicit request for a timestep-vs-rotation-angle
# graph. See cuboid_stillwater_realgravity981_3xsphere_hires_d050_resume_floaterlog_gas.prm's own
# header for why this rerun was needed (neither the completed run nor the original crashed run has
# any usable rotation data -- no db column, VTK never fired) and what it can/can't recover (only the
# t=90000-100000 tail, not the original 0-90000 history, which is genuinely lost).
#
# Read-only diagnostic rerun of an already-completed result -- storeSnapshot=false in the prm, this
# does not need to (and must not) produce a new checkpoint. Only 10000 steps, single floater ->
# logFloaterHistory is cheap here, nowhere near the many-floater/long-run scale that motivated that
# feature's own start/end-step filters.
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

for dir in vtk-out eval floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done
# NOTE: deliberately NOT touching snapshots/ -- this job loads the existing t=90000 checkpoint
# (restored from $HPCVAULT/antidunes-archive/d050_t90000_checkpoint_backup_20260918/ before this
# job was submitted).

if [[ ! -f "$BIN_DIR/snapshots/simState.file" ]]; then
    echo "FATAL: no snapshots/simState.file found in $BIN_DIR -- refusing to run with loadSnapshot=true" >&2
    exit 1
fi

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires, densityRatio=0.5 -- ROTATION LOGGING resume, t=90000 to 100000"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d050_resume_floaterlog_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use), not at job start -- see
# run_140_cuboid_stillwater_hires_d050.sbatch's own 2026-09-16 comment for the full rationale.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_d050_resume_floaterlog"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

for d in vtk-out eval floater-history; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/d050_resume_floaterlog" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/d050_resume_floaterlog/$d"
done
# snapshots/ deliberately left in place (storeSnapshot=false, nothing new to archive, and we want the
# restored t=90000 checkpoint to remain available for any further diagnosis without re-copying it).

echo "## DONE. Archived under $ARCHIVE_ROOT/d050_resume_floaterlog/{eval,vtk-out,floater-history}"
