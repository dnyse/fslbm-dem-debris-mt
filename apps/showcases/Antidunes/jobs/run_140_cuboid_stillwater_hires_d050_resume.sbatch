#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d050-resume
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=140
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d050_resume_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d050_resume_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-18] Thesis item 1a (THESIS_PLAN.md) -- resumes job 4240983's own t=90000
# checkpoint through to the original production target (timesteps=100000), as a normal completion
# run WITHOUT diagnostic verbosity, per explicit request ("run the rest of the crashed 0.5 density
# cuboid still water equilibrium run without the detailed output so the simulation is complete").
# See cuboid_stillwater_realgravity981_3xsphere_hires_d050_resume_gas.prm's own header for the full
# root-cause investigation this follows (7/7 diagnostic resume attempts came back clean, no confirmed
# mechanism for why -- this does NOT claim to have fixed anything, just completes the run).
#
# Only 10000 steps remain (90000->100000) vs. the original job's full 100000 -- 1h requested (vs. the
# original production job's 3h) for margin, matching this project's own convention of scaling time
# budgets to actual remaining work rather than reusing the full-run estimate.
#
# Deliberately does NOT touch snapshots/ at start (unlike a fresh-start job) -- this job LOADS the
# existing t=90000 checkpoint. storeSnapshot=true in the prm means this run WILL advance/overwrite
# snapshots/ as it progresses -- the original t=90000 checkpoint has already been backed up to
# $HPCVAULT/antidunes-archive/d050_t90000_checkpoint_backup_20260918/ before this job was submitted,
# so overwriting the live copy here is safe.
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

for dir in vtk-out eval; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done
# NOTE: deliberately NOT touching snapshots/ -- this job loads the existing t=90000 checkpoint.

if [[ ! -f "$BIN_DIR/snapshots/simState.file" ]]; then
    echo "FATAL: no snapshots/simState.file found in $BIN_DIR -- refusing to run with loadSnapshot=true" >&2
    exit 1
fi

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires, densityRatio=0.5 -- RESUME from t=90000 to completion (t=100000)"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d050_resume_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use), not at job start -- see
# run_140_cuboid_stillwater_hires_d050.sbatch's own 2026-09-16 comment for the full rationale
# (autofs staleness, local-only is a last resort not a default).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_d050_resume"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/d050_resume" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/d050_resume/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/d050_resume/{eval,vtk-out,snapshots}"
