#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-all
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_all_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3 (THESIS_PLAN.md) -- validates the new ClogMeasurementEvaluator (Evaluators.h) and the
# new PIDParameters.enabled constant-forcing flag (Forcing.h) on a cheap two-pillar sanity scenario,
# NOT yet item 4's real clogging scenarios. Runs all 3 prms back-to-back in ONE job (same run_case()
# combined-job pattern as item 1a/1b's own sbatch scripts) so the allocation is only requested once:
#   1. no-clog reference (no floating spheres at all)
#   2. sanity-clog, PID enabled (controlled-flume mode)
#   3. sanity-clog, PID disabled (constant-forcing/fixed-hydraulic-gradient mode)
# 72 processes (cellsPerBlock <64,20,15>, numBlocks <6,2,6>=72) -- matches this project's standard
# single-node rank count. An earlier version of this file used only 8 ranks (cellsPerBlock
# <48,40,90>, numBlocks <8,1,1>) on the theory that a small "sanity check" run should also request
# few cores -- wrong: --partition=singlenode --nodes=1 already reserves a whole 72-core node
# regardless of --ntasks, so 8 ranks just left 64 cores idle for zero benefit while making each rank
# own ~9x more cells (slower wall-clock, not faster). Fixed to 72 ranks here.
# Time limit 06:00:00: based on a REAL measured throughput from a comparable existing run
# (cuboid_stillwater_gravity_comparison_3xsphere_d025, job 4035436, archived log: domain
# <96,135,120>=1,555,200 cells / 72 ranks = 21600 cells/rank, steady-state 56.4ms/timestep --
# excluding one-time setup) rather than a hand extrapolation -- this scenario's ~19200 cells/rank is
# the same order of magnitude, so ~56500 steps * 56.4ms ~= 53 min/run * 3 runs ~= 2.6-2.9h, +~2x
# buffer for this scenario's extra pillar/BodyBC and 12-floater DEM overhead (not present in that
# reference run).
# [Claude, 2026-09-14]: prms' release ladder pushed from t=5000 to t=20000 and timesteps 35000->56500
# (see each prm's own header) -- flow wasn't spun up (meanVelocity_x still converging) when the first
# sphere used to arrive at t=5000; time limit bumped from 04:00:00 accordingly.
set -euo pipefail

# [Claude, 2026-09-14]: this REPO derivation previously used 4 "../" (assuming sbatch was invoked
# from inside jobs/), inconsistent with the rest of this project's scripts (e.g.
# run_72_two_cylinder_floater.sbatch), which are always submitted with `sbatch` from
# apps/showcases/Antidunes itself (same dir as antidunes.sh) -- so SLURM_SUBMIT_DIR sits only 3
# levels below the repo root. Fixed to 3 "../" to match; `sbatch jobs/run_....sbatch` invoked from
# Antidunes/ now works, consistent with every other script in this directory.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$1" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$1/$d"
    done
}

run_case() {  # $1 = banner, $2 = prm file, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "no-clog reference"        clog_measurement_validation_noclog_reference_gas.prm    noclog_reference
run_case "sanity-clog, PID on"      clog_measurement_validation_sanity_clog_pidon_gas.prm   sanity_clog_pidon
run_case "sanity-clog, PID off"     clog_measurement_validation_sanity_clog_pidoff_gas.prm  sanity_clog_pidoff

echo "## DONE. Archived under $ARCHIVE_ROOT/{noclog_reference,sanity_clog_pidon,sanity_clog_pidoff}_{vtk-out,eval,snapshots}"
