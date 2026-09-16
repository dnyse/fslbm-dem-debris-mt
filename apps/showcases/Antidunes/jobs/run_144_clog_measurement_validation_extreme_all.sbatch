#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-extreme-all
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=16:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_extreme_all_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_extreme_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3 (THESIS_PLAN.md) -- the "make it more extreme" follow-up to
# run_72_clog_measurement_validation_all.sbatch's small sanity-clog trio, requested once that trio's
# tool validation and flow-spinup-timing fix (release ladder pushed from t=5000 to t=20000) were
# confirmed correct. Same run_case() combined-job pattern, 3 prms back-to-back in ONE job:
#   1. no-clog reference, extreme scale  (clog_measurement_validation_noclog_reference_extreme_gas.prm)
#   2. extreme clog, PID enabled          (clog_measurement_validation_extreme_pidon_gas.prm)
#   3. extreme clog, PID disabled         (clog_measurement_validation_extreme_pidoff_gas.prm)
#
# Scale-up vs. the small trio: domain X doubled (384->768, every pillar/probe/release-zone X position
# is simply the small scenario's own X doubled -- a clean geometric scale-up, not a redesign; Y/Z and
# the pier gap itself are UNCHANGED, so the "one sphere physically cannot pass" guarantee still holds),
# and floater count tripled (12->36, same 12-entry X/density/Z cycle repeated 3x, same release ladder
# convention, just continued for 3x as many entries). See each new prm's own header for the exact
# derivation of every scaled number.
#
# 144 ranks (2x the small trio's 72): cellsPerBlock <64,20,15> unchanged, numBlocks
# (768/64, 40/20, 90/15) = (12, 2, 6) = 144 -- needs 2 nodes, hence --partition=multinode (the small
# trio's --partition=singlenode only ever reserves 1 node/72 cores, not enough here).
#
# Time limit 16:00:00: this scenario's own small-trio sibling measured ~19200 cells/rank at 72 ranks;
# this one is ~19200 cells/rank at 144 ranks too (same cellsPerBlock, so identical per-rank cell count
# by construction) -- so per-step cost should be similar to the small trio's own ~56ms/step baseline,
# but timesteps are 97500 (1.73x the small trio's 56500) AND floater count is 3x (more DEM/contact
# resolution cost per step, not accounted for in the per-cell throughput estimate) AND this is the
# project's first run at 2-node scale for this scenario (untested inter-node MPI overhead for the
# particle/DEM sync). Given that much extra uncertainty stacked on top of the small trio's own already
#-applied ~2x buffer, doubling further to a generous 16h ceiling (well under the multinode partition's
# 24h cap) rather than trying to extrapolate a tight number from three unverified multipliers at once.
# Adjust down once real per-step timing is known from this run's own log.
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
    "$SH" --run --np 144 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "no-clog reference, extreme"   clog_measurement_validation_noclog_reference_extreme_gas.prm  noclog_reference_extreme
run_case "extreme clog, PID on"         clog_measurement_validation_extreme_pidon_gas.prm              extreme_pidon
run_case "extreme clog, PID off"        clog_measurement_validation_extreme_pidoff_gas.prm             extreme_pidoff

echo "## DONE. Archived under $ARCHIVE_ROOT/{noclog_reference_extreme,extreme_pidon,extreme_pidoff}_{vtk-out,eval,snapshots}"
