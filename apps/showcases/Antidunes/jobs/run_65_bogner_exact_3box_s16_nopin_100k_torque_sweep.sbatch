#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-exact-s16-nopin-100k-torque-sweep
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=3:00:00
#SBATCH --output=sbatch-output/output/bogner_exact_3box_s16_nopin_100k_torque_sweep_%j.out
#SBATCH --error=sbatch-output/error/bogner_exact_3box_s16_nopin_100k_torque_sweep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-22] Thesis item 1a (THESIS_PLAN.md): follow-up to job 4255587 (the 50000-step nopin damping sweep, which found the FINAL
# resting magnitude ~30-32 deg was the SAME at every damping level tested -- damping only changes escape TIMING, not where it lands). Per
# explicit instruction, SIX runs this time, all rho=0.75, fully free (--nopin), 100000 steps (2x longer -- does it get closer to 26.565 deg
# given more time?), and evaluationFrequency 1 (--evalfreq 1) so [BoxForces]/[BoxTorque] log EVERY timestep instead of every 500 (BoxTorque
# already existed in ParticleKernels.h -- hydroTorque, buoyLeverTorque_y, totalTorque -- just throttled; no new code needed, just this prm
# key). Runs:
#   1-3. omega 1.9/1.5/1.0, seed 5 deg (same as 4255587, just longer + per-step torque)
#   4-5. omega 1.9 (Bogner's own damping), seed 3 deg and seed 4 deg (NEW: does a SMALLER seed -> smaller initial torque -> no overshoot ->
#        land AT 26.565 deg instead of ~30-32?)
#   6.   omega 1.9, seed AT the equilibrium heave (26.565 deg, z0 16.9223, same start as the PINNED run 4253988) but this time FULLY FREE
#        (nopin) -- 4253988 was position/rotation-pinned (roll-only) and stayed flat near 26.8-27 deg; this tests whether the FREE box also
#        stays near equilibrium (torque ~0, genuinely stable) or drifts away toward the same ~30-32 deg band the seeded-5-deg runs reach
#        (torque nonzero but too weak to move it within 4253988's shorter 50000-step window) -- directly answers "why doesn't it turn away
#        from equilibrium" with per-step torque data instead of just tiltY.
# All SIX in ONE job/ONE directory (matches 4255587's own precedent), sequential, build once, archived to per-run subfolders, continues past
# a crash. Time budget: 4255587's own 3 runs at 50000 steps took 34m58s total (~11-12 min each); 100000 steps ~2x that -> ~23-24 min each x 6
# = ~140 min + build/archive overhead -- 3h limit for margin.
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

"$SH" --build --jobs "$(nproc)"

STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_exact_s16_nopin_100k_torque_sweep"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

OVERALL_RC=0

run_one() {
    local tag="$1" prm="$2"

    local RESCUE="$BIN_DIR/archive/rescued_before_${tag}_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
    for dir in vtk-out eval snapshots; do
        if [[ -d "$BIN_DIR/$dir" ]]; then
            mkdir -p "$RESCUE"
            echo "WARNING: leftover $BIN_DIR/$dir before $tag -- moving to $RESCUE (NOT deleting)" >&2
            mv "$BIN_DIR/$dir" "$RESCUE/$dir"
        fi
    done

    echo "############################################################"
    echo "## Bogner exact box equilibrium, s=16 nopin 100k torque, $tag"
    echo "############################################################"
    set +e
    "$SH" --run --np 65 --prm "$prm"
    local rc=$?
    set -e
    echo "## $tag exited with status $rc -- archiving whatever output exists"
    [[ $rc -ne 0 ]] && OVERALL_RC=1

    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$tag" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$tag/$d"
    done
    echo "## DONE $tag. Archived under $ARCHIVE_ROOT/$tag/{eval,vtk-out,snapshots}"
}

run_one omega1p9_seed5 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p9_100k_torque_gas.prm"
run_one omega1p5_seed5 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p5_100k_torque_gas.prm"
run_one omega1p0_seed5 "$PRM_DIR/bogner_exact_3box_s16_seed5_rho075_nopin_omega1p0_100k_torque_gas.prm"
run_one omega1p9_seed3 "$PRM_DIR/bogner_exact_3box_s16_seed3_rho075_nopin_omega1p9_100k_torque_gas.prm"
run_one omega1p9_seed4 "$PRM_DIR/bogner_exact_3box_s16_seed4_rho075_nopin_omega1p9_100k_torque_gas.prm"
run_one omega1p9_seedEq "$PRM_DIR/bogner_exact_3box_s16_seedEq_rho075_nopin_omega1p9_100k_torque_gas.prm"

echo "## ALL SIX DONE. Archive root: $ARCHIVE_ROOT"
exit "$OVERALL_RC"
