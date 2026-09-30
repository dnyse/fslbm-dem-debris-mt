#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d050-diag-repeat
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=140
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d050_diag_repeat_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d050_diag_repeat_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-17] Thesis item 1a (THESIS_PLAN.md) -- repeats the single-shot diagnostic
# (run_140_cuboid_stillwater_hires_d050_diag.sbatch) N times in a row, ALL from the SAME job 4240983
# checkpoint (t=90000), to test reproducibility of the quaternion-normalization crash
# (Quaternion.h:241, step ~91960) after job 4242078's single dense-scan attempt came back completely
# clean (zero precursor firings, zero NaN, full completion to step 92500) -- the opposite of the
# original production run's own behaviour from the identical checkpoint.
#
# IMPORTANT CAVEAT (found while investigating why 4242078 didn't reproduce it): checkpoint resume is
# NOT a bit-for-bit re-creation of the original run's own microstate at t=90000. Antidunes.cpp:110
# unconditionally calls initFlagsFromFillLevel() on every setup, fresh or resumed -- this does a
# domain-wide "clear all liquid/interface/gas flags, then reclassify purely from the current fill
# field" pass (src/lbm/free_surface/InitFunctions.h:74-117), using handling->forceFlag()/removeFlag()
# directly rather than the setDomain()/setBoundary() API that BoundaryHandling.h's own
# getNearBoundaryFlag() comment says is the ONLY sanctioned way to touch the near-boundary flag
# ("Never (ever) set near boundary flags manually outside the boundary handler!", BoundaryHandling.h:105).
# This is exactly why 4242078's log showed a BURST of 302 [STALE-NEARBOUNDARY] corrections clustered
# right at step 90002 (setup time) rather than the trickle a continuously-running simulation would
# generate as cells convert one at a time. The existing mitigation (AntidunesBoundaryHandling.h:180ish)
# correctly clears every case it catches -- no NaN resulted -- but the resumed run's exact starting
# flag-field microstate is demonstrably NOT identical to what continuous execution would have had at
# that same instant. So this repeat test can only show whether the CHECKPOINT-RESUME path itself is
# flaky across repeated attempts (same nodes, same checkpoint, same reclassification pass each time)
# -- it cannot prove the original continuous-execution bug is gone, since it never recreates the exact
# conditions that triggered it.
#
# All repeats run on the SAME node allocation (fixed for the whole job) -- if MPI/LBM execution is
# bit-for-bit deterministic on fixed hardware (the expected case), all N repeats should produce an
# IDENTICAL outcome (all clean or all crash); divergence between repeats within this one job would
# itself be a notable finding (real non-determinism beyond node-to-node rounding differences). A
# genuine test of "does node allocation change the outcome" would need separate job submissions
# instead -- left as a natural follow-up if this comes back "always identical".
#
# 140 ranks (matches the checkpoint's own block layout, same as the single-shot diag). Each repeat is
# ~2500 steps at ~3.5min (per 4242078's own measured timing) -- 6 repeats budgeted, 2h requested for
# margin (build + 6x runs + summary).
set -uo pipefail   # NOTE: no -e -- a crashing repeat must not abort the remaining repeats

NUM_REPEATS=6

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
LOG_DIR="$BIN_DIR/diag_repeat_logs"
mkdir -p "$LOG_DIR"

if [[ ! -f "$BIN_DIR/snapshots/simState.file" ]]; then
    echo "FATAL: no snapshots/simState.file found in $BIN_DIR -- refusing to run with loadSnapshot=true" >&2
    exit 1
fi

set -e
"$SH" --build --jobs "$(nproc)"
set +e

declare -a RESULTS
for i in $(seq 1 "$NUM_REPEATS"); do
    echo "############################################################"
    echo "## REPEAT $i/$NUM_REPEATS -- cuboid stillwater hires d050, resume from t=90000 checkpoint"
    echo "############################################################"

    # Clear leftover per-run outputs (VTK/eval) between repeats -- NEVER touch snapshots/, every
    # repeat reads the SAME job-4240983 checkpoint (storeSnapshot=false in the diag prm).
    for dir in vtk-out eval; do
        [[ -d "$BIN_DIR/$dir" ]] && rm -rf "${BIN_DIR:?}/$dir"
    done

    LOGFILE="$LOG_DIR/repeat_${i}.log"
    "$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d050_diag_gas.prm" \
        2>&1 | tee "$LOGFILE"
    RC=${PIPESTATUS[0]}

    if [[ "$RC" -eq 0 ]] && ! grep -qE 'DRY-LINK|LOW-RHO|quaternion|Assertion|WALBERLA_ABORT|non-finite' "$LOGFILE"; then
        RESULTS[$i]="CLEAN"
    else
        RESULTS[$i]="CRASHED (exit=$RC)"
    fi
    echo "## REPEAT $i result: ${RESULTS[$i]}"
done

echo "############################################################"
echo "## SUMMARY: $NUM_REPEATS repeats from the same t=90000 checkpoint"
echo "############################################################"
for i in $(seq 1 "$NUM_REPEATS"); do
    echo "  repeat $i: ${RESULTS[$i]}  (log: $LOG_DIR/repeat_${i}.log)"
done
