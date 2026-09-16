#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-debug-extreme-pidon-phase2
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/clog_debug_extreme_pidon_phase2_%j.out
#SBATCH --error=sbatch-output/error/clog_debug_extreme_pidon_phase2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3 (THESIS_PLAN.md) -- PHASE 2: dense-diagnostic resume from phase 1's own
# checkpoint, to determine whether phase 1.5's near-instant `reynoldsNumber is nan!` crash
# (job 4235729, ~96 steps after resume) is a checkpoint save/load bug or a genuine (if fast)
# post-restart numerical divergence. beginTimeStep is 49501 (phase 1's last checkpoint write,
# an off-by-one relative to its own timesteps=50000 -- see clog_debug_extreme_pidon_phase1_5_gas.prm's
# own header), and this prm's timesteps=50000, so this run is naturally SHORT (~499 steps) --
# cheap to run and analyze, not a repeat of phase 1's own full-length run.
#
# denseNanScanEnabled=true from step 0 (i.e. from the very first post-load step, t=49501) with
# both precursor thresholds (lowRho/highVelocity) and tight per-step VTK (writeFrequency=10,
# velocity+density+fill_level+mapped_flag) means: if the FIRST post-load step already shows a
# precursor warning or the crash fires within just 1-2 steps, that points at the loaded state
# itself being bad (checkpoint bug); if precursors build up gradually over tens of steps of
# genuine new computation before the crash, that points at a real post-restart numerical
# divergence instead (not a checkpoint bug, but a hypersensitive state).
#
# storeSnapshot=false in the prm -- this run must NOT overwrite phase 1's own still-intact
# checkpoint (snapshots/, timestamped 06:48 from phase 1, untouched by phase 1.5's crash since
# phase 1.5 never reached its own next scheduled checkpoint write at t=50000).
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
# NOTE: deliberately NOT touching snapshots/ -- this job loads phase 1's checkpoint and (with
# storeSnapshot=false in the prm) must not overwrite it.

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## clog debug: extreme, PID on -- PHASE 2 (dense-diagnostic resume, short window)"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_debug_extreme_pidon_phase2_gas.prm"

echo "## PHASE 2 finished -- check output for [LOW-RHO-PRECURSOR]/[HIGH-VELOCITY-PRECURSOR] onset step and tight VTK series."
