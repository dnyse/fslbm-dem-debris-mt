#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-81floaters-pass3-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_81floaters_pass3_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_81floaters_pass3_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# A/B comparison of PdfReconstructionManager's new Pass-3 local re-resolution
# sweep (usePass3Reresolution): same 81-sphere drop (plain useHardcodedUnits
# gravity, proven clean in job 3856318), 30000 timesteps, run TWICE in this
# one job -- once with Pass 3 off (sphere_only_81floaters_pass3off_30k_gas.prm,
# the baseline/default behaviour), once with it on
# (sphere_only_81floaters_pass3on_30k_gas.prm). The binary is built once (the
# toggle is a runtime prm flag, not compile-time); outputs from run A are
# archived before run B starts so their vtk-out/eval/snapshots don't collide
# (both use the same timestep range 0-30000, so file names WOULD clash
# otherwise). [FORMERMO-CONVERT] logs (evaluationFrequency 1000) in each run's
# .out give the per-branch FormerMO-conversion breakdown; run B's additionally
# report pass3ToLiquid/pass3ToGas/pass3StillInterface/pass3UnresolvedIsland --
# i.e. which sweep (Pass 2 directly vs. Pass 3) resolved each cell.
#
# set -e so a fatal in-simulation crash actually fails this SLURM job instead
# of being masked by the archiving step's own clean exit.
set -euo pipefail

# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

archive_outputs() {
    local label="$1"
    local stamp
    stamp="$(date +%Y%m%d_%H%M%S)"
    local archive_dir
    if [[ -n "${HPCVAULT:-}" ]]; then
        archive_dir="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${label}_${stamp}"
    else
        archive_dir="$BIN_DIR/archive/${label}_${stamp}"
    fi
    mkdir -p "$archive_dir"
    for dir in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$archive_dir/$dir"
    done
    echo "Archived $label outputs to $archive_dir"
}

# Build once -- shared binary for both runs.
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)"

# Run A: Pass 3 OFF (baseline).
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_only_81floaters_pass3off_30k_gas.prm"
archive_outputs "pass3off"

# Run B: Pass 3 ON.
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_only_81floaters_pass3on_30k_gas.prm"
archive_outputs "pass3on"
