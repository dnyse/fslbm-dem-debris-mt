#!/bin/bash -l
#SBATCH --job-name=antidunes-droptest-sat24
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/droptest_sat24_%j.out
#SBATCH --error=sbatch-output/error/droptest_sat24_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# NOTE: --constraint=el9 removed 2026-09-11 -- sbatch was rejecting ANY job (including pre-existing,
# previously-working scripts) with "allocation failure: Invalid feature specification" while this
# constraint was present, a cluster-side scheduler issue unrelated to this script. Every node in the
# singlenode partition already carries the el9 feature (confirmed via sinfo), so omitting the
# constraint is a safe no-op for node selection, not a behavior change -- re-add it once the scheduler
# issue is resolved cluster-side, if ever needed again.
#
# CRASH-REPRODUCTION run (per explicit request): job 4216120 (36 processes, useSpheroBoxContact=true)
# crashed at step 6506 with an MPI-level error -- "[ParseMessage] Received invalid notification type:
# e from sender: 30" on rank 18 -- a corrupted inter-rank message, not a physics/logic abort. This
# job runs the SAME 1000-cuboid scenario at 24 processes (72/3, not 36) with useSpheroBoxContact=false
# (plain analytic boxBoxSAT instead of GJK/EPA-with-margin) to see whether the crash reproduces under
# a different process count and a different Box x Box contact algorithm -- see
# prms/DropTest_sat_24proc_gas.prm's own header for full context.
#
# Fluid-free, DEM-only cuboid drop test (DropTest.cpp) at real block-forest scale. 1000 cuboids, 10
# released every 1500 steps, run continues 10000 steps past the last release (timesteps=158500
# total). Time limit kept low (1h) since this is a small, fluid-free DEM-only run -- resubmit with a
# longer limit if it doesn't finish in time. Built directly via cmake/mpirun rather than
# antidunes.sh, since that script is hardcoded to the Antidunes/BedGeneration/tsunami target names
# and is left untouched per this session's own constraint (existing Antidunes files stay as-is).
# Run from a SEPARATE cluster checkout so it can't race any other run's build/output.
set -euo pipefail

# antidunes.sh normally handles this module loading internally -- since this script bypasses it
# (see header comment), the same modules it would have loaded must be loaded here explicitly.
module load gcc/11.5.0
module load openmpi/4.1.3-gcc11.5.0
module load cmake/3.31.8
module load python/3.13-miniforge

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_DIR="$REPO/build"
BIN_DIR="$BUILD_DIR/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

cmake --build "$BUILD_DIR" --target DropTest -- -j"$(nproc)"

cd "$BIN_DIR"
mkdir -p eval snapshots
mpirun -n 24 ./DropTest "$REPO/apps/showcases/Antidunes/prms/DropTest_sat_24proc_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_DIR"
fi
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
    fi
done
