#!/bin/bash -l
#SBATCH --job-name=antidunes-droptest-fraige
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/droptest_fraige_%j.out
#SBATCH --error=sbatch-output/error/droptest_fraige_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# NOTE: --constraint=el9 omitted, see run_36_droptest_all.sbatch's own comment for why.
#
# Fluid-free, DEM-only cuboid drop test (DropTest.cpp), SAME 1000-cuboid scenario as
# run_36_droptest_all.sbatch's own "cuboids" case, but using prms/DropTest_fraige_gas.prm instead of
# DropTest_gas.prm -- box-box contact routed through boxBoxSAT + the new
# paper_analytical::boxBoxFraigeMultiPointContacts (useFraigeBoxContactModel=true,
# useSpheroBoxContact=false) instead of the GJK/EPA-margin spherosimplex path
# (useSpheroBoxContact=true) job 4226400 used. See DropTest_fraige_gas.prm's own header for the full
# motivation: job 4226400 crashed at step 7928/158500 with "Received invalid notification type: ...
# from sender: 12" (mesa_pd/mpi/notifications/ParseMessage.h:273, a corrupted/desynchronized
# inter-rank particle-sync message) -- this run tests whether the different box-contact code path
# avoids or reproduces that same failure. Cylinders are NOT re-run here (the crash was in the cuboids
# phase, which never even got far enough to hand off to cylinders in job 4226400) -- run
# run_36_droptest_all.sbatch again separately if the cylinders phase itself needs re-checking.
set -euo pipefail

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
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

cmake --build "$BUILD_DIR" --target DropTest -- -j"$(nproc)"

cd "$BIN_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

mkdir -p eval snapshots
mpirun -n 36 ./DropTest "$PRM_DIR/DropTest_fraige_gas.prm"

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/cuboids_fraige_$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/cuboids_fraige_{vtk-out,eval,snapshots}"
