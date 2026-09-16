#!/bin/bash -l
#SBATCH --job-name=antidunes-droptest-fraige-250
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/droptest_fraige_250_%j.out
#SBATCH --error=sbatch-output/error/droptest_fraige_250_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# NOTE: --constraint=el9 omitted, see run_36_droptest_all.sbatch's own comment for why.
#
# Same as run_36_droptest_fraige.sbatch (SAT + boxBoxFraigeMultiPointContacts box-box contact,
# useSpheroBoxContact=false/useFraigeBoxContactModel=true), but with the floater count reduced
# 1000 -> 250 (prms/DropTest_fraige_250_gas.prm instead of DropTest_fraige_gas.prm, per explicit
# request 2026-09-13) -- same 10-column spawn grid truncated to the first 25 waves/column, timesteps
# reduced 158500 -> 46000 accordingly. --time reduced to 4h (generous guess for the shorter run, not
# a measured budget -- neither this reduced scale nor the full 1000-floater Fraige run has completed
# yet). See DropTest_fraige_250_gas.prm's own header, and DropTest_fraige_gas.prm's, for the full
# motivation (job 4226400's crash under the OLD useSpheroBoxContact=true contact model).
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
mpirun -n 36 ./DropTest "$PRM_DIR/DropTest_fraige_250_gas.prm"

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/cuboids_fraige_250_$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/cuboids_fraige_250_{vtk-out,eval,snapshots}"
