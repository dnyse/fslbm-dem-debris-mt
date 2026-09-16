#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-calibration-r1-sphere-mem
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/buoyancy_calibration_r1_sphere_mem_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_calibration_r1_sphere_mem_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1b (THESIS_PLAN.md) -- calibration run to measure the sphere's actual rise rate before
# fixing `timesteps` for the real 18-job comparison sweep. See the prm's own header for full context.
# --constraint=el9 deliberately omitted (see run_72_bedgen_buoyancy_r1.sbatch's own header comment).
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

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)" --run --np 72 --prm "$PRM_DIR/buoyancy_calibration_r1_sphere_mem_gas.prm"
