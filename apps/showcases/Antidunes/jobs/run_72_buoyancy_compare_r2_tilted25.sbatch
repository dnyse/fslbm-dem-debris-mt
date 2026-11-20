#!/bin/bash -l
#SBATCH --job-name=antidunes-buoyancy-r2-compare-tilted25
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/buoyancy_compare_r2_tilted25_%j.out
#SBATCH --error=sbatch-output/error/buoyancy_compare_r2_tilted25_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# MEM vs. analytical buoyancy at R = 2: sphere, cylinder and box, each with MEM-only and analytical buoyancy,
# run back-to-back in one job so the allocation is requested only once.
# The box prms are the tilted25 variants (initialTiltAngleDeg 25, scripts/gen_buoyancy_tilted25_box_prm.py):
# upright, the box sits on a mapping tie (face planes through cell centres) and can map flush onto the floor with
# no fluid cell below it, which starves MEM-only buoyancy. The sphere and cylinder prms are upright.
# R1, R2 and R4 can be chained (--dependency=afterany) in one checkout that holds all three bed files
# spheres_out_buoyancy_r{1,2,4}.dat.

# R = 2: 72 ranks / 1 node, cellsPerBlock <32,32,24>; 1.5 h requested (~1 h 09 min measured for the 6 runs).
# --constraint=el9 is omitted on purpose (see run_72_bedgen_buoyancy_r2.sbatch).
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

run_case "R=2 sphere MEM"        buoyancy_compare_sphere_r2_mem_gas.prm        sphere_mem
run_case "R=2 sphere analytic"   buoyancy_compare_sphere_r2_analytic_gas.prm   sphere_analytic
run_case "R=2 cylinder MEM"      buoyancy_compare_cylinder_r2_mem_gas.prm      cylinder_mem
run_case "R=2 cylinder analytic" buoyancy_compare_cylinder_r2_analytic_gas.prm cylinder_analytic
run_case "R=2 box MEM tilted25"           buoyancy_compare_box_r2_mem_tilted25_gas.prm      box_mem
run_case "R=2 box analytic tilted25"      buoyancy_compare_box_r2_analytic_tilted25_gas.prm box_analytic

echo "## DONE. Archived under $ARCHIVE_ROOT/{sphere,cylinder,box}_{mem,analytic}_{vtk-out,eval,snapshots}"
