#!/bin/bash -l
#SBATCH --job-name=antidunes-droptest-250-three-shapes
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/droptest_250_three_shapes_%j.out
#SBATCH --error=sbatch-output/error/droptest_250_three_shapes_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Fluid-free 250-floater drop tests (DropTest), one case per job, each in its own archive directory so outputs
# cannot mix. The three main cases share container, release positions and timing:
#   cylinders  DropTest_250_cylinders_gas.prm  250 cylinders (R 1, L 12), closed form with GJK/EPA fallback
#   cuboids    DropTest_250_sat_gas.prm        250 cuboids (6 x 1 x 1.5 half-extents), plain SAT
#   mixed      DropTest_250_mixed_gas.prm      84 cuboids, 83 cylinders, 83 spheres (R 2), cycled in release order
# Usage, from apps/showcases/Antidunes: build once in its own job (so the run jobs never build concurrently in
# the same build directory), then one job per case:
#   b=$(sbatch --parsable --time=00:30:00 --ntasks=8 jobs/run_36_droptest_250_three_shapes.sbatch build)
#   for c in cylinders cuboids mixed; do sbatch --dependency=afterok:$b --job-name=antidunes-droptest-250-$c \
#       jobs/run_36_droptest_250_three_shapes.sbatch $c; done
# mpirun uses $SLURM_NTASKS, so pass each case's process count with --ntasks: 36 by default, 24 for mixed
# (2x2x6 blocks); the other cases' counts are noted below. The cuboid case took ~55 min on 36 ranks.
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
BIN="$BUILD_DIR/apps/showcases/Antidunes/DropTest"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"

CASE="${1:?usage: sbatch run_36_droptest_250_three_shapes.sbatch build|cylinders|cuboids|mixed}"
if [[ "$CASE" == build ]]; then
    cmake --build "$BUILD_DIR" --target DropTest -- -j"$(nproc)"
    exit 0
fi
[[ -x "$BIN" ]] || { echo "DropTest not built: $BIN" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BUILD_DIR/apps/showcases/Antidunes/archive/droptest_250_three_shapes_$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi
echo "## Archive: $ARCHIVE_ROOT"

run_case() {  # $1 = label, $2 = prm file
    echo "############################################################"
    echo "## DropTest 250: $1 ($2), start $(date +%T)"
    echo "############################################################"
    local dir="$ARCHIVE_ROOT/$1"
    mkdir -p "$dir/eval" "$dir/snapshots"
    cp "$PRM_DIR/$2" "$dir/"
    ( cd "$dir" && mpirun -n "${SLURM_NTASKS:-36}" "$BIN" "$2" > run.log 2>&1 ) && status=0 || status=$?
    echo "## $1 finished $(date +%T) with exit $status"
    return $status
}

case "$CASE" in
    cylinders) run_case cylinders DropTest_250_cylinders_gas.prm ;;
    cuboids)   run_case cuboids   DropTest_250_sat_gas.prm ;;
    mixed)     run_case mixed     DropTest_250_mixed_gas.prm ;;
    # randomised orientations/offsets, wood-like contact (prms/gen_droptest_250_random.py); 48 tasks
    random_cylinders) run_case random_cylinders DropTest_250_random_cylinders_gas.prm ;;
    random_cuboids)   run_case random_cuboids   DropTest_250_random_cuboids_gas.prm ;;
    random_mixed)     run_case random_mixed     DropTest_250_random_mixed_gas.prm ;;
    # Kodam et al. (2010) packing in a round container, one-by-one release (scripts/gen_kodam_round_prm.py); 16 tasks
    kodam_round)      run_case kodam_round      DropTest_kodam_round_gas.prm ;;
    # the same with kodamBandEdgeAxialCheck true (checks the energy spike of kodam_round); 16 tasks
    kodam_round_axial) run_case kodam_round_axial DropTest_kodam_round_axial_gas.prm ;;
    # Kodam packing comparison, the four contact-detection variants (9 tasks each), run on the same code so that
    # their wall-clock times are comparable
    kodam_fallback)   run_case kodam_fallback   DropTest_kodam_packing_gas.prm ;;
    kodam_gjk)        run_case kodam_gjk        DropTest_kodam_packing_gjk_gas.prm ;;
    kodam_pure)       run_case kodam_pure       DropTest_kodam_packing_pureanalytical_gas.prm ;;
    kodam_refined)    run_case kodam_refined    DropTest_kodam_packing_refined_gas.prm ;;
    *) echo "unknown case: $CASE" >&2; exit 1 ;;
esac
echo "## DONE. Output under $ARCHIVE_ROOT/$CASE"
