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

# [Claude, 2026-10-01] Rerun of the fluid-free 250-floater drop tests for thesis Section 5.2 (author's request), on the
# current code (DropTest no-fluid fixes, sphere VTK output). Three cases, one after another, same container, release
# positions and timing (prms/DropTest_250_*_gas.prm, all derived from DropTest_250_sat_gas.prm):
#   cylinders  DropTest_250_cylinders_gas.prm  250 cylinders (R 1, L 12), closed form with GJK/EPA fallback
#   cuboids    DropTest_250_sat_gas.prm        250 cuboids (6 x 1 x 1.5 half-extents), plain SAT
#   mixed      DropTest_250_mixed_gas.prm      84 cuboids, 83 cylinders, 83 spheres (R 2), cycled in release order
# Each case runs in its own directory, so outputs cannot mix. The previous cuboid run took ~55 min on 36 ranks
# (DropTest_250_sat_output_20260923).
# [Claude, 2026-10-01, later] One case per job (author: 1 h limit, three separate jobs). Usage, from
# apps/showcases/Antidunes of ~/mt-dennys-mirror2 (~/mt-dennys was in use by job 4282267):
#   b=$(sbatch --parsable --time=00:30:00 --ntasks=8 jobs/run_36_droptest_250_three_shapes.sbatch build)
#   for c in cylinders cuboids mixed; do sbatch --dependency=afterok:$b --job-name=antidunes-droptest-250-$c \
#       jobs/run_36_droptest_250_three_shapes.sbatch $c; done
# The mixed case needs --ntasks=24 (DropTest_250_mixed_gas.prm uses 2x2x6 blocks, see its header); mpirun uses
# $SLURM_NTASKS.
# The build runs once, in its own job, so the three run jobs never build concurrently in the same build directory.
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
    # [Claude, 2026-10-01] randomised orientations/offsets, wood-like contact (prms/gen_droptest_250_random.py); 24 tasks
    random_cylinders) run_case random_cylinders DropTest_250_random_cylinders_gas.prm ;;
    random_cuboids)   run_case random_cuboids   DropTest_250_random_cuboids_gas.prm ;;
    random_mixed)     run_case random_mixed     DropTest_250_random_mixed_gas.prm ;;
    # [Claude, 2026-10-03] Kodam et al. (2010b) packing with a round container, one-by-one release (gen_kodam_round_prm.py); 16 tasks
    kodam_round)      run_case kodam_round      DropTest_kodam_round_gas.prm ;;
    # [Claude, 2026-10-10] the same with kodamBandEdgeAxialCheck true (energy spike of 4292501); 16 tasks as 4292501
    kodam_round_axial) run_case kodam_round_axial DropTest_kodam_round_axial_gas.prm ;;
    # [Claude, 2026-10-03] Kodam packing comparison, the four detection variants of thesis Sec. 5.2 (9 tasks each),
    # rerun on the current code so that their wall-clock times are measured together
    kodam_fallback)   run_case kodam_fallback   DropTest_kodam_packing_gas.prm ;;
    kodam_gjk)        run_case kodam_gjk        DropTest_kodam_packing_gjk_gas.prm ;;
    kodam_pure)       run_case kodam_pure       DropTest_kodam_packing_pureanalytical_gas.prm ;;
    kodam_refined)    run_case kodam_refined    DropTest_kodam_packing_refined_gas.prm ;;
    *) echo "unknown case: $CASE" >&2; exit 1 ;;
esac
echo "## DONE. Output under $ARCHIVE_ROOT/$CASE"
