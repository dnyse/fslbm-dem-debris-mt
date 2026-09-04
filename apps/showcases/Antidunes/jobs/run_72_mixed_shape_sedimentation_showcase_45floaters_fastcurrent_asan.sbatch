#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-shape-sedimentation-showcase-45floaters-fastcurrent-asan
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_sedimentation_showcase_45floaters_fastcurrent_asan_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_sedimentation_showcase_45floaters_fastcurrent_asan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# ASan (AddressSanitizer) instrumented rerun of the 5-column/45-floater fastcurrent
# variant (mixed_shape_sedimentation_showcase_45floaters_fastcurrent_gas.prm), which
# crashed via the same quaternion-NaN assertion (Quaternion.h:241) as the earlier
# 27-floater crash investigation, job 3902147. Log-forensics on that run traced the
# corruption to a single block-local cell (26,7,23) on rank 18 going to rho=NaN out of
# nowhere at step 16793 -- no warning the step before -- then spreading via ordinary PDF
# streaming over ~34 steps until it reached two floaters' MEM coupling (identical,
# non-varying garbage |F_hydro| across all 10 DEM sub-cycles, the same "reading
# stale/wrong memory" fingerprint as the already-fixed ghost-layer OOB read) and
# finally the quaternion-NaN abort at step ~16828. Local X=26 (cellsPerBlock.x=32) is
# only 6 cells from that block's right edge -- the same geometric risk zone (near a
# ghost-layer/block boundary) as the bug already fixed in PdfReconstructionManager.h,
# but NOT proven to be the identical code path (that fix only covered one specific
# neighbour-read lambda in classifyFormerObstacleCell). Since ASan pinpointed the exact
# bad memory access for the 27-floater crash in a single run after months of log-based
# hypothesis testing, this script does the same here rather than digging further
# through logs.
#
# Also notable (not yet confirmed as causal, just flagged): the corrupted floater
# (uid ...748, a box) started at X=28 -- the column with the TIGHTEST seam clearance
# (4 LU from the seam at X=32) in this variant's layout, see the prm's own header for
# the full clearance breakdown. If ASan's stack trace points at mesa_pd::mpi::ContactFilter
# rather than the Bogner reconstruction path, that would implicate the tighter X
# clearance directly.
#
# Builds into the SAME build-asan/ directory the 27-floater ASan run used (wiped and
# reconfigured with identical flags below -- only the prm differs, not the source, so
# reusing the directory is safe and avoids a redundant fresh clone/configure).
#
# Uses the fastcurrent 45-floater config as-is (not a cheaper isolation variant) --
# it crashed at step ~16828 in ~800s without ASan, so even at ASan's typical 2-3x
# slowdown this should complete well within the 2h walltime requested above.
#
# Resolve the repo root programmatically (mirrors antidunes.sh's own SCRIPT_DIR/REPO_ROOT
# derivation) -- no hardcoded path, works across mirror clones (mt-dennys, mt-dennys-mirror,
# mt-dennys-mirror2, ...).
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_ASAN="$REPO/build-asan"
BIN_DIR="$BUILD_ASAN/apps/showcases/Antidunes"
PRM="$REPO/apps/showcases/Antidunes/prms/mixed_shape_sedimentation_showcase_45floaters_fastcurrent_gas.prm"

echo "==> Loading modules"
# gcc/11.5.0's libasan.so linker stub points at /usr/lib64/libasan.so.6.0.0, which isn't
# installed on this system (confirmed via `rpm -qa | grep asan` -> empty) -- ASan compiles
# but fails to link. gcc/14.2.0 (Spack-provided) bundles its own working libasan runtime,
# confirmed working with openmpi/4.1.3-gcc11.5.0 (there is no openmpi module built against
# gcc14/15) in the earlier 27-floater ASan run (job 3893620). `mpicxx -show` hardcodes
# /usr/bin/g++ (system gcc11) regardless of which gcc module is loaded, so OMPI_CXX/OMPI_CC
# are set explicitly below as a safety net -- though this script never actually invokes
# mpicxx as the compiler (CMake's FindMPI only scrapes flags from the wrapper;
# CMAKE_CXX_COMPILER below points straight at the loaded g++).
module load gcc/14.2.0
module load openmpi/4.1.3-gcc11.5.0
module load cmake/3.31.8
module load python/3.13-miniforge

export OMPI_CC="$(which gcc)"
export OMPI_CXX="$(which g++)"

echo "==> Activating conda environment: walberla"
CONDA_BASE="$(conda info --base 2>/dev/null)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate walberla

PYTHON_ROOT="${CONDA_PREFIX:-}"

# CMake refuses an in-place compiler swap on an existing build directory -- wipe any stale
# build-asan/ (e.g. from the earlier 27-floater ASan run, or a prior gcc/11.5.0 configure
# attempt) before reconfiguring.
echo "==> Clearing any stale $BUILD_ASAN"
rm -rf "$BUILD_ASAN"

echo "==> Configuring ASan build (build dir: $BUILD_ASAN)"
cmake -B "$BUILD_ASAN" \
    -DCMAKE_C_COMPILER="$(which gcc)" \
    -DCMAKE_CXX_COMPILER="$(which g++)" \
    -DWALBERLA_BUILD_SHOWCASES=ON \
    -DWALBERLA_BUILD_WITH_CODEGEN=ON \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_C_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
    -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address" \
    -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=address" \
    ${PYTHON_ROOT:+-DPython_ROOT_DIR="$PYTHON_ROOT"} \
    "$REPO"

echo "==> Building Antidunes (ASan)"
"$REPO/apps/showcases/Antidunes/antidunes.sh" --local \
    --build --build-dir "$BUILD_ASAN" --jobs "$(nproc)"

echo "==> Running Antidunes under ASan (np=72, prm=$PRM)"
# detect_leaks=0: confirmed necessary on-cluster -- OpenMPI leaks ~120KB/3400 allocations
# internally at shutdown per rank; with leak detection on, that report drowns out the actual
# corruption report and makes the job exit nonzero even on an otherwise-clean run.
# abort_on_error=1 + halt_on_error=1: stop immediately at the first real violation with a
# full stack trace, rather than continuing and potentially cascading into a second, harder
# to interpret failure.
# log_path=$BIN_DIR/asan_log (not stdout): at 72 ranks, only the rank(s) that actually hit a
# violation write a log_path.<pid> file -- keeps the report isolated instead of interleaved
# into 72 ranks' worth of stdout.
mkdir -p "$BIN_DIR"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=$BIN_DIR/asan_log"

cd "$BIN_DIR"
mpirun -n 72 ./Antidunes "$PRM"

echo "==> ASan report files (if any):"
ls -la "$BIN_DIR"/asan_log.* 2>/dev/null || echo "   (none written -- no violation caught, or run didn't reach one)"
