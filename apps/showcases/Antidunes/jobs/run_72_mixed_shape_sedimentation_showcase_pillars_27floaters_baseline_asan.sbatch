#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-shape-sedimentation-showcase-pillars-27floaters-baseline-asan
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_sedimentation_showcase_pillars_27floaters_baseline_asan_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_sedimentation_showcase_pillars_27floaters_baseline_asan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# ASan (AddressSanitizer) instrumented rerun of the 27-floater crash investigation's
# earliest/cleanest failure: mixed_shape_sedimentation_showcase_pillars_27floaters_baseline_gas.prm
# (27 floaters, radius 4.5, targetMeanVelocityMagnitude 0.005 -- pure floater-count isolation)
# crashed via a glibc/quaternion-NaN heap-corruption-class bug at t~16518 in job 3885722, and
# the same class of bug recurred at the geometrically-equivalent step (~33838) in the 1.5x
# finer-resolution rerun (job 3887984) -- ruling out a discretization/resolution artifact and
# pointing at a genuine memory-corruption bug (most likely mesa_pd::mpi::ContactFilter's
# rank-boundary contact-ownership race, this file's own header's documented earlier crash
# class, jobs 3854234/3854357). Since log forensics couldn't pin the exact bad write, this
# script rebuilds Antidunes with AddressSanitizer to catch it directly with an exact stack
# trace instead of more circumstantial timing correlations.
#
# Uses the ORIGINAL (non-finegrid) 27floaters_baseline config deliberately -- it's the
# cheapest reproducer (crashed at step 16518 in ~700s without ASan), so even at ASan's
# typical 2-3x slowdown this should complete well within the 2h walltime requested above.
#
# Builds into a SEPARATE build-asan/ directory (sibling of the normal build/) so this
# never disturbs the regular (non-instrumented) build antidunes.sh maintains. antidunes.sh
# itself has no option to inject extra compiler flags, so this script inlines the same
# module-load/conda/CMake logic antidunes.sh uses for --configure, adding the sanitizer
# flags, then reuses antidunes.sh only for the --build/--run steps against --build-dir
# build-asan (which just invokes `cmake --build`/`mpirun` and doesn't touch configure flags).
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
PRM="$REPO/apps/showcases/Antidunes/prms/mixed_shape_sedimentation_showcase_pillars_27floaters_baseline_gas.prm"

echo "==> Loading modules"
# gcc/11.5.0's libasan.so linker stub points at /usr/lib64/libasan.so.6.0.0, which isn't
# installed on this system (confirmed via `rpm -qa | grep asan` -> empty) -- ASan compiles
# but fails to link. gcc/14.2.0 (Spack-provided) bundles its own working libasan runtime.
# Verified on-cluster (compute node f0355): a trivial MPI+C++ program built with gcc/14.2.0
# + openmpi/4.1.3-gcc11.5.0 (there is no openmpi module built against gcc14/15) links and
# runs cleanly across ranks, and correctly catches a deliberate heap-buffer-overflow with a
# full stack trace -- libstdc++ from gcc-14.2.0 is a backward-compatible superset of what the
# gcc-11-built OpenMPI needs. `mpicxx -show` hardcodes /usr/bin/g++ (system gcc11) regardless
# of which gcc module is loaded, so OMPI_CXX/OMPI_CC are set explicitly below as a safety net
# -- though this script never actually invokes mpicxx as the compiler (CMake's FindMPI only
# scrapes flags from the wrapper; CMAKE_CXX_COMPILER below points straight at the loaded g++).
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
# build-asan/ from a prior (e.g. gcc/11.5.0) configure attempt before reconfiguring.
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
