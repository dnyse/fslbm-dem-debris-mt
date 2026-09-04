#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-cuboids-slitpass-asan-debug
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_cuboids_slitpass_asan_debug_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_cuboids_slitpass_asan_debug_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Debug (NOT RelWithDebInfo) + ASan rerun of thousand_cuboids_two_pillars_slitpass_gas.prm, after
# fixing the heap-buffer-overflow the prior RelWithDebInfo+ASan run (job 4099923) pinpointed exactly:
# SimpleBB.h's treatDirection() looked up a particle index via uidToIdx() and used it to read
# position/velocity through mesa_pd::getVelocityAtWFPoint() WITHOUT checking it was valid in a
# release build -- the WALBERLA_ASSERT_UNEQUAL guarding that read is compiled out under NDEBUG
# (RelWithDebInfo defines it: CMAKE_CXX_FLAGS_RELWITHDEBINFO=-O2 -g -DNDEBUG), so a stale UID -> idx
# mapping silently read 8 bytes before a ParticleStorage-backed vector's allocation instead of
# aborting cleanly. Fixed in SimpleBB.h to skip the link (same pattern as the existing dry-link
# guard right above it in that file) instead of reading through an invalid index either way.
#
# THIS run switches CMAKE_BUILD_TYPE to Debug specifically (CMAKE_CXX_FLAGS_DEBUG=-g, no -DNDEBUG,
# confirmed from this project's own CMakeCache.txt) so every WALBERLA_ASSERT_* in the codebase is
# actually active, not just the ones already promoted to WALBERLA_CHECK/real runtime checks --
# the whole point being to surface any OTHER assertion that would otherwise fail silently in
# production (RelWithDebInfo/Release), the same way this exact one did before today's fix. Kept ASan
# on top since it's a different, complementary class of bug detector and costs nothing extra to run
# alongside Debug's own checking. Time limit raised to 8h (from the ASan-only run's 3h) since
# Debug's -O0 (no optimization at all) on top of ASan's own overhead is expected to be considerably
# slower than the RelWithDebInfo+ASan run that found the original bug in under an hour.
#
# Builds into its own build-debug-asan-cuboids/ directory, independent of build-asan-cuboids/ (the
# RelWithDebInfo+ASan build from the prior run) and any other build in progress.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_DIR="$REPO/build-debug-asan-cuboids"
BIN_DIR="$BUILD_DIR/apps/showcases/Antidunes"
PRM="$REPO/apps/showcases/Antidunes/prms/thousand_cuboids_two_pillars_slitpass_gas.prm"

echo "==> Loading modules"
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

echo "==> Clearing any stale $BUILD_DIR"
rm -rf "$BUILD_DIR"

echo "==> Configuring Debug+ASan build (build dir: $BUILD_DIR)"
cmake -B "$BUILD_DIR" \
    -DCMAKE_C_COMPILER="$(which gcc)" \
    -DCMAKE_CXX_COMPILER="$(which g++)" \
    -DWALBERLA_BUILD_SHOWCASES=ON \
    -DWALBERLA_BUILD_WITH_CODEGEN=ON \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
    -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address" \
    -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=address" \
    ${PYTHON_ROOT:+-DPython_ROOT_DIR="$PYTHON_ROOT"} \
    "$REPO"

echo "==> Building Antidunes (Debug+ASan)"
"$REPO/apps/showcases/Antidunes/antidunes.sh" --local \
    --build --build-dir "$BUILD_DIR" --jobs "$(nproc)"

echo "==> Copying spheres_out.dat from the regular build dir"
mkdir -p "$BIN_DIR"
cp "$REPO/build/apps/showcases/Antidunes/spheres_out.dat" "$BIN_DIR/spheres_out.dat"

echo "==> Running Antidunes under Debug+ASan (np=72, prm=$PRM)"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=$BIN_DIR/asan_log"

cd "$BIN_DIR"
mpirun -n 72 ./Antidunes "$PRM"

echo "==> ASan report files (if any):"
ls -la "$BIN_DIR"/asan_log.* 2>/dev/null || echo "   (none written -- no violation caught, or run didn't reach one)"
