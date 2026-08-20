#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-asan
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_asan_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_asan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# ASan (AddressSanitizer) instrumented rerun of test_bridge_v1_gas.prm (now with piers/arches/floaters
# wired up: CylinderPillars+Pillars for the piers, ArchTunnels+ComputeTunnelContacts for the arches,
# 6 FloatingSpheres incl. a startInactive dormant pair), which crashed with a plain SIGSEGV (rank 31,
# no NaN-SCAN trigger) at ~step 5070 of job 4064993 -- right after the first floater batch
# (settlingSteps=2000) had been active for ~3000 steps, during a burst of "Bogner state selection ...
# unclassified FormerMO neighbour(s) ... provisionally interface pending Pass-3 re-resolution"
# activity -- the same signature as two earlier crashes this session in the same reconstruction path
# (both already fixed: the ghost-layer OOB read, commit d8c0a265; and this session's own
# seedConvertedCellFill obstacle-exclusion fix). Could be yet another instance of that same bug class,
# or something new from the piers/arches/floaters/dormant-floater code added this session -- following
# this project's own precedent of using ASan to get a direct answer instead of further log-based
# guessing.
#
# Builds into its own build-asan-bridge/ directory (not build-asan-v6/, which may be reused by that
# other scenario's own ASan work independently).
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_ASAN="$REPO/build-asan-bridge"
BIN_DIR="$BUILD_ASAN/apps/showcases/Antidunes"
PRM="$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_gas.prm"

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

echo "==> Copying spheres_out.dat from the regular build dir (Utility.h reads its header
# unconditionally even though useParticleBed=false here -- this fresh build-asan-bridge/ dir has
# never had BedGeneration run in it, unlike build/ which already has one from earlier this session)"
mkdir -p "$BIN_DIR"
cp "$REPO/build/apps/showcases/Antidunes/spheres_out.dat" "$BIN_DIR/spheres_out.dat"

echo "==> Running Antidunes under ASan (np=72, prm=$PRM)"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=$BIN_DIR/asan_log"

cd "$BIN_DIR"
mpirun -n 72 ./Antidunes "$PRM"

echo "==> ASan report files (if any):"
ls -la "$BIN_DIR"/asan_log.* 2>/dev/null || echo "   (none written -- no violation caught, or run didn't reach one)"
