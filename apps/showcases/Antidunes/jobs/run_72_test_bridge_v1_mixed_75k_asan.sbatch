#!/bin/bash -l
#SBATCH --job-name=antidunes-test-bridge-v1-mixed-75k-asan
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_test_bridge_v1_mixed_75k_asan_%j.out
#SBATCH --error=sbatch-output/error/antidunes_test_bridge_v1_mixed_75k_asan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# ASan-instrumented rerun of test_bridge_v1_mixed_75k_gas.prm, which crashed with a plain SIGSEGV
# (rank 19, exit 139) at ~step 3986 of job 4065181 -- same "Bogner state selection ... unclassified
# FormerMO neighbour(s) ... provisionally interface pending Pass-3 re-resolution" flood immediately
# preceding the crash as job 4064993 (the smaller radius-2-floater run), which was traced via ASan
# (job 4065016) to a heap-buffer-overflow caused by an unbounded overlap in ComputeTunnelContacts --
# already fixed with a radial gate (verified clean on job 4065061, 30k steps, radius-2 floaters). This
# run uses much larger floaters (radius/halfExtents 4.5/<3.5,7,3.5>, close to or bigger than the
# side-arch radius of 5) and crashes much earlier (step ~3986 vs ~5070), so this is either a residual
# instance of the same bug class at larger scale, or something new triggered by a floater whose
# effective radius approaches/exceeds the tunnel or pillar geometry itself. Reuses build-asan-bridge/
# (already exists from the earlier ASan run; rebuilt fresh here regardless).
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_ASAN="$REPO/build-asan-bridge"
BIN_DIR="$BUILD_ASAN/apps/showcases/Antidunes"
PRM="$REPO/apps/showcases/Antidunes/prms/test_bridge_v1_mixed_75k_gas.prm"

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

echo "==> Copying spheres_out.dat from the regular build dir"
mkdir -p "$BIN_DIR"
cp "$REPO/build/apps/showcases/Antidunes/spheres_out.dat" "$BIN_DIR/spheres_out.dat"

echo "==> Running Antidunes under ASan (np=72, prm=$PRM)"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=$BIN_DIR/asan_log"

cd "$BIN_DIR"
mpirun -n 72 ./Antidunes "$PRM"

echo "==> ASan report files (if any):"
ls -la "$BIN_DIR"/asan_log.* 2>/dev/null || echo "   (none written -- no violation caught, or run didn't reach one)"
