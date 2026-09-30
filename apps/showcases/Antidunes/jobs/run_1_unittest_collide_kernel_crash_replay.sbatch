#!/bin/bash -l
#SBATCH --job-name=antidunes-unittest-collide-crash-replay
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --output=sbatch-output/output/unittest_collide_crash_replay_%j.out
#SBATCH --error=sbatch-output/error/unittest_collide_crash_replay_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-25] Runs the new ANTIDUNES_CollideKernelCrashReplay unit test (see
# tests/showcases_antidunes/CollideKernelCrashReplayTest.cpp) on the cluster itself, on the SAME
# compiler/toolchain/optimization flags used for the real production hopper-clogging runs, as a
# cross-check of the local (workstation) result. The test replays job 4264997's exact captured
# pre-collide PDF/force state for the hopper crash cell through the REAL compiled generated Cumulant
# collision kernel -- locally it came out fully finite (rho=1.0069218, allFinite=true), matching a
# separate hand-ported Python replay of the same arithmetic. Running it here checks whether a
# cluster-specific compiler/flag difference (e.g. -march=native, -ffast-math, or any other
# optimization the production cluster build uses that the local workstation build might not) changes
# that outcome -- if it does NOT, that closes off "compiler-specific floating point behaviour" as an
# explanation too, further narrowing the real root cause toward genuine memory corruption/stale reads
# needing an ASan/Valgrind pass on the full simulation binary (see THESIS_PLAN.md's 2026-09-25 entries).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_DIR="$REPO/build"

module load gcc/11.5.0
module load openmpi/4.1.3-gcc11.5.0
module load cmake/3.31.8
module load python/3.13-miniforge

CONDA_BASE="$(conda info --base 2>/dev/null)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate walberla 2>/dev/null || true

echo "############################################################"
echo "## Reconfiguring (new test target added to CMakeLists.txt)"
echo "############################################################"
cmake "$BUILD_DIR"

echo "############################################################"
echo "## Building ANTIDUNES_CollideKernelCrashReplay"
echo "############################################################"
cmake --build "$BUILD_DIR" --target ANTIDUNES_CollideKernelCrashReplay -- -j"$(nproc)"

echo "############################################################"
echo "## Running ANTIDUNES_CollideKernelCrashReplay"
echo "############################################################"
ctest --test-dir "$BUILD_DIR" -R ANTIDUNES_CollideKernelCrashReplay --output-on-failure
