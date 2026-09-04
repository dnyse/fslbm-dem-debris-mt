#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-cuboids-slitpass-asan
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_cuboids_slitpass_asan_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_cuboids_slitpass_asan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# ASan (AddressSanitizer) instrumented rerun of thousand_cuboids_two_pillars_slitpass_gas.prm
# (100 floating boxes, single un-staggered release at step 10000, useHardcodedUnits false), which
# has crashed with a DIFFERENT signature on every plain rerun of the identical config: glibc heap
# corruption ("free(): invalid next size (fast)", job 4094949; "corrupted size vs. prev_size" on
# the cylinder sibling job 4094948) on some runs, and an invalid-quaternion NaN assertion on others
# (job 4095440) -- that inconsistency across identical reruns is the classic fingerprint of a real
# out-of-bounds write, not a numerical/physics instability. Box and cylinder floaters share the
# same V_sub sampling/clip code path (ParticleKernels.h's "combined box + cylinder 'clip floater'
# list", the plane-fit + clip buoyancy pipeline, [BoxBuoyancy] log tag) -- spheres use a completely
# separate implementation and have never shown this failure family, which is why this targets boxes
# specifically. Following this project's own precedent
# (run_72_mixed_shape_sedimentation_showcase_45floaters_fastcurrent_asan.sbatch and
# run_72_thousand_spheres_two_pillars_slitpass_v6_asan.sbatch, both of which pinpointed a real bad
# memory access in a single ASan run instead of further log-based guessing) rather than continuing
# to hypothesis-test through plain-log reruns.
#
# Builds into its own build-asan-cuboids/ directory so it can run independently of any other ASan
# build in progress.
#
# The plain (non-ASan) reruns of this exact scenario have crashed within ~9-53 minutes; even at
# ASan's typical 2-3x slowdown this should complete well within the 3h walltime requested above.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_ASAN="$REPO/build-asan-cuboids"
BIN_DIR="$BUILD_ASAN/apps/showcases/Antidunes"
PRM="$REPO/apps/showcases/Antidunes/prms/thousand_cuboids_two_pillars_slitpass_gas.prm"

echo "==> Loading modules"
# See run_72_mixed_shape_sedimentation_showcase_45floaters_fastcurrent_asan.sbatch's own comment for
# why gcc/14.2.0 (not the system gcc/11.5.0) is required for a working libasan runtime here.
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
# unconditionally for the avgParticleDiameter/domain-scaling derivation even though
# useParticleBed=false here -- this fresh build-asan-cuboids/ dir has never had BedGeneration run
# in it, unlike build/ which already has one from earlier this session)"
mkdir -p "$BIN_DIR"
cp "$REPO/build/apps/showcases/Antidunes/spheres_out.dat" "$BIN_DIR/spheres_out.dat"

echo "==> Running Antidunes under ASan (np=72, prm=$PRM)"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=$BIN_DIR/asan_log"

cd "$BIN_DIR"
mpirun -n 72 ./Antidunes "$PRM"

echo "==> ASan report files (if any):"
ls -la "$BIN_DIR"/asan_log.* 2>/dev/null || echo "   (none written -- no violation caught, or run didn't reach one)"
