#!/bin/bash -l
#SBATCH --job-name=antidunes-hopper-clog-cyl-halfrows-pass3-precursor-asan
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/hopper_clog_cyl_halfrows_pass3_precursor_asan_%j.out
#SBATCH --error=sbatch-output/error/hopper_clog_cyl_halfrows_pass3_precursor_asan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-25] AddressSanitizer instrumented rerun of hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm
# (the scenario behind jobs 4262624/4264769/4264868/4264997), to catch the collide-stage NaN's real
# mechanism directly. THESIS_PLAN.md's 2026-09-25 entries establish that the crash cell's captured
# pre-collide PDF/force state, replayed through the REAL compiled generated collision kernel three
# independent ways (a hand-ported Python arithmetic check, a local C++ unit test, and the same unit
# test rerun on this cluster's own production toolchain -- ANTIDUNES_CollideKernelCrashReplay, job
# 4265384), ALWAYS comes out fully finite. PDFs, the force-density field, omega, and the collide
# Sweep call site's field-indexing/ghost-layer wiring have all been individually ruled out as the
# source. The only remaining explanation is memory corruption or a stale read landing in this exact
# PDF field slot from somewhere else in the real multi-rank simulation, between the moment our own
# log probe reads it and the moment collide's own kernel invocation reads it -- something no amount
# of further log-based probing or input replay can distinguish. Following this project's own
# established precedent (run_72_test_bridge_v1_asan.sbatch) of using ASan to get a direct answer
# instead of further indirect guessing.
#
# Builds into its own build-asan-hopper/ directory (does not touch build/ or any other build-asan-*/
# directory another scenario's ASan work may be using concurrently).
#
# NOTE: this scenario requires domainSize<768,40,180> / cellsPerBlock<64,20,30> => exactly 144
# processes (12*2*6) -- setupBlockForest aborts on any other count, so this cannot be shrunk to a
# cheaper single-node ASan run without changing the domain decomposition (and possibly losing the
# crash reproduction). useParticleBed=false in this prm -- no BedGeneration/spheres_out.dat needed.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BUILD_ASAN="$REPO/build-asan-hopper"
BIN_DIR="$BUILD_ASAN/apps/showcases/Antidunes"
PRM="$REPO/apps/showcases/Antidunes/prms/hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm"

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

echo "==> Reusing existing $BUILD_ASAN if present (already built successfully in job 4265406 -- that
# run's failure was a UCX/InfiniBand RDMA registration error at step 400, 'ibv_reg_mr(...) failed: Bad
# address', a well-known ASan-vs-RDMA incompatibility, NOT a build problem or the real crash we're
# after -- see the UCX_TLS/OMPI_MCA fix below. Only wipe it if it doesn't exist yet.)"
[[ -d "$BUILD_ASAN" ]] || echo "   (no existing build dir -- will configure+build from scratch)"

echo "==> Configuring ASan build (build dir: $BUILD_ASAN)"
cmake -B "$BUILD_ASAN" \
    -DCMAKE_C_COMPILER="$(which gcc)" \
    -DCMAKE_CXX_COMPILER="$(which g++)" \
    -DWALBERLA_BUILD_SHOWCASES=ON \
    -DWALBERLA_BUILD_WITH_CODEGEN=ON \
    -DWALBERLA_BUILD_WITH_MPI=ON \
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

mkdir -p "$BIN_DIR"

# Move aside any leftover vtk-out/eval/snapshots from the previous (RDMA-failed) attempt rather than
# let the new run collide with them.
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

echo "==> Running Antidunes under ASan (np=144, prm=$PRM)"
export ASAN_OPTIONS="detect_leaks=0:abort_on_error=1:halt_on_error=1:log_path=$BIN_DIR/asan_log"

# Job 4265406's first attempt died at step 400 with a UCX/InfiniBand RDMA memory-registration error
# ("ibv_reg_mr(...) failed: Bad address" / "failed to register address ... on md[4]=mlx5_0"), not a
# real crash or ASan violation -- a known incompatibility between ASan's shadow-memory/redzone layout
# and RDMA pinned-memory registration. Force UCX/OpenMPI onto TCP+shared-memory transports only,
# skipping InfiniBand RDMA registration entirely, at the cost of slower inter-node communication
# (acceptable for a debug run).
export UCX_TLS="tcp,self,sm"
export UCX_NET_DEVICES="all"
export OMPI_MCA_pml="ob1"
export OMPI_MCA_btl="^openib"

# a crashed run must still reach the ASan-report listing below, so capture the exit status
set +e
cd "$BIN_DIR"
mpirun -n 144 ./Antidunes "$PRM"
RUN_RC=$?
set -e
echo "## run exited with status $RUN_RC"

echo "==> ASan report files (if any):"
ls -la "$BIN_DIR"/asan_log.* 2>/dev/null || echo "   (none written -- no violation caught, or run didn't reach one)"

exit "$RUN_RC"
