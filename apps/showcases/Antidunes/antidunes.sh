#!/bin/bash -l
# antidunes.sh — build and run helper for the Antidunes showcase
# See docs/HOWTO_BUILD_AND_RUN.md for full documentation.
set -euo pipefail

# ── paths derived from script location ────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# ── defaults ──────────────────────────────────────────────────────────────────
LOAD_MODULES=true
DO_CONFIGURE=false
DO_BUILD=false
DO_TESTS=false
DO_BED=false
DO_RUN=false
FRESH=false
BUILD_DIR="$REPO_ROOT/build"
JOBS="$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)"
NP=72
BED_NP=72
PRM="$SCRIPT_DIR/Antidunes.prm"
BED_PRM="$SCRIPT_DIR/BedGeneration.prm"
CONDA_ENV="walberla"
CONDA_PY_VERSION="3.11"
CONDA_PKGS=("pystencils~=1.3" "lbmpy~=1.3" "jinja2")
DO_SYNC=false
REMOTE="fritz.iwia"
REMOTE_PATH="~/mt-dennys/"

# ── usage ─────────────────────────────────────────────────────────────────────
usage() {
cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Steps (pick one or more, or use --all):
  --configure       CMake configure
  --build           Build Antidunes + BedGeneration targets
  --tests           Build + run the Antidunes contact/geometry unit tests (ctest, single-process,
                    no mpirun needed) -- ANTIDUNES_ContactDetection + ANTIDUNES_ClipAnalyticGeometry
  --bed             Run BedGeneration (generates spheres_out.dat)
  --run             Run Antidunes main simulation
  --all             configure + build + tests + bed + run
  --sync            rsync repo to remote (excludes build/ and dotfiles)

Environment:
  --local           Skip module loading — use this on a local machine
  --conda-env NAME  Conda environment to activate (default: $CONDA_ENV)
  --fresh           Pass --fresh to CMake (force clean reconfigure)

Paths:
  --build-dir DIR   CMake build directory   (default: $BUILD_DIR)
  --prm FILE        Antidunes parameter file (default: Antidunes.prm in script dir)
  --bed-prm FILE    BedGeneration prm file   (default: BedGeneration.prm in script dir)
  --remote HOST     SSH host for --sync      (default: $REMOTE)
  --remote-path P   Remote destination path  (default: $REMOTE_PATH)

Parallelism:
  --jobs N          Build parallelism        (default: $JOBS)
  --np N            MPI ranks for Antidunes  (default: $NP)
  --bed-np N        MPI ranks for BedGeneration (default: $BED_NP)

  -h, --help        Show this help

Examples:
  # Sync changes to supercomputer, then build and run there
  ./antidunes.sh --sync

  # Local machine — full pipeline
  ./antidunes.sh --local --all --np 8 --bed-np 1

  # Supercomputer — build only, 32 jobs
  ./antidunes.sh --configure --build --jobs 32

  # Supercomputer — build and run the contact/geometry unit tests (fast, no allocation needed)
  ./antidunes.sh --configure --build --tests --jobs 32

  # Supercomputer — rerun with different prm, skip rebuild
  ./antidunes.sh --run --prm myparams.prm --np 72
EOF
}

# ── argument parsing ──────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && { usage; exit 1; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --local)       LOAD_MODULES=false ;;
        --configure)   DO_CONFIGURE=true ;;
        --build)       DO_BUILD=true ;;
        --tests)       DO_TESTS=true ;;
        --bed)         DO_BED=true ;;
        --run)         DO_RUN=true ;;
        --all)         DO_CONFIGURE=true; DO_BUILD=true; DO_TESTS=true; DO_BED=true; DO_RUN=true ;;
        --fresh)       FRESH=true ;;
        --build-dir)   BUILD_DIR="$2"; shift ;;
        --jobs)        JOBS="$2"; shift ;;
        --np)          NP="$2"; shift ;;
        --bed-np)      BED_NP="$2"; shift ;;
        --prm)         PRM="$2"; shift ;;
        --bed-prm)     BED_PRM="$2"; shift ;;
        --conda-env)   CONDA_ENV="$2"; shift ;;
        --sync)        DO_SYNC=true ;;
        --remote)      REMOTE="$2"; shift ;;
        --remote-path) REMOTE_PATH="$2"; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown option: $1"; echo; usage; exit 1 ;;
    esac
    shift
done

BIN_DIR="$BUILD_DIR/apps/showcases/Antidunes"

# Resolve prm paths to absolute so they survive the 'cd $BIN_DIR' before mpirun
PRM="$(realpath "$PRM")"
BED_PRM="$(realpath "$BED_PRM")"

# ── helpers ───────────────────────────────────────────────────────────────────
log() { echo; echo "==> $*"; }

# ── sync to remote ────────────────────────────────────────────────────────────
if $DO_SYNC; then
    log "Syncing to $REMOTE:$REMOTE_PATH"
    rsync -avz \
        --exclude='build/' \
        --exclude='.*' \
        "$REPO_ROOT/" "$REMOTE:$REMOTE_PATH"
fi

# ── module loading (supercomputer only) ───────────────────────────────────────
if $LOAD_MODULES; then
    log "Loading modules"
    module load gcc/11.5.0
    module load openmpi/4.1.3-gcc11.5.0
    module load cmake/3.31.8
    module load python/3.13-miniforge

    log "Activating conda environment: $CONDA_ENV"
    # conda activate does not work in non-interactive shells without this
    CONDA_BASE="$(conda info --base 2>/dev/null)"
    # shellcheck source=/dev/null
    source "$CONDA_BASE/etc/profile.d/conda.sh"

    # Create the env on first use if it doesn't exist yet (activation is the test).
    if ! conda activate "$CONDA_ENV" 2>/dev/null; then
        log "Conda env '$CONDA_ENV' not found — creating it (python=$CONDA_PY_VERSION)"
        # Shared module installs (e.g. NHR miniforge) keep their pkgs/envs dirs
        # read-only; point conda at writable dirs in $HOME so create can proceed.
        conda config --add pkgs_dirs "$HOME/.conda/pkgs" 2>/dev/null || true
        conda config --add envs_dirs "$HOME/.conda/envs" 2>/dev/null || true
        conda create -y -n "$CONDA_ENV" "python=$CONDA_PY_VERSION"
        conda activate "$CONDA_ENV"
        log "Installing Python packages: ${CONDA_PKGS[*]}"
        pip install "${CONDA_PKGS[@]}"
    fi
else
    log "Skipping module load (--local)"
    # If a conda env is already active, honour it; otherwise use the system Python.
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
        echo "    No conda env active. Using system Python: $(which python3 2>/dev/null || echo 'not found')"
    else
        echo "    Using active conda env: $CONDA_PREFIX"
    fi
fi

# Resolve Python root for CMake (after potential conda activation)
PYTHON_ROOT="${CONDA_PREFIX:-$(python3 -c 'import sys; print(sys.prefix)' 2>/dev/null || echo '')}"

# ── configure ─────────────────────────────────────────────────────────────────
if $DO_CONFIGURE; then
    log "Configuring (repo root: $REPO_ROOT)"
    CMAKE_CMD=(
        cmake
        $($FRESH && echo "--fresh" || true)
        -B "$BUILD_DIR"
        -DCMAKE_C_COMPILER="$(which gcc)"
        -DCMAKE_CXX_COMPILER="$(which g++)"
        -DWALBERLA_BUILD_SHOWCASES=ON
        -DWALBERLA_BUILD_WITH_CODEGEN=ON
        -DCMAKE_BUILD_TYPE=Release
    )
    [[ -n "$PYTHON_ROOT" ]] && CMAKE_CMD+=(-DPython_ROOT_DIR="$PYTHON_ROOT")
    CMAKE_CMD+=("$REPO_ROOT")

    cd "$REPO_ROOT"
    "${CMAKE_CMD[@]}"
fi

# ── build ─────────────────────────────────────────────────────────────────────
if $DO_BUILD; then
    log "Building Antidunes (jobs: $JOBS)"
    cmake --build "$BUILD_DIR" --target Antidunes     -- -j"$JOBS"

    log "Building BedGeneration (jobs: $JOBS)"
    cmake --build "$BUILD_DIR" --target BedGeneration -- -j"$JOBS"
fi

# ── contact/geometry unit tests ───────────────────────────────────────────────
if $DO_TESTS; then
    log "Building Antidunes unit tests (jobs: $JOBS)"
    cmake --build "$BUILD_DIR" --target ANTIDUNES_ContactDetection ANTIDUNES_ClipAnalyticGeometry -- -j"$JOBS"

    log "Running Antidunes unit tests (ctest, single-process)"
    ctest --test-dir "$BUILD_DIR" -R ANTIDUNES_ --output-on-failure
fi

# ── bed generation ────────────────────────────────────────────────────────────
if $DO_BED; then
    log "Running BedGeneration (np=$BED_NP, prm=$BED_PRM)"
    [[ -f "$BIN_DIR/BedGeneration" ]] || { echo "BedGeneration binary not found in $BIN_DIR — did you --build?"; exit 1; }
    cd "$BIN_DIR"
    mpirun -n "$BED_NP" ./BedGeneration "$BED_PRM"
fi

# ── main simulation ───────────────────────────────────────────────────────────
if $DO_RUN; then
    log "Running Antidunes (np=$NP, prm=$PRM)"
    [[ -f "$BIN_DIR/Antidunes" ]] || { echo "Antidunes binary not found in $BIN_DIR — did you --build?"; exit 1; }
    [[ -f "$BIN_DIR/spheres_out.dat" ]] || {
        echo "WARNING: spheres_out.dat not found in $BIN_DIR."
        echo "         Run --bed first, or copy it there manually."
    }
    cd "$BIN_DIR"
    mpirun -n "$NP" ./Antidunes "$PRM"
fi

log "Done."
