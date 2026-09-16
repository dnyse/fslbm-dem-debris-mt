#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-stillwater-d075-reproducibility
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_stillwater_d075_reproducibility_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_stillwater_d075_reproducibility_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Reproducibility + root-cause probe for the cylinder d=0.75 real-gravity run that diverged
# (rho=nan at step 11000) in the original 3x cylinder gravity-comparison sweep. Three cases,
# one build, run back to back:
#   1. asis        -- cylinder_stillwater_realgravity981_3xsphere_d075_gas.prm, unchanged.
#                      Literal rerun of the exact crashed config: is the divergence deterministic
#                      (reproduces at the same step) or a one-off?
#   2. reducedz     -- same but domainSize Z 120->90 (headroom above the liquid ~60->~30 LU).
#                      Tests whether Z domain size affects the divergence. 6x3x3=54 processes
#                      (fewer blocks in Z), launched with --np 54 within this job's 72-task
#                      allocation.
#   3. dyngravity   -- same but useHardcodedUnits-only gravity (no override): gravity=7.02618e-06
#                      (g_SI~=1.36 m/s^2), ~7.2x smaller than the solved real-9.81 value. Tests
#                      whether gravity magnitude itself (not just density-dependent settling
#                      transient) drives the divergence.
#
# Each case's own stdout is also teed to a per-case log and grepped for [LOW-RHO] so the job's
# final summary states DIVERGED/completed without having to dig through the combined .out file.

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
STAMP="$(date +%Y%m%d_%H%M%S)"

if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_ROOT"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$1" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$1/$d"
    done
}

run_case() {  # $1 = banner, $2 = prm, $3 = label, $4 = np
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    local LOG="$ARCHIVE_ROOT/${3}_run.log"
    "$SH" --run --np "$4" --prm "$PRM_DIR/$2" 2>&1 | tee "$LOG"
    local N_LOWRHO
    N_LOWRHO=$(grep -c '\[LOW-RHO\]' "$LOG" || true)
    if [[ "$N_LOWRHO" -gt 0 ]]; then
        echo "## RESULT [$3]: DIVERGED -- $N_LOWRHO [LOW-RHO] warnings logged"
    else
        echo "## RESULT [$3]: completed cleanly (no [LOW-RHO] warnings)"
    fi
    archive_outputs "$3"
}

run_case "RUN asis (exact crashed config, real gravity g_SI=9.81, domain Z=120)" \
    cylinder_stillwater_realgravity981_3xsphere_d075_gas.prm asis 72
run_case "RUN reducedz (domain Z=120->90, same gravity)" \
    cylinder_stillwater_realgravity981_3xsphere_d075_reducedz_gas.prm reducedz 54
run_case "RUN dyngravity (native gravity=7.02618e-06, domain Z=120 unchanged)" \
    cylinder_stillwater_realgravity981_3xsphere_d075_dyngravity_gas.prm dyngravity 72

echo "## DONE. Archived under $ARCHIVE_ROOT/{asis,reducedz,dyngravity}_{vtk-out,eval,snapshots,run.log}"
