#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-d025-reproducibility
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cuboid_stillwater_d025_reproducibility_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cuboid_stillwater_d025_reproducibility_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Reproducibility + root-cause probe for the cuboid d=0.25 Bogner-gravity run that diverged
# (rho=nan at step 2500) in the original 3x cuboid gravity-comparison sweep. Three cases,
# one build, run back to back:
#   1. asis        -- cuboid_stillwater_bognerboxgravity_3xsphere_d025_gas.prm, unchanged (this
#                      now uses the NEW s x 2s x s = 30x60x30 cuboid shape, not the original
#                      48x48x60 box that actually crashed -- see reproducibility caveat below).
#                      Bogner gravity g=7.5e-4, tau=1.9, domain Z=120.
#   2. reducedz     -- same but domainSize Z 120->90 (headroom ~60->~30 LU). 6x3x3=54 processes,
#                      launched with --np 54 within this job's 72-task allocation.
#   3. dyngravity   -- same but native gravity (no override): gravity=7.02618e-06 (g_SI~=1.36
#                      m/s^2) instead of Bogner's own g=7.5e-4 -- tests whether it's specifically
#                      Bogner's large gravity magnitude (not the shape or density) driving the
#                      divergence.
#
# CAVEAT: because the cuboid's halfExtents were changed to the s x 2s x s ratio (30x60x30)
# after this crash was found, "asis" here is NOT a byte-identical rerun of the original crashed
# setup (48x48x60) -- it tests whether the NEW shape still diverges under otherwise-identical
# conditions, which is the more useful question at this point, but keep the shape change in mind
# when reading "reproduced" vs "not reproduced" below.
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

run_case "RUN asis (new 30x60x30 shape, exact crashed density/gravity, Bogner gravity, domain Z=120)" \
    cuboid_stillwater_bognerboxgravity_3xsphere_d025_gas.prm asis 72
run_case "RUN reducedz (domain Z=120->90, same gravity)" \
    cuboid_stillwater_bognerboxgravity_3xsphere_d025_reducedz_gas.prm reducedz 54
run_case "RUN dyngravity (native gravity=7.02618e-06, domain Z=120 unchanged)" \
    cuboid_stillwater_bognerboxgravity_3xsphere_d025_dyngravity_gas.prm dyngravity 72

echo "## DONE. Archived under $ARCHIVE_ROOT/{asis,reducedz,dyngravity}_{vtk-out,eval,snapshots,run.log}"
