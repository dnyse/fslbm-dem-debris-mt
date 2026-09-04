#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-d075-pass3
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cuboid_stillwater_d075_pass3_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cuboid_stillwater_d075_pass3_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Pass-3 re-resolution check for cuboid d=0.75, real gravity g_SI=9.81, smaller 30x60x30
# (s x 2s x s) shape -- the same config as the reproducibility "asis" run (job 4038344),
# which already completed cleanly with usePass3Reresolution OFF (default). This run flips
# that one flag on to see whether re-resolving FormerMO cells that fell back to interface
# on Pass 2 (PdfReconstructionManager.h, usePass3Reresolution, SimulationParams.h:439)
# changes the settling trajectory or reconstruction diagnostics relative to that baseline.
# Not expected to "fix a crash" -- there isn't one to fix in this shape/gravity combination
# any more -- this is a single-variable check of Pass 3's effect on its own.

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
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_ROOT"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/${1}_$d"
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

run_case "RUN pass3 (30x60x30 shape, real gravity g_SI=9.81, domain Z=120, usePass3Reresolution=true)" \
    cuboid_stillwater_realgravity981_3xsphere_d075_pass3_gas.prm pass3 72

echo "## DONE. Archived under $ARCHIVE_ROOT/pass3_{vtk-out,eval,snapshots,run.log}"
