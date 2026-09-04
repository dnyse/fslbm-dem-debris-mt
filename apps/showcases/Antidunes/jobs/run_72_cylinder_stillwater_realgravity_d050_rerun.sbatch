#!/bin/bash -l
#SBATCH --job-name=antidunes-cylinder-stillwater-realgravity-d050-rerun
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cylinder_stillwater_realgravity_d050_rerun_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cylinder_stillwater_realgravity_d050_rerun_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Plain rerun of cylinder_stillwater_realgravity981_3xsphere_gas.prm (density ratio 0.5,
# real gravity g_SI=9.81). The original sweep run (job 4035437) completed but was still
# visibly rising at step 50,000 (60.18 -> 60.31 LU, not yet flat like the Bogner-gravity
# sibling) -- this reruns the identical, unmodified config.

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

run_case "RUN real gravity d=0.50 rerun (g_SI=9.81, domain Z=120)" \
    cylinder_stillwater_realgravity981_3xsphere_gas.prm rerun 72

echo "## DONE. Archived under $ARCHIVE_ROOT/rerun_{vtk-out,eval,snapshots,run.log}"
