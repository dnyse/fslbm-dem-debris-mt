#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-smallshape-missing
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cuboid_stillwater_smallshape_missing_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cuboid_stillwater_smallshape_missing_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Fills the gap discovered while auditing the cuboid gravity-comparison sweep artifact: the
# ORIGINAL sweep jobs (4035436/4035438/4035440) ran before the cuboid shape was redesigned
# from 48x60x48 (halfExtents<24,30,24>) to the current 30x60x30 s x 2s x s box
# (halfExtents<15,30,15>) -- confirmed by grepping "halfExtents" out of the archived .out log,
# which shows <24,30,24> even though every current cuboid prm on disk already has <15,30,15>.
# Combined with one gravity crashing per density in that old-shape sweep, the artifact ended
# up with exactly one (old-shape) line per density panel, never a real-vs-Bogner pair.
#
# The reproducibility experiment already produced three clean SMALL-shape runs that go unused
# in the current chart: bogner d=0.25 and d=0.50 (job 4038342/4038343 "asis"), and real d=0.75
# (job 4038344 "asis"). This job runs exactly the three MISSING small-shape/gravity
# combinations -- real d=0.25, real d=0.50, bogner d=0.75 -- so each density can show a
# genuine same-shape two-gravity comparison.

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

run_case "RUN real d025 (small 30x60x30 shape, g_SI=9.81, domain Z=120)" \
    cuboid_stillwater_realgravity981_3xsphere_d025_gas.prm real_d025 72
run_case "RUN real d050 (small 30x60x30 shape, g_SI=9.81, domain Z=120)" \
    cuboid_stillwater_realgravity981_3xsphere_gas.prm real_d050 72
run_case "RUN bogner d075 (small 30x60x30 shape, g=7.5e-4 tau=1.9, domain Z=120)" \
    cuboid_stillwater_bognerboxgravity_3xsphere_d075_gas.prm bogner_d075 72

echo "## DONE. Archived under $ARCHIVE_ROOT/{real_d025,real_d050,bogner_d075}_{vtk-out,eval,snapshots,run.log}"
