#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d075-resume
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=140
#SBATCH --time=07:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d075_resume_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d075_resume_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-18] Thesis item 1a (THESIS_PLAN.md) -- resumes job 4246210's own t=90000
# checkpoint and extends the run from timesteps=100000 to timesteps=300000 (i.e. 210000 further
# steps from the loaded checkpoint), per explicit request ("put the files in place on the cluster
# and run it for another 200k for a total of 300k"). See
# cuboid_stillwater_realgravity981_3xsphere_hires_d075_resume_gas.prm's own header for the full
# motivation (densityRatio=0.75 still visibly relaxing at t=100000 -- ~12deg drift in the last 20000
# steps alone, vs 0.25/0.5 both essentially flat over the same window) and that same header's own
# retraction of an earlier (incorrect) attempt to explain WHY via Bogner's Scribanti formula -- that
# formula is for a different validation case in the thesis and doesn't actually apply here. The
# slower relaxation is real and measured; the mechanism is not confirmed.
#
# Time: the original 100000-step production run measured ~1h46m (job 4246210, see its own sacct
# timing) -- roughly 63ms/step at 140 ranks. 210000 further steps at that rate is ~3.7h; 7h requested
# for margin (this resume also re-does the already-completed 90000->100000 tail once more before
# reaching genuinely new territory, since the checkpoint is 10000 steps behind the prior run's own
# final state -- cheap, not worth avoiding via more complex step-accounting).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done
# NOTE: deliberately NOT touching snapshots/ -- this job loads the existing t=90000 checkpoint.

if [[ ! -f "$BIN_DIR/snapshots/simState.file" ]]; then
    echo "FATAL: no snapshots/simState.file found in $BIN_DIR -- refusing to run with loadSnapshot=true" >&2
    exit 1
fi

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires, densityRatio=0.75 -- RESUME from t=90000, extended to t=300000"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_d075_resume_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use) -- see
# run_140_cuboid_stillwater_hires_d050.sbatch's own 2026-09-16 comment for the full rationale
# (autofs staleness, local-only is a last resort not a default).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_d075_resume"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/d075_resume" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/d075_resume/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/d075_resume/{eval,vtk-out,snapshots}"
