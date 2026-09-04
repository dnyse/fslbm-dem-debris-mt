#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-125floaters-realgravity-twopillars-flow
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=05:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_125floaters_realgravity_twopillars_flow_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_125floaters_realgravity_twopillars_flow_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Two-pillars/flow sibling of run_144_box_sedimentation_showcase_analytical_125floaters_realgravity.sbatch
# (job 4104656, completed cleanly, still water). Same 125-box grid/domain/gravity, but: two CylinderPillars
# (radius 12, X=138, Y=50/110) added with N/S walls (periodicity Y turned off); PID target raised
# 0 -> 0.01 (real streamwise current, matching this project's other "two pillars slitpass" scenarios);
# floaters raised +30 LU in Z and their X-columns compressed to 12/34/56/78/100 so all 125 sit entirely
# upstream (left) of the pillar footprint (X=126-150); all floaters startInactive+settlingSteps 20000
# (dormant until the current has 20000 steps to develop), timesteps raised 100000 -> 120000. See
# prms/box_sedimentation_showcase_analytical_125floaters_realgravity_twopillars_flow_gas.prm's header for
# the full derivation. Domain/decomposition unchanged: <160,160,180> / cellsPerBlock <40,40,20> = 144 procs.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 144 \
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_125floaters_realgravity_twopillars_flow_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
# `if` (not `[[ ]] &&`) so a missing directory's false test doesn't leak into this being the
# script's own final exit status -- bash's `for` returns its last command's exit code, and a bare
# `[[ -d ]] && mv` as the LAST statement of the LAST iteration (floater-history, only created when
# logFloaterHistory is on) would make the whole job report FAILED/ExitCode 1 to SLURM even when the
# simulation itself completed cleanly (confirmed happening exactly this way for job 4104656).
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"; fi
done
