#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-1000floaters-realgravity-density05-75k-collapsed-v3-lastbatch-floaterlog
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3_lastbatch_floaterlog_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3_lastbatch_floaterlog_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# logFloaterHistory sibling of job 4127217 (plain v3), restricted via the new
# ModelParameters.logFloaterHistoryReleaseStepFilter=9000 to just the 100 floaters released at
# step 9000 (the last of 10 release batches) -- that batch alone is the single heaviest anomaly
# contributor (17/81) in job 4127217's own statistics. Every diagnostic in that run only ever
# caught a floater's launch AFTER it was already airborne (max Z-range within any single
# MEM-WATCH/DEM-WATCH burst, across all ~1000 floaters: 1.36 LU) -- this run instead gets a
# magnitude-independent, per-subcycle position+velocity(+force/torque/hydroForce/density)
# snapshot for the whole run, for just that one batch, cheap enough (~1/10th the floater count of
# the 126GB all-floaters attempt that blew quota) to actually show a climb happening instead of
# just its aftermath. Classic VTK output (fluid_field/particles/floaterBox) is disabled in the prm
# (writeFrequency 1e10) since it isn't needed and the point is to save quota for the CSV log.
# See the prm's own header for the full rationale.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3_lastbatch_floaterlog_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
# `if` (not `[[ ]] &&`) so a missing directory's false test doesn't leak into this being the
# script's own final exit status -- see the 125-floater job's own comment (job 4104656) for the
# full mechanism this avoids.
for dir in vtk-out eval snapshots floater-history; do
    if [[ -d "$BIN_DIR/$dir" ]]; then mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"; fi
done
