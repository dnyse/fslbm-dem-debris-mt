#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-1000floaters-realgravity-floaterhistory
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_floaterhistory_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_floaterhistory_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Diagnostic rerun of run_144_box_sedimentation_showcase_analytical_1000floaters_realgravity.sbatch
# with ModelParameters.logFloaterHistory turned on and VTK "density" added to fluid_field writers
# (see the prm's own header note). This scenario has crashed with the same [LOW-RHO]/invalid-
# quaternion signature in jobs 4106364 and 4112811, always within ~650-840s (well under 1h), so
# the time limit here is trimmed to 2h -- this run is expected to crash on its own, same as those
# two, and the floater-history CSVs + density-enabled VTK dumps are what we actually want out of it.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_1000floaters_realgravity_floaterhistory_gas.prm"

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
