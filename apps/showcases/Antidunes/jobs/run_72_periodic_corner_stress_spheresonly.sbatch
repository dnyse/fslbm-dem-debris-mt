#!/bin/bash -l
#SBATCH --job-name=antidunes-periodic-corner-stress-spheresonly
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_periodic_corner_stress_spheresonly_%j.out
#SBATCH --error=sbatch-output/error/antidunes_periodic_corner_stress_spheresonly_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Shape-independence stress test for the mesa_pd sync desync chased in
# eighteen_stacked_release_pillar_bogner_gas.prm (see that file and
# periodic_corner_stress_spheresonly_gas.prm's own header for the full history). Three
# free-floating spheres only, no boxes/cylinders, no pillars, deliberately pushed to lap the
# periodic-X domain repeatedly through the exact corner (periodic seam x=0/128 + Z-block boundary
# z=60) implicated in two of four mixed-shape crash repros. If this crashes with the same "invalid
# notification type" abort, the bug is confirmed shape-independent.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

# --bed included even though useParticleBed=false for the floaters themselves: Antidunes reads
# spheres_out.dat's header unconditionally at startup (getAvgDiameterScalingFactor(), Utility.h)
# to derive the SI-unit scaling for the WHOLE run, so a fresh bed must exist in this clone's
# build dir regardless of whether the bed is physically used.
"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --bed --bed-np 72 \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/periodic_corner_stress_spheresonly_gas.prm"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
