#!/bin/bash -l
#SBATCH --job-name=antidunes-18stacked-bogner-spheresonly
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_eighteen_stacked_bogner_spheresonly_%j.out
#SBATCH --error=sbatch-output/error/antidunes_eighteen_stacked_bogner_spheresonly_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Shape-isolation test for the eighteen-stack mesa_pd sync crash (ParseMessage.h "Received
# invalid notification type", first seen in job 3956749, reproduced again in 3963620 after
# raising the floaters out of the water -- see eighteen_stacked_release_pillar_bogner_gas.prm's
# header). This run uses eighteen_stacked_release_pillar_bogner_gas_spheresonly.prm: identical
# 18-position layout, but every floater is a FloatingSphere instead of the mixed sphere/box/
# cylinder Latin square. If this run completes cleanly while the mixed-shape run (running
# concurrently on another clone, instrumented with ParseMessage.h's new [PARSE-HISTORY] ring
# buffer) crashes again, that implicates the box/cylinder shape-specific code path rather than a
# generic particle-count/contact-density race.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/eighteen_stacked_release_pillar_bogner_gas_spheresonly.prm"

# archive this run's own output so the next job doesn't overwrite it
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
