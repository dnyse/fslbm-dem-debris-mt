#!/bin/bash -l
#SBATCH --job-name=antidunes-18stacked-bogner-resume
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/antidunes_eighteen_stacked_bogner_resume_%j.out
#SBATCH --error=sbatch-output/error/antidunes_eighteen_stacked_bogner_resume_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Checkpoint-resume rerun of the mixed-shape 18-stack mesa_pd sync crash (see
# eighteen_stacked_release_pillar_bogner_gas_resume.prm's header for the full story): loads the
# checkpoint from job 3976230's last surviving snapshot (step 40001, archived under
# $HPCVAULT/antidunes-archive/antidunes-18stacked-bogner_20260807_011731/snapshots) instead of
# starting from t=0, so only the ~3000-step window up to the crash (~42890) needs to be
# re-simulated. Short --time budget accordingly.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/eighteen_stacked_release_pillar_bogner_gas_resume.prm"

# No output archiving needed here -- storeSnapshot is off, and vtk-out/eval from this short
# diagnostic window aren't needed for anything downstream.
