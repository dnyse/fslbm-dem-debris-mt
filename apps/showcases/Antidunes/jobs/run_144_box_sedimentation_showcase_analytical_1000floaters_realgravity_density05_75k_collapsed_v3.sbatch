#!/bin/bash -l
#SBATCH --job-name=antidunes-box-sedimentation-showcase-analytical-1000floaters-realgravity-density05-75k-collapsed-v3
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3_%j.out
#SBATCH --error=sbatch-output/error/antidunes_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# V3 of run_144_box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_125k_collapsed_v2.sbatch
# (job 4124938, cancelled at step 15743 -- ran clean, zero [LOW-RHO], but classic VTK output alone
# was on track to blow the 200GB hard home quota well before 125000 steps). V3 keeps V2's spawn
# height (Z=118) and release interval (1000 steps, all released by step 9000), drops "mapped_flag"
# from the VTK writers, raises VTK writeFrequency 1000 -> 2000, and shortens the run to 75000
# timesteps -- landing around ~140GB total, comfortably under quota. logFloaterHistory stays OFF
# (job 4124682 blew the 200GB hard quota with it on at step ~18254). See the prm's own header.
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
    --prm "$REPO/apps/showcases/Antidunes/prms/box_sedimentation_showcase_analytical_1000floaters_realgravity_density05_75k_collapsed_v3_gas.prm"

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
