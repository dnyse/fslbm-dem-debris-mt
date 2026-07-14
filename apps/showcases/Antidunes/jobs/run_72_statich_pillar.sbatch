#!/bin/bash -l
#SBATCH --job-name=antidunes-statich-pillar
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_statich_pillar_%j.out
#SBATCH --error=sbatch-output/error/antidunes_statich_pillar_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# build + run: two-pillar clogging scenario with a FLAT static MEM reference height H pinned to the
# empirically observed pillar pile-up level (H=64, vs the domain-wide flat absoluteLiquidHeight
# ~56.02) instead of live dynamic-H tracking. Single-variable test: does simply picking the right
# CONSTANT already fix most of the problem, before spending more effort on dynamic-H's self-pollution
# or the wide-net plane-fit's own blending risk. usePhantomForceCorrection and
# useDynamicReferenceHeight both stay OFF. timesteps 100000 (30000 held + 70000 free), matching the
# wide-net dynamic-H run for a fair comparison.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_drop_two_pillars_statich_pillar_gas.prm"

# archive this run's own output so the next job doesn't overwrite it
STAMP="$(date +%Y%m%d_%H%M%S)"
# if $HPCVAULT is set, archive to persistent vault storage; otherwise the local build dir
if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_DIR="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_$STAMP"
else
    ARCHIVE_DIR="$BIN_DIR/archive/$STAMP"
fi
mkdir -p "$ARCHIVE_DIR"
for dir in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$dir" ]] && mv "$BIN_DIR/$dir" "$ARCHIVE_DIR/$dir"
done
