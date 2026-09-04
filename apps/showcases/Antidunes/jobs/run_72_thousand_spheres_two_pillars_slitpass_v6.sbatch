#!/bin/bash -l
#SBATCH --job-name=antidunes-thousand-spheres-two-pillars-slitpass-v6
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=12:00:00
#SBATCH --output=sbatch-output/output/antidunes_thousand_spheres_two_pillars_slitpass_v6_%j.out
#SBATCH --error=sbatch-output/error/antidunes_thousand_spheres_two_pillars_slitpass_v6_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# v6: same 1000-sphere/two-pillar slit-pass scenario, but testing the SIZE hypothesis for why
# v5 (radius=2.0, job 4060371) ended with 983/993 floaters sunk below the water layer entirely.
# radius 2.0->4.0 (4->8 lattice cells across the diameter, matching this project's own validated
# MEM-only-buoyancy showcases much more closely), domain grown <160,75,180> -> <256,120,240> to fit
# the bigger spheres with adequate clearance (avgParticleDiameter/floorHeight/absoluteLiquidHeight
# rescale accordingly -- hand-estimated in the prm header, to be cross-checked against the real
# runtime-printed eval/computedParams.txt), release layers raised further above the water and spaced
# further apart, first release delayed from settlingSteps=300 to 2000, pillars scaled up
# proportionally with the domain. See prms/thousand_spheres_two_pillars_slitpass_v6_gas.prm's header
# for the full derivation of every number.
# domainSize <256, 120, 240> / cellsPerBlock <64, 40, 40> = 4x3x6 = 72 processes (same convention)
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

"$REPO/apps/showcases/Antidunes/antidunes.sh" \
    --build --jobs "$(nproc)" \
    --run --np 72 \
    --prm "$REPO/apps/showcases/Antidunes/prms/thousand_spheres_two_pillars_slitpass_v6_gas.prm"

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
