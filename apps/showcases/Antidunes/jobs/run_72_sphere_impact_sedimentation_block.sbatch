#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-impact-sedimentation-block
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_impact_sedimentation_block_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_impact_sedimentation_block_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Replication of Bogner (2017 dissertation) Fig 5.11: a densely packed
# 10x10x10 (1000-sphere) cubic block, three densities (2/1/0.5) arranged in a
# diagonal-plane weave -- (i+j+k)%3 -- matching the figure's visible banding,
# released together into a walled tank of quiescent liquid, real 9.81 m/s^2
# gravity, 100000 timesteps. Full run following the successful 20000-step
# exploratory pass (job 3999350, completed clean in 1h13m50s at ~4.15
# steps/sec once settled) -- 10h walltime budgeted off that observed rate
# (100000 steps * ~0.2215 sec/step =~ 6.15h, plus margin for possibly slower
# throughput while the pile stays contact-dense on the floor for longer).
# See prm file header for full geometry/gravity derivation.
#
# set -e so a fatal in-simulation crash actually fails this SLURM job instead
# of being masked by the archiving step's own clean exit.
set -euo pipefail

# domainSize <150, 141, 240> / cellsPerBlock <50, 47, 30> = 3x3x8 = 72 blocks
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
    --prm "$REPO/apps/showcases/Antidunes/prms/sphere_impact_sedimentation_block_gas.prm"

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
