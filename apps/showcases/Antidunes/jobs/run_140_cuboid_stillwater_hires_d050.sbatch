#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-hires-d050
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --exclude=f0220,f0258
#SBATCH --ntasks=140
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_hires_d050_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_hires_d050_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 1a (THESIS_PLAN.md) -- higher-resolution, longer-duration rerun of the cuboid
# equilibrium (buoyancy) sweep, densityRatio=0.5 case. See
# cuboid_stillwater_realgravity981_3xsphere_hires_gas.prm's own header for the full resolution
# derivation (domain <119,168,150>, cellsPerBlock <119,7,25> = 144 ranks (full 2-node allocation), ~1.24x finer than the
# original <96,135,120>/<16,45,30>=72-rank run, timesteps 50000->100000).
#
# One of 3 density cases (0.25/0.5/0.75), each its own job (not combined) because all 3 sequentially
# would exceed the multinode partition's 24h MaxTime (~9-9.5h/case measured-derived estimate * 3 ~=
# 28h). This one (d050) runs first; d050 and d075 are submitted with --dependency=afterok on this
# job's ID so they queue safely behind it on the SAME checkout without a build-directory collision,
# rather than requesting 3 separate mirrors (only one of the 3 -- this "mt-dennys" checkout -- was
# free when submitted; mirror/mirror2 both had other jobs in flight).
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP" 2>/dev/null; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
    mkdir -p "$ARCHIVE_ROOT"
fi

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## cuboid stillwater hires, densityRatio=0.5"
echo "############################################################"
"$SH" --run --np 140 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_hires_gas.prm"

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/d050" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/d050/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/d050/{eval,vtk-out,snapshots}"
