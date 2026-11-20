#!/bin/bash -l
#SBATCH --job-name=antidunes-bedgen-buoyancy-r2
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=00:30:00
#SBATCH --output=sbatch-output/output/bedgen_buoyancy_r2_%j.out
#SBATCH --error=sbatch-output/error/bedgen_buoyancy_r2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Scaling-reference bed for the R = 2 level of the MEM vs. analytical buoyancy comparison
# (prms/BedGeneration_buoyancy_r2.prm, writes spheres_out_buoyancy_r2.dat). The R = 2 buoyancy prms read only the
# particle diameter and dx_SI from it (useParticleBed false), so its particles never enter the physics. Small
# and cheap: 72 ranks, 30 min. Uses antidunes.sh, which loads the modules itself.
# --constraint=el9 is omitted on purpose: every node of the singlenode partition carries that feature anyway.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72 --bed-prm "$PRM_DIR/BedGeneration_buoyancy_r2.prm"
