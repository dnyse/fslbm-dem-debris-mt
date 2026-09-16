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

# Thesis item 1b (THESIS_PLAN.md) -- scaling-reference bed generation for the R=1 resolution level.
# See prms/BedGeneration_buoyancy_r2.prm's own header for the full derivation. Cheap, small run (this
# bed is a pure SI-scaling reference -- useParticleBed=false in every Antidunes prm here means its
# particles are never loaded into the actual physics). Uses antidunes.sh directly (not bypassed --
# it handles module loading internally), unlike DropTest's own sbatch scripts.
# --constraint=el9 deliberately omitted -- was rejecting ALL jobs cluster-wide with "Invalid feature
# specification" earlier this session (a scheduler-side issue, unrelated to this script); every node
# in the singlenode partition already carries that feature anyway (confirmed via sinfo).
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
