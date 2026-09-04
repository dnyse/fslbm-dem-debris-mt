#!/bin/bash -l
#SBATCH --job-name=antidunes-pont-ferroviere-realscale
#SBATCH --partition=multinode
#SBATCH --constraint=el9
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=72
#SBATCH --time=12:00:00
#SBATCH --output=sbatch-output/output/antidunes_pont_ferroviere_realscale_%j.out
#SBATCH --error=sbatch-output/error/antidunes_pont_ferroviere_realscale_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Real-scale reconstruction of Pont ferroviaire, Rue Goffontaine, Trooz (river Vesdre, EPFL/ULiège
# bridge-clogging database #30071). Two steps: (1) BedGeneration against
# BedGeneration_pont_ferroviere.prm (produces spheres_out_pont_ferroviere.dat -- required even
# though useParticleBed=false, see that file's own header for why), then (2) the main run against
# prms/pont_ferroviere_rue_goffontaine_realscale_gas.prm (domain <2688,1056,432>,
# cellsPerBlock <112,44,54> = 24x24x8 = 4608 processes, 64 nodes -- ~33 LU/m resolution, 115.8
# cells across the water depth, exceeding even the original `antidunes` branch's own 81-cell
# reference. Time limit set to 12h; expected actual runtime ~6-7h for 50000 steps. See that prm's
# header for the full resolution-vs-node-count derivation).
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
    --bed --bed-np 72 --bed-prm "$REPO/apps/showcases/Antidunes/BedGeneration_pont_ferroviere.prm" \
    --run --np 4608 \
    --prm "$REPO/apps/showcases/Antidunes/prms/pont_ferroviere_rue_goffontaine_realscale_gas.prm"

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
