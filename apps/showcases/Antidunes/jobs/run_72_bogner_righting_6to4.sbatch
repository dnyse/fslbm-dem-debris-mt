#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-righting-6to4
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=12:00:00
#SBATCH --output=sbatch-output/output/antidunes_bogner_righting_6to4_%j.out
#SBATCH --error=sbatch-output/error/antidunes_bogner_righting_6to4_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# BOGNER & RÜDE 2013 Fig 5.6 — righting-moment (stability-curve) validation, ratio 6:4, COARSE resolution.
# A 6:4 cuboid (b=24 x, height h=16 z, length 2h=32 y) is held FIXED, half-immersed, at a sequence of heel
# angles a ∈ {0,5,10,15,20,25,30}° about its longitudinal y-axis, in still water (U=0, hardcoded units).
# Each run logs [BoxRightingMoment] heelDeg / M_y / GZ / V_sub — one point on the stability curve M(a).
# Analytic buoyancy ON, so the clip's exact submerged centroid drives the moment (converges to the ideal).
# Collect M(a) across the 7 runs and plot vs the ideal curve (and vs Bogner's 24dx line).
# NOTE: this is the COARSE (24 dx) resolution only; the 48 dx / 96 dx rows in Bogner need larger custom
# domains (box length 2h = 64/128 exceeds the y=60 of this domain) — set up separately.
# grep -E 'RUN a|BoxRightingMoment|BoxBuoyancy' on the .out.
# domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 blocks
# Resolve the repo root programmatically (mirrors antidunes.sh's own SCRIPT_DIR/REPO_ROOT
# derivation) -- no hardcoded path, so this works across mirror clones (mt-dennys,
# mt-dennys-mirror, mt-dennys-mirror2, ...) regardless of which one submitted the job. Jobs
# are always submitted with `sbatch` from apps/showcases/Antidunes (same dir as antidunes.sh),
# so SLURM_SUBMIT_DIR (guaranteed set by SLURM to that submission directory) sits 3 levels
# below the repo root; fall back to deriving it from this script's own location (four levels
# up from apps/showcases/Antidunes/jobs/) for direct/non-sbatch execution.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
STAMP="$(date +%Y%m%d_%H%M%S)"

"$SH" --build --jobs "$(nproc)"

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$BIN_DIR/archive/$STAMP" && mv "$BIN_DIR/$d" "$BIN_DIR/archive/$STAMP/${1}_$d"
    done
}

run_case() {  # $1 = heel angle, $2 = prm, $3 = label
    echo "############################################################"
    echo "## RUN a=$1 deg  (6:4 coarse, fixed half-immersed)"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

for a in 00 05 10 15 20 25 30; do
    run_case "$a" "bogner_righting_6to4_a${a}_gas.prm" "a${a}"
done

echo "## DONE. grep -E 'RUN a|BoxRightingMoment' on this .out; collect M_y vs heelDeg into the stability curve."
