#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-gravity-comparison-3x-d075
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_cuboid_stillwater_gravity_comparison_3xsphere_d075_%j.out
#SBATCH --error=sbatch-output/error/antidunes_cuboid_stillwater_gravity_comparison_3xsphere_d075_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Still-water two-way gravity comparison for a single free-floating densityRatio=0.5 cuboid
# (halfExtents <24,30,24>, 48x60x48, Y-elongated -- 48-LU cross-section matches the 3x sphere's
# diameter, 60-LU length matches the cylinder sibling), started at initialTiltAngleDeg=5 as a
# small-angle stability probe (tiltY logged every evaluationFrequency by [BoxBuoyancy]). Bogner
# three-phase MEM buoyancy throughout, run back to back in a single job on one build:
#   1. cuboid_stillwater_realgravity981_3xsphere_d075_gas.prm -- gravityOverride=5.06585e-05,
#      solved to give true g_SI=9.81 m/s^2 (same value as the sphere family, same dx_SI/dt_SI)
#   2. cuboid_stillwater_bognerboxgravity_3xsphere_d075_gas.prm -- Bogner & Rude (2013)'s own
#      g=7.5e-4, tau=1.9, unscaled
# Domain <96,135,120>, same as the 3x sphere family, 50000 timesteps.
# domainSize <96, 135, 120> / cellsPerBlock <16, 45, 30> = 6x3x4 = 72 processes

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

if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi

# --bed: spheres_out.dat's header is read unconditionally at startup (getAvgDiameterScalingFactor,
# Utility.h) for the SI-unit scaling, even though useParticleBed=false for these floater-only runs.
"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

archive_outputs() {  # $1 = label
    for d in vtk-out eval snapshots; do
        [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/$1" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$1/$d"
    done
}

run_case() {  # $1 = banner, $2 = prm, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "RUN real gravity g_SI=9.81 (gravityOverride=5.06585e-05)" \
    cuboid_stillwater_realgravity981_3xsphere_d075_gas.prm real_gravity_981
run_case "RUN Bogner box gravity (g=7.5e-4, tau=1.9)" \
    cuboid_stillwater_bognerboxgravity_3xsphere_d075_gas.prm bogner_box_gravity

echo "## DONE. Archived under $ARCHIVE_ROOT/{real_gravity_981,bogner_box_gravity}_{vtk-out,eval,snapshots}"
