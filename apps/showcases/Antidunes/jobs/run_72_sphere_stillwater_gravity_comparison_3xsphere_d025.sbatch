#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-stillwater-gravity-comparison-3x-d025
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=08:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_stillwater_gravity_comparison_3xsphere_d025_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_stillwater_gravity_comparison_3xsphere_d025_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Still-water two-way gravity comparison for a single free-floating densityRatio=0.5 sphere,
# radius 24 (3x the original 8-radius family), Bogner three-phase MEM buoyancy throughout, run
# back to back in a single job on one build:
#   1. sphere_stillwater_realgravity981_3xsphere_d025_gas.prm -- corrected gravityOverride=5.06585e-05,
#      solved to give true g_SI=9.81 m/s^2 against the current (rescaled) BedGeneration.prm's
#      dx_SI/dt_SI (the original "fixed dynamic" g=1.16888e-05 only converted to ~2.26 m/s^2 --
#      see that file's header for the derivation)
#   2. sphere_stillwater_bognerboxgravity_3xsphere_d025_gas.prm -- Bogner & Rude (2013)'s own
#      g=7.5e-4, tau=1.9, unscaled
# Domain <96,135,120> (9x cells vs the original 8-radius family -- X kept unchanged, only Y/Z
# grown, see either prm's header for the full derivation), 50000 timesteps (up from 30000).
# domainSize <96, 135, 120> / cellsPerBlock <16, 45, 30> = 6x3x4 = 72 processes

# Resolve the repo root programmatically (mirrors antidunes.sh's own SCRIPT_DIR/REPO_ROOT
# derivation) -- no hardcoded path, so this works across mirror clones (mt-dennys,
# mt-dennys-mirror, mt-dennys-mirror2, ...) regardless of which one submitted the job.
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
    sphere_stillwater_realgravity981_3xsphere_d025_gas.prm real_gravity_981
run_case "RUN Bogner box gravity (g=7.5e-4, tau=1.9)" \
    sphere_stillwater_bognerboxgravity_3xsphere_d025_gas.prm bogner_box_gravity

echo "## DONE. Archived under $ARCHIVE_ROOT/{real_gravity_981,bogner_box_gravity}_{vtk-out,eval,snapshots}"
