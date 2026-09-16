#!/bin/bash -l
#SBATCH --job-name=antidunes-sphere-stillwater-fixeddyn-density-sweep
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=02:00:00
#SBATCH --output=sbatch-output/output/antidunes_sphere_stillwater_fixeddyn_density_sweep_%j.out
#SBATCH --error=sbatch-output/error/antidunes_sphere_stillwater_fixeddyn_density_sweep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Still-water densityRatio sweep (0.25 / 0.5 / 0.75) for a single free-floating sphere, all three
# using the "fixed dynamic-value gravity" regime (gravityOverride=1.16888e-05, the lattice-unit
# value job 3947080's dynamic Fr/Re/liquidHeight computation produced for this domain -- see
# sphere_stillwater_fixed_dynamicgravity_gas.prm for the full rationale). Companion to the
# corresponding fix applied to nine_stacked_release_pillar_bogner_gas.prm: BedGeneration.prm's
# diameter_SI/domainSize_SI/generationSpacing_SI were rescaled (by 6.89286, matching the original
# antidunes branch's own 2.9mm-sphere physical scale) so spheres_out.dat's avgParticleDiameter --
# and hence liquidHeight, and hence the physical (SI) meaning of dx_SI/dt_SI -- now converts this
# same lattice gravity to a physically real ~9-ish m/s^2 instead of the ~0.02 m/s^2 the previous
# (20mm-sphere) BedGeneration values implied. IMPORTANT: this rescaling is exact-ratio-preserving
# (see this session's derivation), so avgParticleDiameter_lattice, liquidHeight_lattice, and hence
# every lattice-unit number in these three runs -- gravity_lattice, the whole submerged-fraction
# trajectory -- come out BIT-IDENTICAL to the earlier (uncorrected) runs of these same three prm
# files; only the logged dx_SI/dt_SI/g_SI (the SI-unit interpretation) change. This run exists to
# capture that corrected SI-unit context alongside a fresh, from-a-known-good-BedGeneration log.
#
# Requires spheres_out.dat to already be regenerated (from the rescaled BedGeneration.prm) in this
# clone's build directory -- run jobs/bed.sbatch first if not already done.
#
# domainSize <96, 45, 40> / cellsPerBlock <16, 15, 10> = 6x3x4 = 72 processes

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

if [[ -n "$HPCVAULT" ]]; then
    ARCHIVE_ROOT="$HPCVAULT/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$STAMP"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/$STAMP"
fi

[[ -f "$BIN_DIR/spheres_out.dat" ]] || { echo "spheres_out.dat not found in $BIN_DIR -- run jobs/bed.sbatch first"; exit 1; }

"$SH" --build --jobs "$(nproc)"

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

run_case "RUN densityRatio 0.25, fixed dynamic-value gravity (rescaled BedGeneration)" \
    sphere_stillwater_fixed_dynamicgravity_d025_gas.prm d025
run_case "RUN densityRatio 0.50, fixed dynamic-value gravity (rescaled BedGeneration)" \
    sphere_stillwater_fixed_dynamicgravity_gas.prm d050
run_case "RUN densityRatio 0.75, fixed dynamic-value gravity (rescaled BedGeneration)" \
    sphere_stillwater_fixed_dynamicgravity_d075_gas.prm d075

echo "## DONE. Archived under $ARCHIVE_ROOT/{d025,d050,d075}_{vtk-out,eval,snapshots}"
