#!/bin/bash -l
#SBATCH --job-name=antidunes-mixed-shape-density-arrangement-sweep
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/antidunes_mixed_shape_density_arrangement_sweep_%j.out
#SBATCH --error=sbatch-output/error/antidunes_mixed_shape_density_arrangement_sweep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# DENSITY-ARRANGEMENT DIAGNOSTIC SWEEP -- four short (20000-step) runs on ONE build, back to
# back, to isolate exactly which density arrangement within a vertically-stacked floater group
# (sphere+box+cylinder sharing an X,Y position, Z=64/76/88) causes the DEM/MEM blowup seen in
# job 3858669 (mix_sed_pil.log, 54-floater SCRAMBLED grid: avgParticleVelocity non-finite
# ~step 1000, both the pillars and non-pillars variants crashed identically).
#
# Background: the ORIGINAL 9/81-floater "clean-factorial" grids (mixed_sed2.log job 3854507,
# mixed_sed_worked.log) stack all 3 shapes per (X,Y) too, but keep densityRatio CONSTANT across
# the whole stack (set by X-anchor only) -- both ran clean (20000 and 100000 steps respectively).
# The crashing 54-floater SCRAMBLED grid decorrelated density from shape+position, so a stack
# could contain 3 DIFFERENT densities, with the lightest (most buoyant/reboundy) object
# systematically at the BOTTOM of the stack, closest to the water, with denser objects still
# falling directly onto it once released. Since stacking itself is proven safe, this per-stack
# density mixing is the prime suspect -- these 4 runs test that directly:
#
#   uniform                    CONTROL: density constant per stack (the proven-working design).
#                               Expect: no crash.
#   mixed_light                REPRODUCTION: density mixed per stack, lightest at the bottom
#                               (the exact crashing pattern). Expect: reproduces the crash.
#   mixed_heavy                ORDERING TEST: same density mix, reversed (heaviest at bottom).
#                               Distinguishes "light object rebounds into falling heavy ones"
#                               from "any density mismatch in a stack is enough".
#   same_shape_mixed_density   SHAPE-ISOLATION TEST: all spheres (best-conditioned DEM contact
#                               case, no cross-shape geometry), same light-at-bottom density mix.
#                               Distinguishes density-mismatch-alone from a needed shape-mismatch
#                               contribution.
#
# ParticleKernels.h also gained unconditional [DEM-WATCH]/[MEM-WATCH] force-magnitude logs (every
# DEM sub-cycle, threshold ~0.05 vs. the ~1e-3..1e-2 magnitudes seen in stable runs) so that IF a
# variant reproduces the crash, the .out can be grepped to tell a DEM contact-resolution blowup
# (F_contact+lub growing across subCycle=0..9 within one step) apart from a fluid momentum-
# exchange (MEM) one (F_hydro already large at subCycle=0, before any DEM contact that step could
# have caused it).
#
# No pillars, no flow forcing (still water) in any of the 4 -- isolates the density-arrangement
# variable only. domainSize <128, 60, 180> / cellsPerBlock <32, 20, 30> = 4x3x6 = 72 processes.
#
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

run_case() {  # $1 = banner, $2 = prm, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "RUN uniform: CONTROL, density constant per stack (expect: no crash)" \
    mixed_shape_density_sweep_uniform_gas.prm uniform
run_case "RUN mixed_light: REPRODUCTION, lightest object at bottom of stack (expect: crash)" \
    mixed_shape_density_sweep_mixed_light_gas.prm mixed_light
run_case "RUN mixed_heavy: ORDERING TEST, heaviest object at bottom of stack" \
    mixed_shape_density_sweep_mixed_heavy_gas.prm mixed_heavy
run_case "RUN same_shape_mixed_density: SHAPE-ISOLATION TEST, all spheres, mixed density per stack" \
    mixed_shape_density_sweep_same_shape_mixed_density_gas.prm same_shape_mixed_density

echo "## DONE. grep -E 'RUN |DEM-WATCH|MEM-WATCH|non-finite|corrupted|invalid next size' on this .out/.err"
echo "## to see which variant(s) crashed and, if so, whether DEM-WATCH or MEM-WATCH fired first."
