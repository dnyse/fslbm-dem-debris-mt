#!/bin/bash -l
#SBATCH --job-name=antidunes-dynamich-ab
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=10:00:00
#SBATCH --output=sbatch-output/output/antidunes_dynamich_ab_%j.out
#SBATCH --error=sbatch-output/error/antidunes_dynamich_ab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# DYNAMIC-H re-run, now explicitly framed as a CLOGGING test (not a free-float-equilibrium test): the
# pillar gap (14 LU clear) is narrower than the sphere (20 LU diameter) — it geometrically cannot pass
# through, only wedge/jam against both pillars. Two-pillar floater (r=10 sphere, densityRatio 0.5,
# released at step 30000 into developed U=0.02 flow with the pile-up), run twice on ONE build:
#   A) BASELINE   static flat H (sphere_drop_two_pillars_gas.prm)
#   B) DYNAMIC-H  useDynamicReferenceHeight ON (sphere_drop_two_pillars_dynamich_gas.prm)
# B's MEM reference height is EMA-tracked live from the floaters' fitted local waterlines (force-side
# only, confirmed no PDF/fluid feedback) instead of held at the flat absoluteLiquidHeight — this cancels
# the uncancelled-pile-up phantom exactly, with no A_proj/taper geometry, unlike Option A (F_pred).
# usePhantomForceCorrection is OFF in both runs, to isolate dynamic-H from Option A.
#
# What to look at (this is a diagnostic run, not a pass/fail test):
#   - [DynamicH] in RUN B: H(EMA) should climb off the flat ~56 baseline toward the pillar pile-up
#     level (~62-65) as the sphere approaches, then track it.
#   - [FloaterForces] F_hydro=(...) z-component in A vs B at matched submersion: B should have
#     LESS of the open-cap upward phantom than A.
#   - [FloaterForces] F_contact+lub=(...) z-component, BOTH runs: from the phantom.log review, this is
#     nonzero and persistently DOWNWARD once the sphere wedges against the pillars (separation ~0 —
#     sphere too big for the gap) — likely friction from being jammed, not hydrodynamics. Compare its
#     magnitude/sign between A and B to see whether cancelling the hydrodynamic phantom changes how
#     firmly/where the sphere clogs.
#   - [phantomCorr] is absent in both (usePhantomForceCorrection off) — the new F_pred_raw / F_pred_applied
#     saturation-flagging log only fires when that toggle is on; not exercised by this job.
# grep -E 'RUN A|RUN B|FloaterForces|FloaterBuoyancy|DynamicH|MEM reference height' on the .out.
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

run_case() {  # $1 = banner, $2 = prm, $3 = label
    echo "############################################################"
    echo "## $1"
    echo "############################################################"
    "$SH" --run --np 72 --prm "$PRM_DIR/$2"
    archive_outputs "$3"
}

run_case "RUN A: BASELINE  (static flat H)"        sphere_drop_two_pillars_gas.prm             baseline
run_case "RUN B: DYNAMIC-H (EMA-tracked local H)"  sphere_drop_two_pillars_dynamich_gas.prm    dynamich

echo "## DONE. grep -E 'RUN A|RUN B|FloaterForces|FloaterBuoyancy|DynamicH' on this .out."
