#!/usr/bin/env python3
"""Bogner-like FREE floating box, three densities (0.25 / 0.5 / 0.75), DEEP basin, 150k steps, checkpointing ON (THESIS_PLAN.md 1a, 2026-09-28).
Follows the earlier gen_bognerlike_box_scan.py set-up (real g_SI=9.81, box s=18 -> 18x36x18, long axis Y, heave + X/Y translation FREE, pitch/yaw pinned via
pinRotationAxes <1,0,1>, roll about Y free) with three changes: (1) basin <96,45,60> and liquid depth 30 (was 40 / 20) so a sinking rho 0.75 box cannot touch the
floor within the observed ~5-cell drop; (2) rho 0.5 added, ALL seeds 2.86 deg (Bogner's); (3) 150000 steps, storeSnapshot every 25000 steps (single rolling checkpoint, keepLastTwoSnapshots would need a per-generation resume prm) into a
run-specific folder, plus a matching *_resume prm (loadSnapshot true, timesteps 400000). The eval sqlite (totalFluidMass + maxFluidPos_z every 250 steps) is KEPT by
the job to test the water-loss hypothesis. Generates prms/bogner3_r{025,050,075}_deep{,_resume}_gas.prm."""
import math, re, sys
sys.path.insert(0, ".")
from gen_bognerlike_box_scan import square, clip_below, area
RHOS = [0.25, 0.5, 0.75]; SEED = 2.86; S = 18.0; AVGD = 6.00146; FLOOR = 0.0600146; LIQ_FACTOR = 5.0
WL0 = LIQ_FACTOR*AVGD + FLOOR   # 30.067 (absoluteLiquidHeight; the ref set had 3.3227*6.00146+0.06 = 20.0011)
STEPS = 150000; RESUME_STEPS = 400000
def heave(rho, alpha):
    lo, hi = WL0 - 2*S, WL0 + 2*S
    for _ in range(80):
        m = 0.5*(lo+hi); lo, hi = (m, hi) if area(clip_below(square(m, alpha), WL0)) > rho*S*S else (lo, m)
    return round(0.5*(lo+hi), 3)
BASE = open("prms/bognerlike_r075_h286_gas.prm", encoding="utf-8").read()
def prm(rho, resume, variant=""):
    tag = f"r{int(round(rho*100)):03d}"; zc = heave(rho, SEED); t = BASE
    def rep(pat, new):
        nonlocal t; t, n = re.subn(pat, new, t, count=1, flags=re.S); assert n == 1, pat
    rep(r"cellsPerBlock\s+<[^>]*>;", "cellsPerBlock                 < 16, 15, 15 >;")
    rep(r"loadSnapshot\s+\w+;", f"loadSnapshot                  {'true' if resume else 'false'};")
    rep(r"storeSnapshot\s+\w+;", "storeSnapshot                 true;")
    rep(r"snapshotFrequency\s+\d+;", "snapshotFrequency             25000;")
    rep(r"snapshotBaseFolder\s+\w+;", f"snapshotBaseFolder            snapshots_b3_{tag}{variant};")
    rep(r"domainSize\s+<[^>]*>;", "domainSize         <96, 45, 60>;")
    rep(r"liquidHeightFactor\s+[\d.]+;", f"liquidHeightFactor {LIQ_FACTOR};")
    if variant == "_legacyrec":
        rep(r"useBognerReconstruction\s+\w+;", "useBognerReconstruction       false;   // LEGACY uniform all-interface FormerMO reconstruction, no phi seeding (control)")
    if variant == "_nobogner":
        rep(r"useBognerReconstruction\s+\w+;", "useBognerReconstruction       false;   // ALL Bogner options OFF (control): legacy uniform all-interface FormerMO reconstruction")
        rep(r"useFillWeighting\s+\w+;", "useFillWeighting              false;")
        rep(r"useBognerThreePhaseMEM\s+\w+;", "useBognerThreePhaseMEM        false;   // rhoHydStat subtracted from the MEM force (pre-Bogner Antidunes behaviour)")
        rep(r"useMemOnlyBuoyancy\s+\w+;", "useMemOnlyBuoyancy            false;   // analytic vSub*g buoyancy + analytic righting torque")
    rep(r"timesteps\s+\d+;", f"timesteps              {RESUME_STEPS if resume else STEPS};")
    rep(r"logDensityColumnStartStep\s+\d+;", f"logDensityColumnStartStep     {STEPS-2000};")
    rep(r"useMemOnlyBuoyancy\s+\w+;", "useMemOnlyBuoyancy            true;\n   repairDegenerateObstaclePdfs  false;   // Option B (obstacle-cell repair guard) OFF\n   seedStreamPdfCloneOnFirstUse  false;   // Option A (scratch PDF clone seeding) OFF")
    if variant == "_bognergrav":
        rep(r"gravityOverride\s+[\d.e+-]+;[^\n]*", "gravityOverride        7.5e-4;        // Bogner and Rude 2013's own g (thesis 5.3.3)\n   relaxationRateOverride 1.9;           // Bogner and Rude 2013's tau = 1/1.9 = 0.526 (thesis 5.3.3)")
    rep(r"logDensityColumnEndStep\s+\d+;", f"logDensityColumnEndStep       {STEPS-1};\n\n   // fill-integrated water height on a 12 x 3 grid of columns (x = 4, 12, ..., 92; y = 7.5, 22.5, 37.5) every 250 steps over the WHOLE run\n   logWaterHeight                true;\n   waterHeightBaseFolder         water-height;\n   waterHeightGridX              12;\n   waterHeightGridY              3;\n   waterHeightEveryNSteps        250;")
    rep(r"CellBB_filter_0 \{[^}]*\}", "CellBB_filter_0 { min < 0, 22, 0 >;   max < 96, 23, 60 >; }")
    rep(r"CellBB_filter_1 \{[^}]*\}", "CellBB_filter_1 { min < 0, 0, 0 >;    max < 1, 45, 60 >; }")
    rep(r"CellBB_filter_2 \{[^}]*\}", "CellBB_filter_2 { min < 0, 0, 0 >;    max < 96, 45, 60 >; }")
    blk = f"""FloatingBoxes
{{
   FloatingBox
   {{
      // Bogner-like free box: s=18 (18 x 36 x 18), roll about Y only (pitch/yaw pinned), heave + X/Y translation FREE, heel seed {SEED} deg
      position             < 48, 22.5, {zc} >;
      halfExtents          < 9, 18, 9 >;
      densityRatio         {rho};
      settlingSteps        0;
      initialTiltAngleDeg  {SEED};
      comOffsetZ           0;
      pinPosition          false;
      pinRotation          true;
      pinRotationAxes      < 1, 0, 1 >;
   }}
}}"""
    rep(r"FloatingBoxes\n\{.*?\n\}\n", blk + "\n")
    hdr = (f"// bogner3_{tag}_deep{variant}{'_resume' if resume else ''}_gas.prm -- generated by gen_bogner3_deep_scan.py (do not hand-edit)\n"
           f"// [Claude, 2026-09-28] Bogner-like FREE box rho {rho}, seed {SEED} deg, real g_SI=9.81, DEEP basin <96,45,60> (liquid depth {WL0:.3f}), heave {zc}, "
           f"{'RESUME from the newest snapshots_b3_'+tag+' checkpoint up to '+str(RESUME_STEPS) if resume else str(STEPS)+' steps'}; single rolling checkpoint every 25000 steps.\n"
           f"// Free heave + X/Y translation, roll about Y only (pinRotationAxes <1,0,1>). Eval sqlite (totalFluidMass, maxFluidPos_z) every 250 steps is kept by the job.\n"
           f"// ---- original header follows ----\n")
    return hdr + t, zc, tag
SBATCH = """#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner3-deep
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=05:00:00
#SBATCH --output=sbatch-output/output/bogner3_deep_%j.out
#SBATCH --error=sbatch-output/error/bogner3_deep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-28] Bogner-like free box, deep basin, 150k steps, checkpointing (gen_bogner3_deep_scan.py). ONE density per job/checkout:
#   fresh : cd <checkout>/apps/showcases/Antidunes && sbatch --export=ALL,B3_RHO=075 jobs/run_72_bogner3_deep.sbatch      (B3_RHO = 025 | 050 | 075)
#   resume: cd <checkout>/apps/showcases/Antidunes && sbatch --export=ALL,B3_RHO=075,B3_RESUME=1 jobs/run_72_bogner3_deep.sbatch  (continues to 400000 from the newest checkpoint)
# Results in pinscan_results/<jobid>/ (never overwritten, backup in ~/pinscan_backup), snapshots kept in build/.../snapshots_b3_r<RHO>/ (NOT deleted), eval sqlite kept.
set -uo pipefail
RHO="${B3_RHO:?set B3_RHO=025|050|075 via --export=ALL,B3_RHO=...}"; RESUME="${B3_RESUME:-0}"; VAR="${B3_VARIANT:-}"   # VAR="" | "_legacyrec" | "_nobogner"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
OUT="$REPO/apps/showcases/Antidunes/pinscan_results/${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
SNAP="$BIN_DIR/snapshots_b3_r${RHO}${VAR}"
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }
BACKUP="$HOME/pinscan_backup/$(basename "$OUT")_bogner3_r${RHO}${VAR}"
if [[ "$RESUME" == "1" ]]; then
    PRM="bogner3_r${RHO}_deep${VAR}_resume_gas.prm"; LAB="b3_r${RHO}${VAR}_resume"
    [[ -f "$SNAP/simState.file" ]] || { echo "FATAL: no $SNAP/simState.file to resume from"; exit 1; }
    cp -r "$SNAP" "${SNAP}_backup_before_resume_${SLURM_JOB_ID:-manual}" || { echo "FATAL: checkpoint backup failed"; exit 1; }
else
    PRM="bogner3_r${RHO}_deep${VAR}_gas.prm"; LAB="b3_r${RHO}${VAR}"
    [[ -e "$SNAP" ]] && { echo "FATAL: $SNAP already exists -- refusing to start a fresh run over an existing checkpoint (move it away or resume)"; exit 1; }
fi
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

"$SH" --build --jobs "$(nproc)" --bed --bed-np 72

TMPLOG="$OUT/$LAB.fulllog.tmp"
"$SH" --run --np 72 --prm "$PRM_DIR/$PRM" > "$TMPLOG" 2>&1
echo "exit=$?" > "$OUT/$LAB.status"
grep -E "\\[BoxForces\\]|\\[BoxTorque\\]" "$TMPLOG" | cut -c1-460 > "$OUT/$LAB.forces.txt"
grep -E "Assertion|FATAL" "$TMPLOG" | head -5 >> "$OUT/$LAB.status"
grep -c "WARNING" "$TMPLOG" | sed 's/^/warnings=/' >> "$OUT/$LAB.status"
cat "$BIN_DIR"/density-column/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { s[$1]+=$8; r[$5]+=$7; c[$5]++ } END { for (k in s) { n++; t+=s[k] } if (n) printf "waterline_mean %.5f nsteps %d\\n", t/n, n; for (z in r) printf "rho %s %.6f\\n", z, r[z]/c[z] }' | sort -k1,1 -k2,2n > "$OUT/$LAB.column.txt"
cat "$BIN_DIR"/water-height/rank_*.csv 2>/dev/null | awk -F, '$1 ~ /^[0-9]+$/ { k=$1" "$3" "$4; h[k]+=$5; o[k]+=$6 } END { for (k in h) print k, h[k], o[k] }' | sort -k1,1n -k2,2n -k3,3n > "$OUT/$LAB.waterheight.txt"   # step gx gy height nObstacleCells
cp "$BIN_DIR/eval/computedParams.txt" "$OUT/$LAB.computedParams.txt" 2>/dev/null
cp "$BIN_DIR/eval/antidunes.sqlite" "$OUT/$LAB.antidunes.sqlite" 2>/dev/null
gzip -c "$TMPLOG" > "$OUT/$LAB.fulllog.gz" && rm -f "$TMPLOG"
rm -rf "$BIN_DIR/density-column" "$BIN_DIR/eval" "$BIN_DIR/vtk-out" "$BIN_DIR/water-height"

mkdir -p "$BACKUP" && cp -n "$OUT"/* "$BACKUP"/ && chmod -R a-w "$OUT" "$BACKUP"
echo "## DONE. Results in $OUT (backup: $BACKUP); checkpoints kept in $SNAP"
"""
if __name__ == "__main__":
    for r in RHOS:
        for resume in (False, True):
            text, zc, tag = prm(r, resume)
            open(f"prms/bogner3_{tag}_deep{'_resume' if resume else ''}_gas.prm", "w", encoding="utf-8").write(text)
        print(f"rho {r}: heave {zc} (level {WL0:.3f})")
    for var in ("_legacyrec", "_nobogner", "_bognergrav"):   # control variants, rho 0.75 only
        for resume in (False, True):
            text, zc, tag = prm(0.75, resume, var)
            open(f"prms/bogner3_{tag}_deep{var}{'_resume' if resume else ''}_gas.prm", "w", encoding="utf-8").write(text)
    open("jobs/run_72_bogner3_deep.sbatch", "w", encoding="utf-8").write(SBATCH)
