#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-righting-scan2-s16
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/bogner_righting_scan2_s16_%j.out
#SBATCH --error=sbatch-output/error/bogner_righting_scan2_s16_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-19] Thesis item 1a (THESIS_PLAN.md): TWO-BOX fixed-heel righting-moment scan (rho 0.75 and 0.25 only), raw MEM torque, NO analytic buoyancy, s = 16.
# Follow-up to job 4248900 (three boxes), where heel 0/5/10 ran but 15 died twice (NaN ~1150 steps in) and 20 was failing: the rho=0.5 box (largest torque, largest disturbance) is dropped
# to test whether fewer boxes make the start-up survive. DIAGNOSTIC FIRST: a15 and a20 run first; if BOTH fail twice the job stops (the start-up problem is not fixed) instead of burning the
# rest of the angle list. Otherwise it continues with the fine scan around the equilibrium: 16.5, 18, 19.5, 21, 22.5, 24, 25.5, 26.565, 27.5, 29, 30.5, 32 deg.
# Prms: prms/bogner_righting_scan2_s16_a*_gas.prm (gen_bogner_righting_scan_prm.py --boxes 0.75,0.25 --prefix bogner_righting_scan2), ideal reference in
# prms/bogner_righting_scan2_s16_manifest.json, analysis: analyze_bogner_righting_scan.py <archive_dir> --manifest prms/bogner_righting_scan2_s16_manifest.json.
# Ideal theory: the rho=0.25 and rho=0.75 torque curves are IDENTICAL. The three-box scan showed the measured MEM torque scatters by ~0.35-0.45 LU around the ideal (2-5x the ideal
# torque near the equilibrium): this scan resolves that staircase. Output per run is MOVED to the archive right after the run (history CSVs are big); a dead run is retried once.
# Time: UNVERIFIED. 14 runs x (20000 steps ~4.5 min + startup) ~ 1.2 h, plus retries; 3 h requested.
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
SH="$REPO/apps/showcases/Antidunes/antidunes.sh"
PRM_DIR="$REPO/apps/showcases/Antidunes/prms"
BIN_DIR="$REPO/build/apps/showcases/Antidunes"
OUTDIRS=(vtk-out eval snapshots floater-history)

# Never DELETE a previous run's output (chained jobs share this checkout): move leftovers aside.
RESCUE="$BIN_DIR/archive/rescued_before_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)"
for dir in "${OUTDIRS[@]}"; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        mkdir -p "$RESCUE"
        echo "WARNING: leftover $BIN_DIR/$dir from a previous run -- moving to $RESCUE (NOT deleting)" >&2
        mv "$BIN_DIR/$dir" "$RESCUE/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

# Archive root: prefer $HPCVAULT, decided and verified once at first use (right before the first archive), else a LOCAL archive as a last resort.
STAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_ROOT=""
choose_archive_root() {
    [[ -n "$ARCHIVE_ROOT" ]] && return 0
    local cand="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
    if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$cand" 2>/dev/null && touch "$cand/.write_test" 2>/dev/null; then
        rm -f "$cand/.write_test"; ARCHIVE_ROOT="$cand"
    else
        ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_righting_scan2_s16"; mkdir -p "$ARCHIVE_ROOT"
        echo "WARNING: \$HPCVAULT unavailable -- LOCAL archive at $ARCHIVE_ROOT (move off manually)" >&2
    fi
}
archive_run() {  # $1 = tag, $2 = attempt
    choose_archive_root
    local dest="$ARCHIVE_ROOT/$1"; [[ "$2" != "1" ]] && dest="${dest}_attempt$2"
    mkdir -p "$dest"
    du -sh "$BIN_DIR/floater-history" 2>/dev/null | sed "s/^/floater-history size for $1: /" || true
    for d in "${OUTDIRS[@]}"; do [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$dest/$d"; done
    return 0
}

TAGS=(a15 a20 a16p5 a18 a19p5 a21 a22p5 a24 a25p5 a26p565 a27p5 a29 a30p5 a32)
DIAG=(a15 a20)
FAILED=()
for tag in "${TAGS[@]}"; do
    ok=0
    for attempt in 1 2; do
        echo "############################################################"
        echo "## righting scan s=16, ${tag}  (attempt ${attempt})"
        echo "############################################################"
        set +e
        "$SH" --run --np 65 --prm "$PRM_DIR/bogner_righting_scan2_s16_${tag}_gas.prm"
        rc=$?
        set -e
        archive_run "$tag" "$attempt"
        if [[ $rc -eq 0 ]]; then ok=1; break; fi
        echo "## ${tag} attempt ${attempt} exited with status ${rc}" >&2
    done
    [[ $ok -eq 1 ]] || FAILED+=("$tag")
    # diagnostic gate: after the two diagnostic angles, stop if BOTH failed twice (start-up problem not fixed)
    if [[ "$tag" == "${DIAG[1]}" ]]; then
        nd=0; for d in "${DIAG[@]}"; do for f in "${FAILED[@]:-}"; do [[ "$f" == "$d" ]] && nd=$((nd+1)); done; done
        if [[ $nd -ge ${#DIAG[@]} ]]; then echo "## DIAGNOSTIC GATE: both ${DIAG[*]} failed twice -- stopping the scan" >&2; break; fi
    fi
done

choose_archive_root
{ echo "scan finished $(date)"; echo "angles: ${TAGS[*]}"; echo "failed after 2 attempts: ${FAILED[*]:-none}"; } | tee "$ARCHIVE_ROOT/SCAN_SUMMARY.txt"
echo "## DONE. Archived under $ARCHIVE_ROOT/<tag>/{eval,vtk-out,snapshots,floater-history}"
[[ ${#FAILED[@]} -eq 0 ]]
