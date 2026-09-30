#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-righting-jitter-s16
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=65
#SBATCH --time=02:30:00
#SBATCH --output=sbatch-output/output/bogner_righting_jitter_s16_%j.out
#SBATCH --error=sbatch-output/error/bogner_righting_jitter_s16_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-20] Thesis item 1a (THESIS_PLAN.md): SUB-CELL POSITION JITTER of the two-box fixed-heel MEM-torque scan (s=16). The two-box scan (4249061) found the held-box MEM torque is a static jagged
# staircase (rms 0.14-0.17 LU vs an ideal signal of 0.04-0.19 LU). Here the same two boxes (rho 0.75/0.25) at four heel angles (25.5, 26.565, 29, 30.5 deg) are shifted by (dz,dx) in {(-0.25,0),(+0.25,0),(0,+0.25),
# (+0.25,+0.25)} cells (same shift for both boxes): 16 runs of 20000 steps (~4.5 min each => ~1.3 h), each archived right after it runs, a dead run retried once. Prms: prms/bogner_righting_jitter_s16_<tag>_gas.prm
# (gen_bogner_righting_jitter_prm.py; ideal reference at the SHIFTED heave in prms/bogner_righting_jitter_s16_manifest.json). MEM-only buoyancy, NO analytic buoyancy. 2.5 h requested.
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
        ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_righting_jitter_s16"; mkdir -p "$ARCHIVE_ROOT"
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

TAGS=(a25p5_zm25 a25p5_zp25 a25p5_xp25 a25p5_zp25_xp25 a26p565_zm25 a26p565_zp25 a26p565_xp25 a26p565_zp25_xp25 a29_zm25 a29_zp25 a29_xp25 a29_zp25_xp25 a30p5_zm25 a30p5_zp25 a30p5_xp25 a30p5_zp25_xp25)
FAILED=()
for tag in "${TAGS[@]}"; do
    ok=0
    for attempt in 1 2; do
        echo "############################################################"
        echo "## righting scan JITTER s=16, ${tag}  (attempt ${attempt})"
        echo "############################################################"
        set +e
        "$SH" --run --np 65 --prm "$PRM_DIR/bogner_righting_jitter_s16_${tag}_gas.prm"
        rc=$?
        set -e
        archive_run "$tag" "$attempt"
        if [[ $rc -eq 0 ]]; then ok=1; break; fi
        echo "## ${tag} attempt ${attempt} exited with status ${rc}" >&2
    done
    [[ $ok -eq 1 ]] || FAILED+=("$tag")
done

choose_archive_root
{ echo "scan finished $(date)"; echo "angles: ${TAGS[*]}"; echo "failed after 2 attempts: ${FAILED[*]:-none}"; } | tee "$ARCHIVE_ROOT/SCAN_SUMMARY.txt"
echo "## DONE. Archived under $ARCHIVE_ROOT/<tag>/{eval,vtk-out,snapshots,floater-history}"
[[ ${#FAILED[@]} -eq 0 ]]
