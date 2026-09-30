#!/bin/bash -l
#SBATCH --job-name=antidunes-bogner-righting-scan2-s32
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=130
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/bogner_righting_scan2_s32_%j.out
#SBATCH --error=sbatch-output/error/bogner_righting_scan2_s32_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-20] Thesis item 1a (THESIS_PLAN.md): the two-box fixed-heel MEM-torque scan at TWICE the resolution (s=32, basin 260x80x80, g scaled to 3.75e-4, tau=1/1.9), reduced angle set
# {22.5, 25, 26.565, 28, 30} deg, rho 0.75 and 0.25 held at their ideal heaves, raw MEM torque from the floater-history CSV (window 39000..39300 of 40000 steps), NO analytic buoyancy. Question: does the static torque
# scatter (rms 0.14-0.17 LU at s=16, vs an ideal signal 0.04-0.19 LU) shrink with resolution? 130 ranks on 2 nodes (cellsPerBlock <20,16,40> = 13x5x2). Time: UNVERIFIED -- the free s=32 run took 64 min for 100000 steps
# on 130 ranks (~0.038 s/step), so 40000 steps ~ 26 min per angle, ~2.2 h for five, more with retries; 6 h requested. Each run's output is archived right after it (history CSVs are big); a dead run is retried once.
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
        ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_bogner_righting_scan2_s32"; mkdir -p "$ARCHIVE_ROOT"
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

TAGS=(a22p5 a25 a26p565 a28 a30)
FAILED=()
for tag in "${TAGS[@]}"; do
    ok=0
    for attempt in 1 2; do
        echo "############################################################"
        echo "## righting scan s=32, ${tag}  (attempt ${attempt})"
        echo "############################################################"
        set +e
        "$SH" --run --np 130 --prm "$PRM_DIR/bogner_righting_scan2_s32_${tag}_gas.prm"
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
