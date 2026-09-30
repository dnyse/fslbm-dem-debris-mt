#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-extreme-pidon-fresh
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks=288
#SBATCH --time=04:00:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_extreme_pidon_fresh_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_extreme_pidon_fresh_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3/4 (THESIS_PLAN.md) -- clean, SINGLE-JOB fresh rerun of the extreme clog/PID-on case,
# superseding the earlier phase1+phase1.5 checkpoint-resume reconstruction (jobs 4235135+4237896).
# That split was only needed because the checkpoint-resume floater-mapping bug, the ghost-sync gap,
# and the framework logging slowdown were all still unfixed/undiscovered at the time. All three are
# now fixed, so a single fresh run of the real validation prm (not a capped/debug variant) should
# complete in one shot: 97500 steps at the confirmed post-fix rate (~66ms/step on this scenario
# family, see R3's own post-fix timing) is ~1.8h; 4h requested for margin.
#
# Also rerun (not just PID-off) because clog_measurement_validation_extreme_pidon_gas.prm's own
# release-timing ladder was widened today (200 -> 250 step gap) to match the fix applied to the
# PID-off sibling file after ITS release-overlap crash -- see both files' own 2026-09-16 header
# notes. Keeping both prms byte-identical except PIDParameters.enabled is the whole point of this
# comparison pair, so PID-on is rerun with the same corrected timing even though its original
# 200-step schedule never actually failed.
#
# storeSnapshot enabled (snapshotFrequency=5000, prm-level change, matching PID-off) as a safety net.
#
# [Claude, 2026-09-18]: RENAMED from run_144_..._fresh.sbatch (144 ranks / 2 nodes) to run_288_...
# (288 ranks / 4 nodes) -- the completed run under the OLD prm (job 4242919) showed
# clogUpstreamElevation pinned at the domain's own old ceiling (Z=90) for ~40% of the run, a hard
# box-height limitation rather than a converged value. domainSize's Z was doubled to 180 to fix this
# (see the prm's own header), which at the unchanged cellsPerBlock<64,20,15> doubles numBlocks in Z
# (6 -> 12) and therefore total ranks (144 -> 288). Per-rank work is unchanged (same cellsPerBlock),
# so the per-step rate and the ~1.8h estimate above should hold; 4h kept for the same margin as before.
# A 2-node/144-rank fallback (coarser Z resolution, same physical domain) exists as
# run_144_clog_measurement_validation_extreme_pidon_fresh_2node.sbatch if 4 nodes aren't available.
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

for dir in vtk-out eval snapshots; do
    if [[ -d "$BIN_DIR/$dir" ]]; then
        echo "WARNING: removing leftover $BIN_DIR/$dir from an unarchived previous run" >&2
        rm -rf "${BIN_DIR:?}/$dir"
    fi
done

"$SH" --build --jobs "$(nproc)"

echo "############################################################"
echo "## clog measurement validation: extreme clog, PID on (fresh, single-job, 4-node/288-rank, Z=180)"
echo "############################################################"
"$SH" --run --np 288 --prm "$PRM_DIR/clog_measurement_validation_extreme_pidon_gas.prm"

# [Claude, 2026-09-16]: prefer $HPCVAULT, decided and verified HERE (right before use), not at job
# start -- $HPCVAULT is autofs-mounted, and a prior run this session (jobs 4238846/4239167) showed a
# job-start mkdir can go stale by the time a multi-hour run's end-of-script mv actually happens
# ("mv: cannot move ... /home/vault/... No such file or directory" despite $HPCVAULT being valid and
# the earlier mkdir having reported success). Falling back to local-only risks the WORKING DIRECTORY's
# quota on an unattended overnight run, so local is a last resort, not a default.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_extreme_pidon"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
