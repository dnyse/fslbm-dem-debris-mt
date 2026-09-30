#!/bin/bash -l
#SBATCH --job-name=antidunes-cuboid-stillwater-realgravity981-d050
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=03:00:00
#SBATCH --output=sbatch-output/output/cuboid_stillwater_realgravity981_d050_%j.out
#SBATCH --error=sbatch-output/error/cuboid_stillwater_realgravity981_d050_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# [Claude, 2026-09-18] Thesis item 1a (THESIS_PLAN.md) -- rerun of the ORIGINAL, lower-resolution
# d050 equilibrium scenario (cuboid_stillwater_realgravity981_3xsphere_gas.prm, 72 ranks, domain
# <96,135,120>, 50000 timesteps), per explicit request. This exact config already completed cleanly
# once before (job 4038343, COMPLETED, exit 0, 2h26m48s) -- this rerun is specifically to capture the
# tiltY log field just added to ParticleKernels.h's [BoxForces] line (job 4038343 predates that
# change and has no rotation-angle history at all: only 3 unrelated one-off mentions of
# [BoxBuoyancy]/[BoxRightingMoment] in its whole log, not a per-step series -- that evaluator does
# not fire for this free-floating scenario type).
#
# Also directly informative for the still-open hires-resolution crash investigation
# (cuboid_stillwater_realgravity981_3xsphere_hires_gas.prm, jobs 4240983/4244969, both crashed):
# every ModelParameters flag between this file and the hires one is byte-for-byte identical -- the
# only differences are resolution/rank-decomposition (72 ranks/6x3x4 blocks vs 140 ranks/7x4x5
# blocks) and the derived timesteps/gravityOverride/box-size values that scale with it. A second
# clean completion here (expected, given job 4038343's own precedent) would reinforce that the
# crash is resolution/decomposition-specific, not a general property of this physical setup.
#
# Requesting 3h (vs. job 4038343's own measured 2h27m) for margin. singlenode partition, matching
# this scenario family's own established pattern (run_72_cuboid_stillwater_d0{25,50,75}
# _reproducibility.sbatch) -- but WITHOUT that family's own --constraint=el9. Isolated precisely via
# sbatch --test-only: --partition=singlenode alone queues fine (confirmed, would start on node
# f0330); adding --constraint=el9 specifically causes "Invalid feature specification" outright
# rejection (not a queue/pending state), even though sinfo still lists el9 as an AVAIL_FEATURE for
# this partition. Cluster-side issue with that specific constraint right now, not with singlenode
# itself or with this file -- dropping just the constraint keeps the requested partition intact.
#
# [Claude, 2026-09-18]: deliberately does NOT run `--bed --bed-np 72` before `--build`, unlike this
# scenario family's older reference scripts. Checked directly rather than assuming: the existing
# spheres_out.dat already on this checkout has header "0.0464248483188 0.0218704558877
# 0.145077650996", an exact match to BedGeneration.prm's own domainSize_SI and to what both this
# prm and its hires sibling expect. Since useParticleBed=false here, the actual settled particle
# positions in that file are never read at all -- only this header value (a fixed property of
# BedGeneration.prm's own domainSize_SI, not of bed-np or the settling dynamics) matters, and it
# already matches. Regenerating would be redundant, not incorrect, but isn't needed.
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
echo "## cuboid stillwater realgravity981, densityRatio=0.5, ORIGINAL 72-rank/50000-step config"
echo "############################################################"
"$SH" --run --np 72 --prm "$PRM_DIR/cuboid_stillwater_realgravity981_3xsphere_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use), not at job start -- see
# run_140_cuboid_stillwater_hires_d050.sbatch's own 2026-09-16 comment for the full rationale.
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_d050_realgravity981_72rank"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi

for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mkdir -p "$ARCHIVE_ROOT/d050_72rank" && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/d050_72rank/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT/d050_72rank/{eval,vtk-out,snapshots}"
