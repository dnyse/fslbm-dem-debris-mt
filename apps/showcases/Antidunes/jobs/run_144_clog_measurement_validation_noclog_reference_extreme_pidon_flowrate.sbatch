#!/bin/bash -l
#SBATCH --job-name=antidunes-clog-measurement-validation-noclog-reference-extreme-pidon-flowrate
#SBATCH --partition=multinode
#SBATCH --nodes=2
#SBATCH --ntasks=144
#SBATCH --time=01:45:00
#SBATCH --output=sbatch-output/output/clog_measurement_validation_noclog_reference_extreme_pidon_flowrate_%j.out
#SBATCH --error=sbatch-output/error/clog_measurement_validation_noclog_reference_extreme_pidon_flowrate_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# Thesis item 3/4 (THESIS_PLAN.md) -- standalone run of the extreme no-clog reference (PID-on), needed
# to properly baseline the completed extreme_pidon_gas.prm run (job 4240982, COMPLETED, 115000 steps).
# Found directly by inspecting that run's own sqlite results (2026-09-17): the pre-release phase of the
# elevation-difference metric is dominated by a large flow-spinup transient (up to +-11 LU swing
# between the upstream/downstream probes) that has nothing to do with clogging, so a naive
# pre-release-vs-final comparison on the with-debris run alone cannot distinguish "clog signal" from
# "spin-up wave noise." This run reproduces that exact same transient with zero debris, giving the fair
# wave-controlled baseline the afflux/baseline-subtraction methodology
# (clogging/probe-and-forcing.tex's own "Baseline-Referenced and Complementary Metrics") already calls
# for. See clog_measurement_validation_noclog_reference_extreme_flowrate_gas.prm's own header for the full
# derivation (timesteps=115000 to match the with-debris run's own extended duration exactly).
#
# 144 ranks / 2 nodes, matching every other job in this "extreme" clog-measurement family. No debris at
# all, so no [RELEASE-OVERLAP] risk; 115000 steps at the with-debris pair's own measured rate
# (~1h15-1h20m) should be similar or somewhat faster (no floater/contact overhead) -- 3h requested for
# margin.
# [Claude, 2026-09-22] RERUN of the no-clog PID-ON reference (job 4242918, source of the 'Clog Flow-Forcing Comparison' artifact's no-clog baseline) with the
# ClogMeasurement flow-rate evaluator turned on (flowRateEnabled true, Evaluators.h flowRateAt -- didn't exist when 4242918 ran). Everything else identical: 144
# ranks, no debris, PID on, extreme-scenario domain/pillar geometry. Time limit 1h45m (4242918 itself COMPLETED in 1h14m52s -- tightened from the original 3h budget).
# Queue: afterok:4256467 (the currently RUNNING 5-run torque sweep), SAME checkout ~/mt-dennys-mirror2. See prms/clog_measurement_validation_noclog_reference_extreme_flowrate_gas.prm.
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
echo "## clog measurement validation: no-clog reference, extreme, PID on"
echo "############################################################"
"$SH" --run --np 144 --prm "$PRM_DIR/clog_measurement_validation_noclog_reference_extreme_flowrate_gas.prm"

# Prefer $HPCVAULT, decided and verified HERE (right before use), not at job start -- see
# run_144_clog_measurement_validation_extreme_pidoff.sbatch's own 2026-09-16 comment for the full
# rationale (autofs staleness, local-only is a last resort not a default).
STAMP="$(date +%Y%m%d_%H%M%S)"
VAULT_CANDIDATE="${HPCVAULT:-}/antidunes-archive/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${STAMP}"
if [[ -n "${HPCVAULT:-}" ]] && mkdir -p "$VAULT_CANDIDATE" 2>/dev/null \
   && touch "$VAULT_CANDIDATE/.write_test" 2>/dev/null; then
    rm -f "$VAULT_CANDIDATE/.write_test"
    ARCHIVE_ROOT="$VAULT_CANDIDATE"
else
    ARCHIVE_ROOT="$BIN_DIR/archive/${STAMP}_noclog_reference_extreme_pidon"
    mkdir -p "$ARCHIVE_ROOT"
    echo "WARNING: \$HPCVAULT unavailable at archive time -- falling back to LOCAL archive at $ARCHIVE_ROOT (move this off manually to avoid quota pressure)" >&2
fi
for d in vtk-out eval snapshots; do
    [[ -d "$BIN_DIR/$d" ]] && mv "$BIN_DIR/$d" "$ARCHIVE_ROOT/$d"
done

echo "## DONE. Archived under $ARCHIVE_ROOT"
