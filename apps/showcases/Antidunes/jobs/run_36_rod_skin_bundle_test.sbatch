#!/bin/bash -l
#SBATCH --job-name=antidunes-rod-skin-bundle
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=01:30:00
#SBATCH --output=sbatch-output/output/rod_skin_bundle_%j.out
#SBATCH --error=sbatch-output/error/rod_skin_bundle_%j.err

# [Claude, 2026-09-30] Rod contact skin test (ModelParameters.cylinderContactSkin, SimulationParams.h), single node.
# Four independent 9-rank runs of the local three-rod bundle channel (gen_rod_bundle_channel_test_prm.py), started
# in parallel (4 x 9 = 36 of 72 cores), each in its own empty result directory:
#   skinA_gapm0p03      prescribed bundle motion, surface gap -0.03 LU (baseline: the sealed-channel geometry)
#   skinA_gap1p0        prescribed bundle motion, gap 1.0 LU (what a skin of 0.5 enforces)
#   skinA_gap2p0        prescribed bundle motion, gap 2.0 LU (what a skin of 1.0 enforces)
#   skinB_free_skin0p5  free rods (no prescribed motion), start touching, cylinderContactSkin 0.5, floater history on
# A answers "does a resolved fluid gap remove the density collapse?" (DEM-independent); B checks that the skin
# actually separates touching rods to the expected gap of 2m = 1.0 LU.
# Submit from Antidunes/ (project convention):  cd .../apps/showcases/Antidunes && sbatch jobs/run_36_rod_skin_bundle_test.sbatch
set -uo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO="$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
fi
AD="$REPO/apps/showcases/Antidunes"
SH="$AD/antidunes.sh"
BIN="$REPO/build/apps/showcases/Antidunes/Antidunes"
OUT="$AD/skin_results/${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT" || { echo "results dir $OUT already exists -- refusing to overwrite"; exit 1; }

"$SH" --build --jobs "$(nproc)" || { echo "build failed"; exit 1; }
module load gcc/11.5.0 openmpi/4.1.3-gcc11.5.0

CASES="skinA_gapm0p03 skinA_gap1p0 skinA_gap2p0 skinB_free_skin0p5"
for n in $CASES; do
    mkdir -p "$OUT/$n"
    cp "$AD/prms/${n}_gas.prm" "$OUT/$n/"
    ( cd "$OUT/$n" && mpirun -n 9 --bind-to none "$BIN" "${n}_gas.prm" > run.log 2>&1; echo "done $n exit=$?" >> run.log ) &
done
wait
for n in $CASES; do tail -1 "$OUT/$n/run.log"; done
echo "## DONE. Results in $OUT"
