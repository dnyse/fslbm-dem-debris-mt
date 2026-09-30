#!/bin/bash
set -uo pipefail
BASE=/home/dhuber/University/Thesis/mt-dennys/apps/showcases/Antidunes
BIN=/home/dhuber/University/Thesis/mt-dennys/build-release/apps/showcases/Antidunes/Antidunes
for name in sphere_pair box_pair; do
  echo "=== [$(date '+%F %T')] starting $name ==="
  cd "$BASE/${name}_output"
  rm -rf vtk-out eval snapshots
  mpirun -np 9 "$BIN" "$BASE/prms/${name}_repro_gas.prm" > run.log 2>&1
  RC=$?
  echo "=== [$(date '+%F %T')] $name finished, exit=$RC ==="
  cd "$BASE"
done
echo "=== SHAPE COMPARISON DONE ==="
