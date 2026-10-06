#!/usr/bin/env python3
"""Relabels a fluid_field_snapshots VTK output directory so its folder/file names and .pvd
timestep attributes match the REAL simulation step, instead of a sequential write-count.

WHY THIS EXISTS: the "fluid_field_snapshots" mechanism (VTKParams::snapshotsEnabled, gated in
src/lbm/free_surface/VtkWriter.h's addVTKOutput) only invokes the underlying VTKOutput's write()
call on the exact listed steps -- every other simulation tick, write() is never called at all. Its
internal executionCounter_ (src/vtk/VTKOutput.cpp) increments once per CALL, not once per
simulation tick, so after N snapshots it reads 0..N-1 regardless of which real steps they were --
giving folders/files named simulation_step_0, simulation_step_1, ... and a .pvd with
timestep="0","1",.... The particle/floater point-data VTK writers (createVTKOutput_PointData,
wired via BuildTimeloop.h's addSnapshotAwareVTKWrite) don't have this problem: they're called
every tick (just conditionally FORCED to write via VTKOutput::forceWriteNextStep() on the matching
steps), so their own executionCounter_ naturally tracks the real step count throughout.

A real fix needs either widening vtk::initializeVTKOutput's core, project-wide OutputFunction API
(std::function<void()>, no handle to the underlying VTKOutput) to expose the raw object, or
hand-rolling a parallel VTKOutput construction bypassing that generic factory just for this one
named block (risking silent drift from its filter/writer setup) -- neither is a good trade for a
purely cosmetic ParaView-timeline issue. This script fixes the SYMPTOM instead, safely and
losslessly: it only renames directories/files and rewrites path references/timestep attributes; it
never touches the actual field data (the base64 velocity/mapped_flag arrays inside each .vtu are
byte-for-byte untouched).

Usage:
    python3 fix_fluid_snapshot_naming.py <vtk-out-dir> <step1> <step2> ... <stepN>

The steps must be given in the SAME order as the prm's own "steps" list for this run (ascending,
matching VTKParams::snapshotSteps, which is itself a std::set<uint_t> and therefore always sorted
ascending regardless of the order written in the prm) -- i.e. just pass them sorted ascending.
If the run crashed partway through and only the first K of the full list actually fired, pass only
those first K real steps (this script asserts the snapshot count it finds on disk matches the
count of steps you pass, so a mismatch is caught immediately rather than silently mislabeling).

Example (the full 8-step production schedule):
    python3 fix_fluid_snapshot_naming.py \\
        /path/to/vtk-out 500 1000 2000 4000 8000 16000 32000 49000

Example (a run that crashed after only the first 6 steps fired):
    python3 fix_fluid_snapshot_naming.py /path/to/vtk-out 500 1000 2000 4000 8000 16000
"""
import os
import re
import sys


def fix(vtk_out_dir: str, steps: list[int]) -> None:
    base_name = "fluid_field_snapshots"
    base = os.path.join(vtk_out_dir, base_name)
    pvd_path = base + ".pvd"
    if not os.path.isdir(base) or not os.path.isfile(pvd_path):
        raise SystemExit(f"'{base}' / '{pvd_path}' not found under {vtk_out_dir} -- wrong path?")

    n = len(steps)
    for i in range(n):
        d = os.path.join(base, f"simulation_step_{i}")
        f = os.path.join(base, f"simulation_step_{i}.pvtu")
        if not os.path.isdir(d) or not os.path.isfile(f):
            raise SystemExit(
                f"expected {n} snapshots (indices 0..{n-1}) matching the {n} steps you passed, "
                f"but '{d}' or '{f}' is missing -- recount how many steps actually fired (check "
                f"the .pvd's own <DataSet> count) and pass exactly that many real step numbers."
            )

    # Reverse order: avoids any risk of a rename colliding with a not-yet-processed index (even
    # though real step numbers never collide with small sequential indices in practice).
    for i in reversed(range(n)):
        real = steps[i]
        old_dir = os.path.join(base, f"simulation_step_{i}")
        new_dir = os.path.join(base, f"simulation_step_{real}")
        old_pvtu = os.path.join(base, f"simulation_step_{i}.pvtu")
        new_pvtu = os.path.join(base, f"simulation_step_{real}.pvtu")

        with open(old_pvtu) as fh:
            content = fh.read()
        patched, nsub = re.subn(rf'Source="simulation_step_{i}/', f'Source="simulation_step_{real}/', content)
        if nsub == 0:
            raise SystemExit(f"no Source= references found to patch in {old_pvtu} -- unexpected file format")
        with open(old_pvtu, "w") as fh:
            fh.write(patched)

        os.rename(old_dir, new_dir)
        os.rename(old_pvtu, new_pvtu)
        print(f"index {i} -> step {real}: renamed dir+pvtu, patched {nsub} Source= refs")

    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        ' <Collection>',
    ]
    for real in steps:
        lines.append(f'  <DataSet timestep="{real}" file="{base_name}/simulation_step_{real}.pvtu"/>')
    lines += [' </Collection>', '</VTKFile>', '']
    with open(pvd_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"rewrote {pvd_path} with real timesteps {steps}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        raise SystemExit(1)
    fix(sys.argv[1], [int(s) for s in sys.argv[2:]])
