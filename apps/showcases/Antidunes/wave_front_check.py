#!/usr/bin/env python3
"""Track free-surface wave-front position and liquid-volume conservation over an Antidunes-format
fluid_field VTK run (dam-break/tsunami scenarios in particular, but works for any run using the
shared free_surface::addVTKOutput writer with fill_level + mapped_flag enabled).

Reuses find_degenerate_cells.py's dependency-free VTU parser (plain Python 3 + numpy, no
pyvista/vtk/paraview needed -- runs on a login node) rather than re-implementing it; only adds a
loader that also pulls fill_level (find_degenerate_cells.load_step only loads mapped_flag +
velocity, since that script doesn't need fill level).

Checks performed, per requested step:
  1. wave-front position: max x-coordinate among liquid(1)/interface(2) cells with
     fill_level > --fill-threshold. Advancing monotonically and settling near a domain boundary
     is the expected dam-break/tsunami signature; a stalled, teleporting, or oscillating front
     signals a setup or numerical problem.
  2. liquid volume conservation: sum of fill_level over liquid+interface cells, reported as a
     percent drift relative to the first requested step. FSLBM's excess-mass redistribution
     should keep this within a fraction of a percent over the whole run; large drift indicates a
     broken mass-conservation path.

Usage:
  python3 wave_front_check.py <run-dir-or-vtk-out-or-fluid_field-dir> [--steps N [N ...]] \\
      [--fill-threshold F]

Example:
  python3 wave_front_check.py \\
      ~/University/Thesis/animations/tsunami-pole-forest-driftwood_20260906 \\
      --steps 0 2500 5000 7500 9500 12000 14750
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from find_degenerate_cells import parse_vtu, find_step_dirs, LIQUID, INTERFACE  # noqa: E402


def load_step_with_fill(step_dir):
    """Like find_degenerate_cells.load_step, but also stitches together fill_level."""
    files = sorted(glob.glob(os.path.join(step_dir, '*.vtu')))
    if not files:
        raise FileNotFoundError(f'no .vtu block files in {step_dir}')
    centers_list, flags_list, fill_list = [], [], []
    for f in files:
        arrs = parse_vtu(f)
        points = arrs['_unnamed_0']
        conn = arrs['connectivity'].reshape(-1, 8)
        centers_list.append(points[conn].mean(axis=1))
        flags_list.append(arrs['mapped_flag'])
        fill_list.append(arrs['fill_level'])
    centers = np.concatenate(centers_list, axis=0)
    flags = np.concatenate(flags_list, axis=0)
    fill = np.concatenate(fill_list, axis=0)
    return centers, flags, fill


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path')
    ap.add_argument('--steps', type=int, nargs='+', default=None)
    ap.add_argument('--fill-threshold', type=float, default=0.01)
    args = ap.parse_args()

    all_steps = find_step_dirs(args.path)
    if args.steps is not None:
        wanted = set(args.steps)
        all_steps = [(s, d) for s, d in all_steps if s in wanted]
        if not all_steps:
            print('no matching steps found', file=sys.stderr)
            sys.exit(1)

    print(f"{'step':>8} {'ncells':>10} {'front_x':>10} {'liq_vol':>12} {'drift_%':>10}")
    base_liq_vol = None
    for step, step_dir in all_steps:
        centers, flags, fill = load_step_with_fill(step_dir)
        liquid_or_interface = np.isin(flags, [LIQUID, INTERFACE])
        liq_fill = np.where(liquid_or_interface, fill, 0.0)

        liq_vol = float(liq_fill.sum())
        if base_liq_vol is None:
            base_liq_vol = liq_vol
        drift = 100.0 * (liq_vol - base_liq_vol) / base_liq_vol if base_liq_vol else 0.0

        front_mask = liquid_or_interface & (fill > args.fill_threshold)
        front_x = float(centers[front_mask, 0].max()) if front_mask.any() else float('nan')

        print(f"{step:>8} {len(flags):>10} {front_x:>10.2f} {liq_vol:>12.2f} {drift:>+10.4f}")


if __name__ == '__main__':
    main()
