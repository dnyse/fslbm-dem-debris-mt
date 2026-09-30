#!/usr/bin/env python3
"""Scan Antidunes fluid_field VTK output for degenerate free-surface cell configurations.

Parses the per-block "*.vtu" files written by walberla::free_surface::addVTKOutput directly
(base64-encoded legacy VTK XML) -- no pyvista/vtk dependency, so it runs in a plain Python 3 +
numpy environment (e.g. on a login node without the simulation's build environment).

Requires the "mapped_flag" writer enabled in the run's .prm (VTK.fluid_field.writers) -- see
VtkWriter.h for the flag->integer mapping reproduced in FLAG_NAMES below. "velocity" is optional;
if absent, the velocity-based checks are skipped for that run.

Checks performed, all against the FSLBM invariant that liquid must never be a stencil neighbour
of gas (see PdfReconstructionManager.h's bognerTargetFlag for the production-code statement of
this invariant):
  1. liquid<->gas direct face adjacency (no interface cell in between)          -- hard violation
  2. orphaned interface cells: no liquid AND no gas among the 6 face neighbours -- can never
     resolve into a normal liquid/gas transition without help; a precursor pattern, not itself
     a violation of (1)
  3. unmapped/anomalous mapped_flag values (not in the known set) -- catches bitwise-OR flag
     collisions (field/vtk/FlagFieldMapping.h ORs every matching mapped flag together)
  4. non-finite velocity (nan/inf) at a liquid/interface/gas cell -- already-corrupted PDFs
  5. velocity magnitude above --vel-threshold -- precursor to the kind of blow-up documented in
     changes.md (an acute force/velocity spike overflowing the LBM's stability limit a few steps
     later)

Usage:
  python3 find_degenerate_cells.py <run-dir-or-vtk-out-or-fluid_field-dir> [options]

  <path> may be:
    - a top-level run directory containing vtk-out/fluid_field/simulation_step_*
    - the vtk-out/fluid_field directory itself
    - a single simulation_step_NNNN directory (checks just that one step)

Options:
  --steps N [N ...]   only check these step numbers (default: all available)
  --max-examples N    world-coordinate examples to print per category (default 15)
  --vel-threshold F   velocity-magnitude |u| flag threshold, lattice units (default 0.17 = lattice
                       Mach number 0.3 with cs = 1/sqrt(3); Kruger et al. 2017 Sect. 7.2 p. 279: Ma < 0.3 is
                       'small', larger Ma gives significant compressibility errors. An accuracy warning, not a
                       collapse threshold. Before 2026-09-21 the default was |u| = 0.3, which is Ma = 0.52.)
  --non-periodic-x    domain is NOT periodic in X (default: periodic, matches this showcase)
  --non-periodic-y    domain is NOT periodic in Y (default: periodic, matches this showcase)
  --periodic-z        domain IS periodic in Z (default: not periodic)
  --quiet             only print the summary timeline table, suppress per-step example dumps

Example:
  python3 find_degenerate_cells.py \\
      ~/University/Thesis/animations/antidunes-cylinder-sedimentation-showcase_20260723_153223
"""
import argparse
import base64
import glob
import os
import re
import sys

import numpy as np

DTYPE_MAP = {
    'Float32': np.float32, 'Float64': np.float64,
    'Int32': np.int32, 'UInt32': np.uint32,
    'Int64': np.int64, 'UInt64': np.uint64,
    'UInt8': np.uint8, 'Int8': np.int8,
}

# VtkWriter.h's flagMapper->addMapping(...) calls -- the only values a well-formed mapped_flag
# cell should ever take (0 = no mapped flag set at all, e.g. a moving-obstacle/particle-interior
# cell, which is expected and not itself an anomaly).
KNOWN_FLAGS = {
    0: 'unmapped (obstacle/particle interior)',
    1: 'liquid', 2: 'interface', 3: 'gas',
    4: 'noSlip', 6: 'freeSlip/UBB', 7: 'UBB-inflow',
    8: 'pressure', 9: 'pressure-outflow', 10: 'outlet',
}
LIQUID, INTERFACE, GAS = 1, 2, 3


def parse_vtu(fname):
    """Return {DataArray name -> numpy array} for one legacy-XML base64 .vtu file."""
    with open(fname, 'rb') as fh:
        text = fh.read().decode('latin1')
    arrs = {}
    unnamed = 0
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', text, re.DOTALL):
        attrs = m.group(1)
        namem = re.search(r'Name="([^"]*)"', attrs)
        if namem:
            name = namem.group(1)
        else:
            name = f'_unnamed_{unnamed}'
            unnamed += 1
        vtype_m = re.search(r'type="([^"]*)"', attrs)
        assert vtype_m is not None
        vtype = vtype_m.group(1)
        ncompm = re.search(r'NumberOfComponents="([^"]*)"', attrs)
        ncomp = int(ncompm.group(1)) if ncompm else 1
        raw = base64.b64decode(re.sub(r'\s+', '', m.group(2)))
        nbytes = int(np.frombuffer(raw[:4], dtype=np.uint32)[0])
        arr = np.frombuffer(raw[4:4 + nbytes], dtype=DTYPE_MAP[vtype])
        if ncomp > 1:
            arr = arr.reshape(-1, ncomp)
        arrs[name] = arr
    return arrs


def load_step(step_dir):
    """Stitch every block file in a simulation_step_NNNN dir into one global (centers, flags,
    velocity|None) triple, in world coordinates (each block's Points are already world coords,
    so no manual block-offset bookkeeping is needed -- robust to the Hilbert-curve rank<->block
    assignment used by blockforest::createUniformBlockGrid)."""
    files = sorted(glob.glob(os.path.join(step_dir, '*.vtu')))
    if not files:
        raise FileNotFoundError(f'no .vtu block files in {step_dir}')
    centers_list, flags_list, vel_list = [], [], []
    have_vel = True
    for f in files:
        arrs = parse_vtu(f)
        points = arrs['_unnamed_0']
        conn = arrs['connectivity'].reshape(-1, 8)
        centers_list.append(points[conn].mean(axis=1))
        flags_list.append(arrs['mapped_flag'])
        if 'velocity' in arrs:
            vel_list.append(arrs['velocity'])
        else:
            have_vel = False
    centers = np.concatenate(centers_list, axis=0)
    flags = np.concatenate(flags_list, axis=0)
    vel = np.concatenate(vel_list, axis=0) if have_vel else None
    return centers, flags, vel


def to_grid(centers, flags, vel):
    """Bin (centers, flags, vel) into dense 3D arrays indexed by integer cell coordinate.
    Assumes unit cell spacing and cell-centered sampling at (i+0.5, j+0.5, k+0.5), which is what
    walberla's VTK writer always produces for a uniform block-structured field."""
    ix = np.round(centers[:, 0] - 0.5).astype(np.int64)
    iy = np.round(centers[:, 1] - 0.5).astype(np.int64)
    iz = np.round(centers[:, 2] - 0.5).astype(np.int64)
    nx, ny, nz = int(ix.max()) + 1, int(iy.max()) + 1, int(iz.max()) + 1
    if nx * ny * nz != len(flags):
        print(f'  WARNING: grid dims {nx}x{ny}x{nz}={nx*ny*nz} != {len(flags)} cells loaded '
              f'-- domain may not be a single uniform level, or cells are missing/duplicated.',
              file=sys.stderr)
    grid = np.zeros((nx, ny, nz), dtype=np.uint32)
    grid[ix, iy, iz] = flags
    velgrid = None
    if vel is not None:
        velgrid = np.zeros((nx, ny, nz, 3), dtype=np.float64)
        velgrid[ix, iy, iz] = vel
    return grid, velgrid


def check_axis_adjacency(grid, axis, periodic):
    g0 = grid
    g1 = np.roll(grid, -1, axis=axis)
    bad = ((g0 == LIQUID) & (g1 == GAS)) | ((g0 == GAS) & (g1 == LIQUID))
    if not periodic:
        sl: list = [slice(None)] * 3
        sl[axis] = -1
        bad[tuple(sl)] = False
    return bad


def analyze_step(grid, velgrid, periodicity, vel_threshold):
    px, py, pz = periodicity
    results = {}

    # 1. direct liquid<->gas adjacency
    bad_x = check_axis_adjacency(grid, 0, px)
    bad_y = check_axis_adjacency(grid, 1, py)
    bad_z = check_axis_adjacency(grid, 2, pz)
    bad_cells = np.zeros(grid.shape, dtype=bool)
    for bad, axis in ((bad_x, 0), (bad_y, 1), (bad_z, 2)):
        bad_cells |= bad
        bad_cells |= np.roll(bad, 1, axis=axis)
    results['liquid_gas_adjacency'] = (int(bad_x.sum() + bad_y.sum() + bad_z.sum()),
                                       np.argwhere(bad_cells))

    # 2. orphaned interface cells (no liquid AND no gas among 6 face neighbours)
    is_intf = grid == INTERFACE
    has_liquid_nbr = np.zeros(grid.shape, dtype=bool)
    has_gas_nbr = np.zeros(grid.shape, dtype=bool)
    for axis, periodic in ((0, px), (1, py), (2, pz)):
        for shift in (1, -1):
            rolled = np.roll(grid, shift, axis=axis)
            if not periodic:
                sl: list = [slice(None)] * 3
                sl[axis] = 0 if shift == 1 else -1
                rolled[tuple(sl)] = 255  # sentinel: no such neighbour, matches nothing
            has_liquid_nbr |= (rolled == LIQUID)
            has_gas_nbr |= (rolled == GAS)
    orphaned = is_intf & ~has_liquid_nbr & ~has_gas_nbr
    results['orphaned_interface'] = (int(orphaned.sum()), np.argwhere(orphaned))

    # 3. anomalous / unmapped flag values
    known_mask = np.isin(grid, list(KNOWN_FLAGS.keys()))
    anomalous = ~known_mask
    results['anomalous_flag'] = (int(anomalous.sum()), np.argwhere(anomalous))

    # 4/5. velocity checks (only for liquid/interface/gas cells)
    if velgrid is not None:
        fluid_mask = (grid == LIQUID) | (grid == INTERFACE) | (grid == GAS)
        vmag = np.linalg.norm(velgrid, axis=-1)
        nonfinite = fluid_mask & ~np.isfinite(vmag)
        results['nonfinite_velocity'] = (int(nonfinite.sum()), np.argwhere(nonfinite))
        finite_fluid = fluid_mask & np.isfinite(vmag)
        toofast = finite_fluid & (vmag > vel_threshold)
        results['high_velocity'] = (int(toofast.sum()), np.argwhere(toofast))
        results['_vmag_max'] = float(np.nanmax(np.where(finite_fluid, vmag, np.nan))) \
            if finite_fluid.any() else 0.0
    else:
        results['nonfinite_velocity'] = (0, np.empty((0, 3), dtype=int))
        results['high_velocity'] = (0, np.empty((0, 3), dtype=int))
        results['_vmag_max'] = None

    return results


def print_examples(label, idxs, grid, velgrid, max_examples):
    print(f'    {label}: {len(idxs)} cell(s)' + (' (showing up to '
          f'{max_examples})' if len(idxs) > max_examples else ''))
    for x, y, z in idxs[:max_examples]:
        flag = grid[x, y, z]
        flagname = KNOWN_FLAGS.get(int(flag), f'UNKNOWN({flag})')
        extra = ''
        if velgrid is not None:
            v = velgrid[x, y, z]
            vmag = np.linalg.norm(v)
            extra = f'  vel=({v[0]:.4g},{v[1]:.4g},{v[2]:.4g}) |v|={vmag:.4g}'
        print(f'      world≈({x + 0.5},{y + 0.5},{z + 0.5})  flag={flag}({flagname}){extra}')


def find_step_dirs(path):
    """Resolve `path` (run dir / vtk-out / fluid_field / single step dir) to a sorted list of
    (step_number, step_dir) pairs."""
    path = os.path.abspath(os.path.expanduser(path))
    base = os.path.basename(path.rstrip('/'))
    if base.startswith('simulation_step_'):
        step = int(base[len('simulation_step_'):])
        return [(step, path)]
    candidates = [
        path,
        os.path.join(path, 'fluid_field'),
        os.path.join(path, 'vtk-out', 'fluid_field'),
    ]
    for c in candidates:
        matches = sorted(glob.glob(os.path.join(c, 'simulation_step_*')))
        matches = [m for m in matches if os.path.isdir(m)]
        if matches:
            out = []
            for m in matches:
                b = os.path.basename(m)
                out.append((int(b[len('simulation_step_'):]), m))
            out.sort(key=lambda t: t[0])
            return out
    raise FileNotFoundError(
        f'could not find simulation_step_* directories under any of: {candidates}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path')
    ap.add_argument('--steps', type=int, nargs='+', default=None)
    ap.add_argument('--max-examples', type=int, default=15)
    ap.add_argument('--vel-threshold', type=float, default=0.17)
    ap.add_argument('--non-periodic-x', action='store_true')
    ap.add_argument('--non-periodic-y', action='store_true')
    ap.add_argument('--periodic-z', action='store_true')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()

    periodicity = (not args.non_periodic_x, not args.non_periodic_y, args.periodic_z)

    all_steps = find_step_dirs(args.path)
    if args.steps is not None:
        wanted = set(args.steps)
        all_steps = [(s, d) for s, d in all_steps if s in wanted]
        if not all_steps:
            print('no matching steps found', file=sys.stderr)
            sys.exit(1)

    timeline = []
    for step, step_dir in all_steps:
        print(f'=== step {step} ({step_dir}) ===')
        try:
            centers, flags, vel = load_step(step_dir)
        except Exception as e:
            print(f'  FAILED to load: {e}', file=sys.stderr)
            continue
        grid, velgrid = to_grid(centers, flags, vel)
        results = analyze_step(grid, velgrid, periodicity, args.vel_threshold)

        n_adj, idx_adj = results['liquid_gas_adjacency']
        n_orph, idx_orph = results['orphaned_interface']
        n_anom, idx_anom = results['anomalous_flag']
        n_nfv, idx_nfv = results['nonfinite_velocity']
        n_hv, idx_hv = results['high_velocity']
        vmax = results['_vmag_max']

        print(f'  liquid<->gas direct adjacency (no interface between): {n_adj} link(s)')
        print(f'  orphaned interface cells (no liquid/gas face-neighbour): {n_orph}')
        print(f'  anomalous/unmapped mapped_flag values: {n_anom}')
        if velgrid is not None:
            print(f'  non-finite (nan/inf) velocity at fluid cells: {n_nfv}')
            print(f'  |velocity| > {args.vel_threshold}: {n_hv}  (max |v| this step: {vmax:.4g})')
        else:
            print('  (no "velocity" writer found -- velocity checks skipped)')

        if not args.quiet:
            if n_adj:
                print_examples('liquid<->gas adjacency', idx_adj, grid, velgrid,
                                args.max_examples)
            if n_orph:
                print_examples('orphaned interface', idx_orph, grid, velgrid, args.max_examples)
            if n_anom:
                print_examples('anomalous flag', idx_anom, grid, velgrid, args.max_examples)
            if n_nfv:
                print_examples('non-finite velocity', idx_nfv, grid, velgrid, args.max_examples)
            if n_hv:
                print_examples('high velocity', idx_hv, grid, velgrid, args.max_examples)

        timeline.append((step, n_adj, n_orph, n_anom, n_nfv, n_hv, vmax))
        print()

    if len(timeline) > 1:
        print('=== timeline summary ===')
        header = f'{"step":>8} {"adj":>6} {"orph":>6} {"anom":>6} {"nfv":>6} {"hivel":>6} {"maxV":>10}'
        print(header)
        for step, n_adj, n_orph, n_anom, n_nfv, n_hv, vmax in timeline:
            vmax_s = f'{vmax:.4g}' if vmax is not None else '-'
            print(f'{step:>8} {n_adj:>6} {n_orph:>6} {n_anom:>6} {n_nfv:>6} {n_hv:>6} {vmax_s:>10}')


if __name__ == '__main__':
    main()
