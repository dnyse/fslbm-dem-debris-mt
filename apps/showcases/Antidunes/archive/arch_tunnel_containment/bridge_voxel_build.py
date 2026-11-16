#!/usr/bin/env python3
"""Stage 1/2: compute the bridge's NoSlip occupancy grid at TRUE simulation resolution
(1 LU/cell, matching the real prm exactly) and save it to a file -- separated from
rendering because the full-bridge region is ~140M cells, too large for matplotlib's
ax.voxels() to render directly (see bridge_voxel_render_full.py for the actual render,
which only needs the much smaller SURFACE of this grid).

Classifies each cell with the same boolean composition as the real BodyBC blocks (see
bridge_voxel_view.py's docstring for the full pier-box/pier-cylinder/arch breakdown), but
stores a single uint8 category code per cell instead of an RGBA color string (much
cheaper to hold in memory / write to disk for a grid this size):
    0 = empty, 1 = pier box (Pillar), 2 = pier end-cap (CylinderPillar), 3 = arch collar

Processed in Y-slices so peak memory stays bounded regardless of total grid size, with a
progress print after each slice batch (this is the "which stage is it in" the run asked
for -- elapsed/estimated-remaining time, not just a silent multi-minute wait).

Usage:
    python3 bridge_voxel_build.py
Writes bridge_voxels_full.npz (category grid + region/voxel_size metadata).
"""
import time
import numpy as np

from bridge_layout_optimizer import BRIDGE, BridgeLayout


def build_occupancy(layout: BridgeLayout, deck_top: float, region, voxel_size: float,
                     y_chunk: int = 32, out_path: str = "bridge_voxels_full.npz") -> None:
    x0, x1, y0, y1, z0, z1 = region
    xs = np.arange(x0, x1, voxel_size)
    ys = np.arange(y0, y1, voxel_size)
    zs = np.arange(z0, z1, voxel_size)
    nx, ny, nz = len(xs), len(ys), len(zs)
    total_cells = nx * ny * nz
    print(f"[build] grid {nx} x {ny} x {nz} = {total_cells:,} cells at {voxel_size:g} LU/cell "
          f"(region X=[{x0:g},{x1:g}] Y=[{y0:g},{y1:g}] Z=[{z0:g},{z1:g}])")
    est_bytes = total_cells  # uint8, one byte/cell
    print(f"[build] estimated in-memory size: {est_bytes/1e9:.2f} GB (uint8 category grid)")

    category = np.zeros((nx, ny, nz), dtype=np.uint8)
    cx = xs + voxel_size / 2
    cz = zs + voxel_size / 2
    X2d, Z2d = np.meshgrid(cx, cz, indexing="ij")  # reused every Y-slice-batch, shape (nx, nz)

    n_batches = (ny + y_chunk - 1) // y_chunk
    t_start = time.time()
    for b in range(n_batches):
        j0 = b * y_chunk
        j1 = min(j0 + y_chunk, ny)
        cy = ys[j0:j1] + voxel_size / 2  # shape (bj,)
        # Broadcast to (nx, bj, nz): X2d[:,None,:], cy[None,:,None], Z2d[:,None,:]
        X = X2d[:, None, :]
        Y = cy[None, :, None]
        Z = Z2d[:, None, :]

        cat = np.zeros((nx, j1 - j0, nz), dtype=np.uint8)
        for p in layout.piers:
            in_box = ((p.x_center - p.x_half <= X) & (X <= p.x_center + p.x_half) &
                      (p.y_min <= Y) & (Y <= p.y_max) & (0.0 <= Z) & (Z <= p.z_top))
            cat[in_box] = 1
            in_cyl_a = ((X - p.x_cyl_a) ** 2 + (Y - p.y_center) ** 2 <= p.radius ** 2) & (0.0 <= Z) & (Z <= p.z_top)
            in_cyl_b = ((X - p.x_cyl_b) ** 2 + (Y - p.y_center) ** 2 <= p.radius ** 2) & (0.0 <= Z) & (Z <= p.z_top)
            cat[in_cyl_a | in_cyl_b] = 2
        for a in layout.arches:
            in_box = ((a.x_min <= X) & (X <= a.x_max) & (a.y_min <= Y) & (Y <= a.y_max) &
                      (a.z_center <= Z) & (Z <= deck_top))
            outside_bore = (Y - a.y_center) ** 2 + (Z - a.z_center) ** 2 > a.radius ** 2
            arch_mask = in_box & outside_bore
            cat[arch_mask & (cat == 0)] = 3  # don't overwrite a pier cell already set this slice

        category[:, j0:j1, :] = cat

        elapsed = time.time() - t_start
        frac = (b + 1) / n_batches
        eta = elapsed / frac - elapsed if frac > 0 else 0.0
        print(f"[build] Y-slice batch {b+1}/{n_batches} (Y={ys[j0]:.0f}..{ys[j1-1]+voxel_size:.0f}) "
              f"-- {elapsed:5.1f}s elapsed, ~{eta:5.1f}s remaining", flush=True)

    n_solid = int(np.count_nonzero(category))
    print(f"[build] done classifying: {n_solid:,} solid cells ({n_solid/total_cells*100:.2f}%), "
          f"{time.time()-t_start:.1f}s total")

    print(f"[build] writing {out_path} (compressed) ...", flush=True)
    t_save = time.time()
    np.savez_compressed(out_path, category=category,
                         x0=x0, x1=x1, y0=y0, y1=y1, z0=z0, z1=z1,
                         voxel_size=voxel_size, deck_top=deck_top,
                         waterline_z=layout.waterline_z)
    print(f"[build] saved in {time.time()-t_save:.1f}s")

    import os
    size_mb = os.path.getsize(out_path) / 1e6
    print(f"[build] {out_path}: {size_mb:.1f} MB on disk")


if __name__ == "__main__":
    VOXEL_SIZE = 1.0  # true simulation resolution
    deck_top = max((p.z_top for p in BRIDGE.piers), default=BRIDGE.domain[2] * 0.4)
    arch_x0 = min((a.x_min for a in BRIDGE.arches), default=0.0)
    arch_x1 = max((a.x_max for a in BRIDGE.arches), default=BRIDGE.domain[0])
    margin = 40.0
    REGION = (arch_x0 - margin, arch_x1 + margin, 0.0, BRIDGE.domain[1], 0.0, deck_top + margin)

    build_occupancy(BRIDGE, deck_top, REGION, VOXEL_SIZE)
