#!/usr/bin/env python3
"""Stage 2/2: render the full-resolution occupancy grid saved by bridge_voxel_build.py.

Does NOT hand the dense ~140M-cell grid to ax.voxels() (infeasible at that size) -- instead
extracts only the SURFACE faces (a cell face is visible only where a solid cell borders an
empty one, or the region edge). Two things matter for keeping this actually interactive at
2M+ solid cells:

  1. Visibility is computed ONCE against the union of all solid categories, not per-category
     separately -- otherwise a pier-box cell touching a pier-cylinder cell (different
     category, same material) gets a spurious "internal" face drawn between them, which is
     both visually wrong and inflates the face count for no reason.
  2. Each flat run of same-category, same-depth faces is merged into ONE big quad instead of
     one 1x1 LU quad per cell (run-length merging along one in-plane axis per row) -- a
     pier's flat side wall is physically two rectangles, not tens of thousands of unit
     squares, and matplotlib's Poly3DCollection genuinely cannot redraw millions of polygons
     interactively.

Usage:
    python3 bridge_voxel_render_full.py
Reads bridge_voxels_full.npz (run bridge_voxel_build.py first if it doesn't exist).
"""
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

CATEGORY_COLORS = {
    1: "#8a7f76",  # pier box (Pillar)
    2: "#c1552b",  # pier end-caps (CylinderPillar)
    3: "#5b6b82",  # arch collar
}
CATEGORY_NAMES = {1: "pier box", 2: "pier end-caps", 3: "arch collar"}


def visible_faces(occ: np.ndarray):
    """occ: boolean array, True = solid (any category). Returns 6 boolean masks (same shape
    as occ), one per face direction, True where that cell's face on that side is visible
    (borders empty space or the grid edge)."""
    padded = np.pad(occ, 1, mode="constant", constant_values=False)
    nx, ny, nz = occ.shape
    return {
        "-x": occ & ~padded[0:nx, 1:-1, 1:-1],
        "+x": occ & ~padded[2:nx + 2, 1:-1, 1:-1],
        "-y": occ & ~padded[1:-1, 0:ny, 1:-1],
        "+y": occ & ~padded[1:-1, 2:ny + 2, 1:-1],
        "-z": occ & ~padded[1:-1, 1:-1, 0:nz],
        "+z": occ & ~padded[1:-1, 1:-1, 2:nz + 2],
    }


def _runs_1d_by_category(row_cat: np.ndarray, max_run: int = 24):
    """row_cat: 1D uint8 array, 0 = not visible/empty, else the category code. Yields
    (start, end, category) [start,end) index triples for each maximal run of constant
    NONZERO category -- splitting a run at a category change (e.g. pier box -> pier
    end-cap) is essential: those are geometrically adjacent and visibility-contiguous
    (same face plane) but must not be merged into one quad with one color, or one category
    silently swallows the other (this is exactly the bug that made the pier box vanish in
    an earlier version of this script).

    Also caps each yielded run at `max_run` cells, splitting a longer one into several --
    NOT for correctness (either way covers the exact same solid region) but because
    matplotlib's mplot3d has no real depth buffer: it sorts whole polygons back-to-front by
    centroid distance, which visibly fails for a few giant polygons spanning a huge depth
    range (a large arch wall was painting over a nearer pier cap it should've been behind
    in parts). Many medium polygons give that painter's-algorithm sort far less room to get
    the order wrong, at a modest cost in polygon count versus fully unbounded merging."""
    if not row_cat.any():
        return
    n = len(row_cat)
    change = np.flatnonzero(np.diff(row_cat.astype(np.int16))) + 1
    boundaries = np.concatenate(([0], change, [n]))
    for s, e in zip(boundaries[:-1], boundaries[1:]):
        v = int(row_cat[s])
        if v == 0:
            continue
        for sub_s in range(s, e, max_run):
            yield int(sub_s), int(min(sub_s + max_run, e)), v


def merged_quads_and_colors(mask: np.ndarray, category: np.ndarray, direction: str,
                             const_axis: int, merge_axis: int, other_axis: int,
                             xe, ye, ze, edges):
    """mask: boolean visibility for this direction. const_axis is the face-normal axis
    (fixed per layer), merge_axis is the in-plane axis runs get merged along, other_axis is
    the remaining in-plane axis (kept at unit resolution, looped row by row). Returns
    (quads (N,4,3), colors (N,) object array of hex strings)."""
    quads = []
    colors = []
    # Layers with ANY visible face on this direction -- for our extruded shapes this is a
    # small set for -y/+y/-z/+z (piers/arches are uniform along most of their own extent),
    # but can be large for -x/+x where the circular pier caps live (expected, not a bug).
    axes_order = [0, 1, 2]
    layer_idx = [np.moveaxis(mask, const_axis, 0).any(axis=(1, 2))]
    active_layers = np.flatnonzero(layer_idx[0])
    for c in active_layers:
        idx = [slice(None)] * 3
        idx[const_axis] = c
        layer_mask = mask[tuple(idx)]    # 2D, axes = the two remaining dims in original order
        layer_cat = category[tuple(idx)]  # same shape, the category at each of those cells
        remaining_axes = [a for a in axes_order if a != const_axis]
        merge_pos = remaining_axes.index(merge_axis)
        if merge_pos != 0:
            layer_mask = layer_mask.T  # make merge_axis the first axis of both...
            layer_cat = layer_cat.T    # ...so row indexing below lines up
        n_other = layer_mask.shape[1]
        cat_idx = [slice(None)] * 3
        cat_idx[const_axis] = c
        for o in range(n_other):
            row_mask = layer_mask[:, o]
            if not row_mask.any():
                continue
            row_cat = np.where(row_mask, layer_cat[:, o], 0)
            oc = cat_idx.copy()
            oc[other_axis] = o
            for s, e, cat_val in _runs_1d_by_category(row_cat):
                color = CATEGORY_COLORS.get(cat_val)
                if color is None:
                    continue
                lo = [0, 0, 0]
                hi = [0, 0, 0]
                lo[const_axis] = hi[const_axis] = c
                lo[other_axis] = hi[other_axis] = o
                lo[merge_axis], hi[merge_axis] = s, e
                quads.append(_quad(direction, const_axis, merge_axis, other_axis,
                                    lo, hi, xe, ye, ze))
                colors.append(color)
    if not quads:
        return np.empty((0, 4, 3)), np.empty((0,), dtype=object)
    return np.stack(quads, axis=0), np.array(colors, dtype=object)


def _quad(direction, const_axis, merge_axis, other_axis, lo, hi, xe, ye, ze):
    edges = [xe, ye, ze]
    is_hi_face = direction[0] == "+"
    c_edge = edges[const_axis][lo[const_axis] + (1 if is_hi_face else 0)]
    o_edge0 = edges[other_axis][lo[other_axis]]
    o_edge1 = edges[other_axis][lo[other_axis] + 1]
    m_edge0 = edges[merge_axis][lo[merge_axis]]
    m_edge1 = edges[merge_axis][hi[merge_axis]]  # hi is already exclusive-end cell index
    corners_axes = {const_axis: c_edge}
    pts = []
    for (m, o) in [(m_edge0, o_edge0), (m_edge1, o_edge0), (m_edge1, o_edge1), (m_edge0, o_edge1)]:
        p = [0.0, 0.0, 0.0]
        p[const_axis] = c_edge
        p[merge_axis] = m
        p[other_axis] = o
        pts.append(p)
    # Winding order doesn't matter for matplotlib's unlit Poly3DCollection fill.
    return np.array(pts)


def load_and_render(path: str = "bridge_voxels_full.npz"):
    print(f"[render] loading {path} ...", flush=True)
    t0 = time.time()
    data = np.load(path)
    category = data["category"]
    x0, x1, y0, y1, z0, z1 = (float(data[k]) for k in ("x0", "x1", "y0", "y1", "z0", "z1"))
    voxel_size = float(data["voxel_size"])
    nx, ny, nz = category.shape
    xe = np.arange(x0, x1 + voxel_size, voxel_size)
    ye = np.arange(y0, y1 + voxel_size, voxel_size)
    ze = np.arange(z0, z1 + voxel_size, voxel_size)
    print(f"[render] loaded {nx}x{ny}x{nz} grid in {time.time()-t0:.1f}s")

    occ_all = category > 0
    print(f"[render] computing visible faces against the COMBINED solid region "
          f"({int(occ_all.sum()):,} solid cells total) ...", flush=True)
    t1 = time.time()
    faces = visible_faces(occ_all)
    print(f"[render] face-visibility masks computed in {time.time()-t1:.1f}s")

    # direction -> (const_axis, merge_axis, other_axis). merge_axis picked as the axis each
    # direction's flat runs are typically longest along, for max reduction.
    dir_axes = {
        "-x": (0, 2, 1), "+x": (0, 2, 1),   # merge along Z, loop Y
        "-y": (1, 0, 2), "+y": (1, 0, 2),   # merge along X, loop Z
        "-z": (2, 0, 1), "+z": (2, 0, 1),   # merge along X, loop Y
    }

    all_quads, all_colors = [], []
    t2 = time.time()
    for direction, mask in faces.items():
        const_axis, merge_axis, other_axis = dir_axes[direction]
        n_before = int(mask.sum())
        if n_before == 0:
            continue
        q, c = merged_quads_and_colors(mask, category, direction, const_axis, merge_axis,
                                        other_axis, xe, ye, ze, (xe, ye, ze))
        print(f"[render]   {direction}: {n_before:,} unit faces -> {len(q):,} merged quads", flush=True)
        all_quads.append(q)
        all_colors.append(c)
    quads = np.concatenate(all_quads, axis=0)
    colors = np.concatenate(all_colors, axis=0)
    print(f"[render] total: {len(quads):,} merged quads (from {sum(int(m.sum()) for m in faces.values()):,} "
          f"unit faces), {time.time()-t2:.1f}s to merge")

    fig = plt.figure(figsize=(12, 8.5))
    ax = fig.add_subplot(projection="3d")
    print("[render] building Poly3DCollection ...", flush=True)
    t3 = time.time()
    coll = Poly3DCollection(quads, facecolor=colors, edgecolor="#241f1a15", linewidths=0.1, alpha=1.0)
    ax.add_collection3d(coll)
    print(f"[render] collection built in {time.time()-t3:.1f}s")

    ax.set_xlabel("X (flow direction)")
    ax.set_ylabel("Y (cross-channel)")
    ax.set_zlabel("Z (height)")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_zlim(z0, z1)
    ax.set_box_aspect((x1 - x0, y1 - y0, z1 - z0))
    ax.set_title("Bridge NoSlip cells at true simulation resolution (1 LU/cell)")

    legend_handles = [
        Patch(facecolor=CATEGORY_COLORS[1], edgecolor="none", label="Pier box (Pillar)"),
        Patch(facecolor=CATEGORY_COLORS[2], edgecolor="none",
              label="Pier end-caps (CylinderPillar, Z-infinite for floaters)"),
        Patch(facecolor=CATEGORY_COLORS[3], edgecolor="none", label="Deck / arch collar"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(0.0, 0.95), fontsize=9)

    fig.tight_layout()
    return fig


if __name__ == "__main__":
    print(f"[render] backend: {matplotlib.get_backend()}")
    fig = load_and_render()
    print("[render] saving static snapshot to bridge_voxels_full.png ...", flush=True)
    fig.savefig("bridge_voxels_full.png", dpi=150)
    print("[render] saved. Opening interactive window ...", flush=True)
    plt.show()
