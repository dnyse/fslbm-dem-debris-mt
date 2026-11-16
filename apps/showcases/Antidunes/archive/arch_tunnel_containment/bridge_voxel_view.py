#!/usr/bin/env python3
"""3D voxel view of the bridge's NoSlip cells -- shows what actually gets painted onto the
lattice grid (piers via CylinderPillars/Pillars, arch collars via ArchTunnels), not just the
abstract analytic shape. Reuses the exact geometry from bridge_layout_optimizer.py (same
Pier/ArchTunnel classes, same BRIDGE definition) and classifies each voxel with the SAME
boolean composition the real BodyBC blocks use:

    pier box:      Box                                            -> grey   (Pillar)
    pier end-caps: Cylinder_a OR Cylinder_b                        -> terracotta (CylinderPillar --
                    Z-infinite for floater collision in the real DEM model, see Utility.h's
                    createCylinder(); only capped here on the LBM/BodyBC side)
    arch collar:   Box[xMin..xMax,...] AND NOT Cylinder             -> slate  (ArchTunnel)

Rendered with matplotlib's mplot3d ax.voxels() in a real interactive window (drag to orbit,
scroll to zoom) -- no new dependencies beyond what bridge_layout_optimizer.py already uses.

The real simulation resolution is 1 LU/cell; rendering every cell in the bridge's own
bounding region (~330 x 1056 x 350 LU = ~120M cells) is not feasible for an interactive
voxel plot, so this renders at a coarser VOXEL_SIZE (LU per display-voxel, default 8 --
tune it, and REGION below, to zoom into a specific seam at finer resolution if you want to
check a particular join closely).

Usage:
    python3 bridge_voxel_view.py
Edit VOXEL_SIZE / REGION near the bottom and rerun.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from bridge_layout_optimizer import BRIDGE, LU_PER_METER, BridgeLayout


def classify_voxels(layout: BridgeLayout, deck_top: float, region, voxel_size: float):
    """region: (x0,x1,y0,y1,z0,z1) in LU. Returns (occupied, colors, edges) -- boolean/RGBA
    arrays shaped for ax.voxels(), plus the coordinate arrays used to build them."""
    x0, x1, y0, y1, z0, z1 = region
    xs = np.arange(x0, x1, voxel_size)
    ys = np.arange(y0, y1, voxel_size)
    zs = np.arange(z0, z1, voxel_size)
    # Sample at each voxel's CENTER (voxel_size/2 offset) against the exact analytic shape --
    # faithful to the real BodyBC test, just coarser than the true 1-LU cell grid.
    cx, cy, cz = (xs + voxel_size / 2, ys + voxel_size / 2, zs + voxel_size / 2)
    X, Y, Z = np.meshgrid(cx, cy, cz, indexing="ij")

    # Tracked separately (not just "pier" vs "arch") so the two DEM shapes that make up a
    # pier -- CylinderPillar (the round end-caps, uncapped/infinite in Z for floater
    # collision -- see createCylinder()'s own comment in Utility.h) and Pillar (the
    # connecting box, genuinely capped) -- can be told apart visually, since they behave
    # differently for a floater even though both look the same on the LBM/BodyBC side.
    pier_box = np.zeros(X.shape, dtype=bool)
    pier_cyl = np.zeros(X.shape, dtype=bool)
    for p in layout.piers:
        pier_box |= ((p.x_center - p.x_half <= X) & (X <= p.x_center + p.x_half) &
                     (p.y_min <= Y) & (Y <= p.y_max) & (0.0 <= Z) & (Z <= p.z_top))
        in_cyl_a = ((X - p.x_cyl_a) ** 2 + (Y - p.y_center) ** 2 <= p.radius ** 2) & (0.0 <= Z) & (Z <= p.z_top)
        in_cyl_b = ((X - p.x_cyl_b) ** 2 + (Y - p.y_center) ** 2 <= p.radius ** 2) & (0.0 <= Z) & (Z <= p.z_top)
        pier_cyl |= in_cyl_a | in_cyl_b

    arch_solid = np.zeros(X.shape, dtype=bool)
    for a in layout.arches:
        in_box = ((a.x_min <= X) & (X <= a.x_max) & (a.y_min <= Y) & (Y <= a.y_max) &
                  (a.z_center <= Z) & (Z <= deck_top))
        outside_bore = (Y - a.y_center) ** 2 + (Z - a.z_center) ** 2 > a.radius ** 2
        arch_solid |= in_box & outside_bore

    occupied = pier_box | pier_cyl | arch_solid
    colors = np.empty(occupied.shape, dtype=object)
    # Priority: cylinder end-caps drawn over box (both belong to the same pier, this just
    # picks which color wins where they overlap), box over arch.
    colors[arch_solid] = "#5b6b82"   # deck/arch collar: cool slate
    colors[pier_box] = "#8a7f76"     # pier box (Pillar): warm grey
    colors[pier_cyl] = "#c1552b"     # pier end-caps (CylinderPillar): terracotta -- the
                                        # one DEM shape here that's actually Z-infinite
    return occupied, colors, (xs, ys, zs)


def plot_voxels(layout: BridgeLayout, deck_top: float, region, voxel_size: float, lu_per_meter: float):
    occupied, colors, (xs, ys, zs) = classify_voxels(layout, deck_top, region, voxel_size)
    n_filled = int(occupied.sum())
    print(f"Voxel grid: {occupied.shape[0]}x{occupied.shape[1]}x{occupied.shape[2]} "
          f"({occupied.size} total), {n_filled} solid ({n_filled/occupied.size*100:.1f}%), "
          f"voxel size {voxel_size:g} LU ({voxel_size/lu_per_meter:.3f} m)")

    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(projection="3d")
    # ax.voxels needs grid EDGES shaped (nx+1, ny+1, nz+1) to bound occupied's (nx,ny,nz)
    # cells -- xs/ys/zs are the LOWER edge of each cell (len == nx/ny/nz), so append one
    # more edge at the end of each axis.
    xe = np.append(xs, xs[-1] + voxel_size)
    ye = np.append(ys, ys[-1] + voxel_size)
    ze = np.append(zs, zs[-1] + voxel_size)
    X, Y, Z = np.meshgrid(xe, ye, ze, indexing="ij")
    ax.voxels(X, Y, Z, occupied, facecolors=colors, edgecolor="#241f1a20", linewidth=0.15)

    x0, x1, y0, y1, z0, z1 = region
    ax.set_xlabel("X (flow direction)")
    ax.set_ylabel("Y (cross-channel)")
    ax.set_zlabel("Z (height)")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_zlim(z0, z1)
    ax.set_box_aspect((x1 - x0, y1 - y0, z1 - z0))  # true relative proportions, not a cube
    ax.set_title(f"Bridge NoSlip cells (voxel size {voxel_size:g} LU = "
                 f"{voxel_size/lu_per_meter:.3f} m/voxel, {n_filled:,} voxels shown)")

    legend_handles = [
        Patch(facecolor="#8a7f76", edgecolor="none", label="Pier box (Pillar)"),
        Patch(facecolor="#c1552b", edgecolor="none",
              label="Pier end-caps (CylinderPillar, Z-infinite for floaters)"),
        Patch(facecolor="#5b6b82", edgecolor="none", label="Deck / arch collar"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(0.0, 0.95), fontsize=9)

    fig.tight_layout()
    return fig


if __name__ == "__main__":
    # Tune these to explore: VOXEL_SIZE trades render responsiveness for how finely you can
    # see a specific seam/join; REGION restricts what's rendered.
    #
    # Default: zoomed-OUT full-bridge view, all 3 arches + both piers, at 8 LU/voxel --
    # confirmed by a 2 LU zoom into the pier1/opening1 seam (see conversation/git history)
    # that this stays visually flush without the coarse-voxelization gap 12 LU showed, while
    # staying light enough to orbit smoothly (~94k solid voxels vs ~29k at 12 LU).
    #
    # For a close-up seam check instead, set e.g.:
    #   VOXEL_SIZE = 1.0  # true simulation resolution
    #   REGION = (1170.0, 1260.0, 305.0, 395.0, 330.0, 352.0)  # pier1/opening1 corner
    VOXEL_SIZE = 8.0
    deck_top = max((p.z_top for p in BRIDGE.piers), default=BRIDGE.domain[2] * 0.4)
    arch_x0 = min((a.x_min for a in BRIDGE.arches), default=0.0)
    arch_x1 = max((a.x_max for a in BRIDGE.arches), default=BRIDGE.domain[0])
    margin = 40.0
    REGION = (arch_x0 - margin, arch_x1 + margin, 0.0, BRIDGE.domain[1], 0.0, deck_top + margin)

    print(f"backend: {matplotlib.get_backend()}")
    fig = plot_voxels(BRIDGE, deck_top, REGION, VOXEL_SIZE, LU_PER_METER)
    fig.savefig("bridge_voxels.png", dpi=150)
    print("Saved static snapshot to bridge_voxels.png")
    plt.show()
