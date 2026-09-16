#!/usr/bin/env python3
"""Visualize the spherosimplex ("rounded box") geometry behind useSpheroBoxContact.

A rounded box is the Minkowski sum of a sharp "core" box (half-extents b) with a sphere
of radius r: outer half-extents H = b + r. Its surface is built here from three exact
analytic patch families rather than a voxelised isosurface, so there is no approximation
error and no extra dependency beyond numpy/matplotlib:

  - 6 flat rectangular faces (the box's faces, inset by r on the two other axes)
  - 12 quarter-cylinder edges of radius r (where two faces would otherwise meet at a
    sharp edge)
  - 8 spherical-octant corners of radius r (where three faces would otherwise meet at a
    sharp vertex)

This is the exact same construction ParticleKernels.h's useSpheroBoxContact relies on
(mesa_pd's GJK/EPA run with margin=r on the sharp core box) -- this script just draws the
resulting outer surface directly, it does not call any C++ code.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

FACE_COLOR = "#4C8BF5"
CORE_COLOR = "#1B1B1B"
EDGE_ALPHA = 0.95


def rounded_box_patches(half_extents, radius, n=14):
    """Return a list of (X, Y, Z) meshgrid triples, one per surface patch, describing the
    rounded box with given outer half-extents and corner/edge radius. n controls the
    angular resolution of the curved (edge/corner) patches."""
    hx, hy, hz = half_extents
    bx, by, bz = hx - radius, hy - radius, hz - radius
    assert min(bx, by, bz) > 0, "radius must be smaller than every half-extent"
    r = radius
    patches = []

    # --- 6 flat faces ---
    for axis in range(3):
        for sign in (-1.0, 1.0):
            b = [bx, by, bz]
            u = np.linspace(-b[(axis + 1) % 3], b[(axis + 1) % 3], 2)
            v = np.linspace(-b[(axis + 2) % 3], b[(axis + 2) % 3], 2)
            U, V = np.meshgrid(u, v)
            coord = [None, None, None]
            coord[axis] = np.full_like(U, sign * half_extents[axis])
            coord[(axis + 1) % 3] = U
            coord[(axis + 2) % 3] = V
            patches.append(tuple(coord))

    # --- 12 quarter-cylinder edges (one per axis-pair x sign-combination) ---
    theta = np.linspace(0, np.pi / 2, n)
    axis_pairs = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]  # (axisA, axisB, axisAlong)
    b = [bx, by, bz]
    for axisA, axisB, axisAlong in axis_pairs:
        length_half = b[axisAlong]
        t = np.linspace(-length_half, length_half, 2)
        Theta, T = np.meshgrid(theta, t)
        for sA in (-1.0, 1.0):
            for sB in (-1.0, 1.0):
                coord = [None, None, None]
                coord[axisA] = sA * (b[axisA] + r * np.cos(Theta))
                coord[axisB] = sB * (b[axisB] + r * np.sin(Theta))
                coord[axisAlong] = T
                patches.append(tuple(coord))

    # --- 8 spherical-octant corners ---
    u = np.linspace(0, np.pi / 2, n)  # polar angle within the octant
    v = np.linspace(0, np.pi / 2, n)  # azimuthal angle within the octant
    U, V = np.meshgrid(u, v)
    nx, ny, nz = np.sin(U) * np.cos(V), np.sin(U) * np.sin(V), np.cos(U)
    for sx in (-1.0, 1.0):
        for sy in (-1.0, 1.0):
            for sz in (-1.0, 1.0):
                X = sx * (bx + r * nx)
                Y = sy * (by + r * ny)
                Z = sz * (bz + r * nz)
                patches.append((X, Y, Z))

    return patches


def core_box_edges(half_extents):
    """12 line segments (as (xs, ys, zs) triples) tracing the sharp core box's wireframe."""
    hx, hy, hz = half_extents
    verts = np.array([[sx * hx, sy * hy, sz * hz]
                       for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
    edges = []
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(verts[i] != verts[j]) == 1:  # differ in exactly one coordinate
                edges.append((verts[i], verts[j]))
    return edges


def draw_rounded_box(ax, half_extents, radius, core_half_extents=None, face_color=FACE_COLOR,
                      alpha=0.55, show_core=True, n=16):
    for X, Y, Z in rounded_box_patches(half_extents, radius, n=n):
        ax.plot_surface(X, Y, Z, color=face_color, alpha=alpha, linewidth=0,
                         antialiased=True, shade=True)
    if show_core:
        core = core_half_extents if core_half_extents is not None else \
            tuple(h - radius for h in half_extents)
        for p, q in core_box_edges(core):
            ax.plot(*zip(p, q), color=CORE_COLOR, linewidth=1.2, linestyle="--", alpha=EDGE_ALPHA)


def set_equal_aspect(ax, limit):
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_box_aspect((1, 1, 1))


def unit_sphere_mesh(center, radius, n=24):
    u = np.linspace(0, 2 * np.pi, n)
    v = np.linspace(0, np.pi, n)
    U, V = np.meshgrid(u, v)
    X = center[0] + radius * np.cos(U) * np.sin(V)
    Y = center[1] + radius * np.sin(U) * np.sin(V)
    Z = center[2] + radius * np.cos(V)
    return X, Y, Z


def core_box_surface_samples(half, n=25):
    """Dense sample of points lying ON the sharp core box's surface (all 6 faces)."""
    b = np.array(half)
    pts = []
    lin = [np.linspace(-b[i], b[i], n) for i in range(3)]
    for axis in range(3):
        for sign in (-1.0, 1.0):
            other = [i for i in range(3) if i != axis]
            U, V = np.meshgrid(lin[other[0]], lin[other[1]])
            P = np.zeros((U.size, 3))
            P[:, axis] = sign * b[axis]
            P[:, other[0]] = U.ravel()
            P[:, other[1]] = V.ravel()
            pts.append(P)
    return np.concatenate(pts, axis=0)


def rot_z(deg):
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def transform_patches(patches, R, pos):
    out = []
    for X, Y, Z in patches:
        pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=0)
        pts = R @ pts + pos.reshape(3, 1)
        out.append((pts[0].reshape(X.shape), pts[1].reshape(Y.shape), pts[2].reshape(Z.shape)))
    return out


def main():
    fig = plt.figure(figsize=(12, 11))

    half_extents = (3.0, 2.0, 1.5)   # illustrative production-shaped floater (12x8x6 full size)
    radius = 0.6                      # exaggerated vs. the production default (0.05) for visibility
    limit = max(half_extents) + radius + 0.3
    core_half = tuple(h - radius for h in half_extents)

    # Panel 1: the sharp "core" box alone.
    ax1 = fig.add_subplot(2, 2, 1, projection="3d")
    for p, q in core_box_edges(core_half):
        ax1.plot(*zip(p, q), color=CORE_COLOR, linewidth=2)
    ax1.set_title("1. Sharp core box\n(what boxBoxSAT sees)")
    set_equal_aspect(ax1, limit)

    # Panel 2: ONE sphere, radius r, swept over the box's surface -- three example positions
    # (a face centre, an edge midpoint, a vertex) shown as translucent copies of the exact
    # SAME sphere, to make clear this is a single sphere swept everywhere over the box, not
    # a separate/different sphere glued on at each corner. The final rounded surface (drawn
    # faintly) is the outer envelope of the union of ALL such copies, one centred at every
    # point of the box's surface -- these three are just illustrative samples of that
    # continuum.
    ax2 = fig.add_subplot(2, 2, 2, projection="3d")
    bx, by, bz = core_half
    sample_centers = {
        "face centre": np.array([0.0, 0.0, bz]),     # middle of the TOP face
        "edge midpoint": np.array([0.0, -by, -bz]),   # a BOTTOM-front edge, far from the above
        "vertex": np.array([-bx, by, -bz]),           # a corner far from both of the above
    }
    for p, q in core_box_edges(core_half):
        ax2.plot(*zip(p, q), color=CORE_COLOR, linewidth=1.5)
    for label, c in sample_centers.items():
        X, Y, Z = unit_sphere_mesh(c, radius, n=22)
        ax2.plot_surface(X, Y, Z, color="#E24A4A", alpha=0.35, linewidth=0, shade=True)
        tip = c + radius * c / np.linalg.norm(c) * 0.0  # label at sphere centre
        ax2.text(c[0], c[1], c[2], f"  {label}", fontsize=8, color="#8A1E1E")
    ax2.set_title("2. ONE sphere (radius r), copied to\n3 example points on the box surface\n"
                  "-- same sphere everywhere, not one per corner")
    set_equal_aspect(ax2, limit)

    # Panel 3: the rounded box (spherosimplex) -- the envelope of copies of that SAME sphere
    # centred at EVERY point of the core box's surface (not just the 3 samples above).
    ax3 = fig.add_subplot(2, 2, 3, projection="3d")
    draw_rounded_box(ax3, half_extents, radius, core_half_extents=core_half)
    ax3.set_title(f"3. Spherosimplex = envelope of that ONE\nsphere swept over the WHOLE core surface\n"
                  f"(core ⊕ sphere(r={radius}), what useSpheroBoxContact sees)")
    set_equal_aspect(ax3, limit)

    # Panel 4: two rounded boxes in GENUINE contact (rounded surfaces actually overlapping),
    # illustrating the continuous closest-point normal -- the vector between the two closest
    # points on the CORE boxes is exactly the reported contact normal, continuous through any
    # face/edge/vertex transition.
    ax4 = fig.add_subplot(2, 2, 4, projection="3d")
    halfA = (2.2, 0.9, 0.7)
    halfB = (2.2, 0.9, 0.7)
    rA = rB = 0.35

    posA = np.array([-3.2, 0.0, 0.0])
    posB = np.array([3.4, 0.55, 0.15])
    RA = rot_z(8)
    RB = rot_z(-14)

    # Solve for the actual gap along the line connecting the two CORE boxes' closest points,
    # then shift box B so the ROUNDED (outer) surfaces overlap by a small, deliberate amount
    # instead of just eyeballing a position -- guarantees the picture shows genuine contact.
    # Iterated a few times since moving the box can change WHICH points are closest.
    sampleA_core = (RA @ core_box_surface_samples(halfA).T).T + posA
    desired_overlap = 0.15  # how deep into each other the ROUNDED surfaces should sit
    for _ in range(6):
        sampleB0 = (RB @ core_box_surface_samples(halfB).T).T + posB
        d0 = np.linalg.norm(sampleA_core[:, None, :] - sampleB0[None, :, :], axis=2)
        ia0, ib0 = np.unravel_index(np.argmin(d0), d0.shape)
        direction = sampleB0[ib0] - sampleA_core[ia0]
        direction /= np.linalg.norm(direction)
        current_core_gap = d0[ia0, ib0]
        shift = current_core_gap - (rA + rB) + desired_overlap
        posB = posB - direction * shift

    for X, Y, Z in transform_patches(rounded_box_patches(halfA, rA, n=14), RA, posA):
        ax4.plot_surface(X, Y, Z, color="#4C8BF5", alpha=0.6, linewidth=0, shade=True)
    for X, Y, Z in transform_patches(rounded_box_patches(halfB, rB, n=14), RB, posB):
        ax4.plot_surface(X, Y, Z, color="#F5824C", alpha=0.6, linewidth=0, shade=True)

    sampleA = (RA @ core_box_surface_samples(halfA).T).T + posA
    sampleB = (RB @ core_box_surface_samples(halfB).T).T + posB
    dists = np.linalg.norm(sampleA[:, None, :] - sampleB[None, :, :], axis=2)
    ia, ib = np.unravel_index(np.argmin(dists), dists.shape)
    pA, pB = sampleA[ia], sampleB[ib]
    coreGap = dists[ia, ib]
    print(f"panel 4: core-box gap={coreGap:.3f}, rA+rB={rA + rB:.3f} "
          f"-> rounded surfaces {'OVERLAP' if coreGap < rA + rB else 'separated'} by "
          f"{(rA + rB) - coreGap:.3f}")

    n = (pB - pA) / np.linalg.norm(pB - pA)
    outerA = pA + n * rA   # point on rounded surface A along the contact normal
    outerB = pB - n * rB   # point on rounded surface B along the contact normal
    ax4.plot(*zip(pA, pB), color="black", linewidth=1, linestyle=":")
    ax4.plot(*zip(outerA, outerB), color="black", linewidth=2.5, marker="o", markersize=5)
    mid = (outerA + outerB) / 2
    ax4.text(mid[0], mid[1], mid[2], "  normal", fontsize=9, fontweight="bold")

    ax4.set_title("4. Two spherosimplices IN CONTACT\nsolid line = contact normal (between the\n"
                  "rounded surfaces); dotted = between the CORE boxes")
    set_equal_aspect(ax4, 3.6)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.xaxis.pane.fill = False; ax.yaxis.pane.fill = False; ax.zaxis.pane.fill = False

    fig.suptitle("useSpheroBoxContact: the spherosimplex (\"rounded box\") geometry", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_path = "spherobox_visualization.png"
    fig.savefig(out_path, dpi=160)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
