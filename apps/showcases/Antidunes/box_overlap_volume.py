#!/usr/bin/env python3
"""
Proof of concept: exact (analytical) overlap volume of two oriented boxes (OBBs).

Idea
----
Each box is the intersection of 6 half-spaces (one per face).  The overlap of
two boxes is therefore the convex polytope defined by all 12 half-spaces.
We enumerate every candidate vertex of that polytope as the intersection of a
triple of planes, keep the ones that satisfy all 12 half-space constraints,
and take the convex hull of the survivors.  The hull's volume is the exact
overlap volume -- no sampling, no discretization.

This works for boxes of any position, orientation and size.

Requires: numpy, scipy.
"""

from itertools import combinations

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation


class Box:
    """An oriented box (OBB) defined by center, full side lengths and rotation."""

    def __init__(self, center, size, rotation=None):
        """
        center   : (3,) world-space center
        size     : (3,) full side lengths (Lx, Ly, Lz) along the box's local axes
        rotation : scipy Rotation, 3x3 matrix, or None for axis-aligned
        """
        self.center = np.asarray(center, dtype=float)
        self.half = np.asarray(size, dtype=float) / 2.0

        if rotation is None:
            self.R = np.eye(3)
        elif isinstance(rotation, Rotation):
            self.R = rotation.as_matrix()
        else:
            self.R = np.asarray(rotation, dtype=float)

    def half_spaces(self):
        """
        Return the box as 6 half-spaces (normal, offset) with the convention
            normal . x <= offset.
        The columns of R are the box's local axes.
        """
        planes = []
        for i in range(3):
            axis = self.R[:, i]
            h = self.half[i]
            for s in (+1.0, -1.0):
                n = s * axis
                # s*axis.(x - c) <= h   ->   n.x <= h + n.c
                offset = h + n @ self.center
                planes.append((n, offset))
        return planes

    def corners(self):
        """The 8 world-space corners (handy for plotting / sanity checks)."""
        signs = np.array([[sx, sy, sz]
                          for sx in (-1, 1)
                          for sy in (-1, 1)
                          for sz in (-1, 1)], dtype=float)
        local = signs * self.half
        return (self.R @ local.T).T + self.center

    def volume(self):
        return float(np.prod(2.0 * self.half))


def overlap_volume(box_a, box_b, tol=1e-9):
    """
    Exact overlap volume of two boxes.

    Returns (volume, vertices) where `vertices` are the corners of the
    intersection polytope (empty if the boxes do not overlap in 3D).
    """
    planes = box_a.half_spaces() + box_b.half_spaces()
    normals = np.array([p[0] for p in planes])
    offsets = np.array([p[1] for p in planes])

    # A vertex of the intersection lies on 3 of the planes simultaneously.
    verts = []
    for i, j, k in combinations(range(len(planes)), 3):
        A = normals[[i, j, k]]
        b = offsets[[i, j, k]]
        det = np.linalg.det(A)
        if abs(det) < 1e-12:
            continue  # planes are (near) parallel -> no unique intersection
        x = np.linalg.solve(A, b)
        # Keep it only if it satisfies every half-space (i.e. is inside both boxes).
        if np.all(normals @ x <= offsets + tol):
            verts.append(x)

    if len(verts) < 4:
        return 0.0, np.empty((0, 3))

    verts = np.array(verts)
    # Deduplicate near-identical points before hulling.
    verts = _dedupe(verts, tol=1e-7)
    if len(verts) < 4:
        return 0.0, verts

    try:
        hull = ConvexHull(verts)
    except Exception:
        # Degenerate (coplanar) intersection -> zero volume.
        return 0.0, verts

    return float(hull.volume), verts[hull.vertices]


def _dedupe(points, tol=1e-7):
    kept = []
    for p in points:
        if not any(np.linalg.norm(p - q) < tol for q in kept):
            kept.append(p)
    return np.array(kept)


def monte_carlo_overlap(box_a, box_b, n=2_000_000, seed=0):
    """Independent stochastic estimate, used only to validate the analytical result."""
    rng = np.random.default_rng(seed)

    # Sample inside box_a (uniform in its local frame), then test membership in box_b.
    local = (rng.random((n, 3)) - 0.5) * (2.0 * box_a.half)
    pts = (box_a.R @ local.T).T + box_a.center

    # Point-in-box test for box_b: transform into box_b's local frame.
    rel = (box_b.R.T @ (pts - box_b.center).T).T
    inside = np.all(np.abs(rel) <= box_b.half, axis=1)

    frac = inside.mean()
    stderr = np.sqrt(frac * (1 - frac) / n) * box_a.volume()
    return frac * box_a.volume(), stderr


def _hull_mesh_and_edges(points):
    """Triangulated faces (i,j,k) and unique edges of the convex hull of `points`."""
    hull = ConvexHull(points)
    i, j, k = hull.simplices.T
    edges = set()
    for a, b, c in hull.simplices:
        for u, v in ((a, b), (b, c), (c, a)):
            edges.add((min(u, v), max(u, v)))
    return (i, j, k), sorted(edges)


def _add_solid(fig, points, color, name, opacity):
    import plotly.graph_objects as go

    (i, j, k), edges = _hull_mesh_and_edges(points)
    fig.add_trace(go.Mesh3d(
        x=points[:, 0], y=points[:, 1], z=points[:, 2],
        i=i, j=j, k=k,
        color=color, opacity=opacity, flatshading=True,
        name=name, showlegend=True, hoverinfo="name",
    ))
    # Wireframe edges (one polyline with None breaks between segments).
    ex, ey, ez = [], [], []
    for a, b in edges:
        ex += [points[a, 0], points[b, 0], None]
        ey += [points[a, 1], points[b, 1], None]
        ez += [points[a, 2], points[b, 2], None]
    fig.add_trace(go.Scatter3d(
        x=ex, y=ey, z=ez, mode="lines",
        line=dict(color=color, width=4),
        name=f"{name} edges", showlegend=False, hoverinfo="skip",
    ))


def build_figure(box_a, box_b, title=None, scene_range=None):
    """Build the 3D plotly figure of both boxes and their exact overlap polytope."""
    import plotly.graph_objects as go

    vol, verts = overlap_volume(box_a, box_b)
    fig = go.Figure()
    _add_solid(fig, box_a.corners(), "royalblue", "Box A", opacity=0.25)
    _add_solid(fig, box_b.corners(), "firebrick", "Box B", opacity=0.25)
    if len(verts) >= 4:
        _add_solid(fig, verts, "gold", "Overlap", opacity=0.85)

    if title is None:
        title = f"overlap volume = {vol:.6f}"
    scene = dict(aspectmode="cube" if scene_range else "data",
                 xaxis_title="x", yaxis_title="y", zaxis_title="z")
    if scene_range:
        lo, hi = scene_range
        scene["xaxis"] = dict(range=[lo, hi])
        scene["yaxis"] = dict(range=[lo, hi])
        scene["zaxis"] = dict(range=[lo, hi])
    fig.update_layout(
        title=title,
        scene=scene,
        legend=dict(itemsizing="constant"),
        margin=dict(l=0, r=0, t=40, b=0),
        uirevision="keep",  # keep camera pose across recomputes
    )
    return fig, vol


def visualize(name, box_a, box_b, filename=None, show=True):
    """Write/show a static interactive 3D view of two boxes and their overlap."""
    fig, vol = build_figure(box_a, box_b, title=f"{name} — overlap volume = {vol_fmt(box_a, box_b)}")
    if filename:
        fig.write_html(filename, include_plotlyjs="cdn", auto_open=False)
        print(f"  wrote interactive view -> {filename}")
    if show:
        fig.show()
    return fig


def vol_fmt(box_a, box_b):
    return f"{overlap_volume(box_a, box_b)[0]:.6f}"


def run_app(box_a=None, box_b=None, port=8050):
    """
    Launch an interactive Dash app: sliders move/rotate both boxes and the
    overlap volume + polytope are recomputed live on every change.
    """
    import dash
    from dash import dcc, html, Input, Output
    import plotly.graph_objects as go  # noqa: F401 (used indirectly via build_figure)

    if box_a is None:
        # Water body: top surface (water line) at z = 0, depth 2.
        box_a = Box(center=[0.0, 0.0, -1.0], size=[6.0, 6.0, 2.0])
    if box_b is None:
        # Small floater centered on the water line -> exactly half submerged.
        box_b = Box(center=[0.0, 0.0, 0.0], size=[1.0, 1.0, 1.0])

    # Fixed sizes for this session; sliders control pose only.
    size_a = 2.0 * box_a.half
    size_b = 2.0 * box_b.half
    rng = (-4.0, 4.0)

    def slider(sid, label, lo, hi, val, step):
        return html.Div([
            html.Label(label, style={"fontSize": "12px"}),
            dcc.Slider(lo, hi, step, value=val, id=sid,
                       marks=None, tooltip={"placement": "bottom",
                                            "always_visible": False}),
        ], style={"margin": "2px 0"})

    def pose_controls(prefix, box, color, with_size=False):
        c = box.center
        s = 2.0 * box.half
        children = [
            html.H4(f"{prefix.upper()} box", style={"color": color, "margin": "4px 0"}),
            slider(f"{prefix}-cx", "center x", -3, 3, c[0], 0.05),
            slider(f"{prefix}-cy", "center y", -3, 3, c[1], 0.05),
            slider(f"{prefix}-cz", "center z", -3, 3, c[2], 0.05),
            slider(f"{prefix}-rx", "rot x (deg)", -180, 180, 0, 1),
            slider(f"{prefix}-ry", "rot y (deg)", -180, 180, 0, 1),
            slider(f"{prefix}-rz", "rot z (deg)", -180, 180, 0, 1),
        ]
        if with_size:
            children += [
                slider(f"{prefix}-sx", "size x", 0.1, 4, s[0], 0.05),
                slider(f"{prefix}-sy", "size y", 0.1, 4, s[1], 0.05),
                slider(f"{prefix}-sz", "size z", 0.1, 4, s[2], 0.05),
            ]
        return html.Div(children,
                        style={"padding": "6px", "border": f"2px solid {color}",
                               "borderRadius": "8px", "margin": "6px"})

    app = dash.Dash(__name__)
    app.layout = html.Div([
        html.Div([
            pose_controls("a", box_a, "royalblue"),
            pose_controls("b", box_b, "firebrick", with_size=True),
            html.Div(id="readout",
                     style={"fontSize": "16px", "fontWeight": "bold",
                            "padding": "8px", "margin": "6px"}),
        ], style={"width": "300px", "float": "left"}),
        dcc.Graph(id="scene", style={"height": "90vh", "marginLeft": "310px"}),
    ])

    slider_ids = ([f"a-{ax}" for ax in ("cx", "cy", "cz", "rx", "ry", "rz")]
                  + [f"b-{ax}" for ax in ("cx", "cy", "cz", "rx", "ry", "rz",
                                          "sx", "sy", "sz")])

    @app.callback(
        Output("scene", "figure"), Output("readout", "children"),
        [Input(sid, "value") for sid in slider_ids],
    )
    def _update(*values):
        v = dict(zip(slider_ids, values))
        ba = Box([v["a-cx"], v["a-cy"], v["a-cz"]], size_a,
                 Rotation.from_euler("xyz", [v["a-rx"], v["a-ry"], v["a-rz"]],
                                     degrees=True))
        bb = Box([v["b-cx"], v["b-cy"], v["b-cz"]],
                 [v["b-sx"], v["b-sy"], v["b-sz"]],
                 Rotation.from_euler("xyz", [v["b-rx"], v["b-ry"], v["b-rz"]],
                                     degrees=True))
        fig, vol = build_figure(ba, bb, scene_range=rng)
        pct_a = 100.0 * vol / ba.volume()
        pct_b = 100.0 * vol / bb.volume()
        text = (f"Overlap volume = {vol:.6f}  "
                f"({pct_a:.1f}% of A, {pct_b:.1f}% of B)")
        return fig, text

    print(f"Interactive app running at http://127.0.0.1:{port}  (Ctrl-C to stop)")
    app.run(debug=False, port=port)


def _demo(name, box_a, box_b):
    vol, verts = overlap_volume(box_a, box_b)
    mc, err = monte_carlo_overlap(box_a, box_b)
    print(f"=== {name} ===")
    print(f"  box A volume        : {box_a.volume():.6f}")
    print(f"  box B volume        : {box_b.volume():.6f}")
    print(f"  analytical overlap  : {vol:.6f}   ({len(verts)} polytope vertices)")
    print(f"  monte-carlo overlap : {mc:.6f} +/- {err:.6f}")
    print()


if __name__ == "__main__":
    import sys

    if "--app" in sys.argv:
        run_app()
        sys.exit(0)

    # 1) Two axis-aligned unit cubes offset by 0.5 in each axis -> overlap 0.5^3.
    _demo(
        "axis-aligned, partial overlap (expect 0.125)",
        Box(center=[0, 0, 0], size=[1, 1, 1]),
        Box(center=[0.5, 0.5, 0.5], size=[1, 1, 1]),
    )

    # 2) A rotated cube fully containing a smaller one -> overlap = smaller volume.
    _demo(
        "small cube fully inside a big rotated cube (expect 1.0)",
        Box(center=[0, 0, 0], size=[1, 1, 1]),
        Box(center=[0, 0, 0], size=[5, 5, 5],
            rotation=Rotation.from_euler("xyz", [30, 45, 60], degrees=True)),
    )

    # 3) Two arbitrarily oriented, differently sized boxes.
    box_a = Box(center=[0.2, -0.1, 0.3], size=[2.0, 1.0, 1.5],
                rotation=Rotation.from_euler("xyz", [20, 35, 10], degrees=True))
    box_b = Box(center=[0.5, 0.4, -0.2], size=[1.2, 2.5, 0.8],
                rotation=Rotation.from_euler("zyx", [50, 15, 70], degrees=True))
    _demo("arbitrary orientations and sizes", box_a, box_b)

    # 4) Disjoint boxes -> zero overlap.
    _demo(
        "disjoint (expect 0.0)",
        Box(center=[0, 0, 0], size=[1, 1, 1]),
        Box(center=[10, 0, 0], size=[1, 1, 1]),
    )

    # 5) Small floater half-submerged in a water body (expect 50% of floater).
    water = Box(center=[0.0, 0.0, -1.0], size=[6.0, 6.0, 2.0])
    floater = Box(center=[0.0, 0.0, 0.0], size=[1.0, 1.0, 1.0])
    _demo("floater 50% submerged in water (expect 0.5)", water, floater)

    # Interactive 3D visualization of the half-submerged case.
    visualize("floater 50% submerged in water", water, floater,
              filename="box_overlap.html", show=False)
