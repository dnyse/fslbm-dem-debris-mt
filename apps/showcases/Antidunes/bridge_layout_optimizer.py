#!/usr/bin/env python3
"""Offline geometric model of the Antidunes real-scale bridge scenario (piers + arch
tunnels), based on Pont ferroviaire, Rue Goffontaine, Trooz (river Vesdre, Belgium) --
EPFL/ULiege bridge-clogging database record #30071, July 2021 flood event.

Recreates the bridge's SHAPE ONLY, purely from the SAME primitives the real simulation
uses -- CylinderPillars + Pillars for piers, ArchTunnels for the arch openings (see
SimulationParams.h's ArchTunnelParams) -- as plain numpy/matplotlib geometry, no waLBerla
involved. Lets you see how the structure sits in the domain and tune pier/arch placement
and dimensions (in lattice units AND real metres) before touching the real simulation,
then export the tuned numbers as prm-ready text to paste into
prms/pont_ferroviere_rue_goffontaine_realscale_gas.prm.

Defaults below are the exact numbers currently in that prm (33 LU/m resolution); its own
header has the full real-value -> LU derivation for every one of them.

Usage:
    python3 bridge_layout_optimizer.py
Edit the BRIDGE block near the bottom and rerun -- no simulation needed.
"""
from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ─────────────────────────────────────────────────────────────────────────────────────────
# Geometry model
# ─────────────────────────────────────────────────────────────────────────────────────────

@dataclass
class Pier:
    """Stadium-shaped pier: a box union two end-cap cylinders, matching CylinderPillars +
    Pillars in the real prm exactly (and the BodyBC LogicalOR that mirrors it on the LBM
    side). Cylinder axes run along Z (vertical); the box/cylinders share one X-center and
    are offset by x_half along X."""
    name: str
    x_center: float
    y_center: float
    x_half: float      # box halfExtents[0]; also the cylinder offset from x_center
    y_half: float       # box halfExtents[1]
    z_top: float          # deck cap height (piers stop here, not at domain top)
    radius: float           # end-cap cylinder radius

    @property
    def x_cyl_a(self) -> float: return self.x_center - self.x_half
    @property
    def x_cyl_b(self) -> float: return self.x_center + self.x_half
    @property
    def y_min(self) -> float: return self.y_center - self.y_half
    @property
    def y_max(self) -> float: return self.y_center + self.y_half
    @property
    def width_y(self) -> float:
        """Overall pier width across the channel (Y): box half-width == cylinder radius
        in every scenario seen so far, so this is just 2*y_half, but computed from both
        so a future asymmetric pier still reports its true footprint."""
        return 2 * max(self.y_half, self.radius)


@dataclass
class ArchTunnel:
    """Cylindrical bore, axis along X (world-X ignored for center per the real
    ArchTunnelParams), clamped to [x_min, x_max]."""
    name: str
    y_center: float
    z_center: float
    radius: float
    x_min: float
    x_max: float

    @property
    def y_min(self) -> float: return self.y_center - self.radius
    @property
    def y_max(self) -> float: return self.y_center + self.radius
    @property
    def crown_z(self) -> float: return self.z_center + self.radius

    def deck_solid_contains(self, x: float, y: float, z: float, deck_top: float) -> bool:
        """True where this arch's own NoSlip collar sits -- mirrors the real BodyBC exactly:
        LogicalAND(Box[xMin..xMax, y_min..y_max, z_center..deck_top], LogicalNOT(Cylinder
        along X, this radius, centered at (y_center,z_center))). I.e. solid deck above the
        tunnel, empty inside the semicircular bore."""
        if not (self.x_min <= x <= self.x_max and self.y_min <= y <= self.y_max
                and self.z_center <= z <= deck_top):
            return False
        return (y - self.y_center) ** 2 + (z - self.z_center) ** 2 > self.radius ** 2


@dataclass
class BridgeLayout:
    domain: Tuple[float, float, float]      # (X, Y, Z)
    waterline_z: float
    piers: List[Pier] = field(default_factory=list)
    arches: List[ArchTunnel] = field(default_factory=list)

    # ---- shape optimization checks --------------------------------------------------------
    def check_channel_tiling(self) -> Tuple[List[str], List[str]]:
        """Piers and arches should tile the channel with no INTERNAL gaps/overlaps between
        consecutive pieces (that's the real structural constraint -- every floater
        approaching the bridge hits either a pier or an opening, never a hole). A gap
        *before the first* or *after the last* band is not itself a bug -- it's bank/
        abutment margin (the real record explicitly reserves ~1 m for this: 31.7 m span -
        3x8.9 m openings - 2x2 m piers = 1 m of margin). Returns (problems, notes)."""
        bands = sorted(
            [(p.y_min, p.y_max, p.name) for p in self.piers] +
            [(a.y_min, a.y_max, a.name) for a in self.arches]
        )
        problems: List[str] = []
        notes: List[str] = []
        if bands[0][0] != 0.0:
            notes.append(f"Y=0..{bands[0][0]:g} left as bank/abutment margin before '{bands[0][2]}'")
        for (_a0, a1, na), (b0, _b1, nb) in zip(bands, bands[1:]):
            if abs(a1 - b0) > 1e-9:
                kind = "overlap" if b0 < a1 else "internal gap"
                problems.append(f"{kind} between '{na}' (ends {a1}) and '{nb}' (starts {b0})")
        if abs(bands[-1][1] - self.domain[1]) > 1e-9:
            notes.append(f"Y={bands[-1][1]:g}..{self.domain[1]:g} left as bank/abutment "
                         f"margin after '{bands[-1][2]}'")
        return problems, notes

    def check_deck_span_consistent(self) -> List[str]:
        """All piers/arches should share the same outer X-footprint (they sit under one
        deck) -- flags any that don't, which is almost always a typo rather than an
        intended design. A pier's OUTER footprint includes its rounded end-caps
        (x_half + radius on each side), not just the connecting box -- the real bridge's
        piers are deliberately inset so their end-caps land flush with the deck edge."""
        spans = {(p.x_center - p.x_half - p.radius, p.x_center + p.x_half + p.radius, p.name)
                 for p in self.piers}
        arch_spans = {(a.x_min, a.x_max, a.name) for a in self.arches}
        all_spans = {(a, b) for a, b, _ in spans | arch_spans}
        if len(all_spans) <= 1:
            return []
        return [f"inconsistent outer X-span across pieces: {sorted(all_spans)}"]

    # ---- dimension summary (LU and real metres) --------------------------------------------
    def dimension_summary(self, lu_per_meter: float) -> None:
        def m(v: float) -> float: return v / lu_per_meter
        print(f"\n--- Dimensions (resolution: {lu_per_meter:g} LU/m) ---")
        print(f"  domain: {self.domain[0]:g} x {self.domain[1]:g} x {self.domain[2]:g} LU  "
              f"= {m(self.domain[0]):.2f} x {m(self.domain[1]):.2f} x {m(self.domain[2]):.2f} m")
        print(f"  waterline Z = {self.waterline_z:g} LU = {m(self.waterline_z):.2f} m")
        for p in self.piers:
            print(f"  {p.name:>10s}: width={p.width_y:7.2f} LU ({m(p.width_y):.2f} m)  "
                  f"Y-span=[{p.y_min:g},{p.y_max:g}]  deck-top Z={p.z_top:g} LU ({m(p.z_top):.2f} m)")
        for a in self.arches:
            print(f"  {a.name:>22s}: opening width={2*a.radius:7.2f} LU ({m(2*a.radius):.2f} m)  "
                  f"Y-span=[{a.y_min:g},{a.y_max:g}]  crown Z={a.crown_z:g} LU ({m(a.crown_z):.2f} m)"
                  f"  springing Z={a.z_center:g} LU ({m(a.z_center):.2f} m)")

    # ---- prm export ------------------------------------------------------------------------
    def to_prm_text(self) -> str:
        lines = ["CylinderPillars", "{"]
        for p in self.piers:
            lines.append(f"   CylinderPillar {{ center < {p.x_cyl_a:g}, {p.y_center:g}, "
                          f"{p.z_top / 2:g} >; radius {p.radius:g}; }}")
            lines.append(f"   CylinderPillar {{ center < {p.x_cyl_b:g}, {p.y_center:g}, "
                          f"{p.z_top / 2:g} >; radius {p.radius:g}; }}")
        lines += ["}", "", "Pillars", "{"]
        for p in self.piers:
            lines.append(f"   Pillar {{ center < {p.x_center:g}, {p.y_center:g}, {p.z_top / 2:g} >; "
                          f"halfExtents < {p.x_half:g}, {p.y_half:g}, {p.z_top / 2:g} >; }}")
        lines += ["}", "", "ArchTunnels", "{"]
        for a in self.arches:
            lines.append(f"   ArchTunnel {{ center < 0, {a.y_center:g}, {a.z_center:g} >;  "
                          f"radius {a.radius:g};  xMin {a.x_min:g}; xMax {a.x_max:g}; }}  // {a.name}")
        lines.append("}")
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────────────────
# Visualization -- bridge shape only, no floaters
# ─────────────────────────────────────────────────────────────────────────────────────────

def plot_layout(layout: BridgeLayout):
    fig, (ax_elev, ax_plan) = plt.subplots(1, 2, figsize=(13, 5.5))
    dom_x, dom_y, dom_z = layout.domain

    # -- Elevation view (Y-Z cross-section, constant along the bridge's X-span) --
    ax_elev.set_title("Elevation (Y-Z cross-section through the bridge)")
    deck_top = max((p.z_top for p in layout.piers), default=dom_z * 0.4)
    for p in layout.piers:
        ax_elev.add_patch(mpatches.Rectangle((p.y_min, 0), p.y_half * 2, deck_top,
                                              facecolor="#8a7f76", edgecolor="none", zorder=1))
    for a in layout.arches:
        y0, y1 = a.y_min, a.y_max
        ax_elev.add_patch(mpatches.Rectangle((y0, a.z_center), y1 - y0, deck_top - a.z_center,
                                              facecolor="#8a7f76", edgecolor="none", zorder=1))
        theta = np.linspace(0, np.pi, 100)
        ax_elev.fill(a.y_center + a.radius * np.cos(theta), a.z_center + a.radius * np.sin(theta),
                     facecolor="white", edgecolor="#241f1a", linewidth=1.2, zorder=2)
    ax_elev.axhline(layout.waterline_z, color="#4a7ba6", linestyle="--", linewidth=1, zorder=3,
                     label=f"waterline Z={layout.waterline_z:g}")
    ax_elev.set_xlim(0, dom_y)
    ax_elev.set_ylim(0, deck_top * 1.1)
    ax_elev.set_xlabel("Y"); ax_elev.set_ylabel("Z"); ax_elev.set_aspect("equal")
    ax_elev.legend(loc="upper right", fontsize=8)

    # -- Plan view (X-Y top-down, full domain) --
    ax_plan.set_title("Plan (X-Y top-down)")
    for a in layout.arches:
        ax_plan.add_patch(mpatches.Rectangle((a.x_min, a.y_min), a.x_max - a.x_min, 2 * a.radius,
                                              facecolor="#4a7ba6", alpha=0.15, edgecolor="none"))
    for p in layout.piers:
        ax_plan.add_patch(mpatches.Rectangle((p.x_center - p.x_half, p.y_min),
                                              2 * p.x_half, 2 * p.y_half,
                                              facecolor="#8a7f76", edgecolor="none"))
        for xc in (p.x_cyl_a, p.x_cyl_b):
            ax_plan.add_patch(mpatches.Circle((xc, p.y_center), p.radius, facecolor="#8a7f76", edgecolor="none"))
    ax_plan.set_xlim(0, dom_x); ax_plan.set_ylim(0, dom_y)
    ax_plan.set_xlabel("X"); ax_plan.set_ylabel("Y"); ax_plan.set_aspect("equal")

    fig.tight_layout()
    return fig


def _nice_ticks(lo: float, hi: float, target: int = 6) -> List[float]:
    """Round-number tick positions covering [lo, hi] -- plain stdlib/math, no matplotlib
    ticker machinery needed for a couple of SVG axes."""
    span = hi - lo
    if span <= 0:
        return [lo]
    raw_step = span / target
    mag = 10 ** np.floor(np.log10(raw_step))
    for mult in (1, 2, 2.5, 5, 10):
        step = mag * mult
        if step >= raw_step:
            break
    start = np.ceil(lo / step) * step
    ticks = []
    v = start
    while v <= hi + 1e-9:
        ticks.append(round(v, 6))
        v += step
    return ticks or [lo, hi]


def to_interactive_html(layout: BridgeLayout, path: str, lu_per_meter: float) -> None:
    """Self-contained, dependency-free interactive view (pan: drag, zoom: mouse wheel,
    hover: native SVG tooltips) of the same elevation + plan drawings plot_layout() makes
    as a static PNG. Pure Python string templating + numpy for the arc points -- no
    plotly/bokeh/mpld3, so it needs nothing beyond what this script already imports."""
    dom_x, dom_y, dom_z = layout.domain
    deck_top = max((p.z_top for p in layout.piers), default=dom_z * 0.4)

    def m(v: float) -> str:
        return f"{v / lu_per_meter:.2f}"

    # -- Elevation panel (Y horizontal, Z vertical, SVG y-axis flipped vs. data) --
    ez_max = deck_top * 1.1
    def ey(y: float) -> float: return y            # SVG x <- data Y directly
    def ez(z: float) -> float: return ez_max - z    # SVG y <- flipped data Z

    elev_shapes = []
    for p in layout.piers:
        elev_shapes.append(
            f'<rect x="{ey(p.y_min):.2f}" y="{ez(deck_top):.2f}" '
            f'width="{p.y_half*2:.2f}" height="{deck_top:.2f}" class="solid">'
            f'<title>{p.name}: width {p.width_y:g} LU ({m(p.width_y)} m), '
            f'deck-top Z={p.z_top:g} LU ({m(p.z_top)} m)</title></rect>')
    for a in layout.arches:
        elev_shapes.append(
            f'<rect x="{ey(a.y_min):.2f}" y="{ez(deck_top):.2f}" '
            f'width="{2*a.radius:.2f}" height="{deck_top - a.z_center:.2f}" class="solid">'
            f'<title>{a.name}: opening width {2*a.radius:g} LU ({m(2*a.radius)} m)</title></rect>')
        theta = np.linspace(0, np.pi, 80)
        pts = " ".join(f"{ey(a.y_center + a.radius*np.cos(t)):.2f},{ez(a.z_center + a.radius*np.sin(t)):.2f}"
                        for t in theta)
        elev_shapes.append(
            f'<polygon points="{pts}" class="void">'
            f'<title>{a.name}: crown Z={a.crown_z:g} LU ({m(a.crown_z)} m), '
            f'springing Z={a.z_center:g} LU ({m(a.z_center)} m)</title></polygon>')
        elev_shapes.append(f'<text x="{ey(a.y_center):.2f}" y="{ez(a.z_center) - 8:.2f}" '
                            f'class="label" text-anchor="middle">{a.name}</text>')
    for p in layout.piers:
        elev_shapes.append(f'<text x="{ey(p.y_center):.2f}" y="{ez(deck_top) + 16:.2f}" '
                            f'class="label label-on-dark" text-anchor="middle">{p.name}</text>')
    elev_shapes.append(f'<line x1="0" y1="{ez(layout.waterline_z):.2f}" x2="{dom_y:.2f}" '
                        f'y2="{ez(layout.waterline_z):.2f}" class="waterline">'
                        f'<title>waterline Z={layout.waterline_z:g} LU ({m(layout.waterline_z)} m)</title></line>')
    for ty in _nice_ticks(0, dom_y):
        elev_shapes.append(f'<line x1="{ey(ty):.2f}" y1="{ez(0):.2f}" x2="{ey(ty):.2f}" '
                            f'y2="{ez(0)+6:.2f}" class="tick"/>')
        elev_shapes.append(f'<text x="{ey(ty):.2f}" y="{ez(0)+20:.2f}" class="axis-label" '
                            f'text-anchor="middle">{ty:g}</text>')
    for tz in _nice_ticks(0, deck_top):
        elev_shapes.append(f'<line x1="0" y1="{ez(tz):.2f}" x2="-6" y2="{ez(tz):.2f}" class="tick"/>')
        elev_shapes.append(f'<text x="-10" y="{ez(tz)+4:.2f}" class="axis-label" '
                            f'text-anchor="end">{tz:g}</text>')

    # -- Plan panel (X horizontal, Y vertical, no flip needed but keep consistent helper) --
    def px(x: float) -> float: return x
    def py(y: float) -> float: return y

    plan_shapes = []
    for a in layout.arches:
        plan_shapes.append(
            f'<rect x="{px(a.x_min):.2f}" y="{py(a.y_min):.2f}" '
            f'width="{a.x_max-a.x_min:.2f}" height="{2*a.radius:.2f}" class="band">'
            f'<title>{a.name}: X-span [{a.x_min:g},{a.x_max:g}] LU</title></rect>')
    for p in layout.piers:
        plan_shapes.append(
            f'<rect x="{px(p.x_center-p.x_half):.2f}" y="{py(p.y_min):.2f}" '
            f'width="{2*p.x_half:.2f}" height="{2*p.y_half:.2f}" class="solid">'
            f'<title>{p.name} box: halfExtents &lt;{p.x_half:g},{p.y_half:g}&gt; LU</title></rect>')
        for xc in (p.x_cyl_a, p.x_cyl_b):
            plan_shapes.append(f'<circle cx="{px(xc):.2f}" cy="{py(p.y_center):.2f}" r="{p.radius:.2f}" '
                                f'class="solid"><title>{p.name} end-cap: radius {p.radius:g} LU</title></circle>')
        plan_shapes.append(f'<text x="{px(p.x_center):.2f}" y="{py(p.y_center)+4:.2f}" '
                            f'class="label label-on-dark" text-anchor="middle">{p.name}</text>')
    for ty in _nice_ticks(0, dom_x):
        plan_shapes.append(f'<line x1="{px(ty):.2f}" y1="{py(0):.2f}" x2="{px(ty):.2f}" '
                            f'y2="{py(0)-6:.2f}" class="tick"/>')
        plan_shapes.append(f'<text x="{px(ty):.2f}" y="{py(0)-10:.2f}" class="axis-label" '
                            f'text-anchor="middle">{ty:g}</text>')
    for tz in _nice_ticks(0, dom_y):
        plan_shapes.append(f'<line x1="{px(0):.2f}" y1="{py(tz):.2f}" x2="{px(0)-6:.2f}" y2="{py(tz):.2f}" class="tick"/>')
        plan_shapes.append(f'<text x="{px(0)-10:.2f}" y="{py(tz)+4:.2f}" class="axis-label" '
                            f'text-anchor="end">{tz:g}</text>')

    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Bridge layout (interactive)</title>
<style>
  body {{ margin:0; padding:20px; font-family:-apple-system,BlinkMacSystemFont,sans-serif;
          background:#f4f2ee; color:#241f1a; }}
  h1 {{ font-size:18px; margin:0 0 4px; }}
  p.hint {{ color:#5a5049; font-size:13px; margin:0 0 16px; }}
  .row {{ display:flex; gap:20px; flex-wrap:wrap; }}
  .panel {{ background:#fff; border:1px solid #ddd5cc; border-radius:6px; padding:10px;
            flex:1 1 480px; min-width:360px; }}
  .panel h2 {{ font-size:14px; margin:0 0 6px; font-weight:600; }}
  .viewport-wrap {{ overflow:hidden; touch-action:none; cursor:grab; }}
  .viewport-wrap:active {{ cursor:grabbing; }}
  svg {{ width:100%; height:auto; display:block; }}
  .solid {{ fill:#8a7f76; stroke:none; }}
  .void {{ fill:#ffffff; stroke:#241f1a; stroke-width:1.2; }}
  .band {{ fill:#4a7ba6; opacity:0.15; stroke:none; }}
  .waterline {{ stroke:#4a7ba6; stroke-width:1.5; stroke-dasharray:6,4; }}
  .tick {{ stroke:#8a7f76; stroke-width:1; }}
  .axis-label {{ font-size:11px; fill:#5a5049; }}
  .label {{ font-size:11px; fill:#241f1a; }}
  .label-on-dark {{ fill:#fff; }}
  button.reset {{ font-size:12px; margin-left:8px; padding:2px 8px; cursor:pointer; }}
</style></head>
<body>
<h1>Bridge layout -- Pont ferroviaire, Rue Goffontaine, Trooz</h1>
<p class="hint">Drag to pan, scroll/pinch to zoom. Hover a pier or opening for its dimensions.
<button class="reset" onclick="resetAll()">Reset view</button></p>
<div class="row">
  <div class="panel">
    <h2>Elevation (Y-Z cross-section through the bridge)</h2>
    <div class="viewport-wrap pz">
      <svg viewBox="-40 -10 {dom_y+60:.0f} {ez_max+40:.0f}" xmlns="http://www.w3.org/2000/svg">
        <g class="viewport">{''.join(elev_shapes)}</g>
      </svg>
    </div>
  </div>
  <div class="panel">
    <h2>Plan (X-Y top-down)</h2>
    <div class="viewport-wrap pz">
      <svg viewBox="-40 -10 {dom_x+60:.0f} {dom_y+40:.0f}" xmlns="http://www.w3.org/2000/svg">
        <g class="viewport">{''.join(plan_shapes)}</g>
      </svg>
    </div>
  </div>
</div>
<script>
  const states = [];
  function resetAll() {{ states.forEach(s => {{ s.x=0; s.y=0; s.k=1; apply(s); }}); }}
  function apply(s) {{ s.g.setAttribute('transform', `translate(${{s.x}},${{s.y}}) scale(${{s.k}})`); }}
  document.querySelectorAll('.pz').forEach(wrap => {{
    const svg = wrap.querySelector('svg');
    const g = wrap.querySelector('.viewport');
    const s = {{ x:0, y:0, k:1, g }};
    states.push(s);
    let dragging = false, lastX = 0, lastY = 0;
    wrap.addEventListener('pointerdown', e => {{ dragging = true; lastX = e.clientX; lastY = e.clientY; wrap.setPointerCapture(e.pointerId); }});
    wrap.addEventListener('pointerup', () => dragging = false);
    wrap.addEventListener('pointermove', e => {{
      if (!dragging) return;
      const vb = svg.viewBox.baseVal; const rect = svg.getBoundingClientRect();
      const scale = vb.width / rect.width;
      s.x += (e.clientX - lastX) * scale; s.y += (e.clientY - lastY) * scale;
      lastX = e.clientX; lastY = e.clientY; apply(s);
    }});
    wrap.addEventListener('wheel', e => {{
      e.preventDefault();
      const factor = e.deltaY < 0 ? 1.1 : 1/1.1;
      s.k *= factor; apply(s);
    }}, {{ passive:false }});
  }});
</script>
</body></html>"""
    with open(path, "w") as f:
        f.write(html)


# ─────────────────────────────────────────────────────────────────────────────────────────
# BRIDGE definition -- edit these numbers to explore alternatives, then rerun.
# Defaults = exact current values from prms/pont_ferroviere_rue_goffontaine_realscale_gas.prm
# (real Pont ferroviaire, Rue Goffontaine, Trooz; see that file's header for the full
# real-value -> LU derivation of every number below).
# ─────────────────────────────────────────────────────────────────────────────────────────

LU_PER_METER = 33.0  # this prm's resolution

BRIDGE = BridgeLayout(
    domain=(2688.0, 1056.0, 432.0),
    waterline_z=115.17,  # 3.49 m flood depth (max water elev. 116.49 m - riverbed 113.0 m)
    piers=[
        # Pier width 2 m = 66 LU (radius 33); center X=1344 (deck X-span [1179,1509])
        Pier(name="pier1", x_center=1344, y_center=348.15, x_half=132, y_half=33,
             z_top=347.16, radius=33),
        Pier(name="pier2", x_center=1344, y_center=707.85, x_half=132, y_half=33,
             z_top=347.16, radius=33),
    ],
    arches=[
        # All 3 openings measured equal in the real record (8.9 m each -> radius 146.85 LU,
        # assumed semicircular -- no arch rise given in the source record).
        #
        # xMin/xMax = the piers' own straight-box span (1212..1476, i.e. pier x_center +-
        # x_half), NOT their outer reach including the rounded end-caps (1179..1509, the
        # deck's full 10 m streamwise width). Recesses each arch collar's flat leading/
        # trailing face behind the piers' rounded nose, so the nose visibly sticks out past
        # the deck edge instead of the deck being flush with the nose tip -- a deliberate
        # design change from the real record (which reports "Pier protrusion: 0", i.e.
        # flush), decided during 3D-voxel review and ported to BOTH the DEM ArchTunnels AND
        # the LBM BodyBC arch boxes in the real prm, so DEM and fluid geometry stay
        # consistent with each other (both use this recessed span, not the old flush one).
        ArchTunnel(name="opening1 (left bank)",  y_center=168.3, z_center=111.21,
                   radius=146.85, x_min=1212, x_max=1476),
        ArchTunnel(name="opening2 (centre)",     y_center=528.0, z_center=111.21,
                   radius=146.85, x_min=1212, x_max=1476),
        ArchTunnel(name="opening3 (right bank)", y_center=887.7, z_center=111.21,
                   radius=146.85, x_min=1212, x_max=1476),
    ],
)


if __name__ == "__main__":
    tiling_problems, tiling_notes = BRIDGE.check_channel_tiling()
    problems = tiling_problems + BRIDGE.check_deck_span_consistent()
    if problems:
        print("!! Shape problems found:")
        for p in problems:
            print("   -", p)
    else:
        print("Shape OK: no internal gaps/overlaps between piers and arches, "
              "deck X-span consistent across all pieces.")
    for n in tiling_notes:
        print("   (note)", n)

    BRIDGE.dimension_summary(LU_PER_METER)

    print("\n--- prm-ready text (paste into CylinderPillars/Pillars/ArchTunnels blocks) ---\n")
    print(BRIDGE.to_prm_text())

    fig = plot_layout(BRIDGE)
    out_path = "bridge_layout.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nSaved static visualization to {out_path}")

    html_path = "bridge_layout.html"
    to_interactive_html(BRIDGE, html_path, LU_PER_METER)
    print(f"Saved interactive visualization to {html_path} (open in a browser: drag to pan, "
          f"scroll to zoom, hover a pier/opening for its dimensions)")

    plt.show()
