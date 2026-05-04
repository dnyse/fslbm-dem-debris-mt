"""
Antidunes simulation — force & layer diagram (XZ side view).
Run:  python3 forces-diagram.py
Output: forces-diagram.pdf  (and .png at 200 dpi)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

matplotlib.rcParams.update({"font.family": "sans-serif", "font.size": 11})

# ── Canvas ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 8))
ax.set_xlim(-10, 430)
ax.set_ylim(-14, 128)
ax.set_aspect("equal")
ax.axis("off")

# ── Colour palette ─────────────────────────────────────────────────────────────
C = dict(
    floor="#5d4037",
    fixed_bed="#8d6e63",
    mobile_bed="#c8a97e",
    liquid="#aed6f1",
    gas="#eaf4fb",
    surface="#1a6ea8",
    pillar="#424242",
    particle="#78909c",
    fixed_particle="#37474f",
    gravity="#c0392b",
    buoyancy="#27ae60",
    drag="#e67e22",
    pid="#8e44ad",
    contact="#2c3e50",
    lubrication="#16a085",
    surface_tension="#f39c12",
    lift="#1e8bc3",
)

# ── Domain geometry ────────────────────────────────────────────────────────────
DX, DZ = 400, 120
x_bed = np.linspace(0, DX, 800)

floor_h = 10.0
fixed_h = 17.0
bed_amp  = 9.0
bed_off  = 30.0
bed_surf = bed_off + bed_amp * np.sin(2 * np.pi * x_bed / DX + np.pi)

surf_off = 70.0
surf_amp = 3.5
free_surf = surf_off + surf_amp * np.sin(2 * np.pi * x_bed / DX)

# ── Layer fills ────────────────────────────────────────────────────────────────
ax.fill_between(x_bed, free_surf, DZ,        color=C["gas"],        zorder=1)
ax.fill_between(x_bed, bed_surf,  free_surf,  color=C["liquid"],     zorder=2)
ax.fill_between(x_bed, fixed_h,   bed_surf,   color=C["mobile_bed"], zorder=3)
ax.fill_between(x_bed, floor_h,   fixed_h,    color=C["fixed_bed"],  zorder=4)
ax.fill_between(x_bed, 0,         floor_h,    color=C["floor"],      zorder=5)
ax.fill_between([0, DX], DZ - 2,  DZ,         color=C["floor"],      zorder=5)

# Hatching
ax.fill_between(x_bed, 0,       floor_h, color="none",
                hatch="////", edgecolor="#3e2723", lw=0.3, zorder=6)
ax.fill_between(x_bed, floor_h, fixed_h, color="none",
                hatch="....", edgecolor="#6d4c41", lw=0.3, zorder=6)

# ── Layer boundary lines ───────────────────────────────────────────────────────
ax.plot(x_bed, free_surf, color=C["surface"], lw=2.2, zorder=10)
ax.plot(x_bed, bed_surf,  color="#7a5230",    lw=1.0, ls="--", alpha=0.5, zorder=10)

# ── Pillar ─────────────────────────────────────────────────────────────────────
pcx, pcz, pr = 200, 38, 10
ax.add_patch(Circle((pcx, pcz), pr, color=C["pillar"], zorder=20))
# 24-HalfSpace polygon outline (visual hint)
for i in range(24):
    th = 2 * np.pi * i / 24
    ax.plot([pcx + pr * np.cos(th), pcx + pr * np.cos(th + 2*np.pi/24)],
            [pcz + pr * np.sin(th), pcz + pr * np.sin(th + 2*np.pi/24)],
            color="#bdbdbd", lw=0.6, zorder=21)
ax.text(pcx, pcz + 1, "Pillar\n(NoSlip +\n24 HalfSpaces)", ha="center", va="center",
        color="white", fontsize=7.5, fontweight="bold", zorder=22)

# ── Particle helper ────────────────────────────────────────────────────────────
p_r = 2.2

def sphere(cx, cz, color=C["particle"]):
    ax.add_patch(Circle((cx, cz), p_r, color=color, ec="black", lw=0.6, zorder=15))

# Fixed particles
for x in [40, 100, 160, 260, 330, 380]:
    sphere(x, fixed_h - 1, color=C["fixed_particle"])

# Mobile bed particles along sinusoid
bed_particles = []
for x in [50, 80, 130, 240, 290, 350, 390]:
    z = bed_off + bed_amp * np.sin(2 * np.pi * x / DX + np.pi) + p_r + 0.5
    bed_particles.append((x, z))
    sphere(x, z)

# Suspended particles in liquid
for x, z in [(110, 48), (175, 52), (310, 56)]:
    sphere(x, z)

# ── Arrow helper ───────────────────────────────────────────────────────────────
def arrow(x0, z0, dx, dz, color, label=None, lw=2.0, hw=2.0, hl=2.5,
          label_offset=None, zorder=30):
    ax.annotate(
        "", xy=(x0 + dx, z0 + dz), xytext=(x0, z0),
        arrowprops=dict(
            arrowstyle=f"->,head_width={hw},head_length={hl}",
            color=color, lw=lw),
        zorder=zorder,
    )
    if label:
        if label_offset:
            lx, lz = x0 + dx + label_offset[0], z0 + dz + label_offset[1]
        else:
            ang  = np.arctan2(dz, dx)
            perp = ang + np.pi / 2
            off  = 3.5
            lx   = x0 + dx + off * np.cos(perp)
            lz   = z0 + dz + off * np.sin(perp)
        ax.text(lx, lz, label, color=color, fontsize=9, fontweight="bold",
                ha="center", va="center", zorder=zorder + 1,
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])

# ── Forces on a focal mobile particle (x=130) ─────────────────────────────────
fx, fz = 130, bed_off + bed_amp * np.sin(2 * np.pi * 130 / DX + np.pi) + p_r + 0.5

arrow(fx, fz,  0,  -13, C["gravity"],         r"$F_g$",       label_offset=(-5, -2))
arrow(fx, fz,  0,  +13, C["buoyancy"],        r"$F_b$",       label_offset=(+5, +2))
arrow(fx, fz, +16,   0, C["drag"],            r"$F_{drag}$",  label_offset=(+1, +4))
arrow(fx, fz,  0,   +8, C["lift"],            r"$F_{lift}$",  lw=1.6, hw=1.4, hl=1.8,
      label_offset=(+7, +1))
# contact with neighbor to the left (x=80)
arrow(fx - p_r, fz, -18, 0, C["contact"],    r"$F_{contact}$", label_offset=(-2, +4))
# lubrication toward right neighbor
arrow(fx + p_r, fz, +10, 0, C["lubrication"],r"$F_{lub}$",
      lw=1.4, hw=1.2, hl=1.5, label_offset=(0, +4))
# wall reaction from floor for a fixed particle
arrow(40, fixed_h - 1 + p_r, 0, +10, C["contact"], r"$F_{wall}$", label_offset=(+6, 0))

# ── PID / body force in the liquid ────────────────────────────────────────────
arrow(280, 52, +25, 0, C["pid"], r"$F_{PID}$", lw=2.4, hw=2.4, hl=3.0,
      label_offset=(0, +5))
ax.text(292, 44, "(body force on fluid\n& all particles)", color=C["pid"],
        fontsize=7.5, ha="center", va="top",
        path_effects=[pe.withStroke(linewidth=2, foreground="white")])

# ── Gravity on fluid ──────────────────────────────────────────────────────────
arrow(350, 56, 0, -11, C["gravity"], label_offset=(+5, 0))
ax.text(357, 52, r"$F_g$", color=C["gravity"], fontsize=9, fontweight="bold",
        path_effects=[pe.withStroke(linewidth=2, foreground="white")])

# ── Surface tension at free surface ──────────────────────────────────────────
st_x = 370
st_z = surf_off + surf_amp * np.sin(2 * np.pi * st_x / DX)
arrow(st_x, st_z, -10, +4, C["surface_tension"], r"$\sigma$",
      lw=1.6, hw=1.4, hl=1.8, label_offset=(-2, +4))

# ── Flow direction ────────────────────────────────────────────────────────────
ax.annotate("", xy=(75, DZ + 4), xytext=(30, DZ + 4),
            arrowprops=dict(arrowstyle="->,head_width=2.5,head_length=3.5",
                            color="black", lw=2), zorder=35)
ax.text(52, DZ + 8, "flow direction", ha="center", fontsize=9, fontweight="bold")

# ── Axis ticks on left and bottom ─────────────────────────────────────────────
for z in [0, 20, 40, 60, 80, 100, 120]:
    ax.plot([-2, 0], [z, z], color="gray", lw=0.8)
    ax.text(-4, z, str(z), ha="right", va="center", fontsize=7.5, color="gray")
for x in [0, 100, 200, 300, 400]:
    ax.plot([x, x], [-2, 0], color="gray", lw=0.8)
    ax.text(x, -4, str(x), ha="center", va="top", fontsize=7.5, color="gray")

ax.text(-10, DZ / 2, "z [lu]", rotation=90, ha="center", va="center",
        fontsize=9, color="gray")
ax.text(DX / 2, -10, "x [lu]", ha="center", va="center",
        fontsize=9, color="gray")

# ── Layer labels (floating, right side) ───────────────────────────────────────
def layer_label(z, text, color):
    ax.text(403, z, text, ha="left", va="center", fontsize=8.5,
            color=color, fontweight="bold",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")])

layer_label(DZ - 8,                         "Gas\n(passive)",         "#607d8b")
layer_label(surf_off - 3,                   "Free surface\n(FSLBM)", C["surface"])
layer_label((fixed_h + surf_off) / 2 + 3,  "Liquid\n(LBM)",         "#1565c0")
layer_label(bed_off - 2,                    "Mobile bed\n(DEM)",     "#7a5230")
layer_label((floor_h + fixed_h) / 2,       "Fixed\nparticles",      C["fixed_particle"])
layer_label(floor_h / 2,                   "Floor\n(NoSlip)",        C["floor"])

# ── Title ─────────────────────────────────────────────────────────────────────
ax.set_title(
    "Antidunes simulation — forces and domain layers  (XZ mid-plane view)",
    fontsize=13, fontweight="bold", pad=10,
)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = [
    Line2D([0],[0], color=C["gravity"],         lw=2.5, label=r"Gravity  $F_g$  (fluid & particles)"),
    Line2D([0],[0], color=C["buoyancy"],        lw=2.5, label=r"Buoyancy  $F_b$  (particles)"),
    Line2D([0],[0], color=C["drag"],            lw=2.5, label=r"Hydrodynamic drag  $F_{drag}$  (MEM)"),
    Line2D([0],[0], color=C["lift"],            lw=2.5, label=r"Hydrodynamic lift  $F_{lift}$  (MEM)"),
    Line2D([0],[0], color=C["contact"],         lw=2.5, label=r"Contact / wall force  $F_{contact}$  (DEM)"),
    Line2D([0],[0], color=C["lubrication"],     lw=2.5, label=r"Lubrication  $F_{lub}$  (DEM)"),
    Line2D([0],[0], color=C["pid"],             lw=2.5, label=r"PID body force  $F_{PID}$  (fluid + particles)"),
    Line2D([0],[0], color=C["surface_tension"], lw=2.5, label=r"Surface tension  $\sigma$  (FSLBM)"),
    mpatches.Patch(color=C["gas"],        label="Gas phase  (passive)"),
    mpatches.Patch(color=C["liquid"],     label="Liquid phase  (LBM)"),
    mpatches.Patch(color=C["mobile_bed"], label="Mobile particle bed  (DEM)"),
    mpatches.Patch(color=C["fixed_bed"],  label="Fixed particle band"),
    mpatches.Patch(color=C["floor"],      label="Floor / top wall  (NoSlip)"),
    mpatches.Patch(color=C["pillar"],     label="Circular pillar  (NoSlip + 24 HalfSpaces)"),
]
ax.legend(handles=legend_elements, loc="upper left",
          fontsize=8, framealpha=0.95, ncol=2,
          bbox_to_anchor=(0.0, 1.0), borderaxespad=0,
          title="Force arrows & region colours", title_fontsize=8.5)

plt.tight_layout(pad=0.5)
plt.savefig("forces-diagram.pdf", bbox_inches="tight")
plt.savefig("forces-diagram.png", dpi=200, bbox_inches="tight")
print("Saved: forces-diagram.pdf  and  forces-diagram.png")
