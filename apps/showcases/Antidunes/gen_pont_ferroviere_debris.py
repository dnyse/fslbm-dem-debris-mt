#!/usr/bin/env python3
"""Generate FloatingCylinders/FloatingBoxes prm blocks for the pont_ferroviere_rue_goffontaine
debris-clogging scenarios, derived from record #30071's Accumulation/Main-debris-content fields
(pont_ferroviere_rue_goffontaine.md): total volume 60 m^3, composition natural wood 40% / building
rubble 40% / anthropogenic wood 20%, a distinguished "main trunk", carpet 10m(streamwise) x
12m(cross-stream) x 1m(tall), located on the left bank (<80% of span).

LU/m = 33 (matches pont_ferroviere_rue_goffontaine_realscale_gas.prm's own resolution derivation).
Piece sizes are NOT independently sourced from the record (which gives only aggregate volume/
composition, no individual piece dimensions) -- they are engineering judgment calls sized to (a)
land close to each category's target volume with a tractable piece count (~50-100/category, matching
the project's existing scenario scale) and (b) stay chunky enough to avoid thin-body DEM contact
degeneracy (same lesson as the ARCH-ENVELOP postmortem). Documented per-class below.
"""
import random

LU_PER_M = 33.0
M3_TO_LU3 = LU_PER_M ** 3

def cyl_volume_m3(radius_m, length_m):
    import math
    return math.pi * radius_m * radius_m * length_m

def box_volume_m3(hx_m, hy_m, hz_m):
    return (2 * hx_m) * (2 * hy_m) * (2 * hz_m)

# ---------------------------------------------------------------------------------------------
# Material classes
# ---------------------------------------------------------------------------------------------
# name, shape, geometry(m), densityRatio, count -> volume computed
classes = []

# 1. Main trunk (natural wood) -- single "key log", the record's own "Main trunk presence: Yes".
#    Sized to span/lodge across one 8.9 m opening: length 9 m, diameter 0.5 m (typical large tree).
trunk = dict(name="trunk", shape="cyl", radius_m=0.25, length_m=9.0, densityRatio=0.6, count=1)
trunk["vol_each"] = cyl_volume_m3(trunk["radius_m"], trunk["length_m"])
classes.append(trunk)

# 2. Natural wood, secondary (branches/smaller logs) -- fills the rest of the 40% natural-wood share.
nw = dict(name="nat_wood", shape="cyl", radius_m=0.18, length_m=2.2, densityRatio=0.6, count=99)
nw["vol_each"] = cyl_volume_m3(nw["radius_m"], nw["length_m"])
classes.append(nw)

# 3. Anthropogenic wood (planks/furniture/pallets -- processed, slightly denser than raw wood but
#    still buoyant). Boxes (rectilinear) to visually/physically distinguish from natural round wood.
aw = dict(name="anthro_wood", shape="box", hx_m=0.6, hy_m=0.3, hz_m=0.2, densityRatio=0.7, count=42)
aw["vol_each"] = box_volume_m3(aw["hx_m"], aw["hy_m"], aw["hz_m"])
classes.append(aw)

# 4. Building rubble (broken concrete/brick chunks) -- real density ~2200-2400 kg/m^3 -> sinks.
#    SIMPLIFICATION: the record's "carpet" is a 1m-tall SURFACE mat, but rigid-body rubble at real
#    density will physically settle toward the bed in this model (no wood-rubble entrainment/binding
#    mechanism exists in this codebase) rather than staying suspended in the surface raft. Kept at
#    real density anyway (rather than fudged toward neutral buoyancy) so the run honestly shows
#    whether sunken rubble still contributes to blockage under the wood mat, or washes downstream --
#    see the scenario prm's own header for this same note.
rb = dict(name="rubble", shape="box", hx_m=0.5, hy_m=0.4, hz_m=0.3, densityRatio=2.3, count=50)
rb["vol_each"] = box_volume_m3(rb["hx_m"], rb["hy_m"], rb["hz_m"])
classes.append(rb)

def report():
    print(f"{'class':14s} {'shape':4s} {'count':>5s} {'vol_each(m3)':>13s} {'total(m3)':>10s}")
    grand = 0.0
    nat_wood_total = 0.0
    for c in classes:
        tot = c["vol_each"] * c["count"]
        grand += tot
        if c["name"] in ("trunk", "nat_wood"):
            nat_wood_total += tot
        print(f"{c['name']:14s} {c['shape']:4s} {c['count']:5d} {c['vol_each']:13.4f} {tot:10.3f}")
    print(f"{'TOTAL':14s} {'':4s} {sum(c['count'] for c in classes):5d} {'':13s} {grand:10.3f}")
    print(f"  natural wood (trunk+nat_wood) = {nat_wood_total:.3f} m3  (target 24.0, 40%)")
    print(f"  anthropogenic wood            = {aw['vol_each']*aw['count']:.3f} m3  (target 12.0, 20%)")
    print(f"  rubble                        = {rb['vol_each']*rb['count']:.3f} m3  (target 24.0, 40%)")
    print(f"  grand total                   = {grand:.3f} m3  (target 60.0)")

# ---------------------------------------------------------------------------------------------
# Geometry constants (LU), matching pont_ferroviere_rue_goffontaine_realscale_gas.prm exactly
# ---------------------------------------------------------------------------------------------
PIER_X_MIN = 1212.0   # piers' / arch tunnels' upstream face
OPENING1_Y = 168.3     # left-bank opening centre (record: "Location at structure: Left bank")
OPENING2_Y = 528.0
OPENING3_Y = 887.7
FREE_SURFACE_Z = 115.17  # target still-water level (floorHeightFactor+liquidHeightFactor)*avgDiam

def m(v):
    return v * LU_PER_M

def flat_pieces():
    """Expand classes into a flat list of individual piece dicts (one per physical debris object)."""
    out = []
    for c in classes:
        for i in range(c["count"]):
            out.append(dict(c))
    return out

# ---------------------------------------------------------------------------------------------
# Scenario 1: RECONSTRUCTION -- pre-formed jam at the recorded footprint (opening 1, left bank),
# carpet 10m(streamwise, X) x 12m(cross-stream, Y) x 1m(tall, Z), centred so it covers the full
# opening-1 span and spills slightly onto pier 1's upstream nose (same as a real raft piling wider
# than the single gap it ultimately blocks). All pieces get distinct positions (grid + jitter) so no
# shared-spawn-point staggering is needed; a short uniform settlingSteps just lets FSLBM's fill-level
# field stabilise around each piece before release.
# ---------------------------------------------------------------------------------------------
def gen_reconstruction():
    random.seed(30071)  # record ID, for reproducibility
    pieces = flat_pieces()
    random.shuffle(pieces)

    carpet_x_max = PIER_X_MIN - m(0.3)     # small standoff so the raft doesn't spawn overlapping the pier
    carpet_x_min = carpet_x_max - m(10.0)  # 10 m streamwise
    carpet_y_min = 0.0
    carpet_y_max = m(12.0)                 # 12 m cross-stream, flush against the left abutment
    z_center = FREE_SURFACE_Z

    lines_cyl, lines_box = [], []
    for idx, p in enumerate(pieces):
        x = random.uniform(carpet_x_min, carpet_x_max)
        y = random.uniform(carpet_y_min, carpet_y_max)
        # rubble rides lower (denser, settles toward the base of the mat / bed); wood rides at the
        # surface +-0.5m within the 1m carpet height.
        if p["name"] == "rubble":
            z = z_center - m(0.6) + random.uniform(-m(0.2), m(0.2))
        else:
            z = z_center + random.uniform(-m(0.4), m(0.4))
        settling = 2000 + (idx % 6) * 1000  # 2000..7000, just de-synchronises the "materialisation"

        if p["shape"] == "cyl":
            axis = random.choice(["< 1, 0, 0 >", "< 0, 1, 0 >", "< 0.7, 0.7, 0 >"])
            lines_cyl.append(
                f"   FloatingCylinder {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; "
                f"radius {m(p['radius_m']):.2f}; length {m(p['length_m']):.2f}; axis {axis}; "
                f"densityRatio {p['densityRatio']}; settlingSteps {settling}; }}"
            )
        else:
            lines_box.append(
                f"   FloatingBox {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; "
                f"halfExtents < {m(p['hx_m']):.2f}, {m(p['hy_m']):.2f}, {m(p['hz_m']):.2f} >; "
                f"densityRatio {p['densityRatio']}; settlingSteps {settling}; }}"
            )
    return lines_cyl, lines_box

# ---------------------------------------------------------------------------------------------
# Scenario 2/3: GENESIS -- spawned upstream, drift downstream under the flow, and are free to jam
# wherever the hydraulics/geometry actually push them (NOT pre-biased toward opening 1) -- tests
# whether the recorded left-bank concentration is an emergent hydraulic result or needs seeding.
# 3 aim points (one per opening's own Y-centre), each releasing its whole share of pieces in staged
# waves (startInactive + staggered settlingSteps, same convention as every prior bridge scenario in
# this project). include_trunk=False produces the "no main trunk" sensitivity control (scenario 3).
# ---------------------------------------------------------------------------------------------
def gen_genesis(include_trunk: bool):
    random.seed(20210714)  # July 2021 flood date, for reproducibility
    pieces = [p for p in flat_pieces() if include_trunk or p["name"] != "trunk"]
    random.shuffle(pieces)

    aim_x = PIER_X_MIN - m(9.5)  # ~313 LU upstream of the pier face -> ~15,600-step travel at
                                  # targetMeanVelocityMagnitude=0.02 LU/step
    aim_points = [OPENING1_Y, OPENING2_Y, OPENING3_Y]
    wave_gap = 400  # steps between successive releases at the same aim point
    first_release = 2000

    lines_cyl, lines_box = [], []
    for i, p in enumerate(pieces):
        aim_y = aim_points[i % 3]
        wave = i // 3
        settling = first_release + wave * wave_gap
        x, y = aim_x, aim_y
        z = FREE_SURFACE_Z - (m(0.6) if p["name"] == "rubble" else 0.0)

        if p["shape"] == "cyl":
            axis = "< 0, 1, 0 >"  # spawned broadside to the flow, matching prior scenarios' convention
            lines_cyl.append(
                f"   FloatingCylinder {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; "
                f"radius {m(p['radius_m']):.2f}; length {m(p['length_m']):.2f}; axis {axis}; "
                f"densityRatio {p['densityRatio']}; settlingSteps {settling}; startInactive true; }}"
            )
        else:
            lines_box.append(
                f"   FloatingBox {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; "
                f"halfExtents < {m(p['hx_m']):.2f}, {m(p['hy_m']):.2f}, {m(p['hz_m']):.2f} >; "
                f"densityRatio {p['densityRatio']}; settlingSteps {settling}; startInactive true; }}"
            )
    last_settling = first_release + (len(pieces) // 3) * wave_gap
    return lines_cyl, lines_box, last_settling

if __name__ == "__main__":
    report()
    print()

    rc_cyl, rc_box = gen_reconstruction()
    with open("/tmp/pf_reconstruction_cyl.txt", "w") as f: f.write("\n".join(rc_cyl) + "\n")
    with open("/tmp/pf_reconstruction_box.txt", "w") as f: f.write("\n".join(rc_box) + "\n")
    print(f"reconstruction: {len(rc_cyl)} cylinders, {len(rc_box)} boxes")

    ge_cyl, ge_box, last_t = gen_genesis(include_trunk=True)
    with open("/tmp/pf_genesis_cyl.txt", "w") as f: f.write("\n".join(ge_cyl) + "\n")
    with open("/tmp/pf_genesis_box.txt", "w") as f: f.write("\n".join(ge_box) + "\n")
    print(f"genesis: {len(ge_cyl)} cylinders, {len(ge_box)} boxes, last release wave at t={last_t}")

    nt_cyl, nt_box, nt_last_t = gen_genesis(include_trunk=False)
    with open("/tmp/pf_genesis_notrunk_cyl.txt", "w") as f: f.write("\n".join(nt_cyl) + "\n")
    with open("/tmp/pf_genesis_notrunk_box.txt", "w") as f: f.write("\n".join(nt_box) + "\n")
    print(f"genesis_notrunk: {len(nt_cyl)} cylinders, {len(nt_box)} boxes, last release wave at t={nt_last_t}")
