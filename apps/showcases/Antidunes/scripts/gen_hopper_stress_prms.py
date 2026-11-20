#!/usr/bin/env python3
"""Variants of the 406-rod capped-PID hopper, written by text substitution on the base prm
BASE = prms/hopper_clog_cyl_frontblock_pidon_2node_skin0p5_curvguard_cappid_gas.prm (block pressed against the pillars).

Changes common to every variant (common()):
 - settlingSteps 2000 -> 500 (the floaters are held fixed for 500 steps);
 - VTK snapshot list gets step 0 (written before the first step) and 120000, 149000;
 - wallOuterGhostLayer true (the outer ghost layer becomes wall too; removes the side-wall NaN normals);
 - the four ClogMeasurement probe points lowered from z = 40 to z = 20 (pressure cell stays submerged downstream);
 - timesteps 100000 -> 150000;
 - ClogMeasurement.blockageRatioEnabled true (blockage ratio at the probe planes, Evaluators.h blockageAt).

Without options, two prms are written:
1) CUBOIDS -> OUT_BOX
   Every FloatingCylinder (R 3, L 18, axis y) is replaced 1:1 by a FloatingBox with the same bounding footprint,
   halfExtents <3, 9, 3> (6 x 18 x 6, long axis y), same position and density (the convention of
   gen_hopper_clog_box_prm.py), plus boxContactSkin 0.5 (box analog of cylinderContactSkin 0.5). Box x Box contacts
   use boxBoxSAT with the closed-form contact point (boxBoxContactPoint, ClosedFormContactsBoxPoint.h). Clearances with
   the skin: box-box 2.0 (x) / 1.52 (z) / 2.0 (y) >= 2 m = 1, box-wall 1.0 and box-floor 1.0 >= m, front row to the
   pillars 1.04 >= m (asserted in make_boxes).
2) CYLINDERS, BLOCK AT THE START -> OUT_CYL ("drifting rods")
   The whole 406-rod block is shifted in x so that its first row sits at x = 8 (centres 144..368 -> 8..232; the
   channel is periodic in x). The front row is then 140 LU upstream of the pillar faces (x = 372) and the upstream
   probe pair (x = 336) starts in clean water. Most --*-only variants below start from this layout.
Each --<variant>-only option writes one variant instead (see main()).
Usage: python3 gen_hopper_stress_prms.py [--<variant>-only]
"""
import math
import re

BASE = "prms/hopper_clog_cyl_frontblock_pidon_2node_skin0p5_curvguard_cappid_gas.prm"
OUT_BOX = "prms/hopper_clog_box_frontblock_skin0p5_cappid_sat_gas.prm"
OUT_CYL = "prms/hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm"
OUT_MIX = "prms/hopper_clog_mixed_startblock_skin0p5_cappid_gas.prm"
OUT_BOXDRIFT = "prms/hopper_clog_box_startblock_skin0p5_cappid_gas.prm"
OUT_SPHDRIFT = "prms/hopper_clog_sph_startblock_skin0p5_cappid_gas.prm"
OUT_CYLNOSKIN = "prms/hopper_clog_cyl_startblock_noskin_cappid_gas.prm"
OUT_CYLCAP = "prms/hopper_clog_cyl_startblock_skin0p5_cap2p5_gas.prm"
OUT_CYLCAP15 = "prms/hopper_clog_cyl_startblock_skin0p5_cap1p5_gas.prm"
OUT_GIANT = "prms/hopper_clog_giantrods100_startblock_cappid_gas.prm"
OUT_GIANT_Z36 = "prms/hopper_clog_giantrods100_startblock_cappid_z36_gas.prm"
OUT_TWOPIL = "prms/twopillar_clog_cyl812_startblock_cappid_gas.prm"
OUT_TWOPIL_CLEAN = "prms/twopillar_clean_reference_gas.prm"
OUT_TWOPIL_CLEAN_LONG = "prms/twopillar_clean_reference_150k_gas.prm"
OUT_HOPPER_CLEAN_LONG = "prms/hopper_clean_reference_150k_gas.prm"
OUT_CENTRED_SB = "prms/clog_centred_pillar_startblock_cappid_gas.prm"
OUT_CYLREP2 = "prms/hopper_clog_cyl_startblock_skin0p5_cappid_rep2_gas.prm"
OUT_CYLMASS = "prms/hopper_clog_cyl_startblock_skin0p5_cappid_drivebymass_gas.prm"
OUT_CYLSKIN = "prms/hopper_clog_cyl_startblock_skin{tag}_cappid_gas.prm"
OUT_POCKET = "prms/skin_pocket_test_m{tag}_gas.prm"
OUT_SIGMA = "prms/hopper_clog_cyl_startblock_skin0p5_cappid_sigma{tag}_gas.prm"
OUT_SPH4 = "prms/hopper_clog_sph4_startblock_cappid_gas.prm"
SPH4_STEPS = ("0,500,1000,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,22000,24000,28000,32000,36000,40000,"
              "48000,56000,64000,80000,99000,120000,149000")
BASE_CENTRED ="prms/clog_measurement_validation_extreme_cappid_2node_flowrate_gas.prm"
OUT_CENTRED = "prms/clog_centred_pillar_cappid_current_gas.prm"
CENTRED_STEPS = "0,500,4000,8000,16000,20000,24000,28000,32000,36000,40000,44000,48000,52000,56000,64000,80000,99000,120000,149000"
SPHERE_YS = (4.0, 12.0, 20.0, 28.0, 36.0)   # 5 spheres R 3 across the channel, as in the 1015-sphere hopper prms
SKIN = 0.5
X_START = 8.0
CYL_RE = re.compile(r"FloatingCylinder \{ position < ([\d.]+), ([\d.]+), ([\d.]+) >; radius 3; length 18; axis < 0, 1, 0 >; "
                    r"densityRatio ([\d.]+); settlingSteps 2000; \}")


def common(s, n_floaters=406):
    s, n1 = re.subn(r"settlingSteps 2000;", "settlingSteps 500;", s)
    s, n2 = re.subn(r'steps "500,1000,([0-9,]*,99000)";[^\n]*',
                    r'steps "0,500,1000,\1,120000,149000";  // denser around the jam', s)
    s, n5 = re.subn(r"(\n\s*timesteps\s+)100000;[^\n]*",
                    r"\g<1>150000;  // base prm: 100000", s)
    s, n3 = re.subn(r"(\n(\s*)useMemOnlyBuoyancy\s+true;[^\n]*)",
                    r"\1\n\2wallOuterGhostLayer           true;   // outer ghost layer becomes wall too (side-wall NaN-normal fix)", s)
    # pressure probes lowered from z = 40 to z = 20: the probe pressure is read at the single cell of the probe position
    # (Evaluators.h ClogMeasurementEvaluator), and behind a jam the downstream level can fall to ~26-32, leaving a z = 40
    # cell in air. Elevation and discharge are column/plane sums and do not depend on z.
    s, n4 = re.subn(r"((?:upstreamPos|downstreamPos|pillarUpstreamPos|pillarDownstreamPos)\s+< *\d+, *\d+, *)40( *>);", r"\g<1>20\2;  // z = 20: pressure cell stays submerged downstream", s)
    s, n6 = re.subn(r"\n(\s*)(flowRateEnabled\s+true;[^\n]*)",
                    r"\n\1\2\n\1blockageRatioEnabled true;   // pillar + jam blockage at the probe planes (clog*BlockageRatio columns)", s)
    assert n1 == n_floaters and n2 == 1 and n3 == 1 and n4 == 4 and n5 == 1 and n6 == 1, (n1, n2, n3, n4, n5, n6)
    return s


def rods(s):
    r = [tuple(float(v) for v in m.groups()) for m in CYL_RE.finditer(s)]
    assert len(r) == 406, len(r)
    return r


def make_boxes(s):
    r = rods(s)
    # clearances with the skin (box half-extents 3, 9, 3; channel y 0..40, floor z 0; pillars R 12 at (384, 0) and (384, 40))
    xs, ys, zs = sorted({p[0] for p in r}), sorted({p[1] for p in r}), sorted({p[2] for p in r})
    gx = min(b - a for a, b in zip(xs, xs[1:])) - 6
    gz = min(b - a for a, b in zip(zs, zs[1:])) - 6
    gy = (ys[1] - 9) - (ys[0] + 9)
    gwall = min(ys[0] - 9, 40 - (ys[-1] + 9))
    gfloor = zs[0] - 3
    xf = xs[-1] + 3
    gpil = min(math.hypot(384 - xf, max(0.0, min(abs(py - c) for py in (ys[0] - 9, ys[0] + 9, ys[1] - 9, ys[1] + 9)))) - 12
               for c in (0, 40))
    assert min(gx, gz, gy) >= 2 * SKIN and min(gwall, gfloor, gpil) >= SKIN, (gx, gz, gy, gwall, gfloor, gpil)
    print(f"  box clearances: x {gx:.2f}, z {gz:.2f}, y {gy:.2f}, wall {gwall:.2f}, floor {gfloor:.2f}, pillar {gpil:.2f}")
    s = CYL_RE.sub(lambda m: f"FloatingBox {{ position < {m.group(1)}, {m.group(2)}, {m.group(3)} >; halfExtents < 3, 9, 3 >; "
                             f"densityRatio {m.group(4)}; settlingSteps 2000; }}", s)
    s, n = re.subn(r"\nFloatingCylinders\n", "\nFloatingBoxes\n", s)
    assert n == 1
    s, n = re.subn(r"(\n(\s*)cylinderContactSkin\s+0\.5;[^\n]*)",
                   r"\1\n\2boxContactSkin                0.5;   // box analog of cylinderContactSkin (each half-extent + 0.5)", s)
    assert n == 1
    return s


def make_startblock(s):
    r = rods(s)
    dx = X_START - min(p[0] for p in r)
    s = CYL_RE.sub(lambda m: m.group(0).replace(f"< {m.group(1)},", f"< {float(m.group(1)) + dx:.2f},", 1), s)
    r2 = rods(s)
    print(f"  shift dx = {dx:+.2f}: block x {min(p[0] for p in r2):.2f}..{max(p[0] for p in r2):.2f}")
    return s


def make_mixed(s):
    """Mixed floaters: the start-block layout (x 8..232, 29 rows x 7 layers, densities 1.2 / 1.0 / 0.5 by layer as in
    the base), rows cycling cylinder -> box -> sphere from the front row: a cylinder row keeps its 2 rods (R 3, L 18, axis y),
    a box row has 2 boxes 6 x 18 x 6 in the same places, a sphere row has 5 spheres R 3 at y = 4, 12, ..., 36 (the grid of
    the 1015-sphere hopper prms). All shapes share the 6 LU cross-section. Skins: cylinderContactSkin 0.5 (already in the
    base) and boxContactSkin 0.5; mixed pairs (box-cylinder, box-sphere, cylinder-sphere) have no skin in the code."""
    s = make_startblock(s)
    r = rods(s)
    xs = sorted({p[0] for p in r}, reverse=True)          # front row first
    kind_of_x = {x: ("cyl", "box", "sph")[i % 3] for i, x in enumerate(xs)}
    cyl, box, sph, aabbs = [], [], [], []
    for (x, y, z, d) in r:
        k = kind_of_x[x]
        if k == "cyl":
            cyl.append(f"   FloatingCylinder {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; radius 3; length 18; axis < 0, 1, 0 >; "
                       f"densityRatio {d:g}; settlingSteps 2000; }}")
            aabbs.append(((x - 3, y - 9, z - 3), (x + 3, y + 9, z + 3), "cyl"))
        elif k == "box":
            box.append(f"   FloatingBox {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; halfExtents < 3, 9, 3 >; densityRatio {d:g}; "
                       f"settlingSteps 2000; }}")
            aabbs.append(((x - 3, y - 9, z - 3), (x + 3, y + 9, z + 3), "box"))
        elif y == min(p[1] for p in r):                    # one sphere row per (x, z): emit the 5 spheres once
            for sy in SPHERE_YS:
                sph.append(f"   FloatingSphere {{ position < {x:.2f}, {sy:.2f}, {z:.2f} >; radius 3; densityRatio {d:g}; "
                           f"settlingSteps 2000; }}")
                aabbs.append(((x - 3, sy - 3, z - 3), (x + 3, sy + 3, z + 3), "sph"))
    # clearances on the axis-aligned bounding boxes (conservative: the true surface gap is never smaller)
    def gap(a, b):
        return math.sqrt(sum(max(0.0, max(a[0][i] - b[1][i], b[0][i] - a[1][i])) ** 2 for i in range(3)))
    gmin = min(gap(a, b) for i, a in enumerate(aabbs) for b in aabbs[i + 1:]
               if abs(a[0][0] - b[0][0]) < 20 and abs(a[0][2] - b[0][2]) < 20)
    gwall = min(min(a[0][1], 40 - a[1][1]) for a in aabbs)
    gfloor = min(a[0][2] for a in aabbs)
    assert gmin >= 2 * SKIN and gwall >= SKIN and gfloor >= SKIN, (gmin, gwall, gfloor)
    print(f"  mixed: {len(cyl)} cylinders, {len(box)} boxes, {len(sph)} spheres; min AABB gap {gmin:.2f}, wall {gwall:.2f}, "
          f"floor {gfloor:.2f}")
    a0 = s.index("\nFloatingCylinders\n")
    a1 = s.index("\n}", a0) + 2
    blocks = ("\nFloatingCylinders\n{\n" + "\n".join(cyl) + "\n}\n\nFloatingBoxes\n{\n" + "\n".join(box)
              + "\n}\n\nFloatingSpheres\n{\n" + "\n".join(sph) + "\n}\n")
    s = s[:a0] + blocks + s[a1:]
    s, n = re.subn(r"(\n(\s*)cylinderContactSkin\s+0\.5;[^\n]*)",
                   r"\1\n\2boxContactSkin                0.5;   // box analog of cylinderContactSkin (each half-extent + 0.5)", s)
    assert n == 1
    return s, len(cyl) + len(box) + len(sph)


def make_spheres(s):
    """Every rod position of the (shifted) block becomes a row of 5 spheres R 3 across the channel
    at y = 4, 12, ..., 36 (the grid of the 1015-sphere hopper prms): 29 x 7 x 5 = 1015 spheres, densities by layer as the
    rods. Spheres have no contact skin (sphere-sphere pairs get the lubrication correction instead)."""
    r = rods(s)
    ymin = min(p[1] for p in r)
    sph = [f"   FloatingSphere {{ position < {x:.2f}, {sy:.2f}, {z:.2f} >; radius 3; densityRatio {d:g}; settlingSteps 2000; }}"
           for (x, y, z, d) in r if y == ymin for sy in SPHERE_YS]
    xs = sorted({p[0] for p in r}); zs = sorted({p[2] for p in r})
    gx = min(b - a for a, b in zip(xs, xs[1:])) - 6; gz = min(b - a for a, b in zip(zs, zs[1:])) - 6
    gy = min(b - a for a, b in zip(SPHERE_YS, SPHERE_YS[1:])) - 6
    gw = min(SPHERE_YS[0] - 3, 40 - SPHERE_YS[-1] - 3); gf = zs[0] - 3
    assert min(gx, gz, gy) >= 1.0 and min(gw, gf) >= 0.5, (gx, gz, gy, gw, gf)
    print(f"  spheres: {len(sph)}; gaps x {gx:.2f}, z {gz:.2f}, y {gy:.2f}, wall {gw:.2f}, floor {gf:.2f}")
    a0 = s.index("\nFloatingCylinders\n")
    a1 = s.index("\n}", a0) + 2
    return s[:a0] + "\nFloatingSpheres\n{\n" + "\n".join(sph) + "\n}\n" + s[a1:], len(sph)


def make_spheres_r4(s):
    """Spheres that can bridge: R 4 instead of R 3, so the
    16 LU gap is 2 sphere diameters (R 3: 2.7). Same block extent as the rods (x 5..235 -> sphere rows x 9..225, pitch 9),
    4 spheres across (y 5.5, 15, 24.5, 34), 6 layers z 5..50 (pitch 9), densities by layer 1.2 / 1.2 / 1.0 / 1.0 / 0.5 /
    0.5 (rods: 2 / 3 / 2 layers of 7). 25 x 4 x 6 = 600 spheres, solid volume 1.61e5 (R 3 spheres 1.15e5, rods 2.07e5).
    Surface gaps >= 1 LU between spheres, >= 1.5 to the walls; spheres have no contact skin."""
    xs = [9.0 + 9.0 * i for i in range(25)]
    ys = (5.5, 15.0, 24.5, 34.0)
    zd = ((5.0, 1.2), (14.0, 1.2), (23.0, 1.0), (32.0, 1.0), (41.0, 0.5), (50.0, 0.5))
    sph = [f"   FloatingSphere {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; radius 4; densityRatio {d:g}; settlingSteps 2000; }}"
           for x in xs for (z, d) in zd for y in ys]
    gx = xs[1] - xs[0] - 8; gy = min(b - a for a, b in zip(ys, ys[1:])) - 8; gz = zd[1][0] - zd[0][0] - 8
    gw = min(ys[0] - 4, 40 - ys[-1] - 4); gf = zd[0][0] - 4
    assert min(gx, gy, gz) >= 1.0 and min(gw, gf) >= 0.5 and xs[-1] + 4 < 372, (gx, gy, gz, gw, gf)
    print(f"  R4 spheres: {len(sph)}; gaps x {gx:.2f}, y {gy:.2f}, z {gz:.2f}, wall {gw:.2f}, floor {gf:.2f}, "
          f"volume {len(sph) * 4 / 3 * math.pi * 64:.0f}")
    a0 = s.index("\nFloatingCylinders\n")
    a1 = s.index("\n}", a0) + 2
    return s[:a0] + "\nFloatingSpheres\n{\n" + "\n".join(sph) + "\n}\n" + s[a1:], len(sph)


def jitter_rods(s, amp=0.25, seed=1):
    """Second realisation of the drifting rods: every rod moved by a seeded uniform random offset in
    [-amp, amp] in x and z (y unchanged). The rod grid has surface gaps of 2.0 (x) and 1.52 (z), so the smallest gap stays
    >= 1.02 > 2 x skin; floor clearance 1.0 - amp >= skin."""
    import random
    rng = random.Random(seed)

    def mv(m):
        x, z = float(m.group(1)), float(m.group(3))
        return m.group(0).replace(f"< {m.group(1)}, {m.group(2)}, {m.group(3)} >",
                                  f"< {x + rng.uniform(-amp, amp):.3f}, {m.group(2)}, {z + rng.uniform(-amp, amp):.3f} >", 1)
    s2 = re.sub(r"FloatingCylinder \{ position < ([\d.]+), ([\d.]+), ([\d.]+) >;[^}]*\}", mv, s)
    r = [tuple(float(v) for v in m.groups()) for m in re.finditer(
        r"FloatingCylinder \{ position < ([\d.]+), ([\d.]+), ([\d.]+) >; radius 3; length 18; axis < 0, 1, 0 >; "
        r"densityRatio ([\d.]+); settlingSteps 2000; \}", s2)]
    assert len(r) == 406
    gmin = min(math.hypot(a[0] - b[0], a[2] - b[2]) - 6 for i, a in enumerate(r) for b in r[i + 1:]
               if a[1] == b[1] and abs(a[0] - b[0]) < 12 and abs(a[2] - b[2]) < 12)
    gfloor = min(p[2] for p in r) - 3
    assert gmin >= 2 * SKIN and gfloor >= SKIN, (gmin, gfloor)
    print(f"  jitter +-{amp} (seed {seed}): min rod-rod gap {gmin:.3f}, floor {gfloor:.3f}, x {min(p[0] for p in r):.2f}.."
          f"{max(p[0] for p in r):.2f}")
    return s2


def make_giant_rods(s):
    """100 giant rods drifting from upstream: R 6, L 36, axis y at y = 20 (rod spans y 2..38, wall gap 2),
    4 layers z 7.5 / 21 / 34.5 / 48 (pitch 13.5, gap 1.5), densities 1.2 / 1.0 / 1.0 / 0.5 by layer, 25 rows at pitch 13.5
    with the front row at x = 222.75 (front face 228.75, ~143 LU upstream of the pillar faces at 372, as the 406-rod start
    block). The block (336 LU long) does not fit between x = 0 and the front, so the rear rows wrap through the periodic
    boundary: rows at x = 222.75 - 13.5 k mod 768, k = 0..24 (rear row x = 666.75, 265 LU behind the pillars); the row
    positions are chosen so that no rod straddles x = 0 (rows at 6.75 and 761.25). NB for the analysis: rods starting at
    x > 396 are UPSTREAM of the pillars through the periodic boundary (unwrap their start x by -768)."""
    xs = [(222.75 - 13.5 * k) % 768 for k in range(25)]
    zd = ((7.5, 1.2), (21.0, 1.0), (34.5, 1.0), (48.0, 0.5))
    assert all(min(x - 6, 768 - (x + 6)) >= 0.5 for x in xs), xs
    rods = [f"   FloatingCylinder {{ position < {x:.2f}, 20.00, {z:.2f} >; radius 6; length 36; axis < 0, 1, 0 >; "
            f"densityRatio {d:g}; settlingSteps 2000; }}" for x in xs for (z, d) in zd]
    gx = 13.5 - 12; gz = zd[1][0] - zd[0][0] - 12; gw = 20 - 18; gf = zd[0][0] - 6
    assert min(gx, gz) > 2 * SKIN and min(gw, gf) >= SKIN, (gx, gz, gw, gf)
    print(f"  giant rods: {len(rods)}; gaps x {gx:.2f}, z {gz:.2f}, wall {gw:.2f}, floor {gf:.2f}; rows x "
          f"{sorted(xs)[0]:.2f}..{sorted(xs)[-1]:.2f}; volume {len(rods) * math.pi * 36 * 36:.0f}")
    a0 = s.index("\nFloatingCylinders\n")
    a1 = s.index("\n}", a0) + 2
    return s[:a0] + "\nFloatingCylinders\n{\n" + "\n".join(rods) + "\n}\n" + s[a1:], len(rods)


def make_two_pillars(s, with_rods=True):
    """Two-pillar layout in which a rod fits through all three slots (pillars R 5): channel width 40 -> 80
    (domainSize y, cellsPerBlock unchanged -> 12 x 4 x 6 = 288 blocks = 288 ranks), the two half-pillars R 12 at the walls
    replaced by two full pillars R 5 at (384, 25) and (384, 55): three slots of 20 LU (y 0..20, 30..50, 60..80), each
    longer than a rod (L 18; 19 with the skin). The 406-rod start block is mirrored to y + 40 (rods at y 10 / 30 / 50 / 70,
    812 rods, wall gap 1, rod-rod gap in y 2). Probes: centre pair in the middle slot (y = 40), pillar pair in line with
    the first pillar (y = 25). VTK cell filters widened to y 80. with_rods=False: clean reference (no floaters)."""
    s, n1 = re.subn(r"domainSize\s+<768, 40, 180>;", "domainSize         <768, 80, 180>;  // two pillars: width 80", s)
    s, n2 = re.subn(r"\n(\s*)CylinderPillar \{ center < 384, 0, 45 >; radius 12; \}[^\n]*\n\s*CylinderPillar \{ center < 384, 40, 45 >; radius 12; \}[^\n]*",
                    r"\n\1CylinderPillar { center < 384, 25, 45 >; radius 5; }  // two pillars R 5: three 20 LU slots"
                    r"\n\1CylinderPillar { center < 384, 55, 45 >; radius 5; }", s)
    s, n3 = re.subn(r"\n(\s*)BodyBC \{ shape Cylinder; min < 384, 0, 0 >;  max < 384, 0, 180 >;  radius 12; NoSlip\{\} \}[^\n]*\n\s*BodyBC \{ shape Cylinder; min < 384, 40, 0 >; max < 384, 40, 180 >; radius 12; NoSlip\{\} \}[^\n]*",
                    r"\n\1BodyBC { shape Cylinder; min < 384, 25, 0 >;  max < 384, 25, 180 >;  radius 5; NoSlip{} }  // two pillars R 5"
                    r"\n\1BodyBC { shape Cylinder; min < 384, 55, 0 >;  max < 384, 55, 180 >;  radius 5; NoSlip{} }", s)
    s, n4 = re.subn(r"((?:upstreamPos|downstreamPos)\s+< *\d+, *)20(, *\d+ *>;)", r"\g<1>40\2  // middle slot", s)
    s, n5 = re.subn(r"((?:pillarUpstreamPos|pillarDownstreamPos)\s+< *\d+, *)6(, *\d+ *>;)", r"\g<1>25\2  // in line with pillar 1", s)
    s, n6 = re.subn(r"(CellBB_filter_0 \{ min < 0, 0, 0 >;  max < 768, )40(, 180 >; \})", r"\g<1>80\2", s)
    assert (n1, n2, n3, n4, n5, n6) == (1, 1, 1, 2, 2, 2), (n1, n2, n3, n4, n5, n6)
    a0 = s.index("\nFloatingCylinders\n")
    a1 = s.index("\n}", a0) + 2
    if not with_rods:
        print("  two pillars: clean reference, no floaters")
        return s[:a0] + "\n" + s[a1:], 0
    r = rods(s)
    lines = []
    for (x, y, z, d) in r:
        for yy in (y, y + 40):
            lines.append(f"   FloatingCylinder {{ position < {x:.2f}, {yy:.2f}, {z:.2f} >; radius 3; length 18; axis < 0, 1, 0 >; "
                         f"densityRatio {d:g}; settlingSteps 2000; }}")
    ys = sorted({y for (_, y, _, _) in r} | {y + 40 for (_, y, _, _) in r})
    assert ys == [10.0, 30.0, 50.0, 70.0] and min(ys[0] - 9, 80 - (ys[-1] + 9)) >= SKIN, ys
    print(f"  two pillars: {len(lines)} rods at y {ys}, slots 20 LU (pillars R 5 at y 25 / 55)")
    return s[:a0] + "\nFloatingCylinders\n{\n" + "\n".join(lines) + "\n}\n" + s[a1:], len(lines)


def make_centred_startblock(s):
    """Centred pillar with the 100 spheres R 5 as a start block: the staggered release of OUT_CENTRED (x 40..120,
    z 73..75, released one by one between steps 20000 and 56000) replaced by a block at the start of the domain, released together at step 500 like
    the hopper start blocks. 11 rows at x = 13 + 11.5 k (k = 0..10) x 3 across (y 6.5 / 20 / 33.5) x 3 layers by density
    (z 6: 1.2 -> here 1.5, z 18: 1.0, z 40: 0.5), plus one light sphere at (139.5, 20, 40): 34 x 0.5, 33 x 1.0, 33 x 1.5 as before.
    Gaps: x 1.5, y 3.5, z >= 2, walls 1.5, floor 1 (spheres have no skin). VTK at the hopper snapshot steps."""
    rows = [13.0 + 11.5 * k for k in range(11)]
    ys = (6.5, 20.0, 33.5)
    zd = ((6.0, 1.5), (18.0, 1.0), (40.0, 0.5))
    sph = [f"   FloatingSphere {{ position < {x:.2f}, {y:.2f}, {z:.2f} >; radius 5; densityRatio {d:g}; settlingSteps 500; }}"
           for x in rows for (z, d) in zd for y in ys]
    sph.append("   FloatingSphere { position < 139.50, 20.00, 40.00 >; radius 5; densityRatio 0.5; settlingSteps 500; }")
    gx = 11.5 - 10; gy = 13.5 - 10; gz = min(b[0] - a[0] for a, b in zip(zd, zd[1:])) - 10
    assert len(sph) == 100 and min(gx, gy, gz) >= 1.0 and ys[0] - 5 >= 0.5 and 40 - ys[-1] - 5 >= 0.5 and zd[0][0] - 5 >= 0.5
    assert 139.5 + 5 < 384 - 12
    print(f"  centred start block: {len(sph)} spheres R 5, gaps x {gx}, y {gy}, z {gz}, walls 1.5, floor 1, x 8..144.5")
    a0 = s.index("\nFloatingSpheres\n")
    a1 = s.index("\n}", a0) + 2
    s = s[:a0] + "\nFloatingSpheres\n{\n" + "\n".join(sph) + "\n}\n" + s[a1:]
    s, k = re.subn(r'steps "0,500,4000,[0-9,]*,149000";[^\n]*',
                   'steps "0,500,1000,2000,4000,8000,12000,16000,20000,24000,32000,48000,64000,80000,99000,120000,149000";  // hopper snapshot steps (start block)', s)
    assert k == 1 and "startInactive" not in s.split("FloatingSpheres",1)[1].split("}\n}",1)[0]
    return s


def make_centred(s):
    """Centred-pillar layout (BASE_CENTRED: 100 spheres, staggered release 20000..56000, cap 2.65e-6) with the
    settings of the current hopper runs: 150000 steps, wallOuterGhostLayer, probes at z = 20, blockage ratio, and VTK
    at fixed snapshot steps (fluid_field periodic output off, fluid_field_snapshots added; floater/pillar VTK follow it)."""
    s, n1 = re.subn(r"(\n\s*timesteps\s+)124750;", r"\g<1>150000;  // base prm: 124750", s)
    s, n2 = re.subn(r"(\n(\s*)useMemOnlyBuoyancy\s+true;[^\n]*)",
                    r"\1\n\2wallOuterGhostLayer           true;   // outer ghost layer becomes wall too (side-wall NaN-normal fix)", s)
    s, n3 = re.subn(r"((?:upstreamPos|downstreamPos|pillarUpstreamPos|pillarDownstreamPos)\s+< *\d+, *\d+, *)40( *>);",
                    r"\g<1>20\2;  // z = 20: pressure cell stays submerged downstream", s)
    s, n4 = re.subn(r"\n(\s*)(flowRateEnabled\s+true;[^\n]*)",
                    r"\n\1\2\n\1blockageRatioEnabled true;   // pillar + jam blockage at the probe planes", s)
    s, n5 = re.subn(r"(\n\s*fluid_field\s*\n\s*\{\s*\n\s*writeFrequency\s+)24950;",
                    r"\g<1>0;   // periodic output off; see fluid_field_snapshots", s)
    snap = ("\n   // fixed snapshot steps, denser during the staggered release (20000..56000)\n"
            "   fluid_field_snapshots\n   {\n      writeFrequency       1;\n      ghostLayers          0;\n"
            "      baseFolder           vtk-out;\n      samplingResolution   1;\n      useMPIIO false;\n"
            f'      steps "{CENTRED_STEPS}";\n'
            "      writers\n      {\n         velocity;\n         mapped_flag;\n      }\n"
            "      before_functions { gas_cell_zero_setter; }\n   }\n")
    s, n6 = re.subn(r"(\n\s*domain_decomposition\s*\n)", lambda m: snap + m.group(1), s, count=1)
    assert (n1, n2, n3, n4, n5, n6) == (1, 1, 4, 1, 1, 1), (n1, n2, n3, n4, n5, n6)
    return s


def main():
    import sys
    base = open(BASE).read()
    if "--centred-only" in sys.argv:
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--centred-only) -- CENTRED PILLAR with the current hopper settings:\n"
               f"// {BASE_CENTRED} + 150000 steps, wallOuterGhostLayer, probes z = 20, blockage ratio,\n"
               "// VTK snapshots from step 0. Same cap (2.65e-6), so the same driving force as the hopper runs.\n")
        print(OUT_CENTRED)
        open(OUT_CENTRED, "w").write(hdr + make_centred(open(BASE_CENTRED).read()))
        return
    if "--cyl-startblock-cap2p5-only" in sys.argv:
        # force-level control: the drifting rods with the cap halved, 2.65e-6 -> 1.325e-6 = 2.5 a_ref (head limit
        # 17.5 LU instead of 34.9). A second force level on the same jam type separates the jam's resistance from the
        # cap (leakage against level difference at two forces).
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-cap2p5-only) -- HALF CAP: the drifting-rods run\n"
               "// with PIDParameters.maxActuatingVariable 1.325e-6 (2.5 a_ref, head limit 17.5 LU);\n"
               "// everything else identical.\n")
        out = common(make_startblock(base))
        out, k = re.subn(r"(\n\s*maxActuatingVariable\s+)2\.65e-6;[^\n]*",
                         r"\g<1>1.325e-6;  // half cap: 2.5 a_ref (drifting rods: 2.65e-6 = 5.0 a_ref)", out)
        assert k == 1
        print(OUT_CYLCAP)
        open(OUT_CYLCAP, "w").write(hdr + out)
        return
    if "--giant-rods-startblock-only" in sys.argv:
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--giant-rods-startblock-only) -- 100 GIANT RODS:\n"
               "// R 6, L 36 drifting from upstream (front row ~143 LU upstream of the pillars, rear rows wrapped through the periodic\n"
               "// boundary, see make_giant_rods); otherwise the drifting-rods run (hopper, cap 5 a_ref).\n")
        out = common(*make_giant_rods(make_startblock(base)))
        print(OUT_GIANT)
        open(OUT_GIANT, "w").write(hdr + out)
        return
    if "--giant-rods-startblock-z36-only" in sys.argv:
        # giant rods with 36-cell blocks in z: a giant-rod run with <64, 20, 30> blew up at step ~33400 with a NaN that
        # started in a free-surface cell on the block (rank) boundary z = 60, where the backed-up water surface had
        # settled. Blocks of 36 cells put the boundaries at 36 / 72 / 108 / 144, away from the surface (49..60).
        # Only cellsPerBlock changes; 120 ranks.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--giant-rods-startblock-z36-only) -- 100 GIANT RODS,\n"
               "// 36-CELL BLOCKS: the giant-rod prm with cellsPerBlock <64, 20, 36> (12 x 2 x 5 = 120 ranks) instead of\n"
               "// <64, 20, 30>, so that no block boundary lies at z = 60, where the backed-up water surface settles.\n")
        out = common(*make_giant_rods(make_startblock(base)))
        out, k = re.subn(r"(\n\s*cellsPerBlock\s+)< 64, 20, 30 >;[^\n]*",
                         r"\g<1>< 64, 20, 36 >;  // z 30 -> 36: no block boundary at the water surface (z = 60)", out)
        assert k == 1
        print(OUT_GIANT_Z36)
        open(OUT_GIANT_Z36, "w").write(hdr + out)
        return
    if "--two-pillars-only" in sys.argv:
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--two-pillars-only) -- TWO PILLARS, 812 DRIFTING RODS: the drifting-rods\n"
               "// start block mirrored across a channel of twice the width (see make_two_pillars in the generator).\n"
               "// Release at step 500, VTK at step 0, wallOuterGhostLayer true, probes at z = 20.\n")
        out = common(*make_two_pillars(make_startblock(base)))
        print(OUT_TWOPIL)
        open(OUT_TWOPIL, "w").write(hdr + out)
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--two-pillars-only) -- TWO PILLARS, CLEAN REFERENCE: the two-pillar channel\n"
               "// of twopillar_clog_cyl812_startblock_cappid_gas.prm without floaters (a_ref, Q_ref, clean levels of this layout).\n")
        out = common(*make_two_pillars(make_startblock(base)))       # same settings as the rod run, then rods removed
        a0 = out.index("\nFloatingCylinders\n")
        a1 = out.index("\n}", a0) + 2
        out = out[:a0] + "\n" + out[a1:]
        out, k = re.subn(r"(\n\s*timesteps\s+)150000;[^\n]*", r"\g<1>40000;  // clean reference", out)
        assert k == 1 and "FloatingCylinder" not in out.split("\n", 3)[3].replace("FloatingCylinders", "")
        print(OUT_TWOPIL_CLEAN)
        open(OUT_TWOPIL_CLEAN, "w").write(hdr + out)
        return
    if "--centred-startblock-only" in sys.argv:
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--centred-startblock-only) -- CENTRED PILLAR, START BLOCK:\n"
               f"// {OUT_CENTRED} with the 100 spheres R 5 as a block at the start of the domain, released at step 500\n"
               "// (instead of one by one between 20000 and 56000); see make_centred_startblock.\n")
        out = make_centred_startblock(open(OUT_CENTRED).read())
        print(OUT_CENTRED_SB)
        open(OUT_CENTRED_SB, "w").write(hdr + out)
        return
    if "--hopper-clean-long-only" in sys.argv:
        # converged clean reference for all hopper-gap runs: the drifting-rods prm (same channel, pillars, probes, PID,
        # cap, wallOuterGhostLayer, VTK fields) with the FloatingCylinders block removed, 150000 steps, fluid VTK at 0
        # and 149000 only -- the hopper analogue of --two-pillars-clean-long-only. A 30000-step reference was too short:
        # its discharge still swung by about +-12 % over 5000-step windows.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--hopper-clean-long-only) -- HOPPER GAP, CLEAN REFERENCE: the drifting-rods\n"
               "// run (hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm) without floaters; clean a_ref, Q_ref and levels of the\n"
               "// hopper-gap layout for all hopper runs. 150000 steps, VTK at steps 0 and 149000 only.\n")
        out = common(make_startblock(base))
        a0 = out.index("\nFloatingCylinders\n")
        a1 = out.index("\n}", a0) + 2
        out = out[:a0] + "\n" + out[a1:]
        out, k = re.subn(r'steps "0,500,1000,[0-9,]*,149000";[^\n]*', 'steps "0,149000";  // clean reference: start and end only', out)
        assert k == 1 and "FloatingCylinder {" not in out
        print(OUT_HOPPER_CLEAN_LONG)
        open(OUT_HOPPER_CLEAN_LONG, "w").write(hdr + out)
        return
    if "--two-pillars-clean-long-only" in sys.argv:
        # longer clean two-pillar reference: 40000 steps did not settle -- in the wide, low-friction channel the flow
        # starts above the target and decays slowly, with the controller at 0 most of the time. Same channel, 150000
        # steps, fluid VTK only at 0 and 149000.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--two-pillars-clean-long-only) -- TWO PILLARS, CLEAN REFERENCE, LONG: as\n"
               "// twopillar_clean_reference_gas.prm but 150000 steps (40000 did not settle) and VTK at steps 0 and 149000 only.\n")
        out = common(*make_two_pillars(make_startblock(base)))
        a0 = out.index("\nFloatingCylinders\n")
        a1 = out.index("\n}", a0) + 2
        out = out[:a0] + "\n" + out[a1:]
        out, k = re.subn(r'steps "0,500,1000,[0-9,]*,149000";[^\n]*', 'steps "0,149000";  // clean reference: start and end only', out)
        assert k == 1 and "FloatingCylinder {" not in out
        print(OUT_TWOPIL_CLEAN_LONG)
        open(OUT_TWOPIL_CLEAN_LONG, "w").write(hdr + out)
        return
    if "--cyl-startblock-cap1p5-only" in sys.argv:
        # third force level: the drifting rods with the cap at 1.5 a_ref = 7.95e-7 (a_ref 5.30e-7 of the clean hopper;
        # head limit 10.5 LU). Does the block still jam when the drive is barely above the clean one? Also a third point
        # on leakage against head (with 2.5 and 5 a_ref).
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-cap1p5-only) -- CAP 1.5 a_ref: the drifting-rods run\n"
               "// with PIDParameters.maxActuatingVariable 7.95e-7 (1.5 a_ref, head limit 10.5 LU);\n"
               "// everything else identical.\n")
        out = common(make_startblock(base))
        out, k = re.subn(r"(\n\s*maxActuatingVariable\s+)2\.65e-6;[^\n]*",
                         r"\g<1>7.95e-7;  // cap 1.5 a_ref (drifting rods: 2.65e-6 = 5.0 a_ref)", out)
        assert k == 1
        print(OUT_CYLCAP15)
        open(OUT_CYLCAP15, "w").write(hdr + out)
        return
    if "--cyl-startblock-sigma-only" in sys.argv:
        # surface-tension scale test: the drifting-rods prm as written, with surfaceTensionOverride = sigma/10 and
        # sigma/100 (sigma = 1.63945e-3 from We = 15.6188, eval/computedParams.txt of the drifting-rods run). Bond
        # number of the rods g R^2/sigma: 0.32 -> 3.2 -> 32. Re, Fr unchanged (a larger system at the same Froude
        # number but the same Re, like a Froude-scaled flume); everything else identical.
        for fac, tag in ((10, "div10"), (100, "div100")):
            sig = 1.63945e-3 / fac
            hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-sigma-only) -- SURFACE TENSION / %d:\n"
                   "// the drifting-rods run with surfaceTensionOverride %.6g (We-derived 1.63945e-3 / %d; rod Bond number %.2g);\n"
                   "// everything else identical.\n" % (fac, sig, fac, 5.82575e-5 * 9 / sig))
            out, k = re.subn(r"(\n(\s*)useHardcodedUnits\s+false;[^\n]*)",
                             r"\g<1>\n\g<2>surfaceTensionOverride %.6g;  // sigma / %d (Bond number of the rods %.2g)"
                             % (sig, fac, 5.82575e-5 * 9 / sig), open(OUT_CYL).read())
            assert k == 1
            fn = OUT_SIGMA.format(tag=tag)
            print(fn)
            open(fn, "w").write(hdr + out)
        return
    if "--skin-pocket-test-only" in sys.argv:
        # cheap test of the sub-cell pocket mechanism (floaters only, small domain): the drifting-rods start block cut
        # to its first 6 columns (x 8..48, 84 rods: 2 across, 7 layers, same densities, row pitch, orientation, rod ends
        # 1 LU from the side walls), periodic channel 128 LU long, no pillars, same flow/drive/cap; skins 0.05 / 0.1 /
        # 0.25 / 0.5; 8000 steps; floater + fluid VTK every 250 steps; separate output folders per skin so that the four
        # runs can share one checkout and one job.
        base = open(OUT_CYL).read()
        for skin, tag in ((0.05, "0p05"), (0.1, "0p1"), (0.25, "0p25"), (0.5, "0p5")):
            out = base
            def sub(pat, rep, cnt=1, flags=0):
                nonlocal out
                out, k = re.subn(pat, rep, out, flags=flags)
                assert k == cnt, (pat, k)
            sub(r"(\n\s*cylinderContactSkin\s+)0\.5;[^\n]*", r"\g<1>%s;   // pocket test" % skin)
            sub(r"cellsPerBlock(\s+)< 64, 20, 30 >;", r"cellsPerBlock\g<1>< 32, 20, 48 >;   // 4 x 2 x 2 = 16 blocks")
            sub(r"domainSize(\s+)<768, 40, 180>;", r"domainSize\g<1><128, 40, 96>;")
            sub(r"(\n\s*storeSnapshot\s+)true;", r"\g<1>false;")
            sub(r"(\n\s*snapshotBaseFolder\s+)snapshots;", r"\g<1>snapshots_m%s;" % tag)
            sub(r"(\n\s*timesteps\s+)150000;[^\n]*", r"\g<1>8000;")
            sub(r"\n\s*CylinderPillar \{[^\n]*", "", cnt=2)
            sub(r"\n\s*BodyBC \{ shape Cylinder;[^\n]*", "", cnt=2)
            # keep only the first 6 columns of the block (x <= 48)
            def keep(m):
                return m.group(0) if float(m.group(1)) <= 48.0 else ""
            out, k = re.subn(r"\n\s*FloatingCylinder \{ position < ([\d.]+),[^\n]*", keep, out)
            assert out.count("FloatingCylinder {") == 84, out.count("FloatingCylinder {")
            sub(r"(\n\s*baseFolderName\s+)eval;", r"\g<1>eval_m%s;" % tag)
            sub(r"(upstreamPos\s+)< 336, 20, 20 >", r"\g<1>< 96, 20, 20 >")
            sub(r"(downstreamPos\s+)< 432, 20, 20 >", r"\g<1>< 120, 20, 20 >")
            sub(r"(pillarUpstreamPos\s+)< 336, 6, 20 >", r"\g<1>< 96, 6, 20 >")
            sub(r"(pillarDownstreamPos\s+)< 432, 6, 20 >", r"\g<1>< 120, 6, 20 >")
            sub(r"max < 768, 40, 180 >", "max < 128, 40, 96 >", cnt=2)
            sub(r"(\n\s*baseFolder\s+)vtk-out;", r"\g<1>vtk-out_m%s;" % tag, cnt=3)
            sub(r'steps "0,500,1000[^"]*";[^\n]*', 'steps "%s";' % ",".join(str(i) for i in range(0, 8001, 250)))
            sub(r"(\n\s*diagnoseNaNNormals\s+)false;", r"\g<1>true;", cnt=0) if "diagnoseNaNNormals" in out else None
            out = out.replace("   wallOuterGhostLayer           true;", "   wallOuterGhostLayer           true;\n   diagnoseNaNNormals            true;   // pocket test: world coordinates of NaN normals")
            hdr = ("// GENERATED by gen_hopper_stress_prms.py (--skin-pocket-test-only) -- SUB-CELL POCKET TEST, skin %s:\n"
                   "// first 6 columns (84 rods) of the drifting-rods start block in a 128 LU periodic channel, no pillars,\n"
                   "// 8000 steps; rod geometry, densities, flow, drive and cap as in the drifting rods. Output folders *_m%s.\n" % (skin, tag))
            fn = OUT_POCKET.format(tag=tag)
            print(fn)
            open(fn, "w").write(hdr + out)
        return
    if "--cyl-startblock-skin-sweep-only" in sys.argv:
        # skin sensitivity: the drifting-rods prm (skin 0.5) as written, with cylinderContactSkin 0.25 / 0.1 / 0.05;
        # rod geometry (row pitch 8) unchanged.
        for skin, tag in ((0.25, "0p25"), (0.1, "0p1"), (0.05, "0p05")):
            hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-skin-sweep-only) -- SKIN %s:\n"
                   "// the drifting-rods run (hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm, skin 0.5) with cylinderContactSkin %s;\n"
                   "// everything else identical (rod positions, row pitch 8, drive, cap).\n" % (skin, skin))
            out, k = re.subn(r"(\n\s*cylinderContactSkin\s+)0\.5;[^\n]*",
                             r"\g<1>%s;   // skin sensitivity run (drifting rods: 0.5)" % skin, open(OUT_CYL).read())
            assert k == 1
            fn = OUT_CYLSKIN.format(tag=tag)
            print(fn)
            open(fn, "w").write(hdr + out)
        return
    if "--cyl-startblock-drivebymass-only" in sys.argv:
        # drive by mass: the drifting-rods prm as written, with driveFloatersByMass true, so floaters are pushed
        # streamwise with m_p * a_x instead of V * a_x; nothing else changes.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-drivebymass-only) -- DRIVE BY MASS: the drifting-rods run\n"
               "// (hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm) with driveFloatersByMass true, i.e. the streamwise floater drive\n"
               "// is m_p * a_x instead of V * a_x.\n")
        out, k = re.subn(r"(\n(\s*)useMemOnlyBuoyancy\s+true;[^\n]*)",
                         r"\g<1>\n\g<2>driveFloatersByMass           true;  // streamwise floater drive m_p * a_x instead of V * a_x",
                         open(OUT_CYL).read())
        assert k == 1
        print(OUT_CYLMASS)
        open(OUT_CYLMASS, "w").write(hdr + out)
        return
    if "--cyl-startblock-rep2-only" in sys.argv:
        # second realisation of the drifting rods: every rod jittered by +-0.25 LU in x and z (seed 1); everything
        # else identical.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-rep2-only) -- REALISATION 2:\n"
               "// the drifting-rods run with every rod moved by a seeded random offset in [-0.25, 0.25] LU in x and z\n"
               "// (seed 1, jitter_rods); everything else identical.\n")
        print(OUT_CYLREP2)
        open(OUT_CYLREP2, "w").write(hdr + common(jitter_rods(make_startblock(base))))
        return
    if "--sphere4-startblock-only" in sys.argv:
        # spheres that can bridge: R 4, gap / diameter = 2; denser VTK because the passing spheres circulate (snapshots
        # 20-30k apart make the passing counts ambiguous).
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--sphere4-startblock-only) -- DRIFTING SPHERES R 4:\n"
               "// 600 spheres R 4 in the start-block extent of the drifting rods (x 5..229), gap / diameter = 2 (R 3: 2.7); see\n"
               "// make_spheres_r4. Release at step 500, wallOuterGhostLayer true, probes at z = 20; VTK every 2-4k steps to 40k.\n")
        sp, n = make_spheres_r4(make_startblock(base))
        out = common(sp, n)
        out, k = re.subn(r'steps "0,500,1000,[0-9,]*,149000";[^\n]*', f'steps "{SPH4_STEPS}";  // denser: circulating spheres', out)
        assert k == 1
        print(OUT_SPH4)
        open(OUT_SPH4, "w").write(hdr + out)
        return
    if "--cyl-startblock-noskin-only" in sys.argv:
        # skin control: the drifting rods with cylinderContactSkin 0 (does the skin prevent the breakdown?), nothing
        # else changed.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--cyl-startblock-noskin-only) -- SKIN CONTROL: the drifting-rods run\n"
               "// (hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm) with cylinderContactSkin 0 (contacts at the true rod surfaces).\n"
               "// Everything else identical.\n")
        out = common(make_startblock(base))
        out, k = re.subn(r"(\n\s*cylinderContactSkin\s+)0\.5;[^\n]*",
                         r"\g<1>0;     // skin control: contacts at the true rod surfaces", out)
        assert k == 1
        print(OUT_CYLNOSKIN)
        open(OUT_CYLNOSKIN, "w").write(hdr + out)
        return
    if "--sphere-startblock-only" in sys.argv:
        # drifting spheres: the sphere twin of the drifting rods / cuboids, same block position x 8..232, so that all
        # three shapes have the same upstream setup.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--sphere-startblock-only) -- HOPPER GAP, DRIFTING SPHERES: the drifting-rods\n"
               f"// start block ({BASE} shifted to x 8..232), each rod replaced by 5 spheres R 3 across\n"
               "// (1015 spheres, y = 4, 12, ..., 36), release at step 500, VTK at step 0, wallOuterGhostLayer true, probes at z = 20.\n")
        print(OUT_SPHDRIFT)
        sp, n = make_spheres(make_startblock(base))
        open(OUT_SPHDRIFT, "w").write(hdr + common(sp, n))
        return
    if "--box-startblock-only" in sys.argv:
        # drifting cuboids: the cuboid twin of the drifting rods -- block shifted to x 8..232, every rod replaced by a
        # box 6 x 18 x 6, boxContactSkin 0.5.
        hdr = ("// GENERATED by gen_hopper_stress_prms.py (--box-startblock-only) -- HOPPER GAP, DRIFTING CUBOIDS: the drifting-rods\n"
               f"// start block ({BASE} shifted to x 8..232) with 406 boxes 6 x 18 x 6 instead of the rods,\n"
               "// halfExtents <3, 9, 3> (long axis y), boxContactSkin 0.5, release at step 500, VTK at step 0, probes at z = 20.\n")
        print(OUT_BOXDRIFT)
        open(OUT_BOXDRIFT, "w").write(hdr + common(make_boxes(make_startblock(base))))
        return
    if "--mixed-only" in sys.argv:
        hdr_mix = ("// GENERATED by gen_hopper_stress_prms.py (--mixed-only) -- HOPPER GAP, MIXED FLOATERS: the drifting-rods start block\n"
                   f"// ({BASE} shifted to x 8..232) with the rows cycling cylinder / box / sphere from the front row,\n"
                   "// release at step 500, VTK at step 0, wallOuterGhostLayer true, probes at z = 20 (see make_mixed in the generator).\n")
        print(OUT_MIX)
        mixed, n = make_mixed(base)
        open(OUT_MIX, "w").write(hdr_mix + common(mixed, n))
        return
    assert "boxContactSkin" not in base and "wallOuterGhostLayer" not in base
    hdr_box = ("// GENERATED by gen_hopper_stress_prms.py -- CUBOID version of the 406-rod capped-PID hopper\n"
               f"// ({BASE}): boxes 6 x 18 x 6, boxContactSkin 0.5, closed-form SAT contact point, release at\n"
               "// step 500, VTK at step 0, wallOuterGhostLayer true, probes at z = 20. See the generator's docstring.\n")
    hdr_cyl = ("// GENERATED by gen_hopper_stress_prms.py -- HOPPER GAP, DRIFTING RODS: the 406-rod capped-PID hopper\n"
               f"// ({BASE}) with the whole block moved to the start of the domain (first row x = {X_START:g}),\n"
               "// release at step 500, VTK at step 0, wallOuterGhostLayer true. See the generator's docstring.\n")
    print(OUT_BOX)
    box = common(make_boxes(base))
    open(OUT_BOX, "w").write(hdr_box + box)
    print(OUT_CYL)
    cyl = common(make_startblock(base))
    open(OUT_CYL, "w").write(hdr_cyl + cyl)


if __name__ == "__main__":
    main()
