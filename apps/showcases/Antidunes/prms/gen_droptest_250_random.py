#!/usr/bin/env python3
"""[Claude, 2026-10-01] Generates the randomised 250-floater drop tests (author's request after the first reruns ended as
ten vertical towers: axis-aligned floaters released at ten fixed points fall straight down and stack face on face).

   DropTest_250_random_cylinders_gas.prm   250 cylinders (R 1, L 12)
   DropTest_250_random_cuboids_gas.prm     250 cuboids (half-extents 6 x 1 x 1.5)
   DropTest_250_random_mixed_gas.prm       84 cuboids, 83 cylinders, 83 spheres (R 2), cycled in release order

Base: DropTest_250_sat_gas.prm (container 32 x 30 x 180, gravity, contact detection settings). Changes:
 - Orientation uniformly random (uniform unit quaternion -> axis-angle; cylinders: axis = rotated body axis), seeded.
 - Release: 4 points per wave, (8, 7.5), (24, 7.5), (8, 22.5), (24, 22.5) at z = 170, each with a random lateral offset of
   up to +-0.75; centres of one wave are >= 13.5 apart, more than the largest circumscribed diameter (12.5), and >= 6.75
   from the side walls (largest circumscribed radius 6.26). 63 waves every 1500 steps [revised 2026-10-02: 1200, see below] (a released floater falls
   0.5 g t^2 = 24 cells in that time), last release at step 93000, then 10000 steps to settle: timesteps 103000.
 - Wood-like contact: collisionTimeFactor 1/8.0019 -> collision time 1 LBM step = 20 DEM sub-steps (Rettinger 2023,
   Sec. 4.3.4: >= 10 steps per collision for an accurate restitution; was 32 LBM steps), restitution 0.25, friction 0.3
   (driftwood DEM of Watanabe et al. 2023, Sec. 3.2).
 - Blocks 16 x 15 x 30 (2 x 2 x 6 = 24 processes): every block wider than the particle sync reach of 12.5 (see the
   SyncNextNeighbors fix in SetupStages*.h, THESIS_PLAN.md 2026-10-01).
 - [2026-10-02] Revised after the first runs timed out (3 h; cuboids at step 60006, mixed at 93370; no errors):
   blocks 16 x 15 x 15 (2 x 2 x 12 = 48 processes) -- the pile (~25-35 LU high) lay in the bottom block layer alone, so
   4 of 24 processes did all contact work and the rest waited in the particle sync (~75 % of all process time); now
   it spans two layers. Every block is still wider than the sync reach of 12.5. Waves every 1200 steps instead of 1500
   (a floater falls 0.5 g t^2 = 15.4 LU in that time, more than the 12.5 needed before the next release at the same
   point): last release at step 74400, then 10000 steps to settle -> timesteps 84400.
Usage: python3 gen_droptest_250_random.py  (from this directory)
"""
import math
import re

import numpy as np

SEED = 20261001
BASE = "DropTest_250_sat_gas.prm"
N, PER_WAVE, WAVE_DT, Z0, JITTER = 250, 4, 1200, 170.0, 0.75
POINTS = [(8.0, 7.5), (24.0, 7.5), (8.0, 22.5), (24.0, 22.5)]
SETTLE = 10000


def random_rotation(rng):
    q = rng.standard_normal(4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    if w < 0:
        w, x, y, z = -w, -x, -y, -z
    angle = 2 * math.degrees(math.acos(min(1.0, w)))
    s = math.sqrt(max(1e-16, 1 - w * w))
    axis = np.array([x, y, z]) / s
    R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    return axis, angle, R


def fmt(v):
    return ", ".join(f"{c:.6f}" for c in v)


def main():
    base = open(BASE).read()
    dens = [float(m) for m in re.findall(r"FloatingBox \{[^}]*densityRatio ([0-9.]+);", base)]
    assert len(dens) == N
    start = base.index("\nFloatingBoxes")
    end = base.index("\n}", start) + 2
    for case in ("cylinders", "cuboids", "mixed"):
        rng = np.random.default_rng(SEED)
        boxes, cyls, sphs = [], [], []
        for k in range(N):
            wave, slot = divmod(k, PER_WAVE)
            px, py = POINTS[slot]
            pos = (px + rng.uniform(-JITTER, JITTER), py + rng.uniform(-JITTER, JITTER), Z0)
            axis, angle, R = random_rotation(rng)
            st = wave * WAVE_DT
            tail = f"densityRatio {dens[k]}; settlingSteps {st};" + (" startInactive true;" if st > 0 else "")
            shape = {"cylinders": "cyl", "cuboids": "box"}.get(case) or ("box", "cyl", "sph")[k % 3]
            if shape == "box":
                boxes.append(f"   FloatingBox {{ position < {fmt(pos)} >; halfExtents < 6, 1, 1.5 >; "
                             f"initialTiltAxis < {fmt(axis)} >; initialTiltAngleDeg {angle:.6f}; {tail} }}")
            elif shape == "cyl":
                # FloatingCylinderParams: the axis is set directly (spin about the axis is irrelevant for a cylinder)
                a = R @ np.array([1.0, 0.0, 0.0])  # same body direction as the base prm's axis <1,0,0>
                cyls.append(f"   FloatingCylinder {{ position < {fmt(pos)} >; radius 1.0; length 12; axis < {fmt(a)} >; {tail} }}")
            else:
                sphs.append(f"   FloatingSphere {{ position < {fmt(pos)} >; radius 2.0; {tail} }}")
        blocks = ""
        for name, items in (("FloatingBoxes", boxes), ("FloatingCylinders", cyls), ("FloatingSpheres", sphs)):
            if items:
                blocks += f"\n{name}\n{{\n" + "\n".join(items) + "\n}\n"
        s = base[:start] + blocks + base[end:]
        last = ((N - 1) // PER_WAVE) * WAVE_DT
        s = re.sub(r"(\n\s*timesteps\s+)\d+;[^\n]*", rf"\g<1>{last + SETTLE};  // last release {last} + {SETTLE}", s, 1)
        s = re.sub(r"(\n\s*cellsPerBlock\s+)<[^>]*>;[^\n]*", r"\1< 16, 15, 15 >;  // 2 x 2 x 12 = 48 processes", s, 1)
        s = re.sub(r"(\n\s*frictionCoefficient\s+)[0-9.]+;", r"\g<1>0.3;", s, 1)
        s = re.sub(r"(\n\s*restitutionCoefficient\s+)[0-9.]+;",
                   r"\g<1>0.25;\n   collisionTimeFactor         0.124969648;  // collision time 1 LBM step (= 1/avgParticleDiameterOverride)", s, 1)
        if case != "cuboids":
            s = re.sub(r"(\n\s*useSpheroBoxContact\s+)true;", r"\1false;", s, 1)
        hdr = (f"// [Claude, 2026-10-01] GENERATED by gen_droptest_250_random.py (seed {SEED}) -- randomised 250-floater "
               f"drop test, case '{case}'.\n// See the generator's docstring for every change relative to {BASE}.\n")
        open(f"DropTest_250_random_{case}_gas.prm", "w").write(hdr + s)
        print(case, len(boxes), len(cyls), len(sphs))


if __name__ == "__main__":
    main()
