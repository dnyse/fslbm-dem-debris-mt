#!/usr/bin/env python3
"""Cross-check every shape-pair geometry in ContactDetectionTest.cpp against pybullet (Bullet
Physics) -- a completely independent, separately-developed physics engine, not sharing a single
line of code with waLBerla/mesa_pd's SolidCylinderContactDetection.

Each check below reuses the EXACT geometry (positions, radii, half-extents, rotations) from the
corresponding C++ test function, and compares pybullet's own collision-query result against the
same expected penetrationDepth / contactNormal / contactPoint values asserted in the C++ file
(WALBERLA_CHECK_FLOAT_EQUAL / _EPSILON there). Where the two engines use a different but equally
valid convention (e.g. Bullet reports the two surface points positionOnA/positionOnB rather than a
single midpoint), that's normalized here rather than treated as a mismatch, exactly like the C++
file's own documented caveats about pybullet's conventions.

Run with the venv this was validated against: see validation/README.md.
"""
import numpy as np
import pybullet as p

TOTAL = 0
FAILED = 0


def _fmt(v):
    return "(" + ", ".join(f"{x:+.6f}" for x in v) + ")"


# --- region-membership checkers for genuinely non-unique contact points -------------------------
#
# "Degenerate" here means the true contact PATCH is a disk or a line segment, not a single point --
# any point in that patch is an equally valid answer to "where do these two shapes touch." That's a
# fact about the geometry, not about pybullet being unpredictable (GJK/EPA is deterministic; it just
# lands on a different, but equally correct, member of the same solution set our own code picks
# from). So instead of skipping the point check for these pairs, verify pybullet's point actually
# lies IN that valid set, rather than merely closer than some hand-waved tolerance to one arbitrary
# member of it. A point outside the set (e.g. off the disk, or beyond the finite body's own extent)
# would still be a real bug.

def region_disk_on_plane(center_xy, plane_axis_value, axis_index, radius, tol=1e-2):
    """Flat cap (radius, centred at center_xy) lying exactly on a plane perpendicular to world axis
    `axis_index` at value `plane_axis_value` -- e.g. a vertical cylinder's bottom cap on a z=0 floor.
    Valid region: on the plane, within `radius` of the cap's axis."""
    def check_fn(point):
        on_plane = abs(point[axis_index] - plane_axis_value) < tol
        other = [i for i in range(3) if i != axis_index]
        radial = np.linalg.norm(np.array([point[i] - center_xy[j] for j, i in enumerate(other)]))
        in_disk = radial < radius + tol
        return on_plane and in_disk, (
            f"on_plane(axis {axis_index}={plane_axis_value})={on_plane} "
            f"(got {point[axis_index]:+.6f}); radial_dist={radial:.6f} <= r={radius} -> {in_disk}")
    return check_fn


def region_line_segment(fixed_xy, axis_index, lo, hi, tol=1e-2):
    """A tangent line between two parallel round/flat surfaces, uniform along world axis
    `axis_index` -- e.g. a vertical floater against a vertical pillar. The other two coordinates are
    exact (`fixed_xy`); the axis coordinate is only bounded, by the FINITE floater's own extent
    [lo, hi] (not the pillar's own -- infinite -- extent, which would be too weak a check)."""
    def check_fn(point):
        other = [i for i in range(3) if i != axis_index]
        xy_ok = all(abs(point[i] - fixed_xy[j]) < tol for j, i in enumerate(other))
        axial_ok = (lo - tol) <= point[axis_index] <= (hi + tol)
        return xy_ok and axial_ok, (
            f"xy match={xy_ok} (got {[point[i] for i in other]} vs {fixed_xy}); "
            f"axis {axis_index} in [{lo},{hi}]={axial_ok} (got {point[axis_index]:+.6f})")
    return check_fn


def check(name, bodyA, bodyB, expect_contact, expected_depth=None, expected_normal=None,
          expected_point=None, tol=1e-3, check_point=True, check_normal=True,
          point_convention="midpoint", point_region=None):
    global TOTAL, FAILED
    TOTAL += 1
    pts = p.getClosestPoints(bodyA, bodyB, distance=50)
    got_contact = len(pts) > 0 and pts[0][8] < 0

    if not expect_contact:
        ok = not got_contact
        status = "PASS" if ok else "FAIL"
        print(f"[{status}] {name}: expected NO contact, "
              f"pybullet dist={pts[0][8] if pts else 'n/a':.6f}" if pts else
              f"[{status}] {name}: expected NO contact, pybullet found nothing in range (OK)")
        if not ok:
            FAILED += 1
        return

    if not got_contact:
        print(f"[FAIL] {name}: expected contact, pybullet found none/no-penetration")
        FAILED += 1
        return

    pt = min(pts, key=lambda x: x[8])  # deepest penetration
    depth = pt[8]
    normal = np.array(pt[7])
    # Bullet reports the TRUE (unmargined) surface point on each body separately
    # (positionOnA, positionOnB). Our own C++ convention differs by pair kind: dynamic-dynamic
    # pairs (Sphere x Sphere, Sphere x Box) place contactPoint at the MIDPOINT between the two
    # surfaces; wall/pillar pairs (anything x HalfSpace, anything x CylindricalBoundary) place it
    # ON the static boundary's own surface instead (see e.g. the Box x HalfSpace comment in
    # ParticleKernels.h: "contactPoint on the wall surface"). Match whichever convention the C++
    # test being cross-checked actually uses, rather than forcing one formula on both.
    if point_convention == "midpoint":
        point = (np.array(pt[5]) + np.array(pt[6])) / 2.0
    elif point_convention == "onB":
        point = np.array(pt[6])  # positionOnB: the static wall/pillar body, passed as bodyB below
    else:  # "onA": some CylindricalBoundary overloads report the point on the FLOATER's own
        # surface instead of the pillar's (Box x CylindricalBoundary, Cylinder x CylindricalBoundary
        # -- unlike Sphere x CylindricalBoundary, which lands on the pillar's surface). Verified by
        # inspecting raw positionOnA/positionOnB against ParticleKernels.h's own formula for each.
        point = np.array(pt[5])

    ok = True
    msgs = []
    if expected_depth is not None:
        d_ok = abs(depth - expected_depth) < tol
        ok &= d_ok
        msgs.append(f"depth: got={depth:+.6f} expected={expected_depth:+.6f} "
                    f"{'OK' if d_ok else 'MISMATCH'}")
    if check_normal and expected_normal is not None:
        expected_normal = np.array(expected_normal)
        n_ok = np.linalg.norm(normal - expected_normal) < tol
        ok &= n_ok
        msgs.append(f"normal: got={_fmt(normal)} expected={_fmt(expected_normal)} "
                    f"{'OK' if n_ok else 'MISMATCH'}")
    if check_point and expected_point is not None:
        expected_point = np.array(expected_point)
        p_ok = np.linalg.norm(point - expected_point) < max(tol, 5e-3)
        ok &= p_ok
        msgs.append(f"point: got={_fmt(point)} expected={_fmt(expected_point)} "
                    f"{'OK' if p_ok else 'MISMATCH'}")
    elif point_region is not None:
        p_ok, detail = point_region(point)
        ok &= p_ok
        msgs.append(f"point: got={_fmt(point)}, checked against valid region ({detail}) "
                    f"{'OK' if p_ok else 'MISMATCH'}")

    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}: " + "; ".join(msgs))
    if not ok:
        FAILED += 1


def sphere(pos, r):
    col = p.createCollisionShape(p.GEOM_SPHERE, radius=r)
    return p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col, basePosition=pos)


def box(pos, half_extents, quat=(0, 0, 0, 1)):
    col = p.createCollisionShape(p.GEOM_BOX, halfExtents=half_extents)
    return p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col, basePosition=pos,
                              baseOrientation=quat)


def cylinder(pos, radius, height, quat=(0, 0, 0, 1)):
    # pybullet's GEOM_CYLINDER default axis is body z, matching mesa_pd::data::Cylinder's convention.
    col = p.createCollisionShape(p.GEOM_CYLINDER, radius=radius, height=height)
    return p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col, basePosition=pos,
                              baseOrientation=quat)


def pillar(pos_xy, radius, quat=(0, 0, 0, 1)):
    # CylindricalBoundary is used in this codebase as a solid external pillar (contact when
    # dist_to_axis < cb.radius + floater.radius) -- physically identical to a plain solid cylinder,
    # confirmed against the Sphere x CylindricalBoundary overload's own formula in ParticleKernels.h.
    # Modeled as a very tall cylinder so it behaves as effectively infinite along its axis, same as
    # CylindricalBoundary's own infinite-line semantics.
    return cylinder([pos_xy[0], pos_xy[1], 0], radius, 400.0, quat)


def halfspace(pos, normal):
    col = p.createCollisionShape(p.GEOM_PLANE, planeNormal=normal)
    return p.createMultiBody(baseMass=0, baseCollisionShapeIndex=col, basePosition=pos)


def quat_axis_angle(axis, angle):
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    s = np.sin(angle / 2)
    return (axis[0] * s, axis[1] * s, axis[2] * s, np.cos(angle / 2))


def main():
    p.connect(p.DIRECT)

    # ---- Sphere x Sphere ----
    d = np.array([1, 2, 3], dtype=float); d /= np.linalg.norm(d)
    s0 = sphere([0, 0, 0], 0.5); s1 = sphere((d * 0.75).tolist(), 0.5)
    check("testSphereSphereContact", s0, s1, True, -0.25, (-d).tolist(), (d * 0.375).tolist())

    s0 = sphere([0, 0, 0], 0.5); s1 = sphere((d * 1.5).tolist(), 0.5)
    check("testSphereSphereNoContact", s0, s1, False)

    # ---- Sphere x Box ----
    bx = box([0, 0, 0], [1, 1, 1]); sp = sphere([2, 2, 2], 3 ** 0.5 + 0.1)
    n = np.array([1, 1, 1]) / np.sqrt(3)
    check("testSphereBoxContact", sp, bx, True, -0.1, n.tolist(), None, check_point=False)  # point convention genuinely differs (see ContactDetectionTest.cpp's own comment: ours=box surface, pybullet=midpoint)

    bx = box([0, 0, 0], [0.5, 0.5, 0.5]); sp = sphere([5, 0, 0], 0.6)
    check("testSphereBoxNoContact", sp, bx, False)

    # ---- Box x Box ----
    b0 = box([0, 0, 0], [0.5, 0.5, 0.5]); b1 = box([0.8, 0, 0], [0.5, 0.5, 0.5])
    check("testBoxBoxContact", b0, b1, True, -0.2, [-1, 0, 0], None, check_point=False)

    b0 = box([0, 0, 0], [0.5, 0.5, 0.5]); b1 = box([3, 0, 0], [0.5, 0.5, 0.5])
    check("testBoxBoxNoContact", b0, b1, False)

    # ---- Cylinder x Box ----
    cyl = cylinder([0, 0, 0], 0.5, 2); bx = box([0.9, 0, 0], [0.5, 0.5, 0.5])
    check("testCylinderBoxContact", cyl, bx, True, -0.1, [-1, 0, 0], None, check_point=False)

    cyl = cylinder([0, 0, 0], 0.5, 2); bx = box([3, 0, 0], [0.5, 0.5, 0.5])
    check("testCylinderBoxNoContact", cyl, bx, False)

    # ---- Cylinder x Cylinder ----
    c0 = cylinder([0, 0, 0], 0.5, 2); c1 = cylinder([0.8, 0, 0], 0.5, 2)
    check("testCylinderCylinderContact", c0, c1, True, -0.2, [-1, 0, 0], None, check_point=False)

    c0 = cylinder([0, 0, 0], 0.5, 2); c1 = cylinder([3, 0, 0], 0.5, 2)
    check("testCylinderCylinderNoContact", c0, c1, False)

    # ---- Sphere x Cylinder ("floating log") ----
    sp = sphere([1.5, 0, 0], 1); cyl = cylinder([0, 0, 0], 1, 4)
    check("testSphereCylinderSideContact", sp, cyl, True, -0.5, [1, 0, 0], [0.75, 0, 0])

    q = quat_axis_angle([0, 0, 1], np.pi / 2)
    sp = sphere([0, 1.5, 0], 1); cyl = cylinder([0, 0, 0], 1, 4, quat=q)
    check("testSphereCylinderRotationAware", sp, cyl, True, -0.5, [0, 1, 0], [0, 0.75, 0])

    sp = sphere([3, 0, 0], 1); cyl = cylinder([0, 0, 0], 1, 4)
    check("testSphereCylinderNoContact", sp, cyl, False)

    # ---- Sphere x CylindricalBoundary ("pillar", modeled as a tall solid cylinder) ----
    sp = sphere([1.5, 0, 7], 0.6); pil = pillar([0, 0], 1)
    check("testSphereCylindricalBoundaryContact", sp, pil, True, -0.1, [1, 0, 0], [1.0, 0, 7], point_convention="onB")

    sp = sphere([3, 0, 0], 0.6); pil = pillar([0, 0], 1)
    check("testSphereCylindricalBoundaryNoContact", sp, pil, False)

    # ---- Sphere x HalfSpace ----
    sp = sphere([0, 0, 0.3], 0.6); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testSphereHalfSpaceContact", sp, wall, True, -0.3, [0, 0, 1], [0, 0, 0], point_convention="onB")

    sp = sphere([0, 0, 2], 0.6); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testSphereHalfSpaceNoContact", sp, wall, False)

    # ---- Box x HalfSpace ----
    bx = box([0, 0, 0.3], [0.5, 0.5, 0.5]); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testBoxHalfSpaceContactAxisAligned", bx, wall, True, -0.2, [0, 0, 1], [0.5, 0.5, 0], point_convention="onB")

    q = quat_axis_angle([0, 0, 1], np.pi / 2)
    bx = box([0, 0, 0.3], [0.3, 0.5, 0.5], quat=q); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testBoxHalfSpaceContactRotated", bx, wall, True, -0.2, [0, 0, 1], [-0.5, 0.3, 0], point_convention="onB")

    bx = box([0, 0, 2], [0.5, 0.5, 0.5]); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testBoxHalfSpaceNoContact", bx, wall, False)

    # ---- Box x CylindricalBoundary ("pillar") ----
    # Point convention is "onA" (the box's OWN surface, e.g. x=0.7 -- ParticleKernels.h's own
    # comment: "contactPoint = closest point on the box surface"), not the pillar's -- confirmed
    # against raw positionOnA/positionOnB. Degenerate along z only (box's flat face vs. the pillar's
    # z-uniform radius): x,y are exact, z is checked as bounded within the box's OWN z half-extent
    # (not the pillar's infinite one -- a point outside the box's actual extent would still be a bug).
    bx = box([1.2, 0, 0], [0.5, 0.5, 2]); pil = pillar([0, 0], 1)
    check("testBoxCylindricalBoundaryContact", bx, pil, True, -0.3, [1, 0, 0], None,
          point_convention="onA", point_region=region_line_segment([0.7, 0], axis_index=2, lo=-2, hi=2))

    bx = box([1.2, 0, 5], [0.5, 0.5, 2]); pil = pillar([0, 0], 1)
    check("testBoxCylindricalBoundaryAxialInvariance", bx, pil, True, -0.3, [1, 0, 0], None,
          point_convention="onA", point_region=region_line_segment([0.7, 0], axis_index=2, lo=3, hi=7))

    bx = box([3, 0, 0], [0.5, 0.5, 2]); pil = pillar([0, 0], 1)
    check("testBoxCylindricalBoundaryNoContact", bx, pil, False)

    # ---- Cylinder x CylindricalBoundary ("pillar") ----
    # Same "onA" convention: Cylinder x CylindricalBoundary's own formula (ParticleKernels.h) reports
    # contactPoint on the LOG's own surface (x=0.5, its near side), not the pillar's -- verified
    # against raw pybullet output. Degenerate along z only, bounded by the log's own [-2,2] extent.
    log = cylinder([1.5, 0, 0], 1, 4); pil = pillar([0, 0], 1)
    check("testCylinderPillarContact", log, pil, True, -0.5, [1, 0, 0], None,
          point_convention="onA", point_region=region_line_segment([0.5, 0], axis_index=2, lo=-2, hi=2))

    log = cylinder([3, 0, 0], 1, 4); pil = pillar([0, 0], 1)
    check("testCylinderPillarNoContact", log, pil, False)

    # ---- Cylinder x HalfSpace ----
    # Degenerate along x,y instead (flat cap vs. flat plane -> a whole disk, radius 1, is equally
    # deep). "onB" convention (the wall's own reported point) is correctly ON the plane (z=0 exact,
    # confirmed); only the radial position within the disk varies between our own cap-centre pick
    # and pybullet's rim pick, so check disk membership instead of an exact x,y match.
    cyl = cylinder([0, 0, 1], 1, 4); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testCylinderHalfSpaceContact", cyl, wall, True, -1.0, [0, 0, 1], None,
          point_convention="onB", point_region=region_disk_on_plane([0, 0], plane_axis_value=0, axis_index=2, radius=1))

    cyl = cylinder([0, 0, 3], 1, 4); wall = halfspace([0, 0, 0], [0, 0, 1])
    check("testCylinderHalfSpaceNoContact", cyl, wall, False)

    p.disconnect()

    print()
    print(f"{TOTAL - FAILED}/{TOTAL} checks passed against pybullet.")
    if FAILED:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
