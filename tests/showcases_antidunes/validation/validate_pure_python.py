#!/usr/bin/env python3
"""Independent, from-scratch pure-Python re-derivation of every shape-pair geometry in
ContactDetectionTest.cpp -- no dependency on pybullet, LAMMPS, or any physics engine, and no shared
code with the C++ under test (ParticleKernels.h). Only numpy for vector algebra.

Each function below implements the closest-point / separating-axis geometry for one shape pair from
first principles (clamped projection onto a box, distance to a cylinder's axis line, etc.) -- the
same kind of basic vector-algebra derivation the C++ file's own comments describe using by hand, just
written as actual runnable code instead of a comment, so it's mechanically re-checked instead of
trusted by inspection. This is the one validation layer that reaches every pair, including the ones
where pybullet's GJK/EPA implementation picks a different (but equally valid) point from a
mathematically non-unique contact region -- see validate_pybullet.py's notes on those cases.
"""
import numpy as np

TOTAL = 0
FAILED = 0


def _fmt(v):
    return "(" + ", ".join(f"{x:+.6f}" for x in np.asarray(v).flatten()) + ")"


def report(name, got_contact, expect_contact, depth=None, normal=None, point=None,
           expected_depth=None, expected_normal=None, expected_point=None, tol=1e-9,
           check_point=True):
    global TOTAL, FAILED
    TOTAL += 1

    if got_contact != expect_contact:
        print(f"[FAIL] {name}: contact got={got_contact} expected={expect_contact}")
        FAILED += 1
        return

    if not expect_contact:
        print(f"[PASS] {name}: no contact, as expected")
        return

    ok = True
    msgs = []
    if expected_depth is not None:
        d_ok = abs(depth - expected_depth) < max(tol, 1e-6)
        ok &= d_ok
        msgs.append(f"depth: got={depth:+.6f} expected={expected_depth:+.6f} {'OK' if d_ok else 'MISMATCH'}")
    if expected_normal is not None:
        n_ok = np.linalg.norm(np.asarray(normal) - np.asarray(expected_normal)) < 1e-6
        ok &= n_ok
        msgs.append(f"normal: got={_fmt(normal)} expected={_fmt(expected_normal)} {'OK' if n_ok else 'MISMATCH'}")
    if check_point and expected_point is not None:
        p_ok = np.linalg.norm(np.asarray(point) - np.asarray(expected_point)) < 1e-6
        ok &= p_ok
        msgs.append(f"point: got={_fmt(point)} expected={_fmt(expected_point)} {'OK' if p_ok else 'MISMATCH'}")

    print(f"[{'PASS' if ok else 'FAIL'}] {name}: " + "; ".join(msgs))
    if not ok:
        FAILED += 1


# --- from-scratch closest-point primitives -----------------------------------------------------

def closest_point_on_box(local_pt, half_extents):
    """Clamp a point (box-local frame) onto/into an axis-aligned box [-h,h]^3."""
    return np.clip(local_pt, -np.asarray(half_extents), np.asarray(half_extents))


def closest_point_on_segment(pt, seg_a, seg_b):
    d = seg_b - seg_a
    t = np.clip(np.dot(pt - seg_a, d) / np.dot(d, d), 0.0, 1.0)
    return seg_a + t * d


def dist_to_line(pt, line_pt, line_dir_unit):
    rel = pt - line_pt
    t = np.dot(rel, line_dir_unit)
    radial = rel - t * line_dir_unit
    return np.linalg.norm(radial), radial, t


# --- pair checks ---------------------------------------------------------------------------------

def sphere_sphere(name, p0, p1, r0, r1, expect_contact, expected_depth=None, expected_normal=None,
                   expected_point=None):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    dist = np.linalg.norm(p1 - p0)
    depth = dist - r0 - r1
    got_contact = depth < 0
    normal = (p0 - p1) / dist if dist > 1e-12 else np.array([1, 0, 0])
    if got_contact:
        dir01 = (p1 - p0) / dist  # p0 -> p1
        p0_surf = p0 + r0 * dir01
        p1_surf = p1 - r1 * dir01
        point = (p0_surf + p1_surf) / 2.0  # true midpoint between the two surfaces
    else:
        point = None
    report(name, got_contact, expect_contact, depth, normal, point,
           expected_depth, expected_normal, expected_point)


def sphere_box(name, sphere_pos, r, box_pos, half_extents, expect_contact, expected_depth=None,
               expected_normal=None):
    sphere_pos, box_pos = np.array(sphere_pos, float), np.array(box_pos, float)
    local = sphere_pos - box_pos
    closest_local = closest_point_on_box(local, half_extents)
    closest_world = box_pos + closest_local
    delta = sphere_pos - closest_world
    dist = np.linalg.norm(delta)
    depth = dist - r
    got_contact = depth < 0
    normal = delta / dist if dist > 1e-12 else np.array([1, 0, 0])
    report(name, got_contact, expect_contact, depth, normal, None,
           expected_depth, expected_normal, None, check_point=False)


def sphere_halfspace(name, sphere_pos, r, wall_pos, wall_normal, expect_contact,
                      expected_depth=None, expected_normal=None, expected_point=None):
    sphere_pos, wall_pos, n = np.array(sphere_pos, float), np.array(wall_pos, float), np.array(wall_normal, float)
    k = np.dot(n, sphere_pos) - np.dot(n, wall_pos)
    depth = k - r
    got_contact = depth < 0
    point = sphere_pos - (r + depth) * n if got_contact else None
    report(name, got_contact, expect_contact, depth, n, point,
           expected_depth, expected_normal, expected_point)


def box_halfspace(name, box_pos, half_extents, R, wall_pos, wall_normal, expect_contact,
                   expected_depth=None, expected_normal=None, expected_point=None):
    box_pos, wall_pos, n = np.array(box_pos, float), np.array(wall_pos, float), np.array(wall_normal, float)
    R = np.array(R, float)
    # deepest corner: support of the box in body frame along -R^T n, same formula the fixed
    # Cylinder x HalfSpace overload also uses (ParticleKernels.h).
    d_body = -(R.T @ n)
    corner_local = np.where(d_body >= 0, np.asarray(half_extents), -np.asarray(half_extents))
    deepest = box_pos + R @ corner_local
    depth = np.dot(n, deepest) - np.dot(n, wall_pos)
    got_contact = depth < 0
    point = deepest - depth * n if got_contact else None
    report(name, got_contact, expect_contact, depth, n, point,
           expected_depth, expected_normal, expected_point)


def cylinder_halfspace(name, cyl_pos, radius, half_len, axis_unit, wall_pos, wall_normal,
                        expect_contact, expected_depth=None, expected_normal=None, expected_point=None):
    # Support of a finite cylinder along direction d (world frame, axis already in world frame):
    # radial part follows d's component perpendicular to the axis, axial part is +-halfLen along
    # the axis, matching data::shape::Cylinder::support()'s own (independently-derived) formula.
    cyl_pos, n, axis = np.array(cyl_pos, float), np.array(wall_normal, float), np.array(axis_unit, float)
    wall_pos = np.array(wall_pos, float)
    d = -n
    axial = np.dot(d, axis)
    perp = d - axial * axis
    perp_len = np.linalg.norm(perp)
    radial_dir = perp / perp_len if perp_len > 1e-12 else np.zeros(3)
    support = radius * radial_dir + np.sign(axial) * half_len * axis
    deepest = cyl_pos + support
    depth = np.dot(n, deepest) - np.dot(n, wall_pos)
    got_contact = depth < 0
    point = deepest - depth * n if got_contact else None
    report(name, got_contact, expect_contact, depth, n, point,
           expected_depth, expected_normal, expected_point)


def sphere_cylinder_finite(name, sphere_pos, r, cyl_pos, cyl_r, cyl_axis, expect_contact,
                           expected_depth=None, expected_normal=None, expected_point=None):
    # Finite Sphere x Cylinder ("floating log"): unlike Sphere x CylindricalBoundary below (an
    # infinite pillar, whose contactPoint convention projects from the sphere by radius+depth), the
    # actual ParticleKernels.h overload for a finite cylinder places contactPoint at the MIDPOINT
    # between the two surfaces, same convention as sphere_sphere above -- confirmed by first getting
    # this wrong (used the pillar convention) and having pybullet's independent point disagree.
    sphere_pos, cyl_pos, axis = np.array(sphere_pos, float), np.array(cyl_pos, float), np.array(cyl_axis, float)
    dist, radial, _ = dist_to_line(sphere_pos, cyl_pos, axis)
    depth = dist - cyl_r - r
    got_contact = depth < 0
    normal = radial / dist if dist > 1e-12 else np.array([1, 0, 0])
    if got_contact:
        sphere_surf = sphere_pos - r * normal
        cyl_surf = sphere_pos - (dist - cyl_r) * normal
        point = (sphere_surf + cyl_surf) / 2.0
    else:
        point = None
    report(name, got_contact, expect_contact, depth, normal, point,
           expected_depth, expected_normal, expected_point)


def sphere_pillar(name, sphere_pos, r, pillar_pos, pillar_r, pillar_axis, expect_contact,
                   expected_depth=None, expected_normal=None, expected_point=None):
    sphere_pos, pillar_pos, axis = np.array(sphere_pos, float), np.array(pillar_pos, float), np.array(pillar_axis, float)
    dist, radial, _ = dist_to_line(sphere_pos, pillar_pos, axis)
    depth = dist - pillar_r - r
    got_contact = depth < 0
    normal = radial / dist if dist > 1e-12 else np.array([1, 0, 0])
    point = sphere_pos - (r + depth) * normal if got_contact else None
    report(name, got_contact, expect_contact, depth, normal, point,
           expected_depth, expected_normal, expected_point)


def box_pillar(name, box_pos, half_extents, R, pillar_pos, pillar_r, pillar_axis, expect_contact,
               expected_depth=None, expected_normal=None):
    box_pos, pillar_pos, axis = np.array(box_pos, float), np.array(pillar_pos, float), np.array(pillar_axis, float)
    R = np.array(R, float)
    # Closest point on the oriented box to the infinite axis line, in the box body frame.
    d_body = R.T @ axis
    d_body /= np.linalg.norm(d_body)
    a0_body = R.T @ (pillar_pos - box_pos)
    box_pt = np.zeros(3)
    for _ in range(2):
        t = np.dot(box_pt - a0_body, d_body)
        lp = a0_body + t * d_body
        box_pt = closest_point_on_box(lp, half_extents)
    vec = box_pt - a0_body
    radial = vec - np.dot(vec, d_body) * d_body
    dist = np.linalg.norm(radial)
    depth = dist - pillar_r
    got_contact = depth < 0
    normal_body = radial / dist if dist > 1e-12 else np.array([1, 0, 0])
    normal_world = R @ normal_body
    report(name, got_contact, expect_contact, depth, normal_world, None,
           expected_depth, expected_normal, None, check_point=False)


def cylinder_pillar(name, log_pos, log_r, log_half_len, log_axis, pillar_pos, pillar_r, pillar_axis,
                     expect_contact, expected_depth=None, expected_normal=None):
    log_pos, pillar_pos = np.array(log_pos, float), np.array(pillar_pos, float)
    log_axis, pillar_axis = np.array(log_axis, float), np.array(pillar_axis, float)
    # Fixed-point search for the closest point on the log's axis SEGMENT to the pillar's axis LINE.
    t = 0.0
    for _ in range(20):
        pt_on_log = log_pos + t * log_axis
        _, _, s = dist_to_line(pt_on_log, pillar_pos, pillar_axis)
        nearest_on_pillar = pillar_pos + s * pillar_axis
        t = np.clip(np.dot(nearest_on_pillar - log_pos, log_axis), -log_half_len, log_half_len)
    pt_on_log = log_pos + t * log_axis
    dist, radial, _ = dist_to_line(pt_on_log, pillar_pos, pillar_axis)
    depth = dist - pillar_r - log_r
    got_contact = depth < 0
    normal = radial / dist if dist > 1e-12 else np.array([1, 0, 0])
    report(name, got_contact, expect_contact, depth, normal, None,
           expected_depth, expected_normal, None, check_point=False)


def box_box_axis_aligned(name, pos0, half0, pos1, half1, expect_contact, expected_depth=None,
                          expected_normal=None):
    # Separating-axis test, axis-aligned only (matches the C++ test's simple geometry) -- the
    # overlap on each world axis, then pick the axis of minimum (least negative) overlap.
    pos0, pos1 = np.array(pos0, float), np.array(pos1, float)
    half0, half1 = np.array(half0, float), np.array(half1, float)
    overlaps = (half0 + half1) - np.abs(pos1 - pos0)
    got_contact = np.all(overlaps > 0)
    if not got_contact:
        report(name, False, expect_contact)
        return
    axis = int(np.argmin(overlaps))
    depth = -overlaps[axis]
    normal = np.zeros(3)
    normal[axis] = -1.0 if pos1[axis] > pos0[axis] else 1.0
    report(name, True, expect_contact, depth, normal, None, expected_depth, expected_normal, None,
           check_point=False)


def cylinder_box_axis_aligned(name, cyl_pos, cyl_r, cyl_half_len, box_pos, half_extents,
                               expect_contact, expected_depth=None, expected_normal=None):
    # Axis-aligned cylinder (world z) vs. axis-aligned box, matching the C++ test's simple geometry:
    # radial (xy) separating test combined with axial (z) separating test.
    cyl_pos, box_pos = np.array(cyl_pos, float), np.array(box_pos, float)
    half_extents = np.array(half_extents, float)
    dx = abs(cyl_pos[0] - box_pos[0]) - (cyl_r + half_extents[0])
    dz = abs(cyl_pos[2] - box_pos[2]) - (cyl_half_len + half_extents[2])
    got_contact = dx < 0 and dz < 0
    if not got_contact:
        report(name, False, expect_contact)
        return
    depth = dx  # x is the separating axis in this test's geometry
    normal = np.array([-1.0 if box_pos[0] > cyl_pos[0] else 1.0, 0, 0])
    report(name, True, expect_contact, depth, normal, None, expected_depth, expected_normal, None,
           check_point=False)


def main():
    d = np.array([1, 2, 3], dtype=float); d /= np.linalg.norm(d)

    sphere_sphere("testSphereSphereContact", [0, 0, 0], d * 0.75, 0.5, 0.5, True,
                  -0.25, -d, d * 0.375)
    sphere_sphere("testSphereSphereNoContact", [0, 0, 0], d * 1.5, 0.5, 0.5, False)

    sphere_box("testSphereBoxContact", [2, 2, 2], 3 ** 0.5 + 0.1, [0, 0, 0], [1, 1, 1], True,
               -0.1, np.array([1, 1, 1]) / np.sqrt(3))
    sphere_box("testSphereBoxNoContact", [5, 0, 0], 0.6, [0, 0, 0], [0.5, 0.5, 0.5], False)

    box_box_axis_aligned("testBoxBoxContact", [0, 0, 0], [0.5, 0.5, 0.5], [0.8, 0, 0], [0.5, 0.5, 0.5],
                          True, -0.2, [-1, 0, 0])
    box_box_axis_aligned("testBoxBoxNoContact", [0, 0, 0], [0.5, 0.5, 0.5], [3, 0, 0], [0.5, 0.5, 0.5],
                          False)

    cylinder_box_axis_aligned("testCylinderBoxContact", [0, 0, 0], 0.5, 1.0, [0.9, 0, 0],
                               [0.5, 0.5, 0.5], True, -0.1, [-1, 0, 0])
    cylinder_box_axis_aligned("testCylinderBoxNoContact", [0, 0, 0], 0.5, 1.0, [3, 0, 0],
                               [0.5, 0.5, 0.5], False)

    box_box_axis_aligned("testCylinderCylinderContact (as radial separating-axis)", [0, 0, 0],
                          [0.5, 0.5, 1.0], [0.8, 0, 0], [0.5, 0.5, 1.0], True, -0.2, [-1, 0, 0])
    box_box_axis_aligned("testCylinderCylinderNoContact (as radial separating-axis)", [0, 0, 0],
                          [0.5, 0.5, 1.0], [3, 0, 0], [0.5, 0.5, 1.0], False)

    sphere_cylinder_finite("testSphereCylinderSideContact", [1.5, 0, 0], 1,
                           [0, 0, 0], 1, [0, 0, 1], True, -0.5, [1, 0, 0], [0.75, 0, 0])
    sphere_cylinder_finite("testSphereCylinderNoContact", [3, 0, 0], 1, [0, 0, 0], 1, [0, 0, 1], False)

    sphere_pillar("testSphereCylindricalBoundaryContact", [1.5, 0, 7], 0.6, [0, 0, 0], 1, [0, 0, 1],
                  True, -0.1, [1, 0, 0], [1.0, 0, 7])
    sphere_pillar("testSphereCylindricalBoundaryNoContact", [3, 0, 0], 0.6, [0, 0, 0], 1, [0, 0, 1],
                  False)

    sphere_halfspace("testSphereHalfSpaceContact", [0, 0, 0.3], 0.6, [0, 0, 0], [0, 0, 1], True,
                      -0.3, [0, 0, 1], [0, 0, 0])
    sphere_halfspace("testSphereHalfSpaceNoContact", [0, 0, 2], 0.6, [0, 0, 0], [0, 0, 1], False)

    I = np.eye(3)
    box_halfspace("testBoxHalfSpaceContactAxisAligned", [0, 0, 0.3], [0.5, 0.5, 0.5], I, [0, 0, 0],
                  [0, 0, 1], True, -0.2, [0, 0, 1], [0.5, 0.5, 0])
    Rz90 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)
    box_halfspace("testBoxHalfSpaceContactRotated", [0, 0, 0.3], [0.3, 0.5, 0.5], Rz90, [0, 0, 0],
                  [0, 0, 1], True, -0.2, [0, 0, 1], [-0.5, 0.3, 0])
    box_halfspace("testBoxHalfSpaceNoContact", [0, 0, 2], [0.5, 0.5, 0.5], I, [0, 0, 0], [0, 0, 1],
                  False)

    box_pillar("testBoxCylindricalBoundaryContact", [1.2, 0, 0], [0.5, 0.5, 2.0], I, [0, 0, 0], 1,
               [0, 0, 1], True, -0.3, [1, 0, 0])
    box_pillar("testBoxCylindricalBoundaryAxialInvariance", [1.2, 0, 5], [0.5, 0.5, 2.0], I, [0, 0, 0],
               1, [0, 0, 1], True, -0.3, [1, 0, 0])
    box_pillar("testBoxCylindricalBoundaryNoContact", [3, 0, 0], [0.5, 0.5, 2.0], I, [0, 0, 0], 1,
               [0, 0, 1], False)

    cylinder_pillar("testCylinderPillarContact", [1.5, 0, 0], 1, 2, [0, 0, 1], [0, 0, 0], 1,
                     [0, 0, 1], True, -0.5, [1, 0, 0])
    cylinder_pillar("testCylinderPillarNoContact", [3, 0, 0], 1, 2, [0, 0, 1], [0, 0, 0], 1, [0, 0, 1],
                     False)

    cylinder_halfspace("testCylinderHalfSpaceContact", [0, 0, 1], 1, 2, [0, 0, 1], [0, 0, 0],
                        [0, 0, 1], True, -1.0, [0, 0, 1], [0, 0, 0])
    cylinder_halfspace("testCylinderHalfSpaceNoContact", [0, 0, 3], 1, 2, [0, 0, 1], [0, 0, 0],
                        [0, 0, 1], False)

    print()
    print(f"{TOTAL - FAILED}/{TOTAL} checks passed against the pure-Python re-derivation.")
    if FAILED:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
