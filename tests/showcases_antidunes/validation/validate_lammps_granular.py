#!/usr/bin/env python3
"""Cross-check the Sphere x Sphere and Sphere x HalfSpace contact cases in
ContactDetectionTest.cpp against LAMMPS's GRANULAR package -- a real, independent DEM/granular
contact-mechanics implementation (the actual physics LIGGGHTS forked from; PyGran itself has no
physics of its own, it's just a ctypes wrapper around a LIGGGHTS/LAMMPS shared library -- see
validation/README.md for why this script calls LAMMPS directly instead).

LAMMPS's GRANULAR package is fundamentally sphere-based (pair_style granular computes
sphere-sphere contact; fix wall/gran computes sphere-vs-flat-wall contact) -- it has no box or
cylinder shape, so it can only cross-check the two Sphere x {Sphere, HalfSpace} pairs, not the
full matrix (pybullet's validate_pybullet.py covers all 12 pairs instead).

What's actually being checked: given the exact same geometry as the corresponding C++ test, does
LAMMPS's own Hookean contact law (1) detect contact/no-contact at the same threshold our
penetrationDepth sign does, and (2) apply the resulting repulsive force along the same direction as
our contactNormal? Force magnitude is a further, free arithmetic check: F = kn * |overlap|, and
LAMMPS's own |overlap| must equal our |penetrationDepth| exactly (both are simple distance-minus-radii
geometry) -- an independent confirmation of the geometric quantity our code reports.

Requires: LD_LIBRARY_PATH pointing at an MPICH-ABI-compatible libmpi.so.12 (the pip lammps wheel's
own ABI) -- see validation/README.md.
"""
import numpy as np
from lammps import lammps

TOTAL = 0
FAILED = 0

KN = 1000.0
KT = 1000.0
GAMMA_N = 50.0
GAMMA_T = 50.0
XMU = 0.5


def new_lammps():
    L = lammps(cmdargs=["-log", "none", "-screen", "none"])
    L.commands_string("""
units si
atom_style sphere
boundary f f f
comm_modify vel yes
region box block -20 20 -20 20 -20 20
create_box 1 box
""")
    return L


def report(name, expect_contact, overlap_expected, force, normal_expected, tol=1e-6):
    global TOTAL, FAILED
    TOTAL += 1
    got_contact = np.linalg.norm(force) > 1e-9

    if not expect_contact:
        ok = not got_contact
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: expected NO contact, "
              f"LAMMPS force magnitude={np.linalg.norm(force):.6f}")
        if not ok:
            FAILED += 1
        return

    if not got_contact:
        print(f"[FAIL] {name}: expected contact, LAMMPS reported zero force")
        FAILED += 1
        return

    f_mag = np.linalg.norm(force)
    overlap_from_force = f_mag / KN  # Hooke's law: F = kn * overlap
    overlap_ok = abs(overlap_from_force - overlap_expected) < tol
    f_dir = force / f_mag
    normal_expected = np.array(normal_expected, dtype=float)
    dir_ok = np.linalg.norm(f_dir - normal_expected) < 1e-6

    ok = overlap_ok and dir_ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: "
          f"overlap (from F=kn*x): got={overlap_from_force:.6f} expected={overlap_expected:.6f} "
          f"{'OK' if overlap_ok else 'MISMATCH'}; "
          f"force direction: got={tuple(f_dir)} expected={tuple(normal_expected)} "
          f"{'OK' if dir_ok else 'MISMATCH'}")
    if not ok:
        FAILED += 1


def sphere_sphere(name, shift, expect_contact):
    # Mirrors testSphereSphereContact / testSphereSphereNoContact: two r=0.5 spheres, one at the
    # origin, one at distance `shift` along a fixed direction. Overlap = 2r - shift = 1 - shift.
    L = new_lammps()
    L.commands_string(f"""
create_atoms 1 single 0 0 0
create_atoms 1 single {shift} 0 0
set atom 1 diameter 1.0 density 2500
set atom 2 diameter 1.0 density 2500
pair_style granular
pair_coeff * * hooke {KN} {GAMMA_N} tangential linear_nohistory 0 0
fix 1 all nve/sphere
run 0
""")
    f = np.array(L.numpy.extract_atom("f"))
    L.close()
    overlap_expected = max(0.0, 1.0 - shift)
    # Force on atom 0 points away from atom 1, i.e. -x; matches our contactNormal=-dir convention
    # for testSphereSphereContact (dir here reduces to +x).
    report(name, expect_contact, overlap_expected, f[0], [-1, 0, 0])


def sphere_halfspace(name, z, radius, expect_contact):
    # Mirrors testSphereHalfSpaceContact / testSphereHalfSpaceNoContact: sphere centred at height z
    # above a wall at z=0, normal (0,0,1). Overlap = radius - z.
    L = new_lammps()
    L.commands_string(f"""
create_atoms 1 single 0 0 {z}
set atom 1 diameter {2*radius} density 2500
pair_style granular
pair_coeff * * hooke {KN} {GAMMA_N} tangential linear_nohistory 0 0
fix 1 all nve/sphere
fix 2 all wall/gran hooke {KN} {GAMMA_N} {KT} {GAMMA_T} {XMU} 1 zplane 0.0 NULL
run 0
""")
    f = np.array(L.numpy.extract_atom("f"))
    L.close()
    overlap_expected = max(0.0, radius - z)
    # Force pushes the sphere away from the wall, i.e. +z; matches our contactNormal=(0,0,1)
    # convention for testSphereHalfSpaceContact / testBoxHalfSpaceContactAxisAligned.
    report(name, expect_contact, overlap_expected, f[0], [0, 0, 1])


def main():
    sphere_sphere("testSphereSphereContact", shift=0.75, expect_contact=True)
    sphere_sphere("testSphereSphereNoContact", shift=1.5, expect_contact=False)
    sphere_halfspace("testSphereHalfSpaceContact", z=0.3, radius=0.6, expect_contact=True)
    sphere_halfspace("testSphereHalfSpaceNoContact", z=2.0, radius=0.6, expect_contact=False)

    print()
    print(f"{TOTAL - FAILED}/{TOTAL} checks passed against LAMMPS GRANULAR.")
    if FAILED:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
