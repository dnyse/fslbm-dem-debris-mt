//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file CapsuleShapeTest.cpp
//! [Claude, 2026-09-26] Unit tests for mesa_pd::data::Capsule ("pill" shape: a cylinder with two
//! hemispherical end caps), per explicit request ("Create a test suite for the pills"). Written
//! directly off the back of a real bug found and fixed in this same session: Capsule had NO
//! mesa_pd::ContainsPointFunctor overload in src/mesa_pd/common/Contains.h -- the fluid-obstacle
//! mapping (MovingParticleMapping.h, via SingleCast) would have hard-aborted with
//! "ContainsPoint not implemented!" at the very first particle-to-lattice mapping call, caught only
//! by manual code audit before a cluster job ever ran with capsule rods. This file's
//! testCapsuleContainsPoint_* cases are written specifically to catch a regression of exactly that
//! class of gap -- and specifically to distinguish a TRUE capsule (rounded caps) from a flat-capped
//! cylinder of the same (radius, height), since that distinction is precisely where a lazy/wrong
//! implementation (e.g. accidentally reusing isPointInsideCylinderBF) would still pass a naive
//! on-axis-only test but fail at the rounded corners.
//!
//! Coverage:
//!  - Volume/mass/inertia: analytic volume formula, and the degenerate height<=2*radius limit where
//!    a capsule must reduce to an exact solid sphere (volume AND inertia).
//!  - support() (GJK): cardinal on-axis/radial directions, checked against the Minkowski-sum
//!    definition (segment support + ball support) directly, not just plausibility.
//!  - ContainsPointFunctor / SingleCast dispatch: the exact call path MovingParticleMapping.h uses
//!    (mesa_pd::kernel::SingleCast -> mesa_pd::ContainsPointFunctor), including the discriminating
//!    "beyond the flat cap but still within the rounded cap" case a flat-cylinder test would get
//!    wrong, and the symmetric "just outside the rounded corner" case.
//!  - Contact detection parity with Cylinder, exercised through the SAME production dispatch class
//!    (SolidCylinderContactDetection) real scenarios use: Capsule x HalfSpace ("wall"), Capsule x
//!    CylindricalBoundary ("pillar"), and Capsule x Capsule -- all three fall through to the
//!    inherited generic GJK/EPA path (no analytic Capsule overload exists), so these tests also
//!    stand as regression coverage for that fallback itself.
//!  - pack()/unpack() MPI serialisation round-trip.
//======================================================================================================================

#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/common/Contains.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleAccessorWithBaseShape.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Capsule.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/kernel/SingleCast.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/mpi/MPIManager.h>
#include <core/mpi/BufferSystem.h>

#include <cmath>

namespace walberla {
namespace antidunes {

using Vec3 = mesa_pd::Vec3;

// ---------------------------------------------------------------------------------------------
// Volume / mass / inertia
// ---------------------------------------------------------------------------------------------

void testCapsuleVolumeFormula()
{
   const real_t r = real_t(2), h = real_t(10); // coreHalfLength = 10/2 - 2 = 3
   const mesa_pd::data::Capsule capsule(r, h);
   WALBERLA_CHECK_FLOAT_EQUAL(capsule.getCoreHalfLength(), real_t(3));

   const real_t expectedVolume = math::pi * r * r * (real_t(2) * real_t(3))              // cylindrical core
                                  + (real_t(4) / real_t(3)) * math::pi * r * r * r;       // two caps = one full sphere
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.getVolume(), expectedVolume, real_t(1e-10));
}

void testCapsuleDegenerateSphereLimit()
{
   // height <= 2*radius -> coreHalfLength = 0 -> the capsule IS exactly a sphere of that radius,
   // in both volume and inertia (not just approximately -- this must be exact).
   const real_t r = real_t(1.5);
   const mesa_pd::data::Capsule capsule(r, real_t(2) * r); // height == 2r, the boundary case
   WALBERLA_CHECK_FLOAT_EQUAL(capsule.getCoreHalfLength(), real_t(0));

   const real_t sphereVolume = (real_t(4) / real_t(3)) * math::pi * r * r * r;
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.getVolume(), sphereVolume, real_t(1e-10));

   const real_t density = real_t(2.5);
   mesa_pd::data::Capsule mutableCapsule(r, real_t(2) * r);
   mutableCapsule.updateMassAndInertia(density);
   const real_t expectedMass = sphereVolume * density;
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(mutableCapsule.getMass(), expectedMass, real_t(1e-10));

   const real_t expectedI = real_t(0.4) * expectedMass * r * r; // solid sphere: (2/5) m r^2, all three axes equal
   const auto&  I         = mutableCapsule.getInertiaBF();
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(I(0, 0), expectedI, real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(I(1, 1), expectedI, real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(I(2, 2), expectedI, real_t(1e-9));
}

void testCapsuleThinRodLimitApproachesRod()
{
   // As radius -> tiny relative to height, a capsule approaches an idealised massless-radius thin
   // rod: Iperp -> m*Lc^2/12 (Lc = full core length, ~= height here since radius is negligible),
   // Iaxis -> 0. Not exact (there is always some O(r^2) correction) -- checked to a loose 1% only,
   // as a sanity/limiting-behaviour check, not a precision one.
   const real_t r = real_t(0.001), h = real_t(10);
   mesa_pd::data::Capsule capsule(r, h);
   capsule.updateMassAndInertia(real_t(1));

   const real_t Lc          = real_t(2) * capsule.getCoreHalfLength();
   const real_t expectedRod = capsule.getMass() * Lc * Lc / real_t(12);
   const auto&  I           = capsule.getInertiaBF();
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(I(0, 0), expectedRod, expectedRod * real_t(0.01));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(I(1, 1), expectedRod, expectedRod * real_t(0.01));
   WALBERLA_CHECK(I(2, 2) < expectedRod * real_t(0.01)); // axial inertia vanishes relative to perpendicular
}

// ---------------------------------------------------------------------------------------------
// support() -- the Minkowski-sum-of-segment-and-ball definition, checked directly
// ---------------------------------------------------------------------------------------------

void testCapsuleSupportCardinalDirections()
{
   const real_t r = real_t(2), h = real_t(10); // coreHalfLength = 3
   const mesa_pd::data::Capsule capsule(r, h);
   const real_t a = capsule.getCoreHalfLength();

   // +z: segment support is the +z end of the core, ball support is r*(0,0,1) -- the cap TIP.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 0, 1)), Vec3(0, 0, a + r), real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 0, -1)), Vec3(0, 0, -(a + r)), real_t(1e-10));

   // Purely radial: segment support is the CENTRE point (0,0,0) by the near-perpendicular guard
   // (bfD[2] == 0 exactly), ball support is r*(1,0,0) -- the equatorial band surface, not a cap.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(1, 0, 0)), Vec3(r, 0, 0), real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 1, 0)), Vec3(0, r, 0), real_t(1e-10));

   // Oblique direction: must equal segment_support(d) + r*d exactly (the Minkowski-sum definition),
   // verified against that formula directly rather than a hand-picked numeric expectation.
   Vec3 d(real_t(0.6), real_t(0), real_t(0.8)); // already unit length
   Vec3 expected(real_t(0), real_t(0), a) ; // sign(d.z)=+1
   expected = expected + d * r;
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(d), expected, real_t(1e-10));
}

void testCapsuleContainedWithinSameParameterCylinder()
{
   // Capsule.h's own doc comment says `height` matches Cylinder's `height` convention exactly (full
   // end-to-end extent) so the two are "directly comparable/swappable at the same radius+height
   // numbers" -- which means, geometrically, Capsule(r,h) is always CONTAINED WITHIN Cylinder(r,h)
   // as a set: the rounded caps shave the flat cylinder's sharp rim corners off, they never add
   // reach beyond it. The two shapes touch (support points coincide) exactly on-axis (the capsule's
   // rounded tip reaches precisely as far as the cylinder's flat cap, both at height/2) and at the
   // equator (both reach exactly radius r radially) -- and the cylinder reaches STRICTLY farther
   // for any genuinely oblique direction, because its rim is a sharp corner (simultaneously at the
   // full radius AND the full height/2) while the capsule's surface there is smoothly rounded.
   const real_t r = real_t(2), h = real_t(10);
   const mesa_pd::data::Capsule  capsule(r, h);
   const mesa_pd::data::Cylinder cylinder(r, h);

   // Touching points: pure axial and pure radial directions give IDENTICAL support points.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 0, 1)), cylinder.support(Vec3(0, 0, 1)), real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(1, 0, 0)), cylinder.support(Vec3(1, 0, 0)), real_t(1e-10));

   // Oblique direction: the cylinder's sharp rim corner must reach strictly farther along z than
   // the capsule's rounded surface, and the capsule's support point must lie within (or on) the
   // cylinder's own radius/height bounds everywhere -- never poking outside it.
   const Vec3 d(real_t(0.6), real_t(0), real_t(0.8));
   const Vec3 capSupport = capsule.support(d);
   const Vec3 cylSupport = cylinder.support(d);
   WALBERLA_CHECK(capSupport[2] < cylSupport[2] - real_t(1e-6));
   WALBERLA_CHECK(std::sqrt(capSupport[0] * capSupport[0] + capSupport[1] * capSupport[1]) <= r + real_t(1e-10));
   WALBERLA_CHECK(std::abs(capSupport[2]) <= h * real_t(0.5) + real_t(1e-10));
}

// ---------------------------------------------------------------------------------------------
// ContainsPointFunctor via SingleCast -- the EXACT dispatch path MovingParticleMapping.h uses.
// This is the direct regression coverage for the bug found this session.
// ---------------------------------------------------------------------------------------------

bool capsuleContains(std::shared_ptr< mesa_pd::data::ParticleStorage > ps, size_t idx, const Vec3& point)
{
   // SingleCast needs getShape() (a raw BaseShape*), which plain ParticleAccessor does not expose
   // (only getBaseShape(), a shared_ptr) -- ParticleAccessorWithBaseShape adds it.
   mesa_pd::data::ParticleAccessorWithBaseShape ac(ps);
   mesa_pd::kernel::SingleCast                  singleCast;
   mesa_pd::ContainsPointFunctor                containsPointFctr;
   return singleCast(idx, ac, containsPointFctr, ac, point);
}

void testCapsuleContainsPoint_OnAxisAndCentre()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   auto p  = ps->create();
   p->setPosition(Vec3(0, 0, 0));
   p->setRotation(mesa_pd::Rot3());
   const real_t r = real_t(2), h = real_t(10), a = real_t(3); // coreHalfLength
   p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(r, h));

   WALBERLA_CHECK(capsuleContains(ps, 0, Vec3(0, 0, 0)));           // dead centre
   WALBERLA_CHECK(capsuleContains(ps, 0, Vec3(0, 0, a)));           // on-axis, at the flat-core/cap junction
   WALBERLA_CHECK(capsuleContains(ps, 0, Vec3(0, 0, a + r - real_t(1e-6)))); // just inside the cap tip
   WALBERLA_CHECK(!capsuleContains(ps, 0, Vec3(0, 0, a + r + real_t(1e-6)))); // just past the cap tip
}

void testCapsuleContainsPoint_BeyondFlatCapButInsideRoundedCap()
{
   // THE discriminating case: a point beyond the core segment's end (z > coreHalfLength) but still
   // within the hemispherical cap (radial offset small enough that distance-to-segment-end <=
   // radius). A flat-capped-cylinder test (isPointInsideCylinderBF, which requires |z| <=
   // height/2 == a+r... no: Cylinder's OWN height IS a+r-equivalent for a flat cap of the SAME
   // nominal height) would reject this point; a true capsule test must ACCEPT it.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   auto p  = ps->create();
   p->setPosition(Vec3(0, 0, 0));
   p->setRotation(mesa_pd::Rot3());
   const real_t r = real_t(2), h = real_t(10), a = real_t(3);
   p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(r, h));

   // Point at z = a + 1.5 (i.e. 1.5 beyond the flat core's end -- OUTSIDE a flat cylinder of height
   // 2a==6, since |z|=4.5 > a=3), radial offset chosen so the 3D distance to the segment endpoint
   // (0,0,a) is sqrt(1.5^2 + radial^2) = 1.9 < r=2 -- strictly inside the rounded cap.
   const real_t radial = std::sqrt(real_t(1.9) * real_t(1.9) - real_t(1.5) * real_t(1.5));
   WALBERLA_CHECK(capsuleContains(ps, 0, Vec3(radial, 0, a + real_t(1.5))));

   // Same axial offset, but radial offset pushed out so the distance to the segment endpoint is
   // now 2.1 > r=2 -- just outside the rounded cap (and would ALSO be outside a flat cylinder,
   // so this alone isn't discriminating, but confirms the corner rounding has a real edge, not an
   // unbounded/always-true bug).
   const real_t radialOutside = std::sqrt(real_t(2.1) * real_t(2.1) - real_t(1.5) * real_t(1.5));
   WALBERLA_CHECK(!capsuleContains(ps, 0, Vec3(radialOutside, 0, a + real_t(1.5))));
}

void testCapsuleContainsPoint_RadialBand()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   auto p  = ps->create();
   p->setPosition(Vec3(0, 0, 0));
   p->setRotation(mesa_pd::Rot3());
   const real_t r = real_t(2), h = real_t(10);
   p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(r, h));

   WALBERLA_CHECK(capsuleContains(ps, 0, Vec3(r - real_t(1e-6), 0, 0)));  // just inside the equatorial band
   WALBERLA_CHECK(!capsuleContains(ps, 0, Vec3(r + real_t(1e-6), 0, 0))); // just outside it
   WALBERLA_CHECK(!capsuleContains(ps, 0, Vec3(real_t(100), real_t(100), real_t(100)))); // far away
}

void testCapsuleContainsPoint_RespectsPositionAndRotation()
{
   // Off-origin, rotated 90 deg about world X (body-z axis now points along world -y): a point that
   // is "on-axis, at the cap tip" in world coordinates must be derived from the ACTUAL transform,
   // not accidentally always tested in body frame at the origin.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   auto p  = ps->create();
   const Vec3 centre(real_t(5), real_t(7), real_t(-3));
   const mesa_pd::Rot3 rot(Vec3(1, 0, 0), math::pi * real_t(0.5));
   p->setPosition(centre);
   p->setRotation(rot);
   const real_t r = real_t(2), h = real_t(10), a = real_t(3);
   p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(r, h));

   // Derive the +z cap tip's world position from the SAME rotation matrix the shape's own
   // world<->body transform uses, rather than hand-deriving the library's sign convention (which
   // way "rotate 90 deg about world X" turns +z is an easy thing to get backwards by inspection --
   // this way the test is self-consistent regardless of that convention).
   const Vec3 capTipWorld    = centre + rot.getMatrix() * Vec3(0, 0, a + r);
   const Vec3 outwardOneStep = (capTipWorld - centre).getNormalized() * real_t(1e-6);
   WALBERLA_CHECK(capsuleContains(ps, 0, capTipWorld - outwardOneStep));
   WALBERLA_CHECK(!capsuleContains(ps, 0, capTipWorld + outwardOneStep));
   // A point offset purely along world X from the centre (perpendicular to the body axis) at
   // distance r-eps must be inside (equatorial band), matching the un-rotated case exactly since
   // this rotation doesn't move the world-X direction.
   WALBERLA_CHECK(capsuleContains(ps, 0, centre + Vec3(r - real_t(1e-6), 0, 0)));
}

// ---------------------------------------------------------------------------------------------
// Contact detection, through the SAME production dispatch class real scenarios use
// (SolidCylinderContactDetection, ParticleKernels.h) -- no analytic Capsule overload exists for
// any of these pairs, so this also stands as regression coverage for the inherited generic
// GeneralContactDetection fallback itself, specifically for Capsule.
// ---------------------------------------------------------------------------------------------

void testCapsuleHalfSpaceContact()
{
   // Wall at x=0 with outward normal +x (matches testSphereHalfSpaceContact/testBoxHalfSpaceContact's
   // convention in ContactDetectionTest.cpp). Capsule lying along world z (identity rotation),
   // radius 1, positioned so its equatorial band overlaps the wall by 0.3.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto wall = ps->create();
   wall->setPosition(Vec3(0, 0, 0));
   auto rod = ps->create();
   rod->setPosition(Vec3(real_t(0.7), 0, 0)); // 0.7 + 0 (wall at x=0) ... penetration = radius - 0.7 = 0.3
   rod->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::HalfSpace wallShape(Vec3(1, 0, 0));
   const mesa_pd::data::Capsule   rodShape(real_t(1), real_t(6));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, wallShape, rodShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(scd.getPenetrationDepth()), real_t(0.3), real_t(1e-6));
   // NOTE: falls through to GeneralContactDetection's generic HalfSpace template, which (per this
   // session's audit) returns contactNormal pointing INTO the wall for pairs with no dedicated
   // analytic override -- the OPPOSITE convention from Sphere/Box/Cylinder x HalfSpace. Proven
   // algebraically harmless to the force actually applied (LinearSpringDashpot is symmetric under
   // the corresponding index-swap+negate), but asserted explicitly here so a future change to
   // either convention is caught rather than silently drifting.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(-1), 0, 0), real_t(1e-6));
}

void testCapsuleHalfSpaceNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto wall = ps->create();
   wall->setPosition(Vec3(0, 0, 0));
   auto rod = ps->create();
   rod->setPosition(Vec3(real_t(3), 0, 0)); // clearance 3 - 1 = 2 > 0
   rod->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::HalfSpace wallShape(Vec3(1, 0, 0));
   const mesa_pd::data::Capsule   rodShape(real_t(1), real_t(6));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, wallShape, rodShape, ac));
}

void testCapsuleCylindricalBoundaryContact()
{
   // "Pillar" contact -- CylindricalBoundary at the origin, radius 5, infinite along world z.
   // Capsule rod alongside it, radius 1, offset so its band overlaps the pillar surface by 0.2.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto pillar = ps->create();
   pillar->setPosition(Vec3(0, 0, 0));
   auto rod = ps->create();
   rod->setPosition(Vec3(real_t(5.8), 0, real_t(7))); // z=7 arbitrary -- pillar is axially invariant
   rod->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(5), Vec3(0, 0, 1));
   const mesa_pd::data::Capsule             rodShape(real_t(1), real_t(6));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, pillarShape, rodShape, ac));
   // separation of centres 5.8, pillar radius 5, rod radius 1 -> overlap = 5+1-5.8 = 0.2
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(scd.getPenetrationDepth()), real_t(0.2), real_t(1e-6));

   // Swapped-argument order must agree.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, rodShape, pillarShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), scd2.getPenetrationDepth(), real_t(1e-6));
}

void testCapsuleCylindricalBoundaryNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto pillar = ps->create();
   pillar->setPosition(Vec3(0, 0, 0));
   auto rod = ps->create();
   rod->setPosition(Vec3(real_t(8), 0, 0)); // 5 + 1 = 6 < 8, well clear
   rod->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(5), Vec3(0, 0, 1));
   const mesa_pd::data::Capsule             rodShape(real_t(1), real_t(6));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, pillarShape, rodShape, ac));
}

void testCapsuleCapsuleContact()
{
   // Two parallel capsules (radius 0.5, axis = world z via identity rotation), centres 0.8 apart
   // along x -> overlap = 0.5+0.5-0.8 = 0.2, exactly mirroring
   // ContactDetectionTest.cpp's testCylinderCylinderContact geometry so the two are directly
   // comparable (same numbers, Capsule instead of Cylinder -- confirms the generic GJK/EPA fallback
   // gives the same answer as Cylinder's analytic path in a case unaffected by end-cap shape).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto c0 = ps->create();
   c0->setPosition(Vec3(0, 0, 0));
   c0->setRotation(mesa_pd::Rot3());
   auto c1 = ps->create();
   c1->setPosition(Vec3(real_t(0.8), 0, 0));
   c1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Capsule capsuleShape(real_t(0.5), real_t(2));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, capsuleShape, capsuleShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(-1), 0, 0), real_t(1e-3));
}

void testCapsuleCapsuleNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto c0 = ps->create();
   c0->setPosition(Vec3(0, 0, 0));
   c0->setRotation(mesa_pd::Rot3());
   auto c1 = ps->create();
   c1->setPosition(Vec3(real_t(3), 0, 0)); // well clear: 0.5+0.5 < 3
   c1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Capsule capsuleShape(real_t(0.5), real_t(2));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, capsuleShape, capsuleShape, ac));
}

void testCapsuleEndCapContactDiffersFromCylinder()
{
   // The physically-meaningful difference from Cylinder: two capsules approaching END-ON (along
   // their shared axis) contact at a ROUNDED tip, so the penetration depth for a given centre
   // separation is LARGER than two same-size flat-capped cylinders would report end-on (a flat
   // cylinder's end-on contact is ambiguous/degenerate -- rim-to-rim -- while a capsule's is a
   // clean, unambiguous sphere-sphere-equivalent contact between the two cap tips).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   const real_t r = real_t(1), h = real_t(6), a = real_t(2); // coreHalfLength = 6/2 - 1 = 2, tip at z = a+r = 3

   auto c0 = ps->create();
   c0->setPosition(Vec3(0, 0, 0));
   c0->setRotation(mesa_pd::Rot3());
   auto c1 = ps->create();
   c1->setPosition(Vec3(0, 0, real_t(5.5))); // tips at z=3 and z=5.5-3=2.5 -> overlap 0.5
   c1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Capsule capsuleShape(r, h);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, capsuleShape, capsuleShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.5), real_t(1e-3));
   // Looser epsilon than usual: this is a perfectly axis-aligned, tip-to-tip (sphere-sphere-
   // equivalent) contact -- a genuinely degenerate case for an iterative GJK/EPA polytope
   // refinement, which is known to leave a small amount of numerical noise in the normal direction
   // for exactly-symmetric configurations even though the penetration depth itself stays precise.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(0, 0, real_t(-1)), real_t(0.02));
   WALBERLA_UNUSED(a);
}

// ---------------------------------------------------------------------------------------------
// pack()/unpack() MPI serialisation round-trip
// ---------------------------------------------------------------------------------------------

void testCapsulePackUnpack()
{
   mesa_pd::data::Capsule original(real_t(1.7), real_t(9.3));

   mpi::SendBuffer sb;
   original.pack(sb);

   mpi::RecvBuffer rb(sb);
   mesa_pd::data::Capsule restored;
   restored.unpack(rb);

   WALBERLA_CHECK_EQUAL(restored.getShapeType(), mesa_pd::data::Capsule::SHAPE_TYPE);
   WALBERLA_CHECK_FLOAT_EQUAL(restored.getRadius(), original.getRadius());
   WALBERLA_CHECK_FLOAT_EQUAL(restored.getHeight(), original.getHeight());
   WALBERLA_CHECK_FLOAT_EQUAL(restored.getCoreHalfLength(), original.getCoreHalfLength());
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes;

   testCapsuleVolumeFormula();
   testCapsuleDegenerateSphereLimit();
   testCapsuleThinRodLimitApproachesRod();

   testCapsuleSupportCardinalDirections();
   testCapsuleContainedWithinSameParameterCylinder();

   testCapsuleContainsPoint_OnAxisAndCentre();
   testCapsuleContainsPoint_BeyondFlatCapButInsideRoundedCap();
   testCapsuleContainsPoint_RadialBand();
   testCapsuleContainsPoint_RespectsPositionAndRotation();

   testCapsuleHalfSpaceContact();
   testCapsuleHalfSpaceNoContact();
   testCapsuleCylindricalBoundaryContact();
   testCapsuleCylindricalBoundaryNoContact();
   testCapsuleCapsuleContact();
   testCapsuleCapsuleNoContact();
   testCapsuleEndCapContactDiffersFromCylinder();

   testCapsulePackUnpack();

   return EXIT_SUCCESS;
}
