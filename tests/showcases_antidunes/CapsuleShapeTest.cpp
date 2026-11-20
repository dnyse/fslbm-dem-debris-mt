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
//! Unit tests for mesa_pd::data::Capsule (a cylinder with two hemispherical end caps).
//!
//! The ContainsPoint cases guard the ContainsPointFunctor overload in src/mesa_pd/common/Contains.h; without it
//! the fluid-obstacle mapping (MovingParticleMapping.h, via SingleCast) aborts with "ContainsPoint not implemented!".
//! They distinguish rounded caps from flat caps, which an on-axis-only test would not.
//!
//! Coverage:
//!  - Volume/mass/inertia: analytic volume, and the height <= 2*radius limit where the capsule is an exact sphere.
//!  - support() (GJK): cardinal and oblique directions against the Minkowski-sum definition (segment + ball).
//!  - ContainsPointFunctor through SingleCast, the call path of MovingParticleMapping.h, including points
//!    beyond the core end but inside the rounded cap, and just outside the rounded corner.
//!  - Contact detection through the production class SolidCylinderContactDetection: Capsule x HalfSpace,
//!    Capsule x CylindricalBoundary and Capsule x Capsule. No analytic Capsule overload exists, so all three
//!    use the generic GJK/EPA fallback.
//!  - pack()/unpack() MPI round-trip.
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
   // height <= 2*radius -> coreHalfLength = 0 -> exactly a sphere of that radius (volume and inertia).
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
   // For r << h the capsule approaches a thin rod: Iperp -> m*Lc^2/12 (Lc = core length ~= h), Iaxis -> 0.
   // Not exact (O(r^2) correction), so checked to 1% as a limiting-behaviour check.
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

   // +z: segment support is the +z end of the core, ball support is r*(0,0,1): the cap tip.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 0, 1)), Vec3(0, 0, a + r), real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 0, -1)), Vec3(0, 0, -(a + r)), real_t(1e-10));

   // Purely radial: segment support is the centre (0,0,0) via the near-perpendicular guard (bfD[2] == 0),
   // ball support is r*(1,0,0), i.e. the equatorial band.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(1, 0, 0)), Vec3(r, 0, 0), real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 1, 0)), Vec3(0, r, 0), real_t(1e-10));

   // Oblique direction: must equal segment_support(d) + r*d (Minkowski-sum definition).
   Vec3 d(real_t(0.6), real_t(0), real_t(0.8)); // already unit length
   Vec3 expected(real_t(0), real_t(0), a) ; // sign(d.z)=+1
   expected = expected + d * r;
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(d), expected, real_t(1e-10));
}

void testCapsuleContainedWithinSameParameterCylinder()
{
   // Capsule uses Cylinder's height convention (full end-to-end extent), so Capsule(r,h) is contained in
   // Cylinder(r,h): the caps round off the rim. Support points coincide on-axis (both reach h/2) and radially
   // (both reach r); for an oblique direction the cylinder's sharp rim reaches strictly farther.
   const real_t r = real_t(2), h = real_t(10);
   const mesa_pd::data::Capsule  capsule(r, h);
   const mesa_pd::data::Cylinder cylinder(r, h);

   // Axial and radial directions give identical support points.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(0, 0, 1)), cylinder.support(Vec3(0, 0, 1)), real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(capsule.support(Vec3(1, 0, 0)), cylinder.support(Vec3(1, 0, 0)), real_t(1e-10));

   // Oblique direction: the cylinder rim reaches strictly farther along z, and the capsule's support point
   // stays within the cylinder's radius/height bounds.
   const Vec3 d(real_t(0.6), real_t(0), real_t(0.8));
   const Vec3 capSupport = capsule.support(d);
   const Vec3 cylSupport = cylinder.support(d);
   WALBERLA_CHECK(capSupport[2] < cylSupport[2] - real_t(1e-6));
   WALBERLA_CHECK(std::sqrt(capSupport[0] * capSupport[0] + capSupport[1] * capSupport[1]) <= r + real_t(1e-10));
   WALBERLA_CHECK(std::abs(capSupport[2]) <= h * real_t(0.5) + real_t(1e-10));
}

// ---------------------------------------------------------------------------------------------
// ContainsPointFunctor via SingleCast, the dispatch path MovingParticleMapping.h uses
// ---------------------------------------------------------------------------------------------

bool capsuleContains(std::shared_ptr< mesa_pd::data::ParticleStorage > ps, size_t idx, const Vec3& point)
{
   // SingleCast needs getShape() (raw BaseShape*); plain ParticleAccessor only has getBaseShape() (shared_ptr),
   // ParticleAccessorWithBaseShape adds it.
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
   // A point beyond the core segment's end (z > coreHalfLength) but within the hemispherical cap
   // (distance to the segment end <= radius) must be accepted.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   auto p  = ps->create();
   p->setPosition(Vec3(0, 0, 0));
   p->setRotation(mesa_pd::Rot3());
   const real_t r = real_t(2), h = real_t(10), a = real_t(3);
   p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(r, h));

   // z = a + 1.5 = 4.5 (outside a flat cylinder of height 2a = 6); radial offset such that the distance to the
   // segment end (0,0,a) is sqrt(1.5^2 + radial^2) = 1.9 < r = 2: inside the rounded cap.
   const real_t radial = std::sqrt(real_t(1.9) * real_t(1.9) - real_t(1.5) * real_t(1.5));
   WALBERLA_CHECK(capsuleContains(ps, 0, Vec3(radial, 0, a + real_t(1.5))));

   // Same axial offset, radial offset such that the distance to the segment end is 2.1 > r: just outside the
   // rounded cap (also outside the height-2a flat cylinder); rules out an always-true test.
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
   // Off-origin, rotated 90 deg about world X: checks that the test uses the particle's actual transform,
   // not the body frame at the origin.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   auto p  = ps->create();
   const Vec3 centre(real_t(5), real_t(7), real_t(-3));
   const mesa_pd::Rot3 rot(Vec3(1, 0, 0), math::pi * real_t(0.5));
   p->setPosition(centre);
   p->setRotation(rot);
   const real_t r = real_t(2), h = real_t(10), a = real_t(3);
   p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(r, h));

   // The cap tip's world position comes from the same rotation matrix, so the test does not depend on the
   // library's rotation sign convention.
   const Vec3 capTipWorld    = centre + rot.getMatrix() * Vec3(0, 0, a + r);
   const Vec3 outwardOneStep = (capTipWorld - centre).getNormalized() * real_t(1e-6);
   WALBERLA_CHECK(capsuleContains(ps, 0, capTipWorld - outwardOneStep));
   WALBERLA_CHECK(!capsuleContains(ps, 0, capTipWorld + outwardOneStep));
   // A point at distance r - eps along world X (unchanged by this rotation, perpendicular to the axis) is inside.
   WALBERLA_CHECK(capsuleContains(ps, 0, centre + Vec3(r - real_t(1e-6), 0, 0)));
}

// ---------------------------------------------------------------------------------------------
// Contact detection through the production class SolidCylinderContactDetection (ParticleKernels.h). No analytic
// Capsule overload exists, so these also cover the generic GeneralContactDetection fallback.
// ---------------------------------------------------------------------------------------------

void testCapsuleHalfSpaceContact()
{
   // Wall at x=0 with outward normal +x (as in ContactDetectionTest.cpp). Capsule along world z, radius 1,
   // band overlapping the wall by 0.3.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto wall = ps->create();
   wall->setPosition(Vec3(0, 0, 0));
   auto rod = ps->create();
   rod->setPosition(Vec3(real_t(0.7), 0, 0)); // penetration = radius - 0.7 = 0.3
   rod->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::HalfSpace wallShape(Vec3(1, 0, 0));
   const mesa_pd::data::Capsule   rodShape(real_t(1), real_t(6));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, wallShape, rodShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(scd.getPenetrationDepth()), real_t(0.3), real_t(1e-6));
   // The generic HalfSpace template of GeneralContactDetection returns a normal pointing into the wall, the
   // opposite convention from the analytic Sphere/Box/Cylinder x HalfSpace overloads. The applied force is
   // unaffected (LinearSpringDashpot is symmetric under index swap + negated normal); asserted so that a change
   // of either convention is noticed.
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
   // Pillar: CylindricalBoundary at the origin, radius 5, infinite along world z. Capsule of radius 1 alongside,
   // band overlapping the pillar surface by 0.2.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto pillar = ps->create();
   pillar->setPosition(Vec3(0, 0, 0));
   auto rod = ps->create();
   rod->setPosition(Vec3(real_t(5.8), 0, real_t(7))); // z=7 arbitrary: pillar is axially invariant
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
   // Two parallel capsules (radius 0.5, axis world z), centres 0.8 apart along x -> overlap 0.5+0.5-0.8 = 0.2.
   // Same geometry as testCylinderCylinderContact in ContactDetectionTest.cpp; the end caps play no role, so
   // the GJK/EPA fallback must match the cylinder result.
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
   // Two coaxial capsules end-on: the rounded tips touch like two spheres, giving a well-defined contact,
   // whereas the end-on contact of flat-capped cylinders is degenerate (face to face).
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
   // Looser epsilon: an exactly symmetric tip-to-tip contact is degenerate for the iterative EPA refinement,
   // which leaves some noise in the normal while the depth stays precise.
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
