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
//! \file ContactDetectionTest.cpp
//! Unit tests for walberla::antidunes::SolidCylinderContactDetection. Covers every floater shape (Sphere, Box,
//! Cylinder) against every floater, HalfSpace (wall/floor) and CylindricalBoundary (pillar): 12 pairs, since the two
//! static shapes never contact each other. Each pair has a contact case and a no-contact control; mixed-type pairs
//! other than Sphere x Box also check argument-order invariance.
//!
//! By resolution path:
//! - 6 custom overloads (hand-traced expected values): Sphere x Cylinder, Sphere x CylindricalBoundary,
//!   Box x HalfSpace, Box x CylindricalBoundary, Cylinder x CylindricalBoundary, Cylinder x HalfSpace. Cylinder x
//!   HalfSpace has its own overload (same support()-based formula as Box x HalfSpace) because the generic
//!   HalfSpace template of GeneralContactDetection returns contactNormal = -wallNormal, the opposite sign from
//!   the Sphere and Box pairs.
//! - 3 analytic overloads inherited from AnalyticContactDetection (Sphere x Sphere, Sphere x Box,
//!   Sphere x HalfSpace), tested through this class to confirm `using GeneralContactDetection::operator();` does
//!   not shadow them.
//! - 3 pairs on the GJK/EPA path: Box x Box, Cylinder x Box, Cylinder x Cylinder (Box x Box and Cylinder x Cylinder
//!   have closed forms only when useAnalyticalContactDetection is set, default false). Simple axis-aligned
//!   geometry whose result follows from the separating-axis theorem.
//!
//! Further sections: CylindricalBoundary containment mode, secondary contact points, and the cylinder/box contact
//! skins. testCylinderPillarBruteForceCrossCheck verifies the pillar closest-point search by dense surface sampling.
//======================================================================================================================

// SimulationObjects.h (included via ParticleKernels.h) uses blockforest::SimpleCommunication without including it;
// Antidunes.cpp gets it earlier through BuildTimeloop.h. Included explicitly so this file does not rely on that order.
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/Cylinder.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/mpi/MPIManager.h>

#include <cmath>
#include <limits>

namespace walberla {
namespace antidunes {

using Vec3 = mesa_pd::Vec3;

// ---------------------------------------------------------------------------------------------
// Analytic overloads inherited from AnalyticContactDetection (Sphere x Sphere, Sphere x Box), tested through
// SolidCylinderContactDetection to confirm `using GeneralContactDetection::operator();` does not shadow them.
// ---------------------------------------------------------------------------------------------

void testSphereSphereContact()
{
   // Geometry and expected values from checkSphereSphereCollision() in
   // tests/mesa_pd/collision_detection/AnalyticContactDetection.cpp. pybullet agrees to ~1e-6 on depth, normal and
   // contact point (for two spheres the contact point is unique: on the line between the centres).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto s0 = ps->create();
   s0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   const Vec3   dir   = Vec3(real_t(1), real_t(2), real_t(3)).getNormalized();
   const real_t shift = real_t(0.75);
   auto s1 = ps->create();
   s1->setPosition(dir * shift);

   const mesa_pd::data::Sphere sphereShape(real_t(0.5));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, sphereShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), shift - real_t(1));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), -dir);
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), dir * shift * real_t(0.5));
}

void testSphereSphereNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto s0 = ps->create();
   s0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   const Vec3 dir = Vec3(real_t(1), real_t(2), real_t(3)).getNormalized();
   auto s1 = ps->create();
   s1->setPosition(dir * real_t(1.5)); // separation: 1.5 - (0.5+0.5) = 0.5 > 0

   const mesa_pd::data::Sphere sphereShape(real_t(0.5));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, sphereShape, ac));
}

void testSphereBoxContact()
{
   // Hand-traced from analytic::detectSphereBoxCollision (AnalyticCollisionFunctions.h): sphere outside the box,
   // clamped projection point p=(0.5,0,0) (box half-extent), n=d-q=(0.5,0,0), dist=0.5-0.6=-0.1.
   // pybullet agrees on depth and normal; its contact point is the mid-surface point (0.45,0,0), whereas
   // detectSphereBoxCollision returns the box surface point (contactPoint=pos2+q), which is what is checked.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(1.0), real_t(0), real_t(0)));
   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   const mesa_pd::data::Box    boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.1));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0.5), real_t(0), real_t(0)));
}

void testSphereBoxNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(5), real_t(0), real_t(0))); // well clear
   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   const mesa_pd::data::Box    boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, boxShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Pairs resolved by the inherited GJK/EPA path (GeneralContactDetection::collideGJKEPA, margin=1e-4); with the
// default useAnalyticalContactDetection=false this includes Box x Box and Cylinder x Cylinder. No mesa_pd test
// (e.g. tests/mesa_pd/collision_detection/EPA.cpp) covers these three pairs.
//
// Depth and normal match pybullet to ~1e-6 for the same geometry, and the separating-axis prediction for these
// single-axis overlaps.
//
// Only the overlap-axis component of contactPoint is asserted. For convex polytopes the contact point is not
// unique: identical aligned cubes share a whole face of valid points, and EPA and pybullet pick different ones.
// pybullet also approximates cylinders by a discretised hull, so its off-axis components are not a reference.
// ---------------------------------------------------------------------------------------------

void testBoxBoxContact()
{
   // Two unit cubes (half-extent 0.5), centres 0.8 apart along x, y/z aligned -> overlap = 0.5+0.5-0.8 = 0.2 along x.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto b0 = ps->create();
   b0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   b0->setRotation(mesa_pd::Rot3());
   auto b1 = ps->create();
   b1->setPosition(Vec3(real_t(0.8), real_t(0), real_t(0)));
   b1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Box boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   // Only x is asserted: any (y,z) in [-0.5,0.5]^2 on the overlap face is valid; x=0.4 is the overlap midpoint.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint()[0], real_t(0.4), real_t(1e-3));
}

void testBoxBoxNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto b0 = ps->create();
   b0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   b0->setRotation(mesa_pd::Rot3());
   auto b1 = ps->create();
   b1->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // well clear: 0.5+0.5 < 3
   b1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Box boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, boxShape, boxShape, ac));
}

void testCylinderBoxContact()
{
   // Cylinder radius 0.5 (axis world z) vs. unit cube (half-extent 0.5), centres 0.9 apart along x
   // -> overlap = 0.5+0.5-0.9 = 0.1 along x (z-extents overlap: cylinder [-1,1], box [-0.5,0.5]).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3());
   auto box = ps->create();
   box->setPosition(Vec3(real_t(0.9), real_t(0), real_t(0)));
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Cylinder cylShape(real_t(0.5), real_t(2));
   const mesa_pd::data::Box      boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, cylShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.1), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   // Only x asserted, see the section header.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint()[0], real_t(0.45), real_t(1e-3));
}

void testCylinderBoxNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3());
   auto box = ps->create();
   box->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // well clear: 0.5+0.5 < 3
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Cylinder cylShape(real_t(0.5), real_t(2));
   const mesa_pd::data::Box      boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, cylShape, boxShape, ac));
}

void testCylinderBoxSwappedOrder()
{
   // GeneralContactDetection's GJK/EPA path reorders by particle UID before computing (its "idx2's UID > idx1's
   // UID" check), so both argument orders must give exactly the same result. Box x Box and Cylinder x Cylinder
   // have identical shapes on both sides, so this mixed pair is the one that exercises the reordering.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3());
   auto box = ps->create();
   box->setPosition(Vec3(real_t(0.9), real_t(0), real_t(0)));
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Cylinder cylShape(real_t(0.5), real_t(2));
   const mesa_pd::data::Box      boxShape(Vec3(real_t(1), real_t(1), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, cylShape, boxShape, ac));

   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, boxShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testCylinderCylinderContact()
{
   // Two parallel cylinders (radius 0.5, axis world z), centres 0.8 apart along x -> overlap 0.5+0.5-0.8 = 0.2.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto c0 = ps->create();
   c0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   c0->setRotation(mesa_pd::Rot3());
   auto c1 = ps->create();
   c1->setPosition(Vec3(real_t(0.8), real_t(0), real_t(0)));
   c1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Cylinder cylShape(real_t(0.5), real_t(2));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, cylShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   // Only x asserted, see the section header.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint()[0], real_t(0.4), real_t(1e-3));
}

void testCylinderCylinderNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto c0 = ps->create();
   c0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   c0->setRotation(mesa_pd::Rot3());
   auto c1 = ps->create();
   c1->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // well clear: 0.5+0.5 < 3
   c1->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Cylinder cylShape(real_t(0.5), real_t(2));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, cylShape, cylShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Sphere (idx0) x Cylinder (idx1). For the side-contact geometry (r=1 sphere at (1.5,0,0), r=1 h=4 cylinder at the
// origin) pybullet agrees to ~1e-4 on depth, normal and contact point: both use the mid-overlap contact point.
// ---------------------------------------------------------------------------------------------

void testSphereCylinderSideContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(1.5), real_t(0), real_t(0)));
   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3()); // identity

   const mesa_pd::data::Sphere   sphereShape(real_t(1));
   const mesa_pd::data::Cylinder cylShape(real_t(1), real_t(4)); // radius 1, height 4 (halfLen 2)

   SolidCylinderContactDetection scd;
   const bool contact = scd(0, 1, sphereShape, cylShape, ac);

   WALBERLA_CHECK(contact);
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.5));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0.75), real_t(0), real_t(0)));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, cylShape, sphereShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testSphereCylinderRotationAware()
{
   // Cylinder rotated 90 deg about its own z axis: the axis is unchanged, but the radial normal and point must be
   // rotated (R applied to nBody), which exercises the body-frame transform.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(0), real_t(1.5), real_t(0))); // = R * (1.5,0,0)
   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3(Vec3(real_t(0), real_t(0), real_t(1)), math::pi * real_t(0.5)));

   const mesa_pd::data::Sphere   sphereShape(real_t(1));
   const mesa_pd::data::Cylinder cylShape(real_t(1), real_t(4));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.5));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(1), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0), real_t(0.75), real_t(0)));
}

void testSphereInsideCylinderEnvelopBranch()
{
   // Sphere centre inside the cylinder's radial extent, deep overlap -> pushed out radially
   // because dSide < dCap in this configuration.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(0.2), real_t(0), real_t(0)));
   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Sphere   sphereShape(real_t(1.5));
   const mesa_pd::data::Cylinder cylShape(real_t(1), real_t(4));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-2.3));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(-0.15), real_t(0), real_t(0)));
}

void testSphereCylinderNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // separation: 3-(1+1)=1 > 0
   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   cyl->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Sphere   sphereShape(real_t(1));
   const mesa_pd::data::Cylinder cylShape(real_t(1), real_t(4));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, cylShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Sphere (idx0, "sphere floater") x CylindricalBoundary (idx1, "pillar")
// ---------------------------------------------------------------------------------------------

void testSphereCylindricalBoundaryContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(1.5), real_t(0), real_t(7))); // z=7 arbitrary, see the axial-invariance note
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere              sphereShape(real_t(0.6));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, pillarShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.1));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   // contactPoint z = 7 (the sphere's z): the radial projection (delta - (delta.axis)*axis) discards the axial
   // component instead of biasing toward the pillar's z=0.
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(1.0), real_t(0), real_t(7)));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, pillarShape, sphereShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testSphereCylindricalBoundaryNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // separation = 3 - 1 - 0.6 = 1.4 > 0
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere              sphereShape(real_t(0.6));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, pillarShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Sphere (idx0, "sphere floater") x HalfSpace (idx1, "wall/floor"): inherited from AnalyticContactDetection
// (detectSphereHalfSpaceCollision), tested here to confirm the custom overloads do not shadow it.
// ---------------------------------------------------------------------------------------------

void testSphereHalfSpaceContact()
{
   // Same wall as testBoxHalfSpaceContactAxisAligned. Expected values hand-traced from detectSphereHalfSpaceCollision
   // (AnalyticCollisionFunctions.h): d = n.pos2 = 0, k = n.pos1 = 0.3,
   // depth = k - r - d = 0.3 - 0.6 = -0.3, contactPoint = pos1 - (r+depth)*n = (0,0,0.3) - 0.3*(0,0,1)
   // = (0,0,0), contactNormal = n (unchanged sign, "points away from the wall").
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(0), real_t(0), real_t(0.3)));
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0))); // plane at z=0

   const mesa_pd::data::Sphere    sphereShape(real_t(0.6));
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1))); // normal points up, away from the wall

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, wallShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.3));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0), real_t(0), real_t(0)));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, wallShape, sphereShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testSphereHalfSpaceNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(0), real_t(0), real_t(2))); // well clear of the wall
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere    sphereShape(real_t(0.6));
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, wallShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Box (idx0, "box floater") x HalfSpace (idx1, "wall/floor")
// ---------------------------------------------------------------------------------------------

void testBoxHalfSpaceContactAxisAligned()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(0.3)));
   box->setRotation(mesa_pd::Rot3()); // identity
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0))); // plane at z=0

   const mesa_pd::data::Box       boxShape(Vec3(real_t(1), real_t(1), real_t(1))); // edge length -> half-extent 0.5
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1))); // normal points up, away from the wall

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, wallShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.2));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)));
   // Flush face-on: the contact region is the whole bottom face; Box::support() returns zero offset on axes with a
   // zero direction component, so the contact point is the face centre (0,0,0), not an arbitrary corner.
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0), real_t(0), real_t(0)));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, wallShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testBoxHalfSpaceContactRotated()
{
   // Rotated about X, not Z: a rotation about the wall normal leaves the body-frame wall direction unchanged and
   // gives the same face-centre point as the unrotated case. Tilting about X puts the box edge-first against the
   // wall and checks that R is applied when finding the deepest point.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   const real_t tilt = math::pi / real_t(6); // 30 deg about X

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(0.3)));
   box->setRotation(mesa_pd::Rot3(Vec3(real_t(1), real_t(0), real_t(0)), tilt));
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box       boxShape(Vec3(real_t(0.6), real_t(1.0), real_t(1.0))); // half-extents (0.3,0.5,0.5)
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, wallShape, ac));

   // Body-frame direction toward the wall: R^T*(0,0,-1) = (0, -sin(tilt), -cos(tilt)). The X component is 0, so the
   // deepest feature is an edge along X and support() returns its midpoint (X offset 0); in the Y-Z section it is
   // the corner: body-frame support point (0, -0.5, -0.5).
   const real_t s = std::sin(tilt), c = std::cos(tilt);
   const Vec3   worldOffset(real_t(0), real_t(0.5) * (s - c), real_t(-0.5) * (s + c)); // R * (0,-0.5,-0.5)
   const Vec3   deepest            = box->getPosition() + worldOffset;
   const real_t expectedPenetration = deepest[2]; // wall plane through the origin, normal (0,0,1)
   const Vec3   expectedContactPoint(deepest[0], deepest[1], real_t(0));

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), expectedPenetration, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint(), expectedContactPoint, real_t(1e-6));
}

void testBoxHalfSpaceNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(2))); // well clear of the wall
   box->setRotation(mesa_pd::Rot3());
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box       boxShape(Vec3(real_t(1), real_t(1), real_t(1)));
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, boxShape, wallShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Box (idx0, "box floater") x CylindricalBoundary (idx1, "pillar")
// ---------------------------------------------------------------------------------------------

void testBoxCylindricalBoundaryContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(1.2), real_t(0), real_t(0)));
   box->setRotation(mesa_pd::Rot3()); // identity
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box                 boxShape(Vec3(real_t(1), real_t(1), real_t(4))); // half-extents (0.5,0.5,2.0)
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, pillarShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.3));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0.7), real_t(0), real_t(0)));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, pillarShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testBoxCylindricalBoundaryAxialInvariance()
{
   // As testBoxCylindricalBoundaryContact, box shifted by z=5 along the pillar axis: depth and normal unchanged,
   // contact point z follows the box (not the pillar's z=0).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(1.2), real_t(0), real_t(5)));
   box->setRotation(mesa_pd::Rot3());
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box                 boxShape(Vec3(real_t(1), real_t(1), real_t(4)));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, pillarShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.3));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0.7), real_t(0), real_t(5)));
}

void testBoxCylindricalBoundaryNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // well clear of the pillar
   box->setRotation(mesa_pd::Rot3());
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box                 boxShape(Vec3(real_t(1), real_t(1), real_t(4)));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, boxShape, pillarShape, ac));
}

// ---------------------------------------------------------------------------------------------
// Cylinder (idx0, "floating log") x CylindricalBoundary (idx1, "pillar")
// ---------------------------------------------------------------------------------------------

void testCylinderPillarContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto log = ps->create();
   log->setPosition(Vec3(real_t(1.5), real_t(0), real_t(0)));
   log->setRotation(mesa_pd::Rot3()); // identity: log axis = world z
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Cylinder            logShape(real_t(1), real_t(4));         // radius 1, height 4
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1))); // radius 1, axis = world z

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, logShape, pillarShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.5));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0.5), real_t(0), real_t(0)));

   // Swap order.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, pillarShape, logShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testCylinderPillarNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto log = ps->create();
   log->setPosition(Vec3(real_t(3), real_t(0), real_t(0))); // far enough apart: separation = 2 - 1 = 1 > 0
   log->setRotation(mesa_pd::Rot3());
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Cylinder            logShape(real_t(1), real_t(4));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, logShape, pillarShape, ac));
}

void testCylinderPillarBruteForceCrossCheck()
{
   // Independent check of the closest-point search: densely sample the log's curved surface, take the minimum
   // perpendicular distance to the pillar axis, and compare with the implementation (geometry of
   // testCylinderPillarContact).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   const Vec3 logPos(real_t(1.5), real_t(0), real_t(0));
   const Vec3 pillarPos(real_t(0), real_t(0), real_t(0));
   const Vec3 pillarAxis(real_t(0), real_t(0), real_t(1));
   const real_t logRadius = real_t(1), logHalfLen = real_t(2), pillarRadius = real_t(1);

   auto log = ps->create();
   log->setPosition(logPos);
   log->setRotation(mesa_pd::Rot3()); // identity: log axis = world z
   auto pillar = ps->create();
   pillar->setPosition(pillarPos);

   const mesa_pd::data::Cylinder            logShape(logRadius, real_t(2) * logHalfLen);
   const mesa_pd::data::CylindricalBoundary pillarShape(pillarRadius, pillarAxis);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, logShape, pillarShape, ac));
   const real_t analyticDistToAxis = scd.getPenetrationDepth() + pillarRadius; // undo the "- pillarRadius" step

   // Sampling of the curved surface (identity rotation, body == world): theta in [0,2*pi), z in [-halfLen,halfLen].
   constexpr uint_t nTheta = 720;
   constexpr uint_t nZ     = 400;
   real_t minDist = std::numeric_limits< real_t >::max();
   for (uint_t i = 0; i < nTheta; ++i)
   {
      const real_t theta = real_t(2) * math::pi * real_t(i) / real_t(nTheta);
      const real_t px = logPos[0] + logRadius * std::cos(theta);
      const real_t py = logPos[1] + logRadius * std::sin(theta);
      for (uint_t j = 0; j < nZ; ++j)
      {
         const real_t z  = logPos[2] - logHalfLen + real_t(2) * logHalfLen * real_t(j) / real_t(nZ - 1);
         // Perpendicular distance from (px,py,z) to the pillar axis line (through pillarPos, direction pillarAxis).
         const Vec3   rel = Vec3(px, py, z) - pillarPos;
         const real_t t   = rel * pillarAxis;
         const Vec3   radial = rel - t * pillarAxis;
         minDist = std::min(minDist, radial.length());
      }
   }

   // Tolerance set by sampling: 720 angles on a radius-1 circle (0.5 deg spacing) give a sagitta well under 1e-4.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(minDist, analyticDistToAxis, real_t(1e-3));
}

// ---------------------------------------------------------------------------------------------
// Cylinder (idx0, "floating log") x HalfSpace (idx1, "wall/floor"): own overload in SolidCylinderContactDetection
// (ParticleKernels.h), same support()-based formula as Box x HalfSpace, normal = +wallNormal as for Sphere and Box
// (the generic HalfSpace template would give -wallNormal). The flat end cap faces the wall; Cylinder::support()
// returns the cap centre for a purely axial direction, so the case is well defined.
// ---------------------------------------------------------------------------------------------

void testCylinderHalfSpaceContact()
{
   // Cylinder axis = world z (identity rotation), radius 1, height 4 (halfLen 2), centred at
   // (0,0,1) -> bottom cap centre at world (0,0,-1), 1 LU below the wall plane at z=0.
   // Hand-traced from the overload:
   // deepest = c + R*cyl.support(-(R^T*n)) = (0,0,1) + cyl.support((0,0,-1));
   // Cylinder::support((0,0,-1)) in body frame = (0,0,-2) (cap centre)
   // -> deepest = (0,0,-1); penetrationDepth = n.deepest - d = -1-0 = -1;
   // contactNormal = n = (0,0,1); contactPoint = deepest - penetrationDepth*n = (0,0,-1)+(0,0,1)
   // = (0,0,0) (on the wall plane, same convention as the Box case).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(1)));
   cyl->setRotation(mesa_pd::Rot3()); // identity: axis = world z
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0))); // plane at z=0

   const mesa_pd::data::Cylinder  cylShape(real_t(1), real_t(4)); // radius 1, height 4
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1))); // normal points up, away from the wall

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, cylShape, wallShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-1));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0), real_t(0), real_t(0)));

   // Swapped-argument overload must give the identical result (the HalfSpace x Cylinder overload delegates to
   // Cylinder x HalfSpace with the indices swapped).
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, wallShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testCylinderHalfSpaceNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(3))); // bottom cap at z=1, well clear of z=0
   cyl->setRotation(mesa_pd::Rot3());
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Cylinder  cylShape(real_t(1), real_t(4));
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, cylShape, wallShape, ac));
}

// ---------------------------------------------------------------------------------------------
// CylindricalBoundary containment mode (cb.setContainment(true)), the mirror of the pillar case: a body poking
// outside the tube is pushed back toward the axis. Covers the farthest-point search (Sphere/Box/Cylinder), the
// optional axial bounds and the optional half-clip plane (e.g. an arch solid only above its springing line).
// isContainment() defaults to false, so the tests above are unaffected.
// ---------------------------------------------------------------------------------------------

void testSphereCylindricalBoundaryContainmentContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(2), real_t(4.7), real_t(0)));
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   const real_t expectedRadial = real_t(4.7) + real_t(0.6); // sphere is isotropic: exact, no basis needed
   const real_t expectedOverlap = expectedRadial - real_t(5);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, tubeShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), -expectedOverlap, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(0), real_t(-1), real_t(0)), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint(), Vec3(real_t(2), expectedRadial, real_t(0)), real_t(1e-6));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, tubeShape, sphereShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testSphereCylindricalBoundaryContainmentNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(2), real_t(3), real_t(0))); // farthest reach 3+0.6=3.6 < radius 5
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, tubeShape, ac));
}

void testSphereCylindricalBoundaryContainmentAxialBounds()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   // Radially the sphere pokes outside (farthest reach 4.7+0.6=5.3 > radius 5), but at X=10 it is outside the
   // tube's axial span [-3,3] (offsets from the tube position X=0): the axial gate must reject it.
   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(10), real_t(4.7), real_t(0)));
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);
   tubeShape.setAxialBounds(real_t(-3), real_t(3));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, tubeShape, ac));
}

void testSphereCylindricalBoundaryContainmentHalfClipExcludesBelow()
{
   // Sphere entirely on the non-solid side of the half-clip plane (Z < 0, e.g. below an arch's springing line):
   // although it pokes outside the full tube radius, there is no wall there, so no contact.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(2), real_t(0), real_t(-4.7))); // top of sphere at z=-4.1, still < 0
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);
   tubeShape.setHalfClipNormal(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, tubeShape, ac));
}

void testSphereCylindricalBoundaryContainmentHalfClipStraddling()
{
   // Sphere straddling the clip plane (Z=0): the contact point must stay on the solid side. Checked via invariants
   // (contact, point on the solid side, genuine radial overlap); the exact arc-endpoint position depends on the
   // quartic solver, which is tested elsewhere.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(2), real_t(4.9), real_t(0.3)));
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere sphereShape(real_t(0.6));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);
   tubeShape.setHalfClipNormal(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, tubeShape, ac));
   WALBERLA_CHECK(scd.getContactPoint()[2] >= real_t(-1e-9)); // stayed on the solid side of the clip
   const real_t radialDist = std::sqrt(scd.getContactPoint()[1] * scd.getContactPoint()[1] +
                                        scd.getContactPoint()[2] * scd.getContactPoint()[2]);
   WALBERLA_CHECK(radialDist > real_t(5)); // genuine overlap, not a clamped-to-boundary non-contact
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), -(radialDist - real_t(5)), real_t(1e-6));
}

void testBoxCylindricalBoundaryContainmentContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(4.5), real_t(0.3))); // z-shifted to avoid a Z-symmetric tie
   box->setRotation(mesa_pd::Rot3());
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box boxShape(Vec3(real_t(2), real_t(2), real_t(2))); // half-extent 1 each way
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   // Farthest corner in Y-Z: (4.5+1, 0.3+1).
   const real_t expectedY = real_t(5.5), expectedZ = real_t(1.3);
   const real_t expectedRadial = std::sqrt(expectedY * expectedY + expectedZ * expectedZ);
   const real_t expectedOverlap = expectedRadial - real_t(5);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, tubeShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), -expectedOverlap, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint(), Vec3(real_t(-1), expectedY, expectedZ), real_t(1e-6));
}

void testBoxCylindricalBoundaryContainmentNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(3), real_t(0))); // farthest corner reach sqrt(4^2+1^2) < 5
   box->setRotation(mesa_pd::Rot3());
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box boxShape(Vec3(real_t(2), real_t(2), real_t(2)));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, boxShape, tubeShape, ac));
}

void testBoxCylindricalBoundaryContainmentRadialPreGate()
{
   // Several tubes can share one axial span (e.g. adjacent arches). A box inside the span but far from this tube's
   // axis must be rejected by the radial pre-gate; the farthest-corner search alone would report a large spurious
   // overlap. Axial span X in [-3,3]; the box is 30 LU away in Y.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(1), real_t(30), real_t(0))); // inside the axial span, far away radially
   box->setRotation(mesa_pd::Rot3());
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box boxShape(Vec3(real_t(2), real_t(2), real_t(2))); // half-extent 1 each way
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);
   tubeShape.setAxialBounds(real_t(-3), real_t(3));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, boxShape, tubeShape, ac));
}

void testCylinderCylindricalBoundaryContainmentContact()
{
   // Generic tilt on purpose: a cylinder exactly along a world axis perpendicular to the tube axis is a degenerate
   // quartic in farthestThetaOnCircleToLine (ClosedFormContactsCylinderPlaneGeometry.h; c4, c2, c0 vanish together
   // and the Durand-Kerner solve may miss the root). Validated via invariants, not a closed form.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(4.7), real_t(0.3)));
   cyl->setRotation(mesa_pd::Rot3(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.17)));
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Cylinder cylShape(real_t(0.4), real_t(2)); // radius 0.4, height 2 (hh=1)
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, cylShape, tubeShape, ac));

   const auto&  cp = scd.getContactPoint();
   const real_t radialDist = std::sqrt(cp[1] * cp[1] + cp[2] * cp[2]);
   WALBERLA_CHECK(radialDist > real_t(5)); // genuine overlap, not a boundary/no-contact case
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), -(radialDist - real_t(5)), real_t(1e-6));
   // The point must lie on the floater's rim: cylRadius from its axis line and at one end (+-hh), never at an
   // interior axial offset (the farthest point is always on an end rim).
   const mesa_pd::Vec3 logAxis = cyl->getRotation().getMatrix() * mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
   const mesa_pd::Vec3 rel     = cp - cyl->getPosition();
   const real_t         axialT  = rel * logAxis;
   const mesa_pd::Vec3 radialOnCyl = rel - axialT * logAxis;
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(radialOnCyl.length(), real_t(0.4), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(axialT), real_t(1), real_t(1e-6));
   // Normal must point inward (back toward the tube axis).
   WALBERLA_CHECK(scd.getContactNormal() * mesa_pd::Vec3(real_t(0), cp[1], cp[2]) < real_t(0));
}

void testCylinderCylindricalBoundaryContainmentNoContact()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(3), real_t(0))); // farthest reach 3+0.4=3.4 < radius 5
   cyl->setRotation(mesa_pd::Rot3());
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Cylinder cylShape(real_t(0.4), real_t(2));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, cylShape, tubeShape, ac));
}

// Vertical round container (DropTest ContainerCylinders, Utility.h createContainerCylinder): containment tube with
// axis z and radius 25.3 (Kodam et al. 2010b: 50.6 mm container), Kodam's cylinders (radius 4, length 5.3). Three
// orientations: disc lying flat (axis parallel to the tube axis), disc upright (axis exactly perpendicular, the
// degenerate-quartic case above), and generic tilts. Expected overlaps are hand-derived (flat, upright) and always
// checked against a dense sampling of both rims.
real_t containerSampledOverlap(const Vec3& pos, const mesa_pd::Rot3& rot, real_t rc, real_t hh, real_t R)
{
   // farthest radial distance (from the z axis through the origin) over a dense sampling of both rims
   const mesa_pd::Vec3 a = rot.getMatrix() * Vec3(0, 0, 1);
   mesa_pd::Vec3       u, v;
   closed_form_contacts::detail::perpendicularBasis(a, u, v);
   real_t best = 0;
   for (const real_t sgn : { real_t(-1), real_t(1) })
      for (int k = 0; k < 200000; ++k)
      {
         const real_t th = real_t(2) * math::pi * real_t(k) / real_t(200000);
         const Vec3   p  = pos + sgn * hh * a + rc * (std::cos(th) * u + std::sin(th) * v);
         best            = std::max(best, std::sqrt(p[0] * p[0] + p[1] * p[1]));
      }
   return best - R;
}

void checkVerticalContainer(const char* name, const Vec3& pos, const mesa_pd::Rot3& rot, real_t expectedOverlap)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto cyl = ps->create();
   cyl->setPosition(pos);
   cyl->setRotation(rot);
   cyl->setInteractionRadius(std::sqrt(real_t(4) * real_t(4) + real_t(2.65) * real_t(2.65)));
   auto tube = ps->create();
   tube->setPosition(Vec3(0, 0, 0));
   const mesa_pd::data::Cylinder      cylShape(real_t(4), real_t(5.3));
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(25.3), Vec3(0, 0, 1));
   tubeShape.setContainment(true);
   SolidCylinderContactDetection scd;
   const bool   hit     = scd(0, 1, cylShape, tubeShape, ac);
   const real_t sampled = containerSampledOverlap(pos, rot, real_t(4), real_t(2.65), real_t(25.3));
   WALBERLA_LOG_INFO_ON_ROOT("[vertical container] " << name << ": contact=" << hit << " overlap="
                             << (hit ? -scd.getPenetrationDepth() : real_t(0)) << " expected=" << expectedOverlap
                             << " sampled=" << sampled);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sampled, expectedOverlap, real_t(1e-6), name);
   if (expectedOverlap <= real_t(0))
   {
      WALBERLA_CHECK(!hit, name);
      return;
   }
   WALBERLA_CHECK(hit, name);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(-scd.getPenetrationDepth(), expectedOverlap, real_t(1e-6), name);
   const auto& cp = scd.getContactPoint();
   const auto& n  = scd.getContactNormal();
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::sqrt(cp[0] * cp[0] + cp[1] * cp[1]), real_t(25.3) + expectedOverlap,
                                      real_t(1e-6), name);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[2], real_t(0), real_t(1e-9), name); // horizontal
   WALBERLA_CHECK_LESS(n * Vec3(cp[0], cp[1], 0), real_t(0), name);       // inward
}

void testCylinderVerticalContainer()
{
   const real_t R = real_t(25.3), rc = real_t(4), hh = real_t(2.65), ov = real_t(0.05);
   // flat disc (axis || z), centre at distance R - rc + ov from the axis along x and along a diagonal
   checkVerticalContainer("flat, offset +x", Vec3(R - rc + ov, 0, 10), mesa_pd::Rot3(), ov);
   checkVerticalContainer("flat, offset diagonal", Vec3((R - rc + ov) / std::sqrt(real_t(2)), (R - rc + ov) / std::sqrt(real_t(2)), 10),
                          mesa_pd::Rot3(), ov);
   checkVerticalContainer("flat, inside", Vec3(R - rc - real_t(0.5), 0, 10), mesa_pd::Rot3(), real_t(-0.5));
   // upright disc, axis exactly along y (perpendicular to the tube axis), centre on the x axis: its farthest points
   // are the two rim points at (x0 + rc, +-hh): radial distance sqrt((x0 + rc)^2 + hh^2)
   {
      const real_t x0 = std::sqrt((R + ov) * (R + ov) - hh * hh) - rc;
      checkVerticalContainer("upright, axis y, centre on x", Vec3(x0, 0, 10),
                             mesa_pd::Rot3(Vec3(1, 0, 0), math::pi * real_t(0.5)), ov);
   }
   // upright disc, axis exactly along x, centre on the x axis: in plan view the outer face is the segment from
   // (x0 + hh, -rc) to (x0 + hh, rc), farthest radial distance sqrt((x0 + hh)^2 + rc^2)
   {
      const real_t x0 = std::sqrt((R + ov) * (R + ov) - rc * rc) - hh;
      checkVerticalContainer("upright, axis x, centre on x", Vec3(x0, 0, 10),
                             mesa_pd::Rot3(Vec3(0, 1, 0), math::pi * real_t(0.5)), ov);
   }
   // generic tilts: the centre is moved outward along a fixed direction until the sampled overlap is ov (bisection)
   for (const real_t angle : { real_t(0.3), real_t(1.1), real_t(2.2) })
   {
      const mesa_pd::Rot3 rot(Vec3(real_t(0.3), real_t(-0.7), real_t(0.2)).getNormalized(), angle);
      const Vec3          dir = Vec3(real_t(0.89), real_t(0.46), 0).getNormalized();
      real_t              lo = 0, hi = R;
      for (int it = 0; it < 60; ++it)
      {
         const real_t m = real_t(0.5) * (lo + hi);
         (containerSampledOverlap(dir * m + Vec3(0, 0, 10), rot, rc, hh, R) < ov ? lo : hi) = m;
      }
      const Vec3 pos = dir * lo + Vec3(0, 0, 10);
      checkVerticalContainer("generic tilt", pos, rot, containerSampledOverlap(pos, rot, rc, hh, R));
   }
}

// ---------------------------------------------------------------------------------------------
// Secondary contact points (containment, Box/Cylinder floaters only; a sphere has a single point). The
// farthest-point search tracks the top two points and exposes a second breaching point via
// hasSecondaryContact()/getSecondaryContact*(), so a body touching the curve at two points gets a force at both.
// ---------------------------------------------------------------------------------------------

void testBoxCylindricalBoundaryContainmentDualCornerContact()
{
   // Geometry of testBoxCylindricalBoundaryContainmentContact. For an unrotated box the two farthest corners differ
   // only in local X (same Y,Z, same radial distance), so both ends of the edge breach at once. The loop order
   // (sx outer, -1 first) makes X=-1 primary and X=+1 the tied secondary.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(4.5), real_t(0.3)));
   box->setRotation(mesa_pd::Rot3());
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box boxShape(Vec3(real_t(2), real_t(2), real_t(2))); // half-extent 1 each way
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, tubeShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint(), Vec3(real_t(-1), real_t(5.5), real_t(1.3)), real_t(1e-6));

   WALBERLA_CHECK(scd.hasSecondaryContact());
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getSecondaryContactPoint(), Vec3(real_t(1), real_t(5.5), real_t(1.3)), real_t(1e-6));
   // Tied corners -> identical overlap and normal, just the other end of the same edge.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getSecondaryPenetrationDepth(), scd.getPenetrationDepth(), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getSecondaryContactNormal(), scd.getContactNormal(), real_t(1e-6));
}

void testBoxCylindricalBoundaryContainmentSingleCornerNoSecondary()
{
   // Rotated on purpose: an unrotated box always ties its farthest corner pair (see the dual-corner test), so the
   // single-point case needs a generic tilt (same as testCylinderCylindricalBoundaryContainmentContact).
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   // Pose found by corner enumeration with margin on both sides: farthest corner at radial distance ~5.208
   // (overlap ~0.208), runner-up ~4.814 (~0.186 inside the radius).
   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(3.5), real_t(0.6)));
   box->setRotation(mesa_pd::Rot3(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.17)));
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box boxShape(Vec3(real_t(2.4), real_t(2.4), real_t(2.4))); // half-extent 1.2
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, tubeShape, ac));
   WALBERLA_CHECK(!scd.hasSecondaryContact());
}

void testCylinderCylindricalBoundaryContainmentDualRimContact()
{
   // A tilted log across the tube with both end rims past the tube radius. Validated via invariants (point on the
   // log's rim, genuine overlap at both ends), as in testCylinderCylindricalBoundaryContainmentContact.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(4.6), real_t(0.2)));
   cyl->setRotation(mesa_pd::Rot3(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.17)));
   auto tube = ps->create();
   tube->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Cylinder cylShape(real_t(0.4), real_t(6)); // radius 0.4, height 6 (hh=3)
   mesa_pd::data::CylindricalBoundary tubeShape(real_t(5), Vec3(real_t(1), real_t(0), real_t(0)));
   tubeShape.setContainment(true);

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, cylShape, tubeShape, ac));

   const mesa_pd::Vec3 logAxis = cyl->getRotation().getMatrix() * mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
   auto checkOwnRim = [&](const mesa_pd::Vec3& cp, real_t expectedPenetration, const mesa_pd::Vec3& normal) {
      const real_t         radialDist = std::sqrt(cp[1] * cp[1] + cp[2] * cp[2]);
      WALBERLA_CHECK(radialDist > real_t(5));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(expectedPenetration, -(radialDist - real_t(5)), real_t(1e-6));
      const mesa_pd::Vec3 rel        = cp - cyl->getPosition();
      const real_t         axialT     = rel * logAxis;
      const mesa_pd::Vec3 radialOnCyl = rel - axialT * logAxis;
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(radialOnCyl.length(), real_t(0.4), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(axialT), real_t(3), real_t(1e-6)); // hh=3, always at an end
      WALBERLA_CHECK(normal * mesa_pd::Vec3(real_t(0), cp[1], cp[2]) < real_t(0));   // inward
   };

   checkOwnRim(scd.getContactPoint(), scd.getPenetrationDepth(), scd.getContactNormal());
   WALBERLA_CHECK(scd.hasSecondaryContact());
   checkOwnRim(scd.getSecondaryContactPoint(), scd.getSecondaryPenetrationDepth(), scd.getSecondaryContactNormal());
   // The two points are on opposite ends of the log, not a duplicate.
   WALBERLA_CHECK((scd.getContactPoint() - scd.getSecondaryContactPoint()).length() > real_t(1));
}


} // namespace antidunes
} // namespace walberla

namespace walberla {
namespace antidunes {
using Vec3 = mesa_pd::Vec3;

// ---------------------------------------------------------------------------------------------
// cylinderContactSkin (ModelParameters): every floating cylinder is dilated by the skin m for
// contact detection. Expected values follow from the Minkowski dilation: the signed surface distance to another
// dilated rod shrinks by 2m, to a plane or a solid pillar by m; normals are unchanged.
// ---------------------------------------------------------------------------------------------
void testCylinderSkinCylinderCylinder()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto c0 = ps->create();
   c0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   c0->setRotation(mesa_pd::Rot3());
   auto c1 = ps->create();
   c1->setPosition(Vec3(real_t(1.3), real_t(0), real_t(0))); // parallel, surface gap 0.3
   c1->setRotation(mesa_pd::Rot3());
   const mesa_pd::data::Cylinder cylShape(real_t(0.5), real_t(2));

   // Without skin: separated.
   SolidCylinderContactDetection scd0;
   WALBERLA_CHECK(!scd0(0, 1, cylShape, cylShape, ac));

   // Skin 0.25: dilated gap 0.3 - 0.5 = -0.2 -> contact, normal along -x (same sign as testCylinderCylinderContact),
   // contact point in the middle of the dilated overlap, x = (0.5+0.25 + 1.3-0.5-0.25)/2 = 0.65.
   SolidCylinderContactDetection scd;
   scd.cylinderContactSkin = real_t(0.25);
   WALBERLA_CHECK(scd(0, 1, cylShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint()[0], real_t(0.65), real_t(1e-3));

   // Skin 0.1: dilated gap 0.3 - 0.2 = 0.1 > 0 -> still separated.
   SolidCylinderContactDetection scdSmall;
   scdSmall.cylinderContactSkin = real_t(0.1);
   WALBERLA_CHECK(!scdSmall(0, 1, cylShape, cylShape, ac));

   // Already overlapping cores (centres 0.8 apart, overlap 0.2): depth grows by 2m.
   c1->setPosition(Vec3(real_t(0.8), real_t(0), real_t(0)));
   SolidCylinderContactDetection scdOv;
   scdOv.cylinderContactSkin = real_t(0.25);
   WALBERLA_CHECK(scdOv(0, 1, cylShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scdOv.getPenetrationDepth(), real_t(-0.7), real_t(1e-3));
}

void testCylinderSkinHalfSpaceAndPillar()
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   // Wall: cylinder (radius 1, height 4, axis z) with its bottom cap 0.3 above the plane z = 0.
   auto cyl = ps->create();
   cyl->setPosition(Vec3(real_t(0), real_t(0), real_t(2.3)));
   cyl->setRotation(mesa_pd::Rot3());
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   const mesa_pd::data::Cylinder  cylShape(real_t(1), real_t(4));
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd0;
   WALBERLA_CHECK(!scd0(0, 1, cylShape, wallShape, ac));
   SolidCylinderContactDetection scd;
   scd.cylinderContactSkin = real_t(0.5);
   WALBERLA_CHECK(scd(0, 1, cylShape, wallShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.2));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0), real_t(0), real_t(0))); // stays on the plane
   SolidCylinderContactDetection scdSw;
   scdSw.cylinderContactSkin = real_t(0.5);
   WALBERLA_CHECK(scdSw(1, 0, wallShape, cylShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scdSw.getPenetrationDepth(), real_t(-0.2));

   // Pillar: log (radius 1, axis z) at x = 2.3 next to a solid pillar of radius 1 at the origin, gap 0.3.
   auto ps2 = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac2(ps2);
   auto log = ps2->create();
   log->setPosition(Vec3(real_t(2.3), real_t(0), real_t(0)));
   log->setRotation(mesa_pd::Rot3());
   auto pillar = ps2->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scdP0;
   WALBERLA_CHECK(!scdP0(0, 1, cylShape, pillarShape, ac2));
   SolidCylinderContactDetection scdP;
   scdP.cylinderContactSkin = real_t(0.5);
   WALBERLA_CHECK(scdP(0, 1, cylShape, pillarShape, ac2));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scdP.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scdP.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)), real_t(1e-6));
   // Threshold is restored after the call.
   WALBERLA_CHECK_FLOAT_EQUAL(scdP.getContactThreshold(), real_t(0));
   SolidCylinderContactDetection scdPSw;
   scdPSw.cylinderContactSkin = real_t(0.5);
   WALBERLA_CHECK(scdPSw(1, 0, pillarShape, cylShape, ac2));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scdPSw.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
}

// ---------------------------------------------------------------------------------------------
// boxContactSkin (ModelParameters): every floating box is inflated by m (each half-extent + m) for
// contact detection, so the gap to another inflated box shrinks by 2m, to a plane or a solid pillar by m; normals are
// unchanged. Checked on the SAT path (useAnalyticalContactDetection true, as in the hopper prms) and on GJK/EPA.
// ---------------------------------------------------------------------------------------------
void testBoxSkin()
{
   const mesa_pd::data::Box boxShape(Vec3(real_t(2), real_t(2), real_t(2)));

   // Box x Box: axis-aligned, centres 2.3 apart along x -> face gap 0.3.
   for (bool analytic : { true, false })
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto b0 = ps->create();
      b0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
      b0->setRotation(mesa_pd::Rot3());
      auto b1 = ps->create();
      b1->setPosition(Vec3(real_t(2.3), real_t(0), real_t(0)));
      b1->setRotation(mesa_pd::Rot3());

      SolidCylinderContactDetection scd0;
      scd0.useAnalyticalContactDetection = analytic;
      WALBERLA_CHECK(!scd0(0, 1, boxShape, boxShape, ac));

      // skin 0.25: inflated gap 0.3 - 0.5 = -0.2; point in the middle of the inflated overlap, x = (1.25 + 1.05) / 2
      SolidCylinderContactDetection scd;
      scd.useAnalyticalContactDetection = analytic;
      scd.boxContactSkin                = real_t(0.25);
      WALBERLA_CHECK(scd(0, 1, boxShape, boxShape, ac));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(scd.getContactNormal()[0]), real_t(1), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactPoint()[0], real_t(1.15), real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL(scd.boxContactSkin, real_t(0.25)); // restored after the recursive call

      // skin 0.1: inflated gap 0.1 > 0 -> separated
      SolidCylinderContactDetection scdSmall;
      scdSmall.useAnalyticalContactDetection = analytic;
      scdSmall.boxContactSkin                = real_t(0.1);
      WALBERLA_CHECK(!scdSmall(0, 1, boxShape, boxShape, ac));

      // rotated pair: the skin gives exactly the result of the inflated shapes without skin
      b1->setPosition(Vec3(real_t(2.2), real_t(0.4), real_t(-0.3)));
      b1->setRotation(mesa_pd::Rot3(Vec3(real_t(0.3), real_t(-0.5), real_t(0.8)).getNormalized(), real_t(0.6)));
      const mesa_pd::data::Box inflated(Vec3(real_t(2.5), real_t(2.5), real_t(2.5)));
      SolidCylinderContactDetection scdRot;
      scdRot.useAnalyticalContactDetection = analytic;
      scdRot.boxContactSkin                = real_t(0.25);
      SolidCylinderContactDetection scdInfl;
      scdInfl.useAnalyticalContactDetection = analytic;
      const bool hitRot  = scdRot(0, 1, boxShape, boxShape, ac);
      const bool hitInfl = scdInfl(0, 1, inflated, inflated, ac);
      WALBERLA_CHECK_EQUAL(hitRot, hitInfl);
      WALBERLA_CHECK(hitRot);
      WALBERLA_CHECK_FLOAT_EQUAL(scdRot.getPenetrationDepth(), scdInfl.getPenetrationDepth());
      WALBERLA_CHECK_FLOAT_EQUAL(scdRot.getContactNormal(), scdInfl.getContactNormal());
   }

   // Box x wall: box (half-extent 1) with its bottom face 0.3 above the plane z = 0; skin 0.5 -> depth -0.2, normal +z.
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto bx = ps->create();
      bx->setPosition(Vec3(real_t(0), real_t(0), real_t(1.3)));
      bx->setRotation(mesa_pd::Rot3());
      auto wall = ps->create();
      wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
      const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));
      SolidCylinderContactDetection scd0;
      WALBERLA_CHECK(!scd0(0, 1, boxShape, wallShape, ac));
      SolidCylinderContactDetection scd;
      scd.boxContactSkin = real_t(0.5);
      WALBERLA_CHECK(scd(0, 1, boxShape, wallShape, ac));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)), real_t(1e-6));
      SolidCylinderContactDetection scdSw;
      scdSw.boxContactSkin = real_t(0.5);
      WALBERLA_CHECK(scdSw(1, 0, wallShape, boxShape, ac));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scdSw.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
   }

   // Box x solid pillar: box (half-extent 1) centred at x = 2.3 next to a pillar of radius 1 along z -> gap 0.3.
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto bx = ps->create();
      bx->setPosition(Vec3(real_t(2.3), real_t(0), real_t(0)));
      bx->setRotation(mesa_pd::Rot3());
      auto pillar = ps->create();
      pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
      const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));
      SolidCylinderContactDetection scd0;
      WALBERLA_CHECK(!scd0(0, 1, boxShape, pillarShape, ac));
      SolidCylinderContactDetection scd;
      scd.boxContactSkin = real_t(0.5);
      WALBERLA_CHECK(scd(0, 1, boxShape, pillarShape, ac));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs(scd.getContactNormal()[0]), real_t(1), real_t(1e-6));
      SolidCylinderContactDetection scdSw;
      scdSw.boxContactSkin = real_t(0.5);
      WALBERLA_CHECK(scdSw(1, 0, pillarShape, boxShape, ac));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scdSw.getPenetrationDepth(), real_t(-0.2), real_t(1e-6));
   }
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes;
   testSphereSphereContact();
   testSphereSphereNoContact();
   testSphereBoxContact();
   testSphereBoxNoContact();
   testBoxBoxContact();
   testBoxBoxNoContact();
   testCylinderBoxContact();
   testCylinderBoxNoContact();
   testCylinderBoxSwappedOrder();
   testCylinderCylinderContact();
   testCylinderCylinderNoContact();
   testSphereCylinderSideContact();
   testSphereCylinderRotationAware();
   testSphereInsideCylinderEnvelopBranch();
   testSphereCylinderNoContact();
   testSphereCylindricalBoundaryContact();
   testSphereCylindricalBoundaryNoContact();
   testSphereHalfSpaceContact();
   testSphereHalfSpaceNoContact();
   testBoxHalfSpaceContactAxisAligned();
   testBoxHalfSpaceContactRotated();
   testBoxHalfSpaceNoContact();
   testBoxCylindricalBoundaryContact();
   testBoxCylindricalBoundaryAxialInvariance();
   testBoxCylindricalBoundaryNoContact();
   testCylinderPillarContact();
   testCylinderPillarNoContact();
   testCylinderPillarBruteForceCrossCheck();
   testCylinderHalfSpaceContact();
   testCylinderHalfSpaceNoContact();
   testCylinderSkinCylinderCylinder();
   testCylinderSkinHalfSpaceAndPillar();
   testBoxSkin();
   testSphereCylindricalBoundaryContainmentContact();
   testSphereCylindricalBoundaryContainmentNoContact();
   testSphereCylindricalBoundaryContainmentAxialBounds();
   testSphereCylindricalBoundaryContainmentHalfClipExcludesBelow();
   testSphereCylindricalBoundaryContainmentHalfClipStraddling();
   testBoxCylindricalBoundaryContainmentContact();
   testBoxCylindricalBoundaryContainmentNoContact();
   testBoxCylindricalBoundaryContainmentRadialPreGate();
   testCylinderCylindricalBoundaryContainmentContact();
   testCylinderVerticalContainer();
   testCylinderCylindricalBoundaryContainmentNoContact();
   testBoxCylindricalBoundaryContainmentDualCornerContact();
   testBoxCylindricalBoundaryContainmentSingleCornerNoSecondary();
   testCylinderCylindricalBoundaryContainmentDualRimContact();

   return EXIT_SUCCESS;
}
