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
//! Unit tests for walberla::antidunes::SolidCylinderContactDetection. Covers the full matrix every
//! floater shape (Sphere, Box, Cylinder) can be asked to resolve against: every floater x floater
//! pair, every floater x HalfSpace ("wall/floor") pair, and every floater x CylindricalBoundary
//! ("pillar") pair -- 12 physically-meaningful pairs in total (HalfSpace/CylindricalBoundary never
//! contact each other or themselves, since both are static). Each pair gets a positive contact case
//! plus a no-contact negative control; pairs with their own swapped-argument overload (the 5 custom
//! ones below, plus the GJK/EPA and pillar ones) additionally assert order-invariance.
//!
//! Breakdown by resolution path:
//! - 6 CUSTOM overloads (hand-traced expected values): Sphere x Cylinder "floating log", Sphere x
//!   CylindricalBoundary, Box x HalfSpace, Box x CylindricalBoundary, Cylinder x CylindricalBoundary,
//!   Cylinder x HalfSpace. The last one was ADDED after testCylinderHalfSpaceContact's derivation
//!   surfaced that this pair, when left to fall through to GeneralContactDetection's generic
//!   HalfSpace-vs-anything template, returned contactNormal = -wallNormal (into the wall) -- the
//!   opposite sign from every other HalfSpace pair (Sphere and Box both give +wallNormal, matching
//!   detectSphereHalfSpaceCollision's proven-correct-in-production convention). Fixed in
//!   ParticleKernels.h by giving Cylinder x HalfSpace its own overload using the identical
//!   support()-based formula as Box x HalfSpace.
//! - 3 pairs INHERITED as exact analytic overloads (Sphere x Sphere, Sphere x Box, Sphere x
//!   HalfSpace), exercised through THIS class specifically to confirm
//!   `using GeneralContactDetection::operator();` doesn't shadow them. Sphere x HalfSpace was
//!   missing from an earlier version of this coverage list -- an actual gap, not a stylistic
//!   omission, since AnalyticContactDetection provides it the same way as the other two.
//! - 3 pairs with NO custom/analytic override, falling through to the inherited GJK/EPA path: Box x
//!   Box, Cylinder x Box, Cylinder x Cylinder -- simple axis-aligned geometry, unambiguous from
//!   separating-axis theorem, sign convention confirmed by an instrumented run.
//!
//! testCylinderPillarBruteForceCrossCheck additionally verifies the fixed-point closest-point search
//! independently, via dense surface sampling that never calls the algorithm under test.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition itself -- it only compiles in the real Antidunes.cpp because
// BuildTimeloop.h happens to pull it in first via the free-surface dynamics headers. Include it
// explicitly here too (test-side only, no production file touched) so this translation unit doesn't
// depend on that same include-order accident.
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
// Analytic overloads inherited from AnalyticContactDetection (Sphere x Sphere, Sphere x Box) --
// NOT custom to SolidCylinderContactDetection, but exercised here through OUR class specifically
// to confirm the `using GeneralContactDetection::operator();` declaration correctly exposes them
// rather than accidentally shadowing them with our own overloads.
// ---------------------------------------------------------------------------------------------

void testSphereSphereContact()
{
   // Reuses the EXACT geometry already validated in
   // tests/mesa_pd/collision_detection/AnalyticContactDetection.cpp's checkSphereSphereCollision()
   // -- known-correct expected values taken directly from that independently-authored test, not
   // re-derived here. Additionally cross-checked against pybullet (Bullet Physics, a completely
   // separate physics engine) for this same geometry: matched to ~1e-6 on penetration depth,
   // normal, AND contact point -- for two spheres the contact point IS a unique geometric ground
   // truth (the point on the line between centres), unlike the general-polytope pairs below.
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
   // Hand-traced against analytic::detectSphereBoxCollision (AnalyticCollisionFunctions.h): sphere
   // outside the box, clamped projection point p=(0.5,0,0) (box half-extent), n=d-q=(0.5,0,0),
   // dist=0.5-0.6=-0.1.
   // Cross-checked against pybullet for this same geometry: depth and normal matched exactly, but
   // the contact POINT did not (pybullet returned the midpoint between the two surfaces, 0.45,
   // 0,0; waLBerla's detectSphereBoxCollision explicitly returns the BOX surface point, 0.5,0,0 --
   // see its own contactPoint=pos2+q). That's a genuine, documented convention difference between
   // the two engines, not a bug in either -- so the contactPoint assertion below is checked against
   // waLBerla's own documented convention, not against pybullet.
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
// Pairs with NO custom override in this class -- fall through to the inherited GJK/EPA path
// (GeneralContactDetection::collideGJKEPA, margin=1e-4). This is a genuine gap that existed in the
// codebase before these tests: neither tests/mesa_pd/collision_detection/EPA.cpp nor any other
// existing test exercises Box x Box, Cylinder x Box, or Cylinder x Cylinder (a comment on
// SolidCylinderContactDetection claiming otherwise was checked and found to be incorrect).
//
// INDEPENDENT VERIFICATION: penetrationDepth and contactNormal below were cross-checked against
// pybullet (Bullet Physics) -- a completely separate, independently-developed physics engine, not
// a re-derivation of waLBerla's own formulas -- using the SAME geometry. Depth and normal matched
// pybullet's output exactly (to ~1e-6) for every pair, on top of independently matching the
// unambiguous separating-axis-theorem prediction for these simple, single-axis-overlap geometries.
//
// contactPoint is NOT asserted beyond its well-defined axis component, and this is deliberate, not
// an oversight: for two general convex polytopes (unlike two spheres, where the contact point is
// trivially the point on the line between centres), the exact contact-point coordinate is an
// ALGORITHM-SPECIFIC CHOICE, not a unique geometric ground truth -- confirmed directly by testing
// this. A first attempt at a fully symmetric box-box case (identical cubes, offset along one axis
// only) has a whole SHARED FACE of equally-valid contact points; waLBerla's EPA and pybullet's GJK
// picked different points on that same face (both "correct"). Breaking the symmetry (different box
// sizes, offset in y and z too) gives a case with a unique deepest-point pair, but pybullet's own
// cylinder shapes are internally approximated as a discretized convex hull (confirmed separately:
// a probe placed at the exact analytic surface distance read ~0.15 instead of 0.2, a systematic
// discretization artifact) -- so for the cylinder-involving pairs even pybullet's own
// off-axis contact-point components carry some approximation error, not a value to assert against.
// ---------------------------------------------------------------------------------------------

void testBoxBoxContact()
{
   // Two unit cubes (half-extent 0.5 each), centres 0.8 apart along x, y/z fully aligned -> x is
   // the unambiguous minimum-penetration axis: overlap = 0.5+0.5-0.8 = 0.2.
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
   // Only the x-component is asserted: with identical, perfectly y/z-aligned boxes, the entire
   // overlapping face (any y,z within [-0.5,0.5]^2) is an equally valid contact point -- confirmed
   // by cross-checking against pybullet, which picks a DIFFERENT (also valid) point on that same
   // face. x=0.4 (the overlap region's midpoint along the one axis that actually matters) is the
   // one component with a unique, unambiguous answer.
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
   // Cylinder radius 0.5 (axis = world z via identity rotation) vs. a unit cube (half-extent 0.5),
   // centres 0.9 apart along x -> overlap = 0.5+0.5-0.9 = 0.1 along the unambiguous x axis (both
   // shapes' z-extents fully overlap: cylinder z in [-1,1], box z in [-0.5,0.5]).
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
   // Only x asserted -- see the section header: the exact contact point isn't a unique ground
   // truth for general convex shapes, and pybullet's own cylinder approximation adds further noise
   // to any off-axis component specifically for cylinder-involving pairs.
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
   // Unlike the custom overloads (whose "reverse" version hands off to the identical formula with
   // idx1/idx2 swapped by construction), this pair falls through to GeneralContactDetection's
   // GJK/EPA path, which internally re-orders by particle UID before computing anything (see
   // GeneralContactDetection::operator(), the "ensure idx2's UID > idx1's UID" check) -- so calling
   // with either argument order must converge on the exact SAME underlying computation, not just an
   // equivalent one. That's what this test actually verifies, for the one GJK/EPA-fallback pair
   // that didn't have a swap check yet (Box x Box and Cylinder x Cylinder use identical shapes on
   // both sides, so a swap test there wouldn't exercise anything a mixed-shape pair doesn't already
   // cover more thoroughly).
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
   // Two parallel cylinders (radius 0.5, axis = world z via identity rotation), centres 0.8 apart
   // along x -> overlap = 0.5+0.5-0.8 = 0.2 along the unambiguous x axis (z-extents fully overlap).
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
   // Only x asserted -- see the section header comment.
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
// Sphere (idx0) x Cylinder (idx1) -- "floating log". This pair's geometry (r=1 sphere at 1.5,0,0
// vs r=1 h=4 cylinder at origin) was additionally cross-checked against pybullet: depth, normal,
// AND contact point all matched to ~1e-4 -- waLBerla's own documented convention here ("mid-overlap
// point... half the overlap past the sphere surface") happens to be the same midpoint convention
// pybullet uses, so unlike Sphere x Box above, this one's contact point IS fully independently
// confirmed, not just depth/normal.
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
   // Rotate the cylinder 90 deg about its own (z) axis: this leaves the axis direction unchanged
   // but must correctly rotate the RADIAL contact normal/point (R applied to nBody), exercising
   // the body-frame transform pipeline rather than being a no-op.
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
   sphere->setPosition(Vec3(real_t(1.5), real_t(0), real_t(7))); // z=7 arbitrary -- see axial-invariance note below
   auto pillar = ps->create();
   pillar->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Sphere              sphereShape(real_t(0.6));
   const mesa_pd::data::CylindricalBoundary pillarShape(real_t(1), Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, pillarShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.1));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(1), real_t(0), real_t(0)));
   // contactPoint's z matches the sphere's own z=7 -- confirms the pillar's radial-only projection
   // (delta - (delta.axis)*axis) correctly discards the axial component instead of biasing toward
   // the pillar's own z=0, i.e. the infinite-pillar axial invariance is exercised here directly.
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
// Sphere (idx0, "sphere floater") x HalfSpace (idx1, "wall/floor") -- inherited unchanged from
// AnalyticContactDetection (detectSphereHalfSpaceCollision), same rationale as testSphereSphereContact
// / testSphereBoxContact: exercised through THIS class to confirm it isn't shadowed by the custom
// overloads above, not because the math itself is new. This pair was previously untested in this
// file entirely -- see the file header note.
// ---------------------------------------------------------------------------------------------

void testSphereHalfSpaceContact()
{
   // Same wall geometry as testBoxHalfSpaceContactAxisAligned (position, normal) for direct
   // cross-shape comparison. Expected values hand-traced from detectSphereHalfSpaceCollision
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
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(0.5), real_t(0.5), real_t(0)));

   // Swapped-argument overload must give the identical result.
   SolidCylinderContactDetection scd2;
   WALBERLA_CHECK(scd2(1, 0, wallShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), scd2.getPenetrationDepth());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), scd2.getContactNormal());
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), scd2.getContactPoint());
}

void testBoxHalfSpaceContactRotated()
{
   // Same z-offset as the axis-aligned case (so the penetration depth matches), but an ASYMMETRIC
   // box (edge lengths differ in x/y) rotated 90 deg about z -- a cube would look identical after
   // this rotation, so this specifically exercises that R is applied to find the correct deepest
   // corner in the rotated body frame, not just passed through as a no-op.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(0.3)));
   box->setRotation(mesa_pd::Rot3(Vec3(real_t(0), real_t(0), real_t(1)), math::pi * real_t(0.5)));
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box       boxShape(Vec3(real_t(0.6), real_t(1.0), real_t(1.0))); // half-extents (0.3,0.5,0.5)
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, wallShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getPenetrationDepth(), real_t(-0.2)); // same as axis-aligned: only z matters here
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(1)));
   // Different from the axis-aligned case's (0.5,0.5,0) -- confirms the rotation genuinely changes
   // which corner is deepest, not just passing the unrotated corner through.
   WALBERLA_CHECK_FLOAT_EQUAL(scd.getContactPoint(), Vec3(real_t(-0.5), real_t(0.3), real_t(0)));
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
   // Same geometry as testBoxCylindricalBoundaryContact, box shifted along the pillar's (infinite)
   // axis by z=5 -- penetration depth/normal must be unchanged, and the contact point's z must
   // follow the box (not stay at the pillar's own z=0), confirming the fixed-point search's
   // closest-point-to-the-axis-LINE logic is genuinely axis-position-independent.
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
   // Independent check of the fixed-point closest-point search: densely sample the log's curved
   // surface directly (no fixed-point iteration, no clamping logic -- just perpendicular distance
   // to the pillar axis line) and confirm the true minimum distance matches what the analytic
   // implementation computed for the SAME geometry as testCylinderPillarContact.
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

   // Dense brute-force sampling of the log's curved surface (identity rotation, so body == world
   // frame here): theta in [0,2*pi), z in [-halfLen,halfLen].
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

   // Sampling-resolution-limited tolerance: with 720x400 samples over a radius-1 circle, the
   // angular gap between samples is ~0.5 deg, giving a worst-case sagitta well under 1e-4.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(minDist, analyticDistToAxis, real_t(1e-3));
}

// ---------------------------------------------------------------------------------------------
// Cylinder (idx0, "floating log") x HalfSpace (idx1, "wall/floor") -- has its own overload in
// SolidCylinderContactDetection (ParticleKernels.h), identical support()-based formula to Box x
// HalfSpace above. Geometry below is chosen so the cylinder's flat end cap faces the wall (axis
// parallel to the wall normal): Cylinder::support() explicitly returns the CAP CENTRE (not an
// ambiguous rim point) whenever the search direction is purely axial, so this is a clean,
// unambiguous case to hand-trace, not a degenerate one.
//
// This pair used to fall through to GeneralContactDetection's generic HalfSpace-vs-anything
// template, which returns contactNormal = -wallNormal (into the wall) -- the opposite sign from
// every other HalfSpace pair (Sphere and Box both give +wallNormal, matching
// detectSphereHalfSpaceCollision's proven-correct-in-production convention -- see
// testSphereHalfSpaceContact / testBoxHalfSpaceContactAxisAligned). Fixed by adding a dedicated
// Cylinder x HalfSpace overload using the same formula as Box's; the expected values below reflect
// the fixed, +wallNormal convention.
// ---------------------------------------------------------------------------------------------

void testCylinderHalfSpaceContact()
{
   // Cylinder axis = world z (identity rotation), radius 1, height 4 (halfLen 2), centred at
   // (0,0,1) -> bottom cap centre at world (0,0,-1), 1 LU below the wall plane at z=0.
   // Hand-traced from the fixed overload (same formula as Box x HalfSpace):
   // deepest = c + R*cyl.support(-(R^T*n)) = (0,0,1) + cyl.support((0,0,-1));
   // Cylinder::support((0,0,-1)) in body frame = (0,0,-2) (pure axial -> cap CENTRE, not a rim
   // point) -> deepest = (0,0,-1); penetrationDepth = n.deepest - d = -1-0 = -1;
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

   // Swapped-argument overload must give the identical result (the generic HalfSpace overloads
   // always internally reorder so the half-space is geo1, regardless of call order).
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

   return EXIT_SUCCESS;
}
