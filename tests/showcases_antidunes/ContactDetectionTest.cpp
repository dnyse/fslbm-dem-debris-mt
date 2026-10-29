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
   // Flush face-on (query direction exactly axis-aligned, both lateral components exactly zero): the
   // TRUE contact region is the whole bottom face, not one arbitrary corner -- (0,0,0), the face
   // centre, is the physically-symmetric representative point (Box::support()'s own fix, see its
   // header comment). Was (0.5,0.5,0) before that fix (an arbitrary corner among 4 equally valid
   // ones) -- confirmed via mesa_pd::data::Cylinder/Box::support()'s degenerate-direction guard.
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
   // NOTE (updated alongside Box::support()'s degenerate-direction fix): a rotation purely about Z --
   // the SAME axis as the wall normal -- never changes the wall-normal direction expressed in the
   // box's own body frame (Z stays Z under any Z-rotation), so it stays exactly as degenerate/flush
   // as the unrotated axis-aligned case above, regardless of angle. That no longer exercises rotated
   // corner-finding at all post-fix (both now correctly return the SAME face-centre point) -- so this
   // test instead rotates about X (an axis that genuinely changes the wall-facing direction in body
   // frame), which puts the box corner-first against the wall, a well-defined, non-degenerate single
   // contact point, and DOES exercise that R is applied to find the correct deepest corner in the
   // rotated body frame rather than passed through as a no-op.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   const real_t tilt = math::pi / real_t(6); // 30 deg about X -- clearly non-degenerate, not a special angle

   auto box = ps->create();
   box->setPosition(Vec3(real_t(0), real_t(0), real_t(0.3)));
   box->setRotation(mesa_pd::Rot3(Vec3(real_t(1), real_t(0), real_t(0)), tilt));
   auto wall = ps->create();
   wall->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));

   const mesa_pd::data::Box       boxShape(Vec3(real_t(0.6), real_t(1.0), real_t(1.0))); // half-extents (0.3,0.5,0.5)
   const mesa_pd::data::HalfSpace wallShape(Vec3(real_t(0), real_t(0), real_t(1)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, boxShape, wallShape, ac));

   // Body-frame direction toward the wall: R^T*(0,0,-1) = (0, -sin(tilt), -cos(tilt)) for a rotation
   // about X -- X stays exactly 0 (rotating about X can never introduce an X-component), so that axis
   // is still degenerate here (the deepest CORNER is really a deepest EDGE spanning the box's full
   // X-extent) and the fixed support() correctly keeps the X-offset at 0 (the edge's own midpoint) --
   // but Y and Z are both clearly, unambiguously nonzero, so the deepest point is otherwise a genuine,
   // well-defined corner in the Y-Z cross-section: body-frame support point (0, -0.5, -0.5).
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

// ---------------------------------------------------------------------------------------------
// CylindricalBoundary CONTAINMENT mode (cb.setContainment(true)) -- the mirror image of the
// pillar (outside-solid) overloads above: a body poking OUTSIDE the tube gets pushed back toward
// the axis. Exercises the new farthest-point search (Sphere/Box/Cylinder), the optional axial
// bounds, and the optional half-clip plane (e.g. an arch tunnel is only solid above its own
// springing line) added to SolidCylinderContactDetection alongside the pillar overloads it never
// touches (isContainment() defaults to false, so every test above this section is unaffected).
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

   // Radially this sphere clearly pokes outside (farthest reach 4.7+0.6=5.3 > radius 5), but it
   // sits at world X=10, well outside the tube's bounded [xMin=-3,xMax=3] (offsets from the tube's
   // own position at X=0) axial span -- the axial gate must reject it before the radial search
   // even runs, exactly like ComputeTunnelContacts' own X-extent gate (ParticleKernels.h).
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
   // Sphere entirely on the non-solid side of the half-clip plane (e.g. below an arch's own
   // springing line, Z < 0 here): even though it pokes radially outside the FULL tube's radius,
   // there is no real wall there, so this must report no contact -- mirrors ComputeTunnelContacts'
   // own springing-line exclusion (ParticleKernels.h, "below the springing line -- no real wall
   // here").
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
   // Sphere straddling the clip plane (Z=0): its unconstrained farthest-from-axis point would dip
   // below Z=0 (non-solid), so the resulting contact point must be pulled back onto the clip arc's
   // boundary -- verified via invariants (contact fires, the reported point stays on the solid
   // side, and it is a genuine radial overlap) rather than a hand-derived closed form, since the
   // exact arc-endpoint position depends on the quartic solver already validated elsewhere.
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

   // Farthest corner is unambiguously (Y=4.5+1, Z=0.3+1) -- both offsets add, no tie to break.
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
   // job 4141165 postmortem: two arch tunnels can share the same axial (X) span (this bridge's own
   // 3 arches all span X=[61,67]), and a box merely within that shared span but nowhere near THIS
   // tube's own axis in Y-Z used to still get a (huge, spurious) farthest-corner overlap computed
   // against it -- confirmed in production as overlap=15.75 (251% of interaction radius) for a box
   // ~26 LU away from the arch it was wrongly "contained" by. This tube's axial span is world-X
   // [-3,3] (offsets from the tube's own position at X=0); the box sits well within that span but
   // 30 LU away in Y -- the radial pre-gate must reject it before the corner search ever runs.
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
   // Tilted (not axis-aligned) on purpose: a floating cylinder oriented EXACTLY along a world axis
   // that is itself exactly perpendicular to the tube's own axis hits a genuine pre-existing
   // degenerate-quartic edge case in farthestThetaOnCircleToLine (PaperAnalyticalContacts.h) --
   // found while writing this test (c4/c2/c0 all vanish simultaneously for that special alignment,
   // and the Durand-Kerner solve doesn't reliably recover the true root then). That function is
   // shared, already exercised elsewhere, and out of scope for the containment feature itself, so
   // this test instead uses a generic tilt (representative of a real tumbling log) and validates
   // via invariants rather than a hand-derived closed form.
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
   // The reported point must genuinely lie on the floater's own rim: exactly cylRadius from the
   // floater's OWN axis line, and exactly at one end (+-hh) of it (never a fractional/interior
   // axial offset -- see the "always at an endpoint" derivation in the calling code's comment).
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

// [Claude, 2026-10-03] VERTICAL round container (DropTest ContainerCylinders, Utility.h createContainerCylinder): a
// containment tube with its axis along z, radius 25.3 (Kodam et al. 2010b: 50.6 mm container), Kodam's cylinders
// (radius 4, length 5.3). Three orientation classes that occur in a packing: a disc lying flat (axis PARALLEL to the
// tube axis), a disc standing upright (axis exactly PERPENDICULAR to the tube axis -- the degenerate-quartic case noted
// above), and a generic tilt. Expected overlaps are derived by hand (flat, upright) or checked against a dense sampling
// of the rim (all cases).
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
   // upright disc, axis exactly along x, centre on the x axis: farthest point is the rim face at x0 + hh, the whole
   // flat face reaches radial distance sqrt((x0 + hh)^2 + rc^2) at its top/bottom... in plan view the face is a
   // segment from (x0 + hh, -rc) to (x0 + hh, rc): farthest radial sqrt((x0 + hh)^2 + rc^2)
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
// Secondary contact points (arch/tunnel containment, Box/Cylinder floaters only -- Sphere is
// inherently single-point). The farthest-point search's top-2 tracking exposes a SECOND
// independently-breaching point via hasSecondaryContact()/getSecondaryContact*() instead of
// silently discarding it, so an object touching the curve at two points genuinely gets a force
// response at both, not just whichever one happens to be worse.
// ---------------------------------------------------------------------------------------------

void testBoxCylindricalBoundaryContainmentDualCornerContact()
{
   // Same geometry as testBoxCylindricalBoundaryContainmentContact -- but an UNROTATED box's farthest
   // (Y,Z) corner pair is always tied along local X (two corners differing only in local X share the
   // same world Y,Z, hence the identical radial distance from the tube axis). That means BOTH ends of
   // this box's top edge genuinely poke past the tube radius at once, not just one -- exactly the
   // "does an object touching the arch at two points get modeled at both" case. Loop order (sx outer,
   // so sx=-1 is visited before sx=+1) makes the X=-1 corner win the tie and become primary; X=+1 is
   // the tied runner-up.
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
   // Rotated (not axis-aligned) on purpose: an unrotated box always ties its farthest corner pair
   // along local X (see the dual-corner test above), so honestly testing the "only one point, no
   // secondary" case requires breaking that degeneracy. Same tilt convention as
   // testCylinderCylindricalBoundaryContainmentContact and for the same reason: a generic,
   // non-axis-aligned rotation has no reason to reintroduce an exact tie between two different
   // corners.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   // Numerically searched (not hand-derived) for a comfortable margin on both sides: worst corner
   // radial dist ~5.208 (overlap ~0.208), runner-up ~4.814 (safely ~0.186 short of the radius) --
   // confirmed via a standalone corner-enumeration script against this exact rotation before writing
   // this test, not assumed.
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
   // A log lying across the tube (tilted, same convention as the primary containment test above) with
   // both end-cap rims independently poking past the tube radius at once -- the "touches the arch at
   // two points" case for a Cylinder floater specifically. Validated via invariants (own-rim
   // membership, genuine overlap on both ends) rather than a hand-derived closed form, matching
   // testCylinderCylindricalBoundaryContainmentContact's own reasoning: the exact farthest-angle
   // solve is shared, already exercised elsewhere, and out of scope for this feature.
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
   // The two rims are genuinely different points (opposite ends of the log), not a duplicate report.
   WALBERLA_CHECK((scd.getContactPoint() - scd.getSecondaryContactPoint()).length() > real_t(1));
}

// ---------------------------------------------------------------------------------------------
// Arch deck-fill / haunch solid geometry (SimulationParams.h: archDeckFillBoxGeometry,
// archHaunchBoxGeometry). These are pure center/edgeLength math with no mesa_pd particle involved
// -- SetupStages.h::createArchSolidElements is the only caller that turns the result into a real
// FixedBox via createBox(). Covers: (1) the real bridge's own three arches produce the exact
// deck-fill box the prm files' `deckTopZ 76;` was chosen to match, (2) the pre-existing "no opt-in"
// default stays a no-op (backward compatibility), (3) this bridge's own haunch boxes are correctly
// degenerate (its BodyBC abutments have zero width beyond the circular opening), (4) a synthetic
// wider abutment produces real, correctly-placed haunch boxes, and (5) the resulting box actually
// collides with a floater above the crown, and does NOT collide with one still inside the open bore.
// ---------------------------------------------------------------------------------------------

void testArchDeckFillBoxGeometryRealBridge()
{
   // Left/central/right arches, exactly as defined in prms/test_bridge_v1_*_gas.prm's ArchTunnels
   // block, all with deckTopZ 76 (matching the BodyBC arch boxes' own max Z -- see the artifact's
   // own derivation). yMin/yMax left at their defaults (center.y -+ radius) for all three.
   const ArchTunnelParams left(Vector3< real_t >(0, 5, 56), real_t(5), real_t(61), real_t(67),
                                std::numeric_limits< real_t >::quiet_NaN(), std::numeric_limits< real_t >::quiet_NaN(), real_t(76));
   const ArchTunnelParams central(Vector3< real_t >(0, 30, 56), real_t(10), real_t(61), real_t(67),
                                   std::numeric_limits< real_t >::quiet_NaN(), std::numeric_limits< real_t >::quiet_NaN(), real_t(76));
   const ArchTunnelParams right(Vector3< real_t >(0, 55, 56), real_t(5), real_t(61), real_t(67),
                                 std::numeric_limits< real_t >::quiet_NaN(), std::numeric_limits< real_t >::quiet_NaN(), real_t(76));

   const auto leftBox = archDeckFillBoxGeometry(left);
   WALBERLA_CHECK(leftBox.has_value());
   // crownZ = 56+5 = 61; center Z = (61+76)/2 = 68.5; edgeLength Z = 76-61 = 15.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(leftBox->center, Vec3(real_t(64), real_t(5), real_t(68.5)), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(leftBox->edgeLength, Vec3(real_t(6), real_t(10), real_t(15)), real_t(1e-9));

   const auto centralBox = archDeckFillBoxGeometry(central);
   WALBERLA_CHECK(centralBox.has_value());
   // crownZ = 56+10 = 66; center Z = (66+76)/2 = 71; edgeLength Z = 76-66 = 10.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(centralBox->center, Vec3(real_t(64), real_t(30), real_t(71)), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(centralBox->edgeLength, Vec3(real_t(6), real_t(20), real_t(10)), real_t(1e-9));

   const auto rightBox = archDeckFillBoxGeometry(right);
   WALBERLA_CHECK(rightBox.has_value());
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rightBox->center, Vec3(real_t(64), real_t(55), real_t(68.5)), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rightBox->edgeLength, Vec3(real_t(6), real_t(10), real_t(15)), real_t(1e-9));
}

void testArchDeckFillBoxGeometryDefaultIsNoOp()
{
   // The 4-arg constructor (no deckTopZ given) is exactly what every ArchTunnels block that hasn't
   // opted into deck-fill still uses -- deckTopZ defaults to the crown itself, so the box would have
   // zero height. Confirms every pre-existing scenario keeps its exact prior behaviour (no box at all).
   const ArchTunnelParams tunnel(Vector3< real_t >(0, 5, 56), real_t(5), real_t(61), real_t(67));
   WALBERLA_CHECK(!archDeckFillBoxGeometry(tunnel).has_value());
}

void testArchBlockParamsDimsConversionMatchesRealBridge()
{
   // createArchBlock's dimension-based builder (SetupStages.h) takes base/xDim/yDim/zDim/radius
   // instead of xMin/xMax/yMin/yMax/deckTopZ -- archBlockParams (SimulationParams.h) is the pure
   // conversion between the two, checked here against the left arch's own real numbers: springing
   // line at Z=56 (base), deck at Z=76 -> zDim=20 (NOT 15 -- zDim is measured from the springing
   // line, not the crown, unlike deckTopZ which IS measured from world zero the same way center.z
   // is), full opening width yDim=10 (2*radius, matching this bridge's own zero-haunch geometry),
   // block thickness xDim=6 (61 to 67).
   const auto tunnel = archBlockParams(Vector3< real_t >(64, 5, 56), real_t(6), real_t(10), real_t(20), real_t(5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(tunnel.center, Vector3< real_t >(64, 5, 56), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL(tunnel.radius, real_t(5));
   WALBERLA_CHECK_FLOAT_EQUAL(tunnel.xMin, real_t(61));
   WALBERLA_CHECK_FLOAT_EQUAL(tunnel.xMax, real_t(67));
   WALBERLA_CHECK_FLOAT_EQUAL(tunnel.yMin, real_t(0));
   WALBERLA_CHECK_FLOAT_EQUAL(tunnel.yMax, real_t(10));
   WALBERLA_CHECK_FLOAT_EQUAL(tunnel.deckTopZ, real_t(76));

   // Feeding this straight into the already-tested deck-fill geometry function must reproduce
   // exactly the same box as testArchDeckFillBoxGeometryRealBridge's own left-arch case -- confirms
   // the dims-based builder and the prm-based path agree, not just that each half works alone.
   const auto deckFill = archDeckFillBoxGeometry(tunnel);
   WALBERLA_CHECK(deckFill.has_value());
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(deckFill->center, Vec3(real_t(64), real_t(5), real_t(68.5)), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(deckFill->edgeLength, Vec3(real_t(6), real_t(10), real_t(15)), real_t(1e-9));
}

void testArchHaunchBoxGeometryRealBridgeIsDegenerate()
{
   // This bridge's own BodyBC arch boxes have Y-span EXACTLY equal to each arch's own diameter (no
   // separate flat haunch wall beside the curve) -- confirmed directly against the BodyBC coordinates
   // in prms/test_bridge_v1_*_gas.prm. yMin/yMax at their defaults must therefore make BOTH haunch
   // sides degenerate for all three arches, for every one of this bridge's scenario prms.
   const ArchTunnelParams left(Vector3< real_t >(0, 5, 56), real_t(5), real_t(61), real_t(67));
   const ArchTunnelParams central(Vector3< real_t >(0, 30, 56), real_t(10), real_t(61), real_t(67));
   const ArchTunnelParams right(Vector3< real_t >(0, 55, 56), real_t(5), real_t(61), real_t(67));

   for (const auto& tunnel : { left, central, right })
   {
      WALBERLA_CHECK(!archHaunchBoxGeometry(tunnel, ArchHaunchSide::Left).has_value());
      WALBERLA_CHECK(!archHaunchBoxGeometry(tunnel, ArchHaunchSide::Right).has_value());
   }
}

void testArchHaunchBoxGeometryWidenedAbutment()
{
   // Synthetic scenario for an abutment wider than its own opening (not true of today's bridge, but
   // a real future case per ArchTunnelParams's own comment) -- radius 5 opening at Y=5, abutment
   // running Y=[-2,12] (3 LU wider on each side than the circular opening's own Y=[0,10] span).
   const ArchTunnelParams tunnel(Vector3< real_t >(0, 5, 56), real_t(5), real_t(61), real_t(67), real_t(-2), real_t(12));

   const auto leftHaunch = archHaunchBoxGeometry(tunnel, ArchHaunchSide::Left);
   WALBERLA_CHECK(leftHaunch.has_value());
   // width = (5-5) - (-2) = 2; center Y = -2 + 2/2 = -1; center Z = 56 + 5/2 = 58.5.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(leftHaunch->center, Vec3(real_t(64), real_t(-1), real_t(58.5)), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(leftHaunch->edgeLength, Vec3(real_t(6), real_t(2), real_t(5)), real_t(1e-9));

   const auto rightHaunch = archHaunchBoxGeometry(tunnel, ArchHaunchSide::Right);
   WALBERLA_CHECK(rightHaunch.has_value());
   // width = 12 - (5+5) = 2; center Y = 12 - 2/2 = 11.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rightHaunch->center, Vec3(real_t(64), real_t(11), real_t(58.5)), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rightHaunch->edgeLength, Vec3(real_t(6), real_t(2), real_t(5)), real_t(1e-9));
}

void testArchDeckFillBoxCollidesSphereAboveCrown()
{
   // Left arch's own deck-fill box (see testArchDeckFillBoxGeometryRealBridge): center (64,5,68.5),
   // edgeLength (6,10,15) -> half-extents (3,5,7.5) -> Z in [61,76]. A sphere centred just below the
   // crown (Z=60.5, radius 1) reaches up to Z=61.5 -- 0.5 LU into the deck-fill box. This is exactly
   // the case SetupStages.h::createArchSolidElements exists to catch: a floater kicked hard enough to
   // reach above the crown while still inside the arch's own axial/Y span.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(64), real_t(5), real_t(60.5)));
   auto box = ps->create();
   box->setPosition(Vec3(real_t(64), real_t(5), real_t(68.5)));
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Sphere sphereShape(real_t(1));
   const mesa_pd::data::Box    boxShape(Vec3(real_t(6), real_t(10), real_t(15)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(scd(0, 1, sphereShape, boxShape, ac));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getPenetrationDepth(), real_t(-0.5), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(scd.getContactNormal(), Vec3(real_t(0), real_t(0), real_t(-1)), real_t(1e-9));
}

void testArchDeckFillBoxNoContactInsideOpenBore()
{
   // Same deck-fill box as above; sphere safely below the crown (top of sphere at Z=60.5 < crownZ=61)
   // -- must NOT register a contact against the deck-fill slab while merely passing through the open
   // tunnel bore below it.
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   auto sphere = ps->create();
   sphere->setPosition(Vec3(real_t(64), real_t(5), real_t(60)));
   auto box = ps->create();
   box->setPosition(Vec3(real_t(64), real_t(5), real_t(68.5)));
   box->setRotation(mesa_pd::Rot3());

   const mesa_pd::data::Sphere sphereShape(real_t(0.5));
   const mesa_pd::data::Box    boxShape(Vec3(real_t(6), real_t(10), real_t(15)));

   SolidCylinderContactDetection scd;
   WALBERLA_CHECK(!scd(0, 1, sphereShape, boxShape, ac));
}

} // namespace antidunes
} // namespace walberla

namespace walberla {
namespace antidunes {
using Vec3 = mesa_pd::Vec3;

// ---------------------------------------------------------------------------------------------
// [Claude, 2026-09-30] cylinderContactSkin (ModelParameters): every floating cylinder is dilated by the skin m for
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
// [Claude, 2026-10-03] boxContactSkin (ModelParameters): every floating box is inflated by m (each half-extent + m) for
// contact detection, so the gap to another inflated box shrinks by 2m and to a plane or a solid pillar by m; normals are
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
   testBoxSkin(); // [Claude, 2026-10-03]
   testSphereCylindricalBoundaryContainmentContact();
   testSphereCylindricalBoundaryContainmentNoContact();
   testSphereCylindricalBoundaryContainmentAxialBounds();
   testSphereCylindricalBoundaryContainmentHalfClipExcludesBelow();
   testSphereCylindricalBoundaryContainmentHalfClipStraddling();
   testBoxCylindricalBoundaryContainmentContact();
   testBoxCylindricalBoundaryContainmentNoContact();
   testBoxCylindricalBoundaryContainmentRadialPreGate();
   testCylinderCylindricalBoundaryContainmentContact();
   testCylinderVerticalContainer(); // [Claude, 2026-10-03]
   testCylinderCylindricalBoundaryContainmentNoContact();
   testBoxCylindricalBoundaryContainmentDualCornerContact();
   testBoxCylindricalBoundaryContainmentSingleCornerNoSecondary();
   testCylinderCylindricalBoundaryContainmentDualRimContact();
   testArchDeckFillBoxGeometryRealBridge();
   testArchDeckFillBoxGeometryDefaultIsNoOp();
   testArchBlockParamsDimsConversionMatchesRealBridge();
   testArchHaunchBoxGeometryRealBridgeIsDegenerate();
   testArchHaunchBoxGeometryWidenedAbutment();
   testArchDeckFillBoxCollidesSphereAboveCrown();
   testArchDeckFillBoxNoContactInsideOpenBore();

   return EXIT_SUCCESS;
}
