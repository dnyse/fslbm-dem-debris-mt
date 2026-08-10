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
//! \file PaperAnalyticalContactTest.cpp
//! Compares PaperAnalyticalContacts.h (closed-form Box x Box / Cylinder x Cylinder / Cylinder x
//! Plane algorithms taken directly from published papers) against SolidCylinderContactDetection's
//! EXISTING results (ParticleKernels.h, the production code actually used by the simulation) for
//! identical geometry.
//!
//! NOT wired into the simulation -- this file only exercises PaperAnalyticalContacts.h and, for
//! comparison, the already-existing SolidCylinderContactDetection, neither of which this test
//! modifies.
//!
//! Sign-convention reminder: PaperAnalyticalContacts.h reports `overlap` using the PAPERS' OWN
//! convention (positive = penetrating); SolidCylinderContactDetection reports `penetrationDepth`
//! using the OPPOSITE convention (negative = penetrating). Every comparison below negates one side
//! to match the other -- this is a deliberate, documented sign flip, not a bug in either file.
//!
//! Testing strategy: for geometries simple enough to hand-derive an exact expected value (axis-
//! aligned separations, single-axis rotations with a clean resulting geometry), both a hand-derived
//! literal AND a GJK/EPA cross-check are asserted. For geometries too intricate to hand-verify
//! reliably (general rotated-box orientations; the face-edge/band-edge/edge-edge cylinder scenarios,
//! which were flagged in PaperAnalyticalContacts.h as implemented-but-not-yet-tested and are the
//! least rigorously founded part of that file), ONLY the GJK/EPA cross-check is asserted -- treating
//! the already-proven-in-production GJK/EPA path as the ground truth oracle rather than risking a
//! second, independently-wrong hand derivation on already-shaky ground. The scenario that actually
//! fired is logged either way, for inspection.
//======================================================================================================================

#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"
#include "PaperAnalyticalContacts.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/Cylinder.h>
#include <mesa_pd/data/shape/HalfSpace.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/mpi/MPIManager.h>

namespace walberla
{
namespace antidunes
{
namespace paper_analytical
{
namespace test
{
using Vec3 = mesa_pd::Vec3;
using Mat3 = mesa_pd::Mat3;
using Rot3 = mesa_pd::Rot3;

//======================================================================================================================
// Helpers: run the existing, already-in-production SolidCylinderContactDetection for the same
// geometry, so every test below can cross-check against it without repeating the particle-storage
// boilerplate each time.
//======================================================================================================================
struct GjkResult
{
   bool   contact;
   real_t penetrationDepth;
   Vec3   normal;
   Vec3   point;
};

GjkResult gjkBoxBox(const Vec3& posA, const Rot3& rotA, const Vec3& edgeA, const Vec3& posB, const Rot3& rotB,
                     const Vec3& edgeB)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto body0 = ps->create();
   body0->setPosition(posA);
   body0->setRotation(rotA);
   auto body1 = ps->create();
   body1->setPosition(posB);
   body1->setRotation(rotB);
   const mesa_pd::data::Box boxA(edgeA), boxB(edgeB);

   SolidCylinderContactDetection gjk;
   const bool                    contactFound = gjk(0, 1, boxA, boxB, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

GjkResult gjkCylCyl(const Vec3& posI, const Rot3& rotI, real_t radiusI, real_t lengthI, const Vec3& posJ,
                     const Rot3& rotJ, real_t radiusJ, real_t lengthJ)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto body0 = ps->create();
   body0->setPosition(posI);
   body0->setRotation(rotI);
   auto body1 = ps->create();
   body1->setPosition(posJ);
   body1->setRotation(rotJ);
   const mesa_pd::data::Cylinder cylI(radiusI, lengthI), cylJ(radiusJ, lengthJ);

   SolidCylinderContactDetection gjk;
   const bool                    contactFound = gjk(0, 1, cylI, cylJ, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

GjkResult gjkCylPlane(const Vec3& posCyl, const Rot3& rotCyl, real_t radius, real_t length, const Vec3& planePos,
                       const Vec3& planeNormal)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto body = ps->create();
   body->setPosition(posCyl);
   body->setRotation(rotCyl);
   auto wall = ps->create();
   wall->setPosition(planePos);
   const mesa_pd::data::Cylinder  cylShape(radius, length);
   const mesa_pd::data::HalfSpace wallShape(planeNormal);

   SolidCylinderContactDetection gjk;
   const bool                    contactFound = gjk(0, 1, cylShape, wallShape, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

//======================================================================================================================
// Hybrid dispatcher: PaperAnalyticalContacts.h's cylinderCylinderContact(), with
// cylinderCylinderResultIsUncertain() deciding whether to trust it or fall back to the existing
// GJK/EPA path. Lives HERE (not in PaperAnalyticalContacts.h) specifically so that header stays free
// of any dependency on ParticleKernels.h / SolidCylinderContactDetection -- the detection LOGIC
// (cylinderCylinderResultIsUncertain) is pure and lives in the header; only the actual GJK/EPA CALL,
// which needs production types, lives here.
//======================================================================================================================
struct HybridResult
{
   bool        contact;
   real_t      overlap;  // paper convention throughout (positive = penetrating), even when the
                          // fallback fired -- GJK/EPA's penetrationDepth is negated on the way in, so
                          // callers never have to know which path actually answered.
   Vec3        normal;
   Vec3        point;
   const char* scenario; // analytical scenario name, or "GJK-EPA-fallback"
   bool        usedFallback;
};

HybridResult cylinderCylinderContactHybrid(const CylinderState& ci, const Rot3& rotI, const CylinderState& cj,
                                            const Rot3& rotJ)
{
   const ContactResult analytical = cylinderCylinderContact(ci, cj);

   if (!cylinderCylinderResultIsUncertain(ci, cj, analytical))
   {
      return HybridResult{ analytical.contact, analytical.overlap, analytical.normal,
                            analytical.point,   analytical.scenario, false };
   }

   const GjkResult gjk = gjkCylCyl(ci.pos, rotI, ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   return HybridResult{ gjk.contact, -gjk.penetrationDepth, gjk.normal, gjk.point, "GJK-EPA-fallback", true };
}

//======================================================================================================================
// Box x Box -- OBBTree SAT vs. the existing GJK/EPA fallback
//======================================================================================================================

//! Identical geometry to ContactDetectionTest.cpp's testBoxBoxContact: two unit cubes (half-extent
//! 0.5), centres 0.8 apart along x -> the unambiguous minimum-penetration axis is x, overlap 0.2.
void testBoxBoxSATMatchesExisting()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(0.8), real_t(0), real_t(0));
   const Mat3 R; // identity
   const Vec3 half(real_t(0.5), real_t(0.5), real_t(0.5));

   const ContactResult sat = boxBoxSAT(posA, R, half, posB, R, half);
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.point[0], real_t(0.4), real_t(1e-3));

   const GjkResult gjk = gjkBoxBox(posA, Rot3(), Vec3(real_t(1), real_t(1), real_t(1)), posB, Rot3(),
                                    Vec3(real_t(1), real_t(1), real_t(1)));
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.point[0], gjk.point[0], real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Box x Box (face-face, x): SAT overlap=" << sat.overlap << " (agree)");
}

void testBoxBoxSATNoContact()
{
   const Mat3           R;
   const Vec3           half(real_t(0.5), real_t(0.5), real_t(0.5));
   const ContactResult  sat = boxBoxSAT(Vec3(real_t(0), real_t(0), real_t(0)), R, half,
                                        Vec3(real_t(3), real_t(0), real_t(0)), R, half);
   WALBERLA_CHECK(!sat.contact);
}

//! Same as the x-axis case, offset along y instead -- catches an axis-specific indexing mistake the
//! x-only version couldn't.
void testBoxBoxSATFaceFaceYAxis()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(0), real_t(0.8), real_t(0));
   const Mat3 R;
   const Vec3 half(real_t(0.5), real_t(0.5), real_t(0.5));

   const ContactResult sat = boxBoxSAT(posA, R, half, posB, R, half);
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, Vec3(real_t(0), real_t(-1), real_t(0)), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.point[1], real_t(0.4), real_t(1e-3));

   const GjkResult gjk = gjkBoxBox(posA, Rot3(), Vec3(real_t(1), real_t(1), real_t(1)), posB, Rot3(),
                                    Vec3(real_t(1), real_t(1), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
}

//! Same again, offset along z.
void testBoxBoxSATFaceFaceZAxis()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(0), real_t(0), real_t(0.8));
   const Mat3 R;
   const Vec3 half(real_t(0.5), real_t(0.5), real_t(0.5));

   const ContactResult sat = boxBoxSAT(posA, R, half, posB, R, half);
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, Vec3(real_t(0), real_t(0), real_t(-1)), real_t(1e-3));

   const GjkResult gjk = gjkBoxBox(posA, Rot3(), Vec3(real_t(1), real_t(1), real_t(1)), posB, Rot3(),
                                    Vec3(real_t(1), real_t(1), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
}

//! Different-sized boxes (non-cubic half-extents), still an unambiguous single-axis contact.
void testBoxBoxSATAsymmetricSizes()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(0.7), real_t(0), real_t(0));
   const Mat3 R;
   const Vec3 halfA(real_t(0.3), real_t(0.4), real_t(0.5));
   const Vec3 halfB(real_t(0.6), real_t(0.2), real_t(0.3));
   // overlap along x = 0.3+0.6-0.7 = 0.2; along y = 0.4+0.2-0 = 0.6; along z = 0.5+0.3-0 = 0.8 -> x wins.

   const ContactResult sat = boxBoxSAT(posA, R, halfA, posB, R, halfB);
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.point[0], real_t(0.2), real_t(1e-3)); // midpoint(0.3, 0.7-0.6=0.1)=0.2

   const GjkResult gjk = gjkBoxBox(posA, Rot3(), real_t(2) * halfA, posB, Rot3(), real_t(2) * halfB);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.point[0], gjk.point[0], real_t(1e-3));
}

//! Corner of a 45-degree-rotated box poking into an axis-aligned box's flat face -- a genuinely
//! different SAT branch (an edge-cross axis or a face axis depending on exact geometry) than the
//! flush face-face cases above. Too intricate to fully hand-verify which of the 15 axes wins by
//! inspection alone, so this relies on the GJK/EPA cross-check as the actual verification.
void testBoxBoxSATRotatedCornerIntoFace()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(1.15), real_t(0), real_t(0));
   const Mat3 identity;
   const Vec3 half(real_t(0.5), real_t(0.5), real_t(0.5));

   Rot3 rotB;
   rotB.rotate(Vec3(real_t(0), real_t(0), real_t(1)), math::pi * real_t(0.25)); // 45 degrees about z

   const ContactResult sat = boxBoxSAT(posA, identity, half, posB, rotB.getMatrix(), half);
   const GjkResult      gjk = gjkBoxBox(posA, Rot3(), real_t(2) * half, posB, rotB, real_t(2) * half);

   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact)
   {
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-2));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-2));
      WALBERLA_LOG_INFO_ON_ROOT("Box x Box (rotated corner into face): SAT overlap=" << sat.overlap
                                                                                      << " normal=" << sat.normal
                                                                                      << " (agrees with GJK/EPA)");
   }
}

//! Two GENERALLY oriented boxes (both rotated, about different axes) -- the case SAT's 15-axis test
//! is actually FOR. No hand-derived expectation; GJK/EPA cross-check only.
void testBoxBoxSATBothRotatedGeneralOrientation()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(0.9), real_t(0.2), real_t(0.1));
   const Vec3 halfA(real_t(0.5), real_t(0.4), real_t(0.6));
   const Vec3 halfB(real_t(0.4), real_t(0.5), real_t(0.3));

   Rot3 rotA;
   rotA.rotate(Vec3(real_t(1), real_t(0), real_t(0)), real_t(0.5)); // ~28.6 degrees about x
   Rot3 rotB;
   rotB.rotate(Vec3(real_t(0), real_t(1), real_t(0)), real_t(-0.35)); // ~-20 degrees about y

   const ContactResult sat = boxBoxSAT(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(), halfB);
   const GjkResult      gjk = gjkBoxBox(posA, rotA, real_t(2) * halfA, posB, rotB, real_t(2) * halfB);

   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact)
   {
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-2));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-2));
   }
   WALBERLA_LOG_INFO_ON_ROOT("Box x Box (both rotated, general orientation): contact="
                              << sat.contact << " overlap=" << sat.overlap << " (agrees with GJK/EPA)");
}

//! Deep interpenetration (box B's centre lies well inside box A) -- SAT should handle an arbitrarily
//! large overlap the same uniform way as a shallow one, no special-casing.
void testBoxBoxSATDeepOverlap()
{
   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB(real_t(0.5), real_t(0), real_t(0));
   const Mat3 R;
   const Vec3 half(real_t(1), real_t(1), real_t(1));

   const ContactResult sat = boxBoxSAT(posA, R, half, posB, R, half);
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(1.5), real_t(1e-3)); // 1+1-0.5
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.point[0], real_t(0.25), real_t(1e-3));

   const GjkResult gjk = gjkBoxBox(posA, Rot3(), real_t(2) * half, posB, Rot3(), real_t(2) * half);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-2));
}

//! Boundary cases: exactly touching (gap == 0, must still register as contact) and just barely
//! separated (must NOT).
void testBoxBoxSATBoundaryTouching()
{
   const Mat3 R;
   const Vec3 half(real_t(0.5), real_t(0.5), real_t(0.5));

   const ContactResult touching =
      boxBoxSAT(Vec3(real_t(0), real_t(0), real_t(0)), R, half, Vec3(real_t(1), real_t(0), real_t(0)), R, half);
   WALBERLA_CHECK(touching.contact); // gap == 0 exactly, not < 0 -- must still count as contact
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(touching.overlap, real_t(0), real_t(1e-6));

   const ContactResult separated = boxBoxSAT(Vec3(real_t(0), real_t(0), real_t(0)), R, half,
                                              Vec3(real_t(1.0001), real_t(0), real_t(0)), R, half);
   WALBERLA_CHECK(!separated.contact);
}

//======================================================================================================================
// Cylinder x Cylinder -- Kodam et al. (2010) vs. the existing GJK/EPA fallback
//======================================================================================================================

//! Identical geometry to ContactDetectionTest.cpp's testCylinderCylinderContact.
void testCylinderCylinderParallelMatchesExisting()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   CylinderState cj{ Vec3(real_t(0.8), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "band-band-parallel");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[0], real_t(0.4), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, Rot3(), cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[0], gjk.point[0], real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (band-band-parallel): overlap=" << kodam.overlap << " (agree)");
}

void testCylinderCylinderParallelNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   CylinderState cj{ Vec3(real_t(3), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   WALBERLA_CHECK(!cylinderCylinderContact(ci, cj).contact);
}

//! Asymmetric radii/lengths, still band-band-parallel -- catches a hard-coded-equal-radius mistake.
void testCylinderCylinderParallelAsymmetric()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.3), real_t(3) };
   CylinderState cj{ Vec3(real_t(0.65), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(1) };
   // radial overlap = 0.3+0.5-0.65 = 0.15; axial overlap = 1.5+0.5-0 = 2.0 -> radial (band-band) wins.

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "band-band-parallel");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.15), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, Rot3(), cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Two cylinders stacked along their shared axis -> face-face branch.
void testCylinderCylinderFaceFace()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.5)), Mat3(), real_t(0.5), real_t(2) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-face");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.5), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(0), real_t(0), real_t(-1)), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, Rot3(), cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

void testCylinderCylinderFaceFaceNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(3)), Mat3(), real_t(0.5), real_t(2) };
   WALBERLA_CHECK(!cylinderCylinderContact(ci, cj).contact);
}

//! Skewed axes (i along z, j along x), sides touching. No equivalent in ContactDetectionTest.cpp at
//! all -- GJK/EPA was never previously exercised with non-parallel cylinder axes.
void testCylinderCylinderBandBandSkewed90Degrees()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };

   Rot3 rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5)); // body z -> world x
   CylinderState cj{ Vec3(real_t(0), real_t(0.8), real_t(1.5)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cj.axis(), Vec3(real_t(1), real_t(0), real_t(0)), real_t(1e-6));

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "band-band-skewed");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.2), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(0), real_t(-1), real_t(0)), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-2));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-2));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (skewed, 90deg): overlap=" << kodam.overlap << " (agree)");
}

//! Same idea at a genuinely oblique angle (60 degrees, not the "clean" 90) -- exercises the general
//! Vega & Lago closest-line-points formula rather than a right-angle special case.
void testCylinderCylinderBandBandSkewed60Degrees()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.4), real_t(4) };

   Rot3 rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(3)); // 60 degrees about y
   const Vec3 axisJ = rotJ.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1));

   // Place cj so its centreline passes near ci's axis at a small offset, in the plane containing
   // both cross(axisI, axisJ) directions -- offset purely along y (perpendicular to both the z-axis
   // and the xz-plane axisJ lies in), analogous to the 90-degree case above.
   CylinderState cj{ Vec3(real_t(0), real_t(0.6), real_t(1)), rotJ.getMatrix(), real_t(0.4), real_t(4) };
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(axisJ, cj.axis(), real_t(1e-6));

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk =
      gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK_EQUAL(kodam.contact, gjk.contact);
   if (kodam.contact)
   {
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-2));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-2));
   }
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (skewed, 60deg): scenario=" << kodam.scenario << " overlap="
                                                                               << kodam.overlap << " (agrees with GJK/EPA)");
}

void testCylinderCylinderBandBandSkewedNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
   CylinderState cj{ Vec3(real_t(0), real_t(3), real_t(1.5)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
   WALBERLA_CHECK(!cylinderCylinderContact(ci, cj).contact);
}

//! Cylinder j's curved side pressed into cylinder i's flat top face, with j's axis perpendicular to
//! i's (eq. 34's precondition) and the contact line lying fully within i's face (eq. 42's case).
//! Hand-derived: overlap = (halfLi + radiusJ) - |zIJ| = (1+0.5) - 1.4 = 0.1.
void testCylinderCylinderFaceBand()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5)); // body z -> world x
   CylinderState cj{ Vec3(real_t(0), real_t(0.45), real_t(1.4)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cj.axis(), Vec3(real_t(1), real_t(0), real_t(0)), real_t(1e-6));

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-band");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.1), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(0), real_t(0), real_t(-1)), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-2));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-2));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (face-band): overlap=" << kodam.overlap << " (agree)");
}

void testCylinderCylinderFaceBandNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
   CylinderState cj{ Vec3(real_t(0), real_t(0.45), real_t(2.5)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
   WALBERLA_CHECK(!cylinderCylinderContact(ci, cj).contact);
}

//! Cylinder j tilted 45 degrees to cylinder i, positioned so j's bottom rim region is near i's flat
//! face. This test EXISTS TO DOCUMENT A REAL, DISCOVERED LIMITATION, not to demonstrate agreement:
//! kodam.overlap=0.299 vs. GJK/EPA's 0.235 (27% off), and GJK/EPA's contact normal comes out
//! PERPENDICULAR to cylinder j's own axis -- a band-surface signature, not a rim-edge one.
//!
//! RETESTED after implementing Guo, Wassgren, Ketterhagen, Hancock & Curtis (2012)'s corrections
//! (see PaperAnalyticalContacts.h): UNCHANGED, still 0.299 vs 0.235. This confirms the diagnosis --
//! Guo et al.'s corrections fix MISCLASSIFICATION between rim-based candidate points (their E for
//! face-edge, P for band-edge/special-contact), not the deeper issue here, where the true closest
//! feature isn't a rim point AT ALL but a point on j's curved band surface away from either rim.
//! Neither Kodam et al.'s original paper nor Guo et al.'s corrections search that surface as a
//! candidate -- this is a genuine gap beyond what either published source addresses. Only CONTACT
//! EXISTENCE is asserted here for that reason: a boolean mismatch would mean a real missed/spurious
//! contact, but this specific magnitude disagreement is a documented, unresolved limitation, not this
//! test's job to hide.
void testCylinderCylinderTiltedRimIntoFace()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.6), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.25)); // 45 degrees about y
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.55)), rotJ.getMatrix(), real_t(0.2), real_t(2) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK_EQUAL(kodam.contact, gjk.contact);
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (tilted, oblique -- KNOWN LIMITATION CASE): kodam scenario="
                              << kodam.scenario << " overlap=" << kodam.overlap
                              << " vs GJK/EPA penetrationDepth=" << gjk.penetrationDepth
                              << " -- contact existence agrees; magnitude does not (see comment above)");
}

//! z=2.5 (only ~1 unit from ci's face) turned out to still be in actual contact according to
//! GJK/EPA (on j's BAND, not its rim -- the same discovered limitation as the test above), which
//! kodam's rim-only search incorrectly reported as separated: a genuine false negative, not just a
//! magnitude error. Moved much further away here so this specific test is unambiguous regardless of
//! that limitation; the false-negative failure mode itself is exactly what the comment above (and
//! PaperAnalyticalContacts.h's face-edge/band-edge/edge-edge section comment) documents.
void testCylinderCylinderTiltedRimIntoFaceNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.6), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.25));
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(6)), rotJ.getMatrix(), real_t(0.2), real_t(2) };
   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(!kodam.contact);
   WALBERLA_CHECK(!gjk.contact);
}

//! Cylinder j tilted obliquely toward cylinder i's curved side -> intended to exercise band-edge,
//! but hits the SAME rim-vs-band limitation as the test above: kodam's rim-only search reported NO
//! contact here while GJK/EPA (searching j's full surface) found one. RETESTED after implementing
//! Guo et al. (2012)'s corrections (bidirectional rimVsFaceOrBand, the corrected face-edge/
//! special-contact/band-edge dispatch): UNCHANGED, still a complete miss. Same conclusion as the
//! test above -- this is the true-closest-feature-is-on-the-band case neither Kodam et al. nor Guo
//! et al. search for, not a scenario-boundary classification bug their corrections target. Only
//! contact existence is logged, not asserted equal, since it's already known and documented to
//! disagree for this class of geometry.
void testCylinderCylinderTiltedRimIntoBand()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(0.3)); // ~54 degrees about x
   CylinderState cj{ Vec3(real_t(0.75), real_t(0), real_t(1.2)), rotJ.getMatrix(), real_t(0.25), real_t(1.5) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (tilted rim toward band -- KNOWN LIMITATION CASE): kodam.contact="
                              << kodam.contact << " scenario=" << kodam.scenario
                              << " | gjk.contact=" << gjk.contact
                              << " penetrationDepth=" << gjk.penetrationDepth);
}

void testCylinderCylinderTiltedRimIntoBandNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(0.3));
   CylinderState cj{ Vec3(real_t(1.5), real_t(0), real_t(1.2)), rotJ.getMatrix(), real_t(0.25), real_t(1.5) };
   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(!kodam.contact);
   WALBERLA_CHECK(!gjk.contact); // confirm this is a TRUE negative, not another instance of the
                                  // rim-vs-band limitation masking a real contact
}

//! Both cylinders' rim EDGES facing each other (neither a flat face nor a curved side is the closest
//! feature on either body) -> edge-edge-ish (actually resolved as face-edge, see below). Originally
//! a "KNOWN LIMITATION" case: with Kodam et al.'s ORIGINAL (uncorrected) dispatch logic and a
//! one-directional (j-rim-vs-i only) search, this reported NO CONTACT (scenario="none") for geometry
//! GJK/EPA confirms is genuinely touching (penetrationDepth=-0.05) -- a false negative.
//!
//! FIXED by applying Guo et al. (2012)'s corrections: rimVsFaceOrBand is now checked in BOTH
//! directions (i's rim against j's face/band, not just j's rim against i's), and this geometry turns
//! out to resolve cleanly as a face-edge contact from the i-rim-vs-j-face direction once that second
//! direction is actually checked -- overlap now matches GJK/EPA EXACTLY (0.05 vs 0.05). Restored to a
//! full assertion accordingly. (Not every case in this file's history was this fully resolved --
//! testCylinderCylinderTiltedRimIntoFace/IntoBand below are still open; see their comments for why
//! this fix doesn't reach them.)
void testCylinderCylinderEdgeEdge()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.4), real_t(1.6) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5)); // body z -> world x
   CylinderState cj{ Vec3(real_t(1.15), real_t(0), real_t(1.1)), rotJ.getMatrix(), real_t(0.4), real_t(1.6) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-2));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (edge-edge-ish, FIXED by Guo et al. corrections): scenario="
                              << kodam.scenario << " overlap=" << kodam.overlap
                              << " (was a false negative before the bidirectional check + corrections)");
}

//======================================================================================================================
// Hybrid dispatcher: does cylinderCylinderResultIsUncertain() correctly route the two still-open
// "known limitation" geometries to GJK/EPA (recovering the right answer), while leaving every
// already-verified scenario family on the fast analytical path (no unnecessary fallback calls)?
//======================================================================================================================

//! TIER 1: identical geometry to testCylinderCylinderTiltedRimIntoBand, where the analytical
//! dispatch finds NOTHING (scenario="none") despite the bodies genuinely touching. The bounding-sum
//! gate (centre distance vs. combined radius+half-length "reach") should recognise these bodies are
//! close enough to warrant a GJK/EPA double-check, and the fallback should recover the correct answer.
void testHybridTier1TriggersOnMissedContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(0.3));
   CylinderState cj{ Vec3(real_t(0.75), real_t(0), real_t(1.2)), rotJ.getMatrix(), real_t(0.25), real_t(1.5) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
   const GjkResult     gjk    = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(hybrid.usedFallback);
   WALBERLA_CHECK_EQUAL(hybrid.contact, gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 1 (missed contact): usedFallback=" << hybrid.usedFallback
                              << " contact=" << hybrid.contact << " overlap=" << hybrid.overlap
                              << " -- recovers the contact the pure-analytical path missed");
}

//! TIER 1 must NOT fire for a genuinely, unambiguously separated pair -- otherwise the fast path
//! would be worthless (every non-contact would pay for a GJK/EPA call). Centres 10 units apart,
//! reach only ~3 -- nowhere near the bounding-sum gate.
void testHybridTier1DoesNotTriggerWhenClearlySeparated()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
   CylinderState cj{ Vec3(real_t(10), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, Rot3());
   WALBERLA_CHECK(!hybrid.usedFallback);
   WALBERLA_CHECK(!hybrid.contact);
}

//! TIER 1 at the exact boundary of the bounding-sum gate: testCylinderCylinderParallelNoContact's
//! geometry has centreDist == reach EXACTLY (3.0 == 0.5+0.5+1+1), i.e. the gate's `<` (not `<=`)
//! should leave it untriggered -- worth pinning down explicitly since it's the closest a genuinely
//! separated pair gets to accidentally tripping Tier 1 in this test suite.
void testHybridTier1BoundaryCase()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   CylinderState cj{ Vec3(real_t(3), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, Rot3());
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 1 boundary: usedFallback=" << hybrid.usedFallback
                              << " (centreDist == reach exactly)");
   WALBERLA_CHECK(!hybrid.contact);
}

//! TIER 2: identical geometry to testCylinderCylinderTiltedRimIntoFace, where the analytical
//! dispatch DOES find a scenario (face-edge) but with the wrong magnitude (0.299 vs GJK/EPA's 0.235).
//! The circumscribing-sphere cross-check should flag this and recover the correct value via fallback.
void testHybridTier2TriggersOnWrongMagnitude()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.6), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.25));
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.55)), rotJ.getMatrix(), real_t(0.2), real_t(2) };

   const ContactResult analyticalOnly = cylinderCylinderContact(ci, cj); // for the "before" log line
   const HybridResult   hybrid        = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
   const GjkResult       gjk           = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(hybrid.usedFallback);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 2 (wrong magnitude): analytical-only overlap="
                              << analyticalOnly.overlap << ", hybrid (with fallback) overlap=" << hybrid.overlap
                              << " -- corrected from the wrong analytical value to GJK/EPA's " << -gjk.penetrationDepth);
}

//! TIER 2 must NOT fire for the scenario families already fully verified against GJK/EPA (face-face,
//! band-band-parallel, band-band-skewed, face-band), even though they involve tilted/oblique
//! geometry too -- these are excluded from the circumscribing-sphere check by scenario name, not by
//! re-deriving correctness each time. Reuses testCylinderCylinderFaceFace/BandBandSkewed60Degrees'
//! geometries as representative cases.
void testHybridTier2DoesNotTriggerOnVerifiedScenarios()
{
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.5)), Mat3(), real_t(0.5), real_t(2) };
      const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, Rot3());
      WALBERLA_CHECK(!hybrid.usedFallback);
      WALBERLA_CHECK_EQUAL(std::string(hybrid.scenario), "face-face");
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      CylinderState cj{ Vec3(real_t(0.8), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, Rot3());
      WALBERLA_CHECK(!hybrid.usedFallback);
      WALBERLA_CHECK_EQUAL(std::string(hybrid.scenario), "band-band-parallel");
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.4), real_t(4) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(3));
      CylinderState cj{ Vec3(real_t(0), real_t(0.6), real_t(1)), rotJ.getMatrix(), real_t(0.4), real_t(4) };
      const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
      WALBERLA_CHECK(!hybrid.usedFallback);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
      CylinderState cj{ Vec3(real_t(0), real_t(0.45), real_t(1.4)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
      const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
      WALBERLA_CHECK(!hybrid.usedFallback);
      WALBERLA_CHECK_EQUAL(std::string(hybrid.scenario), "face-band");
   }
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 2: none of face-face/band-band-parallel/band-band-skewed/"
                              "face-band triggered an unnecessary fallback call");
}

//! TIER 2, the empirically interesting case: the NOW-FIXED edge-edge-ish geometry resolves as
//! face-edge with an EXACT match to GJK/EPA (0.05 == 0.05) -- but face-edge IS one of the three
//! scenario names subject to the circumscribing-sphere check. Does the (necessarily loose) bound
//! over-trigger here anyway, spending an unnecessary GJK/EPA call on an already-correct answer? This
//! is reported, not asserted either way -- both outcomes are "correct" (the final number is right
//! regardless), this just characterises how tight the Tier 2 heuristic actually is in practice.
void testHybridTier2OnAlreadyCorrectFaceEdgeCase()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.4), real_t(1.6) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
   CylinderState cj{ Vec3(real_t(1.15), real_t(0), real_t(1.1)), rotJ.getMatrix(), real_t(0.4), real_t(1.6) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
   const GjkResult     gjk    = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(hybrid.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-2));
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 2 on already-correct face-edge case: usedFallback="
                              << hybrid.usedFallback
                              << " (either answer is 'correct' here -- this just characterises the "
                                 "heuristic's tightness, see comment)");
}

//======================================================================================================================
// Cylinder x Plane -- Kodam et al. face/band/edge-plane vs. the existing Cylinder x HalfSpace overload
//======================================================================================================================

//! Identical geometry to ContactDetectionTest.cpp's testCylinderHalfSpaceContact.
void testCylinderPlaneFaceMatchesExisting()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(1)), Mat3(), real_t(1), real_t(4) };
   const Vec3    planePos(real_t(0), real_t(0), real_t(0));
   const Vec3    planeNormal(real_t(0), real_t(0), real_t(1));

   const ContactResult kodam = cylinderPlaneContact(cyl, planePos, planeNormal);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-plane");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(1), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, planeNormal, real_t(1e-6));

   const GjkResult gjk = gjkCylPlane(cyl.pos, Rot3(), cyl.radius, cyl.length, planePos, planeNormal);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-6));

   // Contact POINT deliberately NOT asserted equal -- genuine convention difference, not a bug. See
   // the original version of this test / conversation history: Kodam's eq. 6 places the point at the
   // MIDPOINT of the overlap region; SolidCylinderContactDetection places it ON the wall surface.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(0), real_t(0), real_t(-0.5)), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(gjk.point, Vec3(real_t(0), real_t(0), real_t(0)), real_t(1e-6));
}

void testCylinderPlaneFaceNoContact()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(3)), Mat3(), real_t(1), real_t(4) };
   WALBERLA_CHECK(
      !cylinderPlaneContact(cyl, Vec3(real_t(0), real_t(0), real_t(0)), Vec3(real_t(0), real_t(0), real_t(1)))
          .contact);
}

//! Band-plane: cylinder lying on its side on the floor.
void testCylinderPlaneBand()
{
   Rot3 rot;
   rot.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(0.5)); // body z -> world y (lying flat)
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0.9)), rot.getMatrix(), real_t(1), real_t(4) };
   const Vec3    planePos(real_t(0), real_t(0), real_t(0));
   const Vec3    planeNormal(real_t(0), real_t(0), real_t(1));

   const ContactResult kodam = cylinderPlaneContact(cyl, planePos, planeNormal);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "band-plane");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.1), real_t(1e-3)); // radius 1 - |0.9| = 0.1
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, planeNormal, real_t(1e-3));

   const GjkResult gjk = gjkCylPlane(cyl.pos, rot, cyl.radius, cyl.length, planePos, planeNormal);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

void testCylinderPlaneBandNoContact()
{
   Rot3 rot;
   rot.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(0.5));
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(2)), rot.getMatrix(), real_t(1), real_t(4) };
   WALBERLA_CHECK(
      !cylinderPlaneContact(cyl, Vec3(real_t(0), real_t(0), real_t(0)), Vec3(real_t(0), real_t(0), real_t(1)))
          .contact);
}

//! Edge-plane: cylinder tilted 30 degrees from vertical, touching the floor at its bottom rim edge.
//! Hand-derived using the CORRECTED eq. 16-17 formula (see PaperAnalyticalContacts.h's edge-plane
//! comment -- an earlier version of both this derivation and the code it checks had the radial
//! direction and a sign wrong, caught when this test's original expected value of 0.116 disagreed
//! with GJK/EPA's actual 0.416): radialAwayFromPlane normalizes to (0.8660, 0, -0.5), giving rim
//! point E = (-0.2402, 0, 0.484 - 0.9) = (-0.2402, 0, -0.416) at zc = 0.6, so zE = -0.416, opposite
//! sign from zCyl = 0.6 -> contact, overlap = 0.416.
void testCylinderPlaneEdge()
{
   Rot3 rot;
   rot.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(6)); // 30 degrees about y
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0.6)), rot.getMatrix(), real_t(0.3), real_t(2) };
   const Vec3    planePos(real_t(0), real_t(0), real_t(0));
   const Vec3    planeNormal(real_t(0), real_t(0), real_t(1));

   const ContactResult kodam = cylinderPlaneContact(cyl, planePos, planeNormal);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "edge-plane");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.416), real_t(2e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, planeNormal, real_t(1e-3));

   const GjkResult gjk = gjkCylPlane(cyl.pos, rot, cyl.radius, cyl.length, planePos, planeNormal);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(2e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Plane (edge, NEW coverage): overlap=" << kodam.overlap << " (agree)");
}

void testCylinderPlaneEdgeNoContact()
{
   Rot3 rot;
   rot.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(6));
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(1.5)), rot.getMatrix(), real_t(0.3), real_t(2) };
   WALBERLA_CHECK(
      !cylinderPlaneContact(cyl, Vec3(real_t(0), real_t(0), real_t(0)), Vec3(real_t(0), real_t(0), real_t(1)))
          .contact);
}

} // namespace test
} // namespace paper_analytical
} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes::paper_analytical::test;

   // Box x Box
   testBoxBoxSATMatchesExisting();
   testBoxBoxSATNoContact();
   testBoxBoxSATFaceFaceYAxis();
   testBoxBoxSATFaceFaceZAxis();
   testBoxBoxSATAsymmetricSizes();
   testBoxBoxSATRotatedCornerIntoFace();
   testBoxBoxSATBothRotatedGeneralOrientation();
   testBoxBoxSATDeepOverlap();
   testBoxBoxSATBoundaryTouching();

   // Cylinder x Cylinder
   testCylinderCylinderParallelMatchesExisting();
   testCylinderCylinderParallelNoContact();
   testCylinderCylinderParallelAsymmetric();
   testCylinderCylinderFaceFace();
   testCylinderCylinderFaceFaceNoContact();
   testCylinderCylinderBandBandSkewed90Degrees();
   testCylinderCylinderBandBandSkewed60Degrees();
   testCylinderCylinderBandBandSkewedNoContact();
   testCylinderCylinderFaceBand();
   testCylinderCylinderFaceBandNoContact();
   testCylinderCylinderTiltedRimIntoFace();
   testCylinderCylinderTiltedRimIntoFaceNoContact();
   testCylinderCylinderTiltedRimIntoBand();
   testCylinderCylinderTiltedRimIntoBandNoContact();
   testCylinderCylinderEdgeEdge();

   // Hybrid dispatcher (analytical + GJK/EPA fallback via cylinderCylinderResultIsUncertain)
   testHybridTier1TriggersOnMissedContact();
   testHybridTier1DoesNotTriggerWhenClearlySeparated();
   testHybridTier1BoundaryCase();
   testHybridTier2TriggersOnWrongMagnitude();
   testHybridTier2DoesNotTriggerOnVerifiedScenarios();
   testHybridTier2OnAlreadyCorrectFaceEdgeCase();

   // Cylinder x Plane
   testCylinderPlaneFaceMatchesExisting();
   testCylinderPlaneFaceNoContact();
   testCylinderPlaneBand();
   testCylinderPlaneBandNoContact();
   testCylinderPlaneEdge();
   testCylinderPlaneEdgeNoContact();

   return EXIT_SUCCESS;
}
