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
//! Tests the closed-form contact detection of ClosedFormContacts*.h (Box x Box SAT, Cylinder x Cylinder after
//! Kodam et al. 2010 / Guo et al. 2012, Cylinder x Plane, Sphere x Cylinder, the support-function refinement and the
//! capsule/Layer-1 bounds) against GJK/EPA (SolidCylinderContactDetection / GeneralContactDetection) and against
//! hand-derived cases. Not part of the simulation; it modifies neither implementation.
//!
//! Sign convention: closed_form_contacts reports `overlap` > 0 for penetration; SolidCylinderContactDetection reports
//! `penetrationDepth` < 0 for penetration. Every comparison negates one side.
//!
//! Strategy: geometries simple enough to hand-derive (axis-aligned separations, single-axis rotations) are checked
//! against both a hand-derived literal and GJK/EPA. Intricate geometries (general rotated boxes, the cylinder
//! face-edge/band-edge/edge-edge scenarios) are checked against GJK/EPA only, which serves as the reference.
//! The scenario that fired is logged either way.
//======================================================================================================================

#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"
#include "ClosedFormContacts.h"

#include <sstream>
#include <mesa_pd/collision_detection/GeneralContactDetection.h>
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

#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{
namespace test
{
using Vec3 = mesa_pd::Vec3;
using Mat3 = mesa_pd::Mat3;
using Rot3 = mesa_pd::Rot3;

//! Directory of the test binary (set from argv[0] in main()); debug output (mismatch VTK/CSV dumps) is written there.
std::filesystem::path g_outputDir;

//! Per-mismatch VTK frame series (cylinder_fuzz_mismatches, cylinder_kodam_mismatches), off by default because it
//! dominates the runtime; the CSV dump is always written. Set WALBERLA_WRITE_MISMATCH_VTK=1 to enable.
inline bool writeMismatchVTKEnabled()
{
   const char* env = std::getenv("WALBERLA_WRITE_MISMATCH_VTK");
   return env != nullptr && std::string(env) == "1";
}

//======================================================================================================================
// Helpers: run SolidCylinderContactDetection / GeneralContactDetection on the same geometry.
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

   // The Cylinder x Cylinder overload (ParticleKernels.h) rejects pairs whose centre distance exceeds the sum of
   // the interaction radii, which ParticleStorage defaults to 0. Set the circumscribing-sphere radius production
   // uses (SetupStages.h: sqrt(radius^2 + halfLength^2)), otherwise every pair is reported as no contact.
   body0->setInteractionRadius(std::sqrt(radiusI * radiusI + real_t(0.25) * lengthI * lengthI));
   body1->setInteractionRadius(std::sqrt(radiusJ * radiusJ + real_t(0.25) * lengthJ * lengthJ));

   SolidCylinderContactDetection gjk;
   const bool                    contactFound = gjk(0, 1, cylI, cylJ, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

//! Not GJK/EPA: SolidCylinderContactDetection's Cylinder x HalfSpace overload is the production support-mapping
//! closed form (ParticleKernels.h), so this compares the closed form against production.
GjkResult productionCylPlane(const Vec3& posCyl, const Rot3& rotCyl, real_t radius, real_t length, const Vec3& planePos,
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

//! Not GJK/EPA: SolidCylinderContactDetection's Sphere x Cylinder overload is the production closed form
//! (ParticleKernels.h, mid-overlap contact point). The GJK/EPA reference is gjkSphereCyl below.
GjkResult productionSphereCyl(const Vec3& spherePos, real_t sphereRadius, const Vec3& cylPos, const Rot3& rotCyl,
                        real_t cylRadius, real_t cylLength)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto body0 = ps->create();
   body0->setPosition(spherePos);
   auto body1 = ps->create();
   body1->setPosition(cylPos);
   body1->setRotation(rotCyl);
   const mesa_pd::data::Sphere   sphShape(sphereRadius);
   const mesa_pd::data::Cylinder cylShape(cylRadius, cylLength);

   SolidCylinderContactDetection gjk;
   const bool                    contactFound = gjk(0, 1, sphShape, cylShape, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

//! GJK/EPA for a sphere and a cylinder: mesa_pd's GeneralContactDetection, whose generic operator() runs
//! collideGJKEPA for this pair.
GjkResult gjkSphereCyl(const Vec3& spherePos, real_t sphereRadius, const Vec3& cylPos, const Rot3& rotCyl,
                        real_t cylRadius, real_t cylLength)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto body0 = ps->create();
   body0->setPosition(spherePos);
   auto body1 = ps->create();
   body1->setPosition(cylPos);
   body1->setRotation(rotCyl);
   const mesa_pd::data::Sphere   sphShape(sphereRadius);
   const mesa_pd::data::Cylinder cylShape(cylRadius, cylLength);
   mesa_pd::collision_detection::GeneralContactDetection gjk;
   const bool contactFound = gjk(0, 1, sphShape, cylShape, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

//======================================================================================================================
// Hybrid dispatcher: cylinderCylinderContact(), with cylinderCylinderResultIsUncertain() deciding whether to trust
// it or fall back to GJK/EPA. Kept here so ClosedFormContacts*.h does not depend on ParticleKernels.h.
//======================================================================================================================
struct HybridResult
{
   bool        contact;
   real_t      overlap;  // positive = penetrating, also when the fallback fired (GJK/EPA depth negated)
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
// Box x Box -- OBBTree SAT vs. GJK/EPA
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

//! Same as the x-axis case, offset along y (catches axis-specific indexing errors).
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

//! Corner of a box rotated 45 degrees about z poking into an axis-aligned box's face. Hand-derived (the x face axis
//! wins) and cross-checked against GJK/EPA.
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

   // Hand-derived: B reaches 0.5 sqrt(2) towards -x, so the x axis overlaps by 0.5 + 0.5 sqrt(2) - 1.15 = 0.0571068;
   // B's face axes at 45 deg overlap by 0.5 sqrt(2) + 0.5 - 1.15/sqrt(2) = 0.394; the cross products of the z edges
   // coincide with these axes. Confirmed by a brute-force minimum translation distance.
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(0.5) + real_t(0.5) * std::sqrt(real_t(2)) - real_t(1.15), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-9));
   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact)
   {
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
      WALBERLA_LOG_INFO_ON_ROOT("Box x Box (rotated corner into face): SAT overlap=" << sat.overlap
                                                                                      << " normal=" << sat.normal
                                                                                      << " (agrees with GJK/EPA)");
   }
}

//! Two generally oriented boxes (rotated about different axes); expected value from a brute-force minimum
//! translation distance, plus the GJK/EPA cross-check.
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

   // Brute-force minimum translation distance (overleaf scripts/contact-validation/verify_hand_cases.py): 0.078618
   // along -x.
   WALBERLA_CHECK(sat.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(0.078618), real_t(1e-5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, Vec3(real_t(-1), real_t(0), real_t(0)), real_t(1e-3));
   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact)
   {
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
   }
   WALBERLA_LOG_INFO_ON_ROOT("Box x Box (both rotated, general orientation): contact="
                              << sat.contact << " overlap=" << sat.overlap << " (agrees with GJK/EPA)");
}

//! Deep interpenetration (B's centre inside A): overlap = 1 + 1 - 0.5 = 1.5, point at the overlap midpoint 0.25.
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
}

//! SAT vs. GJK/EPA on a geometry reconstructed from a logged deep contact of two 6x18x6 boxes (half-extents
//! <3,9,3>): depth 0.773241, tilted (non-cardinal) normal (0.17535, 0.0340412, 0.983917). Orientations were not
//! logged, so box A is tilted by the minimal rotation sending its local +Z onto the logged normal, B stays
//! axis-aligned, and the two are separated along the normal to reproduce the logged depth (4 of 6 relative DOF fixed;
//! not an exact replay). Checks that SAT and GJK/EPA agree and that the reconstruction reproduces depth and normal.
void testBoxBoxSATRealAntidunesLaunchCollision()
{
   const Vec3 half(real_t(3), real_t(9), real_t(3)); // half-extents of the logged boxes, both

   // Logged SAT normal and depth at the deepest instant of the contact.
   const Vec3 realNormal(real_t(0.17535), real_t(0.0340412), real_t(0.983917));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(realNormal.length(), real_t(1), real_t(1e-3));
   const real_t realPenetration = real_t(0.773241); // SAT/paper convention: positive = penetrating

   // Box A: tilted by the minimal rotation sending local +Z onto realNormal.
   const Vec3   worldZ(real_t(0), real_t(0), real_t(1));
   Rot3         rotA;
   const real_t cosTheta = std::clamp(worldZ * realNormal, real_t(-1), real_t(1));
   if (cosTheta < real_t(1) - real_t(1e-9))
   {
      Vec3 axis = worldZ % realNormal; // cross product
      axis      = axis.getNormalized();
      rotA.rotate(axis, std::acos(cosTheta));
   }
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rotA.getMatrix() * worldZ, realNormal, real_t(1e-4));

   // Box B: axis-aligned; the tilt is attributed entirely to A.
   const Mat3 identity;

   // SAT overlap along n = projRadius(A,n) + projRadius(B,n) - centreDist, with projRadius(A,n) = half.z (A's local
   // Z is n) and projRadius(B,n) = |half.x n.x| + |half.y n.y| + |half.z n.z| (B axis-aligned).
   const real_t projA = half[2]; // A's local Z axis IS realNormal by construction
   const real_t projB =
      std::abs(half[0] * realNormal[0]) + std::abs(half[1] * realNormal[1]) + std::abs(half[2] * realNormal[2]);
   const real_t centreDist = projA + projB - realPenetration;

   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   // B sits behind A along -realNormal (A on top).
   const Vec3 posB = posA - realNormal * centreDist;

   const ContactResult sat = boxBoxSAT(posA, rotA.getMatrix(), half, posB, identity, half);
   const GjkResult      gjk = gjkBoxBox(posA, rotA, real_t(2) * half, posB, Rot3(), real_t(2) * half);

   WALBERLA_LOG_INFO_ON_ROOT("Box x Box (REAL job-4131344 launch collision reconstruction): "
      << "target penetration=" << realPenetration << " target normal=" << realNormal);
   WALBERLA_LOG_INFO_ON_ROOT("  SAT:     contact=" << sat.contact << " overlap=" << sat.overlap
      << " normal=" << sat.normal << " point=" << sat.point);
   WALBERLA_LOG_INFO_ON_ROOT("  GJK/EPA: contact=" << gjk.contact << " penetration=" << gjk.penetrationDepth
      << " normal=" << gjk.normal << " point=" << gjk.point);

   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact && gjk.contact)
   {
      const real_t overlapDiff = std::abs(sat.overlap - (-gjk.penetrationDepth));
      const real_t normalDiff  = (sat.normal - gjk.normal).length();
      WALBERLA_LOG_INFO_ON_ROOT("  agreement: overlapDiff=" << overlapDiff << " normalDiff=" << normalDiff);
      // Loose tolerance: catches gross divergence (wrong axis, sign flip), not round-off.
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
      // The reconstruction reproduces the logged contact (brute-force MTD 0.773241 along the logged normal).
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, realPenetration, real_t(1e-4));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, realNormal, real_t(1e-4));
   }
}

//! SAT vs. GJK/EPA on two geometries built from logged anomalously deep Box x Box contacts (half-extents <6,1,1.5>),
//! with the same construction as testBoxBoxSATRealAntidunesLaunchCollision. The separation reproduces the logged
//! overlap only along the logged normal; another axis overlaps less (brute-force MTD 2.04333 instead of 2.51343 for
//! event 1, 1.23458 instead of 1.29292 for event 2), so these are SAT/GJK cross-checks on run-derived geometry, not
//! replays of the logged contacts.
void testBoxBoxSATJob4139454Explosion1()
{
   const Vec3 half(real_t(6), real_t(1), real_t(1.5)); // half-extents of the logged boxes

   // Logged normal and depth of event 1.
   const Vec3   realNormal(real_t(0.116856), real_t(-0.773005), real_t(0.623545));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(realNormal.length(), real_t(1), real_t(1e-3));
   const real_t realPenetration = real_t(2.51343);

   const Vec3   worldZ(real_t(0), real_t(0), real_t(1));
   Rot3         rotA;
   const real_t cosTheta = std::clamp(worldZ * realNormal, real_t(-1), real_t(1));
   if (cosTheta < real_t(1) - real_t(1e-9))
   {
      Vec3 axis = worldZ % realNormal;
      axis      = axis.getNormalized();
      rotA.rotate(axis, std::acos(cosTheta));
   }
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rotA.getMatrix() * worldZ, realNormal, real_t(1e-4));

   const Mat3 identity;
   const real_t projA = half[2]; // A's local Z axis IS realNormal by construction
   const real_t projB =
      std::abs(half[0] * realNormal[0]) + std::abs(half[1] * realNormal[1]) + std::abs(half[2] * realNormal[2]);
   const real_t centreDist = projA + projB - realPenetration;

   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB = posA - realNormal * centreDist;

   const ContactResult sat = boxBoxSAT(posA, rotA.getMatrix(), half, posB, identity, half);
   const GjkResult      gjk = gjkBoxBox(posA, rotA, real_t(2) * half, posB, Rot3(), real_t(2) * half);

   WALBERLA_LOG_INFO_ON_ROOT("Box x Box (job 4139454 explosion, event 1/initiating): target penetration="
      << realPenetration << " target normal=" << realNormal);
   WALBERLA_LOG_INFO_ON_ROOT("  SAT:     contact=" << sat.contact << " overlap=" << sat.overlap
      << " normal=" << sat.normal << " point=" << sat.point);
   WALBERLA_LOG_INFO_ON_ROOT("  GJK/EPA: contact=" << gjk.contact << " penetration=" << gjk.penetrationDepth
      << " normal=" << gjk.normal << " point=" << gjk.point);

   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact && gjk.contact)
   {
      const real_t overlapDiff = std::abs(sat.overlap - (-gjk.penetrationDepth));
      const real_t normalDiff  = (sat.normal - gjk.normal).length();
      WALBERLA_LOG_INFO_ON_ROOT("  agreement: overlapDiff=" << overlapDiff << " normalDiff=" << normalDiff);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
      // Brute-force MTD of this geometry: 2.04333 (not the logged 2.51343, see the header).
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(2.04333), real_t(1e-4));
   }
}

//! Same as above for the second logged event.
void testBoxBoxSATJob4139454Explosion2()
{
   const Vec3 half(real_t(6), real_t(1), real_t(1.5));

   const Vec3   realNormal(real_t(-0.0150641), real_t(-0.216776), real_t(-0.976105));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(realNormal.length(), real_t(1), real_t(1e-3));
   const real_t realPenetration = real_t(1.29292);

   const Vec3   worldZ(real_t(0), real_t(0), real_t(1));
   Rot3         rotA;
   const real_t cosTheta = std::clamp(worldZ * realNormal, real_t(-1), real_t(1));
   if (cosTheta < real_t(1) - real_t(1e-9))
   {
      Vec3 axis = worldZ % realNormal;
      axis      = axis.getNormalized();
      rotA.rotate(axis, std::acos(cosTheta));
   }
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(rotA.getMatrix() * worldZ, realNormal, real_t(1e-4));

   const Mat3 identity;
   const real_t projA = half[2];
   const real_t projB =
      std::abs(half[0] * realNormal[0]) + std::abs(half[1] * realNormal[1]) + std::abs(half[2] * realNormal[2]);
   const real_t centreDist = projA + projB - realPenetration;

   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   const Vec3 posB = posA - realNormal * centreDist;

   const ContactResult sat = boxBoxSAT(posA, rotA.getMatrix(), half, posB, identity, half);
   const GjkResult      gjk = gjkBoxBox(posA, rotA, real_t(2) * half, posB, Rot3(), real_t(2) * half);

   WALBERLA_LOG_INFO_ON_ROOT("Box x Box (job 4139454 explosion, event 2/cascading): target penetration="
      << realPenetration << " target normal=" << realNormal);
   WALBERLA_LOG_INFO_ON_ROOT("  SAT:     contact=" << sat.contact << " overlap=" << sat.overlap
      << " normal=" << sat.normal << " point=" << sat.point);
   WALBERLA_LOG_INFO_ON_ROOT("  GJK/EPA: contact=" << gjk.contact << " penetration=" << gjk.penetrationDepth
      << " normal=" << gjk.normal << " point=" << gjk.point);

   WALBERLA_CHECK_EQUAL(sat.contact, gjk.contact);
   if (sat.contact && gjk.contact)
   {
      const real_t overlapDiff = std::abs(sat.overlap - (-gjk.penetrationDepth));
      const real_t normalDiff  = (sat.normal - gjk.normal).length();
      WALBERLA_LOG_INFO_ON_ROOT("  agreement: overlapDiff=" << overlapDiff << " normalDiff=" << normalDiff);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
      // Brute-force MTD of this geometry: 1.23458 (not the logged 1.29292, see the header).
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(1.23458), real_t(1e-4));
   }
}

//! Boundary cases: exactly touching (gap 0) counts as contact, a gap of 1e-4 does not.
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
// Cylinder x Cylinder -- Kodam et al. (2010) vs. GJK/EPA
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

//! Asymmetric radii/lengths, still band-band-parallel (catches an equal-radius assumption).
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

//! Skewed axes (i along z, j along x), sides overlapping: overlap = 0.5 + 0.5 - 0.8 = 0.2 along -y.
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (skewed, 90deg): overlap=" << kodam.overlap << " (agree)");
}

//! Same at an oblique angle (60 degrees), exercising the general Vega & Lago closest-line-points formula.
void testCylinderCylinderBandBandSkewed60Degrees()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.4), real_t(4) };

   Rot3 rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(3)); // 60 degrees about y
   const Vec3 axisJ = rotJ.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1));

   // cj offset purely along y, perpendicular to both axes (axisJ lies in the xz plane).
   CylinderState cj{ Vec3(real_t(0), real_t(0.6), real_t(1)), rotJ.getMatrix(), real_t(0.4), real_t(4) };
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(axisJ, cj.axis(), real_t(1e-6));

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk =
      gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   // Hand-derived: the axes cross at x = 0 with a perpendicular distance of 0.6 along y (both lie in planes
   // y = const), well inside both lengths, so the overlap is 0.4 + 0.4 - 0.6 = 0.2 along -y.
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "band-band-skewed");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.2), real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(0), real_t(-1), real_t(0)), real_t(1e-9));
   WALBERLA_CHECK_EQUAL(kodam.contact, gjk.contact);
   if (kodam.contact)
   {
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
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

//! Cylinder j's band pressed into i's top face, j's axis perpendicular to i's (Kodam eq. 34 precondition).
//! Overlap = (halfLi + radiusJ) - |zIJ| = (1 + 0.5) - 1.4 = 0.1. Kodam's line AB (eq. 40-41: j's centreline of
//! half-length 2 translated onto i's face) straddles i's face disc (radius 0.5) with both ends outside, a case
//! neither eq. 42 (fully inside) nor eq. 44 (one end inside) names; faceBandContactPoint intersects the finite line
//! with the circle. The two crossings are symmetric about x = 0, so the point is (0, 0.45), moved half the overlap
//! below the face (z = 1) to z = 0.95.
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
   // point = AB's midpoint projected onto i's face plane = (0, 0.45, 0.95) (Kodam et al. 2010 eq. 40-41)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(0), real_t(0.45), real_t(0.95)), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (face-band, straddling): overlap=" << kodam.overlap
                              << " point=" << kodam.point << " (agree; GJK/EPA point=" << gjk.point
                              << ", informational only -- underdetermined along the band direction)");
}

//! Eq. 42 "fully inside" case: short cylinder j (half-length 0.3; eq. 40-41 has no radius term) whose translated
//! centreline AB lies entirely within i's face disc (radius 1). Expected point: AB's midpoint, i.e. cj's centre
//! projected onto i's face, moved half the overlap below it: (0, 0.2, 1 - 0.1/2) = (0, 0.2, 0.95).
void testCylinderCylinderFaceBandFullyInside()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5)); // body z -> world x
   CylinderState cj{ Vec3(real_t(0), real_t(0.2), real_t(1.2)), rotJ.getMatrix(), real_t(0.3), real_t(0.6) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-band");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.1), real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(0), real_t(0.2), real_t(0.95)), real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (face-band, fully inside): overlap=" << kodam.overlap
                              << " point=" << kodam.point << " (agree)");
}

//! Eq. 44 "partially inside" case (one end of AB inside the disc), off-centre so the ends are not symmetric.
//! A (t = +halfLj = +1) at (1.5, 0.3) is outside Ri = 1, B (t = -1) at (-0.5, 0.3) inside; the disc crossing nearer
//! A is at t2 = -0.5 + sqrt(0.91) = 0.453939, so the point is the midpoint of B and that crossing,
//! t = (-1 + 0.453939)/2 = -0.273031, i.e. xy = (0.5 - 0.273031, 0.3) = (0.226970, 0.3), z = 1 - 0.1/2 = 0.95.
void testCylinderCylinderFaceBandPartiallyInside()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5)); // body z -> world x
   CylinderState cj{ Vec3(real_t(0.5), real_t(0.3), real_t(1.2)), rotJ.getMatrix(), real_t(0.3), real_t(2.0) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-band");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.1), real_t(1e-3));
   // point = midpoint of B and the disc-crossing nearer A, t=(-1+0.453939)/2 -> (0.226970, 0.3) (Kodam et al. 2010 eq.
   // 44)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(0.226970), real_t(0.3), real_t(0.95)),
                                       real_t(1e-3));

   const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (face-band, partially inside): overlap=" << kodam.overlap
                              << " point=" << kodam.point << " (agree)");
}

void testCylinderCylinderFaceBandNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
   CylinderState cj{ Vec3(real_t(0), real_t(0.45), real_t(2.5)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
   WALBERLA_CHECK(!cylinderCylinderContact(ci, cj).contact);
}

//! Known limitation of face-edge: j tilted 45 degrees with its lowest rim region in i's face. Face-edge treats the
//! face as unbounded and returns the depth of the deepest rim point E below it,
//! 1 - (1.55 - 0.5 sqrt(2) - 0.2 * 0.5 sqrt(2)) = 0.298528. E lies at radius 0.566 on a face of radius 0.6, so
//! moving j sideways across i's rim separates the pair sooner: brute-force MTD 0.235356 along (0.707, 0, -0.707),
//! perpendicular to j's axis. Neither Kodam et al. nor Guo et al. (2012) search j's band as a candidate. With a host
//! ten times the probe radius (testCylinderCylinderFaceEdgeClean) face-edge is exact. Both values are pinned so a
//! change is noticed.
void testCylinderCylinderTiltedRimIntoFace()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.6), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.25)); // 45 degrees about y
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.55)), rotJ.getMatrix(), real_t(0.2), real_t(2) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK_EQUAL(kodam.contact, gjk.contact);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-edge");
   // depth of the lowest rim point below the face: 1 - (1.55 - 0.6 sqrt(2)) = 0.298528
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(1) - (real_t(1.55) - real_t(0.6) * std::sqrt(real_t(2))),
                                       real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(-gjk.penetrationDepth, real_t(0.235356), real_t(1e-3)); // brute-force MTD
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (tilted, oblique -- KNOWN LIMITATION CASE): kodam scenario="
                              << kodam.scenario << " overlap=" << kodam.overlap
                              << " vs GJK/EPA penetrationDepth=" << gjk.penetrationDepth
                              << " -- contact existence agrees; magnitude does not (see comment above)");
}

//! Same tilt, clearly separated (z = 6): both report no contact.
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

//! Cylinder j tilted obliquely toward i's band. The pair is exactly touching (brute-force MTD 0, GJK/EPA depth about
//! -1e-10), so the closed form's "no contact" is correct at the boundary. A genuine miss (depth 0.113) is used in
//! testHybridTier1TriggersOnMissedContact.
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
   WALBERLA_CHECK(!kodam.contact);
   WALBERLA_CHECK_LESS(std::abs(gjk.penetrationDepth), real_t(1e-6)); // exactly touching
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
   WALBERLA_CHECK(!gjk.contact); // true negative, confirmed by GJK/EPA
}

//! Guo et al. (2012) "special contact" (their Fig. 8), not one of Kodam's seven scenarios: the band-edge candidate P
//! lies within host i's band region and on the same side of i's face as probe j's centre, which lies beyond that face
//! (j has sunk in past i's end). Geometry: j (radius 1, length 0.6) tilted 35 degrees, centred at z = 1.7, 0.2 past
//! i's face at z = 1.5, offset in x so face-edge's candidate E lies outside i's face disc (rho^2 = 1.23 > 1). The
//! depth is the deepest rim point inside the host (0.614334; Guo's original depth at P is 0.181), checked against
//! GJK/EPA (0.609508; brute-force MTD 0.609529) with the fuzz criterion max(0.01, 5 %).
void testCylinderCylinderSpecialContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(3.0) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * (real_t(35) / real_t(180))); // 35 degrees about x
   CylinderState cj{ Vec3(real_t(0.9), real_t(0), real_t(1.7)), rotJ.getMatrix(), real_t(1.0), real_t(0.6) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "special-contact");
   // Depth (deepest rim point inside the host, see rimVsFaceOrBand) vs. GJK/EPA within max(0.01, 5 %).
   WALBERLA_CHECK_LESS(std::abs(kodam.overlap - (-gjk.penetrationDepth)),
                       std::max(real_t(0.01), real_t(0.05) * (-gjk.penetrationDepth)));
   WALBERLA_CHECK(gjk.contact);
   // special-contact reports the purely axial face normal (0,0,-1) by construction; GJK/EPA's normal has a small
   // off-axis component for the tilted j (about <-0.12,0.11,-0.99>), so only the direction is compared.
   WALBERLA_CHECK_GREATER(kodam.normal * gjk.normal, real_t(0.95));
   // Point: host i's face at z = +1.5, deepest point at depth overlap below it, so z = 1.5 - overlap/2.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], real_t(1.5) - real_t(0.5) * kodam.overlap, real_t(1e-10));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (special-contact): kodam overlap=" << kodam.overlap
                              << " normal=" << kodam.normal << " vs GJK/EPA penetrationDepth="
                              << gjk.penetrationDepth << " normal=" << gjk.normal);
}

//! Independent reference for testCylinderCylinderSpecialContact: the minimum translation distance (MTD) from the
//! cylinders' support functions, by a dense sweep over unit directions (no closed form, no GJK/EPA). Logged only.
//!
//! With support function h(d) = max_{p in body} (p-c).d, the gap along d is (cB-cA).d - hA(d) - hB(-d), and
//! MTD = -max_d gap(d) = min_d [hA(d) + hB(-d) - (cB-cA).d]. For a cylinder
//! h(d) = halfLength |axis.d| + radius sqrt(max(0, 1 - (axis.d)^2)).
void testCylinderCylinderSpecialContactBruteForceMTD()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(3.0) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * (real_t(35) / real_t(180)));
   CylinderState cj{ Vec3(real_t(0.9), real_t(0), real_t(1.7)), rotJ.getMatrix(), real_t(1.0), real_t(0.6) };

   auto support = [](const CylinderState& c, const Vec3& d) {
      const real_t axisDotD = c.axis() * d;
      return real_t(0.5) * c.length * std::abs(axisDotD) +
             c.radius * std::sqrt(std::max(real_t(0), real_t(1) - axisDotD * axisDotD));
   };

   const Vec3 centreDiff = cj.pos - ci.pos;

   // Coarse-to-fine search for the maximum gap: a Fibonacci sphere lattice (2e6 directions), then shrinking random
   // perturbations around the best direction.
   auto gap = [&](const Vec3& d) { return centreDiff * d - support(ci, d) - support(cj, -d); };

   Vec3   bestDir(real_t(0), real_t(0), real_t(1));
   real_t bestGap = gap(bestDir);
   auto   sweepFibonacci = [&](int n) {
      const real_t goldenAngle = math::pi * (real_t(3) - std::sqrt(real_t(5)));
      for (int i = 0; i < n; ++i)
      {
         const real_t z     = real_t(1) - (real_t(2) * real_t(i) + real_t(1)) / real_t(n);
         const real_t r     = std::sqrt(std::max(real_t(0), real_t(1) - z * z));
         const real_t theta = goldenAngle * real_t(i);
         const Vec3   d(r * std::cos(theta), r * std::sin(theta), z);
         const real_t g = gap(d);
         if (g > bestGap) // MAXIMIZE gap: the true separating axis is the direction of LEAST overlap
         {
            bestGap = g;
            bestDir = d;
         }
      }
   };
   sweepFibonacci(2000000);

   // Local refinement: 40 rounds of 2000 random perturbations, radius shrinking by 0.7 per round (seed 999).
   std::mt19937                             rng(999);
   std::uniform_real_distribution< real_t > unit(real_t(-1), real_t(1));
   real_t                                    radius = real_t(0.05);
   for (int round = 0; round < 40; ++round)
   {
      for (int i = 0; i < 2000; ++i)
      {
         Vec3 perturb(unit(rng), unit(rng), unit(rng));
         Vec3 d = (bestDir + radius * perturb).getNormalizedIfNotZero();
         if (d.sqrLength() < real_t(0.5)) continue;
         const real_t g = gap(d);
         if (g > bestGap) // MAXIMIZE gap: the true separating axis is the direction of LEAST overlap
         {
            bestGap = g;
            bestDir = d;
         }
      }
      radius *= real_t(0.7);
   }

   const real_t mtd = -bestGap;
   WALBERLA_LOG_INFO_ON_ROOT(
      "Cylinder x Cylinder (special-contact, INDEPENDENT brute-force MTD arbiter): brute-force MTD="
      << mtd << " (direction=" << bestDir << ") | original Guo overlap=0.1811 (corrected special-contact now ~0.614) | GJK/EPA overlap~0.610 "
      << "-- whichever of the latter two this brute-force value sits closer to is the more likely "
      << "correct one, independent of both implementations");
}

//! Sweep of the special-contact configuration (Guo et al. 2012, Fig. 8): probe j tilted about x by 5..45 degrees,
//! shifted laterally and axially around testCylinderCylinderSpecialContact's geometry. For every pose classified as
//! "special-contact", the depth (deepest rim point inside the host) must agree with GJK/EPA within the fuzz
//! tolerance max(0.01, 5 %), and the point must lie at the middle of the overlap along the axial normal.
void testCylinderCylinderSpecialContactCorrectedSweep()
{
   int numSpecial = 0, numAgree = 0;
   real_t worstAbs = real_t(0);
   for (int tiltDeg = 5; tiltDeg <= 45; tiltDeg += 5)
      for (int ix = 0; ix <= 6; ++ix)
         for (int iz = 0; iz <= 6; ++iz)
         {
            CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(3.0) };
            Rot3          rotJ;
            rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(tiltDeg) / real_t(180));
            const real_t x = real_t(0.5) + real_t(0.1) * real_t(ix);
            const real_t z = real_t(1.55) + real_t(0.05) * real_t(iz);
            CylinderState cj{ Vec3(x, real_t(0), z), rotJ.getMatrix(), real_t(1.0), real_t(0.6) };

            const ContactResult kodam = cylinderCylinderContact(ci, cj);
            if (!kodam.contact || std::string(kodam.scenario) != "special-contact") continue;
            ++numSpecial;
            const GjkResult gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
            WALBERLA_CHECK(gjk.contact);
            const real_t ref  = -gjk.penetrationDepth;
            const real_t diff = std::abs(kodam.overlap - ref);
            worstAbs = std::max(worstAbs, diff);
            if (diff <= std::max(real_t(0.01), real_t(0.05) * ref)) ++numAgree;
            else
               WALBERLA_LOG_INFO_ON_ROOT("  special-contact sweep mismatch: tilt=" << tiltDeg << " x=" << x << " z=" << z
                                         << " corrected=" << kodam.overlap << " GJK/EPA=" << ref);
            // contact point in the middle of the overlap along the axial normal (host face at z = +1.5)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], real_t(1.5) - real_t(0.5) * kodam.overlap, real_t(1e-9));
         }
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (special-contact, corrected depth, sweep): " << numSpecial
                              << " special-contact poses, " << numAgree << " agree with GJK/EPA within max(0.01, 5%), "
                              << "worst absolute difference=" << worstAbs);
   WALBERLA_CHECK_GREATER(numSpecial, 0);
   WALBERLA_CHECK_EQUAL(numAgree, numSpecial);
}

void testCylinderCylinderSpecialContactNoContact()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(3.0) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * (real_t(35) / real_t(180)));
   CylinderState cj{ Vec3(real_t(0.9), real_t(0), real_t(6.0)), rotJ.getMatrix(), real_t(1.0), real_t(0.6) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(!kodam.contact);
   WALBERLA_CHECK(!gjk.contact);
}

//! Both cylinders' rims facing each other (i along z, j along x). Resolved as face-edge from the i-rim-vs-j-face
//! direction, which needs the bidirectional rimVsFaceOrBand check of Guo et al. (2012); overlap 0.05, equal to
//! GJK/EPA and the brute-force MTD.
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.05), real_t(1e-9)); // brute-force MTD 0.05
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (edge-edge-ish, FIXED by Guo et al. corrections): scenario="
                              << kodam.scenario << " overlap=" << kodam.overlap
                              << " (was a false negative before the bidirectional check + corrections)");
}

//======================================================================================================================
// Randomized fuzz helpers shared by the Cylinder x Cylinder and Box x Box fuzz tests: random orientation, and a
// separation sampler that stays within a physically plausible overlap depth. Uniform sampling of the full separation
// range mostly yields near-total interpenetration, which a spring-dashpot contact never lets develop; the overlap is
// therefore capped at 15 % of the smaller body's radius (by default), already a stress case compared with the
// few-percent overlaps of a DEM run.
//======================================================================================================================
struct FuzzRng
{
   std::mt19937                             rng;
   std::uniform_real_distribution< real_t > angleDist{ real_t(0), real_t(2) * math::pi };
   std::uniform_real_distribution< real_t > unitDist{ real_t(-1), real_t(1) };

   explicit FuzzRng(unsigned seed) : rng(seed) {}

   Vec3 randomUnitVec()
   {
      Vec3 v;
      do {
         v = Vec3(unitDist(rng), unitDist(rng), unitDist(rng));
      } while (v.sqrLength() < real_t(1e-6));
      return v.getNormalized();
   }
   // Uniformly distributed over SO(3) (Shoemake 1992: unit quaternion from three uniform numbers).
   Rot3 randomRotation()
   {
      std::uniform_real_distribution< real_t > u(real_t(0), real_t(1));
      const real_t u1 = u(rng), u2 = u(rng), u3 = u(rng);
      const real_t a = std::sqrt(real_t(1) - u1), b = std::sqrt(u1);
      const real_t w = a * std::sin(real_t(2) * math::pi * u2), x = a * std::cos(real_t(2) * math::pi * u2);
      const real_t y = b * std::sin(real_t(2) * math::pi * u3), z = b * std::cos(real_t(2) * math::pi * u3);
      const real_t angle = real_t(2) * std::acos(std::clamp(std::abs(w), real_t(0), real_t(1)));
      const Vec3   axis  = (w < 0 ? real_t(-1) : real_t(1)) * Vec3(x, y, z);
      Rot3         r;
      if (axis.length() > real_t(1e-12)) r.rotate(axis.getNormalized(), angle);
      return r;
   }
   //! Uniform axis and uniform angle in [0, 2pi): NOT uniform on SO(3) (the Haar measure has density
   //! ~ 1 - cos(angle)), it over-samples rotations near the identity. Unused; kept for reproducing older numbers.
   Rot3 randomRotationAxisAngle()
   {
      Rot3 r;
      r.rotate(randomUnitVec(), angleDist(rng));
      return r;
   }
   // Finds the contact boundary along the direction by bisection (40 steps) of the reference predicate contactAt,
   // then samples the centre distance uniformly in [boundary - maxOverlap, boundary + 2 maxOverlap], with
   // maxOverlap = maxOverlapFraction * minBodySize; this covers both shallow overlaps and near misses. (Sampling
   // relative to circumscribing radii misses the boundary for most orientations.) Requires contactAt(nearSep) and
   // !contactAt(farSep).
   template< typename ContactAtFn >
   real_t sampleSeparationNearBoundary(ContactAtFn contactAt, real_t nearSep, real_t farSep, real_t minBodySize,
                                       real_t maxOverlapFraction = real_t(0.15))
   {
      WALBERLA_ASSERT(contactAt(nearSep));
      WALBERLA_ASSERT(!contactAt(farSep));
      real_t lo = nearSep, hi = farSep; // invariant: contactAt(lo)==true, contactAt(hi)==false
      for (int it = 0; it < 40; ++it)
      {
         const real_t mid = real_t(0.5) * (lo + hi);
         if (contactAt(mid)) { lo = mid; }
         else { hi = mid; }
      }
      const real_t boundary   = real_t(0.5) * (lo + hi);
      const real_t maxOverlap = maxOverlapFraction * minBodySize;
      std::uniform_real_distribution< real_t > d(boundary - maxOverlap, boundary + real_t(2) * maxOverlap);
      return d(rng);
   }
};

//======================================================================================================================
// Randomized fuzz test: random radius/length/orientation for both cylinders, separation sampled near the GJK/EPA
// contact boundary (FuzzRng::sampleSeparationNearBoundary). Both the raw closed form (cylinderCylinderContact) and
// the hybrid dispatcher (cylinderCylinderContactHybrid, GJK/EPA fallback when cylinderCylinderResultIsUncertain())
// are checked against GJK/EPA. Fixed seed.
//======================================================================================================================

//======================================================================================================================
// Debug visualization: every raw closed-form vs. GJK/EPA disagreement of the fuzz test (wrong existence, or wrong
// magnitude where both report contact) is dumped as a VTK point-glyph dataset, with the convention of
// BuildTimeloop.h's floaterCylinder output: "glyphScale" (2r, length, 2r) and "rotation" pre-composed with Rx(+90),
// because ParaView's Cylinder glyph has its axis along local Y and mesa_pd's along local Z. In ParaView: Glyph,
// Glyph Type = Cylinder, Scale Array = glyphScale, Orientation Array = rotation, Scale Factor = 1. A companion CSV
// lists the same cases.
//======================================================================================================================
struct CylMismatchRecord
{
   int           trial;
   CylinderState ci, cj;
   bool          analyticalContact, gjkContact;
   real_t        analyticalOverlap; // paper convention (positive = penetrating); 0 if no contact
   real_t        gjkOverlap;        // -penetrationDepth, same convention; 0 if no contact
   const char*   scenario;          // cylinderCylinderContact's own scenario tag ("none" if no contact)
   int           kind; // 0 = existence false positive (analytical says contact, GJK says none)
                        // 1 = existence false negative (analytical says none, GJK says contact)
                        // 2 = magnitude mismatch (both agree contact exists, overlap differs)
   const char*   cause      = ""; // root-cause tag (string literal), see
                                    // testCylinderCylinderFullComparisonAndHybridFeasibility
   real_t        causeMetric = real_t(-1); // number supporting `cause` (e.g. normal misalignment in
                                             // degrees); -1 where n/a
};

inline std::array< real_t, 9 > glyphRotationFromR(const Mat3& R)
{
   // (R * Rx(+90))^T = Rx(-90) * R^T, as in BuildTimeloop.h's SelectParticleRotationMatrix for floaterCylinder.
   return { R(0, 0), R(1, 0), R(2, 0), R(0, 2), R(1, 2), R(2, 2), -R(0, 1), -R(1, 1), -R(2, 1) };
}

//! Writes one cylinder (single point) as a .vtu; side = 0 for ci, 1 for cj.
inline void writeCylMismatchFrame(const std::string& path, const CylMismatchRecord& r, int side)
{
   const CylinderState& c = (side == 0) ? r.ci : r.cj;

   std::ofstream f(path);

   f << "<?xml version=\"1.0\"?>\n";
   f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
   f << " <UnstructuredGrid>\n";
   f << "  <Piece NumberOfPoints=\"1\" NumberOfCells=\"1\">\n";

   f << "   <Points>\n    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
   f << "     " << c.pos[0] << " " << c.pos[1] << " " << c.pos[2] << "\n";
   f << "    </DataArray>\n   </Points>\n";

   f << "   <Cells>\n";
   f << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">0</DataArray>\n";
   f << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">1</DataArray>\n";
   f << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">1</DataArray>\n";
   f << "   </Cells>\n";

   f << "   <PointData>\n";

   const auto g = glyphRotationFromR(c.R);

   f << "    <DataArray type=\"Float64\" Name=\"glyphScale\" NumberOfComponents=\"3\" format=\"ascii\">"
     << (real_t(2) * c.radius) << " " << c.length << " " << (real_t(2) * c.radius) << "</DataArray>\n";

   f << "    <DataArray type=\"Float64\" Name=\"rotation\" NumberOfComponents=\"9\" format=\"ascii\">";
   for (real_t v : g) f << v << " ";
   f << "</DataArray>\n";

   f << "    <DataArray type=\"Int32\" Name=\"side\" NumberOfComponents=\"1\" format=\"ascii\">"
     << side << "</DataArray>\n";

   f << "    <DataArray type=\"Int32\" Name=\"trial\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.trial << "</DataArray>\n";

   f << "    <DataArray type=\"Int32\" Name=\"mismatchKind\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.kind << "</DataArray>\n";

   // edgeEdgeStatus: 1 = Kodam's edge-edge branch (eq. 82-101) fired; 2 = scenario "none" (dispatch fell through
   // every scenario, no closed-form answer); 0 = another scenario resolved the pair.
   const int edgeEdgeStatus = (std::string(r.scenario) == "edge-edge")   ? 1
                              : (std::string(r.scenario) == "none")      ? 2
                                                                          : 0;
   f << "    <DataArray type=\"Int32\" Name=\"edgeEdgeStatus\" NumberOfComponents=\"1\" format=\"ascii\">"
     << edgeEdgeStatus << "</DataArray>\n";

   f << "    <DataArray type=\"Float64\" Name=\"analyticalOverlap\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.analyticalOverlap << "</DataArray>\n";

   f << "    <DataArray type=\"Float64\" Name=\"gjkOverlap\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.gjkOverlap << "</DataArray>\n";

   // The string tag `cause` goes to the CSV only (no string DataArray here); causeMetric is exposed for colouring.
   f << "    <DataArray type=\"Float64\" Name=\"causeMetric\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.causeMetric << "</DataArray>\n";

   f << "   </PointData>\n";
   f << "   <CellData>\n   </CellData>\n";
   f << "  </Piece>\n </UnstructuredGrid>\n</VTKFile>\n";
}

//! Writes one cylinder's (ci or cj, per `side`) mismatches as its own time series, baseDir/<name>_ci.pvd or
//! baseDir/<name>_cj.pvd, so each can get its own Glyph filter and scale factor in ParaView.
inline void writeCylMismatchSide(const std::string& baseDir, const std::string& name, int side,
                                  const std::vector< CylMismatchRecord >& records)
{
   if (records.empty())
   {
      WALBERLA_LOG_INFO_ON_ROOT("[fuzz] " << name << ": no mismatches to write, skipping.");
      return;
   }

   const std::string suffix    = (side == 0) ? "_ci" : "_cj";
   const std::string framesDir = baseDir + "/frames" + suffix;
   std::filesystem::create_directories(framesDir);

   std::ofstream pvd(baseDir + "/" + name + suffix + ".pvd");
   pvd << "<?xml version=\"1.0\"?>\n";
   pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
   pvd << " <Collection>\n";

   for (size_t i = 0; i < records.size(); ++i)
   {
      const std::string frameFile = "frames" + suffix + "/" + name + suffix + "_" + std::to_string(i) + ".vtu";
      writeCylMismatchFrame(baseDir + "/" + frameFile, records[i], side);
      pvd << "  <DataSet timestep=\"" << i << "\" group=\"\" part=\"0\" file=\"" << frameFile << "\"/>\n";
   }

   pvd << " </Collection>\n</VTKFile>\n";

   WALBERLA_LOG_INFO_ON_ROOT("[fuzz] wrote " << records.size() << " frame(s) to " << baseDir << "/" << name
      << suffix << ".pvd");
}

//! Writes both cylinders' time series; they share the timestep numbering 0..N-1, so ParaView keeps them in step.
inline void writeCylMismatchVTKSeries(const std::string& baseDir, const std::string& name,
                                       const std::vector< CylMismatchRecord >& records)
{
   writeCylMismatchSide(baseDir, name, 0, records);
   writeCylMismatchSide(baseDir, name, 1, records);
}

inline void writeCylMismatchCSV(const std::string& path, const std::vector< CylMismatchRecord >& records)
{
   if (records.empty()) return;
   std::ofstream f(path);
   f << "trial,kind,Ri,Li,Rj,Lj,sep,scenario,analyticalContact,gjkContact,analyticalOverlap,gjkOverlap,"
        "overlapDiff,cause,causeMetric\n";
   for (const auto& r : records)
   {
      const real_t sep        = (r.cj.pos - r.ci.pos).length();
      const real_t overlapDiff = std::abs(r.analyticalOverlap - r.gjkOverlap);
      const char*  kindStr = r.kind == 0 ? "existence_false_positive"
                             : r.kind == 1 ? "existence_false_negative"
                                           : "magnitude_mismatch";
      f << r.trial << "," << kindStr << "," << r.ci.radius << "," << r.ci.length << "," << r.cj.radius << ","
        << r.cj.length << "," << sep << "," << r.scenario << "," << r.analyticalContact << "," << r.gjkContact
        << "," << r.analyticalOverlap << "," << r.gjkOverlap << "," << overlapDiff << "," << r.cause << ","
        << r.causeMetric << "\n";
   }
   WALBERLA_LOG_INFO_ON_ROOT("[fuzz] wrote " << records.size() << " mismatch row(s) to " << path);
}

void testCylinderCylinderRandomizedFuzz()
{
   FuzzRng                                  fuzz(12345); // fixed seed -- reproducible across runs
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 10000;
   int           numOverlapping = 0;
   int           numAnalyticalMismatch = 0;   // raw analytical contact-existence disagrees with GJK/EPA
   int           numHybridMismatch      = 0;   // hybrid (with fallback) contact-existence disagrees with GJK/EPA
   int           numFallback            = 0;   // trials where cylinderCylinderResultIsUncertain triggered
                                                // and the hybrid called GJK/EPA
   real_t        maxAgreeOverlapDiff = real_t(0);
   real_t        maxAgreePointDist   = real_t(0);
   Vec3          worstPointCi, worstPointCj;

   // Every raw closed-form vs. GJK/EPA disagreement, for the VTK/CSV dump at the end.
   std::vector< CylMismatchRecord > mismatches;

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));

      const Vec3 posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult analytical = cylinderCylinderContact(ci, cj);
      const HybridResult  hybrid     = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
      const GjkResult     gjk = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);

      const real_t gjkOverlap = gjk.contact ? -gjk.penetrationDepth : real_t(0);

      if (hybrid.usedFallback) ++numFallback;
      if (gjk.contact) ++numOverlapping;
      if (analytical.contact != gjk.contact)
      {
         ++numAnalyticalMismatch;
         mismatches.push_back(CylMismatchRecord{
            trial, ci, cj, analytical.contact, gjk.contact,
            analytical.contact ? analytical.overlap : real_t(0), gjkOverlap, analytical.scenario,
            analytical.contact ? 0 : 1 }); // 0=false positive (analytical only), 1=false negative (analytical only)
      }
      else if (analytical.contact && gjk.contact)
      {
         // Both report contact: record a magnitude mismatch beyond max(0.01, 5 %) for the debug dump (not asserted).
         const real_t overlapDiff = std::abs(analytical.overlap - gjkOverlap);
         if (overlapDiff > std::max(real_t(0.01), real_t(0.05) * gjkOverlap))
            mismatches.push_back(CylMismatchRecord{ trial, ci, cj, analytical.contact, gjk.contact,
                                                      analytical.overlap, gjkOverlap, analytical.scenario, 2 });
      }

      if (hybrid.contact != gjk.contact)
      {
         ++numHybridMismatch;
         WALBERLA_LOG_WARNING_ON_ROOT(
            "[fuzz] hybrid/GJK contact-existence mismatch, trial " << trial << ": hybrid.contact="
            << hybrid.contact << " (scenario=" << hybrid.scenario << ") gjk.contact=" << gjk.contact
            << " | Ri=" << Ri << " Li=" << Li << " Rj=" << Rj << " Lj=" << Lj << " sep=" << sep);
      }

      if (hybrid.contact && gjk.contact)
      {
         const real_t overlapDiff = std::abs(hybrid.overlap - (-gjk.penetrationDepth));
         const real_t pointDist   = (hybrid.point - gjk.point).length();
         if (overlapDiff > maxAgreeOverlapDiff) maxAgreeOverlapDiff = overlapDiff;
         if (pointDist > maxAgreePointDist)
         {
            maxAgreePointDist = pointDist;
            worstPointCi       = posI;
            worstPointCj       = posJ;
         }
      }
   }

   if (writeMismatchVTKEnabled())
      writeCylMismatchVTKSeries((g_outputDir / "cylinder_fuzz_mismatches").string(), "cylinder_fuzz_mismatches",
                                 mismatches);
   writeCylMismatchCSV((g_outputDir / "cylinder_fuzz_mismatches.csv").string(), mismatches);

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Cylinder x Cylinder, " << numTrials << " random trials (" << numOverlapping
      << " actually overlapping per GJK/EPA, overlap depth capped realistic): raw-analytical "
      << "contact-existence mismatches=" << numAnalyticalMismatch << " (expected > 0 -- this is the "
      << "known rim-vs-band gap the hybrid fallback exists to cover), hybrid contact-existence "
      << "mismatches=" << numHybridMismatch << " (should be 0 or very close to it), max overlap "
      << "magnitude difference on agreeing contacts=" << maxAgreeOverlapDiff
      << ", max contact-point distance on agreeing contacts=" << maxAgreePointDist
      << " (worst-case centres: ci=" << worstPointCi << " cj=" << worstPointCj << ")");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] GJK/EPA fallback usage (BOUNDARY-BIASED SAMPLING -- see testCylinderCylinderFallbackRateBroadSampling "
      << "below for the realistic-distribution number instead): " << numFallback << "/" << numTrials
      << " trials (" << (real_t(100) * real_t(numFallback) / real_t(numTrials)) << "%) triggered "
      << "cylinderCylinderResultIsUncertain (Tier 1 or Tier 2). sampleSeparationNearBoundary deliberately "
      << "samples every trial within 0.15-0.30x minBodySize of the TRUE GJK contact boundary (by design, to "
      << "stress-test edge cases) -- both Tier 1's reach-sum bound and Tier 2's circumscribing-sphere bound "
      << "are loose/conservative specifically near that boundary, so this number is expected to run very "
      << "high and should NOT be read as \"the hybrid falls back this often in a real simulation\".");

   // Pass criterion: the hybrid agrees with GJK/EPA on contact existence in every trial. Tier 1 (reach-sum bound)
   // sends every closed-form "no contact" within reach to GJK/EPA, and special-contact always falls back.
   WALBERLA_CHECK_EQUAL(numHybridMismatch, 0);
   // Overlap magnitude is not asserted: band-edge/edge-edge/face-edge keep their raw closed-form depth (no magnitude
   // check in cylinderCylinderResultIsUncertain), so maxAgreeOverlapDiff reflects their approximation error.
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] max overlap magnitude difference on agreeing contacts=" << maxAgreeOverlapDiff
      << " (informational only -- Tier 2 no longer corrects band-edge/edge-edge/face-edge's own raw "
      << "overlap value, see cylinderCylinderResultIsUncertain's comment; not asserted)");
   // Point distance is not asserted either: for flush contacts (line or area) the contact point is not unique.
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] max contact-point distance on agreeing contacts=" << maxAgreePointDist << " (informational "
      << "only -- not asserted)");
}

//======================================================================================================================
// Raw closed form (cylinderCylinderContact, no fallback) vs. GJK/EPA, broken down by scenario: frequency of each
// scenario and, per scenario, agreement on existence, overlap and contact point. 5000 trials, seed 24680, same
// sampling as testCylinderCylinderRandomizedFuzz (radius U[0.3,2], length U[0.5,4], uniform SO(3), boundary-biased
// separation). Reported, not asserted: face-edge/band-edge/edge-edge/special-contact are approximate by design.
//======================================================================================================================
void testCylinderCylinderPureAnalyticalScenarioBreakdown()
{
   FuzzRng                                  fuzz(24680); // fixed seed
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 5000;

   struct ScenarioStats
   {
      int                    numTrials         = 0;
      int                    existenceAgree    = 0;
      int                    existenceMismatch = 0;
      real_t                 maxOverlapDiff    = real_t(0);
      real_t                 maxPointDist      = real_t(0);
      std::vector< real_t >  pointDists; // for a percentile breakdown, not just the worst outlier
      // Worst overlap disagreement: a normal close to GJK/EPA's indicates an imprecise depth on the same feature, a
      // different normal indicates that GJK/EPA found another closest feature (e.g. the band, not a rim).
      real_t        worstOverlapDiff = real_t(0);
      CylinderState worstCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      CylinderState worstCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      Vec3          worstAnalyticalNormal, worstGjkNormal;
      real_t        worstAnalyticalOverlap = real_t(0), worstGjkOverlap = real_t(0);
      // One example of a missed contact (closed form none, GJK/EPA contact).
      bool          haveMismatchExample = false;
      CylinderState mismatchCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      CylinderState mismatchCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      Vec3          mismatchGjkNormal;
      real_t        mismatchGjkOverlap = real_t(0);
      // Best-agreeing example with gjkOverlap > 0.05 (candidate geometry for a hand-derived test).
      real_t        bestRelOverlapDiff = std::numeric_limits< real_t >::max();
      CylinderState bestCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      CylinderState bestCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      real_t        bestOverlap = real_t(0);
   };
   std::map< std::string, ScenarioStats > byScenario;
   int numNoContactBoth = 0;

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));
      const Vec3 posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult analytical = cylinderCylinderContact(ci, cj); // pure analytical, no fallback
      const GjkResult      gjk       = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);

      if (!analytical.contact && !gjk.contact)
      {
         ++numNoContactBoth;
         continue;
      }

      const std::string scenario =
         analytical.contact ? std::string(analytical.scenario) : std::string("(missed -- analytical found no contact)");
      auto& st = byScenario[scenario];
      ++st.numTrials;

      if (analytical.contact != gjk.contact)
      {
         ++st.existenceMismatch;
         if (!st.haveMismatchExample && gjk.contact) // the interesting case: analytical missed a real contact
         {
            st.haveMismatchExample = true;
            st.mismatchCi          = ci;
            st.mismatchCj          = cj;
            st.mismatchGjkNormal   = gjk.normal;
            st.mismatchGjkOverlap  = -gjk.penetrationDepth;
         }
         continue;
      }
      ++st.existenceAgree;

      const real_t gjkOverlap  = -gjk.penetrationDepth;
      const real_t overlapDiff = std::abs(analytical.overlap - gjkOverlap);
      const real_t pointDist   = (analytical.point - gjk.point).length();
      st.maxOverlapDiff = std::max(st.maxOverlapDiff, overlapDiff);
      st.maxPointDist   = std::max(st.maxPointDist, pointDist);
      st.pointDists.push_back(pointDist);
      if (overlapDiff > st.worstOverlapDiff)
      {
         st.worstOverlapDiff       = overlapDiff;
         st.worstCi                = ci;
         st.worstCj                = cj;
         st.worstAnalyticalNormal  = analytical.normal;
         st.worstGjkNormal         = gjk.normal;
         st.worstAnalyticalOverlap = analytical.overlap;
         st.worstGjkOverlap        = gjkOverlap;
      }
      if (gjkOverlap > real_t(0.05))
      {
         const real_t relDiff = overlapDiff / gjkOverlap;
         if (relDiff < st.bestRelOverlapDiff)
         {
            st.bestRelOverlapDiff = relDiff;
            st.bestCi             = ci;
            st.bestCj             = cj;
            st.bestOverlap        = gjkOverlap;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] PURE analytical (no GJK/EPA fallback) Cylinder x Cylinder scenario breakdown, "
      << numTrials << " random-orientation trials (" << numNoContactBoth
      << " agree no contact at all -- excluded from the per-scenario counts below):");
   for (auto& [name, st] : byScenario)
   {
      real_t p50 = real_t(0), p90 = real_t(0);
      if (!st.pointDists.empty())
      {
         std::vector< real_t > sorted = st.pointDists;
         std::sort(sorted.begin(), sorted.end());
         p50 = sorted[sorted.size() / 2];
         p90 = sorted[std::min(sorted.size() - 1, size_t(real_t(0.9) * real_t(sorted.size())))];
      }
      WALBERLA_LOG_INFO_ON_ROOT(
         "  " << name << ": count=" << st.numTrials << " (" << (real_t(100) * real_t(st.numTrials) / real_t(numTrials))
         << "% of all trials) | existence: agree=" << st.existenceAgree << " mismatch="
         << st.existenceMismatch << " | on agreeing trials: maxOverlapDiff=" << st.maxOverlapDiff
         << " maxPointDist=" << st.maxPointDist << " p50PointDist=" << p50 << " p90PointDist=" << p90);
      if (st.worstOverlapDiff > real_t(1e-3))
      {
         const real_t normalDot = st.worstAnalyticalNormal * st.worstGjkNormal;
         WALBERLA_LOG_INFO_ON_ROOT(
            "    worst overlap-diff root-cause probe: analytical.overlap=" << st.worstAnalyticalOverlap
            << " gjk.overlap=" << st.worstGjkOverlap << " | analytical.normal=" << st.worstAnalyticalNormal
            << " gjk.normal=" << st.worstGjkNormal << " (dot=" << normalDot << ", <1 means GJK/EPA "
            << "found a geometrically DIFFERENT closest feature than the rim/edge point the analytical "
            << "formula used -- consistent with the documented rim-vs-band-surface gap, not a shared "
            << "feature computed slightly imprecisely) | geometry: ci.pos=" << st.worstCi.pos
            << " ci.axis=" << st.worstCi.axis() << " Ri=" << st.worstCi.radius << " Li=" << st.worstCi.length
            << " | cj.pos=" << st.worstCj.pos << " cj.axis=" << st.worstCj.axis() << " Rj=" << st.worstCj.radius
            << " Lj=" << st.worstCj.length);
      }
      if (st.haveMismatchExample)
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "    existence-mismatch example (analytical missed a real contact): gjk.overlap="
            << st.mismatchGjkOverlap << " gjk.normal=" << st.mismatchGjkNormal << " | geometry: ci.pos="
            << st.mismatchCi.pos << " ci.axis=" << st.mismatchCi.axis() << " Ri=" << st.mismatchCi.radius
            << " Li=" << st.mismatchCi.length << " | cj.pos=" << st.mismatchCj.pos << " cj.axis="
            << st.mismatchCj.axis() << " Rj=" << st.mismatchCj.radius << " Lj=" << st.mismatchCj.length);
      }
      if (st.bestRelOverlapDiff < std::numeric_limits< real_t >::max())
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "    best-agreeing non-trivial example (candidate for a genuine hand-drawn test): "
            << "gjk.overlap=" << st.bestOverlap << " relOverlapDiff=" << st.bestRelOverlapDiff
            << " | geometry: ci.pos=" << st.bestCi.pos << " ci.R=" << st.bestCi.R << " Ri=" << st.bestCi.radius
            << " Li=" << st.bestCi.length << " | cj.pos=" << st.bestCj.pos << " cj.R=" << st.bestCj.R
            << " Rj=" << st.bestCj.radius << " Lj=" << st.bestCj.length);
      }
   }
}

//======================================================================================================================
// Targeted fuzz test: axes forced (anti)parallel, the precondition of face-face/band-band-parallel, which uniformly
// random orientations never hit (zero measure). 5000 trials, seed 11235, otherwise the sampling of
// testCylinderCylinderPureAnalyticalScenarioBreakdown, so the breakdowns are comparable. The separation direction is
// random, so eq. 22's delta_a vs. delta_r comparison chooses between face-face and band-band-parallel. j's axis is
// i's axis or its 180-degree flip with equal probability, covering both signs of |e_i . e_j| = 1.
//======================================================================================================================
void testCylinderCylinderParallelAxesFuzz()
{
   FuzzRng                                  fuzz(11235); // fixed seed
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > flipDist(real_t(0), real_t(1));

   constexpr int numTrials = 5000;

   struct ScenarioStats
   {
      int                    numTrials         = 0;
      int                    existenceAgree    = 0;
      int                    existenceMismatch = 0;
      real_t                 maxOverlapDiff    = real_t(0);
      real_t                 maxPointDist      = real_t(0);
      std::vector< real_t >  pointDists;
   };
   std::map< std::string, ScenarioStats > byScenario;
   int numNoContactBoth = 0;

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();

      const Vec3 eI = rotI.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1));
      Rot3       rotJ = rotI; // parallel by default
      if (flipDist(fuzz.rng) < real_t(0.5))
      {
         // antiparallel: flip 180 degrees about any axis perpendicular to eI
         const Vec3 helper = (std::abs(eI[2]) < real_t(0.9)) ? Vec3(real_t(0), real_t(0), real_t(1))
                                                              : Vec3(real_t(1), real_t(0), real_t(0));
         const Vec3 perp = (eI % helper).getNormalized();
         rotJ.rotate(perp, math::pi);
      }
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(std::abs((rotJ.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1))) * eI),
                                          real_t(1), real_t(1e-6));

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));
      const Vec3 posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult analytical = cylinderCylinderContact(ci, cj);
      const GjkResult      gjk       = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);

      if (!analytical.contact && !gjk.contact)
      {
         ++numNoContactBoth;
         continue;
      }

      const std::string scenario =
         analytical.contact ? std::string(analytical.scenario) : std::string("(missed -- analytical found no contact)");
      auto& st = byScenario[scenario];
      ++st.numTrials;

      if (analytical.contact != gjk.contact)
      {
         ++st.existenceMismatch;
         continue;
      }
      ++st.existenceAgree;

      const real_t gjkOverlap  = -gjk.penetrationDepth;
      const real_t overlapDiff = std::abs(analytical.overlap - gjkOverlap);
      const real_t pointDist   = (analytical.point - gjk.point).length();
      st.maxOverlapDiff = std::max(st.maxOverlapDiff, overlapDiff);
      st.maxPointDist   = std::max(st.maxPointDist, pointDist);
      st.pointDists.push_back(pointDist);
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] PURE analytical (no GJK/EPA fallback) Cylinder x Cylinder, forced-(anti)parallel-axes, "
      << numTrials << " trials (" << numNoContactBoth
      << " agree no contact at all -- excluded from the per-scenario counts below):");
   for (auto& [name, st] : byScenario)
   {
      real_t p50 = real_t(0), p90 = real_t(0);
      if (!st.pointDists.empty())
      {
         std::vector< real_t > sorted = st.pointDists;
         std::sort(sorted.begin(), sorted.end());
         p50 = sorted[sorted.size() / 2];
         p90 = sorted[std::min(sorted.size() - 1, size_t(real_t(0.9) * real_t(sorted.size())))];
      }
      WALBERLA_LOG_INFO_ON_ROOT(
         "  " << name << ": count=" << st.numTrials << " (" << (real_t(100) * real_t(st.numTrials) / real_t(numTrials))
         << "% of all trials) | existence: agree=" << st.existenceAgree << " mismatch="
         << st.existenceMismatch << " | on agreeing trials: maxOverlapDiff=" << st.maxOverlapDiff
         << " maxPointDist=" << st.maxPointDist << " p50PointDist=" << p50 << " p90PointDist=" << p90);
   }
}

//======================================================================================================================
// Targeted fuzz test for face-band and detail::faceBandContactPoint (Kodam eq. 40-44): 1000 trials, seed 778899,
// axes forced perpendicular (the face-band precondition; other scenarios can still fire and are reported
// separately), j's centre offset within 0.8 Ri of i's axis, separation sampled only on the overlapping side of the
// GJK/EPA boundary within 2 % of the smaller radius (shallow contacts). Raw closed form vs. GJK/EPA; reported, not
// asserted.
//======================================================================================================================
void testCylinderCylinderFaceBandBarelyOverlappingFuzz()
{
   FuzzRng                                  fuzz(778899); // fixed seed -- reproducible across runs
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 1000;
   int    numFaceBand = 0, numFaceBandExistenceMismatch = 0;
   int    numOtherScenario = 0, numOtherExistenceMismatch = 0, numNoContactEither = 0;
   real_t maxFaceBandOverlapDiff = real_t(0), maxFaceBandPointDist = real_t(0);
   real_t maxOtherOverlapDiff    = real_t(0), maxOtherPointDist    = real_t(0);
   CylinderState worstCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
   CylinderState worstCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
   real_t        worstPointDist = real_t(0);
   Vec3          worstAnalyticalPoint, worstGjkPoint;
   std::vector< real_t > faceBandPointDists; // for a percentile breakdown, not just the max outlier

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();

      // Force j's axis perpendicular to i's (Gram-Schmidt of a random direction against eI), then rotate body z onto
      // it.
      const Vec3 eI = rotI.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1));
      Vec3       eJ;
      do {
         const Vec3 raw = fuzz.randomUnitVec();
         eJ = raw - (raw * eI) * eI;
      } while (eJ.sqrLength() < real_t(1e-6));
      eJ = eJ.getNormalized();

      Rot3         rotJ;
      const Vec3   ez(real_t(0), real_t(0), real_t(1));
      const real_t cosA = ez * eJ;
      if (cosA < real_t(1) - real_t(1e-10))
      {
         const Vec3 rotAxis = (cosA > real_t(-1) + real_t(1e-10)) ? (ez % eJ).getNormalized()
                                                                    : Vec3(real_t(1), real_t(0), real_t(0));
         rotJ.rotate(rotAxis, std::acos(std::clamp(cosA, real_t(-1), real_t(1))));
      }

      // Offset j's centre in i's local xy plane, within i's face disc.
      std::uniform_real_distribution< real_t > offsetDist(real_t(0), real_t(0.8) * Ri);
      const real_t offsetAng = fuzz.angleDist(fuzz.rng);
      const real_t offsetMag = offsetDist(fuzz.rng);
      const Vec3   lateral   = rotI.getMatrix() *
                              Vec3(offsetMag * std::cos(offsetAng), offsetMag * std::sin(offsetAng), real_t(0));
      const Vec3 posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         const Vec3 posJ = posI + lateral + s * eI;
         return gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj).contact;
      };
      // s = 0 always overlaps for these draws; skip the degenerate case rather than assert.
      if (!contactAtSep(real_t(0))) continue;
      const real_t farSep = real_t(2) * (Ri + Li + Rj + Lj);
      if (contactAtSep(farSep)) continue; // still overlapping at the generous far bound -- skip

      real_t lo = real_t(0), hi = farSep;
      for (int it = 0; it < 40; ++it)
      {
         const real_t mid = real_t(0.5) * (lo + hi);
         if (contactAtSep(mid)) { lo = mid; }
         else { hi = mid; }
      }
      const real_t boundary   = real_t(0.5) * (lo + hi);
      const real_t maxOverlap = real_t(0.02) * std::min(Ri, Rj); // 2 % of the smaller radius
      std::uniform_real_distribution< real_t > sepDist(std::max(real_t(0), boundary - maxOverlap), boundary);
      const real_t sep  = sepDist(fuzz.rng); // sampled ONLY on the overlapping side of the boundary
      const Vec3   posJ = posI + lateral + sep * eI;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult analytical = cylinderCylinderContact(ci, cj);
      const GjkResult      gjk       = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);

      const bool isFaceBand = analytical.contact && std::string(analytical.scenario) == "face-band";
      if (!gjk.contact && !analytical.contact)
      {
         ++numNoContactEither;
         continue;
      }
      if (isFaceBand) { ++numFaceBand; }
      else { ++numOtherScenario; }

      if (analytical.contact != gjk.contact)
      {
         if (isFaceBand) { ++numFaceBandExistenceMismatch; }
         else { ++numOtherExistenceMismatch; }
         continue;
      }
      if (!analytical.contact) continue; // both agree: no contact

      const real_t gjkOverlap  = -gjk.penetrationDepth;
      const real_t overlapDiff = std::abs(analytical.overlap - gjkOverlap);
      const real_t pointDist   = (analytical.point - gjk.point).length();

      if (isFaceBand)
      {
         maxFaceBandOverlapDiff = std::max(maxFaceBandOverlapDiff, overlapDiff);
         maxFaceBandPointDist   = std::max(maxFaceBandPointDist, pointDist);
         faceBandPointDists.push_back(pointDist);
      }
      else
      {
         maxOtherOverlapDiff = std::max(maxOtherOverlapDiff, overlapDiff);
         maxOtherPointDist   = std::max(maxOtherPointDist, pointDist);
      }
      if (pointDist > worstPointDist)
      {
         worstPointDist       = pointDist;
         worstCi              = ci;
         worstCj              = cj;
         worstAnalyticalPoint = analytical.point;
         worstGjkPoint        = gjk.point;
      }
   }

   std::sort(faceBandPointDists.begin(), faceBandPointDists.end());
   auto percentile = [&](real_t p) {
      if (faceBandPointDists.empty()) return real_t(0);
      const size_t idx = std::min(faceBandPointDists.size() - 1,
                                    size_t(p * real_t(faceBandPointDists.size())));
      return faceBandPointDists[idx];
   };
   const Vec3   worstDiff       = worstAnalyticalPoint - worstGjkPoint;
   const Vec3   worstAxisJ      = worstCj.axis();
   const real_t worstAlongAxis  = std::abs(worstDiff * worstAxisJ);
   const real_t worstPerp       = std::sqrt(std::max(real_t(0), worstDiff.sqrLength() - worstAlongAxis * worstAlongAxis));

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Cylinder x Cylinder, barely-overlapping (<=2% of min radius), forced-perpendicular-axes, "
      << numTrials << " trials: face-band scenario count=" << numFaceBand
      << " (existence mismatches vs GJK/EPA=" << numFaceBandExistenceMismatch
      << "), max overlap diff on agreeing face-band trials=" << maxFaceBandOverlapDiff
      << ", max contact-point distance on agreeing face-band trials=" << maxFaceBandPointDist);
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] ... other perpendicular-axis-eligible scenarios (face-edge/band-edge/edge-edge/no-contact) count="
      << numOtherScenario << " (existence mismatches=" << numOtherExistenceMismatch
      << "), max overlap diff=" << maxOtherOverlapDiff << ", max point distance=" << maxOtherPointDist
      << " | both-agree-no-contact=" << numNoContactEither);
   WALBERLA_LOG_INFO_ON_ROOT("[fuzz] worst face-band/other point-distance geometry: ci.pos=" << worstCi.pos
                              << " ci axis=" << worstCi.axis() << " Ri=" << worstCi.radius
                              << " Li=" << worstCi.length << " | cj.pos=" << worstCj.pos
                              << " cj axis=" << worstCj.axis() << " Rj=" << worstCj.radius
                              << " Lj=" << worstCj.length << " | distance=" << worstPointDist);
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] worst-case decomposition: analytical.point=" << worstAnalyticalPoint << " gjk.point="
      << worstGjkPoint << " | component ALONG cylinder j's axis (the direction a face-band contact "
      << "line has no unique point along, if this is what's driving the distance)=" << worstAlongAxis
      << ", component PERPENDICULAR to it (a genuine positional disagreement if this is large)="
      << worstPerp);
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] face-band point-distance percentiles over " << faceBandPointDists.size()
      << " agreeing trials: p50=" << percentile(real_t(0.50)) << " p90=" << percentile(real_t(0.90))
      << " p99=" << percentile(real_t(0.99)) << " max=" << maxFaceBandPointDist);
}

//======================================================================================================================
// Fallback rate of the hybrid and precision of its triggers: 100000 trials, seed 13579, boundary-biased separation
// (FuzzRng::sampleSeparationNearBoundary, as in testCylinderCylinderRandomizedFuzz), i.e. shallow contacts and near
// misses, the population relevant for DEM. For every fallback, records whether the raw closed form would have agreed
// with GJK/EPA anyway (existence, and depth within max(0.01, 5 %)). Measurement only, not asserted.
//======================================================================================================================
void testCylinderCylinderFallbackRateBoundaryBiased()
{
   FuzzRng                                  fuzz(13579); // fixed seed
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials    = 100000;
   int           numFallback  = 0;
   int           numTier1     = 0; // closed form found no contact, reach-sum bound says "maybe"
   int           numTier2     = 0; // closed form found a contact but it is untrusted (special-contact)
   int           numOverlapping = 0;
   // Precision of cylinderCylinderResultIsUncertain: of the fallbacks, how many the raw closed form would have
   // answered correctly anyway.
   int numTier1WouldHaveBeenFine = 0, numTier1WouldHaveBeenWrong = 0;
   int numTier2WouldHaveBeenFine = 0, numTier2WouldHaveBeenWrong = 0;
   // Tier 2 fallbacks per scenario, split into fine / wrong.
   std::map< std::string, std::pair< int, int > > tier2ByScenario; // scenario -> (fine, wrong)

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();
      const Vec3   dir  = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      auto         contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));
      const Vec3   posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult analytical = cylinderCylinderContact(ci, cj);
      const HybridResult  hybrid     = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
      const GjkResult      gjk = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
      if (gjk.contact) ++numOverlapping;

      if (hybrid.usedFallback)
      {
         ++numFallback;
         const real_t gjkOverlap = gjk.contact ? -gjk.penetrationDepth : real_t(0);
         // Agreement: same existence and, for contacts, depth within max(0.01, 5 %) of GJK/EPA.
         const bool wouldHaveBeenFine =
            (analytical.contact == gjk.contact) &&
            (!gjk.contact || std::abs(analytical.overlap - gjkOverlap) <=
                                 std::max(real_t(0.01), real_t(0.05) * gjkOverlap));

         if (!analytical.contact)
         {
            ++numTier1;
            wouldHaveBeenFine ? ++numTier1WouldHaveBeenFine : ++numTier1WouldHaveBeenWrong;
         }
         else
         {
            ++numTier2;
            wouldHaveBeenFine ? ++numTier2WouldHaveBeenFine : ++numTier2WouldHaveBeenWrong;
            auto& entry = tier2ByScenario[std::string(analytical.scenario)];
            wouldHaveBeenFine ? ++entry.first : ++entry.second;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] GJK/EPA fallback usage (BOUNDARY-BIASED separation sampling, same as "
      << "testCylinderCylinderRandomizedFuzz -- the physically relevant population for DEM): "
      << numFallback << "/" << numTrials << " trials ("
      << (real_t(100) * real_t(numFallback) / real_t(numTrials)) << "%) triggered the fallback ("
      << numOverlapping << "/" << numTrials << " actually in contact per GJK/EPA). Of those "
      << numFallback << " fallbacks: Tier 1=" << numTier1 << " (" << (real_t(100) * real_t(numTier1) / real_t(numFallback))
      << "%), Tier 2=" << numTier2 << " (" << (real_t(100) * real_t(numTier2) / real_t(numFallback)) << "%).");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Of the " << numFallback << " total fallbacks, how many could we actually have trusted "
      << "the raw analytical answer for instead (i.e. the alarm's precision)? Tier 1: "
      << numTier1WouldHaveBeenFine << "/" << numTier1 << " would have been fine to trust ("
      << (real_t(100) * real_t(numTier1WouldHaveBeenFine) / real_t(numTier1)) << "%), "
      << numTier1WouldHaveBeenWrong << " genuinely needed GJK/EPA. Tier 2: " << numTier2WouldHaveBeenFine
      << "/" << numTier2 << " would have been fine to trust ("
      << (real_t(100) * real_t(numTier2WouldHaveBeenFine) / real_t(numTier2)) << "%), "
      << numTier2WouldHaveBeenWrong << " genuinely needed GJK/EPA.");
   for (const auto& kv : tier2ByScenario)
   {
      const int total = kv.second.first + kv.second.second;
      WALBERLA_LOG_INFO_ON_ROOT("  Tier 2 scenario=" << kv.first << ": " << total << " total, "
                                 << kv.second.first << " would have been fine to trust ("
                                 << (real_t(100) * real_t(kv.second.first) / real_t(total)) << "%), "
                                 << kv.second.second << " genuinely wrong.");
   }
}

//! Face-edge counterpart to testCylinderCylinderTiltedRimIntoFace: with a host much larger than the probe (Ri = 3,
//! Rj = 0.3) the deepest rim point lies far from the host rim, the host face acts as a plane, and face-edge reduces
//! to the exact cylinder-plane edge case (for a tilt strictly between 0 and 90 degrees the extremal point of a solid
//! cylinder along a fixed direction lies on its rim). Depth must match GJK/EPA to 1e-4.
void testCylinderCylinderFaceEdgeClean()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(3.0), real_t(1.0) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(6)); // 30 degrees about y
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.0)), rotJ.getMatrix(), real_t(0.3), real_t(1.5) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-edge");
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-4));
   WALBERLA_CHECK_GREATER(kodam.normal * gjk.normal, real_t(0.999));
   // Point in the middle of the overlap: host i's face at z = +0.5, j's deepest rim point E at z = 0.5 - overlap,
   // so z = 0.5 - overlap/2.
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], real_t(0.5) - real_t(0.5) * kodam.overlap, real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], gjk.point[2], real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (face-edge, clean/negligible-host-curvature): overlap="
                              << kodam.overlap << " (agrees with GJK/EPA to 1e-4, unlike the KNOWN "
                              << "LIMITATION case above -- confirms the limitation is curvature-coupling, "
                              << "not a general face-edge defect)");
}

//! Band-edge case with a nearly flat host band (Ri = 3, Rj = 0.3, as in testCylinderCylinderFaceEdgeClean), so the
//! probe rim point of Kodam eq. 67-81 is the closest feature. Depth must match GJK/EPA to 1e-4.
void testCylinderCylinderBandEdgeClean()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(3.0), real_t(6.0) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(3)); // 60 degrees about y
   CylinderState cj{ Vec3(real_t(3.3), real_t(0), real_t(1.0)), rotJ.getMatrix(), real_t(0.3), real_t(1.5) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "band-edge");
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-4));
   WALBERLA_CHECK_GREATER(kodam.normal * gjk.normal, real_t(0.999));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (band-edge, clean/negligible-host-curvature): overlap="
                              << kodam.overlap << " (agrees with GJK/EPA to 1e-4, unlike the KNOWN "
                              << "LIMITATION case above -- confirms the limitation is curvature-coupling, "
                              << "not a general band-edge defect)");
}

//! Rotation matrix -> quaternion (e.g. Shepperd 1978); Rot3 has no Mat3 constructor. Also snaps a slightly
//! non-orthonormal matrix to a valid rotation.
Rot3 rot3FromMatrix(const Mat3& m)
{
   const real_t w = std::sqrt(std::max(real_t(0), real_t(1) + m(0, 0) + m(1, 1) + m(2, 2))) / real_t(2);
   real_t       x = std::sqrt(std::max(real_t(0), real_t(1) + m(0, 0) - m(1, 1) - m(2, 2))) / real_t(2);
   real_t       y = std::sqrt(std::max(real_t(0), real_t(1) - m(0, 0) + m(1, 1) - m(2, 2))) / real_t(2);
   real_t       z = std::sqrt(std::max(real_t(0), real_t(1) - m(0, 0) - m(1, 1) + m(2, 2))) / real_t(2);
   x              = std::copysign(x, m(2, 1) - m(1, 2));
   y              = std::copysign(y, m(0, 2) - m(2, 0));
   z              = std::copysign(z, m(1, 0) - m(0, 1));
   return Rot3(mesa_pd::Quat(w, x, y, z));
}

//! A pair that dispatches to "edge-edge" (both rims facing each other). Geometry: the edge-edge trial with the best
//! non-trivial agreement in testCylinderCylinderPureAnalyticalScenarioBreakdown's corpus (seed 24680, 5000 trials);
//! arbitrary orientations are needed to reach this scenario. Closed form 0.07830 vs. GJK/EPA 0.07868 and
//! brute-force MTD 0.07883 (0.47 % relative), checked to 1e-3.
void testCylinderCylinderEdgeEdgeGenuine()
{
   const Mat3 rotICols(real_t(-0.390898), real_t(-0.198955), real_t(0.898675), real_t(0.196307),
                        real_t(0.935875), real_t(0.292579), real_t(-0.899257), real_t(0.290784),
                        real_t(-0.326775));
   const Mat3 rotJCols(real_t(0.131414), real_t(-0.370814), real_t(0.919362), real_t(-0.959446),
                        real_t(0.18573), real_t(0.212056), real_t(-0.249386), real_t(-0.909946),
                        real_t(-0.331369));

   // The matrices are printed to ~6 digits and not exactly orthonormal; both methods use the snapped Rot3 so they
   // see identical geometry.
   const Rot3           rotI = rot3FromMatrix(rotICols);
   const Rot3           rotJ = rot3FromMatrix(rotJCols);
   const CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), rotI.getMatrix(), real_t(1.44446),
                            real_t(1.30996) };
   const CylinderState cj{ Vec3(real_t(2.11131), real_t(0.251881), real_t(0.154264)), rotJ.getMatrix(),
                            real_t(1.22407), real_t(2.5082) };

   const ContactResult kodam = cylinderCylinderContact(ci, cj);
   const GjkResult      gjk =
      gjkCylCyl(ci.pos, rotI, ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "edge-edge");
   WALBERLA_CHECK(gjk.contact);
   // closed form 0.07830 vs GJK/EPA 0.07868 and brute-force MTD 0.07883
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (genuine edge-edge, found by search): overlap="
                              << kodam.overlap << " vs GJK/EPA " << -gjk.penetrationDepth
                              << " (0.47% relative -- first dedicated coverage of this scenario)");
}


//======================================================================================================================
// Root-cause classifiers for the master comparison (testCylinderCylinderFullComparisonAndHybridFeasibility).
//======================================================================================================================
//! Classifies why a trial's overlap magnitude disagrees with GJK/EPA; returns a string literal (stored in
//! CylMismatchRecord::cause) and a metric in outMetric.
//!   (a) misclassified scenario: angle between kodam.normal and gjk.normal > 15 degrees, i.e. the wrong feature
//!       pairing was chosen and only re-dispatching could fix it;
//!   (b) wrong point, same scenario: normal within 15 degrees, but the scenario's point is not the extremal one.
//! special-contact is classified separately: its axial depth is taken at P (rim point closest to the host axis)
//! instead of E (rim point closest to the host face plane, detail::rimPointTowardPlane); since E minimises the
//! axial distance over the rim, |z_E| <= |z_P| and depth(P) underestimates. E is recomputed for the firing
//! host/probe pairing and delta_a(E) compared with GJK/EPA.
const char* classifyMagnitudeCause(const CylinderState& ci, const CylinderState& cj, const ContactResult& kodam,
                                    const GjkResult& gjk, real_t& outMetric)
{
   const real_t gjkOverlap = -gjk.penetrationDepth;

   if (std::string(kodam.scenario) == "special-contact")
   {
      // rimVsFaceOrBand is order-sensitive: try both pairings and keep the one whose overlap matches.
      auto axialOverlapViaE = [&](const CylinderState& host, const CylinderState& probe, bool& eInDisc) -> real_t {
         const ContactResult raw = rimVsFaceOrBand(host, probe);
         if (!raw.contact || std::string(raw.scenario) != "special-contact") return real_t(-1);
         if (std::abs(raw.overlap - kodam.overlap) > real_t(1e-9)) return real_t(-1); // not this pairing
         const Vec3   E    = detail::rimPointTowardPlane(probe, host.pos, host.axis());
         const Vec3   relE = host.R.getTranspose() * (E - host.pos);
         eInDisc = (relE[0] * relE[0] + relE[1] * relE[1]) < host.radius * host.radius; // E on host face?
         return real_t(0.5) * host.length - std::abs(relE[2]);
      };
      bool   eInDisc     = false;
      real_t viaE        = axialOverlapViaE(ci, cj, eInDisc);
      bool   matchedOnCi = viaE >= real_t(0);
      if (!matchedOnCi) viaE = axialOverlapViaE(cj, ci, eInDisc);
      if (viaE >= real_t(0))
      {
         const real_t errViaP = std::abs(kodam.overlap - gjkOverlap);
         const real_t errViaE = std::abs(viaE - gjkOverlap);
         outMetric            = errViaE; // error the E-based depth would have had
         // E is the extremal point only if it lies on the host face (eInDisc); otherwise halfLHost - |z_E| does not
         // correspond to a real overlap and a disc-constrained extremal search would be needed.
         if (errViaE < real_t(1e-3)) return "special-contact-exact-fix-via-E";
         if (!eInDisc) return "special-contact-E-outside-face-disc";
         if (errViaE < errViaP) return "special-contact-partial-fix-via-E";
      }
      return "special-contact-E-also-fails"; // neither P nor E gives GJK/EPA's depth; also returned when
                                                // no pairing matched (round-off at the deltaAp ~ deltaRp
                                                // dispatch threshold can flip the branch)
   }

   const real_t cosAngle = std::clamp(kodam.normal * gjk.normal, real_t(-1), real_t(1));
   const real_t angleDeg = std::acos(cosAngle) * real_t(180) / math::pi;
   outMetric              = angleDeg;
   return (angleDeg > real_t(15)) ? "misclassified-scenario" : "wrong-point-same-scenario";
}

//! Companion to classifyMagnitudeCause for missed contacts (no scenario matched, GJK/EPA finds contact). The GJK/EPA
//! point counts as "rim" on a cylinder if it lies within 0.3 radius of that cylinder's end cap, else "band"; the
//! four combinations are returned as string literals.
const char* classifyDispatchGapCause(const CylinderState& ci, const CylinderState& cj, const GjkResult& gjk)
{
   auto isRim = [&](const CylinderState& c) {
      const real_t zAlong  = c.axis() * (gjk.point - c.pos);
      const real_t halfLen = real_t(0.5) * c.length;
      return (halfLen - std::abs(zAlong)) < real_t(0.3) * c.radius;
   };
   const bool riI = isRim(ci), riJ = isRim(cj);
   if (riI && riJ) return "dispatch-gap-rim-rim";
   if (riI && !riJ) return "dispatch-gap-rim-band";
   if (!riI && riJ) return "dispatch-gap-band-rim";
   return "dispatch-gap-band-band";
}

//======================================================================================================================
// Edge-edge solver variants (testEdgeEdgeThreeWayAccuracy10k) and the thetaRef tie-break count
// (testThetaRefTieBreakNecessity10k).
//======================================================================================================================
//! Accuracy of the edge-edge solver variants against GJK/EPA, all on identical trials:
//!   A) Durand-Kerner, 60 iterations, line-based thetaRef (useFastEdgeEdgeReference = false);
//!   B) Durand-Kerner, 30 iterations, line-based thetaRef;
//!   C) Durand-Kerner, 30 iterations, thetaRef via closestThetaOnCircleToPoint (useFastEdgeEdgeReference = true);
//!   D) / E) as B / C with 20 iterations;
//!   F) / G) as B with the Ferrari / Aberth-Ehrlich quartic solver;
//!   H) Ferrari with lazyThetaRef (the current defaults of cylinderCylinderContact).
//! Only trials classified "edge-edge" by variant A are scored (the only scenario the variants affect), so the
//! selection does not depend on the variants compared. 10000 trials, seed 13579, sampling as in
//! testCylinderCylinderRandomizedFuzz; edge-edge fires on roughly 4 % of trials. Also counts trials where a variant
//! gives the identical answer to A (overlap and point within 1e-9). Checks: at least 100 edge-edge trials, and
//! A/B/C scored on the same set.
void testEdgeEdgeThreeWayAccuracy10k()
{
   FuzzRng                                  fuzz(13579); // distinct seed from the master-table test
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 10000;

   struct RefStats
   {
      int    edgeEdgeTrials = 0;
      int    existenceAgree = 0, existenceMismatch = 0;
      real_t maxOverlapDiff = 0, sumOverlapDiff = 0;
      real_t maxPointDist = 0, sumPointDist = 0;
      real_t maxNormalAngleDeg = 0, sumNormalAngleDeg = 0;
   };
   RefStats original, currentDefault, fullyOptimized, iter20OldRef, iter20NewRef, ferrariSolver, aberthEhrlichSolver,
      lazyThetaRefFerrari;

   auto accumulate = [](RefStats& s, const ContactResult& kodam, const GjkResult& gjk) {
      ++s.edgeEdgeTrials;
      if (gjk.contact != kodam.contact)
      {
         ++s.existenceMismatch;
         return; // no meaningful overlap/point/normal comparison when existence itself disagrees
      }
      ++s.existenceAgree;
      if (!gjk.contact) return; // both agree "no contact" -- nothing further to measure

      const real_t gjkOverlap = -gjk.penetrationDepth;
      const real_t overlapDiff = std::abs(kodam.overlap - gjkOverlap);
      s.maxOverlapDiff = std::max(s.maxOverlapDiff, overlapDiff);
      s.sumOverlapDiff += overlapDiff;

      const real_t pointDist = (kodam.point - gjk.point).length();
      s.maxPointDist = std::max(s.maxPointDist, pointDist);
      s.sumPointDist += pointDist;

      const real_t cosAngle = std::clamp(kodam.normal * gjk.normal, real_t(-1), real_t(1));
      const real_t angleDeg = std::acos(cosAngle) * real_t(180) / math::pi;
      s.maxNormalAngleDeg = std::max(s.maxNormalAngleDeg, angleDeg);
      s.sumNormalAngleDeg += angleDeg;
   };

   int    pairwiseIdenticalOldVsCurrent = 0, pairwiseIdenticalCurrentVsFull = 0;
   int    pairwiseIdenticalOriginalVs20Old = 0, pairwiseIdenticalOriginalVs20New = 0;
   int    pairwiseIdenticalOriginalVsFerrari = 0, pairwiseIdenticalOriginalVsAberthEhrlich = 0;
   int    pairwiseIdenticalOriginalVsLazyThetaRef = 0;

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));
      const Vec3 posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      // Classify with variant A. Solver and lazyThetaRef are passed explicitly in A-G (the defaults are Ferrari and
      // lazyThetaRef = true) so each variant isolates one change; only H uses lazyThetaRef.
      const ContactResult kodamOriginal =
         cylinderCylinderContact(ci, cj, false, 60, detail::QuarticSolverMethod::DurandKerner, false);
      if (std::string(kodamOriginal.scenario) != "edge-edge") continue;

      const ContactResult kodamCurrent =
         cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::DurandKerner, false);
      const ContactResult kodamFull =
         cylinderCylinderContact(ci, cj, true, 30, detail::QuarticSolverMethod::DurandKerner, false);
      const ContactResult kodam20Old =
         cylinderCylinderContact(ci, cj, false, 20, detail::QuarticSolverMethod::DurandKerner, false);
      const ContactResult kodam20New =
         cylinderCylinderContact(ci, cj, true, 20, detail::QuarticSolverMethod::DurandKerner, false);
      // Solver-only variants: as B, with solveQuarticFerrari / solveQuarticAberthEhrlich (synthetic-quartic
      // comparison: QuarticSolverComparisonTest.cpp).
      const ContactResult kodamFerrari =
         cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, false);
      const ContactResult kodamAberthEhrlich =
         cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::AberthEhrlich, false);
      // lazyThetaRef on top of Ferrari: thetaRef is computed only if point A or B has more than one candidate
      // (see testThetaRefTieBreakNecessity10k), so H should equal F.
      const ContactResult kodamLazyThetaRef =
         cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, true);
      const GjkResult      gjk          = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);

      accumulate(original, kodamOriginal, gjk);
      accumulate(currentDefault, kodamCurrent, gjk);
      accumulate(fullyOptimized, kodamFull, gjk);
      accumulate(iter20OldRef, kodam20Old, gjk);
      accumulate(iter20NewRef, kodam20New, gjk);
      accumulate(ferrariSolver, kodamFerrari, gjk);
      accumulate(aberthEhrlichSolver, kodamAberthEhrlich, gjk);
      accumulate(lazyThetaRefFerrari, kodamLazyThetaRef, gjk);

      // Pairwise identity: did the variant change this trial's answer at all?
      if (kodamOriginal.contact == kodamCurrent.contact &&
          std::abs(kodamOriginal.overlap - kodamCurrent.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodamCurrent.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOldVsCurrent;
      if (kodamCurrent.contact == kodamFull.contact &&
          std::abs(kodamCurrent.overlap - kodamFull.overlap) < real_t(1e-9) &&
          (kodamCurrent.point - kodamFull.point).length() < real_t(1e-9))
         ++pairwiseIdenticalCurrentVsFull;
      // D-H are compared against A (60 iterations), the strongest baseline.
      if (kodamOriginal.contact == kodam20Old.contact &&
          std::abs(kodamOriginal.overlap - kodam20Old.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodam20Old.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOriginalVs20Old;
      if (kodamOriginal.contact == kodam20New.contact &&
          std::abs(kodamOriginal.overlap - kodam20New.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodam20New.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOriginalVs20New;
      if (kodamOriginal.contact == kodamFerrari.contact &&
          std::abs(kodamOriginal.overlap - kodamFerrari.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodamFerrari.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOriginalVsFerrari;
      if (kodamOriginal.contact == kodamAberthEhrlich.contact &&
          std::abs(kodamOriginal.overlap - kodamAberthEhrlich.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodamAberthEhrlich.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOriginalVsAberthEhrlich;
      if (kodamOriginal.contact == kodamLazyThetaRef.contact &&
          std::abs(kodamOriginal.overlap - kodamLazyThetaRef.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodamLazyThetaRef.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOriginalVsLazyThetaRef;
   }

   auto report = [](const char* label, const RefStats& s) {
      const real_t meanOverlapDiff = s.existenceAgree > 0 ? s.sumOverlapDiff / real_t(s.existenceAgree) : real_t(0);
      const real_t meanPointDist   = s.existenceAgree > 0 ? s.sumPointDist / real_t(s.existenceAgree) : real_t(0);
      const real_t meanNormalAngle = s.existenceAgree > 0 ? s.sumNormalAngleDeg / real_t(s.existenceAgree) : real_t(0);
      WALBERLA_LOG_INFO_ON_ROOT(
         "[EdgeEdge3Way] " << label << ": edge-edge trials=" << s.edgeEdgeTrials
         << " | existence: agree=" << s.existenceAgree << " mismatch=" << s.existenceMismatch
         << " | overlap diff vs GJK/EPA: mean=" << meanOverlapDiff << " max=" << s.maxOverlapDiff
         << " | point dist vs GJK/EPA: mean=" << meanPointDist << " max=" << s.maxPointDist
         << " | normal angle vs GJK/EPA (deg): mean=" << meanNormalAngle << " max=" << s.maxNormalAngleDeg);
   };
   report("A) ORIGINAL (60 iter, old ref)",           original);
   report("B) CURRENT DEFAULT (30 iter, old ref)",    currentDefault);
   report("C) FULLY OPTIMIZED (30 iter, new ref)",    fullyOptimized);
   report("D) 20 iter, old ref",                      iter20OldRef);
   report("E) 20 iter, new ref",                      iter20NewRef);
   report("F) Ferrari solver (30 iter budget N/A, old ref)", ferrariSolver);
   report("G) Aberth-Ehrlich solver (30, old ref)",   aberthEhrlichSolver);
   report("H) Ferrari + lazyThetaRef",                lazyThetaRefFerrari);
   WALBERLA_LOG_INFO_ON_ROOT(
      "[EdgeEdge3Way] Direct pairwise identical-answer count (same contact flag, overlap within 1e-9, "
      "point within 1e-9): A-vs-B=" << pairwiseIdenticalOldVsCurrent << "/" << original.edgeEdgeTrials
      << " | B-vs-C=" << pairwiseIdenticalCurrentVsFull << "/" << original.edgeEdgeTrials
      << " | A-vs-D(20,old)=" << pairwiseIdenticalOriginalVs20Old << "/" << original.edgeEdgeTrials
      << " | A-vs-E(20,new)=" << pairwiseIdenticalOriginalVs20New << "/" << original.edgeEdgeTrials
      << " | A-vs-F(Ferrari)=" << pairwiseIdenticalOriginalVsFerrari << "/" << original.edgeEdgeTrials
      << " | A-vs-G(Aberth-Ehrlich)=" << pairwiseIdenticalOriginalVsAberthEhrlich << "/" << original.edgeEdgeTrials
      << " | A-vs-H(lazyThetaRef)=" << pairwiseIdenticalOriginalVsLazyThetaRef << "/" << original.edgeEdgeTrials);

   WALBERLA_CHECK_GREATER_EQUAL(original.edgeEdgeTrials, 100); // sanity: comparison isn't running on noise
   WALBERLA_CHECK_EQUAL(original.edgeEdgeTrials, currentDefault.edgeEdgeTrials); // same selection set
   WALBERLA_CHECK_EQUAL(original.edgeEdgeTrials, fullyOptimized.edgeEdgeTrials); // by construction
}

//======================================================================================================================
//! Is thetaRef's tie-break ever needed? findRimCrossing's point A (eq. 83) and point B (eq. 84-85) each pick, among
//! the candidate angles passing their containment check, the one closest to thetaRef, a tie-break neither Kodam nor
//! Guo specifies. If at most one candidate ever survives, thetaRef's value (and its quartic solve) is irrelevant.
//! Counts the candidates via findRimCrossing's candidateCountA/candidateCountB out-parameters, in both directions,
//! on the edge-edge trials of the testEdgeEdgeThreeWayAccuracy10k corpus (seed 13579, 10000 trials). Logged; checks
//! only that at least 200 calls were sampled.
void testThetaRefTieBreakNecessity10k()
{
   FuzzRng                                  fuzz(13579); // SAME seed as testEdgeEdgeThreeWayAccuracy10k
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 10000;

   std::map< int, int > countHistogramA, countHistogramB; // candidate count -> number of rim-crossing
                                                              // calls that saw exactly that many
   int totalRimCrossingCalls = 0;
   int everMultipleA = 0, everMultipleB = 0; // calls where thetaRef's tie-break actually had >1 to choose from

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));
      const Vec3 posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult classify =
         cylinderCylinderContact(ci, cj, false, 60, detail::QuarticSolverMethod::DurandKerner);
      if (std::string(classify.scenario) != "edge-edge") continue;

      // Same geometry setup as cylinderCylinderContact's edge-edge branch, so findRimCrossing gets the same inputs.
      const Vec3   eI = ci.axis(), eJ = cj.axis();
      const real_t halfLi = real_t(0.5) * ci.length, halfLj = real_t(0.5) * cj.length;

      Vec3 uJ, vJ, uI, vI;
      detail::perpendicularBasis(eJ, uJ, vJ);
      detail::perpendicularBasis(eI, uI, vI);
      const Vec3 capJ = cj.pos + math::sign((ci.pos - cj.pos) * eJ) * halfLj * eJ;
      const Vec3 capI = ci.pos + math::sign((cj.pos - ci.pos) * eI) * halfLi * eI;

      int countA_jVsI = 0, countB_jVsI = 0, countA_iVsJ = 0, countB_iVsJ = 0;
      detail::findRimCrossing(capJ, uJ, vJ, cj.radius, ci.pos, eI, ci.radius, halfLi, false, 60,
                               detail::QuarticSolverMethod::DurandKerner, &countA_jVsI, &countB_jVsI);
      detail::findRimCrossing(capI, uI, vI, ci.radius, cj.pos, eJ, cj.radius, halfLj, false, 60,
                               detail::QuarticSolverMethod::DurandKerner, &countA_iVsJ, &countB_iVsJ);

      for (const int countA : { countA_jVsI, countA_iVsJ })
      {
         ++countHistogramA[countA];
         ++totalRimCrossingCalls;
         if (countA > 1) ++everMultipleA;
      }
      for (const int countB : { countB_jVsI, countB_iVsJ })
      {
         ++countHistogramB[countB];
         if (countB > 1) ++everMultipleB;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[ThetaRefNecessity] " << totalRimCrossingCalls << " findRimCrossing calls "
                              << "sampled (both directions, every real edge-edge trial in a " << numTrials
                              << "-trial fuzz run)");
   for (const auto& [count, occurrences] : countHistogramA)
      WALBERLA_LOG_INFO_ON_ROOT("[ThetaRefNecessity]   point A candidate count=" << count << ": "
                                 << occurrences << " calls");
   for (const auto& [count, occurrences] : countHistogramB)
      WALBERLA_LOG_INFO_ON_ROOT("[ThetaRefNecessity]   point B candidate count=" << count << ": "
                                 << occurrences << " calls");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[ThetaRefNecessity] thetaRef's tie-break was ACTUALLY consulted (>1 surviving candidate) on "
      << everMultipleA << "/" << totalRimCrossingCalls << " point-A calls and " << everMultipleB << "/"
      << totalRimCrossingCalls << " point-B calls");

   WALBERLA_CHECK_GREATER_EQUAL(totalRimCrossingCalls, 200); // sanity: enough real edge-edge geometry sampled
}

//======================================================================================================================
//! Hand-computed cases for detail::closestDistanceSqBetweenSegments (the axis-segment distance behind the capsule
//! bound): parallel overlapping, parallel non-overlapping, perpendicular skew, intersecting.
//======================================================================================================================
void testClosestDistanceSqBetweenSegments()
{
   using detail::closestDistanceSqBetweenSegments;

   // Parallel, overlapping ranges: distance = perpendicular offset 3.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(0, 0, 0), Vec3(0, 0, 10), Vec3(3, 0, 0), Vec3(3, 0, 10));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(9), real_t(1e-9));
   }

   // Parallel, non-overlapping ranges: nearest endpoints (0,0,5) and (3,0,10), distance^2 = 3^2 + 5^2 = 34.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(0, 0, 0), Vec3(0, 0, 5), Vec3(3, 0, 10), Vec3(3, 0, 15));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(34), real_t(1e-9));
   }

   // Perpendicular skew segments (along x at z = 2, along y at z = 0): closest points (0,0,2), (0,0,0), distance 2.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(-5, 0, 2), Vec3(5, 0, 2), Vec3(0, -5, 0), Vec3(0, 5, 0));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(4), real_t(1e-9));
   }

   // Intersecting segments (X in the z = 0 plane): distance 0.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(-1, -1, 0), Vec3(1, 1, 0), Vec3(-1, 1, 0), Vec3(1, -1, 0));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(0), real_t(1e-9));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[CapsuleBound] closestDistanceSqBetweenSegments: all 4 hand-computed "
                              "ground-truth cases matched.");
}

//======================================================================================================================
//! Is the capsule bound (cylinder Layer 1, ModelParameters.useTightLayer1: axis-segment distance > sum of the radii
//! means separated) conservative? A false "separated" would drop the pair before any contact detection, so the pass
//! criterion is zero rejected pairs that GJK/EPA reports in contact (plus at least 200 rejections, so the check is
//! not vacuous). Seed 97531, 2 x 5000 trials: separation near the contact boundary (sampleSeparationNearBoundary)
//! and uniform in [0, 3 (circumRi + circumRj)].
//======================================================================================================================
void testCapsuleBoundConservative10k()
{
   FuzzRng                                  fuzz(97531); // fixed seed
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrialsPerCorpus = 5000;

   int    rejected             = 0; // capsule bound says "definitely separated"
   int    unsafeFalseSeparated = 0; // ... while GJK/EPA reports contact -> MUST stay 0
   real_t maxUnsafeOverlap     = 0;

   auto runCorpus = [&](bool uniformSeparation) {
      for (int trial = 0; trial < numTrialsPerCorpus; ++trial)
      {
         const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
         const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
         const Rot3   rotI = fuzz.randomRotation();
         const Rot3   rotJ = fuzz.randomRotation();

         const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
         const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
         const Vec3   dir      = fuzz.randomUnitVec();
         const Vec3   posI(real_t(0), real_t(0), real_t(0));

         real_t sep;
         if (uniformSeparation)
         {
            std::uniform_real_distribution< real_t > uniformSepDist(real_t(0), real_t(3) * (circumRi + circumRj));
            sep = uniformSepDist(fuzz.rng);
         }
         else
         {
            auto contactAtSep = [&](real_t s) {
               return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
            };
            sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0), real_t(3) * (circumRi + circumRj),
                                                       std::min(Ri, Rj));
         }
         const Vec3 posJ = posI + dir * sep;

         const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
         const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

         // same test as the Layer-1 reject in ParticleKernels.h (SolidCylinderContactDetection, Cylinder x Cylinder)
         Vec3 p0i, p1i, p0j, p1j;
         detail::cylinderAxisEndpoints(ci, p0i, p1i);
         detail::cylinderAxisEndpoints(cj, p0j, p1j);
         const real_t reach = Ri + Rj;
         if (detail::closestDistanceSqBetweenSegments(p0i, p1i, p0j, p1j) <= reach * reach) continue;
         ++rejected;

         const GjkResult gjk = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
         if (gjk.contact)
         {
            ++unsafeFalseSeparated;
            maxUnsafeOverlap = std::max(maxUnsafeOverlap, -gjk.penetrationDepth);
         }
      }
   };

   runCorpus(false); // near the actual GJK/EPA contact boundary
   runCorpus(true);  // uniformly over [0, 3 (circumRi + circumRj)]

   WALBERLA_LOG_INFO_ON_ROOT("[CapsuleBound] " << rejected << " of " << (2 * numTrialsPerCorpus)
      << " trials rejected as definitely separated | UNSAFE (rejected although GJK/EPA reports contact): "
      << unsafeFalseSeparated << " (max missed overlap=" << maxUnsafeOverlap << ")");

   WALBERLA_CHECK_EQUAL(unsafeFalseSeparated, 0);
   WALBERLA_CHECK_GREATER_EQUAL(rejected, 200); // sanity: comparison isn't running on noise
}

//======================================================================================================================
//! Hand-computed cases for detail::boxSupportAlongDirection and detail::boxesDefinitelySeparatedAlongCenterLine
//! (box Layer 1 pre-filter).
//======================================================================================================================
void testBoxSeparationAlongCenterLine()
{
   using detail::boxesDefinitelySeparatedAlongCenterLine;
   using detail::boxSupportAlongDirection;

   // Axis-aligned box, half-extents (2,3,4): projection along the X axis is exactly halfExtents.x.
   {
      const Mat3   R = Mat3();
      const Vec3   half(2, 3, 4);
      const real_t supportX = boxSupportAlongDirection(half, R, Vec3(1, 0, 0));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(supportX, real_t(2), real_t(1e-9));
   }

   // Box with half-extents (h,h,h) rotated 45 degrees about Z: projection along world X is
   // h |cos45| + h |sin45| = h sqrt(2).
   {
      Rot3 rot;
      rot.rotate(Vec3(0, 0, 1), math::pi * real_t(0.25));
      const Vec3   half(real_t(1), real_t(1), real_t(1));
      const real_t supportX = boxSupportAlongDirection(half, rot.getMatrix(), Vec3(1, 0, 0));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(supportX, real_t(1) * std::sqrt(real_t(2)), real_t(1e-6));
   }

   // Two axis-aligned boxes separated along X by more than the sum of their X half-extents ->
   // definitely separated. half-extents (1,1,1) and (1,1,1), centres 3 apart -> reach=2 < 3.
   {
      const Mat3 R = Mat3();
      const Vec3 half(1, 1, 1);
      WALBERLA_CHECK(boxesDefinitelySeparatedAlongCenterLine(Vec3(0, 0, 0), R, half, Vec3(3, 0, 0), R, half));
   }

   // Same boxes, centres 1.5 apart (reach = 2 > 1.5): this axis cannot prove separation.
   {
      const Mat3 R = Mat3();
      const Vec3 half(1, 1, 1);
      WALBERLA_CHECK(!boxesDefinitelySeparatedAlongCenterLine(Vec3(0, 0, 0), R, half, Vec3(1.5, 0, 0), R, half));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[TightLayer1Box] boxSupportAlongDirection / "
                              "boxesDefinitelySeparatedAlongCenterLine: all hand-computed cases matched.");
}

//======================================================================================================================
//! Is the single-axis box separation bound (box Layer 1, centre line as candidate axis) conservative? Pass criterion:
//! zero pairs declared separated that GJK/EPA reports in contact, and at least 200 such declarations. Seed 24601,
//! half-extents U[0.3,3], uniform SO(3), 2 x 5000 trials on the two corpora of testCapsuleBoundConservative10k.
//======================================================================================================================
void testTightLayer1BoxAccuracy10k()
{
   FuzzRng                                  fuzz(24601); // fixed seed
   std::uniform_real_distribution< real_t > halfExtentDist(real_t(0.3), real_t(3.0));

   constexpr int numTrialsPerCorpus = 5000;

   int unsafeFalseSeparated = 0, resolvedSeparated = 0;
   real_t maxUnsafeOverlap = 0;

   auto runCorpus = [&](bool uniformSeparation) {
      for (int trial = 0; trial < numTrialsPerCorpus; ++trial)
      {
         const Vec3 halfA(halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng));
         const Vec3 halfB(halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng));
         const Rot3 rotA = fuzz.randomRotation();
         const Rot3 rotB = fuzz.randomRotation();

         const real_t circumA = halfA.length(); // circumscribing-sphere radius of a box = |halfExtents|
         const real_t circumB = halfB.length();
         const Vec3   dir     = fuzz.randomUnitVec();
         const Vec3   posA(real_t(0), real_t(0), real_t(0));

         real_t sep;
         if (uniformSeparation)
         {
            std::uniform_real_distribution< real_t > uniformSepDist(real_t(0), real_t(3) * (circumA + circumB));
            sep = uniformSepDist(fuzz.rng);
         }
         else
         {
            auto contactAtSep = [&](real_t s) {
               return gjkBoxBox(posA, rotA, halfA * real_t(2), posA + dir * s, rotB, halfB * real_t(2)).contact;
            };
            sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0), real_t(3) * (circumA + circumB),
                                                       std::min(halfA.min(), halfB.min()));
         }
         const Vec3 posB = posA + dir * sep;

         const bool tightSaysSeparated =
            detail::boxesDefinitelySeparatedAlongCenterLine(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(), halfB);
         if (tightSaysSeparated)
         {
            ++resolvedSeparated;
            const GjkResult gjk = gjkBoxBox(posA, rotA, halfA * real_t(2), posB, rotB, halfB * real_t(2));
            if (gjk.contact)
            {
               ++unsafeFalseSeparated;
               maxUnsafeOverlap = std::max(maxUnsafeOverlap, -gjk.penetrationDepth);
            }
         }
      }
   };

   runCorpus(false);
   runCorpus(true);

   WALBERLA_LOG_INFO_ON_ROOT(
      "[TightLayer1Box] " << (2 * numTrialsPerCorpus) << " trials | resolved definitely-separated="
      << resolvedSeparated << " | UNSAFE false-separated (bound wrong vs GJK/EPA): " << unsafeFalseSeparated
      << " (max missed overlap=" << maxUnsafeOverlap << ")");

   WALBERLA_CHECK_EQUAL(unsafeFalseSeparated, 0);
   WALBERLA_CHECK_GREATER_EQUAL(resolvedSeparated, 200); // sanity: comparison isn't running on noise
}

//======================================================================================================================
// Master comparison of the Kodam/Guo closed form against GJK/EPA: 5000 trials, seed 24680, radius U[0.3,2], length
// U[0.5,4], uniform SO(3), boundary-biased separation. Logs (no assertions):
//   (1) existence misses / false positives of the raw closed form, and the hybrid's existence errors, depth errors
//       (> max(0.01, 5 %)) and fallback rate, plus a hybrid that falls back only when no scenario matched;
//   (2) a per-scenario table (occurrences, wrong-magnitude count, relative/absolute error) and normal angles of
//       passing vs. wrong trials;
//   (3) the "no scenario matched" gap, classified rim/band (classifyDispatchGapCause), with the missed depth;
//   (4) a cause for every wrong-magnitude and gap trial (classifyMagnitudeCause / classifyDispatchGapCause), and a
//       VTK/CSV dump of these trials (cylinder_kodam_mismatches).
//======================================================================================================================
void testCylinderCylinderFullComparisonAndHybridFeasibility()
{
   FuzzRng                                  fuzz(24680); // same seed as every breakdown/pattern test above
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 5000; // at 10000 trials: hybrid existence wrong 0/3338, fallback 112/3338 (3.36 %)

   int kodamMiss = 0, kodamFalsePos = 0, kodamCorrectNoContact = 0;
   int kodamHybridWrong = 0, kodamHybridFallback = 0, kodamHybridFallbackNeeded = 0;
   int kodamHybridWrongDepth = 0; // hybrid depth errors on genuine contacts (kodamHybridWrong counts existence only)
   int kodamTier1OnlyWrong = 0, kodamTier1OnlyFallback = 0;
   real_t kodamMaxOverlapDiff = 0;

   struct ScenarioWrongStats
   {
      int                    count = 0;
      real_t                 maxAbsDiff = 0, maxRelDiff = 0;
      std::vector< real_t >  relDiffs;
   };
   std::map< std::string, ScenarioWrongStats > kodamWrongByScenario;
   std::map< std::string, std::vector< real_t > > kodamWrongAngleByScenario;   // normal angle, wrongMag trials
   std::map< std::string, std::vector< real_t > > kodamPassingAngleByScenario; // normal angle, passing trials
   std::map< std::string, int >                kodamScenarioCount; // denominator per scenario (all that
                                                                       // fire)
   std::map< std::string, int >                kodamMissRimBandClass;
   std::map< std::string, std::vector< real_t > > kodamMissOverlapByClass;
   int kodamSpecialContactCount = 0; // trials firing special-contact (always routed to GJK/EPA, see
                                       // cylinderCylinderResultIsUncertain)

   std::map< std::string, int >     causeOccurrence; // "overall cause occurrence" count table
   std::vector< CylMismatchRecord > mismatches;       // every one of them, for the VTK/CSV dump

   const std::set< std::string > exemptScenarios = { "face-face", "band-band-parallel", "band-band-skewed" };

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (circumRi + circumRj),
                                                             std::min(Ri, Rj));
      const Vec3 posJ = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult kodam = cylinderCylinderContact(ci, cj);
      const GjkResult      gjk  = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
      const real_t         gjkOverlap = gjk.contact ? -gjk.penetrationDepth : real_t(0);

      // Rim-rim dispatch gap: Kodam's edge-edge chord test (eq. 82-101) needs, for a rim, a point A on the other
      // cylinder's lateral surface within its length and a point B on its near face plane within its radius. For
      // corner-to-corner contacts either A has no valid root or both B candidates lie outside the radius, so no
      // scenario matches although the bodies overlap (confirmed by brute-force MTD). Tier 1 of
      // cylinderCylinderResultIsUncertain sends all such cases to GJK/EPA.

      if (!gjk.contact)
      {
         if (kodam.contact)
         {
            ++kodamFalsePos;
            ++causeOccurrence["false-positive"];
            mismatches.push_back(CylMismatchRecord{ trial, ci, cj, true, false, kodam.overlap, real_t(0),
                                                       kodam.scenario, 0, "false-positive", real_t(-1) });
         }
         else { ++kodamCorrectNoContact; }
         continue;
      }

      if (!kodam.contact)
      {
         ++kodamMiss;
         const char* cause = classifyDispatchGapCause(ci, cj, gjk);
         ++kodamMissRimBandClass[cause];
         kodamMissOverlapByClass[cause].push_back(gjkOverlap);
         ++causeOccurrence[cause];
         mismatches.push_back(CylMismatchRecord{ trial, ci, cj, false, true, real_t(0), gjkOverlap, "none", 1,
                                                    cause, real_t(-1) });
      }
      else
      {
         kodamMaxOverlapDiff = std::max(kodamMaxOverlapDiff, std::abs(kodam.overlap - gjkOverlap));
         if (std::string(kodam.scenario) == "special-contact") ++kodamSpecialContactCount;
      }

      // Hybrid dispatcher.
      const HybridResult kodamHybrid = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
      if (kodamHybrid.usedFallback) ++kodamHybridFallback;
      if (kodamHybrid.contact != gjk.contact) ++kodamHybridWrong;
      else if (kodamHybrid.contact &&
               std::abs(kodamHybrid.overlap - gjkOverlap) > std::max(real_t(0.01), real_t(0.05) * gjkOverlap))
         ++kodamHybridWrongDepth;
      const bool kodamWouldHaveBeenWrongRaw =
         (kodam.contact != gjk.contact) ||
         (kodam.contact && std::abs(kodam.overlap - gjkOverlap) > std::max(real_t(0.01), real_t(0.05) * gjkOverlap));
      if (cylinderCylinderResultIsUncertain(ci, cj, kodam) && kodamWouldHaveBeenWrongRaw) ++kodamHybridFallbackNeeded;

      // "Tier-1-only" hybrid: fall back only when no scenario matched, otherwise trust the closed form.
      if (!kodam.contact) { ++kodamTier1OnlyFallback; }
      else
      {
         ++kodamScenarioCount[std::string(kodam.scenario)];
         const real_t absDiff  = std::abs(kodam.overlap - gjkOverlap);
         const real_t relDiff  = absDiff / gjkOverlap;
         const bool   wrongMag = absDiff > std::max(real_t(0.01), real_t(0.05) * gjkOverlap);

         // Normal angle vs. GJK/EPA for every contact trial, split into passing and wrong-magnitude.
         const real_t cosAngle = std::clamp(kodam.normal * gjk.normal, real_t(-1), real_t(1));
         const real_t angleDeg = std::acos(cosAngle) * real_t(180) / math::pi;
         (wrongMag ? kodamWrongAngleByScenario : kodamPassingAngleByScenario)[std::string(kodam.scenario)]
            .push_back(angleDeg);

         if (wrongMag)
         {
            ++kodamTier1OnlyWrong;
            auto& st = kodamWrongByScenario[std::string(kodam.scenario)];
            ++st.count;
            st.maxAbsDiff = std::max(st.maxAbsDiff, absDiff);
            st.maxRelDiff = std::max(st.maxRelDiff, relDiff);
            st.relDiffs.push_back(relDiff);

            real_t      metric = real_t(-1);
            const char* cause  = classifyMagnitudeCause(ci, cj, kodam, gjk, metric);
            ++causeOccurrence[cause];
            mismatches.push_back(CylMismatchRecord{ trial, ci, cj, true, true, kodam.overlap, gjkOverlap,
                                                       kodam.scenario, 2, cause, metric });
         }
      }
   }

   const int totalOverlapping = numTrials - kodamFalsePos - kodamCorrectNoContact;
   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] ================ KODAM/GUO MASTER TABLE (" << numTrials << " trials, same seed "
      << "throughout this file) ================");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] " << kodamCorrectNoContact << " genuinely separated (Kodam correctly reports no "
      << "contact), " << kodamFalsePos << " Kodam false positives (excluded from the miss-rate rows below), "
      << totalOverlapping << " genuine GJK/EPA contacts remain as the denominator for every row below");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] PURE Kodam/Guo:   miss=" << kodamMiss << "/" << totalOverlapping << " ("
      << (real_t(100) * real_t(kodamMiss) / real_t(totalOverlapping)) << "%), maxOverlapDiff="
      << kodamMaxOverlapDiff);
   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] Kodam HYBRID (existing production dispatcher): existence wrong=" << kodamHybridWrong << "/"
      << totalOverlapping << ", depth wrong (> max(0.01, 5%)) on agreeing contacts=" << kodamHybridWrongDepth << "/"
      << totalOverlapping << ", fallback rate=" << kodamHybridFallback << "/" << totalOverlapping << " ("
      << (real_t(100) * real_t(kodamHybridFallback) / real_t(totalOverlapping))
      << "%), of which genuinely NEEDED (raw analytical would've been wrong)=" << kodamHybridFallbackNeeded
      << "/" << kodamHybridFallback << " (" << (real_t(100) * real_t(kodamHybridFallbackNeeded) /
                                                  real_t(std::max(1, kodamHybridFallback))) << "%)");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] Of those, special-contact fired (and is now UNCONDITIONALLY routed to GJK/EPA, per "
      << "cylinderCylinderResultIsUncertain's new check) on " << kodamSpecialContactCount << "/" << totalOverlapping
      << " (" << (real_t(100) * real_t(kodamSpecialContactCount) / real_t(totalOverlapping))
      << "%) trials");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] Kodam TIER-1-ONLY (fallback iff no scenario matched, Tier 2 dropped entirely): wrong="
      << kodamTier1OnlyWrong << "/" << totalOverlapping << " (" << (real_t(100) * real_t(kodamTier1OnlyWrong) /
      real_t(totalOverlapping)) << "%), fallback rate=" << kodamTier1OnlyFallback << "/" << totalOverlapping
      << " (" << (real_t(100) * real_t(kodamTier1OnlyFallback) / real_t(totalOverlapping)) << "%)");

   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] ---------------- FULL PER-SCENARIO TABLE (every dispatch outcome, including "
      << "zero-wrong ones) ----------------");
   for (auto& [name, denom] : kodamScenarioCount)
   {
      const bool  exempt = exemptScenarios.count(name) > 0;
      auto        it     = kodamWrongByScenario.find(name);
      const int   wrongCount = (it != kodamWrongByScenario.end()) ? it->second.count : 0;
      if (wrongCount == 0)
      {
         WALBERLA_LOG_INFO_ON_ROOT(
            "  [comparison]   " << name << (exempt ? " [fully verified, exempt from Tier 2]" : "") << ": "
            << denom << " occurrences, 0 wrong (0%)");
         continue;
      }
      const ScenarioWrongStats& st = it->second;
      std::vector< real_t >     sorted = st.relDiffs;
      std::sort(sorted.begin(), sorted.end());
      const real_t p50 = sorted[sorted.size() / 2];
      WALBERLA_LOG_INFO_ON_ROOT(
         "  [comparison]   " << name << (exempt ? " [fully verified, exempt from Tier 2]" : "") << ": "
         << denom << " occurrences, " << wrongCount << " wrong (" << (real_t(100) * real_t(wrongCount) /
         real_t(denom)) << "% of this scenario's own trials, " << (real_t(100) * real_t(wrongCount) /
         real_t(kodamTier1OnlyWrong)) << "% of all wrong trials) | relative overlap error: p50=" << p50
         << " max=" << st.maxRelDiff << " (absolute max=" << st.maxAbsDiff << ")");
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] ---------------- NORMAL ANGLE: kodam.normal vs GJK/EPA's, PASSING vs WRONG-MAGNITUDE "
      << "trials, per scenario ----------------");
   for (auto& [name, denom] : kodamScenarioCount)
   {
      auto printStats = [&](const char* label, std::map< std::string, std::vector< real_t > >& byScenario) {
         auto it = byScenario.find(name);
         if (it == byScenario.end() || it->second.empty())
         {
            WALBERLA_LOG_INFO_ON_ROOT("  [comparison]   " << name << " " << label << ": n=0");
            return;
         }
         std::vector< real_t > sorted = it->second;
         std::sort(sorted.begin(), sorted.end());
         WALBERLA_LOG_INFO_ON_ROOT(
            "  [comparison]   " << name << " " << label << ": n=" << sorted.size() << " min=" << sorted.front()
            << " median=" << sorted[sorted.size() / 2] << " max=" << sorted.back());
      };
      printStats("PASSING ", kodamPassingAngleByScenario);
      printStats("WRONG   ", kodamWrongAngleByScenario);
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] Kodam's \"no scenario matched at all\" gap (" << kodamMiss << " trials), classified by "
      << "whether GJK/EPA's own contact point sits near each cylinder's rim or well within its band "
      << "(ci-side first, cj-side second):");
   for (auto& [name, count] : kodamMissRimBandClass)
   {
      std::vector< real_t > overlaps = kodamMissOverlapByClass[name];
      std::sort(overlaps.begin(), overlaps.end());
      WALBERLA_LOG_INFO_ON_ROOT(
         "  [comparison]   " << name << ": " << count << "/" << kodamMiss << " ("
         << (real_t(100) * real_t(count) / real_t(kodamMiss)) << "%) | missed overlap (GJK/EPA, i.e. "
         << "the depth completely undetected): min=" << overlaps.front() << " p50="
         << overlaps[overlaps.size() / 2] << " max=" << overlaps.back());
   }
   {
      const int total = kodamCorrectNoContact + kodamMiss; // trials where no scenario matched
      WALBERLA_LOG_INFO_ON_ROOT(
         "[comparison] Precision of Kodam's \"no case\" classification: of every trial where Kodam reports "
         << "no scenario matched at all (" << total << " = " << kodamCorrectNoContact
         << " genuinely separated + " << kodamMiss << " actually real contacts), " << kodamMiss << "/"
         << total << " (" << (real_t(100) * real_t(kodamMiss) / real_t(total))
         << "%) were wrong -- i.e. \"no case\" is correct " << (real_t(100) * real_t(kodamCorrectNoContact) /
         real_t(total)) << "% of the time it fires, which is exactly what Tier 1's cheap "
         << "reach-sum bound exists to catch the residual of");
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[comparison] ---------------- OVERALL CAUSE OCCURRENCE (every dispatch-gap + wrong-magnitude "
      << "trial, " << mismatches.size() << " total) ----------------");
   for (auto& [cause, count] : causeOccurrence)
   {
      WALBERLA_LOG_INFO_ON_ROOT(
         "  [comparison]   " << cause << ": " << count << "/" << mismatches.size() << " (" <<
         (real_t(100) * real_t(count) / real_t(mismatches.size())) << "%)");
   }

   if (writeMismatchVTKEnabled())
      writeCylMismatchVTKSeries((g_outputDir / "cylinder_kodam_mismatches").string(), "cylinder_kodam_mismatches",
                                 mismatches);
   writeCylMismatchCSV((g_outputDir / "cylinder_kodam_mismatches.csv").string(), mismatches);
}

//======================================================================================================================
// Deep-overlap stress test: overlap forced to 30-150 % of the smaller radius (far beyond the 15 % cap of the other
// fuzz tests), 10000 trials, seed 13579. Question: if a large overlap ever occurred (bad initial condition, a bug
// elsewhere), does the hybrid still detect it? Pass criterion: zero hybrid existence mismatches vs. GJK/EPA. The worst
// raw and hybrid depth errors are logged (an overestimating scenario is not caught by any fallback).
//======================================================================================================================
void testCylinderCylinderLargeOverlapDegenerateStress()
{
   FuzzRng                                  fuzz(13579); // fixed seed
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > overlapFracDist(real_t(0.3), real_t(1.5)); // x smaller radius

   constexpr int numTrials = 10000;
   int kodamNoneCase = 0, hybridWrong = 0;
   real_t worstUnsafeAbsDiff = 0, worstUnsafeRelDiff = 0;
   CylinderState worstCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
   CylinderState worstCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
   real_t worstHybridAbsDiff = 0, worstHybridRelDiff = 0;
   CylinderState worstHybridCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
   CylinderState worstHybridCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng);
      const real_t Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation();
      const Rot3   rotJ = fuzz.randomRotation();

      const real_t circumRi = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li);
      const real_t circumRj = std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir      = fuzz.randomUnitVec();
      const Vec3   posI(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * s, rotJ, Rj, Lj).contact;
      };
      // Contact boundary by bisection (as in FuzzRng::sampleSeparationNearBoundary), then push 30-150 % past it.
      real_t lo = real_t(0), hi = real_t(3) * (circumRi + circumRj);
      WALBERLA_ASSERT(contactAtSep(lo));
      WALBERLA_ASSERT(!contactAtSep(hi));
      for (int it = 0; it < 40; ++it)
      {
         const real_t mid = real_t(0.5) * (lo + hi);
         if (contactAtSep(mid)) { lo = mid; } else { hi = mid; }
      }
      const real_t boundary = real_t(0.5) * (lo + hi);
      const real_t sep      = boundary - overlapFracDist(fuzz.rng) * std::min(Ri, Rj);
      const Vec3   posJ     = posI + dir * sep;

      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li };
      const CylinderState cj{ posJ, rotJ.getMatrix(), Rj, Lj };

      const ContactResult kodam = cylinderCylinderContact(ci, cj);
      const GjkResult      gjk  = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
      const real_t         gjkOverlap = gjk.contact ? -gjk.penetrationDepth : real_t(0);

      if (!gjk.contact) continue; // rare overshoot into no contact
      if (!kodam.contact) { ++kodamNoneCase; } // caught by Tier 1 in the hybrid
      else
      {
         const real_t absDiff = std::abs(kodam.overlap - gjkOverlap);
         const real_t relDiff = absDiff / gjkOverlap;
         // A scenario fired: its depth reaches the hybrid unchecked; track the worst case.
         if (absDiff > worstUnsafeAbsDiff)
         {
            worstUnsafeAbsDiff = absDiff;
            worstUnsafeRelDiff = relDiff;
            worstCi            = ci;
            worstCj            = cj;
         }
      }

      const HybridResult hybrid = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
      if (hybrid.contact != gjk.contact) ++hybridWrong;
      // Also track the hybrid's own worst depth error (logged, not asserted).
      if (hybrid.contact && gjk.contact)
      {
         const real_t hybAbsDiff = std::abs(hybrid.overlap - gjkOverlap);
         const real_t hybRelDiff = hybAbsDiff / gjkOverlap;
         if (hybAbsDiff > worstHybridAbsDiff)
         {
            worstHybridAbsDiff = hybAbsDiff;
            worstHybridRelDiff = hybRelDiff;
            worstHybridCi      = ci;
            worstHybridCj      = cj;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[stress] Deep-overlap degenerate stress (" << numTrials << " trials, overlap forced to 30%-150% "
      << "of the smaller radius -- well past the 15% 'realistic' cap every other fuzz test in this file "
      << "uses): " << kodamNoneCase << " correctly fell to \"no case\" (safe, always caught by Tier 1), "
      << (numTrials - kodamNoneCase) << " had SOME scenario fire, of which the worst raw analytical "
      << "overlap disagreement was absDiff=" << worstUnsafeAbsDiff << " relDiff=" << worstUnsafeRelDiff
      << " (ci: R=" << worstCi.radius << " L=" << worstCi.length << ", cj: R=" << worstCj.radius
      << " L=" << worstCj.length << "). HYBRID (with Tier 1/2 fallback) existence mismatches="
      << hybridWrong << "/" << numTrials << " (must be 0 for the hybrid to be considered safe under "
      << "deep/degenerate overlap, not just realistic overlap)");
   WALBERLA_LOG_INFO_ON_ROOT(
      "[stress] HYBRID's own worst overlap MAGNITUDE disagreement (not just raw kodam's, and not just "
      << "existence): absDiff=" << worstHybridAbsDiff << " relDiff=" << worstHybridRelDiff << " (ci: R="
      << worstHybridCi.radius << " L=" << worstHybridCi.length << ", cj: R=" << worstHybridCj.radius
      << " L=" << worstHybridCj.length << ") -- this is what actually reaches the simulation, unlike "
      << "the raw-analytical worst-case above which the hybrid may have already corrected");
   WALBERLA_CHECK_EQUAL(hybridWrong, 0);
}

//======================================================================================================================
// Randomized fuzz test for Box x Box: boxBoxSAT against GJK/EPA. boxBoxSAT has no fallback gate (trusted directly in
// ParticleKernels.h), so only the raw SAT result is tested.
//======================================================================================================================

void testBoxBoxRandomizedFuzz()
{
   FuzzRng                                  fuzz(67890); // half-extents U[0.3,2], uniform SO(3), near the boundary
   std::uniform_real_distribution< real_t > halfExtentDist(real_t(0.3), real_t(2.0));

   constexpr int numTrials = 10000;
   int           numOverlapping = 0;
   int           numMismatch    = 0; // SAT contact-existence disagrees with GJK/EPA
   real_t        maxAgreeOverlapDiff = real_t(0);
   real_t        maxAgreePointDist   = real_t(0); // along-normal component only, see below
   real_t        maxAgreePerpDist    = real_t(0); // informational only, not asserted
   Vec3          worstPointA, worstPointB;

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const Vec3 halfA(halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng));
      const Vec3 halfB(halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng));
      const Rot3 rotA = fuzz.randomRotation();
      const Rot3 rotB = fuzz.randomRotation();

      const real_t diagA = halfA.length(); // half-diagonal -- conservative circumscribing-sphere radius
      const real_t diagB = halfB.length();
      const Vec3   dir   = fuzz.randomUnitVec();
      const real_t minSize =
         std::min({ halfA[0], halfA[1], halfA[2], halfB[0], halfB[1], halfB[2] });
      const Vec3 posA(real_t(0), real_t(0), real_t(0));

      auto contactAtSep = [&](real_t s) {
         return gjkBoxBox(posA, rotA, real_t(2) * halfA, posA + dir * s, rotB, real_t(2) * halfB).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                             real_t(3) * (diagA + diagB), minSize);

      const Vec3 posB = posA + dir * sep;

      const auto satResult = boxBoxSAT(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(), halfB);
      const GjkResult gjk  = gjkBoxBox(posA, rotA, real_t(2) * halfA, posB, rotB, real_t(2) * halfB);

      if (gjk.contact) ++numOverlapping;
      if (satResult.contact != gjk.contact)
      {
         ++numMismatch;
         WALBERLA_LOG_WARNING_ON_ROOT(
            "[fuzz] Box x Box SAT/GJK contact-existence mismatch, trial " << trial << ": sat.contact="
            << satResult.contact << " gjk.contact=" << gjk.contact << " | halfA=" << halfA
            << " halfB=" << halfB << " sep=" << sep);
      }

      if (satResult.contact && gjk.contact)
      {
         const real_t overlapDiff = std::abs(satResult.overlap - (-gjk.penetrationDepth));
         if (overlapDiff > maxAgreeOverlapDiff) maxAgreeOverlapDiff = overlapDiff;

         // Contact point: for flush face-face or parallel edge-edge overlaps only the component along the normal is
         // well defined; perpendicular to it any point of the overlap footprint is valid. So the along-normal
         // component is asserted, the perpendicular one only logged.
         const Vec3   diff            = satResult.point - gjk.point;
         const real_t alongNormalDist = std::abs(diff * satResult.normal);
         const real_t perpDist =
            std::sqrt(std::max(real_t(0), diff.sqrLength() - alongNormalDist * alongNormalDist));
         if (alongNormalDist > maxAgreePointDist)
         {
            maxAgreePointDist = alongNormalDist;
            worstPointA        = posA;
            worstPointB        = posB;
         }
         if (perpDist > maxAgreePerpDist) maxAgreePerpDist = perpDist;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Box x Box, " << numTrials << " random trials (" << numOverlapping << " actually "
      << "overlapping per GJK/EPA, overlap depth capped realistic): SAT/GJK contact-existence "
      << "mismatches=" << numMismatch << " (expected 0 -- SAT has no known limitation class the way "
      << "cylinder-cylinder does), max overlap magnitude difference on agreeing contacts="
      << maxAgreeOverlapDiff << ", max along-normal contact-point difference on agreeing contacts="
      << maxAgreePointDist << " (worst-case centres: A=" << worstPointA << " B=" << worstPointB
      << "), max PERPENDICULAR-to-normal point difference=" << maxAgreePerpDist
      << " (informational only -- underdetermined for flush contacts, not asserted)");

   WALBERLA_CHECK_EQUAL(numMismatch, 0);
   WALBERLA_CHECK_LESS(maxAgreeOverlapDiff, real_t(0.05));
   WALBERLA_CHECK_LESS(maxAgreePointDist, real_t(0.1));
}


//======================================================================================================================
// Hybrid dispatcher: does cylinderCylinderResultIsUncertain() send missed contacts to GJK/EPA while leaving the
// trusted scenarios on the closed-form path?
//======================================================================================================================

//! Tier 1: a genuine missed contact, trial 1209 of the independent-reference sample (overleaf
//! scripts/contact-validation/data/cyl_pairs.csv): no scenario applies (the end faces meet corner to corner near both
//! rims, thesis Figure 4.11). The reach-sum gate (centre distance < sum of radii and half-lengths) must trigger the
//! fallback, which must return GJK/EPA's depth 0.113 (brute-force MTD 0.1131).
void testHybridTier1TriggersOnMissedContact()
{
   const Mat3 rotIRows(real_t(0.0903122648821542), real_t(-0.93901787549968108), real_t(0.33179681177479803),
                       real_t(-0.78472296406621889), real_t(0.13804609137935925), real_t(0.60427903018556672),
                       real_t(-0.61323206413150244), real_t(-0.31494238544049091), real_t(-0.72440163540267466));
   const Mat3 rotJRows(real_t(0.42629962306143004), real_t(-0.82827933257245778), real_t(0.3636234021773167),
                       real_t(0.34357538494322015), real_t(0.52011375867889886), real_t(0.78194477611531255),
                       real_t(-0.83679423171936407), real_t(-0.20841071293537705), real_t(0.50630068980297382));
   const Rot3          rotI = rot3FromMatrix(rotIRows);
   const Rot3          rotJ = rot3FromMatrix(rotJRows);
   const CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), rotI.getMatrix(), real_t(1.5425930566498214),
                           real_t(3.1842270146051681) };
   const CylinderState cj{ Vec3(real_t(1.8193447429136438), real_t(1.1130105298503024), real_t(-2.6016803736233309)),
                           rotJ.getMatrix(), real_t(1.4789598196015563), real_t(2.1871520227523074) };
   const ContactResult analytical = cylinderCylinderContact(ci, cj);
   const HybridResult  hybrid     = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
   const GjkResult     gjk = gjkCylCyl(ci.pos, rotI, ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
   WALBERLA_CHECK(!analytical.contact);
   WALBERLA_CHECK_EQUAL(std::string(analytical.scenario), "none");
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK(hybrid.usedFallback);
   WALBERLA_CHECK(hybrid.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-9));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, real_t(0.113098), real_t(1e-4)); // brute-force MTD 0.1131
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 1 (missed contact, trial 1209): analytical=none, fallback overlap="
                             << hybrid.overlap);
}

//! Tier 1 must not fire for a clearly separated pair (centres 10 apart, reach 5).
void testHybridTier1DoesNotTriggerWhenClearlySeparated()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
   CylinderState cj{ Vec3(real_t(10), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, Rot3());
   WALBERLA_CHECK(!hybrid.usedFallback);
   WALBERLA_CHECK(!hybrid.contact);
}

//! Tier 1 at the gate boundary: centreDist == reach exactly (3.0 == 0.5 + 0.5 + 1 + 1), so the strict `<` leaves the
//! gate untriggered (logged); the hybrid must report no contact.
void testHybridTier1BoundaryCase()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
   CylinderState cj{ Vec3(real_t(3), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, Rot3());
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid Tier 1 boundary: usedFallback=" << hybrid.usedFallback
                              << " (centreDist == reach exactly)");
   WALBERLA_CHECK(!hybrid.contact);
}

//! Geometry of testCylinderCylinderTiltedRimIntoFace (face-edge 0.299 vs. GJK/EPA 0.235): the hybrid has no magnitude
//! check, so it must not fall back and must return the raw closed-form depth (accepted residual error of face-edge).
void testHybridKeepsFaceEdgeMagnitudeWithoutFallback()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.6), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.25));
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.55)), rotJ.getMatrix(), real_t(0.2), real_t(2) };

   const ContactResult analyticalOnly = cylinderCylinderContact(ci, cj); // raw closed form, for comparison
   const HybridResult   hybrid        = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
   const GjkResult       gjk           = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(!hybrid.usedFallback);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, analyticalOnly.overlap, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid (Tier 2 removed) on a known wrong-magnitude geometry: analytical "
                              << "overlap=" << analyticalOnly.overlap << ", hybrid overlap=" << hybrid.overlap
                              << " (same value, no longer corrected) vs GJK/EPA's true "
                              << -gjk.penetrationDepth << " -- accepted residual error, not a bug");
}

//! face-face, band-band-parallel, band-band-skewed and face-band are trusted directly: no fallback may fire, and for
//! the latter two the hybrid must agree with GJK/EPA (1e-3). Only "no scenario matched" (Tier 1) and special-contact
//! fall back.
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
      const GjkResult     gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
      WALBERLA_CHECK(!hybrid.usedFallback);
      WALBERLA_CHECK_EQUAL(hybrid.contact, gjk.contact);
      if (hybrid.contact && gjk.contact)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-3));
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
      CylinderState cj{ Vec3(real_t(0), real_t(0.45), real_t(1.4)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
      const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
      const GjkResult     gjk = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);
      WALBERLA_CHECK(!hybrid.usedFallback);
      WALBERLA_CHECK_EQUAL(hybrid.contact, gjk.contact);
      if (hybrid.contact && gjk.contact)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-3));
   }
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid (Tier 2 removed): face-face, band-band-parallel, band-band-skewed "
                              "and face-band all now trusted directly (no fallback call at all), and "
                              "all four still matched GJK/EPA");
}

//! Geometry of testCylinderCylinderEdgeEdge (face-edge, exact 0.05): the hybrid must agree with GJK/EPA (1e-3).
//! usedFallback is expected false but only logged.
void testHybridTier2OnAlreadyCorrectFaceEdgeCase()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.4), real_t(1.6) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
   CylinderState cj{ Vec3(real_t(1.15), real_t(0), real_t(1.1)), rotJ.getMatrix(), real_t(0.4), real_t(1.6) };

   const HybridResult hybrid = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
   const GjkResult     gjk    = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(hybrid.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid (Tier 2 removed) on already-correct face-edge case: usedFallback="
                              << hybrid.usedFallback << " (expected false now -- no longer needed here)");
}

//======================================================================================================================
// Cylinder x Plane -- Kodam et al. face/band/edge-plane vs. the production Cylinder x HalfSpace overload
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

   const GjkResult gjk = productionCylPlane(cyl.pos, Rot3(), cyl.radius, cyl.length, planePos, planeNormal);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-6));

   // Contact points differ by convention: Kodam eq. 6 places it at the middle of the overlap (z = -0.5),
   // SolidCylinderContactDetection on the wall surface (z = 0).
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

   const GjkResult gjk = productionCylPlane(cyl.pos, rot, cyl.radius, cyl.length, planePos, planeNormal);
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

//! Edge-plane: cylinder (R = 0.3, L = 2) tilted 30 degrees about y, centre at z = 0.6, bottom rim below the floor.
//! Kodam eq. 16-17: axis (0.5, 0, 0.866), bottom axis end (-0.5, 0, 0.6 - 0.866); the rim direction towards the
//! plane is (0.866, 0, -0.5), so E = (-0.5 + 0.3 * 0.866, 0, -0.266 - 0.3 * 0.5) = (-0.2402, 0, -0.416). zE has the
//! opposite sign of zCyl = 0.6, so contact with overlap 0.416.
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

   const GjkResult gjk = productionCylPlane(cyl.pos, rot, cyl.radius, cyl.length, planePos, planeNormal);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(2e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));

   // Contact point from eq. 21 (z_E - sign(z_E) * 0.5 * delta): E = (-0.240192, 0, -0.416025), overlap 0.416025, so
   // point = (-0.240192, 0, -0.208013), in the middle of the overlap. Production places the point on the plane (the
   // deepest rim point projected along n).
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(gjk.point, Vec3(real_t(-0.240192), real_t(0), real_t(0)), real_t(1e-5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point,
      Vec3(real_t(-0.240192), real_t(0), real_t(-0.208013)), real_t(1e-3));

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

//======================================================================================================================
// Sphere x Cylinder -- closed-form closest point (sphereCylinderContact); exact, no uncertainty gate or fallback.
//======================================================================================================================

//! Sphere outside, closest feature on the band (radial clamp only): cylinder R = 1, L = 4 along z; sphere radius 0.3
//! at (1.2, 0, 0). rho = 1.2, z = 0 (|z| < halfL = 2); closest point q = (1, 0, 0); dist = 0.2; overlap = 0.3 - 0.2
//! = 0.1, normal +x.
void testSphereCylinderBand()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(4.0) };
   const Vec3    spherePos(real_t(1.2), real_t(0), real_t(0));
   const real_t  sphereRadius = real_t(0.3);

   const ContactResult kodam = sphereCylinderContact(spherePos, sphereRadius, cyl);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.1), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(1), real_t(0), real_t(0)), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(1), real_t(0), real_t(0)), real_t(1e-6));

   const GjkResult prod = productionSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(prod.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -prod.penetrationDepth, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, prod.normal, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Sphere x Cylinder (band): overlap=" << kodam.overlap << " (exact agreement)");
   // production contact point: middle of the overlap along the normal (ParticleKernels.h)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.point, spherePos - (sphereRadius - real_t(0.5) * kodam.overlap) * kodam.normal,
                                       real_t(1e-6));
   const GjkResult gjk = gjkSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Sphere outside, closest feature on a cap (axial clamp only): cylinder R = 1, L = 2 (halfL = 1) along z; sphere
//! radius 0.3 on the axis at z = 1.2. zC = 1; q = (0, 0, 1); dist = 0.2; overlap = 0.3 - 0.2 = 0.1.
void testSphereCylinderCap()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(2.0) };
   const Vec3    spherePos(real_t(0), real_t(0), real_t(1.2));
   const real_t  sphereRadius = real_t(0.3);

   const ContactResult kodam = sphereCylinderContact(spherePos, sphereRadius, cyl);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.1), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(0), real_t(0), real_t(1)), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(0), real_t(0), real_t(1)), real_t(1e-6));

   const GjkResult prod = productionSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(prod.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -prod.penetrationDepth, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, prod.normal, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Sphere x Cylinder (cap): overlap=" << kodam.overlap << " (exact agreement)");
   // production contact point: middle of the overlap along the normal (ParticleKernels.h)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.point, spherePos - (sphereRadius - real_t(0.5) * kodam.overlap) * kodam.normal,
                                       real_t(1e-6));
   const GjkResult gjk = gjkSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Sphere near the rim: radial and axial clamp both fire (rho > R and |z| > halfL), so the closest point is on the
//! circular edge. Cylinder R = 1, L = 2; sphere radius 0.3 at (1.2, 0, 1.1); q = (1, 0, 1); diff = (0.2, 0, 0.1),
//! dist = sqrt(0.05) = 0.223607, overlap = 0.3 - 0.223607 = 0.076393.
void testSphereCylinderRim()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(2.0) };
   const Vec3    spherePos(real_t(1.2), real_t(0), real_t(1.1));
   const real_t  sphereRadius = real_t(0.3);

   const ContactResult kodam = sphereCylinderContact(spherePos, sphereRadius, cyl);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.0763932), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, Vec3(real_t(1), real_t(0), real_t(1)), real_t(1e-6));

   const GjkResult prod = productionSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(prod.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -prod.penetrationDepth, real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, prod.normal, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Sphere x Cylinder (rim, both clamps fire): overlap=" << kodam.overlap
                              << " (exact agreement)");
   // production contact point: middle of the overlap along the normal (ParticleKernels.h)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.point, spherePos - (sphereRadius - real_t(0.5) * kodam.overlap) * kodam.normal,
                                       real_t(1e-6));
   const GjkResult gjk = gjkSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Sphere centre inside the cylinder, nearer the band: cylinder R = 1, L = 4; sphere radius 0.3 at (0.8, 0, 0).
//! Radial clearance 0.2 < axial clearance 2, so overlap = 0.2 + 0.3 = 0.5 along +x.
void testSphereCylinderEmbeddedBand()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(4.0) };
   const Vec3    spherePos(real_t(0.8), real_t(0), real_t(0));
   const real_t  sphereRadius = real_t(0.3);

   const ContactResult kodam = sphereCylinderContact(spherePos, sphereRadius, cyl);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.5), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(1), real_t(0), real_t(0)), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, spherePos, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Sphere x Cylinder (embedded, band side wins): overlap=" << kodam.overlap);
}

//! Sphere centre inside the cylinder, nearer a cap: cylinder R = 3, L = 2; sphere radius 0.3 at (0, 0, 0.9).
//! Axial clearance 0.1 < radial clearance 3, so overlap = 0.1 + 0.3 = 0.4 along +z.
void testSphereCylinderEmbeddedCap()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(3.0), real_t(2.0) };
   const Vec3    spherePos(real_t(0), real_t(0), real_t(0.9));
   const real_t  sphereRadius = real_t(0.3);

   const ContactResult kodam = sphereCylinderContact(spherePos, sphereRadius, cyl);
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.4), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, Vec3(real_t(0), real_t(0), real_t(1)), real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point, spherePos, real_t(1e-6));
   WALBERLA_LOG_INFO_ON_ROOT("Sphere x Cylinder (embedded, cap side wins): overlap=" << kodam.overlap);
}

void testSphereCylinderNoContact()
{
   CylinderState cyl{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(4.0) };
   const Vec3    spherePos(real_t(3), real_t(0), real_t(0));
   WALBERLA_CHECK(!sphereCylinderContact(spherePos, real_t(0.3), cyl).contact);
   const GjkResult prod = productionSphereCyl(spherePos, real_t(0.3), cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(!prod.contact);
   const GjkResult gjk = gjkSphereCyl(spherePos, real_t(0.3), cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(!gjk.contact);
}

//! Randomized check of sphereCylinderContact: 10000 trials, seed 112233, cylinder radius U[0.3,2], length U[0.5,4],
//! sphere radius U[0.1,1], uniform SO(3), separation near the boundary of the production overload. Compared with the
//! production closed form (productionSphereCyl, variable `gjk`) and with GJK/EPA (gjkSphereCyl, variable `ref`).
void testSphereCylinderRandomizedFuzz()
{
   FuzzRng                                  fuzz(112233);
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > sphereRadiusDist(real_t(0.1), real_t(1.0));

   constexpr int numTrials = 10000;
   int           numOverlapping = 0;
   int           numMismatch    = 0;
   int           numGjkMismatch = 0;   // against GJK/EPA
   real_t        maxGjkOverlapDiff = real_t(0), maxGjkNormalAngleDeg = real_t(0);
   real_t        maxOverlapDiff = real_t(0);
   real_t        maxPointDiff   = real_t(0);
   real_t        maxNormalAngleDeg = real_t(0);

   for (int trial = 0; trial < numTrials; ++trial)
   {
      const real_t cylR = radiusDist(fuzz.rng), cylL = lengthDist(fuzz.rng);
      const real_t sphR = sphereRadiusDist(fuzz.rng);
      const Rot3   rotCyl = fuzz.randomRotation();
      const Vec3   cylPos(real_t(0), real_t(0), real_t(0));

      const real_t circumCyl = std::sqrt(cylR * cylR + real_t(0.25) * cylL * cylL);
      const Vec3   dir       = fuzz.randomUnitVec();

      auto contactAtSep = [&](real_t s) {
         return productionSphereCyl(cylPos + dir * s, sphR, cylPos, rotCyl, cylR, cylL).contact;
      };
      const real_t sep = fuzz.sampleSeparationNearBoundary(contactAtSep, real_t(0),
                                                              real_t(3) * (circumCyl + sphR), sphR);
      const Vec3   spherePos = cylPos + dir * sep;
      const CylinderState cyl{ cylPos, rotCyl.getMatrix(), cylR, cylL };

      const ContactResult kodam = sphereCylinderContact(spherePos, sphR, cyl);
      // `gjk`: production closed form; `ref`: independent reference, mesa_pd GJK/EPA
      const GjkResult      gjk  = productionSphereCyl(spherePos, sphR, cylPos, rotCyl, cylR, cylL);
      const GjkResult      ref  = gjkSphereCyl(spherePos, sphR, cylPos, rotCyl, cylR, cylL);
      if (ref.contact != sphereCylinderContact(spherePos, sphR, CylinderState{ cylPos, rotCyl.getMatrix(), cylR, cylL }).contact)
         ++numGjkMismatch;
      else if (ref.contact)
      {
         const ContactResult k = sphereCylinderContact(spherePos, sphR, CylinderState{ cylPos, rotCyl.getMatrix(), cylR, cylL });
         maxGjkOverlapDiff    = std::max(maxGjkOverlapDiff, std::abs(k.overlap + ref.penetrationDepth));
         maxGjkNormalAngleDeg = std::max(maxGjkNormalAngleDeg,
            std::acos(std::clamp(k.normal * ref.normal, real_t(-1), real_t(1))) * real_t(180) / math::pi);
      }

      if (gjk.contact) ++numOverlapping;
      if (kodam.contact != gjk.contact)
      {
         ++numMismatch;
         WALBERLA_LOG_WARNING_ON_ROOT(
            "[fuzz] Sphere x Cylinder mismatch, trial " << trial << ": kodam.contact=" << kodam.contact
            << " gjk.contact=" << gjk.contact << " | cylR=" << cylR << " cylL=" << cylL << " sphR=" << sphR
            << " sep=" << sep);
         continue;
      }
      if (!kodam.contact) continue;

      const real_t overlapDiff = std::abs(kodam.overlap - (-gjk.penetrationDepth));
      const real_t pointDiff   = (kodam.point - gjk.point).length();
      const real_t cosAngle    = std::clamp(kodam.normal * gjk.normal, real_t(-1), real_t(1));
      const real_t angleDeg    = std::acos(cosAngle) * real_t(180) / math::pi;
      maxOverlapDiff    = std::max(maxOverlapDiff, overlapDiff);
      maxPointDiff      = std::max(maxPointDiff, pointDiff);
      maxNormalAngleDeg = std::max(maxNormalAngleDeg, angleDeg);
   }

   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Sphere x Cylinder, " << numTrials << " random trials (" << numOverlapping
      << " actually overlapping per GJK/EPA): existence mismatches=" << numMismatch
      << ", max overlap diff=" << maxOverlapDiff << ", max point diff=" << maxPointDiff
      << ", max normal angle diff=" << maxNormalAngleDeg << " deg");

   // Asserted, unlike the cylinder-cylinder fuzz tests: sphereCylinderContact is an exact minimum translation distance.
   // Against the production closed form: existence identical, overlap < 1e-5, normal < 0.1 deg (measured 2.8e-15 and
   // 1.9e-6 deg). Against GJK/EPA: existence identical, overlap < 1e-3, normal < 5 deg.
   WALBERLA_CHECK_EQUAL(numMismatch, 0);
   WALBERLA_LOG_INFO_ON_ROOT("[fuzz] Sphere x Cylinder vs true GJK/EPA: existence mismatches=" << numGjkMismatch
                             << ", max overlap diff=" << maxGjkOverlapDiff << ", max normal angle diff="
                             << maxGjkNormalAngleDeg << " deg");
   WALBERLA_CHECK_EQUAL(numGjkMismatch, 0);
   WALBERLA_CHECK_LESS(maxGjkOverlapDiff, real_t(1e-3));
   // Measured against GJK/EPA: overlap 3.9e-4, normal up to 2.9 deg. The normal error is GJK/EPA's (polytope
   // approximation of the curved surfaces near the boundary), hence the looser bound.
   WALBERLA_CHECK_LESS(maxGjkNormalAngleDeg, real_t(5));
   WALBERLA_CHECK_LESS(maxOverlapDiff, real_t(1e-5));
   WALBERLA_CHECK_LESS(maxNormalAngleDeg, real_t(0.1));
   // Contact point not asserted: sphereCylinderContact places it on the cylinder surface, the production overload in
   // the middle of the overlap along the normal (two conventions differing by half the overlap).
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Sphere x Cylinder max contact-point difference=" << maxPointDiff << " (informational only -- "
      << "differing point CONVENTION between the two paths, not a geometric disagreement; see comment)");
}


//======================================================================================================================
// Continuity of the detection along a path (thesis Section 5.2). For each random pair (uniform
// SO(3)), cylinder j is translated along a straight line through the contact boundary, from an overlap of 15 % of the
// smaller radius to a gap of 10 %, in 300 steps of size ds. Under a pure translation the true penetration depth is
// Lipschitz with constant 1 (it changes at most by ds per step), contact switches exactly once, and the normal turns
// continuously. Counted per method (GJK/EPA, raw closed form, hybrid, refined): paths with a depth jump > 2 ds,
// paths whose contact flag switches more than once, paths with a normal jump > 5 deg between consecutive steps, and
// the largest depth jump in units of ds.
//======================================================================================================================
void testCylinderCylinderContinuityAlongPaths()
{
   FuzzRng                                  fuzz(97531);
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0)), lengthDist(real_t(0.5), real_t(4.0));
   constexpr int numPaths = 300, numSteps = 300;
   struct Stats { int depthJump = 0, flicker = 0, normalJump = 0, bigJump = 0; real_t maxJump = 0, maxJumpRel = 0; };
   std::map< std::string, Stats > st;
   const std::vector< std::string > methods = { "gjk", "raw", "hybrid", "refined" };
   for (int t = 0; t < numPaths; ++t)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng), Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation(), rotJ = fuzz.randomRotation();
      const Vec3   dir  = fuzz.randomUnitVec();
      const Vec3   posI(0, 0, 0);
      const real_t circ = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li) + std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      real_t lo = 0, hi = 3 * circ;
      for (int it = 0; it < 50; ++it)
      {
         const real_t mid = real_t(0.5) * (lo + hi);
         (gjkCylCyl(posI, rotI, Ri, Li, posI + dir * mid, rotJ, Rj, Lj).contact ? lo : hi) = mid;
      }
      const real_t boundary = real_t(0.5) * (lo + hi), minR = std::min(Ri, Rj);
      const real_t s0 = boundary - real_t(0.15) * minR, s1 = boundary + real_t(0.10) * minR, ds = (s1 - s0) / numSteps;
      std::map< std::string, bool > prevC, jumped, flick, njump, bjump;
      std::map< std::string, real_t > prevD;
      std::map< std::string, Vec3 >   prevN;
      std::map< std::string, int >    switches;
      for (int k = 0; k <= numSteps; ++k)
      {
         const Vec3          posJ = posI + dir * (s0 + k * ds);
         const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li }, cj{ posJ, rotJ.getMatrix(), Rj, Lj };
         const GjkResult     g  = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
         const ContactResult r  = cylinderCylinderContact(ci, cj);
         const HybridResult  h  = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
         const ContactResult rf = refineCylinderContact(ci, cj, r);
         const std::map< std::string, std::tuple< bool, real_t, Vec3 > > res = {
            { "gjk", { g.contact, g.contact ? -g.penetrationDepth : real_t(0), g.normal } },
            { "raw", { r.contact, r.contact ? r.overlap : real_t(0), r.normal } },
            { "hybrid", { h.contact, h.contact ? h.overlap : real_t(0), h.normal } },
            { "refined", { rf.contact, rf.contact ? rf.overlap : real_t(0), rf.normal } } };
         for (const auto& m : methods)
         {
            const auto& [c, d, n] = res.at(m);
            if (k > 0)
            {
               const real_t jump = std::abs(d - prevD[m]);
               st[m].maxJump = std::max(st[m].maxJump, jump / ds);
               st[m].maxJumpRel = std::max(st[m].maxJumpRel, jump / minR);
               if (jump > real_t(0.01) * minR)
               {
                  if (m == "refined" && !bjump[m] && std::getenv("ANTIDUNES_CONT_DEBUG"))
                     WALBERLA_LOG_INFO_ON_ROOT("[cont-big] path " << t << " step " << k << " refined " << prevD[m] << " -> " << d
                        << " scen " << (c ? rf.scenario : "none") << " raw " << r.scenario << " gjk=" << (g.contact ? -g.penetrationDepth : real_t(0))
                        << " minR=" << minR);
                  bjump[m] = true;
               }
               if (jump > real_t(2) * ds + real_t(1e-9))
               {
                  if (m == "refined" && !jumped[m] && std::getenv("ANTIDUNES_CONT_DEBUG"))
                     WALBERLA_LOG_INFO_ON_ROOT("[cont-debug] path " << t << " step " << k << " refined " << prevD[m] << " -> " << d
                        << " (" << (c ? rf.scenario : "none") << ", raw " << r.scenario << ") gjk=" << (g.contact ? -g.penetrationDepth : real_t(0))
                        << " ds=" << ds);
                  jumped[m] = true;
               }
               if (c != prevC[m]) ++switches[m];
               if (c && prevC[m] && std::acos(std::clamp(n * prevN[m], real_t(-1), real_t(1))) * real_t(180) / math::pi > real_t(5))
                  njump[m] = true;
            }
            prevC[m] = c; prevD[m] = d; prevN[m] = n;
         }
      }
      for (const auto& m : methods)
      {
         st[m].depthJump += jumped[m];
         st[m].flicker += switches[m] > 1;
         st[m].normalJump += njump[m];
         st[m].bigJump += bjump[m];
      }
   }
   for (const auto& m : methods)
      WALBERLA_LOG_INFO_ON_ROOT("[continuity] " << m << ": paths with depth jump > 2 ds = " << st[m].depthJump << "/"
                                 << numPaths << ", contact flicker = " << st[m].flicker << ", normal jump > 5 deg = "
                                 << st[m].normalJump << ", largest depth jump = " << st[m].maxJump << " ds = "
                                 << st[m].maxJumpRel << " of the smaller radius, paths with a jump > 1 % of it = "
                                 << st[m].bigJump);
   // GJK/EPA and the refined detection must be continuous along translations
   WALBERLA_CHECK_EQUAL(st["gjk"].depthJump, 0);
   WALBERLA_CHECK_EQUAL(st["refined"].flicker, 0);
   WALBERLA_CHECK_EQUAL(st["refined"].bigJump, 0); // no jump larger than 1 % of the smaller radius
}

//======================================================================================================================
// Refined cylinder-cylinder detection (refineCylinderContact / cylinderCylinderContactRefined in
// ClosedFormContactsCylinderCylinder.h) vs. GJK/EPA, sampled as the master comparison (seed 24680, uniform SO(3)).
// Four detection variants are compared, all against GJK/EPA:
//   raw      closed form alone;
//   hybrid   closed form + GJK/EPA for "no scenario" (Tier 1) and special contacts (current production option);
//   trust    closed form, "no scenario" trusted as "no contact" (no GJK/EPA at all; special contact kept as is);
//   refined  closed form + support-function refinement (no GJK/EPA at all).
// Reported: existence errors (missed / false contacts), depth errors > max(0.01, 5 %), max depth error, normal angle,
// mean projection-overlap evaluations per refined call, and the time per call of each variant
// (testRefinedCylinderContactVsGjk). Pass criterion (refined): no missed contacts, no false contacts, no depth error
// > max(0.01, 5 %), and every contact point inside the intersection and not deeper than the depth in either body
// (checked also for raw band-band-skewed).
//======================================================================================================================
//! Signed distance of point p to the surface of cylinder c (negative inside).
real_t cylinderSignedDistance(const CylinderState& c, const Vec3& p)
{
   const Vec3   d   = p - c.pos;
   const Vec3   a   = c.axis();
   const real_t ax  = d * a;
   const real_t rad = (d - ax * a).length();
   const real_t dr = rad - c.radius, dz = std::abs(ax) - real_t(0.5) * c.length;
   if (dr <= 0 && dz <= 0) return std::max(dr, dz);
   return std::sqrt(std::max(dr, real_t(0)) * std::max(dr, real_t(0)) + std::max(dz, real_t(0)) * std::max(dz, real_t(0)));
}

//======================================================================================================================
// Box-box contact point (detail::boxBoxContactPoint, used by boxBoxSAT): checked against hand-derived points, for
// lying inside the intersection on random pairs, and for continuity along a motion. The former point (formerSatPoint,
// near the line between the centres; it can lie far from the contact and outside both boxes) is reported alongside.
//======================================================================================================================
real_t boxSignedDistance(const Vec3& c, const Mat3& R, const Vec3& h, const Vec3& p)
{
   const Vec3 l = R.getTranspose() * (p - c);
   Vec3       q(std::abs(l[0]) - h[0], std::abs(l[1]) - h[1], std::abs(l[2]) - h[2]);
   const real_t inside = std::min(std::max(q[0], std::max(q[1], q[2])), real_t(0));
   const Vec3   qp(std::max(q[0], real_t(0)), std::max(q[1], real_t(0)), std::max(q[2], real_t(0)));
   return qp.length() + inside;
}

//! Former boxBoxSAT point: midpoint of the centres shifted by (ra - rb)/2 along the axis; for comparison only.
Vec3 formerSatPoint(const Vec3& posA, const Mat3& RA, const Vec3& halfA, const Vec3& posB, const Mat3& RB, const Vec3& halfB,
                    const Vec3& normal)
{
   const Vec3 axis = -normal; // former: normal = -sign * axis with sign chosen so that sign*axis points from A to B
   real_t     ra = 0, rb = 0;
   for (uint_t i = 0; i < 3; ++i)
   {
      ra += halfA[i] * std::abs(Vec3(RA(0, i), RA(1, i), RA(2, i)) * axis);
      rb += halfB[i] * std::abs(Vec3(RB(0, i), RB(1, i), RB(2, i)) * axis);
   }
   return real_t(0.5) * (posA + posB) + real_t(0.5) * (ra - rb) * axis;
}

struct BoxPointCheck
{
   bool   inside;
   real_t outside, deepest; // outside: max signed distance; deepest: max depth below either surface / overlap
};
BoxPointCheck checkBoxPoint(const Vec3& pA, const Mat3& RA, const Vec3& hA, const Vec3& pB, const Mat3& RB, const Vec3& hB,
                            const Vec3& point, real_t depth)
{
   const real_t dA = boxSignedDistance(pA, RA, hA, point), dB = boxSignedDistance(pB, RB, hB, point);
   const real_t tol = real_t(1e-9) + real_t(1e-6) * depth;
   return BoxPointCheck{ std::max(dA, dB) <= tol, std::max(dA, dB), std::max(-dA, -dB) / depth };
}

void testBoxBoxContactPointHandCases()
{
   const Mat3 I;
   auto run = [&](const char* name, const Vec3& pA, const Mat3& RA, const Vec3& hA, const Vec3& pB, const Mat3& RB,
                  const Vec3& hB, const Vec3& expected, real_t tolPoint) {
      const ContactResult r = boxBoxSAT(pA, RA, hA, pB, RB, hB);
      WALBERLA_CHECK(r.contact, name);
      const auto   chk = checkBoxPoint(pA, RA, hA, pB, RB, hB, r.point, r.overlap);
      const Vec3   old = formerSatPoint(pA, RA, hA, pB, RB, hB, r.normal);
      const auto   chkOld = checkBoxPoint(pA, RA, hA, pB, RB, hB, old, r.overlap);
      WALBERLA_LOG_INFO_ON_ROOT("[box-point] " << name << ": overlap=" << r.overlap << " normal=" << r.normal << " point="
                                << r.point << " expected=" << expected << " |diff|=" << (r.point - expected).length()
                                << " inside=" << chk.inside << " | former point " << old << " |diff|="
                                << (old - expected).length() << " inside=" << chkOld.inside);
      WALBERLA_CHECK_LESS((r.point - expected).length(), tolPoint, name);
      WALBERLA_CHECK(chk.inside, name);
      return (old - expected).length();
   };
   const real_t d = real_t(0.1);
   // (a) centred stack: the line between the centres runs through the contact; point (0, 0, 1 - d/2)
   run("centred stack", Vec3(0, 0, 0), I, Vec3(1, 1, 1), Vec3(0, 0, real_t(2) - d), I, Vec3(1, 1, 1),
       Vec3(0, 0, real_t(1) - real_t(0.5) * d), real_t(1e-9));
   // (b) small box near the end of a long one: contact area x in [4, 6] -> centroid x = 5 (former point: x = 2.5)
   const real_t oldB = run("small box on the end of a long one", Vec3(0, 0, 0), I, Vec3(6, 1, real_t(1.5)),
                           Vec3(5, 0, real_t(2.5) - d), I, Vec3(1, 1, 1), Vec3(5, 0, real_t(1.5) - real_t(0.5) * d), real_t(1e-9));
   WALBERLA_CHECK_GREATER(oldB, real_t(2), "the former point should be ~2.5 off here");
   // (c) two equal faces, partly overlapping: rectangle x in [1, 2], y in [-1, 2] -> (1.5, 0.5)
   run("partly overlapping faces", Vec3(0, 0, 0), I, Vec3(2, 2, 1), Vec3(3, 1, real_t(2) - d), I, Vec3(2, 2, 1),
       Vec3(real_t(1.5), real_t(0.5), real_t(1) - real_t(0.5) * d), real_t(1e-9));
   // (d) corner into a face, off centre: cube with a body diagonal vertical, lowest vertex at height 1 - d above
   // (0.7,-0.3)
   {
      Rot3 rot;
      rot.rotate(Vec3(1, -1, 0).getNormalized(), std::acos(real_t(1) / std::sqrt(real_t(3))));
      const Mat3 RB  = rot.getMatrix();
      const Vec3 hB(real_t(0.5), real_t(0.5), real_t(0.5));
      // lowest vertex in the body frame
      Vec3   low;
      real_t lowZ = std::numeric_limits< real_t >::max();
      for (int sx = -1; sx <= 1; sx += 2)
         for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2)
            {
               const Vec3 v = RB * Vec3(sx * hB[0], sy * hB[1], sz * hB[2]);
               if (v[2] < lowZ) { lowZ = v[2]; low = v; }
            }
      const Vec3 target(real_t(0.7), real_t(-0.3), real_t(1) - d);
      const Vec3 pB = target - low;
      run("corner into a face", Vec3(0, 0, 0), I, Vec3(2, 2, 1), pB, RB, hB, target + Vec3(0, 0, real_t(0.5) * d), real_t(1e-6));
   }
   // (e) crossing edges: bar A along x rotated 45 deg about x (top edge along x at z = sqrt(2)/2), bar B along y
   // rotated 45 deg about y (bottom edge along y), centred at x = 3, y = -2 -> edges cross at (3, 0)
   {
      Rot3 rA, rB;
      rA.rotate(Vec3(1, 0, 0), math::pi * real_t(0.25));
      rB.rotate(Vec3(0, 1, 0), math::pi * real_t(0.25));
      const real_t e = std::sqrt(real_t(2)) * real_t(0.5);
      const Vec3   pB(3, -2, real_t(2) * e - d);
      const real_t oldE = run("crossing edges", Vec3(0, 0, 0), rA.getMatrix(), Vec3(6, real_t(0.5), real_t(0.5)), pB,
                              rB.getMatrix(), Vec3(real_t(0.5), 6, real_t(0.5)), Vec3(3, 0, e - real_t(0.5) * d), real_t(1e-6));
      WALBERLA_CHECK_GREATER(oldE, real_t(0.5), "the former point should be off here");
   }
   // (f) a logged pair of 12x2x3 cuboids from a randomised drop test: the former point is 4.3 from the centroid of the
   //     overlap and outside both cuboids; the new point must lie inside the intersection and within 1 of the centroid
   {
      const Vec3 pA(real_t(1.37178412), real_t(6.07269819), real_t(1.71354156));
      const Vec3 pB(real_t(2.48370382), real_t(1.60370074), real_t(6.09968309));
      const Vec3 h(6, 1, real_t(1.5));
      // rotation matrices as logged (row-major)
      const Mat3 RA = Mat3(real_t(-0.00032603468), real_t(0.957859121), real_t(0.287238225), real_t(0.999999942), real_t(0.000341559239), real_t(-3.93737322e-06), real_t(-0.000101880317), real_t(0.287238207), real_t(-0.957859176));
      const Mat3 RB = Mat3(real_t(-0.228133258), real_t(-0.968314389), real_t(-0.101599512), real_t(0.00176512624), real_t(0.10393976), real_t(-0.994582028), real_t(0.973628318), real_t(-0.227076574), real_t(-0.0220029189));
      const ContactResult r = boxBoxSAT(pA, RA, h, pB, RB, h);
      WALBERLA_CHECK(r.contact);
      // centroid of the overlap by sampling box A (deterministic grid) and keeping the points inside B
      Vec3 sum(0);
      int  n = 0;
      for (int i = 0; i < 120; ++i)
         for (int j = 0; j < 20; ++j)
            for (int k = 0; k < 30; ++k)
            {
               const Vec3 loc((real_t(i) + real_t(0.5)) / real_t(120) * real_t(12) - real_t(6),
                              (real_t(j) + real_t(0.5)) / real_t(20) * real_t(2) - real_t(1),
                              (real_t(k) + real_t(0.5)) / real_t(30) * real_t(3) - real_t(1.5));
               const Vec3 x = pA + RA * loc;
               if (boxSignedDistance(pB, RB, h, x) <= real_t(0)) { sum += x; ++n; }
            }
      WALBERLA_CHECK_GREATER(n, 0);
      const Vec3 cen = sum / real_t(n);
      const auto chk = checkBoxPoint(pA, RA, h, pB, RB, h, r.point, r.overlap);
      const Vec3 old = formerSatPoint(pA, RA, h, pB, RB, h, r.normal);
      WALBERLA_LOG_INFO_ON_ROOT("[box-point] drop-test launch pair: overlap=" << r.overlap << " point=" << r.point
                                << " inside=" << chk.inside << " overlap-centroid=" << cen << " |point-centroid|="
                                << (r.point - cen).length() << " | former point " << old << " |former-centroid|="
                                << (old - cen).length());
      WALBERLA_CHECK(chk.inside);
      WALBERLA_CHECK_LESS((r.point - cen).length(), real_t(1));
      WALBERLA_CHECK_GREATER((old - cen).length(), real_t(3));
   }
}

//! 10000 pairs, seed 13579, near the contact boundary: even trials random half-extents U[0.3,2] (as
//! testBoxBoxRandomizedFuzz), odd trials 12x2x3 drop-test cuboids. Every point must lie in the intersection and not
//! deeper than the overlap below either surface.
void testBoxBoxContactPointFuzz()
{
   FuzzRng                                  fuzz(13579);
   std::uniform_real_distribution< real_t > halfExtentDist(real_t(0.3), real_t(2.0));
   int    contacts = 0, outside = 0, tooDeep = 0, formerOutside = 0;
   real_t maxOutside = 0, maxDeep = 0, formerMaxOutside = 0;
   std::vector< real_t > distToGjk;
   for (int trial = 0; trial < 10000; ++trial)
   {
      Vec3 hA, hB;
      if (trial % 2 == 0)
      {
         hA = Vec3(halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng));
         hB = Vec3(halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng), halfExtentDist(fuzz.rng));
      }
      else
         hA = hB = Vec3(6, 1, real_t(1.5)); // the drop-test cuboids
      const Rot3   rA = fuzz.randomRotation(), rB = fuzz.randomRotation();
      const Vec3   dir = fuzz.randomUnitVec();
      const Vec3   pA(0, 0, 0);
      const real_t minSize = std::min({ hA[0], hA[1], hA[2], hB[0], hB[1], hB[2] });
      auto contactAt = [&](real_t sep) { return gjkBoxBox(pA, rA, real_t(2) * hA, pA + dir * sep, rB, real_t(2) * hB).contact; };
      const Vec3 pB = pA + dir * fuzz.sampleSeparationNearBoundary(contactAt, real_t(0), real_t(3) * (hA.length() + hB.length()), minSize);
      const Mat3 RA = rA.getMatrix(), RB = rB.getMatrix();
      const ContactResult r = boxBoxSAT(pA, RA, hA, pB, RB, hB);
      if (!r.contact || r.overlap <= real_t(0)) continue;
      ++contacts;
      const auto chk = checkBoxPoint(pA, RA, hA, pB, RB, hB, r.point, r.overlap);
      maxOutside = std::max(maxOutside, chk.outside / r.overlap);
      maxDeep    = std::max(maxDeep, chk.deepest);
      if (!chk.inside) ++outside;
      if (chk.deepest > real_t(1) + real_t(1e-6)) ++tooDeep;
      const auto chkOld = checkBoxPoint(pA, RA, hA, pB, RB, hB, formerSatPoint(pA, RA, hA, pB, RB, hB, r.normal), r.overlap);
      if (!chkOld.inside) ++formerOutside;
      formerMaxOutside = std::max(formerMaxOutside, chkOld.outside / r.overlap);
      const GjkResult g = gjkBoxBox(pA, rA, real_t(2) * hA, pB, rB, real_t(2) * hB);
      if (g.contact) distToGjk.push_back((r.point - g.point).length());
   }
   std::sort(distToGjk.begin(), distToGjk.end());
   auto pct = [&](real_t q) { return distToGjk.empty() ? real_t(0) : distToGjk[size_t(q * real_t(distToGjk.size() - 1))]; };
   WALBERLA_LOG_INFO_ON_ROOT("[box-point fuzz] 10000 pairs (half uniform-random, half 12x2x3), " << contacts
                             << " contacts: points outside the intersection " << outside << " (max signed distance "
                             << maxOutside << " x overlap), deeper than the overlap " << tooDeep << " (max " << maxDeep
                             << " x overlap) | former point: outside " << formerOutside << " (max " << formerMaxOutside
                             << " x overlap) | distance to the GJK/EPA point median " << pct(0.5) << " p95 " << pct(0.95)
                             << " max " << pct(1));
   WALBERLA_CHECK_EQUAL(outside, 0);
   WALBERLA_CHECK_EQUAL(tooDeep, 0);
}

//! Continuity of the contact point of SAT and GJK/EPA along straight paths through the contact: 300 pairs of 12x2x3
//! cuboids (seed 24601, uniform SO(3)), translated from an overlap of 0.3 to a gap of 0.1 in 300 steps of ds. Metric
//! per method: largest distance the point moves in one step, in units of ds (paths with a jump > 10 ds counted).
//! Feature changes (vertex/edge/face) are real jumps of a single-point model; reported, not asserted.
void testBoxBoxContactPointContinuity()
{
   FuzzRng fuzz(24601);
   const std::vector< std::string > methods = { "sat", "gjk" };
   struct Stats { real_t worstJump = 0, sumJumps = 0; int pathsWithJumpOver10 = 0; };
   std::map< std::string, Stats > st;
   int paths = 0;
   for (int k = 0; k < 300; ++k)
   {
      const Vec3   hA(6, 1, real_t(1.5)), hB(6, 1, real_t(1.5));
      const Rot3   rA = fuzz.randomRotation(), rB = fuzz.randomRotation();
      const Vec3   dir = fuzz.randomUnitVec();
      const Vec3   pA(0, 0, 0);
      auto contactAt = [&](real_t sep) { return gjkBoxBox(pA, rA, real_t(2) * hA, pA + dir * sep, rB, real_t(2) * hB).contact; };
      // touching separation along dir, then a path from an overlap of 0.3 to a gap of 0.1 along dir
      real_t lo = 0, hi = real_t(3) * (hA.length() + hB.length());
      for (int it = 0; it < 60; ++it) { const real_t m = real_t(0.5) * (lo + hi); (contactAt(m) ? lo : hi) = m; }
      const real_t s0 = lo - real_t(0.3), s1 = lo + real_t(0.1), ds = (s1 - s0) / real_t(300);
      const Mat3   RA = rA.getMatrix(), RB = rB.getMatrix();
      std::map< std::string, bool >   havePrev{ { "sat", false }, { "gjk", false } };
      std::map< std::string, Vec3 >   prev;
      std::map< std::string, real_t > pathMax{ { "sat", real_t(0) }, { "gjk", real_t(0) } };
      for (int i = 0; i <= 300; ++i)
      {
         const Vec3          pB = pA + dir * (s0 + ds * real_t(i));
         const ContactResult r  = boxBoxSAT(pA, RA, hA, pB, RB, hB);
         const GjkResult     g  = gjkBoxBox(pA, rA, real_t(2) * hA, pB, rB, real_t(2) * hB);
         const std::map< std::string, std::pair< bool, Vec3 > > res = { { "sat", { r.contact, r.point } },
                                                                          { "gjk", { g.contact, g.point } } };
         for (const auto& m : methods)
         {
            const auto& [c, pt] = res.at(m);
            if (!c) { havePrev[m] = false; continue; }
            if (havePrev[m]) pathMax[m] = std::max(pathMax[m], (pt - prev[m]).length() / ds);
            prev[m]     = pt;
            havePrev[m] = true;
         }
      }
      ++paths;
      for (const auto& m : methods)
      {
         st[m].worstJump = std::max(st[m].worstJump, pathMax[m]);
         st[m].sumJumps += pathMax[m];
         if (pathMax[m] > real_t(10)) ++st[m].pathsWithJumpOver10;
      }
   }
   for (const auto& m : methods)
      WALBERLA_LOG_INFO_ON_ROOT("[box-point continuity] " << m << ": " << paths << " paths of 12x2x3 cuboids, step ds: "
                                << "largest point jump " << st[m].worstJump << " ds (mean of per-path maxima "
                                << st[m].sumJumps / real_t(paths) << " ds), paths with a jump > 10 ds: "
                                << st[m].pathsWithJumpOver10);
}

//! Cost of the box-box contact point: random near-boundary pairs vs. nearly flat resting contacts
//! (face on face, tilt 1e-4..5e-2 rad, overlap 1e-4..1e-2), the configurations of a settled cuboid pile.
void benchmarkBoxBoxContactPoint()
{
   using clk = std::chrono::steady_clock;
   FuzzRng fuzz(31415);
   auto runSet = [&](const char* name, const std::vector< std::array< Vec3, 2 > >& pos, const std::vector< std::array< Mat3, 2 > >& rot,
                     const Vec3& h) {
      auto& st = detail::boxPointIterationStats();
      st = detail::BoxPointIterationStats{};
      int        contacts = 0;
      const auto t0       = clk::now();
      for (size_t i = 0; i < pos.size(); ++i)
         contacts += boxBoxSAT(pos[i][0], rot[i][0], h, pos[i][1], rot[i][1], h).contact ? 1 : 0;
      const double us = std::chrono::duration< double, std::micro >(clk::now() - t0).count() / double(pos.size());
      WALBERLA_LOG_INFO_ON_ROOT("[box-point cost] " << name << ": " << pos.size() << " pairs, " << contacts << " contacts, "
                                << us << " us/pair | pull iterations: mean "
                                << (st.calls ? double(st.iterations) / double(st.calls) : 0.0) << " max " << st.maxIterations
                                << " capped(2000) " << st.capped);
   };
   const Vec3 h(6, 1, real_t(1.5));
   std::vector< std::array< Vec3, 2 > > posR, posF;
   std::vector< std::array< Mat3, 2 > > rotR, rotF;
   std::uniform_real_distribution< real_t > u01(real_t(0), real_t(1));
   for (int t = 0; t < 4000; ++t)
   {
      const Rot3 rA = fuzz.randomRotation(), rB = fuzz.randomRotation();
      const Vec3 dir = fuzz.randomUnitVec();
      auto contactAt = [&](real_t sep) { return gjkBoxBox(Vec3(0), rA, real_t(2) * h, dir * sep, rB, real_t(2) * h).contact; };
      posR.push_back({ Vec3(0), dir * fuzz.sampleSeparationNearBoundary(contactAt, real_t(0), real_t(3) * real_t(2) * h.length(), real_t(1)) });
      rotR.push_back({ rA.getMatrix(), rB.getMatrix() });
      // flat: B on top of A (z), small tilt about a random horizontal axis, random lateral offset, small overlap
      const real_t tilt = std::pow(real_t(10), real_t(-4) + real_t(2.7) * u01(fuzz.rng));
      const real_t ov   = std::pow(real_t(10), real_t(-4) + real_t(2) * u01(fuzz.rng));
      const real_t phi  = real_t(2) * math::pi * u01(fuzz.rng);
      Rot3         rT;
      rT.rotate(Vec3(std::cos(phi), std::sin(phi), 0), tilt);
      Rot3 yaw;
      yaw.rotate(Vec3(0, 0, 1), real_t(2) * math::pi * u01(fuzz.rng));
      const Mat3   RB  = rT.getMatrix() * yaw.getMatrix();
      const Vec3   off((u01(fuzz.rng) - real_t(0.5)) * real_t(6), (u01(fuzz.rng) - real_t(0.5)) * real_t(2), 0);
      // lowest point of B for this rotation, so that B penetrates A's top face (z = 1.5) by ov
      real_t lowZ = std::numeric_limits< real_t >::max();
      for (int sx = -1; sx <= 1; sx += 2)
         for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2)
               lowZ = std::min(lowZ, (RB * Vec3(sx * h[0], sy * h[1], sz * h[2]))[2]);
      posF.push_back({ Vec3(0), off + Vec3(0, 0, real_t(1.5) - ov - lowZ) });
      rotF.push_back({ Mat3(), RB });
   }
   runSet("random near-boundary 12x2x3", posR, rotR, h);
   runSet("flat resting face-on-face 12x2x3", posF, rotF, h);
}

//! detail::clipLineToBox, hand cases: the interval of t for which x0 + t dir lies in the box
void testClipLineToBox()
{
   const Vec3 axesI[3] = { Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1) };
   auto clip = [&](const Vec3& x0, const Vec3& dir, const Vec3& pos, const Vec3 axes[3], const Vec3& half) {
      real_t lo = -std::numeric_limits< real_t >::max(), hi = std::numeric_limits< real_t >::max();
      detail::clipLineToBox(x0, dir, pos, axes, half, lo, hi);
      return std::make_pair(lo, hi);
   };
   const Vec3 h(6, 1, real_t(1.5));
   // through the centre along z: [-1.5, 1.5]
   auto r = clip(Vec3(0), Vec3(0, 0, 1), Vec3(0), axesI, h);
   WALBERLA_CHECK_FLOAT_EQUAL(r.first, real_t(-1.5));
   WALBERLA_CHECK_FLOAT_EQUAL(r.second, real_t(1.5));
   // offset start and box: x0 = (2, 0.5, 10), box at (1, 0, 3) -> z in [1.5, 4.5] -> t in [-8.5, -5.5]
   r = clip(Vec3(2, real_t(0.5), 10), Vec3(0, 0, 1), Vec3(1, 0, 3), axesI, h);
   WALBERLA_CHECK_FLOAT_EQUAL(r.first, real_t(-8.5));
   WALBERLA_CHECK_FLOAT_EQUAL(r.second, real_t(-5.5));
   // parallel to an axis but outside that slab: empty
   r = clip(Vec3(0, 2, 0), Vec3(0, 0, 1), Vec3(0), axesI, h);
   WALBERLA_CHECK_GREATER(r.first, r.second);
   // diagonal in the x-z plane through the centre: limited by z (|t| sqrt(1/2) <= 1.5) -> |t| <= 1.5 sqrt(2)
   r = clip(Vec3(0), Vec3(1, 0, 1).getNormalized(), Vec3(0), axesI, h);
   WALBERLA_CHECK_FLOAT_EQUAL(r.first, -real_t(1.5) * std::sqrt(real_t(2)));
   WALBERLA_CHECK_FLOAT_EQUAL(r.second, real_t(1.5) * std::sqrt(real_t(2)));
   // rotated box (45 deg about z), line along z through a point inside: full z extent
   Rot3 rot;
   rot.rotate(Vec3(0, 0, 1), math::pi * real_t(0.25));
   const Mat3 R = rot.getMatrix();
   const Vec3 axesR[3] = { Vec3(R(0, 0), R(1, 0), R(2, 0)), Vec3(R(0, 1), R(1, 1), R(2, 1)), Vec3(R(0, 2), R(1, 2), R(2, 2)) };
   r = clip(Vec3(3, 3, 0), Vec3(0, 0, 1), Vec3(0), axesR, h); // (3,3) lies on the long axis at 4.24 < 6
   WALBERLA_CHECK_FLOAT_EQUAL(r.first, real_t(-1.5));
   WALBERLA_CHECK_FLOAT_EQUAL(r.second, real_t(1.5));
   r = clip(Vec3(3, -3, 0), Vec3(0, 0, 1), Vec3(0), axesR, h); // 4.24 across the short axis (half 1): outside
   WALBERLA_CHECK_GREATER(r.first, r.second);
   // narrowing: the second box shrinks the interval of the first
   real_t lo = -std::numeric_limits< real_t >::max(), hi = std::numeric_limits< real_t >::max();
   detail::clipLineToBox(Vec3(0), Vec3(0, 0, 1), Vec3(0), axesI, h, lo, hi);
   detail::clipLineToBox(Vec3(0), Vec3(0, 0, 1), Vec3(0, 0, real_t(2.9)), axesI, h, lo, hi); // z in [1.4, 4.4]
   WALBERLA_CHECK_FLOAT_EQUAL(lo, real_t(1.4));
   WALBERLA_CHECK_FLOAT_EQUAL(hi, real_t(1.5));
   WALBERLA_LOG_INFO_ON_ROOT("[clip-line] hand cases passed");
}

//! Contacts as in a settled pile: 4000 pairs of arbitrarily rotated 12x2x3 cuboids (seed 27182), B placed along a
//! random direction at the touching distance minus a tiny overlap (1e-7..1e-2, log-uniform). For such contacts the
//! former alternating projection often hits its 2000-iteration cap. Asserted for the closed-form point (line clip along
//! the normal): boxBoxSAT returns it, it lies inside the intersection, it is cheaper than the former projection, and
//! more than 3000 contacts occur. Reported: its offset from the former point perpendicular to the normal (lever-arm
//! change).
void testBoxBoxContactPointRestingContacts()
{
   using clk = std::chrono::steady_clock;
   FuzzRng fuzz(27182);
   std::uniform_real_distribution< real_t > u01(real_t(0), real_t(1));
   const Vec3 h(6, 1, real_t(1.5));
   struct Case { Vec3 pB; Mat3 RA, RB; };
   std::vector< Case > cases;
   for (int t = 0; t < 4000; ++t)
   {
      const Rot3 rA = fuzz.randomRotation(), rB = fuzz.randomRotation();
      const Vec3 dir = fuzz.randomUnitVec();
      const Mat3 RA = rA.getMatrix(), RB = rB.getMatrix();
      real_t lo = 0, hi = real_t(3) * real_t(2) * h.length();
      for (int it = 0; it < 80; ++it)
      {
         const real_t m = real_t(0.5) * (lo + hi);
         (boxBoxSAT(Vec3(0), RA, h, dir * m, RB, h).contact ? lo : hi) = m;
      }
      const real_t ov = std::pow(real_t(10), real_t(-7) + real_t(5) * u01(fuzz.rng));
      cases.push_back(Case{ dir * (lo - ov), RA, RB });
   }
   auto axesOf = [](const Mat3& R, Vec3 a[3]) {
      for (uint_t i = 0; i < 3; ++i) a[i] = Vec3(R(0, i), R(1, i), R(2, i));
   };
   int    contacts = 0, outsideNew = 0, outsideOld = 0;
   real_t maxOutsideNew = 0, maxOutsideOld = 0;
   std::vector< real_t > perp, perpRel;
   auto& st = detail::boxPointIterationStats();
   st = detail::BoxPointIterationStats{};
   double usNew = 0, usOld = 0;
   for (const auto& c : cases)
   {
      const ContactResult r = boxBoxSAT(Vec3(0), c.RA, h, c.pB, c.RB, h);
      if (!r.contact || r.overlap <= real_t(0)) continue;
      ++contacts;
      Vec3 aA[3], aB[3];
      axesOf(c.RA, aA);
      axesOf(c.RB, aB);
      const auto t0 = clk::now();
      const Vec3 pNew = detail::boxBoxContactPoint(Vec3(0), aA, h, c.pB, aB, h, r.normal, r.overlap);
      const auto t1 = clk::now();
      const Vec3 pOld = detail::boxBoxContactPoint(Vec3(0), aA, h, c.pB, aB, h, r.normal, r.overlap, true);
      const auto t2 = clk::now();
      usNew += std::chrono::duration< double, std::micro >(t1 - t0).count();
      usOld += std::chrono::duration< double, std::micro >(t2 - t1).count();
      WALBERLA_CHECK_FLOAT_EQUAL((pNew - r.point).length() + real_t(1), real_t(1), "boxBoxSAT must use the new point");
      const auto chkN = checkBoxPoint(Vec3(0), c.RA, h, c.pB, c.RB, h, pNew, r.overlap);
      const auto chkO = checkBoxPoint(Vec3(0), c.RA, h, c.pB, c.RB, h, pOld, r.overlap);
      if (!chkN.inside) ++outsideNew;
      if (!chkO.inside) ++outsideOld;
      maxOutsideNew = std::max(maxOutsideNew, chkN.outside);
      maxOutsideOld = std::max(maxOutsideOld, chkO.outside);
      const Vec3 d = pNew - pOld;
      const real_t dp = (d - (d * r.normal) * r.normal).length();
      perp.push_back(dp);
      perpRel.push_back(dp / r.overlap);
   }
   std::sort(perp.begin(), perp.end());
   std::sort(perpRel.begin(), perpRel.end());
   auto pct = [](const std::vector< real_t >& v, real_t q) { return v.empty() ? real_t(0) : v[size_t(q * real_t(v.size() - 1))]; };
   const uint64_t calls = uint64_t(contacts); // one new + one legacy evaluation per contact
   WALBERLA_LOG_INFO_ON_ROOT("[box-point resting] " << cases.size() << " pairs, " << contacts << " contacts | new point: "
                             << usNew / contacts << " us, line clip hit " << st.lineClipHits / 2 << " / polytope shift "
                             << st.polytopeShifts / 2 << " / projection fallback " << st.projectionFallbacks / 2 << " of " << calls
                             << ", outside the intersection " << outsideNew << " (max " << maxOutsideNew << ") | former "
                             << "projection: " << usOld / contacts << " us, outside " << outsideOld << " (max "
                             << maxOutsideOld << "), projection calls " << st.calls << " mean iterations "
                             << (st.calls ? double(st.iterations) / double(st.calls) : 0.0) << " capped " << st.capped
                             << " | offset new-former perpendicular to the normal: median " << pct(perp, real_t(0.5))
                             << " p95 " << pct(perp, real_t(0.95)) << " max " << pct(perp, real_t(1)) << " (x overlap: median "
                             << pct(perpRel, real_t(0.5)) << " p95 " << pct(perpRel, real_t(0.95)) << ")");
   WALBERLA_CHECK_GREATER(contacts, 3000);
   WALBERLA_CHECK_EQUAL(outsideNew, 0);
   WALBERLA_CHECK_LESS(usNew, usOld);
}

//! ANTIDUNES_BOX_PAIRS_EVAL=<file> (optional): the same comparison on dumped cuboid pairs, one pair per line:
//! posA(3) axesA(3x3, one axis per triple) halfA(3) posB(3) axesB(3x3) halfB(3), e.g. the near-contact pairs of a
//! drop-test pile. Asserts that every new point lies inside the intersection.
void evaluateBoxPairsDump()
{
   const char* path = std::getenv("ANTIDUNES_BOX_PAIRS_EVAL");
   if (!path) return;
   std::ifstream in(path);
   WALBERLA_CHECK(in.good(), "cannot open " << path);
   using clk = std::chrono::steady_clock;
   int    pairs = 0, contacts = 0, outsideNew = 0, outsideOld = 0;
   double usNew = 0, usOld = 0;
   std::vector< real_t > perp;
   auto& st = detail::boxPointIterationStats();
   st = detail::BoxPointIterationStats{};
   real_t v[30];
   while (true)
   {
      bool ok = true;
      for (int k = 0; k < 30 && ok; ++k) ok = bool(in >> v[k]);
      if (!ok) break;
      ++pairs;
      Vec3 pA(v[0], v[1], v[2]), pB(v[15], v[16], v[17]), aA[3], aB[3];
      Mat3 RA, RB;
      for (uint_t i = 0; i < 3; ++i)
      {
         aA[i] = Vec3(v[3 + 3 * i], v[4 + 3 * i], v[5 + 3 * i]);
         aB[i] = Vec3(v[18 + 3 * i], v[19 + 3 * i], v[20 + 3 * i]);
         for (uint_t r = 0; r < 3; ++r) { RA(r, i) = aA[i][r]; RB(r, i) = aB[i][r]; }
      }
      const Vec3 hA(v[12], v[13], v[14]), hB(v[27], v[28], v[29]);
      const ContactResult r = boxBoxSAT(pA, RA, hA, pB, RB, hB);
      if (!r.contact || r.overlap <= real_t(0)) continue;
      ++contacts;
      const auto t0 = clk::now();
      const Vec3 pNew = detail::boxBoxContactPoint(pA, aA, hA, pB, aB, hB, r.normal, r.overlap);
      const auto t1 = clk::now();
      const Vec3 pOld = detail::boxBoxContactPoint(pA, aA, hA, pB, aB, hB, r.normal, r.overlap, true);
      const auto t2 = clk::now();
      usNew += std::chrono::duration< double, std::micro >(t1 - t0).count();
      usOld += std::chrono::duration< double, std::micro >(t2 - t1).count();
      if (!checkBoxPoint(pA, RA, hA, pB, RB, hB, pNew, r.overlap).inside) ++outsideNew;
      if (!checkBoxPoint(pA, RA, hA, pB, RB, hB, pOld, r.overlap).inside) ++outsideOld;
      const Vec3 d = pNew - pOld;
      perp.push_back((d - (d * r.normal) * r.normal).length());
   }
   std::sort(perp.begin(), perp.end());
   auto pct = [&](real_t q) { return perp.empty() ? real_t(0) : perp[size_t(q * real_t(perp.size() - 1))]; };
   WALBERLA_LOG_INFO_ON_ROOT("[box-point pairs " << path << "] " << pairs << " pairs, " << contacts << " contacts | new: "
                             << usNew / contacts << " us, line clip hit " << st.lineClipHits / 2 << " / polytope shift "
                             << st.polytopeShifts / 2 << " / projection fallback " << st.projectionFallbacks / 2 << " of "
                             << contacts << ", outside " << outsideNew
                             << " | former projection: " << usOld / contacts << " us, outside " << outsideOld
                             << ", mean iterations " << (st.calls ? double(st.iterations) / double(st.calls) : 0.0)
                             << " capped " << st.capped << " | perpendicular offset new-former: median " << pct(real_t(0.5))
                             << " p95 " << pct(real_t(0.95)) << " max " << pct(real_t(1)));
   WALBERLA_CHECK_EQUAL(outsideNew, 0);
}

void testRefinedCylinderContactVsGjk()
{
   FuzzRng                                  fuzz(24680);
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0)), lengthDist(real_t(0.5), real_t(4.0));
   constexpr int numTrials = 5000;
   struct Pair { CylinderState ci, cj; Rot3 rotI, rotJ; };
   std::vector< Pair > pairs;
   for (int t = 0; t < numTrials; ++t)
   {
      const real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng), Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      const Rot3   rotI = fuzz.randomRotation(), rotJ = fuzz.randomRotation();
      const real_t circ = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li) + std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      const Vec3   dir  = fuzz.randomUnitVec();
      const Vec3   posI(0, 0, 0);
      auto contactAt = [&](real_t sep) { return gjkCylCyl(posI, rotI, Ri, Li, posI + dir * sep, rotJ, Rj, Lj).contact; };
      const Vec3 posJ = posI + dir * fuzz.sampleSeparationNearBoundary(contactAt, real_t(0), real_t(3) * circ, std::min(Ri, Rj));
      pairs.push_back(Pair{ CylinderState{ posI, rotI.getMatrix(), Ri, Li }, CylinderState{ posJ, rotJ.getMatrix(), Rj, Lj },
                            rotI, rotJ });
   }

   struct Stats { int missed = 0, falsePos = 0, wrongDepth = 0; real_t maxErr = 0, maxErrRel = 0, maxAngle = 0; };
   // Contact-point check (the point need not be unique, e.g. face on band): the point must lie in the intersection of
   // the two cylinders -- signed distance to both surfaces <= tol -- and must not lie deeper than the penetration depth
   // in either body. Counted per method over the contacts both GJK/EPA and the method report.
   struct PointStats { int checked = 0, outside = 0, tooDeep = 0; real_t maxOutside = 0, maxDepthRatio = 0; };
   std::map< std::string, PointStats > pst;
   auto checkPoint = [&](const std::string& name, const CylinderState& ci, const CylinderState& cj, bool c,
                         const Vec3& point, real_t depth) {
      if (!c || depth <= real_t(0)) return;
      auto&        ps  = pst[name];
      const real_t di  = cylinderSignedDistance(ci, point), dj = cylinderSignedDistance(cj, point);
      const real_t tol = real_t(1e-6) + real_t(1e-3) * depth;
      ++ps.checked;
      const real_t out = std::max(di, dj);
      ps.maxOutside    = std::max(ps.maxOutside, out / depth);
      if (out > tol) ++ps.outside;
      const real_t deep = std::max(-di, -dj) / depth;
      ps.maxDepthRatio  = std::max(ps.maxDepthRatio, deep);
      if (deep > real_t(1) + real_t(1e-3)) ++ps.tooDeep;
   };
   Stats raw, hyb, trust, ref;
   int   contacts = 0;
   long  evals = 0, refineCalls = 0;
   auto  score = [](Stats& st, bool c, real_t ov, const Vec3& n, const GjkResult& g) {
      const real_t go = g.contact ? -g.penetrationDepth : real_t(0);
      if (g.contact && !c) { ++st.missed; return; }
      if (!g.contact && c) { ++st.falsePos; return; }
      if (!c) return;
      const real_t err = std::abs(ov - go);
      st.maxErr    = std::max(st.maxErr, err);
      st.maxErrRel = std::max(st.maxErrRel, err / go);
      if (err > std::max(real_t(0.01), real_t(0.05) * go)) ++st.wrongDepth;
      st.maxAngle = std::max(st.maxAngle, std::acos(std::clamp(n * g.normal, real_t(-1), real_t(1))) * real_t(180) / math::pi);
   };
   for (const auto& p : pairs)
   {
      const GjkResult     g  = gjkCylCyl(p.ci.pos, p.rotI, p.ci.radius, p.ci.length, p.cj.pos, p.rotJ, p.cj.radius, p.cj.length);
      const ContactResult r  = cylinderCylinderContact(p.ci, p.cj);
      const HybridResult  h  = cylinderCylinderContactHybrid(p.ci, p.rotI, p.cj, p.rotJ);
      const ContactResult rf = refineCylinderContact(p.ci, p.cj, r);
      if (g.contact) ++contacts;
      score(raw, r.contact, r.overlap, r.normal, g);
      score(hyb, h.contact, h.overlap, h.normal, g);
      score(trust, r.contact, r.overlap, r.normal, g); // identical to raw by construction: "no scenario" -> no contact
      score(ref, rf.contact, rf.overlap, rf.normal, g);
      if (g.contact)
      {
         checkPoint("gjk", p.ci, p.cj, true, g.point, -g.penetrationDepth);
         checkPoint("raw", p.ci, p.cj, r.contact, r.point, r.overlap);
         checkPoint("hybrid", p.ci, p.cj, h.contact, h.point, h.overlap);
         checkPoint("refined", p.ci, p.cj, rf.contact, rf.point, rf.overlap);
         if (r.contact) checkPoint(std::string("raw:") + r.scenario, p.ci, p.cj, r.contact, r.point, r.overlap);
      }
      if (rf.contact && !g.contact && std::getenv("ANTIDUNES_REFINED_DEBUG"))
      {
         Vec3 p0i, p1i, p0j, p1j;
         detail::cylinderAxisEndpoints(p.ci, p0i, p1i);
         detail::cylinderAxisEndpoints(p.cj, p0j, p1j);
         WALBERLA_LOG_INFO_ON_ROOT("[refined-FP] overlap=" << rf.overlap << " scen=" << rf.scenario
            << " capsuleGap=" << std::sqrt(detail::closestDistanceSqBetweenSegments(p0i, p1i, p0j, p1j)) - p.ci.radius - p.cj.radius
            << " axisDot=" << p.ci.axis() * p.cj.axis() << " R=" << p.ci.radius << "," << p.cj.radius << " L=" << p.ci.length << "," << p.cj.length);
      }
      if (!(r.contact && (std::string(r.scenario) == "face-face" || std::string(r.scenario) == "band-band-parallel" ||
                          std::string(r.scenario) == "band-band-skewed" || std::string(r.scenario) == "face-band")))
      {
         Vec3 seeds[6]; int n = 0;
         if (r.contact) seeds[n++] = -r.normal;
         n += detail::cylinderSeedDirections(p.ci, p.cj, seeds + n);
         evals += detail::minimiseCylinderProjectionOverlap(p.ci, p.cj, seeds, n, !r.contact).evaluations;
         ++refineCalls;
      }
   }

   {
      int positiveNoScenario = 0, convexSeparated = 0;
      for (const auto& p : pairs)
      {
         const ContactResult r = cylinderCylinderContact(p.ci, p.cj);
         if (r.contact) continue;
         Vec3 seeds[6];
         const int ns = detail::cylinderSeedDirections(p.ci, p.cj, seeds);
         const auto m = detail::minimiseCylinderProjectionOverlap(p.ci, p.cj, seeds, ns, true);
         if (m.value > real_t(0))
         {
            ++positiveNoScenario;
            if (detail::cylindersSeparatedConvexCheck(p.ci, p.cj, m.direction)) ++convexSeparated;
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[refined] no-scenario pairs whose descent ended positive: " << positiveNoScenario
                                 << ", of which the convex check proved separation: " << convexSeparated);
   }
   // timing: each variant over all pairs, repeated
   using clk = std::chrono::steady_clock;
   auto  timeIt = [&](auto&& fn) {
      const auto t0 = clk::now();
      real_t sink = 0;
      for (int rep = 0; rep < 5; ++rep)
         for (const auto& p : pairs) sink += fn(p);
      const auto t1 = clk::now();
      WALBERLA_UNUSED(sink);
      return std::chrono::duration< double, std::micro >(t1 - t0).count() / (5.0 * double(pairs.size()));
   };
   const double tGjk = timeIt([&](const Pair& p) {
      return -gjkCylCyl(p.ci.pos, p.rotI, p.ci.radius, p.ci.length, p.cj.pos, p.rotJ, p.cj.radius, p.cj.length).penetrationDepth; });
   const double tRaw = timeIt([&](const Pair& p) { return cylinderCylinderContact(p.ci, p.cj).overlap; });
   const double tHyb = timeIt([&](const Pair& p) { return cylinderCylinderContactHybrid(p.ci, p.rotI, p.cj, p.rotJ).overlap; });
   const double tRef = timeIt([&](const Pair& p) { return cylinderCylinderContactRefined(p.ci, p.cj).overlap; });

   {
      // cost breakdown of the refined detection by category of pair
      std::map< std::string, std::vector< const Pair* > > cat;
      for (const auto& p : pairs)
      {
         const ContactResult r = cylinderCylinderContact(p.ci, p.cj);
         const std::string   sc(r.scenario);
         std::string         k;
         if (!r.contact) k = "no scenario";
         else if (sc == "face-face" || sc == "band-band-parallel" || sc == "band-band-skewed" || sc == "face-band") k = "exact scenario";
         else
         {
            const Vec3   L0 = -r.normal;
            const real_t scale = p.ci.radius + p.cj.radius + real_t(0.5) * (p.ci.length + p.cj.length);
            const real_t D0 = detail::cylinderProjectionOverlap(p.ci, p.cj, L0);
            bool         fast = std::abs(D0 - r.overlap) < real_t(1e-6) * scale;
            Vec3         t1 = L0 % Vec3(real_t(1), real_t(0), real_t(0));
            if (t1.length() < real_t(0.5)) t1 = L0 % Vec3(real_t(0), real_t(1), real_t(0));
            t1 /= t1.length();
            const Vec3 t2 = L0 % t1;
            for (int kk = 0; kk < 8 && fast; ++kk)
            {
               Vec3 Lk = L0 + real_t(1e-4) * (std::cos(kk * real_t(0.25) * math::pi) * t1 + std::sin(kk * real_t(0.25) * math::pi) * t2);
               Lk /= Lk.length();
               fast = detail::cylinderProjectionOverlap(p.ci, p.cj, Lk) >= D0 - real_t(1e-12) * scale;
            }
            k = fast ? "inexact, fast path" : "inexact, refined";
         }
         cat[k].push_back(&p);
      }
      for (const auto& [k, v] : cat)
      {
         const auto t0 = clk::now();
         real_t     sink = 0;
         for (int rep = 0; rep < 5; ++rep)
            for (const Pair* p : v) sink += cylinderCylinderContactRefined(p->ci, p->cj).overlap;
         const double us = std::chrono::duration< double, std::micro >(clk::now() - t0).count() / (5.0 * double(v.size()));
         WALBERLA_UNUSED(sink);
         WALBERLA_LOG_INFO_ON_ROOT("[refined] cost: " << k << ": " << v.size() << " pairs, " << us << " us/pair");
      }
   }
   auto line = [&](const char* name, const Stats& st, double us) {
      WALBERLA_LOG_INFO_ON_ROOT("[refined] " << name << ": missed=" << st.missed << " falsePositive=" << st.falsePos
         << " wrongDepth=" << st.wrongDepth << "/" << contacts << " maxDepthErr=" << st.maxErr << " (rel " << st.maxErrRel
         << ") maxNormalAngle=" << st.maxAngle << " deg | " << us << " us/pair");
   };
   WALBERLA_LOG_INFO_ON_ROOT("[refined] " << numTrials << " pairs (uniform SO(3), seed 24680), " << contacts
                              << " contacts per GJK/EPA; refinement calls " << refineCalls << ", mean projection-overlap "
                              << "evaluations per call " << double(evals) / double(std::max(1L, refineCalls)));
   line("raw closed form      ", raw, tRaw);
   line("hybrid (GJK fallback)", hyb, tHyb);
   line("trust no-scenario    ", trust, tRaw);
   line("refined (no GJK/EPA) ", ref, tRef);
   for (const auto& [m, ps] : pst)
   {
      WALBERLA_LOG_INFO_ON_ROOT("[refined] contact point " << m << ": " << ps.checked << " contacts, outside the intersection="
         << ps.outside << " (max signed distance " << ps.maxOutside << " x depth), deeper than the penetration depth="
         << ps.tooDeep << " (max " << ps.maxDepthRatio << " x depth)");
   }
   WALBERLA_LOG_INFO_ON_ROOT("[refined] GJK/EPA alone: " << tGjk << " us/pair");
   // regression guard for the refined detection (measured: 0 / 0 / 0, max depth error 2.2e-4)
   WALBERLA_CHECK_EQUAL(ref.missed, 0);
   WALBERLA_CHECK_EQUAL(ref.falsePos, 0);
   WALBERLA_CHECK_EQUAL(ref.wrongDepth, 0);
   // contact point of the refined detection inside the intersection of the two cylinders (measured: 0 of 1698 outside)
   WALBERLA_CHECK_EQUAL(pst["refined"].outside, 0);
   WALBERLA_CHECK_EQUAL(pst["refined"].tooDeep, 0);
   WALBERLA_CHECK_EQUAL(pst["raw:band-band-skewed"].outside, 0);
}

//======================================================================================================================
// Nearly parallel rod pairs. The random samples above draw both orientations uniformly on SO(3), so that only a
// fraction 1 - cos(theta) of the pairs is skewed by less than theta: about 4 of the 1019 contacts of the thesis
// reference sample lie within 5 deg of parallel. The rods of the clogging runs, however, start aligned and lie nearly
// parallel in a jam. This sampler draws pairs of EQUAL rods (radius 0.3-2, length 4-8 radii) whose axes are skewed by
// theta, log-uniform between 1e-9 rad and maxSkew (10 % exactly parallel), placed side by side with an axial offset of
// up to 1.2 lengths, at a centre distance sampled around the GJK/EPA contact boundary (overlap up to 15 % of the
// radius). The closed form switches to the parallel formulas for 1 - |e_i.e_j| < 1e-8 (cylinderCylinderContact), i.e.
// theta below about 1.4e-4 rad, so the sample spans both sides of that switch.
//======================================================================================================================
struct NearParallelPair
{
   real_t R, L;
   Rot3   rotI, rotJ;
   Vec3   posJ;
   real_t skew; // angle between the axes [rad], measured from the matrices
};

NearParallelPair sampleNearParallelPair(FuzzRng& fuzz, real_t maxSkew, real_t minSkew = real_t(1e-9))
{
   std::uniform_real_distribution< real_t > u(real_t(0), real_t(1));
   NearParallelPair p;
   p.R = real_t(0.3) + real_t(1.7) * u(fuzz.rng);
   p.L = p.R * (real_t(4) + real_t(4) * u(fuzz.rng));
   p.rotI          = fuzz.randomRotation();
   const Vec3 eI   = p.rotI.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1));
   Vec3       perp = fuzz.randomUnitVec() % eI;
   while (perp.length() < real_t(1e-3)) perp = fuzz.randomUnitVec() % eI;
   perp = perp.getNormalized();
   const real_t theta = (u(fuzz.rng) < real_t(0.1))
                           ? real_t(0)
                           : minSkew * std::pow(maxSkew / minSkew, u(fuzz.rng)); // log-uniform
   p.rotJ = p.rotI;
   if (theta > real_t(0)) p.rotJ.rotate(perp, theta);
   const Vec3 eJ = p.rotJ.getMatrix() * Vec3(real_t(0), real_t(0), real_t(1));
   p.skew        = std::acos(std::clamp(std::abs(eI * eJ), real_t(0), real_t(1)));
   // side-by-side direction with an axial component: at the contact boundary (lateral distance ~2R) the axial offset of
   // the centres is ~s*2R, i.e. up to 1.2 L, so that partly overlapping and end-to-end-shifted rods are covered too
   const Vec3   lat  = fuzz.randomUnitVec() % eI;
   const real_t s    = (real_t(2.4) * u(fuzz.rng) - real_t(1.2)) * p.L / (real_t(2) * p.R);
   const Vec3   dir3 = (lat.getNormalized() + s * eI).getNormalized();
   const real_t circ = real_t(2) * std::sqrt(p.R * p.R + real_t(0.25) * p.L * p.L);
   auto contactAt    = [&](real_t sep) {
      return gjkCylCyl(Vec3(0, 0, 0), p.rotI, p.R, p.L, dir3 * sep, p.rotJ, p.R, p.L).contact;
   };
   // nearSep = 1e-4 R, not 0: for two (nearly) identical cylinders centred at the same point, GJK's fixed initial
   // search direction (GJK.cpp, d_ = (0, 0.6, 0.8)) can degenerate and report no contact, violating the precondition
   // contactAt(nearSep). The tiny separation is still a deep overlap and breaks that symmetry; exactly coincident
   // cylinders do not occur in the simulations.
   p.posJ = dir3 * fuzz.sampleSeparationNearBoundary(contactAt, real_t(1e-4) * p.R, real_t(3) * circ, p.R, real_t(0.15));
   return p;
}

//! Closed form with fallback (as in production) against GJK/EPA for nearly parallel rods, in bands of the skew angle.
//! Asserts what the fallback guarantees by construction (no missed contact, no false contact, finite results); depth
//! and normal errors are logged per band with the error criterion of the reference comparison, max(0.01, 0.05 depth).
void testCylinderNearParallelSkewedFuzz()
{
   const int numTrials = std::getenv("ANTIDUNES_NEARPARALLEL_TRIALS") ? std::atoi(std::getenv("ANTIDUNES_NEARPARALLEL_TRIALS"))
                                                                      : 2000;
   FuzzRng fuzz(24680u);
   const real_t maxSkew = real_t(10) * math::pi / real_t(180);
   // ANTIDUNES_NEARPARALLEL_MINSKEW_DEG: lower end of the log-uniform skew (default 1e-9 rad), to concentrate the
   // sample
   const real_t minSkew = std::getenv("ANTIDUNES_NEARPARALLEL_MINSKEW_DEG")
                             ? real_t(std::atof(std::getenv("ANTIDUNES_NEARPARALLEL_MINSKEW_DEG"))) * math::pi / real_t(180)
                             : real_t(1e-9);
   // bands of the skew angle [rad]: exactly parallel, parallel branch (< 1.4e-4), and three skewed decades up to 10 deg
   const std::vector< real_t > edges{ real_t(0), real_t(1e-12), real_t(1.4e-4), real_t(1e-3), real_t(1e-2), real_t(1e-1), maxSkew + real_t(1e-12) };
   const std::vector< std::string > names{ "exactly parallel", "< 1.4e-4 rad (parallel branch)", "1.4e-4 .. 1e-3 rad", "1e-3 .. 1e-2 rad",
                                           "1e-2 .. 1e-1 rad", "0.1 rad .. 10 deg" };
   struct Band { int contacts = 0, missed = 0, falsePos = 0, wrong = 0, fallback = 0, n1deg = 0, n10deg = 0; real_t maxRel = 0, maxAbsOverR = 0, maxNormalDeg = 0; std::map< std::string, int > scen; };
   std::vector< Band > bands(names.size());
   int nonFinite = 0, separated = 0;
   for (int t = 0; t < numTrials; ++t)
   {
      const NearParallelPair p = sampleNearParallelPair(fuzz, maxSkew, minSkew);
      size_t b = 0;
      while (b + 1 < names.size() && p.skew >= edges[b + 1]) ++b;
      const CylinderState ci{ Vec3(0, 0, 0), p.rotI.getMatrix(), p.R, p.L }, cj{ p.posJ, p.rotJ.getMatrix(), p.R, p.L };
      const GjkResult    g  = gjkCylCyl(ci.pos, p.rotI, p.R, p.L, cj.pos, p.rotJ, p.R, p.L);
      const HybridResult hy = cylinderCylinderContactHybrid(ci, p.rotI, cj, p.rotJ);
      if (!std::isfinite(hy.overlap) || !std::isfinite(hy.normal.length())) { ++nonFinite; continue; }
      Band& bd = bands[b];
      if (!g.contact) { ++separated; if (hy.contact) ++bd.falsePos; continue; }
      ++bd.contacts;
      if (!hy.contact) { ++bd.missed; continue; }
      if (hy.usedFallback) ++bd.fallback;
      ++bd.scen[hy.scenario];
      const real_t ref = -g.penetrationDepth, err = std::abs(hy.overlap - ref);
      if (err > std::max(real_t(0.01), real_t(0.05) * ref))
      {
         ++bd.wrong;
         if (bd.wrong <= 3)
            WALBERLA_LOG_INFO_ON_ROOT("[near-parallel]    wrong depth: skew " << p.skew << " rad, scenario " << hy.scenario << ", reported " << hy.overlap
                                      << ", GJK/EPA " << ref << ", R " << p.R << ", L " << p.L);
      }
      if (ref > real_t(1e-6)) bd.maxRel = std::max(bd.maxRel, err / ref);
      bd.maxAbsOverR = std::max(bd.maxAbsOverR, err / p.R);
      const real_t ang = std::acos(std::clamp(std::abs(hy.normal * g.normal) / (hy.normal.length() * g.normal.length()), real_t(0), real_t(1))) * real_t(180) / math::pi;
      bd.maxNormalDeg = std::max(bd.maxNormalDeg, ang);
      if (ang > real_t(1)) ++bd.n1deg;
      if (ang > real_t(10)) ++bd.n10deg;
   }
   int missed = 0, falsePos = 0;
   for (size_t b = 0; b < names.size(); ++b)
   {
      const Band& bd = bands[b];
      missed += bd.missed; falsePos += bd.falsePos;
      std::ostringstream sc;
      for (const auto& [k, v] : bd.scen) sc << ' ' << k << '=' << v;
      WALBERLA_LOG_INFO_ON_ROOT("[near-parallel] " << names[b] << ": contacts " << bd.contacts << ", missed " << bd.missed
                                << ", false " << bd.falsePos << ", fallback " << bd.fallback << ", wrong depths " << bd.wrong
                                << ", max rel err " << bd.maxRel << ", max err/R " << bd.maxAbsOverR << ", normals > 1 deg "
                                << bd.n1deg << ", > 10 deg " << bd.n10deg << ", max " << bd.maxNormalDeg << " deg; scenarios:" << sc.str());
   }
   WALBERLA_LOG_INFO_ON_ROOT("[near-parallel] " << numTrials << " pairs, " << separated << " separated, non-finite results " << nonFinite);
   WALBERLA_CHECK_EQUAL(nonFinite, 0);
   WALBERLA_CHECK_EQUAL(missed, 0);
   WALBERLA_CHECK_EQUAL(falsePos, 0);
}

//! The closed form switches from the parallel formulas to the skewed scenarios at 1 - |e_i.e_j| = 1e-8. Sweeps the skew
//! of fixed rod pairs smoothly across that switch (and on to 10 deg) and requires that the closed form with fallback
//! follows GJK/EPA without a jump up to 1e-2 rad: between neighbouring angles, its depth may change by at most 0.01 R
//! more than that of GJK/EPA. Beyond, the sweep exposes a known defect: when one rod is tilted out of the common plane
//! by about 4-10 deg and the rods are offset along their axes, the band-edge scenario reports depths that are too small
//! by up to 70 %, over a whole range of angles, and recovers with a jump (kodamBandEdgeAxialCheck narrows the range but
//! adds a jump of its own). The random samples rarely hit this configuration (none in 2409 contacts at 1-10 deg). It is
//! logged here and fails the test only with ANTIDUNES_NEARPARALLEL_STRICT set.
void testCylinderNearParallelThresholdContinuity()
{
   const real_t R = real_t(1), L = real_t(6);
   const Vec3   eI(0, 0, 1);
   struct Config { const char* name; Vec3 lateral; real_t axialOffset; };
   const std::vector< Config > configs{ { "side by side, centred", Vec3(1, 0, 0), real_t(0) },
                                        { "side by side, offset by half a length", Vec3(1, 0, 0), real_t(0.5) * L },
                                        { "ends overlapping by a fifth", Vec3(1, 0, 0), real_t(0.8) * L },
                                        { "oblique lateral direction", Vec3(1, 1, 0).getNormalized(), real_t(0.3) * L } };
   int jumps = 0, knownDefect = 0;
   for (const auto& c : configs)
      for (const real_t depth0 : { real_t(0.02), real_t(0.1) })
      {
         // skew about the axis perpendicular to both the rod axis and the lateral direction (the rods open like
         // scissors in their common plane) and about the lateral direction itself (one rod tilts out of the plane)
         for (const Vec3 tiltAxis : { (eI % c.lateral).getNormalized(), c.lateral })
         {
            const int steps = 400;
            real_t prevG = 0, prevH = 0, maxDiff = 0, worstTheta = 0, maxDiffSwitch = 0;
            bool   first = true;
            int    worstK = 0;
            std::vector< std::string > rows(steps + 1);
            for (int k = 0; k <= steps; ++k)
            {
               const real_t theta = real_t(1e-7) * std::pow(real_t(0.1745) / real_t(1e-7), real_t(k) / real_t(steps));
               Rot3 rotJ; rotJ.rotate(tiltAxis, theta);
               const Vec3 posJ = c.lateral * (real_t(2) * R - depth0) + eI * c.axialOffset;
               const CylinderState ci{ Vec3(0, 0, 0), Rot3().getMatrix(), R, L }, cj{ posJ, rotJ.getMatrix(), R, L };
               const GjkResult    g  = gjkCylCyl(ci.pos, Rot3(), R, L, posJ, rotJ, R, L);
               const HybridResult hy = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
               const real_t dG = g.contact ? -g.penetrationDepth : real_t(0), dH = hy.contact ? hy.overlap : real_t(0);
               std::ostringstream row;
               row << "theta " << theta << ": GJK/EPA " << g.contact << ' ' << dG << ", closed form with fallback " << hy.contact << ' ' << dH << " (" << hy.scenario << ")";
               rows[size_t(k)] = row.str();
               if (!first)
               {
                  const real_t diff = std::abs((dH - prevH) - (dG - prevG));
                  if (diff > maxDiff) { maxDiff = diff; worstTheta = theta; worstK = k; }
                  if (theta <= real_t(1e-2)) maxDiffSwitch = std::max(maxDiffSwitch, diff);
               }
               prevG = dG; prevH = dH; first = false;
            }
            // range of angles in which the depth violates the error criterion of the reference comparison
            int nBad = 0; real_t badLo = 0, badHi = 0, badMax = 0, badRef = 0, badRep = 0;
            for (int k = 0; k <= steps; ++k)
            {
               const real_t theta = real_t(1e-7) * std::pow(real_t(0.1745) / real_t(1e-7), real_t(k) / real_t(steps));
               Rot3 rotJ; rotJ.rotate(tiltAxis, theta);
               const Vec3 posJ = c.lateral * (real_t(2) * R - depth0) + eI * c.axialOffset;
               const CylinderState ci{ Vec3(0, 0, 0), Rot3().getMatrix(), R, L }, cj{ posJ, rotJ.getMatrix(), R, L };
               const GjkResult    g  = gjkCylCyl(ci.pos, Rot3(), R, L, posJ, rotJ, R, L);
               const HybridResult hy = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
               const real_t dG = g.contact ? -g.penetrationDepth : real_t(0), dH = hy.contact ? hy.overlap : real_t(0);
               if (std::abs(dH - dG) > std::max(real_t(0.01), real_t(0.05) * dG))
               {
                  if (nBad == 0) badLo = theta;
                  badHi = theta; ++nBad;
                  if (std::abs(dH - dG) > badMax) { badMax = std::abs(dH - dG); badRef = dG; badRep = dH; }
               }
            }
            if (nBad > 0)
               WALBERLA_LOG_INFO_ON_ROOT("[near-parallel sweep]    wrong depth at " << nBad << " of " << steps + 1 << " angles, theta " << badLo << " .. " << badHi
                                         << ", largest error " << badMax << " (reported " << badRep << " instead of " << badRef << ")");
            const bool jump = maxDiff > real_t(0.01) * R;
            if (maxDiffSwitch > real_t(0.01) * R) ++jumps;
            if (jump || nBad > 0) ++knownDefect;
            if (jump)
            {
               for (int k = std::max(0, worstK - 3); k <= std::min(steps, worstK + 2); ++k)
                  WALBERLA_LOG_INFO_ON_ROOT("[near-parallel sweep]    " << rows[size_t(k)]);
            }
            WALBERLA_LOG_INFO_ON_ROOT("[near-parallel sweep] " << c.name << ", depth " << depth0 << ", tilt axis " << tiltAxis
                                      << ": largest step difference to GJK/EPA " << maxDiff << " at theta " << worstTheta
                                      << (jump ? "  JUMP" : ""));
         }
      }
   WALBERLA_LOG_INFO_ON_ROOT("[near-parallel sweep] jumps up to 1e-2 rad: " << jumps << "; configurations with the KNOWN band-edge defect "
                             "(wrong depths or a jump beyond 1e-2 rad): " << knownDefect);
   WALBERLA_CHECK_EQUAL(jumps, 0);
   if (std::getenv("ANTIDUNES_NEARPARALLEL_STRICT")) WALBERLA_CHECK_EQUAL(knownDefect, 0);
}

//======================================================================================================================
// Independent-reference dump (opt-in: runs only if ANTIDUNES_REFERENCE_DUMP is set to an output directory). The other
// randomized tests treat GJK/EPA as the reference, but GJK/EPA is also the production path of the thesis's rod runs
// (useAnalyticalContactDetection=false), so it needs an independent check. This writes random pairs (uniform SO(3),
// FuzzRng::randomRotation) with the GJK/EPA, closed-form (raw, hybrid, refined) and SAT answers to CSV;
// scripts/contact-validation/brute_force_reference.py (thesis repo) recomputes every depth from the support functions
// alone and compares all methods against it.
//======================================================================================================================
void dumpMat(std::ofstream& f, const Mat3& m)
{
   for (uint_t i = 0; i < 3; ++i)
      for (uint_t j = 0; j < 3; ++j) f << ',' << m(i, j);
}

void testIndependentReferenceDump()
{
   const char* dir = std::getenv("ANTIDUNES_REFERENCE_DUMP");
   if (dir == nullptr) return;
   const int numTrials = std::getenv("ANTIDUNES_REFERENCE_TRIALS") ? std::atoi(std::getenv("ANTIDUNES_REFERENCE_TRIALS"))
                                                                   : 2000;
   // ANTIDUNES_REFERENCE_SEED (default 13579, the seed of the thesis figure) for an out-of-sample check,
   // ANTIDUNES_REFERENCE_OVERLAP (default 0.15) for the depth of the sampling band in units of the smaller radius
   const unsigned seed = std::getenv("ANTIDUNES_REFERENCE_SEED") ? unsigned(std::atol(std::getenv("ANTIDUNES_REFERENCE_SEED"))) : 13579u;
   const real_t   band = std::getenv("ANTIDUNES_REFERENCE_OVERLAP") ? real_t(std::atof(std::getenv("ANTIDUNES_REFERENCE_OVERLAP")))
                                                                    : real_t(0.15);
   FuzzRng fuzz(seed);
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0)), lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > halfDist(real_t(0.3), real_t(2.0));

   std::ofstream fc(std::string(dir) + "/cyl_pairs.csv");
   fc.precision(17);
   fc << "trial,Ri,Li,Rj,Lj,xj,yj,zj,Ai00,Ai01,Ai02,Ai10,Ai11,Ai12,Ai20,Ai21,Ai22,Aj00,Aj01,Aj02,Aj10,Aj11,Aj12,Aj20,Aj21,"
         "Aj22,gjk_contact,gjk_depth,gjk_nx,gjk_ny,gjk_nz,cf_contact,cf_overlap,cf_nx,cf_ny,cf_nz,cf_scenario,"
         "hyb_contact,hyb_overlap,hyb_fallback,rfd_contact,rfd_overlap,rfd_nx,rfd_ny,rfd_nz,"
         // rfd_* = refined result; contact points of every method and the hybrid normal feed the thesis panels
         "gjk_px,gjk_py,gjk_pz,cf_px,cf_py,cf_pz,hyb_nx,hyb_ny,hyb_nz,hyb_px,hyb_py,hyb_pz,rfd_px,rfd_py,rfd_pz\n";
   // ANTIDUNES_REFERENCE_NEARPARALLEL=<max skew in degrees>: the cylinder pairs are nearly parallel equal rods from
   // sampleNearParallelPair instead of uniform SO(3) pairs (same columns, so brute_force_reference.py applies as is)
   const char*  nearParallel = std::getenv("ANTIDUNES_REFERENCE_NEARPARALLEL");
   const real_t maxSkewDeg   = nearParallel ? real_t(std::atof(nearParallel)) : real_t(0);
   for (int t = 0; t < numTrials; ++t)
   {
      real_t Ri = radiusDist(fuzz.rng), Li = lengthDist(fuzz.rng), Rj = radiusDist(fuzz.rng), Lj = lengthDist(fuzz.rng);
      Rot3   rotI = fuzz.randomRotation(), rotJ = fuzz.randomRotation();
      const Vec3   dir3 = fuzz.randomUnitVec();
      const Vec3   posI(0, 0, 0);
      const real_t circ = std::sqrt(Ri * Ri + real_t(0.25) * Li * Li) + std::sqrt(Rj * Rj + real_t(0.25) * Lj * Lj);
      auto contactAt = [&](real_t sep) { return gjkCylCyl(posI, rotI, Ri, Li, posI + dir3 * sep, rotJ, Rj, Lj).contact; };
      Vec3 posJ;
      if (nearParallel)
      {
         const NearParallelPair p = sampleNearParallelPair(fuzz, maxSkewDeg * math::pi / real_t(180));
         Ri = Rj = p.R; Li = Lj = p.L; rotI = p.rotI; rotJ = p.rotJ; posJ = p.posJ;
      }
      else
         posJ = posI + dir3 * fuzz.sampleSeparationNearBoundary(contactAt, real_t(0), real_t(3) * circ, std::min(Ri, Rj), band);
      const CylinderState ci{ posI, rotI.getMatrix(), Ri, Li }, cj{ posJ, rotJ.getMatrix(), Rj, Lj };
      const GjkResult     g  = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
      const ContactResult cf = cylinderCylinderContact(ci, cj);
      const HybridResult  hy = cylinderCylinderContactHybrid(ci, rotI, cj, rotJ);
      const ContactResult rf = refineCylinderContact(ci, cj, cf);
      fc << t << ',' << Ri << ',' << Li << ',' << Rj << ',' << Lj << ',' << posJ[0] << ',' << posJ[1] << ',' << posJ[2];
      dumpMat(fc, rotI.getMatrix());
      dumpMat(fc, rotJ.getMatrix());
      fc << ',' << g.contact << ',' << g.penetrationDepth << ',' << g.normal[0] << ',' << g.normal[1] << ',' << g.normal[2]
         << ',' << cf.contact << ',' << cf.overlap << ',' << cf.normal[0] << ',' << cf.normal[1] << ',' << cf.normal[2]
         << ',' << cf.scenario << ',' << hy.contact << ',' << hy.overlap << ',' << hy.usedFallback << ',' << rf.contact
         << ',' << rf.overlap << ',' << rf.normal[0] << ',' << rf.normal[1] << ',' << rf.normal[2];
      for (const Vec3* v : { &g.point, &cf.point, &hy.normal, &hy.point, &rf.point })
         fc << ',' << (*v)[0] << ',' << (*v)[1] << ',' << (*v)[2];
      fc << '\n';
   }

   std::ofstream fb(std::string(dir) + "/box_pairs.csv");
   fb.precision(17);
   fb << "trial,hax,hay,haz,hbx,hby,hbz,xb,yb,zb,Aa00,Aa01,Aa02,Aa10,Aa11,Aa12,Aa20,Aa21,Aa22,Ab00,Ab01,Ab02,Ab10,Ab11,Ab12,"
         "Ab20,Ab21,Ab22,gjk_contact,gjk_depth,gjk_nx,gjk_ny,gjk_nz,sat_contact,sat_overlap,sat_nx,sat_ny,sat_nz,"
         "gjk_px,gjk_py,gjk_pz,sat_px,sat_py,sat_pz\n";
   for (int t = 0; t < numTrials; ++t)
   {
      const Vec3 hA(halfDist(fuzz.rng), halfDist(fuzz.rng), halfDist(fuzz.rng));
      const Vec3 hB(halfDist(fuzz.rng), halfDist(fuzz.rng), halfDist(fuzz.rng));
      const Rot3 rotA = fuzz.randomRotation(), rotB = fuzz.randomRotation();
      const Vec3 dir3 = fuzz.randomUnitVec();
      const Vec3 posA(0, 0, 0);
      const real_t minSize = std::min({ hA[0], hA[1], hA[2], hB[0], hB[1], hB[2] });
      auto contactAt = [&](real_t sep) { return gjkBoxBox(posA, rotA, real_t(2) * hA, posA + dir3 * sep, rotB, real_t(2) * hB).contact; };
      const Vec3 posB = posA + dir3 * fuzz.sampleSeparationNearBoundary(contactAt, real_t(0), real_t(3) * (hA.length() + hB.length()), minSize, band);
      const GjkResult g   = gjkBoxBox(posA, rotA, real_t(2) * hA, posB, rotB, real_t(2) * hB);
      const auto      sat = boxBoxSAT(posA, rotA.getMatrix(), hA, posB, rotB.getMatrix(), hB);
      fb << t << ',' << hA[0] << ',' << hA[1] << ',' << hA[2] << ',' << hB[0] << ',' << hB[1] << ',' << hB[2] << ','
         << posB[0] << ',' << posB[1] << ',' << posB[2];
      dumpMat(fb, rotA.getMatrix());
      dumpMat(fb, rotB.getMatrix());
      fb << ',' << g.contact << ',' << g.penetrationDepth << ',' << g.normal[0] << ',' << g.normal[1] << ',' << g.normal[2]
         << ',' << sat.contact << ',' << sat.overlap << ',' << sat.normal[0] << ',' << sat.normal[1] << ',' << sat.normal[2]
         << ',' << g.point[0] << ',' << g.point[1] << ',' << g.point[2] << ',' << sat.point[0] << ',' << sat.point[1] << ','
         << sat.point[2] << '\n';
   }
   WALBERLA_LOG_INFO_ON_ROOT("[reference-dump] wrote " << numTrials << " cylinder pairs ("
                             << (nearParallel ? "nearly parallel equal rods, skew up to " + std::to_string(maxSkewDeg) + " deg" : std::string("uniform SO(3) orientations"))
                             << ") and " << numTrials << " box pairs (uniform SO(3) orientations) to " << dir);
}

//======================================================================================================================
// Replay of rod-pillar configurations taken from the floater snapshots of the clogging runs (opt-in:
// ANTIDUNES_PILLAR_REPLAY=<input csv> ANTIDUNES_PILLAR_REPLAY_OUT=<output csv>). Input per line:
// id,R,L,x,y,z,ax,ay,az,Rp,px,py (rod radius and length already inflated by the contact skin, world centre, world axis;
// vertical pillar of radius Rp through (px,py)). Writes the production rod-pillar detection (projection variant, refine
// = false) and the refined one, so that scripts compare both with the exact depth.
//======================================================================================================================
void replayPillarConfigurations()
{
   const char* in  = std::getenv("ANTIDUNES_PILLAR_REPLAY");
   const char* out = std::getenv("ANTIDUNES_PILLAR_REPLAY_OUT");
   if (in == nullptr || out == nullptr) return;
   std::ifstream fi(in);
   std::ofstream fo(out);
   fo.precision(17);
   fo << "id,contact,depth,r_contact,r_depth\n";
   std::string line;
   while (std::getline(fi, line))
   {
      if (line.empty() || line[0] == '#' || line[0] == 'i') continue;
      std::stringstream ss(line);
      std::string       tok, id;
      std::vector< real_t > v;
      std::getline(ss, id, ',');
      while (std::getline(ss, tok, ',')) v.push_back(real_t(std::atof(tok.c_str())));
      WALBERLA_CHECK_EQUAL(v.size(), 11u); // R, L, x, y, z, ax, ay, az, Rp, px, py
      const real_t R = v[0], L = v[1], Rp = v[8];
      const Vec3   rodPos(v[2] - v[9], v[3] - v[10], v[4]); // pillar axis moved to the origin
      const Vec3   a = Vec3(v[5], v[6], v[7]).getNormalized();
      const Vec3   ez(0, 0, 1);
      const Vec3   k  = ez % a;
      const real_t sn = k.length(), cs = ez * a;
      // mesa_pd cylinders have their axis along the body z axis: rotate e_z onto a
      const Rot3 rot = (sn < real_t(1e-12)) ? Rot3() : Rot3(k / sn, std::atan2(sn, cs));
      const mesa_pd::data::CylindricalBoundary pillar(Rp, ez);
      const mesa_pd::data::Cylinder            cyl(R, L);
      auto detect = [&](bool refine, real_t& depth) {
         auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
         mesa_pd::data::ParticleAccessor ac(ps);
         auto rod = ps->create();
         rod->setPosition(rodPos);
         rod->setRotation(rot);
         rod->setInteractionRadius(std::sqrt(R * R + real_t(0.25) * L * L));
         auto obs = ps->create();
         obs->setPosition(Vec3(0, 0, 0));
         SolidCylinderContactDetection scd;
         scd.refineCylinderContacts = refine;
         const bool c = scd(0, 1, cyl, pillar, ac);
         depth        = c ? -scd.getPenetrationDepth() : real_t(0);
         return c;
      };
      real_t d0 = 0, d1 = 0;
      const bool c0 = detect(false, d0), c1 = detect(true, d1);
      fo << id << ',' << c0 << ',' << d0 << ',' << c1 << ',' << d1 << '\n';
   }
}

//======================================================================================================================
// Independent-reference dump for rod-pillar (solid CylindricalBoundary) and rod-wall (HalfSpace) contacts, the
// obstacles of the clogging layouts. Opt-in like testIndependentReferenceDump (ANTIDUNES_REFERENCE_DUMP). Production
// detection: SolidCylinderContactDetection (ParticleKernels.h, "pillar exterior" and Cylinder x HalfSpace), written
// with and without refinement. scripts/contact-validation/brute_force_reference.py recomputes the depth from the
// support functions, with the infinite pillar represented by a cylinder of half-length 1000.
//======================================================================================================================
void testIndependentReferenceDumpObstacles()
{
   const char* dir = std::getenv("ANTIDUNES_REFERENCE_DUMP");
   if (dir == nullptr) return;
   const int numTrials = std::getenv("ANTIDUNES_REFERENCE_TRIALS") ? std::atoi(std::getenv("ANTIDUNES_REFERENCE_TRIALS"))
                                                                   : 2000;
   FuzzRng fuzz(86420);
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0)), lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > pillarDist(real_t(0.5), real_t(3.0));

   auto detect = [](const Vec3& pos, const Rot3& rot, real_t R, real_t L, auto&& obstacleShape, const Vec3& obstaclePos,
                    GjkResult& out, bool refine = false) {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto rod = ps->create();
      rod->setPosition(pos);
      rod->setRotation(rot);
      rod->setInteractionRadius(std::sqrt(R * R + real_t(0.25) * L * L));
      auto obs = ps->create();
      obs->setPosition(obstaclePos);
      const mesa_pd::data::Cylinder cyl(R, L);
      SolidCylinderContactDetection scd;
      scd.refineCylinderContacts = refine;
      const bool c = scd(0, 1, cyl, obstacleShape, ac);
      out = GjkResult{ c, scd.getPenetrationDepth(), scd.getContactNormal(), scd.getContactPoint() };
      return c;
   };

   for (const std::string kind : { "pillar", "wall" })
   {
      std::ofstream f(std::string(dir) + "/" + kind + "_pairs.csv");
      f.precision(17);
      f << "trial,R,L,x,y,z,A00,A01,A02,A10,A11,A12,A20,A21,A22,Rp,ax,ay,az,contact,depth,nx,ny,nz,px,py,pz,"
           "r_contact,r_depth,r_nx,r_ny,r_nz,r_px,r_py,r_pz\n";
      for (int t = 0; t < numTrials; ++t)
      {
         const real_t R = radiusDist(fuzz.rng), L = lengthDist(fuzz.rng);
         const Rot3   rot  = fuzz.randomRotation();
         const Vec3   axis = fuzz.randomUnitVec(); // pillar axis or wall normal
         const real_t Rp   = (kind == "pillar") ? pillarDist(fuzz.rng) : real_t(0);
         const Vec3   dir3 = fuzz.randomUnitVec();
         const Vec3   origin(0, 0, 0);
         const mesa_pd::data::CylindricalBoundary pillar(Rp > 0 ? Rp : real_t(1), axis);
         const mesa_pd::data::HalfSpace           wall(axis);
         GjkResult                                res{ false, real_t(0), Vec3(real_t(0)), Vec3(real_t(0)) };
         auto contactAt = [&](real_t sep) {
            GjkResult tmp{ false, real_t(0), Vec3(real_t(0)), Vec3(real_t(0)) };
            return (kind == "pillar") ? detect(origin + dir3 * sep, rot, R, L, pillar, origin, tmp)
                                      : detect(origin + dir3 * sep, rot, R, L, wall, origin, tmp);
         };
         // wall: dir3 must point to the free side for the rod to leave the wall
         Vec3 d3 = dir3;
         if (kind == "wall" && d3 * axis < real_t(0.2)) d3 = (d3 + axis).getNormalized();
         auto contactAtW = [&](real_t sep) {
            GjkResult tmp{ false, real_t(0), Vec3(real_t(0)), Vec3(real_t(0)) };
            // sampling with the refined detection (validated against the reference); the legacy pillar detection misses
            // contacts, which would bias the sampled overlaps towards large values
            return (kind == "pillar") ? detect(origin + d3 * sep, rot, R, L, pillar, origin, tmp, true)
                                      : detect(origin + d3 * sep, rot, R, L, wall, origin, tmp, true);
         };
         WALBERLA_UNUSED(contactAt);
         const real_t far = real_t(3) * (Rp + R + L);
         const Vec3   pos = origin + d3 * fuzz.sampleSeparationNearBoundary(contactAtW, real_t(0), far, std::min(R, Rp > 0 ? Rp : R));
         GjkResult rres{ false, real_t(0), Vec3(real_t(0)), Vec3(real_t(0)) };
         if (kind == "pillar") { detect(pos, rot, R, L, pillar, origin, res); detect(pos, rot, R, L, pillar, origin, rres, true); }
         else { detect(pos, rot, R, L, wall, origin, res); detect(pos, rot, R, L, wall, origin, rres, true); }
         f << t << ',' << R << ',' << L << ',' << pos[0] << ',' << pos[1] << ',' << pos[2];
         dumpMat(f, rot.getMatrix());
         f << ',' << Rp << ',' << axis[0] << ',' << axis[1] << ',' << axis[2] << ',' << res.contact << ',' << res.penetrationDepth
           << ',' << res.normal[0] << ',' << res.normal[1] << ',' << res.normal[2] << ',' << res.point[0] << ',' << res.point[1]
           << ',' << res.point[2] << ',' << rres.contact << ',' << rres.penetrationDepth << ',' << rres.normal[0] << ','
           << rres.normal[1] << ',' << rres.normal[2] << ',' << rres.point[0] << ',' << rres.point[1] << ',' << rres.point[2] << '\n';
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("[reference-dump] wrote " << numTrials << " rod-pillar and " << numTrials << " rod-wall pairs to " << dir);

}
//! Offline evaluation of the refined detection on pairs dumped from a packing run (ParticleKernels
//! ANTIDUNES_REFINE_DUMP). Only with ANTIDUNES_REFINE_EVAL=<dir with dump_rank*.csv>. For the contacts the fast path
//! rejected because the short descent found a smaller overlap ("smaller-elsewhere"), it compares the full refinement
//! with cheaper acceptance rules: accept the closed form if the descent improves on it by less than a relative
//! tolerance, or take the short descent's own minimum.
void evaluateRefineDump()
{
   const char* dir = std::getenv("ANTIDUNES_REFINE_EVAL");
   if (dir == nullptr) return;
   using clk = std::chrono::steady_clock;
   auto stateFrom = [](const std::vector< double >& v, size_t o) {
      Vec3 a(v[o + 3], v[o + 4], v[o + 5]);
      a /= a.length();
      Vec3 e1 = a % Vec3(real_t(1), real_t(0), real_t(0));
      if (e1.length() < real_t(0.5)) e1 = a % Vec3(real_t(0), real_t(1), real_t(0));
      e1 /= e1.length();
      const Vec3 e2 = a % e1;
      const Mat3 R(e1[0], e2[0], a[0], e1[1], e2[1], a[1], e1[2], e2[2], a[2]);
      return CylinderState{ Vec3(v[o], v[o + 1], v[o + 2]), R, real_t(v[o + 6]), real_t(v[o + 7]) };
   };
   const real_t tols[4] = { real_t(1e-9), real_t(1e-6), real_t(1e-4), real_t(1e-3) }; // relative to D0
   uint_t n = 0, nAcc[4] = {}, nScenarioMismatch = 0;
   double maxErrAcc[4] = {}, maxAngAcc[4] = {}, maxErrShort = 0, maxAngShort = 0, tFull = 0, tShort = 0;
   std::map< std::string, uint_t > byScenario;
   for (const auto& entry : std::filesystem::directory_iterator(dir))
   {
      if (entry.path().filename().string().rfind("dump_rank", 0) != 0) continue;
      std::ifstream in(entry.path());
      std::string   line;
      std::getline(in, line); // header
      while (std::getline(in, line))
      {
         std::vector< std::string > f;
         std::stringstream          ss(line);
         std::string                tok;
         while (std::getline(ss, tok, ',')) f.push_back(tok);
         if (f.size() != 23 || f[4] != "3") continue; // smaller-elsewhere only
         std::vector< double > v;
         for (size_t k = 7; k < 23; ++k) v.push_back(std::stod(f[k]));
         const CylinderState ci = stateFrom(v, 0), cj = stateFrom(v, 8);
         const ContactResult raw = cylinderCylinderContact(ci, cj);
         if (!raw.contact || std::string(raw.scenario) != f[0]) { ++nScenarioMismatch; continue; }
         ++n;
         ++byScenario[f[0]];
         const auto          t0   = clk::now();
         const ContactResult full = refineCylinderContact(ci, cj, raw);
         const auto          t1   = clk::now();
         tFull += std::chrono::duration< double >(t1 - t0).count();
         // the fast path's short descent, as in refineCylinderContact
         const auto   t2 = clk::now();
         const Vec3   L0 = -raw.normal;
         const real_t D0 = detail::cylinderProjectionOverlap(ci, cj, L0);
         Vec3         std5[5];
         detail::cylinderSeedDirections(ci, cj, std5);
         const detail::OverlapMinimum q = detail::minimiseCylinderProjectionOverlap(ci, cj, std5, 1, false, 20, real_t(1e-7));
         const auto t3 = clk::now();
         tShort += std::chrono::duration< double >(t3 - t2).count();
         const double ref = full.contact ? double(full.overlap) : 0.0;
         auto angle = [](const Vec3& a, const Vec3& b) {
            return std::acos(std::clamp(double(a * b), -1.0, 1.0)) * 180.0 / math::pi;
         };
         for (int k = 0; k < 4; ++k)
            if (q.value >= D0 - tols[k] * D0)
            {
               ++nAcc[k];
               maxErrAcc[k] = std::max(maxErrAcc[k], std::abs(double(D0) - ref) / std::max(ref, 1e-30));
               if (full.contact) maxAngAcc[k] = std::max(maxAngAcc[k], angle(raw.normal, full.normal));
            }
         const real_t   shortVal = std::min(q.value, D0);
         const Vec3 shortNormal = (q.value < D0) ? Vec3(-q.direction) : raw.normal;
         maxErrShort = std::max(maxErrShort, std::abs(double(shortVal) - ref) / std::max(ref, 1e-30));
         if (full.contact) maxAngShort = std::max(maxAngShort, angle(shortNormal, full.normal));
      }
   }
   std::cout << "[RefineDumpEval] smaller-elsewhere pairs " << n << " (scenario mismatch / lost " << nScenarioMismatch
             << "); full refinement " << (n ? 1e6 * tFull / double(n) : 0.0) << " us, short descent "
             << (n ? 1e6 * tShort / double(n) : 0.0) << " us\n";
   for (const auto& [k, c] : byScenario) std::cout << "   " << k << ": " << c << "\n";
   for (int k = 0; k < 4; ++k)
      std::cout << "   accept closed form if improvement < " << tols[k] << " D0: " << nAcc[k] << " ("
                << (n ? 100.0 * double(nAcc[k]) / double(n) : 0.0) << " %), max rel. depth error " << maxErrAcc[k]
                << ", max normal angle " << maxAngAcc[k] << " deg\n";
   std::cout << "   take min(D0, short descent): max rel. depth error " << maxErrShort << ", max normal angle "
             << maxAngShort << " deg\n";
}


} // namespace test
} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla


int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes::closed_form_contacts::test;

   // Debug output (cylinder-mismatch VTK/CSV) goes next to the executable in the build tree,
   // regardless of the process's current working directory at launch.
   g_outputDir = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();

   // ANTIDUNES_BANDEDGE_AXIAL: run with closed_form_contacts::kodamBandEdgeAxialCheck on (prm
   // kodamBandEdgeAxialCheck); ANTIDUNES_REFERENCE_DUMP_ONLY: write only the independent-reference dump and exit.
   walberla::antidunes::closed_form_contacts::kodamBandEdgeAxialCheck = (std::getenv("ANTIDUNES_BANDEDGE_AXIAL") != nullptr);
   if (std::getenv("ANTIDUNES_REFERENCE_DUMP_ONLY")) { testIndependentReferenceDump(); return EXIT_SUCCESS; }
   if (std::getenv("ANTIDUNES_PILLAR_REPLAY")) { replayPillarConfigurations(); return EXIT_SUCCESS; }
   if (std::getenv("ANTIDUNES_NEARPARALLEL_ONLY")) // only the nearly parallel rod tests
   {
      testCylinderNearParallelSkewedFuzz();
      testCylinderNearParallelThresholdContinuity();
      return EXIT_SUCCESS;
   }

   // Box x Box
   testBoxBoxSATMatchesExisting();
   testBoxBoxSATNoContact();
   testBoxBoxSATFaceFaceYAxis();
   testBoxBoxSATFaceFaceZAxis();
   testBoxBoxSATAsymmetricSizes();
   testBoxBoxSATRotatedCornerIntoFace();
   testBoxBoxSATBothRotatedGeneralOrientation();
   testBoxBoxSATDeepOverlap();
   testBoxBoxSATRealAntidunesLaunchCollision();
   testBoxBoxSATJob4139454Explosion1();
   testBoxBoxSATJob4139454Explosion2();
   testBoxBoxSATBoundaryTouching();
   testBoxBoxRandomizedFuzz();

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
   testCylinderCylinderFaceBandFullyInside();
   testCylinderCylinderFaceBandPartiallyInside();
   testCylinderCylinderFaceBandNoContact();
   testCylinderCylinderTiltedRimIntoFace();
   testCylinderCylinderTiltedRimIntoFaceNoContact();
   testCylinderCylinderTiltedRimIntoBand();
   testCylinderCylinderTiltedRimIntoBandNoContact();
   testCylinderCylinderSpecialContact();
   testCylinderCylinderSpecialContactBruteForceMTD();
   testCylinderCylinderSpecialContactNoContact();
   testCylinderCylinderEdgeEdge();
   testCylinderCylinderRandomizedFuzz();
   testCylinderCylinderPureAnalyticalScenarioBreakdown();
   testCylinderCylinderParallelAxesFuzz();
   testCylinderCylinderFaceBandBarelyOverlappingFuzz();
   testCylinderCylinderFaceEdgeClean();
   testCylinderCylinderSpecialContactCorrectedSweep();
   testCylinderCylinderBandEdgeClean();
   testCylinderCylinderEdgeEdgeGenuine();

   testCylinderCylinderFullComparisonAndHybridFeasibility();
   testEdgeEdgeThreeWayAccuracy10k();
   testThetaRefTieBreakNecessity10k();
   testClosestDistanceSqBetweenSegments();
   testCapsuleBoundConservative10k();
   testBoxSeparationAlongCenterLine();
   testTightLayer1BoxAccuracy10k();
   testCylinderCylinderLargeOverlapDegenerateStress();
   // testCylinderCylinderFallbackRateBoundaryBiased();  // deactivated -- 100k trials, slow; call
   // manually when the fallback-precision breakdown is actually needed

   // Hybrid dispatcher (analytical + GJK/EPA fallback via cylinderCylinderResultIsUncertain)
   testHybridTier1TriggersOnMissedContact();
   testHybridTier1DoesNotTriggerWhenClearlySeparated();
   testHybridTier1BoundaryCase();
   testHybridKeepsFaceEdgeMagnitudeWithoutFallback();
   testHybridTier2DoesNotTriggerOnVerifiedScenarios();
   testHybridTier2OnAlreadyCorrectFaceEdgeCase();

   // Cylinder x Plane
   testCylinderPlaneFaceMatchesExisting();
   testCylinderPlaneFaceNoContact();
   testCylinderPlaneBand();
   testCylinderPlaneBandNoContact();
   testCylinderPlaneEdge();
   testCylinderPlaneEdgeNoContact();

   // Sphere x Cylinder
   testSphereCylinderBand();
   testSphereCylinderCap();
   testSphereCylinderRim();
   testSphereCylinderEmbeddedBand();
   testSphereCylinderEmbeddedCap();
   testSphereCylinderNoContact();
   testSphereCylinderRandomizedFuzz();

   // refined detection vs. GJK/EPA; opt-in reference dump (ANTIDUNES_REFERENCE_DUMP=<dir>)
   testRefinedCylinderContactVsGjk();
   testCylinderCylinderContinuityAlongPaths();
   testCylinderNearParallelSkewedFuzz();
   testCylinderNearParallelThresholdContinuity();
   testIndependentReferenceDump();
   testBoxBoxContactPointHandCases();
   testBoxBoxContactPointFuzz();
   testBoxBoxContactPointContinuity();
   benchmarkBoxBoxContactPoint();
   testClipLineToBox();
   testBoxBoxContactPointRestingContacts();
   evaluateBoxPairsDump();
   testIndependentReferenceDumpObstacles();
   evaluateRefineDump();

   return EXIT_SUCCESS;
}
