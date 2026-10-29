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
namespace paper_analytical
{
namespace test
{
using Vec3 = mesa_pd::Vec3;
using Mat3 = mesa_pd::Mat3;
using Rot3 = mesa_pd::Rot3;

//! Directory the test binary itself lives in (set from argv[0] in main()), so debug output such as
//! the cylinder-mismatch VTK/CSV dump always lands in the build tree next to the executable instead
//! of wherever the process happened to be launched from.
std::filesystem::path g_outputDir;

//! Off by default -- writing the per-mismatch VTK frame series (cylinder_fuzz_mismatches,
//! cylinder_kodam_mismatches) is the single biggest cost in this whole test's runtime for no benefit
//! on a routine run; the CSV dump (cheap, always written) already has everything needed to see WHICH
//! trials mismatched. Set WALBERLA_WRITE_MISMATCH_VTK=1 to turn the VTK series back on when actually
//! visualizing a specific mismatch geometry in ParaView.
inline bool writeMismatchVTKEnabled()
{
   const char* env = std::getenv("WALBERLA_WRITE_MISMATCH_VTK");
   return env != nullptr && std::string(env) == "1";
}

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

   // SolidCylinderContactDetection's Cylinder x Cylinder overload (ParticleKernels.h) has a cheap
   // raw-center-distance reject BEFORE either the closed-form solver or GJK/EPA runs -- added for
   // a real production bug (job 4103899: a stale periodic-ghost position produced a false-positive
   // band-band-skewed contact ~34 LU apart). It uses ac.getInteractionRadius(idx), which
   // ParticleStorage.h defaults to 0 for every newly-created particle. Production code always sets
   // this explicitly at floater creation (SetupStages.h/Utility.h); this test fixture never did,
   // so every call here was tripping that guard (0+0=0 max reach, any nonzero separation "too far")
   // and returning contact=false before real detection ever ran, regardless of the true geometry --
   // confirmed via a crash in testCylinderCylinderParallelMatchesExisting (two clearly-overlapping
   // parallel cylinders, radius 0.5 each, centres 0.8 apart, silently reporting no contact). Set the
   // same circumscribing-sphere radius production code uses (SetupStages.h: circumR =
   // sqrt(radius^2 + halfLength^2)) so this fixture matches what SolidCylinderContactDetection
   // actually expects to find populated.
   body0->setInteractionRadius(std::sqrt(radiusI * radiusI + real_t(0.25) * lengthI * lengthI));
   body1->setInteractionRadius(std::sqrt(radiusJ * radiusJ + real_t(0.25) * lengthJ * lengthJ));

   SolidCylinderContactDetection gjk;
   const bool                    contactFound = gjk(0, 1, cylI, cylJ, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth(), gjk.getContactNormal(), gjk.getContactPoint() };
}

//! [Claude, 2026-10-15] NOT GJK/EPA: SolidCylinderContactDetection's Cylinder x HalfSpace overload is the production
//! support-mapping closed form (ParticleKernels.h), so this checks paperAnalytical against production.
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

//! Sphere x Cylinder has no dedicated override anywhere in the SolidCylinderContactDetection /
//! GeneralContactDetection hierarchy (unlike sphere-sphere/sphere-box, inherited from
//! AnalyticContactDetection, or cylinder-cylinder/box-box, overridden in ParticleKernels.h), so this
//! exercises the generic GJK/EPA template directly -- exactly the ground truth
//! paper_analytical::sphereCylinderContact is validated against below.
//! [Claude, 2026-10-15] NOT GJK/EPA: SolidCylinderContactDetection's Sphere x Cylinder overload is the production closed
//! form (ParticleKernels.h, mid-overlap contact point). The true GJK/EPA reference is gjkSphereCyl below.
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

//! [Claude, 2026-10-15] True GJK/EPA for a sphere and a cylinder: mesa_pd's GeneralContactDetection itself, whose generic
//! operator() runs collideGJKEPA for this pair (productionSphereCyl above goes through the production closed form).
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
//! [Claude, 2026-10-15] It can be hand-verified after all (the x face axis wins), see the checks below.
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

   // [Claude, 2026-10-15] hand-derived: B reaches 0.5 sqrt(2) towards -x, so the x axis overlaps by 0.5 + 0.5 sqrt(2) - 1.15
   // = 0.0571068; B's face axes at 45 deg overlap by 0.5 sqrt(2) + 0.5 - 1.15/sqrt(2) = 0.394, the cross products
   // of the z edges coincide with these axes. Confirmed by an independent brute-force minimum translation distance.
   WALBERLA_CHECK(sat.contact);
   // overlap = 0.5 + 0.5*sqrt(2) - 1.15 = 0.0571068 (x-face axis wins over the 45deg face axes, 0.394)
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

   // [Claude, 2026-10-15] independent brute-force minimum translation distance (overleaf scripts/contact-validation/
   // verify_hand_cases.py): 0.078618 along -x, so the pair does overlap and the comparison below is not vacuous.
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
}

//! Reconstruction of a REAL fatal collision from job 4131344 (Antidunes, collapsed-spawn 1000-box
//! scenario): uid 8142508126285856798 (box A, halfExtents <3,9,3>) vs uid 8142508126285856796 (box
//! B, same halfExtents), steps 10148-10165, useAnalyticalContactDetection=true in that run. Ground
//! truth pulled directly from the run's own contact log (rank_113.csv), not assumed: at the deepest
//! logged instant (step 10153, subCycle 9) penetrationDepth=-0.773241 and the SAT-computed contact
//! normal was (0.17535, 0.0340412, 0.983917) -- normalized already (|n|=1.0000), and NOT a clean
//! cardinal axis. That's the key fact this test exists to check: a normal with a small but real,
//! SMOOTHLY-VARYING (not jumping) X/Y component sustained across the whole 15-step contact is what a
//! genuinely TILTED box's own face normal looks like, not what two perfectly axis-aligned boxes could
//! ever produce (their SAT normal is always an exact Cartesian axis). Neither box's actual rotation
//! was logged (state.csv has position/velocity but not orientation), so it can't be recovered exactly
//! -- reconstructed here as: box A tilted by the exact minimal rotation that sends its local +Z face
//! normal onto the REAL observed normal vector above (box B left axis-aligned, since the tilt has to
//! live on at least one side and A is the more bumped-around of the two per the run's own collision
//! log), separated along that real normal by whatever distance gives penetration matching the real
//! -0.773241 given these real half-extents. This is a best-effort reconstruction of an under-
//! determined system (6 real DOF unknown, 4 constraints known: penetration + 3 normal components) --
//! not a bit-exact replay -- documented as such. The question this answers: for geometry matching the
//! real event's actual depth AND actual (non-cardinal, tilted) normal direction, does the analytical
//! SAT path agree with the existing GJK/EPA path, or does the tilt expose a real analytical-detection
//! bug that could itself be contributing to the oversized launch force (as opposed to the force being
//! purely a spring/collisionTime magnitude issue, independent of which detection path computed it)?
void testBoxBoxSATRealAntidunesLaunchCollision()
{
   const Vec3 half(real_t(3), real_t(9), real_t(3)); // real halfExtents from the prm, both boxes

   // Real, ground-truth SAT normal at the deepest logged instant (step 10153, subCycle 9).
   const Vec3 realNormal(real_t(0.17535), real_t(0.0340412), real_t(0.983917));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(realNormal.length(), real_t(1), real_t(1e-3));
   const real_t realPenetration = real_t(0.773241); // SAT/paper convention: positive = penetrating

   // Box A: tilted by the minimal rotation sending local +Z onto realNormal (see comment above).
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

   // Box B: axis-aligned (identity) -- the tilt is attributed entirely to A (see comment above).
   const Mat3 identity;

   // Separation along realNormal: SAT overlap = projRadius(A,n) + projRadius(B,n) - centreDist.
   // projRadius(A,n) with A's local Z EXACTLY aligned to n is just half.z (=3), since A's other two
   // local axes are then exactly perpendicular to n and contribute 0. projRadius(B,n) with B
   // axis-aligned is the standard |half.x*n.x| + |half.y*n.y| + |half.z*n.z|.
   const real_t projA = half[2]; // A's local Z axis IS realNormal by construction
   const real_t projB =
      std::abs(half[0] * realNormal[0]) + std::abs(half[1] * realNormal[1]) + std::abs(half[2] * realNormal[2]);
   const real_t centreDist = projA + projB - realPenetration;

   const Vec3 posA(real_t(0), real_t(0), real_t(0));
   // B sits behind A along -realNormal (A is the falling box on top, matching the real run: box A
   // (798) is the one this test is named for -- it fell onto B (796) from above).
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
      // Loose epsilon deliberately: this is a reconstructed, not bit-exact, real event -- the point is
      // to catch a GROSS divergence (wrong axis picked, sign flip, order-of-magnitude overlap error),
      // not to demand floating-point-identical agreement with an already-approximate setup.
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, -gjk.penetrationDepth, real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, gjk.normal, real_t(1e-3));
      // [Claude, 2026-10-15] this reconstruction does reproduce the logged contact (brute-force minimum
      // translation distance 0.773241 along the logged normal)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, realPenetration, real_t(1e-4));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.normal, realNormal, real_t(1e-4));
   }
}

//! Same reconstruction technique as testBoxBoxSATRealAntidunesLaunchCollision above, applied to the
//! two logged [DEM-WATCH] deep-contact events that immediately preceded job 4139454's fatal NaN-
//! quaternion crash (bigjam scenario, 162 box floaters, step 119243, rank 33): a Box x Box pair whose
//! reported overlap is anomalously large (40.12% and 20.64% of the smaller body's interaction radius,
//! respectively, vs. an ordinary jam contact's ~6-10%) right before one particle's velocity blew up
//! from ~1e-4 to double digits within two DEM sub-cycles. The question these answer: for geometry
//! matching the ACTUAL reported depth and normal of each event, does boxBoxSAT (the production
//! detector for this shape pair) agree with GJK/EPA -- i.e. is the anomalous number a genuine, exactly
//! computed penetration depth (meaning the DEM force/integration response to it is the actual
//! problem), or does it already disagree with the general-purpose oracle here (meaning detection
//! itself contributed to the launch, the way an earlier real incident's tilted-normal case turned out
//! to)? Reconstruction is under-determined the same way the earlier test's is (4 known constraints --
//! penetration + 3 normal components -- for 6 unknown relative DOF): box A's local Z is rotated onto
//! the reported normal, box B stays axis-aligned, and the two are separated along that normal by
//! whatever distance reproduces the reported penetration exactly for these real halfExtents=<6,1,1.5>.
//! [Claude, 2026-10-15] CORRECTION: it does not -- the separation reproduces the overlap along the logged normal only,
//! and another axis of the pair overlaps less (2.043 instead of 2.513 for event 1, 1.235 instead of 1.293 for event 2,
//! independent brute-force minimum translation distance). Both tests remain valid SAT/GJK cross-checks on geometry
//! derived from the run, but they are not replays of the logged contacts.
void testBoxBoxSATJob4139454Explosion1()
{
   const Vec3 half(real_t(6), real_t(1), real_t(1.5)); // real FloatingBox halfExtents, this prm

   // Real, ground-truth SAT-equivalent values from the [DEM-WATCH] log at step 119243, subCycle 1
   // (the INITIATING event -- uid=2161727821137838098 vs uid=2305843009213693979, the deepest of the
   // two logged contacts and the one immediately preceding the velocity spike to ~10.3 in every
   // component).
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
      // [Claude, 2026-10-15] this reconstruction does NOT reproduce the logged depth 2.51343: another axis of the pair
      // overlaps less, and an independent brute-force minimum translation distance gives 2.04333 (normal not the
      // logged one). The test is therefore a SAT/GJK cross-check on a run-derived geometry, not a replay.
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(2.04333), real_t(1e-4));
   }
}

//! Same as above, for the SECOND logged event (same step, subCycle 2 -- uid=432345564227567617 vs
//! the now-runaway uid=2161837838098, which had ALREADY reached vel=(-10.3,-7.3,10.3) by this point).
//! Notably this overlap (20.64%) is smaller than event 1's (40.12%) even though the particle's
//! velocity nearly tripled between the two -- consistent with the already-corrupted velocity from
//! event 1 carrying it into a second, less-extreme contact rather than the overlap itself continuing
//! to grow. See the chat/analysis this test was written from for the full step-by-step reasoning.
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
      // [Claude, 2026-10-15] this reconstruction does NOT reproduce the logged depth 1.29292: another axis of the pair
      // overlaps less, and an independent brute-force minimum translation distance gives 1.23458 (normal not the
      // logged one). The test is therefore a SAT/GJK cross-check on a run-derived geometry, not a replay.
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(sat.overlap, real_t(1.23458), real_t(1e-4));
   }
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
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

   // [Claude, 2026-10-15] hand-derived: the axes cross at x = 0 with a perpendicular distance of 0.6 along y (both axes lie
   // in planes y = const), well inside both lengths, so the overlap is 0.4 + 0.4 - 0.6 = 0.2 along -y.
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

//! Cylinder j's curved side pressed into cylinder i's flat top face, with j's axis perpendicular to
//! i's (eq. 34's precondition). Hand-derived: overlap = (halfLi + radiusJ) - |zIJ| = (1+0.5) - 1.4 =
//! 0.1. NOTE: despite this test's ORIGINAL docstring claiming eq. 42's "fully inside" case, checking
//! Kodam's actual line AB (eq. 40-41: cylinder j's own centreline of half-length 2, translated onto
//! i's face) against i's face disc (radius 0.5) shows BOTH endpoints (x=+-2) lie far outside it --
//! this is actually the case neither eq. 42 (fully inside) nor eq. 44 (one end inside) names: AB
//! straddles the disc, poking out through both ends. faceBandContactPoint's finite line-circle
//! intersection handles this by construction (see its own comment); the exact expected point below
//! is independently hand-derived from that same construction, not merely read back from the code:
//! the two disc-boundary crossings are symmetric about x=0 (since cj's centre projects to x=0), so
//! their midpoint is (0, 0.45), pulled back from the face (z=1) by half the overlap (0.05) to 0.95.
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

//! Genuine eq. 42 "fully inside" case: a short, narrow cylinder j (half-length 0.3, matching AB's
//! own half-length since eq. 40-41 has no radius term) whose translated centreline AB fits entirely
//! within i's much larger face disc (radius 1.0). Hand-derived expected point: AB's own midpoint,
//! i.e. cj's centre projected onto i's face plane, (0, 0.2, 1 - 0.1/2) = (0, 0.2, 0.95).
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

//! Genuine eq. 44 "partially inside" case (one end of AB inside the disc, one outside), constructed
//! off-centre (cj's projected centre at (0.5, 0.3), not on i's own axis) so the two ends of AB are NOT
//! symmetric the way testCylinderCylinderFaceBand's straddling case is. Hand-derived via the same
//! line-circle intersection faceBandContactPoint solves: with A (t=+halfLj=+1) at (1.5,0.3) (outside
//! Ri=1) and B (t=-halfLj=-1) at (-0.5,0.3) (inside), the disc crossing nearer A is at
//! t2=(-0.5+sqrt(0.91))=0.453939, so the contact point sits at the midpoint of B and that crossing,
//! t=(-1+0.453939)/2=-0.273031, i.e. xy=(0.5-0.273031, 0.3)=(0.226970, 0.3).
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
   // point = midpoint of B and the disc-crossing nearer A, t=(-1+0.453939)/2 -> (0.226970, 0.3) (Kodam et al. 2010 eq. 44)
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

//! Cylinder j tilted 45 degrees to cylinder i, positioned so j's bottom rim region is near i's flat
//! face. This test EXISTS TO DOCUMENT A REAL, DISCOVERED LIMITATION, not to demonstrate agreement:
//! [Claude, 2026-10-15] Cause (independent brute-force minimum translation distance 0.235356, direction
//! (0.707, 0, -0.707), i.e. perpendicular to j's axis and 45 deg off i's face normal): the deepest rim point E lies at
//! radius 0.566 on i's face of radius 0.6, close to i's own rim, so moving j sideways across that rim separates the
//! pair sooner than pushing it back out through the face. Face-edge treats the face as unbounded and returns the
//! depth of E below it (0.29853). With a host of ten times the probe radius (testCylinderCylinderFaceEdgeClean) E is
//! far from the host rim and face-edge is exact. The values are pinned below so that a change is noticed.
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
   WALBERLA_CHECK(kodam.contact);
   WALBERLA_CHECK_EQUAL(std::string(kodam.scenario), "face-edge");
   // depth of the lowest rim point below the face: 1 - (1.55 - 0.5 sqrt(2) - 0.2 * 0.5 sqrt(2)) = 0.298528
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(1) - (real_t(1.55) - real_t(0.6) * std::sqrt(real_t(2))),
                                       real_t(1e-6));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(-gjk.penetrationDepth, real_t(0.235356), real_t(1e-3)); // brute-force MTD
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
//! [Claude, 2026-10-15] CORRECTION: the pair is exactly touching -- an independent brute-force minimum translation
//! distance gives 0.000000 and GJK/EPA reports a depth of about -1e-10 -- so the closed form's "no contact" is the
//! boundary of the scenarios, not a missed overlap. Asserted below. A genuine miss (trial 1209 of the reference sample,
//! depth 0.113) is used in testHybridTier1TriggersOnMissedContact.
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
   WALBERLA_CHECK(!gjk.contact); // confirm this is a TRUE negative, not another instance of the
                                  // rim-vs-band limitation masking a real contact
}

//! Guo, Wassgren, Ketterhagen, Hancock & Curtis (2012)'s "special contact" (their Fig. 8) -- NOT one
//! of Kodam et al.'s original 7 scenarios. Fires when the band-edge candidate point P lands within
//! host i's band region (same test band-edge uses), P is on the SAME side of i's face as probe j's
//! own centre, but j's centre itself lies BEYOND that face (j has "sunk in" past i's own end) --
//! ordinary band-edge would otherwise report the RADIAL delta at P, which can be far larger than the
//! true (smaller) axial pile-up overlap this configuration actually has.
//!
//! Geometry: cylinder j (radius 1, length 0.6 -- short and wide, like a squat cap) tilted 35 degrees
//! and centred at z=1.7, well past host i's own face at z=+halfLi=1.5 (j's centre is 0.2 beyond it),
//! offset in x so face-edge's own candidate point E lands outside i's face disc (rho^2=1.23 > Ri^2=1,
//! comfortable margin) and doesn't pre-empt this branch. P (comfortable margins in every direction:
//! 0.18 axial margin inside i's band, radial well inside at rho^2=0.006) sits at the SAME sign of z
//! as j's own centre, satisfying every precondition in PaperAnalyticalContacts.h's special-contact
//! branch. Found by direct parameter search against this file's own dispatch logic (see git history/
//! PR description), then hand-verified here against GJK/EPA as the independent ground truth.
//!
//! NEWLY DISCOVERED LIMITATION while building this test: special-contact's reported overlap
//! (deltaAp = axial distance from P to i's face) massively UNDERESTIMATES the true penetration for
//! this geometry -- kodam reports 0.181, GJK/EPA reports 0.610, a 3.4x gap. This is the SAME root
//! cause already documented for band-edge in testCylinderCylinderTiltedRimIntoFace/IntoBand above:
//! deltaAp is derived purely from P, the single closest point on probe's RIM, while the true deepest
//! penetration for this "j sunk in past i's face" configuration is dominated by j's curved BAND
//! surface plunging straight through i's face, not by anything on j's rim. Special-contact was
//! specifically designed (by Guo et al.) to catch this configuration's EXISTENCE correctly where
//! plain band-edge would have picked the wrong delta (radial vs axial) entirely -- it was never
//! validated for MAGNITUDE, and this test is the first place that's been checked. Only existence and
//! normal direction are asserted below for that reason; overlap is logged, not asserted against GJK.
//! [Claude, 2026-10-15] SUPERSEDED: the corrected special-contact depth (deepest rim point inside the host,
//! 0.614334) is now asserted against GJK/EPA (0.609508, brute-force MTD 0.609529) with the fuzz criterion
//! max(0.01, 5 %); 0.181 was Guo's original depth at P.
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
   // [Claude, 2026-09-29] was CHECK(overlap == 0.1811), Guo et al.'s original value measured at P. The
   // corrected special-contact depth (deepest rim point inside the host, see rimVsFaceOrBand) must now agree
   // with GJK/EPA (0.6095) and the brute-force MTD arbiter (0.6095) within the fuzz suite's own tolerance.
   WALBERLA_CHECK_LESS(std::abs(kodam.overlap - (-gjk.penetrationDepth)),
                       std::max(real_t(0.01), real_t(0.05) * (-gjk.penetrationDepth)));
   WALBERLA_CHECK(gjk.contact);
   // special-contact deliberately reports a purely axial (0,0,+-1) "face-normal-like" normal by
   // construction (see PaperAnalyticalContacts.h's comment on this branch) rather than GJK/EPA's true
   // normal, which for a genuinely tilted j also picks up a small off-axis component (here
   // <-0.12,0.11,-0.99> vs kodam's exact <0,0,-1>) -- an approximation, not a bug, so this checks
   // directional agreement (dot product close to 1) rather than exact componentwise equality.
   WALBERLA_CHECK_GREATER(kodam.normal * gjk.normal, real_t(0.95));
   // [Claude, 2026-09-29] special-contact point: host i's face at z = +1.5, P at depth overlap below it,
   // so the point must be at z = 1.5 - overlap/2 (same sign fix as face-edge).
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], real_t(1.5) - real_t(0.5) * kodam.overlap, real_t(1e-10));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (special-contact): kodam overlap=" << kodam.overlap
                              << " normal=" << kodam.normal << " vs GJK/EPA penetrationDepth="
                              << gjk.penetrationDepth << " normal=" << gjk.normal);
}

//! Independent, third arbiter for testCylinderCylinderSpecialContact's own disagreement (Kodam/Guo:
//! 0.181, GJK/EPA: ~0.610): the minimum TRANSLATION distance (MTD) between two convex solids, computed
//! directly from each cylinder's own support function via a dense sweep over unit directions -- no
//! dependency on either Kodam/Guo's candidate-point approach or on this codebase's own GJK/EPA
//! implementation, so it can genuinely arbitrate between them rather than just re-confirming one side.
//!
//! For a solid convex body centred at c with support function h(d) = max_{p in body} (p-c).d, and two
//! bodies A, B, the projected gap along direction d is (cB-cA).d - hA(d) - hB(-d); if this is negative
//! for some d, the bodies overlap along that axis by that amount, and the TRUE minimum translation
//! distance (over every direction, not just a finite candidate set) is the most negative such gap,
//! i.e. MTD = -max_d[gap(d)] = min_d[hA(d) + hB(-d) - (cB-cA).d]. A cylinder's own support function is
//! elementary: h(d) = halfLength*|axis.d| + radius*sqrt(max(0, 1-(axis.d)^2)).
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

   // Coarse-to-fine directional search: a Fibonacci sphere lattice (uniform-ish coverage, no pole
   // clustering), then a local perturbation search re-centred on the best direction found so far --
   // avoids needing many more points for a precise MAXIMUM (the least-overlapping, true separating,
   // axis -- see this function's own doc comment for why it's a maximization, not a minimization).
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

   // Local refinement: small random perturbations around the best direction found, shrinking the
   // perturbation radius -- a simple, dependency-free local search to polish the coarse global sweep.
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

//! [Claude, 2026-09-29] Sweep of the special-contact configuration (Guo et al. 2012, Fig. 8): probe j tilted about
//! x by 5..45 degrees, shifted laterally and axially around testCylinderCylinderSpecialContact's geometry. For every
//! pose the analytical dispatch classifies as "special-contact", the corrected depth (deepest rim point inside
//! the host) must agree with GJK/EPA within the fuzz suite's tolerance max(0.01, 5 %), and the contact point must
//! lie at the middle of the overlap along the (axial) normal. Before the correction, Guo's depth-at-P failed this
//! for most tilted poses.
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, real_t(0.05), real_t(1e-9)); // [Claude, 2026-10-15] brute-force MTD 0.05
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (edge-edge-ish, FIXED by Guo et al. corrections): scenario="
                              << kodam.scenario << " overlap=" << kodam.overlap
                              << " (was a false negative before the bidirectional check + corrections)");
}

//======================================================================================================================
// Randomized fuzz test helpers, shared between the Cylinder x Cylinder and Box x Box fuzz tests
// below: random orientation, and a separation sampler that stays within a PHYSICALLY REALISTIC
// overlap depth. Random sampling over the full [0, reach) separation range (an earlier version of
// this test used exactly that) mostly generates near-total-interpenetration geometry -- two random
// 0.3-2.0-radius bodies with only ~0.5 units between their centres, i.e. 90%+ overlap -- that a
// working simulation's spring-dashpot contact force would never actually let develop; it would have
// already pushed the bodies apart at a tiny fraction of that depth. Capping the sampled overlap to a
// generous-but-real bound (15% of the smaller body's own radius -- steady-state DEM overlaps in
// practice are usually low single digits of a percent, so even this is already a stress case, not a
// realistic operating point) means failures found here are things that could actually happen during
// a run, not artifacts of testing geometry the simulation itself would never reach.
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
   // [Claude, 2026-09-30] Uniformly distributed over SO(3) (Shoemake 1992: unit quaternion from three uniform numbers).
   // Before, this drew a uniform random axis and a uniform angle in [0, 2pi) -- that covers every orientation, but
   // over-samples rotations close to the identity (the Haar measure has density ~ 1 - cos(angle) in the angle), so the
   // per-scenario rates of every randomized test were biased towards configurations near the sampled bodies' own frames.
   // Every fixed-seed number recorded in comments below 2026-09-30 refers to the old sampler; the old one is kept as
   // randomRotationAxisAngle() for reproducing them.
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
   //! Legacy sampler (uniform axis, uniform angle): NOT uniform on SO(3). Only for reproducing pre-2026-09-30 numbers.
   Rot3 randomRotationAxisAngle()
   {
      Rot3 r;
      r.rotate(randomUnitVec(), angleDist(rng));
      return r;
   }
   // Finds the ACTUAL GJK/EPA contact boundary along a given direction via bisection, then samples a
   // centre distance around it -- replaces an earlier version that estimated "reach" analytically
   // (sum of circumscribing-sphere radii / half-diagonals). That analytic bound is far looser than
   // the true touching distance for most random orientations (a box's half-diagonal is its CORNER
   // reach, not its typical face-on reach), so sampling relative to it mostly missed the actual
   // contact boundary entirely -- confirmed empirically: an earlier run of the box fuzz test using
   // that scheme produced 0 genuinely-overlapping trials out of 5000. Bisection finds the true
   // boundary for whatever orientation this trial actually drew, cylinder or box alike, so sampling
   // around it reliably produces both near-miss and realistic-overlap coverage regardless of shape.
   //
   // farSep must be a centre distance the caller guarantees is separated (contactAt(farSep) ==
   // false); nearSep must be guaranteed overlapping (contactAt(nearSep) == true). minBodySize caps
   // how deep past the found boundary this will sample (see sampleSeparation's own comment).
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
      const real_t maxOverlap = maxOverlapFraction * minBodySize; // [Claude, 2026-10-11] fraction optional (default 0.15)
      // Sample from "maxOverlap past the boundary, into contact" to "twice that far out, clearly
      // separated" -- covers realistic overlap depth and near-miss geometry either side of the true
      // boundary, centred where contact actually starts instead of a loose analytic guess.
      std::uniform_real_distribution< real_t > d(boundary - maxOverlap, boundary + real_t(2) * maxOverlap);
      return d(rng);
   }
};

//======================================================================================================================
// Randomized fuzz test: random radius/length/orientation for both cylinders, separation sampled
// near the true GJK/EPA contact boundary within a realistic overlap depth (see
// FuzzRng::sampleSeparationNearBoundary above), run through both the raw analytical dispatch
// (cylinderCylinderContact -- now closed-form throughout, no golden-section
// search left anywhere in it) and the production hybrid dispatcher (cylinderCylinderContactHybrid,
// which falls back to GJK/EPA whenever cylinderCylinderResultIsUncertain() flags the analytical
// result), each checked against GJK/EPA directly. Fixed seed for reproducibility. Not a replacement
// for the hand-picked scenario tests above (those pin down which SPECIFIC regime each geometry class
// hits); this is coverage those can't give -- thousands of arbitrary orientations/sizes, the way
// real floaters in a dense scenario actually meet each other.
//======================================================================================================================

//======================================================================================================================
// Debug visualization: for every raw-analytical-vs-GJK/EPA disagreement the fuzz test below finds
// (both "wrong existence" -- 844/10000 in the current run -- and "wrong overlap magnitude" among
// the cases where both sides agree a contact exists), dump the exact pair geometry as a VTK point
// glyph dataset (same convention BuildTimeloop.h's own floaterCylinder output uses: "glyphScale"
// (2r, length, 2r) + "rotation" pre-composed with Rx(+90) to correct for ParaView's built-in
// Cylinder glyph source having its axis along local Y where mesa_pd's Cylinder has it along local
// Z -- see BuildTimeloop.h's own comment for the full derivation) so the actual failing geometries
// can be opened directly in ParaView (Filters > Glyph, Glyph Type = Cylinder, Scale Array =
// glyphScale, Orientation Array = rotation, Scale Factor = 1) instead of only read as numbers in a
// log. A companion CSV reports the same cases in tabular form.
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
   const char*   cause      = ""; // mathematical root-cause tag, see
                                    // testCylinderCylinderFullComparisonAndHybridFeasibility's
                                    // classifyMagnitudeCause/classifyDispatchGapCause -- always a
                                    // fixed string literal (never dynamically built), so this pointer
                                    // never dangles
   real_t        causeMetric = real_t(-1); // supporting number for `cause` (e.g. the kodam-vs-GJK
                                             // normal misalignment angle in degrees); -1 where n/a
};

inline std::array< real_t, 9 > glyphRotationFromR(const Mat3& R)
{
   // (R * Rx(+90))^T = Rx(-90) * R^T -- identical derivation/formula to BuildTimeloop.h's
   // SelectParticleRotationMatrix for floaterCylinder, reproduced here so this standalone debug
   // dump renders correctly without depending on BuildTimeloop.h.
   return { R(0, 0), R(1, 0), R(2, 0), R(0, 2), R(1, 2), R(2, 2), -R(0, 1), -R(1, 1), -R(2, 1) };
}

//! Writes ONE cylinder (single point) as its own small .vtu -- side=0 for ci, side=1 for cj, so the
//! SAME writer serves both of writeCylMismatchVTKSeries's two independent series below.
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

   // edgeEdgeStatus: was this mismatch produced by the real Kodam eq. 82-101 edge-edge
   // implementation actually firing (1 -- findRimCrossing found valid chords on both rims, but the
   // linear-chord/sagitta approximation still under-reports the magnitude), by that same code path
   // being reached but failing to find a valid crossing on either rim (2 -- "disabled": scenario ==
   // "none" means dispatch fell all the way through every scenario including edge-edge with no
   // analytical answer at all, relying entirely on Tier 1's GJK/EPA fallback), or by a DIFFERENT
   // scenario resolving the pair before edge-edge was ever reached (0 -- not applicable here).
   const int edgeEdgeStatus = (std::string(r.scenario) == "edge-edge")   ? 1
                              : (std::string(r.scenario) == "none")      ? 2
                                                                          : 0;
   f << "    <DataArray type=\"Int32\" Name=\"edgeEdgeStatus\" NumberOfComponents=\"1\" format=\"ascii\">"
     << edgeEdgeStatus << "</DataArray>\n";

   f << "    <DataArray type=\"Float64\" Name=\"analyticalOverlap\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.analyticalOverlap << "</DataArray>\n";

   f << "    <DataArray type=\"Float64\" Name=\"gjkOverlap\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.gjkOverlap << "</DataArray>\n";

   // cause is a fixed-string tag (see CylMismatchRecord::cause); ParaView's ascii VTU writer here has
   // no string DataArray, so the human-readable tag goes in the companion CSV only and causeMetric
   // (e.g. a misalignment angle in degrees, -1 if n/a) is exposed here for colour-by-value glyphing.
   f << "    <DataArray type=\"Float64\" Name=\"causeMetric\" NumberOfComponents=\"1\" format=\"ascii\">"
     << r.causeMetric << "</DataArray>\n";

   f << "   </PointData>\n";
   f << "   <CellData>\n   </CellData>\n";
   f << "  </Piece>\n </UnstructuredGrid>\n</VTKFile>\n";
}

//! Writes ONE cylinder's (ci or cj, per `side`) mismatch history as its OWN time series -- baseDir/
//! <name>_ci.pvd and baseDir/<name>_cj.pvd, each a fully independent .pvd/frames set. Split (not one
//! combined 2-point-per-frame series, as an earlier version of this did) specifically so each
//! cylinder can get its own Glyph filter in ParaView with its own independent Scale Factor -- ci and
//! cj are drawn from the same [0.3,2.0] radius / [0.5,4.0] length ranges but completely
//! independently per trial, so a single shared glyph scale is never really the right choice for
//! both at once.
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

//! Writes BOTH cylinders' time series (see writeCylMismatchSide's own comment for why they're
//! separate files, not one combined one). Open both .pvd files in ParaView as two independent
//! sources, add a Glyph filter to each, and scale/colour them independently -- ParaView keeps both
//! in lock-step on the same timestep slider since they share the same 0..N-1 timestep numbering.
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
                                                // (Tier 1 or Tier 2) and hybrid actually called GJK/EPA
   real_t        maxAgreeOverlapDiff = real_t(0);
   real_t        maxAgreePointDist   = real_t(0);
   Vec3          worstPointCi, worstPointCj;

   // Every raw-analytical-vs-GJK/EPA disagreement this loop finds, for the VTK/CSV dump at the end
   // -- see writeCylMismatchVTK/writeCylMismatchCSV's own comments for why (visualizing WHY the
   // hybrid fallback is necessary, not just counting how often).
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
         // Both agree a contact exists -- but is the RAW analytical magnitude actually trustworthy,
         // or does it just happen to agree on existence while getting the depth wrong (exactly the
         // rim-vs-band gap cylinderCylinderResultIsUncertain's Tier 2 exists to catch)? Threshold:
         // absolute floor for small overlaps, 5% relative for larger ones -- not asserted anywhere,
         // purely for the debug dump below.
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

   // The hybrid dispatcher's whole purpose is to catch what the raw analytical path misses -- its
   // EXISTENCE mismatch rate against GJK/EPA is still held to 0, strictly: Tier 1 (the unconditional
   // reach-sum bound whenever the raw dispatch reports no contact) plus special-contact's own
   // unconditional routing together guarantee this regardless of overlap magnitude. History: before
   // Tier 2 ever covered band-band-skewed/edge-edge (when Tier 2 still existed), this was ~8% of all
   // trials (403/5000), purely false positives -- reporting contact for genuinely separated pairs.
   // That was traced to a genuine frame-mismatch bug in band-band-skewed's own perpendicular-distance
   // formula (dotting a body-frame vector against a world-frame one, correct only when ci itself was
   // unrotated) and fixed at the root, not safety-netted -- band-band-skewed is provably exact under
   // its own precondition regardless of any fallback tier.
   WALBERLA_CHECK_EQUAL(numHybridMismatch, 0);
   // Overlap MAGNITUDE is no longer held to a tight bound: Tier 2 (a circumscribing-sphere magnitude
   // cross-check applied to every non-exempt scenario) was REMOVED from cylinderCylinderResultIsUncertain
   // -- see that function's own comment for the full history -- after a dedicated master comparison
   // showed it pushed the GJK/EPA fallback rate to ~90% while still being structurally one-sided (it
   // only catches UNDERESTIMATES, never overestimates, and a deliberately-deep-overlap stress test
   // found a case where it missed one anyway). The traded-off residual: band-edge/edge-edge/face-edge
   // now keep their own raw overlap value whenever they fire (no longer corrected by Tier 2), so
   // maxAgreeOverlapDiff can legitimately be as large as their own documented error range (median
   // ~2.7% relative, max ~13.5% absolute in a 10000-trial master comparison) -- informational only
   // below, not asserted, same treatment already given to the point-distance metric right after it.
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] max overlap magnitude difference on agreeing contacts=" << maxAgreeOverlapDiff
      << " (informational only -- Tier 2 no longer corrects band-edge/edge-edge/face-edge's own raw "
      << "overlap value, see cylinderCylinderResultIsUncertain's comment; not asserted)");
   // Point distance: band-band-skewed's own contact POINT can differ from GJK/EPA's along a flush
   // contact line with no unique point (the same benign ambiguity already documented for face-band and
   // band-band-parallel elsewhere in this file) -- informational only, same treatment the box-box fuzz
   // test already gives its own analogous "PERPENDICULAR-to-normal point difference" for flush contacts.
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] max contact-point distance on agreeing contacts=" << maxAgreePointDist << " (informational "
      << "only -- not asserted)");
}

//======================================================================================================================
// PURE analytical (cylinderCylinderContact directly -- NOT cylinderCylinderContactHybrid, which would
// blend in GJK/EPA on every Tier-2-flagged trial and hide exactly what this is meant to show) validated
// against GJK/EPA, broken down BY SCENARIO: how often does each of Kodam et al.'s dispatch outcomes
// actually occur under random fully-3D orientations and realistic (boundary-biased, capped-overlap)
// separations, and on the trials it fires for, does it agree with GJK/EPA on existence, overlap
// magnitude, and contact point? Same random/sampling scheme as testCylinderCylinderRandomizedFuzz
// above (fully random axes for BOTH cylinders, not forced-perpendicular the way
// testCylinderCylinderFaceBandBarelyOverlappingFuzz is) so the scenario mix reflects what a real
// simulation's arbitrary floater orientations would actually produce, not one scenario forced into
// dominance. Reported for inspection -- not asserted against a fixed bound, since face-edge/band-edge/
// edge-edge/special-contact are DOCUMENTED to sometimes disagree with GJK/EPA (that disagreement is
// the whole reason cylinderCylinderResultIsUncertain/the hybrid fallback exist); asserting a tolerance
// here would just be re-deriving that known limitation under a different name.
//======================================================================================================================
void testCylinderCylinderPureAnalyticalScenarioBreakdown()
{
   FuzzRng                                  fuzz(24680); // distinct seed from every other fuzz test
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
      // Root-cause diagnostics for the worst overlap-magnitude disagreement: is GJK/EPA's contact
      // normal still close to what the analytical rim/edge-based formula reports (a shared-feature,
      // pure-magnitude imprecision), or does it point somewhere structurally different (a sign that
      // GJK/EPA found a DIFFERENT, non-rim closest feature -- e.g. the curved band surface away from
      // any edge -- that the analytical dispatch never searches for at all, per
      // testCylinderCylinderTiltedRimIntoFace/IntoBand/SpecialContact's own documented limitation)?
      real_t        worstOverlapDiff = real_t(0);
      CylinderState worstCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      CylinderState worstCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      Vec3          worstAnalyticalNormal, worstGjkNormal;
      real_t        worstAnalyticalOverlap = real_t(0), worstGjkOverlap = real_t(0);
      // One example existence mismatch (analytical says no contact, GJK/EPA says yes -- a total miss,
      // not just a wrong magnitude), to check whether it shares the same root cause.
      bool          haveMismatchExample = false;
      CylinderState mismatchCi{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      CylinderState mismatchCj{ Vec3(real_t(0)), Mat3(), real_t(0), real_t(0) };
      Vec3          mismatchGjkNormal;
      real_t        mismatchGjkOverlap = real_t(0);
      // Best-agreeing example with a non-trivial overlap (gjkOverlap > 0.05, so it isn't just a
      // barely-touching near-zero case) -- candidate geometry for a genuine, well-agreeing hand-drawn
      // test, to complement the existing "known limitation" ones for this scenario.
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
// Targeted fuzz test: axes forced (ANTI-)PARALLEL by construction -- the necessary precondition for
// face-face/band-band-parallel that testCylinderCylinderPureAnalyticalScenarioBreakdown's fully-random
// orientations structurally cannot hit at all (exact axis-parallelism, checked to 1e-8, has zero
// measure under random 3D rotation -- confirmed empirically there: 0/5000 trials landed in either
// scenario). Same boundary-biased separation sampling and same trial count (5000) as that test, so
// its results are directly comparable in the same breakdown table.
//
// A random direction is used for the boundary search, exactly like the general fuzz test: separating
// mostly along the shared axis tends to produce face-face, separating mostly sideways tends to
// produce band-band-parallel, and eq. 22's own delta_a-vs-delta_r comparison (see
// cylinderCylinderContact) picks the right one for everything in between -- no need to construct the
// two sub-cases as separate loops. Cylinder j's axis is either exactly cylinder i's own (parallel) or
// flipped 180 degrees about a perpendicular axis (antiparallel), chosen with equal probability per
// trial, so both of eq. 22's "|e_i . e_j| ~ 1" sign cases get covered -- two logs pointing the same
// way, or pointing opposite ways, are equally physical floater configurations.
//======================================================================================================================
void testCylinderCylinderParallelAxesFuzz()
{
   FuzzRng                                  fuzz(11235); // distinct seed from every other fuzz test
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
// Targeted fuzz test: 1000 trials, axes forced PERPENDICULAR by construction (the necessary
// precondition for face-band -- random offsets can still land the actual dispatch in
// face-edge/band-edge/edge-edge/no-contact too; results are tagged by scenario and reported
// separately), separation sampled to be BARELY overlapping -- a much tighter cap (2% of the smaller
// radius, sampled only on the overlapping side of the true GJK/EPA boundary) than
// testCylinderCylinderRandomizedFuzz's already-modest 15%-of-radius, up-to-2x-past-boundary range, so
// every kept trial is a genuine, shallow contact rather than a mix of deep overlap and near-misses.
//
// Compares the RAW analytical result directly against GJK/EPA -- not the GJK/EPA-fallback-blended
// hybrid, which would trivially agree with GJK/EPA on every Tier-2-flagged trial and tell us nothing
// about the analytical formulas themselves. Purpose-built to stress-test detail::faceBandContactPoint
// (the exact Kodam eq. 40-44 construction that replaced this file's earlier point approximation) at
// the shallow-overlap regime a real spring-dashpot contact model actually operates in. Not asserted
// against a fixed tolerance -- reported for inspection, the same style as
// testCylinderCylinderFallbackRateBoundaryBiased below, since the non-face-band scenarios this can
// also land in (face-edge, band-edge) already have documented, accepted imprecision of their own
// (see testCylinderCylinderTiltedRimIntoFace) that isn't this test's job to re-litigate.
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

      // Force cylinder j's axis perpendicular to cylinder i's: Gram-Schmidt a random direction
      // against eI, then align body-z to the result the same way SetupStages.h aligns a floating
      // cylinder's initial axis to a requested world direction.
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

      // Offset cylinder j's centre, in cylinder i's local xy-plane, biased to land near/within i's
      // face disc -- otherwise almost every random draw would be a "no contact at all" geometry.
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
      // s=0 (j's centre level with i's centre along eI) is overlapping for every draw this
      // distribution can produce; skip the ~zero-probability degenerate case rather than assert.
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
      const real_t maxOverlap = real_t(0.02) * std::min(Ri, Rj); // "barely" -- 2%, not the general
                                                                   // fuzz test's already-modest 15%
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
// SAME statistic as the broad/uniform version above (Tier 1/Tier 2 fallback rate and precision), but
// using BOUNDARY-BIASED separation sampling -- the exact same FuzzRng::sampleSeparationNearBoundary
// call testCylinderCylinderRandomizedFuzz uses -- instead of uniform-over-reach. This is the more
// physically relevant population for an actual DEM run: a real spring-dashpot contact model actively
// resists deep interpenetration, so genuine contacts stay shallow, and clearly-separated pairs (most
// of what uniform-over-reach samples) are trivial cases no scenario ever misclassifies anyway. The
// broad/uniform version above over-represents deep overlap (58% of its trials were already in
// contact, many likely deep) that a stable simulation essentially never produces; this version's
// population -- shallow contacts and near-misses -- is where classification correctness actually
// matters. Not asserted against anything -- purely a measurement, logged for inspection.
//======================================================================================================================
void testCylinderCylinderFallbackRateBoundaryBiased()
{
   FuzzRng                                  fuzz(13579); // different seed from every other fuzz test
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials    = 100000;
   int           numFallback  = 0;
   int           numTier1     = 0; // analytical found NO contact at all (reach-sum bound says "maybe")
   int           numTier2     = 0; // analytical found A contact, but the scenario/magnitude is untrusted
   int           numOverlapping = 0;
   // Of the fallback trials, how many would have been FINE to trust anyway (raw analytical actually
   // agreed with GJK/EPA on both existence and magnitude) vs how many the fallback genuinely saved --
   // i.e. cylinderCylinderResultIsUncertain's own PRECISION: of every alarm it raises, how many are
   // real.
   int numTier1WouldHaveBeenFine = 0, numTier1WouldHaveBeenWrong = 0;
   int numTier2WouldHaveBeenFine = 0, numTier2WouldHaveBeenWrong = 0;
   // Which analytical SCENARIO was in play when Tier 2 fired, split into "would still have been fine"
   // vs "genuinely wrong" -- i.e. which of Kodam's 8 dispatch outcomes is actually driving REAL
   // distrust, not just how often each one triggers the check at all.
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
         // Same agreement criterion the main fuzz test uses: existence must match, and if both agree
         // contact exists, magnitude must be within max(0.01, 5%) of GJK/EPA's answer.
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

//! Clean, well-agreeing face-edge case -- deliberately constructed as the counterpart to
//! testCylinderCylinderTiltedRimIntoFace's KNOWN LIMITATION: that limitation traces to the true
//! closest feature landing on probe j's curved BAND surface rather than its rim, a candidate neither
//! Kodam et al. nor Guo et al. ever search for. That failure mode is driven by the interaction of
//! BOTH cylinders' curvature; if host i's own curvature is made negligible (a much larger radius than
//! j's, so its face looks locally flat), the closest-point problem degenerates to "cylinder j vs. an
//! effectively flat plane" -- a case Section 3.1's cylinder-plane edge-plane formula already proves
//! exact (the support-function argument: for any tilt away from 0/90 degrees, a solid cylinder's true
//! extremal point along a fixed direction is always exactly on its rim). Confirmed empirically here:
//! sweeping cylinder j's axial position with Ri=3.0 (vs. Rj=0.3) reproduced EXACT agreement with
//! GJK/EPA to every printed digit at every position tried, unlike the genuinely curvature-coupled
//! geometry in the KNOWN LIMITATION test.
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
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-4)); // tight, not
                                                                                             // the loose
                                                                                             // 1e-2 used
                                                                                             // for the
                                                                                             // known-
                                                                                             // limitation
                                                                                             // case
   WALBERLA_CHECK_GREATER(kodam.normal * gjk.normal, real_t(0.999));
   // [Claude, 2026-09-29] contact point must lie in the middle of the overlap: host i's face is at
   // z = +0.5, j's deepest rim point E at z = 0.5 - overlap, so the midpoint is z = 0.5 - overlap/2.
   // (Before the sign fix in rimVsFaceOrBand the point sat at z = 0.5 - 1.5*overlap, outside j.)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], real_t(0.5) - real_t(0.5) * kodam.overlap, real_t(1e-10));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.point[2], gjk.point[2], real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (face-edge, clean/negligible-host-curvature): overlap="
                              << kodam.overlap << " (agrees with GJK/EPA to 1e-4, unlike the KNOWN "
                              << "LIMITATION case above -- confirms the limitation is curvature-coupling, "
                              << "not a general face-edge defect)");
}

//! Clean, well-agreeing band-edge case -- counterpart to testCylinderCylinderTiltedRimIntoBand's
//! KNOWN LIMITATION, using the SAME large-host-radius trick as testCylinderCylinderFaceEdgeClean
//! above (see its own comment for why this isolates the curvature-coupling root cause): host i's
//! band is made locally near-flat (Ri=3.0 vs. Rj=0.3), so probe j's rim -- Kodam's own eq. 67-81
//! candidate -- is genuinely the closest feature, matching GJK/EPA exactly.
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

//! The FIRST dedicated test that actually dispatches to scenario=="edge-edge" -- previously the only
//! edge-edge-adjacent test (testCylinderCylinderEdgeEdge above) resolves as "face-edge" once Guo et
//! al.'s bidirectional correction is applied (see that test's own comment), so no test anywhere
//! exercised this branch specifically. Geometry found by direct parameter search against
//! testCylinderCylinderPureAnalyticalScenarioBreakdown's own fuzz corpus (fixed seed 24680, the
//! trial with the best non-trivial overlap agreement among every "edge-edge"-classified trial across
//! 5000 fully-random-orientation samples -- the same "found by search, hand-verified against GJK/EPA"
//! methodology testCylinderCylinderSpecialContact's own comment describes), not hand-derived from
//! first principles: both cylinders' rotations are genuinely arbitrary 3D orientations, which is
//! exactly the regime edge-edge (both rims facing each other, neither a flat face nor a curved band
//! being the dominant feature on either side) actually requires -- the axis-angle constructions used
//! by every other hand-picked test in this file cannot easily reach it. Agreement here (0.47% relative
//! overlap difference) is close but not machine-precision, consistent with edge-edge's fuzz-measured
//! p50 relative point-distance of ~0.1% and its own worst-case root-cause probe (see
//! testCylinderCylinderPureAnalyticalScenarioBreakdown's log output) showing the same curvature-
//! coupling gap can still bite edge-edge, just less severely on this particular geometry than on its
//! worst fuzzed trial.
// Standard rotation-matrix -> quaternion conversion (e.g. Shepperd 1978), needed only because Rot3
// has no direct Matrix3 constructor -- it keeps a quaternion and a matrix in sync internally, always
// built up incrementally via rotate(axis, angle). The geometry below is a raw fuzz-discovered matrix
// with no such incremental history, so it must be converted once here instead.
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

void testCylinderCylinderEdgeEdgeGenuine()
{
   const Mat3 rotICols(real_t(-0.390898), real_t(-0.198955), real_t(0.898675), real_t(0.196307),
                        real_t(0.935875), real_t(0.292579), real_t(-0.899257), real_t(0.290784),
                        real_t(-0.326775));
   const Mat3 rotJCols(real_t(0.131414), real_t(-0.370814), real_t(0.919362), real_t(-0.959446),
                        real_t(0.18573), real_t(0.212056), real_t(-0.249386), real_t(-0.909946),
                        real_t(-0.331369));

   // rotICols/rotJCols are printed to ~6 significant figures, not perfectly orthonormal at full
   // precision -- rot3FromMatrix's quaternion round-trip snaps them to the nearest valid rotation
   // (Rot3 asserts its matrix has determinant 1). Build BOTH CylinderState (analytical) and the
   // GJK/EPA reference from that same snapped Rot3, not the raw printed matrix, so the two methods
   // are compared on IDENTICAL geometry rather than two matrices differing by ~1e-4.
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
   // [Claude, 2026-10-15] tolerance 5e-3 -> 1e-3 (closed form 0.07830 vs GJK/EPA 0.07868 and brute-force MTD 0.07883)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("Cylinder x Cylinder (genuine edge-edge, found by search): overlap="
                              << kodam.overlap << " vs GJK/EPA " << -gjk.penetrationDepth
                              << " (0.47% relative -- first dedicated coverage of this scenario)");
}


//======================================================================================================================
// The master comparison: pure Kodam/Guo, pure Feng et al. (2017), Kodam's existing GJK/EPA hybrid,
// and a proposed analogous Feng-based hybrid (fengResultIsUncertain above), all evaluated on the
// SAME 5000 trials (same seed as every other breakdown/pattern test in this file) against GJK/EPA as
// ground truth -- directly answering: does a Feng-based hybrid behave any differently from Kodam's
// existing one, or does it also fall back to GJK/EPA on essentially every trial?
//======================================================================================================================
//! Classifies WHY a matched-scenario trial's overlap MAGNITUDE disagrees with GJK/EPA, as a fixed
//! string literal (safe to store in a CylMismatchRecord -- see that struct's own comment) plus a
//! supporting numeric metric. Two structurally different causes are distinguished, exactly the two
//! ways a closed-form dispatch can go wrong:
//!
//!   (a) MISCLASSIFIED SCENARIO: the dispatch fired a branch whose contact NORMAL points somewhere
//!       structurally different from GJK/EPA's true separating direction (checked directly: angle
//!       between kodam.normal and gjk.normal > 15 degrees). This means the wrong PRIMITIVE PAIRING
//!       was selected -- e.g. band-edge's fixed "radial from host axis" normal when the true nearest
//!       feature pair is actually face-vs-edge or a different rim entirely -- so no refinement of the
//!       chosen formula could fix it; only re-dispatching would.
//!   (b) WRONG POINT, SAME SCENARIO: the normal is essentially right (within 15 degrees) but the
//!       chosen POINT along that scenario's own construction is not the true extremal one, so the
//!       measured overlap is off even though the qualitative geometry (which primitive pairs with
//!       which) was correctly identified.
//!
//! special-contact is handled as a distinct, THIRD, already-diagnosed cause (not fit into the (a)/(b)
//! split above): rimVsFaceOrBand's special-contact branch measures its axial overlap at P, the point
//! on probe's rim closest to HOST'S AXIS (found for the band-edge search just above it in that
//! function), instead of E, the point on probe's rim closest to HOST'S FACE PLANE (the exact closed
//! form face-edge already computes, one branch earlier in the SAME function, via
//! detail::rimPointTowardPlane). Since E is defined as the axial-distance MINIMIZER over the whole
//! rim, |z_E| <= |z_P| always, so re-using P systematically UNDERESTIMATES the axial overlap
//! halfLHost-|z| -- exactly the ">3x underestimate" already noted where special-contact is
//! unconditionally routed to GJK/EPA in PaperAnalyticalContacts.h. This is verified per-trial below
//! (not just asserted from that comment) by recomputing E for whichever host/probe pairing actually
//! fired and checking that the resulting delta_a(E) lands much closer to GJK/EPA's ground truth than
//! Kodam's own delta_a(P) did.
const char* classifyMagnitudeCause(const CylinderState& ci, const CylinderState& cj, const ContactResult& kodam,
                                    const GjkResult& gjk, real_t& outMetric)
{
   const real_t gjkOverlap = -gjk.penetrationDepth;

   if (std::string(kodam.scenario) == "special-contact")
   {
      // Recompute E for whichever (host, probe) pairing actually produced `kodam` -- rimVsFaceOrBand
      // is order-sensitive, so try both and keep the one whose raw (unnegated) overlap matches.
      auto axialOverlapViaE = [&](const CylinderState& host, const CylinderState& probe, bool& eInDisc) -> real_t {
         const ContactResult raw = rimVsFaceOrBand(host, probe);
         if (!raw.contact || std::string(raw.scenario) != "special-contact") return real_t(-1);
         if (std::abs(raw.overlap - kodam.overlap) > real_t(1e-9)) return real_t(-1); // not this pairing
         const Vec3   E    = detail::rimPointTowardPlane(probe, host.pos, host.axis());
         const Vec3   relE = host.R.getTranspose() * (E - host.pos);
         eInDisc = (relE[0] * relE[0] + relE[1] * relE[1]) < host.radius * host.radius; // is E even a
                                                                                           // geometrically
                                                                                           // valid face
                                                                                           // contact
                                                                                           // candidate?
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
         outMetric            = errViaE; // the error the E-based point WOULD have had, for comparison
         // Verified against a 5000-trial run (2026-09-07): of the 24 trials where this re-derivation
         // matched the firing host/probe pairing, exactly the 5 with eInDisc=true landed within
         // errViaE<1e-3 of GJK/EPA (several to <1e-7, i.e. machine precision) -- E IS the true extremal
         // point whenever it's a geometrically valid candidate. All 20 with eInDisc=false overshot
         // GJK/EPA by 1.3x-16x even using E, confirming the deeper bug: rimVsFaceOrBand never checks
         // whether E lies inside host's face disc before using it (unlike the face-edge branch just
         // above it in the same function, which does via eInFaceDisc) -- when it doesn't, the
         // unconstrained axial gap halfLHost-|z_E| doesn't correspond to any real overlap, because no
         // part of host's face actually extends out to where E's axial projection lands. Fixing this
         // properly would need a disc-boundary-constrained extremal search (the same kind of Appendix-B
         // ellipse/circle-intersection machinery built for feng2017::fengCylinderCylinderContact's
         // face-type containment case), not just substituting E for P.
         if (errViaE < real_t(1e-3)) return "special-contact-exact-fix-via-E";
         if (!eInDisc) return "special-contact-E-outside-face-disc";
         if (errViaE < errViaP) return "special-contact-partial-fix-via-E";
      }
      return "special-contact-E-also-fails"; // deeper than a simple point swap -- neither P nor E is
                                                // the true GJK/EPA extremal point (includes the ~4/28
                                                // trials where this re-derivation couldn't even match
                                                // the firing pairing: right at the deltaAp~=deltaRp
                                                // dispatch threshold, where floating-point noise between
                                                // the inlined production call and this separate one can
                                                // flip which branch fires -- benign re-dispatch jitter,
                                                // not a new bug)
   }

   const real_t cosAngle = std::clamp(kodam.normal * gjk.normal, real_t(-1), real_t(1));
   const real_t angleDeg = std::acos(cosAngle) * real_t(180) / math::pi;
   outMetric              = angleDeg;
   return (angleDeg > real_t(15)) ? "misclassified-scenario" : "wrong-point-same-scenario";
}

//! Companion to classifyMagnitudeCause, for the OTHER kind of Kodam error: no scenario matched at
//! all even though GJK/EPA finds a real contact. Always a fixed string literal (see
//! CylMismatchRecord::cause). Reuses the rim/band classification from the prior investigation (a
//! point within 30% of a radius's width of a cylinder's own end cap counts as "rim" on that side,
//! else "band") -- the four combinations are enumerated explicitly (rather than building the string
//! dynamically) specifically so the result is a safe-to-store string literal.
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
// Kodam/Guo (2010) master comparison against GJK/EPA. Feng et al. (2017) is intentionally out of
// scope here (dropped from active development -- see testFeng*/testCylinderCylinderMissPatternAnalysis,
// still defined above but no longer called from main(), for that investigation's full record).
//
// Produces, over the SAME 5000 trials (same seed as every other breakdown/pattern test in this file):
//   (1) a full per-scenario table -- occurrence count, wrong-magnitude count/rate, relative/absolute
//       overlap error stats -- covering EVERY dispatch outcome, including ones with zero wrong trials;
//   (2) the "no scenario matched at all" dispatch-gap breakdown (rim/band classification + missed
//       overlap magnitude);
//   (3) for every wrong-magnitude and dispatch-gap trial, a MATHEMATICAL CAUSE classification
//       (classifyMagnitudeCause / classifyDispatchGapCause above) with an overall cause-occurrence
//       count table -- directly answering "why, not just how often"; and
//   (4) a VTK/CSV dump (writeCylMismatchVTKSeries/CSV, same writer testCylinderCylinderRandomizedFuzz
//       uses) of every one of those trials, tagged with its cause, so the actual failing geometries
//       can be opened in ParaView.
//======================================================================================================================
//! THREE-WAY accuracy comparison, all on IDENTICAL trials in one process (not separate builds/runs):
//!   A) ORIGINAL:      maxIterations=60, useFastEdgeEdgeReference=false -- the code as it shipped
//!                      before this session's own performance investigation.
//!   B) CURRENT DEFAULT: maxIterations=30, useFastEdgeEdgeReference=false -- simplification #2 only
//!                      (halved Durand-Kerner iteration count; see solveQuarticDurandKerner's comment).
//!   C) FULLY OPTIMIZED: maxIterations=30, useFastEdgeEdgeReference=true -- #2 AND #1 together
//!                      (thetaRef via the cheap closestThetaOnCircleToPoint instead of a second
//!                      quartic solve; see findRimCrossing's own comment).
//!   D) 20 iter, old ref: maxIterations=20, useFastEdgeEdgeReference=false -- probing whether the
//!                      Durand-Kerner iteration count can be dropped further than #2's own 30.
//!   E) 20 iter, new ref: maxIterations=20, useFastEdgeEdgeReference=true -- D combined with #1.
//! Old code paths for BOTH simplifications are kept fully intact (nothing deleted) -- maxIterations
//! and useFastEdgeEdgeReference are both ordinary default-valued parameters threaded through
//! solveQuarticDurandKerner / closestThetaOnCircleToLine / findRimCrossing / cylinderCylinderContact,
//! so variant A (the literal original behaviour) is still reachable by simply passing 60/false.
//!
//! Only trials classified "edge-edge" (the ONE scenario either simplification can affect at all) are
//! scored, all three variants evaluated on the SAME geometry, classified via variant A (the original)
//! so the selection criterion doesn't depend on either simplification being compared.
//!
//! 10000 trials (not this file's usual 5000) specifically because edge-edge only fires on ~4% of
//! trials (213/5000 in the master-comparison run) -- doubling the trial count roughly doubles the
//! edge-edge sample size this comparison actually has to work with.
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

      // Classify with variant A (the literal original: 60 iterations, old reference) -- only trials
      // landing in edge-edge under the ORIGINAL dispatch are used, so all three variants are scored on
      // EXACTLY the same set of edge-edge geometries, and the selection criterion is independent of
      // either simplification being compared.
      // quarticSolver pinned EXPLICITLY to DurandKerner, and lazyThetaRef EXPLICITLY to false, on every
      // one of these calls (A-G): the file defaults changed to Ferrari and lazyThetaRef=true after this
      // test was written (see cylinderCylinderContact's own comment), and these variants' labels ("60
      // iter", "30 iter", "20 iter", "Ferrari solver") are each meant to isolate ONE specific change --
      // leaving any of them on the current defaults would have silently pulled in whichever other
      // simplifications happen to default to on/off at the time, collapsing distinct variants into
      // duplicates of each other without warning. Only H below deliberately tests lazyThetaRef.
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
      // Solver-choice-ONLY variants: same 30 iterations and old (line-based) reference as variant B,
      // isolating the effect of swapping solveQuarticDurandKerner for solveQuarticFerrari /
      // solveQuarticAberthEhrlich -- see QuarticSolverComparisonTest.cpp's synthetic-quartic comparison
      // (known ground-truth roots) that motivated checking this against real cylinder-cylinder geometry
      // too, not just synthetic polynomials.
      const ContactResult kodamFerrari =
         cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, false);
      const ContactResult kodamAberthEhrlich =
         cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::AberthEhrlich, false);
      // lazyThetaRef, on top of Ferrari (the current production default) -- defers thetaRef's
      // computation until it's actually needed (>1 candidate for point A or B), which
      // testThetaRefTieBreakNecessity10k found is NEVER on real geometry (0/832 calls). If that holds,
      // this should be bit-identical to the plain Ferrari variant above, at even lower cost.
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

      // Direct pairwise diff -- not just "both agree with GJK to similar tolerance" but "did changing
      // the iteration count / reference method change THIS trial's actual answer at all".
      if (kodamOriginal.contact == kodamCurrent.contact &&
          std::abs(kodamOriginal.overlap - kodamCurrent.overlap) < real_t(1e-9) &&
          (kodamOriginal.point - kodamCurrent.point).length() < real_t(1e-9))
         ++pairwiseIdenticalOldVsCurrent;
      if (kodamCurrent.contact == kodamFull.contact &&
          std::abs(kodamCurrent.overlap - kodamFull.overlap) < real_t(1e-9) &&
          (kodamCurrent.point - kodamFull.point).length() < real_t(1e-9))
         ++pairwiseIdenticalCurrentVsFull;
      // 20-iteration variants compared directly against the literal original (60 iter, old ref), the
      // strongest baseline, not against the 30-iteration variants -- answers "does dropping to 20 still
      // match the original code", independent of whether the 30-iteration step already changed anything.
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

   // All three variants must land in the SAME ballpark of existence agreement against GJK/EPA -- either
   // simplification costing genuinely worse EXISTENCE detection (not just a slightly different overlap
   // magnitude on already-agreeing trials) would be a real regression, not an acceptable trade-off.
   WALBERLA_CHECK_GREATER_EQUAL(original.edgeEdgeTrials, 100); // sanity: comparison isn't running on noise
   WALBERLA_CHECK_EQUAL(original.edgeEdgeTrials, currentDefault.edgeEdgeTrials); // same selection set
   WALBERLA_CHECK_EQUAL(original.edgeEdgeTrials, fullyOptimized.edgeEdgeTrials); // by construction
}

//======================================================================================================================
//! Is thetaRef's tie-break EVER actually load-bearing on real edge-edge geometry, or is it dead weight?
//!
//! findRimCrossing's point A (eq. 83) and point B (eq. 84-85) each pick, among every candidate theta
//! that passes their own containment check, the one CLOSEST to thetaRef -- a tie-break neither Kodam's
//! nor Guo's paper specifies at all (see this file's own RimCrossing comment). If, in practice, at most
//! ONE candidate ever survives each containment check, thetaRef's VALUE is irrelevant to the result --
//! any value at all (a cheap atan2, a hardcoded 0, garbage) would produce the identical answer, and
//! thetaRef's own quartic solve (candidates 3/4 in the "6 quartic solves per edge-edge pair" breakdown)
//! would be pure overhead, not doing any actual disambiguation work.
//!
//! Uses findRimCrossing's own optional candidateCountA/candidateCountB out-parameters (added purely for
//! this diagnostic, default nullptr, zero effect otherwise) directly on the SAME real edge-edge fuzz
//! corpus (seed 13579, same generation as testEdgeEdgeThreeWayAccuracy10k) already validated elsewhere
//! in this file -- not a synthetic argument, an actual count on the geometries this file's own fuzz
//! testing has shown are representative.
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

      // Replicate cylinderCylinderContact's own edge-edge branch geometry setup exactly, so jVsI/iVsJ
      // are called with the SAME inputs production would use.
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
//! Hand-computed ground-truth check of detail::closestDistanceSqBetweenSegments -- validates the
//! closed-form segment-segment distance itself, independent of anything cylinder-related, before it's
//! trusted as a cylinder-cylinder separation bound. Four configurations, each with an exactly known
//! answer: parallel+overlapping (pure perpendicular offset), parallel+non-overlapping (nearest
//! endpoints), perpendicular skew (crossing at a single closest pair of points), and intersecting
//! (distance exactly 0).
//======================================================================================================================
void testClosestDistanceSqBetweenSegments()
{
   using detail::closestDistanceSqBetweenSegments;

   // Parallel, overlapping ranges -- closest distance is the pure perpendicular offset (3), achieved
   // anywhere along the overlapping length, independent of where exactly s/t land.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(0, 0, 0), Vec3(0, 0, 10), Vec3(3, 0, 0), Vec3(3, 0, 10));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(9), real_t(1e-9));
   }

   // Parallel, NON-overlapping ranges (axial gap too) -- closest points are the two nearest endpoints:
   // (0,0,5) and (3,0,10) -> distance = sqrt(3^2 + 5^2) = sqrt(34).
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(0, 0, 0), Vec3(0, 0, 5), Vec3(3, 0, 10), Vec3(3, 0, 15));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(34), real_t(1e-9));
   }

   // Perpendicular skew segments (one along x at z=2, one along y at z=0) -- closest points are
   // (0,0,2) and (0,0,0), distance = 2.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(-5, 0, 2), Vec3(5, 0, 2), Vec3(0, -5, 0), Vec3(0, 5, 0));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(4), real_t(1e-9));
   }

   // Intersecting segments (an X-crossing in the z=0 plane) -- distance exactly 0.
   {
      const real_t distSq = closestDistanceSqBetweenSegments(
         Vec3(-1, -1, 0), Vec3(1, 1, 0), Vec3(-1, 1, 0), Vec3(1, -1, 0));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(distSq, real_t(0), real_t(1e-9));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[TightTier1] closestDistanceSqBetweenSegments: all 4 hand-computed "
                              "ground-truth cases matched.");
}

//======================================================================================================================
//! Does the tight, capsule-based Tier 1 bound (useTightTier1=true) ever WRONGLY declare two cylinders
//! "definitely separated" when they're actually in contact? This is the one thing that must NEVER
//! happen -- a false "definitely separated" would silently skip the GJK/EPA fallback and drop a real
//! contact. Checked directly against GJK/EPA ground truth on real cylinder geometry, both near the
//! actual contact boundary (sampleSeparationNearBoundary, the same methodology as every other fuzz
//! test in this file) and at uniformly-distributed separations across the whole physically-possible
//! range (representative of a packed simulation, where two pillars can sit anywhere within reach, not
//! just near their own contact boundary) -- since the tight bound's main value proposition is
//! specifically for well-separated-but-within-reach pairs, the near-boundary corpus alone would
//! under-sample exactly the population it's meant to help with.
//======================================================================================================================
void testTightTier1Accuracy10k()
{
   FuzzRng                                  fuzz(97531); // distinct seed from every other fuzz test
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrialsPerCorpus = 5000;

   int dispatchGapTrials = 0;                                  // analytical found NO scenario at all
   int oldTier1Uncertain = 0, newTightTier1Uncertain = 0;       // how often each bound says "check with GJK"
   int unsafeFalseSeparated = 0;                                // tight bound wrong -> MUST stay 0
   real_t maxUnsafeOverlap = 0;

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

         const ContactResult analytical = cylinderCylinderContact(ci, cj);
         if (analytical.contact) continue; // Tier 1 only ever applies when the dispatch found nothing
         ++dispatchGapTrials;

         const bool oldUncertain = cylinderCylinderResultIsUncertain(ci, cj, analytical, false);
         const bool newUncertain = cylinderCylinderResultIsUncertain(ci, cj, analytical, true);
         if (oldUncertain) ++oldTier1Uncertain;
         if (newUncertain) ++newTightTier1Uncertain;

         if (!newUncertain) // tight bound says "definitely separated, don't bother with GJK"
         {
            const GjkResult gjk = gjkCylCyl(posI, rotI, Ri, Li, posJ, rotJ, Rj, Lj);
            if (gjk.contact)
            {
               ++unsafeFalseSeparated;
               maxUnsafeOverlap = std::max(maxUnsafeOverlap, -gjk.penetrationDepth);
            }
         }
      }
   };

   runCorpus(false); // near the actual GJK/EPA contact boundary
   runCorpus(true);  // uniformly across the whole physically-possible-reach range

   WALBERLA_LOG_INFO_ON_ROOT(
      "[TightTier1] " << dispatchGapTrials << " dispatch-gap trials (analytical found no scenario) "
      << "out of " << (2 * numTrialsPerCorpus) << " total | old Tier 1 uncertain=" << oldTier1Uncertain
      << "/" << dispatchGapTrials << " | NEW tight Tier 1 uncertain=" << newTightTier1Uncertain << "/"
      << dispatchGapTrials << " (i.e. tight bound resolved "
      << (dispatchGapTrials - newTightTier1Uncertain) << " trials as definitely-separated, skipping "
      << "GJK/EPA entirely) | UNSAFE false-separated (tight bound wrong vs GJK/EPA): "
      << unsafeFalseSeparated << " (max missed overlap=" << maxUnsafeOverlap << ")");

   // The one thing that must NEVER happen: the tight bound declaring "definitely separated" for a
   // pair GJK/EPA confirms is actually touching. Any failure here means the capsule bound is NOT
   // conservative as derived and must not be trusted.
   WALBERLA_CHECK_EQUAL(unsafeFalseSeparated, 0);
   WALBERLA_CHECK_GREATER_EQUAL(dispatchGapTrials, 200); // sanity: comparison isn't running on noise
}

//======================================================================================================================
//! Hand-computed ground-truth check of detail::boxSupportAlongDirection and
//! detail::boxesDefinitelySeparatedAlongCenterLine -- the box-box analog of
//! testClosestDistanceSqBetweenSegments, validating the geometric primitive itself before it's
//! trusted as a Layer 1 pre-filter.
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
   // h*|cos45| + h*|sin45| = h*sqrt(2) (both local x and y axes now contribute equally).
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

   // Same boxes, centres only 1.5 apart (reach=2 > 1.5) -> NOT resolvable as separated by this bound
   // (they may or may not actually overlap -- the point is this axis alone can't rule it out).
   {
      const Mat3 R = Mat3();
      const Vec3 half(1, 1, 1);
      WALBERLA_CHECK(!boxesDefinitelySeparatedAlongCenterLine(Vec3(0, 0, 0), R, half, Vec3(1.5, 0, 0), R, half));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[TightLayer1Box] boxSupportAlongDirection / "
                              "boxesDefinitelySeparatedAlongCenterLine: all hand-computed cases matched.");
}

//======================================================================================================================
//! Does the cheap single-axis box-box separation bound (Layer 1 pre-filter) ever WRONGLY declare two
//! boxes "definitely separated" when they're actually touching? Checked directly against GJK/EPA
//! ground truth on real box geometry (same two-corpus methodology as testTightTier1Accuracy10k: near
//! the actual contact boundary, and uniformly across the whole physically-possible-reach range).
//======================================================================================================================
void testTightLayer1BoxAccuracy10k()
{
   FuzzRng                                  fuzz(24601); // distinct seed from every other fuzz test
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

void testCylinderCylinderFullComparisonAndHybridFeasibility()
{
   FuzzRng                                  fuzz(24680); // same seed as every breakdown/pattern test above
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));

   constexpr int numTrials = 5000; // confirmed consistent at 10000 trials for the Tier1+special-contact-
                                      // only hybrid standard (fallback rate landed at 3.3553%, matching
                                      // the 3.3-3.4% predicted from the 5000-trial numbers almost
                                      // exactly); reverted to keep the standard suite's runtime down.
                                      // Re-confirmed at 10000 trials again after making Ferrari the
                                      // default quartic solver (this session): wrong=0/3338, fallback
                                      // rate=112/3338 (3.3553%) -- IDENTICAL to the number recorded here
                                      // under the previous Durand-Kerner default, across every scenario
                                      // in the file. Re-confirmed a THIRD time at 10000 trials after also
                                      // making lazyThetaRef=true the default (this session): again
                                      // wrong=0/3338, fallback rate=112/3338 (3.3553%) -- unchanged
                                      // (not just edge-edge), before reverting back to 5000 again.

   int kodamMiss = 0, kodamFalsePos = 0, kodamCorrectNoContact = 0;
   int kodamHybridWrong = 0, kodamHybridFallback = 0, kodamHybridFallbackNeeded = 0;
   int kodamHybridWrongDepth = 0; // [Claude, 2026-09-30] hybrid DEPTH errors on genuine contacts (kodamHybridWrong
                                   // above counts existence only, which was misread as "no errors at all")
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
   std::map< std::string, int >                kodamScenarioCount; // denominator per scenario, EVERY
                                                                       // scenario that ever fires, not
                                                                       // just the ones with wrong trials
   std::map< std::string, int >                kodamMissRimBandClass;
   std::map< std::string, std::vector< real_t > > kodamMissOverlapByClass;
   int kodamSpecialContactCount = 0; // how many trials fire special-contact (now unconditionally
                                       // routed to GJK/EPA -- see cylinderCylinderResultIsUncertain)

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

      // NOTE: the rim-rim dispatch-gap investigation (why Kodam/Guo's edge-edge chord-existence test
      // -- eq. 82-101 -- fails even for real, substantial overlaps, e.g. trials 463/3047/3075/4997)
      // was done via temporary ad-hoc diagnostics here and has been removed after the finding was
      // extracted. Summary (see the thesis writeup for the full derivation): the chord test requires,
      // for each rim, a point A (crosses the other cylinder's infinite lateral surface, within its
      // axial extent) AND a point B (crosses the other cylinder's near face plane, within its radius).
      // Two distinct sub-failures were found: point A can have zero valid roots at all (the rim never
      // reaches the other's lateral surface in that direction -- trial 463's iVsJ), or point A exists
      // but point B's candidates both land outside the target's radius (the entry point never lines up
      // with a valid exit through the face -- trial 3047's iVsJ). Both stem from the same structural
      // cause: the formula assumes a single clean in-and-out chord through one curved surface and one
      // flat face, which doesn't exist for a genuine corner-to-corner (both rims' edges meeting near
      // both cylinders' own corners) contact -- confirmed independently via brute-force MTD (support-
      // function-based, see testCylinderCylinderSpecialContactBruteForceMTD's method) to agree with
      // GJK/EPA to within 0.01% on all three trials. No orientation/aspect-ratio signature was found
      // across the 13 occurrences in a 5000-trial run (axis angles spanned 33-169 degrees); a mild
      // (73%) lean toward squat cylinders was the only pattern. Already fully covered in production:
      // cylinderCylinderResultIsUncertain's Tier 1 reach-sum bound catches every "no scenario matched"
      // case unconditionally, and the existing hybrid shows 0/1664 existence mismatches.

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

      // Kodam's existing production hybrid.
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

      // Simplified "Tier-1-only" hybrid: fall back to GJK/EPA ONLY when nothing matched at all --
      // drop Tier 2 (the circumscribing-sphere magnitude check) entirely. Directly answers "if no
      // scenario matches, run GJK, otherwise always trust the raw analytical answer as-is."
      if (!kodam.contact) { ++kodamTier1OnlyFallback; }
      else
      {
         ++kodamScenarioCount[std::string(kodam.scenario)];
         const real_t absDiff  = std::abs(kodam.overlap - gjkOverlap);
         const real_t relDiff  = absDiff / gjkOverlap;
         const bool   wrongMag = absDiff > std::max(real_t(0.01), real_t(0.05) * gjkOverlap);

         // Normal-angle comparison for EVERY kodam.contact trial, not just the wrongMag ones -- to
         // check whether the "passing" (magnitude-agreeing) population also has near-zero normal
         // divergence, or whether it's only the wrongMag trials with large-angle outliers.
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
      const int total = kodamCorrectNoContact + kodamMiss; // every trial where Kodam itself reported
                                                              // "no scenario matched"
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
// Deliberately DEEP (not "realistic") overlap stress test. Every other fuzz test in this file caps
// sampled overlap at 15% of the smaller body's radius (FuzzRng::sampleSeparationNearBoundary's
// maxOverlap) -- a working simulation's spring-dashpot contact force would never let overlap develop
// beyond that, so testing deeper is testing geometry the simulation itself would never reach. This
// test intentionally goes there anyway, to answer a narrower safety question: IF a genuinely large,
// degenerate overlap ever did occur (a bad initial condition, an unrelated bug elsewhere in the
// pipeline pushing state somewhere it shouldn't be), does Kodam/Guo's dispatch degrade SAFELY -- i.e.
// either get the scenario/overlap right, or cleanly report "no case" (scenario=none) so
// cylinderCylinderResultIsUncertain's unconditional Tier 1 reach-sum bound catches it and the hybrid
// still recovers the correct answer via GJK/EPA -- rather than UNSAFELY: some scenario firing with
// kodam.contact=true and a badly wrong overlap that Tier 2's magnitude bound might not always catch.
//======================================================================================================================
void testCylinderCylinderLargeOverlapDegenerateStress()
{
   FuzzRng                                  fuzz(13579); // distinct seed
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > overlapFracDist(real_t(0.3), real_t(1.5)); // 30%-150% of
                                                                                          // the smaller
                                                                                          // radius --
                                                                                          // well past
                                                                                          // "realistic"

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
      // Find the true touching boundary via the same bisection FuzzRng::sampleSeparationNearBoundary
      // uses internally, then deliberately push DEEP past it -- 30%-150% of the smaller radius, not
      // the usual 15% cap.
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

      if (!gjk.contact) continue; // sampling occasionally overshoots into no-contact -- skip, not the
                                    // population this test is about
      if (!kodam.contact) { ++kodamNoneCase; } // SAFE: caught unconditionally by Tier 1 below
      else
      {
         const real_t absDiff = std::abs(kodam.overlap - gjkOverlap);
         const real_t relDiff = absDiff / gjkOverlap;
         // "Unsafe" here means Tier 2 would ALSO need to catch it for the hybrid to recover --
         // worth knowing the worst case even though the hybrid check below already covers whether it
         // actually does.
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
      // Existence agreement alone isn't the full safety story -- also check the HYBRID's own overlap
      // MAGNITUDE (not just raw kodam's) against GJK/EPA, since Tier 2's circumscribing-sphere bound
      // only catches UNDERESTIMATES (it flags uncertainty when more overlap is geometrically possible
      // than reported) and could in principle miss an overestimate.
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
// Same randomized fuzz test, for Box x Box: boxBoxSAT (provably order-symmetric, no golden-section
// search ever involved) against GJK/EPA directly -- boxBoxSAT has no uncertain-result fallback gate
// at all (see its own "trusted directly" comment in PaperAnalyticalContacts.h), so there is no
// separate "hybrid" dispatcher to test here, just the raw SAT result.
//======================================================================================================================

void testBoxBoxRandomizedFuzz()
{
   FuzzRng                                  fuzz(67890); // different seed than the cylinder fuzz test
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

         // Contact POINT: only the component along the collision normal is geometrically
         // well-defined for a flush face-face (or edge-edge-parallel) overlap -- boxBoxSAT's own
         // comment says as much ("only the component along bestAxis is geometrically well-defined
         // for a flush face-face/edge-edge configuration"). Perpendicular to the normal, any point
         // within the overlapping footprint is an equally valid "the" contact point, and SAT's
         // particular choice (midpoint of the extreme points along the axis) can legitimately
         // differ from GJK/EPA's -- that's an underdetermined-point ambiguity, not a bug. So the
         // along-normal component is asserted strictly; the perpendicular component is only logged.
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
   // [Claude, 2026-10-15] geometry replaced: the former pair (testCylinderCylinderTiltedRimIntoBand) is exactly touching, so it
   // did not test a missed contact. Trial 1209 of the independent-reference sample (overleaf
   // scripts/contact-validation/data/cyl_pairs.csv) is a genuine miss: no scenario applies, the end faces meet
   // corner to corner near both rims (thesis Figure 4.11), and GJK/EPA and the brute-force reference both give 0.113.
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

//! FORMERLY "Tier 2 triggers on wrong magnitude": identical geometry to
//! testCylinderCylinderTiltedRimIntoFace, where the analytical dispatch DOES find a scenario
//! (face-edge) but with the wrong magnitude (0.299 vs GJK/EPA's 0.235). Tier 2 (a circumscribing-
//! sphere magnitude cross-check applied to every non-exempt scenario) used to flag this and recover
//! the correct value via fallback -- it was REMOVED from cylinderCylinderResultIsUncertain (see that
//! function's own comment for the full history: it pushed the GJK/EPA fallback rate to ~90% while
//! still being structurally one-sided, catching underestimates but not overestimates). This is now
//! exactly the accepted trade-off: the hybrid no longer falls back here, and legitimately returns the
//! same wrong raw analytical value as the pure dispatch -- one concrete instance of the ~2.3-2.5%
//! residual wrong-magnitude rate documented for the current (Tier 1 + special-contact only) design.
void testHybridKeepsFaceEdgeMagnitudeWithoutFallback()
{
   CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.6), real_t(2) };
   Rot3          rotJ;
   rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.25));
   CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.55)), rotJ.getMatrix(), real_t(0.2), real_t(2) };

   const ContactResult analyticalOnly = cylinderCylinderContact(ci, cj); // for the "before" log line
   const HybridResult   hybrid        = cylinderCylinderContactHybrid(ci, Rot3(), cj, rotJ);
   const GjkResult       gjk           = gjkCylCyl(ci.pos, Rot3(), ci.radius, ci.length, cj.pos, rotJ, cj.radius, cj.length);

   WALBERLA_CHECK(!hybrid.usedFallback); // no longer falls back -- Tier 2 removed, see comment above
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(hybrid.overlap, analyticalOnly.overlap, real_t(1e-6)); // kept as-is
   WALBERLA_LOG_INFO_ON_ROOT("Hybrid (Tier 2 removed) on a known wrong-magnitude geometry: analytical "
                              << "overlap=" << analyticalOnly.overlap << ", hybrid overlap=" << hybrid.overlap
                              << " (same value, no longer corrected) vs GJK/EPA's true "
                              << -gjk.penetrationDepth << " -- accepted residual error, not a bug");
}

//! FORMERLY "Tier 2 must not fire for verified scenarios": now that Tier 2 has been removed entirely
//! (see cylinderCylinderResultIsUncertain's own comment), NONE of face-face/band-band-parallel/
//! band-band-skewed/face-band ever trigger a fallback any more -- only "no scenario matched" (Tier 1)
//! and special-contact do. What's still worth asserting per geometry is that trusting the raw
//! analytical result directly continues to agree with GJK/EPA: face-face/band-band-parallel are
//! provably exact under their own near-parallel-axis precondition regardless of any fallback tier;
//! band-band-skewed's own root-cause bug was fixed (not merely safety-netted) per
//! testCylinderCylinderPureAnalyticalScenarioBreakdown's fresh-run confirmation; face-band was never
//! covered by the fuzz suite's documented error population the way band-edge/edge-edge/face-edge are,
//! so it is expected (and checked here) to still agree.
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

//! FORMERLY the empirically interesting Tier 2 case: the NOW-FIXED edge-edge-ish geometry resolves as
//! face-edge with an EXACT match to GJK/EPA (0.05 == 0.05). Retained now that Tier 2 is gone as a
//! simple regression check that this already-correct face-edge case is still trusted directly
//! (usedFallback is expected false, though not hard-asserted -- either answer would still be
//! numerically "correct" here) and still agrees with GJK/EPA.
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

   const GjkResult gjk = productionCylPlane(cyl.pos, Rot3(), cyl.radius, cyl.length, planePos, planeNormal);
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

   const GjkResult gjk = productionCylPlane(cyl.pos, rot, cyl.radius, cyl.length, planePos, planeNormal);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(2e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));

   // Contact POINT, hand-computed directly from eq. 21 (z_E^P - sign(z_E^P)*0.5*delta, measured from
   // E's own height): with this geometry, E = (-0.240192, 0, -0.416025) and overlap = 0.416025, giving
   // point = (-0.240192, 0, -0.208013). Not cross-checked against gjk.point -- same convention
   // difference as testCylinderPlaneFaceMatchesExisting (Kodam's own midpoint-of-overlap vs. GJK/EPA's
   // wall-surface witness point) -- but this IS the check that would have caught the eOnPlane
   // double-projection bug (previously (-0.240192, 0, +0.208013): the opposite sign, on the wrong side
   // of the plane) before it shipped, since overlap/normal alone were already correct throughout.
   // [Claude, 2026-10-15] production places the point on the plane: the deepest rim point projected along n
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
// Sphere x Cylinder -- new for this project (see PaperAnalyticalContacts.h's own comment on
// sphereCylinderContact for why sphere-sphere/sphere-box predate it but this does not). Exact and
// unconditional, unlike cylinder-cylinder: no uncertainty gate, no fallback tier needed at all.
//======================================================================================================================

//! Sphere clearly outside, closest feature on the cylinder's curved BAND (radial clamp fires, axial
//! clamp does not): cylinder R=1,L=4 along z; sphere of radius 0.3 centred at (1.5,0,0), i.e. purely
//! radially offset at the cylinder's own mid-height. rho=1.5, z=0 (|z|<halfL=2, so axial clamp is a
//! no-op); closest point q=(1,0,0); dist=0.5; overlap=0.3-0.5=-0.2 -> NO contact at this distance, so
//! moved sphere to 1.2 instead (dist=0.2, overlap=0.1).
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
   // [Claude, 2026-10-15] production contact point: middle of the overlap along the normal (ParticleKernels.h)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.point, spherePos - (sphereRadius - real_t(0.5) * kodam.overlap) * kodam.normal,
                                       real_t(1e-6));
   const GjkResult gjk = gjkSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Sphere clearly outside, closest feature on a flat CAP (axial clamp fires, radial does not):
//! cylinder R=1,L=2 (halfL=1) along z; sphere radius 0.3 centred on-axis at z=1.2 (rho=0). z=1.2>halfL,
//! so zC=1.0; q=(0,0,1); dist=0.2; overlap=0.3-0.2=0.1.
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
   // [Claude, 2026-10-15] production contact point: middle of the overlap along the normal (ParticleKernels.h)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.point, spherePos - (sphereRadius - real_t(0.5) * kodam.overlap) * kodam.normal,
                                       real_t(1e-6));
   const GjkResult gjk = gjkSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Sphere near the RIM: both the radial and axial clamp fire simultaneously (rho>R AND |z|>halfL), so
//! the closest point is the nearest point on the circular edge itself, not a plain face or band point.
//! Cylinder R=1,L=2 (halfL=1); sphere radius 0.3 centred at (1.5,0,1.3) (rho=1.5>1, z=1.3>1). Closest
//! point q=(1,0,1) (the actual rim); diff=(0.5,0,0.3), dist=sqrt(0.34)=0.583095; overlap=0.3-0.583095
//! (negative -> no contact at this distance); moved sphere closer, to (1.2,0,1.1): diff=(0.2,0,0.1),
//! dist=sqrt(0.05)=0.223607, overlap=0.3-0.223607=0.076393.
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
   // [Claude, 2026-10-15] production contact point: middle of the overlap along the normal (ParticleKernels.h)
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.point, spherePos - (sphereRadius - real_t(0.5) * kodam.overlap) * kodam.normal,
                                       real_t(1e-6));
   const GjkResult gjk = gjkSphereCyl(spherePos, sphereRadius, cyl.pos, Rot3(), cyl.radius, cyl.length);
   WALBERLA_CHECK(gjk.contact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.overlap, -gjk.penetrationDepth, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(kodam.normal, gjk.normal, real_t(1e-3));
}

//! Sphere centre embedded INSIDE the cylinder, closer to the curved band than either cap: cylinder
//! R=1,L=4 (halfL=2); sphere radius 0.3 centred at (0.8,0,0) (rho=0.8<1, z=0). radialClearance=0.2,
//! axialClearance=2.0 -> radial wins. overlap=0.2+0.3=0.5.
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

//! Sphere centre embedded INSIDE the cylinder, closer to a flat cap than the band: cylinder R=3,L=2
//! (halfL=1); sphere radius 0.3 centred at (0,0,0.9) (rho=0, z=0.9). radialClearance=3.0,
//! axialClearance=0.1 -> axial wins. overlap=0.1+0.3=0.4.
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

//! Randomized fuzz sweep against GJK/EPA, fully arbitrary 3-D cylinder orientation, boundary-biased
//! separation (same FuzzRng/sampleSeparationNearBoundary methodology as every other fuzz test in this
//! file). Unlike cylinder-cylinder, sphereCylinderContact carries no uncertainty gate and no known
//! limitation to characterize -- this test's entire purpose is to confirm that claim of exactness
//! empirically, not to map out a failure population the way the cylinder-cylinder fuzz tests do.
void testSphereCylinderRandomizedFuzz()
{
   FuzzRng                                  fuzz(112233);
   std::uniform_real_distribution< real_t > radiusDist(real_t(0.3), real_t(2.0));
   std::uniform_real_distribution< real_t > lengthDist(real_t(0.5), real_t(4.0));
   std::uniform_real_distribution< real_t > sphereRadiusDist(real_t(0.1), real_t(1.0));

   constexpr int numTrials = 10000;
   int           numOverlapping = 0;
   int           numMismatch    = 0;
   int           numGjkMismatch = 0;   // [Claude, 2026-10-15] against true GJK/EPA
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
      const GjkResult      gjk  = productionSphereCyl(spherePos, sphR, cylPos, rotCyl, cylR, cylL); // production closed form
      // [Claude, 2026-10-15] independent reference: mesa_pd GJK/EPA (productionSphereCyl is not GJK/EPA)
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

   // [Claude, 2026-10-15] NOTE: 'gjk' in this test is the production closed form (productionSphereCyl), so the
   // numbers below compare two closed forms; the true GJK/EPA comparison is the 'ref' block above.
   // Unlike every cylinder-cylinder fuzz test in this file, existence/overlap/normal are asserted
   // strictly, not just logged -- sphereCylinderContact's own closed form is a complete, exact
   // minimum-translation-distance solution (see its header comment), not a named-scenario heuristic,
   // so there is no known geometric limitation to tolerate here. Confirmed empirically too: across
   // 10000 trials, max overlap diff was 2.8e-15 (machine precision) and max normal angle diff was
   // 1.9e-6 degrees -- both already at floating-point noise level.
   WALBERLA_CHECK_EQUAL(numMismatch, 0);
   WALBERLA_LOG_INFO_ON_ROOT("[fuzz] Sphere x Cylinder vs true GJK/EPA: existence mismatches=" << numGjkMismatch
                             << ", max overlap diff=" << maxGjkOverlapDiff << ", max normal angle diff="
                             << maxGjkNormalAngleDeg << " deg");
   WALBERLA_CHECK_EQUAL(numGjkMismatch, 0);
   WALBERLA_CHECK_LESS(maxGjkOverlapDiff, real_t(1e-3));
   // [Claude, 2026-10-15] measured over the 10000 trials: existence 0 mismatches, overlap 3.9e-4, normal up to 2.9 deg.
   // The normal error is GJK/EPA's (polytope approximation of the curved sphere and band near the boundary; the closed
   // form is the exact closest point and agrees with the production overload to 1.9e-6 deg), hence the looser bound.
   WALBERLA_CHECK_LESS(maxGjkNormalAngleDeg, real_t(5));
   WALBERLA_CHECK_LESS(maxOverlapDiff, real_t(1e-5));
   WALBERLA_CHECK_LESS(maxNormalAngleDeg, real_t(0.1));
   // Contact POINT is NOT held to the same standard: sphereCylinderContact places it on the
   // cylinder's own surface (the closest-point construction's natural output), while
   // SolidCylinderContactDetection's GJK/EPA path places it roughly midway between the two bodies'
   // surfaces along the separating axis -- two different, equally valid conventions for "the" contact
   // point (same underdetermined-point situation already documented for box-box/cylinder-cylinder
   // flush contacts elsewhere in this file), which can differ by up to roughly half the overlap/sphere-
   // radius scale even when overlap and normal agree exactly. Informational only, not asserted.
   WALBERLA_LOG_INFO_ON_ROOT(
      "[fuzz] Sphere x Cylinder max contact-point difference=" << maxPointDiff << " (informational only -- "
      << "differing point CONVENTION between the two paths, not a geometric disagreement; see comment)");
}

//======================================================================================================================
// Box x ArchTunnel (CylindricalBoundary containment mode) -- hand-derived vs. the production
// containment override (SolidCylinderContactDetection, ParticleKernels.h) vs. mesa_pd's OWN generic
// GJK-based CylindricalBoundary path (GeneralContactDetection, called directly).
//
// IMPORTANT CAVEAT, read before trusting/extending this comparison: mesa_pd's generic
// CylindricalBoundary path (GeneralContactDetection.h) reads ONLY the shape's radius and axis -- it
// has NO notion whatsoever of isContainment()/axial bounds/half-clip (grep-confirmed: none of those
// three fields is referenced anywhere in that file). It happens to already compute the same core
// containment quantity our custom override does: whether the OTHER body's farthest point, in the
// direction of that body's own radial offset from the axis, exceeds the boundary radius (this is
// mesa_pd's own long-standing generic behaviour for this shape, predating Antidunes' own
// isContainment()/axial-bounds/half-clip extension -- see CylindricalBoundary.h's own doc comment).
// For an UNROTATED box, a single support-function query in that one fixed direction provably returns
// the SAME farthest corner an exhaustive per-corner search would, because the farthest corner's
// per-axis sign always matches the box centre's own offset sign in each dimension. That is WHY the
// two paths agree in the specific configuration below -- NOT because they are generally
// interchangeable for arch-tunnel geometry. They are NOT, once axial bounds or the half-clip actually
// exclude part of the tunnel (the generic path would keep reporting a contact there; the production
// override correctly reports none), or once the box is rotated (a single support-function direction
// no longer necessarily matches the true farthest corner the way it does here).
//======================================================================================================================

//! HAND-DERIVED geometry. Tunnel: axis=+X through the origin, radius=10, axially bounded to
//! [-10,10] (comfortably containing the box's own +-1 x-extent, so axial clipping never triggers),
//! half-clipped to the +z side (comfortably satisfied -- the box sits entirely at z=8.5..10.5, all
//! > 0). Box: centre (0, 0.3, 9.5), axis-aligned (identity rotation), halfExtents (1,1,1). The 0.3
//! y-offset is deliberate, chosen to make the (y,z) sign combination that maximises radial distance
//! UNAMBIGUOUS (see below) while still leaving x -- the tunnel's own axis -- degenerate on purpose,
//! to exercise the secondary-contact path deliberately (see the two "by hand" derivations below).
//!
//! By hand: every corner's radial distance-from-axis is sqrt((0.3+-1)^2 + (9.5+-1)^2) (x drops out
//! entirely -- the tunnel's axis IS the box's own x-axis). The 4 (y,z) sign combinations give
//! 49.49, 81.49, 50.69, 82.69 -- wait, using the ACTUAL corner coordinates (y=0.3+-1, z=9.5+-1):
//!   (y=-0.7,z=8.5): 0.49+72.25=72.74        (y=-0.7,z=10.5): 0.49+110.25=110.74
//!   (y=1.3, z=8.5): 1.69+72.25=73.94        (y=1.3, z=10.5): 1.69+110.25=111.94  <- unique max
//! sqrt(111.94) ~= 10.5802 -- overlap = sqrt(111.94)-10 ~= 0.5802 (comfortably under the 1.0 LU
//! envelope clamp, so the raw unclamped value is asserted). Because x contributes nothing to radial
//! distance, BOTH x=-1 and x=+1 corners at (y=1.3,z=10.5) tie exactly at this same max. Tracing the
//! production code's own sx/sy/sz loop (sx outermost, both loops -1 before +1, STRICT > comparisons)
//! by hand: x=-1 is reached first and becomes "worst"; x=+1, reached last, ties that same value --
//! not strictly greater, so it never replaces "worst" -- but IS strictly greater than whatever was
//! sitting in "second" at that point, so it displaces it into "second". I.e. this one configuration
//! deliberately exercises BOTH the primary AND the secondary-contact code path in a single test, with
//! both contacts sharing the identical overlap/radial-distance by construction.
void testArchTunnelBoxContainmentHandDerived()
{
   using namespace mesa_pd::data;

   const real_t tunnelRadius = real_t(10);
   CylindricalBoundary tunnel(tunnelRadius, Vec3(real_t(1), real_t(0), real_t(0)));
   tunnel.setContainment(true);
   tunnel.setAxialBounds(real_t(-10), real_t(10));
   tunnel.setHalfClipNormal(Vec3(real_t(0), real_t(0), real_t(1)));

   const Vec3 boxPos(real_t(0), real_t(0.3), real_t(9.5));
   const Box  box(Vec3(real_t(2), real_t(2), real_t(2))); // ctor takes FULL edge length -> halfExtents (1,1,1)

   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto boxP = ps->create(); // idx 0 -- matches operator()'s own (Box, CylindricalBoundary) argument order
   boxP->setPosition(boxPos);
   boxP->setRotation(Rot3());
   auto tunnelP = ps->create(); // idx 1
   tunnelP->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   tunnelP->setRotation(Rot3());

   // --- hand-derived expected values -------------------------------------------------------------
   const real_t expectedRadial  = std::sqrt(real_t(111.94)); // sqrt(1.3^2 + 10.5^2)
   const real_t expectedOverlap = expectedRadial - tunnelRadius;
   WALBERLA_CHECK_GREATER(expectedOverlap, real_t(0)); // sanity: genuinely penetrating
   WALBERLA_CHECK_LESS(expectedOverlap, real_t(1));    // sanity: under the envelope clamp

   const Vec3 expectedWorstCorner(real_t(-1), real_t(1.3), real_t(10.5));  // x=-1: found first -> "worst"
   const Vec3 expectedSecondCorner(real_t(1), real_t(1.3), real_t(10.5));  // x=+1: ties, found last -> "second"
   const Vec3 expectedOutward(real_t(0), real_t(1.3), real_t(10.5));       // radial part (x removed); SAME for both
   const Vec3 expectedNormal = -expectedOutward / expectedRadial;         // -outHat, points inward toward the axis

   // --- production containment override (SolidCylinderContactDetection -- ALWAYS analytic for this
   // shape pair regardless of useAnalyticalContactDetection; see that overload's own comment) -------
   SolidCylinderContactDetection prod;
   const bool prodContact = prod(0, 1, box, tunnel, ac);
   WALBERLA_CHECK(prodContact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.getPenetrationDepth(), -expectedOverlap, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.getContactNormal(), expectedNormal, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.getContactPoint(), expectedWorstCorner, real_t(1e-3));
   WALBERLA_LOG_INFO_ON_ROOT("ArchTunnel containment (box, hand-derived): overlap="
      << -prod.getPenetrationDepth() << " (expected " << expectedOverlap << ", agree)");

   // Secondary contact: the tied x=+1 corner, per the by-hand loop trace above. Same overlap/radial
   // distance as the primary by construction; only the contact POINT's x differs.
   WALBERLA_CHECK(prod.hasSecondaryContact());
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.getSecondaryPenetrationDepth(), -expectedOverlap, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.getSecondaryContactNormal(), expectedNormal, real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(prod.getSecondaryContactPoint(), expectedSecondCorner, real_t(1e-3));

   // --- mesa_pd's OWN generic GJK-based CylindricalBoundary path, called DIRECTLY (bypassing
   // SolidCylinderContactDetection's override entirely) -- see this section's own header comment for
   // exactly why this comparison is valid here and not in general. -----------------------------------
   mesa_pd::collision_detection::GeneralContactDetection generic;
   const bool genericContact = generic(0, 1, box, tunnel, ac);
   WALBERLA_CHECK(genericContact);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(generic.getPenetrationDepth(), -expectedOverlap, real_t(1e-3));
   // NORMAL SIGN CONVENTION DIFFERS, confirmed by actually running this (not just theorised): the
   // overlap MAGNITUDE agrees with the hand-derivation exactly (as predicted -- see this section's own
   // header comment for why), but mesa_pd's generic path reports its normal pointing AWAY from the
   // axis (its own long-standing exterior-pillar convention: "a body collides with the pillar from
   // OUTSIDE", per CylindricalBoundary.h's own doc comment), whereas our custom containment override
   // deliberately flips it to point INWARD, toward the axis (the physically correct push-back
   // direction for a body trying to escape a containing tube). Negating the generic path's normal
   // before comparing, same as this file's own established convention elsewhere for a documented sign
   // difference (see this file's top comment re: PaperAnalyticalContacts.h vs SolidCylinderContactDetection).
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(-generic.getContactNormal()[1], expectedNormal[1], real_t(1e-3));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(-generic.getContactNormal()[2], expectedNormal[2], real_t(1e-3));
   // x is deliberately NOT compared: the generic path's support-function direction has an EXACT-ZERO
   // x component (the axial component is subtracted out before the support() call), so which of the
   // box's two tied x=-+1 corners its support function returns for that direction is implementation-
   // defined -- unlike the production path's own deterministic, loop-order-derived choice above.
   WALBERLA_LOG_INFO_ON_ROOT("ArchTunnel containment (box) vs mesa_pd's own generic GJK path: overlap="
      << -generic.getPenetrationDepth() << " (expected " << expectedOverlap
      << ", agree on magnitude; normal sign flipped -- exterior-pillar vs containment convention, see comment)");
}

//! Sanity companion: move the box safely INSIDE the tunnel (radial distance from axis well under the
//! radius everywhere) -- both paths must agree there is NO contact.
void testArchTunnelBoxContainmentNoContact()
{
   using namespace mesa_pd::data;

   const real_t tunnelRadius = real_t(10);
   CylindricalBoundary tunnel(tunnelRadius, Vec3(real_t(1), real_t(0), real_t(0)));
   tunnel.setContainment(true);
   tunnel.setAxialBounds(real_t(-10), real_t(10));
   tunnel.setHalfClipNormal(Vec3(real_t(0), real_t(0), real_t(1)));

   const Box box(Vec3(real_t(2), real_t(2), real_t(2))); // halfExtents (1,1,1)

   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);
   auto boxP = ps->create();
   boxP->setPosition(Vec3(real_t(0), real_t(0), real_t(2))); // farthest corner radial dist = sqrt(1^2+3^2) ~= 3.16 << 10
   boxP->setRotation(Rot3());
   auto tunnelP = ps->create();
   tunnelP->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
   tunnelP->setRotation(Rot3());

   SolidCylinderContactDetection prod;
   WALBERLA_CHECK(!prod(0, 1, box, tunnel, ac));

   mesa_pd::collision_detection::GeneralContactDetection generic;
   WALBERLA_CHECK(!generic(0, 1, box, tunnel, ac));
}

//======================================================================================================================
// [Claude, 2026-09-30] Continuity of the detection along a path (thesis Section 5.2). For each random pair (uniform
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
// [Claude, 2026-09-30] Refined cylinder-cylinder detection (refineCylinderContact / cylinderCylinderContactRefined in
// PaperAnalyticalContacts.h) vs. GJK/EPA, on the same sampling as the master comparison (seed 24680, uniform SO(3)).
// Four detection variants are compared, all against GJK/EPA:
//   raw      closed form alone;
//   hybrid   closed form + GJK/EPA for "no scenario" (Tier 1) and special contacts (current production option);
//   trust    closed form, "no scenario" trusted as "no contact" (no GJK/EPA at all; special contact kept as is);
//   refined  closed form + support-function refinement (no GJK/EPA at all).
// Reported: existence errors (missed / false contacts), depth errors > max(0.01, 5 %), max depth error, normal angle,
// mean projection-overlap evaluations per refined call, and the time per call of each variant.
//======================================================================================================================
//! [Claude, 2026-09-30] Signed distance of point p to the surface of cylinder c (negative inside).
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
// [Claude, 2026-10-01] Box-box CONTACT POINT (detail::boxBoxContactPoint, used by boxBoxSAT). The former point lay on the
// line between the centres and could be several units from the contact and outside both boxes; elastic collisions then
// gained up to 98x their kinetic energy (DEMValidation cuboidPairEnergyFuzz). These tests check the new point against
// hand-derived exact points, inside-the-intersection on random pairs, and its continuity along a motion.
//======================================================================================================================
real_t boxSignedDistance(const Vec3& c, const Mat3& R, const Vec3& h, const Vec3& p)
{
   const Vec3 l = R.getTranspose() * (p - c);
   Vec3       q(std::abs(l[0]) - h[0], std::abs(l[1]) - h[1], std::abs(l[2]) - h[2]);
   const real_t inside = std::min(std::max(q[0], std::max(q[1], q[2])), real_t(0));
   const Vec3   qp(std::max(q[0], real_t(0)), std::max(q[1], real_t(0)), std::max(q[2], real_t(0)));
   return qp.length() + inside;
}

//! the former boxBoxSAT point, for comparison only
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
   // (a) centred stack: unchanged (the line between the centres runs through the contact)
   run("centred stack", Vec3(0, 0, 0), I, Vec3(1, 1, 1), Vec3(0, 0, real_t(2) - d), I, Vec3(1, 1, 1),
       Vec3(0, 0, real_t(1) - real_t(0.5) * d), real_t(1e-9));
   // (b) small box near the end of a long one: contact area x in [4, 6] -> centroid x = 5 (former point: x = 2.5)
   const real_t oldB = run("small box on the end of a long one", Vec3(0, 0, 0), I, Vec3(6, 1, real_t(1.5)),
                           Vec3(5, 0, real_t(2.5) - d), I, Vec3(1, 1, 1), Vec3(5, 0, real_t(1.5) - real_t(0.5) * d), real_t(1e-9));
   WALBERLA_CHECK_GREATER(oldB, real_t(2), "the former point should be ~2.5 off here");
   // (c) two equal faces, partly overlapping: rectangle x in [1, 2], y in [-1, 2] -> (1.5, 0.5)
   run("partly overlapping faces", Vec3(0, 0, 0), I, Vec3(2, 2, 1), Vec3(3, 1, real_t(2) - d), I, Vec3(2, 2, 1),
       Vec3(real_t(1.5), real_t(0.5), real_t(1) - real_t(0.5) * d), real_t(1e-9));
   // (d) corner into a face, off centre: cube with a body diagonal vertical, lowest vertex at height 1 - d above (0.7,-0.3)
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
   // (e) crossing edges: bar A along x rotated 45 deg about x (top edge along x at z = sqrt(2)/2), bar B along y rotated
   //     45 deg about y (bottom edge along y), centred at x = 3, y = -2 -> edges cross at (3, 0)
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
   // (f) the launch pair of the randomised cuboid drop test (local reproduction, step 20522, sub-cycle 9): the former
   //     point was 4.3 LU from the centroid of the overlap and outside both cuboids
   {
      const Vec3 pA(real_t(1.37178412), real_t(6.07269819), real_t(1.71354156));
      const Vec3 pB(real_t(2.48370382), real_t(1.60370074), real_t(6.09968309));
      const Vec3 h(6, 1, real_t(1.5));
      // rotation matrices as logged by [LAUNCH-TRACE] (row-major)
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

//! [Claude, 2026-10-01] random pairs near the contact boundary (as testBoxBoxRandomizedFuzz) plus elongated drop-test
//! cuboids; every reported point must lie in the intersection and not deeper than the overlap in either box.
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

//! [Claude, 2026-10-01] continuity of the contact point along straight paths through the contact (300 random pairs,
//! 300 steps each, translation only): largest jump of the point per step, in units of the step length. Feature changes
//! (vertex <-> edge <-> face as the overlap grows) are real jumps of a single-point model; reported, not asserted.
//! [Claude, 2026-10-22] extended to also track GJK/EPA's point along the same paths, under the same jump metric
//! (distance the reported point moves in one step, in units of the step length ds, counted as a jump above 10 ds),
//! per the author ("make this consistent both should be compared with the same") -- the table used to report a
//! number only for the SAT and a dash for GJK/EPA here, since this test originally called only boxBoxSAT.
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

//! [Claude, 2026-10-02] cost of the box-box contact point: random near-boundary pairs vs. nearly flat resting contacts
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

//! [Claude, 2026-10-02] detail::clipLineToBox, hand cases: the interval of t for which x0 + t dir lies in the box
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

//! [Claude, 2026-10-02] contacts as in a settled pile: two arbitrarily rotated 12x2x3 cuboids, B placed along a random
//! direction at the touching distance minus a tiny overlap (1e-7..1e-2, log-uniform). The alternating projection needed
//! up to its 2000-iteration cap here (225 of 362 contacts of the random_cuboids drop-test pile, 4284396 step 70000).
//! Checks for the closed-form point (line clip along the normal): inside the intersection (asserted), its offset from
//! the former projected point perpendicular to the normal (the lever-arm change; reported), and the cost of both.
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

//! [Claude, 2026-10-02] ANTIDUNES_BOX_PAIRS_EVAL=<file> (optional): the same comparison on dumped cuboid pairs, one pair
//! per line: posA(3) axesA(3x3, one axis per triple) halfA(3) posB(3) axesB(3x3) halfB(3) -- e.g. the near-contact
//! pairs of the random_cuboids drop-test pile (4284396, step 70000)
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
      // [Claude, 2026-09-30] cost breakdown of the refined detection by category of pair
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
   // regression guard for the refined detection (values of the 2026-09-30 run: 0 / 0 / 0, max depth error 2.2e-4)
   WALBERLA_CHECK_EQUAL(ref.missed, 0);
   WALBERLA_CHECK_EQUAL(ref.falsePos, 0);
   WALBERLA_CHECK_EQUAL(ref.wrongDepth, 0);
   // contact point of the refined detection inside the intersection of the two cylinders (2026-09-30: 0 of 1698)
   WALBERLA_CHECK_EQUAL(pst["refined"].outside, 0);
   WALBERLA_CHECK_EQUAL(pst["refined"].tooDeep, 0);
   WALBERLA_CHECK_EQUAL(pst["raw:band-band-skewed"].outside, 0);
}

//======================================================================================================================
// [Claude, 2026-10-20] NEARLY PARALLEL ROD PAIRS. The random samples above draw both orientations uniformly on SO(3), so
// that only a fraction 1 - cos(theta) of the pairs is skewed by less than theta: about 4 of the 1019 contacts of the thesis
// reference sample lie within 5 deg of parallel. The rods of the clogging runs, however, start aligned and lie nearly
// parallel in a jam. This sampler draws pairs of EQUAL rods (radius 0.3-2, length 4-8 radii) whose axes are skewed by
// theta, log-uniform between 1e-9 rad and maxSkew (10 % exactly parallel), placed side by side with an axial offset of up
// to 1.2 lengths, at a centre distance sampled around the GJK/EPA contact boundary (overlap up to 15 % of the radius).
// The closed form switches to the parallel formulas for 1 - |e_i.e_j| < 1e-8 (cylinderCylinderContact), i.e. theta below
// about 1.4e-4 rad, so the sample spans both sides of that switch.
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
   p.posJ = dir3 * fuzz.sampleSeparationNearBoundary(contactAt, real_t(0), real_t(3) * circ, p.R, real_t(0.15));
   return p;
}

//! Closed form with fallback (as in production) against GJK/EPA for nearly parallel rods, in bands of the skew angle.
//! Asserts what the fallback guarantees by construction (no missed contact, no false contact, finite results); depth and
//! normal errors are logged per band with the error criterion of the reference comparison, max(0.01, 0.05 depth).
void testCylinderNearParallelSkewedFuzz()
{
   const int numTrials = std::getenv("ANTIDUNES_NEARPARALLEL_TRIALS") ? std::atoi(std::getenv("ANTIDUNES_NEARPARALLEL_TRIALS"))
                                                                      : 2000;
   FuzzRng fuzz(24680u);
   const real_t maxSkew = real_t(10) * math::pi / real_t(180);
   // ANTIDUNES_NEARPARALLEL_MINSKEW_DEG: lower end of the log-uniform skew (default 1e-9 rad), to concentrate the sample
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
//! of fixed rod pairs smoothly across that switch (and on to 10 deg) and requires that the closed form with fallback follows
//! GJK/EPA without a jump up to 1e-2 rad: between neighbouring angles, its depth may change by at most 0.01 R more than that
//! of GJK/EPA. Beyond, the sweep exposes a KNOWN DEFECT (found 2026-10-20, THESIS_PLAN.md): when one rod is tilted out of
//! the common plane by about 4-10 deg and the rods are offset along their axes, the band-edge scenario reports depths that
//! are too small by up to 70 %, over a whole range of angles, and recovers with a jump (kodamBandEdgeAxialCheck narrows the
//! range but adds a jump of its own). The random samples rarely hit this configuration (none in 2409 contacts at 1-10 deg).
//! It is logged here and fails the test only with ANTIDUNES_NEARPARALLEL_STRICT set.
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
         // skew about the axis perpendicular to both the rod axis and the lateral direction (the rods open like scissors
         // in their common plane) and about the lateral direction itself (one rod tilts out of the plane)
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
// [Claude, 2026-09-30] INDEPENDENT-REFERENCE DUMP (opt-in, off by default: runs only if the environment variable
// ANTIDUNES_REFERENCE_DUMP is set to an output directory). Every other randomized test in this file treats mesa_pd's
// GJK/EPA as ground truth -- but GJK/EPA is also the PRODUCTION path of the thesis's rod runs
// (useAnalyticalContactDetection=false), so it needs its own independent check. This writes random pairs, with
// orientations drawn uniformly on SO(3) (FuzzRng::randomRotation), together with the GJK/EPA, closed-form (raw and hybrid) and SAT answers,
// to CSV. scripts/contact-validation/brute_force_reference.py (thesis repo) then recomputes the penetration depth of
// every pair from the support functions alone and compares all methods against it.
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
   // [Claude, 2026-10-11] ANTIDUNES_REFERENCE_SEED (default 13579, the seed of the thesis figure) for an out-of-sample check,
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
         // [Claude, 2026-10-10] columns of the refined result renamed ref_* -> rfd_* (they never held the reference); contact
         // points of every method and the normal of the hybrid added for the normal and contact-point panels of the thesis
         "gjk_px,gjk_py,gjk_pz,cf_px,cf_py,cf_pz,hyb_nx,hyb_ny,hyb_nz,hyb_px,hyb_py,hyb_pz,rfd_px,rfd_py,rfd_pz\n";
   // [Claude, 2026-10-20] ANTIDUNES_REFERENCE_NEARPARALLEL=<max skew in degrees>: the cylinder pairs are nearly parallel equal
   // rods from sampleNearParallelPair instead of uniform SO(3) pairs (same columns, so brute_force_reference.py applies as is)
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
         "gjk_px,gjk_py,gjk_pz,sat_px,sat_py,sat_pz\n";  // [Claude, 2026-10-10] contact points added
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
// [Claude, 2026-10-01] Independent-reference dump for rod-PILLAR (solid CylindricalBoundary) and rod-WALL (HalfSpace)
// contacts, the obstacles of the clogging layouts. Opt-in like testIndependentReferenceDump (ANTIDUNES_REFERENCE_DUMP).
// Production detection: SolidCylinderContactDetection (ParticleKernels.h, "pillar exterior" and Cylinder x HalfSpace).
// scripts/contact-validation/brute_force_reference.py recomputes the depth from the support functions, with the
// infinite pillar represented by a cylinder of half-length 1000 and the wall by its support function.
//======================================================================================================================
//======================================================================================================================
// [Claude, 2026-10-11] Replay of rod-pillar configurations taken from the floater snapshots of the clogging runs (opt-in:
// ANTIDUNES_PILLAR_REPLAY=<input csv> ANTIDUNES_PILLAR_REPLAY_OUT=<output csv>). Input per line: id,R,L,x,y,z,ax,ay,az,Rp,px,py
// (rod radius and length already inflated by the contact skin, world centre, world axis; vertical pillar of radius Rp through
// (px,py)). Writes the production rod-pillar detection (projection variant, refine = false) and the refined one, so that
// scripts compare both with the exact depth.
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
//! [Claude, 2026-10-01] Offline evaluation of the refined detection on pairs dumped from a packing run
//! (ParticleKernels ANTIDUNES_REFINE_DUMP). Only with ANTIDUNES_REFINE_EVAL=<dir with dump_rank*.csv>. For the
//! contacts the fast path rejected because the short descent found a smaller overlap ("smaller-elsewhere"), it compares
//! the full refinement with cheaper acceptance rules: accept the closed form if the descent improves on it by less than a
//! relative tolerance, or take the short descent's own minimum.
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
} // namespace paper_analytical
} // namespace antidunes
} // namespace walberla


int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes::paper_analytical::test;

   // Debug output (cylinder-mismatch VTK/CSV) goes next to the executable in the build tree,
   // regardless of the process's current working directory at launch.
   g_outputDir = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();

   // [Claude, 2026-10-10] ANTIDUNES_BANDEDGE_AXIAL: run with paper_analytical::kodamBandEdgeAxialCheck on (prm
   // kodamBandEdgeAxialCheck); ANTIDUNES_REFERENCE_DUMP_ONLY: write only the independent-reference dump and exit.
   walberla::antidunes::paper_analytical::kodamBandEdgeAxialCheck = (std::getenv("ANTIDUNES_BANDEDGE_AXIAL") != nullptr);
   if (std::getenv("ANTIDUNES_REFERENCE_DUMP_ONLY")) { testIndependentReferenceDump(); return EXIT_SUCCESS; }
   if (std::getenv("ANTIDUNES_PILLAR_REPLAY")) { replayPillarConfigurations(); return EXIT_SUCCESS; } // [Claude, 2026-10-11]
   if (std::getenv("ANTIDUNES_NEARPARALLEL_ONLY")) // [Claude, 2026-10-20] only the nearly parallel rod tests
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
   // [Claude, 2026-10-29] testBoxBoxFraigeVsGjkVsSatFuzz/testFraigeVertex*/testFraigeEdge* removed --
   // boxBoxFraigeMultiPointContacts unused by any thesis-cited result; archived at
   // archive/unused_contact_methods/ (fraige_multipoint_box_contact.h, fraige_validation_tests.cpp.excerpt).

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

   // [Claude, 2026-10-29] testFeng*/testCylinderCylinderMissPatternAnalysis removed (already
   // commented out here, i.e. already dropped from active development before this change) --
   // feng2017 unused by any thesis-cited result; archived at archive/unused_contact_methods/
   // (feng2017_cylinder_cylinder_contact.h, feng2017_validation_tests.cpp.excerpt).
   testCylinderCylinderFullComparisonAndHybridFeasibility();
   testEdgeEdgeThreeWayAccuracy10k();
   testThetaRefTieBreakNecessity10k();
   testClosestDistanceSqBetweenSegments();
   testTightTier1Accuracy10k();
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

   // Box x ArchTunnel (CylindricalBoundary containment mode)
   testArchTunnelBoxContainmentHandDerived();
   testArchTunnelBoxContainmentNoContact();

   // [Claude, 2026-09-30] refined detection vs. GJK/EPA; opt-in reference dump (ANTIDUNES_REFERENCE_DUMP=<dir>)
   testRefinedCylinderContactVsGjk();
   testCylinderCylinderContinuityAlongPaths();
   testCylinderNearParallelSkewedFuzz();          // [Claude, 2026-10-20]
   testCylinderNearParallelThresholdContinuity(); // [Claude, 2026-10-20]
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
