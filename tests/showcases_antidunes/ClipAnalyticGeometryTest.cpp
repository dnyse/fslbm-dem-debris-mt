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
//! \file ClipAnalyticGeometryTest.cpp
//! Unit tests for the analytic box/cylinder submerged-volume clip used by the Antidunes
//! ParticleKernels (clipBoxAnalytic, clipCylinderAnalytic), plus a sanity check of the documented
//! spherical-cap formula used inline in ParticleKernels::computeSubmergedVolumes(). See
//! /home/dhuber/.claude/plans/shiny-toasting-wind.md for the ground-truth derivation of every
//! expected value used below -- most rest on general geometric identities (central symmetry,
//! volume conservation, textbook centroid formulas) independent of this implementation.
//!
//! The non-trivial cylinder cases (general tilt, and all the position-varying/partial-submersion
//! cases) were additionally cross-checked against TWO independent methods, neither of which calls
//! clipCylinderAnalytic or re-derives its closed-form antiderivatives: (1) a 60M-sample Monte
//! Carlo integration (numpy, scratchpad script) sampling random points and classifying them
//! directly against the plane equation, and (2) trimesh (an established, separately-developed
//! mesh-processing library) actually slicing a 720-facet polygon mesh with the same plane and
//! integrating the resulting solid's volume/centroid directly. All three sources (this code, the
//! Monte Carlo integration, and trimesh's exact mesh-slicing) agree to ~1e-4, the expected
//! precision limit from trimesh's polygon approximation of the true circular cross-section -- the
//! box cases matched trimesh's exact (non-approximated) polygon-slicing to machine precision.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition itself -- it only compiles in the real Antidunes.cpp because
// BuildTimeloop.h happens to pull it in first via the free-surface dynamics headers. Include it
// explicitly here too (test-side only, no production file touched) so this translation unit doesn't
// depend on that same include-order accident.
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h" // pulls in mesa_pd::Mat3, Vector3<real_t>, math::pi, and the clip*Analytic() functions

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/mpi/MPIManager.h>

#include <cmath>
#include <cstdlib>

namespace walberla {
namespace antidunes {

using Vec3 = Vector3< real_t >;

// Note: WALBERLA_CHECK_FLOAT_EQUAL(/_EPSILON) has a Vector3<real_t> specialization
// (core/math/Vector3.h) that compares component-wise, so it can be used directly on Vec3 values
// below -- same convention as tests/mesa_pd/collision_detection/AnalyticContactDetection.cpp.

// ---------------------------------------------------------------------------------------------
// clipBoxAnalytic -- ground truth: central-symmetry theorem (through-centre -> exact half) and
// the volume-conservation identity (flipping the plane normal gives the exact complement). Both
// are geometric facts true for ANY correct clip implementation, independent of clipBoxAnalytic's
// own internals.
// ---------------------------------------------------------------------------------------------

void testBoxThroughCenterAxisAligned()
{
   const Vec3 h(real_t(0.5), real_t(1.0), real_t(2.0)); // half-extents; asymmetric so x/y/z aren't interchangeable
   const real_t vFull = real_t(8) * h[0] * h[1] * h[2];
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat(real_t(0), real_t(0), real_t(1));
   const Vec3 xBar = pos; // plane through the box centre

   real_t vSub; Vec3 cbOff;
   clipBoxAnalytic(h, R, nHat, xBar, pos, vSub, cbOff);

   WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull * real_t(0.5));
   WALBERLA_CHECK_FLOAT_EQUAL(cbOff, Vec3(real_t(0), real_t(0), -h[2] * real_t(0.5)));
}

void testBoxThroughCenterArbitraryRotation()
{
   // Same box, arbitrary rotation R, with nHat = R * e_z (the cut stays "axis-aligned" in the
   // box's OWN body frame) -- exercises that R is correctly applied both to the input (via R^T
   // acting on nHat internally) and to the output centroid offset, while still resting on the
   // same central-symmetry identity (through-centre -> exact half, for ANY orientation).
   const Vec3 h(real_t(0.5), real_t(1.0), real_t(2.0));
   const real_t vFull = real_t(8) * h[0] * h[1] * h[2];
   const mesa_pd::Rot3 rot(Vec3(real_t(1), real_t(1), real_t(1)).getNormalized(), real_t(0.9));
   const mesa_pd::Mat3 R = rot.getMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat = R * Vec3(real_t(0), real_t(0), real_t(1));
   const Vec3 xBar = pos;

   real_t vSub; Vec3 cbOff;
   clipBoxAnalytic(h, R, nHat, xBar, pos, vSub, cbOff);

   WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull * real_t(0.5));
   const Vec3 expectedCbOff = R * Vec3(real_t(0), real_t(0), -h[2] * real_t(0.5));
   WALBERLA_CHECK_FLOAT_EQUAL(cbOff, expectedCbOff);
}

void testBoxFullyAboveAndFullyBelow()
{
   const Vec3 h(real_t(0.5), real_t(1.0), real_t(2.0));
   const real_t vFull = real_t(8) * h[0] * h[1] * h[2];
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat(real_t(0), real_t(0), real_t(1));

   // Plane far below the box's bottom face (box entirely on the gas side) -> vSub = 0.
   {
      real_t vSub; Vec3 cbOff;
      clipBoxAnalytic(h, R, nHat, Vec3(real_t(0), real_t(0), real_t(-3)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, real_t(0));
   }
   // Plane far above the box's top face (box entirely submerged) -> vSub = vFull, cbOff = 0.
   {
      real_t vSub; Vec3 cbOff;
      clipBoxAnalytic(h, R, nHat, Vec3(real_t(0), real_t(0), real_t(3)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull);
      WALBERLA_CHECK_FLOAT_EQUAL(cbOff, Vec3(real_t(0)));
   }
}

void testBoxConservationIdentity()
{
   // vSub(plane) + vSub(flipped normal, same point) == vFull for an off-centre, tilted plane:
   // the two clips are, by construction, exact geometric complements of the same solid, for ANY
   // plane -- an algebraic identity, not something specific to this implementation.
   const Vec3 h(real_t(0.7), real_t(0.9), real_t(1.3));
   const real_t vFull = real_t(8) * h[0] * h[1] * h[2];
   const mesa_pd::Rot3 rot(Vec3(real_t(1), real_t(2), real_t(-1)).getNormalized(), real_t(0.9));
   const mesa_pd::Mat3 R = rot.getMatrix();
   const Vec3 pos(real_t(1), real_t(-2), real_t(0.5));
   const Vec3 nHat = Vec3(real_t(0.3), real_t(-0.6), real_t(0.7)).getNormalized();
   const Vec3 xBar = pos + Vec3(real_t(0.05), real_t(-0.1), real_t(0.2));

   real_t vSub1, vSub2; Vec3 cbOff1, cbOff2;
   clipBoxAnalytic(h, R, nHat, xBar, pos, vSub1, cbOff1);
   clipBoxAnalytic(h, R, -nHat, xBar, pos, vSub2, cbOff2);

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub1 + vSub2, vFull, real_t(1e-6));
   // First moments about the box centre must cancel too (the full box's own centroid is its centre).
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub1 * cbOff1 + vSub2 * cbOff2, Vec3(real_t(0)), real_t(1e-6));
}

void testBoxPartialSubmersion()
{
   // The floater at different depths relative to the water, not just through-centre/fully-wet/
   // fully-dry. For an axis-aligned box cut by a horizontal plane, the cross-section area is
   // constant with height, so the submerged volume is an EXACT linear function of the plane
   // offset -- a closed form independent of clipBoxAnalytic's own internals, not just another
   // self-consistency check.
   const Vec3 h(real_t(0.5), real_t(1.0), real_t(2.0)); // box spans z in [-2,2]
   const real_t vFull = real_t(8) * h[0] * h[1] * h[2];
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat(real_t(0), real_t(0), real_t(1));

   // Plane at z=-1 -> submerged slab z in [-2,-1], hcap=1 -> 1/(2*hz=4) = 25% submerged.
   {
      real_t vSub; Vec3 cbOff;
      clipBoxAnalytic(h, R, nHat, Vec3(real_t(0), real_t(0), real_t(-1)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, real_t(2));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub / vFull, real_t(0.25), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL(cbOff, Vec3(real_t(0), real_t(0), real_t(-1.5)));
   }
   // Plane at z=+1 -> submerged slab z in [-2,1], hcap=3 -> 3/4 = 75% submerged.
   {
      real_t vSub; Vec3 cbOff;
      clipBoxAnalytic(h, R, nHat, Vec3(real_t(0), real_t(0), real_t(1)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, real_t(6));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub / vFull, real_t(0.75), real_t(1e-6));
      WALBERLA_CHECK_FLOAT_EQUAL(cbOff, Vec3(real_t(0), real_t(0), real_t(-0.5)));
   }
}

// ---------------------------------------------------------------------------------------------
// clipCylinderAnalytic -- one case per internal branch. Ground truth: same central-symmetry /
// conservation identities as the box for the through-centre and flat-cap/vertical-cut cases, PLUS
// the textbook half-disk centroid 4r/(3*pi) for the vertical-cut case (an external reference, not
// derived from this code). The general-tilt case was additionally cross-checked against an
// independent 60M-sample Monte Carlo integration (see the plan file) -- that MC run first
// reproduced the flat-cap and vertical-cut cases to ~1e-4 (validating the method itself), then
// gave cbOff=(0,-0.14443,-0.45829) for the general-tilt case, matching what's asserted here.
// ---------------------------------------------------------------------------------------------

void testCylinderFlatCapBranch()
{
   // rho ~ 0: nHat parallel to the body z (cylinder) axis.
   const real_t r = real_t(1), halfLen = real_t(1);
   const real_t vFull = math::pi * r * r * real_t(2) * halfLen; // 2*pi
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat(real_t(0), real_t(0), real_t(1));

   // Through centre: half volume, offset purely axial = -halfLen/2, zero radial component.
   {
      real_t vSub; Vec3 cbOff;
      clipCylinderAnalytic(r, halfLen, R, nHat, pos, pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull * real_t(0.5));
      WALBERLA_CHECK_FLOAT_EQUAL(cbOff, Vec3(real_t(0), real_t(0), real_t(-0.5)));
   }
   // Fully dry (plane below the bottom cap).
   {
      real_t vSub; Vec3 cbOff;
      clipCylinderAnalytic(r, halfLen, R, nHat, Vec3(real_t(0), real_t(0), real_t(-2)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, real_t(0));
   }
   // Fully submerged (plane above the top cap).
   {
      real_t vSub; Vec3 cbOff;
      clipCylinderAnalytic(r, halfLen, R, nHat, Vec3(real_t(0), real_t(0), real_t(2)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull);
   }
}

void testCylinderVerticalCutBranch()
{
   // nz ~ 0: nHat perpendicular to the axis -> classic circular-segment / half-disk geometry.
   const real_t r = real_t(1), halfLen = real_t(1);
   const real_t vFull = math::pi * r * r * real_t(2) * halfLen;
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat(real_t(1), real_t(0), real_t(0));

   real_t vSub; Vec3 cbOff;
   clipCylinderAnalytic(r, halfLen, R, nHat, pos, pos, vSub, cbOff);

   WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull * real_t(0.5));
   // Textbook half-disk centroid distance from the flat edge: 4r/(3*pi); on the -x (submerged) side.
   const real_t expectedCx = -real_t(4) * r / (real_t(3) * math::pi);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cbOff, Vec3(expectedCx, real_t(0), real_t(0)), real_t(1e-6));
}

void testCylinderGeneralTiltBranch()
{
   // Both rho and nz non-negligible: nHat tilted 30 deg from the axis, plane through centre.
   // Independently verified via Monte Carlo (see file header + plan) -- epsilon set an order of
   // magnitude looser than the MC statistical uncertainty (~1e-4), so this only fails on a real
   // implementation change, not MC noise.
   const real_t r = real_t(1), halfLen = real_t(1);
   const real_t vFull = math::pi * r * r * real_t(2) * halfLen;
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));
   const Vec3 nHat(real_t(0), real_t(0.5), std::sqrt(real_t(0.75))); // unit vector, 30 deg tilt from z

   real_t vSub; Vec3 cbOff;
   clipCylinderAnalytic(r, halfLen, R, nHat, pos, pos, vSub, cbOff);

   // Through-centre identity (finite cylinder is centrally symmetric about its own centre).
   WALBERLA_CHECK_FLOAT_EQUAL(vSub, vFull * real_t(0.5));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cbOff, Vec3(real_t(0), real_t(-0.14443), real_t(-0.45829)), real_t(1e-3));
}

void testCylinderPartialSubmersion()
{
   // The floater at different depths relative to the water, not just through-centre. Unlike the
   // box, a cylinder's cross-section curves, so these offsets are NOT a linear function of depth
   // -- independently verified via the same 60M-sample Monte Carlo integration as the general-tilt
   // case (see the plan file), not re-derived from clipCylinderAnalytic's own formulas.
   const real_t r = real_t(1), halfLen = real_t(1);
   const real_t vFull = math::pi * r * r * real_t(2) * halfLen;
   const mesa_pd::Mat3 R = mesa_pd::Mat3::makeIdentityMatrix();
   const Vec3 pos(real_t(0), real_t(0), real_t(0));

   // Vertical-cut branch (nHat perpendicular to axis), offset +0.4 along nHat -> MC: 74.77% submerged.
   {
      const Vec3 nHat(real_t(1), real_t(0), real_t(0));
      real_t vSub; Vec3 cbOff;
      clipCylinderAnalytic(r, halfLen, R, nHat, Vec3(real_t(0.4), real_t(0), real_t(0)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub / vFull, real_t(0.74773), real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cbOff, Vec3(real_t(-0.2185), real_t(0), real_t(0)), real_t(1e-3));
   }
   // Same branch, offset -0.4 along nHat -> MC: 25.23% submerged (the geometric complement of the
   // above, roughly -- confirms the offset direction/sign convention, not just the magnitude).
   {
      const Vec3 nHat(real_t(1), real_t(0), real_t(0));
      real_t vSub; Vec3 cbOff;
      clipCylinderAnalytic(r, halfLen, R, nHat, Vec3(real_t(-0.4), real_t(0), real_t(0)), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub / vFull, real_t(0.25230), real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cbOff, Vec3(real_t(-0.6474), real_t(0), real_t(0)), real_t(1e-3));
   }
   // General-tilt branch (30 deg from axis), offset +0.3 along nHat -> MC: 67.32% submerged.
   {
      const Vec3 nHat(real_t(0), real_t(0.5), std::sqrt(real_t(0.75)));
      real_t vSub; Vec3 cbOff;
      clipCylinderAnalytic(r, halfLen, R, nHat, nHat * real_t(0.3), pos, vSub, cbOff);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub / vFull, real_t(0.67322), real_t(1e-3));
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(cbOff, Vec3(real_t(0), real_t(-0.10722), real_t(-0.29587)), real_t(1e-3));
   }
}

void testCylinderConservationIdentity()
{
   const real_t r = real_t(1.3), halfLen = real_t(0.8);
   const real_t vFull = math::pi * r * r * real_t(2) * halfLen;
   const mesa_pd::Rot3 rot(Vec3(real_t(0), real_t(1), real_t(1)).getNormalized(), real_t(0.4));
   const mesa_pd::Mat3 R = rot.getMatrix();
   const Vec3 pos(real_t(-0.3), real_t(0.2), real_t(0.1));
   const Vec3 nHat = Vec3(real_t(0.4), real_t(0.2), real_t(0.5)).getNormalized();
   const Vec3 xBar = pos + Vec3(real_t(0.1), real_t(0), real_t(-0.05));

   real_t vSub1, vSub2; Vec3 cbOff1, cbOff2;
   clipCylinderAnalytic(r, halfLen, R, nHat, xBar, pos, vSub1, cbOff1);
   clipCylinderAnalytic(r, halfLen, R, -nHat, xBar, pos, vSub2, cbOff2);

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(vSub1 + vSub2, vFull, real_t(1e-6));
}

// ---------------------------------------------------------------------------------------------
// Sphere submerged-cap formula (documented closed form, ParticleKernels.h ~line 1551) -- sanity
// checks against INDEPENDENTLY known volumes (hemisphere, full sphere, empty), not a re-derivation
// of the formula itself. This is pre-existing math, not new this project, so it's a lighter check.
// ---------------------------------------------------------------------------------------------

void testSphereCapFormula()
{
   const real_t r = real_t(1);
   const real_t fullSphereVol = math::pi * real_t(4) / real_t(3) * r * r * r;
   auto capVolume = [&](real_t hSigned) {
      const real_t hCap = std::clamp(r + hSigned, real_t(0), real_t(2) * r);
      return math::pi * hCap * hCap * (real_t(3) * r - hCap) / real_t(3);
   };

   WALBERLA_CHECK_FLOAT_EQUAL(capVolume(real_t(-1.5)), real_t(0));                 // fully above -> empty
   WALBERLA_CHECK_FLOAT_EQUAL(capVolume(real_t(1.5)), fullSphereVol);              // fully below -> full sphere
   WALBERLA_CHECK_FLOAT_EQUAL(capVolume(real_t(0)), fullSphereVol * real_t(0.5));  // through centre -> hemisphere
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes;
   testBoxThroughCenterAxisAligned();
   testBoxThroughCenterArbitraryRotation();
   testBoxFullyAboveAndFullyBelow();
   testBoxConservationIdentity();
   testBoxPartialSubmersion();
   testCylinderFlatCapBranch();
   testCylinderVerticalCutBranch();
   testCylinderGeneralTiltBranch();
   testCylinderPartialSubmersion();
   testCylinderConservationIdentity();
   testSphereCapFormula();

   return EXIT_SUCCESS;
}
