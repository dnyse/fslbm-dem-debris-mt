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
//! Unit tests for the analytic box/cylinder submerged-volume clip of ParticleKernels (clipBoxAnalytic,
//! clipCylinderAnalytic), plus a check of the spherical-cap formula used in computeSubmergedVolumes(). Most expected
//! values rest on geometric identities independent of the implementation (central symmetry, volume conservation,
//! textbook centroid formulas).
//!
//! The non-trivial cylinder cases (general tilt, position-varying and partial submersion) were cross-checked against
//! two methods that do not use clipCylinderAnalytic: a 60M-sample Monte Carlo integration classifying random points
//! against the plane, and trimesh slicing a 720-facet polygon mesh with the same plane. All three agree to ~1e-4,
//! the limit set by the polygon approximation of the circle; the box cases match trimesh to machine precision.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition; Antidunes.cpp gets it earlier via BuildTimeloop.h (free-surface dynamics
// headers). Included explicitly here.
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

// WALBERLA_CHECK_FLOAT_EQUAL(/_EPSILON) has a component-wise Vector3<real_t> specialization (core/math/Vector3.h),
// used directly on Vec3 values as in tests/mesa_pd/collision_detection/AnalyticContactDetection.cpp.

// ---------------------------------------------------------------------------------------------
// clipBoxAnalytic -- ground truth: central symmetry (through-centre plane -> exact half) and volume conservation
// (flipping the plane normal gives the exact complement), true for any correct clip.
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
   // Same box, arbitrary rotation R, nHat = R * e_z (cut axis-aligned in the body frame): checks that R is applied to
   // the input (R^T nHat) and to the output centroid offset; through-centre -> exact half for any orientation.
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
   // vSub(plane) + vSub(flipped normal, same point) == vFull for an off-centre, tilted plane (exact complements).
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
   // Partial submersion at several depths. For an axis-aligned box cut by a horizontal plane the cross-section is
   // constant with height, so the submerged volume is exactly linear in the plane offset.
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
// clipCylinderAnalytic -- one case per internal branch. Ground truth: the same symmetry/conservation identities as
// for the box, plus the textbook half-disk centroid 4r/(3*pi) for the vertical cut. The general-tilt case was
// cross-checked by the 60M-sample Monte Carlo integration (which reproduces the flat-cap and vertical-cut cases to
// ~1e-4): cbOff = (0,-0.14443,-0.45829), as asserted here.
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
   // Monte Carlo reference (file header); epsilon an order of magnitude above the MC uncertainty (~1e-4).
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
   // Partial submersion at several depths. The cylinder's cross-section curves, so the volume is not linear in
   // depth; expected values from the same Monte Carlo integration as the general-tilt case.
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
   // Same branch, offset -0.4 along nHat -> MC: 25.23% submerged (complement of the above; checks the sign convention).
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
// Sphere submerged-cap formula: calls sphereCapVolume() (ParticleKernels.h, used by computeSubmergedVolumes())
// and checks it against known volumes (empty, hemisphere, full sphere). The clamp of hSigned to hCap in [0,2r]
// mirrors computeSubmergedVolumes(); it is the caller's job, not part of the formula under test.
// ---------------------------------------------------------------------------------------------

void testSphereCapFormula()
{
   const real_t r = real_t(1);
   const real_t fullSphereVol = math::pi * real_t(4) / real_t(3) * r * r * r;
   auto capVolume = [&](real_t hSigned) {
      const real_t hCap = std::clamp(r + hSigned, real_t(0), real_t(2) * r);
      return sphereCapVolume(r, hCap);
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
