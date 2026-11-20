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
//! \file ParticleKernels.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "AntidunesBoundaryHandling.h"

#include "core/timing/TimingPool.h"

#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/collision_detection/GeneralContactDetection.h"
#include "mesa_pd/data/shape/Cylinder.h"
#include "mesa_pd/data/shape/CylindricalBoundary.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"

#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "field/FlagField.h"
#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/shape/Box.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"

#include "ClosedFormContacts.h"
#include "SimulationParams.h"
#include "SimulationObjects.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <numeric>
#include <filesystem>
#include <fstream>
#include <memory>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>


namespace walberla {
namespace antidunes {

// Shared geometry for CylindricalBoundary CONTAINMENT mode (CylindricalBoundary::isContainment()). The
// Sphere/Box/Cylinder x CylindricalBoundary overloads below use it to find a body's FARTHEST point from the
// tube axis (mirror of the pillar overloads' closest-point search), so a body poking outside the tube is
// pushed back toward the axis.
namespace containment_detail {

// Result of circleFarthestFromAxisClipped; valid == false if no point of the circle lies on the solid side.
struct CircleFarthestResult
{
   mesa_pd::Vec3 point;
   real_t        radialDistSq = real_t(0);
   bool          valid        = false;
};

// Farthest point (from a target axis line) on a circle (circleCenter, radius, spanned by orthonormal u/v),
// restricted to the solid half of the plane through clipOrigin with normal clipNormal (zero clipNormal = no
// restriction, i.e. CylindricalBoundary::hasHalfClip()==false). Uses the closed-form solver
// farthestThetaOnCircleToLine and the closed-form circle/plane crossing (amplitude/phase).
inline CircleFarthestResult circleFarthestFromAxisClipped(
   const mesa_pd::Vec3& circleCenter, const mesa_pd::Vec3& u, const mesa_pd::Vec3& v, real_t circleRadius,
   const mesa_pd::Vec3& axisPoint, const mesa_pd::Vec3& axisDir,
   const mesa_pd::Vec3& clipOrigin, const mesa_pd::Vec3& clipNormal)
{
   auto pointAt = [&](real_t theta) {
      return circleCenter + circleRadius * std::cos(theta) * u + circleRadius * std::sin(theta) * v;
   };
   auto radialSq = [&](const mesa_pd::Vec3& P) {
      const mesa_pd::Vec3 rel    = P - axisPoint;
      const mesa_pd::Vec3 radial = rel - (rel * axisDir) * axisDir;
      return radial.sqrLength();
   };
   const bool hasClip = clipNormal.sqrLength() > real_t(1e-20);
   auto onSolidSide = [&](const mesa_pd::Vec3& P) {
      return !hasClip || (P - clipOrigin) * clipNormal >= real_t(0);
   };

   const real_t thetaBest = closed_form_contacts::detail::farthestThetaOnCircleToLine(
      circleCenter, u, v, circleRadius, axisPoint, axisDir);
   const mesa_pd::Vec3 candidate = pointAt(thetaBest);

   if (!hasClip || onSolidSide(candidate)) return { candidate, radialSq(candidate), true };

// Unconstrained farthest point is on the clipped-away side: restrict to the (up to two) points where the
// circle crosses the clip plane.
   const real_t A1 = u * clipNormal, A2 = v * clipNormal;
   const real_t Amag = std::sqrt(A1 * A1 + A2 * A2);
   if (Amag < real_t(1e-12)) return { candidate, real_t(0), false }; // circle plane || clip plane: all on wrong side

   const real_t rhs = ((clipOrigin - circleCenter) * clipNormal) / circleRadius;
   if (std::abs(rhs) > Amag) return { candidate, real_t(0), false }; // whole circle on the wrong side

   const real_t phi   = std::atan2(A2, A1);
   const real_t delta = std::acos(std::clamp(rhs / Amag, real_t(-1), real_t(1)));
   const mesa_pd::Vec3 e1 = pointAt(phi + delta);
   const mesa_pd::Vec3 e2 = pointAt(phi - delta);
   const real_t r1 = radialSq(e1), r2 = radialSq(e2);
   return (r1 > r2) ? CircleFarthestResult{ e1, r1, true } : CircleFarthestResult{ e2, r2, true };
}

// Axial coordinate of a world point relative to the CylindricalBoundary position (the reference of
// CylindricalBoundary::getAxialMin()/getAxialMax()).
inline real_t axialCoord(const mesa_pd::Vec3& worldPoint, const mesa_pd::Vec3& cbPos, const mesa_pd::Vec3& axis)
{
   return (worldPoint - cbPos) * axis;
}

} // namespace containment_detail

// Hides the hollow-cylinder overload of AnalyticContactDetection with exterior (solid cylinder) contact math;
// DoubleCast deduces the concrete functor type, so overload resolution picks the methods here.
//
// Base class is GeneralContactDetection: pairs without an overload here (e.g. Cylinder x Box, Cylinder x
// Ellipsoid; Box x Box and Cylinder x Cylinder by default, see their overloads) use its GJK/EPA path, which
// uses the same contactNormal/penetrationDepth convention as the analytic functions (asserted in
// tests/mesa_pd/collision_detection/EPA.cpp). Analytic overloads (sphere-sphere, sphere-box, sphere-halfspace)
// still win overload resolution. Cylinder x HalfSpace has its own overload below because the base's generic
// HalfSpace template returns contactNormal = -wallNormal, opposite to the Sphere/Box HalfSpace convention
// (+wallNormal, out of the wall); see ContactDetectionTest.cpp testCylinderHalfSpaceContact.
class SolidCylinderContactDetection : public mesa_pd::collision_detection::GeneralContactDetection
{
public:
   using GeneralContactDetection::operator(); // prevent name hiding of all base overloads

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Sphere& s,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   Accessor& ac)
   {
      getIdx1() = idx1;
      getIdx2() = idx2;

      const auto& posSphere = ac.getPosition(idx1);
      const auto& posCyl    = ac.getPosition(idx2);
      const auto& axis      = cb.getAxis();

      if (cb.isContainment())
      {
         // Envelope-safety cap: bounds the force a single degenerate/runaway overlap can inject.
         const real_t maxOverlap = real_t(1.0);

         // Axial gate: skip if the sphere's axial extent misses the tube's span (vacuous for unbounded
         // axialMin/axialMax).
         const real_t axialC = containment_detail::axialCoord(posSphere, posCyl, axis);
         if (axialC + s.getRadius() < cb.getAxialMin() || axialC - s.getRadius() > cb.getAxialMax())
            return false;

         // Radial pre-gate: several tubes can share the same axial span (e.g. the 3 arches of a bridge), so
         // a body near another tube must not get a spurious huge overlap against this one. dist - effR is a
         // lower bound on the body's minimum distance to the axis for any shape.
         {
            const mesa_pd::Vec3 relCenter  = posSphere - posCyl;
            const real_t        distCenter = (relCenter - (relCenter * axis) * axis).length();
            if (distCenter - s.getRadius() > cb.getRadius()) return false;
         }

         mesa_pd::Vec3 u, v;
         closed_form_contacts::detail::perpendicularBasis(axis, u, v);
         const auto far = containment_detail::circleFarthestFromAxisClipped(
            posSphere, u, v, s.getRadius(), posCyl, axis, posCyl, cb.getHalfClipNormal());
         if (!far.valid) return false; // sphere entirely on the non-solid (clipped-away) side

         real_t farthestRadial = std::sqrt(far.radialDistSq);
         real_t overlap = farthestRadial - cb.getRadius(); // > 0 -> pokes outside the tube
         if (overlap <= real_t(0)) return false;           // safely inside

         // SAFETY: clamp the overlap (as the Box/Cylinder pillar-exterior overloads do) so a body wildly
         // outside the tube (already-corrupted state) cannot feed an unbounded spring force into the sim.
         const bool envelopClamped = overlap > maxOverlap;
         if (envelopClamped)
         {
            static uint_t archEnvelopLogCount = uint_t(0);
            if (archEnvelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[ARCH-ENVELOP] sphere uid=" << ac.getUid(idx1)
                  << " deeply overlaps/envelops arch tunnel uid=" << ac.getUid(idx2)
                  << " (raw overlap=" << overlap << ", farthest-point radial dist=" << farthestRadial
                  << ") — clamping overlap to " << maxOverlap << ".");
               ++archEnvelopLogCount;
            }
            overlap = maxOverlap;
            farthestRadial = cb.getRadius() + maxOverlap;
         }

         const mesa_pd::Vec3 rel     = far.point - posCyl;
         const mesa_pd::Vec3 outward = rel - (rel * axis) * axis;
         const real_t        outDist = outward.length();
         const mesa_pd::Vec3 outHat  = (outDist > real_t(1e-10))
                                          ? outward / outDist
                                          : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));

         getPenetrationDepth() = -overlap;  // pokes outside -> penetrating (this class's convention)
         getContactNormal()    = -outHat;   // points inward, toward the axis
         getContactPoint()     = far.point;
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "arch/tunnel containment (Sphere, ENVELOP-clamped)" : "arch/tunnel containment (Sphere)");
         return true;
      }

      // Radial vector from cylinder axis to sphere center (component perpendicular to axis)
      const mesa_pd::Vec3 delta  = posSphere - posCyl;
      const mesa_pd::Vec3 radial = delta - (delta * axis) * axis;
      const real_t dist          = radial.length();

      const real_t separation = dist - cb.getRadius() - s.getRadius();
      if (separation < getContactThreshold())
      {
         getPenetrationDepth() = separation;
         // Normal points radially outward from the cylinder axis toward the sphere
         getContactNormal() = (dist > real_t(1e-10))
                                 ? radial / dist
                                 : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         getContactPoint() = posSphere - (s.getRadius() + getPenetrationDepth()) * getContactNormal();
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), "pillar exterior (Sphere)");
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   const mesa_pd::data::Sphere& s,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, s, cb, ac);
   }

   // Box vs HalfSpace, rotation-aware: the wall normal is transformed to the body frame (R^T) for support()
   // and the result rotated back. Same convention as Sphere-HalfSpace:
   //   penetrationDepth < 0  -> penetrating (contact active)
   //   contactNormal = hs.getNormal() (points into the domain, away from the wall)
   //   contactPoint on the wall surface (deepest box corner projected onto the plane)
   // boxContactSkin (ModelParameters::boxContactSkin): detect with the inflated box (each half-extent + m),
   // as cylinderContactSkin does in the Cylinder x HalfSpace overload below.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Box& bx,
                   const mesa_pd::data::HalfSpace& hs,
                   Accessor& ac)
   {
      if (boxContactSkin > real_t(0))
      {
         const real_t m = boxContactSkin;
         const mesa_pd::data::Box inflated(bx.getEdgeLength() + mesa_pd::Vec3(real_t(2) * m));
         boxContactSkin = real_t(0);
         const bool hit = operator()(idx1, idx2, inflated, hs, ac);
         boxContactSkin = m;
         return hit;
      }
      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Vec3& c = ac.getPosition(idx1);
      const mesa_pd::Vec3& n = hs.getNormal();
      const real_t d = n * ac.getPosition(idx2); // plane displacement: n · pos_plane

      // Rotation-aware deepest corner: transform wall normal to body frame, find the
      // support corner there, then rotate back to world frame.
      const mesa_pd::Mat3& R    = ac.getRotation(idx1).getMatrix();
      const mesa_pd::Vec3 deepest = c + R * bx.support(-(R.getTranspose() * n));
      const real_t penetrationDepth = n * deepest - d;

      if (penetrationDepth < getContactThreshold())
      {
         getPenetrationDepth() = penetrationDepth;
         getContactNormal()    = n;
         getContactPoint()     = deepest - penetrationDepth * n;
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::HalfSpace& hs,
                   const mesa_pd::data::Box& bx,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, bx, hs, ac);
   }

   // Cylinder vs HalfSpace: same formula/convention as Box vs HalfSpace (Cylinder has an exact support()).
   // Explicit override because the base's generic HalfSpace template returns contactNormal = -wallNormal
   // (see the class comment and ContactDetectionTest.cpp testCylinderHalfSpaceContact).
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Cylinder& cyl,
                   const mesa_pd::data::HalfSpace& hs,
                   Accessor& ac)
   {
      // cylinderContactSkin (see ModelParameters): detect with the inflated cylinder (radius + m, half-length + m).
      if (cylinderContactSkin > real_t(0))
      {
         const real_t m = cylinderContactSkin;
         const mesa_pd::data::Cylinder inflated(cyl.getRadius() + m, cyl.getHeight() + real_t(2) * m);
         cylinderContactSkin = real_t(0);
         const bool hit = operator()(idx1, idx2, inflated, hs, ac);
         cylinderContactSkin = m;
         return hit;
      }
      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Vec3& c = ac.getPosition(idx1);
      const mesa_pd::Vec3& n = hs.getNormal();
      const real_t d = n * ac.getPosition(idx2); // plane displacement: n · pos_plane

      const mesa_pd::Mat3& R      = ac.getRotation(idx1).getMatrix();
      const mesa_pd::Vec3 deepest = c + R * cyl.support(-(R.getTranspose() * n));
      const real_t penetrationDepth = n * deepest - d;

      if (penetrationDepth < getContactThreshold())
      {
         getPenetrationDepth() = penetrationDepth;
         getContactNormal()    = n;
         getContactPoint()     = deepest - penetrationDepth * n;
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::HalfSpace& hs,
                   const mesa_pd::data::Cylinder& cyl,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, cyl, hs, ac);
   }

   // Box vs CylindricalBoundary (solid exterior), rotation-aware: works in the box body frame. The cylinder
   // centre/axis are transformed by R^T, the closest box point to the (infinite) axis is found by clamping to
   // the half-extents, the radial distance gives the penetration, and normal/contact point are rotated back.
   //   penetrationDepth = radial_dist - cylinder_radius  (< 0 -> penetrating)
   //   contactNormal points radially outward from the cylinder axis toward the box (world frame)
   //   contactPoint = closest point on the box surface, world frame
   // boxContactSkin (ModelParameters::boxContactSkin): detect with the inflated box (each half-extent + m);
   // containment mode (arch/tunnel) is not inflated, as in the Cylinder x CylindricalBoundary overload below.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Box& bx,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   Accessor& ac)
   {
      const real_t m = boxContactSkin;
      if (m <= real_t(0) || cb.isContainment()) return boxPillarCore(idx1, idx2, bx, cb, ac);
      const mesa_pd::data::Box inflated(bx.getEdgeLength() + mesa_pd::Vec3(real_t(2) * m));
      return boxPillarCore(idx1, idx2, inflated, cb, ac);
   }

   template< typename Accessor >
   bool boxPillarCore(const size_t idx1, const size_t idx2,
                      const mesa_pd::data::Box& bx,
                      const mesa_pd::data::CylindricalBoundary& cb,
                      Accessor& ac)
   {
      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Vec3& c         = ac.getPosition(idx1);
      const mesa_pd::Vec3& cylCenter = ac.getPosition(idx2);
      const mesa_pd::Mat3  R         = ac.getRotation(idx1).getMatrix(); // body → world
      const mesa_pd::Mat3  Rt        = R.getTranspose();                 // world → body

      const mesa_pd::Vec3 half = bx.getEdgeLength() * real_t(0.5);

      if (cb.isContainment())
      {
         const real_t maxOverlap = real_t(1.0); // same envelope-safety scale as the pillar overload

         // Radial pre-gate: several tubes can share the same axial span (e.g. the 3 arches of a bridge), so a
         // body near another tube must not get a spurious huge overlap against this one. Uses the
         // interaction radius as a conservative bound on the box's reach from its centre.
         {
            const real_t        effR      = ac.getInteractionRadius(idx1);
            const auto&         axisCheck = cb.getAxis();
            const mesa_pd::Vec3 relCenter = c - cylCenter;
            const real_t distCenter = (relCenter - (relCenter * axisCheck) * axisCheck).length();
            if (distCenter - effR > cb.getRadius()) return false;
         }

         // The farthest point from the axis is a corner (radial distance to a line is convex, so its max over
         // a convex polytope is at a vertex): check all 8 and keep the largest that satisfies the half-clip,
         // if any (e.g. an arch tunnel is solid only above its springing line).
         const auto&   axis = cb.getAxis();
         const mesa_pd::Vec3& clipN = cb.getHalfClipNormal();
         real_t        farthestRadialSq = real_t(-1);
         mesa_pd::Vec3 worstCorner;
         // Runner-up corner (top-2 tracking, single pass), used when a second corner also breaches the tube
         // radius (an elongated box wedged with two corners out); see hasSecondaryContact().
         real_t        secondRadialSq = real_t(-1);
         mesa_pd::Vec3 secondCorner;
         real_t        axialLo = std::numeric_limits<real_t>::infinity();
         real_t        axialHi = -std::numeric_limits<real_t>::infinity();
         for (real_t sx : { real_t(-1), real_t(1) })
            for (real_t sy : { real_t(-1), real_t(1) })
               for (real_t sz : { real_t(-1), real_t(1) })
               {
                  const mesa_pd::Vec3 corner = c + R * mesa_pd::Vec3(sx * half[0], sy * half[1], sz * half[2]);
                  const real_t axialC = containment_detail::axialCoord(corner, cylCenter, axis);
                  axialLo = std::min(axialLo, axialC);
                  axialHi = std::max(axialHi, axialC);
                  if (clipN.sqrLength() > real_t(1e-20) && (corner - cylCenter) * clipN < real_t(0))
                     continue; // on the non-solid side (e.g. below an arch's springing line)
                  const mesa_pd::Vec3 rel     = corner - cylCenter;
                  const mesa_pd::Vec3 radial  = rel - (rel * axis) * axis;
                  const real_t        radSq   = radial.sqrLength();
                  if (radSq > farthestRadialSq)
                  {
                     secondRadialSq = farthestRadialSq; secondCorner = worstCorner;
                     farthestRadialSq = radSq; worstCorner = corner;
                  }
                  else if (radSq > secondRadialSq)
                  {
                     secondRadialSq = radSq; secondCorner = corner;
                  }
               }
         if (axialHi < cb.getAxialMin() || axialLo > cb.getAxialMax()) return false; // outside the tube's span
         if (farthestRadialSq < real_t(0)) return false; // every corner on the non-solid side

         real_t farthestRadial = std::sqrt(farthestRadialSq);
         real_t overlap = farthestRadial - cb.getRadius();
         if (overlap <= real_t(0)) return false; // safely inside

         // SAFETY: envelope clamp, caps the force a single degenerate/runaway overlap can inject.
         const bool envelopClamped = overlap > maxOverlap;
         if (envelopClamped)
         {
            static uint_t archEnvelopLogCount = uint_t(0);
            if (archEnvelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[ARCH-ENVELOP] box uid=" << ac.getUid(idx1)
                  << " deeply overlaps/envelops arch tunnel uid=" << ac.getUid(idx2)
                  << " (raw overlap=" << overlap << ", farthest-point radial dist=" << farthestRadial
                  << ") — clamping overlap to " << maxOverlap << ".");
               ++archEnvelopLogCount;
            }
            overlap = maxOverlap;
            farthestRadial = cb.getRadius() + maxOverlap;
         }

         const mesa_pd::Vec3 rel     = worstCorner - cylCenter;
         const mesa_pd::Vec3 outward = rel - (rel * axis) * axis;
         const real_t        outDist = outward.length();
         const mesa_pd::Vec3 outHat  = (outDist > real_t(1e-10))
                                          ? outward / outDist
                                          : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));

         getPenetrationDepth() = -overlap;
         getContactNormal()    = -outHat; // points inward, toward the axis
         getContactPoint()     = worstCorner;

         // Secondary contact (see hasSecondaryContact()): only if the runner-up corner itself pokes past the
         // radius, so an ordinary single-corner contact never reports a second one.
         if (secondRadialSq > cb.getRadius() * cb.getRadius())
         {
            real_t secondOverlap = std::sqrt(secondRadialSq) - cb.getRadius();
            if (secondOverlap > maxOverlap) secondOverlap = maxOverlap; // same envelope clamp as the primary
            const mesa_pd::Vec3 relS     = secondCorner - cylCenter;
            const mesa_pd::Vec3 outwardS = relS - (relS * axis) * axis;
            const real_t        outDistS = outwardS.length();
            const mesa_pd::Vec3 outHatS  = (outDistS > real_t(1e-10))
                                              ? outwardS / outDistS
                                              : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
            hasSecondaryContact_       = true;
            secondaryContactPoint_     = secondCorner;
            secondaryContactNormal_    = -outHatS;
            secondaryPenetrationDepth_ = -secondOverlap;
         }

         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "arch/tunnel containment (Box, ENVELOP-clamped)" : "arch/tunnel containment (Box)");
         return true;
      }

      // Cylinder axis as a LINE in the box body frame (cylinder is infinite along its axis).
      mesa_pd::Vec3 d  = Rt * cb.getAxis();
      const real_t  dl = d.length();
      d = (dl > real_t(1e-10)) ? d / dl : mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
      const mesa_pd::Vec3 a0 = Rt * (cylCenter - c); // a point on the axis (the cylinder centre), body frame

      // Closest point on the oriented box to the infinite axis line: project the current estimate onto the
      // line, clamp that line point into the box, refine once. (Clamping the cylinder centre directly would
      // bias toward the wrong corner for a rotated/elongated box.)
      mesa_pd::Vec3 boxPt(real_t(0), real_t(0), real_t(0));
      for (int it = 0; it < 2; ++it)
      {
         const real_t        t  = (boxPt - a0) * d;
         const mesa_pd::Vec3 lp = a0 + t * d; // nearest line point to the current box estimate
         boxPt = mesa_pd::Vec3(std::clamp(lp[0], -half[0], half[0]),
                               std::clamp(lp[1], -half[1], half[1]),
                               std::clamp(lp[2], -half[2], half[2]));
      }

      // Perpendicular distance from boxPt to the axis line.
      const mesa_pd::Vec3 vec    = boxPt - a0;
      const mesa_pd::Vec3 radial = vec - (vec * d) * d;
      const real_t        dist   = radial.length();

      real_t penetrationDepth = dist - cb.getRadius();
      if (penetrationDepth < getContactThreshold())
      {
         // SAFETY: a body enveloping the thin pillar (axis inside the box) yields a near-radius penetration
         // with a degenerate normal, which the spring-dashpot turns into a force that blows up the coupled
         // fluid (global NaN). The overlap is clamped to 1.0 LU and a warning logged. 1.0 LU stays above the
         // deepest real contacts logged by [DEM-WATCH] (0.68-0.77 LU) and the 0.5 LU overlap verified in
         // testCylinderPillarContact, so ordinary deep contacts are never clamped.
         const real_t maxOverlap = real_t(1.0);
         bool envelopClamped = false;
         if (penetrationDepth < -maxOverlap)
         {
            envelopClamped = true;
            static uint_t envelopLogCount = uint_t(0);
            if (envelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[BOX-CYL-ENVELOP] box uid=" << ac.getUid(idx1)
                  << " deeply overlaps/envelops cylinder uid=" << ac.getUid(idx2)
                  << " (raw penetration=" << penetrationDepth << ", dist-to-axis=" << dist
                  << ") — clamping overlap to " << -maxOverlap
                  << "; the box is likely too large for this pillar geometry.");
               ++envelopLogCount;
            }
            penetrationDepth = -maxOverlap;
         }
         getPenetrationDepth() = penetrationDepth;
         const mesa_pd::Vec3 nBody = (dist > real_t(1e-6))
                                        ? radial / dist
                                        : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         getContactNormal() = R * nBody;     // body → world
         getContactPoint()  = c + R * boxPt; // body → world
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "pillar exterior (Box, ENVELOP-clamped)" : "pillar exterior (Box)");
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   const mesa_pd::data::Box& bx,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, bx, cb, ac);
   }

   // Sphere vs finite Cylinder (floating log), analytic: the log lies among many bed spheres and GJK/EPA per
   // pair per sub-cycle is slower and noisier than the closed form. Works in the cylinder body frame
   // (axis = body z). Convention identical to sphere-sphere:
   //   penetrationDepth < 0 -> penetrating; contactNormal points from the CYLINDER (idx2) to the
   //   SPHERE (idx1), world frame; contactPoint = mid-overlap point, halfway between the two surfaces.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Sphere& s,
                   const mesa_pd::data::Cylinder& cyl,
                   Accessor& ac)
   {
      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Mat3 R  = ac.getRotation(idx2).getMatrix(); // cylinder body → world
      const mesa_pd::Mat3 Rt = R.getTranspose();
      const mesa_pd::Vec3 p  = Rt * (ac.getPosition(idx1) - ac.getPosition(idx2)); // sphere centre, body frame

      const real_t rc = cyl.getRadius();
      const real_t hh = cyl.getHeight() * real_t(0.5);
      const real_t rho = std::sqrt(p[0] * p[0] + p[1] * p[1]);

      mesa_pd::Vec3 nBody;      // outward normal at the closest surface point (body frame)
      real_t        separation; // signed surface distance minus sphere radius

      if (rho <= rc && std::abs(p[2]) <= hh)
      {
         // Sphere centre INSIDE the cylinder (deep overlap): push out through the nearest face so
         // the normal never degenerates. separation is negative by construction.
         const real_t dSide = rc - rho;
         const real_t dCap  = hh - std::abs(p[2]);
         if (dSide < dCap)
         {
            nBody = (rho > real_t(1e-10)) ? mesa_pd::Vec3(p[0] / rho, p[1] / rho, real_t(0))
                                          : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
            separation = -dSide - s.getRadius();
         }
         else
         {
            nBody      = mesa_pd::Vec3(real_t(0), real_t(0), math::sign(p[2]));
            separation = -dCap - s.getRadius();
         }
      }
      else
      {
         // Closest point on the (solid) cylinder surface via radial + axial clamp.
         const real_t scale = (rho > real_t(1e-10)) ? std::min(real_t(1), rc / rho) : real_t(0);
         const mesa_pd::Vec3 q(p[0] * scale, p[1] * scale, std::clamp(p[2], -hh, hh));
         const mesa_pd::Vec3 delta = p - q;
         const real_t dist = delta.length();
         nBody      = (dist > real_t(1e-10)) ? delta / dist : mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
         separation = dist - s.getRadius();
      }

      if (separation < getContactThreshold())
      {
         getPenetrationDepth() = separation;
         getContactNormal()    = R * nBody; // cylinder → sphere, world frame
         // Mid-overlap point: from the sphere centre, half the overlap past the sphere surface.
         getContactPoint() = ac.getPosition(idx1) - (s.getRadius() + real_t(0.5) * separation) * getContactNormal();
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Cylinder& cyl,
                   const mesa_pd::data::Sphere& s,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, s, cyl, ac);
   }

   // Finite Cylinder (floating log) vs CylindricalBoundary (solid pillar, infinite along its axis),
   // rotation-aware, analogous to Box x CylindricalBoundary: the pillar axis is transformed into the log's
   // body frame, the closest log point to that line is found by project-and-clamp iteration (radial + axial
   // clamp), and the radial distance to the axis gives the penetration.
   //   penetrationDepth = radial_dist - pillar_radius (< 0 -> penetrating)
   //   contactNormal points radially outward from the pillar axis toward the log (world frame)
   // cylinderContactSkin (ModelParameters::cylinderContactSkin): detect with the inflated cylinder
   // (radius + m, half-length + m); containment mode (arch/tunnel) is not inflated.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Cylinder& cyl,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   Accessor& ac)
   {
      const real_t m = cylinderContactSkin;
      if (m <= real_t(0) || cb.isContainment()) return cylinderPillarCore(idx1, idx2, cyl, cb, ac);
      const mesa_pd::data::Cylinder inflated(cyl.getRadius() + m, cyl.getHeight() + real_t(2) * m);
      return cylinderPillarCore(idx1, idx2, inflated, cb, ac);
   }

   template< typename Accessor >
   bool cylinderPillarCore(const size_t idx1, const size_t idx2,
                           const mesa_pd::data::Cylinder& cyl,
                           const mesa_pd::data::CylindricalBoundary& cb,
                           Accessor& ac)
   {
      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Vec3& c         = ac.getPosition(idx1);
      const mesa_pd::Vec3& cylCenter = ac.getPosition(idx2);
      const mesa_pd::Mat3  R         = ac.getRotation(idx1).getMatrix(); // log body → world
      const mesa_pd::Mat3  Rt        = R.getTranspose();                 // world → log body

      const real_t rc = cyl.getRadius();
      const real_t hh = cyl.getHeight() * real_t(0.5);

      if (cb.isContainment())
      {
         const real_t maxOverlap = real_t(1.0); // same envelope-safety scale as the pillar overload

         // Axial gate with the circumscribing-sphere bound (interactionRadius >= half-length): conservative,
         // never a false negative.
         const auto&   axis  = cb.getAxis();
         const real_t  effR  = ac.getInteractionRadius(idx1);
         const real_t  axialC = containment_detail::axialCoord(c, cylCenter, axis);
         if (axialC + effR < cb.getAxialMin() || axialC - effR > cb.getAxialMax()) return false;

         // Radial pre-gate: several tubes can share the same axial span; skip tubes whose axis is out of
         // reach (see the Box overload).
         {
            const mesa_pd::Vec3 relCenter  = c - cylCenter;
            const real_t        distCenter = (relCenter - (relCenter * axis) * axis).length();
            if (distCenter - effR > cb.getRadius()) return false;
         }

         // A solid cylinder's farthest point from a line lies on one of its two end-cap rims (squared radial
         // distance is convex in the axial parameter): search each rim's farthest point (with the half-clip,
         // if any) and keep the farther one.
         const mesa_pd::Vec3 logAxis = R * mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
         mesa_pd::Vec3 uB, vB;
         closed_form_contacts::detail::perpendicularBasis(logAxis, uB, vB);

         real_t        farthestRadialSq = real_t(-1);
         mesa_pd::Vec3 worstPoint;
         // Runner-up = the other rim, if it also breaches (a log lying across the tube with both end caps
         // out); see hasSecondaryContact().
         real_t        secondRadialSq = real_t(-1);
         mesa_pd::Vec3 secondPoint;
         for (const real_t sgn : { real_t(-1), real_t(1) })
         {
            const mesa_pd::Vec3 rimC = c + sgn * hh * logAxis;
            const auto far = containment_detail::circleFarthestFromAxisClipped(
               rimC, uB, vB, rc, cylCenter, axis, cylCenter, cb.getHalfClipNormal());
            if (!far.valid) continue;
            if (far.radialDistSq > farthestRadialSq)
            {
               secondRadialSq = farthestRadialSq; secondPoint = worstPoint;
               farthestRadialSq = far.radialDistSq; worstPoint = far.point;
            }
            else if (far.radialDistSq > secondRadialSq)
            {
               secondRadialSq = far.radialDistSq; secondPoint = far.point;
            }
         }
         if (farthestRadialSq < real_t(0)) return false; // both rims entirely on the non-solid side

         real_t farthestRadial = std::sqrt(farthestRadialSq);
         real_t overlap = farthestRadial - cb.getRadius();
         if (overlap <= real_t(0)) return false; // safely inside

         // SAFETY: envelope clamp, as in the pillar-exterior overloads.
         const bool envelopClamped = overlap > maxOverlap;
         if (envelopClamped)
         {
            static uint_t archEnvelopLogCount = uint_t(0);
            if (archEnvelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[ARCH-ENVELOP] cylinder uid=" << ac.getUid(idx1)
                  << " deeply overlaps/envelops arch tunnel uid=" << ac.getUid(idx2)
                  << " (raw overlap=" << overlap << ", farthest-point radial dist=" << farthestRadial
                  << ") — clamping overlap to " << maxOverlap << ".");
               ++archEnvelopLogCount;
            }
            overlap = maxOverlap;
            farthestRadial = cb.getRadius() + maxOverlap;
         }

         const mesa_pd::Vec3 rel     = worstPoint - cylCenter;
         const mesa_pd::Vec3 outward = rel - (rel * axis) * axis;
         const real_t        outDist = outward.length();
         const mesa_pd::Vec3 outHat  = (outDist > real_t(1e-10))
                                          ? outward / outDist
                                          : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));

         getPenetrationDepth() = -overlap;
         getContactNormal()    = -outHat; // points inward, toward the axis
         getContactPoint()     = worstPoint;

         // Secondary contact: the other rim, only if it breaches on its own (see hasSecondaryContact()).
         if (secondRadialSq > cb.getRadius() * cb.getRadius())
         {
            real_t secondOverlap = std::sqrt(secondRadialSq) - cb.getRadius();
            if (secondOverlap > maxOverlap) secondOverlap = maxOverlap; // same envelope clamp as the primary
            const mesa_pd::Vec3 relS     = secondPoint - cylCenter;
            const mesa_pd::Vec3 outwardS = relS - (relS * axis) * axis;
            const real_t        outDistS = outwardS.length();
            const mesa_pd::Vec3 outHatS  = (outDistS > real_t(1e-10))
                                              ? outwardS / outDistS
                                              : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
            hasSecondaryContact_       = true;
            secondaryContactPoint_     = secondPoint;
            secondaryContactNormal_    = -outHatS;
            secondaryPenetrationDepth_ = -secondOverlap;
         }

         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "arch/tunnel containment (Cylinder, ENVELOP-clamped)" : "arch/tunnel containment (Cylinder)");
         return true;
      }

      // refineCylinderContacts: the solid pillar is represented locally by a finite cylinder of the pillar radius,
      // centred on the projection of the rod centre onto the pillar axis and long enough that its end faces never
      // matter, and resolved by the refined cylinder-cylinder detection. The default path below (3 alternating
      // projections onto the axis line) misses some contacts and underestimates depths against the independent
      // reference (thesis Section 5.2; scripts/contact-validation/obstacle_reference.py in the thesis repository).
      if (refineCylinderContacts)
      {
         const mesa_pd::Vec3 ap      = cb.getAxis().getNormalizedIfNotZero();
         const mesa_pd::Vec3 foot    = cylCenter + ((c - cylCenter) * ap) * ap;
         // half-length: any overlap along the pillar axis is then >= halfLen > the largest possible radial depth
         // (R_p + R), so the finite end faces never give the minimum
         const real_t        halfLen = real_t(2) * (rc + hh + cb.getRadius());
         mesa_pd::Vec3 e1 = ap % mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         if (e1.length() < real_t(0.5)) e1 = ap % mesa_pd::Vec3(real_t(0), real_t(1), real_t(0));
         e1 = e1.getNormalized();
         const mesa_pd::Vec3 e2 = ap % e1;
         const mesa_pd::Mat3 Rp(e1[0], e2[0], ap[0], e1[1], e2[1], ap[1], e1[2], e2[2], ap[2]); // columns e1, e2, axis
         const closed_form_contacts::CylinderState rodS{ c, R, rc, real_t(2) * hh };
         const closed_form_contacts::CylinderState pilS{ foot, Rp, cb.getRadius(), real_t(2) * halfLen };
         const auto res = closed_form_contacts::cylinderCylinderContactRefined(rodS, pilS);
         if (!res.contact) return false;
         real_t penetrationDepth = -res.overlap;
         const real_t maxOverlap = real_t(1.0); // same envelope safety clamp as the default path below
         const bool envelopClamped = penetrationDepth < -maxOverlap;
         if (envelopClamped) penetrationDepth = -maxOverlap;
         getPenetrationDepth() = penetrationDepth;
         getContactNormal()    = res.normal; // points towards the first argument, the rod: pillar -> log, as below
         getContactPoint()     = res.point;
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "pillar exterior (Cylinder, refined, ENVELOP-clamped)" : "pillar exterior (Cylinder, refined)");
         return true;
      }

      // Pillar axis as a LINE in the log body frame.
      mesa_pd::Vec3 d  = Rt * cb.getAxis();
      const real_t  dl = d.length();
      d = (dl > real_t(1e-10)) ? d / dl : mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
      const mesa_pd::Vec3 a0 = Rt * (cylCenter - c);

      // Closest point of the log to the infinite pillar axis: project the estimate onto the line,
      // clamp the nearest line point into the log (radial in xy, axial in z), refine.
      mesa_pd::Vec3 logPt(real_t(0), real_t(0), real_t(0));
      for (int it = 0; it < 3; ++it)
      {
         const real_t        t  = (logPt - a0) * d;
         const mesa_pd::Vec3 lp = a0 + t * d;
         const real_t rho   = std::sqrt(lp[0] * lp[0] + lp[1] * lp[1]);
         const real_t scale = (rho > real_t(1e-10)) ? std::min(real_t(1), rc / rho) : real_t(0);
         logPt = mesa_pd::Vec3(lp[0] * scale, lp[1] * scale, std::clamp(lp[2], -hh, hh));
      }

      const mesa_pd::Vec3 vec    = logPt - a0;
      const mesa_pd::Vec3 radial = vec - (vec * d) * d;
      const real_t        dist   = radial.length();

      real_t penetrationDepth = dist - cb.getRadius();
      if (penetrationDepth < getContactThreshold())
      {
         // Envelope safety clamp as in Box x CylindricalBoundary: a log enveloping the thin pillar yields a
         // near-radius overlap + degenerate normal -> unbounded spring force -> fluid NaN.
         const real_t maxOverlap = real_t(1.0);
         bool envelopClamped = false;
         if (penetrationDepth < -maxOverlap)
         {
            envelopClamped = true;
            static uint_t envelopLogCount = uint_t(0);
            if (envelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[CYL-CYL-ENVELOP] floating cylinder uid=" << ac.getUid(idx1)
                  << " deeply overlaps/envelops pillar uid=" << ac.getUid(idx2)
                  << " (raw penetration=" << penetrationDepth << ", dist-to-axis=" << dist
                  << ") — clamping overlap to " << -maxOverlap
                  << "; the log is likely too large for this pillar geometry.");
               ++envelopLogCount;
            }
            penetrationDepth = -maxOverlap;
         }
         getPenetrationDepth() = penetrationDepth;
         const mesa_pd::Vec3 nBody = (dist > real_t(1e-6))
                                        ? radial / dist
                                        : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         getContactNormal() = R * nBody;     // body → world, pillar → log
         getContactPoint()  = c + R * logPt; // body → world
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "pillar exterior (Cylinder, ENVELOP-clamped)" : "pillar exterior (Cylinder)");
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   const mesa_pd::data::Cylinder& cyl,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, cyl, cb, ac);
   }

   // Capsule (pill-shaped rod) vs CylindricalBoundary (solid pillar). Needed because the generic
   // GeneralContactDetection::operator()(CylindricalBoundary, GEO2_T) fallback always uses CONTAINMENT
   // semantics (collision when the body's FARTHEST point from the axis exceeds the radius) and ignores
   // CylindricalBoundary::isContainment(); a solid pillar needs collision when the CLOSEST point gets nearer
   // than the radius (checked by testCapsuleCylindricalBoundaryContact). Containment mode is left to the
   // generic path, which implements it correctly.
   //
   // Closed form, not iterative: a capsule's closest surface point to a line is its core segment's closest
   // point moved `radius` toward the line, so only a 1-D clamp of the axial parameter to
   // [-coreHalfLength, +coreHalfLength] is needed. From the two-lines closest-point normal equations (capsule
   // body-z axis vs. infinite pillar axis), then clamped; exact because the squared distance is a convex
   // quadratic in that parameter.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Capsule& capsule,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   Accessor& ac)
   {
      if (cb.isContainment())
         return GeneralContactDetection::operator()(idx1, idx2, capsule, cb, ac);

      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Vec3& c         = ac.getPosition(idx1);
      const mesa_pd::Vec3& cylCenter = ac.getPosition(idx2);
      const mesa_pd::Mat3  R         = ac.getRotation(idx1).getMatrix(); // capsule body → world
      const mesa_pd::Mat3  Rt        = R.getTranspose();                 // world → capsule body

      const real_t rCap = capsule.getRadius();
      const real_t a    = capsule.getCoreHalfLength();

      // Pillar axis as a line in the capsule's body frame.
      mesa_pd::Vec3 d  = Rt * cb.getAxis();
      const real_t  dl = d.length();
      d = (dl > real_t(1e-10)) ? d / dl : mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
      const mesa_pd::Vec3 a0 = Rt * (cylCenter - c);

      // Closed-form closest point of the (infinite) body-z axis to the pillar axis line: from
      // d/dz|P|^2=0 and d/dt|P|^2=0 for P=(0,0,z)-a0-t*d, eliminating t gives
      // z*(1-d_z^2) = a0_z - (a0.d)*d_z. Falls back to z=0 for (near-)parallel axes (no unique closest point).
      const real_t denom  = real_t(1) - d[2] * d[2];
      const real_t zStar  = (std::abs(denom) > real_t(1e-10)) ? (a0[2] - (a0 * d) * d[2]) / denom : real_t(0);
      const real_t zClamp = std::clamp(zStar, -a, a);

      const mesa_pd::Vec3 segPt(real_t(0), real_t(0), zClamp);
      const mesa_pd::Vec3 vec    = segPt - a0;
      const mesa_pd::Vec3 radial = vec - (vec * d) * d;
      // core-segment-to-axis distance; subtract rCap below for the surface
      const real_t        dist   = radial.length();

      real_t penetrationDepth = (dist - rCap) - cb.getRadius();
      if (penetrationDepth < getContactThreshold())
      {
         // Envelope safety clamp as in the Cylinder overload: a capsule enveloping the pillar yields a
         // near-radius overlap + degenerate normal -> unbounded spring force.
         const real_t maxOverlap = real_t(1.0);
         bool envelopClamped = false;
         if (penetrationDepth < -maxOverlap)
         {
            envelopClamped = true;
            static uint_t envelopLogCount = uint_t(0);
            if (envelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[CAPSULE-CYL-ENVELOP] floating capsule uid=" << ac.getUid(idx1)
                  << " deeply overlaps/envelops pillar uid=" << ac.getUid(idx2)
                  << " (raw penetration=" << penetrationDepth << ", core-segment-to-axis dist=" << dist
                  << ") — clamping overlap to " << -maxOverlap
                  << "; the capsule is likely too large for this pillar geometry.");
               ++envelopLogCount;
            }
            penetrationDepth = -maxOverlap;
         }
         getPenetrationDepth() = penetrationDepth;
         const mesa_pd::Vec3 nBody = (dist > real_t(1e-6))
                                        ? radial / dist
                                        : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         getContactNormal() = R * nBody;                        // body → world, pillar → capsule
         getContactPoint()  = c + R * (segPt - rCap * nBody);   // body → world, actual capsule surface point
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            envelopClamped ? "pillar exterior (Capsule, ENVELOP-clamped)" : "pillar exterior (Capsule)");
         return true;
      }
      return false;
   }

   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::CylindricalBoundary& cb,
                   const mesa_pd::data::Capsule& capsule,
                   Accessor& ac)
   {
      return operator()(idx2, idx1, capsule, cb, ac);
   }

   // Box vs Box / Cylinder vs Cylinder: optional closed-form alternative (ClosedFormContacts*.h) to the
   // GJK/EPA overloads inherited from GeneralContactDetection, switched at runtime by
   // useAnalyticalContactDetection (ModelParameters, SimulationParams.h; default false = GJK/EPA only).
   // Cross-checked against GJK/EPA in PaperAnalyticalContactTest.cpp.
   //
   // Box x Box: boxBoxSAT has no known discrepancy vs. GJK/EPA across the test suite, so no fallback gate.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Box& bx1,
                   const mesa_pd::data::Box& bx2,
                   Accessor& ac)
   {
      // boxContactSkin (see ModelParameters): both boxes are replaced by inflated boxes (each half-extent + m)
      // and resolved by the dispatch below; the skin is cleared during the recursive call so the inflation is
      // applied exactly once (as for cylinderContactSkin in Cylinder x Cylinder).
      if (boxContactSkin > real_t(0))
      {
         const real_t m = boxContactSkin;
         const mesa_pd::data::Box inflated1(bx1.getEdgeLength() + mesa_pd::Vec3(real_t(2) * m));
         const mesa_pd::data::Box inflated2(bx2.getEdgeLength() + mesa_pd::Vec3(real_t(2) * m));
         boxContactSkin = real_t(0);
         const bool hit = operator()(idx1, idx2, inflated1, inflated2, ac);
         boxContactSkin = m;
         return hit;
      }

      if (!useAnalyticalContactDetection)
      {
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2),
            "GJK/EPA (Box-Box, useAnalyticalContactDetection=false)");
         return GeneralContactDetection::operator()(idx1, idx2, bx1, bx2, ac);
      }

      // useTightLayer1: cheap pre-filter ahead of the full 15-axis boxBoxSAT test that tests one legitimate
      // separating axis (the line between the centres), so it can reject but never confirm contact (see
      // closed_form_contacts::detail::boxesDefinitelySeparatedAlongCenterLine). Saves a SAT call for
      // obviously separated pairs (a smaller win than the Cylinder x Cylinder gate below, which saves GJK/EPA).
      if (useTightLayer1 &&
          closed_form_contacts::detail::boxesDefinitelySeparatedAlongCenterLine(
             ac.getPosition(idx1), ac.getRotation(idx1).getMatrix(), bx1.getEdgeLength() * real_t(0.5),
             ac.getPosition(idx2), ac.getRotation(idx2).getMatrix(), bx2.getEdgeLength() * real_t(0.5)))
      {
         ++contactStats.boxTightLayer1CheaplyProvedSeparated;
         return false;
      }

      getIdx1() = idx1;
      getIdx2() = idx2;

      const mesa_pd::Vec3 halfA = bx1.getEdgeLength() * real_t(0.5), halfB = bx2.getEdgeLength() * real_t(0.5);
      const mesa_pd::Vec3 posA = ac.getPosition(idx1), posB = ac.getPosition(idx2);
      const mesa_pd::Mat3 RA = ac.getRotation(idx1).getMatrix(), RB = ac.getRotation(idx2).getMatrix();

      const auto result = closed_form_contacts::boxBoxSAT(posA, RA, halfA, posB, RB, halfB);
      if (!result.contact) return false;

      ++contactStats.boxBoxSAT;
      recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), "closed-form analytical (Box-Box SAT)");
      // closed_form_contacts convention: overlap > 0 = penetrating. This class: penetrationDepth < 0.
      getPenetrationDepth() = -result.overlap;
      getContactNormal()    = result.normal;
      getContactPoint()     = result.point;
      return true;
   }

   // Cylinder x Cylinder: cylinderCylinderContact's rim-vs-face/band branches (face-edge, special-contact,
   // band-edge) can report the wrong overlap magnitude for some oblique tilts (the true closest feature lies on
   // a curved band neither Kodam et al. nor Guo et al.'s corrections search). cylinderCylinderResultIsUncertain
   // flags only special-contact and pairs the closed forms found no contact for within the centre-distance bound
   // (Tier 1); those fall back to GJK/EPA, every other scenario is trusted as computed.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Cylinder& cyl1,
                   const mesa_pd::data::Cylinder& cyl2,
                   Accessor& ac)
   {
      // Canonicalise by uid before cylinderCylinderContact, as GeneralContactDetection does for GJK/EPA:
      // the rim-vs-face/band branch is not symmetric in its arguments (rimVsFaceOrBand(ci,cj) is preferred
      // over rimVsFaceOrBand(cj,ci)), so the argument order must not depend on the rank-local storage index,
      // or two ranks holding the same pair could compute two different valid contacts. uid is global, and
      // positions/orientations are already synced before contact detection.
      if (ac.getUid(idx2) < ac.getUid(idx1))
         return operator()(idx2, idx1, cyl2, cyl1, ac);

      // cylinderContactSkin (see ModelParameters): both rods are replaced by inflated cylinders (radius + m,
      // half-length + m) and resolved by the dispatch below; the skin is cleared during the recursive call so
      // the inflation is applied exactly once.
      if (cylinderContactSkin > real_t(0))
      {
         const real_t m = cylinderContactSkin;
         const mesa_pd::data::Cylinder inflated1(cyl1.getRadius() + m, cyl1.getHeight() + real_t(2) * m);
         const mesa_pd::data::Cylinder inflated2(cyl2.getRadius() + m, cyl2.getHeight() + real_t(2) * m);
         cylinderContactSkin = real_t(0);
         const bool hit = operator()(idx1, idx2, inflated1, inflated2, ac);
         cylinderContactSkin = m;
         return hit;
      }

      // Cheap reject before the closed forms or GJK/EPA run. Cylinder-Cylinder has no built-in early bound like
      // detectSphereSphereCollision's |pos1-pos2| < r1+r2, and a bad (e.g. non-minimum-image ghost) position
      // could otherwise yield a false-positive contact for a pair far beyond its combined reach.
      // useTightLayer1=false: circumscribing-sphere bound, getInteractionRadius() = sqrt(radius^2 + halfLength^2)
      // + 2*cylinderContactSkin (set at floater creation in SetupStages.h); skipped if it is unset (0, e.g. unit
      // tests building particles directly), so pairs are not rejected wholesale.
      // useTightLayer1=true: closed_form_contacts::detail::closestDistanceSqBetweenSegments of the two axis
      // segments vs. the sum of radii; equally conservative but orientation-aware, so it rejects far more
      // separated pairs.
      const closed_form_contacts::CylinderState ci{ ac.getPosition(idx1), ac.getRotation(idx1).getMatrix(),
                                                 cyl1.getRadius(), cyl1.getHeight() };
      const closed_form_contacts::CylinderState cj{ ac.getPosition(idx2), ac.getRotation(idx2).getMatrix(),
                                                 cyl2.getRadius(), cyl2.getHeight() };
      if (useTightLayer1)
      {
         mesa_pd::Vec3 p0i, p1i, p0j, p1j;
         closed_form_contacts::detail::cylinderAxisEndpoints(ci, p0i, p1i);
         closed_form_contacts::detail::cylinderAxisEndpoints(cj, p0j, p1j);
         const real_t reach = ci.radius + cj.radius;
         if (closed_form_contacts::detail::closestDistanceSqBetweenSegments(p0i, p1i, p0j, p1j) > reach * reach)
            return false;
      }
      else
      {
         const real_t maxPossibleReach = ac.getInteractionRadius(idx1) + ac.getInteractionRadius(idx2);
         if (maxPossibleReach > real_t(0) &&
             (ac.getPosition(idx1) - ac.getPosition(idx2)).sqrLength() > maxPossibleReach * maxPossibleReach)
            return false;
      }

      if (!useAnalyticalContactDetection)
      {
         const bool gjkHit = GeneralContactDetection::operator()(idx1, idx2, cyl1, cyl2, ac);
         if (shadowAnalyticLog) shadowCompare(idx1, idx2, ci, cj, gjkHit, ac);
         return gjkHit;
      }

      auto result = closed_form_contacts::cylinderCylinderContact(ci, cj);

      // Capsule-bound safety net: a cylinder lies inside the capsule of the same radius around its axis segment,
      // so two cylinders can only touch if their axis segments are within R_i + R_j. An analytic "contact"
      // violating this is a closed-form false positive (cf. CylinderNearParallelSkewedTest.cpp) and must not
      // reach the solver, in particular not with usePureAnalyticalCylinder (no GJK/EPA cross-check). Cheap, and
      // never rejects a real contact (tiny relative slack for round-off only).
      // Skipped when useTightLayer1=true: the pre-filter above already rejected dist(axis segments) > R_i + R_j
      // with the same CylinderState, so this check could not change the outcome.
      if (result.contact && !useTightLayer1)
      {
         mesa_pd::Vec3 p0i, p1i, p0j, p1j;
         closed_form_contacts::detail::cylinderAxisEndpoints(ci, p0i, p1i);
         closed_form_contacts::detail::cylinderAxisEndpoints(cj, p0j, p1j);
         const real_t reach = (ci.radius + cj.radius) * (real_t(1) + real_t(1e-9));
         if (closed_form_contacts::detail::closestDistanceSqBetweenSegments(p0i, p1i, p0j, p1j) > reach * reach)
         {
            ++contactStats.cylCapsuleRejected;
            return false;
         }
      }

      // refineCylinderContacts (opt-in): the closed-form result is refined by minimising the projection overlap
      // over directions (closed_form_contacts::refineCylinderContact). Exact scenarios are unchanged, inexact ones
      // (edge-edge, band-edge, face-edge, special contact) get the minimising depth and normal, and "no scenario"
      // pairs are decided by the same descent (separation proven if a negative overlap is found). No GJK/EPA call;
      // GJK/EPA remains the verification reference in PaperAnalyticalContactTest.cpp.
      if (refineCylinderContacts)
      {
         static const bool shadowTiming = std::getenv("ANTIDUNES_REFINE_SHADOW_TIMING") != nullptr;
         if (shadowTiming)
         {
            // outputs written by the shadow GJK/EPA call are overwritten below before any contact is returned
            const auto s0 = std::chrono::steady_clock::now();
            const auto rawAgain = closed_form_contacts::cylinderCylinderContact(ci, cj);
            const auto s1 = std::chrono::steady_clock::now();
            if (closed_form_contacts::cylinderCylinderResultIsUncertain(ci, cj, rawAgain))
            {
               ++contactStats.shadowFallbackCalls;
               (void) GeneralContactDetection::operator()(idx1, idx2, cyl1, cyl2, ac);
            }
            const auto s2 = std::chrono::steady_clock::now();
            contactStats.shadowClosedFormSeconds += std::chrono::duration< double >(s1 - s0).count();
            contactStats.shadowFallbackSeconds += std::chrono::duration< double >(s2 - s1).count();
         }
         closed_form_contacts::RefineTrace trace;
         const std::string rawScenario = result.contact ? std::string(result.scenario) : std::string();
         const real_t      rawOverlap  = result.overlap;
         const auto        t0          = std::chrono::steady_clock::now();
         result = closed_form_contacts::refineCylinderContact(ci, cj, result, &trace);
         const double dt = std::chrono::duration< double >(std::chrono::steady_clock::now() - t0).count();
         ++contactStats.cylRefinePath[int(trace.path)];
         contactStats.cylRefinePathSeconds[int(trace.path)] += dt;
         ++contactStats.cylRefineFastReject[int(trace.fastReject)];
         contactStats.cylRefineKinks += trace.kinks ? 1 : 0;
         contactStats.cylRefineConvex += trace.convex ? 1 : 0;
         if (trace.path == closed_form_contacts::RefinePath::ContactDescended ||
             trace.path == closed_form_contacts::RefinePath::ContactDescendedSeparated)
         {
            // Profiling dump (ANTIDUNES_REFINE_DUMP=<prefix>): every 499th re-descended pair, at most 20000 per
            // rank, to <prefix>_rank<r>.csv, for offline analysis of the fast-path rejections.
            static const char*   dumpPrefix = std::getenv("ANTIDUNES_REFINE_DUMP");
            static std::ofstream dump;
            static uint_t        dumped = 0, seen = 0;
            if (dumpPrefix && dumped < uint_t(20000) && (seen++ % uint_t(499)) == uint_t(0))
            {
               if (!dump.is_open())
               {
                  dump.open(std::string(dumpPrefix) + "_rank" + std::to_string(mpi::MPIManager::instance()->rank()) + ".csv");
                  dump << "scenario,raw_overlap,ref_overlap,ref_contact,fast_reject,kinks,us,fast_d0,fast_q,"
                          "pix,piy,piz,aix,aiy,aiz,Ri,Li,pjx,pjy,pjz,ajx,ajy,ajz,Rj,Lj\n";
                  dump.precision(17);
               }
               const mesa_pd::Vec3 ai = ci.axis(), aj = cj.axis();
               dump << rawScenario << "," << rawOverlap << "," << result.overlap << "," << int(result.contact) << ","
                    << int(trace.fastReject) << "," << int(trace.kinks) << "," << 1e6 * dt << "," << trace.fastD0 << ","
                    << trace.fastQ << "," << ci.pos[0] << ","
                    << ci.pos[1] << "," << ci.pos[2] << "," << ai[0] << "," << ai[1] << "," << ai[2] << "," << ci.radius
                    << "," << ci.length << "," << cj.pos[0] << "," << cj.pos[1] << "," << cj.pos[2] << "," << aj[0] << ","
                    << aj[1] << "," << aj[2] << "," << cj.radius << "," << cj.length << "\n";
               ++dumped;
            }
            static const char* names[8] = { "face-face", "band-band-parallel", "face-band", "band-band-skewed",
                                            "face-edge", "special-contact", "band-edge", "edge-edge" };
            for (int k = 0; k < 8; ++k)
               if (rawScenario == names[k]) ++contactStats.cylRefineDescendedScenario[k];
         }
         if (!result.contact)
         {
            ++contactStats.cylRefinedSeparated;
            return false;
         }
         ++contactStats.cylRefinedContact;
         const std::string path = std::string("refined closed form (scenario=") + result.scenario + ")";
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), path);
         logCylComputation(ac.getUid(idx1), ac.getUid(idx2), path, ci, cj, result);
         getIdx1() = idx1;
         getIdx2() = idx2;
         getPenetrationDepth() = -result.overlap;
         getContactNormal()    = result.normal;
         getContactPoint()     = result.point;
         return true;
      }

      if (closed_form_contacts::cylinderCylinderResultIsUncertain(ci, cj, result))
      {
         // Tier 1 (closed forms found no contact) vs special-contact (a contact was found, but that scenario is
         // always distrusted); see cylinderCylinderResultIsUncertain in ClosedFormContactsCylinderCylinder.h.
         if (usePureAnalyticalCylinder)
         {
            // "Pure Kodam+Guo" mode: no GJK/EPA safety net; trust the closed-form answer even where
            // cylinderCylinderResultIsUncertain() flags it. See ModelParameters.usePureAnalyticalCylinder.
            if (!result.contact)
            {
               ++contactStats.cylPureTier1Rejected;
               const std::string path = "pure analytical (Kodam+Guo): tier1, no contact found, "
                  "no fallback attempted (usePureAnalyticalCylinder)";
               recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), path);
               logCylComputation(ac.getUid(idx1), ac.getUid(idx2), path, ci, cj, result);
               return false;
            }
            ++contactStats.cylPureSpecialContactTrusted;
            const std::string path = "pure analytical (Kodam+Guo): special-contact result trusted despite "
               "uncertainty flag (usePureAnalyticalCylinder)";
            recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), path);
            logCylComputation(ac.getUid(idx1), ac.getUid(idx2), path, ci, cj, result);
            contactStats.bumpCylScenario(result.scenario);
            getIdx1() = idx1;
            getIdx2() = idx2;
            getPenetrationDepth() = -result.overlap;
            getContactNormal()    = result.normal;
            getContactPoint()     = result.point;
            return true;
         }

         if (result.contact) ++contactStats.cylFallbackSpecialContact;
         else ++contactStats.cylFallbackTier1;
         const std::string path =
            result.contact ? "GJK/EPA fallback (special-contact: unconditionally distrusted scenario)"
                            : "GJK/EPA fallback (tier1: analytical found no contact at all)";
         recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), path);
         logCylComputation(ac.getUid(idx1), ac.getUid(idx2), path, ci, cj, result);
         const bool fallbackContact = GeneralContactDetection::operator()(idx1, idx2, cyl1, cyl2, ac);
         if (result.contact)
         {
            if (fallbackContact) ++contactStats.cylFallbackSpecialContactGjkFoundContact;
            else ++contactStats.cylFallbackSpecialContactGjkFoundNoContact;
         }
         else
         {
            if (fallbackContact) ++contactStats.cylFallbackTier1GjkFoundContact;
            else ++contactStats.cylFallbackTier1GjkFoundNoContact;
         }
         return fallbackContact;
      }

      if (!result.contact)
      {
         // Dispatch-gap pair (no scenario) that cylinderCylinderResultIsUncertain() did NOT flag: its
         // centre-distance bound proves the cylinders cannot touch (rare, since the Layer-1 reject above is
         // tighter), so "no contact" is trusted without GJK/EPA. Counted separately from cylFallbackTier1.
         ++contactStats.cylTier1CheaplyProvedSeparated;
         return false;
      }
      contactStats.bumpCylScenario(result.scenario);
      const std::string path = std::string("closed-form analytical (scenario=") + result.scenario + ")";
      recordDetectionPath(ac.getUid(idx1), ac.getUid(idx2), path);
      logCylComputation(ac.getUid(idx1), ac.getUid(idx2), path, ci, cj, result);
      getIdx1() = idx1;
      getIdx2() = idx2;
      getPenetrationDepth() = -result.overlap;
      getContactNormal()    = result.normal;
      getContactPoint()     = result.point;
      return true;
   }

   // Diagnostic side-channel (rank-local, not synced): per uid pair, the most recently recorded detection path
   // (Box-Box / Cylinder-Cylinder dispatch, pillar/containment overloads), so a later [ContactHistory DUPLICATE]
   // warning (ContactHistoryNotification.h, same reduction cycle) can report it via
   // contactHistoryDuplicateDiagnosticHook.
   void recordDetectionPath(walberla::id_t uidA, walberla::id_t uidB, std::string path)
   {
      const auto key = uidA < uidB ? std::make_pair(uidA, uidB) : std::make_pair(uidB, uidA);
      pairDetectionPath[key] = std::move(path);
   }

   const std::string* lookupDetectionPath(walberla::id_t uidA, walberla::id_t uidB) const
   {
      const auto key = uidA < uidB ? std::make_pair(uidA, uidB) : std::make_pair(uidB, uidA);
      const auto it = pairDetectionPath.find(key);
      return it == pairDetectionPath.end() ? nullptr : &it->second;
   }

   // Log of every Cylinder x Cylinder computation on this rank: rank, both uids, path/tier, both poses and the
   // resulting contact. Ghost-owning ranks may compute a pair's contact that the reducing ranks never see, so
   // grepping a pair's uids across all ranks shows whether every rank got the same scenario/point/overlap.
   void logCylComputation(walberla::id_t uidA, walberla::id_t uidB, const std::string& path,
                           const closed_form_contacts::CylinderState& ci, const closed_form_contacts::CylinderState& cj,
                           const closed_form_contacts::ContactResult& result)
   {
      static uint_t logCount = uint_t(0);
      if (logCount >= uint_t(20000)) return; // cap on log volume
      ++logCount;
      // WARNING, not INFO: MPI_ABORT kills all ranks without flushing buffered stdio, so INFO output from a
      // crashing reduction cycle can be lost; WARNING output survives.
      WALBERLA_LOG_WARNING(
         "[CylComputation] rank=" << walberla::mpi::MPIManager::instance()->rank()
         << " uidA=" << uidA << " uidB=" << uidB << " path=[" << path << "]"
         << " ci.pos=" << ci.pos << " cj.pos=" << cj.pos
         << " contact=" << result.contact << " overlap=" << result.overlap
         << " normal=" << result.normal << " point=" << result.point
         << " scenario=" << result.scenario);
   }

   // Diagnostics: per-rank tallies of which closed-form scenario / GJK/EPA fallback path fired, incremented only
   // while useAnalyticalContactDetection is on. Fixed named counters (not a string-keyed map) so
   // ParticleKernels::logContactDetectionStatistics() can MPI-sum them with scalar allReduceInplace calls.
   struct ContactDetectionStats
   {
      uint64_t boxBoxSAT           = 0;
      uint64_t cylFaceFace         = 0;
      uint64_t cylBandBandParallel = 0;
      uint64_t cylFaceBand         = 0;
      uint64_t cylBandBandSkewed   = 0;
      uint64_t cylFaceEdge         = 0;
      uint64_t cylSpecialContact   = 0;
      uint64_t cylBandEdge         = 0;
      uint64_t cylEdgeEdge         = 0;
      // cylinderCylinderResultIsUncertain() overruled the closed form and GJK/EPA was used instead.
      // cylFallbackTier1: closed forms found no scenario. cylFallbackSpecialContact: a contact was found, but
      // the "special-contact" scenario is always distrusted.
      uint64_t cylFallbackTier1          = 0;
      uint64_t cylFallbackSpecialContact = 0;
      // What GJK/EPA found on each fallback path: Tier 1 fallbacks that GJK/EPA also resolves to "no contact"
      // are the pairs a coarse bounding check should reject cheaply (cf. Kodam et al. 2010, Section 4). Purely a
      // diagnostic, not part of Kodam's contact-detection equations.
      uint64_t cylFallbackTier1GjkFoundContact            = 0;
      uint64_t cylFallbackTier1GjkFoundNoContact           = 0;
      uint64_t cylFallbackSpecialContactGjkFoundContact   = 0;
      uint64_t cylFallbackSpecialContactGjkFoundNoContact = 0;
      // Used instead of the counters above when usePureAnalyticalCylinder is true (no GJK/EPA fallback).
      uint64_t cylPureTier1Rejected          = 0; // analytical found no contact; accepted as "no contact" as-is
      uint64_t cylPureSpecialContactTrusted  = 0; // special-contact flagged uncertain; trusted as-is anyway
      // Dispatch-gap pairs (no scenario) NOT flagged uncertain: a bound proved them separated, so "no contact"
      // is trusted without GJK/EPA. Mutually exclusive with cylFallbackTier1.
      uint64_t cylTier1CheaplyProvedSeparated = 0;
      // Box x Box pairs the useTightLayer1 single-axis pre-filter proved separated before boxBoxSAT ran.
      uint64_t boxTightLayer1CheaplyProvedSeparated = 0;
      // Cylinder x Cylinder analytic "contacts" rejected by the capsule bound (axis segments farther apart than
      // R_i + R_j). Must stay 0 unless the closed-form dispatch produced a false positive.
      uint64_t cylCapsuleRejected = 0;
      // refineCylinderContacts: pairs the support-function refinement reported as touching / as separated
      // (no GJK/EPA call in either case, see closed_form_contacts::refineCylinderContact).
      uint64_t cylRefinedContact   = 0;
      uint64_t cylRefinedSeparated = 0;
      // Profiling of refineCylinderContact: pairs and wall time per exit (closed_form_contacts::RefinePath), why
      // the fast path rejected a closed-form contact (FastReject), how often the kink search / convex check ran,
      // and the closed-form scenario of the re-descended contacts (index as in shadow mode, 0..7).
      static constexpr int kRefinePaths = int(closed_form_contacts::RefinePath::Count);
      uint64_t cylRefinePath[kRefinePaths]        = {};
      double   cylRefinePathSeconds[kRefinePaths] = {};
      uint64_t cylRefineFastReject[int(closed_form_contacts::FastReject::Count)] = {};
      uint64_t cylRefineKinks = 0, cylRefineConvex = 0;
      uint64_t cylRefineDescendedScenario[8] = {};
      // Shadow timing (ANTIDUNES_REFINE_SHADOW_TIMING set, refined mode only): on the same pairs, the time of
      // the closed form and of the GJK/EPA fallback that the closed form with fallback would have called.
      uint64_t shadowFallbackCalls = 0;
      double   shadowClosedFormSeconds = 0, shadowFallbackSeconds = 0;

      void bumpCylScenario(const char* scenario)
      {
         const std::string s(scenario);
         if      (s == "face-face")           ++cylFaceFace;
         else if (s == "band-band-parallel")  ++cylBandBandParallel;
         else if (s == "face-band")           ++cylFaceBand;
         else if (s == "band-band-skewed")    ++cylBandBandSkewed;
         else if (s == "face-edge")           ++cylFaceEdge;
         else if (s == "special-contact")     ++cylSpecialContact;
         else if (s == "band-edge")           ++cylBandEdge;
         else if (s == "edge-edge")           ++cylEdgeEdge;
      }
   };
   ContactDetectionStats contactStats;

   // ---- SHADOW MODE (ModelParameters.shadowAnalyticCylinderLog): closed-form result computed next to the GJK/EPA
   // result that is actually used, never applied. Cumulative (NOT reset with contactStats).
   struct ShadowScenarioStat
   {
      // wrong: |overlap - GJK overlap| > 5 % of the smaller radius (GJK overlap = 0 where GJK sees no contact);
      // over/under = closed-form overlap minus GJK overlap (> 0 over-estimate)
      uint64_t n = 0, wrong = 0, verifierFlagged = 0, normalFlipped = 0, normalOff = 0;
      double   sumAbsErr = 0, maxUnder = 0, maxOver = 0;
   };
   // scenarios: face-face, band-band-parallel, face-band, band-band-skewed, face-edge, special-contact, band-edge,
   // edge-edge, none (closed form found nothing)
   static constexpr size_t kShadowScenarios = 9;
   struct ShadowStats
   {
      std::array< ShadowScenarioStat, kShadowScenarios > byScenario{};
      uint64_t gjkContacts = 0, bothContact = 0, closedFormOnly = 0, gjkOnly = 0;
   };
   static const char* shadowScenarioName(size_t i)
   {
      static const char* names[kShadowScenarios] = { "face-face", "band-band-parallel", "face-band", "band-band-skewed", "face-edge", "special-contact", "band-edge", "edge-edge", "none" };
      return names[i];
   }
   static size_t shadowScenarioIndex(const char* scenario)
   {
      for (size_t i = 0; i + 1 < kShadowScenarios; ++i) if (std::string(scenario) == shadowScenarioName(i)) return i;
      return kShadowScenarios - 1;
   }
   bool  shadowAnalyticLog = false;
   uint_t shadowStartStep = 0, shadowEndStep = 0;
   uint_t shadowMainStep = 0, shadowSubCycle = 0;
   ShadowStats shadowStats;
   std::shared_ptr< std::ofstream > shadowFile_;

   static double shadowSupport(const mesa_pd::Vec3& n, const mesa_pd::Vec3& axis, double radius, double halfLength)
   {
      const double c = std::abs(double(n * axis));
      return halfLength * c + radius * std::sqrt(std::max(0.0, 1.0 - c * c));
   }
   //! distance B must move along +-n (the smaller of the two orientations) to separate from A
   static double shadowNeeded(const mesa_pd::Vec3& n, const closed_form_contacts::CylinderState& a, const closed_form_contacts::CylinderState& b)
   {
      const mesa_pd::Vec3 d = b.pos - a.pos, ea = a.axis(), eb = b.axis();
      auto f = [&](const mesa_pd::Vec3& m) {
         return shadowSupport(m, ea, double(a.radius), 0.5 * double(a.length)) + shadowSupport(m, eb, double(b.radius), 0.5 * double(b.length)) - double(d * m); };
      return std::min(f(n), f(-n));
   }
   //! upper bound of the exact penetration depth: min of shadowNeeded over a few canonical directions
   static double shadowUpperBound(const mesa_pd::Vec3& reported, const closed_form_contacts::CylinderState& a, const closed_form_contacts::CylinderState& b)
   {
      const mesa_pd::Vec3 d = b.pos - a.pos, ea = a.axis(), eb = b.axis();
      std::array< mesa_pd::Vec3, 6 > cand{ reported, ea, eb, d, ea % eb, d - ea * (d * ea) };
      double ub = shadowNeeded(cand[0], a, b);
      for (size_t k = 1; k < cand.size(); ++k)
         if (cand[k].length() > real_t(1e-9)) ub = std::min(ub, shadowNeeded(cand[k] / cand[k].length(), a, b));
      const mesa_pd::Vec3 r2 = d - eb * (d * eb);
      if (r2.length() > real_t(1e-9)) ub = std::min(ub, shadowNeeded(r2 / r2.length(), a, b));
      return ub;
   }

   //! Called right after the GJK/EPA result (the one actually used) is known; reads it from
   //! getPenetrationDepth()/getContactNormal() if gjkHit.
   template< typename Accessor >
   void shadowCompare(const size_t idx1, const size_t idx2, const closed_form_contacts::CylinderState& ci, const closed_form_contacts::CylinderState& cj, const bool gjkHit, Accessor& ac)
   {
      const auto res = closed_form_contacts::cylinderCylinderContact(ci, cj);
      if (!gjkHit && !res.contact) return;
      const double tol = 0.05 * double(std::min(ci.radius, cj.radius));
      const double gjkOverlap = gjkHit ? -double(getPenetrationDepth()) : 0.0;
      const mesa_pd::Vec3 gjkN = gjkHit ? getContactNormal() : mesa_pd::Vec3(real_t(0), real_t(0), real_t(0));
      const size_t si = shadowScenarioIndex(res.contact ? res.scenario : "none");
      auto& st = shadowStats.byScenario[si];
      double err = 0.0, need = -1.0, ub = -1.0; bool flagged = false, flipped = false, normalOff = false;
      if (gjkHit) ++shadowStats.gjkContacts;
      if (gjkHit && res.contact) ++shadowStats.bothContact; else if (res.contact) ++shadowStats.closedFormOnly; else ++shadowStats.gjkOnly;
      if (res.contact)
      {
         err = double(res.overlap) - gjkOverlap;
         need = shadowNeeded(res.normal, ci, cj); ub = shadowUpperBound(res.normal, ci, cj);
         flagged = (need - double(res.overlap) > tol) || (double(res.overlap) - ub > tol);
         if (gjkHit) { const double dot = double(res.normal * gjkN); flipped = dot < 0.0; normalOff = !flipped && dot < 0.94; }
      }
      // closed form missed a contact GJK sees
      else err = -gjkOverlap;
      ++st.n; st.sumAbsErr += std::abs(err); st.maxUnder = std::max(st.maxUnder, -err); st.maxOver = std::max(st.maxOver, err);
      if (std::abs(err) > tol) ++st.wrong;
      if (flagged) ++st.verifierFlagged;
      if (flipped) ++st.normalFlipped;
      if (normalOff) ++st.normalOff;

      if (shadowEndStep != uint_t(0) && shadowMainStep >= shadowStartStep && shadowMainStep <= shadowEndStep)
      {
         if (!shadowFile_)
         {
            const int rank = mpi::MPIManager::instance()->rank();
            std::error_code ec; std::filesystem::create_directories("shadow-analytic", ec);
            shadowFile_ = std::make_shared< std::ofstream >("shadow-analytic/rank_" + std::to_string(rank) + ".csv", std::ios::out | std::ios::trunc);
            *shadowFile_ << "mainStep,subCycle,rank,uidA,uidB,gjkContact,gjkOverlap,gjkNx,gjkNy,gjkNz,closedFormContact,scenario,cfOverlap,cfNx,cfNy,cfNz,errCfMinusGjk,neededAlongCfNormal,upperBound,checkFlags,normalDot\n";
         }
         *shadowFile_ << shadowMainStep << ',' << shadowSubCycle << ',' << mpi::MPIManager::instance()->rank() << ',' << ac.getUid(idx1) << ',' << ac.getUid(idx2) << ','
                      << int(gjkHit) << ',' << gjkOverlap << ',' << gjkN[0] << ',' << gjkN[1] << ',' << gjkN[2] << ',' << int(res.contact) << ',' << (res.contact ? res.scenario : "none") << ','
                      << (res.contact ? double(res.overlap) : 0.0) << ',' << (res.contact ? double(res.normal[0]) : 0.0) << ',' << (res.contact ? double(res.normal[1]) : 0.0) << ',' << (res.contact ? double(res.normal[2]) : 0.0) << ','
                      << err << ',' << need << ',' << ub << ',' << int(flagged) << ',' << ((gjkHit && res.contact) ? double(res.normal * gjkN) : 0.0) << '\n';
         shadowFile_->flush();
      }
   }

   // Runtime switch (default false = GJK/EPA only). Public so ComputeContacts can set it at construction from
   // ModelParameters.useAnalyticalContactDetection.
   bool useAnalyticalContactDetection = false;

   // Runtime switch (default false), only meaningful with useAnalyticalContactDetection: removes the GJK/EPA
   // safety net for Cylinder x Cylinder and trusts the closed-form Kodam et al. (2010) + Guo et al. (2012)
   // dispatch as-is ("pure Kodam+Guo"). See ModelParameters.usePureAnalyticalCylinder (SimulationParams.h).
   bool usePureAnalyticalCylinder = false;
   // See ModelParameters.refineCylinderContacts.
   bool refineCylinderContacts = false;

   // Runtime switch (default true), see ModelParameters.useTightLayer1 (SimulationParams.h): tight,
   // orientation-aware "too far apart to touch" pre-filters. Cylinder x Cylinder: the axis-segment distance bound
   // replaces the circumscribing-sphere sum, regardless of useAnalyticalContactDetection. Box x Box: an extra
   // single-axis (centre line) projection reject ahead of boxBoxSAT, so only when useAnalyticalContactDetection
   // is true.
   bool useTightLayer1 = true;

   // See ModelParameters::useRotationalEffectiveMass (SimulationParams.h). Changes the force-law calibration, so
   // it is read by ComputeContacts where it calls collisionResponse_.setStiffnessAndDamping(...), the point at
   // which the generalized effective mass replaces the translational one.
   bool useRotationalEffectiveMass = false;

   // See ModelParameters::cylinderContactSkin (SimulationParams.h). 0 = off. Temporarily set to 0 around the
   // recursive inflated-shape call in the Cylinder x Cylinder / Cylinder x HalfSpace overloads above.
   real_t cylinderContactSkin = real_t(0);
   // See ModelParameters::boxContactSkin (SimulationParams.h). 0 = off. Box analogue of cylinderContactSkin,
   // toggled the same way in the Box x Box / Box x HalfSpace overloads above.
   real_t boxContactSkin = real_t(0);

   // Secondary contact point (arch/tunnel containment only, Box/Cylinder floaters): the containment
   // farthest-point search can find two independent points poking past the tube radius at once (an elongated
   // floater wedged with a leading and a trailing corner, or both end-cap rims, out). The base interface carries
   // one contact per call, so the worse point takes the normal path (full LinearSpringDashpot response via
   // ComputeContacts); the second is exposed here so ComputeContacts can give it a simplified one-sided
   // response instead of discarding it.
   bool hasSecondaryContact() const { return hasSecondaryContact_; }
   const mesa_pd::Vec3& getSecondaryContactPoint() const { return secondaryContactPoint_; }
   const mesa_pd::Vec3& getSecondaryContactNormal() const { return secondaryContactNormal_; }
   real_t getSecondaryPenetrationDepth() const { return secondaryPenetrationDepth_; }
   // Called once per pair by ComputeContacts before dispatching: this object is reused for every pair in a
   // sub-cycle, so a stale hasSecondaryContact_ must not leak into a pair whose branch (anything other than the
   // two containment branches above) never touches these fields.
   void resetSecondaryContact() { hasSecondaryContact_ = false; }

private:
   // Backing storage for recordDetectionPath()/lookupDetectionPath() above.
   std::map<std::pair<walberla::id_t, walberla::id_t>, std::string> pairDetectionPath;

   // Backing storage for the secondary-contact getters above.
   bool          hasSecondaryContact_       = false;
   mesa_pd::Vec3 secondaryContactPoint_;
   mesa_pd::Vec3 secondaryContactNormal_;
   real_t        secondaryPenetrationDepth_ = real_t(0);
};


// Per-rank CSV writer for ModelParameters.logFloaterHistory (see SimulationParams.h). Owned by ParticleKernels
// and shared by pointer with ComputeContacts and ComputeLubrication, the places where a contact/lubrication pair
// is computed and its force applied. Disabled (enabled_ == false) is the zero-cost default: every call below
// returns early before touching the file.
class FloaterHistoryLogger {
public:
   // Called once, after ModelParameters is known, before the main timeloop. baseFolder must already exist
   // (Antidunes.cpp creates it on root); every rank opens (and truncates) its own file. No barrier needed:
   // collective MPI calls (setupBlockForest etc.) already lie between directory creation and this call.
   void init(const std::string& baseFolder, bool enabled)
   {
      enabled_ = enabled;
      if (!enabled_) return;
      const int rank = mpi::MPIManager::instance()->rank();
      file_.open(baseFolder + "/rank_" + std::to_string(rank) + ".csv", std::ios::out | std::ios::trunc);
      WALBERLA_CHECK(file_.is_open(), "[FloaterHistoryLogger] failed to open CSV for rank " << rank
         << " in \"" << baseFolder << "\" -- does the directory exist?");
      file_ << "mainStep,subCycle,kind,computingRank,uidA,uidB,ownerRankA,ownerRankB,isGhostA,isGhostB,"
               "shapeTypeA,shapeTypeB,detectionPath,contactX,contactY,contactZ,normalX,normalY,normalZ,"
               "penetrationDepth,overlapMagnitude,forceOnA_x,forceOnA_y,forceOnA_z,forceOnB_x,forceOnB_y,forceOnB_z\n";
      file_.flush();

      // Second file, same toggle: a per-floater world-state snapshot (not per-contact) of every non-dormant
      // sphere/box/cylinder floater this rank holds (local and ghost copies), taken right after each sub-cycle's
      // first sync in doDEMSubCycles. Comparing two ranks' rows for the same uid at the same (mainStep,subCycle) shows
      // sync staleness/divergence (e.g. a ghost copy lagging the owner), whether or not a contact was computed.
      stateFile_.open(baseFolder + "/rank_" + std::to_string(rank) + "_state.csv",
                       std::ios::out | std::ios::trunc);
      WALBERLA_CHECK(stateFile_.is_open(), "[FloaterHistoryLogger] failed to open state CSV for rank " << rank
         << " in \"" << baseFolder << "\" -- does the directory exist?");
      // hydroForce/hydroTorque: the raw MEM hydrodynamic force/torque (ac.getHydrodynamicForce/Torque),
      // computed once per main LBM step by reduceAndAvgHydroForces() before the DEM sub-cycles, i.e. the value
      // AddExternalForces checks for NaN/Inf. Unlike position/velocity/rotation it is not part of
      // ParticleGhostCopyNotification, so it is only valid on the owning rank's row (isGhost=0); ghost rows hold
      // stale/default values and must not be compared across ranks.
      stateFile_ << "mainStep,subCycle,observingRank,uid,ownerRank,isGhost,isLocal,shapeType,"
                    "posX,posY,posZ,velX,velY,velZ,angVelX,angVelY,angVelZ,forceX,forceY,forceZ,"
                    "torqueX,torqueY,torqueZ,hydroForceX,hydroForceY,hydroForceZ,"
                    "hydroTorqueX,hydroTorqueY,hydroTorqueZ,fluidDensity,quatR,quatI,quatJ,quatK\n";
      stateFile_.flush();
   }

   bool enabled() const { return enabled_; }

   // nullptr (default) -> every floater is logged. Non-null -> log() and logState() skip rows whose particle(s)
   // are not in this set (ModelParams::logFloaterHistoryReleaseStepFilter: restricts the heavy per-sub-cycle
   // snapshot to one release batch). Not owned; the caller (ParticleKernels) keeps the set alive.
   void setUidFilter(const std::unordered_set< walberla::id_t >* uidFilter) { uidFilter_ = uidFilter; }

   // 0 (default) -> log from step 0. Otherwise skip earlier steps (ModelParams::logFloaterHistoryStartStep), so a
   // targeted rerun logs only the window of interest.
   void setMinStep(uint_t minStep) { minStep_ = minStep; }

   // 0 (default, converted to the max uint_t sentinel) -> log to the end of the run. Nonzero -> upper bound
   // counterpart of setMinStep (ModelParams::logFloaterHistoryEndStep).
   void setMaxStep(uint_t maxStep)
   {
      maxStep_ = (maxStep == uint_t(0)) ? std::numeric_limits< uint_t >::max() : maxStep;
   }

   // kind: "contact" (ComputeContacts) or "lubrication" (ComputeLubrication). detectionPath names the algorithm
   // that produced the row: the string recorded by SolidCylinderContactDetection::recordDetectionPath (Box-Box /
   // Cylinder-Cylinder dispatch, pillar/containment overloads) if any; otherwise, and for lubrication rows, a
   // fixed label from the caller. forceOnA/forceOnB are the force DELTA of this pair in this sub-cycle (the caller
   // reads ac.getForceRef before/after, since the kernels add onto earlier contacts' accumulated force).
   template< typename ParticleAccessor_T >
   void log(uint_t mainStep, uint_t subCycle, const char* kind, const ParticleAccessor_T& ac,
            size_t idx1, size_t idx2, const mesa_pd::Vec3& contactPoint, const mesa_pd::Vec3& normal,
            real_t penetrationDepth, const mesa_pd::Vec3& forceOnA, const mesa_pd::Vec3& forceOnB,
            const std::string& detectionPath)
   {
      if (!enabled_ || mainStep < minStep_ || mainStep > maxStep_) return;
      if (uidFilter_ && !uidFilter_->count(ac.getUid(idx1)) && !uidFilter_->count(ac.getUid(idx2))) return;
      using namespace mesa_pd::data::particle_flags;
      file_ << mainStep << ',' << subCycle << ',' << kind << ',' << mpi::MPIManager::instance()->rank() << ','
            << ac.getUid(idx1) << ',' << ac.getUid(idx2) << ','
            << ac.getOwner(idx1) << ',' << ac.getOwner(idx2) << ','
            << isSet(ac.getFlags(idx1), GHOST) << ',' << isSet(ac.getFlags(idx2), GHOST) << ','
            << (ac.getBaseShape(idx1) ? ac.getBaseShape(idx1)->getShapeType() : uint_t(999)) << ','
            << (ac.getBaseShape(idx2) ? ac.getBaseShape(idx2)->getShapeType() : uint_t(999)) << ','
            << '"' << detectionPath << '"' << ','
            << contactPoint[0] << ',' << contactPoint[1] << ',' << contactPoint[2] << ','
            << normal[0] << ',' << normal[1] << ',' << normal[2] << ',' << penetrationDepth << ','
            << (-penetrationDepth) << ','
            << forceOnA[0] << ',' << forceOnA[1] << ',' << forceOnA[2] << ','
            << forceOnB[0] << ',' << forceOnB[1] << ',' << forceOnB[2] << '\n';
      // Flush every row: the run can abort (WALBERLA_ABORT/MPI_ABORT) at any point, and the final rows before a
      // crash are the most important ones.
      file_.flush();
   }

   // One row per floater particle this rank holds, local or ghost. Called from doDEMSubCycles right after each
   // sub-cycle's first sync; deliberately not restricted to SelectLocal(), so ghost
   // copies are captured for the staleness comparison described above.
   // fluidDensity: rho at the floater position, sampled by the caller (ParticleKernels::sampleFluidDensityAt);
   // the quantity the [LOW-RHO] guard in ReconstructInterfaceCellABB.h checks. Only meaningful for isLocal=1
   // rows; for a position outside this rank's blocks the caller passes NaN.
   template< typename ParticleAccessor_T >
   void logState(uint_t mainStep, uint_t subCycle, const ParticleAccessor_T& ac, size_t idx, real_t fluidDensity)
   {
      if (!enabled_ || mainStep < minStep_ || mainStep > maxStep_) return;
      if (uidFilter_ && !uidFilter_->count(ac.getUid(idx))) return;
      using namespace mesa_pd::data::particle_flags;
      const auto& vel         = ac.getLinearVelocity(idx);
      const auto& angVel      = ac.getAngularVelocity(idx);
      const auto& force       = ac.getForce(idx);
      const auto& torque      = ac.getTorque(idx);
      const auto& pos         = ac.getPosition(idx);
      const auto& hydroForce  = ac.getHydrodynamicForce(idx);
      const auto& hydroTorque = ac.getHydrodynamicTorque(idx);
      // Quaternion (r,i,j,k), body->world: lets an offline check reconstruct each floater's oriented footprint
      // (box corners / cylinder axis) at any logged instant.
      const auto& quat = ac.getRotation(idx).getQuaternion();
      stateFile_ << mainStep << ',' << subCycle << ',' << mpi::MPIManager::instance()->rank() << ','
                 << ac.getUid(idx) << ',' << ac.getOwner(idx) << ','
                 << isSet(ac.getFlags(idx), GHOST) << ',' << !isSet(ac.getFlags(idx), GHOST) << ','
                 << (ac.getBaseShape(idx) ? ac.getBaseShape(idx)->getShapeType() : uint_t(999)) << ','
                 << pos[0] << ',' << pos[1] << ',' << pos[2] << ','
                 << vel[0] << ',' << vel[1] << ',' << vel[2] << ','
                 << angVel[0] << ',' << angVel[1] << ',' << angVel[2] << ','
                 << force[0] << ',' << force[1] << ',' << force[2] << ','
                 << torque[0] << ',' << torque[1] << ',' << torque[2] << ','
                 << hydroForce[0] << ',' << hydroForce[1] << ',' << hydroForce[2] << ','
                 << hydroTorque[0] << ',' << hydroTorque[1] << ',' << hydroTorque[2] << ','
                 << fluidDensity << ','
                 << quat[0] << ',' << quat[1] << ',' << quat[2] << ',' << quat[3] << '\n';
      stateFile_.flush(); // same crash-survivability rationale as log() above
   }

   // One-sided companion to log() for the archived ComputeTunnelContacts (archive/arch_tunnel_containment/; arch
   // tunnels are not mesa_pd particles); no caller in the active code. Same CSV
   // columns (kind="tunnel"): uidB is the sentinel 900000000 + tunnelIndex (real UIDs are far larger), ownerRankB
   // is -1, shapeTypeB is 997 (not a valid mesa_pd SHAPE_TYPE), so "shapeTypeB==997" selects tunnel-wall rows.
   // forceOnB = -forceOnA (reaction on the infinite-mass wall, for readability).
   template< typename ParticleAccessor_T >
   void logTunnelContact(uint_t mainStep, uint_t subCycle, const ParticleAccessor_T& ac, size_t idx,
      size_t tunnelIndex, const mesa_pd::Vec3& contactPoint, const mesa_pd::Vec3& normal,
      real_t penetrationDepth, const mesa_pd::Vec3& forceOnFloater)
   {
      if (!enabled_ || mainStep < minStep_ || mainStep > maxStep_) return;
      if (uidFilter_ && !uidFilter_->count(ac.getUid(idx))) return;
      using namespace mesa_pd::data::particle_flags;
      file_ << mainStep << ',' << subCycle << ',' << "tunnel" << ',' << mpi::MPIManager::instance()->rank() << ','
            << ac.getUid(idx) << ',' << (uint_t(900000000) + tunnelIndex) << ','
            << ac.getOwner(idx) << ',' << -1 << ','
            << isSet(ac.getFlags(idx), GHOST) << ',' << 0 << ','
            << (ac.getBaseShape(idx) ? ac.getBaseShape(idx)->getShapeType() : uint_t(999)) << ',' << 997 << ','
            << '"' << "one-sided arch-tunnel containment (ComputeTunnelContacts, tunnel index "
            << tunnelIndex << ")" << '"' << ','
            << contactPoint[0] << ',' << contactPoint[1] << ',' << contactPoint[2] << ','
            << normal[0] << ',' << normal[1] << ',' << normal[2] << ',' << penetrationDepth << ','
            << (-penetrationDepth) << ','
            << forceOnFloater[0] << ',' << forceOnFloater[1] << ',' << forceOnFloater[2] << ','
            << -forceOnFloater[0] << ',' << -forceOnFloater[1] << ',' << -forceOnFloater[2] << '\n';
      // Flush every row (crash survivability, see log()).
      file_.flush();
   }

private:
   bool enabled_ = false;
   std::ofstream file_;
   std::ofstream stateFile_;
   const std::unordered_set< walberla::id_t >* uidFilter_ = nullptr;
   uint_t minStep_ = uint_t(0);
   uint_t maxStep_ = std::numeric_limits< uint_t >::max();
};


//! Diagnostic artificial caps for cylinder floaters (ModelParameters.rodCap*), per rank; see SimulationParams.h.
//! All zero = off.
struct RodCaps { real_t hydroForce = real_t(0), hydroTorque = real_t(0), velocity = real_t(0), angularVelocity = real_t(0), overlap = real_t(0);
                 bool any() const { return hydroForce > real_t(0) || hydroTorque > real_t(0) || velocity > real_t(0) || angularVelocity > real_t(0) || overlap > real_t(0); } };
struct RodCapStats { uint64_t hydroForce = 0, hydroTorque = 0, velocity = 0, angularVelocity = 0, overlap = 0, nonFinite = 0;
                     real_t maxHydroForce = real_t(0), maxHydroTorque = real_t(0), maxVelocity = real_t(0), maxAngularVelocity = real_t(0), maxOverlap = real_t(0); };
inline RodCaps     g_rodCaps;
inline RodCapStats g_rodCapStats;
//! ModelParameters.rodPrescribedVelocityX/Y/Z (test option, default off): constant rigid velocity imposed on all
//! sphere/box/cylinder floaters after every sub-cycle's velocity update
inline Vector3< real_t > g_rodPrescribedVelocity{ real_t(0) };
inline bool              g_rodPrescribedVelocityOn = false;
inline void logRodCapHit(const char* what, walberla::id_t uid, uint_t step, uint_t sub, real_t before, real_t cap)
{
   static uint_t logCount = uint_t(0);
   if (logCount < uint_t(40)) { ++logCount; WALBERLA_LOG_WARNING("[ROD-CAP " << what << "] uid=" << uid << " step=" << step << " subCycle=" << sub << " value " << before << " -> capped to " << cap); }
}

class ComputeLubrication {
public:
   ComputeLubrication(const std::shared_ptr< mesa_pd::domain::BlockForestDomain >& domain, real_t viscosity,
      bool useContactFilterOwnerTiebreak = false)
      : domain_(domain)
      , lubricationCorrectionKernel_(viscosity,
         [](real_t r) { return (real_c(0.001) + real_c(0.00007) * r) * r; })
   {
      contact_filter_.useOwnerRankTiebreak = useContactFilterOwnerTiebreak;
   }

   // See ComputeContacts::setHistoryLogger/setStepInfo -- same wiring, same rationale.
   void setHistoryLogger(FloaterHistoryLogger* logger) { historyLogger_ = logger; }
   void setStepInfo(uint_t mainStep, uint_t subCycle) { mainStep_ = mainStep; subCycle_ = subCycle; }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx1, const size_t idx2, ParticleAccessor_T& ac){
      // lubrication is only defined for sphere-sphere pairs
      if (ac.getBaseShape(idx1)->getShapeType() != mesa_pd::data::Sphere::SHAPE_TYPE ||
          ac.getBaseShape(idx2)->getShapeType() != mesa_pd::data::Sphere::SHAPE_TYPE) return;

      acd_.getContactThreshold() = lubricationCorrectionKernel_.getNormalCutOffDistance();
      if (double_cast_(idx1, idx2, ac, acd_, ac))
      {
         if (contact_filter_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(), *domain_))
         {
            const bool doLog = historyLogger_ && historyLogger_->enabled();
            const mesa_pd::Vec3 forceBeforeA = doLog ? ac.getForceRef(acd_.getIdx1()) : mesa_pd::Vec3(0_r);
            const mesa_pd::Vec3 forceBeforeB = doLog ? ac.getForceRef(acd_.getIdx2()) : mesa_pd::Vec3(0_r);

            double_cast_(acd_.getIdx1(), acd_.getIdx2(), ac, lubricationCorrectionKernel_, ac,
                        acd_.getContactNormal(), acd_.getPenetrationDepth());

            if (doLog)
            {
               historyLogger_->log(mainStep_, subCycle_, "lubrication", ac, acd_.getIdx1(), acd_.getIdx2(),
                  acd_.getContactPoint(), acd_.getContactNormal(), acd_.getPenetrationDepth(),
                  ac.getForceRef(acd_.getIdx1()) - forceBeforeA, ac.getForceRef(acd_.getIdx2()) - forceBeforeB,
                  "lubrication (LubricationCorrectionKernel, sphere-sphere exact analytic geometry, "
                  "no analytical/GJK choice)");
            }
         }
      }
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain > domain_;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel_;
   mesa_pd::collision_detection::AnalyticContactDetection acd_;
   mesa_pd::kernel::DoubleCast double_cast_;
   mesa_pd::mpi::ContactFilter contact_filter_;

   FloaterHistoryLogger* historyLogger_ = nullptr;
   uint_t mainStep_ = uint_t(0);
   uint_t subCycle_ = uint_t(0);
};

// Generalized (translational + rotational) effective mass for a contact at `contactPoint` with unit normal
// `contactNormal` between bodies idx1/idx2 (see ModelParameters::useRotationalEffectiveMass, SimulationParams.h):
//     1/meff = 1/mA + 1/mB + (rA x n)^T IA_world^-1 (rA x n) + (rB x n)^T IB_world^-1 (rB x n)
// with rA/rB the contact point's offset from each centre of mass. The inverse inertia is stored in the body frame
// (getInvInertiaBF) and rotated to world with math::transformMatrixRART, as in mesa_pd::kernel::VelocityVerlet.
// For a fixed/infinite-mass body (walls, pillars) getInvMass and getInvInertiaBF are zero, so its terms vanish.
// For sphere-sphere contact rA and rB are parallel to the normal, so this reduces exactly to the translational
// formula.
template < typename Accessor >
inline real_t computeRotationalEffectiveMass(Accessor& ac, size_t idx1, size_t idx2,
                                              const mesa_pd::Vec3& contactPoint,
                                              const mesa_pd::Vec3& contactNormal)
{
   const mesa_pd::Vec3 rA = contactPoint - ac.getPosition(idx1);
   const mesa_pd::Vec3 rB = contactPoint - ac.getPosition(idx2);

   const mesa_pd::Vec3 crossA = rA % contactNormal;
   const mesa_pd::Vec3 crossB = rB % contactNormal;

   real_t rotationalTerm = real_t(0);
   if (crossA.sqrLength() > real_t(0))
   {
      const mesa_pd::Mat3 invInertiaWF_A =
         math::transformMatrixRART(ac.getRotation(idx1).getMatrix(), ac.getInvInertiaBF(idx1));
      rotationalTerm += crossA * (invInertiaWF_A * crossA);
   }
   if (crossB.sqrLength() > real_t(0))
   {
      const mesa_pd::Mat3 invInertiaWF_B =
         math::transformMatrixRART(ac.getRotation(idx2).getMatrix(), ac.getInvInertiaBF(idx2));
      rotationalTerm += crossB * (invInertiaWF_B * crossB);
   }

   const real_t invMeff = ac.getInvMass(idx1) + ac.getInvMass(idx2) + rotationalTerm;
   WALBERLA_ASSERT_GREATER(invMeff, real_t(0));
   return real_t(1) / invMeff;
}

class ComputeContacts {
public:
   ComputeContacts(const std::shared_ptr< mesa_pd::domain::BlockForestDomain >& domain,
      real_t restitutionCoefficient, real_t kappa,
      real_t collisionTime, real_t particleTimestepSize, real_t frictionCoefficient,
      uint_t numParticleTypes = uint_t(2), bool useAnalyticalContactDetection = false,
      bool usePureAnalyticalCylinder = false, bool useContactFilterOwnerTiebreak = false,
      bool useTightLayer1 = true, bool useRotationalEffectiveMass = false)
      : domain_(domain)
      , collisionResponse_(numParticleTypes)
      , restitutionCoefficient_(restitutionCoefficient)
      , kappa_(kappa)
      , collisionTime_(collisionTime)
      , particleTimestepSize_(particleTimestepSize)
   {
      for (uint_t t = 0; t < numParticleTypes; ++t)
         for (uint_t u = t; u < numParticleTypes; ++u)
            collisionResponse_.setFrictionCoefficientDynamic(t, u, frictionCoefficient);
      acd_.useAnalyticalContactDetection = useAnalyticalContactDetection;
      acd_.usePureAnalyticalCylinder     = usePureAnalyticalCylinder;
      acd_.useRotationalEffectiveMass    = useRotationalEffectiveMass;
      acd_.useTightLayer1                = useTightLayer1;
      contact_filter_.useOwnerRankTiebreak = useContactFilterOwnerTiebreak;

      // Wire the optional ContactHistoryNotification diagnostic hook to acd_'s per-pair detection-path record,
      // so a [ContactHistory DUPLICATE] warning can report whether the contact came from the closed form or the
      // GJK/EPA fallback. Captures `this`, which outlives the simulation.
      mesa_pd::contactHistoryDuplicateDiagnosticHook = [this](walberla::id_t a, walberla::id_t b) -> std::string {
         const std::string* path = acd_.lookupDetectionPath(a, b);
         std::string result = path ? *path
                                    : std::string("not tracked locally on this rank this cycle "
                                                   "(not a Cylinder x Cylinder pair, or computed on a different rank)");
         // This rank's local view of both particles at the ContactFilter decision, if it computed this pair's
         // contact this cycle (see pairContactInputState_).
         if (const ContactInputState* state = lookupContactInputState(a, b))
            result += " || input state: " + state->format();
         return result;
      };
   }

   // Read-only access to acd_'s per-rank scenario/fallback tallies, for
   // ParticleKernels::logContactDetectionStatistics() to MPI-sum and log.
   const SolidCylinderContactDetection::ContactDetectionStats& contactStats() const { return acd_.contactStats; }

   // Zeroes acd_'s per-rank tallies, for ParticleKernels::resetContactDetectionStatistics(), so a periodic log
   // can report per-interval instead of cumulative counts.
   void resetContactStats() { acd_.contactStats = SolidCylinderContactDetection::ContactDetectionStats{}; }

   // Wires in ModelParameters.logFloaterHistory's CSV writer (nullptr/disabled by default, see
   // FloaterHistoryLogger); called once from the ParticleKernels constructor. setStepInfo is called once per
   // sub-cycle from doDEMSubCycles.
   void setHistoryLogger(FloaterHistoryLogger* logger) { historyLogger_ = logger; }
   void setStepInfo(uint_t mainStep, uint_t subCycle) { mainStep_ = mainStep; subCycle_ = subCycle; acd_.shadowMainStep = mainStep; acd_.shadowSubCycle = subCycle; }

   // See ModelParameters::cylinderContactSkin (SimulationParams.h).
   void setCylinderContactSkin(real_t skin) { acd_.cylinderContactSkin = skin; }
   // See ModelParameters::boxContactSkin (SimulationParams.h).
   void setBoxContactSkin(real_t skin) { acd_.boxContactSkin = skin; }
   void setRefineCylinderContacts(bool on) { acd_.refineCylinderContacts = on; }
   void setUsePureAnalyticalCylinder(bool on) { acd_.usePureAnalyticalCylinder = on; }

   // Shadow mode (ModelParameters.shadowAnalyticCylinderLog): see SolidCylinderContactDetection::shadowCompare.
   void setShadowAnalytic(bool enabled, uint_t startStep, uint_t endStep)
   {
      acd_.shadowAnalyticLog = enabled; acd_.shadowStartStep = startStep; acd_.shadowEndStep = endStep;
   }
   const SolidCylinderContactDetection::ShadowStats& shadowStats() const { return acd_.shadowStats; }

   // Debug hook for release-time overlaps (see SimulationObjects::freshlyActivatedFloaterUIDs); nullptr/disabled
   // by default. A contact resolved below that involves a UID in this set proves a floater was released while
   // already overlapping an active neighbour.
   void setFreshlyActivatedUIDs(std::shared_ptr< std::unordered_set< walberla::id_t > > uids)
   {
      freshlyActivatedUIDs_ = std::move(uids);
   }

   // See ModelParams::logContactAnomalies (SimulationParams.h) for the three watchers this enables
   // ([NORMAL-FLIP]/[PROLONGED-CONTACT]/[SIMULTANEOUS-CONTACTS]). endStep == 0 means unbounded, as in
   // FloaterHistoryLogger::setMaxStep.
   void setLogContactAnomalies(bool enabled, uint_t startStep, uint_t endStep)
   {
      logContactAnomalies_     = enabled;
      contactAnomalyStartStep_ = startStep;
      contactAnomalyEndStep_   = (endStep == uint_t(0)) ? std::numeric_limits< uint_t >::max() : endStep;
   }

   // Called once per sub-cycle by doDEMSubCycles after the contact sweep. [SIMULTANEOUS-CONTACTS] reports any
   // floater resolved against kSimultaneousContactsThreshold (3) or more other bodies in the same sub-cycle:
   // the independently computed spring-dashpot forces are summed, so the effective stiffness may exceed what the
   // explicit integrator's timestep was tuned for (dense packing). The tally is cleared even when logging is
   // off, so the map does not grow.
   void reportAndResetSimultaneousContacts()
   {
      if (logContactAnomalies_ && mainStep_ >= contactAnomalyStartStep_ && mainStep_ <= contactAnomalyEndStep_)
      {
         constexpr uint_t kSimultaneousContactsThreshold = uint_t(3);
         static uint_t simultaneousContactsLogCount = uint_t(0);
         for (const auto& entry : contactCountThisSubCycle_)
         {
            if (entry.second >= kSimultaneousContactsThreshold && simultaneousContactsLogCount < uint_t(200))
            {
               WALBERLA_LOG_WARNING("[SIMULTANEOUS-CONTACTS] step=" << mainStep_ << " subCycle=" << subCycle_
                  << " uid=" << entry.first << " activeContactsThisSubCycle=" << entry.second
                  << " -- this floater was resolved against " << entry.second << " other bodies within "
                     "the SAME sub-cycle; each contact's spring-dashpot force is computed independently "
                     "and summed, so the EFFECTIVE combined stiffness on this body may exceed what the "
                     "explicit integrator's timestep was tuned to handle for a single dominant contact.");
               ++simultaneousContactsLogCount;
            }
         }
      }
      contactCountThisSubCycle_.clear();
   }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx1, const size_t idx2, ParticleAccessor_T& ac){
      // Pairs without an overload in SolidCylinderContactDetection (e.g. Cylinder-Box) use the GJK/EPA path of
      // the GeneralContactDetection base, so floating boxes/cylinders collide with each other.
      acd_.resetSecondaryContact(); // acd_ is reused for every pair, see resetSecondaryContact()
      if (double_cast_(idx1, idx2, ac, acd_, ac))
      {
         // Capture this rank's local input state for the ContactFilter decision, keyed by uid pair (see
         // pairContactInputState_); surfaced only if a [ContactHistory DUPLICATE] fires for this pair. Raw values
         // are stored and formatted only on lookup, since formatting every contact is costly in dense piles.
         {
            const auto i1 = acd_.getIdx1();
            const auto i2 = acd_.getIdx2();
            auto side = [&](size_t i) {
               return ContactInputSide{ ac.getUid(i), ac.getOwner(i),
                                        mesa_pd::data::particle_flags::isSet(ac.getFlags(i), mesa_pd::data::particle_flags::GHOST),
                                        ac.getPosition(i), ac.getRotation(i).getQuaternion(), ac.getLinearVelocity(i),
                                        ac.getAngularVelocity(i) };
            };
            recordContactInputState(ac.getUid(i1), ac.getUid(i2),
                                    ContactInputState{ walberla::mpi::MPIManager::instance()->rank(), side(i1), side(i2),
                                                       acd_.getContactPoint(), acd_.getContactNormal(),
                                                       acd_.getPenetrationDepth() });
         }
         if (contact_filter_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(), *domain_))
         {
            // Detect degenerate contact geometry that would corrupt DEM integration
            const auto& n = acd_.getContactNormal();
            if (std::isnan(n[0]) || std::isnan(n[1]) || std::isnan(n[2]))
            {
               const auto s1 = ac.getBaseShape(acd_.getIdx1())
                                  ? ac.getBaseShape(acd_.getIdx1())->getShapeType() : uint_t(999);
               const auto s2 = ac.getBaseShape(acd_.getIdx2())
                                  ? ac.getBaseShape(acd_.getIdx2())->getShapeType() : uint_t(999);
               WALBERLA_LOG_WARNING_ON_ROOT(
                  "[ComputeContacts] NaN contact normal detected! "
                  "type(" << ac.getType(acd_.getIdx1()) << ",shapeType=" << s1 << ") pos=" << ac.getPosition(acd_.getIdx1()) << " <-> "
                  "type(" << ac.getType(acd_.getIdx2()) << ",shapeType=" << s2 << ") pos=" << ac.getPosition(acd_.getIdx2()) << " | "
                  "normal=" << n << " penetration=" << acd_.getPenetrationDepth());
               return; // skip this contact to avoid NaN propagation
            }

            // Debug hook for release-time overlaps (see SimulationObjects::freshlyActivatedFloaterUIDs). UIDs are
            // in freshlyActivatedUIDs_ only during the first doDEMSubCycles() call after
            // releaseSettledFloatingSpheres() unset FIXED.
            if (freshlyActivatedUIDs_)
            {
               const walberla::id_t uid1 = ac.getUid(acd_.getIdx1());
               const walberla::id_t uid2 = ac.getUid(acd_.getIdx2());
               if (freshlyActivatedUIDs_->count(uid1) > 0 || freshlyActivatedUIDs_->count(uid2) > 0)
               {
                  const auto s1 = ac.getBaseShape(acd_.getIdx1())
                                     ? ac.getBaseShape(acd_.getIdx1())->getShapeType() : uint_t(999);
                  const auto s2 = ac.getBaseShape(acd_.getIdx2())
                                     ? ac.getBaseShape(acd_.getIdx2())->getShapeType() : uint_t(999);
                  // Only a contact at sub-cycle 0 can be an insertion overlap: that detection runs on the
                  // positions the floaters had when released. A contact at a later sub-cycle of that step means
                  // they moved or rotated within it, so log both bodies' full state (throttled) instead and
                  // handle the contact like any other.
                  if (subCycle_ == uint_t(0))
                  {
                     WALBERLA_ABORT(
                        "[RELEASE-OVERLAP] Real contact detected (mainStep=" << mainStep_ << " subCycle="
                        << subCycle_ << ") involving a just-released floater: "
                        "uid=" << uid1 << " (shapeType=" << s1 << ", pos=" << ac.getPosition(acd_.getIdx1())
                        << ", freshlyActivated=" << (freshlyActivatedUIDs_->count(uid1) > 0) << ")  <->  "
                        "uid=" << uid2 << " (shapeType=" << s2 << ", pos=" << ac.getPosition(acd_.getIdx2())
                        << ", freshlyActivated=" << (freshlyActivatedUIDs_->count(uid2) > 0) << "). "
                        "penetration=" << acd_.getPenetrationDepth() << " normal=" << n
                        << " contactPoint=" << acd_.getContactPoint() << ". "
                        "This is the release-time insertion-overlap mechanism behind the 'flying floater' "
                        "launches investigated in this scenario -- increase the release-timing stagger "
                        "between floaters sharing this spawn point.");
                  }
                  else
                  {
                     static uint_t laterContactLogCount = uint_t(0);
                     if (laterContactLogCount < uint_t(40))
                     {
                        ++laterContactLogCount;
                        auto describe = [&ac](const size_t i) {
                           using namespace mesa_pd::data::particle_flags;
                           std::ostringstream o;
                           o << "uid=" << ac.getUid(i) << " shapeType="
                             << (ac.getBaseShape(i) ? ac.getBaseShape(i)->getShapeType() : uint_t(999))
                             << " fixed=" << isSet(ac.getFlags(i), FIXED) << " ghost=" << isSet(ac.getFlags(i), GHOST)
                             << " owner=" << ac.getOwner(i) << " pos=" << ac.getPosition(i)
                             << " quat=" << ac.getRotation(i).getQuaternion() << " vel=" << ac.getLinearVelocity(i)
                             << " angVel=" << ac.getAngularVelocity(i) << " force=" << ac.getForce(i)
                             << " torque=" << ac.getTorque(i) << " hydroForce=" << ac.getHydrodynamicForce(i)
                             << " hydroTorque=" << ac.getHydrodynamicTorque(i);
                           return o.str();
                        };
                        WALBERLA_LOG_WARNING("[RELEASE-CONTACT-LATER] contact at sub-cycle " << subCycle_ << " (mainStep="
                           << mainStep_ << ") of the FIRST step after a release, i.e. NOT an insertion overlap (sub-cycle 0 "
                           "was clean): penetration=" << acd_.getPenetrationDepth() << " normal=" << n << " contactPoint="
                           << acd_.getContactPoint() << " | body1: " << describe(acd_.getIdx1()) << " | body2: "
                           << describe(acd_.getIdx2()));
                     }
                  }
               }
            }

            // [LARGE-OVERLAP]: warns (no clamp) when a finite penetration exceeds half the smaller body's
            // interaction radius, which should never happen for separate convex bodies settling under
            // gravity/contact damping. Aimed at the GJK/EPA path (Box-Box, Cylinder-Cylinder, Cylinder-Box, ...),
            // which unlike the pillar overloads has no envelope clamp; the check itself runs for every contact.
            const real_t overlapMag  = -acd_.getPenetrationDepth();
            const real_t smallerRad  = std::min(ac.getInteractionRadius(acd_.getIdx1()),
                                                 ac.getInteractionRadius(acd_.getIdx2()));
            if (overlapMag > real_t(0.5) * smallerRad)
            {
               static uint_t largeOverlapLogCount = uint_t(0);
               if (largeOverlapLogCount < uint_t(20))
               {
                  const auto s1 = ac.getBaseShape(acd_.getIdx1())
                                     ? ac.getBaseShape(acd_.getIdx1())->getShapeType() : uint_t(999);
                  const auto s2 = ac.getBaseShape(acd_.getIdx2())
                                     ? ac.getBaseShape(acd_.getIdx2())->getShapeType() : uint_t(999);
                  WALBERLA_LOG_WARNING_ON_ROOT(
                     "[LARGE-OVERLAP] anomalous GJK/EPA penetration=" << acd_.getPenetrationDepth()
                     << " (smallerInteractionRadius=" << smallerRad << ") "
                     << "shapeType(" << s1 << ") pos=" << ac.getPosition(acd_.getIdx1()) << " <-> "
                     << "shapeType(" << s2 << ") pos=" << ac.getPosition(acd_.getIdx2())
                     << " normal=" << n);
                  ++largeOverlapLogCount;
               }
            }

            // [DEM-WATCH]: a single contact whose overlap exceeds kDemWatchOverlapRatio (6 %) of the smaller
            // interaction radius, i.e. a plausible launch-strength hit. Overlap-based rather than force-based, so
            // it is comparable across body sizes and tied to one contact. Calibration: glancing/settling contacts
            // sat at ~2-3 % (0.18-0.32 LU on a ~9.95 LU box), launch-causing hits at ~7-8 % (0.68-0.77 LU). Far
            // below [LARGE-OVERLAP]'s 50 %, which flags near-degenerate geometry.
            if (overlapMag > kDemWatchOverlapRatio * smallerRad)
            {
               static uint_t demWatchLogCount = uint_t(0);
               if (demWatchLogCount < uint_t(100))
               {
                  const auto s1 = ac.getBaseShape(acd_.getIdx1())
                                     ? ac.getBaseShape(acd_.getIdx1())->getShapeType() : uint_t(999);
                  const auto s2 = ac.getBaseShape(acd_.getIdx2())
                                     ? ac.getBaseShape(acd_.getIdx2())->getShapeType() : uint_t(999);
                  WALBERLA_LOG_WARNING("[DEM-WATCH] step=" << mainStep_ << " subCycle=" << subCycle_
                     << " deep contact: overlap=" << overlapMag << " (" << (overlapMag / smallerRad * real_t(100))
                     << "% of smallerInteractionRadius=" << smallerRad << ")"
                     << "  A: uid=" << ac.getUid(acd_.getIdx1()) << " shapeType=" << s1
                     << " pos=" << ac.getPosition(acd_.getIdx1()) << " vel=" << ac.getLinearVelocity(acd_.getIdx1())
                     << "  |  B: uid=" << ac.getUid(acd_.getIdx2()) << " shapeType=" << s2
                     << " pos=" << ac.getPosition(acd_.getIdx2()) << " vel=" << ac.getLinearVelocity(acd_.getIdx2())
                     << "  |  normal=" << n << " contactPoint=" << acd_.getContactPoint());
                  ++demWatchLogCount;
               }
            }

            // [NORMAL-FLIP] / [PROLONGED-CONTACT] / [SIMULTANEOUS-CONTACTS]: diagnostics for overlap building up
            // over many sub-cycles before a contact explodes. Off by default (logContactAnomalies_) and only active
            // in [contactAnomalyStartStep_, contactAnomalyEndStep_]; see ModelParams::logContactAnomalies.
            if (logContactAnomalies_ && mainStep_ >= contactAnomalyStartStep_ && mainStep_ <= contactAnomalyEndStep_)
            {
               const walberla::id_t uidA = ac.getUid(acd_.getIdx1());
               const walberla::id_t uidB = ac.getUid(acd_.getIdx2());
               const auto pairKey = uidA < uidB ? std::make_pair(uidA, uidB) : std::make_pair(uidB, uidA);

               // [NORMAL-FLIP]: the contact normal of the same uid pair changed by more than ~45 deg (cos < 0.7)
               // since the pair was last seen in contact. Symptom of the discrete feature switch of SAT (Box-Box)
               // or the Kodam scenario dispatch (Cylinder-Cylinder) near symmetric/degenerate configurations
               // (nearly parallel logs, nearly face-aligned boxes in a jam); the spring-dashpot sees an impulsive
               // kick. Sphere-sphere contact has a single feature (line of centres) and cannot do this.
               {
                  const auto it = lastContactNormal_.find(pairKey);
                  if (it != lastContactNormal_.end())
                  {
                     const real_t cosAngle = std::max(real_t(-1), std::min(real_t(1), n * it->second));
                     constexpr real_t kNormalFlipCosThreshold = real_t(0.7); // ~45 degrees
                     if (cosAngle < kNormalFlipCosThreshold)
                     {
                        static uint_t normalFlipLogCount = uint_t(0);
                        if (normalFlipLogCount < uint_t(200))
                        {
                           const auto fs1 = ac.getBaseShape(acd_.getIdx1())
                                              ? ac.getBaseShape(acd_.getIdx1())->getShapeType() : uint_t(999);
                           const auto fs2 = ac.getBaseShape(acd_.getIdx2())
                                              ? ac.getBaseShape(acd_.getIdx2())->getShapeType() : uint_t(999);
                           WALBERLA_LOG_WARNING("[NORMAL-FLIP] step=" << mainStep_ << " subCycle=" << subCycle_
                              << " uidA=" << uidA << " uidB=" << uidB
                              << " shapeTypeA=" << fs1 << " shapeTypeB=" << fs2
                              << " prevNormal=" << it->second << " newNormal=" << n
                              << " angleChangeDeg=" << (std::acos(cosAngle) * real_t(180) / math::pi)
                              << " overlap=" << overlapMag << " (" << (overlapMag / smallerRad * real_t(100))
                              << "% of smallerInteractionRadius=" << smallerRad << ")"
                              << " posA=" << ac.getPosition(acd_.getIdx1()) << " posB=" << ac.getPosition(acd_.getIdx2())
                              << " velA=" << ac.getLinearVelocity(acd_.getIdx1())
                              << " velB=" << ac.getLinearVelocity(acd_.getIdx2()));
                           ++normalFlipLogCount;
                        }
                     }
                  }
                  lastContactNormal_[pairKey] = n;
               }

               // [PROLONGED-CONTACT]: number of sub-cycles this pair has been resolved in the current episode (an
               // episode survives as long as the pair is seen again in the same or the next main step) and the peak
               // overlap. Reported at 30 sub-cycles and then at every multiple of 50; shows whether overlap crept up
               // gradually or appeared within one step (pointing to a detection/sync gap).
               {
                  auto& episode = contactEpisodes_[pairKey];
                  if (mainStep_ > episode.lastSeenMainStep + uint_t(1))
                  {
                     episode.episodeLength    = uint_t(0);
                     episode.episodeStartStep = mainStep_;
                     episode.peakOverlap      = real_t(0);
                  }
                  episode.lastSeenMainStep = mainStep_;
                  ++episode.episodeLength;
                  episode.peakOverlap = std::max(episode.peakOverlap, overlapMag);

                  constexpr uint_t kProlongedContactFirstReport = uint_t(30);
                  constexpr uint_t kProlongedContactReportEvery = uint_t(50);
                  if (episode.episodeLength == kProlongedContactFirstReport ||
                      (episode.episodeLength > kProlongedContactFirstReport &&
                       episode.episodeLength % kProlongedContactReportEvery == uint_t(0)))
                  {
                     static uint_t prolongedContactLogCount = uint_t(0);
                     if (prolongedContactLogCount < uint_t(200))
                     {
                        WALBERLA_LOG_WARNING("[PROLONGED-CONTACT] step=" << mainStep_ << " subCycle=" << subCycle_
                           << " uidA=" << uidA << " uidB=" << uidB
                           << " episodeLength=" << episode.episodeLength << " subCycles (started step "
                           << episode.episodeStartStep << ")"
                           << " currentOverlap=" << overlapMag << " (" << (overlapMag / smallerRad * real_t(100))
                           << "% of smallerInteractionRadius=" << smallerRad << ")"
                           << " peakOverlapThisEpisode=" << episode.peakOverlap
                           << " (" << (episode.peakOverlap / smallerRad * real_t(100)) << "%)"
                           << " posA=" << ac.getPosition(acd_.getIdx1()) << " posB=" << ac.getPosition(acd_.getIdx2()));
                        ++prolongedContactLogCount;
                     }
                  }
               }

               // [SIMULTANEOUS-CONTACTS]: per-sub-cycle contact count per body, reported and reset by
               // reportAndResetSimultaneousContacts() after the sweep.
               ++contactCountThisSubCycle_[uidA];
               ++contactCountThisSubCycle_[uidB];
            }

            // DIAGNOSTIC artificial cap (ModelParameters.rodCapOverlap): limit the penetration depth fed into the force
            // law for contacts involving a rod (type 4)
            if (g_rodCaps.overlap > real_t(0) && (ac.getType(acd_.getIdx1()) == uint_t(4) || ac.getType(acd_.getIdx2()) == uint_t(4)))
            {
               const real_t pen = -acd_.getPenetrationDepth();
               g_rodCapStats.maxOverlap = std::max(g_rodCapStats.maxOverlap, pen);
               if (pen > g_rodCaps.overlap) { logRodCapHit("overlap", ac.getUid(acd_.getIdx1()), mainStep_, subCycle_, pen, g_rodCaps.overlap); acd_.getPenetrationDepth() = -g_rodCaps.overlap; ++g_rodCapStats.overlap; }
            }
            auto meff = real_c(1) / (ac.getInvMass(idx1) + ac.getInvMass(idx2));
            // useRotationalEffectiveMass (any shape pair, see ModelParameters, SimulationParams.h): use the
            // generalized effective mass at the contact point reported by acd_, so LinearSpringDashpot's
            // stiffness/damping calibration is correct for contacts off the line of centres (identical to the plain
            // formula for sphere-sphere).
            if (acd_.useRotationalEffectiveMass)
            {
               meff = computeRotationalEffectiveMass(ac, idx1, idx2, acd_.getContactPoint(), acd_.getContactNormal());
            }
            // Set stiffness for the actual type pair involved in this contact
            const uint_t t1 = std::min(ac.getType(acd_.getIdx1()), ac.getType(acd_.getIdx2()));
            const uint_t t2 = std::max(ac.getType(acd_.getIdx1()), ac.getType(acd_.getIdx2()));
            collisionResponse_.setStiffnessAndDamping(t1, t2, restitutionCoefficient_, collisionTime_, kappa_, meff);

            // Captured unconditionally: [STRONG-CONTACT-FORCE] below needs the delta too (cheap Vec3 copies).
            const mesa_pd::Vec3 forceBeforeA  = ac.getForceRef(acd_.getIdx1());
            const mesa_pd::Vec3 forceBeforeB  = ac.getForceRef(acd_.getIdx2());
            const mesa_pd::Vec3 torqueBeforeA = ac.getTorqueRef(acd_.getIdx1());
            const mesa_pd::Vec3 torqueBeforeB = ac.getTorqueRef(acd_.getIdx2());

            collisionResponse_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(), acd_.getContactNormal(),
                               acd_.getPenetrationDepth(), particleTimestepSize_);

            // [STRONG-CONTACT-FORCE]: always-on safety-net diagnostic, independent of logFloaterHistory. When a
            // single contact's force delta exceeds STRONG_FORCE_THRESHOLD (0.1, well above a typical gravity-scale
            // weight, deliberately lenient), dumps this contact's force/torque, each particle's accumulated
            // force/torque this sub-cycle and its hydrodynamic force/torque. Throttled per rank (100 lines).
            {
               const mesa_pd::Vec3 forceDeltaA  = ac.getForceRef(acd_.getIdx1())  - forceBeforeA;
               const mesa_pd::Vec3 forceDeltaB  = ac.getForceRef(acd_.getIdx2())  - forceBeforeB;
               const mesa_pd::Vec3 torqueDeltaA = ac.getTorqueRef(acd_.getIdx1()) - torqueBeforeA;
               const mesa_pd::Vec3 torqueDeltaB = ac.getTorqueRef(acd_.getIdx2()) - torqueBeforeB;
               constexpr real_t STRONG_FORCE_THRESHOLD = real_t(0.1);
               if (forceDeltaA.length() > STRONG_FORCE_THRESHOLD || forceDeltaB.length() > STRONG_FORCE_THRESHOLD)
               {
                  static uint_t strongForceLogCount = uint_t(0);
                  if (strongForceLogCount < uint_t(100))
                  {
                     const auto s1 = ac.getBaseShape(acd_.getIdx1())
                                        ? ac.getBaseShape(acd_.getIdx1())->getShapeType() : uint_t(999);
                     const auto s2 = ac.getBaseShape(acd_.getIdx2())
                                        ? ac.getBaseShape(acd_.getIdx2())->getShapeType() : uint_t(999);
                     WALBERLA_LOG_WARNING("[STRONG-CONTACT-FORCE] step=" << mainStep_ << " subCycle=" << subCycle_
                        << "  A: uid=" << ac.getUid(acd_.getIdx1()) << " shapeType=" << s1
                        << " pos=" << ac.getPosition(acd_.getIdx1()) << " vel=" << ac.getLinearVelocity(acd_.getIdx1())
                        << " thisContactForce=" << forceDeltaA << " thisContactTorque=" << torqueDeltaA
                        << " totalForceThisSubCycle=" << ac.getForceRef(acd_.getIdx1())
                        << " totalTorqueThisSubCycle=" << ac.getTorqueRef(acd_.getIdx1())
                        << " hydroForce=" << ac.getHydrodynamicForce(acd_.getIdx1())
                        << " hydroTorque=" << ac.getHydrodynamicTorque(acd_.getIdx1())
                        << "  |  B: uid=" << ac.getUid(acd_.getIdx2()) << " shapeType=" << s2
                        << " pos=" << ac.getPosition(acd_.getIdx2()) << " vel=" << ac.getLinearVelocity(acd_.getIdx2())
                        << " thisContactForce=" << forceDeltaB << " thisContactTorque=" << torqueDeltaB
                        << " totalForceThisSubCycle=" << ac.getForceRef(acd_.getIdx2())
                        << " totalTorqueThisSubCycle=" << ac.getTorqueRef(acd_.getIdx2())
                        << " hydroForce=" << ac.getHydrodynamicForce(acd_.getIdx2())
                        << " hydroTorque=" << ac.getHydrodynamicTorque(acd_.getIdx2())
                        << "  |  penetration=" << acd_.getPenetrationDepth() << " normal=" << acd_.getContactNormal()
                        << " contactPoint=" << acd_.getContactPoint());
                     ++strongForceLogCount;
                  }
               }
            }

            // [LAUNCH-TRACE] (env ANTIDUNES_LAUNCH_TRACE_UID=<uid>, ANTIDUNES_LAUNCH_TRACE_FROM=<step>,
            // ANTIDUNES_LAUNCH_TRACE_TO=<step>): every contact of one floater in a step window, with both bodies' full
            // state (position, rotation matrix, linear and angular velocity, shape type), the contact geometry and
            // this contact's force/torque on body A, to locate a contact that injects energy (e.g. floaters launched
            // from a pile). Unthrottled inside the window; inert unless the variable is set.
            {
               static const uint64_t traceUid  = std::getenv("ANTIDUNES_LAUNCH_TRACE_UID")
                                                    ? std::strtoull(std::getenv("ANTIDUNES_LAUNCH_TRACE_UID"), nullptr, 10) : 0;
               static const uint64_t traceFrom = std::getenv("ANTIDUNES_LAUNCH_TRACE_FROM")
                                                    ? std::strtoull(std::getenv("ANTIDUNES_LAUNCH_TRACE_FROM"), nullptr, 10) : 0;
               static const uint64_t traceTo   = std::getenv("ANTIDUNES_LAUNCH_TRACE_TO")
                                                    ? std::strtoull(std::getenv("ANTIDUNES_LAUNCH_TRACE_TO"), nullptr, 10) : 0;
               // Speed trigger (ANTIDUNES_LAUNCH_SPEED=<LU/step>): trace every contact of a floater faster than
               // this, at most 4000 lines per rank, for runs that do not reproduce a specific launch exactly.
               static const double traceSpeed = std::getenv("ANTIDUNES_LAUNCH_SPEED")
                                                   ? std::strtod(std::getenv("ANTIDUNES_LAUNCH_SPEED"), nullptr) : 0.0;
               static uint_t       speedLines = 0;
               const uint64_t uA = uint64_t(ac.getUid(acd_.getIdx1())), uB = uint64_t(ac.getUid(acd_.getIdx2()));
               const bool byUid = traceUid != 0 && (uA == traceUid || uB == traceUid) && uint64_t(mainStep_) >= traceFrom &&
                                  uint64_t(mainStep_) <= traceTo;
               const bool bySpeed = traceSpeed > 0.0 && speedLines < uint_t(4000) &&
                                    (double(ac.getLinearVelocity(acd_.getIdx1()).length()) > traceSpeed ||
                                     double(ac.getLinearVelocity(acd_.getIdx2()).length()) > traceSpeed);
               if (bySpeed) ++speedLines;
               if (byUid || bySpeed)
               {
                  auto st = [&](size_t i) {
                     std::ostringstream o;
                     o.precision(9);
                     const auto  sh = ac.getBaseShape(i);
                     const auto  R  = ac.getRotation(i).getMatrix();
                     o << "uid=" << ac.getUid(i) << " shape=" << (sh ? int(sh->getShapeType()) : -1)
                       << (isSet(ac.getFlags(i), mesa_pd::data::particle_flags::GHOST) ? " ghost" : " local")
                       << " pos=" << ac.getPosition(i) << " v=" << ac.getLinearVelocity(i) << " w=" << ac.getAngularVelocity(i)
                       << " R=[" << R(0, 0) << "," << R(0, 1) << "," << R(0, 2) << ";" << R(1, 0) << "," << R(1, 1) << ","
                       << R(1, 2) << ";" << R(2, 0) << "," << R(2, 1) << "," << R(2, 2) << "]";
                     return o.str();
                  };
                  std::ostringstream o;
                  o.precision(9);
                  o << "[LAUNCH-TRACE] step=" << mainStep_ << " sub=" << subCycle_ << " depth=" << acd_.getPenetrationDepth()
                    << " n=" << acd_.getContactNormal() << " cp=" << acd_.getContactPoint()
                    << " F_A=" << (ac.getForceRef(acd_.getIdx1()) - forceBeforeA)
                    << " T_A=" << (ac.getTorqueRef(acd_.getIdx1()) - torqueBeforeA) << " | A: " << st(acd_.getIdx1())
                    << " | B: " << st(acd_.getIdx2());
                  WALBERLA_LOG_INFO(o.str());
               }
            }

            const bool doLog = historyLogger_ && historyLogger_->enabled();

            if (doLog)
            {
               // Recorded only by the Box-Box / Cylinder-Cylinder dispatch and the pillar/containment overloads
               // (see recordDetectionPath); other pairs get the fixed label below.
               const std::string* trackedPath = acd_.lookupDetectionPath(ac.getUid(acd_.getIdx1()), ac.getUid(acd_.getIdx2()));
               const std::string detectionPath = trackedPath ? *trackedPath
                  : std::string("not Box-Box/Cylinder-Cylinder -- exact analytic overload or "
                                 "always-GJK/EPA pair (no useAnalyticalContactDetection choice)");
               historyLogger_->log(mainStep_, subCycle_, "contact", ac, acd_.getIdx1(), acd_.getIdx2(),
                  acd_.getContactPoint(), acd_.getContactNormal(), acd_.getPenetrationDepth(),
                  ac.getForceRef(acd_.getIdx1()) - forceBeforeA, ac.getForceRef(acd_.getIdx2()) - forceBeforeB,
                  detectionPath);
            }

            // Secondary contact (arch/tunnel containment only, Box/Cylinder floaters; see
            // SolidCylinderContactDetection::hasSecondaryContact()). Not routed through collisionResponse_:
            // LinearSpringDashpot's tangential-spring history is keyed by the particle pair only, so a second call
            // would overwrite the primary contact's history. Instead a history-free normal spring-dashpot plus a
            // Coulomb-capped tangential dashpot, with the stiffness/damping/friction just set for this type pair.
            if (acd_.hasSecondaryContact())
            {
               const mesa_pd::Vec3& secNormal = acd_.getSecondaryContactNormal();
               const mesa_pd::Vec3& secPoint  = acd_.getSecondaryContactPoint();
               const real_t         secOverlap = -acd_.getSecondaryPenetrationDepth();

               const real_t stiffnessN = collisionResponse_.getStiffnessN(t1, t2);
               const real_t dampingN   = collisionResponse_.getDampingN(t1, t2);
               const real_t dampingT   = collisionResponse_.getDampingT(t1, t2);
               const real_t muDynamic  = collisionResponse_.getFrictionCoefficientDynamic(t1, t2);

               const mesa_pd::Vec3 relVel = mesa_pd::getVelocityAtWFPoint(acd_.getIdx2(), ac, secPoint)
                                           - mesa_pd::getVelocityAtWFPoint(acd_.getIdx1(), ac, secPoint);
               const real_t        velN   = relVel * secNormal;
               const mesa_pd::Vec3 velT   = relVel - velN * secNormal;

               real_t fN = stiffnessN * secOverlap + dampingN * velN;
               if (fN < real_t(0)) fN = real_t(0); // a spring-dashpot must never pull the bodies together

               mesa_pd::Vec3 fT = dampingT * velT;
               const real_t  fTmag = fT.length();
               const real_t  fTmax = muDynamic * fN;
               if (fTmag > fTmax && fTmag > real_t(1e-12)) fT *= (fTmax / fTmag);

               const mesa_pd::Vec3 secondaryForce = fN * secNormal + fT;
               mesa_pd::addForceAtWFPosAtomic(acd_.getIdx1(), ac,  secondaryForce, secPoint);
               mesa_pd::addForceAtWFPosAtomic(acd_.getIdx2(), ac, -secondaryForce, secPoint);
            }
         }
      }
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain >    domain_;
   mesa_pd::kernel::LinearSpringDashpot                     collisionResponse_;
   SolidCylinderContactDetection                            acd_;
   mesa_pd::kernel::DoubleCast                              double_cast_;
   mesa_pd::mpi::ContactFilter                              contact_filter_;

   real_t restitutionCoefficient_;
   real_t kappa_;
   real_t collisionTime_;
   real_t particleTimestepSize_;

   // [DEM-WATCH] threshold, see operator().
   static constexpr real_t kDemWatchOverlapRatio = real_t(0.06);

   FloaterHistoryLogger* historyLogger_ = nullptr;
   uint_t mainStep_ = uint_t(0);
   uint_t subCycle_ = uint_t(0);

   // See setFreshlyActivatedUIDs().
   std::shared_ptr< std::unordered_set< walberla::id_t > > freshlyActivatedUIDs_;

   // Diagnostic side-channel (rank-local, not synced), like SolidCylinderContactDetection::pairDetectionPath:
   // per uid pair, the local input state (both particles' pose/velocity, ghost flag, owner rank) this rank fed
   // into contact_filter_ this cycle. Printed by contactHistoryDuplicateDiagnosticHook, so when two ranks
   // disagree about a contact ([ContactHistory DUPLICATE], e.g. a ContactFilter double-accept on a pair
   // straddling a rank boundary) both ranks' views can be compared, e.g. to detect stale ghost positions.
   // Raw values are stored; format() builds the text on lookup.
   struct ContactInputSide
   {
      walberla::id_t                uid;
      int                           owner;
      bool                          ghost;
      mesa_pd::Vec3                 pos;
      math::Quaternion< real_t >    quat;
      mesa_pd::Vec3                 vel, angVel;
   };
   struct ContactInputState
   {
      int              rank;
      ContactInputSide s1, s2;
      mesa_pd::Vec3    point, normal;
      real_t           penetration;
      std::string format() const
      {
         std::ostringstream oss;
         oss << "rank=" << rank;
         int k = 1;
         for (const auto* s : { &s1, &s2 })
         {
            oss << " | idx" << k++ << ": uid=" << s->uid << " owner=" << s->owner << " ghost=" << (s->ghost ? "yes" : "no")
                << " pos=" << s->pos << " rot(quat)=" << s->quat << " vel=" << s->vel << " angVel=" << s->angVel;
         }
         oss << " | computed contactPoint=" << point << " normal=" << normal << " penetration=" << penetration;
         return oss.str();
      }
   };

   void recordContactInputState(walberla::id_t uidA, walberla::id_t uidB, const ContactInputState& state)
   {
      const auto key = uidA < uidB ? std::make_pair(uidA, uidB) : std::make_pair(uidB, uidA);
      pairContactInputState_[key] = state;
   }

   const ContactInputState* lookupContactInputState(walberla::id_t uidA, walberla::id_t uidB) const
   {
      const auto key = uidA < uidB ? std::make_pair(uidA, uidB) : std::make_pair(uidB, uidA);
      const auto it = pairContactInputState_.find(key);
      return it == pairContactInputState_.end() ? nullptr : &it->second;
   }

   std::map< std::pair< walberla::id_t, walberla::id_t >, ContactInputState > pairContactInputState_;

   // See setLogContactAnomalies().
   bool   logContactAnomalies_     = false;
   uint_t contactAnomalyStartStep_ = uint_t(0);
   uint_t contactAnomalyEndStep_   = std::numeric_limits< uint_t >::max();

   // [NORMAL-FLIP]'s per-pair last contact normal (rank-local, not synced).
   std::map< std::pair< walberla::id_t, walberla::id_t >, mesa_pd::Vec3 > lastContactNormal_;

   // [PROLONGED-CONTACT]'s per-pair episode state.
   struct ContactEpisodeState
   {
      uint_t lastSeenMainStep = uint_t(0);
      uint_t episodeStartStep = uint_t(0);
      uint_t episodeLength    = uint_t(0);
      real_t peakOverlap      = real_t(0);
   };
   std::map< std::pair< walberla::id_t, walberla::id_t >, ContactEpisodeState > contactEpisodes_;

   // [SIMULTANEOUS-CONTACTS]' per-sub-cycle tally: filled in operator(), consumed and cleared by
   // reportAndResetSimultaneousContacts().
   std::unordered_map< walberla::id_t, uint_t > contactCountThisSubCycle_;
};



class AddExternalForces {
public:
   AddExternalForces(const std::shared_ptr< mesa_pd::domain::BlockForestDomain >& domain,
      const std::shared_ptr< Vector3< real_t > >& acceleration,
      real_t particleDensityRatio,
      bool useMemOnlyBuoyancy = true,
      bool driveFloatersByMass = false)
      : domain_(domain)
      , acceleration_(acceleration)
      , particleDensityRatio_(particleDensityRatio)
      , useMemOnlyBuoyancy_(useMemOnlyBuoyancy)
      , driveFloatersByMass_(driveFloatersByMass)
      , submergedVolumes_(nullptr)
      , buoyancyCenterOffset_(nullptr)
   {}

   void setSubmergedVolumes(const std::unordered_map< walberla::id_t, real_t >* sv) {
      submergedVolumes_ = sv;
   }

   void setBuoyancyCenterOffset(const std::unordered_map< walberla::id_t, Vector3< real_t > >* off) {
      buoyancyCenterOffset_ = off;
   }

   void setComOffsetZ(const std::unordered_map< walberla::id_t, real_t >* com) {
      comOffsetZ_ = com;
   }

   // Step counter (owned by ParticleKernels, advances once per step), for the once-per-eval-step force log.
   void setStepInfo(const uint_t* stepCounter) { stepCounter_ = stepCounter; }

   // Set once per sub-cycle by doDEMSubCycles; true only on the first sub-cycle of an eval step, so every
   // floater logs its force breakdown once per eval step.
   void setLogForces(bool logForces) { logForces_ = logForces; }

   // Sub-cycle index within the current LBM step (0..numSubCycles-1), for the [MEM-WATCH] log. hydro_force is
   // fixed for the whole outer step, so a value already anomalous at sub-cycle 0 points at MEM, one that grows
   // across sub-cycles at DEM contact resolution.
   void setSubCycle(uint_t subCycle) { subCycle_ = subCycle; }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac){
      // Buoyancy acts only in Z because the hydrostatic MEM contribution is only in Z.
      // Type 1 = floating sphere, type 3 = floating box, type 4 = floating cylinder; density from invMass x volume.
      const uint_t ptype = ac.getType(idx);
      const bool isFloater = (ptype == uint_t(1) || ptype == uint_t(3) || ptype == uint_t(4));
      const real_t rho = isFloater
         ? real_t(1) / (ac.getInvMass(idx) * ac.getVolume(idx))
         : particleDensityRatio_;

      const real_t vFull = ac.getVolume(idx);
      // For floaters use the precomputed submerged volume so the gas-phase portion does not contribute a
      // spurious buoyancy force. Non-floaters are always fully submerged.
      real_t vSub = vFull;
      if (isFloater && submergedVolumes_) {
         auto it = submergedVolumes_->find(ac.getUid(idx));
         if (it != submergedVolumes_->end()) vSub = it->second;
      }

      // X/Y: full-volume driving force (pressure-gradient body force acts on the whole particle).
      // Z:   useMemOnlyBuoyancy_ (default) -> gravity only (rho * vFull); the MEM force alone supplies buoyant
      //      support (Bogner & Rude 2013, no separate geometric channel).
      //      otherwise -> gravity minus the analytic buoyancy term (vSub x g, fluid density = 1).
      // driveFloatersByMass_ -> floaters are driven with their own mass in X/Y (rho * vFull = m_p; slope
      //      reading of the body force) instead of the fluid-density full-volume term.
      const real_t driveVol = (isFloater && driveFloatersByMass_) ? rho * vFull : vFull;
      Vector3< real_t > bodyForce(
         driveVol * (*acceleration_)[0],
         driveVol * (*acceleration_)[1],
         (useMemOnlyBuoyancy_ ? (rho * vFull) : (rho * vFull - vSub)) * (*acceleration_)[2]);

      auto hydro_force  = ac.getHydrodynamicForce(idx);
      auto hydro_torque = ac.getHydrodynamicTorque(idx);
      const Vector3< real_t > hydroForceRaw = hydro_force; // pre-scaling value, kept for the verification log

      // Submersion-scaled hydrodynamic force (V_sub/V_full). The MEM (SimpleBB) treats near-boundary cells as full
      // liquid, so the emerged part of a floater can see a phantom upward force (refilled vacated cells, raised
      // surface under-counted by the flat hydrostatic reference rho_hyd(z) = 1 + 3 a_z (z - h_liquid)).
      real_t submFrac = real_t(1);
      if (isFloater) {
         submFrac = (vFull > real_t(0)) ? std::clamp(vSub / vFull, real_t(0), real_t(1)) : real_t(0);
         // Deliberately OFF: the phantom lift is corrected at its source by the per-link fill weighting phi(x_b)
         // in SimpleBB (Bogner 2017, Eq. 5.9). submFrac is still computed for the [FloaterForces] log.
         // hydro_force  *= submFrac;
         // hydro_torque *= submFrac;
      }

      if (!std::isfinite(hydro_force[0]) || !std::isfinite(hydro_force[1]) || !std::isfinite(hydro_force[2]))
         WALBERLA_LOG_WARNING("[AddExternalForces] NaN/Inf hydrodynamic force on uid=" << ac.getUid(idx)
            << " pos=" << ac.getPosition(idx));
      if (!std::isfinite(hydro_torque[0]) || !std::isfinite(hydro_torque[1]) || !std::isfinite(hydro_torque[2]))
         WALBERLA_LOG_WARNING("[AddExternalForces] NaN/Inf hydrodynamic torque on uid=" << ac.getUid(idx)
            << " pos=" << ac.getPosition(idx));

      // DIAGNOSTIC artificial caps (ModelParameters.rodCapHydroForce/Torque), rods (type 4) only
      if ((g_rodCaps.hydroForce > real_t(0) || g_rodCaps.hydroTorque > real_t(0)) && ac.getType(idx) == uint_t(4))
      {
         const uint_t stepNow = stepCounter_ ? *stepCounter_ : uint_t(0);
         if (!(std::isfinite(hydro_force[0]) && std::isfinite(hydro_force[1]) && std::isfinite(hydro_force[2]))) { hydro_force = Vector3< real_t >(real_t(0)); ++g_rodCapStats.nonFinite; }
         if (!(std::isfinite(hydro_torque[0]) && std::isfinite(hydro_torque[1]) && std::isfinite(hydro_torque[2]))) { hydro_torque = Vector3< real_t >(real_t(0)); ++g_rodCapStats.nonFinite; }
         const real_t fm = hydro_force.length(), tm = hydro_torque.length();
         g_rodCapStats.maxHydroForce = std::max(g_rodCapStats.maxHydroForce, fm); g_rodCapStats.maxHydroTorque = std::max(g_rodCapStats.maxHydroTorque, tm);
         if (g_rodCaps.hydroForce > real_t(0) && fm > g_rodCaps.hydroForce) { logRodCapHit("hydroForce", ac.getUid(idx), stepNow, subCycle_, fm, g_rodCaps.hydroForce); hydro_force *= g_rodCaps.hydroForce / fm; ++g_rodCapStats.hydroForce; }
         if (g_rodCaps.hydroTorque > real_t(0) && tm > g_rodCaps.hydroTorque) { logRodCapHit("hydroTorque", ac.getUid(idx), stepNow, subCycle_, tm, g_rodCaps.hydroTorque); hydro_torque *= g_rodCaps.hydroTorque / tm; ++g_rodCapStats.hydroTorque; }
      }

      // [MEM-WATCH]: |hydro_force| > kMemWatchThreshold. hydro_force is computed once per LBM step
      // (reduceAndAvgHydroForces, before doDEMSubCycles) and fixed across sub-cycles, so a large value at
      // subCycle_==0 points at the fluid coupling rather than DEM. Checked every sub-cycle to catch the exact step.
      // Gated by kEnableMemDemWatch.
      if (kEnableMemDemWatch)
      {
         const real_t hydroMag = hydro_force.length();
         if (hydroMag > kMemWatchThreshold)
         {
            // Throttled to 100 hits per rank (as [STRONG-CONTACT-FORCE]): in dense piles steady-state hydro forces
            // can exceed the threshold every sub-cycle.
            static uint_t memWatchLogCount = uint_t(0);
            if (memWatchLogCount < uint_t(100))
            {
               WALBERLA_LOG_WARNING("[MEM-WATCH] large hydrodynamic force uid=" << ac.getUid(idx)
                  << " step=" << (stepCounter_ ? *stepCounter_ : uint_t(0)) << " subCycle=" << subCycle_
                  << " |F_hydro|=" << hydroMag << " F_hydro=" << hydro_force << " F_hydro_torque=" << hydro_torque
                  << " pos=" << ac.getPosition(idx) << " vel=" << ac.getLinearVelocity(idx)
                  << " bodyForce(grav/driving)=" << bodyForce
                  << " contactForceThisSubCycle=" << ac.getForce(idx));
               ++memWatchLogCount;
            }
         }
      }

      // [DEM-WATCH] is in ComputeContacts::operator() (per-contact overlap ratio, traceable to one contact).

      // Per-floater force breakdown, once per eval step (doDEMSubCycles sets logForces_ on sub-cycle 0).
      // ac.getForce(idx) = contact + lubrication accumulated this sub-cycle (ComputeContacts runs first);
      // hydro_force and bodyForce are the MEM and buoyancy/gravity/driving parts about to be added. Compares
      // F_buoy+grav (model), F_hydro (flow) and F_contact (pile / pillar / wall) per floater.
      if (logForces_ && isFloater)
      {
         const auto fContact = ac.getForce(idx);
         // Y-tilt angle, same convention as [BoxRightingMoment]'s tiltDeg below (Bogner & Rude 2013, Fig. 5.6
         // validation): for a pure rotation about Y, R(0,2)=sin(theta), R(0,0)=cos(theta). Box/cylinder floaters
         // only (spheres carry no physical tilt).
         std::string tiltStr;
         if (ptype == uint_t(3) || ptype == uint_t(4))
         {
            const auto& R = ac.getRotation(idx).getMatrix();
            const real_t tiltYDeg = std::atan2(R(0, 2), R(0, 0)) * real_t(180) / math::pi;
            tiltStr = "  tiltY=" + std::to_string(tiltYDeg) + "deg";
         }
         WALBERLA_LOG_INFO((ptype == uint_t(3) ? "[BoxForces]"
                            : ptype == uint_t(4) ? "[CylForces]" : "[FloaterForces]")
            << " uid=" << ac.getUid(idx) << " step=" << (stepCounter_ ? *stepCounter_ : uint_t(0))
            << "  subm=" << (vFull > real_t(0) ? real_t(100) * vSub / vFull : real_t(0)) << "%"
            << "  F_buoy+grav_z=" << bodyForce[2] << " (>0 up)"
            << "  F_drive_xy=(" << bodyForce[0] << "," << bodyForce[1] << ")"
            << "  F_hydro=(" << hydro_force[0] << "," << hydro_force[1] << "," << hydro_force[2] << ")"
            << "  [submScale=" << submFrac << " applied to F_hydro_raw=("
            << hydroForceRaw[0] << "," << hydroForceRaw[1] << "," << hydroForceRaw[2] << ")]"
            << "  F_contact+lub=(" << fContact[0] << "," << fContact[1] << "," << fContact[2] << ")"
            << "  vel=" << ac.getLinearVelocity(idx)
            << "  pos=" << ac.getPosition(idx)
            << tiltStr);
      }

      mesa_pd::addForceAtomic(idx, ac, hydro_force + bodyForce);

      // Buoyancy righting / capsizing moment (Bogner & Rude 2013, Fig. 5.4): the buoyancy force (vertical,
      // +vSub g) acts at the centre of buoyancy B (centroid of the submerged volume), not at the COM G. When the
      // body tilts, B shifts off G and the buoyancy produces the moment (B - G) x F_buoy, applied here as an
      // explicit torque (net force unchanged). It self-rights when the metacentre is above G and tips further
      // when below. Spheres store no offset (B on the vertical through G, zero torque).
      Vector3< real_t > totalTorque = hydro_torque;
      if (isFloater && buoyancyCenterOffset_)
      {
         auto itc = buoyancyCenterOffset_->find(ac.getUid(idx));
         if (itc != buoyancyCenterOffset_->end())
         {
            const Vector3< real_t >& lever = itc->second;                                   // B − G (world)
            const Vector3< real_t >  Fbuoy(real_t(0), real_t(0), -vSub * (*acceleration_)[2]); // +vSub·g, up
            totalTorque += Vector3< real_t >(lever[1] * Fbuoy[2] - lever[2] * Fbuoy[1],
                                             lever[2] * Fbuoy[0] - lever[0] * Fbuoy[2],
                                             lever[0] * Fbuoy[1] - lever[1] * Fbuoy[0]);
         }
      }

      // Lowered centre of mass (pendulum righting torque): comOffsetZ < 0 places the COM below the geometric
      // centroid (the rotation centre); when the box tilts, the weight produces the restoring moment
      // (COM - centroid) x F_grav, raising the metacentric height (counters GM < 0 rolling of a uniform log/cube).
      // comOffsetZ = 0 reproduces the uniform box (no torque).
      if (isFloater && comOffsetZ_)
      {
         auto ito = comOffsetZ_->find(ac.getUid(idx)); // map holds only nonzero offsets (guarded at build)
         if (ito != comOffsetZ_->end())
         {
            const mesa_pd::Mat3& R = ac.getRotation(idx).getMatrix();                     // body → world
            // COM − centroid (world)
            const Vector3< real_t > leverG = R * Vector3< real_t >(real_t(0), real_t(0), ito->second);
            // weight (down: a_z<0)
            const Vector3< real_t > Fgrav(real_t(0), real_t(0), rho * vFull * (*acceleration_)[2]);
            totalTorque += Vector3< real_t >(leverG[1] * Fgrav[2] - leverG[2] * Fgrav[1],
                                             leverG[2] * Fgrav[0] - leverG[0] * Fgrav[2],
                                             leverG[0] * Fgrav[1] - leverG[1] * Fgrav[0]);
         }
      }
      // Torque breakdown, once per eval step: hydro_torque is the raw MEM momentum-exchange torque (independent
      // of the buoyancy volume/centroid method) vs. the lever torques added above (buoyancy-centre offset and COM
      // offset; buoyLeverTorque_y is their sum), to tell a hydrodynamic roll torque from the lever torques.
      if (logForces_ && isFloater && (ptype == uint_t(3) || ptype == uint_t(4)))
      {
         WALBERLA_LOG_INFO((ptype == uint_t(3) ? "[BoxTorque]" : "[CylTorque]")
            << " uid=" << ac.getUid(idx) << " step=" << (stepCounter_ ? *stepCounter_ : uint_t(0))
            << "  hydroTorque=(" << hydro_torque[0] << "," << hydro_torque[1] << "," << hydro_torque[2] << ")"
            << "  buoyLeverTorque_y=" << (totalTorque[1] - hydro_torque[1])
            << "  totalTorque=(" << totalTorque[0] << "," << totalTorque[1] << "," << totalTorque[2] << ")");
      }
      mesa_pd::addTorqueAtomic(idx, ac, totalTorque);
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain >   domain_;
   std::shared_ptr< Vector3< real_t > >                    acceleration_;
   real_t                                                   particleDensityRatio_;
   // true  -> (default) vertical body force is gravity only (rho*vFull*a_z); buoyant support comes solely from
   //         the MEM force, as in Bogner & Rude (2013). computeSubmergedVolumes()/computeSubmergedVolumesBox()
   //         are skipped in this mode (doDEMSubCycles): no vSub and no righting torque, consistent with Bogner's
   //         Eq. 18 torque (raw MEM only).
   // false -> gravity minus the analytic vSub*g buoyancy term, plus the vSub-derived righting torque (Antidunes'
   //         extended split, not used for the thesis).
   bool                                                     useMemOnlyBuoyancy_;
   // true → floater X/Y drive is m_p * a_xy instead of vFull * a_xy (see SimulationParams driveFloatersByMass).
   bool                                                     driveFloatersByMass_;
   const std::unordered_map< walberla::id_t, real_t >*     submergedVolumes_;
   // Option A vertical MEM correction; null = disabled
   const std::unordered_map< walberla::id_t, Vector3< real_t > >* buoyancyCenterOffset_;
   const std::unordered_map< walberla::id_t, real_t >*     comOffsetZ_ = nullptr; // per box/cylinder COM offset

   // For the per-floater force-breakdown log.
   const uint_t* stepCounter_ = nullptr; // step number, for the log only
   bool          logForces_   = false;   // set true by doDEMSubCycles on the first sub-cycle of an eval step
   uint_t        subCycle_    = uint_t(0); // current DEM sub-cycle index, for DEM/MEM-WATCH (see setSubCycle)

   // [MEM-WATCH] threshold. Typical per-floater forces in stable runs are ~1e-3 to ~1e-2 ([FloaterForces], e.g.
   // F_buoy+grav_z ~ -0.003, F_hydro ~ 0.002); 0.05 is about an order of magnitude above, so it stays silent
   // normally but fires well before a value goes non-finite. Checked every sub-cycle for every local particle.
   static constexpr real_t kMemWatchThreshold = real_t(0.05);
   // [MEM-WATCH] on/off (throttled, see the call site).
   static constexpr bool kEnableMemDemWatch = true;
};


// ─────────────────────────────────────────────────────────────────────────────────────────────────
// Analytical (exact) submerged volume + centre of buoyancy of an oriented box (OBB) clipped by the fitted
// free-surface plane, the closed-form counterpart to the 12³ point-sampling clip (prototype:
// local_only/debug_and_exploratory_scripts/box_overlap_volume.py). The box is split into 6 tetrahedra (Kuhn
// triangulation), each clipped against the half-space, summing volume and volume-weighted centroid.
//
//   h    : box half-extents (body frame)          nHat : fitted surface normal (world), LIQUID→GAS (up)
//   R    : rotation body→world (orthonormal)       xBar : a point on the fitted plane (world)
//   pos  : box centre (world)
// Returns  vSub  = submerged volume, and  cbOff = world-frame offset of the centre of buoyancy from the
// box centre (the lever arm B−G used for the righting torque), same meaning/frame as the sampling clip.

// Unsigned volume of a tetrahedron (a,b,c,e).
inline real_t tetVolume(const Vector3< real_t >& a, const Vector3< real_t >& b,
                        const Vector3< real_t >& c, const Vector3< real_t >& e)
{
   return std::abs(((b - a) % (c - a)) * (e - a)) / real_t(6); // %=cross, *=dot
}

// Accumulate (sign · volume) and (sign · volume · centroid) of tetrahedron (a,b,c,e).
inline void accumTet(const Vector3< real_t >& a, const Vector3< real_t >& b,
                     const Vector3< real_t >& c, const Vector3< real_t >& e,
                     real_t sign, real_t& volAcc, Vector3< real_t >& cvolAcc)
{
   const real_t v = tetVolume(a, b, c, e);
   volAcc  += sign * v;
   cvolAcc += (sign * v) * ((a + b + c + e) / real_t(4));
}

// Clip one tetrahedron against the half-space { x : s(x) = nb·x − d ≤ 0 } (liquid side) and accumulate.
inline void clipTetBelowPlane(const Vector3< real_t > v[4], const Vector3< real_t >& nb, real_t d,
                              real_t& volAcc, Vector3< real_t >& cvolAcc)
{
   real_t s[4];
   int in[4], out[4];
   int nIn = 0, nOut = 0;
   for (int i = 0; i < 4; ++i)
   {
      s[i] = nb * v[i] - d;
      if (s[i] <= real_t(0)) in[nIn++] = i; else out[nOut++] = i;
   }
   // Point on edge v[i]→v[j] where s crosses zero (i inside, j outside): t = s_i/(s_i − s_j).
   auto cut = [&](int i, int j) {
      const real_t t = s[i] / (s[i] - s[j]);
      return v[i] + t * (v[j] - v[i]);
   };

   if (nIn == 0) return;                                             // fully above
   if (nIn == 4) { accumTet(v[0], v[1], v[2], v[3], +1, volAcc, cvolAcc); return; } // fully below

   if (nIn == 1)
   {
      const Vector3< real_t >& A = v[in[0]];
      accumTet(A, cut(in[0], out[0]), cut(in[0], out[1]), cut(in[0], out[2]), +1, volAcc, cvolAcc);
   }
   else if (nIn == 3)
   {
      // Whole tet minus the small corner tet at the single outside vertex.
      const Vector3< real_t >& O = v[out[0]];
      accumTet(v[0], v[1], v[2], v[3], +1, volAcc, cvolAcc);
      accumTet(O, cut(in[0], out[0]), cut(in[1], out[0]), cut(in[2], out[0]), -1, volAcc, cvolAcc);
   }
   else // nIn == 2: the submerged part is a triangular prism → fan into 3 tets.
   {
      const Vector3< real_t >& A = v[in[0]];
      const Vector3< real_t >& B = v[in[1]];
      const Vector3< real_t > eAC = cut(in[0], out[0]);
      const Vector3< real_t > eAD = cut(in[0], out[1]);
      const Vector3< real_t > eBC = cut(in[1], out[0]);
      const Vector3< real_t > eBD = cut(in[1], out[1]);
      // Prism ends (A,eAC,eAD) and (B,eBC,eBD); fan the 3 tets from A and eBD.
      accumTet(A, eAC, eAD, eBD, +1, volAcc, cvolAcc);
      accumTet(A, eAC, eBD, eBC, +1, volAcc, cvolAcc);
      accumTet(A, eBC, eBD, B,   +1, volAcc, cvolAcc);
   }
}

inline void clipBoxAnalytic(const Vector3< real_t >& h, const mesa_pd::Mat3& R,
                            const Vector3< real_t >& nHat, const Vector3< real_t >& xBar,
                            const Vector3< real_t >& pos, real_t& vSub, Vector3< real_t >& cbOff)
{
   // Move to the box body frame: the box is the axis-aligned cuboid [−h,h]³ centred at the origin.
   // world x = pos + R·xb, so the world half-space nHat·(x−xBar) ≤ 0 becomes  nb·xb ≤ d  with
   //   nb = Rᵀ·nHat  and  d = nb · (Rᵀ·(xBar − pos)).
   const mesa_pd::Mat3 Rt = R.getTranspose();
   const Vector3< real_t > nb = Rt * nHat;
   const Vector3< real_t > p0 = Rt * (xBar - pos);
   const real_t d = nb * p0;

   // 8 body-frame corners, addressed by sign bits.
   auto corner = [&](int sx, int sy, int sz) {
      return Vector3< real_t >(real_t(sx) * h[0], real_t(sy) * h[1], real_t(sz) * h[2]);
   };
   // Kuhn 6-tetrahedron triangulation of the cuboid: the 6 monotone paths (axis permutations) from the
   // corner (−,−,−) to (+,+,+) along cube edges. Exact partition, valid for any rectangular box.
   const Vector3< real_t > base = corner(-1, -1, -1);
   const Vector3< real_t > e[3] = { Vector3< real_t >(real_t(2) * h[0], 0, 0),
                                    Vector3< real_t >(0, real_t(2) * h[1], 0),
                                    Vector3< real_t >(0, 0, real_t(2) * h[2]) };
   const int perms[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };

   real_t            volAcc(0);
   Vector3< real_t > cvolAcc(real_t(0));
   for (const auto& pm : perms)
   {
      Vector3< real_t > tet[4];
      tet[0] = base;
      tet[1] = tet[0] + e[pm[0]];
      tet[2] = tet[1] + e[pm[1]];
      tet[3] = tet[2] + e[pm[2]];
      clipTetBelowPlane(tet, nb, d, volAcc, cvolAcc);
   }

   vSub = std::max(real_t(0), volAcc);
   // Body-frame submerged centroid (offset from box centre) → rotate back to world for the lever arm.
   cbOff = (volAcc > real_t(1e-12)) ? R * (cvolAcc / volAcc) : Vector3< real_t >(real_t(0));
}

// ---------------------------------------------------------------------------------------------
// Analytic cylinder∩half-space submerged volume + centroid — exact for any plane orientation,
// closed form (no sampling), same role as clipBoxAnalytic / the spherical-cap formula.
//
//   r      : cylinder radius (body frame)             nHat : fitted surface normal (world), LIQUID→GAS
//   halfLen: cylinder half-length (body-frame Z axis)  xBar : a point on the fitted plane (world)
//   R      : rotation body→world (orthonormal)         pos  : cylinder centre (world)
// Returns vSub = submerged volume, and cbOff = world-frame offset of the centre of buoyancy from the
// cylinder centre (same meaning/frame as clipBoxAnalytic / the sampling clip).
//
// Method: the circular cross-section is invariant under rotation about the body Z (axis), so the plane normal's
// horizontal component is first aligned azimuthally onto body X. In that frame the plane is rho*x + nz*z = d
// (no y dependence), so the submerged column length at a given x is a CLAMPED AFFINE function of x alone.
// Splitting the disk into the (at most) three regions where that clamp is saturated-low / linear /
// saturated-high, each region's contribution to volume and first moments (∫x dV, ∫z dV; ∫y dV = 0) is a
// closed-form combination of the circular antiderivatives ∫w, ∫x·w, ∫x²·w with w(x) = 2·sqrt(r²−x²). No
// numerical integration or sampling; unlike the 12³-point Cartesian sampling grid (not azimuthally symmetric at
// finite resolution), it produces no spurious roll torque on axisymmetric cylinders.
inline void clipCylinderAnalytic(real_t r, real_t halfLen, const mesa_pd::Mat3& R,
                                 const Vector3< real_t >& nHat, const Vector3< real_t >& xBar,
                                 const Vector3< real_t >& pos, real_t& vSub, Vector3< real_t >& cbOff)
{
   vSub  = real_t(0);
   cbOff = Vector3< real_t >(real_t(0));
   if (r < real_t(1e-9) || halfLen < real_t(1e-9)) return;

   const mesa_pd::Mat3 Rt = R.getTranspose();
   const Vector3< real_t > nb = Rt * nHat;           // plane normal, body frame
   const Vector3< real_t > p0 = Rt * (xBar - pos);   // point on plane, body frame
   const real_t d = nb * p0;

   const real_t vFull = math::pi * r * r * real_t(2) * halfLen;
   const real_t rho   = std::sqrt(nb[0] * nb[0] + nb[1] * nb[1]); // horizontal magnitude of the normal
   const real_t nz    = nb[2];
   constexpr real_t eps = real_t(1e-9);

   // Antiderivatives of w(x)=2·sqrt(r²−x²), x·w(x), x²·w(x) over the disk (clamped to the domain [−r,r]).
   auto clampX = [&](real_t x) { return std::clamp(x, -r, r); };
   auto Warea  = [&](real_t x) { x = clampX(x); return x * std::sqrt(std::max(real_t(0), r * r - x * x)) + r * r * std::asin(x / r); };
   auto Wmom   = [&](real_t x) { x = clampX(x); return -(real_t(2) / real_t(3)) * std::pow(std::max(real_t(0), r * r - x * x), real_t(1.5)); };
   auto Wmom2  = [&](real_t x) {
      x = clampX(x);
      const real_t s = std::sqrt(std::max(real_t(0), r * r - x * x));
      return real_t(0.25) * (x * (real_t(2) * x * x - r * r) * s + r * r * r * r * std::asin(x / r));
   };
   auto I0 = [&](real_t a, real_t b) { return Warea(b) - Warea(a); };
   auto I1 = [&](real_t a, real_t b) { return Wmom(b) - Wmom(a); };
   auto I2 = [&](real_t a, real_t b) { return Wmom2(b) - Wmom2(a); };

   real_t V = real_t(0), Mx = real_t(0), Mz = real_t(0); // moments in the azimuthally-aligned frame

   if (rho < eps)
   {
      // Plane (nearly) perpendicular to the axis: a flat cap, uniform over the whole cross-section.
      const real_t A  = d / nz;
      const real_t zc = std::clamp(A, -halfLen, halfLen);
      real_t len, mid;
      if (nz > real_t(0)) { len = zc + halfLen; mid = (zc - halfLen) * real_t(0.5); }
      else                { len = halfLen - zc; mid = (zc + halfLen) * real_t(0.5); }
      V  = len * (math::pi * r * r);
      Mz = mid * V;
   }
   else if (std::abs(nz) < eps)
   {
      // Plane (nearly) parallel to the axis: a vertical cut, the same segment for the whole length.
      const real_t x0 = clampX(d / rho); // submerged: rho·x ≤ d  ⇔  x ≤ x0
      V  = I0(-r, x0) * (real_t(2) * halfLen);
      Mx = I1(-r, x0) * (real_t(2) * halfLen);
   }
   else
   {
      // General tilt: submerged column length is clamp(A − m·x, −halfLen, halfLen)-derived, affine in x.
      const real_t A = d / nz;
      const real_t m = rho / nz;
      const real_t xHz  = clampX((A - halfLen) / m); // where the column height would hit +halfLen
      const real_t xMHz = clampX((A + halfLen) / m); // where it would hit −halfLen
      const real_t xlo  = std::min(xHz, xMHz);
      const real_t xhi  = std::max(xHz, xMHz);

      // Saturated plateaus outside [xlo,xhi] are provably constant (monotonicity of the affine column
      // height); classify each as fully-submerged (2·halfLen) or fully-dry (0) by one point evaluation.
      auto lenAt = [&](real_t x) {
         const real_t zc = std::clamp(A - m * x, -halfLen, halfLen);
         return (nz > real_t(0)) ? (zc + halfLen) : (halfLen - zc);
      };

      if (xlo > -r)
      {
         const real_t len = (lenAt(-r) > halfLen) ? real_t(2) * halfLen : real_t(0);
         V  += len * I0(-r, xlo);
         Mx += len * I1(-r, xlo);
         // Mz contribution is 0: fully-submerged columns are z-centred (mid=0), dry ones carry no weight.
      }
      if (xhi > xlo)
      {
         real_t cEff, mEff; // len(x) = cEff + mEff·x on the ramp (clamp inactive there by construction)
         if (nz > real_t(0)) { cEff = A + halfLen; mEff = -m; }
         else                { cEff = halfLen - A; mEff = m; }
         const real_t i0 = I0(xlo, xhi), i1 = I1(xlo, xhi), i2 = I2(xlo, xhi);
         const real_t Vramp  = cEff * i0 + mEff * i1;
         const real_t Mxramp = cEff * i1 + mEff * i2;
         const real_t L2     = cEff * cEff * i0 + real_t(2) * cEff * mEff * i1 + mEff * mEff * i2; // ∫len²·w dx
         const real_t Mzramp = (nz > real_t(0)) ? (real_t(0.5) * L2 - halfLen * Vramp)
                                                 : (halfLen * Vramp - real_t(0.5) * L2);
         V  += Vramp;
         Mx += Mxramp;
         Mz += Mzramp;
      }
      if (r > xhi)
      {
         const real_t len = (lenAt(r) > halfLen) ? real_t(2) * halfLen : real_t(0);
         V  += len * I0(xhi, r);
         Mx += len * I1(xhi, r);
      }
   }

   const real_t cx = (V > real_t(1e-12)) ? Mx / V : real_t(0);
   const real_t cz = (V > real_t(1e-12)) ? Mz / V : real_t(0);
   const real_t theta = std::atan2(nb[1], nb[0]); // undo the azimuthal alignment (irrelevant when cx=0)
   const Vector3< real_t > bodyOff(cx * std::cos(theta), cx * std::sin(theta), cz);

   vSub  = std::clamp(V, real_t(0), vFull);
   cbOff = (V > real_t(1e-12)) ? R * bodyOff : Vector3< real_t >(real_t(0));
}

// Volume of a spherical cap of height hCap cut from a sphere of radius r (caller clamps
// hCap to [0, 2r]; hCap=0 -> empty, hCap=r -> exact hemisphere, hCap=2r -> full sphere).
// Used by ParticleKernels::computeSubmergedVolumes() for the sphere floater's submerged volume.
inline real_t sphereCapVolume(real_t r, real_t hCap)
{
   return math::pi * hCap * hCap * (real_t(3) * r - hCap) / real_t(3);
}

template< typename ParticleAccessor_T >
class ParticleKernels{
public:
   ParticleKernels(const std::shared_ptr< const SimulationParams >& simParams,
      const std::shared_ptr< SimulationObjects >& simObjs)
      : particleAccessor(simObjs->particleAccessor)
      , particleStorage(simObjs->particleStorage)
      , syncCall(simObjs->syncCall)
      , blockForest_(simObjs->blockForest)
      , fillFieldID_(simObjs->fillFieldID)
      , flagFieldID_(simObjs->flagFieldID)
      , normalFieldID_(simObjs->normalFieldID)
      , pdfFieldID_(simObjs->pdfFieldID)
      , particleFieldID_(simObjs->particleFieldID)
      // Guarded: the fluid-free DropTest (DropTest.cpp) creates no boundary handling. (The ID is currently not
      // read anywhere in this class.)
      , boundaryHandlingID_(simObjs->antidunesBoundaryHandling ? simObjs->antidunesBoundaryHandling->getHandlingID()
                                                               : BlockDataID())
      , hasFluid_(simObjs->antidunesBoundaryHandling != nullptr)
      , diagnosticSweepEnabled_(simParams->model.diagnosticSweepEnabled)
      , diagnosticEvents_(simObjs->diagnosticEvents)
      , vvIntegratorPreForce(simParams->computed.particleTimestepSize)
      , vvIntegratorPostForce(simParams->computed.particleTimestepSize)
      , computeLubrication(simObjs->mesapdDomain, simParams->computed.viscosity,
         simParams->model.useContactFilterOwnerTiebreak)
      , computeContacts(
         simObjs->mesapdDomain,
         simParams->particle.restitutionCoefficient,
         simParams->computed.particleKappa,
         simParams->computed.particleCollisionTime,
         simParams->computed.particleTimestepSize,
         simParams->particle.frictionCoefficient,
         // 6 types: 0=bed/wall, 1=sphere floater, 2=pillar, 3=box floater, 4=cylinder floater, 5=capsule floater
         uint_t(6),
         simParams->model.useAnalyticalContactDetection,
         simParams->model.usePureAnalyticalCylinder,
         simParams->model.useContactFilterOwnerTiebreak,
         simParams->model.useTightLayer1,
         simParams->model.useRotationalEffectiveMass)
      , dormantFloaterUIDs_(simObjs->dormantFloaterUIDs)
      , freshlyActivatedFloaterUIDs_(simObjs->freshlyActivatedFloaterUIDs)
      , addExternalForces(
         simObjs->mesapdDomain,
         simObjs->acceleration,
         simParams->particle.densityRatio,
         simParams->model.useMemOnlyBuoyancy,
         simParams->model.driveFloatersByMass)
      , useLubricationCorrection(simParams->particle.useLubricationCorrection)
      , evalFrequency_(simParams->evaluation.evaluationFrequency)
      , stepCounter_(simObjs->beginTimeStep) // absolute (resume-safe): step numbers in logs match the original run
      , planeFitRadiusFactor_(simParams->model.floaterPlaneFitRadiusFactor)
      , boxPlaneFitBandFactor_(simParams->model.boxPlaneFitBandFactor)
      , surfaceHeight_(simParams->staticReferenceHeight())
      , acceleration_(simObjs->acceleration)
      , memReferenceHeight_(simObjs->memReferenceHeight)
      , useDynamicReferenceHeight_(simParams->model.useDynamicReferenceHeight)
      , useHydrodynamicForceTimeAveraging_(simParams->model.useHydrodynamicForceTimeAveraging)
   {
      // ModelParameters.logFloaterHistory (SimulationParams.h): opens this rank's CSV (no-op when disabled) and
      // wires it into both places a contact/lubrication pair is computed. The directory is created on root in
      // Antidunes.cpp before this constructor runs.
      floaterHistoryLogger_.init(simParams->model.floaterHistoryBaseFolder, simParams->model.logFloaterHistory);
      floaterHistoryLogger_.setMinStep(simParams->model.logFloaterHistoryStartStep);
      floaterHistoryLogger_.setMaxStep(simParams->model.logFloaterHistoryEndStep);
      computeContacts.setHistoryLogger(&floaterHistoryLogger_);
      computeContacts.setCylinderContactSkin(simParams->model.cylinderContactSkin);
      computeContacts.setRefineCylinderContacts(simParams->model.refineCylinderContacts);
      // ModelParameters.kodamBandEdgeAxialCheck: a switch of the closed form itself (free functions in
      // ClosedFormContacts*.h), hence a namespace-level flag like kodamCapSelectionByAxisDistance.
      closed_form_contacts::kodamBandEdgeAxialCheck = simParams->model.kodamBandEdgeAxialCheck;
      WALBERLA_LOG_INFO_ON_ROOT("kodamBandEdgeAxialCheck = " << (simParams->model.kodamBandEdgeAxialCheck ? "true" : "false"));
      if (simParams->model.cylinderContactSkin > real_t(0))
         WALBERLA_LOG_INFO_ON_ROOT("cylinderContactSkin = " << simParams->model.cylinderContactSkin
                                   << " LU (rod-rod contact below a surface gap of " << real_t(2) * simParams->model.cylinderContactSkin
                                   << ", rod-wall/pillar below " << simParams->model.cylinderContactSkin << ")");
      computeContacts.setBoxContactSkin(simParams->model.boxContactSkin);
      if (simParams->model.boxContactSkin > real_t(0))
         WALBERLA_LOG_INFO_ON_ROOT("boxContactSkin = " << simParams->model.boxContactSkin
                                   << " LU (box-box contact below a surface gap of " << real_t(2) * simParams->model.boxContactSkin
                                   << ", box-wall/pillar below " << simParams->model.boxContactSkin << ")");
      computeContacts.setFreshlyActivatedUIDs(freshlyActivatedFloaterUIDs_);
      computeContacts.setShadowAnalytic(simParams->model.shadowAnalyticCylinderLog, simParams->model.shadowAnalyticStartStep, simParams->model.shadowAnalyticEndStep);
      simParams_shadowEnabled_ = simParams->model.shadowAnalyticCylinderLog;
      g_rodCaps.hydroForce = simParams->model.rodCapHydroForce; g_rodCaps.hydroTorque = simParams->model.rodCapHydroTorque; g_rodCaps.velocity = simParams->model.rodCapVelocity;
      g_rodCaps.angularVelocity = simParams->model.rodCapAngularVelocity; g_rodCaps.overlap = simParams->model.rodCapOverlap;
      g_rodPrescribedVelocity = Vector3< real_t >(simParams->model.rodPrescribedVelocityX, simParams->model.rodPrescribedVelocityY, simParams->model.rodPrescribedVelocityZ); g_rodPrescribedVelocityOn = g_rodPrescribedVelocity.length() > real_t(0);
      if (g_rodPrescribedVelocityOn)
         WALBERLA_LOG_WARNING_ON_ROOT("[RodPrescribedVelocity] TEST OPTION ACTIVE: all rods move with the constant velocity " << g_rodPrescribedVelocity << " and zero angular velocity (forces ignored) -- changes the physics, not for production");
      if (g_rodCaps.any())
         WALBERLA_LOG_WARNING_ON_ROOT("[RodCaps] DIAGNOSTIC artificial caps ACTIVE for rods: hydroForce=" << g_rodCaps.hydroForce << " hydroTorque=" << g_rodCaps.hydroTorque << " velocity=" << g_rodCaps.velocity
            << " angularVelocity=" << g_rodCaps.angularVelocity << " overlap=" << g_rodCaps.overlap << " (0 = that cap off) -- changes the physics, not for production");
      computeContacts.setLogContactAnomalies(simParams->model.logContactAnomalies,
         simParams->model.contactAnomalyStartStep, simParams->model.contactAnomalyEndStep);
      computeLubrication.setHistoryLogger(&floaterHistoryLogger_);
      if (simParams->model.logFloaterHistory)
         WALBERLA_LOG_INFO_ON_ROOT("[FloaterHistoryLogger] enabled -- per-rank CSVs under \""
            << simParams->model.floaterHistoryBaseFolder << "/\""
            << (simParams->model.logFloaterHistoryStartStep != uint_t(0)
                   ? (" starting at step " + std::to_string(simParams->model.logFloaterHistoryStartStep))
                   : std::string())
            << (simParams->model.logFloaterHistoryEndStep != uint_t(0)
                   ? (" ending at step " + std::to_string(simParams->model.logFloaterHistoryEndStep))
                   : std::string()));
      if (simParams->model.logContactAnomalies)
         WALBERLA_LOG_INFO_ON_ROOT("[ComputeContacts] logContactAnomalies enabled -- "
            "[NORMAL-FLIP]/[PROLONGED-CONTACT]/[SIMULTANEOUS-CONTACTS] watchers active"
            << (simParams->model.contactAnomalyStartStep != uint_t(0)
                   ? (" from step " + std::to_string(simParams->model.contactAnomalyStartStep))
                   : std::string())
            << (simParams->model.contactAnomalyEndStep != uint_t(0)
                   ? (" to step " + std::to_string(simParams->model.contactAnomalyEndStep))
                   : std::string()));

      // ModelParams::logFloaterHistoryReleaseStepFilter (0 = disabled = log all): build the UID allowlist from the
      // FloatingSphere/Box/Cylinder entries with matching settlingSteps; log() and logState() skip everyone else.
      if (simParams->model.logFloaterHistory && simParams->model.logFloaterHistoryReleaseStepFilter != uint_t(0))
      {
         const uint_t targetStep = simParams->model.logFloaterHistoryReleaseStepFilter;
         for (size_t fi = 0; fi < simParams->floatingSpheres.size(); ++fi)
            if (simParams->floatingSpheres[fi].settlingSteps == targetStep)
               loggedFloaterUIDs_.insert(simObjs->floaterUIDs[fi]);
         for (size_t bi = 0; bi < simParams->floatingBoxes.size(); ++bi)
            if (simParams->floatingBoxes[bi].settlingSteps == targetStep)
               loggedFloaterUIDs_.insert(simObjs->floaterBoxUIDs[bi]);
         for (size_t ci = 0; ci < simParams->floatingCylinders.size(); ++ci)
            if (simParams->floatingCylinders[ci].settlingSteps == targetStep)
               loggedFloaterUIDs_.insert(simObjs->floaterCylinderUIDs[ci]);
         floaterHistoryLogger_.setUidFilter(&loggedFloaterUIDs_);
         WALBERLA_LOG_INFO_ON_ROOT("[FloaterHistoryLogger] restricted to " << loggedFloaterUIDs_.size()
            << " floater(s) with settlingSteps==" << targetStep);
      }

      // UIDs are resolved by setupParticles via an indexed creation loop + allReduce(SUM), so every rank already
      // has the correct UID for every floater in simObjs->floaterUIDs / floaterBoxUIDs / floaterCylinderUIDs;
      // copy them instead of looking them up per rank.
      const size_t nF = simParams->floatingSpheres.size();
      floaterUIDs_ = simObjs->floaterUIDs;
      floaterRadii_.resize(nF);
      floaterDensities_.resize(nF);
      for (size_t fi = 0; fi < nF; ++fi)
      {
         floaterRadii_[fi]     = simParams->floatingSpheres[fi].radius;
         floaterDensities_[fi] = simParams->floatingSpheres[fi].densityRatio;
         // Translation-pinned spheres: same pin registration as boxes/cylinders below, keyed by UID.
         if (simParams->floatingSpheres[fi].pinPosition)
         {
            pinnedBoxPositions_[floaterUIDs_[fi]] = simParams->floatingSpheres[fi].position;
            pinnedBoxAxes_[floaterUIDs_[fi]]      = simParams->floatingSpheres[fi].pinAxes;
         }
         if (simParams->floatingSpheres[fi].pinRotation)
            pinnedRotationAxes_[floaterUIDs_[fi]] = simParams->floatingSpheres[fi].pinRotationAxes;
      }

      useMemOnlyBuoyancy_        = simParams->model.useMemOnlyBuoyancy;
      useFullVolumeAnalyticBuoyancy_ = simParams->model.useFullVolumeAnalyticBuoyancy;
      WALBERLA_LOG_INFO_ON_ROOT("[FloaterForces] buoyancy model: "
         << (simParams->model.useMemOnlyBuoyancy
                ? "MEM-ONLY (default; Bogner & Rude 2013 -- MEM/hydro_force alone supplies buoyancy, "
                  "body force is gravity-only, no vSub/righting-torque computation at all)"
                : "SPLIT (gravity minus analytic vSub*g, plus vSub-derived righting torque; "
                  "Antidunes' extended, opt-in, pre-default-flip model)"));
      WALBERLA_LOG_INFO_ON_ROOT("[FloaterForces] floater streamwise drive: "
         << (simParams->model.driveFloatersByMass
                ? "BY MASS (F_xy = rho_p * V * a_xy = m_p * a_xy)"
                : "FLUID FULL VOLUME (default; F_xy = V * a_xy)"));

      // Combined box + cylinder "clip floater" list: both shapes share the plane-fit + clip V_sub
      // pipeline in computeSubmergedVolumesBox(); a cylinder is described by halfExtents = (r, r, L/2)
      // with its axis on the body z (matching the mesa_pd::data::Cylinder convention).
      for (size_t bi = 0; bi < simParams->floatingBoxes.size(); ++bi)
      {
         const auto& fbp = simParams->floatingBoxes[bi];
         clipFloaters_.push_back({ simObjs->floaterBoxUIDs[bi], false, fbp.halfExtents });
         // Per-floater COM offset below the centroid (constant); keyed by UID for the righting torque.
         if (std::abs(fbp.comOffsetZ) > real_t(0))
            boxComOffsetZ_[simObjs->floaterBoxUIDs[bi]] = fbp.comOffsetZ;
         // Translation-pinned floaters: record the target world position + per-axis pin mask; rotation stays free.
         if (fbp.pinPosition)
         {
            pinnedBoxPositions_[simObjs->floaterBoxUIDs[bi]] = fbp.position;
            pinnedBoxAxes_[simObjs->floaterBoxUIDs[bi]]      = fbp.pinAxes;
         }
         if (fbp.pinRotation)
            pinnedRotationAxes_[simObjs->floaterBoxUIDs[bi]] = fbp.pinRotationAxes;
      }
      for (size_t ci = 0; ci < simParams->floatingCylinders.size(); ++ci)
      {
         const auto& fcp = simParams->floatingCylinders[ci];
         clipFloaters_.push_back({ simObjs->floaterCylinderUIDs[ci], true,
                                   Vector3< real_t >(fcp.radius, fcp.radius, real_t(0.5) * fcp.length) });
         if (std::abs(fcp.comOffsetZ) > real_t(0))
            boxComOffsetZ_[simObjs->floaterCylinderUIDs[ci]] = fcp.comOffsetZ;
         if (fcp.pinPosition)
         {
            pinnedBoxPositions_[simObjs->floaterCylinderUIDs[ci]] = fcp.position;
            pinnedBoxAxes_[simObjs->floaterCylinderUIDs[ci]]      = fcp.pinAxes;
         }
         if (fcp.pinRotation)
            pinnedRotationAxes_[simObjs->floaterCylinderUIDs[ci]] = fcp.pinRotationAxes;
      }
      WALBERLA_LOG_INFO_ON_ROOT("[BoxBuoyancy] V_sub method (box + cylinder): "
         << "SAMPLING (12³ points)"
         << "  pinned-position floaters: " << pinnedBoxPositions_.size());

      WALBERLA_LOG_INFO_ON_ROOT("[DynamicH] MEM reference height: "
         << (useDynamicReferenceHeight_ ? "DYNAMIC (EMA-tracked from floaters' fitted local waterlines, alpha=0.2)"
                                        : "STATIC (flat H=" + std::to_string(surfaceHeight_) + ")"));
   }

   // MPI-sums computeContacts' per-rank closed-form scenario / GJK/EPA fallback tallies (all zero unless
   // ModelParameters.useAnalyticalContactDetection is on) and logs one line on root. Called once after the main
   // timeloop (stepLabel="", see Antidunes.cpp) and, with EvaluationParameters.logContactDetectionStatsPeriodically,
   // also every performanceLogFrequency steps (stepLabel = step) followed by resetContactDetectionStatistics(), so
   // each line reports per-interval counts and the sub-case mix can be followed as a jam develops.
   void logContactDetectionStatistics(const std::string& stepLabel = "") const
   {
      const auto& s = computeContacts.contactStats();
      uint64_t boxBoxSAT = s.boxBoxSAT, cylFaceFace = s.cylFaceFace,
               cylBandBandParallel = s.cylBandBandParallel, cylFaceBand = s.cylFaceBand,
               cylBandBandSkewed = s.cylBandBandSkewed, cylFaceEdge = s.cylFaceEdge,
               cylSpecialContact = s.cylSpecialContact, cylBandEdge = s.cylBandEdge,
               cylEdgeEdge = s.cylEdgeEdge, cylFallbackTier1 = s.cylFallbackTier1,
               cylFallbackSpecialContact = s.cylFallbackSpecialContact,
               cylPureTier1Rejected = s.cylPureTier1Rejected,
               cylPureSpecialContactTrusted = s.cylPureSpecialContactTrusted,
               cylFallbackTier1GjkFoundContact            = s.cylFallbackTier1GjkFoundContact,
               cylFallbackTier1GjkFoundNoContact           = s.cylFallbackTier1GjkFoundNoContact,
               cylFallbackSpecialContactGjkFoundContact   = s.cylFallbackSpecialContactGjkFoundContact,
               cylFallbackSpecialContactGjkFoundNoContact = s.cylFallbackSpecialContactGjkFoundNoContact,
               cylTier1CheaplyProvedSeparated = s.cylTier1CheaplyProvedSeparated,
               boxTightLayer1CheaplyProvedSeparated = s.boxTightLayer1CheaplyProvedSeparated,
               cylCapsuleRejected = s.cylCapsuleRejected,
               cylRefinedContact = s.cylRefinedContact, cylRefinedSeparated = s.cylRefinedSeparated;

      mpi::allReduceInplace(boxBoxSAT,           mpi::SUM);
      mpi::allReduceInplace(cylFaceFace,         mpi::SUM);
      mpi::allReduceInplace(cylBandBandParallel, mpi::SUM);
      mpi::allReduceInplace(cylFaceBand,         mpi::SUM);
      mpi::allReduceInplace(cylBandBandSkewed,   mpi::SUM);
      mpi::allReduceInplace(cylFaceEdge,         mpi::SUM);
      mpi::allReduceInplace(cylSpecialContact,   mpi::SUM);
      mpi::allReduceInplace(cylBandEdge,         mpi::SUM);
      mpi::allReduceInplace(cylEdgeEdge,         mpi::SUM);
      mpi::allReduceInplace(cylFallbackTier1,    mpi::SUM);
      mpi::allReduceInplace(cylFallbackSpecialContact, mpi::SUM);
      mpi::allReduceInplace(cylPureTier1Rejected, mpi::SUM);
      mpi::allReduceInplace(cylPureSpecialContactTrusted, mpi::SUM);
      mpi::allReduceInplace(cylFallbackTier1GjkFoundContact,   mpi::SUM);
      mpi::allReduceInplace(cylFallbackTier1GjkFoundNoContact, mpi::SUM);
      mpi::allReduceInplace(cylFallbackSpecialContactGjkFoundContact,   mpi::SUM);
      mpi::allReduceInplace(cylFallbackSpecialContactGjkFoundNoContact, mpi::SUM);
      mpi::allReduceInplace(cylTier1CheaplyProvedSeparated, mpi::SUM);
      mpi::allReduceInplace(boxTightLayer1CheaplyProvedSeparated, mpi::SUM);
      mpi::allReduceInplace(cylCapsuleRejected, mpi::SUM);
      mpi::allReduceInplace(cylRefinedContact, mpi::SUM);
      mpi::allReduceInplace(cylRefinedSeparated, mpi::SUM);
      constexpr int nP = SolidCylinderContactDetection::ContactDetectionStats::kRefinePaths;
      std::vector< uint64_t > rp(s.cylRefinePath, s.cylRefinePath + nP);
      std::vector< double >   rt(s.cylRefinePathSeconds, s.cylRefinePathSeconds + nP);
      std::vector< uint64_t > rf(std::begin(s.cylRefineFastReject), std::end(s.cylRefineFastReject));
      std::vector< uint64_t > rs(std::begin(s.cylRefineDescendedScenario), std::end(s.cylRefineDescendedScenario));
      uint64_t rk = s.cylRefineKinks, rc = s.cylRefineConvex, sfc = s.shadowFallbackCalls;
      double   scf = s.shadowClosedFormSeconds, sfb = s.shadowFallbackSeconds;
      mpi::allReduceInplace(sfc, mpi::SUM);
      mpi::allReduceInplace(scf, mpi::SUM);
      mpi::allReduceInplace(sfb, mpi::SUM);
      mpi::allReduceInplace(rp, mpi::SUM);
      mpi::allReduceInplace(rt, mpi::SUM);
      mpi::allReduceInplace(rf, mpi::SUM);
      mpi::allReduceInplace(rs, mpi::SUM);
      mpi::allReduceInplace(rk, mpi::SUM);
      mpi::allReduceInplace(rc, mpi::SUM);

      // cylPureSpecialContactTrusted is a sub-metric of cylAnalyticalTrusted (that path also calls bumpCylScenario),
      // so it is reported separately, not added again. cylPureTier1Rejected is not folded into any scenario, so it
      // is added to the grand total like a fallback attempt.
      const uint64_t cylAnalyticalTrusted = cylFaceFace + cylBandBandParallel + cylFaceBand +
         cylBandBandSkewed + cylFaceEdge + cylSpecialContact + cylBandEdge + cylEdgeEdge;
      const uint64_t cylFallbackTotal = cylFallbackTier1 + cylFallbackSpecialContact;
      const uint64_t cylTotal = cylAnalyticalTrusted + cylFallbackTotal + cylPureTier1Rejected +
         cylTier1CheaplyProvedSeparated;
      const real_t fallbackRatePercent =
         (cylTotal > uint64_t(0)) ? real_t(100) * real_t(cylFallbackTotal) / real_t(cylTotal) : real_t(0);

      WALBERLA_LOG_INFO_ON_ROOT(
         "[ContactDetectionStats]" << (stepLabel.empty() ? "" : " step=" + stepLabel)
         << " Box x Box: SAT=" << boxBoxSAT
         << " (Layer1-cheaply-proved-separated=" << boxTightLayer1CheaplyProvedSeparated << ")"
         << " | Cylinder x Cylinder analytical (trusted): face-face=" << cylFaceFace
         << " band-band-parallel=" << cylBandBandParallel
         << " face-band=" << cylFaceBand
         << " band-band-skewed=" << cylBandBandSkewed
         << " face-edge=" << cylFaceEdge
         << " special-contact=" << cylSpecialContact
         << " band-edge=" << cylBandEdge
         << " edge-edge=" << cylEdgeEdge
         << " (trusted total=" << cylAnalyticalTrusted
         << ", of which pure-mode-special-contact-trusted-anyway=" << cylPureSpecialContactTrusted << ")"
         << " | GJK/EPA fallback: tier1(analytical-found-nothing)=" << cylFallbackTier1
         << " [of which GJK says: contact=" << cylFallbackTier1GjkFoundContact
         << " no-contact=" << cylFallbackTier1GjkFoundNoContact << "]"
         << " special-contact(unconditionally-distrusted)=" << cylFallbackSpecialContact
         << " [of which GJK says: contact=" << cylFallbackSpecialContactGjkFoundContact
         << " no-contact=" << cylFallbackSpecialContactGjkFoundNoContact << "]"
         << " (fallback total=" << cylFallbackTotal << ")"
         << " | Tier 1 CHEAPLY PROVED SEPARATED (no GJK/EPA call needed at all)="
         << cylTier1CheaplyProvedSeparated
         << " | usePureAnalyticalCylinder rejections (no fallback attempted): tier1="
         << cylPureTier1Rejected
         << " | capsule-bound rejected analytic false positives=" << cylCapsuleRejected
         << " | refined (no GJK/EPA): contact=" << cylRefinedContact << " separated=" << cylRefinedSeparated
         << " | cylinder-cylinder grand total=" << cylTotal
         << " fallback-rate=" << fallbackRatePercent << "%");

      if (std::accumulate(rp.begin(), rp.end(), uint64_t(0)) > uint64_t(0))
      {
         static const char* pathNames[nP] = { "exact", "fast-accepted", "contact-descended", "contact-descended-separated",
                                              "no-scenario-separated", "no-scenario-contact" };
         std::ostringstream os;
         for (int k = 0; k < nP; ++k)
            os << " " << pathNames[k] << "=" << rp[size_t(k)] << " (" << rt[size_t(k)] << " s, "
               << (rp[size_t(k)] ? 1e6 * rt[size_t(k)] / double(rp[size_t(k)]) : 0.0) << " us)";
         os << " | fast-path rejects: depth-mismatch=" << rf[1] << " not-locally-optimal=" << rf[2]
            << " smaller-elsewhere=" << rf[3] << " | kink search=" << rk << " convex check=" << rc
            << " | re-descended contacts by closed-form scenario: face-face=" << rs[0] << " band-band-parallel=" << rs[1]
            << " face-band=" << rs[2] << " band-band-skewed=" << rs[3] << " face-edge=" << rs[4]
            << " special-contact=" << rs[5] << " band-edge=" << rs[6] << " edge-edge=" << rs[7];
         if (scf > 0.0)
            os << " | shadow: closed form " << scf << " s, fallback GJK/EPA " << sfb << " s (" << sfc << " calls), refinement "
               << std::accumulate(rt.begin(), rt.end(), 0.0) << " s";
         WALBERLA_LOG_INFO_ON_ROOT("[RefineProfile]" << (stepLabel.empty() ? "" : " step=" + stepLabel)
                                   << " (time summed over ranks)" << os.str());
      }

      if (g_rodCaps.any())
      {
         uint64_t c[6] = { g_rodCapStats.hydroForce, g_rodCapStats.hydroTorque, g_rodCapStats.velocity, g_rodCapStats.angularVelocity, g_rodCapStats.overlap, g_rodCapStats.nonFinite };
         real_t m[5] = { g_rodCapStats.maxHydroForce, g_rodCapStats.maxHydroTorque, g_rodCapStats.maxVelocity, g_rodCapStats.maxAngularVelocity, g_rodCapStats.maxOverlap };
         for (int i = 0; i < 6; ++i) mpi::allReduceInplace(c[i], mpi::SUM);
         for (int i = 0; i < 5; ++i) mpi::allReduceInplace(m[i], mpi::MAX);
         WALBERLA_LOG_INFO_ON_ROOT("[RodCaps]" << (stepLabel.empty() ? "" : " step=" + stepLabel) << " cap hits (cumulative): hydroForce=" << c[0] << " hydroTorque=" << c[1] << " velocity=" << c[2] << " angularVelocity=" << c[3]
            << " overlap=" << c[4] << " non-finite-zeroed=" << c[5] << " | largest value seen before capping: |F_hydro|=" << m[0] << " |T_hydro|=" << m[1] << " |v|=" << m[2] << " |omega|=" << m[3] << " overlap=" << m[4]);
      }

      if (simParams_shadowEnabled_)
      {
         const auto& sh = computeContacts.shadowStats();
         constexpr size_t K = SolidCylinderContactDetection::kShadowScenarios;
         std::vector< uint64_t > cnt(K * 5 + 4, uint64_t(0)); std::vector< double > sums(K, 0.0), mx1(K, 0.0), mx2(K, 0.0);
         for (size_t i = 0; i < K; ++i)
         {
            const auto& b = sh.byScenario[i];
            cnt[i * 5 + 0] = b.n; cnt[i * 5 + 1] = b.wrong; cnt[i * 5 + 2] = b.verifierFlagged; cnt[i * 5 + 3] = b.normalFlipped; cnt[i * 5 + 4] = b.normalOff;
            sums[i] = b.sumAbsErr; mx1[i] = b.maxUnder; mx2[i] = b.maxOver;
         }
         cnt[K * 5 + 0] = sh.gjkContacts; cnt[K * 5 + 1] = sh.bothContact; cnt[K * 5 + 2] = sh.closedFormOnly; cnt[K * 5 + 3] = sh.gjkOnly;
         mpi::allReduceInplace(cnt, mpi::SUM); mpi::allReduceInplace(sums, mpi::SUM); mpi::allReduceInplace(mx1, mpi::MAX); mpi::allReduceInplace(mx2, mpi::MAX);
         std::ostringstream os;
         os << "[ShadowAnalytic]" << (stepLabel.empty() ? "" : " step=" + stepLabel) << " GJK/EPA is USED; closed form computed only for comparison. GJK contacts=" << cnt[K * 5 + 0]
            << " closed-form also finds contact=" << cnt[K * 5 + 1] << " closed-form finds NOTHING where GJK has a contact=" << cnt[K * 5 + 3] << " closed-form contact where GJK has none=" << cnt[K * 5 + 2]
            << " (counts include duplicate computations of the same pair on neighbouring ranks). Per closed-form scenario [n | wrong(|d|>5%R) | mean|d| | max under | max over | check flags | normal flipped | normal off>20deg]:";
         for (size_t i = 0; i < K; ++i)
            if (cnt[i * 5] > 0)
               os << "\n   " << SolidCylinderContactDetection::shadowScenarioName(i) << ": n=" << cnt[i * 5] << " wrong=" << cnt[i * 5 + 1] << " (" << (100.0 * double(cnt[i * 5 + 1]) / double(cnt[i * 5]))
                  << "%) mean|d|=" << sums[i] / double(cnt[i * 5]) << " maxUnder=" << mx1[i] << " maxOver=" << mx2[i] << " checkFlags=" << cnt[i * 5 + 2] << " normalFlipped=" << cnt[i * 5 + 3] << " normalOff=" << cnt[i * 5 + 4];
         WALBERLA_LOG_INFO_ON_ROOT(os.str());
      }
   }

   // Zeroes computeContacts' per-rank ContactDetectionStats tallies; call right after
   // logContactDetectionStatistics() in periodic reporting so the next line covers only the new interval.
   // Diagnostic counters only, no effect on the physics.
   void resetContactDetectionStatistics() { computeContacts.resetContactStats(); }

   template< typename TP >
   void resetHydroForces(timing::TimingPool<TP> & tp){
      tp["MESAPD: ResetHydrodynamicForceTorqueKernel"].start();
      particleStorage->forEachParticle(useOpenMP,
         mesa_pd::kernel::SelectAll(), *particleAccessor,
         lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel{}, *particleAccessor);
      tp["MESAPD: ResetHydrodynamicForceTorqueKernel"].end();
   }

   template< typename TP >
   void reduceAndAvgHydroForces(timing::TimingPool<TP> & tp, bool initializeOldValues){
      // [PRE-REDUCE]: this rank's local hydrodynamic accumulator for each floater (local and ghost copies) before
      // any cross-rank reduction. [NAN-LINK] (SimpleBB.h) shows each individual add is finite; this checks that
      // every rank's pre-reduce total entering ReduceProperty is finite too. Covers all floaters this rank knows,
      // ghosts included (a ghost's accumulator also enters the reduction). Once per main step; pushed into
      // diagnosticEvents_ (gated by diagnosticSweepEnabled_).
      {
         particleStorage->forEachParticle(useOpenMP,
            [](size_t idx, auto& ac) {
               const uint_t t = ac.getType(idx);
               return t == uint_t(1) || t == uint_t(3) || t == uint_t(4);
            },
            *particleAccessor,
            [this](size_t idx, auto& ac) {
               using namespace mesa_pd::data::particle_flags;
               const auto& hf = ac.getHydrodynamicForce(idx);
               const auto& ht = ac.getHydrodynamicTorque(idx);
               const bool finite = std::isfinite(hf[0]) && std::isfinite(hf[1]) && std::isfinite(hf[2]) &&
                                    std::isfinite(ht[0]) && std::isfinite(ht[1]) && std::isfinite(ht[2]);
               if (!finite && diagnosticSweepEnabled_)
               {
                  std::ostringstream detail;
                  detail << "isGhost=" << isSet(ac.getFlags(idx), GHOST) << " owner=" << ac.getOwner(idx)
                         << " local pre-reduce hydForce=" << hf << " hydTorque=" << ht;
                  diagnosticEvents_->push_back(
                     { stepCounter_, uint_t(0), "PRE-REDUCE", ac.getUid(idx), ac.getPosition(idx), real_t(0), detail.str() });
               }
            },
            *particleAccessor);
      }

      tp["MESAPD: ReduceProperty(HydrodynamicForceTorqueNotification)"].start();
      mesa_pd::mpi::ReduceProperty{}.operator()< mesa_pd::HydrodynamicForceTorqueNotification >(*particleStorage);
      tp["MESAPD: ReduceProperty(HydrodynamicForceTorqueNotification)"].end();

      // [RAW-HYDRO]: with every per-link contribution finite ([NAN-LINK], SimpleBB.h), a non-finite force seen by
      // AddExternalForces must enter here: the MPI-summed raw force after ReduceProperty, or the blend with a stale
      // "old" value in AverageHydrodynamicForceTorqueKernel. Logs both raw and old values to tell them apart.
      // SelectLocal only, once per main step; pushed into diagnosticEvents_ (gated by diagnosticSweepEnabled_).
      {
         particleStorage->forEachParticle(useOpenMP,
            [](size_t idx, auto& ac) {
               using namespace mesa_pd::data::particle_flags;
               if (isSet(ac.getFlags(idx), GHOST)) return false; // SelectLocal semantics, see doc comment above
               const uint_t t = ac.getType(idx);
               return t == uint_t(1) || t == uint_t(3) || t == uint_t(4); // sphere/box/cylinder floaters
            },
            *particleAccessor,
            [this](size_t idx, auto& ac) {
               const auto& hf  = ac.getHydrodynamicForce(idx);
               const auto& ht  = ac.getHydrodynamicTorque(idx);
               const auto& ohf = ac.getOldHydrodynamicForce(idx);
               const auto& oht = ac.getOldHydrodynamicTorque(idx);
               const bool rawFinite = std::isfinite(hf[0]) && std::isfinite(hf[1]) && std::isfinite(hf[2]) &&
                                       std::isfinite(ht[0]) && std::isfinite(ht[1]) && std::isfinite(ht[2]);
               const bool oldFinite = std::isfinite(ohf[0]) && std::isfinite(ohf[1]) && std::isfinite(ohf[2]) &&
                                       std::isfinite(oht[0]) && std::isfinite(oht[1]) && std::isfinite(oht[2]);
               if ((!rawFinite || !oldFinite) && diagnosticSweepEnabled_)
               {
                  std::ostringstream detail;
                  detail << "rawFinite=" << rawFinite << " oldFinite=" << oldFinite
                         << " raw hydForce=" << hf << " hydTorque=" << ht
                         << " old hydForce=" << ohf << " hydTorque=" << oht;
                  diagnosticEvents_->push_back(
                     { stepCounter_, uint_t(0), "RAW-HYDRO", ac.getUid(idx), ac.getPosition(idx), real_t(0), detail.str() });
               }
            },
            *particleAccessor);
      }

      if (initializeOldValues)
      {
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel{}, *particleAccessor);
      }

      // ModelParameters.useHydrodynamicForceTimeAveraging (SimulationParams.h): AverageHydrodynamicForceTorqueKernel
      // blends this step's reduced force/torque with the previous step's (F = 0.5*(F(t)+F(t-Δt)), Ladd 1994), as
      // Watanabe et al. 2023 (Eq. 48-49, citing Aidun et al. 1998) do "to improve the stability of the fluid-solid
      // coupling with floating objects". Default true; false feeds DEM the raw reduced force (new is still copied
      // into old, so the old values stay valid).
      if (useHydrodynamicForceTimeAveraging_)
      {
         tp["MESAPD: AverageHydrodynamicForceTorqueKernel"].start();
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel{}, *particleAccessor);
         tp["MESAPD: AverageHydrodynamicForceTorqueKernel"].end();
      }
      else
      {
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            [](size_t idx, auto& ac) {
               ac.setOldHydrodynamicForce(idx, ac.getHydrodynamicForce(idx));
               ac.setOldHydrodynamicTorque(idx, ac.getHydrodynamicTorque(idx));
            },
            *particleAccessor);
      }
   }

   template< typename TP >
   void doDEMSubCycles(timing::TimingPool<TP> & tp, uint_t numSubCycles = uint_t(1)){
      // One increment per LBM step, independent of whether spheres/boxes are present (the submerged-volume
      // functions return early for empty lists; gating the counter on them would stall it).
      ++stepCounter_;
      // Submerged volumes, once per LBM step (the fill field is constant during the DEM sub-cycles); both
      // functions write into submergedVolumes_ (keyed by UID). Skipped under useMemOnlyBuoyancy_ (default): its
      // force term does not depend on vSub and its righting torque is the raw MEM torque alone (Bogner & Rude
      // 2013, Eq. 18); running them would add the analytic buoyancy-lever torque via buoyancyCenterOffset_.
      if (!useMemOnlyBuoyancy_)
      {
         computeSubmergedVolumes();
         computeSubmergedVolumesBox();
      }
      addExternalForces.setSubmergedVolumes(&submergedVolumes_);
      addExternalForces.setBuoyancyCenterOffset(&buoyancyCenterOffset_);
      addExternalForces.setComOffsetZ(&boxComOffsetZ_);
      addExternalForces.setStepInfo(&stepCounter_);

      const bool logForcesThisStep = (stepCounter_ % evalFrequency_ == uint_t(0));
      for (auto subCycle = uint_t(0); subCycle < numSubCycles; ++subCycle)
      {
         // [GHOST-VANISHED] snapshot, once per main step (subCycle==0, before this step's sub-cycles/syncs):
         // every floater UID this rank can see (local or ghost), i.e. the population the fluid-obstacle map
         // (particleField_, refreshed once per main step) assumed. A ghost dropped by a later syncCall() before
         // the next refresh leaves particleField_ with an unresolvable UID (the STALE-PARTICLE-IDX precondition),
         // which the local-only [SUBCYCLE-LEFT-RANK] check below cannot see.
         if (subCycle == uint_t(0))
         {
            ghostMappingSnapshot_.clear();
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  const uint_t t = ac.getType(idx);
                  return t == uint_t(1) || t == uint_t(3) || t == uint_t(4); // sphere/box/cylinder floaters
               },
               *particleAccessor,
               [this](size_t idx, auto& ac) { ghostMappingSnapshot_.insert(ac.getUid(idx)); },
               *particleAccessor);
         }

         // [VELOCITY-JUMP]: catches the onset of a runaway between the once-per-evalFrequency_ floater logs. Once
         // per main step (subCycle==0, velocity as left by the previous main step), compares against
         // lastVelocityForJumpCheck_ and flags a change above 0.05 LU/step per main step: about two orders of
         // magnitude above ordinary floater dynamics (<<0.001) and below runaway values (tens to hundreds).
         // Events go into diagnosticEvents_ (gated by diagnosticSweepEnabled_, off by default), drained once per
         // step by BuildTimeloop.h's addAfterTimeStepDiagnosticSweep. The velocity bookkeeping always runs (cheap),
         // so enabling the sweep mid-run has a valid baseline.
         if (subCycle == uint_t(0))
         {
            constexpr real_t velocityJumpThreshold = real_t(0.05);
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  using namespace mesa_pd::data::particle_flags;
                  if (isSet(ac.getFlags(idx), GHOST)) return false;
                  const uint_t t = ac.getType(idx);
                  return t == uint_t(1) || t == uint_t(3) || t == uint_t(4); // sphere/box/cylinder floaters
               },
               *particleAccessor,
               [this, subCycle](size_t idx, auto& ac) {
                  const walberla::id_t uid = ac.getUid(idx);
                  const auto& vel = ac.getLinearVelocity(idx);
                  auto it = lastVelocityForJumpCheck_.find(uid);
                  if (it != lastVelocityForJumpCheck_.end())
                  {
                     const real_t deltaV = (vel - it->second).length();
                     if (diagnosticSweepEnabled_ && deltaV > velocityJumpThreshold)
                     {
                        std::ostringstream detail;
                        detail << "prevVel=" << it->second << " newVel=" << vel << " |deltaV|=" << deltaV
                               << " threshold=" << velocityJumpThreshold
                               << " -- velocity changed by more than the threshold over this ONE main "
                                  "step -- catches the onset of a runaway, not just its aftermath.";
                        diagnosticEvents_->push_back(
                           { stepCounter_, subCycle, "VELOCITY-JUMP", uid, ac.getPosition(idx), deltaV, detail.str() });
                     }
                  }
                  lastVelocityForJumpCheck_[uid] = vel;
               },
               *particleAccessor);
         }

         // Log the per-floater force breakdown once per eval step (only on the first sub-cycle).
         addExternalForces.setLogForces(logForcesThisStep && subCycle == uint_t(0));
         // [MEM-WATCH] sub-cycle index (see AddExternalForces::setSubCycle), every sub-cycle.
         addExternalForces.setSubCycle(subCycle);
         // FloaterHistoryLogger's (mainStep, subCycle) tag for this sub-cycle -- no-op when disabled.
         computeContacts.setStepInfo(stepCounter_, subCycle);
         computeLubrication.setStepInfo(stepCounter_, subCycle);

         // velocity verlet pre-force (update positions) --------------------------------
         tp["MESAPD: VelocityVerletPreForceUpdate"].start();
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            vvIntegratorPreForce, *particleAccessor);
         tp["MESAPD: VelocityVerletPreForceUpdate"].end();

         // [SUBCYCLE-LEFT-RANK]: low-volume safety check (50 hits per rank) whether DEM sub-cycle motion has moved
         // a LOCAL floater outside every block this rank owns before this sub-cycle's syncCall() hands off
         // ownership. The obstacle mapping / ghost sync assume this cannot happen within one sub-cycle; a hit
         // preceding a STALE-PARTICLE-IDX abort shows the assumption broke. Exact (local blocks only).
         {
            static uint_t leftRankLogCount = uint_t(0);
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  using namespace mesa_pd::data::particle_flags;
                  if (isSet(ac.getFlags(idx), GHOST)) return false;
                  const uint_t t = ac.getType(idx);
                  return t == uint_t(1) || t == uint_t(3) || t == uint_t(4); // sphere/box/cylinder floaters
               },
               *particleAccessor,
               [this, subCycle](size_t idx, auto& ac) {
                  const auto& pos = ac.getPosition(idx);
                  real_t minSqDist = std::numeric_limits< real_t >::max();
                  for (auto blockIt = blockForest_->begin(); blockIt != blockForest_->end(); ++blockIt)
                  {
                     const real_t d = blockIt->getAABB().sqDistance(pos);
                     if (d < minSqDist) minSqDist = d;
                     if (minSqDist <= real_t(0)) break;
                  }
                  static uint_t& logCount = leftRankLogCount;
                  if (minSqDist > real_t(0) && logCount < uint_t(50))
                  {
                     WALBERLA_LOG_WARNING("[SUBCYCLE-LEFT-RANK] rank=" << mpi::MPIManager::instance()->rank()
                        << " uid=" << ac.getUid(idx) << " step=" << stepCounter_ << " subCycle=" << subCycle
                        << " pos=" << pos << " vel=" << ac.getLinearVelocity(idx)
                        << " distanceOutsideNearestLocalBlock=" << std::sqrt(minSqDist)
                        << " -- this floater's position is now outside every block this rank owns, "
                           "BEFORE this sub-cycle's syncCall() has had a chance to hand off ownership.");
                     ++logCount;
                  }
               },
               *particleAccessor);
         }

         // synchronize -----------------------------------------------------------------
         tp["MESAPD: syncCall 1"].start();
         syncCall();
         tp["MESAPD: syncCall 1"].end();

         // [GHOST-VANISHED] check (see the snapshot above): after this sub-cycle's syncCall(), a snapshot UID that
         // uidToIdx() reports invalid means this rank lost its only visibility (ghost or local) of a particle that
         // particleField_ may still cite as a cell's obstacle owner. Such UIDs are only erased from the snapshot;
         // no warning or diagnostic event is emitted. Clearing the vanished particle's orphaned obstacle cells here
         // is deliberately not done: the abrupt fill-level discontinuity destabilises the free-surface
         // reconstruction; the SimpleBB.h STALE-PARTICLE-IDX abort stays.
         for (auto it = ghostMappingSnapshot_.begin(); it != ghostMappingSnapshot_.end(); )
         {
            if (particleAccessor->uidToIdx(*it) == particleAccessor->getInvalidIdx())
               it = ghostMappingSnapshot_.erase(it);
            else
               ++it;
         }

         // [MIGRATION-NAN] (see lastKnownOwner_): a LOCAL floater's hydrodynamicForce/oldHydrodynamicForce is
         // non-finite right after this sub-cycle's syncCall(); reports whether its owner rank changed since this
         // rank last saw it (fresh migration vs. a force gone bad otherwise). Pushed into diagnosticEvents_ (gated
         // by diagnosticSweepEnabled_); the lastKnownOwner_ bookkeeping always runs (cheap).
         {
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  using namespace mesa_pd::data::particle_flags;
                  if (isSet(ac.getFlags(idx), GHOST)) return false;
                  const uint_t t = ac.getType(idx);
                  return t == uint_t(1) || t == uint_t(3) || t == uint_t(4);
               },
               *particleAccessor,
               [this, subCycle](size_t idx, auto& ac) {
                  const auto& hf  = ac.getHydrodynamicForce(idx);
                  const auto& ht  = ac.getHydrodynamicTorque(idx);
                  const auto& ohf = ac.getOldHydrodynamicForce(idx);
                  const auto& oht = ac.getOldHydrodynamicTorque(idx);
                  const bool finite = std::isfinite(hf[0]) && std::isfinite(hf[1]) && std::isfinite(hf[2]) &&
                                       std::isfinite(ht[0]) && std::isfinite(ht[1]) && std::isfinite(ht[2]) &&
                                       std::isfinite(ohf[0]) && std::isfinite(ohf[1]) && std::isfinite(ohf[2]) &&
                                       std::isfinite(oht[0]) && std::isfinite(oht[1]) && std::isfinite(oht[2]);
                  const walberla::id_t uid = ac.getUid(idx);
                  const int curOwner = ac.getOwner(idx);
                  const auto it = lastKnownOwner_.find(uid);
                  const bool everSeenOnThisRank = (it != lastKnownOwner_.end());
                  const bool justMigratedHere = everSeenOnThisRank && it->second != curOwner;
                  lastKnownOwner_[uid] = curOwner;
                  if (!finite && diagnosticSweepEnabled_)
                  {
                     std::ostringstream detail;
                     detail << "curOwner=" << curOwner << " everSeenBeforeOnThisRank=" << everSeenOnThisRank
                            << " justMigratedHereSinceLastSeen=" << justMigratedHere
                            << " hydForce=" << hf << " hydTorque=" << ht
                            << " oldHydForce=" << ohf << " oldHydTorque=" << oht;
                     diagnosticEvents_->push_back(
                        { stepCounter_, subCycle, "MIGRATION-NAN", uid, ac.getPosition(idx), real_t(0), detail.str() });
                  }
               },
               *particleAccessor);
         }

         // ModelParameters.logFloaterHistory world-state snapshot (no-op when disabled): this rank's view of every
         // floater right after the sync used by ComputeLubrication/ComputeContacts (see
         // FloaterHistoryLogger::logState). Dormant (startInactive, not yet released) floaters are skipped, as in
         // pairSelector below: they sit motionless at their spawn position.
         if (floaterHistoryLogger_.enabled())
            particleStorage->forEachParticle(useOpenMP,
               [dormant = dormantFloaterUIDs_](size_t idx, auto& ac) {
                  const uint_t t = ac.getType(idx);
                  if (t != uint_t(1) && t != uint_t(3) && t != uint_t(4)) return false; // sphere/box/cylinder floaters
                  return !dormant || dormant->count(ac.getUid(idx)) == 0;
               },
               *particleAccessor,
               [this, subCycle](size_t idx, auto& ac) {
                  const real_t rho = sampleFluidDensityAt(ac.getPosition(idx));
                  floaterHistoryLogger_.logState(stepCounter_, subCycle, ac, idx, rho);
               },
               *particleAccessor);

         // Dormant (startInactive, not yet at settlingSteps) floaters must not participate in any DEM
         // contact -- mesa_pd has no selector-combinator type, so this ANDs ExcludeInfiniteInfinite's
         // check with a dormant-UID lookup on both indices in one local lambda (same pattern as
         // SetupStages::floaterSelector / BuildTimeloop::movingObstacleSelector's fluid-mapping guard).
         auto dormant = dormantFloaterUIDs_;
         auto pairSelector = [dormant](size_t idx, size_t jdx, auto& ac) {
            using namespace mesa_pd::data::particle_flags;
            if (isSet(ac.getFlags(idx), INFINITE) && isSet(ac.getFlags(jdx), INFINITE)) return false;
            if (dormant->count(ac.getUid(idx)) || dormant->count(ac.getUid(jdx))) return false;
            return true;
         };

         // compute lubrication forces --------------------------------------------------
         if (useLubricationCorrection)
         {
            tp["MESAPD: ComputeLubrication"].start();
            particleStorage->forEachParticlePairHalf(useOpenMP,
               pairSelector, *particleAccessor,
               computeLubrication, *particleAccessor);
            tp["MESAPD: ComputeLubrication"].end();
         }

         // compute contact forces ------------------------------------------------------
         tp["MESAPD: ComputeContacts"].start();
         particleStorage->forEachParticlePairHalf(useOpenMP,
            pairSelector, *particleAccessor,
            computeContacts, *particleAccessor);
         // [SIMULTANEOUS-CONTACTS]: report and reset the per-body contact tally after the full sweep.
         computeContacts.reportAndResetSimultaneousContacts();
         tp["MESAPD: ComputeContacts"].end();

         // reduce contact history ------------------------------------------------------
         tp["MESAPD: ReduceContactHistory"].start();
         mesa_pd::mpi::ReduceContactHistory{}(*particleStorage);
         tp["MESAPD: ReduceContactHistory"].end();

         // add external forces ---------------------------------------------------------
         tp["MESAPD: AddExternalForces"].start();
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            addExternalForces, *particleAccessor);
         tp["MESAPD: AddExternalForces"].end();

         // reduce forces ---------------------------------------------------------------
         tp["MESAPD: ReduceProperty(ForceTorque)"].start();
         mesa_pd::mpi::ReduceProperty{}.operator()< mesa_pd::ForceTorqueNotification >(*particleStorage);
         tp["MESAPD: ReduceProperty(ForceTorque)"].end();

         // [TOTAL-FORCE-THRESHOLD]: companion to [STRONG-CONTACT-FORCE] (per-contact only). Several contacts on a
         // squeezed floater can each stay below that threshold but sum to a large force, so this checks the fully
         // accumulated, cross-rank reduced force (contacts + lubrication + external) once per sub-cycle,
         // before VelocityVerletPostForceUpdate consumes and resets it. Threshold 0.03: above the ordinary baseline
         // (O(1e-4..1e-3) for gravity+buoyancy+hydro), below STRONG_FORCE_THRESHOLD (0.1). 100 lines per rank.
         {
            static uint_t totalForceLogCount = uint_t(0);
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  using namespace mesa_pd::data::particle_flags;
                  if (isSet(ac.getFlags(idx), GHOST)) return false;
                  const uint_t t = ac.getType(idx);
                  return t == uint_t(1) || t == uint_t(3) || t == uint_t(4); // sphere/box/cylinder floaters
               },
               *particleAccessor,
               [this, subCycle](size_t idx, auto& ac) {
                  constexpr real_t TOTAL_FORCE_THRESHOLD = real_t(0.03);
                  const mesa_pd::Vec3& totalForce = ac.getForceRef(idx);
                  if (totalForce.length() > TOTAL_FORCE_THRESHOLD && totalForceLogCount < uint_t(100))
                  {
                     WALBERLA_LOG_WARNING("[TOTAL-FORCE-THRESHOLD] rank=" << mpi::MPIManager::instance()->rank()
                        << " uid=" << ac.getUid(idx) << " step=" << stepCounter_ << " subCycle=" << subCycle
                        << " pos=" << ac.getPosition(idx) << " vel=" << ac.getLinearVelocity(idx)
                        << " totalForceThisSubCycle=" << totalForce << " |total|=" << totalForce.length()
                        << " totalTorqueThisSubCycle=" << ac.getTorqueRef(idx)
                        << " hydroForce=" << ac.getHydrodynamicForce(idx)
                        << " -- fully-accumulated per-subcycle force exceeds " << TOTAL_FORCE_THRESHOLD
                        << " even though no single contact may have (see [STRONG-CONTACT-FORCE]'s own"
                           " per-contact-only threshold) -- catches summed-small-overlaps launches.");
                     ++totalForceLogCount;
                  }
               },
               *particleAccessor);
         }

         // velocity verlet post-force (update velocities) ------------------------------
         tp["MESAPD: VelocityVerletPostForceUpdate"].start();
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            vvIntegratorPostForce, *particleAccessor);
         tp["MESAPD: VelocityVerletPostForceUpdate"].end();

         // DIAGNOSTIC artificial velocity caps (ModelParameters.rodCapVelocity/AngularVelocity), rods (type 4),
         // local only
         if (g_rodCaps.velocity > real_t(0) || g_rodCaps.angularVelocity > real_t(0))
         {
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  using namespace mesa_pd::data::particle_flags;
                  if (isSet(ac.getFlags(idx), GHOST)) return false;
                  return ac.getType(idx) == uint_t(4);
               },
               *particleAccessor,
               [this, subCycle](size_t idx, auto& ac) {
                  auto& v = ac.getLinearVelocityRef(idx); auto& w = ac.getAngularVelocityRef(idx);
                  if (!(std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]))) { v = Vector3< real_t >(real_t(0)); ++g_rodCapStats.nonFinite; }
                  if (!(std::isfinite(w[0]) && std::isfinite(w[1]) && std::isfinite(w[2]))) { w = Vector3< real_t >(real_t(0)); ++g_rodCapStats.nonFinite; }
                  const real_t vm = v.length(), wm = w.length();
                  g_rodCapStats.maxVelocity = std::max(g_rodCapStats.maxVelocity, vm); g_rodCapStats.maxAngularVelocity = std::max(g_rodCapStats.maxAngularVelocity, wm);
                  if (g_rodCaps.velocity > real_t(0) && vm > g_rodCaps.velocity) { logRodCapHit("velocity", ac.getUid(idx), stepCounter_, subCycle, vm, g_rodCaps.velocity); v *= g_rodCaps.velocity / vm; ++g_rodCapStats.velocity; }
                  if (g_rodCaps.angularVelocity > real_t(0) && wm > g_rodCaps.angularVelocity) { logRodCapHit("angularVelocity", ac.getUid(idx), stepCounter_, subCycle, wm, g_rodCaps.angularVelocity); w *= g_rodCaps.angularVelocity / wm; ++g_rodCapStats.angularVelocity; }
                  return true;
               },
               *particleAccessor);
         }

         // TEST OPTION (ModelParameters.rodPrescribedVelocityX/Y/Z): rigid prescribed motion of all sphere (1), box (3)
         // and cylinder (4) floaters (capsules excluded), applied after the velocity update and before the sync so
         // ghosts see it.
         if (g_rodPrescribedVelocityOn)
         {
            particleStorage->forEachParticle(useOpenMP,
               [](size_t idx, auto& ac) {
                  using namespace mesa_pd::data::particle_flags;
                  if (isSet(ac.getFlags(idx), GHOST)) return false;
                  const uint_t t = ac.getType(idx);
                  return t == uint_t(1) || t == uint_t(3) || t == uint_t(4);
               },
               *particleAccessor,
               [](size_t idx, auto& ac) { ac.setLinearVelocity(idx, g_rodPrescribedVelocity); ac.setAngularVelocity(idx, Vector3< real_t >(real_t(0))); return true; },
               *particleAccessor);
         }

         // hold the pinned translational/rotational DOFs, before the sync so ghosts see it
         pinBoxTranslations();

         // synchronize -----------------------------------------------------------------
         tp["MESAPD: syncCall 2"].start();
         syncCall();
         tp["MESAPD: syncCall 2"].end();
      }

      // Debug hook (see SimulationObjects::freshlyActivatedFloaterUIDs): covers only the one main step right
      // after a release, so the set is cleared unconditionally every step.
      if (freshlyActivatedFloaterUIDs_ && !freshlyActivatedFloaterUIDs_->empty())
         freshlyActivatedFloaterUIDs_->clear();
   }

   // Holds selected TRANSLATIONAL DOFs of every pinPosition floater (sphere, box or cylinder), rotation free: for
   // each world axis with a nonzero pin mask, the position component is reset to the pinned value and the
   // velocity component zeroed after each sub-cycle; masked-out axes evolve normally. The floater still rotates
   // to its equilibrium tilt in place. Mask <1,1,1> = rigid hold; <1,1,0> pins x,y and lets z find its buoyant
   // depth.
   //
   // Independently, for pinRotation floaters the ANGULAR velocity about each world axis with a nonzero
   // pinRotationAxes component is zeroed every sub-cycle (orientation itself untouched). It is masked before the
   // next VelocityVerletPreForceUpdate integrates rotation from it, but that update also adds 0.5*wdot*dt^2 from
   // the old torque (not masked), so rotation about a pinned axis is suppressed, not exactly zero.
   // Mask <1,1,1> = no rotation; <1,0,1> allows rotation about Y only (roll in the XZ plane).
   void pinBoxTranslations()
   {
      if (pinnedBoxPositions_.empty() && pinnedRotationAxes_.empty()) return;
      for (size_t idx = 0; idx < particleAccessor->size(); ++idx)
      {
         const uint_t pt = particleAccessor->getType(idx);
         if (pt != uint_t(1) && pt != uint_t(3) && pt != uint_t(4)) continue; // pinned spheres, boxes, AND cylinders
         if (mesa_pd::data::particle_flags::isSet(
                particleAccessor->getFlags(idx), mesa_pd::data::particle_flags::GHOST)) continue;
         const auto uid = particleAccessor->getUid(idx);

         const auto it = pinnedBoxPositions_.find(uid);
         if (it != pinnedBoxPositions_.end())
         {
            const Vector3< real_t >& target = it->second;
            const Vector3< real_t >& mask   = pinnedBoxAxes_.at(uid);
            Vector3< real_t > pos = particleAccessor->getPosition(idx);
            Vector3< real_t > vel = particleAccessor->getLinearVelocity(idx);
            for (uint_t d = 0; d < uint_t(3); ++d)
               // pin this axis; free axes untouched
               if (std::abs(mask[d]) > real_t(0.5)) { pos[d] = target[d]; vel[d] = real_t(0); }
            particleAccessor->setPosition(idx, pos);
            particleAccessor->setLinearVelocity(idx, vel);
         }

         const auto rotIt = pinnedRotationAxes_.find(uid);
         if (rotIt != pinnedRotationAxes_.end())
         {
            const Vector3< real_t >& rotMask = rotIt->second;
            mesa_pd::Vec3 angVel = particleAccessor->getAngularVelocity(idx);
            for (uint_t d = 0; d < uint_t(3); ++d)
               if (std::abs(rotMask[d]) > real_t(0.5)) angVel[d] = real_t(0); // pin this axis; free axes untouched
            particleAccessor->setAngularVelocity(idx, angVel);
         }
      }
   }

   // Computes the submerged volume of every floating sphere (type 1) once per LBM step, by plane fitting with the
   // NormalSweep normal field (unit free-surface normal n̂ at every interface cell, pointing FROM LIQUID TO GAS,
   // i.e. upward for a flat surface).
   //
   // For each sphere, interface cells within a lateral disk of radius planeFitRadiusFactor_ * r around the sphere
   // axis are sampled (reaching past any crater/gas moat the sphere digs), keeping only upward normals (n_z > 0)
   // and dropping cells nearer another floater's surface. Weights φ(1−φ) peak at the transition (φ≈0.5) and
   // vanish in bulk liquid/gas:
   //
   //   wSum  = Σ wᵢ
   //   nSum  = Σ wᵢ nᵢ          (weighted-average normal direction)
   //   xSum  = Σ wᵢ xᵢ          (weighted centroid ≈ point on surface)
   //
   // After allReduce(SUM):
   //
   //   n̂     = nSum / |nSum|
   //   x̄     = xSum / wSum
   //   h_signed = n̂ · (x̄ − C)   (> 0 when centre C is on the liquid side, i.e. below surface)
   //   h_cap  = clamp(r + h_signed, 0, 2r)
   //   V_sub  = π · h_cap² · (3r − h_cap) / 3   (spherical-cap formula, exact for any plane)
   //
   // Exact for any free-surface orientation (waves). Uses only existing waLBerla fields (fillFieldID_,
   // flagFieldID_, normalFieldID_); no eigensolver.
   //
   // Fallback: if wSum ≈ 0 (no interface cells near the sphere) the sphere is fully in gas or fully in liquid;
   // the average fill φ̄ of the non-obstacle cells in the upper half of the sampling box decides:
   //   φ̄ ≤ 0.5 (gas)    → h_cap = 0
   //   φ̄ > 0.5 (liquid) → h_cap = 2r
   // (assuming "fully submerged" would launch a sphere starting in the gas phase).
   void computeSubmergedVolumes()
   {
      const size_t nF = floaterUIDs_.size();
      if (nF == 0) { submergedVolumes_.clear(); return; }

      // useFullVolumeAnalyticBuoyancy_: skip the plane-fit/surface-sampling entirely -- vSub=vFull
      // unconditionally for every floater (see SimulationParams.h's own doc comment on this flag).
      if (useFullVolumeAnalyticBuoyancy_)
      {
         submergedVolumes_.clear();
         for (size_t fi = 0; fi < nF; ++fi)
         {
            const real_t r     = floaterRadii_[fi];
            const real_t vFull = (real_t(4) / real_t(3)) * math::pi * r * r * r;
            submergedVolumes_[floaterUIDs_[fi]] = vFull;
         }
         return;
      }

      // Gather sphere centres from owner ranks; allReduce(SUM) broadcasts to all.
      std::vector< real_t > cenX(nF, real_t(0));
      std::vector< real_t > cenY(nF, real_t(0));
      std::vector< real_t > cenZ(nF, real_t(0));

      for (size_t idx = 0; idx < particleAccessor->size(); ++idx)
      {
         if (particleAccessor->getType(idx) != uint_t(1)) continue;
         if (mesa_pd::data::particle_flags::isSet(
                particleAccessor->getFlags(idx),
                mesa_pd::data::particle_flags::GHOST)) continue;
         const auto uid = particleAccessor->getUid(idx);
         for (size_t fi = 0; fi < nF; ++fi)
         {
            if (floaterUIDs_[fi] != uid) continue;
            const auto& pos = particleAccessor->getPosition(idx);
            cenX[fi] = pos[0]; cenY[fi] = pos[1]; cenZ[fi] = pos[2];
            break;
         }
      }
      mpi::allReduceInplace(cenX, mpi::SUM);
      mpi::allReduceInplace(cenY, mpi::SUM);
      mpi::allReduceInplace(cenZ, mpi::SUM);

      // Per-sphere accumulators packed into one flat vector for a single allReduce.
      // Layout per sphere fi:
      //   [0]   wSum      — Σ φ(1−φ)  over interface cells (plane-fitting weight)
      //   [1-3] nSum      — Σ wᵢ nᵢ   (weighted normal)
      //   [4-6] xSum      — Σ wᵢ xᵢ   (weighted centroid)
      //   [7]   fillSum   — Σ φ of non-obstacle cells above sphere equator (fallback)
      //   [8]   fillCount — number of such cells (fallback)
      //   [9]   aboveSum  — Σ φ of non-obstacle cells in the slab just ABOVE the sphere top (cavity probe)
      //   [10]  aboveCount— number of such cells (cavity probe)
      // Fallback: without interface cells (wSum≈0) all valid (gas/liquid) neighbours have the same fill (≈0 or
      // ≈1), so the average is correct even with several contributing ranks; obstacles (pillars, NoSlip walls)
      // are excluded by the fluidMask check.
      // Cavity probe: liquid on top of the sphere means the gas around its surface is an enclosed cavity, not the
      // real free surface, so the sphere is submerged (reference-free, no fixed liquid height).
      constexpr size_t stride = 11;
      std::vector< real_t > accum(nF * stride, real_t(0));

      // De-contamination of touching/clustered floaters: a sampled cell belongs to the floater whose SURFACE is
      // nearest. nearerOtherFloater(p, fi, myDistErr) is true if another floater's surface is closer to p than
      // myDistErr = |distance to fi's surface|, so a neighbour's waterline does not pollute fi's plane fit.
      // Minimum image in X only (assumes X is the only periodic axis; unlike floaterSurfaceDistance for boxes,
      // a periodic Y is not handled).
      const real_t domainXSize = blockForest_->getDomain().xSize();
      auto nearerOtherFloater = [&](const Vector3< real_t >& p, size_t fi, real_t myDistErr) -> bool {
         for (size_t fj = 0; fj < nF; ++fj)
         {
            if (fj == fi) continue;
            real_t ddx = p[0] - cenX[fj];
            ddx -= domainXSize * std::round(ddx / domainXSize); // periodic minimal image in X
            const real_t ddy = p[1] - cenY[fj];
            const real_t ddz = p[2] - cenZ[fj];
            const real_t distFj = std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
            if (std::abs(distFj - floaterRadii_[fj]) < myDistErr) return true;
         }
         return false;
      };

      for (auto& iblock : *blockForest_)
      {
         auto* fillField   = iblock.getData< ScalarField_T >(fillFieldID_);
         auto* normalField = iblock.getData< VectorField_T >(normalFieldID_);
         auto* flagField   = iblock.getData< FlagField_T >  (flagFieldID_);
         const CellInterval& blockBB = blockForest_->getBlockCellBB(iblock);
         const auto interfaceFlag =
            flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
         const auto gasFlag =
            flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
         const auto liquidFlag =
            flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
         const auto fluidMask = interfaceFlag | gasFlag | liquidFlag;

         for (size_t fi = 0; fi < nF; ++fi)
         {
            const real_t r         = floaterRadii_[fi];
            // lateral reach of the surface plane fit (prm floaterPlaneFitRadiusFactor, default 1.8), so it
            // reaches past any crater/gas moat the sphere digs
            const real_t planeLatR = planeFitRadiusFactor_ * r;
            const int    pad       = static_cast< int >(std::ceil(planeLatR)) + 2;

            Vector3< real_t > c(cenX[fi], cenY[fi], cenZ[fi]);
            blockForest_->mapToPeriodicDomain(c);

            const int cxi = static_cast< int >(c[0]);
            const int cyi = static_cast< int >(c[1]);
            const int czi = static_cast< int >(c[2]);

            for (int dz = -pad; dz <= pad; ++dz)
            for (int dy = -pad; dy <= pad; ++dy)
            for (int dx = -pad; dx <= pad; ++dx)
            {
               Vector3< real_t > samplePos(real_t(cxi + dx) + real_t(0.5),
                                           real_t(cyi + dy) + real_t(0.5),
                                           real_t(czi + dz) + real_t(0.5));
               blockForest_->mapToPeriodicDomain(samplePos);

               const Cell globalCell = blockForest_->getCell(samplePos);
               if (!blockBB.contains(globalCell)) continue;

               Cell localCell;
               blockForest_->transformGlobalToBlockLocalCell(localCell, iblock, globalCell);
               if (!fillField->xyzSize().contains(localCell)) continue;

               const auto cellFlag = flagField->get(localCell[0], localCell[1], localCell[2]);

               // Offsets of this cell from the (periodic-mapped) sphere centre (periodic-safe; no seam wrap).
               const real_t offX  = real_t(dx) + real_t(0.5) - (c[0] - real_t(cxi));
               const real_t offY  = real_t(dy) + real_t(0.5) - (c[1] - real_t(cyi));
               const real_t offZ  = real_t(dz) + real_t(0.5) - (c[2] - real_t(czi));
               const real_t distC = std::sqrt(offX * offX + offY * offY + offZ * offZ);

               if (isPartOfMaskSet(cellFlag, fluidMask))
               {
                  // Fallback channel: average fill of non-obstacle cells in the upper half (z ≥ centre).
                  if (samplePos[2] >= c[2])
                  {
                     accum[fi * stride + 7] += fillField->get(localCell[0], localCell[1], localCell[2]);
                     accum[fi * stride + 8] += real_t(1);
                  }
                  // Cavity probe: average fill in the slab just above the sphere top (offZ in (r, r+4], within
                  // the top disk). Liquid here means the gas around the sphere is an enclosed (numerical) cavity,
                  // so it is submerged; a normally floating sphere has its top in gas, so no false trigger.
                  const real_t lateral = std::sqrt(offX * offX + offY * offY);
                  if (offZ > r && offZ <= r + real_t(4) && lateral < r
                      && !nearerOtherFloater(samplePos, fi, offZ - r))
                  {
                     accum[fi * stride + 9]  += fillField->get(localCell[0], localCell[1], localCell[2]);
                     accum[fi * stride + 10] += real_t(1);
                  }
               }

               // Plane-fitting channel: fit the free surface over a wide lateral disk (radius planeLatR) around
               // the sphere axis, past any crater/gas moat, so the fitted waterline tracks the surrounding level.
               // Only interface cells with an upward normal (n_z > 0): crater/cavity walls point sideways or down.
               // (Can bias under strong waves, but correct for a calm floating body.)
               if (!isFlagSet(cellFlag, interfaceFlag)) continue;
               const real_t lateralP = std::sqrt(offX * offX + offY * offY);
               if (lateralP > planeLatR) continue;
               // De-contaminate: drop cells that are nearer a neighbouring floater's surface.
               if (nearerOtherFloater(samplePos, fi, std::abs(distC - r))) continue;

               const real_t phi = fillField->get(localCell[0], localCell[1], localCell[2]);
               const real_t w   = phi * (real_t(1) - phi);
               if (w < real_t(1e-6)) continue;

               const Vector3< real_t >& n = normalField->get(localCell[0], localCell[1], localCell[2]);
               if (n.sqrLength() < real_t(1e-10) || n[2] <= real_t(0)) continue; // upward free-surface only

               accum[fi * stride + 0] += w;
               accum[fi * stride + 1] += w * n[0];
               accum[fi * stride + 2] += w * n[1];
               accum[fi * stride + 3] += w * n[2];
               accum[fi * stride + 4] += w * samplePos[0];
               accum[fi * stride + 5] += w * samplePos[1];
               accum[fi * stride + 6] += w * samplePos[2];
            }
         }
      }

      mpi::allReduceInplace(accum, mpi::SUM);

      submergedVolumes_.clear();
      // Dynamic MEM reference height (ModelParameters.useDynamicReferenceHeight): accumulate a wTot-weighted
      // average of localSurfaceZ over floaters that have a clean fitted ring THIS step (wTot>1e-6). Floaters
      // fully submerged/emerged (no ring) don't contribute — the reference is held, not snapped to a
      // meaningless centre-Z. EMA-blended into memReferenceHeight_ after the loop.
      real_t dynHWeightedSum = real_t(0);
      real_t dynHWeightSum   = real_t(0);
      for (size_t fi = 0; fi < nF; ++fi)
      {
         const real_t r     = floaterRadii_[fi];
         const real_t vFull = (real_t(4) / real_t(3)) * math::pi * r * r * r;

         const Vector3< real_t > C(cenX[fi], cenY[fi], cenZ[fi]);
         const real_t wTot    = accum[fi * stride + 0];
         const real_t fillCnt = accum[fi * stride + 8];
         const real_t avgFill = (fillCnt > real_t(0)) ? (accum[fi * stride + 7] / fillCnt) : real_t(0);
         real_t h_cap;
         real_t localSurfaceZ = C[2]; // reported below; refined to the fitted waterline when one exists
         Vector3< real_t > fittedNormal(real_t(0)); // free-surface normal used for the plane-fit (0 = none)
         // |Σw n̂|/Σw ∈ [0,1]; low ⇒ normals cancel (cavity / neighbour contamination)
         real_t normalCoherence = real_t(0);

         // --- Reference-free submerged depth: the waterline is read locally around THIS sphere ---
         //   • plane fit available (wTot>0): h_signed = n̂·(x̄−C) is the signed depth of the centre below the
         //     fitted plane; the clamp covers the just-submerged / just-emerged ends.
         //   • otherwise the sphere is fully in one phase: the average fill φ̄ of the surrounding non-obstacle
         //     cells decides (φ̄>0.5 liquid -> 2r, else gas -> 0).
         // No flat global surface is assumed, so this holds on deformed surfaces (piles, pillars, wakes).
         // Limitation: a gas pocket hugging a deeply submerged sphere can under-report V_sub.
         if (wTot > real_t(1e-6))
         {
            const Vector3< real_t > nSum(accum[fi*stride+1], accum[fi*stride+2], accum[fi*stride+3]);
            const Vector3< real_t > xBar(accum[fi*stride+4] / wTot,
                                         accum[fi*stride+5] / wTot,
                                         accum[fi*stride+6] / wTot);
            localSurfaceZ   = xBar[2];
            normalCoherence = nSum.length() / wTot; // 1 ⇒ all ring normals aligned; →0 ⇒ they cancel

            // Feed the fitted waterline into the dynamic-H average (ambient surface: the fit reaches out to
            // planeFitRadiusFactor_ * r, past the floater's own crater).
            dynHWeightedSum += wTot * localSurfaceZ;
            dynHWeightSum   += wTot;

            if (nSum.sqrLength() < real_t(1e-10))
            {
               // Weighted normals cancelled (e.g. an enclosing cavity ring) -> no clean waterline
               // direction; decide by the surrounding fill instead.
               h_cap = (avgFill > real_t(0.5)) ? real_t(2) * r : real_t(0);
            }
            else
            {
               // n̂ points liquid->gas; h_signed > 0 when centre C is below the surface.
               const Vector3< real_t > nHat     = nSum / nSum.length();
               fittedNormal = nHat;
               const real_t            h_signed = nHat * (xBar - C);
               h_cap = std::clamp(r + h_signed, real_t(0), real_t(2) * r);
            }
         }
         else
         {
            // No contact ring -> the sphere is fully in one phase -> the surrounding fill decides.
            h_cap = (avgFill > real_t(0.5)) ? real_t(2) * r : real_t(0);
         }

         // CAVITY BLEND (reference-free, smooth): blend h_cap toward fully submerged (2r) as the slab just above
         // the sphere top fills with liquid. avgFillAbove is the mean fill of those "ceiling" cells (gas→0,
         // liquid→1, interface→partial; obstacles excluded); water overhead ⇒ enclosed cavity ⇒ submerged. A
         // smoothstep over [0.3, 0.7] keeps the buoyancy force continuous (a hard switch excites a growing heave
         // oscillation that breaks the free-surface solver).
         const real_t aboveCnt     = accum[fi * stride + 10];           // # ceiling cells checked
         const real_t avgFillAbove = (aboveCnt > real_t(0)) ? (accum[fi * stride + 9] / aboveCnt) : real_t(0);
         const real_t tCav = std::clamp((avgFillAbove - real_t(0.3)) / real_t(0.4), real_t(0), real_t(1));
         const real_t wCav = tCav * tCav * (real_t(3) - real_t(2) * tCav); // smoothstep weight ∈ [0,1]
         h_cap = (real_t(1) - wCav) * h_cap + wCav * (real_t(2) * r);

         const real_t vSub = sphereCapVolume(r, h_cap);

         submergedVolumes_[floaterUIDs_[fi]] = vSub;

         if (stepCounter_ % evalFrequency_ == uint_t(0))
         {
            // Net VERTICAL force from buoyancy + gravity (the same expression AddExternalForces applies),
            // EXCLUDING the hydrodynamic/MEM force: (rho*vFull - vSub) * g_z. >0 = up (sphere rises),
            // <0 = down (sinks). rho is the floater's own density ratio; g_z = acceleration_[2] (<0).
            const real_t rhoP      = (fi < floaterDensities_.size()) ? floaterDensities_[fi] : real_t(1);
            const real_t netForceZ = (rhoP * vFull - vSub) * (*acceleration_)[2];
            // Depth of the sphere centre below the LIVE local waterline (>0 = centre submerged).
            // At equilibrium this is ~0 for a half-density sphere even if centreZ sits far above the
            // undisturbed level absoluteLiquidHeight, because the sphere rides on a raised/dammed surface.
            const real_t centreDepth = localSurfaceZ - C[2];
            WALBERLA_LOG_INFO_ON_ROOT("[FloaterBuoyancy] uid=" << floaterUIDs_[fi]
               << "  h_cap=" << h_cap << " / " << real_t(2)*r << " LU"
               << "  V_sub=" << vSub << " / " << vFull << " LU³"
               << "  submersion=" << real_t(100) * vSub / vFull << "%"
               << "  netF_z(buoy+grav,excl.hydro)=" << netForceZ << " (>0 up, <0 down)"
               << "  centreZ=" << C[2]
               << "  localSurfaceZ=" << localSurfaceZ
               << (wTot > real_t(1e-6) ? " (plane-fit)" : " (no ring)")
               << "  avgFill=" << avgFill
               << "  avgFillAbove=" << avgFillAbove << " (cellsAbove=" << aboveCnt << ", cavityBlend=" << wCav << ")"
               << "  fitNormal=(" << fittedNormal[0] << "," << fittedNormal[1] << "," << fittedNormal[2] << ")"
               << (fittedNormal[2] < real_t(0) ? " (n_z<0 DOWN!)" : "")
               << "  normalCoherence=" << normalCoherence
               << "  centreDepthBelowSurface=" << centreDepth << " LU"
               << "  (oldFixedRef=" << surfaceHeight_ << ", now unused for force)"
               << "  (expect submersion ~densityRatio×100% and centreDepth~0 at equilibrium)");

         }
      }

      // Dynamic MEM reference height: EMA-blend the wTot-weighted ambient waterline average into
      // memReferenceHeight_. Force-side only (SimpleBB's rhoHydStat feeds only forceMEM, never the fluid; see
      // AntidunesBoundaryHandling / SimpleBB.h), so no PDF/reconstruction feedback. Held (not snapped) when no
      // floater has a plane fit this step. A static/offset correction only: it cancels the uncancelled pile-up
      // head but cannot remove a velocity-dependent (dynamic suction) residual.
      if (useDynamicReferenceHeight_ && memReferenceHeight_ && dynHWeightSum > real_t(1e-6))
      {
         constexpr real_t alpha = real_t(0.2);
         const real_t target = dynHWeightedSum / dynHWeightSum;
         *memReferenceHeight_ = (real_t(1) - alpha) * (*memReferenceHeight_) + alpha * target;

         if (stepCounter_ % evalFrequency_ == uint_t(0))
         {
            WALBERLA_LOG_INFO_ON_ROOT("[DynamicH] target(ambient waterline avg)=" << target
               << "  H(EMA)=" << *memReferenceHeight_
               << "  (flatRef=" << surfaceHeight_ << ", nFloatersContributing weight=" << dynHWeightSum << ")");
         }
      }
   }

   // Computes the submerged volume and centre-of-buoyancy offset of every floating box (type 3) and cylinder
   // (type 4) once per LBM step: the local free surface is plane-fitted from interface cells in a shell around
   // each body (as for the spheres), then the rotated body is clipped against that plane (12³-point sampling,
   // used for the force; the analytic clip is computed as a cross-check). Without interface cells nearby, the
   // average shell fill decides fully submerged vs. fully emerged. Followed by the cavity blend.
   void computeSubmergedVolumesBox()
   {
      buoyancyCenterOffset_.clear(); // box/cylinder-only; recomputed every step
      const size_t nB = clipFloaters_.size(); // boxes (type 3) + cylinders (type 4)
      if (nB == 0) return;

      // useFullVolumeAnalyticBuoyancy_: skip the plane fit/sampling; vSub=vFull, no righting lever, for every
      // floater (see SimulationParams.h).
      if (useFullVolumeAnalyticBuoyancy_)
      {
         for (size_t bi = 0; bi < nB; ++bi)
         {
            const real_t hx    = clipFloaters_[bi].halfExtents[0];
            const real_t hy    = clipFloaters_[bi].halfExtents[1];
            const real_t hz    = clipFloaters_[bi].halfExtents[2];
            const bool   isCyl = clipFloaters_[bi].isCylinder;
            const real_t vFull = isCyl ? math::pi * hx * hx * real_t(2) * hz
                                       : real_t(8) * hx * hy * hz;
            submergedVolumes_[clipFloaters_[bi].uid] = vFull;
            // buoyancyCenterOffset_ left empty (cleared above) -- symmetric, no righting lever.
         }
         return;
      }

      // Gather centre positions and rotation matrices from owner ranks.
      // Non-owner ranks contribute 0; allReduce(SUM) distributes the owner's values.
      std::vector< real_t > cenX(nB, real_t(0));
      std::vector< real_t > cenY(nB, real_t(0));
      std::vector< real_t > cenZ(nB, real_t(0));
      std::vector< real_t > rotFlat(nB * 9, real_t(0)); // row-major R (body→world)

      for (size_t idx = 0; idx < particleAccessor->size(); ++idx)
      {
         const uint_t pt = particleAccessor->getType(idx);
         if (pt != uint_t(3) && pt != uint_t(4)) continue;
         if (mesa_pd::data::particle_flags::isSet(
                particleAccessor->getFlags(idx),
                mesa_pd::data::particle_flags::GHOST)) continue;
         const auto uid = particleAccessor->getUid(idx);
         for (size_t bi = 0; bi < nB; ++bi)
         {
            if (clipFloaters_[bi].uid != uid) continue;
            const auto& pos             = particleAccessor->getPosition(idx);
            cenX[bi] = pos[0]; cenY[bi] = pos[1]; cenZ[bi] = pos[2];
            const mesa_pd::Mat3& R      = particleAccessor->getRotation(idx).getMatrix();
            for (uint_t r = 0; r < 3; ++r)
               for (uint_t c = 0; c < 3; ++c)
                  rotFlat[bi * 9 + r * 3 + c] = R(r, c);
            break;
         }
      }
      mpi::allReduceInplace(cenX,    mpi::SUM);
      mpi::allReduceInplace(cenY,    mpi::SUM);
      mpi::allReduceInplace(cenZ,    mpi::SUM);
      mpi::allReduceInplace(rotFlat, mpi::SUM);

      // --- Plane-fit of the local free surface from interface cells in a SHELL around each box ---
      // Independent of the box-interior fill (which becomes obstacle once the box is mapped as a moving
      // obstacle). Mirrors the floating-sphere method: weighted-average free-surface normal n̂ (from the
      // normal field) and weighted centroid x̄ over nearby interface cells; weight φ(1−φ) peaks at the
      // φ≈0.5 waterline. Per box bi in `acc` (stride 10): [0]=wSum [1..3]=Σw·n̂ [4..6]=Σw·x [7]=Σφ (shell)
      // [8]=aboveSum [9]=aboveCount (ceiling/cavity probe -- see below, ported from the sphere method).
      constexpr size_t stride = 10;
      std::vector< real_t > acc(nB * stride, real_t(0));
      std::vector< real_t > fluidCnt(nB, real_t(0)); // # non-obstacle shell cells (fallback denominator)

      // De-contamination of touching/nearby box/cylinder floaters, as nearerOtherFloater() for spheres: with a
      // wide boxPlaneFitBandFactor a shell can reach far past the body and overlap a neighbour's, dragging every
      // fitted plane toward the group's shared centroid.
      // floaterSurfaceDistance(p, idx) is the true distance from world point p to floater idx's ORIENTED surface
      // (world->body via Rᵀ, then clamp to the shape's extent: per axis for a box, radial+axial for a cylinder).
      // A bounding-sphere approximation would be poor for elongated shapes (a radius-8/half-length-12 cylinder
      // reaches 8 LU laterally but has a bounding radius of ~14.4 LU).
      const real_t domainXSize = blockForest_->getDomain().xSize();
      const real_t domainYSize = blockForest_->getDomain().ySize();
      const bool   xPeriodic   = blockForest_->isXPeriodic();
      const bool   yPeriodic   = blockForest_->isYPeriodic();
      auto floaterSurfaceDistance = [&](const Vector3< real_t >& p, size_t idx) -> real_t {
         const auto& cf = clipFloaters_[idx];
         mesa_pd::Mat3 R;
         for (uint_t r = 0; r < 3; ++r)
            for (uint_t c = 0; c < 3; ++c)
               R(r, c) = rotFlat[idx * 9 + r * 3 + c];

         real_t ddx = p[0] - cenX[idx];
         if (xPeriodic) ddx -= domainXSize * std::round(ddx / domainXSize);
         real_t ddy = p[1] - cenY[idx];
         if (yPeriodic) ddy -= domainYSize * std::round(ddy / domainYSize);
         const real_t ddz = p[2] - cenZ[idx];

         const mesa_pd::Vec3 local = R.getTranspose() * mesa_pd::Vec3(ddx, ddy, ddz); // world -> body frame
         const Vector3< real_t >& h = cf.halfExtents;

         mesa_pd::Vec3 clamped;
         if (cf.isCylinder)
         {
            // halfExtents = (r, r, halfLength); clamp radially in xy, axially in z.
            const real_t rho   = std::sqrt(local[0] * local[0] + local[1] * local[1]);
            const real_t scale = (rho > real_t(1e-10)) ? std::min(real_t(1), h[0] / rho) : real_t(0);
            clamped = mesa_pd::Vec3(local[0] * scale, local[1] * scale, std::clamp(local[2], -h[2], h[2]));
         }
         else
         {
            clamped = mesa_pd::Vec3(std::clamp(local[0], -h[0], h[0]), std::clamp(local[1], -h[1], h[1]),
                                    std::clamp(local[2], -h[2], h[2]));
         }
         // 0 if p projects inside the shape, else the distance to the nearest surface point
         return (local - clamped).length();
      };
      auto nearerOtherBoxFloater = [&](const Vector3< real_t >& p, size_t bi, real_t myDistErr) -> bool {
         for (size_t bj = 0; bj < nB; ++bj)
         {
            if (bj == bi) continue;
            if (floaterSurfaceDistance(p, bj) < myDistErr) return true;
         }
         return false;
      };

      for (auto& iblock : *blockForest_)
      {
         auto* fillField   = iblock.getData< ScalarField_T >(fillFieldID_);
         auto* normalField = iblock.getData< VectorField_T >(normalFieldID_);
         auto* flagField   = iblock.getData< FlagField_T >  (flagFieldID_);
         const CellInterval& blockBB = blockForest_->getBlockCellBB(iblock);
         const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
         const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
         const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
         const auto fluidMask     = interfaceFlag | gasFlag | liquidFlag;

         for (size_t bi = 0; bi < nB; ++bi)
         {
            const real_t hx = clipFloaters_[bi].halfExtents[0];
            const real_t hy = clipFloaters_[bi].halfExtents[1];
            const real_t hz = clipFloaters_[bi].halfExtents[2];
            // Wide lateral reach so the plane fit samples past any crater/gas moat the box digs, anchoring the
            // fitted waterline to the surrounding surface (as for the sphere). prm
            // ModelParameters.boxPlaneFitBandFactor (default 0.8).
            const real_t band = boxPlaneFitBandFactor_ * std::max({hx, hy, hz}) + real_t(2);

            mesa_pd::Mat3 R;
            for (uint_t r = 0; r < 3; ++r)
               for (uint_t c = 0; c < 3; ++c)
                  R(r, c) = rotFlat[bi * 9 + r * 3 + c];

            Vector3< real_t > pos(cenX[bi], cenY[bi], cenZ[bi]);
            blockForest_->mapToPeriodicDomain(pos);

            // World-space AABB half-extents of the rotated box, expanded by the shell `band`.
            const Vector3< real_t > axx = R * Vector3< real_t >(hx, real_t(0), real_t(0));
            const Vector3< real_t > ayy = R * Vector3< real_t >(real_t(0), hy, real_t(0));
            const Vector3< real_t > azz = R * Vector3< real_t >(real_t(0), real_t(0), hz);
            const real_t AABB_hx = std::abs(axx[0]) + std::abs(ayy[0]) + std::abs(azz[0]) + band;
            const real_t AABB_hy = std::abs(axx[1]) + std::abs(ayy[1]) + std::abs(azz[1]) + band;
            const real_t AABB_hz = std::abs(axx[2]) + std::abs(ayy[2]) + std::abs(azz[2]) + band;

            // Ceiling/cavity probe setup (as in computeSubmergedVolumes()): world-space top of the shape's true
            // (unpadded) extent and an isotropic lateral cutoff for "directly above" (the sphere's `lateral < r`).
            // Liquid there means the gas pocket around the shape is an enclosed cavity under real liquid; used by
            // the cavity blend below.
            const real_t trueTopZ     = pos[2] + std::abs(axx[2]) + std::abs(ayy[2]) + std::abs(azz[2]);
            const real_t ceilingReach = std::max({ hx, hy, hz });

            const int ixMin = static_cast< int >(std::floor(pos[0] - AABB_hx));
            const int ixMax = static_cast< int >(std::ceil (pos[0] + AABB_hx));
            const int iyMin = static_cast< int >(std::floor(pos[1] - AABB_hy));
            const int iyMax = static_cast< int >(std::ceil (pos[1] + AABB_hy));
            const int izMin = static_cast< int >(std::floor(pos[2] - AABB_hz));
            const int izMax = static_cast< int >(std::ceil (pos[2] + AABB_hz));

            for (int iz = izMin; iz <= izMax; ++iz)
            for (int iy = iyMin; iy <= iyMax; ++iy)
            for (int ix = ixMin; ix <= ixMax; ++ix)
            {
               Vector3< real_t > samplePos(real_t(ix) + real_t(0.5), real_t(iy) + real_t(0.5),
                                           real_t(iz) + real_t(0.5));
               blockForest_->mapToPeriodicDomain(samplePos);

               const Cell globalCell = blockForest_->getCell(samplePos);
               if (!blockBB.contains(globalCell)) continue;

               Cell localCell;
               blockForest_->transformGlobalToBlockLocalCell(localCell, iblock, globalCell);
               if (!fillField->xyzSize().contains(localCell)) continue;

               const auto cellFlag = flagField->get(localCell[0], localCell[1], localCell[2]);

               // Fallback channel: fill of any non-obstacle shell cell (≈0 fully-gas / ≈1 fully-liquid box).
               if (isPartOfMaskSet(cellFlag, fluidMask))
               {
                  acc[bi * stride + 7] += fillField->get(localCell[0], localCell[1], localCell[2]);
                  fluidCnt[bi]         += real_t(1);

                  // Ceiling/cavity probe (as in computeSubmergedVolumes()): average fill of non-obstacle cells in
                  // the slab just above the shape's true top, within its lateral footprint. Liquid here means the
                  // shape is submerged regardless of what the (possibly void-contaminated) plane fit found.
                  const real_t worldOffZ = samplePos[2] - trueTopZ;
                  const real_t lateralXY = std::sqrt((samplePos[0] - pos[0]) * (samplePos[0] - pos[0]) +
                                                      (samplePos[1] - pos[1]) * (samplePos[1] - pos[1]));
                  if (worldOffZ > real_t(0) && worldOffZ <= real_t(4) && lateralXY < ceilingReach
                      && !nearerOtherBoxFloater(samplePos, bi, worldOffZ))
                  {
                     acc[bi * stride + 8] += fillField->get(localCell[0], localCell[1], localCell[2]);
                     acc[bi * stride + 9] += real_t(1);
                  }
               }

               // Plane-fitting channel: only interface cells carry a valid normal.
               if (!isFlagSet(cellFlag, interfaceFlag)) continue;

               // De-contaminate: drop cells nearer a neighbouring box/cylinder floater's own TRUE
               // (oriented) surface than they are to bi's own — see floaterSurfaceDistance() above.
               if (nearerOtherBoxFloater(samplePos, bi, floaterSurfaceDistance(samplePos, bi)))
                  continue;

               const real_t phi = fillField->get(localCell[0], localCell[1], localCell[2]);
               const real_t w   = phi * (real_t(1) - phi);
               if (w < real_t(1e-6)) continue;

               const Vector3< real_t >& n = normalField->get(localCell[0], localCell[1], localCell[2]);
               if (n.sqrLength() < real_t(1e-10) || n[2] <= real_t(0)) continue; // upward free-surface only

               acc[bi * stride + 0] += w;
               acc[bi * stride + 1] += w * n[0];
               acc[bi * stride + 2] += w * n[1];
               acc[bi * stride + 3] += w * n[2];
               acc[bi * stride + 4] += w * samplePos[0];
               acc[bi * stride + 5] += w * samplePos[1];
               acc[bi * stride + 6] += w * samplePos[2];
            }
         }
      }

      mpi::allReduceInplace(acc,      mpi::SUM);
      mpi::allReduceInplace(fluidCnt, mpi::SUM);

      // --- V_sub and centre of buoyancy: clip the rotated box against the fitted free-surface plane ---
      // Sample a regular grid of points inside the box (body frame), rotate to world, and count those
      // BELOW the surface (n̂·(p − x̄) < 0, since n̂ points liquid→gas). V_sub = vFull·(#below/#total); the
      // mean offset of the submerged samples from the box centre is the lever arm (B − G) for the torque.
      // Computed identically on every rank (plane + pose are global after the allReduce), so no further
      // communication is needed.
      for (size_t bi = 0; bi < nB; ++bi)
      {
         const real_t hx    = clipFloaters_[bi].halfExtents[0];
         const real_t hy    = clipFloaters_[bi].halfExtents[1];
         const real_t hz    = clipFloaters_[bi].halfExtents[2];
         const bool   isCyl = clipFloaters_[bi].isCylinder; // body-z axis, radius hx (= hy), halfLength hz
         const real_t vFull = isCyl ? math::pi * hx * hx * real_t(2) * hz
                                    : real_t(8) * hx * hy * hz;

         mesa_pd::Mat3 R;
         for (uint_t r = 0; r < 3; ++r)
            for (uint_t c = 0; c < 3; ++c)
               R(r, c) = rotFlat[bi * 9 + r * 3 + c];
         Vector3< real_t > pos(cenX[bi], cenY[bi], cenZ[bi]);
         blockForest_->mapToPeriodicDomain(pos);

         const real_t wTot = acc[bi * stride + 0];

         // Live local waterline at the box (REFERENCE-FREE, for the log): gradient-weighted mean z of the
         // shell interface cells. The plane-fit + clip below uses the FULL fitted plane (level AND tilt);
         // the clip naturally covers fully-submerged (all box points below the plane → vFull) and
         // fully-emerged (all above → 0), so no fixed absoluteLiquidHeight bound is needed.
         const real_t localSurfaceZ =
            (wTot > real_t(1e-6)) ? (acc[bi * stride + 6] / wTot) : pos[2];
         Vector3< real_t > fittedNormal(real_t(0)); // fitted free-surface normal (for the log)
         real_t            normalCoherence = real_t(0);

         real_t            vSub = real_t(0);
         Vector3< real_t > cbOff(real_t(0));
         real_t            vSubSample   = real_t(-1); // sampling result   (−1 = not computed this step)
         real_t            vSubAnalytic = real_t(-1); // analytical result (−1 = not computed this step)

         if (wTot > real_t(1e-6))
         {
            // Plane-fit the surrounding surface (wide shell, upward normals only), then clip the rotated
            // box against that plane. n̂ points liquid→gas; degenerate (cancelled) normals fall back to
            // vertical (horizontal surface).
            const Vector3< real_t > nSum(acc[bi*stride+1], acc[bi*stride+2], acc[bi*stride+3]);
            normalCoherence = nSum.length() / wTot;
            const Vector3< real_t > nHat = (nSum.sqrLength() > real_t(1e-10))
                                              ? nSum / nSum.length()
                                              : Vector3< real_t >(real_t(0), real_t(0), real_t(1));
            fittedNormal = nHat;
            const Vector3< real_t > xBar(acc[bi*stride+4] / wTot,
                                         acc[bi*stride+5] / wTot,
                                         acc[bi*stride+6] / wTot);

            // ---- Analytical clip (exact OBB∩plane for boxes, exact disk∩plane for cylinders; cheap), logged as a
            //      cross-check of the sampling estimate. ----
            Vector3< real_t > cbAnalytic(real_t(0));
            if (isCyl) clipCylinderAnalytic(hx, hz, R, nHat, xBar, pos, vSubAnalytic, cbAnalytic);
            else       clipBoxAnalytic(Vector3< real_t >(hx, hy, hz), R, nHat, xBar, pos, vSubAnalytic, cbAnalytic);

            // ---- Sampling clip (12³ points), the V_sub used for the force. For cylinders,
            //      grid points outside the circular cross-section (body xy) are rejected, so the kept
            //      samples uniformly fill the cylinder volume. ----
            Vector3< real_t > cbSample(real_t(0));
            {
               const real_t planeOffset = nHat * (xBar - pos); // submerged sample: nHat·off < planeOffset
               constexpr int N = 12; // 12³ = 1728 samples per box (≈ π/4 of them inside a cylinder)
               uint_t            below = uint_t(0);
               uint_t            total = uint_t(0);
               Vector3< real_t > cbSum(real_t(0));
               for (int i = 0; i < N; ++i)
               for (int j = 0; j < N; ++j)
               for (int k = 0; k < N; ++k)
               {
                  const Vector3< real_t > local(-hx + (real_t(i) + real_t(0.5)) * (real_t(2) * hx / real_t(N)),
                                                -hy + (real_t(j) + real_t(0.5)) * (real_t(2) * hy / real_t(N)),
                                                -hz + (real_t(k) + real_t(0.5)) * (real_t(2) * hz / real_t(N)));
                  if (isCyl && (local[0] * local[0] + local[1] * local[1] > hx * hx)) continue;
                  const Vector3< real_t > off = R * local; // offset from box centre (world frame)
                  ++total;
                  if ((nHat * off) < planeOffset) { ++below; cbSum += off; }
               }
               vSubSample = vFull * real_t(below) / real_t(total);
               if (below > uint_t(0)) cbSample = cbSum / real_t(below);
            }

            vSub  = vSubSample;
            cbOff = cbSample;
         }
         else
         {
            // No interface cells near the box: it is entirely in gas or entirely in liquid.
            const real_t cnt     = fluidCnt[bi];
            const real_t avgFill = (cnt > real_t(0)) ? acc[bi * stride + 7] / cnt : real_t(0);
            vSub  = (avgFill > real_t(0.5)) ? vFull : real_t(0);
            cbOff = Vector3< real_t >(real_t(0)); // symmetric → no righting lever
         }

         // CAVITY BLEND (as in computeSubmergedVolumes()): blend vSub toward vFull as the slab just above the
         // shape's top fills with liquid. Deep in a packed pile the plane fit can lock onto a small local void
         // between particles instead of the real (much higher) free surface, under-reporting vSub even for a
         // shape under bulk liquid.
         const real_t aboveCnt     = acc[bi * stride + 9];
         const real_t avgFillAbove = (aboveCnt > real_t(0)) ? (acc[bi * stride + 8] / aboveCnt) : real_t(0);
         const real_t tCav = std::clamp((avgFillAbove - real_t(0.3)) / real_t(0.4), real_t(0), real_t(1));
         const real_t wCav = tCav * tCav * (real_t(3) - real_t(2) * tCav); // smoothstep weight ∈ [0,1]
         vSub = (real_t(1) - wCav) * vSub + wCav * vFull;
         // A fully submerged centrally symmetric shape has its submerged centroid at its centre (cbOff -> 0):
         // blend the lever arm toward that limit with the same weight, so the righting torque does not keep a
         // partial-submersion lever once vSub is corrected toward full submersion.
         cbOff = (real_t(1) - wCav) * cbOff;

         vSub = std::clamp(vSub, real_t(0), vFull);
         submergedVolumes_[clipFloaters_[bi].uid] = vSub;
         if (vSub > real_t(1e-9))
            buoyancyCenterOffset_[clipFloaters_[bi].uid] = cbOff;

         if (stepCounter_ % evalFrequency_ == uint_t(0))
         {
            // Y-tilt angle from rotation matrix: for R = Ry(θ), R(0,2)=sinθ, R(0,0)=cosθ.
            const real_t R00     = rotFlat[bi * 9 + 0];
            const real_t R02     = rotFlat[bi * 9 + 2];
            const real_t tiltDeg = std::atan2(R02, R00) * real_t(180) / math::pi;
            WALBERLA_LOG_INFO_ON_ROOT("[BoxBuoyancy] uid=" << clipFloaters_[bi].uid
               << "  shape=" << (isCyl ? "CYLINDER" : "BOX")
               << "  submersion=" << real_t(100) * vSub / vFull << "%"
               << "  tiltY=" << tiltDeg << "°"
               << "  V_sub=" << vSub << " / " << vFull << " LU³"
               << "  centreZ=" << pos[2]
               << "  localSurfaceZ=" << localSurfaceZ << (wTot > real_t(1e-6) ? " (plane-fit)" : " (no surf)")
               << "  fitNormal=(" << fittedNormal[0] << "," << fittedNormal[1] << "," << fittedNormal[2] << ")"
               << (fittedNormal[2] < real_t(0) ? " (n_z<0 DOWN!)" : "")
               << "  normalCoherence=" << normalCoherence
               << "  avgFillAbove=" << avgFillAbove << " (cellsAbove=" << aboveCnt << ", cavityBlend=" << wCav << ")"
               << "  B-G offset=(" << cbOff[0] << "," << cbOff[1] << "," << cbOff[2] << ")"
               << "  |B-G|_horiz=" << std::sqrt(cbOff[0] * cbOff[0] + cbOff[1] * cbOff[1])
               << "  comOffsetZ=" << (boxComOffsetZ_.count(clipFloaters_[bi].uid)
                                         ? boxComOffsetZ_.at(clipFloaters_[bi].uid) : real_t(0))
               << "  method=SAMPLING"
               // Cross-check both estimators when both were computed this step (both ≥ 0).
               << "  [xcheck V_sub: analytic=" << vSubAnalytic << " sampling=" << vSubSample
               << ((vSubAnalytic >= real_t(0) && vSubSample >= real_t(0) && vFull > real_t(0))
                     ? std::string("  diff=") + std::to_string(real_t(100) * (vSubAnalytic - vSubSample) / vFull) + "%"
                     : std::string(""))
               << "]");

            // --- RIGHTING MOMENT (Bogner & Rüde 2013, Fig. 5.6 validation) ------------------------------
            // Stability curve of a half-immersed box heeled about its longitudinal (y) axis: the buoyancy
            // F_buoy = +V_sub·|g| acts at the centre of buoyancy B = centroid of the submerged volume; the
            // horizontal offset of B from the centre of gravity G (= box centre) gives the restoring moment
            //   M = (B − G) × F_buoy,   heeled about y ⇒ M_y = −(B−G)_x · V_sub·|g|  (>0 restoring here).
            // For a box held at heel α this is one point of the stability curve M(α), for comparison with the
            // ideal analytic curve. B − G is cbOff (from the sampling clip, after the cavity blend). Logs |M|, M_y,
            // the heel angle, the righting arm GZ = |M|/(V_sub·|g|) and V_sub, once per eval step.
            const real_t gMag     = std::abs((*acceleration_)[2]);
            const real_t fBuoyZ   = vSub * gMag;                              // upward buoyant force
            const Vector3< real_t > Mvec(cbOff[1] * fBuoyZ,                   // (B−G) × (0,0,F_buoy)
                                         -cbOff[0] * fBuoyZ,
                                         real_t(0));
            const real_t GZ       = (fBuoyZ > real_t(1e-12)) ? Mvec.length() / fBuoyZ : real_t(0);
            WALBERLA_LOG_INFO_ON_ROOT("[BoxRightingMoment] uid=" << clipFloaters_[bi].uid
               << "  heelDeg=" << tiltDeg
               << "  M_y=" << Mvec[1] << "  |M|=" << Mvec.length() << " (>0 restoring)"
               << "  GZ=" << GZ << " LU  (righting arm)"
               << "  V_sub=" << vSub << "  submersion=" << real_t(100) * vSub / vFull << "%"
               << "  |g|=" << gMag << "  method=SAMPLING");
         }
      }
   }

   // Fluid density rho (the quantity the [LOW-RHO] guard in ReconstructInterfaceCellABB.h checks before u/=rho)
   // at the LBM cell containing the global position, by scanning this rank's blocks (one block per process here,
   // so O(1) in practice). Returns NaN if pos lies in no block of this rank (e.g. a ghost floater owned
   // elsewhere); callers must not treat that as a reading.
   real_t sampleFluidDensityAt(Vector3< real_t > pos) const
   {
      blockForest_->mapToPeriodicDomain(pos);
      const Cell globalCell = blockForest_->getCell(pos);

      for (auto& iblock : *blockForest_)
      {
         const CellInterval& blockBB = blockForest_->getBlockCellBB(iblock);
         if (!blockBB.contains(globalCell)) continue;

         Cell localCell;
         blockForest_->transformGlobalToBlockLocalCell(localCell, iblock, globalCell);
         const auto* pdfField = iblock.getData< const PdfField_T >(pdfFieldID_);
         if (!pdfField->xyzSize().contains(localCell)) continue;

         return pdfField->getDensity(localCell);
      }
      return std::numeric_limits< real_t >::quiet_NaN();
   }

private:
   std::shared_ptr< ParticleAccessor_T >              particleAccessor;
   std::shared_ptr< mesa_pd::data::ParticleStorage >  particleStorage;
   std::function< void() >                            syncCall;
   bool useOpenMP{ false };

   shared_ptr< StructuredBlockForest >                blockForest_;
   BlockDataID                                        fillFieldID_;
   BlockDataID                                        flagFieldID_;
   BlockDataID                                        normalFieldID_;
   BlockDataID                                        pdfFieldID_;
   BlockDataID                                        particleFieldID_;
   BlockDataID                                        boundaryHandlingID_;
   // false in the fluid-free DropTest (no boundary handling); currently not read anywhere (nothing is gated on it)
   bool                                               hasFluid_;
   // Shared buffer (see SimulationObjects::DiagnosticEvent) that per-sub-cycle checks push records into (when
   // diagnosticSweepEnabled_) instead of logging directly; drained once per main step by BuildTimeloop.h's
   // addAfterTimeStepDiagnosticSweep.
   bool                                                diagnosticSweepEnabled_;
   std::shared_ptr< std::vector< SimulationObjects::DiagnosticEvent > > diagnosticEvents_;
   std::vector< walberla::id_t >                      floaterUIDs_;         // sphere floaters, same size on all ranks
   std::vector< real_t >                              floaterRadii_;        // indexed by floaterUIDs_
   // density ratio per sphere floater (indexed by floaterUIDs_)
   std::vector< real_t >                              floaterDensities_;
   // Box (type 3) + cylinder (type 4) floaters share the plane-fit + clip V_sub pipeline.
   // A cylinder is stored as halfExtents = (r, r, L/2) with its axis on the body z-axis.
   struct ClipFloater {
      walberla::id_t    uid;
      bool              isCylinder;  // false = box
      Vector3< real_t > halfExtents; // box: half-extents; cylinder: (r, r, halfLength)
   };
   std::vector< ClipFloater >                         clipFloaters_;          // same order/size on all ranks
   // keyed by UID; populated by both compute* methods
   std::unordered_map< walberla::id_t, real_t >       submergedVolumes_;
   // keyed by UID; tapered F_pred_z per sphere (Option A vertical MEM correction); currently never filled
   // box/cylinder C_b − C_com (world frame), for the righting torque
   std::unordered_map< walberla::id_t, Vector3< real_t > > buoyancyCenterOffset_;
   // per-floater COM offset below the centroid (body z, const), for the pendulum righting torque
   std::unordered_map< walberla::id_t, real_t >       boxComOffsetZ_;
   // Option A: subtract F_pred_z from the vertical MEM force
   // (default) true → skip computeSubmergedVolumes()/computeSubmergedVolumesBox() in doDEMSubCycles: no vSub, no
   // righting torque, as in Bogner & Rude (2013) (buoyant support and righting moment come solely from the raw
   // MEM force/torque). false → run them (Antidunes' extended split).
   bool                                               useMemOnlyBuoyancy_{ true };
   // Only consulted when useMemOnlyBuoyancy_==false. true -> computeSubmergedVolumes()/
   // computeSubmergedVolumesBox() set vSub=vFull (buoyancyCenterOffset_=0) for every floater, skipping the
   // plane-fit/clip machinery (see ModelParams::useFullVolumeAnalyticBuoyancy, SimulationParams.h).
   bool                                               useFullVolumeAnalyticBuoyancy_{ false };
   // floaters whose TRANSLATION is held (rotation free): UID → pinned world position
   std::unordered_map< walberla::id_t, Vector3< real_t > > pinnedBoxPositions_;
   // UID → per-axis pin mask (nonzero = pin that axis, 0 = free)
   std::unordered_map< walberla::id_t, Vector3< real_t > > pinnedBoxAxes_;
   // UID → per-axis rotation pin mask (nonzero = zero that axis's angular velocity every sub-cycle, 0 = free)
   std::unordered_map< walberla::id_t, Vector3< real_t > > pinnedRotationAxes_;

   // [MIGRATION-NAN] support: hydrodynamicForce/oldHydrodynamicForce are not part of
   // ParticleGhostCopyNotification (ghost copies hold stale values) but are carried to a new owner by
   // ParticleMigrationNotification. UID -> last observed owner rank, to tell a fresh migration from an
   // already-owned particle whose force went bad otherwise.
   std::unordered_map< walberla::id_t, int > lastKnownOwner_;

   // [GHOST-VANISHED] support (see doDEMSubCycles): snapshot, rebuilt each main step at subCycle==0, of every
   // floater UID this rank can see (local or ghost), i.e. the population the last fluid-obstacle mapping used.
   std::unordered_set< walberla::id_t > ghostMappingSnapshot_;

   // [VELOCITY-JUMP] support (see doDEMSubCycles): UID -> velocity at the start of the previous main step.
   std::unordered_map< walberla::id_t, Vector3< real_t > > lastVelocityForJumpCheck_;

   // kernels — declaration order must match initializer-list order to avoid -Wreorder
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce;
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce;
   ComputeLubrication computeLubrication;
   ComputeContacts computeContacts;
   // ModelParameters.logFloaterHistory's per-rank CSV writer (no-op unless enabled). Owned here because
   // ParticleKernels outlives computeContacts/computeLubrication, which hold a non-owning pointer.
   FloaterHistoryLogger floaterHistoryLogger_;
   // Shared with SetupStages.h/BuildTimeloop.h (SimulationObjects::dormantFloaterUIDs): excludes dormant
   // (startInactive, not yet at settlingSteps) floaters from DEM contact pairs, as from fluid mapping.
   std::shared_ptr< std::unordered_set< walberla::id_t > > dormantFloaterUIDs_;
   // Shared with SetupStages.h (SimulationObjects::freshlyActivatedFloaterUIDs): wired into computeContacts
   // (setFreshlyActivatedUIDs) and cleared at the end of every doDEMSubCycles() call.
   std::shared_ptr< std::unordered_set< walberla::id_t > > freshlyActivatedFloaterUIDs_;
   // Built in the constructor from ModelParams::logFloaterHistoryReleaseStepFilter (empty when 0);
   // FloaterHistoryLogger keeps a raw pointer to it (same owner, so lifetimes match).
   std::unordered_set< walberla::id_t > loggedFloaterUIDs_;
   AddExternalForces addExternalForces;
   bool useLubricationCorrection;
   uint_t                                             evalFrequency_;
   uint_t                                             stepCounter_;
   bool                                               simParams_shadowEnabled_ = false;
   // sphere surface plane-fit lateral reach, in units of r (prm floaterPlaneFitRadiusFactor)
   real_t                                             planeFitRadiusFactor_;
   // box/cylinder plane-fit shell reach, in units of max(halfExtents) (prm boxPlaneFitBandFactor)
   real_t                                             boxPlaneFitBandFactor_;
   // static free-surface reference height (absoluteLiquidHeight); used in logs and the disabled Option A block
   real_t                                             surfaceHeight_;
   // body acceleration incl. gravity; for net-force logging
   std::shared_ptr< Vector3< real_t > >               acceleration_;

   // Live MEM hydrostatic reference height (ModelParameters.useDynamicReferenceHeight), shared with
   // HydrostaticDensityFtor via SimulationObjects and EMA-updated here from the floaters' fitted waterlines.
   // Force-side only (SimpleBB.h rhoHydStat never touches the fluid), so no PDF/reconstruction feedback.
   std::shared_ptr< real_t >                          memReferenceHeight_;
   bool                                                useDynamicReferenceHeight_ = false;

   // ModelParameters.useHydrodynamicForceTimeAveraging (see its own doc comment, SimulationParams.h):
   // gates the AverageHydrodynamicForceTorqueKernel call in reduceAndAvgHydroForces() below.
   bool                                                useHydrodynamicForceTimeAveraging_ = true;
};

} // namespace antidunes
} // namespace walberla