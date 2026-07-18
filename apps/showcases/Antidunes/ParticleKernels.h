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
#include "core/mpi/Reduce.h"
#include "field/FlagField.h"
#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/shape/Box.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"

#include "SimulationParams.h"
#include "SimulationObjects.h"

#include <algorithm>
#include <unordered_map>


namespace walberla {
namespace antidunes {

// Subclass that hides the hollow-cylinder overload of AnalyticContactDetection with exterior
// (solid cylinder) contact math. DoubleCast deduces the concrete functor type, so C++ overload
// resolution picks this method instead of the base-class version.
//
// Base class is GeneralContactDetection (not AnalyticContactDetection): pairs without an exact
// analytic overload — Box×Box, Cylinder×Box, Cylinder×Cylinder, Cylinder×Ellipsoid — fall through
// to its GJK/EPA path, which works on any convex shape via support() and produces the same
// contactNormal/penetrationDepth convention as the analytic functions (asserted float-equal in
// tests/mesa_pd/collision_detection/EPA.cpp). Analytic overloads (sphere-sphere, sphere-box,
// sphere-halfspace, …) still win overload resolution, so the bed-sphere hot path is unchanged.
// Cylinder×HalfSpace has its OWN overload below (same support()-based formula as Box×HalfSpace) —
// it used to fall through to the base's generic HalfSpace-vs-anything template instead, but that
// template's contactNormal comes out as -wallNormal (into the wall), the OPPOSITE sign from every
// other HalfSpace pair in this class (Sphere and Box both give +wallNormal, out of the wall — see
// detectSphereHalfSpaceCollision and the Box overload below). Discovered while writing
// tests/showcases_antidunes/ContactDetectionTest.cpp testCylinderHalfSpaceContact; fixed here to
// match the proven-correct Sphere/Box convention rather than left as a silent inconsistency.
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

   // Box vs HalfSpace — rotation-aware: normal is transformed to body frame via R^T before
   // calling support() in body frame and rotating the result back to world frame.
   // Same penetrationDepth / contactNormal / contactPoint convention as Sphere-HalfSpace:
   //   penetrationDepth < 0  → penetrating (contact active)
   //   contactNormal = hs.getNormal() (points into the domain, away from the wall)
   //   contactPoint on the wall surface (deepest box corner projected onto the plane)
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Box& bx,
                   const mesa_pd::data::HalfSpace& hs,
                   Accessor& ac)
   {
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

   // Cylinder vs HalfSpace — identical formula/convention to Box vs HalfSpace above (Cylinder has
   // its own exact support() too). Added as an explicit override instead of letting this pair fall
   // through to GeneralContactDetection's generic HalfSpace-vs-anything template: that generic path
   // returns contactNormal = -wallNormal (into the wall), the opposite sign from every other
   // HalfSpace pair here. See the class-level comment above and
   // tests/showcases_antidunes/ContactDetectionTest.cpp testCylinderHalfSpaceContact for the
   // derivation that surfaced the mismatch.
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
                   const mesa_pd::data::Cylinder& cyl,
                   const mesa_pd::data::HalfSpace& hs,
                   Accessor& ac)
   {
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

   // Box vs CylindricalBoundary (solid exterior) — ROTATION-AWARE. Works in the box BODY frame so it
   // is correct for any orientation (the previous version used a world-axis-aligned bounding box and
   // silently let a rotated box pierce the cylinder). The cylinder centre/axis are transformed by Rᵀ,
   // the closest point on the oriented box is found by clamping to the half-extents, the radial
   // distance to the (infinite) cylinder axis gives the penetration, and the normal/contact point are
   // transformed back to world.
   //   penetrationDepth = radial_dist - cylinder_radius  (< 0 → penetrating)
   //   contactNormal points radially outward from the cylinder axis toward the box (world frame)
   //   contactPoint = closest point on the box surface, world frame
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
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

      // Cylinder axis as a LINE in the box body frame (cylinder is infinite along its axis).
      mesa_pd::Vec3 d  = Rt * cb.getAxis();
      const real_t  dl = d.length();
      d = (dl > real_t(1e-10)) ? d / dl : mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
      const mesa_pd::Vec3 a0 = Rt * (cylCenter - c); // a point on the axis (the cylinder centre), body frame

      // Closest point on the oriented box to the INFINITE axis line: project the current estimate onto
      // the line, clamp the nearest line point into the box, refine once. (Using the far cylinder centre
      // directly would bias the clamp toward the wrong box corner for a rotated/elongated box.)
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
         // SAFETY: a body large enough to ENVELOP the thin pillar (axis inside the box) yields a
         // near-radius penetration with a degenerate normal, which the spring-dashpot turns into a
         // catastrophic force that blows up the coupled fluid (global NaN). Bound the reported overlap
         // so the contact stays a finite restoring push, and warn so the pathological geometry is visible.
         const real_t maxOverlap = real_t(1.5);
         if (penetrationDepth < -maxOverlap)
         {
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

   // Sphere vs finite Cylinder (floating log) — ANALYTIC, because the log sits on/next to thousands
   // of bed spheres and running GJK/EPA per bed-sphere pair per sub-cycle is both slower and noisier
   // than the closed form. Works in the cylinder body frame (axis = body z), so it is correct for any
   // orientation. Convention identical to sphere-sphere:
   //   penetrationDepth < 0 → penetrating; contactNormal points from the CYLINDER (idx2) to the
   //   SPHERE (idx1), world frame; contactPoint = mid-overlap point on the cylinder surface.
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

   // Finite Cylinder (floating log) vs CylindricalBoundary (solid pillar, infinite along its axis) —
   // ROTATION-AWARE, mirrors the Box×CylindricalBoundary overload above: the pillar axis is
   // transformed into the log's body frame, the closest point of the log to that line is found by a
   // project-and-clamp fixed-point iteration (radial + axial clamp instead of the box corner clamp),
   // and the radial distance to the axis gives the penetration.
   //   penetrationDepth = radial_dist - pillar_radius (< 0 → penetrating)
   //   contactNormal points radially outward from the pillar axis toward the log (world frame)
   template< typename Accessor >
   bool operator()(const size_t idx1, const size_t idx2,
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
         // Same envelop safety clamp as Box×CylindricalBoundary: a log enveloping the thin pillar
         // yields a near-radius overlap + degenerate normal → unbounded spring force → fluid NaN.
         const real_t maxOverlap = real_t(1.5);
         if (penetrationDepth < -maxOverlap)
         {
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
};


class ComputeLubrication {
public:
   ComputeLubrication(const std::shared_ptr< mesa_pd::domain::BlockForestDomain >& domain, real_t viscosity)
      : domain_(domain)
      , lubricationCorrectionKernel_(viscosity,
         [](real_t r) { return (real_c(0.001) + real_c(0.00007) * r) * r; })
   {}

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
            double_cast_(acd_.getIdx1(), acd_.getIdx2(), ac, lubricationCorrectionKernel_, ac,
                        acd_.getContactNormal(), acd_.getPenetrationDepth());
         }
      }
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain > domain_;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel_;
   mesa_pd::collision_detection::AnalyticContactDetection acd_;
   mesa_pd::kernel::DoubleCast double_cast_;
   mesa_pd::mpi::ContactFilter contact_filter_;
};


class ComputeContacts {
public:
   ComputeContacts(const std::shared_ptr< mesa_pd::domain::BlockForestDomain >& domain,
      real_t restitutionCoefficient, real_t kappa,
      real_t collisionTime, real_t particleTimestepSize, real_t frictionCoefficient,
      uint_t numParticleTypes = uint_t(2))
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
   }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx1, const size_t idx2, ParticleAccessor_T& ac){
      // Box-Box (and every other convex pair without an analytic overload) is handled by the
      // GJK/EPA fallback in SolidCylinderContactDetection's GeneralContactDetection base, so
      // floating boxes/cylinders collide with each other and with FixedBox pillars.
      if (double_cast_(idx1, idx2, ac, acd_, ac))
      {
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

            auto meff = real_c(1) / (ac.getInvMass(idx1) + ac.getInvMass(idx2));
            // Set stiffness for the actual type pair involved in this contact
            const uint_t t1 = std::min(ac.getType(acd_.getIdx1()), ac.getType(acd_.getIdx2()));
            const uint_t t2 = std::max(ac.getType(acd_.getIdx1()), ac.getType(acd_.getIdx2()));
            collisionResponse_.setStiffnessAndDamping(t1, t2, restitutionCoefficient_, collisionTime_, kappa_, meff);
            collisionResponse_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(), acd_.getContactNormal(),
                              acd_.getPenetrationDepth(), particleTimestepSize_);
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
};


class AddExternalForces {
public:
   AddExternalForces(const std::shared_ptr< mesa_pd::domain::BlockForestDomain >& domain,
      const std::shared_ptr< Vector3< real_t > >& acceleration,
      real_t particleDensityRatio)
      : domain_(domain)
      , acceleration_(acceleration)
      , particleDensityRatio_(particleDensityRatio)
      , submergedVolumes_(nullptr)
      , buoyancyCenterOffset_(nullptr)
   {}

   void setSubmergedVolumes(const std::unordered_map< walberla::id_t, real_t >* sv) {
      submergedVolumes_ = sv;
   }

   // Option A: per-sphere tapered phantom force F_pred_z to subtract from the vertical MEM force.
   // Enabled only when usePhantomForceCorrection is true; the map pointer stays null otherwise.
   void setPhantomForceZ(const std::unordered_map< walberla::id_t, real_t >* pf) {
      phantomForceZ_ = pf;
   }

   void setBuoyancyCenterOffset(const std::unordered_map< walberla::id_t, Vector3< real_t > >* off) {
      buoyancyCenterOffset_ = off;
   }

   void setComOffsetZ(const std::unordered_map< walberla::id_t, real_t >* com) {
      comOffsetZ_ = com;
   }

   // Wire up the step counter + eval frequency so the per-box force breakdown can be logged once
   // per eval step (the counter is owned by the enclosing ParticleKernels and advances once per step).
   void setStepInfo(const uint_t* stepCounter) { stepCounter_ = stepCounter; }

   // Set once per sub-cycle by doDEMSubCycles; true only on the first sub-cycle of an eval step, so
   // every floater (any count) logs its force breakdown exactly once per eval step.
   void setLogForces(bool logForces) { logForces_ = logForces; }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac){
      // Buoyancy acts only in Z because the hydrostatic MEM contribution is only in Z.
      // Type 1 = floating sphere, type 3 = floating box, type 4 = floating cylinder;
      // derive density from invMass×volume.
      const uint_t ptype = ac.getType(idx);
      const bool isFloater = (ptype == uint_t(1) || ptype == uint_t(3) || ptype == uint_t(4));
      const real_t rho = isFloater
         ? real_t(1) / (ac.getInvMass(idx) * ac.getVolume(idx))
         : particleDensityRatio_;

      const real_t vFull = ac.getVolume(idx);
      // For floaters use the precomputed submerged volume so the gas-phase portion does not
      // contribute a spurious upward buoyancy force.  Non-floaters are always fully submerged.
      real_t vSub = vFull;
      if (isFloater && submergedVolumes_) {
         auto it = submergedVolumes_->find(ac.getUid(idx));
         if (it != submergedVolumes_->end()) vSub = it->second;
      }

      // X/Y: full-volume driving force (pressure-gradient body force acts on the whole particle).
      // Z:   gravity (rho * vFull) minus correct buoyancy (vSub × g, fluid density = 1).
      Vector3< real_t > bodyForce(
         vFull * (*acceleration_)[0],
         vFull * (*acceleration_)[1],
         (rho * vFull - vSub) * (*acceleration_)[2]);

      auto hydro_force  = ac.getHydrodynamicForce(idx);
      auto hydro_torque = ac.getHydrodynamicTorque(idx);
      const Vector3< real_t > hydroForceRaw = hydro_force; // pre-scaling value, kept for the verification log

      // Option A — cancel the analytic phantom force in the VERTICAL MEM component only. F_pred_z is the
      // double-counted pile-up head (buoyancy the flat-reference subtraction failed to remove and the
      // geometric channel already supplies), which the MEM leaks as a spurious UPWARD lift.
      //
      // ONE-SIDED CLAMP: the phantom is strictly an UPWARD artifact, so cancelling it may only REDUCE a
      // positive (upward) MEM force toward zero — it must never turn it negative and manufacture a downward
      // force that was not there. Subtracting F_pred unclamped over-corrected in the first end-to-end test
      // (phantom.log): near 50% submersion F_pred marginally exceeded the raw MEM force, flipping it
      // negative, which sank the floater past 50% into the submerged-regime down-residual (a trap the
      // taper leaves uncorrected) and down to the bottom. Clamping to max(raw − F_pred, 0) removes the
      // levitation while structurally preventing that sink. Horizontal drag + dynamic fluctuations intact;
      // self-disables where there is no pile-up (F_pred → 0). Particle-side only → no PDF feedback.
      real_t phantomZ  = real_t(0); // clamped value actually subtracted
      real_t rawFPred  = real_t(0); // unclamped F_pred, for the log — reveals when the clamp saturates
      bool   fPredSeen = false;
      if (isFloater && phantomForceZ_) {
         auto itp = phantomForceZ_->find(ac.getUid(idx));
         if (itp != phantomForceZ_->end()) {
            rawFPred  = itp->second;
            fPredSeen = true;
            if (hydro_force[2] > real_t(0)) {
               // F_pred can be NEGATIVE when the floater sits in a local trough (localSurfaceZ < flat
               // reference, e.g. downstream of a pillar or in an antidune trough) — deltaH<0 there. That
               // regime is outside the model (derived for pile-up, deltaH>0) and must not be allowed to
               // ADD to the upward MEM force via hydro_force[2] -= (negative). Clamp to [0, hydro_force[2]].
               phantomZ = std::clamp(rawFPred, real_t(0), hydro_force[2]);
               hydro_force[2] -= phantomZ;
            }
         }
      }

      // Submersion-scaled hydrodynamic force (generalises the old binary vSub==0 airborne guard).
      // The MEM (SimpleBB) builds the force as if every near-boundary cell were full liquid (ρ≈1),
      // but a partially-emerged floater only pushes on a real fluid column over its SUBMERGED part.
      // The emerged part's "force" is a phantom: cells the floater just vacated are reconstructed/
      // refilled to liquid density, and a locally raised/dammed surface (pillar hump, displaced
      // water) is under-counted by the flat hydrostatic reference ρ_hyd(z)=1+3·a_z·(z−h_liquid).
      // Both leak a spurious upward impulse that can exceed the floater's weight and launch it
      // ballistically up into the gas (nine.log: airborne sphere F_hydro_z≈+0.021 vs weight 0.011).
      // Scaling F_hydro and the torque by the geometric submerged fraction V_sub/V_full removes the
      // leak smoothly: fully submerged → ×1 (unchanged), fully emerged → ×0 (the old guard), and
      // partially emerged → proportional. V_sub itself is read off the LIVE local waterline
      // (computeFloaterBuoyancy), so this is correct on raised/dammed/wake surfaces too.
      real_t submFrac = real_t(1);
      if (isFloater) {
         submFrac = (vFull > real_t(0)) ? std::clamp(vSub / vFull, real_t(0), real_t(1)) : real_t(0);
         // Geometric submersion scaling DELIBERATELY OFF: the phantom MEM lift is now corrected at its
         // source by the per-link fill-weighting φ(x_b) in SimpleBB (Bogner & Rüde 2013, Eq. 5.9), so
         // this coarse global V_sub/V_full factor is no longer applied. φ-weighting is the physical fix;
         // this scaling was a band-aid that masked the same leak with a single geometric factor and is
         // intentionally retired. submFrac is still computed for the [FloaterForces] log below.
         // (Two-scenario validation: changes.md "Phantom-force investigation".)
         // hydro_force  *= submFrac;
         // hydro_torque *= submFrac;
      }

      if (!std::isfinite(hydro_force[0]) || !std::isfinite(hydro_force[1]) || !std::isfinite(hydro_force[2]))
         WALBERLA_LOG_WARNING("[AddExternalForces] NaN/Inf hydrodynamic force on uid=" << ac.getUid(idx)
            << " pos=" << ac.getPosition(idx));
      if (!std::isfinite(hydro_torque[0]) || !std::isfinite(hydro_torque[1]) || !std::isfinite(hydro_torque[2]))
         WALBERLA_LOG_WARNING("[AddExternalForces] NaN/Inf hydrodynamic torque on uid=" << ac.getUid(idx)
            << " pos=" << ac.getPosition(idx));

      // Per-floater force breakdown, once per eval step (doDEMSubCycles sets logForces_ on sub-cycle 0,
      // so every sphere AND the box log exactly once). ac.getForce(idx) here = contact + lubrication
      // accumulated THIS sub-cycle (ComputeContacts runs before AddExternalForces); hydro_force and
      // bodyForce are the MEM and buoyancy/gravity/driving parts about to be added. This shows why
      // equal-density floaters end at different depths: compare F_buoy+grav (model), F_hydro (flow),
      // and F_contact (sphere-sphere pile / pillar / wall) per floater.
      if (logForces_ && isFloater)
      {
         const auto fContact = ac.getForce(idx);
         WALBERLA_LOG_INFO((ptype == uint_t(3) ? "[BoxForces]"
                            : ptype == uint_t(4) ? "[CylForces]" : "[FloaterForces]")
            << " uid=" << ac.getUid(idx) << " step=" << (stepCounter_ ? *stepCounter_ : uint_t(0))
            << "  subm=" << (vFull > real_t(0) ? real_t(100) * vSub / vFull : real_t(0)) << "%"
            << "  F_buoy+grav_z=" << bodyForce[2] << " (>0 up)"
            << "  F_drive_xy=(" << bodyForce[0] << "," << bodyForce[1] << ")"
            << "  F_hydro=(" << hydro_force[0] << "," << hydro_force[1] << "," << hydro_force[2] << ")"
            << (fPredSeen ? "  [phantomCorr F_pred_raw=" + std::to_string(rawFPred)
                                + " F_pred_applied(clamped)=" + std::to_string(phantomZ)
                                + (std::abs(rawFPred - phantomZ) > real_t(1e-9) ? " SATURATED" : "")
                                + " subtracted from F_hydro_z (was " + std::to_string(hydroForceRaw[2]) + ")]"
                          : "")
            << "  [submScale=" << submFrac << " applied to F_hydro_raw=("
            << hydroForceRaw[0] << "," << hydroForceRaw[1] << "," << hydroForceRaw[2] << ")]"
            << "  F_contact+lub=(" << fContact[0] << "," << fContact[1] << "," << fContact[2] << ")"
            << "  vel=" << ac.getLinearVelocity(idx)
            << "  pos=" << ac.getPosition(idx));
      }

      mesa_pd::addForceAtomic(idx, ac, hydro_force + bodyForce);

      // Buoyancy righting / capsizing moment (Bogner & Rüde 2013, Fig. 5.4): the buoyancy force
      // (vertical, +vSub·g) acts at the centre of buoyancy B = centroid of the submerged volume,
      // not at the centre of gravity G (the COM). When the box tilts, B shifts horizontally off G,
      // so the vertical buoyancy produces a moment (B − G) × F_buoy. We apply that moment here as an
      // explicit torque (the net force is unchanged, so the equilibrium depth is unaffected). The
      // sign is physical — it self-rights when B's metacentre is above G and tips further when below,
      // exactly as in the figure. Spheres store no offset (B is on the vertical through G → 0 torque).
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

      // Lowered centre of mass (pendulum righting torque). Gravity acts at the true COM, which a
      // comOffsetZ < 0 places BELOW the geometric centroid (= the box's rotation centre). When the
      // box tilts, this offset swings out to the side, so the weight produces a restoring moment
      // (COM − centroid) × F_grav that raises the metacentric height and keeps the box upright —
      // the configurable fix for the GM<0 rolling of a uniform-density log/cube. comOffsetZ = 0
      // reproduces the uniform box exactly (no torque).
      if (isFloater && comOffsetZ_)
      {
         auto ito = comOffsetZ_->find(ac.getUid(idx)); // map holds only nonzero offsets (guarded at build)
         if (ito != comOffsetZ_->end())
         {
            const mesa_pd::Mat3& R = ac.getRotation(idx).getMatrix();                     // body → world
            const Vector3< real_t > leverG = R * Vector3< real_t >(real_t(0), real_t(0), ito->second); // COM − centroid (world)
            const Vector3< real_t > Fgrav(real_t(0), real_t(0), rho * vFull * (*acceleration_)[2]);     // weight (down: a_z<0)
            totalTorque += Vector3< real_t >(leverG[1] * Fgrav[2] - leverG[2] * Fgrav[1],
                                             leverG[2] * Fgrav[0] - leverG[0] * Fgrav[2],
                                             leverG[0] * Fgrav[1] - leverG[1] * Fgrav[0]);
         }
      }
      // Torque breakdown, logged once per eval step: hydro_torque is the raw MEM/PDF momentum-exchange
      // torque (from the real LBM lattice — independent of the buoyancy volume/centroid method, sampling
      // or analytic) vs. the buoyancy-lever (self-righting) contribution added above from cbOff. Needed
      // to tell apart a genuine hydrodynamic roll torque from the buoyancy-lever torque, since the latter
      // is exact/self-righting for cylinders now (clipCylinderAnalytic) but the tiltY drift persisted.
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
   const std::unordered_map< walberla::id_t, real_t >*     submergedVolumes_;
   const std::unordered_map< walberla::id_t, real_t >*     phantomForceZ_ = nullptr; // Option A vertical MEM correction; null = disabled
   const std::unordered_map< walberla::id_t, Vector3< real_t > >* buoyancyCenterOffset_;
   const std::unordered_map< walberla::id_t, real_t >*     comOffsetZ_ = nullptr; // per-box COM offset below centroid

   // For the per-floater force-breakdown log.
   const uint_t* stepCounter_ = nullptr; // step number, for the log only
   bool          logForces_   = false;   // set true by doDEMSubCycles on the first sub-cycle of an eval step
};


// ─────────────────────────────────────────────────────────────────────────────────────────────────
// Analytical (exact) submerged volume + centre of buoyancy of an oriented box (OBB) clipped by the
// fitted free-surface plane — the closed-form counterpart to the 12³ point-sampling clip. Mirrors the
// Python proof of concept in `box_overlap_volume.py` (OBB ∩ half-space), but specialised to a box vs a
// single plane so no convex-hull library is needed: we split the box into 6 tetrahedra (Kuhn
// triangulation) and clip each tet against the half-space, summing volume and volume-weighted centroid.
//
//   h    : box half-extents (body frame)          nHat : fitted surface normal (world), LIQUID→GAS (up)
//   R    : rotation body→world (orthonormal)       xBar : a point on the fitted plane (world)
//   pos  : box centre (world)
// Returns  vSub  = submerged volume, and  cbOff = world-frame offset of the centre of buoyancy from the
// box centre (the lever arm B−G used for the righting torque) — same meaning/frame as the sampling clip.

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
// Method: the circular cross-section is invariant under rotation about the body Z (axis) direction, so
// azimuthally align the plane normal's horizontal component onto the body X axis first — in that frame
// the plane is  rho*x + nz*z = d  (no y dependence), so the submerged column length at a given x is a
// CLAMPED AFFINE function of x alone. Splitting the disk into the (at most) three regions where that
// clamp is saturated-low / linear / saturated-high, each region's contribution to volume and the first
// moments (∫x dV, ∫z dV; ∫y dV = 0 by the x-only dependence) is a closed-form combination of the
// standard circular antiderivatives ∫w, ∫x·w, ∫x²·w with w(x) = 2·sqrt(r²−x²) — the same family used for
// the circular-segment area. No numerical integration or point sampling is involved. This replaces the
// 12³-point Cartesian sampling grid for cylinders, which is not azimuthally symmetric at finite
// resolution and was the source of a spurious roll torque on axisymmetric cylinders (see changes.md).
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
      , vvIntegratorPreForce(simParams->computed.particleTimestepSize)
      , vvIntegratorPostForce(simParams->computed.particleTimestepSize)
      , computeLubrication(simObjs->mesapdDomain, simParams->computed.viscosity)
      , computeContacts(
         simObjs->mesapdDomain,
         simParams->particle.restitutionCoefficient,
         simParams->computed.particleKappa,
         simParams->computed.particleCollisionTime,
         simParams->computed.particleTimestepSize,
         simParams->particle.frictionCoefficient,
         uint_t(5)) // 5 types: 0=bed/wall, 1=sphere floater, 2=pillar, 3=box floater, 4=cylinder floater
      , addExternalForces(
         simObjs->mesapdDomain,
         simObjs->acceleration,
         simParams->particle.densityRatio)
      , useLubricationCorrection(simParams->particle.useLubricationCorrection)
      , evalFrequency_(simParams->evaluation.evaluationFrequency)
      , stepCounter_(uint_t(0))
      , ringBand_(simParams->particle.floaterRingBand)
      , planeFitRadiusFactor_(simParams->model.floaterPlaneFitRadiusFactor)
      , boxPlaneFitBandFactor_(simParams->model.boxPlaneFitBandFactor)
      , surfaceHeight_(simParams->staticReferenceHeight())
      , acceleration_(simObjs->acceleration)
      , memReferenceHeight_(simObjs->memReferenceHeight)
      , useDynamicReferenceHeight_(simParams->model.useDynamicReferenceHeight)
      , floaterPlaneOrigin_(simObjs->floaterPlaneOrigin)
      , floaterPlaneNormal_(simObjs->floaterPlaneNormal)
   {
      // UIDs are resolved by setupParticles via an indexed creation loop + allReduce(SUM),
      // so every rank already has the correct UID for every floater in simObjs->floaterUIDs /
      // floaterBoxUIDs.  Copying avoids the earlier bug where a per-rank assignIdx caused
      // UID collisions at index 0 after allReduce(SUM).
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
      }

      useAnalyticBoxBuoyancy_    = simParams->model.useAnalyticBoxBuoyancy;
      usePhantomForceCorrection_ = simParams->model.usePhantomForceCorrection;
      WALBERLA_LOG_INFO_ON_ROOT("[FloaterForces] phantom-force correction (Option A): "
         << (usePhantomForceCorrection_ ? "ON (F_hydro_z −= tapered F_pred_z)" : "OFF"));

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
      }
      WALBERLA_LOG_INFO_ON_ROOT("[BoxBuoyancy] V_sub method (box + cylinder): "
         << (useAnalyticBoxBuoyancy_ ? "ANALYTIC (box: OBB∩plane 6-tet clip; cylinder: disk∩plane closed form)"
                                     : "SAMPLING (12³ points)")
         << "  pinned-position floaters: " << pinnedBoxPositions_.size());

      // Self-identify the buoyancy sampling mode in the log (answers "did this run use the ring?").
      WALBERLA_LOG_INFO_ON_ROOT("[FloaterBuoyancy] sphere V_sub sampling: floaterRingBand=" << ringBand_
         << (ringBand_ > real_t(1e5) ? "  (OLD full-box plane-fit)" : "  (contact-ring restriction)"));

      WALBERLA_LOG_INFO_ON_ROOT("[DynamicH] MEM reference height: "
         << (useDynamicReferenceHeight_ ? "DYNAMIC (EMA-tracked from floaters' fitted local waterlines, alpha=0.2)"
                                        : "STATIC (flat H=" + std::to_string(surfaceHeight_) + ")"));
   }

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
      tp["MESAPD: ReduceProperty(HydrodynamicForceTorqueNotification)"].start();
      mesa_pd::mpi::ReduceProperty{}.operator()< mesa_pd::HydrodynamicForceTorqueNotification >(*particleStorage);
      tp["MESAPD: ReduceProperty(HydrodynamicForceTorqueNotification)"].end();

      if (initializeOldValues)
      {
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel{}, *particleAccessor);
      }

      tp["MESAPD: AverageHydrodynamicForceTorqueKernel"].start();
      particleStorage->forEachParticle(useOpenMP,
         mesa_pd::kernel::SelectLocal(), *particleAccessor,
         lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel{}, *particleAccessor);
      tp["MESAPD: AverageHydrodynamicForceTorqueKernel"].end();
   }

   template< typename TP >
   void doDEMSubCycles(timing::TimingPool<TP> & tp, uint_t numSubCycles = uint_t(1)){
      // One increment per LBM step, independent of whether spheres/boxes are present. The per-shape
      // submerged-volume functions early-return when their list is empty, so the counter must not be
      // gated on them (otherwise a box-only run never advances and every eval log fires each step).
      ++stepCounter_;
      // Submerged volumes computed once per LBM step (fill field is constant during DEM sub-cycles).
      // Both functions write into the same submergedVolumes_ map (keyed by UID).
      computeSubmergedVolumes();
      computeSubmergedVolumesBox();
      addExternalForces.setSubmergedVolumes(&submergedVolumes_);
      addExternalForces.setBuoyancyCenterOffset(&buoyancyCenterOffset_);
      addExternalForces.setComOffsetZ(&boxComOffsetZ_);
      // Option A: hand the phantom-force map to the force kernel only when the correction is enabled;
      // otherwise leave the pointer null so the subtraction is skipped entirely (pre-fix behaviour).
      addExternalForces.setPhantomForceZ(usePhantomForceCorrection_ ? &phantomForceZ_ : nullptr);
      addExternalForces.setStepInfo(&stepCounter_);

      const bool logForcesThisStep = (stepCounter_ % evalFrequency_ == uint_t(0));
      for (auto subCycle = uint_t(0); subCycle < numSubCycles; ++subCycle)
      {
         // Log the per-floater force breakdown once per eval step (only on the first sub-cycle).
         addExternalForces.setLogForces(logForcesThisStep && subCycle == uint_t(0));

         // velocity verlet pre-force (update positions) --------------------------------
         tp["MESAPD: VelocityVerletPreForceUpdate"].start();
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            vvIntegratorPreForce, *particleAccessor);
         tp["MESAPD: VelocityVerletPreForceUpdate"].end();

         // synchronize -----------------------------------------------------------------
         tp["MESAPD: syncCall 1"].start();
         syncCall();
         tp["MESAPD: syncCall 1"].end();

         // compute lubrication forces --------------------------------------------------
         if (useLubricationCorrection)
         {
            tp["MESAPD: ComputeLubrication"].start();
            particleStorage->forEachParticlePairHalf(useOpenMP,
               mesa_pd::kernel::ExcludeInfiniteInfinite(), *particleAccessor,
               computeLubrication, *particleAccessor);
            tp["MESAPD: ComputeLubrication"].end();
         }

         // compute contact forces ------------------------------------------------------
         tp["MESAPD: ComputeContacts"].start();
         particleStorage->forEachParticlePairHalf(useOpenMP,
            mesa_pd::kernel::ExcludeInfiniteInfinite(), *particleAccessor,
            computeContacts, *particleAccessor);
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

         // velocity verlet post-force (update velocities) ------------------------------
         tp["MESAPD: VelocityVerletPostForceUpdate"].start();
         particleStorage->forEachParticle(useOpenMP,
            mesa_pd::kernel::SelectLocal(), *particleAccessor,
            vvIntegratorPostForce, *particleAccessor);
         tp["MESAPD: VelocityVerletPostForceUpdate"].end();

         // pin translation of pinPosition boxes (rotation left free), before the sync so ghosts see it -
         pinBoxTranslations();

         // synchronize -----------------------------------------------------------------
         tp["MESAPD: syncCall 2"].start();
         syncCall();
         tp["MESAPD: syncCall 2"].end();
      }
   }

   // Holds selected TRANSLATIONAL DOFs of every pinPosition floater (sphere, box, or cylinder) while
   // leaving rotation free: for each world axis whose pin-mask component is nonzero, resets that
   // position component to the pinned value and zeroes that velocity component after each sub-cycle
   // integration; axes with a zero mask evolve normally. The floater still accumulates hydrodynamic +
   // buoyancy torque, so a box/cylinder still rotates to its equilibrium tilt in place (a sphere's own
   // rotation has no visible effect). Mask <1,1,1> = rigid hold; <1,1,0> pins x,y and lets z float
   // (finds its own buoyant depth without lateral drift). No-op when no floater requested pinning.
   void pinBoxTranslations()
   {
      if (pinnedBoxPositions_.empty()) return;
      for (size_t idx = 0; idx < particleAccessor->size(); ++idx)
      {
         const uint_t pt = particleAccessor->getType(idx);
         if (pt != uint_t(1) && pt != uint_t(3) && pt != uint_t(4)) continue; // pinned spheres, boxes, AND cylinders
         if (mesa_pd::data::particle_flags::isSet(
                particleAccessor->getFlags(idx), mesa_pd::data::particle_flags::GHOST)) continue;
         const auto it = pinnedBoxPositions_.find(particleAccessor->getUid(idx));
         if (it == pinnedBoxPositions_.end()) continue;
         const Vector3< real_t >& target = it->second;
         const Vector3< real_t >& mask   = pinnedBoxAxes_.at(particleAccessor->getUid(idx));
         Vector3< real_t > pos = particleAccessor->getPosition(idx);
         Vector3< real_t > vel = particleAccessor->getLinearVelocity(idx);
         for (uint_t d = 0; d < uint_t(3); ++d)
            if (std::abs(mask[d]) > real_t(0.5)) { pos[d] = target[d]; vel[d] = real_t(0); } // pin this axis; free axes untouched
         particleAccessor->setPosition(idx, pos);
         particleAccessor->setLinearVelocity(idx, vel);
      }
   }

   // Estimates the submerged volume of every floating sphere (type 1) once per LBM step.
   // For each Z-slice dz through the sphere, 4 cells just outside the sphere in XY are sampled
   // from the fill field. Their average fill level phi(dz) ≈ 1 below the free surface and ≈ 0
   // above it, so sum_dz phi(dz) approximates h (submersion depth). V_sub is then the spherical
   // cap formula pi*h^2*(3r-h)/3, clamped to [0, V_sphere].
   // Computes the submerged volume of every floating sphere (type 1) once per LBM step.
   //
   // Strategy: plane-fitting using the NormalSweep normal field.
   //
   // NormalSweep (run every LBM step) already stores the unit free-surface normal n̂
   // at every interface cell.  Convention: n̂ points FROM LIQUID TO GAS (i.e. upward
   // for a flat horizontal surface).
   //
   // For each sphere we sample INTERFACE cells in a (r+2)-radius box but keep only those
   // within `ringBand` of the sphere SURFACE (|dist−r| < ringBand) — the waterline contact
   // ring ON the sphere, excluding the surrounding far-field free surface, whose height under
   // waves is at a different phase and would bias the fitted plane. Restricting to the ring also
   // makes x̄ ≈ the circle centre = the perpendicular foot of C on the plane, so (x̄−C) is
   // nearly parallel to the true normal and the fit becomes robust to noise in n̂.
   // We weight by φ(1−φ) which peaks at the transition (φ≈0.5) and vanishes in bulk
   // liquid/gas, and accumulate:
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
   // This is geometrically exact for any free-surface orientation — works correctly
   // even when waves make the local surface non-horizontal.  Uses only pre-existing
   // waLBerla fields (fillFieldID_, flagFieldID_, normalFieldID_); no eigensolver needed.
   //
   // Fallback: if wSum ≈ 0 (no interface cells near the sphere) the sphere is either
   // fully in gas or fully in liquid.  The average fill φ̄ of the non-obstacle
   // cells in the sampling box distinguishes the two cases:
   //   φ̄ ≈ 0 (gas)   → h_cap = 0
   //   φ̄ ≈ 1 (liquid) → h_cap = 2r
   // This prevents the wrong "assume fully submerged" that launched a sphere starting
   // entirely in the gas phase upward with a huge erroneous buoyancy force.
   void computeSubmergedVolumes()
   {
      const size_t nF = floaterUIDs_.size();
      if (nF == 0) { submergedVolumes_.clear(); return; }

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
      // Fallback: when no interface cells exist (wSum≈0), the sphere is entirely in
      // gas or liquid.  In that case ALL valid (gas/liquid) neighbours have the same
      // fill (≈0 or ≈1), so averaging any of them gives the correct answer even when
      // multiple ranks contribute.  Pillars, NoSlip walls, and other obstacles are
      // excluded by the fluidMask check — we keep looking past them automatically.
      //   [9]   aboveSum  — Σ φ of non-obstacle cells in the slab just ABOVE the sphere top (cavity probe)
      //   [10]  aboveCount— number of such cells (cavity probe)
      // Cavity probe: liquid sitting on top of the sphere means the gas around its surface is an
      // enclosed cavity, not the real free surface → it is submerged. This is the reference-free
      // replacement for the old fixed absoluteLiquidHeight bound.
      constexpr size_t stride = 11;
      std::vector< real_t > accum(nF * stride, real_t(0));

      // De-contamination of touching/clustered floaters: a sampled cell belongs to the floater whose
      // SURFACE it is nearest. nearerOtherFloater(p, fi, myDistErr) returns true if some OTHER floater's
      // surface is closer to p than fi's is (by myDistErr = |dist_to_fi_surface|), so the cell is the
      // neighbour's — skip it for fi. Without this, a neighbour's waterline pollutes fi's plane-fit (the
      // ten-sphere "localSurfaceZ scattered 26..80" symptom). X is the only periodic axis (min-image in X).
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
            const real_t planeLatR = planeFitRadiusFactor_ * r; // lateral reach of the surface plane-fit: sample
                                                                // the surrounding surface OUT TO planeFitRadiusFactor_*r
                                                                // (prm: floaterPlaneFitRadiusFactor, default 1.8) so
                                                                // it reaches past any crater/gas-moat the sphere digs.
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
                  // Cavity probe: average fill in the slab JUST ABOVE the sphere top (offZ in (r, r+4],
                  // within the top disk). Liquid here = water sits on top of the sphere = the gas around
                  // its surface is an enclosed (numerical) cavity, not the real free surface → submerged.
                  // A normally floating sphere has its top in the atmosphere, so this stays gas → no false
                  // trigger. Reference-free replacement for the old fixed absoluteLiquidHeight bound.
                  const real_t lateral = std::sqrt(offX * offX + offY * offY);
                  if (offZ > r && offZ <= r + real_t(4) && lateral < r
                      && !nearerOtherFloater(samplePos, fi, offZ - r))
                  {
                     accum[fi * stride + 9]  += fillField->get(localCell[0], localCell[1], localCell[2]);
                     accum[fi * stride + 10] += real_t(1);
                  }
               }

               // Plane-fitting channel: fit the free surface over a WIDE lateral disk (radius planeLatR)
               // around the sphere axis — reaching PAST any crater/gas-moat the sphere digs — so the fitted
               // waterline tracks the undisturbed surrounding level instead of the depressed moat. Keep only
               // interface cells with an UPWARD normal (n_z>0): crater/cavity walls point sideways/down and
               // would drag the plane into the moat. (Wider than the old contact ring → can bias under
               // strong waves, but correct for a calm floating body.)
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
         real_t normalCoherence = real_t(0);        // |Σw n̂|/Σw ∈ [0,1]; low ⇒ normals cancel (cavity / neighbour contamination)

         // --- REFERENCE-FREE submerged depth: read the waterline locally, no absoluteLiquidHeight ---
         // The free surface is wherever the fill field transitions, sampled around THIS sphere:
         //   • contact ring present (wTot>0): the surface crosses the sphere -> plane-fit the local
         //     waterline from the ring interface cells. h_signed = n̂·(x̄−C) is the signed distance of
         //     the centre below that plane; the clamp covers the just-submerged / just-emerged ends too.
         //   • no ring: the sphere is fully in one phase -> the average fill φ̄ of the surrounding
         //     non-obstacle cells decides (φ̄>0.5 liquid -> submerged 2r, else gas -> emerged 0).
         // No flat global surface is assumed, so this holds in deformed/chaotic surfaces (piles,
         // pillars, wakes) where a fixed reference is simply wrong. TRADE-OFF vs the old code: the
         // fixed-surface "cavity-proof" override is gone — a gas pocket hugging a deeply-submerged
         // sphere can again under-report V_sub (we revisit that locally if it shows up).
         if (wTot > real_t(1e-6))
         {
            const Vector3< real_t > nSum(accum[fi*stride+1], accum[fi*stride+2], accum[fi*stride+3]);
            const Vector3< real_t > xBar(accum[fi*stride+4] / wTot,
                                         accum[fi*stride+5] / wTot,
                                         accum[fi*stride+6] / wTot);
            localSurfaceZ   = xBar[2];
            normalCoherence = nSum.length() / wTot; // 1 ⇒ all ring normals aligned; →0 ⇒ they cancel

            // Feed this floater's clean fitted waterline into the dynamic-H average (ambient surface,
            // not this floater's own crater — the ring reaches out to 1.8r, past any dug-in moat).
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

         // CAVITY BLEND (reference-free, SMOOTH): blend the plane-fit h_cap toward fully-submerged (2r)
         // as the slab just ABOVE the sphere top fills with liquid. avgFillAbove is the mean fill of
         // those "ceiling" cells (gas→0, liquid→1, interface→partial; obstacle cells excluded). Water
         // overhead ⇒ the gas around the sphere is an enclosed cavity ⇒ submerged. A smoothstep over
         // [0.3, 0.7] keeps the buoyancy force CONTINUOUS — the old hard >0.5 snap kicked the heave into
         // a growing oscillation that broke the free-surface solver.
         const real_t aboveCnt     = accum[fi * stride + 10];           // # ceiling cells checked
         const real_t avgFillAbove = (aboveCnt > real_t(0)) ? (accum[fi * stride + 9] / aboveCnt) : real_t(0);
         const real_t tCav = std::clamp((avgFillAbove - real_t(0.3)) / real_t(0.4), real_t(0), real_t(1));
         const real_t wCav = tCav * tCav * (real_t(3) - real_t(2) * tCav); // smoothstep weight ∈ [0,1]
         h_cap = (real_t(1) - wCav) * h_cap + wCav * (real_t(2) * r);

         const real_t vSub = sphereCapVolume(r, h_cap);

         submergedVolumes_[floaterUIDs_[fi]] = vSub;

         // Diagnostic plane fit for VTK output (BuildTimeloop.h, "floater"'s planeOrigin/planeNormal) —
         // never read back into the physics. Built from the FINAL, post-cavity-blend h_cap (not the raw
         // plane-fit's own h_cap from above), by sliding the plane along fittedNormal to the signed offset
         // (h_cap−r) that reproduces exactly this vSub — so a ParaView clip against the exported plane
         // always matches the volume actually used, even when the cavity blend has pulled h_cap away from
         // what the raw fitted plane alone would give. No fittedNormal (fully submerged/emerged, or a
         // cancelled-normal cavity ring) -> nothing meaningful to draw.
         if (fittedNormal.sqrLength() > real_t(1e-12))
         {
            const Vector3< real_t > xBarCorrected = C + (h_cap - r) * fittedNormal;
            if (floaterPlaneOrigin_) (*floaterPlaneOrigin_)[floaterUIDs_[fi]] = xBarCorrected;
            if (floaterPlaneNormal_) (*floaterPlaneNormal_)[floaterUIDs_[fi]] = fittedNormal;
         }
         else
         {
            if (floaterPlaneOrigin_) floaterPlaneOrigin_->erase(floaterUIDs_[fi]);
            if (floaterPlaneNormal_) floaterPlaneNormal_->erase(floaterUIDs_[fi]);
         }

         // --- PHANTOM FORCE (Option A) — computed EVERY step and stored for AddExternalForces to subtract.
         // F_pred_z = Δp·A_proj·taper is the double-counted pile-up head on the open wetted cap (see the
         // block below for the full derivation). Computed here unconditionally (not just at eval steps) so
         // the vertical MEM correction is available every sub-cycle; only the LOGGING is gated.
         const real_t submFrac  = std::clamp(vSub / vFull, real_t(0), real_t(1));
         constexpr real_t kneeSub = real_t(0.5);  // taper begins once the body is more than half submerged
         const real_t tS        = std::clamp((submFrac - kneeSub) / (real_t(1) - kneeSub), real_t(0), real_t(1));
         const real_t taper     = real_t(1) - (real_t(3) * tS * tS - real_t(2) * tS * tS * tS); // 1→0 smoothstep
         const real_t az        = (*acceleration_)[2];
         const real_t deltaH    = localSurfaceZ - surfaceHeight_;   // uncancelled head (LU), >0 when piled up
         const real_t deltaP    = -az * deltaH;                     // excess pressure = |a_z|·Δh
         const real_t rContact  = std::sqrt(std::max(real_t(0), h_cap * (real_t(2) * r - h_cap)));
         const real_t aProj     = math::pi * rContact * rContact;   // wetted-cap projected disk area
         const bool   openSurf  = (h_cap > real_t(0) && h_cap < real_t(2) * r); // straddling the surface
         const real_t fPredRaw  = openSurf ? deltaP * aProj : real_t(0);        // untapered (for the log)
         const real_t fPredZ    = fPredRaw * taper;                             // tapered → the Option A correction
         phantomForceZ_[floaterUIDs_[fi]] = fPredZ;

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

            // --- PHANTOM DECOMPOSITION (hypothesis #1: uncancelled pile-up head on the open wetted cap) ---
            // The MEM force subtracts a FLAT hydrostatic reference ρ_hyd(z)=1+3·a_z·(z−oldFixedRef). Where the
            // local surface is raised to localSurfaceZ (pile-up), the extra water head Δh=localSurfaceZ−oldFixedRef
            // is NOT subtracted, leaving a uniform excess pressure Δp=|a_z|·Δh acting on the wetted lower cap.
            // Because the top is gas (open surface), that excess is unbalanced ⇒ net UP force
            //   F_pred_z = Δp · A_proj · taper,   A_proj = π·r_contact² = π·h_cap·(2r−h_cap)  (wetted-cap rim disk).
            // Fully submerged (closed body) ⇒ uniform Δp has zero gradient ⇒ F_pred_z→0 (the observed sign flip).
            // SUBMERSION TAPER (see the every-step block above where these terms are computed): A_proj alone is
            // symmetric about 50% submersion but the phantom is not — taper(s)=1−smoothstep((s−0.5)/0.5) drives
            // F_pred→0 as s→100% (closed body), leaving the validated ≤50% regime exact. Logged here for the
            // live comparison against measured F_hydro_z.
            WALBERLA_LOG_INFO_ON_ROOT("[PhantomDecomp] uid=" << floaterUIDs_[fi]
               << "  deltaH=" << deltaH << " LU (localSurfaceZ=" << localSurfaceZ << " − flatRef=" << surfaceHeight_ << ")"
               << "  |a_z|=" << -az
               << "  deltaP=" << deltaP
               << "  r_contact=" << rContact << "  A_proj=" << aProj
               << "  openSurface=" << (openSurf ? "yes" : "no(closed)")
               << "  submFrac=" << submFrac << "  taper=" << taper
               << "  F_pred_raw=" << fPredRaw
               << "  F_pred_z=" << fPredZ << " (tapered, >0 up)  [compare to measured F_hydro_z above/next]");
         }
      }

      // Dynamic MEM reference height: EMA-blend the (wTot-weighted) ambient waterline average into
      // memReferenceHeight_. Force-side only (SimpleBB's rhoHydStat feeds only forceMEM, never the fluid —
      // see AntidunesBoundaryHandling / SimpleBB.h), so this has no PDF/reconstruction feedback loop. Held
      // (not snapped) when no floater has a clean ring this step. See changes.md (2026-07 dynamic-H
      // history) for why this is a STATIC/offset correction only — it cancels the uncancelled pile-up head
      // exactly, but cannot remove a genuinely velocity-dependent (dynamic suction) residual.
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

   // Computes the submerged volume of every floating box (type 3) once per LBM step.
   //
   // Strategy: integrate the fill field φ over all LBM cells whose centres lie inside
   // the (possibly rotated) box.  Each such cell contributes φ × 1 LU³ to V_sub:
   //   - φ ≈ 1  (liquid cell, below free surface) → full contribution
   //   - φ ≈ 0  (gas cell, above free surface)    → no contribution
   //   - φ ∈ (0,1) (interface cell)               → proportional contribution
   //
   // Cells inside the box are MovingObstacle cells; their fill values are maintained by
   // ObstacleFillLevelSweep (propagated inward from neighbouring fluid cells each step)
   // and by the fill-field state from when the cell was last fluid.  Boundary-layer cells
   // are updated every step; deeper interior cells change slowly and remain bounded [0,1].
   //
   // Rotation is handled exactly: cell centres are transformed to the box body frame
   // before the containment test, so the method is valid for any tilt angle.  No rate-
   // limiting is needed because V_sub ∈ [0, V_full] by construction (φ ∈ [0,1]).
   void computeSubmergedVolumesBox()
   {
      buoyancyCenterOffset_.clear(); // box/cylinder-only; recomputed every step
      const size_t nB = clipFloaters_.size(); // boxes (type 3) + cylinders (type 4)
      if (nB == 0) return;

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
      // φ≈0.5 waterline. Per box bi in `acc` (stride 8): [0]=wSum [1..3]=Σw·n̂ [4..6]=Σw·x [7]=Σφ (shell).
      constexpr size_t stride = 8;
      std::vector< real_t > acc(nB * stride, real_t(0));
      std::vector< real_t > fluidCnt(nB, real_t(0)); // # non-obstacle shell cells (fallback denominator)

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
            // Wide lateral reach so the surface plane-fit samples PAST any crater/gas-moat the box digs,
            // anchoring the fitted waterline to the undisturbed surrounding surface (mirrors the sphere).
            // Configurable via prm ModelParameters.boxPlaneFitBandFactor (default 0.8, matching the old
            // hardcoded value) — raise it if the fit still looks contaminated by the box's own local
            // disturbance once the body has stopped moving.
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
               }

               // Plane-fitting channel: only interface cells carry a valid normal.
               if (!isFlagSet(cellFlag, interfaceFlag)) continue;

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
            const bool logStep = (stepCounter_ % evalFrequency_ == uint_t(0));
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

            // ---- Analytical clip (exact OBB∩plane for boxes, exact disk∩plane for cylinders; cheap) —
            //      used for the force when selected, and always logged as the cross-check against the
            //      sampling estimate. ----
            Vector3< real_t > cbAnalytic(real_t(0));
            if (isCyl) clipCylinderAnalytic(hx, hz, R, nHat, xBar, pos, vSubAnalytic, cbAnalytic);
            else       clipBoxAnalytic(Vector3< real_t >(hx, hy, hz), R, nHat, xBar, pos, vSubAnalytic, cbAnalytic);

            // ---- Sampling clip (12³ points; reference). Computed when it is the selected method, or at
            //      eval steps so the log can compare it against the analytical result. For cylinders,
            //      grid points outside the circular cross-section (body xy) are rejected, so the kept
            //      samples uniformly fill the cylinder volume. ----
            Vector3< real_t > cbSample(real_t(0));
            if (!useAnalyticBoxBuoyancy_ || logStep)
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

            if (useAnalyticBoxBuoyancy_) { vSub = vSubAnalytic; cbOff = cbAnalytic; }
            else                         { vSub = vSubSample;   cbOff = cbSample;   }

            // Diagnostic plane fit for VTK output (BuildTimeloop.h, planeOrigin/planeNormal) — never
            // read back into the physics. (xBar, nHat) is exactly the plane clipCylinderAnalytic()/
            // clipBoxAnalytic() cut against, above — so a ParaView clip against this exported plane
            // reproduces vSubAnalytic exactly, unlike the sphere case there is no separate post-hoc
            // blend here to correct for. CAVEAT: if useAnalyticBoxBuoyancy_ is OFF (12³ sampling
            // selected for the force instead), vSub above is vSubSample, not vSubAnalytic — the
            // visualized clip then matches the analytic cross-check value (typically within a few
            // percent per the logged xcheck diff), not bit-for-bit the sampled value actually used.
            if (floaterPlaneOrigin_) (*floaterPlaneOrigin_)[clipFloaters_[bi].uid] = xBar;
            if (floaterPlaneNormal_) (*floaterPlaneNormal_)[clipFloaters_[bi].uid] = nHat;
         }
         else
         {
            // No interface cells near the box: it is entirely in gas or entirely in liquid.
            const real_t cnt     = fluidCnt[bi];
            const real_t avgFill = (cnt > real_t(0)) ? acc[bi * stride + 7] / cnt : real_t(0);
            vSub  = (avgFill > real_t(0.5)) ? vFull : real_t(0);
            cbOff = Vector3< real_t >(real_t(0)); // symmetric → no righting lever
            // No clean plane fit -> nothing meaningful to draw; clear any stale entry.
            if (floaterPlaneOrigin_) floaterPlaneOrigin_->erase(clipFloaters_[bi].uid);
            if (floaterPlaneNormal_) floaterPlaneNormal_->erase(clipFloaters_[bi].uid);
         }

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
               << "  B-G offset=(" << cbOff[0] << "," << cbOff[1] << "," << cbOff[2] << ")"
               << "  |B-G|_horiz=" << std::sqrt(cbOff[0] * cbOff[0] + cbOff[1] * cbOff[1])
               << "  comOffsetZ=" << (boxComOffsetZ_.count(clipFloaters_[bi].uid)
                                         ? boxComOffsetZ_.at(clipFloaters_[bi].uid) : real_t(0))
               << "  method=" << (useAnalyticBoxBuoyancy_ ? "ANALYTIC" : "SAMPLING")
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
            // For a FIXED box held at heel α this is one point on the stability curve M(α); sweep α ∈ [0,30°]
            // at 3 resolutions and compare to the ideal analytic curve (this uses the analytic clip's exact
            // submerged centroid, so at high resolution it converges to the ideal). We log |M| and M_y, the
            // heel angle, the righting arm GZ = |M|/(V_sub·|g|), and V_sub, one line per eval step.
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
               << "  |g|=" << gMag << "  method=" << (useAnalyticBoxBuoyancy_ ? "ANALYTIC" : "SAMPLING"));
         }
      }
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
   std::vector< walberla::id_t >                      floaterUIDs_;         // sphere floaters, same size on all ranks
   std::vector< real_t >                              floaterRadii_;        // indexed by floaterUIDs_
   std::vector< real_t >                              floaterDensities_;    // density ratio per sphere floater (indexed by floaterUIDs_)
   // Box (type 3) + cylinder (type 4) floaters share the plane-fit + clip V_sub pipeline.
   // A cylinder is stored as halfExtents = (r, r, L/2) with its axis on the body z-axis.
   struct ClipFloater {
      walberla::id_t    uid;
      bool              isCylinder;  // false = box
      Vector3< real_t > halfExtents; // box: half-extents; cylinder: (r, r, halfLength)
   };
   std::vector< ClipFloater >                         clipFloaters_;          // same order/size on all ranks
   std::unordered_map< walberla::id_t, real_t >       submergedVolumes_;    // keyed by UID; populated by both compute* methods
   std::unordered_map< walberla::id_t, real_t >       phantomForceZ_;       // keyed by UID; tapered F_pred_z per sphere (Option A vertical MEM correction)
   std::unordered_map< walberla::id_t, Vector3< real_t > > buoyancyCenterOffset_; // box C_b − C_com (world frame); box-only, for righting torque
   std::unordered_map< walberla::id_t, real_t >       boxComOffsetZ_;       // per-box COM offset below centroid (body z, const); for the pendulum righting torque
   bool                                               useAnalyticBoxBuoyancy_{ false }; // analytical V_sub (box: OBB∩plane clip; cylinder: disk∩plane closed form) instead of the 12³ sampling
   bool                                               usePhantomForceCorrection_{ false }; // Option A: subtract F_pred_z from vertical MEM force
   std::unordered_map< walberla::id_t, Vector3< real_t > > pinnedBoxPositions_;  // boxes whose TRANSLATION is held (rotation free); UID → pinned world position
   std::unordered_map< walberla::id_t, Vector3< real_t > > pinnedBoxAxes_;       // UID → per-axis pin mask (nonzero = pin that axis, 0 = free)

   // kernels — declaration order must match initializer-list order to avoid -Wreorder
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce;
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce;
   ComputeLubrication computeLubrication;
   ComputeContacts computeContacts;
   AddExternalForces addExternalForces;
   bool useLubricationCorrection;
   uint_t                                             evalFrequency_;
   uint_t                                             stepCounter_;
   real_t                                             ringBand_; // sphere V_sub contact-ring half-width (prm: floaterRingBand)
   real_t                                             planeFitRadiusFactor_; // sphere surface plane-fit lateral reach, in units of r (prm: floaterPlaneFitRadiusFactor)
   real_t                                             boxPlaneFitBandFactor_; // box/cylinder plane-fit shell reach, in units of max(halfExtents) (prm: boxPlaneFitBandFactor)
   real_t                                             surfaceHeight_; // mean free-surface height (absoluteLiquidHeight); global V_sub bound
   std::shared_ptr< Vector3< real_t > >               acceleration_;  // body acceleration incl. gravity; for net-force logging

   // Live MEM hydrostatic reference height (ModelParameters.useDynamicReferenceHeight). Shared with
   // HydrostaticDensityFtor via SimulationObjects; EMA-updated here from the floaters' fitted local
   // waterlines each step. Force-side only (SimpleBB.h rhoHydStat never touches the fluid), so updating
   // it has no PDF/reconstruction feedback — see the [DynamicH] log / changes.md for the 2026-07 history.
   std::shared_ptr< real_t >                          memReferenceHeight_;
   bool                                                useDynamicReferenceHeight_ = false;

   // Per-floater buoyancy plane-fit, diagnostic/VTK only (see SimulationObjects.h). Shared with
   // BuildTimeloop.h's "floaterPlane" VTK output; written every step here, never read back into the
   // physics. Null-checked before writing (both are always set in setupBoundaryHandling, but a null
   // guard costs nothing and keeps this class usable standalone in tests).
   std::shared_ptr< std::unordered_map< walberla::id_t, Vector3< real_t > > > floaterPlaneOrigin_;
   std::shared_ptr< std::unordered_map< walberla::id_t, Vector3< real_t > > > floaterPlaneNormal_;
};

} // namespace antidunes
} // namespace walberla