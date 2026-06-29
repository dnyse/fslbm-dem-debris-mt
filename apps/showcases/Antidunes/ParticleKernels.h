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
class SolidCylinderContactDetection : public mesa_pd::collision_detection::AnalyticContactDetection
{
public:
   using AnalyticContactDetection::operator(); // prevent name hiding of all base overloads

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
      // Box-Box contact is not implemented in AnalyticContactDetection (no GJK fallback).
      // Skip defensively — Box-Box pairs only arise if a FixedBox pillar coincides with a
      // FloatingBox, which cannot happen in the current pillar configs (cylinder pillars only).
      if (ac.getBaseShape(idx1)->getShapeType() == mesa_pd::data::Box::SHAPE_TYPE &&
          ac.getBaseShape(idx2)->getShapeType() == mesa_pd::data::Box::SHAPE_TYPE) return;

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
      // Type 1 = floating sphere, type 3 = floating box; derive density from invMass×volume.
      const uint_t ptype = ac.getType(idx);
      const bool isFloater = (ptype == uint_t(1) || ptype == uint_t(3));
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
         hydro_force  *= submFrac;
         hydro_torque *= submFrac;
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
         WALBERLA_LOG_INFO((ptype == uint_t(3) ? "[BoxForces]" : "[FloaterForces]")
            << " uid=" << ac.getUid(idx) << " step=" << (stepCounter_ ? *stepCounter_ : uint_t(0))
            << "  subm=" << (vFull > real_t(0) ? real_t(100) * vSub / vFull : real_t(0)) << "%"
            << "  F_buoy+grav_z=" << bodyForce[2] << " (>0 up)"
            << "  F_drive_xy=(" << bodyForce[0] << "," << bodyForce[1] << ")"
            << "  F_hydro=(" << hydro_force[0] << "," << hydro_force[1] << "," << hydro_force[2] << ")"
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
      mesa_pd::addTorqueAtomic(idx, ac, totalTorque);
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain >   domain_;
   std::shared_ptr< Vector3< real_t > >                    acceleration_;
   real_t                                                   particleDensityRatio_;
   const std::unordered_map< walberla::id_t, real_t >*     submergedVolumes_;
   const std::unordered_map< walberla::id_t, Vector3< real_t > >* buoyancyCenterOffset_;
   const std::unordered_map< walberla::id_t, real_t >*     comOffsetZ_ = nullptr; // per-box COM offset below centroid

   // For the per-floater force-breakdown log.
   const uint_t* stepCounter_ = nullptr; // step number, for the log only
   bool          logForces_   = false;   // set true by doDEMSubCycles on the first sub-cycle of an eval step
};



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
         uint_t(4)) // 4 types: 0=bed/wall, 1=sphere floater, 2=pillar, 3=box floater
      , addExternalForces(
         simObjs->mesapdDomain,
         simObjs->acceleration,
         simParams->particle.densityRatio)
      , useLubricationCorrection(simParams->particle.useLubricationCorrection)
      , evalFrequency_(simParams->evaluation.evaluationFrequency)
      , stepCounter_(uint_t(0))
      , ringBand_(simParams->particle.floaterRingBand)
      , surfaceHeight_(simParams->computed.absoluteLiquidHeight)
      , acceleration_(simObjs->acceleration)
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
      }

      const size_t nB = simParams->floatingBoxes.size();
      floaterBoxUIDs_ = simObjs->floaterBoxUIDs;
      floaterBoxHalfExtents_.resize(nB);
      for (size_t bi = 0; bi < nB; ++bi)
      {
         floaterBoxHalfExtents_[bi] = simParams->floatingBoxes[bi].halfExtents;
         // Per-box COM offset below the centroid (constant); keyed by UID for the righting torque.
         if (std::abs(simParams->floatingBoxes[bi].comOffsetZ) > real_t(0))
            boxComOffsetZ_[floaterBoxUIDs_[bi]] = simParams->floatingBoxes[bi].comOffsetZ;
      }

      // Self-identify the buoyancy sampling mode in the log (answers "did this run use the ring?").
      WALBERLA_LOG_INFO_ON_ROOT("[FloaterBuoyancy] sphere V_sub sampling: floaterRingBand=" << ringBand_
         << (ringBand_ > real_t(1e5) ? "  (OLD full-box plane-fit)" : "  (contact-ring restriction)"));
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

         // synchronize -----------------------------------------------------------------
         tp["MESAPD: syncCall 2"].start();
         syncCall();
         tp["MESAPD: syncCall 2"].end();
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
            const real_t planeLatR = real_t(1.8) * r;   // lateral reach of the surface plane-fit: sample the
                                                        // surrounding surface OUT TO 1.8r so it reaches past
                                                        // any crater/gas-moat the sphere digs.
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

         const real_t vSub = math::pi * h_cap * h_cap * (real_t(3) * r - h_cap) / real_t(3);

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
      buoyancyCenterOffset_.clear(); // box-only; recomputed every step
      const size_t nB = floaterBoxUIDs_.size();
      if (nB == 0) return;

      // Gather centre positions and rotation matrices from owner ranks.
      // Non-owner ranks contribute 0; allReduce(SUM) distributes the owner's values.
      std::vector< real_t > cenX(nB, real_t(0));
      std::vector< real_t > cenY(nB, real_t(0));
      std::vector< real_t > cenZ(nB, real_t(0));
      std::vector< real_t > rotFlat(nB * 9, real_t(0)); // row-major R (body→world)

      for (size_t idx = 0; idx < particleAccessor->size(); ++idx)
      {
         if (particleAccessor->getType(idx) != uint_t(3)) continue;
         if (mesa_pd::data::particle_flags::isSet(
                particleAccessor->getFlags(idx),
                mesa_pd::data::particle_flags::GHOST)) continue;
         const auto uid = particleAccessor->getUid(idx);
         for (size_t bi = 0; bi < nB; ++bi)
         {
            if (floaterBoxUIDs_[bi] != uid) continue;
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
            const real_t hx = floaterBoxHalfExtents_[bi][0];
            const real_t hy = floaterBoxHalfExtents_[bi][1];
            const real_t hz = floaterBoxHalfExtents_[bi][2];
            // Wide lateral reach so the surface plane-fit samples PAST any crater/gas-moat the box digs,
            // anchoring the fitted waterline to the undisturbed surrounding surface (mirrors the sphere).
            const real_t band = real_t(0.8) * std::max({hx, hy, hz}) + real_t(2);

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
         const real_t hx    = floaterBoxHalfExtents_[bi][0];
         const real_t hy    = floaterBoxHalfExtents_[bi][1];
         const real_t hz    = floaterBoxHalfExtents_[bi][2];
         const real_t vFull = real_t(8) * hx * hy * hz;

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
            const real_t planeOffset = nHat * (xBar - pos); // submerged sample: nHat·off < planeOffset

            constexpr int N = 12; // 12³ = 1728 samples per box
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
               const Vector3< real_t > off = R * local; // offset from box centre (world frame)
               ++total;
               if ((nHat * off) < planeOffset) { ++below; cbSum += off; }
            }
            vSub = vFull * real_t(below) / real_t(total);
            if (below > uint_t(0)) cbOff = cbSum / real_t(below);
         }
         else
         {
            // No interface cells near the box: it is entirely in gas or entirely in liquid.
            const real_t cnt     = fluidCnt[bi];
            const real_t avgFill = (cnt > real_t(0)) ? acc[bi * stride + 7] / cnt : real_t(0);
            vSub  = (avgFill > real_t(0.5)) ? vFull : real_t(0);
            cbOff = Vector3< real_t >(real_t(0)); // symmetric → no righting lever
         }

         vSub = std::clamp(vSub, real_t(0), vFull);
         submergedVolumes_[floaterBoxUIDs_[bi]] = vSub;
         if (vSub > real_t(1e-9))
            buoyancyCenterOffset_[floaterBoxUIDs_[bi]] = cbOff;

         if (stepCounter_ % evalFrequency_ == uint_t(0))
         {
            // Y-tilt angle from rotation matrix: for R = Ry(θ), R(0,2)=sinθ, R(0,0)=cosθ.
            const real_t R00     = rotFlat[bi * 9 + 0];
            const real_t R02     = rotFlat[bi * 9 + 2];
            const real_t tiltDeg = std::atan2(R02, R00) * real_t(180) / math::pi;
            WALBERLA_LOG_INFO_ON_ROOT("[BoxBuoyancy] uid=" << floaterBoxUIDs_[bi]
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
               << "  comOffsetZ=" << (boxComOffsetZ_.count(floaterBoxUIDs_[bi])
                                         ? boxComOffsetZ_.at(floaterBoxUIDs_[bi]) : real_t(0)));
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
   std::vector< walberla::id_t >                      floaterBoxUIDs_;        // box floaters, same size on all ranks
   std::vector< Vector3< real_t > >                   floaterBoxHalfExtents_; // indexed by floaterBoxUIDs_
   std::unordered_map< walberla::id_t, real_t >       submergedVolumes_;    // keyed by UID; populated by both compute* methods
   std::unordered_map< walberla::id_t, Vector3< real_t > > buoyancyCenterOffset_; // box C_b − C_com (world frame); box-only, for righting torque
   std::unordered_map< walberla::id_t, real_t >       boxComOffsetZ_;       // per-box COM offset below centroid (body z, const); for the pendulum righting torque

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
   real_t                                             surfaceHeight_; // mean free-surface height (absoluteLiquidHeight); global V_sub bound
   std::shared_ptr< Vector3< real_t > >               acceleration_;  // body acceleration incl. gravity; for net-force logging
};

} // namespace antidunes
} // namespace walberla