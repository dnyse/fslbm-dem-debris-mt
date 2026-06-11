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

   // Box vs CylindricalBoundary (solid exterior) — AABB approach, approximation for small rotations.
   // Finds the closest point on the AABB to the cylinder axis, then checks whether the radial
   // distance from that point to the axis is less than the cylinder radius.
   //   penetrationDepth = radial_dist - cylinder_radius  (< 0 → penetrating)
   //   contactNormal points radially outward from the cylinder axis toward the box
   //   contactPoint = closest point on the AABB surface (on the box, nearest the cylinder)
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
      const mesa_pd::Vec3& axis      = cb.getAxis();

      // AABB extents (valid for identity-rotation FloatingBox)
      const mesa_pd::Vec3 half   = bx.getEdgeLength() * real_t(0.5);
      const mesa_pd::Vec3 boxMin = c - half;
      const mesa_pd::Vec3 boxMax = c + half;

      // Closest point on AABB to the cylinder center (in world frame)
      const mesa_pd::Vec3 closest(
         std::clamp(cylCenter[0], boxMin[0], boxMax[0]),
         std::clamp(cylCenter[1], boxMin[1], boxMax[1]),
         std::clamp(cylCenter[2], boxMin[2], boxMax[2]));

      // Radial vector from cylinder center to closest point, perpendicular to axis
      const mesa_pd::Vec3 vec    = closest - cylCenter;
      const mesa_pd::Vec3 radial = vec - (vec * axis) * axis;
      const real_t dist = radial.length();

      const real_t penetrationDepth = dist - cb.getRadius();
      if (penetrationDepth < getContactThreshold())
      {
         getPenetrationDepth() = penetrationDepth;
         getContactNormal() = (dist > real_t(1e-10))
                                 ? radial / dist
                                 : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         getContactPoint()  = closest;
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

      if (!std::isfinite(hydro_force[0]) || !std::isfinite(hydro_force[1]) || !std::isfinite(hydro_force[2]))
         WALBERLA_LOG_WARNING("[AddExternalForces] NaN/Inf hydrodynamic force on uid=" << ac.getUid(idx)
            << " pos=" << ac.getPosition(idx));
      if (!std::isfinite(hydro_torque[0]) || !std::isfinite(hydro_torque[1]) || !std::isfinite(hydro_torque[2]))
         WALBERLA_LOG_WARNING("[AddExternalForces] NaN/Inf hydrodynamic torque on uid=" << ac.getUid(idx)
            << " pos=" << ac.getPosition(idx));

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
      mesa_pd::addTorqueAtomic(idx, ac, totalTorque);
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain >   domain_;
   std::shared_ptr< Vector3< real_t > >                    acceleration_;
   real_t                                                   particleDensityRatio_;
   const std::unordered_map< walberla::id_t, real_t >*     submergedVolumes_;
   const std::unordered_map< walberla::id_t, Vector3< real_t > >* buoyancyCenterOffset_;
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
         floaterBoxHalfExtents_[bi] = simParams->floatingBoxes[bi].halfExtents;

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
      // Submerged volumes computed once per LBM step (fill field is constant during DEM sub-cycles).
      // Both functions write into the same submergedVolumes_ map (keyed by UID).
      computeSubmergedVolumes();
      computeSubmergedVolumesBox();
      addExternalForces.setSubmergedVolumes(&submergedVolumes_);
      addExternalForces.setBuoyancyCenterOffset(&buoyancyCenterOffset_);

      for (auto subCycle = uint_t(0); subCycle < numSubCycles; ++subCycle)
      {
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
      constexpr size_t stride = 9;
      std::vector< real_t > accum(nF * stride, real_t(0));

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
            const real_t r        = floaterRadii_[fi];
            const int    pad      = static_cast< int >(std::ceil(r)) + 2;
            const real_t ringBand = ringBand_; // contact-ring half-width (prm: floaterRingBand; large ⇒ old full-box)

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

               // Fallback channel: accumulate fill of non-obstacle cells in the upper
               // half of the bounding box (z ≥ sphere centre).  Pillars, NoSlip walls,
               // and any other obstacle cells are skipped by the fluidMask guard, so
               // the search automatically continues past them.  Multiple ranks may each
               // contribute cells; since all valid cells have the same fill (≈0 or ≈1)
               // in the fallback case, averaging across ranks is correct.
               if (isPartOfMaskSet(cellFlag, fluidMask) && samplePos[2] >= c[2])
               {
                  accum[fi * stride + 7] += fillField->get(localCell[0], localCell[1], localCell[2]);
                  accum[fi * stride + 8] += real_t(1);
               }

               // Plane-fitting channel: only interface cells have valid normals.
               if (!isFlagSet(cellFlag, interfaceFlag)) continue;

               // Contact-ring restriction: keep only interface cells on the sphere SURFACE.
               // Distance is computed from the loop offsets relative to the (periodic-mapped)
               // centre c, so it is periodic-safe (no wrap across the seam). Excludes the
               // far-field surface that would otherwise bias the plane under waves.
               const real_t offX  = real_t(dx) + real_t(0.5) - (c[0] - real_t(cxi));
               const real_t offY  = real_t(dy) + real_t(0.5) - (c[1] - real_t(cyi));
               const real_t offZ  = real_t(dz) + real_t(0.5) - (c[2] - real_t(czi));
               const real_t distC = std::sqrt(offX * offX + offY * offY + offZ * offZ);
               if (std::abs(distC - r) > ringBand) continue;

               const real_t phi = fillField->get(localCell[0], localCell[1], localCell[2]);
               const real_t w   = phi * (real_t(1) - phi);
               if (w < real_t(1e-6)) continue;

               const Vector3< real_t >& n = normalField->get(localCell[0], localCell[1], localCell[2]);
               if (n.sqrLength() < real_t(1e-10)) continue;

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
         const real_t wTot  = accum[fi * stride + 0];
         real_t h_cap;

         // --- Live LOCAL free-surface height at the sphere ---
         // accum[..+6]/wTot is the gradient-weighted mean z of the interface cells sitting on the
         // sphere SURFACE (the contact ring), i.e. the actual waterline AROUND the sphere. Because the
         // sphere dams the driven flow and displaces water, this local waterline is a crater/ridge that
         // sits well above (or below) the undisturbed mean level absoluteLiquidHeight. Using it instead
         // of the fixed mean makes submersion (and the bound below) self-consistent with the deformed
         // surface, so a sphere floating in a raised/dimpled surface still reports ~densityRatio×100%.
         // Falls back to the fixed mean only when no contact ring is found (fully submerged/emerged),
         // where that flat global value is the correct cavity-proof bound.
         const real_t localSurfaceZ =
            (wTot > real_t(1e-6)) ? (accum[fi * stride + 6] / wTot) : surfaceHeight_;
         // surfRef = max(fixed, local) keeps the bound CAVITY-PROOF: a hugging gas pocket drags the
         // local estimate down toward the centre, but the max holds the reference at the far-field
         // level so a deeply-submerged sphere still reads fully submerged (the original sink fix),
         // while still following a genuinely risen/displaced surface.
         const real_t surfRef = std::max(surfaceHeight_, localSurfaceZ);

         // --- Local free-surface bound (robust to gas cavities / spurious interface cells) ---
         // BUG FIX: the old local-fill fallback under-reported V_sub when a deeply-submerged sphere sat
         // in an entrained gas cavity (it averaged the surrounding cells, which were partly gas), so the
         // buoyancy collapsed and the sphere stayed sunk. Using the surface height instead:
         //   sphere top  (C_z + r) <= surface  ->  fully submerged (h_cap = 2r), ignore any local gas
         //   sphere bottom (C_z - r) >= surface ->  fully emerged (h_cap = 0)
         // Only a sphere that genuinely STRADDLES the surface uses the local contact-line plane-fit
         // (which captures a tilted / wavy waterline).
         const real_t topZ = C[2] + r;
         const real_t botZ = C[2] - r;

         if (topZ <= surfRef)
         {
            h_cap = real_t(2) * r; // fully below the surface -> fully submerged (cavity-proof)
         }
         else if (botZ >= surfRef)
         {
            h_cap = real_t(0);     // entirely above the surface
         }
         else if (wTot < real_t(1e-6))
         {
            // Straddling but no contact ring found -> geometric cap from the flat global surface.
            h_cap = std::clamp(surfRef - botZ, real_t(0), real_t(2) * r);
         }
         else
         {
            // Straddling WITH a contact ring -> plane-fit (handles a tilted / locally wavy waterline).
            const Vector3< real_t > nSum(accum[fi*stride+1], accum[fi*stride+2], accum[fi*stride+3]);
            const Vector3< real_t > xBar(accum[fi*stride+4] / wTot,
                                         accum[fi*stride+5] / wTot,
                                         accum[fi*stride+6] / wTot);

            if (nSum.sqrLength() < real_t(1e-10))
            {
               // Weighted normals cancelled -> fall back to the vertical distance to the centroid.
               const real_t h_vertical = xBar[2] - C[2];
               h_cap = std::clamp(r + h_vertical, real_t(0), real_t(2) * r);
            }
            else
            {
               // n̂ points liquid->gas; h_signed > 0 when centre C is below the surface.
               const Vector3< real_t > nHat = nSum / nSum.length();
               const real_t h_signed = nHat * (xBar - C);
               h_cap = std::clamp(r + h_signed, real_t(0), real_t(2) * r);
            }
         }

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
               << "  localSurfaceZ=" << localSurfaceZ << " (live; fixed mean=" << surfaceHeight_ << ")"
               << "  centreDepthBelowSurface=" << centreDepth << " LU"
               << "  (expect submersion ~densityRatio×100% and centreDepth~0 at equilibrium)");
         }
      }
      ++stepCounter_;
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

      const real_t band = ringBand_ + real_t(2); // shell thickness (cells) outside the box surface

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
               if (n.sqrLength() < real_t(1e-10)) continue;

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

         // Geometric world-space half-height of the rotated box (NO shell band).
         const Vector3< real_t > axx = R * Vector3< real_t >(hx, real_t(0), real_t(0));
         const Vector3< real_t > ayy = R * Vector3< real_t >(real_t(0), hy, real_t(0));
         const Vector3< real_t > azz = R * Vector3< real_t >(real_t(0), real_t(0), hz);
         const real_t boxHalfZ = std::abs(axx[2]) + std::abs(ayy[2]) + std::abs(azz[2]);
         const real_t topZ = pos[2] + boxHalfZ;
         const real_t botZ = pos[2] - boxHalfZ;

         // Live local waterline at the box: gradient-weighted mean z of the shell interface cells,
         // falling back to the fixed mean surface height. surfRef = max(fixed, local) is CAVITY-PROOF:
         // an entrained gas pocket hugging a deeply-submerged box drags the local estimate down toward
         // the box centre, but the max keeps the reference at the far-field level so the box still reads
         // fully submerged. Mirrors the floating-sphere fix that stopped the sink-to-floor.
         const real_t localSurfaceZ =
            (wTot > real_t(1e-6)) ? (acc[bi * stride + 6] / wTot) : surfaceHeight_;
         const real_t surfRef = std::max(surfaceHeight_, localSurfaceZ);

         real_t            vSub = real_t(0);
         Vector3< real_t > cbOff(real_t(0));

         if (topZ <= surfRef)
         {
            vSub  = vFull;                        // fully below the surface -> fully submerged (cavity-proof)
            cbOff = Vector3< real_t >(real_t(0)); // B at geometric centre -> no righting lever
         }
         else if (botZ >= surfRef)
         {
            vSub  = real_t(0);                    // entirely above the surface
            cbOff = Vector3< real_t >(real_t(0));
         }
         else if (wTot > real_t(1e-6))
         {
            // Straddling WITH interface cells -> plane-fit + clip (handles a tilted / wavy waterline).
            // n̂ points liquid→gas; degenerate (cancelled) normals fall back to vertical (horizontal surface).
            const Vector3< real_t > nSum(acc[bi*stride+1], acc[bi*stride+2], acc[bi*stride+3]);
            const Vector3< real_t > nHat = (nSum.sqrLength() > real_t(1e-10))
                                              ? nSum / nSum.length()
                                              : Vector3< real_t >(real_t(0), real_t(0), real_t(1));
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
               << "  localSurfaceZ=" << localSurfaceZ << " (live; fixed mean=" << surfaceHeight_ << ")"
               << "  B-G offset=(" << cbOff[0] << "," << cbOff[1] << "," << cbOff[2] << ")"
               << "  |B-G|_horiz=" << std::sqrt(cbOff[0] * cbOff[0] + cbOff[1] * cbOff[1]));
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