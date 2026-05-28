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
#include "mesa_pd/data/Flags.h"
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
   {}

   void setSubmergedVolumes(const std::unordered_map< walberla::id_t, real_t >* sv) {
      submergedVolumes_ = sv;
   }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac){
      // Buoyancy acts only in Z because the hydrostatic MEM contribution is only in Z.
      // Type 1 = floating sphere; derive its density ratio from invMass×volume (baked in at creation).
      const real_t rho = (ac.getType(idx) == uint_t(1))
         ? real_t(1) / (ac.getInvMass(idx) * ac.getVolume(idx))
         : particleDensityRatio_;

      const real_t vFull = ac.getVolume(idx);
      // For floating spheres use the precomputed submerged volume so that the gas-phase portion
      // of the sphere does not contribute a spurious upward buoyancy force.
      // Non-floaters are always fully submerged, so vSub = vFull leaves their behaviour unchanged.
      real_t vSub = vFull;
      if (ac.getType(idx) == uint_t(1) && submergedVolumes_) {
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
      mesa_pd::addTorqueAtomic(idx, ac, hydro_torque);
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain >   domain_;
   std::shared_ptr< Vector3< real_t > >                    acceleration_;
   real_t                                                   particleDensityRatio_;
   const std::unordered_map< walberla::id_t, real_t >*     submergedVolumes_;
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
         uint_t(3))
      , addExternalForces(
         simObjs->mesapdDomain,
         simObjs->acceleration,
         simParams->particle.densityRatio)
      , useLubricationCorrection(simParams->particle.useLubricationCorrection)
      , evalFrequency_(simParams->evaluation.evaluationFrequency)
      , stepCounter_(uint_t(0))
   {
      // Build floater UID/radius list with the same size on ALL ranks so that
      // allReduceInplace calls in computeSubmergedVolumes are always MPI-valid.
      // Only the owner rank (non-ghost) contributes each UID; SUM broadcasts it.
      const size_t nF = simParams->floatingSpheres.size();
      floaterUIDs_.assign(nF, walberla::id_t(0));
      floaterRadii_.resize(nF);
      for (size_t fi = 0; fi < nF; ++fi)
         floaterRadii_[fi] = simParams->floatingSpheres[fi].radius;

      uint_t assignIdx = 0;
      for (size_t idx = 0; idx < particleAccessor->size(); ++idx)
      {
         if (particleAccessor->getType(idx) != uint_t(1)) continue;
         if (mesa_pd::data::particle_flags::isSet(
                particleAccessor->getFlags(idx),
                mesa_pd::data::particle_flags::GHOST)) continue;
         if (assignIdx < nF)
            floaterUIDs_[assignIdx++] = particleAccessor->getUid(idx);
      }
      mpi::allReduceInplace(floaterUIDs_, mpi::SUM);
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
      // Submerged volume is computed once per LBM step (fill field does not change during DEM sub-cycles).
      computeSubmergedVolumes();
      addExternalForces.setSubmergedVolumes(&submergedVolumes_);

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
   void computeSubmergedVolumes()
   {
      const size_t nF = floaterUIDs_.size();
      if (nF == 0) { submergedVolumes_.clear(); return; }

      // Gather current positions from owner ranks (skip ghost copies so each
      // sphere contributes exactly once; allReduce SUM then broadcasts to all ranks).
      // Using floaterUIDs_.size() as nF guarantees all 72 ranks call allReduceInplace
      // with the same vector length — avoiding the undefined MPI count mismatch.
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
            cenX[fi] = pos[0];  cenY[fi] = pos[1];  cenZ[fi] = pos[2];
            break;
         }
      }
      mpi::allReduceInplace(cenX, mpi::SUM);
      mpi::allReduceInplace(cenY, mpi::SUM);
      mpi::allReduceInplace(cenZ, mpi::SUM);

      // Build per-floater slice offsets (each floater has 2*r_cells+1 Z-slices).
      // Track raw (numerator, denominator) per slice so that ring cells split across
      // block boundaries are averaged correctly after the global allReduce.
      std::vector< size_t > rCells(nF);
      std::vector< size_t > sliceOffsets(nF + 1, size_t(0));
      for (size_t fi = 0; fi < nF; ++fi)
      {
         rCells[fi]           = static_cast< size_t >(std::ceil(floaterRadii_[fi]));
         sliceOffsets[fi + 1] = sliceOffsets[fi] + size_t(2) * rCells[fi] + size_t(1);
      }
      const size_t totalSlices = sliceOffsets[nF];

      std::vector< real_t > sliceNum(totalSlices, real_t(0)); // sum of phi per slice
      std::vector< real_t > sliceDen(totalSlices, real_t(0)); // cell count per slice

      for (auto& iblock : *blockForest_)
      {
         auto* fillField = iblock.getData< ScalarField_T >(fillFieldID_);
         const CellInterval& blockCellBB = blockForest_->getBlockCellBB(iblock);

         for (size_t fi = 0; fi < nF; ++fi)
         {
            const real_t   r  = floaterRadii_[fi];
            const int      rc = static_cast< int >(rCells[fi]);
            Vector3< real_t > c(cenX[fi], cenY[fi], cenZ[fi]);
            blockForest_->mapToPeriodicDomain(c);

            for (int dz = -rc; dz <= rc; ++dz)
            {
               const real_t dz_r      = static_cast< real_t >(dz);
               const real_t r_ring_sq = r * r - dz_r * dz_r;
               if (r_ring_sq < real_t(0)) continue;

               const int    r_outer = static_cast< int >(std::sqrt(r_ring_sq)) + 1;
               const real_t sz      = c[2] + dz_r;

               if (sz < blockForest_->getDomain().zMin() ||
                   sz >= blockForest_->getDomain().zMax()) continue;

               const size_t sliceIdx = sliceOffsets[fi] + static_cast< size_t >(dz + rc);

               const std::array< std::pair< int, int >, 4 > dirs = {{
                  {r_outer, 0}, {-r_outer, 0}, {0, r_outer}, {0, -r_outer}
               }};

               for (auto [ddx, ddy] : dirs)
               {
                  Vector3< real_t > samplePos(c[0] + static_cast< real_t >(ddx),
                                              c[1] + static_cast< real_t >(ddy), sz);
                  blockForest_->mapToPeriodicDomain(samplePos);

                  const Cell globalCell = blockForest_->getCell(samplePos);
                  if (!blockCellBB.contains(globalCell)) continue;

                  Cell localCell;
                  blockForest_->transformGlobalToBlockLocalCell(localCell, iblock, globalCell);
                  if (!fillField->xyzSize().contains(localCell)) continue;

                  sliceNum[sliceIdx] += fillField->get(localCell[0], localCell[1], localCell[2]);
                  sliceDen[sliceIdx] += real_t(1);
               }
            }
         }
      }

      mpi::allReduceInplace(sliceNum, mpi::SUM);
      mpi::allReduceInplace(sliceDen, mpi::SUM);

      // Compute phiSum = sum_dz( average_phi(dz) ) ≈ submersion depth h in LU
      std::vector< real_t > phiSum(nF, real_t(0));
      for (size_t fi = 0; fi < nF; ++fi)
      {
         const int rc2 = static_cast< int >(rCells[fi]);
         for (int dz = -rc2; dz <= rc2; ++dz)
         {
            const size_t sliceIdx = sliceOffsets[fi] + static_cast< size_t >(dz + rc2);
            if (sliceDen[sliceIdx] > real_t(0))
               phiSum[fi] += sliceNum[sliceIdx] / sliceDen[sliceIdx];
         }
      }

      submergedVolumes_.clear();
      for (size_t fi = 0; fi < nF; ++fi)
      {
         const real_t r     = floaterRadii_[fi];
         const real_t h     = std::clamp(phiSum[fi], real_t(0), real_t(2) * r);
         const real_t vSub  = math::pi * h * h * (real_t(3) * r - h) / real_t(3);
         const real_t vFull = (real_t(4) / real_t(3)) * math::pi * r * r * r;

         WALBERLA_ASSERT_GREATER_EQUAL(h, real_t(0));
         WALBERLA_ASSERT_LESS_EQUAL(h, real_t(2) * r);
         WALBERLA_ASSERT_GREATER_EQUAL(vSub, real_t(0));
         WALBERLA_ASSERT_LESS_EQUAL(vSub, vFull);

         submergedVolumes_[floaterUIDs_[fi]] = vSub;

         if (stepCounter_ % evalFrequency_ == uint_t(0))
         {
            WALBERLA_LOG_INFO_ON_ROOT("[FloaterBuoyancy] uid=" << floaterUIDs_[fi]
               << "  phiSum=" << phiSum[fi]
               << "  h=" << h << " LU"
               << "  V_sub=" << vSub << " / " << vFull << " LU³"
               << "  submersion=" << real_t(100) * vSub / vFull << "%"
               << "  (expect ~densityRatio×100% at equilibrium)");
         }
      }
      ++stepCounter_;
   }

private:
   std::shared_ptr< ParticleAccessor_T >              particleAccessor;
   std::shared_ptr< mesa_pd::data::ParticleStorage >  particleStorage;
   std::function< void() >                            syncCall;
   bool useOpenMP{ false };

   shared_ptr< StructuredBlockForest >                blockForest_;
   BlockDataID                                        fillFieldID_;
   std::vector< walberla::id_t >                      floaterUIDs_;   // same size on all ranks
   std::vector< real_t >                              floaterRadii_;  // indexed by floaterUIDs_
   std::unordered_map< walberla::id_t, real_t >       submergedVolumes_;

   // kernels — declaration order must match initializer-list order to avoid -Wreorder
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce;
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce;
   ComputeLubrication computeLubrication;
   ComputeContacts computeContacts;
   AddExternalForces addExternalForces;
   bool useLubricationCorrection;
   uint_t                                             evalFrequency_;
   uint_t                                             stepCounter_;
};

} // namespace antidunes
} // namespace walberla