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

#include "SimulationParams.h" 
#include "SimulationObjects.h" 


namespace walberla {
namespace antidunes {

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
      real_t collisionTime, real_t particleTimestepSize, real_t frictionCoefficient)
      : domain_(domain)
      , collisionResponse_(1)
      , restitutionCoefficient_(restitutionCoefficient)
      , kappa_(kappa)
      , collisionTime_(collisionTime)
      , particleTimestepSize_(particleTimestepSize)
   {
      collisionResponse_.setFrictionCoefficientDynamic(size_t(0), size_t(0), frictionCoefficient);
   }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx1, const size_t idx2, ParticleAccessor_T& ac){
      if (double_cast_(idx1, idx2, ac, acd_, ac))
      {
         if (contact_filter_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(), *domain_))
         {
            auto meff = real_c(1) / (ac.getInvMass(idx1) + ac.getInvMass(idx2));
            collisionResponse_.setStiffnessAndDamping(0, 0, restitutionCoefficient_, collisionTime_, kappa_, meff);
            collisionResponse_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(), acd_.getContactNormal(),
                              acd_.getPenetrationDepth(), particleTimestepSize_);
         }
      }
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain >    domain_;
   mesa_pd::kernel::LinearSpringDashpot                     collisionResponse_;
   mesa_pd::collision_detection::AnalyticContactDetection   acd_;
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
      const std::shared_ptr< Vector3< real_t > >& acceleration, real_t particleDensityRatio)
      : domain_(domain)
      , acceleration_(acceleration)
      , particleDensityRatio_(particleDensityRatio)
   {}

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac){
      // Buoyancy is removed only from the Z component because the hydrostaic contrbution removed during MEM is only applied in Z 
      Vector3< real_t > bodyForce =  ac.getVolume(idx) * Vector3< real_t >(
         (*acceleration_)[0], (*acceleration_)[1], (particleDensityRatio_ - real_c(1)) * (*acceleration_)[2]);

      mesa_pd::addForceAtomic(idx, ac, ac.getHydrodynamicForce(idx) + bodyForce);
      mesa_pd::addTorqueAtomic(idx, ac, ac.getHydrodynamicTorque(idx));
   }

private:
   std::shared_ptr< mesa_pd::domain::BlockForestDomain > domain_;
   std::shared_ptr< Vector3< real_t > > acceleration_;
   real_t particleDensityRatio_;
};



template< typename ParticleAccessor_T >
class ParticleKernels{
public:
   ParticleKernels(const std::shared_ptr< const SimulationParams >& simParams, 
      const std::shared_ptr< SimulationObjects >& simObjs) 
      : particleAccessor(simObjs->particleAccessor)
      , particleStorage(simObjs->particleStorage)
      , syncCall(simObjs->syncCall)
      , vvIntegratorPreForce(simParams->computed.particleTimestepSize)
      , vvIntegratorPostForce(simParams->computed.particleTimestepSize)
      , computeLubrication(simObjs->mesapdDomain, simParams->computed.viscosity)
      , computeContacts(
         simObjs->mesapdDomain,
         simParams->particle.restitutionCoefficient, 
         simParams->computed.particleKappa, 
         simParams->computed.particleCollisionTime, 
         simParams->computed.particleTimestepSize, 
         simParams->particle.frictionCoefficient)
      , addExternalForces(
         simObjs->mesapdDomain, 
         simObjs->acceleration, 
         simParams->particle.densityRatio)
      , useLubricationCorrection(simParams->particle.useLubricationCorrection)
   {}

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

private:
   std::shared_ptr< ParticleAccessor_T >              particleAccessor;
   std::shared_ptr< mesa_pd::data::ParticleStorage >  particleStorage;
   std::function< void() >                            syncCall;
   bool useOpenMP{ false };
   
   // kernels -------------------------------------------------------------
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce;
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce;
   ComputeLubrication computeLubrication;
   ComputeContacts computeContacts;
   AddExternalForces addExternalForces;
   bool useLubricationCorrection;
};

} // namespace antidunes
} // namespace walberla