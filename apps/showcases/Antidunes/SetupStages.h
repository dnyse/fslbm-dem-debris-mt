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
//! \file SetupStages.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Jonas Plewinski <jonas.plewinski@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/Initialization.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <vector>

#include "lbm/field/AddToStorage.h"

#include "mesa_pd/data/shape/Capsule.h"
#include "mesa_pd/data/shape/Cylinder.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/mpi/ClearGhostOwnerSync.h"
#include "mesa_pd/mpi/ClearNextNeighborSync.h"
#include "mesa_pd/mpi/SyncGhostOwners.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"

#include "DataInitialization.h"
#include "SimulationParams.h"
#include "SimulationObjects.h"
#include "Utility.h"


namespace walberla {
namespace antidunes {

void setupBlockForest(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   const auto& bfParams = simParams->blockForest;
   const auto& domParams  = simParams->domain;

   Vector3< uint_t > numBlocks;
   for (uint_t i = 0; i <= 2; ++i) {
      numBlocks[i] = domParams.domainSize[i] / bfParams.cellsPerBlock[i];
      WALBERLA_CHECK_EQUAL(numBlocks[i] * bfParams.cellsPerBlock[i], domParams.domainSize[i],
         "Domain size in direction " << i << " is not a multiple of cells per block.")

      if(bfParams.periodicity[i]){
         WALBERLA_CHECK_GREATER_EQUAL(numBlocks[i], uint_c(3),
         "When using particles, use at least three blocks per periodic direction.")
      }

   }

   const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   WALBERLA_CHECK_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
      "The number of MPI processes is different from the number of blocks as defined by "
      "\"domainSize/cellsPerBlock\".")

   if (bfParams.loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Loading blockForest from snapshot");
      MPIManager::instance()->useWorldComm();

      simObjs->blockForest = make_shared< StructuredBlockForest >(
            make_shared< BlockForest >(
               uint_c(MPIManager::instance()->rank()),
               (bfParams.snapshotBaseFolder + "/" + simObjs->blockForestFile).c_str(),
               true,
               false
            ),
            bfParams.cellsPerBlock[0], bfParams.cellsPerBlock[1], bfParams.cellsPerBlock[2]);

      simObjs->blockForest->createCellBoundingBoxes();

   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Creating blockForest");
      simObjs->blockForest = blockforest::createUniformBlockGrid(
         numBlocks[0], numBlocks[1], numBlocks[2],
         bfParams.cellsPerBlock[0], bfParams.cellsPerBlock[1], bfParams.cellsPerBlock[2],
         real_c(1.0),   // dx
         true,          // one block per process
         bfParams.periodicity[0], bfParams.periodicity[1], bfParams.periodicity[2]);
   }

   // save block forest to file (but do not overwrite existing file if snapshot is loaded)
   if (bfParams.storeSnapshot && !bfParams.loadSnapshot){
      WALBERLA_LOG_INFO_ON_ROOT("Saving blockForest to snapshot");
      simObjs->blockForest->getBlockForest().saveToFile(bfParams.snapshotBaseFolder + "/" + simObjs->blockForestFile);
   }
}

void setupSimState(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("setupSimState");
   const auto& bfParams = simParams->blockForest;

   simObjs->beginTimeStep = uint_c(0);

   const real_t accelerationX = simParams->pid.accelerationXOverride > real_c(0)
                                 ? simParams->pid.accelerationXOverride
                                 : real_c(3) * simParams->pid.targetMeanVelocityMagnitude * simParams->computed.viscosity /
                                   (simParams->computed.liquidHeight * simParams->computed.liquidHeight);
   simObjs->acceleration = std::make_shared< Vector3< real_t > >(accelerationX, real_c(0.0), -simParams->computed.gravity);

   simObjs->pidController = std::make_shared< PIDController >(
                                    simParams->pid.targetMeanVelocityMagnitude,
                                    accelerationX,
                                    simParams->pid.proportionalGain,
                                    simParams->pid.derivativeGain,
                                    simParams->pid.integralGain,
                                    simParams->pid.maxRamp,
                                    simParams->pid.minActuatingVariable,
                                    simParams->pid.maxActuatingVariable);

   // meanVelocity and totalFluidMass start at zero; they are computed in the evaluation stages
   simObjs->meanVelocity = std::make_shared< Vector3< real_t > >(real_c(0));
   simObjs->totalFluidMass = std::make_shared< real_t >(real_c(0));
   simObjs->dataCols = std::make_shared< DataCols >();

   if (bfParams.loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Loading simState and pidState from snapshot");

      WALBERLA_ROOT_SECTION()
      {
         std::ifstream file;
         file.open(bfParams.snapshotBaseFolder + "/" + simObjs->simStateFile);
         if (file.fail()) WALBERLA_ABORT("Error: " << simObjs->simStateFile << " could not be opened!");
         file >> simObjs->beginTimeStep;
         file >> (*simObjs->acceleration)[0];
         file.close();
      }
      mpi::broadcastObject(simObjs->beginTimeStep);
      mpi::broadcastObject(*simObjs->acceleration);

      simObjs->pidController->readStateFromFile(bfParams.snapshotBaseFolder + "/" + simObjs->pidStateFile);
   }

}

void setupVTKState(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   if (simParams->blockForest.loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Modifying config file to start VTK counter from beginTimeStep rather than from 0");
      std::vector< config::Config::Block* > configVTKBlock;
      simParams->cfgPtr->getWritableGlobalBlock().getWritableBlocks("VTK", configVTKBlock, 1, 1);
      std::vector< config::Config::Block* > configVTKFluidFieldBlock;
      configVTKBlock[0]->getWritableBlocks("fluid_field", configVTKFluidFieldBlock, 1, 1);
      configVTKFluidFieldBlock[0]->setOrAddParameter("initialExecutionCount", std::to_string(simObjs->beginTimeStep));
   }
}

void setupLBMFields(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   auto& bf  = simObjs->blockForest;
   const auto& bfParams = simParams->blockForest;
   const walberla::free_surface::ExcessMassDistributionModel excessMassModel(simParams->model.excessMassDistributionModel);

   simObjs->forceDensityFieldID = field::addToStorage< VectorFieldFlattened_T >(bf, "Force field", real_c(0), field::fzyx, uint_c(1));
   LatticeModel_T latticeModel(simObjs->forceDensityFieldID, simParams->computed.relaxationRate);


   if (bfParams.loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Loading LBM fields from snapshot");
      const std::string base = bfParams.snapshotBaseFolder + "/";
      auto pdfHandling = make_shared< lbm::internal::PdfFieldHandling< LatticeModel_T > >(bf, latticeModel, false, Vector3< real_t >(real_c(0)), real_c(1), uint_c(1), field::fzyx);
      simObjs->pdfFieldID = bf->getBlockStorage().loadBlockData(base + simObjs->pdfFieldFile, pdfHandling, "PDF field");

      auto fillHandling = std::make_shared< ScalarFieldHandling< ScalarField_T > >(bf, uint_c(2));
      simObjs->fillFieldID = bf->getBlockStorage().loadBlockData(base + simObjs->fillFieldFile, fillHandling, "Fill level field");

      if (excessMassModel.isEvenlyAllInterfaceFallbackLiquidType())
      {
         // load additional field that stores excess mass in liquid cells
         auto excessHandling = std::make_shared< ScalarFieldHandling< ScalarField_T > >(bf, uint_c(1));
         simObjs->excessMassFieldID = bf->getBlockStorage().loadBlockData( base + simObjs->excessMassFieldFile, excessHandling, "Excess mass field");
      }
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Creating LBM fields");
      auto initialVelocity = Vector3< real_t >(simParams->pid.targetMeanVelocityMagnitude, 0, 0);
      simObjs->pdfFieldID = lbm::addPdfFieldToStorage( bf, "PDF field", latticeModel, initialVelocity, real_c(1.0), field::fzyx);
      simObjs->fillFieldID =  field::addToStorage< ScalarField_T >(bf, "Fill level field", real_c(0), field::fzyx, uint_c(2));

      if (excessMassModel.isEvenlyAllInterfaceFallbackLiquidType())
      {
         // add additional field for storing excess mass in liquid cells
         simObjs->excessMassFieldID = field::addToStorage< ScalarField_T >(bf, "Excess mass field", real_c(0), field::fzyx, uint_c(1));
      }

      initFillFieldSinus(simParams, simObjs);
      initPdfFieldPoiseuille(simParams, simObjs);
      initPdfFieldHydrostaticPressure(simParams, simObjs);
   }
}

void setupSurfaceGeometryFields(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   auto& bf  = simObjs->blockForest;

   simObjs->curvatureFieldID =      field::addToStorage< ScalarField_T >(bf, "Curvature field", real_c(0), field::fzyx, uint_c(1));

   if (simParams->physics.enableWetting)
   {
      WALBERLA_ABORT("enableWetting cannot be enabled (TODO: use GenerticSurfaceGeometryHandler)")
   }

   // obstacleNormalField is filled by ObstacleNormalSweep every timestep (not gated on enableWetting)
   simObjs->normalFieldID =         field::addToStorage< VectorField_T >(bf, "Normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   simObjs->obstacleNormalFieldID = field::addToStorage< VectorField_T >(bf, "Obstacle normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   simObjs->smoothFillFieldID =     field::addToStorage< ScalarField_T >(bf, "Smooth fill level field", real_c(0), field::fzyx, uint_c(1));
}

void setupParticles(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   const auto& bfParams = simParams->blockForest;
   auto& bf             = simObjs->blockForest;

   simObjs->particleStorage  = std::make_shared< mesa_pd::data::ParticleStorage >(1);
   simObjs->particleAccessor = std::make_shared< ParticleAccessor_T >(simObjs->particleStorage);
   simObjs->mesapdDomain     = std::make_shared< mesa_pd::domain::BlockForestDomain >(bf->getBlockForestPointer());

   uint_t numParticles = uint_c(0);
   uint_t numFixedParticles = uint_c(0);


   if (bfParams.loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Loading particles from snapshot");
      simObjs->particleStorageID = bf->loadBlockData(bfParams.snapshotBaseFolder + "/" + simObjs->particleStorageFile,
         mesa_pd::domain::createBlockForestDataHandling(simObjs->particleStorage), "Particle Storage");

      mesa_pd::mpi::ClearNextNeighborSync{}(*simObjs->particleAccessor);
      mesa_pd::mpi::ClearGhostOwnerSync{}(*simObjs->particleAccessor);

      // count and validate particles loaded from snapshot ---------------------------------------
      // The Sphere cast and the diameter check only apply to a checkpointed sphere bed (useParticleBed=true).
      // In floater-only runs the storage also holds boxes/cylinders, for which static_cast<Sphere*> is
      // undefined behaviour, so particles are only counted.
      if (simParams->particle.useParticleBed)
      {
         real_t avgParticleDiameterTest = real_c(0);
         simObjs->particleStorage->forEachParticle(
            false, mesa_pd::kernel::SelectLocal(), *simObjs->particleAccessor,
            [&avgParticleDiameterTest, &numParticles, &numFixedParticles](const size_t idx, auto& ac) {
               auto sp = static_cast< mesa_pd::data::Sphere* >(ac.getBaseShape(idx).get());
               avgParticleDiameterTest += real_c(2) * sp->getRadius();
               ++numParticles;
               if (!mesa_pd::data::particle_flags::isSet( ac.getFlags(idx), mesa_pd::data::particle_flags::FIXED))
               {
                  ++numFixedParticles;
               }
            },
            *simObjs->particleAccessor);

         mpi::allReduceInplace(numParticles, mpi::SUM);
         mpi::allReduceInplace(avgParticleDiameterTest, mpi::SUM);
         avgParticleDiameterTest /= real_c(numParticles);

         WALBERLA_LOG_INFO_ON_ROOT("Read particles from check pointing file with avg diameter of " << avgParticleDiameterTest)
         if (std::abs(avgParticleDiameterTest - simParams->computed.avgParticleDiameter) / avgParticleDiameterTest > real_c(0.05))
         {
            WALBERLA_ABORT("Particle diameters not correct.")
         }
      }
      else
      {
         simObjs->particleStorage->forEachParticle(
            false, mesa_pd::kernel::SelectLocal(), *simObjs->particleAccessor,
            [&numParticles, &numFixedParticles](const size_t idx, auto& ac) {
               ++numParticles;
               if (!mesa_pd::data::particle_flags::isSet( ac.getFlags(idx), mesa_pd::data::particle_flags::FIXED))
               {
                  ++numFixedParticles;
               }
            },
            *simObjs->particleAccessor);
         mpi::allReduceInplace(numParticles, mpi::SUM);
         WALBERLA_LOG_INFO_ON_ROOT("Read " << numParticles << " particles from checkpointing file (useParticleBed=false, no diameter check).")
      }
   }
   else
   {
      simObjs->particleStorageID = bf->addBlockData(mesa_pd::domain::createBlockForestDataHandling(simObjs->particleStorage), "Particle Storage");

      // Only create particles below the initialization profile and outside any pillar
      const real_t maxSphereRadius = real_t(0.5) * simParams->computed.avgParticleDiameter * real_t(1.1);

      auto createParticleFct = [simParams, maxSphereRadius](Vector3< real_t > pos) {
         if (pos[2] >= initializationProfile(pos[0], simParams->computed.sinusAmplitude, simParams->computed.sinusOffset, simParams->computed.sinusWavelength))
            return false;
         // exclude particles overlapping the Y walls when Y is non-periodic
         if (!simParams->blockForest.periodicity[1])
         {
            if (pos[1] < maxSphereRadius || pos[1] > real_c(simParams->domain.domainSize[1]) - maxSphereRadius)
               return false;
         }
         // box pillars: nearest-point distance in XY
         for (const auto& pillar : simParams->pillars)
         {
            const real_t dx = std::max(real_t(0), std::abs(pos[0] - pillar.center[0]) - pillar.halfExtents[0]);
            const real_t dy = std::max(real_t(0), std::abs(pos[1] - pillar.center[1]) - pillar.halfExtents[1]);
            if (std::sqrt(dx * dx + dy * dy) < maxSphereRadius)
               return false;
         }
         // cylinder pillars: radial distance in XY
         for (const auto& pillar : simParams->cylinderPillars)
         {
            const real_t dx = pos[0] - pillar.center[0];
            const real_t dy = pos[1] - pillar.center[1];
            if (std::sqrt(dx * dx + dy * dy) < pillar.radius + maxSphereRadius)
               return false;
         }
         return true;
      };

      if (simParams->particle.useParticleBed)
      {
         real_t maxParticleHeight = real_c(0); // top of the sphere (i.e. maxPos + its radius)

         WALBERLA_LOG_INFO_ON_ROOT("Initializing particles from bedgeneration output file : " << simParams->particle.inFileName);
         initSpheresFromFile(simParams->particle.inFileName, *simObjs->particleStorage, *simObjs->mesapdDomain,
                             simParams->particle.densityRatio, simParams->domain.domainSize, createParticleFct,
                             bf->getDomain(), simParams->particle.bedCopiesInX, simParams->particle.bedCopiesInY,
                             numParticles, maxParticleHeight, simParams->computed.particleScalingFactor);

         WALBERLA_LOG_INFO_ON_ROOT("Max particle height " << maxParticleHeight);
         // These checks assume active DEM dynamics: with useNoSlipParticles=true, fixingHeightFactor is set large
         // to fix all particles, so particleFixingHeight >> maxParticleHeight would cause a false abort.
         if (!simParams->particle.useNoSlipParticles)
         {
            if ((simParams->computed.sinusOffset + simParams->computed.sinusAmplitude) > maxParticleHeight)
               WALBERLA_ABORT("Created particle bed is below desired sinus shape!");
            if (real_c(2) * simParams->computed.sinusAmplitude > (maxParticleHeight - simParams->computed.particleFixingHeight))
               WALBERLA_ABORT("Created particle bed is not high enough for desired sinus shape!");
         }
         if (simParams->particle.useNoSlipParticles && (simParams->computed.particleFixingHeight < maxParticleHeight))
            WALBERLA_ABORT("You are using no-slip BCs on particles (which does not set hydrodynamic forces) but do not "
                           "fix all particles - this leads to wrong behavior and is not permitted!")

         // fix particles that are below the particleFixingHeight ------------------------------------------------
         WALBERLA_LOG_INFO_ON_ROOT("Fixing particles below the particleFixingHeight: " << simParams->computed.particleFixingHeight);
         simObjs->particleStorage->forEachParticle(
            false, mesa_pd::kernel::SelectAll(), *simObjs->particleAccessor,
            [&simParams, &numFixedParticles](const size_t idx, auto& ac) {
               if (ac.getPosition(idx)[2] < simParams->computed.particleFixingHeight)
               {
                  mesa_pd::data::particle_flags::set(ac.getFlagsRef(idx), mesa_pd::data::particle_flags::FIXED);
                  ++numFixedParticles;
               }
            }, *simObjs->particleAccessor);
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT("Particle bed disabled (useParticleBed = false); running fluid-only simulation.");
      }
   }

   walberla::mpi::allReduceInplace(numFixedParticles, walberla::mpi::SUM);
   simObjs->numParticles = numParticles;

   WALBERLA_LOG_INFO_ON_ROOT("Total number of particles in the system: " << simObjs->numParticles << " (of which " << numFixedParticles << " are fixed)");


   // create top and bottom planes --------------------------------------------------------
   createPlane(*simObjs->particleStorage, bf->getDomain().minCorner(), Vector3< real_t >(0, 0,  1));
   createPlane(*simObjs->particleStorage, bf->getDomain().maxCorner(), Vector3< real_t >(0, 0, -1));

   // create Y-wall planes when Y is non-periodic
   if (!simParams->blockForest.periodicity[1])
   {
      createPlane(*simObjs->particleStorage, bf->getDomain().minCorner(), Vector3< real_t >(0,  1, 0));
      createPlane(*simObjs->particleStorage, bf->getDomain().maxCorner(), Vector3< real_t >(0, -1, 0));
   }

   // Pillars are GLOBAL/FIXED (createBox/createCylinder in Utility.h). BlockForestDataHandling skips GLOBAL
   // particles in serialize()/deserialize(), so pillars are not checkpointed and are recreated from the prm on
   // every run, including resumes. Their positions are logged for verification.
   for (const auto& pillar : simParams->pillars)
   {
      auto p = createBox(*simObjs->particleStorage, pillar.center, real_t(2) * pillar.halfExtents);
      p->setType(2);
      WALBERLA_LOG_INFO_ON_ROOT("Created box pillar: center=" << pillar.center
         << " halfExtents=" << pillar.halfExtents);
   }

   for (const auto& pillar : simParams->cylinderPillars)
   {
      createCylinder(*simObjs->particleStorage, pillar.center, pillar.radius);
      WALBERLA_LOG_INFO_ON_ROOT("Created cylinder pillar: center=" << pillar.center
         << " radius=" << pillar.radius);
   }

   // Floater creation (spheres/boxes/cylinders below) is skipped when resuming from a checkpoint: floaters are
   // ordinary particles restored from particleStorageFile, and creating them again would duplicate UIDs and abort.
   //
   // KNOWN LIMITATION: floaterUIDs/floaterBoxUIDs/floaterCylinderUIDs/dormantFloaterUIDs are only filled by the
   // creation loops, so on a resumed run releaseSettledFloatingSpheres() has nothing to act on. A floater still
   // dormant (FIXED) at checkpoint time is restored as such but never released. Fixing this needs a persisted
   // per-particle floater index (position cannot identify a particle, since particles move).
   if (bfParams.loadSnapshot)
   {
      // sized to the prm-declared counts but filled with the invalid UID 0, so indexed access stays valid
      simObjs->floaterUIDs.assign(simParams->floatingSpheres.size(), walberla::id_t(0));
      simObjs->floaterBoxUIDs.assign(simParams->floatingBoxes.size(), walberla::id_t(0));
      simObjs->floaterCylinderUIDs.assign(simParams->floatingCylinders.size(), walberla::id_t(0));
      WALBERLA_LOG_INFO_ON_ROOT("loadSnapshot=true: skipping floater (re-)creation -- "
         << simParams->floatingSpheres.size() << " sphere(s), " << simParams->floatingBoxes.size()
         << " box(es), " << simParams->floatingCylinders.size() << " cylinder(s) declared in this prm "
         << "are assumed already present (and, if still startInactive at checkpoint time, already "
         << "correctly FIXED) in the loaded particleStorage instead. floaterUIDs/floaterBoxUIDs/"
         << "floaterCylinderUIDs/dormantFloaterUIDs are left empty -- see this block's own comment for "
         << "why that only matters if some floaters were still dormant at checkpoint time.");
   }
   else
   {
   // Indexed loop so each sphere's UID is stored at its index fi. The rank whose subdomain contains pos creates
   // the particle and records the UID; allReduceInplace(SUM) then broadcasts it (only one rank owns each
   // position). ParticleKernels uses simObjs->floaterUIDs directly, since the local particle ordering does not
   // give this mapping when spheres are distributed over several ranks.
   simObjs->floaterUIDs.assign(simParams->floatingSpheres.size(), walberla::id_t(0));
   for (size_t fi = 0; fi < simParams->floatingSpheres.size(); ++fi)
   {
      const auto& fsp = simParams->floatingSpheres[fi];
      const mesa_pd::Vec3 pos(fsp.position[0], fsp.position[1], fsp.position[2]);
      if (simObjs->mesapdDomain->isContainedInProcessSubdomain(
            uint_c(mpi::MPIManager::instance()->rank()), pos))
      {
         auto p = simObjs->particleStorage->create(false);
         p->setPosition(pos);
         p->setBaseShape(std::make_shared< mesa_pd::data::Sphere >(fsp.radius));
         p->getBaseShapeRef()->updateMassAndInertia(fsp.densityRatio);
         p->setOwner(mpi::MPIManager::instance()->rank());
         p->setType(1);
         p->setInteractionRadius(fsp.radius);
         p->setLinearVelocity(mesa_pd::Vec3(fsp.initialVelocity[0], fsp.initialVelocity[1], fsp.initialVelocity[2]));
         // Hold the sphere FIXED during the settling phase so the FSLBM can establish a stable fill-level
         // distribution around the contact line first. pinPosition holds only the translational DOFs (applied
         // each DEM sub-cycle in ParticleKernels) and must not be combined with FIXED (same for boxes/cylinders).
         if (fsp.settlingSteps > 0 && !fsp.pinPosition)
            mesa_pd::data::particle_flags::set(p->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
         simObjs->floaterUIDs[fi] = p->getUid();
      }
      WALBERLA_LOG_INFO_ON_ROOT("Created floating sphere: radius=" << fsp.radius
         << " densityRatio=" << fsp.densityRatio << " position=" << pos
         << " settlingSteps=" << fsp.settlingSteps
         << " initialVelocity=" << fsp.initialVelocity
         << " pinPosition=" << (fsp.pinPosition ? "true" : "false"));
   }
   mpi::allReduceInplace(simObjs->floaterUIDs, mpi::SUM);
   // Register startInactive spheres as dormant after the UID sync, so every rank inserts the same resolved UID
   // (floatingSpheres is parsed identically on every rank, so no communication is needed).
   for (size_t fi = 0; fi < simParams->floatingSpheres.size(); ++fi)
      if (simParams->floatingSpheres[fi].startInactive)
         simObjs->dormantFloaterUIDs->insert(simObjs->floaterUIDs[fi]);

   // Floating boxes (type 3), same indexed-UID pattern as floating spheres. interactionRadius is the
   // circumscribed half-diagonal, so getParticleCellBB yields a conservative cell bounding box for the
   // shape-generic MEM cell test. FloatingBox rotates freely and uses Box::updateMassAndInertia.
   simObjs->floaterBoxUIDs.assign(simParams->floatingBoxes.size(), walberla::id_t(0));
   for (size_t bi = 0; bi < simParams->floatingBoxes.size(); ++bi)
   {
      const auto& fbp = simParams->floatingBoxes[bi];
      const mesa_pd::Vec3 pos(fbp.position[0], fbp.position[1], fbp.position[2]);
      const mesa_pd::Vec3 edge(real_t(2) * fbp.halfExtents[0],
                               real_t(2) * fbp.halfExtents[1],
                               real_t(2) * fbp.halfExtents[2]);
      const real_t halfDiag = std::sqrt(fbp.halfExtents[0] * fbp.halfExtents[0] +
                                        fbp.halfExtents[1] * fbp.halfExtents[1] +
                                        fbp.halfExtents[2] * fbp.halfExtents[2]);
      if (simObjs->mesapdDomain->isContainedInProcessSubdomain(
            uint_c(mpi::MPIManager::instance()->rank()), pos))
      {
         auto p = simObjs->particleStorage->create(false);
         p->setPosition(pos);
         p->setBaseShape(std::make_shared< FloatingBox >(edge));
         p->getBaseShapeRef()->updateMassAndInertia(fbp.densityRatio);
         p->setOwner(mpi::MPIManager::instance()->rank());
         p->setType(3);
         // + 2 * boxContactSkin bounds the circumscribed radius of the inflated box (see
         // ModelParameters::boxContactSkin); same over-approximation as circumR for cylinders below.
         p->setInteractionRadius(halfDiag + real_t(2) * simParams->model.boxContactSkin);
         if (std::abs(fbp.initialTiltAngleDeg) > real_t(1e-10))
         {
            const real_t angle = fbp.initialTiltAngleDeg * math::pi / real_t(180);
            // Rotate about Y-axis (body axis v in Bogner's diagram) by default; about
            // FloatingBoxParams::initialTiltAxis if given (randomised drop tests).
            mesa_pd::Vec3 tiltAxis(fbp.initialTiltAxis[0], fbp.initialTiltAxis[1], fbp.initialTiltAxis[2]);
            tiltAxis = (tiltAxis.length() > real_t(1e-12)) ? tiltAxis / tiltAxis.length() : mesa_pd::Vec3(0, 1, 0);
            p->setRotation(mesa_pd::Rot3(tiltAxis, angle));
         }
         // pinPosition holds only the translational DOFs (applied each DEM sub-cycle in ParticleKernels); FIXED
         // would also freeze rotation. settlingSteps uses FIXED (full freeze) and excludes pinPosition.
         if (fbp.settlingSteps > 0 && !fbp.pinPosition)
            mesa_pd::data::particle_flags::set(p->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
         simObjs->floaterBoxUIDs[bi] = p->getUid();
      }
      WALBERLA_LOG_INFO_ON_ROOT("Created floating box: halfExtents=" << fbp.halfExtents
         << " densityRatio=" << fbp.densityRatio << " position=" << pos
         << " settlingSteps=" << fbp.settlingSteps
         << " initialTiltAngleDeg=" << fbp.initialTiltAngleDeg
         << " pinPosition=" << (fbp.pinPosition ? "true" : "false"));
   }
   mpi::allReduceInplace(simObjs->floaterBoxUIDs, mpi::SUM);
   for (size_t bi = 0; bi < simParams->floatingBoxes.size(); ++bi)
      if (simParams->floatingBoxes[bi].startInactive)
         simObjs->dormantFloaterUIDs->insert(simObjs->floaterBoxUIDs[bi]);

   // Floating cylinders / logs (type 4), same indexed-UID pattern as spheres and boxes. The mesa_pd Cylinder's
   // axis is its body z-axis; the initial rotation aligns body-z with the world `axis` (default <0,1,0>, a
   // spanwise log) and then applies the optional tilt about the world y-axis.
   // interactionRadius = circumscribed-sphere radius sqrt(r² + (L/2)²).
   simObjs->floaterCylinderUIDs.assign(simParams->floatingCylinders.size(), walberla::id_t(0));
   for (size_t ci = 0; ci < simParams->floatingCylinders.size(); ++ci)
   {
      const auto& fcp = simParams->floatingCylinders[ci];
      const mesa_pd::Vec3 pos(fcp.position[0], fcp.position[1], fcp.position[2]);
      const real_t halfLength = real_t(0.5) * fcp.length;
      // + 2 * cylinderContactSkin: bounds the circumscribing radius of the inflated rod (see ModelParameters)
      const real_t circumR    = std::sqrt(fcp.radius * fcp.radius + halfLength * halfLength) + real_t(2) * simParams->model.cylinderContactSkin;

      // Align body-z with fcp.axis: rotate about (ez × axis) by acos(ez·axis); antiparallel → flip
      // about x. Then compose the world-Y tilt on top.
      mesa_pd::Rot3 rot;
      const mesa_pd::Vec3 ez(real_t(0), real_t(0), real_t(1));
      const real_t cosA = ez * fcp.axis;
      if (cosA < real_t(1) - real_t(1e-10))
      {
         const mesa_pd::Vec3 rotAxis = (cosA > real_t(-1) + real_t(1e-10))
                                          ? (ez % fcp.axis).getNormalized()
                                          : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         rot.rotate(rotAxis, std::acos(std::clamp(cosA, real_t(-1), real_t(1))));
      }
      if (std::abs(fcp.initialTiltAngleDeg) > real_t(1e-10))
         rot.rotate(mesa_pd::Vec3(real_t(0), real_t(1), real_t(0)), fcp.initialTiltAngleDeg * math::pi / real_t(180));

      if (simObjs->mesapdDomain->isContainedInProcessSubdomain(
            uint_c(mpi::MPIManager::instance()->rank()), pos))
      {
         auto p = simObjs->particleStorage->create(false);
         p->setPosition(pos);
         // ModelParameters.useCapsuleShapeForCylinders (experimental, default false): a Capsule (cylinder with
         // hemispherical end caps) has the same (radius, length) constructor and body-z axis convention as
         // Cylinder, so nothing else in this loop (position, rotation, interactionRadius, FIXED flag,
         // settlingSteps) needs to change.
         if (simParams->model.useCapsuleShapeForCylinders)
            p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(fcp.radius, fcp.length));
         else
            p->setBaseShape(std::make_shared< mesa_pd::data::Cylinder >(fcp.radius, fcp.length));
         p->getBaseShapeRef()->updateMassAndInertia(fcp.densityRatio);
         p->setOwner(mpi::MPIManager::instance()->rank());
         p->setType(4);
         p->setInteractionRadius(circumR);
         p->setRotation(rot);
         // same settling semantics as boxes: FIXED (full freeze) unless pinPosition holds the translational DOFs
         if (fcp.settlingSteps > 0 && !fcp.pinPosition)
            mesa_pd::data::particle_flags::set(p->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
         simObjs->floaterCylinderUIDs[ci] = p->getUid();
      }
      WALBERLA_LOG_INFO_ON_ROOT("Created floating cylinder: radius=" << fcp.radius
         << " length=" << fcp.length << " axis=" << fcp.axis
         << " densityRatio=" << fcp.densityRatio << " position=" << pos
         << " settlingSteps=" << fcp.settlingSteps
         << " initialTiltAngleDeg=" << fcp.initialTiltAngleDeg
         << " pinPosition=" << (fcp.pinPosition ? "true" : "false"));
   }
   mpi::allReduceInplace(simObjs->floaterCylinderUIDs, mpi::SUM);
   for (size_t ci = 0; ci < simParams->floatingCylinders.size(); ++ci)
      if (simParams->floatingCylinders[ci].startInactive)
         simObjs->dormantFloaterUIDs->insert(simObjs->floaterCylinderUIDs[ci]);
   } // end of !bfParams.loadSnapshot

   // choose particle sync algorithm -----------------------------------------------------
   // blockSyncExtension must be >= the largest floater extent so that ghost copies reach every block a floater
   // overlaps: sphere radius, box half-diagonal, cylinder circumscribed radius (including the contact skin).
   real_t blockSyncExtension = real_c(2.5);
   for (const auto& fsp : simParams->floatingSpheres)
      blockSyncExtension = std::max(blockSyncExtension, fsp.radius);
   for (const auto& fbp : simParams->floatingBoxes)
   {
      const real_t halfDiag = std::sqrt(fbp.halfExtents[0] * fbp.halfExtents[0] +
                                        fbp.halfExtents[1] * fbp.halfExtents[1] +
                                        fbp.halfExtents[2] * fbp.halfExtents[2]);
      blockSyncExtension = std::max(blockSyncExtension, halfDiag);
   }
   for (const auto& fcp : simParams->floatingCylinders)
   {
      const real_t circumR = std::sqrt(fcp.radius * fcp.radius +
                                       real_t(0.25) * fcp.length * fcp.length) + real_t(2) * simParams->model.cylinderContactSkin;
      blockSyncExtension = std::max(blockSyncExtension, circumR);
   }
   WALBERLA_LOG_INFO_ON_ROOT("blockSyncExtension = " << blockSyncExtension);
   // maxParticleDiameter selects the sync algorithm (SyncNextNeighbors vs SyncGhostOwners); 2x the
   // circumscribed radius of each floater type, so particles spanning two or more blocks get SyncGhostOwners.
   real_t maxParticleDiameter = simParams->computed.avgParticleDiameter * real_c(1.1);
   for (const auto& fsp : simParams->floatingSpheres)
      maxParticleDiameter = std::max(maxParticleDiameter, fsp.radius * real_c(2));
   for (const auto& fbp : simParams->floatingBoxes)
   {
      const real_t halfDiag = std::sqrt(fbp.halfExtents[0] * fbp.halfExtents[0] +
                                        fbp.halfExtents[1] * fbp.halfExtents[1] +
                                        fbp.halfExtents[2] * fbp.halfExtents[2]);
      maxParticleDiameter = std::max(maxParticleDiameter, halfDiag * real_c(2));
   }
   for (const auto& fcp : simParams->floatingCylinders)
   {
      const real_t circumR = std::sqrt(fcp.radius * fcp.radius +
                                       real_t(0.25) * fcp.length * fcp.length) + real_t(2) * simParams->model.cylinderContactSkin;
      maxParticleDiameter = std::max(maxParticleDiameter, circumR * real_c(2));
   }

   // SyncNextNeighbors only exchanges with neighbouring ranks. Each floater is ghosted to every block within its
   // circumscribed radius + blockSyncExtension; if that reach exceeds the block width in a direction with at least
   // three blocks, a floater can hold a ghost two blocks away, and after it migrates its new owner sends to a
   // non-neighbour rank ("Received invalid notification type"). The diameter test below does not catch this,
   // so such cases use SyncGhostOwners.
   const real_t particleSyncReach = real_c(0.5) * maxParticleDiameter + blockSyncExtension;
   bool         ghostsStayWithinNeighbours = true;
   for (uint_t d = 0; d < 3; ++d)
   {
      const uint_t blocksAlongD = simParams->domain.domainSize[d] / simParams->blockForest.cellsPerBlock[d];
      if (blocksAlongD >= uint_t(3) && particleSyncReach > real_c(simParams->blockForest.cellsPerBlock[d]))
         ghostsStayWithinNeighbours = false;
   }
   WALBERLA_LOG_INFO_ON_ROOT("Particle sync reach (circumscribed radius + blockSyncExtension) = " << particleSyncReach
                             << (ghostsStayWithinNeighbours ? "" : " exceeds the block width along a direction with "
                                                                   ">= 3 blocks"));
   if (ghostsStayWithinNeighbours &&
       maxParticleDiameter < real_c(2) * real_c(simParams->blockForest.cellsPerBlock.min()) - blockSyncExtension) {
      WALBERLA_LOG_INFO_ON_ROOT("Using SyncNextNeighbors for particles");
      simObjs->syncCall = [&simObjs, blockSyncExtension]() {
         mesa_pd::mpi::SyncNextNeighbors{}(*simObjs->particleStorage, *simObjs->mesapdDomain, blockSyncExtension);
      };
   }
   else {
      WALBERLA_LOG_INFO_ON_ROOT("Using SyncGhostOwners for particles")
      simObjs->syncCall = [&simObjs, blockSyncExtension]() {
         mesa_pd::mpi::SyncGhostOwners{}(*simObjs->particleStorage, *simObjs->mesapdDomain, blockSyncExtension);
      };
   }

   WALBERLA_LOG_INFO_ON_ROOT("Synchronizing particles");
   simObjs->syncCall();

   // Repeat the sync: SyncGhostOwners::checkAndResolveOverlap() only sends new ghost copies to
   // domain.getNeighborProcesses(), i.e. ranks one hop from where the particle is already known, but a large
   // floater can overlap ranks several hops away. In a running simulation the ghost set grows over successive
   // steps; at a fresh start or a checkpoint resume a single sync can leave ranks two or more hops away without
   // a ghost, which triggers the "[STALE-PARTICLE-IDX]" abort in SimpleBB.h on the first step. Each round seeds
   // the next; 4 rounds cover the 2 hops needed in practice (one-time setup cost).
   for (uint_t round = uint_t(0); round < uint_t(4); ++round) simObjs->syncCall();

   // create particle field (mapping of particle UIDs to lbm cells)
   simObjs->particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>(
      bf, "Particle field", simObjs->particleAccessor->getInvalidUid(), field::fzyx, uint_c(2));
}

void setupBoundaryHandling(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("setupBoundaryHandling");

   // Live MEM reference height (ModelParameters.useDynamicReferenceHeight): seeded with the still-water level;
   // with the toggle on, ParticleKernels' computeSubmergedVolumes() EMA-updates it every step from the floaters'
   // fitted local waterlines. It only enters the force (SimpleBB's rhoHydStat), never the fluid. With the
   // toggle off it is never written and the static-height functor below is used.
   simObjs->memReferenceHeight = std::make_shared< real_t >(simParams->staticReferenceHeight());

   WALBERLA_LOG_INFO_ON_ROOT("[MEM reference height] "
      << (simParams->model.useBognerThreePhaseMEM
             ? "DISABLED -- useBognerThreePhaseMEM=true: no rhoHydStat subtraction at all, reproducing "
               "Bogner & Rude (2013) Eq. 5.9's original three-phase MEM force as published"
             : (simParams->model.useDynamicReferenceHeight
                   ? "DYNAMIC (live, EMA-tracked from floaters' fitted local waterlines)"
                   : "STATIC (flat H=" + std::to_string(simParams->staticReferenceHeight())
                        + (simParams->model.staticReferenceHeightOverride >= real_t(0)
                              ? ", OVERRIDDEN via staticReferenceHeightOverride — flat absoluteLiquidHeight would be "
                                   + std::to_string(simParams->computed.absoluteLiquidHeight)
                              : " = absoluteLiquidHeight") + ")")));

   // useBognerThreePhaseMEM: a null functor selects SimpleBB's `hydrostaticDensityFct_ == nullptr` branch, i.e.
   // no reference-height subtraction. Together with the (default) fill-weighting this is Bogner & Rude's
   // Eq. 5.9: a phi(x_b)-weighted blend of the boundary-link PDFs and the gas equilibrium distribution.
   std::function< real_t(const Vector3< real_t >&) > hydrostaticFtor;
   if (!simParams->model.useBognerThreePhaseMEM)
   {
      hydrostaticFtor = simParams->model.useDynamicReferenceHeight
         ? HydrostaticDensityFtor(simObjs->acceleration, simObjs->memReferenceHeight)
         : HydrostaticDensityFtor(simObjs->acceleration, simParams->staticReferenceHeight());
   }

   simObjs->antidunesBoundaryHandling = std::make_shared< AntidunesBoundaryHandling_T >(
      simObjs->blockForest, simObjs->pdfFieldID, simObjs->fillFieldID, simObjs->particleFieldID,
      simObjs->particleAccessor,
      hydrostaticFtor,
      simParams->model.useFillWeighting,
      simParams->model.useInterpolatedBounceBack);

   // initialize BCs from config (NoSlip top/bottom, CellInterval pillar defined in Antidunes.prm)
   simObjs->antidunesBoundaryHandling->initFromConfig(simParams->cfgPtr->getOneBlock("BoundaryParameters"));

   // NoSlip walls on non-periodic lateral faces (x: W/E, y: N/S), so the walls always match
   // BlockForestParameters.periodicity. A non-periodic face without a boundary condition has undefined PDFs and
   // produces NaN. Periodic faces stay wrapped. Top/bottom (z, never periodic) come from the prm; setting an
   // existing wall again is idempotent.
   {
      auto&       bh = *simObjs->antidunesBoundaryHandling;
      const auto& bf = *simObjs->blockForest;
      const cell_idx_t wd = cell_idx_c(-1);
      if (!bf.isXPeriodic()) { bh.setNoSlipAtBorder(stencil::W, wd); bh.setNoSlipAtBorder(stencil::E, wd); }
      if (!bf.isYPeriodic()) { bh.setNoSlipAtBorder(stencil::N, wd); bh.setNoSlipAtBorder(stencil::S, wd); }
   }

   // store the flagField and flagInfo handles in simObjs (only for convenience)
   simObjs->flagFieldID = simObjs->antidunesBoundaryHandling->getFlagFieldID();
   simObjs->flagInfo    = simObjs->antidunesBoundaryHandling->getFlagInfo();
}

void mapParticlesOntoCells(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("mapParticlesOntoCells");

   // Mark the remaining cells as domain first: the mapping kernels below only claim cells for which isDomain()
   // is true, and at this point only initFromConfig() has run (initFlagsFromFillLevel() follows in
   // Antidunes.cpp). Without this the setup mapping maps no floater cells, which is fatal on a checkpoint resume
   // with active floaters (other users of these kernels call fillWithDomain() for the same reason).
   // fillWithDomain() only touches empty cells and keeps the boundary flags already set.
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      auto* handling = blockIt->getData< AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         simObjs->antidunesBoundaryHandling->getHandlingID());
      handling->fillWithDomain(handling->getFlagField()->nrOfGhostLayers());
   }

   // Select floaters: type 1 = floating sphere, type 3 = floating box, type 4 = floating cylinder. Fixed
   // obstacles (pillars, walls) are mapped separately as static NoSlip/FreeSlip. Dormant floaters (startInactive,
   // before settlingSteps) are excluded so they do not claim fluid cells; this lets several floaters share one
   // spawn position. Must stay consistent with movingObstacleSelector in BuildTimeloop.h (per-step mapping).
   auto dormantUIDs = simObjs->dormantFloaterUIDs;
   auto floaterSelector = [dormantUIDs](const size_t particleIdx, const auto& ac)
   {
      const uint_t t = ac.getType(particleIdx);
      if (t != uint_t(1) && t != uint_t(3) && t != uint_t(4)) return false;
      return dormantUIDs->count(ac.getUid(particleIdx)) == 0;
   };

   // map particles onto the LBM domain (initialisation only; the per-step mapping runs in timeloopAfterDEM)
   if (simParams->particle.useNoSlipParticles) {
      // sets noSlip flag (in flagField through boundaryhandling) for all cells inside the particle
      lbm_mesapd_coupling::ParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
         kernel(simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID());
      simObjs->particleStorage->forEachParticle(false, floaterSelector, *simObjs->particleAccessor,
                                             kernel, *simObjs->particleAccessor,
                                             AntidunesBoundaryHandling_T::noSlipFlagID);
   }
   else {
      // same as above, but also stores the particle UID in particleFieldID (no hydrodynamic forces computed here)
      lbm_mesapd_coupling::MovingParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
         kernel(simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID(), simObjs->particleFieldID,
                simParams->model.mappingTieOffset);
      simObjs->particleStorage->forEachParticle(false, floaterSelector, *simObjs->particleAccessor,
                                             kernel, *simObjs->particleAccessor,
                                             AntidunesBoundaryHandling_T::movingObstacleFlagID);
   }
}

// Called every timestep: releases each floating sphere, box and cylinder whose settlingSteps is reached (clears
// FIXED and removes it from the dormant set). While settling, the floater is held FIXED so the FSLBM can establish
// a stable fill-level / normal field around the contact line before the particle moves.
//
// Each released UID is also inserted into simObjs->freshlyActivatedFloaterUIDs (see SimulationObjects.h).
// ComputeContacts (ParticleKernels.h) aborts with diagnostics if, in the next doDEMSubCycles() call, it finds a
// contact involving such a UID, i.e. a floater released on top of an active neighbour. This uses the real
// contact detection instead of a separate geometric check. The set is cleared at the end of that call.
void releaseSettledFloatingSpheres(const std::shared_ptr< const SimulationParams >& simParams,
                                   const std::shared_ptr< SimulationObjects >&       simObjs,
                                   const uint_t currentStep)
{
   // Release each floater individually by UID (floaterUIDs[i] <-> floatingSpheres[i]), so staggered release
   // schedules work.
   for (size_t i = 0; i < simParams->floatingSpheres.size(); ++i)
   {
      const auto& fsp = simParams->floatingSpheres[i];
      if (fsp.settlingSteps == uint_t(0) || currentStep != fsp.settlingSteps) continue;
      const walberla::id_t uid = simObjs->floaterUIDs[i];

      simObjs->freshlyActivatedFloaterUIDs->insert(uid);

      simObjs->particleStorage->forEachParticle(
         false, mesa_pd::kernel::SelectAll(), *simObjs->particleAccessor,
         [uid](const size_t idx, auto& ac) {
            if (ac.getUid(idx) == uid)
               mesa_pd::data::particle_flags::unset(ac.getFlagsRef(idx),
                                                    mesa_pd::data::particle_flags::FIXED);
         },
         *simObjs->particleAccessor);
      simObjs->dormantFloaterUIDs->erase(uid);

      WALBERLA_LOG_INFO_ON_ROOT("Released floating sphere uid=" << uid << " (index " << i
                                << ") from settling phase at step " << currentStep);
   }

   for (size_t i = 0; i < simParams->floatingBoxes.size(); ++i)
   {
      const auto& fbp = simParams->floatingBoxes[i];
      if (fbp.settlingSteps == uint_t(0) || currentStep != fbp.settlingSteps) continue;
      const walberla::id_t uid = simObjs->floaterBoxUIDs[i];

      simObjs->freshlyActivatedFloaterUIDs->insert(uid);

      simObjs->particleStorage->forEachParticle(
         false, mesa_pd::kernel::SelectAll(), *simObjs->particleAccessor,
         [uid](const size_t idx, auto& ac) {
            if (ac.getUid(idx) == uid)
               mesa_pd::data::particle_flags::unset(ac.getFlagsRef(idx),
                                                    mesa_pd::data::particle_flags::FIXED);
         },
         *simObjs->particleAccessor);
      simObjs->dormantFloaterUIDs->erase(uid);

      WALBERLA_LOG_INFO_ON_ROOT("Released floating box uid=" << uid << " (index " << i
                                << ") from settling phase at step " << currentStep);
   }

   for (size_t i = 0; i < simParams->floatingCylinders.size(); ++i)
   {
      const auto& fcp = simParams->floatingCylinders[i];
      if (fcp.settlingSteps == uint_t(0) || currentStep != fcp.settlingSteps) continue;
      const walberla::id_t uid = simObjs->floaterCylinderUIDs[i];

      simObjs->freshlyActivatedFloaterUIDs->insert(uid);

      simObjs->particleStorage->forEachParticle(
         false, mesa_pd::kernel::SelectAll(), *simObjs->particleAccessor,
         [uid](const size_t idx, auto& ac) {
            if (ac.getUid(idx) == uid)
               mesa_pd::data::particle_flags::unset(ac.getFlagsRef(idx),
                                                    mesa_pd::data::particle_flags::FIXED);
         },
         *simObjs->particleAccessor);
      simObjs->dormantFloaterUIDs->erase(uid);

      WALBERLA_LOG_INFO_ON_ROOT("Released floating cylinder uid=" << uid << " (index " << i
                                << ") from settling phase at step " << currentStep);
   }
}

// Monitor a fixed list of world-coordinate cells every timestep.
// Logs their PDF values (all stencil directions), density, fill level, and flag so that
// we can watch the evolution of specific cells before they turn NaN.
// watchCoords: list of world-space cell *centres* to watch, e.g. {8.5, 29.5, 60.5}.
// Call once per timestep (or per sub-sweep) with a label string identifying the sweep.
void trackWatchCells(const std::shared_ptr< SimulationObjects >& simObjs,
                     const char* const label, const uint_t t,
                     const std::vector< Vector3< real_t > >& watchCoords)
{
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto& aabb      = blockIt->getAABB();
      const auto* pdfField  = blockIt->getData< const PdfField_T  >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      // Force-density field: read directly by the generated cumulant collision kernel, separately from the
      // PDFs, so PDF-based NaN/precursor checks cannot see it.
      const auto* forceField = blockIt->getData< const VectorFieldFlattened_T >(simObjs->forceDensityFieldID);
      // Owning UID from particleField: shows whether a cell's obstacle ownership flips between nearby floaters
      // across steps, which the flag classification alone does not reveal.
      const auto* particleFieldPtr =
         blockIt->getData< const lbm_mesapd_coupling::ParticleField_T >(simObjs->particleFieldID);
      // Curvature/normal fields: the finite-difference curvature method used here
      // (CurvatureSweepFiniteDifferences, the least stable of the three methods Bogner compared) can enter a
      // self-sustaining oscillation when spurious currents perturb the interface layer (Bogner 2017,
      // Section 4.4.3.1, p.99). Logging them shows whether they are already degenerate before a density collapse.
      const auto* curvatureField = blockIt->getData< const ScalarField_T >(simObjs->curvatureFieldID);
      const auto* normalField    = blockIt->getData< const VectorField_T >(simObjs->normalFieldID);
      const auto* obstacleNormalField =
         blockIt->getData< const VectorField_T >(simObjs->obstacleNormalFieldID);

      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      // classification as in checkPdfNaN's classify(), distinguishing gas cells from obstacle cells
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const FlagUID formerMOuid("former moving obstacle");
      const flag_t formerMOFlag = flagField->flagExists(formerMOuid) ? flagField->getFlag(formerMOuid) : flag_t(0);

      for (const auto& wc : watchCoords)
      {
         // cell centre → block-local index:  ix = floor(wc.x - aabb.xMin())
         // only process if this block owns the cell
         if (!aabb.contains(wc)) continue;

         const cell_idx_t ix = cell_idx_c(std::floor(wc[0] - aabb.xMin()));
         const cell_idx_t iy = cell_idx_c(std::floor(wc[1] - aabb.yMin()));
         const cell_idx_t iz = cell_idx_c(std::floor(wc[2] - aabb.zMin()));

         // guard against ghost layer confusion
         if (ix < 0 || iy < 0 || iz < 0) continue;
         if (ix >= cell_idx_c(pdfField->xSize()) ||
             iy >= cell_idx_c(pdfField->ySize()) ||
             iz >= cell_idx_c(pdfField->zSize())) continue;

         // compute density from PDFs
         real_t rho = real_t(0);
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
            rho += pdfField->get(ix, iy, iz, f);

         const real_t fill = fillField->get(ix, iy, iz);
         const flag_t fl = flagField->get(ix, iy, iz);
         const char* flagStr =
            (fl & obstacleMask) != flag_t(0)                        ? "OBST" :
            isFlagSet(fl, interfaceFlag)                            ? "INTF" :
            isFlagSet(fl, liquidFlag)                                ? "LIQ " :
            isFlagSet(fl, gasFlag)                                   ? "GAS " :
            (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag)) ? "fMO " : "??? ";

         // build PDF list string
         std::ostringstream pdfStr;
         pdfStr << std::scientific << std::setprecision(4);
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
         {
            pdfStr << "f[" << f << "]=" << pdfField->get(ix, iy, iz, f);
            if (f + 1 < PdfField_T::Stencil::Size) pdfStr << " ";
         }

         const real_t fx = forceField->get(ix, iy, iz, uint_t(0));
         const real_t fy = forceField->get(ix, iy, iz, uint_t(1));
         const real_t fz = forceField->get(ix, iy, iz, uint_t(2));
         const walberla::id_t owningUid = particleFieldPtr->get(ix, iy, iz);
         const real_t curvature = curvatureField->get(ix, iy, iz);
         const Vector3< real_t > normal         = normalField->get(ix, iy, iz);
         const Vector3< real_t > obstacleNormal = obstacleNormalField->get(ix, iy, iz);

         WALBERLA_LOG_INFO("[WATCH t=" << t << " " << label
            << "] world=(" << wc[0] << "," << wc[1] << "," << wc[2] << ")"
            << " local=(" << ix << "," << iy << "," << iz << ")"
            << " flag=" << flagStr
            << " owningUid=" << owningUid
            << " rho=" << std::scientific << std::setprecision(6) << rho
            << " fill=" << fill
            << " curvature=" << curvature
            << " normal=" << normal
            << " obstacleNormal=" << obstacleNormal
            << " force=(" << fx << "," << fy << "," << fz << ")"
            << " | " << pdfStr.str());
      }
   }
}

// ModelParameters.logFluidBoxHistory (+ Min/Max/StartStep/EndStep, SimulationParams.h): logs every obstacle or
// interface cell inside a fixed world-space box once per main step to a per-rank CSV, to follow the spatial
// evolution of a region (e.g. around a crash). Once per main step suffices since the fluid field only changes
// once per main step; restricting to OBST/INTF cells keeps the row count bounded. Each rank opens and truncates
// its own CSV on the first call (function-local static), as FloaterHistoryLogger does (ParticleKernels.h).
// Called from addAfterTimeStepDiagnosticSweep (BuildTimeloop.h).
void logFluidBoxHistory(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t,
                        const std::string& baseFolder,
                        const Vector3< real_t >& boxMin, const Vector3< real_t >& boxMax)
{
   static std::ofstream file;
   static bool opened = false;
   if (!opened)
   {
      const int rank = mpi::MPIManager::instance()->rank();
      file.open(baseFolder + "/rank_" + std::to_string(rank) + ".csv", std::ios::out | std::ios::trunc);
      WALBERLA_CHECK(file.is_open(), "[FluidBoxHistory] failed to open CSV for rank " << rank
         << " in \"" << baseFolder << "\" -- does the directory exist?");
      file << "mainStep,rank,worldX,worldY,worldZ,flag,owningUid,density,fill,maxAbsPdf\n";
      opened = true;
   }

   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto& aabb = blockIt->getAABB();
      // Skip blocks that don't overlap the box at all -- avoids paying the field-lookup cost below
      // for every block on every rank when the box is small relative to the domain.
      if (aabb.xMax() < boxMin[0] || aabb.xMin() > boxMax[0] ||
          aabb.yMax() < boxMin[1] || aabb.yMin() > boxMax[1] ||
          aabb.zMax() < boxMin[2] || aabb.zMin() > boxMax[2]) continue;

      const auto* pdfField  = blockIt->getData< const PdfField_T  >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto* particleFieldPtr =
         blockIt->getData< const lbm_mesapd_coupling::ParticleField_T >(simObjs->particleFieldID);

      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const FlagUID formerMOuid("former moving obstacle");
      const flag_t formerMOFlag = flagField->flagExists(formerMOuid) ? flagField->getFlag(formerMOuid) : flag_t(0);

      for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()); ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(pdfField->ySize()); ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(pdfField->xSize()); ++ix)
      {
         const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
         const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
         const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
         if (wx < boxMin[0] || wx > boxMax[0] || wy < boxMin[1] || wy > boxMax[1] ||
             wz < boxMin[2] || wz > boxMax[2]) continue;

         const flag_t fl = flagField->get(ix, iy, iz);
         const bool isObstacle  = (fl & obstacleMask) != flag_t(0);
         const bool isInterface = isFlagSet(fl, interfaceFlag);
         if (!isObstacle && !isInterface) continue; // keep row count bounded -- see comment above

         const char* flagStr =
            isObstacle  ? "OBST" :
            isInterface ? "INTF" :
            isFlagSet(fl, liquidFlag)                                     ? "LIQ" :
            isFlagSet(fl, gasFlag)                                        ? "GAS" :
            (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag))    ? "fMO" : "???";

         real_t rho = real_t(0);
         real_t maxAbsPdf = real_t(0);
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
         {
            const real_t v = pdfField->get(ix, iy, iz, f);
            rho += v;
            maxAbsPdf = std::max(maxAbsPdf, std::abs(v));
         }
         const real_t fill = fillField->get(ix, iy, iz);
         const walberla::id_t owningUid = particleFieldPtr->get(ix, iy, iz);

         file << t << ',' << mpi::MPIManager::instance()->rank() << ','
              << wx << ',' << wy << ',' << wz << ','
              << flagStr << ',' << owningUid << ','
              << rho << ',' << fill << ',' << maxAbsPdf << '\n';
      }
   }
   // flush on every call (once per main step) so the data survives an abort, as in FloaterHistoryLogger
   file.flush();
}

// ModelParameters.logDensityColumn (+ X/Y/StartStep/EndStep, SimulationParams.h): logs every cell (all flags) of
// one world (x, y) column over the full domain height once per main step to a per-rank CSV. This gives the
// far-field density profile rho(z), to check whether the fluid is hydrostatically stratified under gravity
// (weakly compressible LBM, p = cs^2*rho) rather than at the uniform rho=1 assumed by a simple Archimedes
// estimate, and a fill-level-based free-surface height. Used for the still-water equilibrium of floaters.
void logDensityColumn(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t,
                      const std::string& baseFolder, const real_t columnX, const real_t columnY)
{
   static std::ofstream file;
   static bool opened = false;
   if (!opened)
   {
      const int rank = mpi::MPIManager::instance()->rank();
      file.open(baseFolder + "/rank_" + std::to_string(rank) + ".csv", std::ios::out | std::ios::trunc);
      WALBERLA_CHECK(file.is_open(), "[DensityColumn] failed to open CSV for rank " << rank
         << " in \"" << baseFolder << "\" -- does the directory exist?");
      file << "mainStep,rank,worldX,worldY,worldZ,flag,density,fill\n";
      opened = true;
   }

   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto& aabb = blockIt->getAABB();
      // Skip blocks that don't overlap the (x, y) column at all.
      if (columnX < aabb.xMin() || columnX > aabb.xMax() ||
          columnY < aabb.yMin() || columnY > aabb.yMax()) continue;

      const auto* pdfField  = blockIt->getData< const PdfField_T  >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);

      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const FlagUID formerMOuid("former moving obstacle");
      const flag_t formerMOFlag = flagField->flagExists(formerMOuid) ? flagField->getFlag(formerMOuid) : flag_t(0);

      const cell_idx_t ix = cell_idx_c(std::floor(columnX - aabb.xMin()));
      const cell_idx_t iy = cell_idx_c(std::floor(columnY - aabb.yMin()));
      if (ix < 0 || iy < 0 || ix >= cell_idx_c(pdfField->xSize()) || iy >= cell_idx_c(pdfField->ySize())) continue;

      for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()); ++iz)
      {
         const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);

         const flag_t fl = flagField->get(ix, iy, iz);
         const bool isObstacle  = (fl & obstacleMask) != flag_t(0);
         const bool isInterface = isFlagSet(fl, interfaceFlag);
         const char* flagStr =
            isObstacle  ? "OBST" :
            isInterface ? "INTF" :
            isFlagSet(fl, liquidFlag)                                  ? "LIQ" :
            isFlagSet(fl, gasFlag)                                     ? "GAS" :
            (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag)) ? "fMO" : "???";

         real_t rho = real_t(0);
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
            rho += pdfField->get(ix, iy, iz, f);
         const real_t fill = fillField->get(ix, iy, iz);

         file << t << ',' << mpi::MPIManager::instance()->rank() << ','
              << columnX << ',' << columnY << ',' << wz << ','
              << flagStr << ',' << rho << ',' << fill << '\n';
      }
   }
   file.flush();
}

// ModelParameters.logWaterHeight: fill-integrated free-surface height on a regular X x Y grid of
// columns (see SimulationParams.h). Every rank appends, for each grid column its blocks intersect, the partial sum of
// fill over its liquid+interface cells (obstacle/gas/former-MO cells contribute 0) plus the number of obstacle cells in
// its part of the column. The analysis sums partialFillSum over ranks per (step, gx, gy).
void logWaterHeights(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t, const std::string& baseFolder,
                     const uint_t gridX, const uint_t gridY)
{
   static std::ofstream file;
   static bool opened = false;
   if (!opened)
   {
      const int rank = mpi::MPIManager::instance()->rank();
      file.open(baseFolder + "/rank_" + std::to_string(rank) + ".csv", std::ios::out | std::ios::trunc);
      WALBERLA_CHECK(file.is_open(), "[WaterHeight] failed to open CSV for rank " << rank << " in \"" << baseFolder << "\"");
      file << "mainStep,rank,gx,gy,partialFillSum,nObstacleCells\n";
      opened = true;
   }
   const auto domain = simObjs->blockForest->getDomain();
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto& aabb = blockIt->getAABB();
      const auto* pdfField  = blockIt->getData< const PdfField_T >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      for (uint_t jy = 0; jy < gridY; ++jy)
      for (uint_t jx = 0; jx < gridX; ++jx)
      {
         const real_t wx = (real_c(jx) + real_t(0.5)) * domain.xSize() / real_c(gridX);
         const real_t wy = (real_c(jy) + real_t(0.5)) * domain.ySize() / real_c(gridY);
         if (wx < aabb.xMin() || wx >= aabb.xMax() || wy < aabb.yMin() || wy >= aabb.yMax()) continue;
         const cell_idx_t ix = cell_idx_c(std::floor(wx - aabb.xMin()));
         const cell_idx_t iy = cell_idx_c(std::floor(wy - aabb.yMin()));
         if (ix < 0 || iy < 0 || ix >= cell_idx_c(pdfField->xSize()) || iy >= cell_idx_c(pdfField->ySize())) continue;
         real_t sum = real_t(0);
         uint_t nObst = 0;
         for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()); ++iz)
         {
            const flag_t fl = flagField->get(ix, iy, iz);
            if ((fl & obstacleMask) != flag_t(0)) { ++nObst; continue; }
            if (isFlagSet(fl, liquidFlag) || isFlagSet(fl, interfaceFlag)) sum += fillField->get(ix, iy, iz);
         }
         file << t << ',' << mpi::MPIManager::instance()->rank() << ',' << cell_idx_t(std::floor(wx)) << ',' << cell_idx_t(std::floor(wy))
              << ',' << sum << ',' << nObst << '\n';
      }
   }
   file.flush();
}

// Scan every liquid and interface cell on every local block for NaN/Inf PDFs.
// Logs the world coordinates, block-local coordinates, and PDF index of the first NaN found
// on each rank, then returns true if any NaN was detected globally.
// Call this after each major sweep to pinpoint which one first corrupts the PDF field.
bool checkPdfNaN(const std::shared_ptr< SimulationObjects >& simObjs,
                 const char* const location, const uint_t t,
                 real_t lowRhoPrecursorThreshold = real_t(0),
                 real_t highVelocityPrecursorThreshold = real_t(0),
                 real_t populationPrecursorThreshold = real_t(0),
                 real_t densityPrecursorLow = real_t(0), real_t densityPrecursorHigh = real_t(0), bool densityPrecursorLiquidOnly = false)
{
   bool localNaN = false;
   // density-window watcher state (see PhysicsParams::densityPrecursorLow/High)
   const bool densityWatch = (densityPrecursorLow > real_t(0)) || (densityPrecursorHigh > real_t(0));
   uint_t densNBelow = 0, densNAbove = 0; real_t densMin = std::numeric_limits< real_t >::max(), densMax = std::numeric_limits< real_t >::lowest();
   std::string densFirst, densLowest, densHighest;

   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* pdfField  = blockIt->getData< const PdfField_T >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto* particleFieldPtr = blockIt->getData< const lbm_mesapd_coupling::ParticleField_T >(simObjs->particleFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const FlagUID formerMOuid("former moving obstacle");
      const flag_t formerMOFlag = flagField->flagExists(formerMOuid) ? flagField->getFlag(formerMOuid) : flag_t(0);
      const auto& aabb = blockIt->getAABB();

      // short label for a cell, for the neighbourhood dump on NaN detection
      auto classify = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> std::string {
         const flag_t fl = flagField->get(cx, cy, cz);
         if ((fl & obstacleMask) != flag_t(0))                        return "OBST";
         if (isFlagSet(fl, interfaceFlag))                            return "INTF";
         if (isFlagSet(fl, liquidFlag))                               return "LIQ ";
         if (isFlagSet(fl, gasFlag))                                  return "GAS ";
         if (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag))return "fMO ";
         return "??? ";
      };

      bool blockDone = false;
      for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()) && !blockDone; ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(pdfField->ySize()) && !blockDone; ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(pdfField->xSize()) && !blockDone; ++ix)
      {
         // Only liquid/interface cells are checked: a NaN there is genuine corruption. NaN PDFs in gas cells at an
         // obstacle's waterline are expected and harmless (SimpleBB skips those dry links).
         if (!flagField->isFlagSet(ix, iy, iz, liquidFlag) &&
             !flagField->isFlagSet(ix, iy, iz, interfaceFlag)) continue;

         // Low-rho precursor (PhysicsParams.lowRhoPrecursorThreshold, 0 = off): flags cells whose density is still
         // finite but below the threshold, i.e. before the isfinite scan below would catch them. Drained interface
         // cells (rho -> 0) at a floater's waterline are expected (see SimpleBB.h); this shows whether a cell that
         // later blows up can be told apart from those. Throttled to 20 hits per rank; meant for targeted reruns
         // (denseNanScanEnabled + From/ToStep around a known crash), not production.
         if (lowRhoPrecursorThreshold > real_t(0))
         {
            real_t rhoPrecheck = real_t(0);
            for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f) rhoPrecheck += pdfField->get(ix, iy, iz, f);

            if (std::isfinite(rhoPrecheck) && rhoPrecheck < lowRhoPrecursorThreshold)
            {
               static uint_t lowRhoPrecursorLogCount = uint_t(0);
               if (lowRhoPrecursorLogCount < uint_t(20))
               {
                  const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
                  const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
                  const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
                  WALBERLA_LOG_WARNING("[LOW-RHO-PRECURSOR t=" << t << " after=" << location << "] "
                     << classify(ix, iy, iz) << " cell block-local=(" << ix << "," << iy << "," << iz
                     << ") world=(" << wx << "," << wy << "," << wz << ") rho=" << rhoPrecheck
                     << " (threshold=" << lowRhoPrecursorThreshold << ") fill=" << fillField->get(ix, iy, iz));
                  ++lowRhoPrecursorLogCount;
               }
            }
         }

         // High-velocity precursor (PhysicsParams.highVelocityPrecursorThreshold, 0 = off): a blow-up can also
         // show as a large but finite velocity (local Mach number beyond the stability bound) before the PDFs
         // overflow. getVelocity() divides by the cell's actual rho (compressible model). Same throttling and
         // intended use as the low-rho check.
         if (highVelocityPrecursorThreshold > real_t(0))
         {
            const Vector3< real_t > velPrecheck = pdfField->getVelocity(ix, iy, iz);
            const real_t velMagSq = velPrecheck.sqrLength();

            if (std::isfinite(velMagSq) && velMagSq > highVelocityPrecursorThreshold * highVelocityPrecursorThreshold)
            {
               static uint_t highVelPrecursorLogCount = uint_t(0);
               if (highVelPrecursorLogCount < uint_t(20))
               {
                  const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
                  const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
                  const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
                  WALBERLA_LOG_WARNING("[HIGH-VELOCITY-PRECURSOR t=" << t << " after=" << location << "] "
                     << classify(ix, iy, iz) << " cell block-local=(" << ix << "," << iy << "," << iz
                     << ") world=(" << wx << "," << wy << "," << wz << ") |u|=" << std::sqrt(velMagSq)
                     << " (threshold=" << highVelocityPrecursorThreshold << ") vel=" << velPrecheck
                     << " fill=" << fillField->get(ix, iy, iz));
                  ++highVelPrecursorLogCount;
               }
            }
         }

         // Density WINDOW watcher (PhysicsParams.densityPrecursorLow/High): rho = sum_i f_i outside [low, high]
         if (densityWatch && (!densityPrecursorLiquidOnly || flagField->isFlagSet(ix, iy, iz, liquidFlag)))
         {
            real_t rhoW = real_t(0);
            for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f) rhoW += pdfField->get(ix, iy, iz, f);
            const bool below = densityPrecursorLow > real_t(0) && rhoW < densityPrecursorLow, above = densityPrecursorHigh > real_t(0) && rhoW > densityPrecursorHigh;
            if (std::isfinite(rhoW) && (below || above))
            {
               densNBelow += below; densNAbove += above;
               auto describe = [&]() {
                  real_t maxAbs = real_t(0); for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f) { const real_t v = pdfField->get(ix, iy, iz, f); if (std::isfinite(v)) maxAbs = std::max(maxAbs, std::abs(v)); }
                  uint_t pairs = 0, obst = 0, fmoN = 0, intfN = 0, gasN = 0; std::set< walberla::id_t > owners;
                  for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
                  {
                     if (dx == 0 && dy == 0 && dz == 0) continue;
                     const flag_t fl = flagField->get(ix + dx, iy + dy, iz + dz); const bool o1 = (fl & obstacleMask) != flag_t(0);
                     obst += o1; if (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag)) ++fmoN; if (isFlagSet(fl, interfaceFlag)) ++intfN; if (isFlagSet(fl, gasFlag)) ++gasN;
                     if (o1 && owners.size() < 4) owners.insert(particleFieldPtr->get(ix + dx, iy + dy, iz + dz));
                     if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0)) pairs += (o1 && ((flagField->get(ix - dx, iy - dy, iz - dz) & obstacleMask) != flag_t(0)));
                  }
                  std::ostringstream os;
                  os << classify(ix, iy, iz) << " cell block-local=(" << ix << "," << iy << "," << iz << ") world=(" << aabb.xMin() + real_c(ix) + real_t(0.5) << "," << aabb.yMin() + real_c(iy) + real_t(0.5) << ","
                     << aabb.zMin() + real_c(iz) + real_t(0.5) << ") rho=" << rhoW << " fill=" << fillField->get(ix, iy, iz) << " max|f|=" << maxAbs << " pinched pairs=" << pairs << "/13, obstacle neighbours=" << obst
                     << "/26, fMO neighbours=" << fmoN << ", interface neighbours=" << intfN << ", gas neighbours=" << gasN << ", obstacle owners=[ ";
                  for (const auto u : owners) os << u << " ";
                  os << "]";
                  return os.str();
               };
               if (densNBelow + densNAbove == uint_t(1)) densFirst = describe();
               if (rhoW < densMin) { densMin = rhoW; densLowest = describe(); }
               if (rhoW > densMax) { densMax = rhoW; densHighest = describe(); }
            }
         }

         // population precursor check (PhysicsParams.populationPrecursorThreshold, 0 = off), see its doc comment
         if (populationPrecursorThreshold > real_t(0))
         {
            real_t maxAbs = real_t(0); uint_t maxIdx = uint_t(0);
            for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
            {
               const real_t v = pdfField->get(ix, iy, iz, f);
               if (std::isfinite(v) && std::abs(v) > maxAbs) { maxAbs = std::abs(v); maxIdx = f; }
            }
            if (maxAbs > populationPrecursorThreshold)
            {
               static uint_t popPrecursorLogCount = uint_t(0);
               if (popPrecursorLogCount < uint_t(40))
               {
                  ++popPrecursorLogCount;
                  uint_t pinchedPairs = 0, obstNeighbours = 0;
                  for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
                  {
                     if (dx == 0 && dy == 0 && dz == 0) continue;
                     const bool o1 = (flagField->get(ix + dx, iy + dy, iz + dz) & obstacleMask) != flag_t(0);
                     obstNeighbours += o1;
                     if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0))
                        pinchedPairs += (o1 && ((flagField->get(ix - dx, iy - dy, iz - dz) & obstacleMask) != flag_t(0)));
                  }
                  // neighbourhood context: former-moving-obstacle cells (just left by a body), interface/gas
                  // neighbours, fastest fluid neighbour, owners of the obstacle neighbours
                  uint_t fmoNeighbours = 0, interfaceNeighbours = 0, gasNeighbours = 0; real_t maxNbrSpeed = real_t(0); std::set< walberla::id_t > owners;
                  for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
                  {
                     if (dx == 0 && dy == 0 && dz == 0) continue;
                     const cell_idx_t cx = ix + dx, cy = iy + dy, cz = iz + dz; const flag_t fl = flagField->get(cx, cy, cz);
                     if (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag)) ++fmoNeighbours;
                     if (isFlagSet(fl, interfaceFlag)) ++interfaceNeighbours;
                     if (isFlagSet(fl, gasFlag)) ++gasNeighbours;
                     if ((fl & obstacleMask) != flag_t(0)) { if (owners.size() < 4) owners.insert(particleFieldPtr->get(cx, cy, cz)); }
                     else if (isFlagSet(fl, liquidFlag) || isFlagSet(fl, interfaceFlag)) { const real_t sp = pdfField->getVelocity(cx, cy, cz).length(); if (std::isfinite(sp)) maxNbrSpeed = std::max(maxNbrSpeed, sp); }
                  }
                  std::ostringstream ownerStr; for (const auto u : owners) ownerStr << u << " ";
                  real_t rho = real_t(0); for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f) rho += pdfField->get(ix, iy, iz, f);
                  const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5), wy = aabb.yMin() + real_c(iy) + real_t(0.5), wz = aabb.zMin() + real_c(iz) + real_t(0.5);
                  WALBERLA_LOG_WARNING("[POP-PRECURSOR t=" << t << " after=" << location << "] " << classify(ix, iy, iz) << " cell block-local=(" << ix << "," << iy << "," << iz
                     << ") world=(" << wx << "," << wy << "," << wz << ") max|f|=" << maxAbs << " at f[" << maxIdx << "] (threshold=" << populationPrecursorThreshold << ") rho=" << rho
                     << " fill=" << fillField->get(ix, iy, iz) << " | obstacle geometry: pinched opposite pairs=" << pinchedPairs << "/13, obstacle neighbours=" << obstNeighbours << "/26"
                     << " | context: former-moving-obstacle neighbours=" << fmoNeighbours << ", interface neighbours=" << interfaceNeighbours << ", gas neighbours=" << gasNeighbours
                     << ", fastest fluid neighbour |u|=" << maxNbrSpeed << ", obstacle owners (uids)=[ " << ownerStr.str() << "]");
               }
            }
         }

         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
         {
            const real_t val = pdfField->get(ix, iy, iz, f);
            if (std::isfinite(val)) continue;

            const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
            const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
            const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
            WALBERLA_LOG_WARNING("[NaN-SCAN t=" << t << " after=" << location << "] FATAL fluid-cell NaN"
               << " block-local=(" << ix << "," << iy << "," << iz << ")"
               << " world=(" << wx << "," << wy << "," << wz << ")"
               << " f[" << f << "]=" << val << " flag=" << classify(ix, iy, iz));

            // Dump the 3x3x3 neighbourhood: per cell  [flag rho fill]  — the config at the failure.
            WALBERLA_LOG_WARNING("[NaN-SCAN t=" << t << "] 3x3x3 neighbourhood  [flag rho fill]:");
            for (int ddz = -1; ddz <= 1; ++ddz)
               for (int ddy = -1; ddy <= 1; ++ddy)
               {
                  std::ostringstream row;
                  row << std::scientific << std::setprecision(3);
                  for (int ddx = -1; ddx <= 1; ++ddx)
                  {
                     const cell_idx_t cx = ix + cell_idx_c(ddx);
                     const cell_idx_t cy = iy + cell_idx_c(ddy);
                     const cell_idx_t cz = iz + cell_idx_c(ddz);
                     real_t r = real_t(0);
                     for (uint_t ff = 0; ff < PdfField_T::Stencil::Size; ++ff) r += pdfField->get(cx, cy, cz, ff);
                     row << "[" << classify(cx, cy, cz) << " " << r << " " << fillField->get(cx, cy, cz) << "]  ";
                  }
                  WALBERLA_LOG_WARNING("   dz=" << ddz << " dy=" << ddy << "  (dx=-1,0,+1):  " << row.str());
               }

            localNaN = true;
            blockDone = true;
            break;
         }

         // Velocity-derived NaN check: getVelocity() = momentum/density gives 0/0 = NaN if all PDFs are finite but
         // sum to zero density, which the isfinite scan above cannot see. MeanVelocityComputer (Forcing.h) sums
         // getVelocity() over all liquid/interface cells, so such a cell would corrupt the mean velocity.
         if (!blockDone)
         {
            const Vector3< real_t > velCheck = pdfField->getVelocity(ix, iy, iz);
            if (!std::isfinite(velCheck[0]) || !std::isfinite(velCheck[1]) || !std::isfinite(velCheck[2]))
            {
               real_t rhoCheck = real_t(0);
               for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f) rhoCheck += pdfField->get(ix, iy, iz, f);

               const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
               const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
               const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
               WALBERLA_LOG_WARNING("[VEL-NaN-SCAN t=" << t << " after=" << location << "] FATAL fluid-cell velocity NaN "
                  "(all raw PDFs individually finite -- this is a 0/0-in-getVelocity() case the plain PDF isfinite scan cannot see)"
                  << " block-local=(" << ix << "," << iy << "," << iz << ")"
                  << " world=(" << wx << "," << wy << "," << wz << ")"
                  << " vel=" << velCheck << " rho=" << rhoCheck << " fill=" << fillField->get(ix, iy, iz)
                  << " flag=" << classify(ix, iy, iz));

               WALBERLA_LOG_WARNING("[VEL-NaN-SCAN t=" << t << "] 3x3x3 neighbourhood  [flag rho fill]:");
               for (int ddz = -1; ddz <= 1; ++ddz)
                  for (int ddy = -1; ddy <= 1; ++ddy)
                  {
                     std::ostringstream row;
                     row << std::scientific << std::setprecision(3);
                     for (int ddx = -1; ddx <= 1; ++ddx)
                     {
                        const cell_idx_t cx = ix + cell_idx_c(ddx);
                        const cell_idx_t cy = iy + cell_idx_c(ddy);
                        const cell_idx_t cz = iz + cell_idx_c(ddz);
                        real_t r = real_t(0);
                        for (uint_t ff = 0; ff < PdfField_T::Stencil::Size; ++ff) r += pdfField->get(cx, cy, cz, ff);
                        row << "[" << classify(cx, cy, cz) << " " << r << " " << fillField->get(cx, cy, cz) << "]  ";
                     }
                     WALBERLA_LOG_WARNING("   dz=" << ddz << " dy=" << ddy << "  (dx=-1,0,+1):  " << row.str());
                  }

               localNaN  = true;
               blockDone = true;
            }
         }
      }
   }

   if (densityWatch)
   {
      if (densNBelow + densNAbove > uint_t(0))
      {
         static uint_t densLogCount = uint_t(0);
         if (densLogCount < uint_t(60))
         {
            ++densLogCount;
            WALBERLA_LOG_WARNING("[DENSITY-PRECURSOR t=" << t << " after=" << location << "] this rank: " << densNBelow << " cells with rho < " << densityPrecursorLow << ", " << densNAbove << " cells with rho > "
               << densityPrecursorHigh << " | FIRST offending: " << densFirst << " | LOWEST: " << densLowest << " | HIGHEST: " << densHighest);
         }
      }
      uint_t glob[2] = { densNBelow, densNAbove };
      mpi::allReduceInplace(glob[0], mpi::SUM); mpi::allReduceInplace(glob[1], mpi::SUM);
      if (glob[0] + glob[1] > uint_t(0))
      {
         static uint_t densRootCount = uint_t(0);
         if (densRootCount < uint_t(400)) { ++densRootCount; WALBERLA_LOG_WARNING_ON_ROOT("[DENSITY-WINDOW t=" << t << " after=" << location << "] cells outside rho in [" << densityPrecursorLow << ", " << densityPrecursorHigh
            << "] over all ranks: " << glob[0] << " below, " << glob[1] << " above" << (densityPrecursorLiquidOnly ? " (liquid cells only)" : " (liquid + interface cells)")); }
      }
   }

   // reduce across all ranks so every rank knows whether any NaN exists globally
   uint_t globalNaN = localNaN ? uint_t(1) : uint_t(0);
   mpi::allReduceInplace< uint_t >(globalNaN, mpi::SUM);
   return globalNaN > uint_t(0);
}

// Repairs obstacle cells with degenerate PDFs (non-finite or non-positive density). streamReconstructAdvect never
// writes PDFs of obstacle cells (it only carries the fill level), and its ping-pong scratch buffer comes from
// cloneUninitialized(), so an obstacle cell that SimpleBB's bounce-back never refreshes (e.g. one facing only other
// obstacle cells) can hold all-zero PDFs. collide then computes m_000=0 -> xi_0=1/m_000=inf and every cumulant
// becomes NaN.
//
// Unlike checkPdfNaN (liquid/interface cells, aborts), this scans obstacle cells and resets a degenerate one to
// the D3Q27 rest equilibrium (rho=1, u=0). Nothing but collide reads an obstacle cell's own PDFs:
// SimpleBB::treatDirection reads the fluid-side cell, PdfReconstructionManager's seedConvertedCellFill averages
// neighbours, and MovingParticleMapping does not touch the PDF field. Healthy cells are never modified.
// Returns the number of cells repaired on this rank.
inline uint_t repairDegenerateObstacleCells(const std::shared_ptr< SimulationObjects >& simObjs, const char* const location, const uint_t t)
{
   uint_t totalRepaired = 0;
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      auto* pdfField         = blockIt->getData< PdfField_T >(simObjs->pdfFieldID);
      const auto* flagField  = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto obstacleMask = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const auto& aabb = blockIt->getAABB();

      for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()); ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(pdfField->ySize()); ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(pdfField->xSize()); ++ix)
      {
         if ((flagField->get(ix, iy, iz) & obstacleMask) == flag_t(0)) continue; // not an obstacle cell -- untouched

         real_t rho = real_t(0);
         for (uint_t q = 0; q < PdfField_T::Stencil::Size; ++q) rho += pdfField->get(ix, iy, iz, q);

         if (std::isfinite(rho) && rho > real_t(0)) continue; // healthy -- leave completely alone

         static uint_t repairLogCount = uint_t(0);
         ++repairLogCount;
         if (repairLogCount <= uint_t(20) || (repairLogCount % uint_t(5000) == uint_t(0)))
         {
            const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
            const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
            const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
            WALBERLA_LOG_WARNING("[OBSTACLE-CELL-REPAIR] (#" << repairLogCount << ") t=" << t << " after=" << location
               << " block-local=(" << ix << "," << iy << "," << iz << ") world=(" << wx << "," << wy << "," << wz
               << ") had rho=" << rho << " -- resetting to rest equilibrium (rho=1, u=0)");
         }

         pdfField->get(ix, iy, iz, uint_t(0)) = real_t(8.0 / 27.0);
         for (uint_t q = uint_t(1); q <= uint_t(6);  ++q) pdfField->get(ix, iy, iz, q) = real_t(2.0 / 27.0);
         for (uint_t q = uint_t(7); q <= uint_t(18); ++q) pdfField->get(ix, iy, iz, q) = real_t(1.0 / 54.0);
         for (uint_t q = uint_t(19); q <= uint_t(26); ++q) pdfField->get(ix, iy, iz, q) = real_t(1.0 / 216.0);

         ++totalRepaired;
      }
   }
   return totalRepaired;
}

// Fluid sanity sweep (PhysicsParams::fluidSanity*, see SimulationParams.h): per-step range check of every liquid
// and interface cell (fluidSanityLiquidOnly drops interface cells; gas, obstacle, former-moving-obstacle and
// boundary cells are never checked). Returns the global number of violating cells.
// Log throttling: the first 30 violating steps are printed, afterwards only steps whose count at least doubles
// relative to the last printed one (or every 100th step), so a steady benign violation cannot use up the budget
// and hide an escalation.
inline bool fluidSanityShouldLog(uint_t& printed, uint_t& lastCount, const uint_t count, const uint_t t, const uint_t cap)
{
   if (printed >= cap) return false;
   if (printed < uint_t(30) || count >= uint_t(2) * lastCount || t % uint_t(100) == uint_t(0)) { ++printed; lastCount = count; return true; }
   return false;
}

inline uint_t fluidSanitySweep(const std::shared_ptr< SimulationObjects >& simObjs, const std::shared_ptr< const SimulationParams >& simParams, const uint_t t)
{
   const auto& P = simParams->physics;
   struct Cat { uint_t nonFinite = 0, rhoLow = 0, rhoHigh = 0, velHigh = 0, popHigh = 0, fillBad = 0, excess = 0, velJump = 0, rhoJump = 0, topology = 0; } cat;
   uint_t violating = 0, checked = 0;
   real_t rhoMin = std::numeric_limits< real_t >::max(), uMax = real_t(0), fMax = real_t(0);            // among violating cells
   real_t sRhoMin = std::numeric_limits< real_t >::max(), sRhoMax = std::numeric_limits< real_t >::lowest(), sUMax = real_t(0), sFMax = real_t(0), sFillMin = std::numeric_limits< real_t >::max(), sFillMax = std::numeric_limits< real_t >::lowest();
   std::string first, lowest, fastest, biggest, statsFastest;
   const bool doJumps = P.fluidSanityMaxNeighbourVelocityJump > real_t(0) || P.fluidSanityMaxNeighbourDensityJump > real_t(0);
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* pdfField   = blockIt->getData< const PdfField_T >(simObjs->pdfFieldID);
      const auto* flagField  = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField  = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto* excessField = blockIt->getData< const ScalarField_T >(simObjs->excessMassFieldID);
      const auto* particleFieldPtr = blockIt->getData< const lbm_mesapd_coupling::ParticleField_T >(simObjs->particleFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const FlagUID formerMOuid("former moving obstacle");
      const flag_t formerMOFlag = flagField->flagExists(formerMOuid) ? flagField->getFlag(formerMOuid) : flag_t(0);
      const auto& aabb = blockIt->getAABB();
      const cell_idx_t nx = cell_idx_c(pdfField->xSize()), ny = cell_idx_c(pdfField->ySize()), nz = cell_idx_c(pdfField->zSize());
      auto klass = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> const char* {
         const flag_t fl = flagField->get(cx, cy, cz);
         if ((fl & obstacleMask) != flag_t(0)) return "OBST"; if (isFlagSet(fl, interfaceFlag)) return "INTF"; if (isFlagSet(fl, liquidFlag)) return "LIQ"; if (isFlagSet(fl, gasFlag)) return "GAS";
         if (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag)) return "fMO"; return "???"; };
      auto describe = [&](cell_idx_t ix, cell_idx_t iy, cell_idx_t iz, real_t rho, real_t maxF, real_t speed) {
         uint_t pairs = 0, obst = 0, fmoN = 0, intfN = 0, gasN = 0; std::set< walberla::id_t > owners;
         for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
         {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const flag_t fl = flagField->get(ix + dx, iy + dy, iz + dz); const bool o1 = (fl & obstacleMask) != flag_t(0);
            obst += o1; if (formerMOFlag != flag_t(0) && isFlagSet(fl, formerMOFlag)) ++fmoN; if (isFlagSet(fl, interfaceFlag)) ++intfN; if (isFlagSet(fl, gasFlag)) ++gasN;
            if (o1 && owners.size() < 4) owners.insert(particleFieldPtr->get(ix + dx, iy + dy, iz + dz));
            if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0)) pairs += (o1 && ((flagField->get(ix - dx, iy - dy, iz - dz) & obstacleMask) != flag_t(0)));
         }
         std::ostringstream os;
         os << klass(ix, iy, iz) << " block-local=(" << ix << "," << iy << "," << iz << ") world=(" << aabb.xMin() + real_c(ix) + real_t(0.5) << "," << aabb.yMin() + real_c(iy) + real_t(0.5) << ","
            << aabb.zMin() + real_c(iz) + real_t(0.5) << ") rho=" << rho << " fill=" << fillField->get(ix, iy, iz) << " excessMass=" << excessField->get(ix, iy, iz) << " max|f|=" << maxF << " |u|=" << speed
            << " pinched pairs=" << pairs << "/13, obstacle nbrs=" << obst << "/26, fMO nbrs=" << fmoN << ", interface nbrs=" << intfN << ", gas nbrs=" << gasN << ", owners=[ ";
         for (const auto u : owners) os << u << " ";
         os << "]"; return os.str(); };
      // pass 1: moments of every checked (liquid / interface) cell of this block's interior
      const uint_t nCells = uint_c(nx) * uint_c(ny) * uint_c(nz);
      std::vector< uint8_t > isChecked(nCells, uint8_t(0));
      std::vector< real_t > rhoV(nCells, real_t(0)), maxFV(nCells, real_t(0));
      std::vector< Vector3< real_t > > uV(nCells, Vector3< real_t >(real_t(0)));
      std::vector< uint8_t > finiteV(nCells, uint8_t(1));
      auto idx = [&](cell_idx_t ix, cell_idx_t iy, cell_idx_t iz) { return (uint_c(iz) * uint_c(ny) + uint_c(iy)) * uint_c(nx) + uint_c(ix); };
      for (cell_idx_t iz = 0; iz < nz; ++iz) for (cell_idx_t iy = 0; iy < ny; ++iy) for (cell_idx_t ix = 0; ix < nx; ++ix)
      {
         const bool isLiquid = flagField->isFlagSet(ix, iy, iz, liquidFlag), isInterface = flagField->isFlagSet(ix, iy, iz, interfaceFlag);
         if (!isLiquid && !(isInterface && !P.fluidSanityLiquidOnly)) continue;
         const uint_t c = idx(ix, iy, iz);
         isChecked[c] = uint8_t(1); ++checked;
         real_t rho = real_t(0), maxF = real_t(0); bool finite = true;
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f) { const real_t v = pdfField->get(ix, iy, iz, f); rho += v; if (!std::isfinite(v)) finite = false; else maxF = std::max(maxF, std::abs(v)); }
         rhoV[c] = rho; maxFV[c] = maxF; finiteV[c] = uint8_t(finite);
         if (finite) uV[c] = pdfField->getVelocity(ix, iy, iz);
      }
      // pass 2: checks
      for (cell_idx_t iz = 0; iz < nz; ++iz) for (cell_idx_t iy = 0; iy < ny; ++iy) for (cell_idx_t ix = 0; ix < nx; ++ix)
      {
         const uint_t c = idx(ix, iy, iz);
         if (!isChecked[c]) continue;
         const bool isInterface = flagField->isFlagSet(ix, iy, iz, interfaceFlag);
         const real_t rho = rhoV[c], maxF = maxFV[c], speed = finiteV[c] ? uV[c].length() : real_t(0);
         const real_t fill = fillField->get(ix, iy, iz), excess = excessField->get(ix, iy, iz), tol = P.fluidSanityFillTolerance;
         if (finiteV[c] && std::isfinite(rho)) { sRhoMin = std::min(sRhoMin, rho); sRhoMax = std::max(sRhoMax, rho); sFMax = std::max(sFMax, maxF); }
         if (std::isfinite(fill)) { sFillMin = std::min(sFillMin, fill); sFillMax = std::max(sFillMax, fill); }
         if (finiteV[c] && std::isfinite(speed) && speed > sUMax) { sUMax = speed; if (P.fluidSanityStatsEveryNSteps > uint_t(0)) statsFastest = describe(ix, iy, iz, rho, maxF, speed); }
         bool bad = false;
         if (!finiteV[c] || !std::isfinite(rho) || !std::isfinite(speed) || !std::isfinite(fill) || !std::isfinite(excess)) { ++cat.nonFinite; bad = true; }
         else
         {
            if (P.fluidSanityDensityMin > real_t(0) && rho < P.fluidSanityDensityMin) { ++cat.rhoLow; bad = true; }
            if (P.fluidSanityDensityMax > real_t(0) && rho > P.fluidSanityDensityMax) { ++cat.rhoHigh; bad = true; }
            if (P.fluidSanityMaxVelocity > real_t(0) && speed > P.fluidSanityMaxVelocity) { ++cat.velHigh; bad = true; }
            if (P.fluidSanityMaxPopulation > real_t(0) && maxF > P.fluidSanityMaxPopulation) { ++cat.popHigh; bad = true; }
            if (tol >= real_t(0) && (isInterface ? (fill < -tol || fill > real_t(1) + tol) : (fill < real_t(1) - tol || fill > real_t(1) + tol))) { ++cat.fillBad; bad = true; }
            if (P.fluidSanityMaxExcessMass > real_t(0) && std::abs(excess) > P.fluidSanityMaxExcessMass) { ++cat.excess; bad = true; }
            if (doJumps)
            {   // largest jump to a checked face neighbour inside this block (ghost-layer PDFs are half a step older)
               real_t dU = real_t(0), dRho = real_t(0);
               static const int nb[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
               for (const auto& d : nb)
               {
                  const cell_idx_t jx = ix + d[0], jy = iy + d[1], jz = iz + d[2];
                  if (jx < 0 || jy < 0 || jz < 0 || jx >= nx || jy >= ny || jz >= nz) continue;
                  const uint_t n = idx(jx, jy, jz);
                  if (!isChecked[n] || !finiteV[n]) continue;
                  dU = std::max(dU, (uV[n] - uV[c]).length()); dRho = std::max(dRho, std::abs(rhoV[n] - rho));
               }
               if (P.fluidSanityMaxNeighbourVelocityJump > real_t(0) && dU > P.fluidSanityMaxNeighbourVelocityJump) { ++cat.velJump; bad = true; }
               if (P.fluidSanityMaxNeighbourDensityJump > real_t(0) && dRho > P.fluidSanityMaxNeighbourDensityJump) { ++cat.rhoJump; bad = true; }
            }
         }
         if (P.fluidSanityCheckTopology && !isInterface)
         {   // a liquid cell must not touch a gas cell (closed interface layer over the D3Q27 stencil)
            bool touchesGas = false;
            for (int dz = -1; dz <= 1 && !touchesGas; ++dz) for (int dy = -1; dy <= 1 && !touchesGas; ++dy) for (int dx = -1; dx <= 1 && !touchesGas; ++dx)
               if ((dx != 0 || dy != 0 || dz != 0) && flagField->isFlagSet(ix + dx, iy + dy, iz + dz, gasFlag)) touchesGas = true;
            if (touchesGas) { ++cat.topology; bad = true; }
         }
         if (!bad) continue;
         ++violating;
         if (violating == uint_t(1)) first = describe(ix, iy, iz, rho, maxF, speed);
         if (std::isfinite(rho) && rho < rhoMin) { rhoMin = rho; lowest = describe(ix, iy, iz, rho, maxF, speed); }
         if (std::isfinite(speed) && speed > uMax) { uMax = speed; fastest = describe(ix, iy, iz, rho, maxF, speed); }
         if (maxF > fMax) { fMax = maxF; biggest = describe(ix, iy, iz, rho, maxF, speed); }
      }
   }
   static uint_t rankPrinted = uint_t(0), rankLast = uint_t(0), rootPrinted = uint_t(0), rootLast = uint_t(0);
   if (violating > uint_t(0) && fluidSanityShouldLog(rankPrinted, rankLast, violating, t, uint_t(400)))
      WALBERLA_LOG_WARNING("[FLUID-SANITY t=" << t << "] this rank: " << violating << " of " << checked << " checked cells violate the ranges (non-finite " << cat.nonFinite << ", rho<" << P.fluidSanityDensityMin << " " << cat.rhoLow
         << ", rho>" << P.fluidSanityDensityMax << " " << cat.rhoHigh << ", |u|>" << P.fluidSanityMaxVelocity << " " << cat.velHigh << ", max|f|>" << P.fluidSanityMaxPopulation << " " << cat.popHigh << ", fill " << cat.fillBad
         << ", |excess|>" << P.fluidSanityMaxExcessMass << " " << cat.excess << ", face-nbr |du|>" << P.fluidSanityMaxNeighbourVelocityJump << " " << cat.velJump << ", face-nbr |drho|>" << P.fluidSanityMaxNeighbourDensityJump << " " << cat.rhoJump
         << ", liquid touching gas " << cat.topology << ") | FIRST: " << first << " | LOWEST rho: " << lowest << " | FASTEST: " << fastest << " | LARGEST |f|: " << biggest);
   uint_t g[11] = { violating, cat.nonFinite, cat.rhoLow, cat.rhoHigh, cat.velHigh, cat.popHigh, cat.fillBad, cat.excess, cat.velJump, cat.rhoJump, cat.topology };
   for (int i = 0; i < 11; ++i) mpi::allReduceInplace(g[i], mpi::SUM);
   if (g[0] > uint_t(0) && fluidSanityShouldLog(rootPrinted, rootLast, g[0], t, uint_t(2000)))
      WALBERLA_LOG_WARNING_ON_ROOT("[FLUID-SANITY-SUMMARY t=" << t << "] " << g[0] << " violating cells over all ranks: non-finite " << g[1] << ", rho<" << P.fluidSanityDensityMin << " " << g[2]
         << ", rho>" << P.fluidSanityDensityMax << " " << g[3] << ", |u|>" << P.fluidSanityMaxVelocity << " " << g[4] << ", max|f|>" << P.fluidSanityMaxPopulation << " " << g[5] << ", fill outside range " << g[6]
         << ", |excess|>" << P.fluidSanityMaxExcessMass << " " << g[7] << ", face-nbr |du|>" << P.fluidSanityMaxNeighbourVelocityJump << " " << g[8] << ", face-nbr |drho|>" << P.fluidSanityMaxNeighbourDensityJump << " " << g[9]
         << ", liquid touching gas " << g[10]);
   if (g[0] > uint_t(0) && P.fluidSanityAbortAbove > uint_t(0) && g[0] > P.fluidSanityAbortAbove)
      WALBERLA_ABORT("[FLUID-SANITY] aborting at t=" << t << ": " << g[0] << " violating cells > fluidSanityAbortAbove=" << P.fluidSanityAbortAbove);
   // global extrema time series, independent of violations: shows the growth of max|u| / rho before a
   // threshold is crossed
   if (P.fluidSanityStatsEveryNSteps > uint_t(0) && t % P.fluidSanityStatsEveryNSteps == uint_t(0))
   {
      uint_t nChecked = checked; mpi::allReduceInplace(nChecked, mpi::SUM);
      real_t gRhoMin = sRhoMin, gRhoMax = sRhoMax, gU = sUMax, gF = sFMax, gFillMin = sFillMin, gFillMax = sFillMax;
      mpi::allReduceInplace(gRhoMin, mpi::MIN); mpi::allReduceInplace(gRhoMax, mpi::MAX); mpi::allReduceInplace(gU, mpi::MAX); mpi::allReduceInplace(gF, mpi::MAX);
      mpi::allReduceInplace(gFillMin, mpi::MIN); mpi::allReduceInplace(gFillMax, mpi::MAX);
      WALBERLA_LOG_INFO_ON_ROOT("[FLUID-STATS t=" << t << "] checked cells " << nChecked << ", rho [" << gRhoMin << ", " << gRhoMax << "], max|u| " << gU << ", max|f_i| " << gF << ", fill [" << gFillMin << ", " << gFillMax << "]");
      if (sUMax >= gU && gU > real_t(0)) WALBERLA_LOG_INFO("[FLUID-STATS t=" << t << "] fastest cell: " << statsFastest);
   }
   return g[0];
}

// Same scan as checkPdfNaN but without its early exit per block: lists every liquid/interface cell with a
// non-finite PDF (world and block-local coordinates only, to keep the output compact). Called before aborting
// on a checkPdfNaN hit, to show the full extent of the corruption.
void listAllBadCellsNoEarlyExit(const std::shared_ptr< SimulationObjects >& simObjs, const char* const location)
{
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* pdfField  = blockIt->getData< const PdfField_T >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto& aabb = blockIt->getAABB();

      for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()); ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(pdfField->ySize()); ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(pdfField->xSize()); ++ix)
      {
         if (!flagField->isFlagSet(ix, iy, iz, liquidFlag) && !flagField->isFlagSet(ix, iy, iz, interfaceFlag)) continue;
         bool bad = false;
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
            if (!std::isfinite(pdfField->get(ix, iy, iz, f))) { bad = true; break; }
         if (!bad) continue;
         const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
         const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
         const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
         WALBERLA_LOG_WARNING("[FULL-BAD-CELL-LIST at=" << location << "] world=(" << wx << "," << wy << "," << wz
            << ") local=(" << ix << "," << iy << "," << iz << ")");
      }
   }
}

// ModelParameters.highAltitudeWaterScanEnabled/highAltitudeWaterThresholdZ (see SimulationParams.h): warns for
// every local liquid or interface cell above thresholdZ. CellConversionSweep logs nothing when a gas cell
// converts, so this gives the step and cell where water first appears at an unexpected height. Local cells only
// (as in checkPdfNaN), so one conversion is not reported twice across a ghost layer. Does not abort.
void checkHighAltitudeWater(const std::shared_ptr< SimulationObjects >& simObjs,
                            const char* const location, const uint_t t, const real_t thresholdZ)
{
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto& aabb = blockIt->getAABB();

      // skip blocks that lie entirely below the threshold
      if (aabb.zMax() <= thresholdZ) continue;

      for (cell_idx_t iz = 0; iz < cell_idx_c(flagField->zSize()); ++iz)
      {
         const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
         if (wz <= thresholdZ) continue;

         for (cell_idx_t iy = 0; iy < cell_idx_c(flagField->ySize()); ++iy)
         for (cell_idx_t ix = 0; ix < cell_idx_c(flagField->xSize()); ++ix)
         {
            const bool isLiquid    = flagField->isFlagSet(ix, iy, iz, liquidFlag);
            const bool isInterface = flagField->isFlagSet(ix, iy, iz, interfaceFlag);
            if (!isLiquid && !isInterface) continue;

            const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
            const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
            WALBERLA_LOG_WARNING("[HIGH-ALTITUDE-WATER t=" << t << " after=" << location << "] "
               << (isLiquid ? "LIQUID" : "INTERFACE") << " cell above threshold Z=" << thresholdZ
               << "  world=(" << wx << "," << wy << "," << wz << ")"
               << "  block-local=(" << ix << "," << iy << "," << iz << ")"
               << "  fill=" << fillField->get(ix, iy, iz));
         }
      }
   }
}

// Check the closed-interface-layer invariant: in FSLBM a liquid cell may never be a (D3Q27) neighbour of a gas
// cell; an interface cell must separate them. A violation is a "hole" in the interface layer (the Fig. 4.12a
// configuration). Counts the holes, logs the first one with its 3x3x3 neighbourhood and returns the global
// count. Does not abort.
uint_t checkClosedInterfaceLayer(const std::shared_ptr< SimulationObjects >& simObjs,
                                 const char* const location, const uint_t t)
{
   uint_t localHoles = uint_t(0);
   bool   loggedFirst = false;

   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const auto& aabb = blockIt->getAABB();

      auto classify = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> std::string {
         const flag_t fl = flagField->get(cx, cy, cz);
         if ((fl & obstacleMask) != flag_t(0)) return "OBST";
         if (isFlagSet(fl, interfaceFlag))     return "INTF";
         if (isFlagSet(fl, liquidFlag))        return "LIQ ";
         if (isFlagSet(fl, gasFlag))           return "GAS ";
         return "??? ";
      };

      for (cell_idx_t iz = 0; iz < cell_idx_c(flagField->zSize()); ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(flagField->ySize()); ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(flagField->xSize()); ++ix)
      {
         if (!flagField->isFlagSet(ix, iy, iz, liquidFlag)) continue; // only liquid cells

         for (int ddz = -1; ddz <= 1; ++ddz)
         for (int ddy = -1; ddy <= 1; ++ddy)
         for (int ddx = -1; ddx <= 1; ++ddx)
         {
            if (ddx == 0 && ddy == 0 && ddz == 0) continue;
            if (!flagField->isFlagSet(ix + cell_idx_c(ddx), iy + cell_idx_c(ddy), iz + cell_idx_c(ddz), gasFlag))
               continue;

            ++localHoles;
            if (!loggedFirst)
            {
               loggedFirst = true;
               const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
               const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
               const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
               WALBERLA_LOG_WARNING("[LAYER-CHECK t=" << t << " after=" << location
                  << "] HOLE: liquid cell block-local=(" << ix << "," << iy << "," << iz << ") world=(" << wx
                  << "," << wy << "," << wz << ") has a GAS neighbour at offset (" << ddx << "," << ddy << ","
                  << ddz << "). Closed-interface-layer invariant broken.");
               for (int ez = -1; ez <= 1; ++ez)
                  for (int ey = -1; ey <= 1; ++ey)
                  {
                     std::ostringstream row;
                     row << std::scientific << std::setprecision(3);
                     for (int ex = -1; ex <= 1; ++ex)
                        row << "[" << classify(ix + cell_idx_c(ex), iy + cell_idx_c(ey), iz + cell_idx_c(ez)) << " fill="
                            << fillField->get(ix + cell_idx_c(ex), iy + cell_idx_c(ey), iz + cell_idx_c(ez)) << "]  ";
                     WALBERLA_LOG_WARNING("   dz=" << ez << " dy=" << ey << "  (dx=-1,0,+1):  " << row.str());
                  }
            }
            ddz = ddy = ddx = 2; // stop scanning neighbours of this liquid cell once a hole is found
         }
      }
   }

   uint_t globalHoles = localHoles;
   mpi::allReduceInplace< uint_t >(globalHoles, mpi::SUM);
   if (globalHoles > uint_t(0))
      WALBERLA_LOG_WARNING_ON_ROOT("[LAYER-CHECK t=" << t << " after=" << location << "] total holes (liquid–gas adjacencies) = "
                                   << globalHoles);
   return globalHoles;
}

// Measures the azimuthal asymmetry of the wetted vs. dry surface of each floating sphere (floaterUIDs), a
// diagnostic for lateral drift caused by the SimpleBB dry-link guard.
//
// The MEM force is integrated only over wet links (fluid-side cell carries a valid PDF); the guard skips dry
// links (gas / NaN fluid-side cell at the waterline) for both the force and the bounce-back. For a floating
// sphere the dry set is the cap above the waterline. If it is azimuthally symmetric, its horizontal components
// cancel; if waterline cells turn into gas asymmetrically (e.g. under streamwise flow), the wetted-surface force
// gains a net horizontal component and the sphere drifts sideways.
//
// For each floater, the outward link offsets from its obstacle cells (particleField == UID) to D3Q27 fluid
// neighbours are summed separately for wet (liquid/interface) and dry (gas) links; the horizontal part of the
// dry sum is the lateral-bias signal. Reduced to root and logged; does not abort.
void logFloaterLinkBalance(const std::shared_ptr< SimulationObjects >& simObjs,
                           const char* const location, const uint_t t)
{
   const uint_t nF = uint_c(simObjs->floaterUIDs.size());
   if (nF == uint_t(0)) return;

   std::vector< uint_t > nWet(nF, uint_t(0)), nDry(nF, uint_t(0));
   std::vector< real_t > sumWet(nF * 3, real_t(0)), sumDry(nF * 3, real_t(0));
   std::vector< real_t > pos(nF * 3, real_t(0)), vel(nF * 3, real_t(0));

   // --- link balance from the fields ---
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* flagField     = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* particleField = blockIt->getData< const lbm_mesapd_coupling::ParticleField_T >(simObjs->particleFieldID);
      const auto gasFlag = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);

      for (cell_idx_t iz = 0; iz < cell_idx_c(flagField->zSize()); ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(flagField->ySize()); ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(flagField->xSize()); ++ix)
      {
         const walberla::id_t uid = particleField->get(ix, iy, iz);
         if (uid == simObjs->particleAccessor->getInvalidUid()) continue; // only obstacle cells

         // map this cell's UID to a floater index (skip bed/plane particles)
         uint_t fi = nF;
         for (uint_t k = 0; k < nF; ++k) if (simObjs->floaterUIDs[k] == uid) { fi = k; break; }
         if (fi == nF) continue;

         for (int ddz = -1; ddz <= 1; ++ddz)
         for (int ddy = -1; ddy <= 1; ++ddy)
         for (int ddx = -1; ddx <= 1; ++ddx)
         {
            if (ddx == 0 && ddy == 0 && ddz == 0) continue;
            const cell_idx_t jx = ix + cell_idx_c(ddx), jy = iy + cell_idx_c(ddy), jz = iz + cell_idx_c(ddz);
            // neighbour must be a fluid cell (not another obstacle cell of the same/other body)
            if (particleField->get(jx, jy, jz) != simObjs->particleAccessor->getInvalidUid()) continue;

            if (flagField->isFlagSet(jx, jy, jz, gasFlag)) // dry link (gas-side; would be skipped by the guard)
            {
               ++nDry[fi];
               sumDry[fi * 3 + 0] += real_c(ddx); sumDry[fi * 3 + 1] += real_c(ddy); sumDry[fi * 3 + 2] += real_c(ddz);
            }
            else // wet link (liquid/interface — carries the force)
            {
               ++nWet[fi];
               sumWet[fi * 3 + 0] += real_c(ddx); sumWet[fi * 3 + 1] += real_c(ddy); sumWet[fi * 3 + 2] += real_c(ddz);
            }
         }
      }
   }

   // --- position / velocity of each floater (owner contributes; ghosts skipped) ---
   for (uint_t i = uint_t(0); i < simObjs->particleAccessor->size(); ++i)
   {
      if (isSet(simObjs->particleAccessor->getFlags(i), mesa_pd::data::particle_flags::GHOST)) continue;
      const walberla::id_t uid = simObjs->particleAccessor->getUid(i);
      for (uint_t k = 0; k < nF; ++k)
         if (simObjs->floaterUIDs[k] == uid)
         {
            const auto p = simObjs->particleAccessor->getPosition(i);
            const auto v = simObjs->particleAccessor->getLinearVelocity(i);
            for (uint_t c = 0; c < 3; ++c) { pos[k * 3 + c] = p[c]; vel[k * 3 + c] = v[c]; }
            break;
         }
   }

   mpi::reduceInplace< uint_t >(nWet, mpi::SUM);
   mpi::reduceInplace< uint_t >(nDry, mpi::SUM);
   mpi::reduceInplace< real_t >(sumWet, mpi::SUM);
   mpi::reduceInplace< real_t >(sumDry, mpi::SUM);
   mpi::reduceInplace< real_t >(pos, mpi::SUM);
   mpi::reduceInplace< real_t >(vel, mpi::SUM);

   WALBERLA_ROOT_SECTION()
   {
      // Per-floater lines only for floaters showing a bias (one line per floater per call would flood the log);
      // a one-line summary is always printed.
      constexpr real_t tiltWarnThreshold = real_t(0.5);
      real_t maxTilt = real_t(0);
      uint_t nFlagged = uint_t(0);
      for (uint_t k = 0; k < nF; ++k)
      {
         const real_t dryH = std::sqrt(sumDry[k * 3] * sumDry[k * 3] + sumDry[k * 3 + 1] * sumDry[k * 3 + 1]);
         const real_t dryZ = std::abs(sumDry[k * 3 + 2]);
         // The cap has thousands of links, so the absolute horizontal residual is meaningless; the
         // tilt of the dry cap is |dryXY| relative to its vertical extent |dryZ|. A symmetric cap
         // has ratio ~0; a strongly one-sided cap (lateral wetted-force bias) climbs toward 1.
         const real_t tilt = (dryZ > real_t(0)) ? dryH / dryZ : real_t(0);
         maxTilt = std::max(maxTilt, tilt);
         if (tilt <= tiltWarnThreshold) continue;
         ++nFlagged;
         WALBERLA_LOG_WARNING("[LINK-BALANCE t=" << t << " after=" << location << "] floater " << k
            << " uid=" << simObjs->floaterUIDs[k]
            << " | pos=(" << pos[k * 3] << "," << pos[k * 3 + 1] << "," << pos[k * 3 + 2] << ")"
            << " velXY=(" << vel[k * 3] << "," << vel[k * 3 + 1] << ")"
            << " | wet=" << nWet[k] << " dry=" << nDry[k]
            << " | dryVec=(" << sumDry[k * 3] << "," << sumDry[k * 3 + 1] << "," << sumDry[k * 3 + 2] << ")"
            << " tilt(|dryXY|/|dryZ|)=" << tilt << "   <-- STRONGLY ONE-SIDED DRY CAP (lateral force bias)");
      }
      WALBERLA_LOG_INFO_ON_ROOT("[LINK-BALANCE t=" << t << " after=" << location << "] " << nF
         << " floaters checked, max tilt=" << maxTilt << ", " << nFlagged << " above " << tiltWarnThreshold);
   }
}

} // namespace antidunes
} // namespace walberla