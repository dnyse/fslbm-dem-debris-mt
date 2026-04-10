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

#include "lbm/field/AddToStorage.h"

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

   const real_t accelerationX = real_c(3) * simParams->pid.targetMeanVelocityMagnitude * simParams->computed.viscosity / 
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

   // Initialize meanVelocity and totalFluidMass as zero 
   // They will be later computed during evalutation stages
   simObjs->meanVelocity = std::make_shared< Vector3< real_t > >(real_c(0));
   simObjs->totalFluidMass = std::make_shared< real_t >(real_c(0));
   
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

   // NOTE: obstacleNormalField is actually never computed, but its not wrong when enablewetting is false. (otherwise it has to be computed before curvature sweep)
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
      simObjs->particleStorageID = bf->addBlockData(mesa_pd::domain::createBlockForestDataHandling(simObjs->particleStorage), "Particle Storage");
      
      // Only create particles if the are below the initialization profile -----------------------------------------------
      auto createParticleFct = [&simParams](Vector3< real_t > pos) {
         return pos[2] < initializationProfile(pos[0], simParams->computed.sinusAmplitude, simParams->computed.sinusOffset, simParams->computed.sinusWavelength);
      };

      real_t maxParticleHeight = real_c(0); // top of the sphere (i.e. maxPos + its radius)

      WALBERLA_LOG_INFO_ON_ROOT("Initializing particles from bedgeneration output file : " << simParams->particle.inFileName);
      initSpheresFromFile(simParams->particle.inFileName, *simObjs->particleStorage, *simObjs->mesapdDomain,
                          simParams->particle.densityRatio, simParams->domain.domainSize, createParticleFct,
                          bf->getDomain(), simParams->particle.bedCopiesInX, simParams->particle.bedCopiesInY,
                          numParticles, maxParticleHeight, simParams->computed.particleScalingFactor);

      WALBERLA_LOG_INFO_ON_ROOT("Max particle height " << maxParticleHeight);
      if ((simParams->computed.sinusOffset + simParams->computed.sinusAmplitude) > maxParticleHeight)
         WALBERLA_ABORT("Created particle bed is below desired sinus shape!");
      if (real_c(2) * simParams->computed.sinusAmplitude > (maxParticleHeight - simParams->computed.particleFixingHeight))
         WALBERLA_ABORT("Created particle bed is not high enough for desired sinus shape!");
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

   walberla::mpi::allReduceInplace(numFixedParticles, walberla::mpi::SUM);
   simObjs->numParticles = numParticles;

   WALBERLA_LOG_INFO_ON_ROOT("Total number of particles in the system: " << simObjs->numParticles << " (of which " << numFixedParticles << " are fixed)");
   

   // create top and bottom planes --------------------------------------------------------
   createPlane(*simObjs->particleStorage, bf->getDomain().minCorner(), Vector3< real_t >(0, 0,  1));
   createPlane(*simObjs->particleStorage, bf->getDomain().maxCorner(), Vector3< real_t >(0, 0, -1));

   // choose particle sync algorithm -----------------------------------------------------
   const real_t blockSyncExtension    = real_c(2.5);
   const real_t maxParticleDiameter   = simParams->computed.avgParticleDiameter * real_c(1.1);  // 1.1 for lubrication correction

   if (maxParticleDiameter < real_c(2) * real_c(simParams->blockForest.cellsPerBlock.min()) - blockSyncExtension) {
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

   // create particle field (mapping of particle UIDs to lbm cells)
   simObjs->particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>(
      bf, "Particle field", simObjs->particleAccessor->getInvalidUid(), field::fzyx, uint_c(2));
}

void setupBoundaryHandling(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("setupBoundaryHandling");
   simObjs->antidunesBoundaryHandling = std::make_shared< AntidunesBoundaryHandling_T >(
      simObjs->blockForest, simObjs->pdfFieldID, simObjs->fillFieldID, simObjs->particleFieldID,
      simObjs->particleAccessor, HydrostaticDensityFtor(simObjs->acceleration, simParams->computed.absoluteLiquidHeight));

   // initialize BCs (noSlip top/bottom) from config file
   simObjs->antidunesBoundaryHandling->initFromConfig(simParams->cfgPtr->getOneBlock("BoundaryParameters"));
   
   // store the flagField and flagInfo handles in simObjs (only for convenience)
   simObjs->flagFieldID = simObjs->antidunesBoundaryHandling->getFlagFieldID();
   simObjs->flagInfo    = simObjs->antidunesBoundaryHandling->getFlagInfo();
}

void mapParticlesOntoCells(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("mapParticlesOntoCells");
   auto sphereSelector = [](const size_t particleIdx, const auto& ac)
   {
      return ac.getBaseShape(particleIdx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE;
   };
   
   // map particles onto LBM domain (this is for initilization. iterations happen in timeloopAfterParticles)
   if (simParams->particle.useNoSlipParticles) {
      // sets noSlip flag (in flagField through boundaryhandling) for all cells inside the particle 
      lbm_mesapd_coupling::ParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
         kernel(simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID());
      simObjs->particleStorage->forEachParticle(false, sphereSelector, *simObjs->particleAccessor,
                                             kernel, *simObjs->particleAccessor,
                                             AntidunesBoundaryHandling_T::noSlipFlagID);
   } 
   else {
      // same as above but it also stores the UID of the particle in particleFieldID (but no hydro forces are computed here)
      lbm_mesapd_coupling::MovingParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
         kernel(simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID(), simObjs->particleFieldID);
      simObjs->particleStorage->forEachParticle(false, sphereSelector, *simObjs->particleAccessor,
                                             kernel, *simObjs->particleAccessor,
                                             AntidunesBoundaryHandling_T::movingObstacleFlagID);
   }
}

} // namespace antidunes
} // namespace walberla