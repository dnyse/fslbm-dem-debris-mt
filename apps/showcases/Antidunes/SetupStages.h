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
#include <iomanip>
#include <sstream>
#include <vector>

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

   // NOTE: obstacleNormalField is populated by ObstacleNormalSweep every timestep (always run, not gated on enableWetting).
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
         // TODO: these checks assume DEM dynamics are active; with useNoSlipParticles=true,
         // fixingHeightFactor is set large to fix all particles, making particleFixingHeight >> maxParticleHeight
         // and causing a false abort. Skip them when useNoSlipParticles=true.
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

   for (const auto& pillar : simParams->pillars)
   {
      auto p = createBox(*simObjs->particleStorage, pillar.center, real_t(2) * pillar.halfExtents);
      p->setType(2);
   }

   for (const auto& pillar : simParams->cylinderPillars)
   {
      createCylinder(*simObjs->particleStorage, pillar.center, pillar.radius);
   }

   for (const auto& fsp : simParams->floatingSpheres)
   {
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
         // Hold the sphere fixed during the settling phase so the FSLBM can establish a stable
         // fill-level distribution around the contact line before the sphere starts moving.
         if (fsp.settlingSteps > 0)
            mesa_pd::data::particle_flags::set(p->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      }
      WALBERLA_LOG_INFO_ON_ROOT("Created floating sphere: radius=" << fsp.radius
         << " densityRatio=" << fsp.densityRatio << " position=" << pos
         << " settlingSteps=" << fsp.settlingSteps);
   }

   // choose particle sync algorithm -----------------------------------------------------
   // blockSyncExtension must be >= the largest floating sphere radius so that SyncNextNeighbors
   // ghost-copies the sphere to every block whose interior it physically overlaps. Without this,
   // cells inside the sphere on a neighboring block are not mapped as MovingObstacle, which creates
   // LBM PDF instabilities (NaN torque / quaternion corruption) for spheres larger than ~2.5 LU.
   real_t blockSyncExtension = real_c(2.5);
   for (const auto& fsp : simParams->floatingSpheres)
      blockSyncExtension = std::max(blockSyncExtension, fsp.radius);
   WALBERLA_LOG_INFO_ON_ROOT("blockSyncExtension = " << blockSyncExtension);
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

   // initialize BCs from config (NoSlip top/bottom, CellInterval pillar defined in Antidunes.prm)
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

// Called every timestep: releases any floating sphere whose settlingSteps have elapsed.
// During the settling phase the sphere is held FIXED so the FSLBM can establish a stable
// fill-level / normal field around the contact line before the sphere starts moving.
void releaseSettledFloatingSpheres(const std::shared_ptr< const SimulationParams >& simParams,
                                   const std::shared_ptr< SimulationObjects >&       simObjs,
                                   const uint_t currentStep)
{
   for (const auto& fsp : simParams->floatingSpheres)
   {
      if (fsp.settlingSteps == uint_t(0) || currentStep != fsp.settlingSteps) continue;

      simObjs->particleStorage->forEachParticle(
         false, mesa_pd::kernel::SelectAll(), *simObjs->particleAccessor,
         [](const size_t idx, auto& ac) {
            if (ac.getType(idx) == uint_t(1)) // type 1 = floating sphere
               mesa_pd::data::particle_flags::unset(ac.getFlagsRef(idx),
                                                    mesa_pd::data::particle_flags::FIXED);
         },
         *simObjs->particleAccessor);

      WALBERLA_LOG_INFO_ON_ROOT("Released floating sphere(s) from settling phase at step " << currentStep);
   }
}

// Debug helper: logs cell-type counts and fill-level / normal statistics in the contact-circle
// band around the first floating sphere (type 1).  The band is defined by:
//   |z_world - sphereZ| <= 2  AND  |sqrt((x-sx)²+(y-sy)²) - radius| <= 2
//
// Runs every 50 steps normally; every step from step 490 onward (to catch the run-up to the
// step-581 crash in the surface-start variant).  No-op when there are no floating spheres.
//
// Call from the top of the main time loop, AFTER releaseSettledFloatingSpheres().
void debugContactLine(const std::shared_ptr< const SimulationParams >& simParams,
                      const std::shared_ptr< SimulationObjects >&       simObjs,
                      const uint_t t)
{
   if (simParams->floatingSpheres.empty()) return;

   // Sparse log: every 50 steps up to step 490, then every step.
   const bool denseLog = (t >= uint_t(490));
   if (!denseLog && (t % uint_t(50) != uint_t(0))) return;

   // ---- Gather floating sphere (type 1) position from whichever rank owns it ----
   Vector3< real_t > spherePos(real_t(0));
   real_t            sphereRadius = real_t(0);

   simObjs->particleStorage->forEachParticle(
      false, mesa_pd::kernel::SelectLocal(), *simObjs->particleAccessor,
      [&spherePos, &sphereRadius](const size_t idx, auto& ac) {
         if (ac.getType(idx) == uint_t(1))
         {
            spherePos    = ac.getPosition(idx);
            sphereRadius = static_cast< const mesa_pd::data::Sphere* >(ac.getBaseShape(idx).get())->getRadius();
         }
      },
      *simObjs->particleAccessor);

   // allReduceInplace(SUM): non-owning ranks contribute 0; owning rank contributes the real value.
   mpi::allReduceInplace< real_t >(spherePos[0],    mpi::SUM);
   mpi::allReduceInplace< real_t >(spherePos[1],    mpi::SUM);
   mpi::allReduceInplace< real_t >(spherePos[2],    mpi::SUM);
   mpi::allReduceInplace< real_t >(sphereRadius,    mpi::SUM);

   if (sphereRadius < real_t(1e-6)) return; // sphere not found

   // ---- Count cell types and fill/normal stats in the contact-circle band ----
   uint_t nLiquid = 0, nInterface = 0, nGas = 0, nObs = 0;
   real_t fillSum  = real_t(0);
   real_t fillMin  = real_t(1e10);
   real_t fillMax  = real_t(-1e10);
   real_t normalMagMin  = real_t(1e10);
   uint_t nZeroNormal   = 0;
   // World coordinates and sphere-relative coords of the locally minimum fill cell
   real_t minFillWx = real_t(0), minFillWy = real_t(0), minFillWz = real_t(0);
   real_t minFillDz = real_t(0), minFillDr = real_t(0);

   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* fillField   = blockIt->getData< ScalarField_T >(simObjs->fillFieldID);
      const auto* flagField   = blockIt->getData< FlagField_T   >(simObjs->flagFieldID);
      const auto* normalField = blockIt->getData< VectorField_T >(simObjs->normalFieldID);

      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto gasFlag       = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      // obstacle flag is a mask (may be multiple bits for static+moving obstacle types)
      const auto obstacleMask  = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const auto& aabb         = blockIt->getAABB();

      for (cell_idx_t iz = 0; iz < cell_idx_c(flagField->zSize()); ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(flagField->ySize()); ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(flagField->xSize()); ++ix)
      {
         const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
         const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
         const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);

         const real_t dz  = wz - spherePos[2];
         const real_t rXY = std::sqrt((wx - spherePos[0]) * (wx - spherePos[0]) +
                                      (wy - spherePos[1]) * (wy - spherePos[1]));
         const real_t dr  = rXY - sphereRadius;

         if (std::abs(dz) > real_t(2) || std::abs(dr) > real_t(2)) continue;

         // Check obstacle FIRST: sphere surface cells may have both obstacle AND interface/gas
         // flags set simultaneously (initFlagsFromFillLevel runs after mapParticlesOntoCells
         // and can add fluid flags to cells that are already obstacle-mapped).
         // Use bitwise AND to check if ANY obstacle flag bit is set.
         if (flagField->get(ix, iy, iz) & obstacleMask)
         {
            ++nObs;
         }
         else if (flagField->isFlagSet(ix, iy, iz, interfaceFlag))
         {
            ++nInterface;
            const real_t f  = fillField->get(ix, iy, iz);
            if (f < fillMin) { fillMin = f; minFillWx = wx; minFillWy = wy; minFillWz = wz; minFillDz = dz; minFillDr = dr; }
            fillMax  = std::max(fillMax, f);
            fillSum += f;
            const real_t nm = normalField->get(ix, iy, iz).length();
            normalMagMin = std::min(normalMagMin, nm);
            if (nm < real_t(1e-6)) ++nZeroNormal;
         }
         else if (flagField->isFlagSet(ix, iy, iz, liquidFlag))    ++nLiquid;
         else if (flagField->isFlagSet(ix, iy, iz, gasFlag))        ++nGas;
         // else: unconverted/unknown flag — not counted
      }
   }

   // ---- MPI reduce across all blocks / ranks ----
   mpi::allReduceInplace< uint_t >(nLiquid,      mpi::SUM);
   mpi::allReduceInplace< uint_t >(nInterface,   mpi::SUM);
   mpi::allReduceInplace< uint_t >(nGas,         mpi::SUM);
   mpi::allReduceInplace< uint_t >(nObs,         mpi::SUM);
   mpi::allReduceInplace< real_t >(fillSum,      mpi::SUM);
   // Save local fillMin before MIN-reduction so we can identify which rank holds the global min.
   const real_t fillMinLocal = fillMin;
   mpi::allReduceInplace< real_t >(fillMin,      mpi::MIN);
   mpi::allReduceInplace< real_t >(fillMax,      mpi::MAX);
   mpi::allReduceInplace< real_t >(normalMagMin, mpi::MIN);
   mpi::allReduceInplace< uint_t >(nZeroNormal,  mpi::SUM);

   // Broadcast position of the global-minimum fill cell: the rank that achieved the local min
   // equal to the global min contributes its (wx,wy,wz,dz,dr); others contribute 0.
   // Correct as long as the exact minimum is achieved by only one rank (overwhelmingly the case).
   const bool hasGlobalMin = (std::abs(fillMinLocal - fillMin) < real_t(1e-12));
   real_t minFillCellWx = hasGlobalMin ? minFillWx : real_t(0);
   real_t minFillCellWy = hasGlobalMin ? minFillWy : real_t(0);
   real_t minFillCellWz = hasGlobalMin ? minFillWz : real_t(0);
   real_t minFillCellDz = hasGlobalMin ? minFillDz : real_t(0);
   real_t minFillCellDr = hasGlobalMin ? minFillDr : real_t(0);
   mpi::allReduceInplace< real_t >(minFillCellWx, mpi::SUM);
   mpi::allReduceInplace< real_t >(minFillCellWy, mpi::SUM);
   mpi::allReduceInplace< real_t >(minFillCellWz, mpi::SUM);
   mpi::allReduceInplace< real_t >(minFillCellDz, mpi::SUM);
   mpi::allReduceInplace< real_t >(minFillCellDr, mpi::SUM);

   const real_t fillAvg = (nInterface > uint_t(0)) ? fillSum / real_c(nInterface) : real_t(0);

   WALBERLA_LOG_INFO_ON_ROOT("[DBG " << t
      << "] sphere=(" << spherePos[0] << "," << spherePos[1] << "," << spherePos[2]
      << ") r=" << sphereRadius
      << " | contact band: liq=" << nLiquid << " int=" << nInterface
      << " gas=" << nGas << " obs=" << nObs
      << " | fill[min=" << fillMin << " avg=" << fillAvg << " max=" << fillMax << "]"
      << " | normal_min=" << normalMagMin << " zero_normals=" << nZeroNormal);
   // When fill_min crosses -0.008 log the cell's world position and sphere-relative coords.
   if (fillMin < real_t(-0.008))
   {
      WALBERLA_LOG_INFO_ON_ROOT("[DBG " << t << "] fill_min_cell world=("
         << minFillCellWx << "," << minFillCellWy << "," << minFillCellWz
         << ") sphere-rel dz=" << minFillCellDz << " dr=" << minFillCellDr);
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

      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);

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
         const bool isInterface = flagField->isFlagSet(ix, iy, iz, interfaceFlag);
         const bool isLiquid    = flagField->isFlagSet(ix, iy, iz, liquidFlag);
         const char* flagStr = isInterface ? "INTERFACE" : (isLiquid ? "LIQUID" : "OTHER");

         // build PDF list string
         std::ostringstream pdfStr;
         pdfStr << std::scientific << std::setprecision(4);
         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
         {
            pdfStr << "f[" << f << "]=" << pdfField->get(ix, iy, iz, f);
            if (f + 1 < PdfField_T::Stencil::Size) pdfStr << " ";
         }

         WALBERLA_LOG_INFO("[WATCH t=" << t << " " << label
            << "] world=(" << wc[0] << "," << wc[1] << "," << wc[2] << ")"
            << " local=(" << ix << "," << iy << "," << iz << ")"
            << " flag=" << flagStr
            << " rho=" << std::scientific << std::setprecision(6) << rho
            << " fill=" << fill
            << " | " << pdfStr.str());
      }
   }
}

// Scan every liquid and interface cell on every local block for NaN/Inf PDFs.
// Logs the world coordinates, block-local coordinates, and PDF index of the first NaN found
// on each rank, then returns true if any NaN was detected globally.
// Call this after each major sweep to pinpoint which one first corrupts the PDF field.
bool checkPdfNaN(const std::shared_ptr< SimulationObjects >& simObjs,
                 const char* const location, const uint_t t)
{
   bool localNaN = false;

   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* pdfField  = blockIt->getData< const PdfField_T >(simObjs->pdfFieldID);
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto& aabb = blockIt->getAABB();

      bool blockDone = false;
      for (cell_idx_t iz = 0; iz < cell_idx_c(pdfField->zSize()) && !blockDone; ++iz)
      for (cell_idx_t iy = 0; iy < cell_idx_c(pdfField->ySize()) && !blockDone; ++iy)
      for (cell_idx_t ix = 0; ix < cell_idx_c(pdfField->xSize()) && !blockDone; ++ix)
      {
         // only check fluid-carrying cells
         if (!flagField->isFlagSet(ix, iy, iz, liquidFlag) &&
             !flagField->isFlagSet(ix, iy, iz, interfaceFlag)) continue;

         for (uint_t f = 0; f < PdfField_T::Stencil::Size; ++f)
         {
            const real_t val = pdfField->get(ix, iy, iz, f);
            if (!std::isfinite(val))
            {
               const real_t wx = aabb.xMin() + real_c(ix) + real_t(0.5);
               const real_t wy = aabb.yMin() + real_c(iy) + real_t(0.5);
               const real_t wz = aabb.zMin() + real_c(iz) + real_t(0.5);
               WALBERLA_LOG_WARNING("[NaN-SCAN t=" << t << " after=" << location
                  << "] block-local=(" << ix << "," << iy << "," << iz << ")"
                  << " world=(" << wx << "," << wy << "," << wz << ")"
                  << " f[" << f << "]=" << val
                  << " flag=" << (flagField->isFlagSet(ix,iy,iz,interfaceFlag) ? "interface" : "liquid"));
               localNaN = true;
               blockDone = true;
               break;
            }
         }
      }
   }

   // reduce across all ranks so every rank knows whether any NaN exists globally
   uint_t globalNaN = localNaN ? uint_t(1) : uint_t(0);
   mpi::allReduceInplace< uint_t >(globalNaN, mpi::SUM);
   return globalNaN > uint_t(0);
}

} // namespace antidunes
} // namespace walberla