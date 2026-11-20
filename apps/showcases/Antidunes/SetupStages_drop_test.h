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
//! \file SetupStages_drop_test.h
//! Copy of SetupStages.h for the fluid-free DropTest executable only (see DropTest.cpp). The two files differ
//! in several places, e.g. this copy adds X-wall planes, container cylinders and capsule floaters, and lacks the
//! accelerationXOverride handling, the checkpoint-resume guard around floater creation, the contact-skin
//! inflation of interaction radii and the repeated initial particle sync; diff against SetupStages.h for the
//! full list. Kept as a separate copy so the executables that include SetupStages.h are not affected.
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
      // In floater-only runs the storage also holds other shapes, for which static_cast<Sphere*> is
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

   // create X-wall planes when X is non-periodic (DropTest only: it uses periodicity <0,0,0> to enclose the
   // box on all six sides; SetupStages.h keeps X periodic as the flow direction and has no such block)
   if (!simParams->blockForest.periodicity[0])
   {
      createPlane(*simObjs->particleStorage, bf->getDomain().minCorner(), Vector3< real_t >( 1, 0, 0));
      createPlane(*simObjs->particleStorage, bf->getDomain().maxCorner(), Vector3< real_t >(-1, 0, 0));
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

   // vertical round containers (ContainerCylinders prm block, see ContainerCylinderParams)
   for (const auto& container : simParams->containerCylinders)
   {
      createContainerCylinder(*simObjs->particleStorage, container.center, container.radius);
      WALBERLA_LOG_INFO_ON_ROOT("Created container cylinder (containment, axis z): center=" << container.center
                                << " radius=" << container.radius);
   }

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
         p->setInteractionRadius(halfDiag);
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
      const real_t circumR    = std::sqrt(fcp.radius * fcp.radius + halfLength * halfLength);

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

   // Floating capsules ("pill" shape, type 5), same indexed-UID pattern and body-z axis alignment as cylinders.
   // interactionRadius = length/2 exactly: `length` is the full end-to-end length including the hemispherical
   // caps, and the farthest surface point from the centre is a cap tip on the axis. There is no closed-form
   // contact for capsules; every capsule pair falls through to GJK/EPA in GeneralContactDetection (see
   // SolidCylinderContactDetection in ParticleKernels.h).
   simObjs->floaterCapsuleUIDs.assign(simParams->floatingCapsules.size(), walberla::id_t(0));
   for (size_t pi = 0; pi < simParams->floatingCapsules.size(); ++pi)
   {
      const auto& fpp = simParams->floatingCapsules[pi];
      const mesa_pd::Vec3 pos(fpp.position[0], fpp.position[1], fpp.position[2]);
      const real_t circumR = real_t(0.5) * fpp.length;

      mesa_pd::Rot3 rot;
      const mesa_pd::Vec3 ez(real_t(0), real_t(0), real_t(1));
      const real_t cosA = ez * fpp.axis;
      if (cosA < real_t(1) - real_t(1e-10))
      {
         const mesa_pd::Vec3 rotAxis = (cosA > real_t(-1) + real_t(1e-10))
                                          ? (ez % fpp.axis).getNormalized()
                                          : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
         rot.rotate(rotAxis, std::acos(std::clamp(cosA, real_t(-1), real_t(1))));
      }
      if (std::abs(fpp.initialTiltAngleDeg) > real_t(1e-10))
         rot.rotate(mesa_pd::Vec3(real_t(0), real_t(1), real_t(0)), fpp.initialTiltAngleDeg * math::pi / real_t(180));

      if (simObjs->mesapdDomain->isContainedInProcessSubdomain(
            uint_c(mpi::MPIManager::instance()->rank()), pos))
      {
         auto p = simObjs->particleStorage->create(false);
         p->setPosition(pos);
         p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(fpp.radius, fpp.length));
         p->getBaseShapeRef()->updateMassAndInertia(fpp.densityRatio);
         p->setOwner(mpi::MPIManager::instance()->rank());
         p->setType(5);
         p->setInteractionRadius(circumR);
         p->setRotation(rot);
         if (fpp.settlingSteps > 0 && !fpp.pinPosition)
            mesa_pd::data::particle_flags::set(p->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
         simObjs->floaterCapsuleUIDs[pi] = p->getUid();
      }
      WALBERLA_LOG_INFO_ON_ROOT("Created floating capsule: radius=" << fpp.radius
         << " length=" << fpp.length << " axis=" << fpp.axis
         << " densityRatio=" << fpp.densityRatio << " position=" << pos
         << " settlingSteps=" << fpp.settlingSteps
         << " initialTiltAngleDeg=" << fpp.initialTiltAngleDeg
         << " pinPosition=" << (fpp.pinPosition ? "true" : "false"));
   }
   mpi::allReduceInplace(simObjs->floaterCapsuleUIDs, mpi::SUM);
   for (size_t pi = 0; pi < simParams->floatingCapsules.size(); ++pi)
      if (simParams->floatingCapsules[pi].startInactive)
         simObjs->dormantFloaterUIDs->insert(simObjs->floaterCapsuleUIDs[pi]);

   // choose particle sync algorithm -----------------------------------------------------
   // blockSyncExtension must be >= the largest floater extent so that ghost copies reach every block a floater
   // overlaps: sphere radius, box half-diagonal, cylinder circumscribed radius, capsule half-length.
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
                                       real_t(0.25) * fcp.length * fcp.length);
      blockSyncExtension = std::max(blockSyncExtension, circumR);
   }
   for (const auto& fpp : simParams->floatingCapsules)
      blockSyncExtension = std::max(blockSyncExtension, real_t(0.5) * fpp.length);
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
                                       real_t(0.25) * fcp.length * fcp.length);
      maxParticleDiameter = std::max(maxParticleDiameter, circumR * real_c(2));
   }
   for (const auto& fpp : simParams->floatingCapsules)
      maxParticleDiameter = std::max(maxParticleDiameter, fpp.length);

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
   // Select floaters: type 1 = floating sphere, 3 = floating box, 4 = floating cylinder, 5 = floating capsule.
   // Fixed obstacles (pillars, walls) are mapped separately as static NoSlip/FreeSlip. Dormant floaters
   // (startInactive, before settlingSteps) are excluded so they do not claim fluid cells; this lets several
   // floaters share one spawn position. Must stay consistent with movingObstacleSelector in BuildTimeloop.h.
   auto dormantUIDs = simObjs->dormantFloaterUIDs;
   auto floaterSelector = [dormantUIDs](const size_t particleIdx, const auto& ac)
   {
      const uint_t t = ac.getType(particleIdx);
      if (t != uint_t(1) && t != uint_t(3) && t != uint_t(4) && t != uint_t(5)) return false;
      return dormantUIDs->count(ac.getUid(particleIdx)) == 0;
   };

   // map particles onto the LBM domain (initialisation only; the time loop does the per-step mapping)
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
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
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
      }
   }

   // reduce across all ranks so every rank knows whether any NaN exists globally
   uint_t globalNaN = localNaN ? uint_t(1) : uint_t(0);
   mpi::allReduceInplace< uint_t >(globalNaN, mpi::SUM);
   return globalNaN > uint_t(0);
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

} // namespace antidunes
} // namespace walberla