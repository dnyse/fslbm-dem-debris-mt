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

// Creates the deck-fill and haunch solid FixedBox pieces for one arch tunnel (both no-ops when the
// corresponding ArchSolidBox geometry comes back std::nullopt -- see ArchTunnelParams's own comment,
// SimulationParams.h, for the exact backward-compatible defaulting rule). Kept separate from the
// geometry math itself (archDeckFillBoxGeometry/archHaunchBoxGeometry, SimulationParams.h) so that
// math is unit-testable without a live mesa_pd::data::ParticleStorage.
void createArchSolidElements(mesa_pd::data::ParticleStorage& ps, const ArchTunnelParams& tunnel)
{
   if (const auto deckFill = archDeckFillBoxGeometry(tunnel))
      createBox(ps, deckFill->center, deckFill->edgeLength);

   if (const auto leftHaunch = archHaunchBoxGeometry(tunnel, ArchHaunchSide::Left))
      createBox(ps, leftHaunch->center, leftHaunch->edgeLength);

   if (const auto rightHaunch = archHaunchBoxGeometry(tunnel, ArchHaunchSide::Right))
      createBox(ps, rightHaunch->center, rightHaunch->edgeLength);
}

// The complete arch block as ONE call: the curved containment wall (createArchTunnel) plus every
// solid flat piece around it (createArchSolidElements) -- everything §07 of the "Arch Tunnel
// Containment" writeup describes as reading like one continuous shape, built in one place instead
// of two calls a caller has to remember to pair up. Equivalent to calling both separately; kept as
// two underlying particles (the curve is not a box, and mesa_pd has no CSG boolean shape), just
// created together.
void createArchBlock(mesa_pd::data::ParticleStorage& ps, const ArchTunnelParams& tunnel)
{
   createArchTunnel(ps, tunnel.center, tunnel.radius, tunnel.xMin, tunnel.xMax);
   createArchSolidElements(ps, tunnel);
}

// Dimension-based builder: the arch's own footprint as X/Y/Z extents plus a bore radius, rather than
// requiring the caller to pre-compute xMin/xMax/yMin/yMax/deckTopZ by hand -- see archBlockParams'
// own comment (SimulationParams.h) for the dims->params conversion and what `base`/xDim/yDim/zDim
// mean. A prm-driven ArchTunnelParams stays the normal path for the real bridge scenarios
// (SimulationParams.h parses it directly from the ArchTunnels block); this overload is for callers
// building geometry in code with no prm block to parse from -- e.g. DEMValidationScenarios.h, or a
// future hand-authored bridge layout script.
void createArchBlock(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& base,
                      real_t xDim, real_t yDim, real_t zDim, real_t radius)
{
   WALBERLA_CHECK_GREATER_EQUAL(yDim, real_t(2) * radius,
      "createArchBlock: yDim (" << yDim << ") must be at least 2*radius (" << real_t(2) * radius
      << ") for the circular bore to fit inside the block's own width.");
   WALBERLA_CHECK_GREATER_EQUAL(zDim, radius,
      "createArchBlock: zDim (" << zDim << ") must be at least radius (" << radius
      << ") for the block to reach the crown at all.");

   createArchBlock(ps, archBlockParams(base, xDim, yDim, zDim, radius));
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
      // The static_cast<Sphere*> + avg-diameter consistency check below only makes sense when the
      // checkpointed particleStorage actually holds a settled BED of spheres (useParticleBed=true)
      // to validate against spheres_out.dat's header-derived avgParticleDiameter -- it was never
      // exercised any other way before. In a floater-only scenario (useParticleBed=false, e.g. the
      // mixed-shape sphere/box/cylinder showcases), particleStorage holds whatever shapes the
      // floaters are; blindly static_cast-ing a Box/Cylinder's BaseShape to Sphere* and calling
      // getRadius() on it is undefined behaviour (reads garbage as a radius), which corrupts
      // avgParticleDiameterTest and made the very first attempt to loadSnapshot a mixed-shape
      // checkpoint abort with "Particle diameters not correct." even though nothing was actually
      // wrong. Guarding on useParticleBed keeps the original bed-consistency check fully intact
      // for bed-based runs while skipping it (just counting particles, any shape, no cast) here.
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

   // Pillars are GLOBAL/FIXED (Utility.h's createBox/createCylinder) and, per BlockForestDataHandling's
   // own `if (GLOBAL) continue;` guard in both serialize() and deserialize(), are explicitly EXCLUDED
   // from checkpoint save/load -- unlike floaters (see the loadSnapshot guard around the floater
   // creation loops below), they are correctly recreated fresh from this same prm-declared config on
   // EVERY run, load or not, since their position never changes and re-deriving it from the (byte-
   // identical, same-prm) config is exactly equivalent to restoring it. Logged here (previously silent,
   // unlike floaters' own "Created floating ..." lines) so a checkpoint-resume run's log directly shows
   // the reconstructed position for the same kind of verification floaters already get.
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

   // ModelParams.useArchTunnelAsMesaPDObject (see its own doc comment): when enabled, each arch
   // tunnel becomes a real, infinite-mass mesa_pd body instead of the one-sided ComputeTunnelContacts
   // kernel (ParticleKernels.h's ParticleKernels constructor skips constructing that kernel whenever
   // this flag is set, so the two paths are never both active for the same tunnel).
   if (simParams->model.useArchTunnelAsMesaPDObject)
   {
      for (const auto& tunnel : simParams->archTunnels)
      {
         createArchBlock(*simObjs->particleStorage, tunnel);
      }
   }

   // Floater creation (spheres/boxes/cylinders below) is SKIPPED entirely when resuming from a
   // checkpoint. Unlike the pillars/planes/arch-tunnels above (all GLOBAL, hence explicitly excluded
   // from checkpoint serialization -- see BlockForestDataHandling.cpp's own `if (GLOBAL) continue;` in
   // both serialize() and deserialize() -- so recreating them fresh every run, load or not, is correct
   // by design), floaters are ordinary non-GLOBAL particles and ARE fully captured in
   // particleStorageFile. Re-running this creation loop on top of an already-loaded checkpoint
   // therefore tries to insert a SECOND particle with the same UID as one already restored, aborting
   // with "particle with the same uid(N) already existing at index(M)" -- confirmed directly (a local
   // round-trip test hit this exact abort before this guard was added).
   //
   // KNOWN LIMITATION: since floaterUIDs/floaterBoxUIDs/floaterCylinderUIDs and dormantFloaterUIDs are
   // only ever populated by the creation loops below, skipping them here means those stay empty on a
   // loaded run -- releaseSettledFloatingSpheres() (called every step) then has nothing to act on. This
   // is harmless when the checkpoint was taken AFTER every floater's own settlingSteps has already
   // passed (the common case for resuming near a late-run crash: all floaters are already active by
   // then, nothing left to release) but means a checkpoint taken WHILE some floaters are still
   // startInactive/dormant will leave them parked (FIXED) forever after this kind of resume -- their
   // dormant state IS correctly restored from the checkpoint (the FIXED flag is ordinary particle
   // state, not excluded like GLOBAL), they just never get released without the bookkeeping arrays.
   // Not fixed here: doing so robustly needs a way to re-identify which restored particle corresponds
   // to which prm-declared floater index without relying on position (particles move) -- e.g. a
   // persisted per-particle index/settlingSteps property -- a larger change than this specific
   // investigation needs.
   if (bfParams.loadSnapshot)
   {
      // Sized (to the prm-declared counts) but left as the invalid UID (0) -- see this block's own
      // opening comment for why: correct indexed access for any downstream code that iterates these
      // arrays by the same prm-declared size, just with no valid UID to act on (matching
      // dormantFloaterUIDs being left empty too).
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
   // Indexed loop so each sphere's UID is captured at the correct index fi.
   // The owning rank (the one whose subdomain contains pos) creates the particle and records its UID.
   // allReduceInplace(SUM) broadcasts each UID to all other ranks (no double-counting because only
   // one rank owns each position).  ParticleKernels uses simObjs->floaterUIDs directly instead of
   // trying to infer the mapping from local particle ordering, which breaks for multiple spheres
   // distributed across different ranks.
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
         // Hold the sphere fixed during the settling phase so the FSLBM can establish a stable
         // fill-level distribution around the contact line before the sphere starts moving.
         // pinPosition holds only the TRANSLATIONAL DOFs (done each DEM sub-cycle in ParticleKernels);
         // it must NOT set the FIXED flag too, same mutual-exclusion rule as boxes/cylinders.
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
   // Register startInactive spheres as dormant -- AFTER the UID sync above, so every rank inserts the
   // same, now-fully-resolved UID (this loop itself needs no communication: floatingSpheres[fi] is
   // parsed identically from the same prm on every rank).
   for (size_t fi = 0; fi < simParams->floatingSpheres.size(); ++fi)
      if (simParams->floatingSpheres[fi].startInactive)
         simObjs->dormantFloaterUIDs->insert(simObjs->floaterUIDs[fi]);

   // Floating boxes (type 3) — same indexed-UID pattern as floating spheres.
   // interactionRadius = circumscribed-sphere half-diagonal so getParticleCellBB produces a
   // conservative but correct cell bounding box for the shape-generic MEM cell test.
   // FloatingBox is freely rotating; the base Box::updateMassAndInertia sets the correct
   // slice-based V_sub estimate in computeSubmergedVolumesBox() stays valid.
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
            // Rotate about Y-axis (body axis v in Bogner's diagram)
            p->setRotation(mesa_pd::Rot3(mesa_pd::Vec3(0, 1, 0), angle));
         }
         // pinPosition holds only the TRANSLATIONAL DOFs (done each DEM sub-cycle in ParticleKernels);
         // it must NOT set the FIXED flag, which would also freeze rotation. settlingSteps still uses
         // FIXED (full freeze) but is mutually exclusive with pinPosition.
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

   // Floating cylinders / logs (type 4) — same indexed-UID pattern as spheres and boxes.
   // The mesa_pd Cylinder's axis is its BODY z-axis; the initial rotation aligns body-z with the
   // requested world `axis` (default <0,1,0> = spanwise log), then applies the optional extra tilt
   // about the world Y-axis. interactionRadius = circumscribed-sphere radius sqrt(r² + (L/2)²).
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
         // [Claude, 2026-09-26] EXPERIMENTAL toggle (ModelParameters.useCapsuleShapeForCylinders,
         // default false): a Capsule (cylinder + hemispherical end caps) has the same (radius, length)
         // constructor signature and the same body-z axis convention as Cylinder, so nothing else in
         // this loop (position, rotation, interactionRadius, FIXED flag, settlingSteps) needs to change.
         if (simParams->model.useCapsuleShapeForCylinders)
            p->setBaseShape(std::make_shared< mesa_pd::data::Capsule >(fcp.radius, fcp.length));
         else
            p->setBaseShape(std::make_shared< mesa_pd::data::Cylinder >(fcp.radius, fcp.length));
         p->getBaseShapeRef()->updateMassAndInertia(fcp.densityRatio);
         p->setOwner(mpi::MPIManager::instance()->rank());
         p->setType(4);
         p->setInteractionRadius(circumR);
         p->setRotation(rot);
         // Same settling semantics as boxes: FIXED (full freeze) unless pinPosition holds the
         // translational DOFs instead (rotation free, applied each DEM sub-cycle in ParticleKernels).
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
   } // end of "if (!bfParams.loadSnapshot)" -- see its own opening comment above

   // choose particle sync algorithm -----------------------------------------------------
   // blockSyncExtension must be >= the largest floater extent so that SyncNextNeighbors
   // ghost-copies each floater to every block its body physically overlaps.
   // For spheres: extent = radius. For boxes: extent = circumscribed half-diagonal.
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
   WALBERLA_LOG_INFO_ON_ROOT("blockSyncExtension = " << blockSyncExtension);
   // maxParticleDiameter determines which sync algorithm is safe (SyncNextNeighbors vs
   // SyncGhostOwners). Use 2×circumscribed radius for each floater type so that particles
   // large enough to span 2+ blocks always get SyncGhostOwners.
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

   // ROOT CAUSE FIX (THESIS_PLAN.md item 1a, checkpoint-resume STALE-PARTICLE-IDX abort): mesa_pd's
   // SyncGhostOwners::checkAndResolveOverlap() only sends a new ghost copy to ranks in
   // domain.getNeighborProcesses() -- i.e. ranks topologically ONE HOP away in the block graph from
   // wherever a particle is currently known (owner or existing ghost), not every rank whose subdomain
   // actually falls within its interactionRadius+dx. A large floater (e.g. a box with
   // interactionRadius comparable to several blocks' width) can physically overlap ranks MULTIPLE
   // hops away from its owning rank. During a continuous run this is never a problem past the first
   // step: ghost state accumulates and is incrementally re-verified every single timestep
   // (doDEMSubCycles' own sync calls), so it has long since converged to the particle's full
   // multi-hop footprint by any given later step. On a checkpoint resume, `ClearNextNeighborSync`/
   // `ClearGhostOwnerSync` (see the loadSnapshot branch above) wipes all ghost state, and this is the
   // ONLY sync call before mapParticlesOntoCells() and the very first FSLBM sweep run -- a single call
   // only reaches one hop out from each particle's loaded position. Confirmed directly: job 4236626 (a
   // resumed cuboid-equilibrium run) had a Y-elongated box with interactionRadius~45 LU, cellsPerBlock
   // <17,42,30> -- diagonally-adjacent (2-hop) ranks reported zero knowledge of the box after exactly
   // one syncCall(), producing SimpleBB.h's own "[STALE-PARTICLE-IDX] ... uid no longer maps to a
   // valid local particle index" hard abort on the very first post-load step.
   //
   // [Claude, 2026-09-18]: the multi-round call below was originally gated to loadSnapshot only, on
   // the assumption that "a fresh start never needs more than one hop yet (freshly-created floaters/
   // beds are local to a handful of blocks)". That assumption is false for a floater whose own
   // multi-hop footprint (governed by blockSyncExtension, not its raw geometric size) already spans
   // more than one hop at creation -- confirmed directly via job 4245028/4245089 (this scenario's own
   // FloatingBox, halfExtents<15,30,15>, cellsPerBlock<16,45,30>): the [GHOST-SYNC-DIAG] diagnostic
   // showed the owning rank at block index (3,1,2) with all 26 of its Moore-neighborhood neighbors
   // correctly ghost-synced after the single call above, but ranks 13/17/21 at block indices
   // (1,0,1)/(1,1,1)/(1,2,1) -- Chebyshev distance 2 from the owner -- showed zero knowledge of the
   // box, and each owns a real fluid cell genuinely adjacent (confirmed against the box's own
   // 5-degree-tilted, hence enlarged, AABB) to the box's true surface. Same abort, same mechanism,
   // just a fresh start instead of a resume -- so the gate is now removed. Calling syncCall()
   // repeatedly lets each round's newly-created ghosts become the "seed" for the next round's
   // neighbor-of-neighbor propagation, converging to the true multi-hop footprint after enough
   // rounds -- 4 is a generous margin over the 2 hops needed in both confirmed cases above, and cheap
   // relative to the rest of setup (a one-time cost, not a per-step one) regardless of which path.
   for (uint_t round = uint_t(0); round < uint_t(4); ++round) simObjs->syncCall();

   // create particle field (mapping of particle UIDs to lbm cells)
   simObjs->particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>(
      bf, "Particle field", simObjs->particleAccessor->getInvalidUid(), field::fzyx, uint_c(2));
}

void setupBoundaryHandling(const std::shared_ptr< const SimulationParams >& simParams, const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("setupBoundaryHandling");

   // Live MEM reference height (ModelParameters.useDynamicReferenceHeight): seed to the flat still-water
   // level; ParticleKernels' computeSubmergedVolumes() EMA-updates it every step from the floaters' fitted
   // local waterlines when the toggle is on. Force-side only (SimpleBB's rhoHydStat never touches the
   // fluid), so this is safe regardless of the toggle — when off, the pointer is simply never written to
   // and the static-height functor overload below is used instead.
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

   // useBognerThreePhaseMEM: pass a null functor so SimpleBB's `hydrostaticDensityFct_ == nullptr`
   // branch is taken -- no static/dynamic H subtraction at all. Combined with the (already-default)
   // fill-weighting below, this is exactly Bogner & Rude's Eq. 5.9: a phi(x_b)-weighted blend of the
   // real boundary-link PDFs and the gas-equilibrium distribution, nothing else.
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

   // Auto NoSlip walls on non-periodic LATERAL faces (X = E/W, Y = N/S), so the walls always stay
   // consistent with BlockForestParameters.periodicity — no manual prm sync. A non-periodic face with
   // NO boundary condition has undefined PDFs → NaN (the exact y-boundary NaN seen when periodicity was
   // <1,0,0> but the N/S walls had been removed). Periodic faces are left wrapped (no wall). The prm
   // still declares top/bottom (Z, never periodic); re-adding an existing wall here is idempotent.
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

   // ROOT CAUSE FIX (THESIS_PLAN.md item 3, checkpoint-resume NaN crash): both ParticleMappingKernel and
   // MovingParticleMappingKernel below only ever claim a cell that boundaryHandling->isDomain() already
   // reports true (see MovingParticleMapping.h's own "no else -> will not overwrite pre-existing
   // boundary cells" comment) -- but isDomain() checks the domain_ flag mask, which nothing has set yet
   // at this point: setupBoundaryHandling() only ran initFromConfig() (Border NoSlip), and
   // initFlagsFromFillLevel() -- the call that actually classifies the rest of the field as
   // liquid/interface/gas (all part of domain_) -- runs AFTER this function, in Antidunes.cpp. Every
   // OTHER caller of these mapping kernels in the codebase (tests/lbm_mesapd_coupling/mapping/
   // ParticleMapping.cpp, MovingParticleMappingMEM.cpp, SettlingSphereMEM.cpp, etc.) calls
   // handling->fillWithDomain() before mapping for exactly this reason -- Antidunes never did. The
   // result: this one-time setup mapping has ALWAYS silently mapped ZERO cells for ANY floater, on
   // every run, fresh or resumed. Harmless on a fresh start (t=0's fluid field is freshly initialized
   // sane data everywhere, and the very next per-step mapping in timeloopAfterDEM -- which runs AFTER
   // initFlagsFromFillLevel() has already executed once during setup -- silently fixes the one-step
   // gap before it matters). Fatal on any checkpoint resume with active (non-dormant) floaters present:
   // fillWithDomain() only touches genuinely empty cells (see its own isEmpty() guard), so this is safe
   // to add unconditionally; it does not overwrite the Border-NoSlip flags initFromConfig() already set.
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      auto* handling = blockIt->getData< AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         simObjs->antidunesBoundaryHandling->getHandlingID());
      handling->fillWithDomain(handling->getFlagField()->nrOfGhostLayers());
   }

   // Select all floaters: type 1 = floating sphere, type 3 = floating box, type 4 = floating cylinder.
   // Fixed obstacles (pillars, walls) are mapped separately as static NoSlip/FreeSlip. Dormant
   // (startInactive, not yet at settlingSteps) floaters are excluded too -- they are not "in the
   // simulation" yet, so must not claim fluid cells (this is what lets several floaters share the
   // exact same spawn position without conflicting). Must stay consistent with movingObstacleSelector
   // in BuildTimeloop.h (the per-step equivalent of this init-time mapping).
   auto dormantUIDs = simObjs->dormantFloaterUIDs;
   auto floaterSelector = [dormantUIDs](const size_t particleIdx, const auto& ac)
   {
      const uint_t t = ac.getType(particleIdx);
      if (t != uint_t(1) && t != uint_t(3) && t != uint_t(4)) return false;
      return dormantUIDs->count(ac.getUid(particleIdx)) == 0;
   };

   // map particles onto LBM domain (this is for initilization. iterations happen in timeloopAfterParticles)
   if (simParams->particle.useNoSlipParticles) {
      // sets noSlip flag (in flagField through boundaryhandling) for all cells inside the particle
      lbm_mesapd_coupling::ParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
         kernel(simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID());
      simObjs->particleStorage->forEachParticle(false, floaterSelector, *simObjs->particleAccessor,
                                             kernel, *simObjs->particleAccessor,
                                             AntidunesBoundaryHandling_T::noSlipFlagID);
   }
   else {
      // same as above but it also stores the UID of the particle in particleFieldID (but no hydro forces are computed here)
      lbm_mesapd_coupling::MovingParticleMappingKernel< AntidunesBoundaryHandling_T::BoundaryHandling_T >
         kernel(simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID(), simObjs->particleFieldID,
                simParams->model.mappingTieOffset);
      simObjs->particleStorage->forEachParticle(false, floaterSelector, *simObjs->particleAccessor,
                                             kernel, *simObjs->particleAccessor,
                                             AntidunesBoundaryHandling_T::movingObstacleFlagID);
   }
}

// Called every timestep: releases any floating sphere or box whose settlingSteps have elapsed.
// During the settling phase the floater is held FIXED so the FSLBM can establish a stable
// fill-level / normal field around the contact line before the particle starts moving.
//
// Each release also inserts the UID into simObjs->freshlyActivatedFloaterUIDs, a debug/validation
// hook (see that member's own doc comment, SimulationObjects.h) for a "flying floater" bug traced
// to floaters being released directly on top of an already-active neighbour: rather than a
// separate, necessarily-approximate geometric check here, this lets the floater fall into the
// REAL contact detection it was excluded from while dormant (ComputeContacts, ParticleKernels.h,
// same GJK/EPA/analytical dispatch used for actual physics) on the very next doDEMSubCycles() call
// -- exact, not approximated, and with zero extra geometry code. ComputeContacts hard-aborts, with
// full diagnostic detail, if it finds a genuine contact involving a UID in that set; the set is
// cleared unconditionally at the end of that same doDEMSubCycles() call either way.
void releaseSettledFloatingSpheres(const std::shared_ptr< const SimulationParams >& simParams,
                                   const std::shared_ptr< SimulationObjects >&       simObjs,
                                   const uint_t currentStep)
{
   // Release each floater INDIVIDUALLY by UID. The old code matched on type only, so the first
   // sphere reaching its settlingSteps unset FIXED on EVERY type-1 sphere — staggered release
   // schedules (rows dropped at 30k/40k/50k) all fired at once. floaterUIDs[i] ↔ floatingSpheres[i].
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
      // [Claude, 2026-09-25] force-density field, printed alongside PDFs below: the generated
      // Cumulant collision kernel reads this DIRECTLY (xi_29/30/31 = xi_3*_data_force[...], folded
      // into velocity and every higher-order cumulant) as an input entirely separate from the PDF
      // field -- structurally invisible to every PDF-based NaN/precursor check in this investigation.
      // THESIS_PLAN.md hopper crash investigation: raw PDFs at the actual crash cell were confirmed
      // clean/ordinary right up to collide, so this field is the one remaining unverified input.
      const auto* forceField = blockIt->getData< const VectorFieldFlattened_T >(simObjs->forceDensityFieldID);
      // [Claude, 2026-09-27] particleField_ owning UID, added per explicit request ("what about
      // obstacle neighbors... contested ownership between the two nearby rods" discussion): lets a
      // cube trace show whether a cell's obstacle-ownership stamp is stable or flip-flopping between
      // two nearby floaters across steps, independent of the flag classification above (a cell can be
      // OBST-flagged by one uid one step and a DIFFERENT uid the next without the flag itself ever
      // changing category).
      const auto* particleFieldPtr =
         blockIt->getData< const lbm_mesapd_coupling::ParticleField_T >(simObjs->particleFieldID);
      // [Claude, 2026-09-27] Curvature/normal fields, added per explicit request to test the
      // "free-surface curvature/normal reconstruction degeneracy" hypothesis (mechanism #3) against
      // Bogner's dissertation (2017, Section 4.4.3.1): the FD curvature method this codebase uses
      // (CurvatureSweepFiniteDifferences, BuildTimeloop.h -- Bogner's least-stable of the three
      // methods he compared) can enter a self-sustaining, non-converging oscillation when spurious
      // currents from curvature error perturb the interface layer, which changes the curvature
      // computation again, "thus explaining the oscillations" (Bogner, p.99). Logging these lets a
      // cube trace show whether curvature/normal values are already degenerate/unstable at this
      // heavily-obstacle-surrounded cell before the density collapse.
      const auto* curvatureField = blockIt->getData< const ScalarField_T >(simObjs->curvatureFieldID);
      const auto* normalField    = blockIt->getData< const VectorField_T >(simObjs->normalFieldID);
      const auto* obstacleNormalField =
         blockIt->getData< const VectorField_T >(simObjs->obstacleNormalFieldID);

      const auto liquidFlag    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto interfaceFlag = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      // [Claude, 2026-09-24] richer classification (was INTERFACE/LIQUID/OTHER, lumping GAS and
      // OBSTACLE together) -- matches checkPdfNaN's own classify() so a cube trace can distinguish
      // a genuinely-gas cell from an obstacle-covered one.
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

// [Claude, 2026-09-27] Per explicit request ("create a box around all the cylinders in the crash
// area, log out their values and see the crash evolution"): logs EVERY OBSTACLE-or-INTERFACE cell
// inside a fixed world-space box, once per MAIN step (the fluid field only changes once per main
// step -- see CLAUDE.md's per-timestep execution order; sub-cycle resolution would multiply row
// count 10x for zero additional fluid-state information), to a per-rank CSV. Complements
// trackWatchCells' handful of hardcoded single cells with the FULL spatial evolution of a region
// known to contain a recurring crash's rods -- restricted to OBST/INTF cells (not every cell in the
// box) to keep row count bounded (a 2-3-rod-sized box is a few thousand cells, not tens of
// thousands). Gated by ModelParameters.logFluidBoxHistory (+ Min/Max/StartStep/EndStep,
// SimulationParams.h); each rank opens (and truncates) its own CSV on first call via a function-
// local static, mirroring FloaterHistoryLogger's per-rank-file convention (ParticleKernels.h) but
// without needing its own class since this is called from exactly one place (BuildTimeloop.h's
// addAfterTimeStepDiagnosticSweep).
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
   // Flush every call (once per main step, not per cell), same crash-survivability rationale as
   // FloaterHistoryLogger -- this run can WALBERLA_ABORT/MPI_ABORT at any point.
   file.flush();
}

// [Claude, 2026-09-28] Per explicit request to investigate the sphere/cylinder/cuboid still-water
// equilibrium position bias (settles consistently HIGHER -- less submerged -- than the naive
// constant-density Archimedes target, growing with density ratio): logs EVERY cell (all flags, not
// just OBST/INTF -- unlike logFluidBoxHistory above) along the full domain height at one fixed
// world (x, y) column, once per MAIN step, to a per-rank CSV. Point of this: get the ACTUAL
// far-field density profile rho(z), to check whether the settled fluid is hydrostatically
// stratified under the gravity body force (weakly-compressible LBM: p = cs^2*rho, so rho should
// decrease with height) rather than the uniform rho=1 the simple Archimedes target assumes -- and,
// independently, to get a precise fill-level-based free-surface height instead of guessing at a
// "waterline". Gated by ModelParameters.logDensityColumn (+ X/Y/StartStep/EndStep, SimulationParams.h).
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

// [Claude, 2026-09-28] ModelParameters.logWaterHeight: fill-integrated free-surface height on a regular X x Y grid of
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
         // Only liquid/interface cells are fatal: a NaN there is genuine corruption. NaN PDFs in gas
         // cells at the obstacle waterline are expected (and harmless — SimpleBB skips those dry links);
         // they are NOT scanned here. (The confirmation that they exist was done; see git history.)
         if (!flagField->isFlagSet(ix, iy, iz, liquidFlag) &&
             !flagField->isFlagSet(ix, iy, iz, interfaceFlag)) continue;

         // Low-rho PRECURSOR check (PhysicsParams.lowRhoPrecursorThreshold, 0 = off): the per-PDF
         // isfinite scan below only catches a cell that has ALREADY gone non-finite -- by definition
         // one step too late to see how it got there. This catches the step BEFORE that: a cell whose
         // aggregate density is still finite but has drained to a suspiciously small value (candidate
         // for the "drained interface cell, rho->0" state SimpleBB.h's own dry-link comment describes
         // as expected/benign at a floater's waterline -- the open question this exists to answer is
         // whether a cell that goes on to actually blow up is distinguishable, in the steps
         // beforehand, from the many that presumably drain harmlessly). Throttled to the first 20 hits
         // per rank so a run where this fires often (many ordinary draining interface cells) doesn't
         // flood the log -- this is meant for a targeted, windowed rerun (denseNanScanEnabled +
         // From/ToStep bracketing an already-known crash), not routine production use.
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

         // High-velocity PRECURSOR check (PhysicsParams.highVelocityPrecursorThreshold, 0 = off): the
         // symmetric counterpart to the low-rho check above -- a blowup does not always show up as a
         // cell draining toward rho=0 first; it can equally show up as an unrealistically large but
         // still-finite velocity for a step or two (the classic LBM divergence precursor: local Mach
         // number pushed far past the scheme's stability bound) before the PDFs actually overflow to
         // Inf/NaN. Uses PdfField_T::getVelocity(), which for this project's compressible cumulant
         // model divides by the cell's own actual rho (not an assumed 1) -- see PdfField.h's own
         // getVelocity() implementation. Same throttling (first 20 hits per rank) and same intended
         // usage (a targeted, windowed rerun) as the low-rho check.
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

         // Population PRECURSOR check (PhysicsParams.populationPrecursorThreshold, 0 = off): see that parameter's doc comment.
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
                  // recent-history / dynamics context of the neighbourhood: cells a rod has just LEFT (former moving obstacle), interface/gas neighbours, fastest neighbouring fluid cell, owners of the obstacle neighbours
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

         // [Claude, 2026-09-24] VELOCITY-derived NaN check, distinct from the raw-PDF isfinite scan
         // just above: getVelocity() = momentum/density. If a cell's individual PDFs are ALL
         // technically finite (e.g. all exactly 0) but happen to sum to exactly zero density, the
         // division produces 0/0 = NaN velocity WITHOUT any single f_i ever being non-finite -- the
         // scan above cannot see this case at all. Found via mean-velocity chain-of-custody analysis
         // (THESIS_PLAN.md, step-444 crash investigation): the ONLY consumer that ever sees the
         // failure is Forcing.h's MeanVelocityComputer, which sums getVelocity() over every
         // liquid/interface cell -- exactly the quantity checked here.
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

// [Claude, 2026-09-25] Narrow, reactive safety net for the hopper-clogging collide-stage NaN
// (THESIS_PLAN.md): the confirmed root cause is that an obstacle-flagged cell's PDF array can end up
// reading as genuinely all-zero -- streamReconstructAdvect never writes PDFs for obstacle cells at all
// (only carries the fill level forward), and the internal FieldClone scratch buffer that sweep ping-
// pongs with the real field is created via cloneUninitialized() (never explicitly initialised, unlike
// the real field's own equilibrium fill at setup) -- so a component that is never independently
// refreshed by SimpleBB's per-link bounce-back (i.e. one facing another obstacle cell, not real fluid)
// can be sitting on genuinely untouched memory, which reads back as exactly 0.0. collide then computes
// m_000=0 -> xi_0=1/m_000=inf -> every downstream cumulant collapses to NaN. Confirmed directly via a
// raw in-kernel diagnostic reading collide's own _data_pdfs[...] pointer arithmetic (job 4266339).
//
// Unlike checkPdfNaN above (which only scans liquid/interface cells and aborts), this scans OBSTACLE
// cells specifically and REPAIRS a degenerate one in place -- resetting it to the D3Q27 rest-equilibrium
// distribution (rho=1, u=0; the same weights every ordinary near-equilibrium cell in this codebase's own
// logs already shows: w0=8/27, w[1-6]=2/27, w[7-18]=1/54, w[19-26]=1/216) -- since nothing downstream
// ever reads an obstacle cell's own stored PDF value except collide itself: SimpleBB::treatDirection
// only reads the FLUID-side cell's PDF (never the obstacle cell's own, confirmed from source);
// PdfReconstructionManager's seedConvertedCellFill derives fill from neighbour averaging, not this
// cell's frozen history; MovingParticleMapping's transition code never reads or writes pdfField at all.
// Only touches a cell that is ALREADY non-finite or non-positive -- the vast majority of obstacle
// cells (including ones actively receiving legitimate SimpleBB bounce-back writes every step) are never
// modified. Returns the number of cells repaired on this rank (for logging by the caller).
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

// FLUID SANITY SWEEP (PhysicsParams::fluidSanity*): per-step range check of every LIQUID/INTERFACE cell (gas, obstacle, former-moving-obstacle and boundary cells are never checked; fluidSanityLiquidOnly also drops the interface
// cells), see the parameters' doc comment (SimulationParams.h). Returns the GLOBAL number of violating cells.
// Log throttling: the first 30 violating steps are always printed, afterwards only steps whose count at least doubles vs. the last printed one (or every 100th step), so a benign steady-state violation cannot fill the budget
// and hide the escalation.
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
            {   // largest jump to a checked (liquid/interface) FACE neighbour inside this block (block-boundary faces are skipped: the ghost-layer PDFs are one half-step older than the interior)
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
         {   // a LIQUID cell must not touch a gas cell (the interface layer is one cell thick over the full D3Q27 stencil)
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
   // global extrema time series (independent of violations): lets the growth of max|u| / rho excursions be read off BEFORE a threshold is crossed
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

// TEMPORARY DIAGNOSTIC (THESIS_PLAN.md item 1a, R3 equilibrium fresh-run NaN crash) -- checkPdfNaN
// stops at the FIRST bad cell per block (blockDone=true; break), which is fine for its own
// abort-and-report purpose but hides how many bad cells actually exist per block. Same scan, no
// early exit, listing every LIQ/INTF cell with a non-finite PDF (world coordinate + rho only, to
// keep output compact across many ranks). Call this right before aborting on a checkPdfNaN==true
// result to see the true extent of the corruption instead of one artifact-of-iteration-order cell
// per block. Remove once the determination is made.
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
// END TEMPORARY DIAGNOSTIC

// ModelParameters.highAltitudeWaterScanEnabled/highAltitudeWaterThresholdZ: scans every LOCAL cell
// whose world-Z exceeds thresholdZ and warns immediately if it is liquid or interface -- see
// ModelParams' own doc comment (SimulationParams.h) for why this exists (neither waLBerla's
// CellConversionSweep nor this project's own wrapping of it logs anything at all when a gas cell
// converts, so a "how did water get up there" crash has no log trail to follow otherwise). Local
// cells only (matches checkPdfNaN's own convention) so a ghost-layer ownership boundary never
// double-counts one real conversion as two warnings. Does NOT abort -- purely observational, meant
// to be read back after the run (or grepped live) to find the exact step/cell a conversion first
// happened at, at a height where it needn't have.
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

      // Skip this block entirely if its own Z-extent never reaches the threshold -- cheap early-out,
      // avoids the per-cell loop for every block below the watched altitude on every step.
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

// Verify the "closed interface layer" invariant: in FSLBM a liquid cell may NEVER be a (D3Q27)
// neighbour of a gas cell — they must always be separated by an interface cell. A violation is a
// "hole" in the interface layer (the Fig. 4.12a configuration). This scans every liquid cell for a
// gas neighbour, counts the holes, and logs the first one's location + its 3x3x3 neighbourhood.
// Does NOT abort — it is meant to reveal whether/when the invariant breaks relative to the NaN.
// Returns the global number of holes found.
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

// Measure the azimuthal asymmetry of each floater's wetted vs. dry surface — the diagnostic for the
// lateral drift introduced by the SimpleBB dry-link guard.
//
// Rationale: the MEM bounce-back integrates the hydrodynamic force only over WET links (fluid-side
// cell carries a valid PDF). The guard SKIPS dry links (gas / NaN fluid-side cell at the waterline)
// — both the force contribution AND the bounce-back reflection. For a floating sphere the dry set is
// the cap above the waterline; if it is azimuthally SYMMETRIC its horizontal components cancel and
// the net hydrodynamic force is purely vertical. If the dry set is LOPSIDED (because waterline cells
// collapse to gas asymmetrically — e.g. driven by the streamwise flow under liquid promotion), the
// wetted-surface force integral acquires a net HORIZONTAL component → the sphere drifts sideways.
// The all-interface reconstruction never starves cells to gas at the contact line, so its dry set
// stays symmetric and no drift appears — exactly the difference the user observed.
//
// For each floater we scan its obstacle cells (particleField == floater UID) and, over the D3Q27
// fluid neighbours, accumulate the outward link offsets split into WET (liquid/interface) and DRY
// (gas) sets. The horizontal magnitude of the DRY-set vector is the lateral-bias signal. Reduced to
// root and logged; does not abort.
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
      for (uint_t k = 0; k < nF; ++k)
      {
         const real_t dryH = std::sqrt(sumDry[k * 3] * sumDry[k * 3] + sumDry[k * 3 + 1] * sumDry[k * 3 + 1]);
         const real_t dryZ = std::abs(sumDry[k * 3 + 2]);
         // The cap has thousands of links, so the absolute horizontal residual is meaningless; the
         // tilt of the dry cap is |dryXY| relative to its vertical extent |dryZ|. A symmetric cap
         // has ratio ~0; a strongly one-sided cap (lateral wetted-force bias) climbs toward 1.
         const real_t tilt = (dryZ > real_t(0)) ? dryH / dryZ : real_t(0);
         WALBERLA_LOG_INFO("[LINK-BALANCE t=" << t << " after=" << location << "] floater " << k
            << " uid=" << simObjs->floaterUIDs[k]
            << " | pos=(" << pos[k * 3] << "," << pos[k * 3 + 1] << "," << pos[k * 3 + 2] << ")"
            << " velXY=(" << vel[k * 3] << "," << vel[k * 3 + 1] << ")"
            << " | wet=" << nWet[k] << " dry=" << nDry[k]
            << " | dryVec=(" << sumDry[k * 3] << "," << sumDry[k * 3 + 1] << "," << sumDry[k * 3 + 2] << ")"
            << " tilt(|dryXY|/|dryZ|)=" << tilt
            << ((tilt > real_t(0.5)) ? "   <-- STRONGLY ONE-SIDED DRY CAP (lateral force bias)" : ""));
      }
   }
}

} // namespace antidunes
} // namespace walberla