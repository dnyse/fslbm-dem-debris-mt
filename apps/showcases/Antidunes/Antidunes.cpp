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
//! \file Antidunes.cpp
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Jonas Plewinski <jonas.plewinski@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// This showcase simulates antidunes, i.e., particulate dunes that travel in opposite stream-wise direction. See
// simulation results published in https://doi.org/10.1017/jfm.2023.262.
//======================================================================================================================

#include "core/Environment.h"
#include "core/waLBerlaBuildInfo.h"

#include "sqlite/SQLite.h"
#include "timeloop/all.h"

#include "BuildTimeloop.h"
#include "WriteData.h"
#include "ParticleKernels.h"
#include "SimulationParams.h"
#include "SimulationObjects.h"
#include "SetupStages.h"

namespace walberla {
namespace antidunes {

int main(int argc, char** argv)
{
   Environment walberlaEnv(argc, argv);
   WALBERLA_LOG_INFO_ON_ROOT("waLBerla Revision: " << std::string(WALBERLA_GIT_SHA1).substr(0, 8));
   if (argc < 2) { WALBERLA_ABORT("Please specify a parameter file as input argument.") }
   auto configPtr = walberlaEnv.config();
   WALBERLA_LOG_INFO_ON_ROOT(*configPtr);

   // Simulation params and objects
   // ----------------------------------------------------------------------------------------------------------------------------------
   const auto simParams = std::make_shared< SimulationParams >(configPtr);
   const auto simObjs   = std::make_shared< SimulationObjects >();

   // Create directories for storing snapshots and evaluations
   // ----------------------------------------------------------------------------------------------------------------------------------
   WALBERLA_ROOT_SECTION()
   {
      std::filesystem::create_directories(simParams->blockForest.snapshotBaseFolder);
      std::filesystem::create_directories(simParams->evaluation.baseFolderName);
      if (simParams->model.logFloaterHistory)
         std::filesystem::create_directories(simParams->model.floaterHistoryBaseFolder);
      if (simParams->model.logFluidBoxHistory)
         std::filesystem::create_directories(simParams->model.fluidBoxHistoryBaseFolder);
      if (simParams->model.logDensityColumn)
         std::filesystem::create_directories(simParams->model.densityColumnBaseFolder);
      if (simParams->model.logWaterHeight)
         std::filesystem::create_directories(simParams->model.waterHeightBaseFolder);
   }

   // Log computed params
   // ----------------------------------------------------------------------------------------------------------------------------------
   simParams->computed.logOnRoot();
   simParams->computed.writeToFile(simParams->evaluation.baseFolderName);

   // Setup simulation (initialize or load from snapshot)
   // ----------------------------------------------------------------------------------------------------------------------------------
   setupBlockForest           (simParams, simObjs);
   setupSimState              (simParams, simObjs);   // beginTimestep + acceleration + pid
   setupVTKState              (simParams, simObjs);   // VTK counter
   setupLBMFields             (simParams, simObjs);   // pdf + fill + excessMass
   setupSurfaceGeometryFields (simParams, simObjs);   // curvature + normal + obstacleNormal + smoothFill
   setupParticles             (simParams, simObjs);   // particleStorage + particleField

   // Boundary handling and flags
   // ----------------------------------------------------------------------------------------------------------------------------------
   setupBoundaryHandling      (simParams, simObjs);
   mapParticlesOntoCells      (simParams, simObjs);   // set solid cells in flagField (+ map UIDs to particleField)

   // Set liquid/interface/gas cells in flagField (must be called after all solid cells are set)
   // ----------------------------------------------------------------------------------------------------------------------------------
   simObjs->antidunesBoundaryHandling->initFlagsFromFillLevel();
   checkOuterGhostLayerFlags(simObjs, "after initFlagsFromFillLevel"); // see BuildTimeloop.h
   if (simParams->model.wallOuterGhostLayer) // see SimulationParams.h
   {
      uint_t n = wallOuterGhostLayer(simObjs);
      mpi::reduceInplace(n, mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("[WALL-OUTER-GHOST-LAYER] " << n << " outer ghost-layer cells behind walls set to wall");
      checkOuterGhostLayerFlags(simObjs, "after wallOuterGhostLayer");
   }

   // Init forceDensity (requires valid fill + pdf + flag fields and acceleration)
   // ----------------------------------------------------------------------------------------------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT("Initializing force density field");
   walberla::free_surface::initForceDensityFieldCodegen< PdfField_T, FlagField_T, VectorFieldFlattened_T, ScalarField_T >(
      simObjs->blockForest, simObjs->forceDensityFieldID, simObjs->fillFieldID, simObjs->pdfFieldID, simObjs->flagFieldID,
      simObjs->flagInfo, *simObjs->acceleration);

   // Communicate after initialization
   // ----------------------------------------------------------------------------------------------------------------------------------
   // particleFieldID is included here to match BuildTimeloop.h's own per-step communication of it
   // alongside flagFieldID -- omitting it leaves a rank that only sees a floater as a GHOST with a
   // correctly-communicated "obstacle" flag but a stale/uncommunicated particleField entry, until
   // the first per-step coupling sweep runs. On a checkpoint resume this can hit the SimpleBB.h MEM
   // force sweep's stale-particle-idx guard on the very first post-load step.
   WALBERLA_LOG_INFO_ON_ROOT("Initial communication");
   Communication_T(simObjs->blockForest, simObjs->flagFieldID, simObjs->fillFieldID, simObjs->forceDensityFieldID,
                    simObjs->particleFieldID)();
   PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID)();

   // Create bubble model
   simObjs->bubbleModel = std::make_shared< walberla::free_surface::bubble_model::BubbleModelConstantPressure >(real_c(1));

   // Create timeloops
   // ----------------------------------------------------------------------------------------------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT("Creating time loops");
   const auto timeloopFSLBM     = std::make_unique< SweepTimeloop >(simObjs->blockForest, simParams->physics.timesteps);
   const auto timeloopAfterDEM  = std::make_unique< SweepTimeloop >(simObjs->blockForest, simParams->physics.timesteps);
   const auto timeloopEval      = std::make_unique< SweepTimeloop >(simObjs->blockForest, simParams->physics.timesteps);

   timeloopFSLBM->setCurrentTimeStep(simObjs->beginTimeStep);
   timeloopAfterDEM->setCurrentTimeStep(simObjs->beginTimeStep);
   timeloopEval->setCurrentTimeStep(simObjs->beginTimeStep);

   // Before timestep -----------------------------------------------------
   addBeforeTimeStepVTKOutput          (simParams, simObjs, timeloopFSLBM);

   // FSLBM ---------------------------------------------------------------
   addSweepSurfaceGeometry             (simParams, simObjs, timeloopFSLBM);
   addSweepSurfaceDynamics             (simParams, simObjs, timeloopFSLBM);

   // LBM-DEM coupling ----------------------------------------------------
   addSweepFluidParticleCoupling       (simParams, simObjs, timeloopAfterDEM);

   // DEM -----------------------------------------------------------------
   ParticleKernels< ParticleAccessor_T > particleKernels(simParams, simObjs);

   // Evaluate ------------------------------------------------------------
   timeloopEval->add() << Sweep(EmptySweep(), "Empty sweep: (evaluators/loggers)");
   addAfterTimeStepStabilityChecker    (simParams, simObjs, timeloopEval);
   addAfterTimeStepDataEvaluators      (simParams, simObjs, timeloopEval);
   addAfterTimeStepPerformanceLogger   (simParams, simObjs, timeloopEval);
   addAfterTimeStepForcingAdjuster     (simParams, simObjs, timeloopEval);
   addAfterTimeStepDiagnosticSweep     (simParams, simObjs, timeloopEval);
   addDensityColumnLoggingSweep        (simParams, simObjs, timeloopEval);
   addWaterHeightLoggingSweep          (simParams, simObjs, timeloopEval);

   WcTimingPool timingPool;

   // Watch cells for per-step PDF/density/flag tracing at a fragile free-surface contact-line
   // configuration (checkPdfNaN dumps the first NaN's neighbourhood and aborts; trackWatchCells
   // logs the watched cells' own history). Sampled every WATCH_FREQ steps; WATCH_DENSE_FROM enables
   // per-step tracking from a given step on, disabled by default (past any realistic step count).
   std::vector< Vector3< real_t > > watchCells;
   for (real_t dz : { real_t(-1), real_t(0), real_t(1) })
      for (real_t dy : { real_t(-1), real_t(0), real_t(1) })
         for (real_t dx : { real_t(-1), real_t(0), real_t(1) })
            watchCells.emplace_back(real_t(20.5) + dx, real_t(8.5) + dy, real_t(25.5) + dz);
   constexpr uint_t WATCH_FREQ      = 1000;   // periodic sampling, cheap enough to run routinely
   constexpr uint_t WATCH_DENSE_FROM = 1000000000;

   // Probe region for the FSLBM/DEM coupling stages around a specific, already-resolved
   // hopper-clogging crash cell. Disabled by default (FROM > TO), same convention as
   // WATCH_DENSE_FROM above; re-enable by setting a real step range if that investigation needs
   // revisiting.
   std::vector< Vector3< real_t > > hopperCrashCells;
   for (real_t y = real_t(5.5); y <= real_t(35.5); y += real_t(1))
      for (real_t z : { real_t(17.5), real_t(18.5), real_t(19.5), real_t(20.5) })
         hopperCrashCells.emplace_back(real_t(267.5), y, z);
   constexpr uint_t HOPPER_PROBE_FROM = 1000000000;
   constexpr uint_t HOPPER_PROBE_TO   = 0;

   // Confirms a checkpoint resume actually loaded real particle/step state, rather than silently
   // starting fresh with an empty particleStorage.
   {
      uint_t localOwnedParticleCount = uint_t(0);
      simObjs->particleStorage->forEachParticle(
         false, mesa_pd::kernel::SelectLocal(), *simObjs->particleAccessor,
         [&localOwnedParticleCount](size_t, auto&) { ++localOwnedParticleCount; }, *simObjs->particleAccessor);
      uint_t globalOwnedParticleCount = localOwnedParticleCount;
      mpi::allReduceInplace(globalOwnedParticleCount, mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("[CHECKPOINT-RESUME] loadSnapshot=" << (simParams->blockForest.loadSnapshot ? "true" : "false")
         << " beginTimeStep=" << simObjs->beginTimeStep << " globalOwnedParticleCount=" << globalOwnedParticleCount
         << " -- confirms whether this run actually resumed real state from a checkpoint (nonzero "
            "beginTimeStep + a particle count matching the expected population) or started fresh.");
   }

   // One-time domain-wide density/NaN scan of the just-loaded checkpoint, before the main loop's
   // first iteration touches anything -- reports the checkpoint's own saved fields.
   {
      const bool checkpointHasNaN = checkPdfNaN(simObjs, "CHECKPOINT-LOAD", simObjs->beginTimeStep,
         real_t(0), real_t(0), real_t(0), simParams->physics.densityPrecursorLow,
         simParams->physics.densityPrecursorHigh, simParams->physics.densityPrecursorLiquidOnly);
      WALBERLA_LOG_INFO_ON_ROOT("[CHECKPOINT-CLEAN-CHECK] beginTimeStep=" << simObjs->beginTimeStep
         << " localNaNFound=" << (checkpointHasNaN ? "true" : "false")
         << " -- one-time scan of the just-loaded checkpoint's own fields, before any new physics runs. "
            "false (with no further [NaN-SCAN]/[DENSITY-WINDOW] warnings above) means clean.");
   }

   WALBERLA_LOG_INFO_ON_ROOT("beginTimeStep: " << simObjs->beginTimeStep << " / " << simParams->physics.timesteps);
   // ----------------------------------------------------------------------------------------------------------------------------------
   for (uint_t t = simObjs->beginTimeStep; t != simParams->physics.timesteps; ++t)
   {
      *simObjs->currentMainStep = t;

      releaseSettledFloatingSpheres(simParams, simObjs, t);

      const bool watchNow = (t % WATCH_FREQ == uint_c(0)) || (t >= WATCH_DENSE_FROM);

      // Independent, prm-toggleable per-step NaN scan window (denseNanScanEnabled/From/ToStep in
      // PhysicsParameters) -- see PhysicsParams' own comment for why this is separate from watchNow.
      const bool denseNanScanNow = simParams->physics.denseNanScanEnabled &&
         t >= simParams->physics.denseNanScanFromStep && t < simParams->physics.denseNanScanToStep;
      const bool nanScanNow = watchNow || denseNanScanNow;

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] FSLBM step");
      timeloopFSLBM->singleStep(timingPool, true);
      if (watchNow) trackWatchCells(simObjs, "post-FSLBM", t, watchCells);
      if (watchNow) checkClosedInterfaceLayer(simObjs, "post-FSLBM", t);
      if (nanScanNow && checkPdfNaN(simObjs, "timeloopFSLBM", t, simParams->physics.lowRhoPrecursorThreshold, simParams->physics.highVelocityPrecursorThreshold, simParams->physics.populationPrecursorThreshold,
         simParams->physics.densityPrecursorLow, simParams->physics.densityPrecursorHigh, simParams->physics.densityPrecursorLiquidOnly))
      { listAllBadCellsNoEarlyExit(simObjs, "timeloopFSLBM"); WALBERLA_ABORT("[NaN-SCAN] aborting after timeloopFSLBM at t=" << t); }
      if (simParams->model.highAltitudeWaterScanEnabled)
         checkHighAltitudeWater(simObjs, "timeloopFSLBM", t, simParams->model.highAltitudeWaterThresholdZ);

      const bool hopperProbeNow = t >= HOPPER_PROBE_FROM && t <= HOPPER_PROBE_TO;

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] DEM sub-cycles");
      particleKernels.reduceAndAvgHydroForces(timingPool, t==0);
      particleKernels.doDEMSubCycles(timingPool, simParams->particle.numSubCycles);
      particleKernels.resetHydroForces(timingPool);
      if (hopperProbeNow) trackWatchCells(simObjs, "HOPPER-after-DEMSubCycles", t, hopperCrashCells);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] LBM-DEM coupling");
      timeloopAfterDEM->singleStep(timingPool, true);
      if (hopperProbeNow) trackWatchCells(simObjs, "HOPPER-after-timeloopAfterDEM", t, hopperCrashCells);
      if (watchNow) trackWatchCells(simObjs, "post-AfterDEM", t, watchCells);
      if (watchNow) checkClosedInterfaceLayer(simObjs, "post-AfterDEM", t);
      if (watchNow) logFloaterLinkBalance(simObjs, "post-AfterDEM", t);
      if (nanScanNow && checkPdfNaN(simObjs, "timeloopAfterDEM", t, simParams->physics.lowRhoPrecursorThreshold, simParams->physics.highVelocityPrecursorThreshold, simParams->physics.populationPrecursorThreshold,
         simParams->physics.densityPrecursorLow, simParams->physics.densityPrecursorHigh, simParams->physics.densityPrecursorLiquidOnly))
      { listAllBadCellsNoEarlyExit(simObjs, "timeloopAfterDEM"); WALBERLA_ABORT("[NaN-SCAN] aborting after timeloopAfterDEM at t=" << t); }
      if (simParams->model.highAltitudeWaterScanEnabled)
         checkHighAltitudeWater(simObjs, "timeloopAfterDEM", t, simParams->model.highAltitudeWaterThresholdZ);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Evaluators");
      timeloopEval->singleStep(timingPool, true);
      if (hopperProbeNow) trackWatchCells(simObjs, "HOPPER-after-timeloopEval", t, hopperCrashCells);

      // extra end-of-step fluid sanity sweep (ModelParameters/PhysicsParameters fluidSanity*, default off)
      if (simParams->physics.fluidSanityEnabled && t >= simParams->physics.fluidSanityFromStep && t < simParams->physics.fluidSanityToStep &&
          (simParams->physics.fluidSanityEveryNSteps <= uint_t(1) || t % simParams->physics.fluidSanityEveryNSteps == uint_t(0)))
         fluidSanitySweep(simObjs, simParams, t);

      if (t % simParams->evaluation.evaluationFrequency == uint_c(0) && t > uint_c(0))
      {
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] acceleration: " << (*simObjs->acceleration)[0] << ", meanVelocity: " << (*simObjs->meanVelocity)[0]);
         writeFluidInfo(simParams, simObjs, t);

         WALBERLA_ROOT_SECTION()
         {
            sqlite::storeRunInSqliteDB(simParams->evaluation.baseFolderName + "/antidunes.sqlite",
               simObjs->dataCols->getIntCols(),
               simObjs->dataCols->getStrCols(),
               simObjs->dataCols->getRealCols());
         }
      }

      if (t % simParams->evaluation.performanceLogFrequency == uint_c(0) && t > uint_c(0))
      {
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Performance log");
         timingPool.logResultOnRoot();

         if (simParams->evaluation.logContactDetectionStatsPeriodically)
         {
            particleKernels.logContactDetectionStatistics(std::to_string(t));
            particleKernels.resetContactDetectionStatistics();
         }
      }

      if (simParams->blockForest.storeSnapshot &&
         t % simParams->blockForest.snapshotFrequency == uint_c(0) && t > uint_c(0))
      {
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Writing checkpoint");
         writeCheckpoint(simParams, simObjs, t, false, simParams->blockForest.keepLastTwoSnapshots);
      }
      else if (simParams->blockForest.storeSnapshot && simParams->blockForest.keepAllSnapshots &&
         t >= simParams->blockForest.keepAllSnapshotsFromStep)
      {
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Writing per-step diagnostic checkpoint");
         writeCheckpoint(simParams, simObjs, t, true);
      }
   }

   particleKernels.logContactDetectionStatistics();

   return EXIT_SUCCESS;
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv) { return walberla::antidunes::main(argc, argv); }