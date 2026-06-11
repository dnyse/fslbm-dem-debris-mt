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

   // Init forceDensity (requires valid fill + pdf + flag fields and acceleration)
   // ----------------------------------------------------------------------------------------------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT("Initializing force density field");
   walberla::free_surface::initForceDensityFieldCodegen< PdfField_T, FlagField_T, VectorFieldFlattened_T, ScalarField_T >(
      simObjs->blockForest, simObjs->forceDensityFieldID, simObjs->fillFieldID, simObjs->pdfFieldID, simObjs->flagFieldID,
      simObjs->flagInfo, *simObjs->acceleration);

   // Communicate after initialization
   // ----------------------------------------------------------------------------------------------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT("Initial communication");
   Communication_T(simObjs->blockForest, simObjs->flagFieldID, simObjs->fillFieldID, simObjs->forceDensityFieldID)();
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

   WcTimingPool timingPool;

   // [DEBUG — disabled after fix confirmed at t=15027+]
   // Watch cells used to track PDF/density drift at the contact-circle cells that produced NaN crashes.
   // All three are at distance 7.9215 from sphere centre (15,30,56), radius=8 — the outermost
   // sphere-interior layer, uncovered first as the sphere oscillates at the free surface.
   //   Run 1 (no sweeps):   NaN at t=581 in cells (8.5,29.5,60.5) and (8.5,34.5,56.5)
   //   Run 2 (with sweeps): NaN at t=577 in cell  (7.5,27.5,56.5)
   // NaN-trace instrumentation (ENABLED).
   //   * checkPdfNaN after each PDF-touching phase: locates the first NaN, dumps its 3x3x3
   //     neighbourhood (flag/rho/fill = the real contact-line config), and aborts cleanly.
   //   * trackWatchCells: per-step rho/fill/flag history of the watched cells. The coords below are
   //     a first guess on the giant-sphere contact ring (surface z=56, sphere drifts in x); after
   //     run 1 reports the exact NaN cell, set watchCells to that world coord and rerun for history.
   const std::vector< Vector3< real_t > > watchCells = {
      Vector3< real_t >(real_t(40.5), real_t(30.5), real_t(56.5)),   // near sphere centre at the waterline
      Vector3< real_t >(real_t(20.5), real_t(30.5), real_t(56.5)),   // contact ring -x
      Vector3< real_t >(real_t(60.5), real_t(30.5), real_t(56.5)),   // contact ring +x
      Vector3< real_t >(real_t(40.5), real_t(10.5), real_t(56.5)),   // contact ring -y
      Vector3< real_t >(real_t(40.5), real_t(50.5), real_t(56.5)),   // contact ring +y
   };
   constexpr uint_t WATCH_FREQ      = 1000;        // periodic sampling (lightweight)
   constexpr uint_t WATCH_DENSE_FROM = 1000000000; // per-step dense tracking effectively off; lower it to debug a crash

   WALBERLA_LOG_INFO_ON_ROOT("beginTimeStep: " << simObjs->beginTimeStep << " / " << simParams->physics.timesteps);
   // ----------------------------------------------------------------------------------------------------------------------------------
   for (uint_t t = simObjs->beginTimeStep; t != simParams->physics.timesteps; ++t)
   {
      releaseSettledFloatingSpheres(simParams, simObjs, t);

      const bool watchNow = (t % WATCH_FREQ == uint_c(0)) || (t >= WATCH_DENSE_FROM);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] FSLBM step");
      timeloopFSLBM->singleStep(timingPool, true);
      if (watchNow) trackWatchCells(simObjs, "post-FSLBM", t, watchCells);
      if (watchNow) checkClosedInterfaceLayer(simObjs, "post-FSLBM", t);
      if (watchNow && checkPdfNaN(simObjs, "timeloopFSLBM", t)) WALBERLA_ABORT("[NaN-SCAN] aborting after timeloopFSLBM at t=" << t);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] DEM sub-cycles");
      particleKernels.reduceAndAvgHydroForces(timingPool, t==0);
      particleKernels.doDEMSubCycles(timingPool, simParams->particle.numSubCycles);
      particleKernels.resetHydroForces(timingPool);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] LBM-DEM coupling");
      timeloopAfterDEM->singleStep(timingPool, true);
      if (watchNow) trackWatchCells(simObjs, "post-AfterDEM", t, watchCells);
      if (watchNow) checkClosedInterfaceLayer(simObjs, "post-AfterDEM", t);
      if (watchNow) logFloaterLinkBalance(simObjs, "post-AfterDEM", t);
      if (watchNow && checkPdfNaN(simObjs, "timeloopAfterDEM", t)) WALBERLA_ABORT("[NaN-SCAN] aborting after timeloopAfterDEM at t=" << t);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Evaluators");
      timeloopEval->singleStep(timingPool, true);

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
      }

      if (simParams->blockForest.storeSnapshot &&
         t % simParams->blockForest.snapshotFrequency == uint_c(0) && t > uint_c(0))
      {
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Writing checkpoint");
         writeCheckpoint(simParams, simObjs, t);
      }
   }

   return EXIT_SUCCESS;
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv) { return walberla::antidunes::main(argc, argv); }