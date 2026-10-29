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
//! \file tsunami.cpp
//
// SEPARATE executable, built alongside (never instead of) Antidunes -- the working Antidunes.cpp
// pipeline is left completely untouched. This showcase recreates the "water flows into a dry domain"
// half of Watanabe et al. 2023's validation (Section 3/4: a raised reservoir behind a gate, released
// under gravity across a dry run-out floor into a pole forest / driftwood field), which the existing
// Antidunes-based toy scenarios (prms/toy_pole_forest_driftwood_gas.prm) explicitly could NOT do --
// that file's own header notes it instead treats the forest as sitting in an already-flowing, already-
// wet steady channel (Antidunes' native PID-driven periodic-channel setup), since a genuine dry-bed
// dam-break needs a different initial condition and no flow-rate controller.
//
// This file is Antidunes.cpp's own structure, reused almost line-for-line (same setup-stage call order,
// same three-timeloop / DEM-sub-cycle main loop, same checkpointing/evaluation/VTK machinery) -- see
// the inline comments below for the small number of places this genuinely differs:
//   1. setupLBMFieldsDamBreak (TsunamiDataInitialization.h) replaces setupLBMFields (SetupStages.h):
//      a dry-domain dam-break fill level + still-water PDF field, instead of the sinusoidal free
//      surface + Poiseuille shear profile + hydrostatic-for-a-mean-flow pressure Antidunes.cpp uses.
//   2. addAfterTimeStepForcingAdjuster (BuildTimeloop.h/Forcing.h) is never added to the eval timeloop
//      -- there is no flow-rate PID controller here; gravity alone drives the whole event, exactly like
//      the real experiment (the gate opens and the reservoir simply falls).
// Every other setup stage (block forest, particles/pillars/arch tunnels/floaters, boundary handling,
// DEM sub-cycling via ParticleKernels, VTK output, checkpointing, SQLite logging) is the EXACT same
// function, called the EXACT same way, reused directly from SetupStages.h/BuildTimeloop.h/
// ParticleKernels.h/WriteData.h -- none of those files were modified for this executable.
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
#include "TsunamiDataInitialization.h"

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

   // Dam-break reservoir geometry (see TsunamiDataInitialization.h) -- the one genuinely new prm block
   // this executable needs, parsed separately since SimulationParams.h itself is left untouched.
   const DamBreakParams damBreakParams(configPtr->getOneBlock("DamBreakParameters"),
                                        real_c(simParams->domain.domainSize[1]));

   // Create directories for storing snapshots and evaluations
   // ----------------------------------------------------------------------------------------------------------------------------------
   WALBERLA_ROOT_SECTION()
   {
      std::filesystem::create_directories(simParams->blockForest.snapshotBaseFolder);
      std::filesystem::create_directories(simParams->evaluation.baseFolderName);
      if (simParams->model.logFloaterHistory)
         std::filesystem::create_directories(simParams->model.floaterHistoryBaseFolder);
   }

   // Log computed params
   // ----------------------------------------------------------------------------------------------------------------------------------
   simParams->computed.logOnRoot();
   simParams->computed.writeToFile(simParams->evaluation.baseFolderName);

   // Setup simulation (initialize or load from snapshot)
   // ----------------------------------------------------------------------------------------------------------------------------------
   setupBlockForest           (simParams, simObjs);
   setupSimState              (simParams, simObjs);   // beginTimestep + acceleration (gravity only, see setupSimState's own accelerationX formula -- 0 here since PIDParameters.targetMeanVelocityMagnitude=0) + pid (constructed but never stepped -- no ForcingAdjuster below)
   setupVTKState              (simParams, simObjs);   // VTK counter
   setupLBMFieldsDamBreak     (simParams, simObjs, damBreakParams);   // dry-domain dam-break IC (see this executable's own header)
   setupSurfaceGeometryFields (simParams, simObjs);   // curvature + normal + obstacleNormal + smoothFill
   setupParticles             (simParams, simObjs);   // particleStorage + particleField (pillars/arch tunnels/floaters -- unchanged from Antidunes)
   addTsunamiXWallPlanes      (simParams, simObjs);   // DEM wall planes at x=0/domainSize[0] -- see this function's own comment (TsunamiDataInitialization.h): setupParticles only ever adds these for Z/Y, never X

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
   // NOTE: addAfterTimeStepForcingAdjuster is deliberately NOT added here -- see this file's own header
   // comment. There is no target flow rate to hold; acceleration stays fixed at setupSimState's initial
   // (0, 0, -gravity) for the whole run, and the dam-break event is driven by gravity alone.

   WcTimingPool timingPool;

   WALBERLA_LOG_INFO_ON_ROOT("beginTimeStep: " << simObjs->beginTimeStep << " / " << simParams->physics.timesteps);
   // ----------------------------------------------------------------------------------------------------------------------------------
   for (uint_t t = simObjs->beginTimeStep; t != simParams->physics.timesteps; ++t)
   {
      *simObjs->currentMainStep = t;

      releaseSettledFloatingSpheres(simParams, simObjs, t);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] FSLBM step");
      timeloopFSLBM->singleStep(timingPool, true);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] DEM sub-cycles");
      particleKernels.reduceAndAvgHydroForces(timingPool, t==0);
      particleKernels.doDEMSubCycles(timingPool, simParams->particle.numSubCycles);
      particleKernels.resetHydroForces(timingPool);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] LBM-DEM coupling");
      timeloopAfterDEM->singleStep(timingPool, true);

      WALBERLA_LOG_PROGRESS_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] Evaluators");
      timeloopEval->singleStep(timingPool, true);

      if (t % simParams->evaluation.evaluationFrequency == uint_c(0) && t > uint_c(0))
      {
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps << "] meanVelocity: " << (*simObjs->meanVelocity)[0]);
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

   particleKernels.logContactDetectionStatistics();

   return EXIT_SUCCESS;
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv) { return walberla::antidunes::main(argc, argv); }
