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

   // [Claude, 2026-10-29] TEMPORARY DIAGNOSTIC (THESIS_PLAN.md item 1a, checkpoint-resume
   // STALE-PARTICLE-IDX abort, job 4245028) commented out -- that investigation is closed (the
   // underlying ghost-sync gap was fixed in SetupStages.h's own 4-round syncCall()), and this was
   // the only genuinely unconditional debug print left in the production codebase after a survey
   // (per explicit request to find/clean up leftover debugging prints). Kept for reference below.
#if 0
   // per-rank dump of what the particle accessor actually knows right after setupParticles() (i.e.
   // after loadBlockData + syncCall()), to see which ranks have a ghost/local copy of each floater
   // and which don't, cross-referenced against the block's own AABB (extended by
   // mesapdDomain-relevant margins) to understand the ghost-sync reach.
   //
   // [Claude, 2026-09-18]: un-gated from loadSnapshot-only. SetupStages.h's own comment right above
   // its 4-round syncCall() fix (search "ROOT CAUSE FIX (THESIS_PLAN.md item 1a...") argues a fresh
   // start "never needs more than one hop yet (freshly-created floaters/beds are local to a handful
   // of blocks)" -- but SyncGhostOwners is a one-hop-per-call BFS propagation (only reaches ranks in
   // domain.getNeighborProcesses() of wherever the particle is ALREADY known), and job 4245028's own
   // box (halfExtents<15,30,15> at <48,67.5,60>, cellsPerBlock<16,45,30>) spans 2x3x2=12 blocks just
   // axis-aligned -- plausibly more than one hop from whichever single rank owns it, on a FRESH
   // start too. That run hit the exact same [STALE-PARTICLE-IDX] abort this diagnostic was built for
   // (job 4245028, rank 13/17/21, t=0, the very first FSLBM sweep -- before any of the 49999 further
   // per-step syncCall()s the old comment was counting on to eventually converge). Enabling this
   // unconditionally to check whether the affected ranks show zero knowledge of the box right after
   // setupParticles()'s single fresh-start syncCall(), confirming (or refuting) that fresh starts can
   // hit the identical multi-hop gap already fixed for resumes.
   {
      const int rank = mpi::MPIManager::instance()->rank();
      uint_t totalKnown = uint_c(simObjs->particleAccessor->size());
      std::ostringstream uidList;
      simObjs->particleStorage->forEachParticle(false, mesa_pd::kernel::SelectAll(), *simObjs->particleAccessor,
         [&uidList](const size_t idx, auto& ac) {
            if (ac.getType(idx) == uint_t(3))
               uidList << ac.getUid(idx) << "(pos=" << ac.getPosition(idx) << ") ";
         }, *simObjs->particleAccessor);
      // Print each of this rank's own local blocks' AABB for comparison against the box's own position.
      std::ostringstream blockAABBs;
      for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
         blockAABBs << blockIt->getAABB() << " ";
      WALBERLA_LOG_WARNING("[GHOST-SYNC-DIAG] rank=" << rank << " totalKnownParticles=" << totalKnown
         << " boxes-known=[ " << uidList.str() << "] blocks=[ " << blockAABBs.str() << "]");
   }
#endif

   // Boundary handling and flags
   // ----------------------------------------------------------------------------------------------------------------------------------
   setupBoundaryHandling      (simParams, simObjs);
   mapParticlesOntoCells      (simParams, simObjs);   // set solid cells in flagField (+ map UIDs to particleField)

   // Set liquid/interface/gas cells in flagField (must be called after all solid cells are set)
   // ----------------------------------------------------------------------------------------------------------------------------------
   simObjs->antidunesBoundaryHandling->initFlagsFromFillLevel();
   checkOuterGhostLayerFlags(simObjs, "after initFlagsFromFillLevel"); // [Claude, 2026-10-03] see BuildTimeloop.h
   if (simParams->model.wallOuterGhostLayer) // [Claude, 2026-10-03] see SimulationParams.h
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
   // alongside flagFieldID (see addSweepFluidParticleCoupling's own "Communication (coupling): after
   // PDF reconstruction (flagField + particleField)" AfterFunction) -- omitting it here left a rank
   // whose block only sees a floater as a GHOST (not locally mapped) with a correctly-communicated
   // "obstacle" flag but a stale/never-communicated particleField entry (default invalid UID=0) for
   // that same cell, until the first per-step coupling sweep runs. On a fresh start this is a
   // harmless one-step gap (same reasoning as the fillWithDomain() fix in mapParticlesOntoCells());
   // on a checkpoint resume, if a floater's cross-rank obstacle footprint at load time happens to
   // include a cell on a rank that only ever sees it as a ghost, the SimpleBB.h MEM force sweep's own
   // stale-particle-idx guard hard-aborts on that very first post-load step (confirmed directly: job
   // 4236616, `[STALE-PARTICLE-IDX] ... uid=0` at t=70001, the first step after loading
   // cuboid_stillwater_realgravity981_3xsphere_hires_d025's own t=70000 checkpoint).
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
   // [Claude, 2026-09-24] Repurposed for the tilt26.565deg/rho=0.7665 step-444 crash: the exact
   // failing cell was pinned down at world (20.5,8.5,25.5) (checkPdfNaN + denseNanScanEnabled, see
   // THESIS_PLAN.md) -- a genuine INTF cell (fill=0.5) at the box's own top corner, surrounded by
   // OBST (the box) and GAS with zero LIQ/INTF neighbours. Watching the full 3x3x3 world-coordinate
   // cube around it (27 cells) to get a per-step history across the run-up to the crash, to settle
   // whether this cell's own bognerTargetFlag classification ever ran (it should have resolved to
   // GAS per classifyFormerObstacleNeighbourhood's hasGas&&!hasLiquid&&!hasUnknown branch) or whether
   // it is a pre-existing interface cell whose neighbourhood only became fully gas-locked later,
   // without ever being re-classified (that reclassification path is formerObstacle-only, not a
   // per-step re-scan of every existing interface cell).
   std::vector< Vector3< real_t > > watchCells;
   for (real_t dz : { real_t(-1), real_t(0), real_t(1) })
      for (real_t dy : { real_t(-1), real_t(0), real_t(1) })
         for (real_t dx : { real_t(-1), real_t(0), real_t(1) })
            watchCells.emplace_back(real_t(20.5) + dx, real_t(8.5) + dy, real_t(25.5) + dz);
   constexpr uint_t WATCH_FREQ      = 1000;   // periodic sampling (lightweight) -- unchanged elsewhere
   // [Claude, 2026-09-25] DISABLED (raised past any realistic timesteps) for the full 60k-step Option-B
   // validation run -- the Bogner-box step-444 crash this dense per-step tracking was built for is long
   // since resolved (THESIS_PLAN.md), and per-step tracking from step 425 onward would flood the log
   // over 60k steps of an unrelated (hopper) scenario. WATCH_FREQ's lightweight periodic sampling above
   // still runs every 1000 steps.
   constexpr uint_t WATCH_DENSE_FROM = 1000000000;

   // [Claude, 2026-09-25] Hopper-clogging crash region (THESIS_PLAN.md): the FSLBM-internal stages
   // (boundarySweep..equilibriumRefilling, wired in BuildTimeloop.h's probeHopperCrashRegion) show this
   // OBST-flagged cell's PDFs completely unchanged by anything except collide -- yet its rho measurably
   // drifts step-to-step (e.g. 1.007356 at one step's end -> 1.006870 at the very next step's
   // boundarySweep, with force staying exactly zero throughout). Nothing in the FSLBM timeloop can cause
   // that: the only code that runs in between is DEM sub-cycling + timeloopAfterDEM (particle mapping +
   // PDF reconstruction, i.e. MovingParticleMapping.h / PdfReconstructionManager.h) -- never probed at
   // this exact location before. Same X/Y/Z coverage as probeHopperCrashRegion for direct comparability.
   std::vector< Vector3< real_t > > hopperCrashCells;
   for (real_t y = real_t(5.5); y <= real_t(35.5); y += real_t(1))
      for (real_t z : { real_t(17.5), real_t(18.5), real_t(19.5), real_t(20.5) })
         hopperCrashCells.emplace_back(real_t(267.5), y, z);
   // [Claude, 2026-09-25] DISABLED (impossible range) for the full 60k-step Option-B validation run --
   // the mechanism this probe traced is fully characterised and the fix is validated (THESIS_PLAN.md).
   constexpr uint_t HOPPER_PROBE_FROM = 1000000000;
   constexpr uint_t HOPPER_PROBE_TO   = 0;

   // [Claude, 2026-09-26] Per explicit request ("add logs to confirm it rerun from the
   // checkpoints"): an explicit, unambiguous confirmation that a checkpoint resume actually loaded
   // real particle/step state, rather than silently starting fresh at step 0 with an empty
   // particleStorage (which loadSnapshot=true + a missing/misconfigured snapshotBaseFolder could do
   // without erroring). Counts LOCAL owned (non-ghost) particles and sums across all ranks for a
   // true global total -- comparing this against the known population size (e.g. 24 for the
   // giant/small rod doubling-series scenarios, plus 2 pillars) directly confirms the checkpoint's
   // particleStorage was actually read, not just its step counter.
   //
   // Re-enabled, 2026-09-27, per explicit request ("check if the initial checkpoint is clean
   // everything looks stable"). Was briefly disabled while isolating the gated-guard hang -- its own
   // mpi::allReduceInplace collective was already confirmed unrelated to that hang (it completed and
   // printed fine in the hung job 4270904's own log).
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

   // [Claude, 2026-09-27] Per explicit request ("check if the initial checkpoint is clean everything
   // looks stable"): a ONE-TIME domain-wide density/NaN scan of the just-loaded state, BEFORE the
   // main loop's first iteration touches anything -- i.e. this reports the checkpoint's OWN saved
   // fields, not state produced by any new physics. Reuses checkPdfNaN entirely (same function the
   // diagnostic sweep and denseNanScan machinery both call): its raw-PDF and velocity-derived NaN
   // checks run unconditionally regardless of the threshold args, so this call answers "is there
   // already a NaN/Inf anywhere in this checkpoint" outright; passing the configured
   // densityPrecursorLow/High/LiquidOnly additionally flags any cell that's finite but already
   // suspiciously out of range. localNaN==false plus no [DENSITY-WINDOW]/[NaN-SCAN] output below
   // means the loaded checkpoint is clean.
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