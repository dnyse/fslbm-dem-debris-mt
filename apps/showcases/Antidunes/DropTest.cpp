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
//! \file DropTest.cpp
//
// Fluid-free, DEM-only cuboid (and, as of 2026-09-13, cylinder) drop test at real 72-process
// block-forest scale. Built to isolate
// whether the "flying floater" / contact-chatter / STALE-PARTICLE-IDX desync bug chased across this
// session's Antidunes bridge-scenario cluster jobs (4201656/4202194/4212104/4212265) is a general
// DEM-contact/parallel-decomposition problem, independent of any fluid coupling, rather than
// something specific to the bridge geometry or the LBM-DEM coupling itself.
//
// Reuses Antidunes' own production DEM machinery completely unmodified: setupBlockForest,
// setupParticles, releaseSettledFloatingSpheres (from SetupStages_drop_test.h -- a copy of
// SetupStages.h with one small addition, X-wall creation; see that file's own header comment for why
// it's a copy rather than a direct edit) and ParticleKernels::doDEMSubCycles (from ParticleKernels.h,
// completely unmodified) for the actual contact detection/resolution/integration. No LBM fields are
// ever created or registered; this is safe because doDEMSubCycles' only fluid-field-dependent code
// (computeSubmergedVolumes()/computeSubmergedVolumesBox()) is gated behind
// ModelParameters.useMemOnlyBuoyancy, which defaults to true (and must stay true here) -- with it
// true, AddExternalForces applies plain density-weighted mass*gravity along Z (the only axis with a
// correct, density-weighted force law in this codebase -- X/Y carry a different, non-density-weighted
// "driving force" term used for the production sim's periodic-flow forcing, so Z must stay the
// gravity/height axis here).
//
// setupSimState() is deliberately NOT called (it unconditionally constructs a PIDController, which
// is both meaningless for a non-flow-driven drop test and an unnecessary link dependency). Instead
// this file does the two lines of state setupSimState would otherwise provide (beginTimeStep,
// acceleration) directly.
//
// VTK output is deliberately minimal and NOT BuildTimeloop.h::addBeforeTimeStepVTKOutput (which
// requires real LBM/FlagField fields and bundles in fluid_field/particles/floater-sphere/
// floater-cylinder writers this fluid-free tool has no use for): just one domain_decomposition dump
// at start (as BedGeneration.cpp already does) plus floaterBox/floaterCylinder-only
// ParticleVtkOutput blocks (the exact same selector/writer fields as BuildTimeloop.h's own "floating
// boxes (type 3)"/"floating cylinders (type 4)" blocks), each written every prm-configured
// VTK.writeFrequency steps, and each a no-op if the prm has no floaters of that shape.
//======================================================================================================================

#include "core/Environment.h"
#include "core/waLBerlaBuildInfo.h"

// Both needed only because SetupStages_drop_test.h (a full copy of SetupStages.h) still defines
// setupLBMFields/setupSurfaceGeometryFields/setupBoundaryHandling/mapParticlesOntoCells -- ordinary
// (non-template) functions this file never CALLS but that must still fully compile. Antidunes.cpp
// gets these two transitively via BuildTimeloop.h, which DropTest.cpp deliberately doesn't include.
#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionModel.h"

#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/VTKOutput.h"

#include "ParticleKernels.h"
#include "SimulationObjects.h"
#include "SimulationParams.h"
#include "SetupStages_drop_test.h"

#include <chrono>
#include <filesystem>

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
   // ----------------------------------------------------------------------------------------------------------------
   const auto simParams = std::make_shared< SimulationParams >(configPtr);
   const auto simObjs   = std::make_shared< SimulationObjects >();

   WALBERLA_ROOT_SECTION()
   {
      std::filesystem::create_directories(simParams->blockForest.snapshotBaseFolder);
      std::filesystem::create_directories(simParams->evaluation.baseFolderName);
      if (simParams->model.logFloaterHistory)
         std::filesystem::create_directories(simParams->model.floaterHistoryBaseFolder);
   }

   simParams->computed.logOnRoot();
   simParams->computed.writeToFile(simParams->evaluation.baseFolderName);

   // Setup: block forest + particles only -- no LBM fields, no boundary handling, no bubble model.
   // ----------------------------------------------------------------------------------------------------------------
   setupBlockForest(simParams, simObjs);

   // Minimal sim-state init -- NOT setupSimState() (SetupStages_drop_test.h), which unconditionally
   // constructs a PIDController: meaningless here (no flow to drive) and an avoidable link dependency.
   simObjs->beginTimeStep = uint_t(0);
   simObjs->acceleration =
      std::make_shared< Vector3< real_t > >(real_t(0), real_t(0), -simParams->computed.gravity);

   setupParticles(simParams, simObjs);

   // DEM ---------------------------------------------------------------------------------------------------------
   ParticleKernels< ParticleAccessor_T > particleKernels(simParams, simObjs);

   // Minimal VTK: domain decomposition (once) + floaterBox only (periodic).
   // ----------------------------------------------------------------------------------------------------------------
   // NOT wrapped in WALBERLA_ROOT_SECTION -- createVTKOutput_DomainDecomposition's write() is a
   // COLLECTIVE call (every rank contributes its own local block's info); an earlier version of this
   // code incorrectly called it root-only, which deadlocked at any process count > 1 (rank 0 blocks
   // forever waiting for contributions from ranks that skip the call entirely, while those ranks race
   // ahead into the main loop and hit a different collective -- classic MPI collective mismatch).
   // Matches BedGeneration.cpp's own reference usage, which also calls this unconditionally.
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition(
      simObjs->blockForest, "domain_decomposition", 1, simParams->vtk.baseFolder, "simulation_step");
   vtkDomainOutput->write();

   std::function< void() > writeFloaterBoxVtk = []() {};
   if (!simParams->floatingBoxes.empty())
   {
      struct SelectFloaterDensityRatio
      {
         using return_type = real_t;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto& shape = p.getBaseShape();
            if (!shape || shape->getVolume() < real_t(1e-12)) return real_t(0);
            return shape->getMass() / shape->getVolume();
         }
      };
      struct SelectBoxEdgeLength
      {
         using return_type = mesa_pd::Vec3;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto& shape = p.getBaseShape();
            if (!shape || shape->getShapeType() != mesa_pd::data::Box::SHAPE_TYPE)
               return mesa_pd::Vec3(real_t(0));
            return static_cast< const mesa_pd::data::Box* >(shape.get())->getEdgeLength();
         }
      };
      // Transpose of the body->world rotation matrix, matching BuildTimeloop.h's own floaterBox
      // output exactly -- see that file's comment for why (ParaView Tensor Glyph column-major order).
      struct SelectParticleRotationMatrix
      {
         using return_type = std::array< real_t, 9 >;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const mesa_pd::Mat3 R = p.getRotation().getMatrix();
            return { R(0, 0), R(1, 0), R(2, 0),
                     R(0, 1), R(1, 1), R(2, 1),
                     R(0, 2), R(1, 2), R(2, 2) };
         }
      };

      auto boxFloaterData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      boxFloaterData->addOutput< SelectBoxEdgeLength >("edgeLength");
      boxFloaterData->addOutput< SelectParticleRotationMatrix >("rotation");
      boxFloaterData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      boxFloaterData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      boxFloaterData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      boxFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(3) && !isSet(pIt->getFlags(), GHOST);
      });
      auto boxFloaterVTKOutput = vtk::createVTKOutput_PointData(
         boxFloaterData, "floaterBox", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      writeFloaterBoxVtk = vtk::writeFiles(boxFloaterVTKOutput);
   }

   // Cylinder floaters (type 4) -- same treatment as the box block above, mirroring
   // BuildTimeloop.h's own "floaterCylinder" VTK block (SelectCylinderGlyphScale/rotation-selector
   // pattern) exactly, so this fluid-free tool renders cylinders identically to production.
   std::function< void() > writeFloaterCylinderVtk = []() {};
   if (!simParams->floatingCylinders.empty())
   {
      struct SelectFloaterDensityRatio
      {
         using return_type = real_t;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto& shape = p.getBaseShape();
            if (!shape || shape->getVolume() < real_t(1e-12)) return real_t(0);
            return shape->getMass() / shape->getVolume();
         }
      };
      struct SelectCylinderGlyphScale
      {
         using return_type = mesa_pd::Vec3;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto* shape = dynamic_cast< const mesa_pd::data::Cylinder* >(p.getBaseShape().get());
            if (!shape) return mesa_pd::Vec3(real_t(0));
            return mesa_pd::Vec3(real_t(2) * shape->getRadius(), shape->getHeight(),
                                 real_t(2) * shape->getRadius());
         }
      };
      // (R * Rx(+90))^T = Rx(-90) * R^T -- glyph-local Y is the cylinder axis, body-local Z is; see
      // BuildTimeloop.h's own identical selector for the full derivation.
      struct SelectParticleRotationMatrix
      {
         using return_type = std::array< real_t, 9 >;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const mesa_pd::Mat3 R = p.getRotation().getMatrix();
            return { R(0, 0),  R(1, 0),  R(2, 0),
                     R(0, 2),  R(1, 2),  R(2, 2),
                    -R(0, 1), -R(1, 1), -R(2, 1) };
         }
      };

      auto cylinderFloaterData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      cylinderFloaterData->addOutput< SelectCylinderGlyphScale >("glyphScale");
      cylinderFloaterData->addOutput< SelectParticleRotationMatrix >("rotation");
      cylinderFloaterData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      cylinderFloaterData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      cylinderFloaterData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      cylinderFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(4) && !isSet(pIt->getFlags(), GHOST);
      });
      auto cylinderFloaterVTKOutput = vtk::createVTKOutput_PointData(
         cylinderFloaterData, "floaterCylinder", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      writeFloaterCylinderVtk = vtk::writeFiles(cylinderFloaterVTKOutput);
   }

   // Main loop -- no LBM, no LBM-DEM coupling, no PID. Just wave release + DEM sub-cycles + VTK.
   // ----------------------------------------------------------------------------------------------------------------
   WcTimingPool timingPool;

   WALBERLA_LOG_INFO_ON_ROOT("beginTimeStep: " << simObjs->beginTimeStep << " / " << simParams->physics.timesteps);
   // Per-step wall-clock timing breakdown (release / DEM sub-cycles / VTK) -- added to diagnose an
   // observed multi-rank slowdown (72-process cluster and 12-process local runs both took many
   // minutes to clear even step 1000 for a 1000-cuboid config, vs. 5.4s single-process for 48
   // cuboids) that the existing every-evaluationFrequency-steps progress log was too coarse to see
   // at all (it never even fired before the slow runs were cancelled). Logged every step,
   // unconditionally -- this is a diagnostic instrumentation pass, not meant to stay this verbose
   // once the bottleneck is identified.
   auto stepTimerStart = std::chrono::steady_clock::now();
   for (uint_t t = simObjs->beginTimeStep; t != simParams->physics.timesteps; ++t)
   {
      *simObjs->currentMainStep = t;

      const auto tStepStart = std::chrono::steady_clock::now();

      releaseSettledFloatingSpheres(simParams, simObjs, t);
      const auto tAfterRelease = std::chrono::steady_clock::now();

      particleKernels.doDEMSubCycles(timingPool, simParams->particle.numSubCycles);
      const auto tAfterDEM = std::chrono::steady_clock::now();

      writeFloaterBoxVtk();
      writeFloaterCylinderVtk();
      const auto tAfterVtk = std::chrono::steady_clock::now();

      const double releaseMs = std::chrono::duration< double, std::milli >(tAfterRelease - tStepStart).count();
      const double demMs     = std::chrono::duration< double, std::milli >(tAfterDEM - tAfterRelease).count();
      const double vtkMs     = std::chrono::duration< double, std::milli >(tAfterVtk - tAfterDEM).count();
      const double totalMs   = std::chrono::duration< double, std::milli >(tAfterVtk - tStepStart).count();
      const double sinceStartS = std::chrono::duration< double >(tAfterVtk - stepTimerStart).count();
      WALBERLA_LOG_INFO_ON_ROOT("[StepTiming] step=" << t << " wallSinceStart=" << sinceStartS << "s"
         << " release=" << releaseMs << "ms dem=" << demMs << "ms vtk=" << vtkMs << "ms total=" << totalMs << "ms");

      if (t % simParams->evaluation.evaluationFrequency == uint_c(0) && t > uint_c(0))
      {
         real_t avgVel = real_t(0), maxVel = real_t(0), maxHeight = real_t(0);
         uint_t numParticles = uint_t(0);
         getAverageVelocity(*simObjs->particleAccessor, avgVel, maxVel, numParticles, maxHeight);
         WALBERLA_LOG_INFO_ON_ROOT("[" << t << "/" << simParams->physics.timesteps
            << "] numParticles=" << numParticles << " avgVel=" << avgVel << " maxVel=" << maxVel
            << " maxHeight=" << maxHeight);
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
   }

   particleKernels.logContactDetectionStatistics();

   return EXIT_SUCCESS;
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv) { return walberla::antidunes::main(argc, argv); }
