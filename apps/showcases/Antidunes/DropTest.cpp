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
// Fluid-free, DEM-only drop test (cuboids, cylinders, spheres, capsules) on a real multi-process block forest: checks
// DEM contact and parallel-decomposition behaviour independently of the fluid coupling.
//
// Uses the production DEM code: setupBlockForest, setupParticles, releaseSettledFloatingSpheres (from
// SetupStages_drop_test.h, a variant of SetupStages.h with X walls, container cylinders and capsules; see its header)
// and ParticleKernels::doDEMSubCycles. No LBM fields are created; this is safe because the only fluid-dependent code
// in doDEMSubCycles (computeSubmergedVolumes()/computeSubmergedVolumesBox()) is skipped with
// ModelParameters.useMemOnlyBuoyancy=true (the default, required here). Gravity then acts as density-weighted
// mass*gravity along Z, while X/Y carry the production driving-force term, so Z must stay the gravity axis.
//
// setupSimState() is not called (it constructs a PIDController); beginTimeStep and acceleration are set directly.
//
// VTK output does not use BuildTimeloop.h::addBeforeTimeStepVTKOutput (which needs LBM fields): one
// domain_decomposition dump at start (as in BedGeneration.cpp) plus a floaterBox/floaterCylinder/floater/
// floaterCapsule ParticleVtkOutput per shape present, with the same selectors as BuildTimeloop.h.
//======================================================================================================================

#include "core/Environment.h"
#include "core/waLBerlaBuildInfo.h"

// Needed only because SetupStages_drop_test.h still defines setupLBMFields/setupSurfaceGeometryFields/
// setupBoundaryHandling/mapParticlesOntoCells, which are never called here but must compile. Antidunes.cpp gets these
// headers through BuildTimeloop.h, which this file does not include.
#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionModel.h"

#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/VTKOutput.h"

#include "ParticleKernels.h"
#include "SimulationObjects.h"
#include "SimulationParams.h"
#include "SetupStages_drop_test.h"

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

   // Minimal sim-state init instead of setupSimState(), which constructs a PIDController (no flow to drive here).
   simObjs->beginTimeStep = uint_t(0);
   simObjs->acceleration =
      std::make_shared< Vector3< real_t > >(real_t(0), real_t(0), -simParams->computed.gravity);

   setupParticles(simParams, simObjs);

   // DEM ---------------------------------------------------------------------------------------------------------
   ParticleKernels< ParticleAccessor_T > particleKernels(simParams, simObjs);

   // Minimal VTK: domain decomposition (once) and the floater outputs (periodic).
   // ----------------------------------------------------------------------------------------------------------------
   // Not inside WALBERLA_ROOT_SECTION: createVTKOutput_DomainDecomposition's write() is collective (every rank
   // contributes its blocks), so calling it on root only deadlocks. Same as BedGeneration.cpp.
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
      // Transposed body->world rotation, as in BuildTimeloop.h's floaterBox output (Tensor Glyph is column-major).
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

   // Cylinder floaters (type 4), same output as BuildTimeloop.h's "floaterCylinder" block.
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
      // (R * Rx(+90))^T = Rx(-90) * R^T: the glyph's axis is its local Y, the cylinder's its local Z (BuildTimeloop.h).
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

   // Sphere floaters (type 1) -- mirrors BuildTimeloop.h's own "floater" VTK block (radius from
   // the interaction radius), so that mixed-shape drop tests also render their spheres.
   std::function< void() > writeFloaterSphereVtk = []() {};
   if (!simParams->floatingSpheres.empty())
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
      auto sphereFloaterData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      sphereFloaterData->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
      sphereFloaterData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      sphereFloaterData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      sphereFloaterData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      sphereFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(1) && !isSet(pIt->getFlags(), GHOST);
      });
      auto sphereFloaterVTKOutput = vtk::createVTKOutput_PointData(
         sphereFloaterData, "floater", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      writeFloaterSphereVtk = vtk::writeFiles(sphereFloaterVTKOutput);
   }

   // Capsule floaters (type 5): like the cylinder block, rendered as same-radius/same-length cylinder glyphs because
   // Tensor Glyph has no capsule source (visual stand-in only; the physics uses the real Capsule shape).
   std::function< void() > writeFloaterCapsuleVtk = []() {};
   if (!simParams->floatingCapsules.empty())
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
      struct SelectCapsuleGlyphScale
      {
         using return_type = mesa_pd::Vec3;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto* shape = dynamic_cast< const mesa_pd::data::Capsule* >(p.getBaseShape().get());
            if (!shape) return mesa_pd::Vec3(real_t(0));
            return mesa_pd::Vec3(real_t(2) * shape->getRadius(), shape->getHeight(),
                                 real_t(2) * shape->getRadius());
         }
      };
      // Same (R * Rx(+90))^T construction as the cylinder block (glyph axis local Y, capsule axis local Z).
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

      auto capsuleFloaterData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      capsuleFloaterData->addOutput< SelectCapsuleGlyphScale >("glyphScale");
      capsuleFloaterData->addOutput< SelectParticleRotationMatrix >("rotation");
      capsuleFloaterData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      capsuleFloaterData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      capsuleFloaterData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      capsuleFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(5) && !isSet(pIt->getFlags(), GHOST);
      });
      auto capsuleFloaterVTKOutput = vtk::createVTKOutput_PointData(
         capsuleFloaterData, "floaterCapsule", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      writeFloaterCapsuleVtk = vtk::writeFiles(capsuleFloaterVTKOutput);
   }

   // Main loop -- no LBM, no LBM-DEM coupling, no PID. Just wave release + DEM sub-cycles + VTK.
   // ----------------------------------------------------------------------------------------------------------------
   WcTimingPool timingPool;

   WALBERLA_LOG_INFO_ON_ROOT("beginTimeStep: " << simObjs->beginTimeStep << " / " << simParams->physics.timesteps);
   for (uint_t t = simObjs->beginTimeStep; t != simParams->physics.timesteps; ++t)
   {
      *simObjs->currentMainStep = t;

      releaseSettledFloatingSpheres(simParams, simObjs, t);

      particleKernels.doDEMSubCycles(timingPool, simParams->particle.numSubCycles);

      writeFloaterBoxVtk();
      writeFloaterCylinderVtk();
      writeFloaterSphereVtk();
      writeFloaterCapsuleVtk();

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
