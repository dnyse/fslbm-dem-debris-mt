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
//! \file BuildTimeloop.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Jonas Plewinski <jonas.plewinski@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// Each function adds one logical group of sweeps/functions to the given timeloop.
// Nothing here creates fields or particles. Those are assumed to already exist in SimObjs
//======================================================================================================================

#include "core/SharedFunctor.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/StabilityChecker.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/free_surface/BlockStateDetectorSweep.h"
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/TotalMassComputer.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/dynamics/CellConversionSweep.h"
#include "lbm/free_surface/dynamics/ConversionFlagsResetSweep.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionSweep.h"
#include "lbm/free_surface/dynamics/PdfRefillingSweep.h"
#include "lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h"
#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/free_surface/surface_geometry/ObstacleFillLevelSweep.h"
#include "lbm/free_surface/surface_geometry/ObstacleNormalSweep.h"
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"
#include "lbm/free_surface/surface_geometry/Utility.h"

#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"

#include "mesa_pd/data/shape/Box.h"
#include "mesa_pd/data/shape/Cylinder.h"
#include "mesa_pd/data/shape/CylindricalBoundary.h"
#include "mesa_pd/vtk/OutputSelector.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/VTKOutput.h"

#include "Evaluators.h"
#include "Forcing.h"
#include "SimulationObjects.h"
#include "SimulationParams.h"

namespace walberla
{
namespace antidunes
{

using StateSweep = walberla::free_surface::BlockStateDetectorSweep< FlagField_T >;

// Empty sweep required for blocks that are not sweeped (i.e., excluded via blockstate selectors)
struct EmptySweep
{
   void operator()(IBlock*) {}
};

void addBeforeTimeStepVTKOutput(const std::shared_ptr< const SimulationParams >& simParams,
                                const std::shared_ptr< SimulationObjects >& simObjs,
                                const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addBeforeTimeStepVTKOutput");

   // LBM/freeSurface fields -----------------------------------------------
   walberla::free_surface::addVTKOutput< LatticeModel_T, AntidunesBoundaryHandling_T, PdfField_T, FlagField_T,
                                         ScalarField_T, VectorField_T >(
      simObjs->blockForest, *timeloop, simParams->cfgPtr, simObjs->flagInfo, simObjs->pdfFieldID, simObjs->flagFieldID,
      simObjs->fillFieldID, BlockDataID{}, simObjs->curvatureFieldID, simObjs->normalFieldID,
      simObjs->obstacleNormalFieldID);
   // triangle mesh of free surface ----------------------------------------
   // walberla::free_surface::SurfaceMeshWriter< ScalarField_T, FlagField_T >
   //    surfaceMeshWriter(simObjs->blockForest, simObjs->fillFieldID, simObjs->flagFieldID,
   //       walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, real_c(0), simParams->cfgPtr);
   // timeloop->addFuncBeforeTimeStep(surfaceMeshWriter, "Writer: surface mesh");

   // particleField --------------------------------------------------------
   // auto particleFieldVTKOutput = vtk::createVTKOutput_BlockData(simObjs->blockForest, "particle_field",
   //    simParams->vtk.writeFrequency, 0, false, simParams->vtk.baseFolder);
   // const AABB sliceAABB(real_t(0), real_c(simParams->domain.domainSize[1]) * real_t(0.5) - real_t(1), real_t(0),
   //                      real_c(simParams->domain.domainSize[0]),
   //                      real_c(simParams->domain.domainSize[1]) * real_t(0.5) + real_t(1),
   //                      real_c(simParams->domain.domainSize[2]));
   // particleFieldVTKOutput->addCellInclusionFilter(vtk::AABBCellFilter(sliceAABB));
   // particleFieldVTKOutput->addCellDataWriter(
   //    make_shared< field::VTKWriter< ParticleField_T > >(simObjs->particleFieldID, "particleField"));
   // particleFieldVTKOutput->setSamplingResolution(simParams->vtk.samplingResolution);
   // timeloop->addFuncBeforeTimeStep(vtk::writeFiles(particleFieldVTKOutput), "VTK: particle field data");

   // densityRatio (ρ_solid/ρ_fluid) for coloring in ParaView. Not stored directly on the particle —
   // recovered exactly from mass/volume, since updateMassAndInertia(density) at creation set
   // mass = getVolume() * density (SimulationParams' densityRatio IS that density, as ρ_fluid = 1 in
   // lattice units). Shape-generic: works unchanged for sphere/box/cylinder floaters alike.
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

   // particles (sediment spheres) -----------------------------------------
   auto particlePointData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
   particlePointData->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
   particlePointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
      using namespace walberla::mesa_pd::data::particle_flags;
      const auto& shape = pIt->getBaseShape();
      return shape && shape->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE && !isSet(pIt->getFlags(), GHOST) &&
             pIt->getType() == uint_t(0);
   });
   auto particleVTKOutput = vtk::createVTKOutput_PointData(
      particlePointData, "particles", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
      std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
   timeloop->addFuncBeforeTimeStep(vtk::writeFiles(particleVTKOutput), "VTK: sphere data");

   // floating spheres (type 1) -------------------------------------------
   if (!simParams->floatingSpheres.empty())
   {
      WALBERLA_LOG_INFO_ON_ROOT("AddOutputFloaters")
      auto floaterPointData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      floaterPointData->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
      floaterPointData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      floaterPointData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      floaterPointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(1) && !isSet(pIt->getFlags(), GHOST);
      });
      auto floaterVTKOutput = vtk::createVTKOutput_PointData(floaterPointData, "floater", simParams->vtk.writeFrequency,
                                                             simParams->vtk.baseFolder, std::string("simulation_step"),
                                                             false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(floaterVTKOutput), "VTK: floater data");
   }

   // floating boxes (type 3) ---------------------------------------------
   // One point per floating box, carrying edgeLength, full orientation, and velocity so the
   // actual ORIENTED solid can be rendered in ParaView (a freely-rotating box needs its
   // rotation, unlike a sphere). Recipe in ParaView:
   //   Filters > Glyph: Glyph Type = Box, Scale Array = edgeLength, Scale Mode = Scale by Components,
   //                    Orientation Mode = ... then Filters > Tensor Glyph using "rotation"
   //   (the 3x3 body→world matrix, 9 components) to apply the box's actual orientation; or color
   //   by "velocity". The unit ParaView box (edge 1) scaled by edgeLength gives the true size.
   if (!simParams->floatingBoxes.empty())
   {
      WALBERLA_LOG_INFO_ON_ROOT("AddOutputFloatingBoxes")

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
      // Full 3x3 orientation for a ParaView Tensor Glyph. getMatrix() is body→world
      // (Support.h: world = getMatrix()*body), but ParaView's Tensor Glyph applies the 9-tuple with
      // the opposite (column-major) ordering, which would render the INVERSE rotation — the box glyph
      // tilts the wrong way and mismatches the mapped obstacle hole. So we emit the TRANSPOSE here
      // (Rᵀ, i.e. column-major of R) so the glyph orientation matches the actual obstacle. This is a
      // visualization-only ordering; the buoyancy/torque math uses getMatrix() directly and is unaffected.
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
      boxFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(3) && !isSet(pIt->getFlags(), GHOST);
      });
      auto boxFloaterVTKOutput = vtk::createVTKOutput_PointData(
         boxFloaterData, "floaterBox", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(boxFloaterVTKOutput), "VTK: floating box data");
   }

   // floating cylinders / logs (type 4) -----------------------------------
   // One point per floating cylinder with glyphScale = (2r, L, 2r), the body→world rotation
   // tensor, and velocity. ParaView recipe: Glyph Type = Cylinder scaled by glyphScale (Scale by
   // Components), then Tensor Glyph using "rotation" (Extract Eigenvalues OFF) — NO manual
   // Transform/Rotate-X-90 step needed, unlike an earlier version of this export.
   //
   // Axis-convention fix: ParaView's built-in Cylinder glyph source has its OWN long axis along
   // its local Y, but mesa_pd's Cylinder shape defines its axis along its own local Z. Exporting
   // the plain body→world rotation R (transposed for ParaView's column-major Tensor Glyph, as the
   // box does) applies R as if the glyph's long axis were already Z, so the glyph's Y-axis (a
   // RADIAL direction of the real cylinder, not its axis) gets mapped to wherever the real axis
   // points — the glyph renders standing on the wrong axis (e.g. a horizontal log rendering as a
   // short vertical post). Confirmed visually: a manual ParaView Transform (Rotate X=90, plus a
   // compensating Translate since Transform rotates about the world origin, not the glyph's own
   // centre) was needed to fix a single frame — not viable across an animation where the glyph
   // moves every step. Fixed here instead by pre-composing R with Rx(+90°) (the fixed rotation
   // that maps "glyph-local Y is the axis" onto "body-local Z is the axis") BEFORE transposing:
   // exported = (R · Rx(+90°))ᵀ = Rx(-90°) · Rᵀ. Tensor Glyph applies this to the glyph's OWN
   // local geometry before translating to the point position, so — unlike the manual ParaView
   // Transform — this rotates correctly about the glyph's own centre at every timestep, no
   // per-frame compensation needed.
   if (!simParams->floatingCylinders.empty())
   {
      WALBERLA_LOG_INFO_ON_ROOT("AddOutputFloatingCylinders")

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
      struct SelectParticleRotationMatrix
      {
         using return_type = std::array< real_t, 9 >;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const mesa_pd::Mat3 R = p.getRotation().getMatrix();
            // (R * Rx(+90))^T = Rx(-90) * R^T — see comment above for the derivation.
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
      cylinderFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(4) && !isSet(pIt->getFlags(), GHOST);
      });
      auto cylinderFloaterVTKOutput = vtk::createVTKOutput_PointData(
         cylinderFloaterData, "floaterCylinder", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(cylinderFloaterVTKOutput), "VTK: floating cylinder data");
   }

   // Box pillar VTK — one point per pillar at box center with edgeLength attribute.
   // In ParaView: Filters > Glyph, glyph type Box, scale by edgeLength to render the full solid shape.
   if (!simParams->pillars.empty())
   {
      WALBERLA_LOG_INFO_ON_ROOT("OutputPillar")
      struct SelectBoxEdgeLength
      {
         using return_type = mesa_pd::Vec3;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto& shape = p.getBaseShape();
            if (!shape) return mesa_pd::Vec3(real_t(0));
            return static_cast< const mesa_pd::data::Box* >(shape.get())->getEdgeLength();
         }
      };
      auto boxPointData = std::make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      boxPointData->addOutput< SelectBoxEdgeLength >("edgeLength");
      boxPointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         const auto& shape = pIt->getBaseShape();
         return pIt->getType() == uint_t(2) && shape && shape->getShapeType() == mesa_pd::data::Box::SHAPE_TYPE;
      });
      auto boxVTKOutput = vtk::createVTKOutput_PointData(boxPointData, "pillar", simParams->vtk.writeFrequency,
                                                         simParams->vtk.baseFolder, std::string("simulation_step"),
                                                         false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(boxVTKOutput), "VTK: pillar data");
   }

   // Cylinder pillar VTK — one point per pillar at the axis centre with a glyphScale attribute.
   // In ParaView: Filters > Glyph, Glyph Type = Cylinder, Scale Mode = Scale by Components,
   // Scale Array = glyphScale, Scale Factor = 1. Then add Filters > Transform, Rotate X = 90
   // to align the cylinder axis (VTK glyph Y-axis) with the simulation Z-axis.
   // glyphScale = (2r, domainZ, 2r): the unit cylinder (r=0.5, h=1) scaled component-wise
   // gives the correct radius and full-domain height without any manual fudge factors.
   if (!simParams->cylinderPillars.empty())
   {
      struct SelectCylinderGlyphScale
      {
         using return_type = mesa_pd::Vec3;
         real_t domainZ;
         explicit SelectCylinderGlyphScale(real_t dz) : domainZ(dz) {}
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto* shape = dynamic_cast< const mesa_pd::data::CylindricalBoundary* >(p.getBaseShape().get());
            if (!shape) return mesa_pd::Vec3(real_t(0));
            const real_t r = shape->getRadius();
            return mesa_pd::Vec3(real_t(2) * r, domainZ, real_t(2) * r);
         }
      };

      auto cylinderPointData = std::make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      cylinderPointData->addOutput(
         "glyphScale",
         std::make_shared< mesa_pd::vtk::OutputSelector< SelectCylinderGlyphScale > >(
            SelectCylinderGlyphScale(real_c(simParams->domain.domainSize[2]))));
      cylinderPointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         const auto& shape = pIt->getBaseShape();
         return pIt->getType() == uint_t(2) && shape &&
                shape->getShapeType() == mesa_pd::data::CylindricalBoundary::SHAPE_TYPE;
      });
      auto cylinderVTKOutput = vtk::createVTKOutput_PointData(
         cylinderPointData, "cylinder_pillar", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(vtk::writeFiles(cylinderVTKOutput), "VTK: cylinder pillar data");
   }
}

void addSweepSurfaceGeometry(const std::shared_ptr< const SimulationParams >& simParams,
                             const std::shared_ptr< SimulationObjects >& simObjs,
                             const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addSweepSurfaceGeometry");

   // ObstacleNormalSweep + ObstacleFillLevelSweep — matches the sanctioned FiniteDifferenceMethod pipeline
   // (FiniteDifferenceMethod::addSweeps in src/lbm/free_surface/surface_geometry/CurvatureModel.impl.h).
   // They are required whenever an obstacle touches the free surface (moving particles, walls at the
   // contact line) REGARDLESS of wetting: ObstacleNormalSweep computes the wall normal (obstacleNormalField);
   // ObstacleFillLevelSweep reflects fill levels into obstacle cells so interface cells adjacent to the
   // obstacle have a COMPLETE fill-level gradient stencil. Without them obstacle cells have fill=0
   // → gradient ≈ 0 → zero normal → NaN cascade. ObstacleFillLevelSweep MUST run BEFORE the SmoothingSweep.
   walberla::free_surface::ObstacleNormalSweep< CommunicationStencil_T, FlagField_T, VectorField_T >
      obstacleNormalSweep(simObjs->obstacleNormalFieldID, simObjs->flagFieldID,
                          walberla::free_surface::flagIDs::interfaceFlagID,
                          walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs,
                          simObjs->flagInfo.getObstacleIDSet(), false, true, true);
   timeloop->add() << Sweep(obstacleNormalSweep, "Sweep: obstacle normal computation")
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->obstacleNormalFieldID),
                                    "Communication: after obstacle normal sweep");

   walberla::free_surface::ObstacleFillLevelSweep< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T >
      obstacleFillLevelSweep(simObjs->smoothFillFieldID, simObjs->fillFieldID, simObjs->flagFieldID,
                             simObjs->obstacleNormalFieldID,
                             walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs,
                             simObjs->flagInfo.getObstacleIDSet());
   timeloop->add() << Sweep(obstacleFillLevelSweep, "Sweep: obstacle fill level computation")
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->smoothFillFieldID),
                                    "Communication: after obstacle fill level sweep");

   // smooth fill level field for decreasing error in finite difference normal and curvature computation (see
   // dissertation of S. Bogner, 2017 (section 4.4.2.1))
   walberla::free_surface::SmoothingSweep< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T >
      smoothingSweep(simObjs->smoothFillFieldID, simObjs->fillFieldID, simObjs->flagFieldID,
                     walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, simObjs->flagInfo.getObstacleIDSet(),
                     simParams->physics.enableWetting);

   // IMPORTANT REMARK: SmoothingSweep must be executed on all blocks, because the algorithm works on all liquid,
   // interface and gas cells. This is necessary since the normals are not only computed in interface cells, but also
   // in the neighborhood of interface cells. Therefore, meaningful values for the fill levels of the second
   // neighbors of interface cells are also required in NormalSweep.
   timeloop->add() << Sweep(smoothingSweep, "Sweep: fill level smoothing")
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->smoothFillFieldID),
                                    "Communication: after smoothing sweep");

   // compute interface normals (using smoothed fill level field). With ObstacleNormalSweep +
   // ObstacleFillLevelSweep running above, obstacle cells now carry meaningful fill, so we match the
   // canonical FiniteDifferenceMethod flags (computeInInterfaceNeighbors, includeObstacleNeighbors,
   // modifyNearObstacles, computeInGhostLayer) = (true, true, true, true): obstacle-adjacent interface
   // cells get a complete gradient stencil instead of a degenerate (zero) normal at the contact line.
   walberla::free_surface::NormalSweep< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T > normalSweep(
      simObjs->normalFieldID, simObjs->smoothFillFieldID, simObjs->flagFieldID,
      walberla::free_surface::flagIDs::interfaceFlagID, walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs,
      simObjs->flagInfo.getObstacleIDSet(), true, true, true, true);

   timeloop->add() << Sweep(normalSweep, "Sweep: normal computation", StateSweep::fullFreeSurface)
                   << Sweep(EmptySweep(), "Empty sweep: normal")
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->normalFieldID),
                                    "Communication: after normal sweep");

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with no surface tension
   const bool computeCurvature = !realIsEqual(simParams->computed.surfaceTension, real_c(0), real_c(1e-14));
   if (computeCurvature)
   {
      // compute interface curvature using finite differences according to Brackbill et al.
      walberla::free_surface::CurvatureSweepFiniteDifferences< CommunicationStencil_T, FlagField_T, ScalarField_T,
                                                               VectorField_T >
         curvSweep(simObjs->curvatureFieldID, simObjs->normalFieldID, simObjs->obstacleNormalFieldID,
                   simObjs->flagFieldID, walberla::free_surface::flagIDs::interfaceFlagID,
                   walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, simObjs->flagInfo.getObstacleIDSet(),
                   false, real_c(0));

      timeloop->add() << Sweep(curvSweep, "Sweep: curvature computation (finite difference method)",
                               StateSweep::fullFreeSurface)
                      << Sweep(EmptySweep(), "Empty sweep: curvature")
                      << AfterFunction(Communication_T(simObjs->blockForest, simObjs->curvatureFieldID),
                                       "Communication: after curvature sweep");
   }
}

void addSweepSurfaceDynamics(const std::shared_ptr< const SimulationParams >& simParams,
                             const std::shared_ptr< SimulationObjects >& simObjs,
                             const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addSweepSurfaceDynamics");
   // add standard waLBerla boundary handling
   // only relavant for onlyGasAndBoundary blocks (i.e. not interface, not inflow, not liquid)
   timeloop->add() << Sweep(simObjs->antidunesBoundaryHandling->getBoundarySweep(), "Sweep: boundary handling",
                            Set< SUID >::emptySet(), StateSweep::onlyGasAndBoundary)
                   << Sweep(EmptySweep(), "Empty sweep: boundary handling", StateSweep::onlyGasAndBoundary);

   // comopute force density field (acceleration applied to fluid and interface cells)
   const ForceDensityCodegenSweep< PdfField_T, FlagField_T, VectorFieldFlattened_T, ScalarField_T >
      forceDensityCodegenSweep(simObjs->forceDensityFieldID, simObjs->pdfFieldID, simObjs->flagFieldID,
                               simObjs->fillFieldID, simObjs->flagInfo, simObjs->acceleration);
   timeloop->add() << Sweep(forceDensityCodegenSweep, "Sweep: force weighting", Set< SUID >::emptySet(),
                            StateSweep::onlyGasAndBoundary)
                   << Sweep(EmptySweep(), "Empty sweep: force weighting", StateSweep::onlyGasAndBoundary)
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->forceDensityFieldID),
                                    "Communication: after force weighting sweep");

   // sweep for
   // - reconstruction of PDFs in interface cells
   // - streaming of PDFs in interface cells (and liquid cells on the same block)
   // - advection of mass
   // - update bubble volumes
   // - marking interface cells for conversion
   const walberla::free_surface::StreamReconstructAdvectSweep<
      LatticeModel_T, typename AntidunesBoundaryHandling_T::BoundaryHandling_T, FlagField_T,
      typename AntidunesBoundaryHandling_T::FlagInfo_T, ScalarField_T, VectorField_T, true >
      streamReconstructAdvectSweep(
         simParams->computed.surfaceTension, simObjs->antidunesBoundaryHandling->getHandlingID(), simObjs->fillFieldID,
         simObjs->flagFieldID, simObjs->pdfFieldID, simObjs->normalFieldID, simObjs->curvatureFieldID,
         simObjs->flagInfo, simObjs->bubbleModel.get(), simParams->model.pdfReconstructionModel,
         simParams->model.useSimpleMassExchange, simParams->model.cellConversionThreshold,
         simParams->model.cellConversionForceThreshold);

   // sweep acts only on blocks with at least one interface cell (due to StateSweep::fullFreeSurface)
   timeloop->add()
      << Sweep(streamReconstructAdvectSweep, "Sweep: StreamReconstructAdvect", StateSweep::fullFreeSurface)
      << Sweep(EmptySweep(), "Empty sweep: StreamReconstructAdvect")
      // do not communicate PDFs here:
      // - stream on blocks with "StateSweep::fullFreeSurface" was performed here using post-collision PDFs
      // - stream on other blocks is performed below and should also use post-collision PDFs
      // => if PDFs were communicated here, the ghost layer of other blocks would have post-stream PDFs
      << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID, simObjs->flagFieldID),
                       "Communication: after StreamReconstructAdvect sweep")
      << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest, simObjs->fillFieldID),
                       "Second ghost layer update: after StreamReconstructAdvect sweep (fill level field)")
      << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(simObjs->blockForest, simObjs->flagFieldID),
                       "Second ghost layer update: after StreamReconstructAdvect sweep (flag field)");

   // Generated sweeps
   auto lbmSweepGenerated = typename LatticeModel_T::Sweep(simObjs->pdfFieldID);

   // temporary functor for calling the generated LBM collision
   class CollideSweep
   {
    public:
      CollideSweep(const typename LatticeModel_T::Sweep& sweep) : sweep_(sweep) {};

      void operator()(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0))
      { sweep_.collide(block, numberOfGhostLayersToInclude); }

    private:
      typename LatticeModel_T::Sweep sweep_;
   };

   // Collision for fullFreeSurface (because stream was done above) + StreamCollide for onlyLBM blocks
   timeloop->add() << Sweep(CollideSweep(lbmSweepGenerated), "Sweep: collision (generated)",
                            StateSweep::fullFreeSurface)
                   << Sweep(lbmSweepGenerated, "Sweep: streamCollide (generated)", StateSweep::onlyLBM)
                   << Sweep(EmptySweep(), "Empty sweep: streamCollide (generated)")
                   << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                                    "Communication: after streamCollide (generated)");

   // convert cells
   // - according to the flags from StreamReconstructAdvectSweep (interface -> gas/liquid)
   // - to ensure a closed layer of interface cells (gas/liquid -> interface)
   // - detect and register bubble merges/splits (bubble volumes are already updated in StreamReconstructAdvectSweep)
   // - convert cells and initialize PDFs near inflow boundaries
   const walberla::free_surface::CellConversionSweep< LatticeModel_T, AntidunesBoundaryHandling_T::BoundaryHandling_T,
                                                      ScalarField_T >
      cellConvSweep(simObjs->antidunesBoundaryHandling->getHandlingID(), simObjs->pdfFieldID, simObjs->flagInfo,
                    simObjs->bubbleModel.get());
   timeloop->add() << Sweep(cellConvSweep, "Sweep: cell conversion", StateSweep::fullFreeSurface)
                   << Sweep(EmptySweep(), "Empty sweep: cell conversion")
                   //<< AfterFunction(PdfCommunication_T(blockForest, pdfFieldID),
                   //
                   //                 "Communication: after cell conversion sweep (PDF field)")
                   // communicate the flag field also in corner directions
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->flagFieldID),
                                    "Communication: after cell conversion sweep (flag field)")
                   << AfterFunction(
                         blockforest::UpdateSecondGhostLayer< FlagField_T >(simObjs->blockForest, simObjs->flagFieldID),
                         "Second ghost layer update: after cell conversion sweep (flag field)");

   // reinitialize PDFs, i.e., refill cells that were converted from gas to interface
   // - when the flag "convertedFromGasToInterface" has been set (by CellConversionSweep)
   // - according to the method specified with pdfRefillingModel_
   const walberla::free_surface::EquilibriumRefillingSweep< LatticeModel_T, FlagField_T > equilibriumRefillingSweep(
      simObjs->pdfFieldID, simObjs->flagFieldID, simObjs->flagInfo, true);
   timeloop->add() << Sweep(equilibriumRefillingSweep, "Sweep: EquilibriumRefilling", StateSweep::fullFreeSurface)
                   << Sweep(EmptySweep(), "Empty sweep: EquilibriumRefilling")
                   << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                                    "Communication: after EquilibriumRefilling sweep");

   const walberla::free_surface::ExcessMassDistributionModel excessMassModel(
      simParams->model.excessMassDistributionModel);
   // distribute excess mass:
   // - excess mass: mass that is free after conversion from interface to gas/liquid cells
   // - update the bubble model
   // IMPORTANT REMARK: this sweep computes the mass via the density, i.e., the PDF field must be up-to-date and the
   // PdfRefillingSweep must have been performed
   if (excessMassModel.isEvenlyType())
   {
      const walberla::free_surface::ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T,
                                                                                ScalarField_T, VectorField_T >
         distributeMassSweep(excessMassModel, simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID,
                             simObjs->flagInfo);
      timeloop->add() << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                      << Sweep(EmptySweep(), "Empty sweep: distribute excess mass")
                      << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID),
                                       "Communication: after excess mass distribution sweep")
                      << AfterFunction(
                            blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest,
                                                                                 simObjs->fillFieldID),
                            "Second ghost layer update: after excess mass distribution sweep (fill level field)")
                      // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits are
                      // already detected and registered by CellConversionSweep
                      << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); }, "Sweep: bubble model update");
   }
   else
   {
      if (excessMassModel.isWeightedType())
      {
         const walberla::free_surface::ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T,
                                                                                     ScalarField_T, VectorField_T >
            distributeMassSweep(excessMassModel, simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID,
                                simObjs->flagInfo, simObjs->normalFieldID);
         timeloop->add() << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
                         << Sweep(EmptySweep(), "Empty sweep: distribute excess mass")
                         << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID),
                                          "Communication: after excess mass distribution sweep")
                         << AfterFunction(
                               blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest,
                                                                                    simObjs->fillFieldID),
                               "Second ghost layer update: after excess mass distribution sweep (fill level field)")
                         // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits
                         // are already detected and registered by CellConversionSweep
                         << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); },
                                          "Sweep: bubble model update");
      }
      else
      {
         if (excessMassModel.isEvenlyAllInterfaceFallbackLiquidType())
         {
            const walberla::free_surface::ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T,
                                                                                         ScalarField_T, VectorField_T >
               distributeMassSweep(excessMassModel, simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID,
                                   simObjs->flagInfo, simObjs->excessMassFieldID);
            timeloop->add()
               // perform this sweep also on "onlyLBM" blocks because liquid cells also exchange excess mass here
               << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
               << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::onlyLBM)
               << Sweep(EmptySweep(), "Empty sweep: distribute excess mass")
               << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID, simObjs->excessMassFieldID),
                                "Communication: after excess mass distribution sweep")
               << AfterFunction(
                     blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest, simObjs->fillFieldID),
                     "Second ghost layer update: after excess mass distribution sweep (fill level field)")
               // update bubble model, i.e., perform registered bubble merges/splits; bubble
               // merges/splits are already detected and registered by CellConversionSweep
               << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); }, "Sweep: bubble model update");
         }
      }
   }

   // reset all flags that signal cell conversions (except "keepInterfaceForWettingFlag")
   walberla::free_surface::ConversionFlagsResetSweep< FlagField_T > resetConversionFlagsSweep(simObjs->flagFieldID,
                                                                                              simObjs->flagInfo);
   timeloop->add() << Sweep(resetConversionFlagsSweep, "Sweep: conversion flag reset", StateSweep::fullFreeSurface)
                   << Sweep(EmptySweep(), "Empty sweep: conversion flag reset")
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->flagFieldID),
                                    "Communication: after excess mass distribution sweep")
                   << AfterFunction(
                         blockforest::UpdateSecondGhostLayer< FlagField_T >(simObjs->blockForest, simObjs->flagFieldID),
                         "Second ghost layer update: after excess mass distribution sweep (flag field)");

   // update block states
   auto blockStateUpdate = StateSweep(simObjs->blockForest, simObjs->flagInfo, simObjs->flagFieldID);
   timeloop->add() << Sweep(blockStateUpdate, "Sweep: block state update");
}

void addSweepFluidParticleCoupling(const std::shared_ptr< const SimulationParams >& simParams,
                                   const std::shared_ptr< SimulationObjects >& simObjs,
                                   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addSweepFluidParticleCoupling");

   // TODO: maybe put this FlagUID in the antidunesBoundaryHandling
   const FlagUID FormerMO_Flag("former moving obstacle");

   // Map every sphere (bed type 0 + floating type 1) AND every floating box (type 3) / floating
   // cylinder (type 4) as a moving obstacle, so each is a full two-way-coupled floater — blocks the
   // fluid, receives a MEM hydrodynamic force via SimpleBB, and has its vacated cells reconstructed.
   // Fixed pillars (type 2) are deliberately excluded here; they are mapped separately as static
   // NoSlip/FreeSlip (see SetupStages). Must stay consistent with `floaterSelector` in SetupStages.
   // Dormant (startInactive, not yet at settlingSteps) floaters are excluded too -- this is the
   // per-step equivalent of SetupStages::mapParticlesOntoCells's init-time dormant guard; both must
   // agree or a dormant floater would appear in the fluid the very first step after initialization.
   auto dormantUIDs = simObjs->dormantFloaterUIDs;
   auto movingObstacleSelector = [dormantUIDs](const size_t particleIdx, const auto& ac) {
      const bool isSphere           = ac.getBaseShape(particleIdx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE;
      const bool isFloatingBox      = ac.getType(particleIdx) == uint_t(3);
      const bool isFloatingCylinder = ac.getType(particleIdx) == uint_t(4);
      if (!(isSphere || isFloatingBox || isFloatingCylinder)) return false;
      return dormantUIDs->count(ac.getUid(particleIdx)) == 0;
   };

   // sweep for mapping particles into the LBM cells
   const bool strictlyConserveMomentum = false;
   timeloop->add() << Sweep(
      lbm_mesapd_coupling::makeMovingParticleMapping< PdfField_T, AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         simObjs->blockForest, simObjs->pdfFieldID, simObjs->antidunesBoundaryHandling->getHandlingID(),
         simObjs->particleFieldID, simObjs->particleAccessor, AntidunesBoundaryHandling_T::movingObstacleFlagID,
         FormerMO_Flag,
         movingObstacleSelector, // all spheres + floating boxes (not fixed pillars)
         strictlyConserveMomentum),
      "Sweep (coupling): MovingParticleMapping");

   // sweep for restoring PDFs in cells previously occupied by particles
   //
   // FIX: switch from GradsMomentApproximationReconstructor to EquilibriumReconstructor.
   //
   // GradsMomentApproximation constructs f_i = f_i^eq(ρ_avg, u_particle) + Pneq_correction.
   // The Pneq term involves velocity gradients computed from neighboring FLUID cells.
   // At the contact circle (z ≈ absoluteLiquidHeight, just inside the sphere surface), the
   // neighbors include:
   //   - interface cells with fill-drift-corrupted PDFs (contact-line mass imbalance)
   //   - gas cells (fill=0, PDF=0 — above the free surface)
   //   - obstacle cells (still inside the sphere)
   // The Pneq correction applied on top of an already-low average density (only ~⅔ of the
   // D3Q27 stencil has valid fluid neighbours) pushes individual f_i negative.  After one
   // FSLBM streaming step those negative values propagate into adjacent liquid cells and the
   // Cumulant collision (τ ≈ 0.501) diverges → NaN.
   //
   // EquilibriumReconstructor sets ALL 27 PDFs to f_i^eq(ρ_avg, u_particle) ≥ 0.
   // Streaming from gas neighbours (PDF=0) then reduces ρ, but keeps every f_i ≥ 0, so the
   // collision cannot produce NaN from a negative-density state.
   //
   // Previous (now-commented) GradsMoment version, kept for reference:
   // auto gradReconstructor =
   //    lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor< AntidunesBoundaryHandling_T::BoundaryHandling_T >(
   //       simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID(), simParams->computed.relaxationRate,
   //       false,  // recomputeTargetDensity
   //       true    // useCentralDifferences
   //    );

   // Bogner & Rüde (2013) state-selection + EquilibriumReconstructor.
   //
   // Each FormerMO cell is classified by its D3Q27 non-obstacle neighbourhood B(x):
   //   no liquid            →  GAS        (φ = 0)
   //   no gas, no unknown   →  LIQUID     (φ = 1, clearly submerged)
   //   everything else      →  INTERFACE  (φ = mean fill of interface cells in B(x))
   // The φ seeding (§5.2) is done in the manager's pass 2 via the fill field.
   //
   // Held as a named shared_ptr so the AfterFunction below can flush its per-branch FormerMO→state
   // conversion counters (see logAndResetConversionCounts / contact_ring_crash_analysis.md §5.3).
   // useBognerReconstruction (prm: ModelParameters.useBognerReconstruction):
   //   true  → per-cell gas/liquid/interface state selection (+ φ seeding) as above.
   //   false → every FormerMO cell is reconstructed uniformly as INTERFACE (target = fluid below),
   //           with no φ seeding — the legacy all-interface path. The manager type is unchanged
   //           (the toggle is a runtime flag), so the conversion-count AfterFunction still compiles;
   //           with Bogner off the counters simply stay zero.
   auto bognerReconstructionManager =
      lbm_mesapd_coupling::makePdfReconstructionManagerBogner< PdfField_T,
                                                               AntidunesBoundaryHandling_T::BoundaryHandling_T,
                                                               ScalarField_T >(
         simObjs->blockForest, simObjs->pdfFieldID, simObjs->antidunesBoundaryHandling->getHandlingID(),
         simObjs->particleFieldID, simObjs->particleAccessor, FormerMO_Flag,
         walberla::free_surface::flagIDs::interfaceFlagID,
         strictlyConserveMomentum,
         walberla::free_surface::flagIDs::liquidFlagID,
         walberla::free_surface::flagIDs::gasFlagID,
         walberla::free_surface::flagIDs::interfaceFlagID,
         AntidunesBoundaryHandling_T::movingObstacleFlagID,
         simObjs->fillFieldID,
         simParams->model.useBognerReconstruction,
         /* optimizeForSmallObstacleFraction */ false,
         simParams->model.useHydrostaticReconstruction,
         (*simObjs->acceleration)[2],                    // a_z (= -gravity), same value HydrostaticDensityFtor uses
         simParams->staticReferenceHeight(),             // flat static reference H (same as SimpleBB)
         simParams->model.usePass3Reresolution);

   const uint_t convLogFreq = simParams->evaluation.evaluationFrequency;

   timeloop->add()
      << Sweep(makeSharedSweep(bognerReconstructionManager),
               "Sweep (coupling): PDF Reconstruction (moving particles) [Bogner2way+EquilibriumReconstructor]")
      << AfterFunction(
            [bognerReconstructionManager, convLogFreq, step = uint_t(0)]() mutable {
               ++step;
               if (convLogFreq > uint_t(0) && step % convLogFreq == uint_t(0))
                  bognerReconstructionManager->logAndResetConversionCounts(step);
            },
            "AfterFunction (coupling): log FormerMO→state conversion counts")
      << AfterFunction(
            Communication_T(simObjs->blockForest, simObjs->flagFieldID, simObjs->particleFieldID),
            "Communication (coupling): after PDF reconstruction (flagField + particleField)") // unsure if necessary but
                                                                                              // added for consistency
      << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                       "Communication (coupling): after PDF reconstruction (pdfField)");
}

void addAfterTimeStepPerformanceLogger(const std::shared_ptr< const SimulationParams >& simParams,
                                       const std::shared_ptr< SimulationObjects >& simObjs,
                                       const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepPerformanceLogger");

   lbm::PerformanceLogger< FlagField_T > performanceLogger(simObjs->blockForest, simObjs->flagFieldID,
                                                           walberla::free_surface::flagIDs::liquidInterfaceFlagIDs,
                                                           simParams->evaluation.performanceLogFrequency);
   timeloop->addFuncAfterTimeStep(performanceLogger, "Evaluator: performance logging");
}

void addAfterTimeStepStabilityChecker(const std::shared_ptr< const SimulationParams >& simParams,
                                      const std::shared_ptr< SimulationObjects >& simObjs,
                                      const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepStabilityChecker");

   // logs data and aborts if it detects non-finite values in liquid or interfce cells
   timeloop->addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                     simParams->cfgPtr, simObjs->blockForest, simObjs->pdfFieldID, simObjs->flagFieldID,
                                     walberla::free_surface::flagIDs::liquidInterfaceFlagIDs)),
                                  "Evaluator: LBM stability check");
}

void addAfterTimeStepDataEvaluators(const std::shared_ptr< const SimulationParams >& simParams,
                                    const std::shared_ptr< SimulationObjects >& simObjs,
                                    const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepDataEvaluators");

   // data slice evaluator --------------------------------------------------------------------------------------
   // writes data to file (solidVolumeFraction + fill + velocityX) averaged over X-Z slices in the domain
   AverageDataSliceEvaluator< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
      averageDataSliceEvaluator(simObjs->blockForest, simObjs->flagFieldID, simObjs->fillFieldID, simObjs->pdfFieldID,
                                simParams->evaluation.baseFolderName, simObjs->dataCols,
                                simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(averageDataSliceEvaluator, "Evaluator: AverageDataSlice");

   // total mass evaluator --------------------------------------------------------------------------------------
   walberla::free_surface::TotalMassComputer< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
      totalFluidMassEvaluator(simObjs->blockForest, simObjs->antidunesBoundaryHandling, simObjs->pdfFieldID,
                              simObjs->fillFieldID, simParams->evaluation.evaluationFrequency, simObjs->totalFluidMass);
   timeloop->addFuncAfterTimeStep(totalFluidMassEvaluator, "Evaluator: TotalFluidMass");

   // bed load transport evaluator ------------------------------------------------------------------------------
   const real_t bedLoadTransportAveragingFactor =
      real_c(1) / real_c(simParams->domain.domainSize[0] * simParams->domain.domainSize[1]);
   BedloadTransportEvaluator< ParticleAccessor_T > bedloadTransportEvaluator(
      simObjs->particleAccessor, bedLoadTransportAveragingFactor, simObjs->numParticles, simObjs->dataCols,
      simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(bedloadTransportEvaluator, "Evaluator: BedloadTransport");

   // total hydro force evaluator ------------------------------------------------------------------------------
   TotalHydrodynamicForceEvaluator< ParticleAccessor_T > totalHydrodynamicForceEvaluator(
      simObjs->particleAccessor, simObjs->dataCols, simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(totalHydrodynamicForceEvaluator, "Evaluator: TotalHydrodynamicForce");
}

void addAfterTimeStepForcingAdjuster(const std::shared_ptr< const SimulationParams >& simParams,
                                     const std::shared_ptr< SimulationObjects >& simObjs,
                                     const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepForcingAdjuster");

   // compute fluid's mean velocity -----------------------------------------------------------------------------
   const real_t velocityAveragingFactor =
      real_c(1) /
      (simParams->computed.liquidHeight * real_c(simParams->domain.domainSize[0] * simParams->domain.domainSize[1]));
   MeanVelocityComputer< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T > meanVelocityComputer(
      simObjs->blockForest, simObjs->antidunesBoundaryHandling, simObjs->pdfFieldID, simObjs->meanVelocity,
      velocityAveragingFactor, simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(meanVelocityComputer, "Evaluator: MeanVelocityComputer");

   // update acceleration ---------------------------------------------------------------------------------------
   // using the computed mean velocity, adjust the acceleration (in X dir) such that the target mean velocity is
   // maintained
   ForcingAdjuster forcingAdjuster(simObjs->meanVelocity, simObjs->acceleration, simObjs->pidController,
                                   simParams->evaluation.evaluationFrequency);

   timeloop->addFuncAfterTimeStep(forcingAdjuster, "Evaluator: ForcingAdjuster");
}

} // namespace antidunes
} // namespace walberla