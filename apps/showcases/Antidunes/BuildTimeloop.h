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
// Each add* function adds one logical group of sweeps/functions to the given timeloop; the other functions are
// diagnostics and the wallOuterGhostLayer flag fix-up called from Antidunes.cpp.
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

#include <sstream>

namespace walberla
{
namespace antidunes
{

// Forward declarations of diagnostics defined in SetupStages.h, which Antidunes.cpp includes after this file.
// trackWatchCells: logs the state of a fixed list of watch cells.
void trackWatchCells(const std::shared_ptr< SimulationObjects >& simObjs, const char* const label, const uint_t t,
                     const std::vector< Vector3< real_t > >& watchCoords);

void logFluidBoxHistory(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t,
                        const std::string& baseFolder,
                        const Vector3< real_t >& boxMin, const Vector3< real_t >& boxMax);

void logDensityColumn(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t,
                      const std::string& baseFolder, const real_t columnX, const real_t columnY);

// Fill-integrated water height on a grid of columns (ModelParameters.logWaterHeight).
void logWaterHeights(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t, const std::string& baseFolder,
                     const uint_t gridX, const uint_t gridY);

// Scans every liquid/interface cell, logs non-finite ones with a 3x3x3 neighbourhood dump, and issues the precursor
// warnings configured in PhysicsParameters.
bool checkPdfNaN(const std::shared_ptr< SimulationObjects >& simObjs, const char* const location, const uint_t t,
                 real_t lowRhoPrecursorThreshold, real_t highVelocityPrecursorThreshold,
                 real_t populationPrecursorThreshold, real_t densityPrecursorLow,
                 real_t densityPrecursorHigh, bool densityPrecursorLiquidOnly);

// Resets obstacle cells with non-finite or non-positive PDF sums (ModelParameters.repairDegenerateObstaclePdfs).
uint_t repairDegenerateObstacleCells(const std::shared_ptr< SimulationObjects >& simObjs, const char* const location, const uint_t t);

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

   // The particle/floater point-data writers below also force a write at every step of VTKParams::snapshotSteps
   // (VTKOutput::forceWriteNextStep), in addition to their own writeFrequency. No-op if snapshotsEnabled is false.
   auto addSnapshotAwareVTKWrite = [simObjs, enabled = simParams->vtk.snapshotsEnabled,
                                     steps = simParams->vtk.snapshotSteps](const std::shared_ptr< vtk::VTKOutput >& output) {
      return [output, simObjs, enabled, steps]() {
         if (enabled && steps.count(*simObjs->currentMainStep)) output->forceWriteNextStep();
         vtk::writeFiles(output)();
      };
   };

   // LBM/freeSurface fields -----------------------------------------------
   walberla::free_surface::addVTKOutput< LatticeModel_T, AntidunesBoundaryHandling_T, PdfField_T, FlagField_T,
                                         ScalarField_T, VectorField_T >(
      simObjs->blockForest, *timeloop, simParams->cfgPtr, simObjs->flagInfo, simObjs->pdfFieldID, simObjs->flagFieldID,
      simObjs->fillFieldID, BlockDataID{}, simObjs->curvatureFieldID, simObjs->normalFieldID,
      simObjs->obstacleNormalFieldID,
      simParams->vtk.denseEnabled ? simObjs->currentMainStep : nullptr,
      simParams->vtk.denseFromStep, simParams->vtk.denseToStep,
      simParams->vtk.snapshotsEnabled ? simObjs->currentMainStep : nullptr,
      simParams->vtk.snapshotSteps);
   if (simParams->vtk.denseEnabled)
      WALBERLA_LOG_INFO_ON_ROOT("[VTK] fluid_field_dense active for steps ["
         << simParams->vtk.denseFromStep << ", " << simParams->vtk.denseToStep << ")");
   if (simParams->vtk.snapshotsEnabled)
   {
      std::ostringstream steps;
      for (uint_t s : simParams->vtk.snapshotSteps)
         steps << s << " ";
      WALBERLA_LOG_INFO_ON_ROOT("[VTK] fluid_field_snapshots active at steps [ " << steps.str() << "]");
   }
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
   timeloop->addFuncBeforeTimeStep(addSnapshotAwareVTKWrite(particleVTKOutput), "VTK: sphere data");

   // floating spheres (type 1) -------------------------------------------
   if (!simParams->floatingSpheres.empty())
   {
      WALBERLA_LOG_INFO_ON_ROOT("AddOutputFloaters")
      auto floaterPointData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      floaterPointData->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
      floaterPointData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      floaterPointData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      // uid: lets a floater picked in ParaView be matched against log/CSV output keyed by uid.
      floaterPointData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      floaterPointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(1) && !isSet(pIt->getFlags(), GHOST);
      });
      auto floaterVTKOutput = vtk::createVTKOutput_PointData(floaterPointData, "floater", simParams->vtk.writeFrequency,
                                                             simParams->vtk.baseFolder, std::string("simulation_step"),
                                                             false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(addSnapshotAwareVTKWrite(floaterVTKOutput), "VTK: floater data");
   }

   // floating boxes (type 3) ---------------------------------------------
   // One point per floating box with edgeLength, rotation, a combined "sizeOrientation" tensor, velocity, densityRatio
   // and uid, so oriented boxes of different sizes can be rendered in ParaView. Tensor Glyph takes both size and
   // orientation from one tensor (a pure rotation gives unit-size glyphs), hence "sizeOrientation": the rotation with
   // each column scaled by that axis's edge length.
   // Recipe: Filters > Tensor Glyph, Glyph Type = Box, Extract Eigenvalues = OFF, Tensors = "sizeOrientation",
   // Scale Factor = 1; colour by "velocity" or "densityRatio". "edgeLength" and "rotation" serve a plain Glyph +
   // Scale-by-Components workflow (orientation vector only, no roll).
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
      // Full 3x3 orientation for a ParaView Tensor Glyph. getMatrix() is body->world, but Tensor Glyph reads the
      // 9-tuple column-major, so the transpose R^T is emitted to get the rotation (not its inverse). Visualisation
      // only.
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
      // Same R^T convention, each column (the box's local axis in world space) scaled by that axis's edge length.
      struct SelectBoxSizeOrientationTensor
      {
         using return_type = std::array< real_t, 9 >;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto& shape = p.getBaseShape();
            mesa_pd::Vec3 edgeLength(real_t(0));
            if (shape && shape->getShapeType() == mesa_pd::data::Box::SHAPE_TYPE)
               edgeLength = static_cast< const mesa_pd::data::Box* >(shape.get())->getEdgeLength();
            const mesa_pd::Mat3 R = p.getRotation().getMatrix();
            return { R(0, 0) * edgeLength[0], R(1, 0) * edgeLength[0], R(2, 0) * edgeLength[0],
                     R(0, 1) * edgeLength[1], R(1, 1) * edgeLength[1], R(2, 1) * edgeLength[1],
                     R(0, 2) * edgeLength[2], R(1, 2) * edgeLength[2], R(2, 2) * edgeLength[2] };
         }
      };

      auto boxFloaterData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      boxFloaterData->addOutput< SelectBoxEdgeLength >("edgeLength");
      boxFloaterData->addOutput< SelectParticleRotationMatrix >("rotation");
      boxFloaterData->addOutput< SelectBoxSizeOrientationTensor >("sizeOrientation");
      boxFloaterData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      boxFloaterData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      // uid: see floaterPointData's own comment above for why.
      boxFloaterData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      boxFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(3) && !isSet(pIt->getFlags(), GHOST);
      });
      auto boxFloaterVTKOutput = vtk::createVTKOutput_PointData(
         boxFloaterData, "floaterBox", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(addSnapshotAwareVTKWrite(boxFloaterVTKOutput), "VTK: floating box data");
   }

   // floating cylinders / logs (type 4) -----------------------------------
   // One point per floating cylinder with glyphScale = (2r, L, 2r), rotation, "sizeOrientation", velocity,
   // densityRatio and uid. Recipe: Filters > Tensor Glyph, Glyph Type = Cylinder (radius 0.5, height 1), Extract
   // Eigenvalues = OFF, Tensors = "sizeOrientation", Scale Factor = 1 (see the box block for why a pure rotation is
   // not enough). ParaView's cylinder glyph has its axis along local Y, mesa_pd's Cylinder along local Z, so R is
   // pre-composed with Rx(+90) before transposing: exported = (R Rx(+90))^T = Rx(-90) R^T.
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
      // Same (R Rx(+90))^T convention, each emitted group scaled by the matching glyphScale component (2r, L, 2r).
      struct SelectCylinderSizeOrientationTensor
      {
         using return_type = std::array< real_t, 9 >;
         return_type operator()(const mesa_pd::data::Particle& p) const
         {
            const auto* shape = dynamic_cast< const mesa_pd::data::Cylinder* >(p.getBaseShape().get());
            mesa_pd::Vec3 glyphScale(real_t(0));
            if (shape)
               glyphScale = mesa_pd::Vec3(real_t(2) * shape->getRadius(), shape->getHeight(),
                                          real_t(2) * shape->getRadius());
            const mesa_pd::Mat3 R = p.getRotation().getMatrix();
            return { R(0, 0) * glyphScale[0],  R(1, 0) * glyphScale[0],  R(2, 0) * glyphScale[0],
                     R(0, 2) * glyphScale[1],  R(1, 2) * glyphScale[1],  R(2, 2) * glyphScale[1],
                    -R(0, 1) * glyphScale[2], -R(1, 1) * glyphScale[2], -R(2, 1) * glyphScale[2] };
         }
      };

      auto cylinderFloaterData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
      cylinderFloaterData->addOutput< SelectCylinderGlyphScale >("glyphScale");
      cylinderFloaterData->addOutput< SelectParticleRotationMatrix >("rotation");
      cylinderFloaterData->addOutput< SelectCylinderSizeOrientationTensor >("sizeOrientation");
      cylinderFloaterData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      cylinderFloaterData->addOutput< SelectFloaterDensityRatio >("densityRatio");
      // uid: see floaterPointData's own comment above for why.
      cylinderFloaterData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
      cylinderFloaterData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getType() == uint_t(4) && !isSet(pIt->getFlags(), GHOST);
      });
      auto cylinderFloaterVTKOutput = vtk::createVTKOutput_PointData(
         cylinderFloaterData, "floaterCylinder", simParams->vtk.writeFrequency, simParams->vtk.baseFolder,
         std::string("simulation_step"), false, false, true, false, simObjs->beginTimeStep);
      timeloop->addFuncBeforeTimeStep(addSnapshotAwareVTKWrite(cylinderFloaterVTKOutput), "VTK: floating cylinder data");
   }

   // Pillars get no glyph output: their cells are visible through the fluid field's "mapped_flag" writer
   // (FlagFieldMapping, VtkWriter.h).
}

// Called after flag initialisation (Antidunes.cpp) and by diagnoseNaNNormals every 1000 steps: counts the cells of
// the outer ghost layer (layer 2) beyond non-periodic domain borders that carry liquid/interface/gas flags. A Border
// wall with walldistance -1 occupies only ghost layer 1; layer 2 behind it gets the domain flag from
// fillWithDomain(nrOfGhostLayers) in mapParticlesOntoCells, is not reset by initFlagsFromFillLevels (interior only)
// and is never communicated, so ObstacleFillLevelSweep (which includes ghost layer 1) treats it as fluid next to the
// wall. Logs [GHOST2-FLUID label]: number of such cells, of those with non-finite fill, and their fill range.
inline void checkOuterGhostLayerFlags(const std::shared_ptr< SimulationObjects >& simObjs, const std::string& label)
{
   uint_t nCells = 0, nFluid = 0, nNonFinite = 0;
   real_t fMin = std::numeric_limits< real_t >::max(), fMax = -std::numeric_limits< real_t >::max();
   const auto& bf     = *simObjs->blockForest;
   const auto& domain = bf.getDomain();
   const bool  periodic[3] = { bf.isXPeriodic(), bf.isYPeriodic(), bf.isZPeriodic() };
   for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
   {
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      if (flagField->nrOfGhostLayers() < uint_t(2) || fillField->nrOfGhostLayers() < uint_t(2)) continue;
      const auto fluidMask = flag_t(flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID) |
                                    flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID) |
                                    flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID));
      const auto& aabb = blockIt->getAABB();
      const cell_idx_t size[3] = { cell_idx_c(flagField->xSize()), cell_idx_c(flagField->ySize()), cell_idx_c(flagField->zSize()) };
      for (uint_t d = 0; d < 3; ++d)
      {
         if (periodic[d]) continue;
         for (int side = 0; side < 2; ++side)
         {
            const bool atBorder = side == 0 ? realIsEqual(aabb.min()[d], domain.min()[d]) : realIsEqual(aabb.max()[d], domain.max()[d]);
            if (!atBorder) continue;
            const cell_idx_t layer = side == 0 ? cell_idx_c(-2) : size[d] + cell_idx_c(1);
            CellInterval ci(-2, -2, -2, size[0] + 1, size[1] + 1, size[2] + 1);
            ci.min()[d] = layer;
            ci.max()[d] = layer;
            for (auto c = ci.begin(); c != ci.end(); ++c)
            {
               ++nCells;
               if ((flagField->get(*c) & fluidMask) == flag_t(0)) continue;
               ++nFluid;
               const real_t f = fillField->get(*c);
               if (!std::isfinite(f)) { ++nNonFinite; continue; }
               fMin = std::min(fMin, f);
               fMax = std::max(fMax, f);
            }
         }
      }
   }
   std::vector< uint_t > c{ nCells, nFluid, nNonFinite };
   mpi::reduceInplace(c, mpi::SUM);
   mpi::reduceInplace(fMin, mpi::MIN);
   mpi::reduceInplace(fMax, mpi::MAX);
   WALBERLA_LOG_INFO_ON_ROOT("[GHOST2-FLUID " << label << "] outer ghost layer beyond non-periodic borders: " << c[0]
                             << " cells, " << c[1] << " with liquid/interface/gas flags, " << c[2]
                             << " of those with non-finite fill; fill range of the flagged ones ["
                             << (c[1] > c[2] ? fMin : real_t(0)) << ", " << (c[1] > c[2] ? fMax : real_t(0)) << "]");
}

// ModelParameters.wallOuterGhostLayer (see SimulationParams.h): makes the outer ghost layer beyond
// every non-periodic domain border a wall, copying the NoSlip/FreeSlip flag of the layer-1 cell between it and the
// domain, with fill 1 as initFillLevelsInBoundaries sets for walls. Other boundary types in layer 1 (inflow, pressure,
// outlet) are left alone. Must run after initFlagsFromFillLevel. Returns the number of cells changed on this rank.
inline uint_t wallOuterGhostLayer(const std::shared_ptr< SimulationObjects >& simObjs)
{
   uint_t changed = 0;
   const auto& bf     = *simObjs->blockForest;
   const auto& domain = bf.getDomain();
   const bool  periodic[3] = { bf.isXPeriodic(), bf.isYPeriodic(), bf.isZPeriodic() };
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      auto* handling = blockIt->getData< AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         simObjs->antidunesBoundaryHandling->getHandlingID());
      auto* flagField = handling->getFlagField();
      auto* fillField = blockIt->getData< ScalarField_T >(simObjs->fillFieldID);
      if (flagField->nrOfGhostLayers() < uint_t(2) || fillField->nrOfGhostLayers() < uint_t(2)) continue;
      const flag_t noSlip   = flagField->getFlag(AntidunesBoundaryHandling_T::noSlipFlagID);
      const flag_t freeSlip = flagField->getFlag(AntidunesBoundaryHandling_T::freeSlipFlagID);
      const auto&  aabb     = blockIt->getAABB();
      const cell_idx_t size[3] = { cell_idx_c(flagField->xSize()), cell_idx_c(flagField->ySize()), cell_idx_c(flagField->zSize()) };
      for (uint_t d = 0; d < 3; ++d)
      {
         if (periodic[d]) continue;
         for (int side = 0; side < 2; ++side)
         {
            const bool atBorder = side == 0 ? realIsEqual(aabb.min()[d], domain.min()[d]) : realIsEqual(aabb.max()[d], domain.max()[d]);
            if (!atBorder) continue;
            const cell_idx_t outer = side == 0 ? cell_idx_c(-2) : size[d] + cell_idx_c(1);
            const cell_idx_t inner = side == 0 ? cell_idx_c(-1) : size[d];
            CellInterval ci(-2, -2, -2, size[0] + 1, size[1] + 1, size[2] + 1);
            ci.min()[d] = outer;
            ci.max()[d] = outer;
            for (auto c = ci.begin(); c != ci.end(); ++c)
            {
               Cell w = *c;
               w[d]   = inner;
               const flag_t wf = flagField->get(w);
               const flag_t f  = (wf & noSlip) ? noSlip : ((wf & freeSlip) ? freeSlip : flag_t(0));
               if (f == flag_t(0) || (flagField->get(*c) & f)) continue;
               handling->forceBoundary(f, c->x(), c->y(), c->z());
               fillField->get(*c) = real_c(1);
               ++changed;
            }
         }
      }
   }
   return changed;
}

// ModelParameters.diagnoseNaNNormals: after the normal sweep and its communication, scans inner cells and the first
// ghost layer for non-finite interface normals, smoothed fill levels and obstacle normals and logs [NAN-NORMAL] /
// [NAN-SMOOTHFILL] / [NAN-OBSTACLENORMAL] with the inputs of the normal computation (flags, fill, smoothed fill,
// obstacle neighbours' smoothed fill and normal, neighbours with non-finite smoothed fill, L/I/G neighbour counts).
// Throttled per kind and rank (first 30, then every 20000th); a global count per kind every 1000 steps.
inline void diagnoseNaNNormals(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t)
{
   static uint_t logged[3]  = { 0, 0, 0 };
   static uint_t counted[3] = { 0, 0, 0 };
   static const char* names[3] = { "NAN-NORMAL", "NAN-SMOOTHFILL", "NAN-OBSTACLENORMAL" };
   auto finite3 = [](const Vector3< real_t >& v) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); };
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* flagField   = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField   = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto* smoothField = blockIt->getData< const ScalarField_T >(simObjs->smoothFillFieldID);
      const auto* normalField = blockIt->getData< const VectorField_T >(simObjs->normalFieldID);
      const auto* obstNField  = blockIt->getData< const VectorField_T >(simObjs->obstacleNormalFieldID);
      const auto  liquid      = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto  interface   = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const auto  gas         = flagField->getFlag(walberla::free_surface::flagIDs::gasFlagID);
      const auto  obstacle    = flagField->getMask(simObjs->flagInfo.getObstacleIDSet());
      const auto& aabb        = blockIt->getAABB();
      const real_t dx         = aabb.xSize() / real_c(flagField->xSize());
      const cell_idx_t X = cell_idx_c(flagField->xSize()), Y = cell_idx_c(flagField->ySize()), Z = cell_idx_c(flagField->zSize());
      for (cell_idx_t z = -1; z <= Z; ++z)
         for (cell_idx_t y = -1; y <= Y; ++y)
            for (cell_idx_t x = -1; x <= X; ++x)
            {
               const auto f      = flagField->get(x, y, z);
               const bool isObst = (f & obstacle) != 0;
               const bool ghost  = (x < 0 || y < 0 || z < 0 || x >= X || y >= Y || z >= Z);
               int        kind   = -1;
               if (!isObst && (f & (liquid | interface | gas)) && !finite3(normalField->get(x, y, z))) kind = 0;
               else if (!std::isfinite(smoothField->get(x, y, z))) kind = 1;
               else if (isObst && !finite3(obstNField->get(x, y, z))) kind = 2;
               if (kind < 0) continue;
               if (!ghost) ++counted[kind];
               const uint_t k = ++logged[kind];
               if (!(k <= uint_t(30) || k % uint_t(20000) == uint_t(0))) continue;
               std::ostringstream o;
               o << "[" << names[kind] << " t=" << t << "] (#" << k << ") cell (" << x << "," << y << "," << z << ")"
                 << (ghost ? " [ghost]" : "") << " world~(" << aabb.xMin() + (real_c(x) + real_t(0.5)) * dx << ","
                 << aabb.yMin() + (real_c(y) + real_t(0.5)) * dx << "," << aabb.zMin() + (real_c(z) + real_t(0.5)) * dx
                 << ") flags=" << uint64_t(f) << ((f & interface) ? " I" : "") << ((f & liquid) ? " L" : "")
                 << ((f & gas) ? " G" : "") << (isObst ? " OBST" : "") << " fill=" << fillField->get(x, y, z)
                 << " smoothFill=" << smoothField->get(x, y, z) << " normal=" << normalField->get(x, y, z)
                 << " obstacleNormal=" << obstNField->get(x, y, z);
               int nL = 0, nI = 0, nG = 0, nO = 0;
               std::ostringstream obs, bad;
               for (int dz = -1; dz <= 1; ++dz)
                  for (int dy = -1; dy <= 1; ++dy)
                     for (int ddx = -1; ddx <= 1; ++ddx)
                     {
                        if (!dz && !dy && !ddx) continue;
                        const cell_idx_t nx = x + ddx, ny = y + dy, nz = z + dz;
                        if (nx < -1 || ny < -1 || nz < -1 || nx > X || ny > Y || nz > Z) continue; // outside ghost layer 1
                        const auto nf = flagField->get(nx, ny, nz);
                        if (nf & obstacle)
                        {
                           ++nO;
                           obs << " (" << ddx << "," << dy << "," << dz << "):sf=" << smoothField->get(nx, ny, nz)
                               << ",on=" << obstNField->get(nx, ny, nz);
                        }
                        else if (nf & liquid) ++nL;
                        else if (nf & interface) ++nI;
                        else if (nf & gas) ++nG;
                        if (!std::isfinite(smoothField->get(nx, ny, nz)))
                           bad << " (" << ddx << "," << dy << "," << dz << ")";
                     }
               o << " | neighbours L/I/G/obst=" << nL << "/" << nI << "/" << nG << "/" << nO;
               if (nO) o << " | obstacle neighbours:" << obs.str();
               if (!bad.str().empty()) o << " | non-finite smoothFill at:" << bad.str();
               WALBERLA_LOG_WARNING(o.str());
            }
   }
   if (t % uint_t(1000) == uint_t(0))
   {
      std::vector< uint_t > c(counted, counted + 3);
      mpi::reduceInplace(c, mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("[NAN-NORMAL-SUMMARY t=" << t << "] inner cells with non-finite normal / smoothFill / "
                                "obstacleNormal since the last summary: " << c[0] << " / " << c[1] << " / " << c[2]);
      counted[0] = counted[1] = counted[2] = uint_t(0);
      checkOuterGhostLayerFlags(simObjs, "t=" + std::to_string(t));
   }
}

// Momentum budget diagnostics (ModelParameters.logMomentumBudget).
// Rank-local {fluid Px, Py, Pz, fluid mass, bodies Px, Pz}: liquid + fill-weighted interface inner cells, and m*v of
// the local non-ghost bodies with finite mass (invMass > 0).
inline std::array< real_t, 6 > localFluidBodyMomentum(const std::shared_ptr< SimulationObjects >& simObjs)
{
   real_t px = 0, py = 0, pz = 0, mass = 0;
   for (auto blockIt = simObjs->blockForest->begin(); blockIt != simObjs->blockForest->end(); ++blockIt)
   {
      const auto* flagField = blockIt->getData< const FlagField_T >(simObjs->flagFieldID);
      const auto* fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);
      const auto* pdfField  = blockIt->getData< const PdfField_T >(simObjs->pdfFieldID);
      const auto  liquid    = flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID);
      const auto  interface = flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID);
      const cell_idx_t X = cell_idx_c(flagField->xSize()), Y = cell_idx_c(flagField->ySize()), Z = cell_idx_c(flagField->zSize());
      for (cell_idx_t z = 0; z < Z; ++z)
         for (cell_idx_t y = 0; y < Y; ++y)
            for (cell_idx_t x = 0; x < X; ++x)
      {
         const auto   f = flagField->get(x, y, z);
         const real_t w = isFlagSet(f, liquid) ? real_t(1)
                          : isFlagSet(f, interface) ? std::clamp(fillField->get(x, y, z), real_t(0), real_t(1)) : real_t(0);
         if (w > real_t(0))
         {
            const Vector3< real_t > m = pdfField->getMomentumDensity(x, y, z);
            px += w * m[0]; py += w * m[1]; pz += w * m[2];
            mass += w * pdfField->getDensity(x, y, z);
         }
      }
   }
   real_t bx = 0, bz = 0;
   const auto& ac = *simObjs->particleAccessor;
   for (size_t i = 0; i < ac.size(); ++i)
   {
      if (isSet(ac.getFlags(i), mesa_pd::data::particle_flags::GHOST) || ac.getInvMass(i) <= real_t(0)) continue;
      bx += ac.getLinearVelocity(i)[0] / ac.getInvMass(i);
      bz += ac.getLinearVelocity(i)[2] / ac.getInvMass(i);
   }
   return { px, py, pz, mass, bx, bz };
}

// Change of {fluid Px, fluid Pz, bodies Px, bodies Pz} attributed to the timeloop stage that ends at each probe, summed
// since the last report (rank-local). A probe charges everything since the previous probe to its own stage; the stage
// list is the order in which the probes run within one time step.
struct MomentumStageBudget
{
   static constexpr const char* names[] = { "boundary(gas blocks)", "forceWeighting", "streamReconstructAdvect", "collide",
                                            "cellConversion", "refilling", "excessMass+flagReset", "DEM", "mapping",
                                            "reconstruction" };
   static constexpr size_t nStages = sizeof(names) / sizeof(names[0]);
   std::array< real_t, 4 > last{};
   bool init = false;
   std::array< real_t, 4 * nStages > acc{};
};
inline MomentumStageBudget& momentumStageBudget()
{
   static MomentumStageBudget b;
   return b;
}
inline void momentumStageProbe(const std::shared_ptr< SimulationObjects >& simObjs, const size_t stage)
{
   const auto m = localFluidBodyMomentum(simObjs);
   const std::array< real_t, 4 > p{ m[0], m[2], m[4], m[5] };
   auto& b = momentumStageBudget();
   if (b.init)
      for (size_t k = 0; k < 4; ++k) b.acc[4 * stage + k] += p[k] - b.last[k];
   b.last = p;
   b.init = true;
}

// [MOMENTUM-BUDGET] line: momentum the freed cells carried into the fluid and the covered cells took out of it since
// the last report (each {Px, Py, Pz, L_y about the body centre, mass, cells}, rank-local, see
// PdfReconstructionManager::takeFreedMomentum), the total fluid momentum now and the bodies' momentum m*v.
// [MOMENTUM-STAGES] line: the stage budget above and the MEM link sums of SimpleBB.h (MemLinkMomentumBudget).
// Resets `covered`, the stage budget and the link sums.
inline void logMomentumBudget(const std::shared_ptr< SimulationObjects >& simObjs, const uint_t t,
                              const std::array< real_t, 6 >& freed, std::array< real_t, 6 >& covered)
{
   std::vector< real_t > v(freed.begin(), freed.end());
   v.insert(v.end(), covered.begin(), covered.end());
   covered.fill(real_t(0));
   const auto now = localFluidBodyMomentum(simObjs);
   v.insert(v.end(), { now[0], now[1], now[2], now[3], now[4], now[5] });
   auto& stages = momentumStageBudget();
   v.insert(v.end(), stages.acc.begin(), stages.acc.end());
   stages.acc.fill(real_t(0));
   auto& links = lbm_mesapd_coupling::MemLinkMomentumBudget::sums();
   v.insert(v.end(), links.begin(), links.end());
   links.fill(real_t(0));
   walberla::mpi::reduceInplace(v, walberla::mpi::SUM);
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_INFO("[MOMENTUM-BUDGET t=" << t << "] since last report: freed cells n=" << v[5] << " P=(" << v[0] << ","
         << v[1] << "," << v[2] << ") Ly=" << v[3] << " mass=" << v[4] << " | covered cells n=" << v[11] << " P=(" << v[6] << ","
         << v[7] << "," << v[8] << ") Ly=" << v[9] << " mass=" << v[10] << " | now: fluid P=(" << v[12] << "," << v[13] << ","
         << v[14] << ") mass=" << v[15] << " bodies P=(" << v[16] << ",.," << v[17] << ")");
      std::ostringstream os;
      os << "[MOMENTUM-STAGES t=" << t << "] since last report, per stage dP fluid (x,z) / bodies (x,z):";
      const size_t o = 18;
      for (size_t k = 0; k < MomentumStageBudget::nStages; ++k)
         os << " " << MomentumStageBudget::names[k] << "=(" << v[o + 4 * k] << "," << v[o + 4 * k + 1] << ")/(" << v[o + 4 * k + 2]
            << "," << v[o + 4 * k + 3] << ")";
      const size_t l = o + 4 * MomentumStageBudget::nStages;
      os << " | MEM links (x,z): fluidLoss=(" << v[l] << "," << v[l + 1] << ") fluidLossFillWeighted=(" << v[l + 2] << ","
         << v[l + 3] << ") refSubtraction=(" << v[l + 4] << "," << v[l + 5] << ") wenCorrection=(" << v[l + 6] << "," << v[l + 7]
         << ") givenToBodies=(" << v[l + 8] << "," << v[l + 9] << ")";
      WALBERLA_LOG_INFO(os.str());
   }
}

void addSweepSurfaceGeometry(const std::shared_ptr< const SimulationParams >& simParams,
                             const std::shared_ptr< SimulationObjects >& simObjs,
                             const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addSweepSurfaceGeometry");

   // Debug probe: fixed watch cells (trackWatchCells) at the start of the FSLBM chain during steps 425-444, and the
   // dense NaN scan (PhysicsParameters.denseNanScan*).
   timeloop->add() << Sweep(EmptySweep(), "Empty sweep: probe anchor (surface geometry start)")
                   << AfterFunction(
                         [simObjs]() {
                            const uint_t t = *simObjs->currentMainStep;
                            if (t < uint_t(425) || t > uint_t(444)) return;
                            trackWatchCells(simObjs, "PROBE-start-of-FSLBM-timestep", t,
                                           { Vector3< real_t >(real_t(20.5), real_t(8.5), real_t(25.5)),
                                             Vector3< real_t >(real_t(20.5), real_t(6.5), real_t(25.5)),
                                             Vector3< real_t >(real_t(20.5), real_t(5.5), real_t(25.5)),
                                             Vector3< real_t >(real_t(18.5), real_t(4.5), real_t(26.5)),
                                             Vector3< real_t >(real_t(18.5), real_t(34.5), real_t(26.5)) });
                         },
                         "Probe: start of timeloopFSLBM")
                   << AfterFunction(
                         [simParams, simObjs]() {
                            if (!simParams->physics.denseNanScanEnabled) return;
                            const uint_t t = *simObjs->currentMainStep;
                            if (t < simParams->physics.denseNanScanFromStep || t >= simParams->physics.denseNanScanToStep) return;
                            checkPdfNaN(simObjs, "INTRA-start-of-FSLBM-timestep", t,
                                       simParams->physics.lowRhoPrecursorThreshold,
                                       simParams->physics.highVelocityPrecursorThreshold,
                                       simParams->physics.populationPrecursorThreshold,
                                       simParams->physics.densityPrecursorLow, simParams->physics.densityPrecursorHigh,
                                       simParams->physics.densityPrecursorLiquidOnly);
                         },
                         "Intra-step NaN scan: start of timeloopFSLBM");

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
   // ObstacleFillLevelSweep running above, obstacle cells carry meaningful fill, so we match the
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
                                    "Communication: after normal sweep")
                   // ModelParameters.diagnoseNaNNormals
                   << AfterFunction(
                         [simObjs, on = simParams->model.diagnoseNaNNormals]() {
                            if (on) diagnoseNaNNormals(simObjs, *simObjs->currentMainStep);
                         },
                         "Diagnostic: non-finite normals");

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
      // ModelParameters.guardCurvatureInputs (default false = upstream behaviour)
      curvSweep.setGuardInvalidNeighbours(simParams->model.guardCurvatureInputs);
      WALBERLA_LOG_INFO_ON_ROOT("Curvature input guard (guardCurvatureInputs): "
                                << (simParams->model.guardCurvatureInputs ? "on" : "off"));

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

   // Debug probe: logs fixed watch cells (trackWatchCells) after the sweeps below, during steps 425-444 (not after the
   // weighted / interface-and-liquid excess mass variants).
   auto probeCrashCell = [simObjs](const char* label) {
      const uint_t t = *simObjs->currentMainStep;
      if (t < uint_t(425) || t > uint_t(444)) return;
      trackWatchCells(simObjs, label, t, { Vector3< real_t >(real_t(20.5), real_t(8.5), real_t(25.5)),
                                            Vector3< real_t >(real_t(20.5), real_t(6.5), real_t(25.5)),
                                            Vector3< real_t >(real_t(20.5), real_t(5.5), real_t(25.5)),
                                            Vector3< real_t >(real_t(18.5), real_t(4.5), real_t(26.5)),
                                            Vector3< real_t >(real_t(18.5), real_t(34.5), real_t(26.5)) });
   };

   // Domain-wide NaN scan at the same points as probeCrashCell (checkPdfNaN with its 3x3x3 neighbourhood dump), in the
   // PhysicsParameters.denseNanScan* window. Silent unless a cell is non-finite or a precursor threshold fires; the
   // result is ignored (no abort), so the stage-by-stage progression within the step stays visible.
   auto scanIntraStepNaN = [simParams, simObjs](const char* label) {
      if (!simParams->physics.denseNanScanEnabled) return;
      const uint_t t = *simObjs->currentMainStep;
      if (t < simParams->physics.denseNanScanFromStep || t >= simParams->physics.denseNanScanToStep) return;
      checkPdfNaN(simObjs, label, t, simParams->physics.lowRhoPrecursorThreshold,
                 simParams->physics.highVelocityPrecursorThreshold, simParams->physics.populationPrecursorThreshold,
                 simParams->physics.densityPrecursorLow, simParams->physics.densityPrecursorHigh,
                 simParams->physics.densityPrecursorLiquidOnly);
   };

   // add standard waLBerla boundary handling
   // runs on every block except onlyGasAndBoundary blocks, which get the empty sweep
   timeloop->add() << Sweep(simObjs->antidunesBoundaryHandling->getBoundarySweep(), "Sweep: boundary handling",
                            Set< SUID >::emptySet(), StateSweep::onlyGasAndBoundary)
                   << Sweep(EmptySweep(), "Empty sweep: boundary handling", StateSweep::onlyGasAndBoundary)
                   << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-boundarySweep"); },
                                    "Probe: after boundary handling")
                   << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-boundarySweep"); },
                                    "Intra-step NaN scan: after boundary handling")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 0); },
                                    "Momentum budget probe: after boundary handling");

   // compute force density field (acceleration applied to fluid and interface cells)
   const ForceDensityCodegenSweep< PdfField_T, FlagField_T, VectorFieldFlattened_T, ScalarField_T >
      forceDensityCodegenSweep(simObjs->forceDensityFieldID, simObjs->pdfFieldID, simObjs->flagFieldID,
                               simObjs->fillFieldID, simObjs->flagInfo, simObjs->acceleration);
   timeloop->add() << Sweep(forceDensityCodegenSweep, "Sweep: force weighting", Set< SUID >::emptySet(),
                            StateSweep::onlyGasAndBoundary)
                   << Sweep(EmptySweep(), "Empty sweep: force weighting", StateSweep::onlyGasAndBoundary)
                   << AfterFunction(Communication_T(simObjs->blockForest, simObjs->forceDensityFieldID),
                                    "Communication: after force weighting sweep")
                   << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-forceDensitySweep"); },
                                    "Probe: after force weighting")
                   << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-forceDensitySweep"); },
                                    "Intra-step NaN scan: after force weighting")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 1); },
                                    "Momentum budget probe: after force weighting");

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
         simParams->model.cellConversionForceThreshold, simParams->model.seedStreamPdfCloneOnFirstUse);

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
                       "Second ghost layer update: after StreamReconstructAdvect sweep (flag field)")
      << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-streamReconstructAdvect"); },
                       "Probe: after StreamReconstructAdvect")
      << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-streamReconstructAdvect"); },
                       "Intra-step NaN scan: after StreamReconstructAdvect")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 2); },
                                    "Momentum budget probe: after StreamReconstructAdvect");

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

   // Debug probe: state (PDFs, density, fill, flag) of two fixed watch cells immediately before and after the generated
   // collision, every step
   // (shows whether collide() also iterates obstacle cells, whose PDFs MovingParticleMapping leaves untouched).
   static const std::vector< Vector3< real_t > > collideProbeCells = {
      Vector3< real_t >(real_t(608.5), real_t(18.5), real_t(30.5)),
      Vector3< real_t >(real_t(609.5), real_t(18.5), real_t(30.5))
   };

   timeloop->add() << BeforeFunction(
                         [simObjs]() { trackWatchCells(simObjs, "COLLIDE-BEFORE", *simObjs->currentMainStep, collideProbeCells); },
                         "Probe: before collide (generated)")
                   << Sweep(CollideSweep(lbmSweepGenerated), "Sweep: collision (generated)",
                            StateSweep::fullFreeSurface)
                   << Sweep(lbmSweepGenerated, "Sweep: streamCollide (generated)", StateSweep::onlyLBM)
                   << Sweep(EmptySweep(), "Empty sweep: streamCollide (generated)")
                   << AfterFunction(
                         [simObjs]() { trackWatchCells(simObjs, "COLLIDE-AFTER", *simObjs->currentMainStep, collideProbeCells); },
                         "Probe: after collide (generated)")
                   << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                                    "Communication: after streamCollide (generated)")
                   << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-collide"); },
                                    "Probe: after collide")
                   << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-collide"); },
                                    "Intra-step NaN scan: after collide")
                   << AfterFunction([simObjs, simParams]() {
                         if (!simParams->model.repairDegenerateObstaclePdfs) return; // default off
                         const uint_t t = *simObjs->currentMainStep;
                         repairDegenerateObstacleCells(simObjs, "after-collide", t);
                      }, "Repair: degenerate obstacle-cell PDFs after collide")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 3); },
                                    "Momentum budget probe: after collide");

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
                         "Second ghost layer update: after cell conversion sweep (flag field)")
                   << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-cellConversion"); },
                                    "Probe: after cell conversion")
                   << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-cellConversion"); },
                                    "Intra-step NaN scan: after cell conversion")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 4); },
                                    "Momentum budget probe: after cell conversion");

   // reinitialize PDFs, i.e., refill cells that were converted from gas to interface
   // - when the flag "convertedFromGasToInterface" has been set (by CellConversionSweep)
   // - always with EquilibriumRefillingSweep; ModelParameters.pdfRefillingModel is parsed but not used
   const walberla::free_surface::EquilibriumRefillingSweep< LatticeModel_T, FlagField_T > equilibriumRefillingSweep(
      simObjs->pdfFieldID, simObjs->flagFieldID, simObjs->flagInfo, true);
   timeloop->add() << Sweep(equilibriumRefillingSweep, "Sweep: EquilibriumRefilling", StateSweep::fullFreeSurface)
                   << Sweep(EmptySweep(), "Empty sweep: EquilibriumRefilling")
                   << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                                    "Communication: after EquilibriumRefilling sweep")
                   << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-equilibriumRefilling"); },
                                    "Probe: after EquilibriumRefilling")
                   << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-equilibriumRefilling"); },
                                    "Intra-step NaN scan: after EquilibriumRefilling")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 5); },
                                    "Momentum budget probe: after EquilibriumRefilling");

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
                      << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); }, "Sweep: bubble model update")
                      << AfterFunction([probeCrashCell]() { probeCrashCell("PROBE-after-excessMassEvenly"); },
                                       "Probe: after excess mass distribution (evenly)")
                      << AfterFunction([scanIntraStepNaN]() { scanIntraStepNaN("INTRA-after-excessMassEvenly"); },
                                       "Intra-step NaN scan: after excess mass distribution (evenly)");
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
                         "Second ghost layer update: after excess mass distribution sweep (flag field)")
                   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 6); },
                                    "Momentum budget probe: after excess mass distribution and conversion flag reset");

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
   // Fixed pillars (box pillars type 2, cylinder pillars) are excluded here; they are mapped separately as
   // static NoSlip/FreeSlip (see SetupStages). Dormant (startInactive, not yet at settlingSteps) floaters
   // are excluded too -- the per-step equivalent of the init-time dormant guard in
   // SetupStages::mapParticlesOntoCells (`floaterSelector`); both must agree or a dormant floater would
   // appear in the fluid the very first step after initialization. Note: `floaterSelector` selects types
   // 1/3/4 only, i.e. bed spheres (type 0) are not mapped at initialisation, only from the first step on.
   auto dormantUIDs = simObjs->dormantFloaterUIDs;
   auto movingObstacleSelector = [dormantUIDs](const size_t particleIdx, const auto& ac) {
      const bool isSphere           = ac.getBaseShape(particleIdx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE;
      const bool isFloatingBox      = ac.getType(particleIdx) == uint_t(3);
      const bool isFloatingCylinder = ac.getType(particleIdx) == uint_t(4);
      if (!(isSphere || isFloatingBox || isFloatingCylinder)) return false;
      return dormantUIDs->count(ac.getUid(particleIdx)) == 0;
   };

   // Debug probe: the same two watch cells around MovingParticleMapping and PdfReconstructionManager, the only
   // operations of timeloopAfterDEM that touch PDFs.
   static const std::vector< Vector3< real_t > > couplingProbeCells = {
      Vector3< real_t >(real_t(608.5), real_t(18.5), real_t(30.5)),
      Vector3< real_t >(real_t(609.5), real_t(18.5), real_t(30.5))
   };

   // ModelParameters.useOwnershipHandoffPdfReset: called by MovingParticleMapping once, right before a cell's obstacle
   // ownership passes from one particle to another; resets the cell's PDFs to rest equilibrium (rho 1, u 0) so it does
   // not inherit the previous owner's stale history. The D3Q27 weights are assigned by hand, as in SetupStages.h's
   // repairDegenerateObstacleCells.
   auto onOwnershipHandoff = [simParams, simObjs](IBlock& block, cell_idx_t x, cell_idx_t y, cell_idx_t z,
                                        walberla::id_t formerUid, walberla::id_t newUid) {
      if (!simParams->model.useOwnershipHandoffPdfReset) return; // off: the cell keeps its PDFs
      // MovingParticleMapping.h calls this only for cells in pdfField's inner part (flag and PDF fields have different
      // ghost-layer counts), so (x, y, z) is always safe to write.
      auto* pdfField = block.getData< PdfField_T >(simObjs->pdfFieldID);
      pdfField->get(x, y, z, uint_t(0)) = real_t(8.0 / 27.0);
      for (uint_t q = uint_t(1); q <= uint_t(6); ++q) pdfField->get(x, y, z, q) = real_t(2.0 / 27.0);
      for (uint_t q = uint_t(7); q <= uint_t(18); ++q) pdfField->get(x, y, z, q) = real_t(1.0 / 54.0);
      for (uint_t q = uint_t(19); q <= uint_t(26); ++q) pdfField->get(x, y, z, q) = real_t(1.0 / 216.0);

      static uint_t handoffLogCount = uint_t(0);
      ++handoffLogCount;
      if (handoffLogCount <= uint_t(50) || (handoffLogCount % uint_t(500) == uint_t(0)))
      {
         WALBERLA_LOG_WARNING("[OWNERSHIP-HANDOFF-RESET] (#" << handoffLogCount << ") t="
            << *simObjs->currentMainStep << " block-local=(" << x << "," << y << "," << z
            << ") ownership handed off uid=" << formerUid << " -> uid=" << newUid
            << " -- PDF reset to rest equilibrium (rho=1, u=0)");
      }
   };

   // ModelParameters.logMomentumBudget: fluid momentum lost to cells the bodies cover, weighted like fluid mass (liquid
   // 1, interface its fill, gas 0) -- {Px, Py, Pz, L_y about the covering particle's centre, mass, cells}, rank-local,
   // reported and reset with the conversion counts below.
   auto coveredMomentum = std::make_shared< std::array< real_t, 6 > >();
   coveredMomentum->fill(real_t(0));
   std::function< void(IBlock&, cell_idx_t, cell_idx_t, cell_idx_t, size_t) > onNewlyCovered; // see MovingParticleMapping.h
   if (simParams->model.logMomentumBudget)
      onNewlyCovered = [simObjs, coveredMomentum](IBlock& block, cell_idx_t x, cell_idx_t y, cell_idx_t z, size_t idx) {
         const auto* flagField = block.getData< const FlagField_T >(simObjs->flagFieldID);
         const auto* fillField = block.getData< const ScalarField_T >(simObjs->fillFieldID);
         const auto* pdfField  = block.getData< const PdfField_T >(simObjs->pdfFieldID);
         const auto  f         = flagField->get(x, y, z);
         const real_t w = isFlagSet(f, flagField->getFlag(walberla::free_surface::flagIDs::liquidFlagID)) ? real_t(1)
                          : isFlagSet(f, flagField->getFlag(walberla::free_surface::flagIDs::interfaceFlagID))
                             ? std::clamp(fillField->get(x, y, z), real_t(0), real_t(1))
                             : real_t(0);
         if (w <= real_t(0)) return;
         const Vector3< real_t > m = pdfField->getMomentumDensity(x, y, z);
         const Vector3< real_t > r =
            simObjs->blockForest->getBlockLocalCellCenter(block, Cell(x, y, z)) - simObjs->particleAccessor->getPosition(idx);
         auto& c = *coveredMomentum;
         c[0] += w * m[0]; c[1] += w * m[1]; c[2] += w * m[2];
         c[3] += w * (r[2] * m[0] - r[0] * m[2]);
         c[4] += w * pdfField->getDensity(x, y, z);
         c[5] += real_t(1);
      };

   // sweep for mapping particles into the LBM cells
   const bool strictlyConserveMomentum = false;
   lbm_mesapd_coupling::MemLinkMomentumBudget::enabled() = simParams->model.logMomentumBudget;
   timeloop->add() << BeforeFunction(
      [simObjs]() { trackWatchCells(simObjs, "BEFORE-MAPPING", *simObjs->currentMainStep, couplingProbeCells); },
      "Probe: before MovingParticleMapping")
   << BeforeFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 7); },
                     "Momentum budget probe: after DEM")
   << Sweep(
      lbm_mesapd_coupling::makeMovingParticleMapping< PdfField_T, AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         simObjs->blockForest, simObjs->pdfFieldID, simObjs->antidunesBoundaryHandling->getHandlingID(),
         simObjs->particleFieldID, simObjs->particleAccessor, AntidunesBoundaryHandling_T::movingObstacleFlagID,
         FormerMO_Flag,
         movingObstacleSelector, // all spheres + floating boxes/cylinders (not fixed pillars)
         strictlyConserveMomentum, simParams->model.mappingTieOffset, onOwnershipHandoff, onNewlyCovered),
      "Sweep (coupling): MovingParticleMapping")
   << AfterFunction(
      [simObjs]() { trackWatchCells(simObjs, "AFTER-MAPPING", *simObjs->currentMainStep, couplingProbeCells); },
      "Probe: after MovingParticleMapping, before PdfReconstructionManager")
   << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 8); },
                    "Momentum budget probe: after MovingParticleMapping");

   // Restore PDFs in cells previously occupied by particles: Bogner & Rude (2013) state selection with the
   // EquilibriumReconstructor, which sets all 27 PDFs to f_i^eq(rho_avg, u_particle) >= 0. (Grad's moment
   // approximation adds a non-equilibrium term from neighbouring fluid velocities; near the contact line, with gas and
   // obstacle neighbours, that term can make f_i negative and the collision diverge.)
   // Each FormerMO cell is classified by its D3Q27 non-obstacle neighbourhood B(x) (default rule; see
   // strictBognerClassification for Bogner's literal rule):
   //   gas, no liquid, no unresolved same-pass neighbour     ->  GAS        (fill 0)
   //   liquid, no gas, no unresolved same-pass neighbour     ->  LIQUID     (fill 1)
   //   everything else                                       ->  INTERFACE  (fill seeded in pass 2: mean fill of the
   //                                                             interface cells in B(x), 0.5 if none; see
   //                                                             seedFillFromAllNeighbours)
   // useBognerReconstruction=false reconstructs every FormerMO cell as interface without fill seeding; the conversion
   // counters then stay zero. Held as a shared_ptr so the AfterFunction below can log and reset those counters.
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
         simParams->staticReferenceHeight(),             // flat static H (SimpleBB's, unless useDynamicReferenceHeight)
         simParams->model.usePass3Reresolution);
   // ModelParameters.seedFillFromAllNeighbours (default false)
   bognerReconstructionManager->setSeedFillFromAllNeighbours(simParams->model.seedFillFromAllNeighbours);
   // ModelParameters.seedSkipObstacleNeighbours (default true)
   bognerReconstructionManager->setSeedSkipObstacleNeighbours(simParams->model.seedSkipObstacleNeighbours);
   WALBERLA_LOG_INFO_ON_ROOT("FormerMO interface fill seeding: "
                             << (simParams->model.seedFillFromAllNeighbours ? "all resolved neighbours" : "interface neighbours (legacy)")
                             << ", obstacle-flagged neighbours " << (simParams->model.seedSkipObstacleNeighbours ? "skipped" : "counted if interface (pre-2026-10-01)"));
   // ModelParameters.strictBognerClassification (default false = same-step interface fallback)
   bognerReconstructionManager->setStrictBognerClassification(simParams->model.strictBognerClassification);
   // ModelParameters.strictBognerNoLiquidAtGas (default false)
   bognerReconstructionManager->setStrictBognerNoLiquidAtGas(simParams->model.strictBognerNoLiquidAtGas);
   WALBERLA_LOG_INFO_ON_ROOT("strictBognerNoLiquidAtGas = " << (simParams->model.strictBognerNoLiquidAtGas ? "true" : "false"));
   // ModelParameters.strictBognerRelaxInterface (default false); no-op unless strictBognerClassification and
   // strictBognerNoLiquidAtGas are also on.
   bognerReconstructionManager->setStrictBognerRelaxInterface(simParams->model.strictBognerRelaxInterface);
   WALBERLA_LOG_INFO_ON_ROOT("strictBognerRelaxInterface = " << (simParams->model.strictBognerRelaxInterface ? "true" : "false"));
   WALBERLA_LOG_INFO_ON_ROOT("FormerMO classification: "
                             << (simParams->model.strictBognerClassification ? "strict Bogner rule (same-step freed cells not in B(x))"
                                                                             : "same-step interface fallback"));
   WALBERLA_LOG_INFO_ON_ROOT("Pass-3 local re-resolution: " << (simParams->model.usePass3Reresolution ? "on" : "off"));

   const uint_t convLogFreq = simParams->evaluation.evaluationFrequency;

   timeloop->add()
      << BeforeFunction(
            [simObjs]() { trackWatchCells(simObjs, "BEFORE-RECONSTRUCTION", *simObjs->currentMainStep, couplingProbeCells); },
            "Probe: before PdfReconstructionManager")
      << Sweep(makeSharedSweep(bognerReconstructionManager),
               "Sweep (coupling): PDF Reconstruction (moving particles) [Bogner2way+EquilibriumReconstructor]")
      << AfterFunction(
            [simObjs]() { trackWatchCells(simObjs, "AFTER-RECONSTRUCTION", *simObjs->currentMainStep, couplingProbeCells); },
            "Probe: after PdfReconstructionManager, before its own communication")
      << AfterFunction([simObjs, simParams]() { if (simParams->model.logMomentumBudget) momentumStageProbe(simObjs, 9); },
                       "Momentum budget probe: after PdfReconstructionManager")
      << AfterFunction(
            [bognerReconstructionManager, convLogFreq, simObjs, simParams, coveredMomentum, step = uint_t(0)]() mutable {
               ++step;
               if (convLogFreq > uint_t(0) && step % convLogFreq == uint_t(0))
               {
                  bognerReconstructionManager->logAndResetConversionCounts(step);
                  if (simParams->model.logMomentumBudget) logMomentumBudget(simObjs, step, bognerReconstructionManager->takeFreedMomentum(),
                                                                            *coveredMomentum);
               }
            },
            "AfterFunction (coupling): log FormerMO→state conversion counts")
      << AfterFunction(
            Communication_T(simObjs->blockForest, simObjs->flagFieldID, simObjs->particleFieldID),
            "Communication (coupling): after PDF reconstruction (flagField + particleField)") // unsure if necessary but
                                                                                              // added for consistency
      << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                       "Communication (coupling): after PDF reconstruction (pdfField)")
      << AfterFunction(
            [simObjs]() { trackWatchCells(simObjs, "AFTER-RECONSTRUCTION-COMM", *simObjs->currentMainStep, couplingProbeCells); },
            "Probe: after PdfReconstructionManager's own communication");
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

   // logs data and aborts if it detects non-finite values in liquid or interface cells
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
                                simParams->evaluation.dataSliceFrequency);
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

   // clog measurement evaluator (optional "ClogMeasurement" prm block) -----------------------------------------
   if (simParams->evaluation.clogMeasurementEnabled)
   {
      // Blockage-ratio inputs: the cylinder pillar geometry is passed in here so that ClogMeasurementEvaluator stays
      // independent of CylinderPillarParams. Box pillars (PillarParams) are not passed.
      std::vector< Vector3< real_t > > pillarCenters;
      std::vector< real_t >            pillarRadii;
      for (const auto& pillar : simParams->cylinderPillars)
      {
         pillarCenters.push_back(pillar.center);
         pillarRadii.push_back(pillar.radius);
      }
      const real_t channelWidthY = real_c(simParams->domain.domainSize[1]);

      ClogMeasurementEvaluator< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
         clogMeasurementEvaluator(simObjs->blockForest, simObjs->flagFieldID, simObjs->fillFieldID, simObjs->pdfFieldID,
                                  simParams->evaluation.clogUpstreamPos, simParams->evaluation.clogDownstreamPos,
                                  simParams->evaluation.clogProbeHalfWidthXY, simObjs->acceleration,
                                  simParams->computed.absoluteLiquidHeight, simObjs->dataCols,
                                  simParams->evaluation.evaluationFrequency, "clog",
                                  simParams->evaluation.clogFlowRateEnabled, pillarCenters, pillarRadii,
                                  channelWidthY, simParams->evaluation.clogBlockageRatioEnabled);
      timeloop->addFuncAfterTimeStep(clogMeasurementEvaluator, "Evaluator: ClogMeasurement");
   }

   // Optional second probe pair at the pillar's y (see EvaluationParams, SimulationParams.h): same evaluator class,
   // writing "clogPillar*" columns.
   if (simParams->evaluation.clogPillarMeasurementEnabled)
   {
      ClogMeasurementEvaluator< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
         clogPillarMeasurementEvaluator(simObjs->blockForest, simObjs->flagFieldID, simObjs->fillFieldID, simObjs->pdfFieldID,
                                        simParams->evaluation.clogPillarUpstreamPos, simParams->evaluation.clogPillarDownstreamPos,
                                        simParams->evaluation.clogProbeHalfWidthXY, simObjs->acceleration,
                                        simParams->computed.absoluteLiquidHeight, simObjs->dataCols,
                                        simParams->evaluation.evaluationFrequency, "clogPillar");
      timeloop->addFuncAfterTimeStep(clogPillarMeasurementEvaluator, "Evaluator: ClogPillarMeasurement");
   }
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
                                   simParams->evaluation.evaluationFrequency, simParams->pid.enabled,
                                   simObjs->beginTimeStep, simParams->pid.switchOffStep,
                                   simParams->pid.switchOffAveragingWindow);

   timeloop->addFuncAfterTimeStep(forcingAdjuster, "Evaluator: ForcingAdjuster");
}

// Optional sweep (ModelParameters.diagnosticSweepEnabled), once per main step on timeloopEval, i.e. after all LBM/DEM
// and force computations of the step, within [diagnosticSweepFromStep, diagnosticSweepToStep) (ToStep 0 = no limit):
//   1. drains and logs SimulationObjects::diagnosticEvents, filled by checks in ParticleKernels.h (PRE-REDUCE,
//      RAW-HYDRO, VELOCITY-JUMP, MIGRATION-NAN) instead of logging at the point of detection;
//   2. runs checkPdfNaN (SetupStages.h) once over all liquid/interface cells with the density window
//      [densityPrecursorLow, densityPrecursorHigh], independent of the denseNanScan window;
//   3. logs every local floater's position, quaternion and contact/hydrodynamic force and torque;
//   4. writes the logFluidBoxHistory CSV, if enabled.
void addAfterTimeStepDiagnosticSweep(const std::shared_ptr< const SimulationParams >& simParams,
                                     const std::shared_ptr< SimulationObjects >& simObjs,
                                     const std::unique_ptr< SweepTimeloop >& timeloop)
{
   if (!simParams->model.diagnosticSweepEnabled) return; // not even registered when off -- zero cost
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepDiagnosticSweep");

   auto diagnosticSweep = [simParams, simObjs]() {
      const uint_t t = *simObjs->currentMainStep;
      if (t < simParams->model.diagnosticSweepFromStep) return;
      if (simParams->model.diagnosticSweepToStep != uint_t(0) && t >= simParams->model.diagnosticSweepToStep) return;

      // (1) Drain + log this step's collected per-sub-cycle diagnostic events, then clear for next step.
      auto& events = *simObjs->diagnosticEvents;
      for (const auto& ev : events)
      {
         WALBERLA_LOG_WARNING("[DIAG-SWEEP][" << ev.tag << "] rank=" << mpi::MPIManager::instance()->rank()
            << " step=" << ev.mainStep << " subCycle=" << ev.subCycle << " uid=" << ev.uid
            << " pos=" << ev.worldPos << " value=" << ev.value << " -- " << ev.detail);
      }
      events.clear();

      // (2) Domain-wide density scan -- reuses checkPdfNaN with only the density window of the dense NaN scan.
      checkPdfNaN(simObjs, "DIAG-SWEEP-DENSITY", t, real_t(0), real_t(0), real_t(0),
                 simParams->physics.densityPrecursorLow, simParams->physics.densityPrecursorHigh,
                 simParams->physics.densityPrecursorLiquidOnly);

      // (3) Per-floater position/quaternion/force/torque snapshot.
      simObjs->particleStorage->forEachParticle(false,
         [](size_t idx, auto& ac) {
            using namespace mesa_pd::data::particle_flags;
            if (isSet(ac.getFlags(idx), GHOST)) return false;
            const uint_t pt = ac.getType(idx);
            return pt == uint_t(1) || pt == uint_t(3) || pt == uint_t(4); // sphere/box/cylinder floaters
         },
         *simObjs->particleAccessor,
         [t](size_t idx, auto& ac) {
            const auto& quat = ac.getRotation(idx).getQuaternion();
            WALBERLA_LOG_INFO("[DIAG-SWEEP-FLOATER] rank=" << mpi::MPIManager::instance()->rank()
               << " step=" << t << " uid=" << ac.getUid(idx)
               << " pos=" << ac.getPosition(idx) << " quat=(" << quat[0] << "," << quat[1] << ","
               << quat[2] << "," << quat[3] << ")"
               << " force=" << ac.getForce(idx) << " torque=" << ac.getTorque(idx)
               << " hydroForce=" << ac.getHydrodynamicForce(idx)
               << " hydroTorque=" << ac.getHydrodynamicTorque(idx));
         },
         *simObjs->particleAccessor);

      // ModelParameters.logFluidBoxHistory: obstacle/interface cells of a fixed box, in its own step window.
      if (simParams->model.logFluidBoxHistory &&
          t >= simParams->model.logFluidBoxHistoryStartStep &&
          (simParams->model.logFluidBoxHistoryEndStep == uint_t(0) || t <= simParams->model.logFluidBoxHistoryEndStep))
      {
         logFluidBoxHistory(simObjs, t, simParams->model.fluidBoxHistoryBaseFolder,
                            simParams->model.fluidBoxMin, simParams->model.fluidBoxMax);
      }
   };

   timeloop->addFuncAfterTimeStep(diagnosticSweep, "Diagnostic: toggleable density+floater sweep");
}

// ModelParameters.logDensityColumn, independent of the diagnostic sweep above.
void addDensityColumnLoggingSweep(const std::shared_ptr< const SimulationParams >& simParams,
                                  const std::shared_ptr< SimulationObjects >& simObjs,
                                  const std::unique_ptr< SweepTimeloop >& timeloop)
{
   if (!simParams->model.logDensityColumn) return; // not even registered when off -- zero cost

   auto densityColumnSweep = [simParams, simObjs]() {
      const uint_t t = *simObjs->currentMainStep;
      if (t < simParams->model.logDensityColumnStartStep) return;
      if (simParams->model.logDensityColumnEndStep != uint_t(0) && t > simParams->model.logDensityColumnEndStep) return;

      logDensityColumn(simObjs, t, simParams->model.densityColumnBaseFolder,
                       simParams->model.densityColumnX, simParams->model.densityColumnY);
   };

   timeloop->addFuncAfterTimeStep(densityColumnSweep, "Diagnostic: far-field density column log");
}

// ModelParameters.logWaterHeight, independent of every other diagnostic switch, over the whole run.
void addWaterHeightLoggingSweep(const std::shared_ptr< const SimulationParams >& simParams,
                                const std::shared_ptr< SimulationObjects >& simObjs,
                                const std::unique_ptr< SweepTimeloop >& timeloop)
{
   if (!simParams->model.logWaterHeight) return; // not even registered when off -- zero cost

   auto waterHeightSweep = [simParams, simObjs]() {
      const uint_t t = *simObjs->currentMainStep;
      const uint_t every = simParams->model.waterHeightEveryNSteps == uint_t(0) ? uint_t(1) : simParams->model.waterHeightEveryNSteps;
      if (t % every != uint_t(0)) return;
      logWaterHeights(simObjs, t, simParams->model.waterHeightBaseFolder, simParams->model.waterHeightGridX,
                      simParams->model.waterHeightGridY);
   };

   timeloop->addFuncAfterTimeStep(waterHeightSweep, "Diagnostic: multi-column water height log");
}

} // namespace antidunes
} // namespace walberla