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
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"
#include "lbm/free_surface/surface_geometry/Utility.h"


#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"


#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/VTKOutput.h"

#include "Evaluators.h"
#include "Forcing.h"
#include "SimulationParams.h" 
#include "SimulationObjects.h" 

namespace walberla {
namespace antidunes {

using StateSweep = walberla::free_surface::BlockStateDetectorSweep< FlagField_T >;

// Empty sweep required for blocks that are not sweeped (i.e., excluded via blockstate selectors)
struct EmptySweep { void operator()(IBlock*) {} };

void addBeforeTimeStepVTKOutput(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addBeforeTimeStepVTKOutput");

   // LBM/freeSurface fields -----------------------------------------------
   walberla::free_surface::addVTKOutput< LatticeModel_T, AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T, VectorField_T >(
      simObjs->blockForest, *timeloop, simParams->cfgPtr, simObjs->flagInfo,
      simObjs->pdfFieldID, simObjs->flagFieldID, simObjs->fillFieldID,
      BlockDataID{}, simObjs->curvatureFieldID,
      simObjs->normalFieldID, simObjs->obstacleNormalFieldID);
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


   // particles ------------------------------------------------------------
   auto particlePointData = make_shared< mesa_pd::vtk::ParticleVtkOutput >(simObjs->particleStorage);
   // particlePointData->addOutput< mesa_pd::data::SelectParticleUid >("uid");
   // particlePointData->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
   particlePointData->addOutput< mesa_pd::data::SelectParticleInteractionRadius >("radius");
   particlePointData->setParticleSelector([](const mesa_pd::data::ParticleStorage::iterator& pIt) {
      using namespace walberla::mesa_pd::data::particle_flags;
      return (pIt->getBaseShape()->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE) &&
               !isSet(pIt->getFlags(), GHOST);
   });
   auto particleVTKOutput = vtk::createVTKOutput_PointData(particlePointData, "particles",
      simParams->vtk.writeFrequency, simParams->vtk.baseFolder, std::string("simulation_step"),
      false, false, true, false, simObjs->beginTimeStep);
   timeloop->addFuncBeforeTimeStep(vtk::writeFiles(particleVTKOutput), "VTK: sphere data");
}

void addSweepSurfaceGeometry(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{

   WALBERLA_LOG_INFO_ON_ROOT("addSweepSurfaceGeometry");
   
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
                  << AfterFunction(Communication_T(simObjs->blockForest, simObjs->smoothFillFieldID), "Communication: after smoothing sweep");

   // compute interface normals (using smoothed fill level field)
   walberla::free_surface::NormalSweep< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T >
      normalSweep(simObjs->normalFieldID, simObjs->smoothFillFieldID, simObjs->flagFieldID, 
                  walberla::free_surface::flagIDs::interfaceFlagID, walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, 
                  simObjs->flagInfo.getObstacleIDSet(), true, false, true, false);

   timeloop->add() << Sweep(normalSweep, "Sweep: normal computation", StateSweep::fullFreeSurface)
                  << Sweep(EmptySweep(), "Empty sweep: normal")
                  << AfterFunction(Communication_T(simObjs->blockForest, simObjs->normalFieldID), "Communication: after normal sweep");

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with no surface tension
   const bool computeCurvature = !realIsEqual(simParams->computed.surfaceTension, real_c(0), real_c(1e-14));
   if (computeCurvature)
   {
      // compute interface curvature using finite differences according to Brackbill et al.
      walberla::free_surface::CurvatureSweepFiniteDifferences< CommunicationStencil_T, FlagField_T, ScalarField_T, VectorField_T >
         curvSweep(simObjs->curvatureFieldID, simObjs->normalFieldID, simObjs->obstacleNormalFieldID, simObjs->flagFieldID,
                   walberla::free_surface::flagIDs::interfaceFlagID,
                   walberla::free_surface::flagIDs::liquidInterfaceGasFlagIDs, 
                   simObjs->flagInfo.getObstacleIDSet(), false, real_c(0));

      timeloop->add() << Sweep(curvSweep, "Sweep: curvature computation (finite difference method)", StateSweep::fullFreeSurface)
                     << Sweep(EmptySweep(), "Empty sweep: curvature")
                     << AfterFunction(Communication_T(simObjs->blockForest, simObjs->curvatureFieldID), "Communication: after curvature sweep");
   }
}

void addSweepSurfaceDynamics(
   const std::shared_ptr< const SimulationParams >& simParams, 
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
   timeloop->add() << Sweep(forceDensityCodegenSweep, "Sweep: force weighting", Set< SUID >::emptySet(), StateSweep::onlyGasAndBoundary)
                  << Sweep(EmptySweep(), "Empty sweep: force weighting", StateSweep::onlyGasAndBoundary)
                  << AfterFunction(Communication_T(simObjs->blockForest, simObjs->forceDensityFieldID), "Communication: after force weighting sweep");

   // sweep for
   // - reconstruction of PDFs in interface cells
   // - streaming of PDFs in interface cells (and liquid cells on the same block)
   // - advection of mass
   // - update bubble volumes
   // - marking interface cells for conversion
   const walberla::free_surface::StreamReconstructAdvectSweep<
      LatticeModel_T, typename AntidunesBoundaryHandling_T::BoundaryHandling_T, FlagField_T,
      typename AntidunesBoundaryHandling_T::FlagInfo_T, ScalarField_T, VectorField_T, true >
      streamReconstructAdvectSweep(simParams->computed.surfaceTension, simObjs->antidunesBoundaryHandling->getHandlingID(), 
                                   simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID, simObjs->normalFieldID, 
                                   simObjs->curvatureFieldID, simObjs->flagInfo, simObjs->bubbleModel.get(), simParams->model.pdfReconstructionModel, 
                                   simParams->model.useSimpleMassExchange, simParams->model.cellConversionThreshold,
                                   simParams->model.cellConversionForceThreshold);

   // sweep acts only on blocks with at least one interface cell (due to StateSweep::fullFreeSurface)
   timeloop->add() << Sweep(streamReconstructAdvectSweep, "Sweep: StreamReconstructAdvect", StateSweep::fullFreeSurface)
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
      CollideSweep(const typename LatticeModel_T::Sweep& sweep) : sweep_(sweep){};

      void operator()(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0))
      {
         sweep_.collide(block, numberOfGhostLayersToInclude);
      }

    private:
      typename LatticeModel_T::Sweep sweep_;
   };

   // Collision for fullFreeSurface (because stream was done above) + StreamCollide for onlyLBM blocks
   timeloop->add() << Sweep(CollideSweep(lbmSweepGenerated), "Sweep: collision (generated)", StateSweep::fullFreeSurface)
                  << Sweep(lbmSweepGenerated, "Sweep: streamCollide (generated)", StateSweep::onlyLBM)
                  << Sweep(EmptySweep(), "Empty sweep: streamCollide (generated)")
                  << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID), "Communication: after streamCollide (generated)");

   // convert cells
   // - according to the flags from StreamReconstructAdvectSweep (interface -> gas/liquid)
   // - to ensure a closed layer of interface cells (gas/liquid -> interface)
   // - detect and register bubble merges/splits (bubble volumes are already updated in StreamReconstructAdvectSweep)
   // - convert cells and initialize PDFs near inflow boundaries
   const walberla::free_surface::CellConversionSweep< LatticeModel_T, AntidunesBoundaryHandling_T::BoundaryHandling_T,
                                                      ScalarField_T >
      cellConvSweep(simObjs->antidunesBoundaryHandling->getHandlingID(), simObjs->pdfFieldID, simObjs->flagInfo, simObjs->bubbleModel.get());
   timeloop->add() << Sweep(cellConvSweep, "Sweep: cell conversion", StateSweep::fullFreeSurface)
                  << Sweep(EmptySweep(), "Empty sweep: cell conversion")
                  //<< AfterFunction(PdfCommunication_T(blockForest, pdfFieldID),
                  //
                  //                 "Communication: after cell conversion sweep (PDF field)")
                  // communicate the flag field also in corner directions
                  << AfterFunction(Communication_T(simObjs->blockForest, simObjs->flagFieldID),
                                   "Communication: after cell conversion sweep (flag field)")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(simObjs->blockForest, simObjs->flagFieldID),
                                   "Second ghost layer update: after cell conversion sweep (flag field)");

   // reinitialize PDFs, i.e., refill cells that were converted from gas to interface
   // - when the flag "convertedFromGasToInterface" has been set (by CellConversionSweep)
   // - according to the method specified with pdfRefillingModel_
   const walberla::free_surface::EquilibriumRefillingSweep< LatticeModel_T, FlagField_T > 
      equilibriumRefillingSweep(simObjs->pdfFieldID, simObjs->flagFieldID, simObjs->flagInfo, true);
   timeloop->add() << Sweep(equilibriumRefillingSweep, "Sweep: EquilibriumRefilling", StateSweep::fullFreeSurface)
                  << Sweep(EmptySweep(), "Empty sweep: EquilibriumRefilling")
                  << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID), "Communication: after EquilibriumRefilling sweep");

   const walberla::free_surface::ExcessMassDistributionModel excessMassModel(simParams->model.excessMassDistributionModel);
   // distribute excess mass:
   // - excess mass: mass that is free after conversion from interface to gas/liquid cells
   // - update the bubble model
   // IMPORTANT REMARK: this sweep computes the mass via the density, i.e., the PDF field must be up-to-date and the
   // PdfRefillingSweep must have been performed
   if (excessMassModel.isEvenlyType())
   {
      const walberla::free_surface::ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
         distributeMassSweep(excessMassModel, simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID, simObjs->flagInfo);
      timeloop->add()
         << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
         << Sweep(EmptySweep(), "Empty sweep: distribute excess mass")
         << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID),
                          "Communication: after excess mass distribution sweep")
         << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest, simObjs->fillFieldID),
                          "Second ghost layer update: after excess mass distribution sweep (fill level field)")
         // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits are
         // already detected and registered by CellConversionSweep
         << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); }, "Sweep: bubble model update");
   }
   else
   {
      if (excessMassModel.isWeightedType())
      {
         const walberla::free_surface::ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            distributeMassSweep(excessMassModel, simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID, simObjs->flagInfo, simObjs->normalFieldID);
         timeloop->add()
            << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
            << Sweep(EmptySweep(), "Empty sweep: distribute excess mass")
            << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID),
                             "Communication: after excess mass distribution sweep")
            << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest, simObjs->fillFieldID),
                             "Second ghost layer update: after excess mass distribution sweep (fill level field)")
            // update bubble model, i.e., perform registered bubble merges/splits; bubble merges/splits
            // are already detected and registered by CellConversionSweep
            << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); }, "Sweep: bubble model update");
      }
      else
      {
         if (excessMassModel.isEvenlyAllInterfaceFallbackLiquidType())
         {
            const walberla::free_surface::ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
               distributeMassSweep(excessMassModel, simObjs->fillFieldID, simObjs->flagFieldID, simObjs->pdfFieldID, simObjs->flagInfo, simObjs->excessMassFieldID);
            timeloop->add()
               // perform this sweep also on "onlyLBM" blocks because liquid cells also exchange excess mass here
               << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::fullFreeSurface)
               << Sweep(distributeMassSweep, "Sweep: excess mass distribution", StateSweep::onlyLBM)
               << Sweep(EmptySweep(), "Empty sweep: distribute excess mass")
               << AfterFunction(Communication_T(simObjs->blockForest, simObjs->fillFieldID, simObjs->excessMassFieldID),
                                "Communication: after excess mass distribution sweep")
               << AfterFunction(blockforest::UpdateSecondGhostLayer< ScalarField_T >(simObjs->blockForest, simObjs->fillFieldID),
                                "Second ghost layer update: after excess mass distribution sweep (fill level field)")
               // update bubble model, i.e., perform registered bubble merges/splits; bubble
               // merges/splits are already detected and registered by CellConversionSweep
               << AfterFunction([&simObjs]() { simObjs->bubbleModel->update(); }, "Sweep: bubble model update");
         }
      }
   }

   // reset all flags that signal cell conversions (except "keepInterfaceForWettingFlag")
   walberla::free_surface::ConversionFlagsResetSweep< FlagField_T > 
      resetConversionFlagsSweep(simObjs->flagFieldID, simObjs->flagInfo);
   timeloop->add() << Sweep(resetConversionFlagsSweep, "Sweep: conversion flag reset", StateSweep::fullFreeSurface)
                  << Sweep(EmptySweep(), "Empty sweep: conversion flag reset")
                  << AfterFunction(Communication_T(simObjs->blockForest, simObjs->flagFieldID),
                                   "Communication: after excess mass distribution sweep")
                  << AfterFunction(blockforest::UpdateSecondGhostLayer< FlagField_T >(simObjs->blockForest, simObjs->flagFieldID),
                                   "Second ghost layer update: after excess mass distribution sweep (flag field)");

   // update block states
   auto blockStateUpdate = StateSweep(simObjs->blockForest, simObjs->flagInfo, simObjs->flagFieldID);
   timeloop->add() << Sweep(blockStateUpdate, "Sweep: block state update");
}

void addSweepFluidParticleCoupling(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addSweepFluidParticleCoupling");

   // TODO: maybe put this FlagUID in the antidunesBoundaryHandling
   const FlagUID FormerMO_Flag("former moving obstacle");

   auto sphereSelector = [](const size_t particleIdx, const auto& ac)
   {
      return ac.getBaseShape(particleIdx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE;
   };

   // sweep for mapping particles into the LBM cells
   const bool strictlyConserveMomentum = false;
   timeloop->add() << Sweep(
      lbm_mesapd_coupling::makeMovingParticleMapping< PdfField_T, AntidunesBoundaryHandling_T::BoundaryHandling_T >(
         simObjs->blockForest, simObjs->pdfFieldID, simObjs->antidunesBoundaryHandling->getHandlingID(),
         simObjs->particleFieldID, simObjs->particleAccessor,
         AntidunesBoundaryHandling_T::movingObstacleFlagID,
         FormerMO_Flag,
         sphereSelector,   // map all spheres, regardless if they are fixed or not
         strictlyConserveMomentum),
      "Sweep (coupling): MovingParticleMapping");

   // sweep for restoring PDFs in cells previously occupied by particles
   auto gradReconstructor =
      lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<AntidunesBoundaryHandling_T::BoundaryHandling_T>(
         simObjs->blockForest, simObjs->antidunesBoundaryHandling->getHandlingID(),
         simParams->computed.relaxationRate,
         false,   // recomputeTargetDensity
         true     // useCentralDifferences
      );

   timeloop->add() << Sweep(makeSharedSweep(
            lbm_mesapd_coupling::makePdfReconstructionManager< PdfField_T, AntidunesBoundaryHandling_T::BoundaryHandling_T >(
               simObjs->blockForest, simObjs->pdfFieldID, simObjs->antidunesBoundaryHandling->getHandlingID(),
               simObjs->particleFieldID, simObjs->particleAccessor,
               FormerMO_Flag, walberla::free_surface::flagIDs::liquidFlagID,
               gradReconstructor, strictlyConserveMomentum)),
            "Sweep (coupling): PDF Reconstruction (moving particles)")
      << AfterFunction(Communication_T(simObjs->blockForest, simObjs->flagFieldID, simObjs->particleFieldID),
                       "Communication (coupling): after PDF reconstruction (flagField + particleField)") // unsure if necessary but added for consistency
      << AfterFunction(PdfCommunication_T(simObjs->blockForest, simObjs->pdfFieldID),
                       "Communication (coupling): after PDF reconstruction (pdfField)");
}

void addAfterTimeStepPerformanceLogger(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepPerformanceLogger");

   lbm::PerformanceLogger< FlagField_T > 
      performanceLogger(
         simObjs->blockForest, 
         simObjs->flagFieldID,
         walberla::free_surface::flagIDs::liquidInterfaceFlagIDs,
         simParams->evaluation.performanceLogFrequency);
   timeloop->addFuncAfterTimeStep(performanceLogger, "Evaluator: performance logging");
}

void addAfterTimeStepStabilityChecker(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepStabilityChecker");
   
   // logs data and aborts if it detects non-finite values in liquid or interfce cells
   timeloop->addFuncAfterTimeStep(makeSharedFunctor(
      field::makeStabilityChecker<PdfField_T, FlagField_T>(
         simParams->cfgPtr, 
         simObjs->blockForest, simObjs->pdfFieldID, simObjs->flagFieldID,
         walberla::free_surface::flagIDs::liquidInterfaceFlagIDs)),
      "Evaluator: LBM stability check");
}

void addAfterTimeStepDataEvaluators(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepDataEvaluators");

   // data slice evaluator --------------------------------------------------------------------------------------
   // writes data to file (solidVolumeFraction + fill + velocityX) averaged over X-Z slices in the domain
   AverageDataSliceEvaluator< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
      averageDataSliceEvaluator(
         simObjs->blockForest, 
         simObjs->flagFieldID, 
         simObjs->fillFieldID, 
         simObjs->pdfFieldID,
         simParams->evaluation.baseFolderName, 
         simObjs->dataCols,
         simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(averageDataSliceEvaluator, "Evaluator: AverageDataSlice");

   // total mass evaluator --------------------------------------------------------------------------------------
   walberla::free_surface::TotalMassComputer< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T >
      totalFluidMassEvaluator(
         simObjs->blockForest, 
         simObjs->antidunesBoundaryHandling, 
         simObjs->pdfFieldID,
         simObjs->fillFieldID, 
         simParams->evaluation.evaluationFrequency,
         simObjs->totalFluidMass);
   timeloop->addFuncAfterTimeStep(totalFluidMassEvaluator, "Evaluator: TotalFluidMass");

   // bed load transport evaluator ------------------------------------------------------------------------------
   const real_t bedLoadTransportAveragingFactor = real_c(1) / real_c(simParams->domain.domainSize[0] * simParams->domain.domainSize[1]);
   BedloadTransportEvaluator< ParticleAccessor_T > 
      bedloadTransportEvaluator(
         simObjs->particleAccessor, 
         bedLoadTransportAveragingFactor, 
         simObjs->numParticles,
         simObjs->dataCols,
         simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(bedloadTransportEvaluator, "Evaluator: BedloadTransport");

   // total hydro force evaluator ------------------------------------------------------------------------------
   TotalHydrodynamicForceEvaluator < ParticleAccessor_T > 
      totalHydrodynamicForceEvaluator(
         simObjs->particleAccessor,
         simObjs->dataCols, 
         simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(totalHydrodynamicForceEvaluator, "Evaluator: TotalHydrodynamicForce");
}

void addAfterTimeStepForcingAdjuster(
   const std::shared_ptr< const SimulationParams >& simParams, 
   const std::shared_ptr< SimulationObjects >& simObjs,
   const std::unique_ptr< SweepTimeloop >& timeloop)
{
   WALBERLA_LOG_INFO_ON_ROOT("addAfterTimeStepForcingAdjuster");

   // compute fluid's mean velocity -----------------------------------------------------------------------------
   const real_t velocityAveragingFactor = real_c(1) / (simParams->computed.liquidHeight * 
                                          real_c(simParams->domain.domainSize[0] * simParams->domain.domainSize[1]));
   MeanVelocityComputer< AntidunesBoundaryHandling_T, PdfField_T, FlagField_T > 
      meanVelocityComputer(
         simObjs->blockForest,
         simObjs->antidunesBoundaryHandling,
         simObjs->pdfFieldID,
         simObjs->meanVelocity,
         velocityAveragingFactor,
         simParams->evaluation.evaluationFrequency);
   timeloop->addFuncAfterTimeStep(meanVelocityComputer, "Evaluator: MeanVelocityComputer");

   // update acceleration ---------------------------------------------------------------------------------------
   // using the computed mean velocity, adjust the acceleration (in X dir) such that the target mean velocity is maintained
   ForcingAdjuster forcingAdjuster(
      simObjs->meanVelocity,
      simObjs->acceleration,
      simObjs->pidController,
      simParams->evaluation.evaluationFrequency);

   timeloop->addFuncAfterTimeStep(forcingAdjuster, "Evaluator: ForcingAdjuster");
}


} // namespace antidunes
} // namespace walberla