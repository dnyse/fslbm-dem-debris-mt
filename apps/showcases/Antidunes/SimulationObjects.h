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
//! \file SimulationObjects.h
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "lbm_mesapd_coupling/DataTypes.h"

#include "mesa_pd/data/ParticleAccessorWithBaseShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"

#include "AntidunesBoundaryHandling.h"
#include "AntidunesLatticeModel.h"
#include "DataCols.h"
#include "PIDController.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace walberla {
namespace antidunes {

// Convenience aliases ------------------------------------------------------------------
using ScalarField_T                 = GhostLayerField< real_t, 1 >;
using VectorField_T                 = GhostLayerField< Vector3< real_t >, 1 >;
using VectorFieldFlattened_T        = GhostLayerField< real_t, 3 >;
using LatticeModel_T                = lbm::AntidunesLatticeModel;
using PdfField_T                    = lbm::PdfField< LatticeModel_T >;
using ParticleField_T               = lbm_mesapd_coupling::ParticleField_T;
using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using ParticleAccessor_T            = mesa_pd::data::ParticleAccessorWithBaseShape;
using AntidunesBoundaryHandling_T   = free_surface::AntidunesBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T, ParticleAccessor_T >;
using FlagInfo_T                    = AntidunesBoundaryHandling_T::FlagInfo_T;
using CommunicationStencil_T        = std::conditional_t< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >;
using Communication_T               = blockforest::SimpleCommunication< CommunicationStencil_T >;
using PdfCommunication_T            = blockforest::SimpleCommunication< LatticeModel_T::Stencil >;

// Single struct that holds all simulation objects and travels through every simulation stage
struct SimulationObjects
{
   // BlockForest -----------------------------------------------------------------------------------------------------
   shared_ptr< StructuredBlockForest > blockForest;   const std::string blockForestFile{ "blockForest.file" };

   // LBM fields ------------------------------------------------------------------------------------------------------
   BlockDataID pdfFieldID;                            const std::string pdfFieldFile{ "pdfField.file" };
   BlockDataID fillFieldID;                           const std::string fillFieldFile{ "fillField.file" };
   BlockDataID excessMassFieldID;                     const std::string excessMassFieldFile{ "excessMassField.file" };
   BlockDataID forceDensityFieldID;

   // Surface geometry fields -----------------------------------------------------------------------------------------
   BlockDataID curvatureFieldID;
   BlockDataID normalFieldID;
   BlockDataID obstacleNormalFieldID;
   BlockDataID smoothFillFieldID;

   // Particles -------------------------------------------------------------------------------------------------------
   BlockDataID particleStorageID;                     const std::string particleStorageFile{ "particleStorage.file" };
   BlockDataID particleFieldID;                       // Field for mapping particle UIDs to cells
   shared_ptr< mesa_pd::data::ParticleStorage >                particleStorage;
   shared_ptr< mesa_pd::data::ParticleAccessorWithBaseShape >  particleAccessor;
   shared_ptr< mesa_pd::domain::BlockForestDomain >            mesapdDomain;
   std::function< void() >                                     syncCall;
   uint_t                                                      numParticles;
   // floaterUIDs[fi] = UID of floating sphere fi (matches simParams->floatingSpheres[fi]).
   // floaterBoxUIDs[bi] = UID of floating box bi (matches simParams->floatingBoxes[bi]).
   // floaterCylinderUIDs[ci] = UID of floating cylinder ci (matches simParams->floatingCylinders[ci]).
   // floaterCapsuleUIDs[pi] = UID of floating capsule pi (matches simParams->floatingCapsules[pi]).
   // Filled by setupParticles (indexed creation + allReduce(SUM)), so every rank knows every floater's UID;
   // floaterCapsuleUIDs only by the DropTest variant (SetupStages_drop_test.h).
   std::vector< walberla::id_t >                               floaterUIDs;
   std::vector< walberla::id_t >                               floaterBoxUIDs;
   std::vector< walberla::id_t >                               floaterCylinderUIDs;
   std::vector< walberla::id_t >                               floaterCapsuleUIDs;

   // UIDs of dormant floaters (startInactive, before settlingSteps): excluded from fluid mapping
   // (SetupStages::floaterSelector, BuildTimeloop::movingObstacleSelector) and DEM contact until
   // releaseSettledFloatingSpheres removes them. Shared because SetupStages.h and BuildTimeloop.h need it before
   // ParticleKernels is constructed.
   std::shared_ptr< std::unordered_set< walberla::id_t > >     dormantFloaterUIDs =
      std::make_shared< std::unordered_set< walberla::id_t > >();

   // UIDs of floaters released in this main step (SetupStages::releaseSettledFloatingSpheres). ComputeContacts
   // (ParticleKernels.h) aborts with diagnostics if one of them is in genuine contact at sub-cycle 0, i.e. was
   // released on top of an active neighbour; contacts at later sub-cycles are only logged. Cleared at the end of
   // every doDEMSubCycles(). Debug/validation aid.
   std::shared_ptr< std::unordered_set< walberla::id_t > >     freshlyActivatedFloaterUIDs =
      std::make_shared< std::unordered_set< walberla::id_t > >();

   // Boundary and flags ----------------------------------------------------------------------------------------------
   shared_ptr< AntidunesBoundaryHandling_T > antidunesBoundaryHandling;
   BlockDataID flagFieldID;
   FlagInfo_T flagInfo;

   // Diagnostic sweep event collector ---------------------------------------------------------------------------------
   // Append-only buffer for per-sub-cycle diagnostics in ParticleKernels.h (one push_back instead of logging at the
   // point of detection). BuildTimeloop.h's addAfterTimeStepDiagnosticSweep drains, logs and clears it once per main
   // step. Every push site is gated by ModelParameters::diagnosticSweepEnabled, so it stays empty when that is off.
   struct DiagnosticEvent
   {
      uint_t             mainStep{ uint_t(0) };
      uint_t             subCycle{ uint_t(0) };
      std::string        tag;                              // e.g. "VELOCITY-JUMP", "MIGRATION-NAN"
      walberla::id_t     uid{ 0 };                          // 0 if not tied to a specific floater
      Vector3< real_t >  worldPos{ real_t(-1) };            // (-1,-1,-1) sentinel if not applicable
      real_t             value{ real_t(0) };                // e.g. |deltaV|, density -- meaning is tag-specific
      std::string        detail;                            // free-form extra context for the log line
   };
   std::shared_ptr< std::vector< DiagnosticEvent > > diagnosticEvents =
      std::make_shared< std::vector< DiagnosticEvent > >();

   // bubbleModel -----------------------------------------------------------------------------------------------------
   shared_ptr< walberla::free_surface::bubble_model::BubbleModelBase > bubbleModel;

   // Simulation state ------------------------------------------------------------------------------------------------
   const std::string                simStateFile{ "simState.file" };
   const std::string                pidStateFile{ "pidState.file" };

   uint_t                           beginTimeStep;
   // Live main-loop step counter, set at the top of each iteration in Antidunes.cpp (before timeloopFSLBM runs), so
   // callbacks registered at setup time (e.g. windowed VTK output) can read the current step.
   shared_ptr< uint_t >             currentMainStep = std::make_shared< uint_t >(uint_t(0));
   shared_ptr< PIDController >      pidController;
   shared_ptr< Vector3< real_t > >  acceleration;
   shared_ptr< Vector3< real_t > >  meanVelocity;
   shared_ptr< real_t >             totalFluidMass;

   // Live MEM hydrostatic reference height (ModelParameters.useDynamicReferenceHeight): seeded to the static
   // reference height and EMA-updated from the floaters' fitted local waterlines in ParticleKernels'
   // computeSubmergedVolumes(). Read per boundary link by HydrostaticDensityFtor; force side only (SimpleBB.h applies
   // rhoHydStat to local PDF copies used for forceMEM). Unused when the toggle is off.
   shared_ptr< real_t >             memReferenceHeight;

   shared_ptr< DataCols >           dataCols;
};

} // namespace antidunes
} // namespace walberla