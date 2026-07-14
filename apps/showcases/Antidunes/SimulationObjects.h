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

#include <unordered_map>

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

   // Surface geomtery fields -----------------------------------------------------------------------------------------
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
   // All populated by setupParticles via an indexed creation loop + allReduce(SUM) so that every
   // rank holds the correct UID for every floater regardless of which rank owns each particle.
   std::vector< walberla::id_t >                               floaterUIDs;
   std::vector< walberla::id_t >                               floaterBoxUIDs;
   std::vector< walberla::id_t >                               floaterCylinderUIDs;
    
   // Boundary and flags ----------------------------------------------------------------------------------------------
   shared_ptr< AntidunesBoundaryHandling_T > antidunesBoundaryHandling;
   BlockDataID flagFieldID;
   FlagInfo_T flagInfo;

   // bubbleModel -----------------------------------------------------------------------------------------------------
   shared_ptr< walberla::free_surface::bubble_model::BubbleModelBase > bubbleModel;

   // Simulation state ------------------------------------------------------------------------------------------------
   const std::string                simStateFile{ "simState.file" };
   const std::string                pidStateFile{ "pidState.file" };

   uint_t                           beginTimeStep;
   shared_ptr< PIDController >      pidController;
   shared_ptr< Vector3< real_t > >  acceleration;
   shared_ptr< Vector3< real_t > >  meanVelocity;
   shared_ptr< real_t >             totalFluidMass;

   // Live MEM hydrostatic reference height (ModelParameters.useDynamicReferenceHeight). Seeded to the flat
   // absoluteLiquidHeight and EMA-updated from floaters' fitted local waterlines in ParticleKernels'
   // computeSubmergedVolumes(). Read fresh per boundary link by HydrostaticDensityFtor — force-side only,
   // never touches the fluid/PDF fields (see SimpleBB.h: rhoHydStat only modifies local pdf copies used for
   // forceMEM, after the field write). Null / unused when the toggle is off (static flat H).
   shared_ptr< real_t >             memReferenceHeight;

   // Per-floater free-surface plane fit used for the buoyancy/submerged-volume calculation, keyed by
   // particle UID. Populated every step by ParticleKernels' computeSubmergedVolumes() (spheres) and
   // computeSubmergedVolumesBox() (box/cylinder floaters). Diagnostic/visualization only (VTK output
   // in BuildTimeloop.h, identifier "floaterPlane") — never read back into the physics. floaterPlaneOrigin
   // is a point on the fitted plane (xBar); floaterPlaneNormal is the fitted free-surface normal
   // (liquid->gas). A floater with no entry, or a zero-length normal, had no clean plane fit this step
   // (fully submerged/emerged -> nothing meaningful to draw).
   shared_ptr< std::unordered_map< walberla::id_t, Vector3< real_t > > > floaterPlaneOrigin;
   shared_ptr< std::unordered_map< walberla::id_t, Vector3< real_t > > > floaterPlaneNormal;

   shared_ptr< DataCols >           dataCols;
};

} // namespace antidunes
} // namespace walberla