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
//! \file DEMValidationScenarios.h
//! \brief Fluid-free DEM validation experiments for DEMValidation.cpp.
//!
//! Split into 4 files by scenario family (this header is now a thin aggregator,
//! kept so DEMValidation.cpp's own #include "DEMValidationScenarios.h" needs no change):
//!   - DEMValidationInfra.h: DEMEngine, PassFailAccumulator, shared logging/loop helpers.
//!   - DEMValidationWallImpact.h: single-body wall-impact scenarios.
//!   - DEMValidationPairEnergy.h: pairwise central head-on restitution/energy scenarios.
//!   - DEMValidationPackingAndDrum.h: packing/angle-of-repose/rotating-drum scenarios.
//! (DEMValidationCornerAndMultiPoint.h and the collision replay were archived on 2026-11-20, see
//! archive/unused_contact_methods/.)
//!
//! Every scenario runs the PRODUCTION contact/integration code (SolidCylinderContactDetection,
//! ComputeContacts, mesa_pd::kernel::VelocityVerletPre/PostForceUpdate -- ParticleKernels.h) through a
//! minimal, hand-written sub-cycle loop with no fluid coupling (no lubrication, no hydrodynamic force,
//! no buoyancy -- gravity, where used, is applied directly via addForceAtomic, same as BedGeneration.cpp).
//! Existing tests (ContactDetectionTest.cpp, PaperAnalyticalContactTest.cpp, validation/*.py) only check
//! STATIC contact geometry; this file is the first thing in the repo that runs real dynamics and checks
//! the result against a closed-form or literature reference -- see AnalyticalSolutions.h for the physics.
//======================================================================================================================

#pragma once

#include "DEMValidationInfra.h"
#include "DEMValidationWallImpact.h"
#include "DEMValidationPairEnergy.h"
#include "DEMValidationPackingAndDrum.h"
