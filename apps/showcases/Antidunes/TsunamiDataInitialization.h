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
//! \file TsunamiDataInitialization.h
//
// Dry-domain / dam-break initial-condition machinery for tsunami.cpp -- the ONE piece of physics that
// is genuinely different from Antidunes.cpp's own initialization (see that file's DataInitialization.h):
// instead of a domain that is everywhere wet (sinusoidal free-surface height, Poiseuille shear profile,
// hydrostatic hydro pressure tuned for a nonzero periodic mean flow), here most of the domain starts
// completely DRY (gas), with liquid confined to a reservoir behind a virtual "gate" plane at x=gateX.
// Releasing the domain (there is no actual gate body -- the fill-level step itself IS the gate; FSLBM
// collapses it under gravity from the very first step, exactly like every dam-break-over-dry-bed
// benchmark in the free-surface-LBM literature) reproduces Watanabe et al. 2023's tsunami-flow
// experiments (Section 3/4 therein: a raised reservoir separated by a gate from a dry run-out floor,
// released under gravity toward a pole forest / driftwood field) MUCH more directly than treating the
// forest/driftwood as sitting in an already-flowing, already-wet steady channel (see this project's own
// prms/toy_pole_forest_driftwood_gas.prm header, which explicitly flagged that simplification and the
// gap this file closes).
//
// Everything else this scenario needs (pillars/cylinder pillars, arch tunnels, floating spheres/boxes/
// cylinders, DEM sub-cycling, MEM coupling, boundary handling, evaluators, checkpointing) is reused
// UNCHANGED from Antidunes' own SimulationParams.h/SimulationObjects.h/SetupStages.h/BuildTimeloop.h/
// ParticleKernels.h -- see tsunami.cpp's own header comment for exactly which functions are reused
// as-is vs. replaced.
//======================================================================================================================

#pragma once

#include "core/config/Config.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionModel.h"

// SimulationParams.h/SimulationObjects.h must come before DataInitialization.h -- that file (like
// SetupStages.h) uses SimulationParams/SimulationObjects/ScalarField_T/PdfField_T in its own function
// signatures without including either header itself, relying on whichever translation unit includes it
// to have already done so first (the same pre-existing ordering constraint Antidunes.cpp/SetupStages.h
// satisfy implicitly; see tsunami.cpp's own include order, which follows the same convention).
#include "SimulationParams.h"
#include "SimulationObjects.h"
#include "DataInitialization.h"

namespace walberla {
namespace antidunes {

// Parsed from a new, dedicated "DamBreakParameters" prm block (NOT part of SimulationParams -- that
// file is Antidunes.cpp's own and is left completely untouched, per this executable's own "existing
// working code must not change" constraint). Defines the initially-wet reservoir as an axis-aligned
// box: x in [0, gateX], y in [yMin, yMax] (default: the full domain width), z in [0, waterDepth].
// Everywhere else in the domain starts completely dry (fill level 0, gas).
struct DamBreakParams {
   real_t gateX;        // reservoir's downstream face (the "gate"); x >= gateX starts dry
   real_t waterDepth;    // reservoir's initial, flat water depth (z in [0, waterDepth) is liquid)
   real_t yMin;          // reservoir's spanwise extent; defaults to the full domain width (see below)
   real_t yMax;

   DamBreakParams(const config::Config::BlockHandle& b, real_t domainSizeY)
      : gateX      { b.getParameter< real_t >("gateX") }
      , waterDepth { b.getParameter< real_t >("waterDepth") }
      , yMin       { b.getParameter< real_t >("yMin", real_t(0)) }
      , yMax       { b.getParameter< real_t >("yMax", domainSizeY) }
   {
      WALBERLA_CHECK_GREATER(gateX, real_t(0));
      WALBERLA_CHECK_GREATER(waterDepth, real_t(0));
      WALBERLA_CHECK_GREATER(yMax, yMin);
   }

   // Direct-value constructor for callers with no prm block to parse from (unit tests) -- same
   // defaulting rule as the prm-parsing constructor above (yMin=0, yMax=domainSizeY).
   DamBreakParams(real_t gateXIn, real_t waterDepthIn, real_t yMinIn, real_t yMaxIn)
      : gateX(gateXIn), waterDepth(waterDepthIn), yMin(yMinIn), yMax(yMaxIn)
   {}
};

// fraction of the unit interval [lo, lo+1) that lies below `hi` -- e.g. a cell's overlap with the
// dry/wet step at the gate (x) or the free surface (z). Pure function (no waLBerla field/block-forest
// dependency) so it is directly unit-testable -- see tests/showcases_antidunes/TsunamiDamBreakTest.cpp.
inline real_t overlapBelow(real_t lo, real_t hi) { return std::clamp(hi - lo, real_t(0), real_t(1)); }

// fraction of the unit interval [lo, lo+1) that lies within [rangeLo, rangeHi) -- e.g. a cell's overlap
// with the reservoir's spanwise (y) extent. Same pure-function rationale as overlapBelow above.
inline real_t overlapWithin(real_t lo, real_t rangeLo, real_t rangeHi)
{
   const real_t coveredLo = std::max(lo, rangeLo);
   const real_t coveredHi = std::min(lo + real_t(1), rangeHi);
   return std::clamp(coveredHi - coveredLo, real_t(0), real_t(1));
}

// Exact fractional overlap of the unit cell whose lower corner is (globalX, globalY, globalZ) with the
// reservoir box (x in [0, gateX), y in [yMin, yMax), z in [0, waterDepth)) -- a closed-form product of
// three per-axis clamped-overlap fractions (exact because the reservoir is itself an axis-aligned box,
// unlike Antidunes' own sinusoidal free surface, which needs Monte-Carlo sampling). This is the pure,
// pull-out-of-the-field-loop core of initFillFieldDamBreak below, kept separate specifically so it is
// unit-testable without a live blockforest/ScalarField_T (same rationale as e.g. SimulationParams.h's
// archDeckFillBoxGeometry/archBlockParams) -- see TsunamiDamBreakTest.cpp.
inline real_t damBreakFillFraction(real_t globalX, real_t globalY, real_t globalZ, const DamBreakParams& damBreak)
{
   const real_t fracX = overlapBelow(globalX, damBreak.gateX);
   const real_t fracY = overlapWithin(globalY, damBreak.yMin, damBreak.yMax);
   const real_t fracZ = overlapBelow(globalZ, damBreak.waterDepth);
   return fracX * fracY * fracZ;
}

// Analytic (not Monte-Carlo -- see damBreakFillFraction above) dam-break fill-level initializer. Cells
// fully inside the reservoir box get fill=1, cells fully outside get fill=0, and cells straddling the
// gate (x) or the free surface (z) get the exact partial fill of a step function, giving FSLBM a clean,
// alias-free initial interface on both the vertical gate face and the horizontal free surface.
void initFillFieldDamBreak(const std::shared_ptr< SimulationObjects >& simObjs, const DamBreakParams& damBreak)
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing dam-break fill field: reservoir x in [0, " << damBreak.gateX
      << "), y in [" << damBreak.yMin << ", " << damBreak.yMax << "), z in [0, " << damBreak.waterDepth
      << ") -- everywhere else starts dry.");

   auto& bf = simObjs->blockForest;

   for (auto blockIt = bf->begin(); blockIt != bf->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(simObjs->fillFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         const Cell localCell = fillFieldIt.cell();
         Cell globalCell = localCell;
         bf->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

         fillField->get(localCell) =
            damBreakFillFraction(real_c(globalCell[0]), real_c(globalCell[1]), real_c(globalCell[2]), damBreak);
      }) // WALBERLA_FOR_ALL_CELLS
   }
}

// Still water at rest: rho=1 (uniform pressure) equilibrium, zero velocity, in EVERY cell (wet or dry
// -- a dry/gas cell's PDF values are physically meaningless until FSLBM converts it, but must still
// hold a valid equilibrium so streaming/collision never reads garbage). Deliberately NOT the hydrostatic
// pressure profile Antidunes.cpp uses (initPdfFieldHydrostaticPressure, DataInitialization.h): that
// profile is tuned to pre-condition a channel that already has a nonzero periodic mean flow, which
// would inject a spurious initial pressure gradient into a domain meant to start from genuine rest. A
// uniform-density start is the standard dam-break-over-dry-bed initial condition in the free-surface-LBM
// literature; the correct hydrostatic gradient develops naturally within the first few steps as gravity
// acts on the reservoir's own step-function free surface.
void initPdfFieldStillWater(const std::shared_ptr< SimulationObjects >& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing still water at rest (rho=1, u=0 everywhere)");
   auto& bf = simObjs->blockForest;

   for (auto blockIt = bf->begin(); blockIt != bf->end(); ++blockIt)
   {
      PdfField_T* const pdfField = blockIt->getData< PdfField_T >(simObjs->pdfFieldID);
      WALBERLA_FOR_ALL_CELLS_XYZ(pdfField, pdfField->setToEquilibrium(x, y, z, Vector3< real_t >(real_t(0)), real_t(1));)
   }
}

// Same field-allocation/snapshot-loading structure as SetupStages.h::setupLBMFields, but with the
// non-snapshot (fresh-start) branch calling initFillFieldDamBreak/initPdfFieldStillWater above instead
// of Antidunes.cpp's own initFillFieldSinus/initPdfFieldPoiseuille/initPdfFieldHydrostaticPressure --
// see this file's own header comment for why those two are the wrong initial condition for a dry-domain
// dam break. Duplicated (rather than parameterizing the original) specifically so
// SetupStages.h::setupLBMFields itself never has to change -- this executable's whole point is that the
// already-working Antidunes code path stays untouched.
void setupLBMFieldsDamBreak(const std::shared_ptr< const SimulationParams >& simParams,
                             const std::shared_ptr< SimulationObjects >& simObjs,
                             const DamBreakParams& damBreak)
{
   auto& bf              = simObjs->blockForest;
   const auto& bfParams  = simParams->blockForest;
   const walberla::free_surface::ExcessMassDistributionModel excessMassModel(simParams->model.excessMassDistributionModel);

   simObjs->forceDensityFieldID = field::addToStorage< VectorFieldFlattened_T >(bf, "Force field", real_c(0), field::fzyx, uint_c(1));
   LatticeModel_T latticeModel(simObjs->forceDensityFieldID, simParams->computed.relaxationRate);

   if (bfParams.loadSnapshot)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Loading LBM fields from snapshot");
      const std::string base = bfParams.snapshotBaseFolder + "/";
      auto pdfHandling = make_shared< lbm::internal::PdfFieldHandling< LatticeModel_T > >(bf, latticeModel, false, Vector3< real_t >(real_c(0)), real_c(1), uint_c(1), field::fzyx);
      simObjs->pdfFieldID = bf->getBlockStorage().loadBlockData(base + simObjs->pdfFieldFile, pdfHandling, "PDF field");

      auto fillHandling = std::make_shared< ScalarFieldHandling< ScalarField_T > >(bf, uint_c(2));
      simObjs->fillFieldID = bf->getBlockStorage().loadBlockData(base + simObjs->fillFieldFile, fillHandling, "Fill level field");

      if (excessMassModel.isEvenlyAllInterfaceFallbackLiquidType())
      {
         auto excessHandling = std::make_shared< ScalarFieldHandling< ScalarField_T > >(bf, uint_c(1));
         simObjs->excessMassFieldID = bf->getBlockStorage().loadBlockData(base + simObjs->excessMassFieldFile, excessHandling, "Excess mass field");
      }
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Creating LBM fields");
      simObjs->pdfFieldID = lbm::addPdfFieldToStorage(bf, "PDF field", latticeModel, Vector3< real_t >(real_t(0)), real_c(1.0), field::fzyx);
      simObjs->fillFieldID = field::addToStorage< ScalarField_T >(bf, "Fill level field", real_c(0), field::fzyx, uint_c(2));

      if (excessMassModel.isEvenlyAllInterfaceFallbackLiquidType())
      {
         simObjs->excessMassFieldID = field::addToStorage< ScalarField_T >(bf, "Excess mass field", real_c(0), field::fzyx, uint_c(1));
      }

      initFillFieldDamBreak(simObjs, damBreak);
      initPdfFieldStillWater(simObjs);
   }
}

// SetupStages.h::setupParticles (UNCHANGED, per this executable's own constraint) only ever creates
// mesa_pd wall planes for Z (always) and Y (conditional on !periodicity[1]) -- never for X, regardless
// of periodicity[0]. setupBoundaryHandling's own auto-added NoSlip wall for non-periodic X (SetupStages.h,
// right after initFromConfig) is FLUID-side only (the flag field); it gives no DEM collision response,
// so a free-floating DEM body (a FloatingBox log, here) has nothing physically stopping it at x=0 or
// x=domainSize[0] the way it already does at the non-periodic Y walls. Every existing Antidunes scenario
// either keeps periodicity[0]=true (the antidune channel family) or keeps its floaters well clear of the
// X boundaries, so this gap was never actually exercised before -- tsunami_pole_forest_driftwood_gas.prm
// is the first scenario where a wave-driven floater could plausibly reach x=0 or x=128. Call this once,
// right after setupParticles(simParams, simObjs), to add the missing DEM planes for THIS executable only
// -- SetupStages.h itself stays untouched.
void addTsunamiXWallPlanes(const std::shared_ptr< const SimulationParams >& simParams,
                           const std::shared_ptr< SimulationObjects >& simObjs)
{
   if (simParams->blockForest.periodicity[0]) return;

   WALBERLA_LOG_INFO_ON_ROOT("Adding DEM wall planes at x=0 and x=domainSize[0] (periodicity[0]=false, "
      "see addTsunamiXWallPlanes's own comment for why SetupStages.h::setupParticles doesn't already do this)");
   auto& bf = simObjs->blockForest;
   createPlane(*simObjs->particleStorage, bf->getDomain().minCorner(), Vector3< real_t >( 1, 0, 0));
   createPlane(*simObjs->particleStorage, bf->getDomain().maxCorner(), Vector3< real_t >(-1, 0, 0));
}

} // namespace antidunes
} // namespace walberla
