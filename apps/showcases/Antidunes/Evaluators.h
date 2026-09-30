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
//! \file   Evaluators.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "Utility.h"
#include "DataCols.h"

namespace walberla
{
namespace antidunes
{

// evaluates X-Z slices of solid volume fraction, fill level and X-velocity
template< typename AntidunesBoundaryHandling_T, typename PdfField_T, typename FlagField_T, typename ScalarField_T >
class AverageDataSliceEvaluator
{
 public:
   AverageDataSliceEvaluator(const shared_ptr< StructuredBlockStorage >& blocks, const ConstBlockDataID& flagFieldID,
                             const ConstBlockDataID& fillFieldID, const ConstBlockDataID& pdfFieldID, 
                             const std::string& baseFolderName, const shared_ptr< DataCols >& dataCols, uint_t frequency)
      : blocks_(blocks), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID), pdfFieldID_(pdfFieldID), 
      baseFolderName_(baseFolderName), dataCols_(dataCols), frequency_(frequency), executionCounter_(uint_c(0))
   {
      xlen_ = blocks_->getNumberOfXCells();
      ylen_ = blocks_->getNumberOfYCells();
      zlen_ = blocks_->getNumberOfZCells();

      x_z_SolidVolumeFraction_ = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_FillLevel_           = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_VelocityX_           = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_FluidCellCount_      = std::vector< uint_t >(xlen_ * zlen_, 0);
      maxFluidZPos_            = uint_t(0);
   }

   void operator()(){
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {  
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   void compute()
   {
      // erase data
      std::fill(x_z_SolidVolumeFraction_.begin(), x_z_SolidVolumeFraction_.end(), real_t(0));
      std::fill(x_z_FillLevel_.begin(), x_z_FillLevel_.end(), real_t(0));
      std::fill(x_z_VelocityX_.begin(), x_z_VelocityX_.end(), real_t(0));
      std::fill(x_z_FluidCellCount_.begin(), x_z_FluidCellCount_.end(), 0);

      // fill contributions
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const PdfField_T* const pdfField     = block->getData< const PdfField_T >(pdfFieldID_);
         const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);
         const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);

         const auto solidMO     = flagField->getFlag(AntidunesBoundaryHandling_T::movingObstacleFlagID);
         const auto solidNoSlip = flagField->getFlag(AntidunesBoundaryHandling_T::noSlipFlagID);

         CellInterval xyz = flagField->xyzSize();
         Cell globalCell;

         maxFluidZPos_ = uint_t(0);

         // iterate all (inner) cells in the field
         for (auto cell : xyz)
         {
            blocks_->transformBlockLocalToGlobalCell(globalCell, *block, cell);
            auto entryIdx = uint_c(globalCell.x()) + uint_c(globalCell.z()) * xlen_;
            if (flagField->isFlagSet(cell, solidMO) || flagField->isFlagSet(cell, solidNoSlip))
            {
               x_z_SolidVolumeFraction_[entryIdx] += real_t(1);
               x_z_FillLevel_[entryIdx] += real_t(1);
            }
            else
            {
               auto fillLevel = fillField->get(cell);
               x_z_FillLevel_[entryIdx] += fillLevel;
               if (fillLevel > 0_r)
               {
                  x_z_VelocityX_[entryIdx] += pdfField->getVelocity(cell)[0];
                  ++x_z_FluidCellCount_[entryIdx];

                  maxFluidZPos_ = std::max(uint_t(globalCell.z()), maxFluidZPos_);
               }
            }
         }
      }

      // reduce this information to the root process
      mpi::reduceInplace(x_z_SolidVolumeFraction_, mpi::SUM);
      mpi::reduceInplace(x_z_FillLevel_, mpi::SUM);
      mpi::reduceInplace(x_z_VelocityX_, mpi::SUM);
      mpi::reduceInplace(x_z_FluidCellCount_, mpi::SUM);
      mpi::reduceInplace(maxFluidZPos_, mpi::MAX);

      WALBERLA_ROOT_SECTION()
      {
         // normalize
         for (uint_t i = 0; i < x_z_VelocityX_.size(); ++i)
         {
            if (x_z_FluidCellCount_[i] > 0) x_z_VelocityX_[i] /= real_c(x_z_FluidCellCount_[i]);
         }
         real_t invNumYCells = 1_r / real_c(ylen_);
         std::for_each(x_z_SolidVolumeFraction_.begin(), x_z_SolidVolumeFraction_.end(),
                     [invNumYCells](real_t& n) { n *= invNumYCells; });
         std::for_each(x_z_FillLevel_.begin(), x_z_FillLevel_.end(), [invNumYCells](real_t& n) { n *= invNumYCells; });
      }
      // note: only root process has the complete information!
   }

   void writeToFile(){
      WALBERLA_ROOT_SECTION()
      {
         const std::string base = baseFolderName_ + "/avgDataSlices";
         std::filesystem::create_directories(base);
         write2DVectorToFile(x_z_SolidVolumeFraction_, xlen_, zlen_,
                              base + "/svfSlice_" + std::to_string(executionCounter_) + ".txt");
         write2DVectorToFile(x_z_FillLevel_, xlen_, zlen_,
                              base + "/fillSlice_" + std::to_string(executionCounter_) + ".txt");
         write2DVectorToFile(x_z_VelocityX_, xlen_, zlen_,
                              base + "/velXSlice_" + std::to_string(executionCounter_) + ".txt");

         dataCols_->add("maxFluidPos_z", int(maxFluidZPos_));
      }
   }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const ConstBlockDataID flagFieldID_;
   const ConstBlockDataID fillFieldID_;
   const ConstBlockDataID pdfFieldID_;

   uint_t xlen_;
   uint_t ylen_;
   uint_t zlen_;
   uint_t maxFluidZPos_;
   std::vector< real_t > x_z_SolidVolumeFraction_;
   std::vector< real_t > x_z_FillLevel_;
   std::vector< real_t > x_z_VelocityX_;
   std::vector< uint_t > x_z_FluidCellCount_;

   std::string baseFolderName_;
   shared_ptr< DataCols > dataCols_;
   uint_t frequency_;
   uint_t executionCounter_;
};


template< typename ParticleAccessor_T >
class BedloadTransportEvaluator
{
 public:
   BedloadTransportEvaluator(const shared_ptr< ParticleAccessor_T >& ac, real_t normalizationFactor,
                             uint_t numParticles, const shared_ptr< DataCols >& dataCols, uint_t frequency)
      : ac_(ac), normalizationFactor_(normalizationFactor), numParticles_(numParticles),
        dataCols_(dataCols), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {  
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   void compute(){
      real_t transportRate(real_t(0));
      real_t velocity(real_t(0));

      for (uint_t i = uint_t(0); i < ac_->size(); ++i)
      {
         if (!isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GHOST) &&
             !isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GLOBAL))
         {
            auto velX = ac_->getLinearVelocity(i)[0];
            transportRate += velX * ac_->getVolume(i);
            velocity += velX;
         }
      }

      // only reduce to root
      mpi::reduceInplace(transportRate, mpi::SUM);
      mpi::reduceInplace(velocity, mpi::SUM);

      avgTransportRate_ = transportRate * normalizationFactor_;
      averageVelocity_  = velocity / real_c(numParticles_);
      // note: only root process has the complete information!
   }

   void writeToFile(){
      WALBERLA_ROOT_SECTION()
      {
         dataCols_->add("avgParticleTransportRate", avgTransportRate_);
         dataCols_->add("avgParticleVelocity", averageVelocity_);
      }
   }

   // sum_i V_p,i * u_p,i / (L*W)
   real_t getTransportRate() const { return avgTransportRate_; }
   real_t getAverageVelocity() const { return averageVelocity_; }

 private:
   shared_ptr< ParticleAccessor_T > ac_;
   real_t normalizationFactor_;
   uint_t numParticles_;
   real_t averageVelocity_;
   real_t avgTransportRate_;

   shared_ptr< DataCols > dataCols_;
   uint_t frequency_;
   uint_t executionCounter_;
};


// Thesis item 3 (THESIS_PLAN.md): measures a developing clog's backwater effect directly, rather than
// only observing it visually -- tracks local water elevation and pressure at two explicit probe
// points (e.g. upstream/downstream of a pier or other obstruction; the caller decides where these are
// via prm, this evaluator has no notion of "obstruction" itself, so it works equally for a no-clog
// open-channel reference run with no obstruction at all).
//
// `columnPrefix` (2026-09-14, per Schalko, Follett & Nepf 2023, Water Resources Research -- a
// lateral-gap logjam study directly analogous to a two-pillar-with-a-gap geometry): that paper found
// real lateral variation in water depth (~7.5%+-1.8%) between a partial-width obstruction's own
// blocked-flow region and the open gap beside it, and took SEPARATE longitudinal profiles at each --
// the gap centerline alone is not representative of the pillar-face region where the actual backwater
// pile-up is expected to be strongest. Two instances of this evaluator, each with its own
// `columnPrefix` (e.g. "clog" for the existing gap-centerline probe pair, "clogPillar" for a second
// pair placed at the pillar's own Y), can therefore coexist and write to distinct dataCols columns
// without collision -- see BuildTimeloop.h's own wiring.
//
// Elevation: sum of the fill level over Z (skipping solid/moving-obstacle cells) for every (x,y) in a
// small XY window around the probe, averaged over the window -- the same "integrate fill over the
// column" technique AverageDataSliceEvaluator already uses above, which for a flat, non-solid column
// gives exactly the local liquid depth (fully-liquid cells contribute 1 each, the single interface
// cell contributes its fractional fill -- no separate free-surface search/reconstruction needed).
//
// Pressure: lattice pressure p=(rho-1)/3 (c_s^2=1/3 convention, matching DataInitialization.h's own
// HydrostaticDensityFtor) read at the probe's EXACT cell (the caller is expected to place probes at a
// sensible submerged position via prm) -- reported both raw and as "excess" over the plain
// hydrostatic reference rho_hyd(z)=1+3*a_z*(z-H) at that same position (same formula/inputs as
// HydrostaticDensityFtor, inlined here rather than depending on DataInitialization.h to avoid an
// otherwise-unnecessary include), since the raw number is dominated by local water depth while the
// excess isolates the obstruction/backwater-driven component.
template< typename AntidunesBoundaryHandling_T, typename PdfField_T, typename FlagField_T, typename ScalarField_T >
class ClogMeasurementEvaluator
{
 public:
   ClogMeasurementEvaluator(const shared_ptr< StructuredBlockStorage >& blocks, const ConstBlockDataID& flagFieldID,
                            const ConstBlockDataID& fillFieldID, const ConstBlockDataID& pdfFieldID,
                            const Vector3< real_t >& upstreamPos, const Vector3< real_t >& downstreamPos,
                            real_t probeHalfWidthXY, const shared_ptr< Vector3< real_t > >& acceleration,
                            real_t absoluteLiquidHeight, const shared_ptr< DataCols >& dataCols, uint_t frequency,
                            std::string columnPrefix = "clog", bool measureFlowRate = false,
                            std::vector< Vector3< real_t > > pillarCenters = {},
                            std::vector< real_t > pillarRadii = {}, real_t channelWidthY = real_t(0),
                            bool measureBlockageRatio = false)
      : blocks_(blocks), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID), pdfFieldID_(pdfFieldID),
        upstreamPos_(upstreamPos), downstreamPos_(downstreamPos), probeHalfWidthXY_(probeHalfWidthXY),
        acceleration_(acceleration), absoluteLiquidHeight_(absoluteLiquidHeight), dataCols_(dataCols),
        frequency_(frequency), executionCounter_(uint_c(0)), columnPrefix_(std::move(columnPrefix)),
        measureFlowRate_(measureFlowRate), pillarCenters_(std::move(pillarCenters)),
        pillarRadii_(std::move(pillarRadii)), channelWidthY_(channelWidthY),
        measureBlockageRatio_(measureBlockageRatio)
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      if (executionCounter_ % frequency_ == uint_c(0))
      {
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   real_t hydrostaticPressure(const Vector3< real_t >& pos) const
   {
      // p = (rho_hyd - 1) / 3, rho_hyd(z) = 1 + 3*a_z*(z-H)  =>  p_hyd(z) = a_z*(z-H)
      return (*acceleration_)[2] * (pos[2] - absoluteLiquidHeight_);
   }

   // Sum of the fill level over Z (skipping solid/moving-obstacle cells) for every (x,y) column whose
   // center falls within probeHalfWidthXY_ of probePos, averaged over the number of such columns --
   // gives the local liquid depth at the probe, per the class-level comment above.
   real_t elevationAt(const Vector3< real_t >& probePos) const
   {
      real_t sum = real_t(0);

      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);
         const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);

         const auto solidMO     = flagField->getFlag(AntidunesBoundaryHandling_T::movingObstacleFlagID);
         const auto solidNoSlip = flagField->getFlag(AntidunesBoundaryHandling_T::noSlipFlagID);

         CellInterval xyz = flagField->xyzSize();
         Cell globalCell;

         for (auto cell : xyz)
         {
            blocks_->transformBlockLocalToGlobalCell(globalCell, *block, cell);
            if (std::abs(real_c(globalCell.x()) - probePos[0]) > probeHalfWidthXY_) continue;
            if (std::abs(real_c(globalCell.y()) - probePos[1]) > probeHalfWidthXY_) continue;
            if (flagField->isFlagSet(cell, solidMO) || flagField->isFlagSet(cell, solidNoSlip)) continue;

            sum += fillField->get(cell);
         }
      }

      // number of (x,y) columns in the probe window -- computed analytically (not counted during the
      // scan above) so it is independent of block decomposition/rank ownership.
      const int halfW          = int_c(std::floor(probeHalfWidthXY_));
      const uint_t columnCount = uint_c((2 * halfW + 1) * (2 * halfW + 1));

      mpi::reduceInplace(sum, mpi::SUM);
      return columnCount > 0 ? sum / real_c(columnCount) : real_t(0);
   }

   // Volumetric flow rate Q = sum over the full y-z cross-section at global cell column x = floor(planeX) of
   // (fill level * u_x), and the wet cross-section area A = sum of the fill level, both over the FLUID cells
   // of that plane only, averaged over the (2*floor(probeHalfWidthXY_)+1) neighbouring y-z planes
   // x-halfW..x+halfW (a single 1-cell plane can fall exactly into the 1-LU slit between two rows of bodies and
   // then says nothing about the flow; Q is conserved from plane to plane in steady flow, so the average is
   // the same quantity with less sampling noise). The cell selection is (fill > 0, not a moving-obstacle cell, not a no-slip cell) -- the same cell selection
   // AverageDataSliceEvaluator uses for its slice velocity. Cells covered by floaters, pillars or walls are
   // excluded, so a body standing in the plane reduces A (its cross-section) and adds nothing to Q; the flux
   // through the plane is therefore the fluid flux, independent of what shape the body has, and Q/A is the mean
   // velocity of the water actually crossing the plane. Interface cells count with their fill level as weight,
   // gas cells (fill 0) not at all. Velocity is pdfField->getVelocity, as in MeanVelocityComputer.
   // Lattice units: Q in cells^3 per time step, A in cells^2. Signed: back-flow reduces Q.
   void flowRateAt(real_t planeX, real_t& flowRate, real_t& wetArea) const
   {
      flowRate = real_t(0);
      wetArea  = real_t(0);
      const cell_idx_t planeCell = cell_idx_c(std::floor(planeX));
      const cell_idx_t halfW     = cell_idx_c(std::floor(probeHalfWidthXY_));

      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const PdfField_T* const pdfField     = block->getData< const PdfField_T >(pdfFieldID_);
         const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);
         const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);

         const auto solidMO     = flagField->getFlag(AntidunesBoundaryHandling_T::movingObstacleFlagID);
         const auto solidNoSlip = flagField->getFlag(AntidunesBoundaryHandling_T::noSlipFlagID);

         // block-local x range of the planes planeCell-halfW .. planeCell+halfW that lie inside this block
         Cell localLo, localHi;
         blocks_->transformGlobalToBlockLocalCell(localLo, *block, Cell(planeCell - halfW, cell_idx_t(0), cell_idx_t(0)));
         blocks_->transformGlobalToBlockLocalCell(localHi, *block, Cell(planeCell + halfW, cell_idx_t(0), cell_idx_t(0)));
         CellInterval xyz = flagField->xyzSize();
         const cell_idx_t xLo = std::max(localLo.x(), xyz.xMin());
         const cell_idx_t xHi = std::min(localHi.x(), xyz.xMax());

         for (cell_idx_t x = xLo; x <= xHi; ++x)
            for (cell_idx_t z = xyz.zMin(); z <= xyz.zMax(); ++z)
               for (cell_idx_t y = xyz.yMin(); y <= xyz.yMax(); ++y)
               {
                  const Cell cell(x, y, z);
                  if (flagField->isFlagSet(cell, solidMO) || flagField->isFlagSet(cell, solidNoSlip)) continue;

                  const real_t fill = fillField->get(cell);
                  if (fill > real_t(0))
                  {
                     flowRate += fill * pdfField->getVelocity(cell)[0];
                     wetArea += fill;
                  }
               }
      }

      const real_t planeCount = real_c(2 * halfW + 1);
      flowRate /= planeCount;
      wetArea /= planeCount;

      mpi::reduceInplace(flowRate, mpi::SUM);
      mpi::reduceInplace(wetArea, mpi::SUM);
   }

   // Fraction of the REFERENCE (undisturbed) wetted cross-section occupied by the pillar(s) and the
   // jam of floaters, at the y-z plane through a probe's x (averaged over the same 2w+1 planes as
   // flowRateAt -- same review motivation: a single plane can fall exactly into a gap between two
   // floaters and read as unblocked). This is `sec:clog-secondary-metrics`'s "blockage ratio",
   // finally computed rather than only defined: \parencite{panici_formation_2018} relate afflux to
   // blockage ratio measured against the UNDISTURBED (reference) flow area, not the locally
   // backed-up one, which is why the denominator here is the fixed reference cross-section
   // `channelWidthY_ * absoluteLiquidHeight_`, not a per-step remeasured wetted area -- this makes
   // the ratio a scenario-describing number comparable across runs, matching the existing text's own
   // "natural axis for comparing scenarios" framing, rather than a probe-placement-sensitive local
   // reading that would move every time the water backs up.
   //
   // Pillar contribution is obtained ANALYTICALLY from the known, fixed pillar geometry rather than by
   // scanning flags: each `CylinderPillar` is a vertical prism of circular footprint (radius R, centre
   // (cx,cy)) spanning the full domain height (see SetupStages.h's `createCylinder` call), so at a
   // given x its intersection with the channel's valid y-range is a vertical strip whose width is the
   // circle's y-chord at that x, clipped to [0, channelWidthY_]; its contribution to the wetted area up
   // to the reference depth is simply that width times the reference depth (constant along z, since
   // the prism's cross-section does not vary with height). Jam contribution counts moving-obstacle
   // (floater) cells at or below the same reference depth, EXCLUDING any cell that geometrically falls
   // inside a pillar's own footprint -- pillars are themselves ordinary (fixed, infinite-mass) mesa_pd
   // bodies and so may also carry the moving-obstacle flag; this exclusion prevents double-counting a
   // floater resting directly against a pillar as extra pillar area, and is correct regardless of
   // whether pillar cells actually carry that flag or not (a purely geometric, not flag-based, test).
   // Units: cells^2 for the area outputs, dimensionless for the ratio.
   void blockageAt(real_t planeX, real_t& blockageRatio, real_t& pillarArea, real_t& jamArea) const
   {
      const cell_idx_t planeCell  = cell_idx_c(std::floor(planeX));
      const cell_idx_t halfW      = cell_idx_c(std::floor(probeHalfWidthXY_));
      const real_t referenceDepth = absoluteLiquidHeight_;
      const cell_idx_t zRefCell   = cell_idx_c(std::floor(referenceDepth));

      pillarArea = real_t(0);
      for (size_t p = 0; p < pillarCenters_.size(); ++p)
      {
         const real_t dx = planeX - pillarCenters_[p][0];
         const real_t R  = pillarRadii_[p];
         if (std::abs(dx) >= R) continue;
         const real_t chord = std::sqrt(R * R - dx * dx);
         const real_t cy    = pillarCenters_[p][1];
         const real_t yLo   = std::max(real_t(0), cy - chord);
         const real_t yHi   = std::min(channelWidthY_, cy + chord);
         const real_t width = std::max(real_t(0), yHi - yLo);
         pillarArea += width * referenceDepth;
      }

      jamArea = real_t(0);
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const FlagField_T* const flagField = block->getData< const FlagField_T >(flagFieldID_);
         const auto solidMO = flagField->getFlag(AntidunesBoundaryHandling_T::movingObstacleFlagID);

         Cell localLo, localHi;
         blocks_->transformGlobalToBlockLocalCell(localLo, *block, Cell(planeCell - halfW, cell_idx_t(0), cell_idx_t(0)));
         blocks_->transformGlobalToBlockLocalCell(localHi, *block, Cell(planeCell + halfW, cell_idx_t(0), cell_idx_t(0)));
         CellInterval xyz = flagField->xyzSize();
         const cell_idx_t xLo = std::max(localLo.x(), xyz.xMin());
         const cell_idx_t xHi = std::min(localHi.x(), xyz.xMax());

         Cell globalCell;
         for (cell_idx_t x = xLo; x <= xHi; ++x)
            for (cell_idx_t z = xyz.zMin(); z <= xyz.zMax(); ++z)
               for (cell_idx_t y = xyz.yMin(); y <= xyz.yMax(); ++y)
               {
                  const Cell cell(x, y, z);
                  if (!flagField->isFlagSet(cell, solidMO)) continue;

                  blocks_->transformBlockLocalToGlobalCell(globalCell, *block, cell);
                  if (globalCell.z() > zRefCell) continue;

                  bool insidePillar = false;
                  for (size_t p = 0; p < pillarCenters_.size(); ++p)
                  {
                     const real_t ddx = real_c(globalCell.x()) - pillarCenters_[p][0];
                     const real_t ddy = real_c(globalCell.y()) - pillarCenters_[p][1];
                     if (ddx * ddx + ddy * ddy < pillarRadii_[p] * pillarRadii_[p]) { insidePillar = true; break; }
                  }
                  if (insidePillar) continue;

                  jamArea += real_t(1);
               }
      }

      const real_t planeCount = real_c(2 * halfW + 1);
      jamArea /= planeCount;
      mpi::reduceInplace(jamArea, mpi::SUM);

      const real_t referenceArea = channelWidthY_ * referenceDepth;
      blockageRatio = referenceArea > real_t(0) ? (pillarArea + jamArea) / referenceArea : real_t(0);
   }

   void compute()
   {
      upstreamElevation_   = elevationAt(upstreamPos_);
      downstreamElevation_ = elevationAt(downstreamPos_);

      if (measureFlowRate_)
      {
         flowRateAt(upstreamPos_[0], upstreamFlowRate_, upstreamWetArea_);
         flowRateAt(downstreamPos_[0], downstreamFlowRate_, downstreamWetArea_);
      }

      if (measureBlockageRatio_)
      {
         blockageAt(upstreamPos_[0], upstreamBlockageRatio_, upstreamPillarArea_, upstreamJamArea_);
         blockageAt(downstreamPos_[0], downstreamBlockageRatio_, downstreamPillarArea_, downstreamJamArea_);
      }

      upstreamPressure_   = real_t(0);
      downstreamPressure_ = real_t(0);

      const Cell upstreamGlobalCell(cell_idx_c(upstreamPos_[0]), cell_idx_c(upstreamPos_[1]), cell_idx_c(upstreamPos_[2]));
      const Cell downstreamGlobalCell(cell_idx_c(downstreamPos_[0]), cell_idx_c(downstreamPos_[1]), cell_idx_c(downstreamPos_[2]));

      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const FlagField_T* const flagField = block->getData< const FlagField_T >(flagFieldID_);
         const PdfField_T* const pdfField   = block->getData< const PdfField_T >(pdfFieldID_);
         const CellInterval xyz              = flagField->xyzSize();

         Cell localCell;
         blocks_->transformGlobalToBlockLocalCell(localCell, *block, upstreamGlobalCell);
         if (xyz.contains(localCell)) { upstreamPressure_ = (pdfField->getDensity(localCell) - real_t(1)) / real_t(3); }

         blocks_->transformGlobalToBlockLocalCell(localCell, *block, downstreamGlobalCell);
         if (xyz.contains(localCell)) { downstreamPressure_ = (pdfField->getDensity(localCell) - real_t(1)) / real_t(3); }
      }

      mpi::reduceInplace(upstreamPressure_, mpi::SUM);
      mpi::reduceInplace(downstreamPressure_, mpi::SUM);
      // note: only root process has the complete information!
   }

   void writeToFile()
   {
      WALBERLA_ROOT_SECTION()
      {
         dataCols_->add(columnPrefix_ + "UpstreamElevation", upstreamElevation_);
         dataCols_->add(columnPrefix_ + "DownstreamElevation", downstreamElevation_);
         dataCols_->add(columnPrefix_ + "UpstreamPressure", upstreamPressure_);
         dataCols_->add(columnPrefix_ + "DownstreamPressure", downstreamPressure_);
         dataCols_->add(columnPrefix_ + "UpstreamExcessPressure", upstreamPressure_ - hydrostaticPressure(upstreamPos_));
         dataCols_->add(columnPrefix_ + "DownstreamExcessPressure", downstreamPressure_ - hydrostaticPressure(downstreamPos_));
         if (measureFlowRate_)
         {
            dataCols_->add(columnPrefix_ + "UpstreamFlowRate", upstreamFlowRate_);
            dataCols_->add(columnPrefix_ + "DownstreamFlowRate", downstreamFlowRate_);
            dataCols_->add(columnPrefix_ + "UpstreamWetArea", upstreamWetArea_);
            dataCols_->add(columnPrefix_ + "DownstreamWetArea", downstreamWetArea_);
         }
         if (measureBlockageRatio_)
         {
            dataCols_->add(columnPrefix_ + "UpstreamBlockageRatio", upstreamBlockageRatio_);
            dataCols_->add(columnPrefix_ + "DownstreamBlockageRatio", downstreamBlockageRatio_);
            dataCols_->add(columnPrefix_ + "UpstreamPillarArea", upstreamPillarArea_);
            dataCols_->add(columnPrefix_ + "DownstreamPillarArea", downstreamPillarArea_);
            dataCols_->add(columnPrefix_ + "UpstreamJamArea", upstreamJamArea_);
            dataCols_->add(columnPrefix_ + "DownstreamJamArea", downstreamJamArea_);
         }
      }
   }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const ConstBlockDataID flagFieldID_;
   const ConstBlockDataID fillFieldID_;
   const ConstBlockDataID pdfFieldID_;

   Vector3< real_t > upstreamPos_;
   Vector3< real_t > downstreamPos_;
   real_t            probeHalfWidthXY_;

   shared_ptr< Vector3< real_t > > acceleration_;
   real_t absoluteLiquidHeight_;

   real_t upstreamElevation_;
   real_t downstreamElevation_;
   real_t upstreamPressure_;
   real_t downstreamPressure_;

   // flow-rate measurement (optional, see flowRateAt): Q and wet area of the y-z plane at each probe's x
   real_t upstreamFlowRate_   = real_t(0);
   real_t downstreamFlowRate_ = real_t(0);
   real_t upstreamWetArea_    = real_t(0);
   real_t downstreamWetArea_  = real_t(0);

   shared_ptr< DataCols > dataCols_;
   uint_t frequency_;
   uint_t executionCounter_;
   std::string columnPrefix_;
   bool        measureFlowRate_;

   // blockage-ratio measurement (optional, see blockageAt): fixed pillar geometry + channel width,
   // and the resulting ratio/area outputs at each probe's x
   std::vector< Vector3< real_t > > pillarCenters_;
   std::vector< real_t >            pillarRadii_;
   real_t                           channelWidthY_;
   bool                             measureBlockageRatio_ = false;

   real_t upstreamBlockageRatio_   = real_t(0);
   real_t downstreamBlockageRatio_ = real_t(0);
   real_t upstreamPillarArea_      = real_t(0);
   real_t downstreamPillarArea_    = real_t(0);
   real_t upstreamJamArea_         = real_t(0);
   real_t downstreamJamArea_       = real_t(0);
};


template< typename ParticleAccessor_T >
class TotalHydrodynamicForceEvaluator
{
 public:
   TotalHydrodynamicForceEvaluator(const shared_ptr< ParticleAccessor_T >& ac,
       const shared_ptr< DataCols >& dataCols, uint_t frequency)
      : ac_(ac), dataCols_(dataCols),
      frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {  
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   void compute(){
      totalHydrodynamicForce_ = Vector3< real_t >(real_c(0));
      maxParticlePosZ_ = real_c(0);

      for (uint_t i = uint_t(0); i < ac_->size(); ++i)
      {
         if (!isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GHOST) &&
            !isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GLOBAL))
         {
            // NOT getHydrodynamicForce(i): this evaluator runs in timeloopEval, which executes
            // AFTER ParticleKernels::resetHydroForces() has already zeroed HydrodynamicForce/Torque
            // for the step (Antidunes.cpp's per-timestep order: reduceAndAvgHydroForces ->
            // doDEMSubCycles -> resetHydroForces -> timeloopAfterDEM -> timeloopEval). Reading it
            // here always summed to exactly zero, for every particle, every step, regardless of
            // buoyancy model (ResetHydrodynamicForceTorqueKernel.h resets unconditionally).
            // OldHydrodynamicForce is set to this same step's real (optionally time-averaged, see
            // useHydrodynamicForceTimeAveraging) value inside reduceAndAvgHydroForces() and is never
            // touched by the reset kernel, so it still holds the correct value here.
            totalHydrodynamicForce_ += ac_->getOldHydrodynamicForce(i);
            maxParticlePosZ_ = std::max(ac_->getPosition(i)[2], maxParticlePosZ_);
         }
      }

      // only reduce to root
      mpi::reduceInplace(totalHydrodynamicForce_, mpi::SUM);
      mpi::reduceInplace(maxParticlePosZ_, mpi::MAX);
   }

   void writeToFile(){
      WALBERLA_ROOT_SECTION()
      {
         dataCols_->add("totalHydroForce_x", totalHydrodynamicForce_[0]);
         dataCols_->add("totalHydroForce_y", totalHydrodynamicForce_[1]);
         dataCols_->add("totalHydroForce_z", totalHydrodynamicForce_[2]);
         dataCols_->add("maxParticlePos_z", maxParticlePosZ_);
      }
   }

   Vector3< real_t > getTotalHydrodynamicForce() const { return totalHydrodynamicForce_; }

 private:
   shared_ptr< ParticleAccessor_T > ac_;
   Vector3< real_t > totalHydrodynamicForce_;
   real_t maxParticlePosZ_;

   shared_ptr< DataCols > dataCols_;
   uint_t frequency_;
   uint_t executionCounter_;
};


} // namespace antidunes
} // namespace walberla
