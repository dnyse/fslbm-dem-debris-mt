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
//! \file StreamReconstructAdvectSweep.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Sweep for reconstruction of PDFs, streaming of PDFs (only in interface cells), advection of mass, update of
//!        bubble volumes and marking interface cells for conversion.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/math/Vector3.h"
#include "core/mpi/Reduce.h"
#include "core/timing/TimingPool.h"

#include "field/FieldClone.h"
#include "field/FlagField.h"

#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"
#include "lbm/sweeps/StreamPull.h"

#include "stencil/D3Q27.h"

#include "PdfReconstructionModel.h"
#include "functionality/AdvectMass.h"
#include "functionality/FindInterfaceCellConversion.h"
#include "functionality/GetOredNeighborhood.h"
#include "functionality/ReconstructInterfaceCellABB.h"

namespace walberla
{
namespace free_surface
{
template< typename LatticeModel_T, typename BoundaryHandling_T, typename FlagField_T, typename FlagInfo_T,
          typename ScalarField_T, typename VectorField_T, bool useCodegen >
class StreamReconstructAdvectSweep
{
 public:
   using flag_t     = typename FlagInfo_T::flag_t;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   StreamReconstructAdvectSweep(real_t sigma, BlockDataID handlingID, BlockDataID fillFieldID, BlockDataID flagFieldID,
                                BlockDataID pdfField, ConstBlockDataID normalFieldID, ConstBlockDataID curvatureFieldID,
                                const FlagInfo_T& flagInfo, BubbleModelBase* bubbleModel,
                                const PdfReconstructionModel& pdfReconstructionModel, bool useSimpleMassExchange,
                                real_t cellConversionThreshold, real_t cellConversionForceThreshold,
                                bool seedPdfCloneOnFirstUse = false)
      : sigma_(sigma), handlingID_(handlingID), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        pdfFieldID_(pdfField), normalFieldID_(normalFieldID), curvatureFieldID_(curvatureFieldID), flagInfo_(flagInfo),
        bubbleModel_(bubbleModel), neighborhoodFlagFieldClone_(flagFieldID), fillFieldClone_(fillFieldID),
        pdfFieldClone_(pdfField), pdfReconstructionModel_(pdfReconstructionModel),
        useSimpleMassExchange_(useSimpleMassExchange), cellConversionThreshold_(cellConversionThreshold),
        cellConversionForceThreshold_(cellConversionForceThreshold), seedPdfCloneOnFirstUse_(seedPdfCloneOnFirstUse)
   {}

   void operator()(IBlock* const block);

 protected:
   real_t sigma_; // surface tension

   BlockDataID handlingID_;
   BlockDataID fillFieldID_;
   BlockDataID flagFieldID_;
   BlockDataID pdfFieldID_;

   ConstBlockDataID normalFieldID_;
   ConstBlockDataID curvatureFieldID_;

   FlagInfo_T flagInfo_;
   bubble_model::BubbleModelBase* const bubbleModel_;

   // efficient clones of fields to provide temporary fields (for writing to)
   field::FieldClone< FlagField_T, true > neighborhoodFlagFieldClone_;
   field::FieldClone< ScalarField_T, true > fillFieldClone_;
   field::FieldClone< PdfField_T, true > pdfFieldClone_;

   PdfReconstructionModel pdfReconstructionModel_;
   bool useSimpleMassExchange_;
   real_t cellConversionThreshold_;
   real_t cellConversionForceThreshold_;
   // [Claude, 2026-09-28] false (DEFAULT, upstream behaviour): the scratch PDF clone keeps its uninitialised first-creation memory.
   // true: copy the real PDF field into the clone once on its first creation (see the 2026-09-25 comment in operator()).
   bool seedPdfCloneOnFirstUse_;
}; // class StreamReconstructAdvectSweep

template< typename LatticeModel_T, typename BoundaryHandling_T, typename FlagField_T, typename FlagInfo_T,
          typename ScalarField_T, typename VectorField_T, bool useCodegen >
void StreamReconstructAdvectSweep< LatticeModel_T, BoundaryHandling_T, FlagField_T, FlagInfo_T, ScalarField_T,
                                   VectorField_T, useCodegen >::operator()(IBlock* const block)
{
   // fetch data
   ScalarField_T* const fillSrcField = block->getData< ScalarField_T >(fillFieldID_);
   PdfField_T* const pdfSrcField     = block->getData< PdfField_T >(pdfFieldID_);

   const ScalarField_T* const curvatureField = block->getData< const ScalarField_T >(curvatureFieldID_);
   const VectorField_T* const normalField    = block->getData< const VectorField_T >(normalFieldID_);
   FlagField_T* const flagField              = block->getData< FlagField_T >(flagFieldID_);

   // temporary fields that act as destination fields
   bool pdfCloneJustCreated = false;
   PdfField_T* const pdfDstField            = pdfFieldClone_.get(block, pdfCloneJustCreated);
   FlagField_T* const neighborhoodFlagField = neighborhoodFlagFieldClone_.get(block);
   ScalarField_T* const fillDstField        = fillFieldClone_.get(block);

   // [Claude, 2026-09-25] FIX (THESIS_PLAN.md, Antidunes hopper-clogging investigation): pdfFieldClone_
   // allocates its buffer via cloneUninitialized() (FieldClone.h) on its ONE-TIME first creation -- this
   // does NOT copy the real field's actual values, genuinely uninitialized memory. The obstacle branch
   // below never writes into pdfDstFieldIt at all (only the fill level is carried forward), so any PDF
   // component of an obstacle-flagged cell that is ALSO never independently refreshed by a boundary
   // condition's per-link bounce-back (i.e. one facing another obstacle cell, not real fluid) can remain
   // on this never-written memory indefinitely -- which typical OS zero-fill-on-first-touch behaviour
   // makes read back as exactly 0.0, and a subsequent 1/rho collision-kernel division then turns that
   // into a full NaN once this buffer becomes "current" via the src/dst swap below. This only matters on
   // the FIRST-ever creation of this clone (every later call reuses the same persistent object), so fix
   // it at the source: seed the fresh clone with the real field's own current, properly-initialised
   // contents once, so it starts as a genuine mirror rather than blank memory.
   if (pdfCloneJustCreated && seedPdfCloneOnFirstUse_)
   {
      WALBERLA_FOR_ALL_CELLS_XYZ(pdfSrcField,
         for (uint_t f = uint_t(0); f < PdfField_T::F_SIZE; ++f)
            pdfDstField->get(x, y, z, f) = pdfSrcField->get(x, y, z, f);
      )
   }

   // combine all neighbor flags using bitwise OR and write them to the neighborhood field
   // this is simply a pre-computation of often required values
   // IMPORTANT REMARK: the "OredNeighborhood" is also required for each cell in the first ghost layer; this requires
   // access to all first ghost layer cell's neighbors (i.e., to the second ghost layer)
   WALBERLA_CHECK_GREATER_EQUAL(flagField->nrOfGhostLayers(), uint_c(2));
   getOredNeighborhood< typename LatticeModel_T::Stencil >(flagField, neighborhoodFlagField);

   // explicitly avoid OpenMP due to bubble model update (reportFillLevelChange)
   WALBERLA_FOR_ALL_CELLS_OMP(
      pdfDstFieldIt, pdfDstField, pdfSrcFieldIt, pdfSrcField, fillSrcFieldIt, fillSrcField, fillDstFieldIt,
      fillDstField, flagFieldIt, flagField, neighborhoodFlagFieldIt, neighborhoodFlagField, normalFieldIt, normalField,
      curvatureFieldIt, curvatureField, omp critical, {
         if (flagInfo_.isInterface(flagFieldIt))
         {
            // get density (rhoGas) for interface PDF reconstruction
            const real_t bubbleDensity = bubbleModel_->getDensity(block, pdfSrcFieldIt.cell());
            const real_t rhoGas        = computeDeltaRhoLaplacePressure(sigma_, *curvatureFieldIt) + bubbleDensity;

            // reconstruct missing PDFs coming from gas neighbors according to the specified model (see dissertation of
            // N. Thuerey, 2007, section 4.2); reconstruction includes streaming of PDFs to interface cells (no LBM
            // stream required here)
            (reconstructInterfaceCellLegacy< LatticeModel_T, FlagField_T >) (flagField, pdfSrcFieldIt, flagFieldIt,
                                                                             normalFieldIt, flagInfo_, rhoGas,
                                                                             pdfDstFieldIt, pdfReconstructionModel_);

            // density before LBM stream (post-collision)
            const real_t oldRho = lbm::getDensity(pdfSrcField->latticeModel(), pdfSrcFieldIt);

            // density after LBM stream in reconstruction
            const real_t newRho = lbm::getDensity(pdfDstField->latticeModel(), pdfDstFieldIt);

            // compute mass advection using post-collision PDFs (explicitly not PDFs updated by stream above)
            const real_t deltaMass =
               (advectMass< LatticeModel_T, FlagField_T, typename ScalarField_T::iterator,
                            typename PdfField_T::iterator, typename FlagField_T::iterator,
                            typename FlagField_T::iterator, FlagInfo_T >) (flagField, fillSrcFieldIt, pdfSrcFieldIt,
                                                                           flagFieldIt, neighborhoodFlagFieldIt,
                                                                           flagInfo_, useSimpleMassExchange_);

            // update fill level after LBM stream and mass exchange
            *fillDstFieldIt        = (*fillSrcFieldIt * oldRho + deltaMass) / newRho;
            const real_t deltaFill = *fillDstFieldIt - *fillSrcFieldIt;

            // [Claude, 2026-09-25] DIAGNOSTIC: prove/disprove the hypothesis that a near-total obstacle
            // enclosure degrades reconstructInterfaceCellLegacy's curvature/normal-based reconstruction
            // enough to make newRho itself degenerate (zero or non-finite), which this division then
            // turns directly into a non-finite fill level with no guard at all (THESIS_PLAN.md,
            // hopper clogging investigation). Fires only when the fill actually IS non-finite, so a
            // healthy run stays silent; throttled the same way as PdfReconstructionManager.h's own
            // [FALLBACK-ISOLATED-INTERFACE]/[PASS3-UNRESOLVED-ISLAND] warnings (first 20, then every
            // 5000th) since this is a generic library file with no Antidunes-specific step-window config
            // available here.
            if (!std::isfinite(*fillDstFieldIt))
            {
               static uint_t degenerateFillLogCount = uint_t(0);
               ++degenerateFillLogCount;
               if (degenerateFillLogCount <= uint_t(20) || (degenerateFillLogCount % uint_t(5000) == uint_t(0)))
               {
                  const Cell cell = pdfSrcFieldIt.cell();
                  uint_t obstacleNbrCount = uint_t(0);
                  for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
                  {
                     const flag_t nbrFlag = flagField->get(cell.x() + dir.cx(), cell.y() + dir.cy(), cell.z() + dir.cz());
                     if (flagInfo_.isObstacle(nbrFlag)) ++obstacleNbrCount;
                  }
                  const auto& aabb = block->getAABB();
                  const real_t cdx = aabb.xSize() / real_c(flagField->xSize());
                  WALBERLA_LOG_WARNING("[DEGENERATE-FILL-RECONSTRUCTION] (#" << degenerateFillLogCount << ") cell ("
                     << cell.x() << "," << cell.y() << "," << cell.z() << ") world~("
                     << (aabb.xMin() + (real_c(cell.x()) + real_t(0.5)) * cdx) << ","
                     << (aabb.yMin() + (real_c(cell.y()) + real_t(0.5)) * cdx) << ","
                     << (aabb.zMin() + (real_c(cell.z()) + real_t(0.5)) * cdx) << ") fillNew=" << *fillDstFieldIt
                     << " fillOld=" << *fillSrcFieldIt << " oldRho=" << oldRho << " newRho=" << newRho
                     << " rhoGas=" << rhoGas << " bubbleDensity=" << bubbleDensity << " deltaMass=" << deltaMass
                     << " curvature=" << *curvatureFieldIt << " normal=" << *normalFieldIt
                     << " obstacleNeighbours=" << obstacleNbrCount << "/26 -- degenerate newRho from"
                        " reconstructInterfaceCellLegacy suspected if obstacleNeighbours is large (near-total"
                        " enclosure starves the curvature/normal computation of real fluid geometry).");
               }
            }

            // update the volume of bubbles
            bubbleModel_->reportFillLevelChange(block, fillSrcFieldIt.cell(), deltaFill);
         }
         else // treat non-interface cells
         {
            // manually adjust the fill level to avoid outdated fill levels being copied from fillSrcField
            if (flagInfo_.isGas(flagFieldIt)) { *fillDstFieldIt = real_c(0); }
            else
            {
               if (flagInfo_.isLiquid(flagFieldIt))
               {
                  const Cell cell = pdfSrcFieldIt.cell();
                  if constexpr (useCodegen)
                  {
                     auto lbmSweepGenerated = typename LatticeModel_T::Sweep(pdfFieldID_);
                     const CellInterval ci(cell, cell);
                     lbmSweepGenerated.streamInCellInterval(pdfSrcField, pdfDstField, ci);
                  }
                  else
                  {
                     lbm::StreamPull< LatticeModel_T >::execute(pdfSrcField, pdfDstField, cell[0], cell[1], cell[2]);
                  }

                  *fillDstFieldIt = real_c(1);
               }
               else // flag is e.g. obstacle or outflow
               {
                  *fillDstFieldIt = *fillSrcFieldIt;
               }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS_XYZ_OMP

   pdfSrcField->swapDataPointers(pdfDstField);
   fillSrcField->swapDataPointers(fillDstField);

   BoundaryHandling_T* const handling = block->getData< BoundaryHandling_T >(handlingID_);

   // mark interface cell for conversion
   findInterfaceCellConversions< LatticeModel_T >(handling, fillSrcField, flagField, neighborhoodFlagField, flagInfo_,
                                                  cellConversionThreshold_, cellConversionForceThreshold_);
}

} // namespace free_surface
} // namespace walberla
