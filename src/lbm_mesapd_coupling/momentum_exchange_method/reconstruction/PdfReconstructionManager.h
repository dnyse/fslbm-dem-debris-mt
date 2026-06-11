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
//! \file PdfReconstructionManager.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "core/mpi/Reduce.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/Field.h"
#include "field/FlagField.h"
#include "field/iterators/IteratorMacros.h"

#include "lbm/field/PdfField.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/data/IAccessor.h"

#include <functional>

namespace walberla
{
namespace lbm_mesapd_coupling
{

//**************************************************************************************************************************************
/*!\brief Class to manage the reconstruction of PDFs that is needed when cells are becoming uncovered by moving
 * obstacles.
 *
 * Due to the explicit mapping of particles into the fluid domain via flags, the PDFs of cells that turned from obstacle
 * to fluid are missing and must be reconstructed in order to continue with the simulation. This class is to be used as
 * a sweep in a LBM simulation with moving obstacles and calls for each cell that is tagged as 'formerObstacle' the
 * specified reconstructor (see Reconstructor.h for the available variants). After the successful reconstruction of all
 * PDFs, the flags are updated to 'fluid'. For small obstacle fractions, an optimized variant is available that only
 * looks for 'formerObstacle' cells in the vicinity of available particles. It is activated via the
 * 'optimizeForSmallObstacleFraction' argument in the constructor.
 *
 *
 * Note:
 * A potential problem arises when the cell that is to be reconstructed, i.e. filled with fluid, is surrounded (given by
 * the LBM stencil) by only boundary cells. Consequently, the successive calls to LBM will simulate a local system
 * containing a single fluid cell surrounded by boundary cells with different velocity boundary conditions. This system
 * is not divergence free (as required by the mass conservation equation) and will thus lead to mass loss/gain depending
 * on the boundary conditions. This is a very rare case but if it happens, it will almost always lead to malicious
 * oscillations in the density. Note that this can also happen if several cells are locally "trapped" by surrounding
 * boundary cells. Such a case can not be detected easily and also not be curred easily. It has to be inspected whether
 * the behavior arising from this issue leads to problems in real simulations.
 */
//**************************************************************************************************************************************

// ScalarField_T (defaulted to void) is the fill-level field type, used only by the Bogner
// state-selection path to seed the volume fraction φ of freshly converted cells
// (Bogner & Rüde 2013, §5.2). For the non-Bogner path it stays void and all fill code is
// compiled out via `if constexpr`.
template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T,
          typename ScalarField_T = void >
class PdfReconstructionManager
{
   static_assert(std::is_base_of_v< mesa_pd::data::IAccessor, ParticleAccessor_T >,
                 "Provide a valid accessor as template");

 public:
   using FlagField_T = typename BoundaryHandling_T::FlagField;
   using flag_t      = typename BoundaryHandling_T::flag_t;

   // Original constructor — all FormerMO cells are uniformly set to `fluid`.
   inline PdfReconstructionManager(const shared_ptr< StructuredBlockStorage >& blockStorage,
                                   const BlockDataID& pdfFieldID, const BlockDataID& boundaryHandlingID,
                                   const BlockDataID& particleFieldID, const shared_ptr< ParticleAccessor_T >& ac,
                                   const FlagUID& formerObstacle, const FlagUID& fluid,
                                   const shared_ptr< Reconstructor_T >& reconstructor,
                                   bool conserveMomentumWhenRestoring,
                                   const bool optimizeForSmallObstacleFraction = false)
      : blockStorage_(blockStorage), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID),
        particleFieldID_(particleFieldID), ac_(ac), reconstructor_(reconstructor), formerObstacle_(formerObstacle),
        fluid_(fluid), conserveMomentumWhenRestoring_(conserveMomentumWhenRestoring),
        optimizeForSmallObstacleFraction_(optimizeForSmallObstacleFraction), useBognerStateSelection_(false),
        liquidFlagUID_(fluid), gasFlagUID_(fluid), interfaceFlagUID_(fluid), fillFieldID_(BlockDataID())
   {}

   // Bogner state-selection constructor [Bogner, 2009].
   // Each FormerMO cell x is converted to the FSLBM state determined by its
   // non-obstacle neighbourhood B(x) = {x + e_q : x + e_q is not an obstacle}:
   //   • B(x) has no liquid cells  → gas       (φ = 0)
   //   • B(x) has no gas cells     → liquid    (φ = 1, PDFs from EquilibriumReconstructor)
   //   • B(x) has liquid AND gas   → interface (φ = mean fill of interface cells in B(x))
   // The fill-level seeding of converted cells (Bogner & Rüde 2013, §5.2) is performed here in
   // pass 2 via `fillFieldID`; this replaces the role that ObstacleFillLevelSweep would play.
   inline PdfReconstructionManager(const shared_ptr< StructuredBlockStorage >& blockStorage,
                                   const BlockDataID& pdfFieldID, const BlockDataID& boundaryHandlingID,
                                   const BlockDataID& particleFieldID, const shared_ptr< ParticleAccessor_T >& ac,
                                   const FlagUID& formerObstacle, const FlagUID& fluid,
                                   const shared_ptr< Reconstructor_T >& reconstructor,
                                   bool conserveMomentumWhenRestoring, const FlagUID& liquidFlagUID,
                                   const FlagUID& gasFlagUID, const FlagUID& interfaceFlagUID,
                                   const BlockDataID& fillFieldID, const bool useBognerStateSelection = true,
                                   const bool optimizeForSmallObstacleFraction = false)
      : blockStorage_(blockStorage), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID),
        particleFieldID_(particleFieldID), ac_(ac), reconstructor_(reconstructor), formerObstacle_(formerObstacle),
        fluid_(fluid), conserveMomentumWhenRestoring_(conserveMomentumWhenRestoring),
        optimizeForSmallObstacleFraction_(optimizeForSmallObstacleFraction),
        useBognerStateSelection_(useBognerStateSelection),
        liquidFlagUID_(liquidFlagUID), gasFlagUID_(gasFlagUID), interfaceFlagUID_(interfaceFlagUID),
        fillFieldID_(fillFieldID)
   {}

   void inline operator()(IBlock* const block);

   // Diagnostic: MPI-reduce (to root), log, and reset the FormerMO→state conversion counters
   // accumulated since the last call. One line per call → invoke periodically (e.g. once per
   // evaluationFrequency steps via an AfterFunction). Counts INNER cells only, so the totals are
   // the true number of physical cells converted, not duplicated across block ghost layers.
   void logAndResetConversionCounts(const uint_t t)
   {
      std::vector< uint_t > c = { convGas_, convLiquid_, convIntfUnknown_, convIntfSurface_ };
      walberla::mpi::reduceInplace(c, walberla::mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-CONVERT t=" << t << "] cells converted since last report: "
         << "gas(A)=" << c[0] << " liquid(B)=" << c[1]
         << " interface-unknown(C)=" << c[2] << " interface-surface(D)=" << c[3]
         << " | total=" << (c[0] + c[1] + c[2] + c[3]));
      convGas_ = convLiquid_ = convIntfUnknown_ = convIntfSurface_ = uint_t(0);
   }

 private:
   // Seed the volume fraction φ of a just-converted cell (Bogner & Rüde 2013, §5.2):
   //   gas → 0, liquid → 1, interface → mean fill of the interface cells in the D3Q27
   //   neighbourhood B(x) (fallback 0.5 if none). Only instantiated for non-void ScalarField_T.
   template< typename ScalarFieldPtr_T >
   void seedConvertedCellFill(ScalarFieldPtr_T fillField, FlagField_T* flagField, cell_idx_t x, cell_idx_t y,
                              cell_idx_t z, flag_t targetFlag, flag_t gasFlag, flag_t liquidFlag,
                              flag_t interfaceFlag) const
   {
      real_t seededFill;
      if (isFlagSet(targetFlag, gasFlag)) { seededFill = real_t(0); }
      else if (isFlagSet(targetFlag, liquidFlag)) { seededFill = real_t(1); }
      else // interface
      {
         real_t sum = real_t(0);
         uint_t cnt = uint_t(0);
         for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
               for (int dx = -1; dx <= 1; ++dx)
               {
                  if (dx == 0 && dy == 0 && dz == 0) continue;
                  if (isFlagSet(flagField->get(x + dx, y + dy, z + dz), interfaceFlag))
                  {
                     sum += fillField->get(x + dx, y + dy, z + dz);
                     ++cnt;
                  }
               }
         seededFill = (cnt > uint_t(0)) ? sum / real_c(cnt) : real_t(0.5);
      }
      fillField->get(x, y, z) = seededFill;
   }

   shared_ptr< StructuredBlockStorage > blockStorage_;

   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID particleFieldID_;

   shared_ptr< ParticleAccessor_T > ac_;

   shared_ptr< Reconstructor_T > reconstructor_;

   const FlagUID formerObstacle_;
   const FlagUID fluid_;

   bool conserveMomentumWhenRestoring_;
   const bool optimizeForSmallObstacleFraction_;

   // Bogner state selection
   const bool useBognerStateSelection_;
   const FlagUID liquidFlagUID_;
   const FlagUID gasFlagUID_;
   const FlagUID interfaceFlagUID_;
   const BlockDataID fillFieldID_; // fill-level field for φ seeding (Bogner path only)

   // FormerMO → state conversion counters, one per branch of bognerTargetFlag (see the table in
   // contact_ring_crash_analysis.md §5.3). Accumulated over inner cells across blocks; flushed by
   // logAndResetConversionCounts. Mutable-free: operator() is non-const, the counting lambda
   // captures `this`.
   uint_t convGas_         = uint_t(0); // branch A: !hasLiquid              → gas
   uint_t convLiquid_      = uint_t(0); // branch B: !hasGas && !hasUnknown  → liquid
   uint_t convIntfUnknown_ = uint_t(0); // branch C: !hasGas &&  hasUnknown  → interface (fallback)
   uint_t convIntfSurface_ = uint_t(0); // branch D:  hasGas                 → interface (real surface)
};

template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T,
          typename ScalarField_T >
void PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T, Reconstructor_T,
                               ScalarField_T >::operator()(IBlock* const block)
{
   WALBERLA_ASSERT_NOT_NULLPTR(block);

   auto* pdfField         = block->getData< PdfField_T >(pdfFieldID_);
   auto* boundaryHandling = block->getData< BoundaryHandling_T >(boundaryHandlingID_);
   auto* flagField        = boundaryHandling->getFlagField();
   auto* particleField    = block->getData< ParticleField_T >(particleFieldID_);

   // Fill-level field for Bogner φ seeding (null / unused unless ScalarField_T is provided).
   [[maybe_unused]] ScalarField_T* fillField = nullptr;
   if constexpr (!std::is_void_v< ScalarField_T >)
   {
      if (useBognerStateSelection_) fillField = block->getData< ScalarField_T >(fillFieldID_);
   }

   WALBERLA_ASSERT_NOT_NULLPTR(pdfField);
   WALBERLA_ASSERT_NOT_NULLPTR(boundaryHandling);
   WALBERLA_ASSERT_NOT_NULLPTR(flagField);
   WALBERLA_ASSERT_NOT_NULLPTR(particleField);

   WALBERLA_ASSERT_EQUAL(particleField->xyzSize(), flagField->xyzSize());

   WALBERLA_ASSERT(flagField->flagExists(formerObstacle_));
   WALBERLA_ASSERT(flagField->flagExists(fluid_));

   const flag_t formerObstacle = flagField->getFlag(formerObstacle_);
   const flag_t fluid          = flagField->getFlag(fluid_);

   // Pre-compute Bogner flag values at function scope so they are accessible in
   // both the reconstruction pass and the flag-update pass below.
   // Initialised to 0; only populated (and used) when useBognerStateSelection_ is true.
   flag_t bognerLiquid    = flag_t(0);
   flag_t bognerGas       = flag_t(0);
   flag_t bognerInterface = flag_t(0);

   if (useBognerStateSelection_)
   {
      WALBERLA_ASSERT(flagField->flagExists(liquidFlagUID_));
      WALBERLA_ASSERT(flagField->flagExists(gasFlagUID_));
      WALBERLA_ASSERT(flagField->flagExists(interfaceFlagUID_));

      bognerLiquid    = flagField->getFlag(liquidFlagUID_);
      bognerGas       = flagField->getFlag(gasFlagUID_);
      bognerInterface = flagField->getFlag(interfaceFlagUID_);
   }

   // Helper lambda: per-cell Bogner target flag.
   // Must be evaluated inside each cell loop — the neighbourhood state varies per cell.
   //
   // The neighbourhood is scanned over the full D3Q27 stencil (all 26 neighbours), not just
   // the 6 face neighbours: the LBM stencil couples to edge/corner neighbours as well, so a gas
   // cell sitting on a diagonal direction must be seen — otherwise a liquid cell could be placed
   // directly adjacent to gas, violating the FSLBM invariant (liquid must never be a stencil
   // neighbour of gas) and producing NaNs.
   //
   // A cell is only promoted to `liquid` when the neighbourhood provably contains no gas AND no
   // still-unclassified FormerMO neighbour (those are converted in the same pass and may yet
   // become gas). If either is present we fall back to the safe `interface` buffer, which the
   // free-surface sweeps promote to liquid on their own once it is genuinely enclosed.
   auto bognerTargetFlag = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> flag_t {
      bool hasLiquid  = false;
      bool hasGas     = false;
      bool hasUnknown = false;
      for (int dz = -1; dz <= 1; ++dz)
         for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
            {
               if (dx == 0 && dy == 0 && dz == 0) continue; // skip centre
               const flag_t nf = flagField->get(cx + dx, cy + dy, cz + dz);
               if (isFlagSet(nf, bognerLiquid) || isFlagSet(nf, bognerInterface)) hasLiquid = true;
               if (isFlagSet(nf, bognerGas) || isFlagSet(nf, bognerInterface)) hasGas = true;
               if (isFlagSet(nf, formerObstacle)) hasUnknown = true;
            }
      // inner-cell guard for the conversion counters: a FormerMO cell in a ghost layer is a
      // duplicate of another rank's inner cell, so only inner cells are tallied.
      const cell_idx_t xs = cell_idx_c(flagField->xSize());
      const cell_idx_t ys = cell_idx_c(flagField->ySize());
      const cell_idx_t zs = cell_idx_c(flagField->zSize());
      const bool inner = (cx >= 0 && cx < xs && cy >= 0 && cy < ys && cz >= 0 && cz < zs);

      if (!hasLiquid) { if (inner) ++convGas_; return bognerGas; }
      if (!hasGas && !hasUnknown)
      {
         // VERIFICATION of the under-resolved-waterline over-promotion hypothesis:
         // this cell qualifies as bulk liquid only on the 1-cell (D3Q27) stencil. Re-scan a 2-cell
         // radius — if the free surface (gas or interface) sits at distance 2, the cell is physically
         // at the waterline and is being OVER-PROMOTED to liquid (fill forced to 1) even though it is
         // within two cells of air. With adequate resolution an interface cell would sit at distance 1
         // (→ hasGas true → classified interface), so this print firing IS the resolution-driven
         // mismatch. It is invariant-legal (gas is >1 cell away → no liquid–gas adjacency, hence
         // checkClosedInterfaceLayer reports 0 holes) but physically wrong.
         const cell_idx_t gl = cell_idx_c(flagField->nrOfGhostLayers());
         bool gasAt2 = false, intfAt2 = false;
         for (int dz = -2; dz <= 2; ++dz)
            for (int dy = -2; dy <= 2; ++dy)
               for (int dx = -2; dx <= 2; ++dx)
               {
                  if (std::abs(dx) <= 1 && std::abs(dy) <= 1 && std::abs(dz) <= 1) continue; // D3Q27 already scanned
                  const cell_idx_t nx = cx + cell_idx_c(dx), ny = cy + cell_idx_c(dy), nz = cz + cell_idx_c(dz);
                  if (nx < -gl || nx >= xs + gl || ny < -gl || ny >= ys + gl || nz < -gl || nz >= zs + gl) continue;
                  const flag_t nf = flagField->get(nx, ny, nz);
                  if (isFlagSet(nf, bognerGas))       gasAt2  = true;
                  if (isFlagSet(nf, bognerInterface)) intfAt2 = true;
               }
         if (gasAt2 || intfAt2)
         {
            const auto&  aabb = block->getAABB();
            const real_t cdx  = aabb.xSize() / real_c(xs);
            WALBERLA_LOG_WARNING("[OVERPROMOTE-CHECK] FormerMO cell block-local (" << cx << ", " << cy << ", " << cz
               << ") world≈(" << (aabb.xMin() + (real_c(cx) + real_t(0.5)) * cdx) << ", "
               << (aabb.yMin() + (real_c(cy) + real_t(0.5)) * cdx) << ", "
               << (aabb.zMin() + (real_c(cz) + real_t(0.5)) * cdx)
               << ") promoted to LIQUID (fill→1), but the free surface is at distance 2:"
               << (gasAt2 ? " GAS" : "") << (intfAt2 ? " INTERFACE" : "")
               << "  --> under-resolved waterline OVER-PROMOTION (should be interface).");
         }
         if (inner) ++convLiquid_;
         return bognerInterface;
      }
      if (!hasGas && hasUnknown)
      {
         // Would qualify as bulk liquid, but an unclassified FormerMO neighbour could still
         // become gas this pass — downgrade to the safe interface buffer instead of liquid.
         WALBERLA_LOG_WARNING("Bogner state selection: cell ("
                              << cx << ", " << cy << ", " << cz
                              << ") has an unclassified FormerMO neighbour; falling back to interface "
                                 "instead of liquid to preserve the liquid/gas separation invariant.");
         if (inner) ++convIntfUnknown_;
         return bognerInterface;
      }
      if (inner) ++convIntfSurface_;
      return bognerInterface;
   };

   // ---- Pass 1: reconstruct missing PDFs (inner domain only) ------------------

   if (optimizeForSmallObstacleFraction_)
   {
      const uint_t numberOfGhostLayersToInclude = uint_t(0);

      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         auto cellInterval = getParticleCellBB(idx, *ac_, *block, *blockStorage_, numberOfGhostLayersToInclude);
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(
            cellInterval, if (isFlagSet(flagField->get(x, y, z), formerObstacle) &&
                              ac_->uidToIdx(particleField->get(x, y, z)) == idx) {
               (*reconstructor_)(block, x, y, z, pdfField, idx, *ac_);
               if (conserveMomentumWhenRestoring_ && pdfField->isInInnerPart(Cell(x, y, z)))
               {
                  Vector3< real_t > force      = pdfField->getMomentumDensity(x, y, z);
                  Vector3< real_t > cellCenter = blockStorage_->getBlockLocalCellCenter(*block, Cell(x, y, z));
                  lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic(idx, *ac_, -force, cellCenter);
               }
            })
      }
   }
   else
   {
      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(
         flagField->xyzSize(), if (isFlagSet(flagField->get(x, y, z), formerObstacle)) {
            size_t particleIdx = ac_->uidToIdx(particleField->get(x, y, z));
            WALBERLA_ASSERT_UNEQUAL(particleIdx, ac_->getInvalidIdx(), "Index of particle is invalid!");
            (*reconstructor_)(block, x, y, z, pdfField, particleIdx, *ac_);
            if (conserveMomentumWhenRestoring_ && pdfField->isInInnerPart(Cell(x, y, z)))
            {
               Vector3< real_t > force      = pdfField->getMomentumDensity(x, y, z);
               Vector3< real_t > cellCenter = blockStorage_->getBlockLocalCellCenter(*block, Cell(x, y, z));
               lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic(particleIdx, *ac_, -force, cellCenter);
            }
         })
   }

   // ---- Pass 2: update flags formerObstacle → target (inner domain + ghost layers),
   //             invalidate particle field entries --------------------------------

   if (optimizeForSmallObstacleFraction_)
   {
      const uint_t numberOfGhostLayersToInclude = flagField->nrOfGhostLayers();

      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         auto cellIntervalToCleanUp =
            getParticleCellBB(idx, *ac_, *block, *blockStorage_, numberOfGhostLayersToInclude);
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(
            cellIntervalToCleanUp, if (isFlagSet(flagField->get(x, y, z), formerObstacle)) {
               const flag_t targetFlag = useBognerStateSelection_ ? bognerTargetFlag(x, y, z) : fluid;
               boundaryHandling->setDomain(targetFlag, x, y, z);
               if constexpr (!std::is_void_v< ScalarField_T >)
                  if (fillField != nullptr)
                     seedConvertedCellFill(fillField, flagField, x, y, z, targetFlag, bognerGas, bognerLiquid,
                                           bognerInterface);
               removeFlag(flagField->get(x, y, z), formerObstacle);
               particleField->get(x, y, z) = ac_->getInvalidUid();
            })
      }
   }
   else
   {
      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(
         flagField->xyzSizeWithGhostLayer(), if (isFlagSet(flagField->get(x, y, z), formerObstacle)) {
            const flag_t targetFlag = useBognerStateSelection_ ? bognerTargetFlag(x, y, z) : fluid;
            boundaryHandling->setDomain(targetFlag, x, y, z);
            if constexpr (!std::is_void_v< ScalarField_T >)
               if (fillField != nullptr)
                  seedConvertedCellFill(fillField, flagField, x, y, z, targetFlag, bognerGas, bognerLiquid,
                                        bognerInterface);
            removeFlag(flagField->get(x, y, z), formerObstacle);
            particleField->get(x, y, z) = ac_->getInvalidUid();
         })
   }
}

// ---- Free factory functions ------------------------------------------------

template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T >
shared_ptr< PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T, Reconstructor_T > >
   makePdfReconstructionManager(const shared_ptr< StructuredBlockStorage >& blockStorage, const BlockDataID& pdfFieldID,
                                const BlockDataID& boundaryHandlingID, const BlockDataID& particleFieldID,
                                const shared_ptr< ParticleAccessor_T >& ac, const FlagUID& formerObstacle,
                                const FlagUID& fluid, const shared_ptr< Reconstructor_T >& reconstructor,
                                bool conserveMomentumWhenRestoring, const bool optimizeForSmallObstacleFraction = false)
{
   using RM_T = PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T, Reconstructor_T >;
   return make_shared< RM_T >(RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac, formerObstacle,
                                   fluid, reconstructor, conserveMomentumWhenRestoring,
                                   optimizeForSmallObstacleFraction));
}

template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T >
shared_ptr< PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                      EquilibriumReconstructor< BoundaryHandling_T > > >
   makePdfReconstructionManager(const shared_ptr< StructuredBlockStorage >& blockStorage, const BlockDataID& pdfFieldID,
                                const BlockDataID& boundaryHandlingID, const BlockDataID& particleFieldID,
                                const shared_ptr< ParticleAccessor_T >& ac, const FlagUID& formerObstacle,
                                const FlagUID& fluid, bool conserveMomentumWhenRestoring,
                                const bool optimizeForSmallObstacleFraction = false)
{
   using RM_T = PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                          EquilibriumReconstructor< BoundaryHandling_T > >;
   return make_shared< RM_T >(RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac, formerObstacle,
                                   fluid,
                                   makeEquilibriumReconstructor< BoundaryHandling_T >(blockStorage, boundaryHandlingID),
                                   conserveMomentumWhenRestoring, optimizeForSmallObstacleFraction));
}

// Bogner state-selection overload: determines gas/liquid/interface per cell from
// the non-obstacle neighbourhood B(x) rather than assigning a uniform target flag, and seeds
// the volume fraction φ of converted cells from `fillFieldID` (Bogner & Rüde 2013, §5.2).
// ScalarField_T (the fill-level field type) is given explicitly; ParticleAccessor_T is deduced.
template< typename PdfField_T, typename BoundaryHandling_T, typename ScalarField_T, typename ParticleAccessor_T >
shared_ptr< PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                      EquilibriumReconstructor< BoundaryHandling_T >, ScalarField_T > >
   makePdfReconstructionManagerBogner(const shared_ptr< StructuredBlockStorage >& blockStorage,
                                      const BlockDataID& pdfFieldID, const BlockDataID& boundaryHandlingID,
                                      const BlockDataID& particleFieldID, const shared_ptr< ParticleAccessor_T >& ac,
                                      const FlagUID& formerObstacle, const FlagUID& fluid,
                                      bool conserveMomentumWhenRestoring, const FlagUID& liquidFlagUID,
                                      const FlagUID& gasFlagUID, const FlagUID& interfaceFlagUID,
                                      const BlockDataID& fillFieldID, const bool useBognerStateSelection = true,
                                      const bool optimizeForSmallObstacleFraction = false)
{
   using RM_T = PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                          EquilibriumReconstructor< BoundaryHandling_T >, ScalarField_T >;
   return make_shared< RM_T >(RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac, formerObstacle,
                                   fluid,
                                   makeEquilibriumReconstructor< BoundaryHandling_T >(blockStorage, boundaryHandlingID),
                                   conserveMomentumWhenRestoring, liquidFlagUID, gasFlagUID, interfaceFlagUID,
                                   fillFieldID, useBognerStateSelection, optimizeForSmallObstacleFraction));
}

} // namespace lbm_mesapd_coupling
} // namespace walberla