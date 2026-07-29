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

#include <array>
#include <functional>
#include <sstream>

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

//! Result of classifying a freshly-uncovered ("formerObstacle") cell's D3Q27 neighbourhood in the
//! Bogner & Rüde state-selection reconstruction.
enum class BognerNeighbourhoodState { Gas, Liquid, Interface };

//**************************************************************************************************************************************
/*!\brief Pure decision logic for the Bogner state-selection classifier (bognerTargetFlag).
 *
 * Given only whether the D3Q27 neighbourhood contains a REAL gas cell, a REAL liquid cell, and/or
 * an unresolved ("formerObstacle") neighbour, decide the target state for the cell itself:
 *  - hasGas && !hasLiquid && !hasUnknown -> Gas       (provably bulk gas, fully resolved)
 *  - hasGas && !hasLiquid &&  hasUnknown -> Interface (an unresolved neighbour could still become
 *                                                       liquid later in the same reconstruction pass)
 *  - hasLiquid && !hasGas && !hasUnknown -> Liquid   (provably bulk liquid, fully resolved)
 *  - hasLiquid && !hasGas &&  hasUnknown -> Interface (an unresolved neighbour could still become
 *                                                       gas later in the same reconstruction pass)
 *  - hasGas && hasLiquid                -> Interface (genuine mixed/free-surface neighbourhood)
 *  - !hasGas && !hasLiquid               -> Interface (neither real phase reachable within one
 *                                                       cell -- e.g. only interface/unresolved
 *                                                       neighbours; safest fallback)
 *
 * The gas and liquid branches are symmetric on purpose: an unresolved same-pass neighbour is exactly
 * as dangerous from the gas side as from the liquid side -- promoting to gas now, only for the
 * neighbour to resolve to liquid a few cells later in the same pass, produces the same illegal
 * direct liquid/gas contact as the mirror-image case, just discovered from the other side.
 *
 * hasGas/hasLiquid must be computed from REAL bognerGas/bognerLiquid flags only -- an interface
 * neighbour must NOT be counted toward either. Interface cells are the legal buffer between liquid
 * and gas (a liquid cell next to an interface cell is an ordinary free surface), so counting them
 * as gas-equivalent enforces nothing beyond what checking real gas alone already guarantees, while
 * causing a self-propagating misclassification in densely packed, actively-jittering particle beds
 * (one interface cell seeded near a contact point makes its neighbour's hasGas true next pass, and
 * so on outward through the narrow gaps between particles, with no real gas anywhere nearby).
 */
//**************************************************************************************************************************************
inline BognerNeighbourhoodState classifyFormerObstacleNeighbourhood(bool hasGas, bool hasLiquid, bool hasUnknown)
{
   if (hasGas && !hasLiquid) return hasUnknown ? BognerNeighbourhoodState::Interface : BognerNeighbourhoodState::Gas;
   if (hasLiquid && !hasGas) return hasUnknown ? BognerNeighbourhoodState::Interface : BognerNeighbourhoodState::Liquid;
   // hasGas && hasLiquid, or neither -> interface in both cases.
   return BognerNeighbourhoodState::Interface;
}

//**************************************************************************************************************************************
/*!\brief Scans a cell's D3Q27 neighbourhood and classifies it (bognerTargetFlag's actual per-cell
 * logic, factored out so it can be driven by a literal neighbourhood -- e.g. in tests -- as well as
 * by a real FlagField).
 *
 * getNeighbourFlag(dx, dy, dz) must return the flag at the given offset from the centre cell for
 * every (dx, dy, dz) in {-1,0,1}^3 except (0,0,0) (the centre itself is never queried). Flag_T only
 * needs to support isFlagSet(value, testFlag) (field::isFlagSet, an integer bitmask test) -- the
 * real production flag_t (BoundaryHandling_T::flag_t) and a small test-local flag type both work.
 *
 * hasGas/hasLiquid/hasUnknown/unknownNbr are written back to the caller (needed for the per-branch
 * conversion counters and the same-pass-race warning, which need the raw booleans, not just the
 * final Gas/Liquid/Interface classification).
 */
//**************************************************************************************************************************************
template< typename Flag_T, typename NeighbourFlagGetter >
inline BognerNeighbourhoodState classifyFormerObstacleCell(NeighbourFlagGetter getNeighbourFlag,
                                                             Flag_T bognerLiquid, Flag_T bognerGas,
                                                             Flag_T formerObstacle, bool& hasGasOut,
                                                             bool& hasLiquidOut, bool& hasUnknownOut,
                                                             uint_t& unknownNbrOut)
{
   hasGasOut = false; hasLiquidOut = false; hasUnknownOut = false; unknownNbrOut = uint_t(0);
   for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
         for (int dx = -1; dx <= 1; ++dx)
         {
            if (dx == 0 && dy == 0 && dz == 0) continue; // skip centre
            const Flag_T nf = getNeighbourFlag(dx, dy, dz);
            if (isFlagSet(nf, bognerLiquid))   hasLiquidOut = true;
            if (isFlagSet(nf, bognerGas))      hasGasOut    = true;
            if (isFlagSet(nf, formerObstacle)) { hasUnknownOut = true; ++unknownNbrOut; }
         }
   return classifyFormerObstacleNeighbourhood(hasGasOut, hasLiquidOut, hasUnknownOut);
}

//**************************************************************************************************************************************
/*!\brief Result of computeRealPhaseComponents() (see PdfReconstructionManager's "[WATERLINE-CHECK]"
 * usage), answering -- for every cell in a scanned region -- whether a real bognerGas / bognerLiquid
 * cell is reachable from it without crossing the OPPOSITE real phase.
 *
 * Used by bognerTargetFlag to resolve cells that would otherwise fall back to interface (an
 * unresolved same-pass neighbour, or no real phase within 1 cell) whenever it can prove -- from
 * real flags only, via genuine connectivity, never a guess -- that the dangerous opposite phase is
 * not reachable at all. Initially landed as a log-only shadow check (see the commit immediately
 * preceding this one) with the intent to validate it against a real cluster run before wiring it
 * into an actual decision; wired in directly instead on the reasoning that a log-only pass, with
 * every cell still forced to interface, wouldn't produce a meaningful comparison anyway -- not yet
 * confirmed against a real cluster run.
 */
//**************************************************************************************************************************************
struct PhaseReachability
{
   cell_idx_t xMin_ = cell_idx_t(0), yMin_ = cell_idx_t(0), zMin_ = cell_idx_t(0);
   cell_idx_t xSize_ = cell_idx_t(0), ySize_ = cell_idx_t(0), zSize_ = cell_idx_t(0);

   std::vector< uint8_t > gasConfidentlyReachable_;
   std::vector< uint8_t > gasInconclusive_;
   std::vector< uint8_t > liquidConfidentlyReachable_;
   std::vector< uint8_t > liquidInconclusive_;

   uint_t totalGasComponents    = uint_t(0); // # of distinct components that contain real gas
   uint_t totalLiquidComponents = uint_t(0); // # of distinct components that contain real liquid

   inline size_t flatIndex(cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      return (size_t(z - zMin_) * size_t(ySize_) + size_t(y - yMin_)) * size_t(xSize_) + size_t(x - xMin_);
   }
   inline bool gasConfidentlyReachable(cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      return gasConfidentlyReachable_[flatIndex(x, y, z)] != uint8_t(0);
   }
   inline bool gasInconclusive(cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      return gasInconclusive_[flatIndex(x, y, z)] != uint8_t(0);
   }
   inline bool liquidConfidentlyReachable(cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      return liquidConfidentlyReachable_[flatIndex(x, y, z)] != uint8_t(0);
   }
   inline bool liquidInconclusive(cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      return liquidInconclusive_[flatIndex(x, y, z)] != uint8_t(0);
   }
};

//**************************************************************************************************************************************
/*!\brief Computes, for every cell in [xMin,xMax]x[yMin,yMax]x[zMin,zMax] (all
 * inclusive, matching the same absolute-coordinate convention as FlagField::get(x,y,z) so
 * production can pass a bare field-access lambda), whether real gas / real liquid is reachable
 * from it without crossing the opposite real phase.
 *
 * The flood fill is free to pass through any cell that is NEITHER real bognerGas NOR real
 * bognerLiquid (formerObstacle, interface, particle-obstacle cells, walls, ...) -- such cells carry
 * no evidence of either phase, so: (a) a floating object sitting in open liquid does not split that
 * liquid into two separate reachability regions (there is a permeable path through the object's own
 * cells connecting both sides), and (b) a cluster of merely-uncovered formerObstacle cells does not
 * itself look like a waterline (it contributes nothing until the flood fill reaches an actual real
 * gas or real liquid cell on the far side).
 *
 * "Confidently reachable" means a real cell of that phase was actually found, connected via only
 * permeable cells. "Inconclusive" means the flood fill for that phase ran all the way to the outer
 * boundary of the scanned region without finding one -- i.e. the true situation might extend further
 * than this region reaches, so "not reachable" cannot be trusted there.
 *
 * getFlag(x, y, z) must return the flag at the given ABSOLUTE coordinates for every cell in the
 * scanned region.
 */
//**************************************************************************************************************************************
template< typename Flag_T, typename FlagGetter >
inline PhaseReachability computeRealPhaseComponents(FlagGetter getFlag, cell_idx_t xMin, cell_idx_t xMax,
                                                     cell_idx_t yMin, cell_idx_t yMax, cell_idx_t zMin,
                                                     cell_idx_t zMax, Flag_T bognerLiquid, Flag_T bognerGas)
{
   PhaseReachability result;
   result.xMin_  = xMin;
   result.yMin_  = yMin;
   result.zMin_  = zMin;
   result.xSize_ = xMax - xMin + cell_idx_t(1);
   result.ySize_ = yMax - yMin + cell_idx_t(1);
   result.zSize_ = zMax - zMin + cell_idx_t(1);

   const size_t nx = size_t(result.xSize_);
   const size_t ny = size_t(result.ySize_);
   const size_t nz = size_t(result.zSize_);
   const size_t n  = nx * ny * nz;

   result.gasConfidentlyReachable_.assign(n, uint8_t(0));
   result.gasInconclusive_.assign(n, uint8_t(0));
   result.liquidConfidentlyReachable_.assign(n, uint8_t(0));
   result.liquidInconclusive_.assign(n, uint8_t(0));

   auto toAbs = [&](size_t flat, cell_idx_t& x, cell_idx_t& y, cell_idx_t& z) {
      x = xMin + cell_idx_c(flat % nx);
      y = yMin + cell_idx_c((flat / nx) % ny);
      z = zMin + cell_idx_c(flat / (nx * ny));
   };
   auto flatOf = [&](cell_idx_t x, cell_idx_t y, cell_idx_t z) -> size_t {
      return (size_t(z - zMin) * ny + size_t(y - yMin)) * nx + size_t(x - xMin);
   };

   // One "graph": nodes are every cell that is NOT `blockingFlag` (i.e. the opposite real phase);
   // edges are 26-connectivity. A resulting component is "productive" if it contains at least one
   // real `productiveFlag` cell, and "touches the edge" if it reaches the outer boundary of the
   // scanned region. Writes confidence/inconclusive results for every cell that IS part of this
   // graph (cells that ARE `blockingFlag` are simply left at 0/0 -- irrelevant for this phase).
   auto runGraph = [&](Flag_T blockingFlag, Flag_T productiveFlag, std::vector< uint8_t >& confidentOut,
                       std::vector< uint8_t >& inconclusiveOut) -> uint_t {
      std::vector< int32_t > label(n, int32_t(-1));
      std::vector< uint8_t > productive;
      std::vector< uint8_t > touchesEdge;
      std::vector< size_t > stack;

      for (size_t startFlat = 0; startFlat < n; ++startFlat)
      {
         if (label[startFlat] != int32_t(-1)) continue;
         cell_idx_t sx, sy, sz;
         toAbs(startFlat, sx, sy, sz);
         if (isFlagSet(getFlag(sx, sy, sz), blockingFlag)) continue; // not part of this graph

         const int32_t myLabel = int32_t(productive.size());
         productive.push_back(uint8_t(0));
         touchesEdge.push_back(uint8_t(0));
         label[startFlat] = myLabel;
         stack.clear();
         stack.push_back(startFlat);

         while (!stack.empty())
         {
            const size_t curFlat = stack.back();
            stack.pop_back();
            cell_idx_t cx, cy, cz;
            toAbs(curFlat, cx, cy, cz);
            const Flag_T cf = getFlag(cx, cy, cz);
            if (isFlagSet(cf, productiveFlag)) productive[size_t(myLabel)] = uint8_t(1);
            if (cx == xMin || cx == xMax || cy == yMin || cy == yMax || cz == zMin || cz == zMax)
               touchesEdge[size_t(myLabel)] = uint8_t(1);

            for (int dz = -1; dz <= 1; ++dz)
               for (int dy = -1; dy <= 1; ++dy)
                  for (int dx = -1; dx <= 1; ++dx)
                  {
                     if (dx == 0 && dy == 0 && dz == 0) continue;
                     const cell_idx_t nx_ = cx + cell_idx_c(dx);
                     const cell_idx_t ny_ = cy + cell_idx_c(dy);
                     const cell_idx_t nz_ = cz + cell_idx_c(dz);
                     if (nx_ < xMin || nx_ > xMax || ny_ < yMin || ny_ > yMax || nz_ < zMin || nz_ > zMax) continue;
                     const size_t nFlat = flatOf(nx_, ny_, nz_);
                     if (label[nFlat] != int32_t(-1)) continue;
                     if (isFlagSet(getFlag(nx_, ny_, nz_), blockingFlag)) continue;
                     label[nFlat] = myLabel;
                     stack.push_back(nFlat);
                  }
         }
      }

      uint_t totalProductive = uint_t(0);
      for (uint8_t p : productive)
         if (p) ++totalProductive;

      for (size_t flat = 0; flat < n; ++flat)
      {
         if (label[flat] == int32_t(-1)) continue; // this cell IS the blocking phase -- not part of this graph
         const size_t c = size_t(label[flat]);
         if (productive[c]) confidentOut[flat] = uint8_t(1);
         else if (touchesEdge[c]) inconclusiveOut[flat] = uint8_t(1);
      }
      return totalProductive;
   };

   result.totalGasComponents = runGraph(bognerLiquid, bognerGas, result.gasConfidentlyReachable_,
                                        result.gasInconclusive_);
   result.totalLiquidComponents = runGraph(bognerGas, bognerLiquid, result.liquidConfidentlyReachable_,
                                           result.liquidInconclusive_);

   return result;
}

//! Verdict of the waterline-reachability check for a single cell. WouldBeLiquid/WouldBeGas are
//! committed directly by bognerTargetFlag; the other three still fall back to interface.
enum class WaterlineVerdict
{
   WouldBeLiquid,
   WouldBeGas,
   StillAmbiguousSurfaceNearby,  // both real gas and real liquid confidently reachable -- a true local surface
   StillAmbiguousNoEvidence,     // neither reachable -- no real evidence of either phase within the scanned region
   ExcludedNearUnresolvedEdge    // the scanned region ended before ruling out gas and/or liquid; not trustworthy
};

//**************************************************************************************************************************************
/*!\brief The waterline-reachability decision rule, given one cell's reachability results from
 * computeRealPhaseComponents().
 */
//**************************************************************************************************************************************
inline WaterlineVerdict classifyViaWaterlineReachability(bool gasConfidentlyReachable, bool gasInconclusive,
                                                          bool liquidConfidentlyReachable, bool liquidInconclusive)
{
   if (gasInconclusive || liquidInconclusive) return WaterlineVerdict::ExcludedNearUnresolvedEdge;
   if (liquidConfidentlyReachable && !gasConfidentlyReachable) return WaterlineVerdict::WouldBeLiquid;
   if (gasConfidentlyReachable && !liquidConfidentlyReachable) return WaterlineVerdict::WouldBeGas;
   if (gasConfidentlyReachable && liquidConfidentlyReachable) return WaterlineVerdict::StillAmbiguousSurfaceNearby;
   return WaterlineVerdict::StillAmbiguousNoEvidence;
}

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
                                   const bool optimizeForSmallObstacleFraction = false,
                                   const bool useHydrostaticReconstruction = false,
                                   const real_t hydrostaticAccelerationZ = real_t(0),
                                   const real_t absoluteLiquidHeight = real_t(0))
      : blockStorage_(blockStorage), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID),
        particleFieldID_(particleFieldID), ac_(ac), reconstructor_(reconstructor), formerObstacle_(formerObstacle),
        fluid_(fluid), conserveMomentumWhenRestoring_(conserveMomentumWhenRestoring),
        optimizeForSmallObstacleFraction_(optimizeForSmallObstacleFraction),
        useBognerStateSelection_(useBognerStateSelection),
        useHydrostaticReconstruction_(useHydrostaticReconstruction),
        hydrostaticAccelerationZ_(hydrostaticAccelerationZ), absoluteLiquidHeight_(absoluteLiquidHeight),
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
      std::vector< uint_t > c = { convGas_, convLiquid_, convIntfUnknown_, convIntfUnknownGas_, convIntfSurface_, convIntfUnresolved_ };
      walberla::mpi::reduceInplace(c, walberla::mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-CONVERT t=" << t << "] cells converted since last report: "
         << "gas(A)=" << c[0] << " liquid(B)=" << c[1]
         << " interface-unknown(C)=" << c[2] << " interface-unknown-gas(A')=" << c[3] << " interface-surface(D)=" << c[4]
         << " interface-unresolved(E)=" << c[5]
         << " | total=" << (c[0] + c[1] + c[2] + c[3] + c[4] + c[5]));

      // Phantom-force investigation: how many unknowns, and how clustered are they?
      std::vector< uint_t > u    = { formerMoTotal_, formerMoWithUnknownNbr_, unknownNbrSum_ };
      std::vector< uint_t > umax = { unknownNbrMax_ };
      std::vector< uint_t > hist(unknownNbrHist_.begin(), unknownNbrHist_.end());
      walberla::mpi::reduceInplace(u, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(umax, walberla::mpi::MAX);
      walberla::mpi::reduceInplace(hist, walberla::mpi::SUM);
      WALBERLA_ROOT_SECTION()
      {
         const uint_t tot     = u[0];
         const real_t pctWith = (tot > uint_t(0)) ? real_t(100) * real_c(u[1]) / real_c(tot) : real_t(0);
         const real_t meanNbr = (tot > uint_t(0)) ? real_c(u[2]) / real_c(tot) : real_t(0);
         std::ostringstream histStr;
         for (uint_t k = uint_t(0); k < hist.size(); ++k)
            if (hist[k] > uint_t(0)) histStr << " " << k << ":" << hist[k];
         WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-UNKNOWNS t=" << t << "] total unknown (FormerMO) cells=" << tot
            << " | with >=1 unknown neighbour=" << u[1] << " (" << pctWith << "%)"
            << " | mean unknown-neighbours/cell=" << meanNbr
            << " | max=" << umax[0]
            << " | hist[#unknownNbr:count]=" << histStr.str());
      }

      convGas_ = convLiquid_ = convIntfUnknown_ = convIntfUnknownGas_ = convIntfSurface_ = convIntfUnresolved_ = uint_t(0);
      formerMoTotal_ = formerMoWithUnknownNbr_ = unknownNbrSum_ = unknownNbrMax_ = uint_t(0);
      unknownNbrHist_.fill(uint_t(0));

      // [WATERLINE-CHECK] aggregate report. waterlineResolvedLiquid_/waterlineResolvedGas_ are cells
      // ACTUALLY committed to liquid/gas by the phase-connectivity check (see bognerTargetFlag)
      // instead of falling back to interface; the remaining three counters are cells that still
      // fell back to interface, broken down by why the check couldn't resolve them.
      std::vector< uint_t > s = { waterlineResolvedLiquid_,       waterlineResolvedGas_, shadowStillAmbiguousSurface_,
                                  shadowStillAmbiguousNoEvidence_, shadowExcludedNearEdge_,
                                  blockGasComponentCountSum_, blockLiquidComponentCountSum_, blockScanCount_ };
      walberla::mpi::reduceInplace(s, walberla::mpi::SUM);
      WALBERLA_ROOT_SECTION()
      {
         const uint_t blocks   = s[7];
         const real_t meanGasComp = (blocks > uint_t(0)) ? real_c(s[5]) / real_c(blocks) : real_t(0);
         const real_t meanLiqComp = (blocks > uint_t(0)) ? real_c(s[6]) / real_c(blocks) : real_t(0);
         WALBERLA_LOG_INFO_ON_ROOT("[WATERLINE-CHECK-SUMMARY t=" << t << "] resolvedLiquid=" << s[0]
            << " resolvedGas=" << s[1] << " stillAmbiguousSurface=" << s[2] << " stillAmbiguousNoEvidence=" << s[3]
            << " excludedNearEdge=" << s[4] << " | mean gas components/block=" << meanGasComp
            << " mean liquid components/block=" << meanLiqComp << " (over " << blocks << " block scans)");
      }
      waterlineResolvedLiquid_ = waterlineResolvedGas_ = shadowStillAmbiguousSurface_ = shadowStillAmbiguousNoEvidence_ =
         shadowExcludedNearEdge_ = uint_t(0);
      blockGasComponentCountSum_ = blockLiquidComponentCountSum_ = blockScanCount_ = uint_t(0);
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
   // Depth-consistent (hydrostatic) reconstruction density. The EquilibriumReconstructor seeds a freed
   // cell with the UNWEIGHTED mean of its domain neighbours' density, which is fill-blind: interface
   // cells report ρ≈1 however empty they are, so a cell uncovered AT or ABOVE the waterline is assigned
   // a full-liquid ρ≈1 while the correct hydrostatic value ρ_hyd(z) = 1 + 3·a_z·(z − H) is < 1 there.
   // The MEM then reads ρ − ρ_hyd > 0 → spurious upward (phantom) lift. When enabled we re-seed the cell
   // to ρ_hyd(z) (a_z = hydrostaticAccelerationZ_ < 0, H = absoluteLiquidHeight_). This uses no fill
   // field, so it is robust to the negative contact-line fill drift that makes Σφ-weighting degenerate.
   // Off by default → identical to the legacy path when unset.
   const bool useHydrostaticReconstruction_ = false;
   const real_t hydrostaticAccelerationZ_   = real_t(0); // a_z (≤ 0); same value fed to HydrostaticDensityFtor
   const real_t absoluteLiquidHeight_       = real_t(0); // H: undisturbed free-surface height in lattice units
   const FlagUID liquidFlagUID_;
   const FlagUID gasFlagUID_;
   const FlagUID interfaceFlagUID_;
   const BlockDataID fillFieldID_; // fill-level field for φ seeding (Bogner path only)

   // FormerMO → state conversion counters, one per branch of bognerTargetFlag (see the table in
   // contact_ring_crash_analysis.md §5.3). Accumulated over inner cells across blocks; flushed by
   // logAndResetConversionCounts. Mutable-free: operator() is non-const, the counting lambda
   // captures `this`.
   uint_t convGas_            = uint_t(0); // branch A:  hasGas && !hasLiquid && !hasUnknown → gas
   uint_t convIntfUnknownGas_ = uint_t(0); // branch A': hasGas && !hasLiquid &&  hasUnknown → interface (fallback, symmetric to C)
   uint_t convLiquid_         = uint_t(0); // branch B:  hasLiquid && !hasGas && !hasUnknown → liquid
   uint_t convIntfUnknown_    = uint_t(0); // branch C:  hasLiquid && !hasGas &&  hasUnknown → interface (fallback)
   uint_t convIntfSurface_    = uint_t(0); // branch D:  hasGas &&  hasLiquid                → interface (real surface)
   uint_t convIntfUnresolved_ = uint_t(0); // branch E: !hasGas && !hasLiquid                → interface (fallback, no real fluid neighbour yet)

   // Unknown-neighbour diagnostics (phantom-force investigation). For every inner FormerMO cell we
   // count how many of its 26 D3Q27 neighbours are themselves still-unknown FormerMO cells — i.e.
   // cells whose PDFs are *also* being reconstructed this same pass. This tells us whether the
   // unknowns are isolated (a reconstructed cell leans on well-defined fluid/gas neighbours) or
   // clustered (it leans on neighbours that are themselves unknown → reconstruction error compounds,
   // the suspected phantom-force amplifier). Flushed together with the conversion counts.
   uint_t formerMoTotal_          = uint_t(0); // total inner FormerMO cells classified this period
   uint_t formerMoWithUnknownNbr_ = uint_t(0); // of those, how many have >=1 unknown FormerMO neighbour
   uint_t unknownNbrSum_          = uint_t(0); // sum of unknown-neighbour counts (→ mean per cell)
   uint_t unknownNbrMax_          = uint_t(0); // largest unknown-neighbour count seen on any cell
   std::array< uint_t, 27 > unknownNbrHist_ = {}; // histogram indexed by #unknown neighbours (0..26)

   // [WATERLINE-CHECK] counters (see PhaseReachability/computeRealPhaseComponents above), tallied
   // for every branch-C/A'/E cell (the ones that used to fall back to interface unconditionally).
   // waterlineResolvedLiquid_/waterlineResolvedGas_ count cells the phase-connectivity check
   // resolved and ACTUALLY committed (see bognerTargetFlag) -- a real behaviour change, not a
   // diagnostic. The remaining three count cells the check could NOT resolve (still fell back to
   // interface exactly as before), broken down by why.
   uint_t waterlineResolvedLiquid_            = uint_t(0);
   uint_t waterlineResolvedGas_               = uint_t(0);
   uint_t shadowStillAmbiguousSurface_    = uint_t(0);
   uint_t shadowStillAmbiguousNoEvidence_ = uint_t(0);
   uint_t shadowExcludedNearEdge_         = uint_t(0);
   uint_t blockGasComponentCountSum_      = uint_t(0); // Σ totalGasComponents over every block scanned
   uint_t blockLiquidComponentCountSum_   = uint_t(0); // Σ totalLiquidComponents over every block scanned
   uint_t blockScanCount_                 = uint_t(0); // # of blocks scanned (→ mean components/block)
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

   // Depth-consistent reconstruction (useHydrostaticReconstruction_): after the reconstructor sets a
   // freed cell to PURE equilibrium at the fill-blind mean neighbour density, re-seed it to the
   // hydrostatic density ρ_hyd(z) = 1 + 3·a_z·(z − H) at the cell centre. f_i^eq ∝ ρ at fixed u, so
   // scaling every PDF by ρ_hyd/ρ_assigned yields f_i^eq(ρ_hyd, u_particle) with the velocity preserved
   // exactly. In the bulk ρ_assigned ≈ ρ_hyd → near-no-op; at/above the waterline it removes the ρ≈1
   // over-density that the MEM reads as phantom lift. Uses no fill field → robust to the negative
   // contact-line fill drift that made Σφ-weighting degenerate. Call right after the reconstructor.
   auto reseedHydrostatic = [&](const cell_idx_t x, const cell_idx_t y, const cell_idx_t z) {
      if (!useHydrostaticReconstruction_) return;
      const real_t rhoAssigned = pdfField->getDensity(x, y, z);
      if (rhoAssigned <= real_t(1e-6)) return;
      const Vector3< real_t > c = blockStorage_->getBlockLocalCellCenter(*block, Cell(x, y, z));
      const real_t rhoHyd = real_t(1) + real_t(3) * hydrostaticAccelerationZ_ * (c[2] - absoluteLiquidHeight_);
      if (rhoHyd <= real_t(0)) return; // above the column where the linear profile underflows: leave as-is
      const real_t scale = rhoHyd / rhoAssigned;
      for (uint_t fIdx = uint_t(0); fIdx < pdfField->fSize(); ++fIdx)
         pdfField->get(x, y, z, fIdx) *= scale;
   };

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

   // [WATERLINE-SHADOW] diagnostic-only phase-connectivity check (PhaseReachability /
   // computeRealPhaseComponents, above). Computed ONCE per block, over the block's full
   // ghost-inclusive extent (2 ghost layers already exist for this flag field -- see
   // AntidunesBoundaryHandling.impl.h -- so this needs no new communication). Consulted below by
   // bognerTargetFlag purely to log what a wider, connectivity-aware look would have concluded;
   // never changes what a cell actually commits to.
   PhaseReachability waterlineShadow;
   if (useBognerStateSelection_)
   {
      const auto ci = flagField->xyzSizeWithGhostLayer();
      waterlineShadow = computeRealPhaseComponents(
         [&](cell_idx_t x, cell_idx_t y, cell_idx_t z) { return flagField->get(x, y, z); }, ci.xMin(), ci.xMax(),
         ci.yMin(), ci.yMax(), ci.zMin(), ci.zMax(), bognerLiquid, bognerGas);
      blockGasComponentCountSum_ += waterlineShadow.totalGasComponents;
      blockLiquidComponentCountSum_ += waterlineShadow.totalLiquidComponents;
      ++blockScanCount_;
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
   // hasGas/hasLiquid count ONLY real bognerGas/bognerLiquid neighbours — an interface neighbour
   // counts toward NEITHER. Interface cells are the legal buffer between liquid and gas (a liquid
   // cell sitting next to an interface cell is exactly what a normal free surface looks like), so
   // there is nothing to protect by treating "interface nearby" as gas-equivalent. Counting it that
   // way (as an earlier version of this code did) causes a self-propagating misclassification in
   // densely packed, actively-jittering particle beds: one interface cell seeded near a contact
   // point makes its neighbour's hasGas true on the next FormerMO pass, which makes THAT cell
   // interface too, and so on outward through the narrow gaps between packed particles — with no
   // real gas anywhere nearby. Confirmed empirically: buried floating cylinders in
   // cylinder_sedimentation_showcase_gas.prm showed persistent locally-fitted "interface" cells and
   // oscillating (17-86%) submersion 30+ LU below the true free surface (~55 LU), long after
   // settling into the packed bed (see changes.md).
   //
   // A cell is only promoted to `liquid` when the neighbourhood provably contains real liquid, no
   // real gas, AND no still-unclassified FormerMO neighbour (those are converted in the same pass
   // and may yet become gas). If neither real gas nor real liquid is found at all (only interface
   // and/or unresolved FormerMO neighbours — e.g. a cell whose entire 1-cell neighbourhood has
   // itself just been through this same classifier), neither promotion is provably safe, so it
   // falls back to the safe `interface` buffer; the ordinary free-surface conversion sweep resolves
   // it later once a real liquid/gas neighbour actually becomes reachable.
   //
   // formerObstacleGhostBounds: the field's full ghost-inclusive extent, computed once here so the
   // neighbour-read lambda below can bounds-check against it. Needed because the Pass-2 update loop
   // (below) iterates x/y/z over this SAME full ghost-inclusive range, but classifyFormerObstacleCell
   // reads a further +/-1 cell around whichever cell it's given -- a formerObstacle cell sitting at
   // the OUTERMOST ghost layer therefore has a neighbour one cell beyond the field's actual allocated
   // buffer. Confirmed with AddressSanitizer on a real cluster run (27-floater mixed-shape pillars
   // showcase, job 3893620): SEGV, invalid READ, exactly at the flagField->get() call directly below.
   const auto formerObstacleGhostBounds = flagField->xyzSizeWithGhostLayer();
   auto bognerTargetFlag = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> flag_t {
      // Scan + decision both live in classifyFormerObstacleCell()/classifyFormerObstacleNeighbourhood()
      // (top of this file) as pure, directly unit-tested functions
      // (tests/showcases_antidunes/BognerStateClassificationTest.cpp) — this lambda only supplies
      // the real field access and attaches the side effects (per-branch counters, the same-pass-race
      // warning) that need cx/cy/cz/inner and can't live in the pure functions themselves.
      bool hasLiquid, hasGas, hasUnknown;
      uint_t unknownNbr;
      const BognerNeighbourhoodState state = classifyFormerObstacleCell(
         [&](int dx, int dy, int dz) -> flag_t {
            const cell_idx_t nx = cx + cell_idx_c(dx);
            const cell_idx_t ny = cy + cell_idx_c(dy);
            const cell_idx_t nz = cz + cell_idx_c(dz);
            // Neighbour falls outside the field's allocated (ghost-inclusive) range: we cannot read
            // it safely. Report it as an unresolved ("formerObstacle") neighbour -- the SAME
            // conservative fallback this classifier already uses for a same-pass-unresolved neighbour
            // (falls back to interface, the always-legal state) -- rather than reading invalid memory.
            if (nx < formerObstacleGhostBounds.xMin() || nx > formerObstacleGhostBounds.xMax() ||
                ny < formerObstacleGhostBounds.yMin() || ny > formerObstacleGhostBounds.yMax() ||
                nz < formerObstacleGhostBounds.zMin() || nz > formerObstacleGhostBounds.zMax())
               return formerObstacle;
            return flagField->get(nx, ny, nz);
         },
         bognerLiquid, bognerGas, formerObstacle, hasGas, hasLiquid, hasUnknown, unknownNbr);

      // inner-cell guard for the conversion counters: a FormerMO cell in a ghost layer is a
      // duplicate of another rank's inner cell, so only inner cells are tallied.
      const cell_idx_t xs = cell_idx_c(flagField->xSize());
      const cell_idx_t ys = cell_idx_c(flagField->ySize());
      const cell_idx_t zs = cell_idx_c(flagField->zSize());
      const bool inner = (cx >= 0 && cx < xs && cy >= 0 && cy < ys && cz >= 0 && cz < zs);

      // Unknown-neighbour tally (independent of the gas/liquid/interface branch chosen below).
      if (inner)
      {
         ++formerMoTotal_;
         unknownNbrSum_ += unknownNbr;
         if (unknownNbr > uint_t(0)) ++formerMoWithUnknownNbr_;
         if (unknownNbr > unknownNbrMax_) unknownNbrMax_ = unknownNbr;
         ++unknownNbrHist_[unknownNbr]; // unknownNbr ∈ [0,26]
      }

      switch (state)
      {
      case BognerNeighbourhoodState::Gas:
         if (inner) ++convGas_;
         return bognerGas;
      case BognerNeighbourhoodState::Liquid:
         // Provably real liquid nearby, no real gas, no unresolved neighbour -> genuinely bulk
         // liquid. (Diagnostic under-resolved-waterline check removed with the old interface-only
         // branch B: promoting to liquid here is now the intended, invariant-safe behaviour, not an
         // over-promotion — the real gas/interface flags right next to this cell are what enforce
         // the liquid-never-touches-gas invariant, not the choice made here.)
         if (inner) ++convLiquid_;
         return bognerLiquid;
      case BognerNeighbourhoodState::Interface:
      {
         if (hasGas && hasLiquid)
         {
            // Genuine mixed/free-surface neighbourhood -- always interface, independent of the
            // waterline-reachability check below (there is nothing to resolve; this is already a
            // correct, real surface cell).
            if (inner) ++convIntfSurface_;
            return bognerInterface;
         }

         // [WATERLINE-CHECK] the remaining three sub-cases (would-be-liquid/gas blocked by an
         // unresolved neighbour, and branch E with neither real phase reachable within 1 cell) all
         // go through the phase-connectivity reachability check (PhaseReachability /
         // computeRealPhaseComponents, above) before falling back to interface. If it can prove --
         // from REAL bognerGas/bognerLiquid flags only, via genuine connectivity, never a guess --
         // that the dangerous opposite phase is not reachable at all, it commits to liquid/gas
         // directly instead of the conservative interface fallback. Validated as a log-only shadow
         // check first (see git history immediately before this commit) before being wired in here.
         const bool gasReach = waterlineShadow.gasConfidentlyReachable(cx, cy, cz);
         const bool gasIncon = waterlineShadow.gasInconclusive(cx, cy, cz);
         const bool liqReach = waterlineShadow.liquidConfidentlyReachable(cx, cy, cz);
         const bool liqIncon = waterlineShadow.liquidInconclusive(cx, cy, cz);
         const WaterlineVerdict verdict = classifyViaWaterlineReachability(gasReach, gasIncon, liqReach, liqIncon);

         if (verdict == WaterlineVerdict::WouldBeLiquid || verdict == WaterlineVerdict::WouldBeGas)
         {
            const bool toLiquid = (verdict == WaterlineVerdict::WouldBeLiquid);
            if (inner)
            {
               if (toLiquid) { ++waterlineResolvedLiquid_; }
               else { ++waterlineResolvedGas_; }
            }
            static uint_t waterlineResolvedLogCount = uint_t(0);
            ++waterlineResolvedLogCount;
            if (waterlineResolvedLogCount <= uint_t(20) || (waterlineResolvedLogCount % uint_t(5000) == uint_t(0)))
            {
               const auto&  aabb = block->getAABB();
               const real_t cdx  = aabb.xSize() / real_c(flagField->xSize());
               WALBERLA_LOG_WARNING("[WATERLINE-RESOLVED] (#" << waterlineResolvedLogCount << ") cell ("
                  << cx << ", " << cy << ", " << cz << ") world≈("
                  << (aabb.xMin() + (real_c(cx) + real_t(0.5)) * cdx) << ","
                  << (aabb.yMin() + (real_c(cy) + real_t(0.5)) * cdx) << ","
                  << (aabb.zMin() + (real_c(cz) + real_t(0.5)) * cdx) << ") had "
                  << (hasUnknown ? unknownNbr : uint_t(0)) << " unclassified FormerMO neighbour(s) and/or no real "
                     "phase within 1 cell, but the phase-connectivity check confirms no real "
                  << (toLiquid ? "gas" : "liquid") << " is reachable -- resolving to "
                  << (toLiquid ? "liquid" : "gas") << " instead of falling back to interface.");
            }
            return toLiquid ? bognerLiquid : bognerGas;
         }

         // Not resolved by the waterline check -- fall back to interface exactly as before, with
         // the existing branch-specific diagnostics.
         if (hasLiquid && !hasGas && hasUnknown)
         {
            // Would qualify as bulk liquid, but an unclassified FormerMO neighbour could still
            // become gas this pass — downgrade to the safe interface buffer instead of liquid.
            WALBERLA_LOG_WARNING("Bogner state selection: cell ("
                                 << cx << ", " << cy << ", " << cz
                                 << ") has " << unknownNbr << " unclassified FormerMO neighbour(s); falling back to "
                                    "interface instead of liquid to preserve the liquid/gas separation invariant.");
            if (inner) ++convIntfUnknown_;
         }
         else if (hasGas && !hasLiquid && hasUnknown)
         {
            // Symmetric case: would qualify as bulk gas, but an unclassified FormerMO neighbour could
            // still become liquid this pass — downgrade to the safe interface buffer instead of gas.
            WALBERLA_LOG_WARNING("Bogner state selection: cell ("
                                 << cx << ", " << cy << ", " << cz
                                 << ") has " << unknownNbr << " unclassified FormerMO neighbour(s); falling back to "
                                    "interface instead of gas to preserve the liquid/gas separation invariant.");
            if (inner) ++convIntfUnknownGas_;
         }
         else
         {
            // Neither a real gas nor a real liquid neighbour within one cell (branch E).
            if (inner) ++convIntfUnresolved_;

            // DIAGNOSTIC (never changes behaviour, throttled like the FormerMO-neighbour warning
            // above): look a few cells further up than the 1-cell stencil this classifier is
            // limited to. If there is real LIQUID directly overhead with no GAS in between, this
            // cell is sitting under a genuine water column -- i.e. very likely actually fully
            // submerged -- and only landed in this ambiguous fallback because its immediate
            // neighbourhood happens to be entirely interface/unresolved cells. That is exactly the
            // interface-contagion pattern the hasGas/hasLiquid fix above targets (see the class
            // comment): confirms at a glance whether THIS instance of the fallback is one of those
            // "actually submerged, just locally under-resolved" cells, e.g. deep inside a densely
            // packed, actively-jittering particle bed (cylinder_sedimentation_showcase_gas.prm).
            // (The waterline check above is the more thorough, connectivity-aware version of this
            // same idea; this narrower column probe still runs for whatever it doesn't resolve.)
            constexpr cell_idx_t columnCheckDepth = cell_idx_t(8);
            const cell_idx_t     gl = cell_idx_c(flagField->nrOfGhostLayers());
            bool columnIsLiquid = false, columnHitsGas = false;
            for (cell_idx_t k = 1; k <= columnCheckDepth; ++k)
            {
               const cell_idx_t nz = cz + k;
               if (nz >= zs + gl) break;
               const flag_t nf = flagField->get(cx, cy, nz);
               if (isFlagSet(nf, bognerGas)) { columnHitsGas = true; break; }
               if (isFlagSet(nf, bognerLiquid)) { columnIsLiquid = true; break; }
            }
            if (columnIsLiquid && !columnHitsGas)
            {
               static uint_t fallbackSubmergedLogCount = uint_t(0);
               ++fallbackSubmergedLogCount;
               if (fallbackSubmergedLogCount <= uint_t(20) || (fallbackSubmergedLogCount % uint_t(5000) == uint_t(0)))
               {
                  WALBERLA_LOG_WARNING("[FORMERMO-FALLBACK-SUBMERGED] (#" << fallbackSubmergedLogCount
                     << ") cell (" << cx << ", " << cy << ", " << cz << ") reconstructed as INTERFACE via the "
                        "ambiguous fallback (no real gas/liquid within 1 cell), but a real LIQUID column was "
                        "found directly overhead within " << columnCheckDepth << " cells with no GAS in "
                        "between -- this cell is very likely genuinely fully submerged, not near any real "
                        "free surface; only its immediate neighbourhood is itself interface/unresolved.");
               }
            }
         }

         if (inner)
         {
            switch (verdict)
            {
            case WaterlineVerdict::StillAmbiguousSurfaceNearby: ++shadowStillAmbiguousSurface_; break;
            case WaterlineVerdict::StillAmbiguousNoEvidence: ++shadowStillAmbiguousNoEvidence_; break;
            case WaterlineVerdict::ExcludedNearUnresolvedEdge: ++shadowExcludedNearEdge_; break;
            default: break; // WouldBeLiquid/WouldBeGas already returned above
            }
         }

         return bognerInterface;
      }
      }
      WALBERLA_ABORT("unreachable: classifyFormerObstacleNeighbourhood returned an unhandled state");
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
               reseedHydrostatic(x, y, z);
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
            reseedHydrostatic(x, y, z);

            // [RECON-DENSITY-CHECK] VERIFY the fill-blind over-density hypothesis with real numbers.
            // EquilibriumReconstructor just set this freed cell to ρ = unweighted mean of its domain
            // (liquid+interface) neighbours' getDensity — and interface cells report ρ≈1 no matter how
            // EMPTY they are (their partial volume lives in the fill field φ, which the reconstructor
            // ignores). Here we recompute that unweighted ρ AND the fill-weighted ρ = Σφρ/Σφ that the
            // proposed fix would use. If the cell sits at the waterline (low-φ interface neighbours) the
            // two diverge sharply → the assigned ρ≈1 is the phantom-force source the MEM then reads.
            if constexpr (!std::is_void_v< ScalarField_T >)
            if (fillField != nullptr)
            {
               const real_t rhoAssigned = pdfField->getDensity(x, y, z);
               real_t sumRho     = real_t(0);
               real_t sumPhiRho  = real_t(0);
               real_t sumPhi     = real_t(0);
               real_t minIntfPhi = real_t(1);
               uint_t nDom  = uint_t(0);
               uint_t nIntf = uint_t(0);
               uint_t nLiq  = uint_t(0);
               for (int dz = -1; dz <= 1; ++dz)
                  for (int dy = -1; dy <= 1; ++dy)
                     for (int dx = -1; dx <= 1; ++dx)
                     {
                        if (dx == 0 && dy == 0 && dz == 0) continue;
                        const cell_idx_t nx = x + cell_idx_c(dx);
                        const cell_idx_t ny = y + cell_idx_c(dy);
                        const cell_idx_t nz = z + cell_idx_c(dz);
                        if (!boundaryHandling->isDomain(Cell(nx, ny, nz))) continue;
                        const real_t rn = pdfField->getDensity(nx, ny, nz);
                        const real_t pn = fillField->get(nx, ny, nz);
                        sumRho += rn;
                        sumPhiRho += pn * rn;
                        sumPhi += pn;
                        ++nDom;
                        if (isFlagSet(flagField->get(nx, ny, nz), bognerInterface)) { ++nIntf; if (pn < minIntfPhi) minIntfPhi = pn; }
                        if (isFlagSet(flagField->get(nx, ny, nz), bognerLiquid)) ++nLiq;
                     }
               const real_t rhoUnw = (nDom > uint_t(0)) ? sumRho / real_c(nDom) : real_t(1);
               const real_t rhoW   = (sumPhi > real_t(1e-6)) ? sumPhiRho / sumPhi : real_t(0);

               // Smoking gun: high density was assigned, but the fill-weighted value is much lower AND
               // the neighbourhood is dominated by partially-empty interface cells (i.e. the waterline).
               if (rhoAssigned > real_t(0.8) && rhoW < real_t(0.6) && nIntf > nLiq)
               {
                  static uint_t reconDensLogCount = uint_t(0);
                  ++reconDensLogCount;
                  if (reconDensLogCount <= uint_t(20) || (reconDensLogCount % uint_t(5000) == uint_t(0)))
                  {
                     const auto&  aabb = block->getAABB();
                     const real_t cdx  = aabb.xSize() / real_c(flagField->xSize());
                     WALBERLA_LOG_WARNING("[RECON-DENSITY-CHECK] (#" << reconDensLogCount << ") world≈("
                        << (aabb.xMin() + (real_c(x) + real_t(0.5)) * cdx) << ","
                        << (aabb.yMin() + (real_c(y) + real_t(0.5)) * cdx) << ","
                        << (aabb.zMin() + (real_c(z) + real_t(0.5)) * cdx) << ")"
                        << "  rhoAssigned=" << rhoAssigned << "  rhoUnweighted=" << rhoUnw
                        << "  rhoFillWeighted=" << rhoW << "  domainNbrs=" << nDom
                        << "  interfaceNbrs=" << nIntf << " (minPhi=" << minIntfPhi << ")"
                        << "  liquidNbrs=" << nLiq
                        << "  --> FILL-BLIND OVER-DENSITY: reconstructor assigned ρ≈1 at the waterline;"
                           " the fill-weighted ρ would be " << rhoW << " (the phantom-force source).");
                  }
               }
            }

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
                                      const bool optimizeForSmallObstacleFraction = false,
                                      const bool useHydrostaticReconstruction = false,
                                      const real_t hydrostaticAccelerationZ = real_t(0),
                                      const real_t absoluteLiquidHeight = real_t(0))
{
   using RM_T = PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                          EquilibriumReconstructor< BoundaryHandling_T >, ScalarField_T >;
   return make_shared< RM_T >(RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac, formerObstacle,
                                   fluid,
                                   makeEquilibriumReconstructor< BoundaryHandling_T >(blockStorage, boundaryHandlingID),
                                   conserveMomentumWhenRestoring, liquidFlagUID, gasFlagUID, interfaceFlagUID,
                                   fillFieldID, useBognerStateSelection, optimizeForSmallObstacleFraction,
                                   useHydrostaticReconstruction, hydrostaticAccelerationZ, absoluteLiquidHeight));
}

} // namespace lbm_mesapd_coupling
} // namespace walberla