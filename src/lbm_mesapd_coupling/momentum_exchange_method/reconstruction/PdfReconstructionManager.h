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
 * `obstacle` excludes a neighbour from ALL evidence (hasGas/hasLiquid/hasUnknown) the moment it's
 * set, checked before any other flag on that neighbour -- a cell currently covered by a particle is
 * definitively solid right now, unlike `formerObstacle` (genuinely mid-reclassification, still
 * counted as evidence of ambiguity) or real liquid/gas (genuine evidence). This matters because
 * MovingParticleMapping's obstacle-covering step only clears/sets the `domain_` mask (liquid |
 * interface for Antidunes) -- a cell that was GAS when a particle's geometry swept over it is never
 * routed through that step at all (gas isn't part of `domain_`), so it can end up simultaneously
 * `obstacle`-flagged (from some other mapping path) or simply left mid-transition while still
 * carrying its old gas bit. Without this exclusion, such a stale bit would be read as genuine gas
 * evidence by any FormerMO neighbour scanning it. Confirmed via direct VTK inspection of a real run:
 * mapped_flag=GAS cells sitting >=1 LU inside 42-54 of 81 spheres, persisting across many timesteps.
 *
 * hasGas/hasLiquid/hasUnknown/unknownNbr are written back to the caller (needed for the per-branch
 * conversion counters and the same-pass-race warning, which need the raw booleans, not just the
 * final Gas/Liquid/Interface classification).
 */
//**************************************************************************************************************************************
template< typename Flag_T, typename NeighbourFlagGetter >
inline BognerNeighbourhoodState classifyFormerObstacleCell(NeighbourFlagGetter getNeighbourFlag,
                                                             Flag_T bognerLiquid, Flag_T bognerGas,
                                                             Flag_T formerObstacle, Flag_T obstacle,
                                                             bool& hasGasOut,
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
            if (isFlagSet(nf, obstacle)) continue; // currently solid -- not evidence, not ambiguity, pure exclusion
            if (isFlagSet(nf, bognerLiquid))   hasLiquidOut = true;
            if (isFlagSet(nf, bognerGas))      hasGasOut    = true;
            if (isFlagSet(nf, formerObstacle)) { hasUnknownOut = true; ++unknownNbrOut; }
         }
   return classifyFormerObstacleNeighbourhood(hasGasOut, hasLiquidOut, hasUnknownOut);
}

//**************************************************************************************************************************************
/*!\brief Pass-3 candidacy rule: should this FormerMO cell's fallback-to-interface outcome from
 * Pass 2 be re-examined once every same-pass FormerMO cell in the block has a final state?
 *
 * True only for the two same-pass-ambiguity fallback branches (A': hasGas&&!hasLiquid&&hasUnknown,
 * and C: hasLiquid&&!hasGas&&hasUnknown) plus branch E (!hasGas&&!hasLiquid) when it ALSO had an
 * unresolved neighbour (hasUnknown) -- a branch-E cell with zero unresolved neighbours at Pass-2
 * time cannot have anything change in its neighbourhood by Pass-3 time, so it is structurally
 * unrescuable and not worth queuing. Never true for branch D (hasGas&&hasLiquid): that is already a
 * genuine, correct free-surface cell, and re-scanning the single largest bucket for zero benefit
 * would be pure waste.
 */
//**************************************************************************************************************************************
inline bool isPass3Candidate(BognerNeighbourhoodState state, bool hasGas, bool hasLiquid, bool hasUnknown)
{
   if (state != BognerNeighbourhoodState::Interface) return false;
   if (hasGas && hasLiquid) return false; // branch D -- genuine surface, nothing to resolve
   return hasUnknown;
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
        liquidFlagUID_(fluid), gasFlagUID_(fluid), interfaceFlagUID_(fluid), movingObstacleFlagUID_(fluid),
        fillFieldID_(BlockDataID())
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
                                   const FlagUID& movingObstacleFlagUID,
                                   const BlockDataID& fillFieldID, const bool useBognerStateSelection = true,
                                   const bool optimizeForSmallObstacleFraction = false,
                                   const bool useHydrostaticReconstruction = false,
                                   const real_t hydrostaticAccelerationZ = real_t(0),
                                   const real_t absoluteLiquidHeight = real_t(0),
                                   const bool usePass3Reresolution = false)
      : blockStorage_(blockStorage), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID),
        particleFieldID_(particleFieldID), ac_(ac), reconstructor_(reconstructor), formerObstacle_(formerObstacle),
        fluid_(fluid), conserveMomentumWhenRestoring_(conserveMomentumWhenRestoring),
        optimizeForSmallObstacleFraction_(optimizeForSmallObstacleFraction),
        useBognerStateSelection_(useBognerStateSelection),
        useHydrostaticReconstruction_(useHydrostaticReconstruction),
        hydrostaticAccelerationZ_(hydrostaticAccelerationZ), absoluteLiquidHeight_(absoluteLiquidHeight),
        liquidFlagUID_(liquidFlagUID), gasFlagUID_(gasFlagUID), interfaceFlagUID_(interfaceFlagUID),
        movingObstacleFlagUID_(movingObstacleFlagUID),
        fillFieldID_(fillFieldID), usePass3Reresolution_(usePass3Reresolution)
   {}

   void inline operator()(IBlock* const block);

   // Diagnostic: MPI-reduce (to root), log, and reset the FormerMO→state conversion counters
   // accumulated since the last call. One line per call → invoke periodically (e.g. once per
   // evaluationFrequency steps via an AfterFunction). Counts INNER cells only, so the totals are
   // the true number of physical cells converted, not duplicated across block ghost layers.
   void logAndResetConversionCounts(const uint_t t)
   {
      std::vector< uint_t > c = { convGas_, convLiquid_, convIntfUnknown_, convIntfUnknownGas_, convIntfSurface_, convIntfUnresolved_,
                                  convPass3ToLiquid_, convPass3ToGas_, convPass3StillInterface_, convPass3UnresolvedIsland_ };
      walberla::mpi::reduceInplace(c, walberla::mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-CONVERT t=" << t << "] cells converted since last report: "
         << "gas(A)=" << c[0] << " liquid(B)=" << c[1]
         << " interface-unknown(C)=" << c[2] << " interface-unknown-gas(A')=" << c[3] << " interface-surface(D)=" << c[4]
         << " interface-unresolved(E)=" << c[5]
         << " | total=" << (c[0] + c[1] + c[2] + c[3] + c[4] + c[5])
         << " | pass3ToLiquid=" << c[6] << " pass3ToGas=" << c[7] << " pass3StillInterface=" << c[8]
         << " pass3UnresolvedIsland=" << c[9]);

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
      convPass3ToLiquid_ = convPass3ToGas_ = convPass3StillInterface_ = convPass3UnresolvedIsland_ = uint_t(0);
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
         // Bounds-checked D3Q27 read -- same fix as getBoundedNeighbourFlag above (commit d8c0a265
         // fixed the analogous unguarded read in bognerTargetFlag/classifyFormerObstacleCell; this
         // function has its OWN, separate unguarded 3x3x3 read that fix never covered). A cell
         // converted at the OUTERMOST ghost layer has a D3Q27 neighbour one cell beyond flagField's/
         // fillField's allocated (ghost-inclusive) extent. Confirmed with AddressSanitizer:
         // heap-buffer-overflow READ at this exact loop (job 4063969, 1000-sphere v6 scenario).
         const auto ghostBounds = flagField->xyzSizeWithGhostLayer();
         real_t sum = real_t(0);
         uint_t cnt = uint_t(0);
         for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
               for (int dx = -1; dx <= 1; ++dx)
               {
                  if (dx == 0 && dy == 0 && dz == 0) continue;
                  const cell_idx_t nx = x + cell_idx_c(dx);
                  const cell_idx_t ny = y + cell_idx_c(dy);
                  const cell_idx_t nz = z + cell_idx_c(dz);
                  if (nx < ghostBounds.xMin() || nx > ghostBounds.xMax() ||
                      ny < ghostBounds.yMin() || ny > ghostBounds.yMax() ||
                      nz < ghostBounds.zMin() || nz > ghostBounds.zMax())
                     continue; // neighbour outside the allocated field -- no evidence, skip (same
                               // conservative stance as getBoundedNeighbourFlag's out-of-bounds case)
                  if (isFlagSet(flagField->get(nx, ny, nz), interfaceFlag))
                  {
                     sum += fillField->get(nx, ny, nz);
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
   const FlagUID movingObstacleFlagUID_; // currently-covered flag, excluded entirely from Bogner classifier evidence
   const BlockDataID fillFieldID_; // fill-level field for φ seeding (Bogner path only)

   // Pass-3 local re-resolution (prm: ModelParameters.usePass3Reresolution): after Pass 2 finishes and
   // every FormerMO cell has a final state, re-examine cells that fell back to interface purely because
   // of same-pass ambiguity (see isPass3Candidate) and resolve them directly. Off by default -> identical
   // to the plain Pass-2-only fallback behaviour.
   const bool usePass3Reresolution_ = false;

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

   // Pass-3 conversion counters (inner cells only, matching the Pass-2 convention above). A Pass-3
   // decision is one of: promote to liquid/gas, or stay interface -- convPass3StillInterface_ counts
   // every "stay interface" outcome, convPass3UnresolvedIsland_ is the subset of those that are
   // genuinely !hasGas && !hasLiquid even after full same-pass resolution (see isPass3Candidate's
   // doc and the [PASS3-UNRESOLVED-ISLAND] log below) -- a strict subset, not an alternative bucket.
   uint_t convPass3ToLiquid_        = uint_t(0);
   uint_t convPass3ToGas_           = uint_t(0);
   uint_t convPass3StillInterface_  = uint_t(0);
   uint_t convPass3UnresolvedIsland_ = uint_t(0);
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
   flag_t movingObstacle  = flag_t(0);

   if (useBognerStateSelection_)
   {
      WALBERLA_ASSERT(flagField->flagExists(liquidFlagUID_));
      WALBERLA_ASSERT(flagField->flagExists(gasFlagUID_));
      WALBERLA_ASSERT(flagField->flagExists(interfaceFlagUID_));
      WALBERLA_ASSERT(flagField->flagExists(movingObstacleFlagUID_));

      bognerLiquid    = flagField->getFlag(liquidFlagUID_);
      bognerGas       = flagField->getFlag(gasFlagUID_);
      bognerInterface = flagField->getFlag(interfaceFlagUID_);
      movingObstacle  = flagField->getFlag(movingObstacleFlagUID_);
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

   // Bounds-checked D3Q27 neighbour read, shared by bognerTargetFlag (Pass 2) and the Pass-3 decide
   // loop below -- factored out so there is only one copy of the logic whose original hand-inlined
   // version caused the ASAN-fixed OOB read (commit d8c0a265). A neighbour falling outside the
   // field's allocated (ghost-inclusive) range is reported as unresolved ("formerObstacle") -- the
   // same conservative fallback used for a genuine same-pass-unresolved neighbour -- rather than
   // reading invalid memory.
   auto getBoundedNeighbourFlag = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, int dx, int dy,
                                      int dz) -> flag_t {
      const cell_idx_t nx = cx + cell_idx_c(dx);
      const cell_idx_t ny = cy + cell_idx_c(dy);
      const cell_idx_t nz = cz + cell_idx_c(dz);
      if (nx < formerObstacleGhostBounds.xMin() || nx > formerObstacleGhostBounds.xMax() ||
          ny < formerObstacleGhostBounds.yMin() || ny > formerObstacleGhostBounds.yMax() ||
          nz < formerObstacleGhostBounds.zMin() || nz > formerObstacleGhostBounds.zMax())
         return formerObstacle;
      return flagField->get(nx, ny, nz);
   };

   // Pass-3 candidates collected during Pass 2 below (see isPass3Candidate) -- a plain local vector
   // is sufficient since Pass 2 and Pass 3 run back-to-back within this single operator() call, no
   // persistent BlockDataID/flag bit needed. Always collected (cheap push_backs); only the Pass-3
   // processing loop further down is gated by usePass3Reresolution_.
   struct Pass3Candidate { Cell cell; bool inner; };
   std::vector< Pass3Candidate > pass3Candidates;

   auto bognerTargetFlag = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> flag_t {
      // Scan + decision both live in classifyFormerObstacleCell()/classifyFormerObstacleNeighbourhood()
      // (top of this file) as pure, directly unit-tested functions
      // (tests/showcases_antidunes/BognerStateClassificationTest.cpp) — this lambda only supplies
      // the real field access and attaches the side effects (per-branch counters, the same-pass-race
      // warning) that need cx/cy/cz/inner and can't live in the pure functions themselves.
      bool hasLiquid, hasGas, hasUnknown;
      uint_t unknownNbr;
      const BognerNeighbourhoodState state = classifyFormerObstacleCell(
         [&](int dx, int dy, int dz) -> flag_t { return getBoundedNeighbourFlag(cx, cy, cz, dx, dy, dz); },
         bognerLiquid, bognerGas, formerObstacle, movingObstacle, hasGas, hasLiquid, hasUnknown, unknownNbr);

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
            // Genuine mixed/free-surface neighbourhood -- always interface, nothing to resolve.
            if (inner) ++convIntfSurface_;
            return bognerInterface;
         }

         // The remaining sub-cases (would-be-liquid/gas blocked by an unresolved same-pass
         // neighbour, and branch E with neither real phase within 1 cell) fall back to interface
         // here, exactly as Pass 2 always has -- Pass 3 (usePass3Reresolution_, below, after every
         // FormerMO cell in the block has a final state) is what actually resolves most of these,
         // not this fallback itself.
         if (hasLiquid && !hasGas && hasUnknown)
         {
            if (inner) ++convIntfUnknown_;
            static uint_t provisionalLogCount = uint_t(0);
            ++provisionalLogCount;
            if (provisionalLogCount <= uint_t(20) || (provisionalLogCount % uint_t(5000) == uint_t(0)))
            {
               WALBERLA_LOG_INFO("Bogner state selection: cell (" << cx << ", " << cy << ", " << cz << ") has "
                  << unknownNbr << " unclassified FormerMO neighbour(s); "
                  << (usePass3Reresolution_
                         ? "provisionally interface pending Pass-3 re-resolution."
                         : "falling back to interface instead of liquid to preserve the liquid/gas separation "
                           "invariant."));
            }
         }
         else if (hasGas && !hasLiquid && hasUnknown)
         {
            // Symmetric case.
            if (inner) ++convIntfUnknownGas_;
            static uint_t provisionalGasLogCount = uint_t(0);
            ++provisionalGasLogCount;
            if (provisionalGasLogCount <= uint_t(20) || (provisionalGasLogCount % uint_t(5000) == uint_t(0)))
            {
               WALBERLA_LOG_INFO("Bogner state selection: cell (" << cx << ", " << cy << ", " << cz << ") has "
                  << unknownNbr << " unclassified FormerMO neighbour(s); "
                  << (usePass3Reresolution_
                         ? "provisionally interface pending Pass-3 re-resolution."
                         : "falling back to interface instead of gas to preserve the liquid/gas separation "
                           "invariant."));
            }
         }
         else
         {
            // Neither a real gas nor a real liquid neighbour within one cell (branch E).
            if (inner) ++convIntfUnresolved_;
         }

         if (isPass3Candidate(state, hasGas, hasLiquid, hasUnknown))
            pass3Candidates.push_back({ Cell(cx, cy, cz), inner });

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

   // ---- Pass 3: local re-resolution (usePass3Reresolution_ only) ------------------------------
   //
   // By the time Pass 2 above has finished, every FormerMO cell in this block has a final,
   // definite state -- there is no more "unknown" anywhere. pass3Candidates (collected during
   // bognerTargetFlag above) holds every cell that fell back to interface purely because of
   // same-pass ambiguity (isPass3Candidate). Re-running the same classifier on them now, with
   // hasUnknown guaranteed false (except the literal outermost ghost layer, which degrades safely
   // to "stays interface" via getBoundedNeighbourFlag), resolves most of them directly.
   //
   // MUST be decide-then-apply (two sub-loops), not decide-and-apply-immediately: two mutual-
   // neighbour candidates must both be classified against the SAME frozen post-Pass-2 snapshot,
   // otherwise the one processed second would see the first's fresh promotion and could cascade
   // into an incorrect classification -- the exact same class of same-pass-order bug Pass 3 exists
   // to fix in Pass 2, reappearing one level in if this isn't respected.
   if (usePass3Reresolution_ && !pass3Candidates.empty())
   {
      struct Pass3Decision { Cell cell; bool inner; BognerNeighbourhoodState state; bool hasGas; bool hasLiquid; };
      std::vector< Pass3Decision > pass3Decisions;
      pass3Decisions.reserve(pass3Candidates.size());

      // Sub-loop 1: decide (read-only against the frozen post-Pass-2 snapshot).
      for (const auto& candidate : pass3Candidates)
      {
         const cell_idx_t cx = candidate.cell.x(), cy = candidate.cell.y(), cz = candidate.cell.z();
         bool hasG, hasL, hasU;
         uint_t nbr;
         const BognerNeighbourhoodState state = classifyFormerObstacleCell(
            [&](int dx, int dy, int dz) -> flag_t { return getBoundedNeighbourFlag(cx, cy, cz, dx, dy, dz); },
            bognerLiquid, bognerGas, formerObstacle, movingObstacle, hasG, hasL, hasU, nbr);
         pass3Decisions.push_back({ candidate.cell, candidate.inner, state, hasG, hasL });

         if (state == BognerNeighbourhoodState::Interface)
         {
            if (candidate.inner) ++convPass3StillInterface_;
            if (!hasG && !hasL)
            {
               // Genuinely unresolvable: even after every same-pass FormerMO cell in the block has
               // a final state, this cell still has no real liquid and no real gas within 1 cell in
               // any direction -- unlike an ordinary Pass-2 branch-E cell (which might simply not
               // have had time to resolve yet), this is an anomalous signal worth surfacing
               // distinctly rather than folding into the generic still-interface count. It must stay
               // interface either way -- there is no safe basis to pick liquid or gas here.
               if (candidate.inner) ++convPass3UnresolvedIsland_;
               static uint_t unresolvedIslandLogCount = uint_t(0);
               ++unresolvedIslandLogCount;
               if (unresolvedIslandLogCount <= uint_t(20) || (unresolvedIslandLogCount % uint_t(5000) == uint_t(0)))
               {
                  const auto&  aabb = block->getAABB();
                  const real_t cdx  = aabb.xSize() / real_c(flagField->xSize());
                  WALBERLA_LOG_WARNING("[PASS3-UNRESOLVED-ISLAND] (#" << unresolvedIslandLogCount << ") cell ("
                     << cx << ", " << cy << ", " << cz << ") world≈("
                     << (aabb.xMin() + (real_c(cx) + real_t(0.5)) * cdx) << ","
                     << (aabb.yMin() + (real_c(cy) + real_t(0.5)) * cdx) << ","
                     << (aabb.zMin() + (real_c(cz) + real_t(0.5)) * cdx) << ") still has no real liquid and no "
                        "real gas within 1 cell even after every same-pass FormerMO neighbour resolved -- "
                        "genuinely isolated from both phases, staying interface (no safe basis to pick "
                        "either side).");
               }
            }
         }
      }

      // Sub-loop 2: apply the buffered decisions.
      for (const auto& decision : pass3Decisions)
      {
         if (decision.state != BognerNeighbourhoodState::Gas && decision.state != BognerNeighbourhoodState::Liquid)
            continue; // stays interface -- already correct from Pass 2, already tallied above

         const flag_t newFlag = (decision.state == BognerNeighbourhoodState::Gas) ? bognerGas : bognerLiquid;
         boundaryHandling->forceDomain(newFlag, decision.cell.x(), decision.cell.y(), decision.cell.z());
         if constexpr (!std::is_void_v< ScalarField_T >)
            if (fillField != nullptr)
               seedConvertedCellFill(fillField, flagField, decision.cell.x(), decision.cell.y(), decision.cell.z(),
                                     newFlag, bognerGas, bognerLiquid, bognerInterface);
         if (decision.inner)
         {
            if (decision.state == BognerNeighbourhoodState::Gas) { ++convPass3ToGas_; }
            else { ++convPass3ToLiquid_; }
         }
      }
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
                                      const FlagUID& movingObstacleFlagUID,
                                      const BlockDataID& fillFieldID, const bool useBognerStateSelection = true,
                                      const bool optimizeForSmallObstacleFraction = false,
                                      const bool useHydrostaticReconstruction = false,
                                      const real_t hydrostaticAccelerationZ = real_t(0),
                                      const real_t absoluteLiquidHeight = real_t(0),
                                      const bool usePass3Reresolution = false)
{
   using RM_T = PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                          EquilibriumReconstructor< BoundaryHandling_T >, ScalarField_T >;
   return make_shared< RM_T >(RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac, formerObstacle,
                                   fluid,
                                   makeEquilibriumReconstructor< BoundaryHandling_T >(blockStorage, boundaryHandlingID),
                                   conserveMomentumWhenRestoring, liquidFlagUID, gasFlagUID, interfaceFlagUID,
                                   movingObstacleFlagUID,
                                   fillFieldID, useBognerStateSelection, optimizeForSmallObstacleFraction,
                                   useHydrostaticReconstruction, hydrostaticAccelerationZ, absoluteLiquidHeight,
                                   usePass3Reresolution));
}

} // namespace lbm_mesapd_coupling
} // namespace walberla