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

#include "core/mpi/MPIManager.h"
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

#include <algorithm>
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
/*!\brief Default (non-strict) decision rule of the Bogner state selection (bognerTargetFlag).
 *
 * Inputs: whether the D3Q27 neighbourhood contains a real gas cell, a real liquid cell and/or a cell freed in the same
 * pass ("unknown", formerObstacle):
 *  - gas,    no liquid, no unknown -> Gas
 *  - liquid, no gas,    no unknown -> Liquid
 *  - gas or liquid only, + unknown -> Interface (the unknown neighbour may still turn into the other phase in this
 *                                     pass; the rule is symmetric for both phases)
 *  - gas and liquid                -> Interface (free surface)
 *  - neither                       -> Interface (only interface or unknown neighbours; safe fallback)
 *
 * hasGas/hasLiquid count real gas/liquid flags only. Interface cells are the legal buffer between the phases; counting
 * them as gas would spread interface cells, pass by pass, through the narrow gaps of a packed particle bed.
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
/*!\brief Scans a cell's D3Q27 neighbourhood and classifies it with classifyFormerObstacleNeighbourhood (the per-cell
 * logic of bognerTargetFlag, separated so that tests can pass a literal neighbourhood).
 *
 * getNeighbourFlag(dx, dy, dz) returns the flag at that offset for every (dx, dy, dz) in {-1,0,1}^3 except (0,0,0).
 * Flag_T only needs isFlagSet(value, testFlag), so the production flag_t and a test flag type both work.
 *
 * A neighbour flagged `obstacle` (currently covered by a particle) is checked first and excluded from all evidence, so
 * liquid/gas bits left on a covered cell are never counted. A `formerObstacle` neighbour counts as unknown.
 *
 * hasGas/hasLiquid/hasUnknown and unknownNbr (number of formerObstacle neighbours) are returned for the per-branch
 * conversion counters.
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
/*!\brief Strict Bogner rule (ModelParameters.strictBognerClassification), Bogner (2017) Sec. 5.2, pp. 113-114, and
 * Bogner & Rüde (2013) Sec. 2.4, applied literally and in this order:
 *  - B(x) contains no liquid -> Gas
 *  - B(x) contains no gas    -> Liquid
 *  - B(x) contains both      -> Interface
 * B(x) is the non-obstacle neighbourhood. A neighbour freed in the same step (formerObstacle) was an obstacle when the
 * step began, so it is not in B(x), whereas classifyFormerObstacleNeighbourhood treats it as unknown. Bogner gives no
 * same-step rule: he assumes the body is slow enough that a freed cell always borders fluid that already has a state.
 */
//**************************************************************************************************************************************
inline BognerNeighbourhoodState classifyStrictBognerNeighbourhood(bool hasGas, bool hasLiquid)
{
   if (!hasLiquid) return BognerNeighbourhoodState::Gas;
   if (!hasGas) return BognerNeighbourhoodState::Liquid;
   return BognerNeighbourhoodState::Interface;
}

//! Scan for classifyStrictBognerNeighbourhood: neighbours flagged obstacle or formerObstacle are not in
//! B(x); hasGas/hasLiquid count real gas/liquid flags only (an interface neighbour counts toward neither, as in
//! classifyFormerObstacleCell).
template< typename Flag_T, typename NeighbourFlagGetter >
inline BognerNeighbourhoodState classifyFormerObstacleCellStrict(NeighbourFlagGetter getNeighbourFlag, Flag_T bognerLiquid,
                                                                   Flag_T bognerGas, Flag_T formerObstacle, Flag_T obstacle,
                                                                   bool& hasGasOut, bool& hasLiquidOut)
{
   hasGasOut = false; hasLiquidOut = false;
   for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
         for (int dx = -1; dx <= 1; ++dx)
         {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const Flag_T nf = getNeighbourFlag(dx, dy, dz);
            if (isFlagSet(nf, obstacle) || isFlagSet(nf, formerObstacle)) continue; // not in B(x)
            if (isFlagSet(nf, bognerLiquid)) hasLiquidOut = true;
            if (isFlagSet(nf, bognerGas)) hasGasOut = true;
         }
   return classifyStrictBognerNeighbourhood(hasGasOut, hasLiquidOut);
}

//**************************************************************************************************************************************
/*!\brief Pass-3 candidacy: should a cell that Pass 2 set to interface be re-examined once every cell freed in the
 * block has a final state?
 *
 * True for the interface fallbacks with an unknown neighbour: A' (gas + unknown), C (liquid + unknown) and E (neither
 * phase) if it has an unknown neighbour. A branch-E cell without one cannot change by Pass 3, and branch D (gas and
 * liquid) is a genuine surface cell.
 */
//**************************************************************************************************************************************
inline bool isPass3Candidate(BognerNeighbourhoodState state, bool hasGas, bool hasLiquid, bool hasUnknown)
{
   if (state != BognerNeighbourhoodState::Interface) return false;
   if (hasGas && hasLiquid) return false; // branch D -- genuine surface, nothing to resolve
   return hasUnknown;
}

//**************************************************************************************************************************************
/*!\brief Fill level (volume fraction) seeded into a freed cell that became interface.
 *
 * getNeighbourFlag(dx, dy, dz) / getNeighbourFill(dx, dy, dz) return flag and fill of the D3Q27 neighbour at that offset
 * (the centre is never queried); a neighbour outside the field should be reported with flag 0 (no evidence).
 *
 * allResolvedNeighbours == false (default; Bogner 2017, Sec. 5.2): mean fill of the interface neighbours.
 * allResolvedNeighbours == true (ModelParameters.seedFillFromAllNeighbours): mean volume fraction of every resolved
 *   neighbour (gas 0, liquid 1, interface its fill). The default ignores gas and liquid evidence: a cell that is
 *   interface only because of an unknown neighbour (branch C or A') usually has no interface neighbour and is seeded
 *   0.5, i.e. half empty inside the water (C) or half full in the air (A'). With all evidence C seeds about 1, A' about
 *   0 and a surface cell (D) the local liquid fraction. Only the fill changes, never the cell state.
 * In both rules the result is clamped to [0,1], and 0.5 is used if no neighbour contributes.
 *
 * skipObstacleNeighbours (ModelParameters.seedSkipObstacleNeighbours, default true) excludes neighbours flagged obstacle
 *   or formerObstacle in both rules. false: the all-neighbour rule still excludes them, but the interface-only rule
 *   counts an interface-flagged neighbour even if it also carries one of these flags.
 *
 * fallbackOut is set when no evidence was found and 0.5 was used.
 */
//**************************************************************************************************************************************
template< typename Flag_T, typename NeighbourFlagGetter, typename NeighbourFillGetter >
inline real_t seedInterfaceFillFromNeighbourhood(NeighbourFlagGetter getNeighbourFlag, NeighbourFillGetter getNeighbourFill,
                                                 Flag_T liquid, Flag_T gas, Flag_T interface, Flag_T formerObstacle,
                                                 Flag_T obstacle, bool allResolvedNeighbours, bool skipObstacleNeighbours,
                                                 bool& fallbackOut)
{
   real_t sum = real_t(0);
   uint_t cnt = uint_t(0);
   for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
         for (int dx = -1; dx <= 1; ++dx)
         {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const Flag_T nf = getNeighbourFlag(dx, dy, dz);
            if ((allResolvedNeighbours || skipObstacleNeighbours) &&
                (isFlagSet(nf, obstacle) || isFlagSet(nf, formerObstacle)))
               continue;
            if (isFlagSet(nf, interface))
            {
               sum += getNeighbourFill(dx, dy, dz);
               ++cnt;
            }
            else if (allResolvedNeighbours && isFlagSet(nf, liquid))
            {
               sum += real_t(1);
               ++cnt;
            }
            else if (allResolvedNeighbours && isFlagSet(nf, gas))
               ++cnt; // contributes 0
         }
   fallbackOut = (cnt == uint_t(0));
   return fallbackOut ? real_t(0.5) : std::clamp(sum / real_c(cnt), real_t(0), real_t(1));
}

// ScalarField_T (defaulted to void) is the fill-level field type, used only by the Bogner
// state-selection path to seed the volume fraction φ of freshly converted cells
// (Bogner 2017, Sec. 5.2) and to weight the freed-cell momentum budget. For the non-Bogner path it
// stays void and all fill code is compiled out via `if constexpr`.
template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T,
          typename ScalarField_T = void >
class PdfReconstructionManager
{
   static_assert(std::is_base_of_v< mesa_pd::data::IAccessor, ParticleAccessor_T >,
                 "Provide a valid accessor as template");

 public:
   using FlagField_T = typename BoundaryHandling_T::FlagField;
   using flag_t      = typename BoundaryHandling_T::flag_t;

   // Plain constructor — all FormerMO cells are uniformly set to `fluid`.
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

   // Bogner state-selection constructor (Bogner 2017, Sec. 5.2, after Bogner 2009): each FormerMO cell is converted to
   // gas, liquid or interface from its non-obstacle neighbourhood B(x) (classifyFormerObstacleCell, or
   // classifyFormerObstacleCellStrict with setStrictBognerClassification), and its fill level is seeded in Pass 2 via
   // fillFieldID (gas 0, liquid 1, interface see seedInterfaceFillFromNeighbourhood).
   // useBognerStateSelection = false sets every FormerMO cell to `fluid`, as the plain constructor does.
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

   //! ModelParameters.seedFillFromAllNeighbours (see seedInterfaceFillFromNeighbourhood)
   void setSeedFillFromAllNeighbours(bool on) { seedFillFromAllNeighbours_ = on; }
   //! ModelParameters.seedSkipObstacleNeighbours (see seedInterfaceFillFromNeighbourhood)
   void setSeedSkipObstacleNeighbours(bool on) { seedSkipObstacleNeighbours_ = on; }
   //! ModelParameters.strictBognerClassification (see classifyStrictBognerNeighbourhood)
   void setStrictBognerClassification(bool on) { strictBognerClassification_ = on; }
   //! ModelParameters.strictBognerNoLiquidAtGas: with the strict rule, a freed cell decided liquid that touches a
   //! (non-obstacle) gas cell once all decisions of the step are applied becomes interface with fill 1 (see Pass 2)
   void setStrictBognerNoLiquidAtGas(bool on) { strictBognerNoLiquidAtGas_ = on; }
   //! ModelParameters.strictBognerRelaxInterface, see the member comment below
   void setStrictBognerRelaxInterface(bool on) { strictBognerRelaxInterface_ = on; }

   //! Momentum budget diagnostic: what the cells freed since the last call carry into the fluid, weighted
   //! like fluid mass (liquid 1, interface its fill, gas 0) after this step's classification -- {Px, Py, Pz, L_y about the
   //! freeing particle's centre, mass, cells}. Rank-local (caller reduces); resets.
   std::array< real_t, 6 > takeFreedMomentum()
   {
      const std::array< real_t, 6 > r = freedMomentum_;
      freedMomentum_.fill(real_t(0));
      return r;
   }

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

      // Unknown-neighbour statistics (see unknownNbrHist_)
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

      // Seeding bookkeeping (see seededInterfaceCells_)
      {
         std::vector< uint_t > sc = { seededInterfaceCells_, seededInterfaceFallback_ };
         real_t                sv = seededInterfaceFillSum_;
         walberla::mpi::reduceInplace(sc, walberla::mpi::SUM);
         walberla::mpi::reduceInplace(sv, walberla::mpi::SUM);
         WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-SEED t=" << t << "] interface cells seeded since last report: " << sc[0]
                                   << " (0.5 fallback, no evidence: " << sc[1] << ") | seeded liquid volume (sum of fill) "
                                   << sv << " cells | rule: "
                                   << (seedFillFromAllNeighbours_ ? "all resolved neighbours" : "interface neighbours"));
         seededInterfaceCells_ = seededInterfaceFallback_ = uint_t(0);
         seededInterfaceFillSum_ = real_t(0);
      }
      // Strict Bogner rule: cells that became gas with no gas in B(x) (only interface neighbours, or none), cells that
      // became liquid but touch gas after the step's decisions, and how many of those were made interface
      if (strictBognerClassification_)
      {
         std::vector< uint_t > sb = { strictGasWithoutGas_, strictLiquidTouchingGas_, strictLiquidToInterface_ };
         walberla::mpi::reduceInplace(sb, walberla::mpi::SUM);
         WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-STRICT t=" << t << "] strict Bogner rule since last report: gas with no gas in B(x)="
                                   << sb[0] << " | liquid touching gas after the step=" << sb[1]
                                   << " | of these made interface (strictBognerNoLiquidAtGas)=" << sb[2]);
         strictGasWithoutGas_ = strictLiquidTouchingGas_ = strictLiquidToInterface_ = uint_t(0);
         if (strictBognerRelaxInterface_)
         {
            std::vector< uint_t > sr = { strictRelaxedToGas_, strictRelaxedToLiquid_, strictRelaxedStillInterface_ };
            walberla::mpi::reduceInplace(sr, walberla::mpi::SUM);
            WALBERLA_LOG_INFO_ON_ROOT("[FORMERMO-RELAX t=" << t << "] strictBognerRelaxInterface since last report: relaxed to gas="
                                      << sr[0] << " relaxed to liquid=" << sr[1] << " still genuinely interface=" << sr[2]);
            strictRelaxedToGas_ = strictRelaxedToLiquid_ = strictRelaxedStillInterface_ = uint_t(0);
         }
      }
      convGas_ = convLiquid_ = convIntfUnknown_ = convIntfUnknownGas_ = convIntfSurface_ = convIntfUnresolved_ = uint_t(0);
      formerMoTotal_ = formerMoWithUnknownNbr_ = unknownNbrSum_ = unknownNbrMax_ = uint_t(0);
      unknownNbrHist_.fill(uint_t(0));
      convPass3ToLiquid_ = convPass3ToGas_ = convPass3StillInterface_ = convPass3UnresolvedIsland_ = uint_t(0);
   }

 private:
   // Seeds the fill level of a just-converted cell: gas 0, liquid 1, interface from its D3Q27 neighbourhood
   // (seedInterfaceFillFromNeighbourhood). Only instantiated for non-void ScalarField_T.
   template< typename ScalarFieldPtr_T >
   void seedConvertedCellFill(ScalarFieldPtr_T fillField, FlagField_T* flagField, cell_idx_t x, cell_idx_t y,
                              cell_idx_t z, flag_t targetFlag, flag_t gasFlag, flag_t liquidFlag,
                              flag_t interfaceFlag, const IBlock* block = nullptr) const
   {
      real_t seededFill;
      if (isFlagSet(targetFlag, gasFlag)) { seededFill = real_t(0); }
      else if (isFlagSet(targetFlag, liquidFlag)) { seededFill = real_t(1); }
      else // interface
      {
         // Bounds-checked neighbour read: a cell in the outermost ghost layer has neighbours outside the allocated
         // field; they are reported with flag 0 (no evidence), so their fill is never read.
         const auto   ghostBounds = flagField->xyzSizeWithGhostLayer();
         auto inRange = [&](int dx, int dy, int dz) {
            const cell_idx_t nx = x + cell_idx_c(dx), ny = y + cell_idx_c(dy), nz = z + cell_idx_c(dz);
            return !(nx < ghostBounds.xMin() || nx > ghostBounds.xMax() || ny < ghostBounds.yMin() ||
                     ny > ghostBounds.yMax() || nz < ghostBounds.zMin() || nz > ghostBounds.zMax());
         };
         const flag_t formerObstacleFlag =
            flagField->flagExists(formerObstacle_) ? flagField->getFlag(formerObstacle_) : flag_t(0);
         const flag_t obstacleFlag =
            flagField->flagExists(movingObstacleFlagUID_) ? flagField->getFlag(movingObstacleFlagUID_) : flag_t(0);
         bool fallback = false;
         seededFill    = seedInterfaceFillFromNeighbourhood< flag_t >(
            [&](int dx, int dy, int dz) {
               return inRange(dx, dy, dz) ? flagField->get(x + cell_idx_c(dx), y + cell_idx_c(dy), z + cell_idx_c(dz))
                                             : flag_t(0);
            },
            [&](int dx, int dy, int dz) {
               return fillField->get(x + cell_idx_c(dx), y + cell_idx_c(dy), z + cell_idx_c(dz));
            },
            liquidFlag, gasFlag, interfaceFlag, formerObstacleFlag, obstacleFlag, seedFillFromAllNeighbours_,
            seedSkipObstacleNeighbours_, fallback);
         const uint_t cnt = fallback ? uint_t(0) : uint_t(1);
         if (x >= 0 && y >= 0 && z >= 0 && x < cell_idx_c(flagField->xSize()) && y < cell_idx_c(flagField->ySize()) &&
             z < cell_idx_c(flagField->zSize()))
         {
            ++seededInterfaceCells_;
            if (fallback) ++seededInterfaceFallback_;
            seededInterfaceFillSum_ += seededFill;
         }

         // Logs interface cells seeded with the 0.5 fallback (no neighbour evidence, an uninformed start); first 20,
         // then every 5000th.
         if (cnt == uint_t(0) && block != nullptr)
         {
            static uint_t fallbackIsolatedLogCount = uint_t(0);
            ++fallbackIsolatedLogCount;
            if (fallbackIsolatedLogCount <= uint_t(20) || (fallbackIsolatedLogCount % uint_t(5000) == uint_t(0)))
            {
               const auto&  aabb = block->getAABB();
               const real_t cdx  = aabb.xSize() / real_c(flagField->xSize());
               WALBERLA_LOG_WARNING("[FALLBACK-ISOLATED-INTERFACE] (#" << fallbackIsolatedLogCount << ") cell ("
                  << x << ", " << y << ", " << z << ") world≈("
                  << (aabb.xMin() + (real_c(x) + real_t(0.5)) * cdx) << ","
                  << (aabb.yMin() + (real_c(y) + real_t(0.5)) * cdx) << ","
                  << (aabb.zMin() + (real_c(z) + real_t(0.5)) * cdx) << ") seeded fill=0.5 with ZERO real "
                  << (seedFillFromAllNeighbours_ ? "resolved (liquid/gas/interface)" : "interface")
                  << " neighbours -- a fragile, uninformed reconstruction that has been traced to a "
                     "real cumulant-collision NaN one step later; worth watching if this cell crashes soon.");
            }
         }
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
   // Hydrostatic reconstruction density (ModelParameters.useHydrostaticReconstruction). The EquilibriumReconstructor
   // uses the unweighted mean density of the domain neighbours, which ignores the fill: interface cells report rho ~ 1
   // however empty they are, so a cell freed at or above the waterline gets rho ~ 1 instead of the hydrostatic
   // rho_hyd(z) = 1 + 3 a_z (z - H) < 1, an over-density at the floater's waterline. When enabled, the cell is
   // re-seeded to rho_hyd(z) (see reseedHydrostatic). Does not use the fill field. Off by default.
   const bool useHydrostaticReconstruction_ = false;
   const real_t hydrostaticAccelerationZ_   = real_t(0); // a_z (≤ 0); same value fed to HydrostaticDensityFtor
   const real_t absoluteLiquidHeight_       = real_t(0); // H: reference liquid height [LU]
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
   bool       seedFillFromAllNeighbours_ = false; // see setSeedFillFromAllNeighbours
   bool       seedSkipObstacleNeighbours_ = true; // see setSeedSkipObstacleNeighbours
   bool       strictBognerClassification_ = false; // see setStrictBognerClassification
   uint_t     strictGasWithoutGas_        = uint_t(0);
   uint_t     strictLiquidTouchingGas_    = uint_t(0);
   bool       strictBognerNoLiquidAtGas_  = false;    // see setStrictBognerNoLiquidAtGas
   uint_t     strictLiquidToInterface_    = uint_t(0);
   // Relax pass (ModelParameters.strictBognerRelaxInterface; acts only together with strictBognerClassification_ and
   // strictBognerNoLiquidAtGas_): every inner freed cell that is now interface is set to gas or liquid if its decided
   // D3Q27 neighbourhood shows only one phase. Decide-then-apply; unlike usePass3Reresolution it reads the corrected
   // flags (LiquidPromoteRelaxClassificationTest.cpp finds no liquid-gas contact in its synthetic scenarios).
   bool       strictBognerRelaxInterface_ = false;
   uint_t     strictRelaxedToGas_         = uint_t(0);
   uint_t     strictRelaxedToLiquid_      = uint_t(0);
   uint_t     strictRelaxedStillInterface_= uint_t(0);
   // Seeding bookkeeping (inner cells, since the last logAndResetConversionCounts): interface cells
   // seeded, how many of them fell back to 0.5, and the sum of their seeded fill (= liquid volume created, in cells)
   mutable uint_t seededInterfaceCells_    = uint_t(0);
   mutable uint_t seededInterfaceFallback_ = uint_t(0);
   mutable real_t seededInterfaceFillSum_  = real_t(0);

   // FormerMO → state conversion counters over inner cells, one per branch of classifyFormerObstacleNeighbourhood (the
   // strict path uses convGas_, convLiquid_ and convIntfSurface_); flushed by logAndResetConversionCounts.
   uint_t convGas_            = uint_t(0); // branch A:  hasGas && !hasLiquid && !hasUnknown → gas
   uint_t convIntfUnknownGas_ = uint_t(0); // branch A': hasGas && !hasLiquid &&  hasUnknown → interface (fallback, symmetric to C)
   uint_t convLiquid_         = uint_t(0); // branch B:  hasLiquid && !hasGas && !hasUnknown → liquid
   uint_t convIntfUnknown_    = uint_t(0); // branch C:  hasLiquid && !hasGas &&  hasUnknown → interface (fallback)
   uint_t convIntfSurface_    = uint_t(0); // branch D:  hasGas &&  hasLiquid                → interface (real surface)
   uint_t convIntfUnresolved_ = uint_t(0); // branch E: !hasGas && !hasLiquid                → interface (fallback, no real fluid neighbour yet)

   // Unknown-neighbour statistics (non-strict path): for every inner FormerMO cell, how many of its 26 D3Q27 neighbours
   // are FormerMO cells of the same pass, i.e. also being reconstructed. Shows whether reconstructed cells lean on
   // resolved fluid or on other reconstructed cells. Flushed together with the conversion counts.
   uint_t formerMoTotal_          = uint_t(0); // total inner FormerMO cells classified this period
   uint_t formerMoWithUnknownNbr_ = uint_t(0); // of those, how many have >=1 unknown FormerMO neighbour
   uint_t unknownNbrSum_          = uint_t(0); // sum of unknown-neighbour counts (→ mean per cell)
   uint_t unknownNbrMax_          = uint_t(0); // largest unknown-neighbour count seen on any cell
   std::array< uint_t, 27 > unknownNbrHist_ = {}; // histogram indexed by #unknown neighbours (0..26)

   // Pass-3 conversion counters (inner cells only, matching the Pass-2 convention above). A Pass-3
   // decision is one of: promote to liquid/gas, or stay interface -- convPass3StillInterface_ counts
   // every "stay interface" outcome, convPass3UnresolvedIsland_ is the subset of those that are
   // genuinely !hasGas && !hasLiquid even after full same-pass resolution (see isPass3Candidate's
   // doc and the [PASS3-UNRESOLVED-ISLAND] log below) -- a subset, not an alternative bucket.
   uint_t convPass3ToLiquid_        = uint_t(0);
   uint_t convPass3ToGas_           = uint_t(0);
   uint_t convPass3StillInterface_  = uint_t(0);
   uint_t convPass3UnresolvedIsland_ = uint_t(0);
   std::array< real_t, 6 > freedMomentum_{}; // see takeFreedMomentum
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

   // Hydrostatic reconstruction (useHydrostaticReconstruction_), called right after the reconstructor: rescales the
   // freed cell's PDFs from the assigned density to rho_hyd(z) = 1 + 3 a_z (z - H) at the cell centre. A compressible
   // equilibrium is proportional to rho at fixed u, so this gives f_i^eq(rho_hyd, u_particle). Nearly a no-op in the
   // bulk; at or above the waterline it removes the rho ~ 1 over-density. Skipped if the assigned density is <= 1e-6
   // or rho_hyd <= 0.
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

   // Per-cell Bogner target flag (bognerTargetFlag below), evaluated inside each cell loop.
   //
   // The neighbourhood is the full D3Q27 stencil, not only the 6 face neighbours: the LBM stencil couples to edge and
   // corner neighbours too, and a liquid cell must never be a stencil neighbour of a gas cell (FSLBM invariant).
   // hasGas/hasLiquid count real gas/liquid neighbours only (see classifyFormerObstacleNeighbourhood). A cell becomes
   // liquid only with real liquid, no real gas and no unknown neighbour; with neither phase in reach it falls back to
   // interface, which the free-surface cell conversion resolves once a real phase is reachable.
   //
   // formerObstacleGhostBounds: Pass 2 iterates over the full ghost-inclusive range and the classification reads one
   // cell further, so a cell in the outermost ghost layer has neighbours outside the allocated field.
   const auto formerObstacleGhostBounds = flagField->xyzSizeWithGhostLayer();

   // Bounds-checked D3Q27 neighbour read, shared by bognerTargetFlag (Pass 2), the strict path and Pass 3. A neighbour
   // outside the allocated (ghost-inclusive) field is reported as formerObstacle, the same conservative fallback as a
   // cell freed in the same pass.
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

   // Pass-3 candidates (see isPass3Candidate), collected by bognerTargetFlag in Pass 2; processed only if
   // usePass3Reresolution_.
   struct Pass3Candidate { Cell cell; bool inner; };
   std::vector< Pass3Candidate > pass3Candidates;

   auto bognerTargetFlag = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) -> flag_t {
      // Scan and decision are the unit-tested pure functions at the top of this file
      // (tests/showcases_antidunes/BognerStateClassificationTest.cpp); this lambda adds the field access and the
      // per-branch counters.
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
         // Real liquid, no real gas, no unknown neighbour: bulk liquid.
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

         // Interface only because of an unknown neighbour (A', C) or with neither phase in reach (E); Pass 3
         // (usePass3Reresolution_) re-examines those with an unknown neighbour once every cell of the pass has a state.
         if (hasLiquid && !hasGas && hasUnknown)
         {
            if (inner) ++convIntfUnknown_;
            static uint_t provisionalLogCount = uint_t(0);
            ++provisionalLogCount;
            // Log disabled: it fires for dozens of cells per step and slows runs (every WALBERLA_LOG_INFO call pays
            // the ignore-regex check, see Logging::isInIgnoreCallerPaths).
            if (provisionalLogCount <= uint_t(20) || (provisionalLogCount % uint_t(5000) == uint_t(0)))
            {
               // WALBERLA_LOG_INFO("Bogner state selection: cell (" << cx << ", " << cy << ", " << cz << ") has "
               //    << unknownNbr << " unclassified FormerMO neighbour(s); "
               //    << (usePass3Reresolution_
               //           ? "provisionally interface pending Pass-3 re-resolution."
               //           : "falling back to interface instead of liquid to preserve the liquid/gas separation "
               //             "invariant."));
            }
         }
         else if (hasGas && !hasLiquid && hasUnknown)
         {
            // Symmetric case.
            if (inner) ++convIntfUnknownGas_;
            static uint_t provisionalGasLogCount = uint_t(0);
            ++provisionalGasLogCount;
            // Log disabled, as above.
            if (provisionalGasLogCount <= uint_t(20) || (provisionalGasLogCount % uint_t(5000) == uint_t(0)))
            {
               // WALBERLA_LOG_INFO("Bogner state selection: cell (" << cx << ", " << cy << ", " << cz << ") has "
               //    << unknownNbr << " unclassified FormerMO neighbour(s); "
               //    << (usePass3Reresolution_
               //           ? "provisionally interface pending Pass-3 re-resolution."
               //           : "falling back to interface instead of gas to preserve the liquid/gas separation "
               //             "invariant."));
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

   struct FreedCell { Cell cell; size_t particleIdx; }; // inner freed cells, for the momentum budget (takeFreedMomentum)
   std::vector< FreedCell > freedInner;

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
               if (pdfField->isInInnerPart(Cell(x, y, z))) freedInner.push_back({ Cell(x, y, z), idx });
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
            if (pdfField->isInInnerPart(Cell(x, y, z))) freedInner.push_back({ Cell(x, y, z), particleIdx });

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

   // Strict Bogner rule (strictBognerClassification_, see classifyStrictBognerNeighbourhood): every freed cell is
   // decided against the flags as they are before this pass (the other cells freed now still carry formerObstacle, so
   // they are not in B(x)), and the decisions are applied only afterwards, so the result does not depend on the loop
   // order. Interface fill: mean fill of the interface cells in B(x) (Bogner), 0.5 if there is none (Bogner gives no
   // rule for that case). Pass 3 has no candidates in this mode.
   if (useBognerStateSelection_ && strictBognerClassification_)
   {
      struct StrictDecision { Cell cell; flag_t target; real_t fill; };
      std::vector< StrictDecision > decisions;
      const cell_idx_t xs = cell_idx_c(flagField->xSize());
      const cell_idx_t ys = cell_idx_c(flagField->ySize());
      const cell_idx_t zs = cell_idx_c(flagField->zSize());
      auto isInner = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) {
         return cx >= 0 && cx < xs && cy >= 0 && cy < ys && cz >= 0 && cz < zs;
      };
      auto decide = [&](cell_idx_t cx, cell_idx_t cy, cell_idx_t cz) {
         auto nbr = [&](int dx, int dy, int dz) -> flag_t { return getBoundedNeighbourFlag(cx, cy, cz, dx, dy, dz); };
         bool hasGas, hasLiquid;
         const BognerNeighbourhoodState state = classifyFormerObstacleCellStrict(nbr, bognerLiquid, bognerGas, formerObstacle,
                                                                                 movingObstacle, hasGas, hasLiquid);
         const bool inner  = isInner(cx, cy, cz);
         flag_t     target = bognerInterface;
         real_t     fill   = real_t(0);
         if (state == BognerNeighbourhoodState::Gas)
         {
            target = bognerGas;
            if (inner) { ++convGas_; if (!hasGas) ++strictGasWithoutGas_; }
         }
         else if (state == BognerNeighbourhoodState::Liquid)
         {
            target = bognerLiquid;
            fill   = real_t(1);
            if (inner) ++convLiquid_;
         }
         else
         {
            if (inner) ++convIntfSurface_;
            bool fallback = false;
            if constexpr (!std::is_void_v< ScalarField_T >)
               if (fillField != nullptr)
                  fill = seedInterfaceFillFromNeighbourhood< flag_t >(
                     nbr,
                     [&](int dx, int dy, int dz) {
                        return fillField->get(cx + cell_idx_c(dx), cy + cell_idx_c(dy), cz + cell_idx_c(dz));
                     },
                     bognerLiquid, bognerGas, bognerInterface, formerObstacle, movingObstacle,
                     /* allResolvedNeighbours */ false, /* skipObstacleNeighbours */ true, fallback);
            if (inner)
            {
               ++seededInterfaceCells_;
               if (fallback) ++seededInterfaceFallback_;
               seededInterfaceFillSum_ += fill;
            }
         }
         decisions.push_back(StrictDecision{ Cell(cx, cy, cz), target, fill });
      };

      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(flagField->xyzSizeWithGhostLayer(),
                                             if (isFlagSet(flagField->get(x, y, z), formerObstacle)) { decide(x, y, z); })

      for (const auto& d : decisions)
      {
         boundaryHandling->setDomain(d.target, d.cell.x(), d.cell.y(), d.cell.z());
         if constexpr (!std::is_void_v< ScalarField_T >)
            if (fillField != nullptr) fillField->get(d.cell) = d.fill;
         removeFlag(flagField->get(d.cell), formerObstacle);
         particleField->get(d.cell) = ac_->getInvalidUid();
      }

      // Rank-local counts of this call only, logged below without MPI reduction (the strict* members accumulate and
      // are reduced in logAndResetConversionCounts).
      uint_t localPromotedThisSweep = uint_t(0), localRelaxedToGasThisSweep = uint_t(0),
             localRelaxedToLiquidThisSweep = uint_t(0);

      // A new liquid cell can touch gas when cells freed in the same step go different ways (Y has no gas in B(x) and
      // becomes liquid, its freed neighbour Z no liquid and becomes gas). Counted; with strictBognerNoLiquidAtGas_ the
      // liquid cell becomes interface with fill 1, since a liquid cell next to gas streams the gas cell's invalid PDFs.
      // Inner cells only: the flag field is communicated right after this sweep, so the owner's decision reaches the
      // ghost layers.
      for (auto& d : decisions)
      {
         if (d.target != bognerLiquid || !isInner(d.cell.x(), d.cell.y(), d.cell.z())) continue;
         bool touches = false;
         for (int dz = -1; dz <= 1 && !touches; ++dz)
            for (int dy = -1; dy <= 1 && !touches; ++dy)
               for (int dx = -1; dx <= 1 && !touches; ++dx)
               {
                  const flag_t nf = getBoundedNeighbourFlag(d.cell.x(), d.cell.y(), d.cell.z(), dx, dy, dz);
                  touches = (dx != 0 || dy != 0 || dz != 0) && isFlagSet(nf, bognerGas) && !isFlagSet(nf, movingObstacle);
               }
         if (touches) ++strictLiquidTouchingGas_;
         if (touches && strictBognerNoLiquidAtGas_)
         {
            boundaryHandling->setDomain(bognerInterface, d.cell.x(), d.cell.y(), d.cell.z());
            if constexpr (!std::is_void_v< ScalarField_T >)
               if (fillField != nullptr) fillField->get(d.cell) = real_t(1);
            ++strictLiquidToInterface_;
            ++localPromotedThisSweep;
            d.target = bognerInterface; // so the relax pass below sees the promoted state
         }
      }

      // Relax pass (strictBognerRelaxInterface_, needs strictBognerNoLiquidAtGas_): every inner freed cell that is now
      // interface is re-scanned against the flag field, which already holds the decisions above, and set to gas or
      // liquid if its neighbourhood shows only one phase; mixed or empty neighbourhoods stay interface.
      // Decide-then-apply, so relax decisions do not see each other.
      if (strictBognerNoLiquidAtGas_ && strictBognerRelaxInterface_)
      {
         struct RelaxDecision { Cell cell; flag_t target; };
         std::vector< RelaxDecision > relaxDecisions;
         for (const auto& d : decisions)
         {
            if (d.target != bognerInterface || !isInner(d.cell.x(), d.cell.y(), d.cell.z())) continue;
            bool hasGas = false, hasLiquid = false;
            for (int dz = -1; dz <= 1; ++dz)
               for (int dy = -1; dy <= 1; ++dy)
                  for (int dx = -1; dx <= 1; ++dx)
                  {
                     if (dx == 0 && dy == 0 && dz == 0) continue;
                     const flag_t nf = getBoundedNeighbourFlag(d.cell.x(), d.cell.y(), d.cell.z(), dx, dy, dz);
                     if (isFlagSet(nf, movingObstacle)) continue; // currently solid -- pure exclusion, as elsewhere
                     if (isFlagSet(nf, bognerLiquid)) hasLiquid = true;
                     if (isFlagSet(nf, bognerGas)) hasGas = true;
                  }
            if (hasGas && !hasLiquid) relaxDecisions.push_back({ d.cell, bognerGas });
            else if (hasLiquid && !hasGas) relaxDecisions.push_back({ d.cell, bognerLiquid });
            else { ++strictRelaxedStillInterface_; }
         }
         for (const auto& r : relaxDecisions)
         {
            boundaryHandling->setDomain(r.target, r.cell.x(), r.cell.y(), r.cell.z());
            if constexpr (!std::is_void_v< ScalarField_T >)
               if (fillField != nullptr) fillField->get(r.cell) = (r.target == bognerLiquid) ? real_t(1) : real_t(0);
            if (r.target == bognerGas) { ++strictRelaxedToGas_; ++localRelaxedToGasThisSweep; }
            else { ++strictRelaxedToLiquid_; ++localRelaxedToLiquidThisSweep; }
         }
      }

      // Logged per rank, only if this call promoted or relaxed a cell.
      const uint_t localRelaxedThisSweep = localRelaxedToGasThisSweep + localRelaxedToLiquidThisSweep;
      if (localPromotedThisSweep > uint_t(0) || localRelaxedThisSweep > uint_t(0))
      {
         WALBERLA_LOG_INFO("[FORMERMO-RELAX-SWEEP] rank=" << walberla::mpi::MPIManager::instance()->rank()
            << ": promoted liquid->interface=" << localPromotedThisSweep
            << " | relaxed interface->gas=" << localRelaxedToGasThisSweep
            << " | relaxed interface->liquid=" << localRelaxedToLiquidThisSweep);
      }
   }
   else if (optimizeForSmallObstacleFraction_)
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
                                           bognerInterface, block);
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
                                        bognerInterface, block);
            removeFlag(flagField->get(x, y, z), formerObstacle);
            particleField->get(x, y, z) = ac_->getInvalidUid();
         })
   }

   // ---- Pass 3: local re-resolution (usePass3Reresolution_ only) ------------------------------
   //
   // After Pass 2 every FormerMO cell of the block has a state, so the candidates (cells set to interface only because
   // of an unknown neighbour) are classified again; hasUnknown is now false except at the outermost ghost layer, where
   // getBoundedNeighbourFlag keeps such cells interface.
   //
   // Decide-then-apply: all candidates are classified against the same post-Pass-2 snapshot, so the loop order cannot
   // cascade a fresh promotion into a neighbour's decision.
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
               // Still no real liquid or gas within one cell after every cell of the pass is resolved: an anomaly,
               // logged separately; the cell stays interface (no basis to pick a phase).
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
                                     newFlag, bognerGas, bognerLiquid, bognerInterface, block);
         if (decision.inner)
         {
            if (decision.state == BognerNeighbourhoodState::Gas) { ++convPass3ToGas_; }
            else { ++convPass3ToLiquid_; }
         }
      }
   }

   // Momentum budget (takeFreedMomentum): the freed cells' final state decides how much of their PDFs'
   // momentum is fluid -- liquid all of it, interface its fill, gas none
   if constexpr (!std::is_void_v< ScalarField_T >)
      if (fillField != nullptr)
         for (const auto& f : freedInner)
         {
            const flag_t fl = flagField->get(f.cell);
            const real_t w  = isFlagSet(fl, bognerLiquid) ? real_t(1)
                              : isFlagSet(fl, bognerInterface) ? std::clamp(fillField->get(f.cell), real_t(0), real_t(1))
                                                                : real_t(0);
            if (w <= real_t(0)) continue;
            const Vector3< real_t > m = pdfField->getMomentumDensity(f.cell.x(), f.cell.y(), f.cell.z());
            const Vector3< real_t > r = blockStorage_->getBlockLocalCellCenter(*block, f.cell) - ac_->getPosition(f.particleIdx);
            freedMomentum_[0] += w * m[0];
            freedMomentum_[1] += w * m[1];
            freedMomentum_[2] += w * m[2];
            freedMomentum_[3] += w * (r[2] * m[0] - r[0] * m[2]);
            freedMomentum_[4] += w * pdfField->getDensity(f.cell.x(), f.cell.y(), f.cell.z());
            freedMomentum_[5] += real_t(1);
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
// the volume fraction φ of converted cells from `fillFieldID` (Bogner 2017, Sec. 5.2).
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