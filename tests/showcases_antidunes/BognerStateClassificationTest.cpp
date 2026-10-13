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
//! \file BognerStateClassificationTest.cpp
//! Unit tests for classifyFormerObstacleNeighbourhood()/isPass3Candidate() (PdfReconstructionManager.h),
//! the pure decision logic behind the Bogner & Rüde state-selection reconstruction's bognerTargetFlag
//! and its optional Pass-3 local re-resolution sweep (ModelParameters.usePass3Reresolution).
//!
//! Background: this function decides what state (gas/liquid/interface) a cell freshly uncovered
//! by a moving mesa_pd particle should get, from three booleans describing its D3Q27 neighbourhood:
//! hasGas/hasLiquid (REAL gas/liquid flags only -- an interface neighbour must count toward
//! NEITHER) and hasUnknown (an adjacent cell freed in the same reconstruction pass, not yet
//! classified). An earlier version of this code counted interface neighbours as evidence of both
//! gas and liquid, which caused a self-propagating misclassification in densely packed, actively
//! jittering particle beds: floating cylinders in cylinder_sedimentation_showcase_gas.prm ended up
//! surrounded by persistent "interface" cells 30+ LU below the true free surface, with oscillating
//! (17-86%) submersion instead of the ~100% a fully submerged particle should read. These tests
//! cover the full truth table of the three input booleans so a future change to this function
//! cannot silently reintroduce that bug or a similar misclassification.
//!
//! The two-round (Pass 2 then Pass 3) tests further down exercise isPass3Candidate() and the same
//! classifier re-run against a mutated neighbourhood, simulating what PdfReconstructionManager's
//! Pass-3 decide sub-loop sees once every same-pass FormerMO cell in the block has a final state.
//======================================================================================================================

#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace lbm_mesapd_coupling {

// ---------------------------------------------------------------------------------------------
// Ground truth: the documented decision table in classifyFormerObstacleNeighbourhood's own
// docstring (PdfReconstructionManager.h). Every one of the 8 combinations of
// (hasGas, hasLiquid, hasUnknown) is exercised explicitly, so a future edit cannot silently change
// a branch's outcome (or reintroduce interface-counts-as-gas/liquid) without a test failing.
// ---------------------------------------------------------------------------------------------

void testGasOnlyFullyResolved()
{
   // hasGas && !hasLiquid && !hasUnknown -> Gas: provably bulk gas, no unresolved neighbour could
   // still turn out to be liquid.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ false, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Gas);
}

void testGasOnlyWithUnresolvedNeighbour()
{
   // hasGas && !hasLiquid && hasUnknown -> Interface: symmetric to the liquid-side case below. An
   // adjacent cell freed in the same pass hasn't been classified yet and could still become liquid,
   // so committing to gas now would risk the liquid-never-touches-gas invariant retroactively, just
   // discovered from the gas side instead of the liquid side.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ false, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

void testLiquidOnlyFullyResolved()
{
   // hasLiquid && !hasGas && !hasUnknown -> Liquid: provably bulk liquid, no unresolved neighbour
   // could still turn out to be gas. This is the case that used to be forced to Interface.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ true, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Liquid);
}

void testLiquidOnlyWithUnresolvedNeighbour()
{
   // hasLiquid && !hasGas && hasUnknown -> Interface: an adjacent cell freed in the same pass
   // hasn't been classified yet and could still become gas, so committing to liquid now would
   // risk the liquid-never-touches-gas invariant retroactively.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ true, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

void testMixedNeighbourhood()
{
   // hasGas && hasLiquid -> Interface: a genuine free-surface cell, independent of hasUnknown.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ true, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ true, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

void testNeitherPhaseReachable()
{
   // !hasGas && !hasLiquid -> Interface (safe fallback): only interface and/or still-unresolved
   // FormerMO neighbours nearby -- neither promotion is provably safe. This is exactly the
   // situation inside a densely packed, actively-jittering particle bed once a few cells nearby
   // have already been (mis)classified interface; the fallback must stay interface here rather
   // than silently defaulting to gas or liquid.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ false, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ false, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

// ---------------------------------------------------------------------------------------------
// Cube-based tests: drive classifyFormerObstacleCell() -- the SAME function bognerTargetFlag
// calls in production -- with a literal 3x3x3 neighbourhood of flags, exactly as a real FlagField
// scan would see it, instead of pre-computed hasGas/hasLiquid/hasUnknown booleans. This exercises
// the neighbour-scanning loop itself, not just the downstream decision table.
//
// A tiny local flag type stands in for the real (BoundaryHandling_T-dependent) flag_t; the only
// requirement on Flag_T is that field::isFlagSet(value, testFlag) works, which holds for any
// unsigned integer used as a bitmask -- exactly how the real flag_t behaves.
// ---------------------------------------------------------------------------------------------

using MockFlag = uint8_t;
constexpr MockFlag LIQUID          = MockFlag(1) << 0;
constexpr MockFlag GAS             = MockFlag(1) << 1;
constexpr MockFlag INTERFACE       = MockFlag(1) << 2;
constexpr MockFlag FORMER_OBSTACLE = MockFlag(1) << 3;
constexpr MockFlag OBSTACLE        = MockFlag(1) << 4; // currently-covered flag, excluded from all evidence

// [dz+1][dy+1][dx+1]; index [1][1][1] (the centre) is never read by classifyFormerObstacleCell,
// so it's left at whatever the cube was filled with.
using Cube = std::array< std::array< std::array< MockFlag, 3 >, 3 >, 3 >;

Cube uniformCube(MockFlag fill)
{
   Cube cube;
   for (auto& plane : cube)
      for (auto& row : plane)
         row.fill(fill);
   return cube;
}

BognerNeighbourhoodState classifyCube(const Cube& cube)
{
   bool hasGas, hasLiquid, hasUnknown;
   uint_t unknownNbr;
   return classifyFormerObstacleCell(
      [&](int dx, int dy, int dz) {
         // dx,dy,dz in {-1,0,1}; casting first (per clang-tidy's own preferred form) then adding 1
         // relies on well-defined unsigned wraparound for -1 -> SIZE_MAX -> (+1, mod 2^64) -> 0.
         const size_t ix = static_cast< size_t >(dx) + size_t(1);
         const size_t iy = static_cast< size_t >(dy) + size_t(1);
         const size_t iz = static_cast< size_t >(dz) + size_t(1);
         return cube[iz][iy][ix];
      },
      LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas, hasLiquid, hasUnknown, unknownNbr);
}

// Same as classifyCube(), but also hands back hasGas/hasLiquid/hasUnknown -- needed by the
// two-round (Pass 2 then Pass 3) tests below, which must feed isPass3Candidate() the same booleans
// production code computes at classification time.
BognerNeighbourhoodState classifyCubeFull(const Cube& cube, bool& hasGas, bool& hasLiquid, bool& hasUnknown)
{
   uint_t unknownNbr;
   return classifyFormerObstacleCell(
      [&](int dx, int dy, int dz) {
         const size_t ix = static_cast< size_t >(dx) + size_t(1);
         const size_t iy = static_cast< size_t >(dy) + size_t(1);
         const size_t iz = static_cast< size_t >(dz) + size_t(1);
         return cube[iz][iy][ix];
      },
      LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas, hasLiquid, hasUnknown, unknownNbr);
}

void testCubeTopLeftInterfaceRestLiquid()
{
   // The user's own example: one corner of the surrounding cube -- (dx,dy,dz) = (-1,-1,-1), "top
   // left" of the neighbourhood -- is interface, the other 25 cells are real liquid, no real gas
   // anywhere. This is the DIRECT regression test for the fix: under the old code (interface
   // counted as both gas- and liquid-equivalent), that one interface neighbour alone would have
   // set hasGas=true, making the whole cell "mixed" -> Interface. Under the fixed code, an
   // interface neighbour counts toward neither, so hasGas stays false and this is unambiguously
   // bulk liquid.
   Cube cube = uniformCube(LIQUID);
   cube[0][0][0] = INTERFACE; // (dx,dy,dz) = (-1,-1,-1)
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Liquid);
}

void testCubeAllLiquid()
{
   WALBERLA_CHECK(classifyCube(uniformCube(LIQUID)) == BognerNeighbourhoodState::Liquid);
}

void testCubeAllGas()
{
   WALBERLA_CHECK(classifyCube(uniformCube(GAS)) == BognerNeighbourhoodState::Gas);
}

void testCubeOneGasRestLiquid()
{
   // A single real gas neighbour amid otherwise-liquid surroundings is a genuine free surface ->
   // Interface, regardless of where in the cube it sits.
   Cube cube = uniformCube(LIQUID);
   cube[2][2][2] = GAS; // (dx,dy,dz) = (+1,+1,+1), the opposite corner from the case above
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Interface);
}

void testCubeAllInterfaceNoRealPhase()
{
   // The entire neighbourhood is interface -- no real gas or liquid reachable within one cell at
   // all. This is exactly the state a cell deep inside a densely packed, actively-jittering
   // particle bed can be in once its immediate surroundings have themselves already been
   // classified interface: neither promotion is provably safe, so the fallback must stay
   // Interface rather than silently defaulting to Gas or Liquid.
   WALBERLA_CHECK(classifyCube(uniformCube(INTERFACE)) == BognerNeighbourhoodState::Interface);
}

void testCubeFormerObstacleNeighbourBlocksLiquidPromotion()
{
   // All-liquid except one still-unresolved FormerMO neighbour, no real gas anywhere -> Interface
   // (branch C): that neighbour could still turn out to be gas later in the same reconstruction
   // pass, so promoting this cell to liquid now would risk the invariant retroactively.
   Cube cube = uniformCube(LIQUID);
   cube[1][0][2] = FORMER_OBSTACLE; // (dx,dy,dz) = (+1,-1,0), an arbitrary non-centre cell
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Interface);
}

void testCubeFormerObstacleNeighbourBlocksGasPromotion()
{
   // Symmetric to the test above: all-gas except one still-unresolved FormerMO neighbour, no real
   // liquid anywhere -> Interface (branch A'): that neighbour could still turn out to be liquid
   // later in the same reconstruction pass, so promoting this cell to gas now would risk the
   // invariant retroactively, same as the liquid-side case, just from the other direction.
   Cube cube = uniformCube(GAS);
   cube[1][0][2] = FORMER_OBSTACLE; // (dx,dy,dz) = (+1,-1,0), same arbitrary non-centre cell
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Interface);
}

// ---------------------------------------------------------------------------------------------
// Stale-flag regression tests: a real cluster run (81-sphere drop, sphere_only_81floaters_*) showed
// via direct VTK inspection that MovingParticleMapping.h can leave a cell simultaneously flagged
// OBSTACLE (currently covered by a particle) and a stale real GAS or LIQUID bit from before it was
// covered (root cause: the mapping only clears/sets flags for cells that were `isDomain()` --
// liquid|interface -- when swept over; a cell that was gas at that moment matches no branch and is
// silently skipped, keeping its old gas bit indefinitely). Without an explicit exclusion, that stale
// bit would be read as genuine phase evidence by any FormerMO neighbour scanning it. These tests
// model that exact corrupted cell directly and confirm the obstacle check in
// classifyFormerObstacleCell neutralises it completely, regardless of what real-phase bit it also
// carries.
// ---------------------------------------------------------------------------------------------

void testObstacleExcludesStaleGasEvidence()
{
   // All-liquid except one neighbour that is BOTH stale-GAS and currently-OBSTACLE. If the obstacle
   // exclusion works, that neighbour contributes no evidence at all -- hasGas must stay false, and
   // the cell must classify as unambiguous bulk Liquid, not Interface.
   Cube cube = uniformCube(LIQUID);
   cube[1][0][2] = GAS | OBSTACLE; // (dx,dy,dz) = (+1,-1,0)
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Liquid);
}

void testObstacleExcludesStaleLiquidEvidence()
{
   // Symmetric case: all-gas except one neighbour that is BOTH stale-LIQUID and currently-OBSTACLE.
   // hasLiquid must stay false -- confirms the exclusion isn't accidentally gas-only.
   Cube cube = uniformCube(GAS);
   cube[1][0][2] = LIQUID | OBSTACLE;
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Gas);
}

// ---------------------------------------------------------------------------------------------
// isPass3Candidate() truth table: should a cell that fell back to Interface in Pass 2 be
// re-examined by Pass 3 once every same-pass FormerMO cell in the block has a final state? True
// only for branches A'/C (blocked purely by an unresolved same-pass neighbour) and branch E when it
// ALSO had an unresolved neighbour -- never for branch D (a genuine, already-correct surface), and
// never for a branch-E cell with zero unresolved neighbours (nothing in its neighbourhood can
// change by Pass-3 time, so it is structurally unrescuable).
// ---------------------------------------------------------------------------------------------

void testPass3CandidateTruthTable()
{
   // Never a candidate for a non-Interface outcome.
   WALBERLA_CHECK(!isPass3Candidate(BognerNeighbourhoodState::Gas, true, false, false));
   WALBERLA_CHECK(!isPass3Candidate(BognerNeighbourhoodState::Liquid, false, true, false));

   // Branch D (hasGas && hasLiquid) -> never a candidate, regardless of hasUnknown.
   WALBERLA_CHECK(!isPass3Candidate(BognerNeighbourhoodState::Interface, true, true, false));
   WALBERLA_CHECK(!isPass3Candidate(BognerNeighbourhoodState::Interface, true, true, true));

   // Branch A' (hasGas && !hasLiquid && hasUnknown) -> candidate.
   WALBERLA_CHECK(isPass3Candidate(BognerNeighbourhoodState::Interface, true, false, true));
   // Branch C (hasLiquid && !hasGas && hasUnknown) -> candidate.
   WALBERLA_CHECK(isPass3Candidate(BognerNeighbourhoodState::Interface, false, true, true));

   // Branch E (!hasGas && !hasLiquid): candidate only when hasUnknown is true.
   WALBERLA_CHECK(isPass3Candidate(BognerNeighbourhoodState::Interface, false, false, true));
   WALBERLA_CHECK(!isPass3Candidate(BognerNeighbourhoodState::Interface, false, false, false));
}

// ---------------------------------------------------------------------------------------------
// Two-round (Pass 2 then Pass 3) tests: call classifyCube[Full]() twice on the same cube, mutating
// the one FormerMO neighbour between rounds to simulate it being resolved elsewhere in the same
// Pass-2 sweep -- exactly what Pass 3 sees when it re-examines a queued candidate against the
// now-fully-known neighbourhood (PdfReconstructionManager.h's decide sub-loop).
// ---------------------------------------------------------------------------------------------

void testCubeTwoRoundLiquidPromotion()
{
   // Round 1 (Pass 2): all-liquid except one still-unresolved FormerMO neighbour -> Interface
   // (branch C), correctly identified as a Pass-3 candidate.
   Cube cube = uniformCube(LIQUID);
   cube[1][0][2] = FORMER_OBSTACLE; // (dx,dy,dz) = (+1,-1,0)
   bool hasGas, hasLiquid, hasUnknown;
   const auto round1 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round1 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(isPass3Candidate(round1, hasGas, hasLiquid, hasUnknown));

   // Round 2 (Pass 3): that neighbour has since resolved to real LIQUID elsewhere in the same
   // Pass-2 sweep -> re-classifying against the now-fully-known neighbourhood correctly promotes to
   // Liquid.
   cube[1][0][2] = LIQUID;
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Liquid);
}

void testCubeTwoRoundGasPromotion()
{
   // Symmetric to the liquid case above (branch A').
   Cube cube = uniformCube(GAS);
   cube[1][0][2] = FORMER_OBSTACLE;
   bool hasGas, hasLiquid, hasUnknown;
   const auto round1 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round1 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(isPass3Candidate(round1, hasGas, hasLiquid, hasUnknown));

   cube[1][0][2] = GAS;
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Gas);
}

void testCubeTwoRoundUnresolvedIsland()
{
   // Round 1 (Pass 2): the entire neighbourhood is interface except one still-unresolved FormerMO
   // cell -- no real gas or liquid anywhere -> Interface (branch E), and (since hasUnknown is true)
   // correctly identified as a Pass-3 candidate.
   Cube cube = uniformCube(INTERFACE);
   cube[1][0][2] = FORMER_OBSTACLE;
   bool hasGas, hasLiquid, hasUnknown;
   const auto round1 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round1 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(isPass3Candidate(round1, hasGas, hasLiquid, hasUnknown));

   // Round 2 (Pass 3): that neighbour has since resolved to a genuine, boring INTERFACE cell
   // elsewhere in the same Pass-2 sweep (not a real phase -- e.g. it turned out to have real gas
   // AND real liquid in its own neighbourhood, branch D). Two things must both hold: (1) the
   // candidate must NOT incorrectly cascade to liquid or gas -- an interface neighbour is still not
   // evidence of either, so it correctly stays Interface; (2) hasGas/hasLiquid must both still be
   // false, which is exactly the [PASS3-UNRESOLVED-ISLAND] trigger condition in production
   // (PdfReconstructionManager.h) -- this cell is genuinely isolated from both real phases even
   // after full same-pass resolution, not simply still-pending.
   cube[1][0][2] = INTERFACE;
   const auto round2 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round2 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(!hasGas);
   WALBERLA_CHECK(!hasLiquid);
   WALBERLA_CHECK(!hasUnknown);
}

// ---------------------------------------------------------------------------------------------
// [Claude, 2026-10-01] seedInterfaceFillFromNeighbourhood (PdfReconstructionManager.h): fill level seeded into a freshly
// uncovered cell that became INTERFACE. Legacy rule: mean fill of the interface neighbours, 0.5 if none. New rule
// (prm seedFillFromAllNeighbours): mean volume fraction of all resolved neighbours (gas 0, liquid 1, interface fill),
// excluding obstacle and still-unclassified (formerObstacle) cells.
// ---------------------------------------------------------------------------------------------
using FillCube = std::array< std::array< std::array< real_t, 3 >, 3 >, 3 >;

FillCube uniformFill(real_t v)
{
   FillCube c;
   for (auto& plane : c)
      for (auto& row : plane)
         row.fill(v);
   return c;
}

real_t seedCube(const Cube& flags, const FillCube& fills, bool allNeighbours, bool& fallback, bool skipObstacle = true)
{
   auto idx = [](int d) { return static_cast< size_t >(d) + size_t(1); };
   return seedInterfaceFillFromNeighbourhood< MockFlag >(
      [&](int dx, int dy, int dz) { return flags[idx(dz)][idx(dy)][idx(dx)]; },
      [&](int dx, int dy, int dz) { return fills[idx(dz)][idx(dy)][idx(dx)]; }, LIQUID, GAS, INTERFACE, FORMER_OBSTACLE,
      OBSTACLE, allNeighbours, skipObstacle, fallback);
}

void checkSeed(const char* name, const Cube& flags, const FillCube& fills, real_t expectLegacy, bool expectLegacyFallback,
               real_t expectNew, bool expectNewFallback)
{
   bool         fbL = false, fbN = false;
   const real_t l = seedCube(flags, fills, false, fbL);
   const real_t n = seedCube(flags, fills, true, fbN);
   WALBERLA_LOG_INFO_ON_ROOT("[seed] " << name << ": legacy " << l << (fbL ? " (fallback)" : "") << ", all-neighbours " << n
                                       << (fbN ? " (fallback)" : ""));
   WALBERLA_CHECK_FLOAT_EQUAL(l, expectLegacy, name);
   WALBERLA_CHECK_EQUAL(fbL, expectLegacyFallback, name);
   WALBERLA_CHECK_FLOAT_EQUAL(n, expectNew, name);
   WALBERLA_CHECK_EQUAL(fbN, expectNewFallback, name);
}

void testSeedBulkLiquidWithUnknownNeighbour()
{
   // branch C: liquid all around, one still-unclassified neighbour -> interface; legacy seeds a half-empty cell
   Cube c = uniformCube(LIQUID);
   c[1][1][2] = FORMER_OBSTACLE;
   checkSeed("bulk liquid + unknown (C)", c, uniformFill(real_t(1)), real_t(0.5), true, real_t(1), false);
}

void testSeedAirWithUnknownNeighbour()
{
   // branch A': gas all around, one unclassified neighbour -> interface; legacy seeds half a cell of water in the air
   Cube c = uniformCube(GAS);
   c[1][1][0] = FORMER_OBSTACLE;
   checkSeed("air + unknown (A')", c, uniformFill(real_t(0)), real_t(0.5), true, real_t(0), false);
}

void testSeedSurfaceLiquidBelowGasAbove()
{
   // branch D without interface neighbours: lower plane liquid, upper plane gas, middle plane obstacle (a floater edge)
   Cube c = uniformCube(OBSTACLE);
   for (auto& row : c[0]) row.fill(LIQUID);
   for (auto& row : c[2]) row.fill(GAS);
   checkSeed("surface, liquid below / gas above", c, uniformFill(real_t(0)), real_t(0.5), true, real_t(0.5), false);
}

void testSeedWithInterfaceNeighbours()
{
   // 9 liquid (below), 9 interface with fill 0.4 (middle plane, 8 + centre ignored), 9 gas (above)
   Cube     c = uniformCube(INTERFACE);
   FillCube f = uniformFill(real_t(0.4));
   for (auto& row : c[0]) row.fill(LIQUID);
   for (auto& row : c[2]) row.fill(GAS);
   // legacy: mean of the 8 interface neighbours = 0.4; new: (9*1 + 8*0.4 + 9*0) / 26
   checkSeed("surface with interface neighbours", c, f, real_t(0.4), false,
             (real_t(9) + real_t(8) * real_t(0.4)) / real_t(26), false);
}

void testSeedNoResolvedNeighbour()
{
   // only obstacle and unclassified neighbours: no evidence either way -> 0.5 in both rules
   Cube c = uniformCube(OBSTACLE);
   c[0][0][0] = FORMER_OBSTACLE;
   checkSeed("no resolved neighbour", c, uniformFill(real_t(0.7)), real_t(0.5), true, real_t(0.5), true);
}

void testSeedIgnoresStaleFillOfExcludedCells()
{
   // obstacle and unclassified cells must not contribute even if a stale fill value is stored there
   Cube     c = uniformCube(OBSTACLE);
   FillCube f = uniformFill(real_t(0.9));
   c[0][1][1] = LIQUID;
   c[2][1][1] = GAS;
   c[1][0][1] = FORMER_OBSTACLE;
   checkSeed("excluded cells with stale fill", c, f, real_t(0.5), true, real_t(0.5), false);
}

void testSeedSkipObstacleNeighboursToggle()
{
   // [Claude, 2026-10-08] an interface neighbour that still carries the formerObstacle flag: skipped by default (since
   // 2026-10-01) -> no evidence, 0.5 fallback; counted with skipObstacle = false (legacy rule before 2026-10-01) -> its fill.
   // The all-neighbour rule always skips it.
   Cube     c = uniformCube(OBSTACLE);
   FillCube f = uniformFill(real_t(0.3));
   c[1][1][2] = MockFlag(INTERFACE | FORMER_OBSTACLE);
   bool fb = false;
   real_t v = seedCube(c, f, false, fb, true);
   WALBERLA_CHECK_FLOAT_EQUAL(v, real_t(0.5), "skip: no evidence");
   WALBERLA_CHECK(fb);
   v = seedCube(c, f, false, fb, false);
   WALBERLA_CHECK_FLOAT_EQUAL(v, real_t(0.3), "no skip: interface fill counted");
   WALBERLA_CHECK(!fb);
   v = seedCube(c, f, true, fb, false);
   WALBERLA_CHECK_FLOAT_EQUAL(v, real_t(0.5), "all-neighbour rule always skips");
   WALBERLA_CHECK(fb);
   WALBERLA_LOG_INFO_ON_ROOT("[seed] skipObstacleNeighbours toggle: ok");
}

// [Claude, 2026-10-08] strict Bogner rule (classifyStrictBognerNeighbourhood / classifyFormerObstacleCellStrict): Bogner's
// three cases in his order, and cells freed in the same step (formerObstacle) are not in B(x).
BognerNeighbourhoodState classifyCubeStrict(const Cube& cube, bool& hasGas, bool& hasLiquid)
{
   return classifyFormerObstacleCellStrict(
      [&](int dx, int dy, int dz) {
         return cube[static_cast< size_t >(dz) + size_t(1)][static_cast< size_t >(dy) + size_t(1)]
                    [static_cast< size_t >(dx) + size_t(1)];
      },
      LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas, hasLiquid);
}

void testStrictBognerTruthTable()
{
   WALBERLA_CHECK(classifyStrictBognerNeighbourhood(true, false) == BognerNeighbourhoodState::Gas);
   WALBERLA_CHECK(classifyStrictBognerNeighbourhood(false, true) == BognerNeighbourhoodState::Liquid);
   WALBERLA_CHECK(classifyStrictBognerNeighbourhood(true, true) == BognerNeighbourhoodState::Interface);
   // neither: "B(x) contains no liquid" is Bogner's first case -> gas
   WALBERLA_CHECK(classifyStrictBognerNeighbourhood(false, false) == BognerNeighbourhoodState::Gas);
}

void testStrictBognerSameStepNeighbours()
{
   bool g, l;
   // dry upper face: gas above, the rest of the freed sheet beside it, body below -> gas (the fallback rule gives interface)
   Cube c = uniformCube(OBSTACLE);
   for (auto& row : c[2]) row.fill(GAS);
   for (auto& row : c[1]) row.fill(FORMER_OBSTACLE);
   WALBERLA_CHECK(classifyCubeStrict(c, g, l) == BognerNeighbourhoodState::Gas);
   WALBERLA_CHECK(g && !l);
   WALBERLA_CHECK(classifyCube(c) == BognerNeighbourhoodState::Interface);
   // submerged face: liquid below, freed sheet beside -> liquid
   c = uniformCube(OBSTACLE);
   for (auto& row : c[0]) row.fill(LIQUID);
   for (auto& row : c[1]) row.fill(FORMER_OBSTACLE);
   WALBERLA_CHECK(classifyCubeStrict(c, g, l) == BognerNeighbourhoodState::Liquid);
   // waterline: liquid below, gas above -> interface
   for (auto& row : c[2]) row.fill(GAS);
   WALBERLA_CHECK(classifyCubeStrict(c, g, l) == BognerNeighbourhoodState::Interface);
   // a stale gas bit on a covered cell is no evidence
   c = uniformCube(LIQUID);
   c[2][1][1] = MockFlag(GAS | OBSTACLE);
   WALBERLA_CHECK(classifyCubeStrict(c, g, l) == BognerNeighbourhoodState::Liquid);
   // a formerObstacle neighbour that still carries a gas bit is not in B(x) either
   c = uniformCube(LIQUID);
   c[2][1][1] = MockFlag(GAS | FORMER_OBSTACLE);
   WALBERLA_CHECK(classifyCubeStrict(c, g, l) == BognerNeighbourhoodState::Liquid);
   WALBERLA_LOG_INFO_ON_ROOT("[strict] Bogner rule: ok");
}

} // namespace lbm_mesapd_coupling
} // namespace walberla

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::lbm_mesapd_coupling;
   testGasOnlyFullyResolved();
   testGasOnlyWithUnresolvedNeighbour();
   testLiquidOnlyFullyResolved();
   testLiquidOnlyWithUnresolvedNeighbour();
   testMixedNeighbourhood();
   testNeitherPhaseReachable();

   testCubeTopLeftInterfaceRestLiquid();
   testCubeAllLiquid();
   testCubeAllGas();
   testCubeOneGasRestLiquid();
   testCubeAllInterfaceNoRealPhase();
   testCubeFormerObstacleNeighbourBlocksLiquidPromotion();
   testCubeFormerObstacleNeighbourBlocksGasPromotion();

   testObstacleExcludesStaleGasEvidence();
   testObstacleExcludesStaleLiquidEvidence();

   testPass3CandidateTruthTable();
   testCubeTwoRoundLiquidPromotion();
   testCubeTwoRoundGasPromotion();
   testCubeTwoRoundUnresolvedIsland();

   testSeedBulkLiquidWithUnknownNeighbour();
   testSeedAirWithUnknownNeighbour();
   testSeedSurfaceLiquidBelowGasAbove();
   testSeedWithInterfaceNeighbours();
   testSeedNoResolvedNeighbour();
   testSeedIgnoresStaleFillOfExcludedCells();
   testSeedSkipObstacleNeighboursToggle();

   testStrictBognerTruthTable();
   testStrictBognerSameStepNeighbours();

   return EXIT_SUCCESS;
}
