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
//! Unit tests for the pure decision logic in PdfReconstructionManager.h behind the Bogner & Rüde state-selection
//! reconstruction: classifyFormerObstacleNeighbourhood()/classifyFormerObstacleCell() (state of a cell freshly
//! uncovered by a moving particle), isPass3Candidate() (optional Pass-3 re-resolution,
//! ModelParameters.usePass3Reresolution), seedInterfaceFillFromNeighbourhood() (initial fill of such a cell) and the
//! strict Bogner rule (classifyFormerObstacleCellStrict()).
//!
//! The classifier takes three booleans from the D3Q27 neighbourhood: hasGas/hasLiquid (real gas/liquid flags only; an
//! interface neighbour counts toward neither) and hasUnknown (a neighbour freed in the same pass, not yet classified).
//! Counting interface neighbours as both gas and liquid evidence makes interface cells propagate through densely
//! packed particle beds far below the free surface; the tests cover the full truth table to guard against that.
//======================================================================================================================

#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace lbm_mesapd_coupling {

// ---------------------------------------------------------------------------------------------
// Decision table of classifyFormerObstacleNeighbourhood (see its docstring): all 8 combinations of
// (hasGas, hasLiquid, hasUnknown).
// ---------------------------------------------------------------------------------------------

void testGasOnlyFullyResolved()
{
   // hasGas && !hasLiquid && !hasUnknown -> Gas: no unresolved neighbour can still turn out liquid.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ false, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Gas);
}

void testGasOnlyWithUnresolvedNeighbour()
{
   // hasGas && !hasLiquid && hasUnknown -> Interface (branch A'): the unresolved neighbour could still become
   // liquid, so committing to gas could break the liquid-never-touches-gas invariant.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ false, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

void testLiquidOnlyFullyResolved()
{
   // hasLiquid && !hasGas && !hasUnknown -> Liquid: no unresolved neighbour can still turn out gas.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ true, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Liquid);
}

void testLiquidOnlyWithUnresolvedNeighbour()
{
   // hasLiquid && !hasGas && hasUnknown -> Interface (branch C): the unresolved neighbour could still become gas.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ true, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

void testMixedNeighbourhood()
{
   // hasGas && hasLiquid -> Interface (branch D): a genuine free-surface cell, independent of hasUnknown.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ true, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ true, /*hasLiquid*/ true, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

void testNeitherPhaseReachable()
{
   // !hasGas && !hasLiquid -> Interface (branch E, safe fallback): only interface and/or unresolved FormerMO
   // neighbours, so neither promotion is provably safe.
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ false, /*hasUnknown*/ false)
                  == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(classifyFormerObstacleNeighbourhood(/*hasGas*/ false, /*hasLiquid*/ false, /*hasUnknown*/ true)
                  == BognerNeighbourhoodState::Interface);
}

// ---------------------------------------------------------------------------------------------
// Cube-based tests: drive classifyFormerObstacleCell() (the function bognerTargetFlag calls) with a literal 3x3x3
// flag neighbourhood, so the neighbour-scanning loop is exercised as well as the decision table.
// MockFlag stands in for the real flag_t; Flag_T only needs field::isFlagSet(value, testFlag), which any unsigned
// bitmask satisfies.
// ---------------------------------------------------------------------------------------------

using MockFlag = uint8_t;
constexpr MockFlag LIQUID          = MockFlag(1) << 0;
constexpr MockFlag GAS             = MockFlag(1) << 1;
constexpr MockFlag INTERFACE       = MockFlag(1) << 2;
constexpr MockFlag FORMER_OBSTACLE = MockFlag(1) << 3;
constexpr MockFlag OBSTACLE        = MockFlag(1) << 4; // currently-covered flag, excluded from all evidence

// [dz+1][dy+1][dx+1]; the centre [1][1][1] is never read by classifyFormerObstacleCell.
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
         // dx,dy,dz in {-1,0,1}; casting first, then adding 1, uses well-defined unsigned wraparound (-1 -> 0).
         const size_t ix = static_cast< size_t >(dx) + size_t(1);
         const size_t iy = static_cast< size_t >(dy) + size_t(1);
         const size_t iz = static_cast< size_t >(dz) + size_t(1);
         return cube[iz][iy][ix];
      },
      LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas, hasLiquid, hasUnknown, unknownNbr);
}

// Same as classifyCube(), but also returns hasGas/hasLiquid/hasUnknown for isPass3Candidate().
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
   // One corner neighbour (-1,-1,-1) is interface, the other 25 are liquid. An interface neighbour counts toward
   // neither phase, so hasGas stays false and the cell is bulk liquid (not "mixed" -> Interface).
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
   // All-interface neighbourhood, no real phase within one cell: the fallback must stay Interface.
   WALBERLA_CHECK(classifyCube(uniformCube(INTERFACE)) == BognerNeighbourhoodState::Interface);
}

void testCubeFormerObstacleNeighbourBlocksLiquidPromotion()
{
   // All-liquid except one unresolved FormerMO neighbour -> Interface (branch C).
   Cube cube = uniformCube(LIQUID);
   cube[1][0][2] = FORMER_OBSTACLE; // (dx,dy,dz) = (+1,-1,0), an arbitrary non-centre cell
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Interface);
}

void testCubeFormerObstacleNeighbourBlocksGasPromotion()
{
   // All-gas except one unresolved FormerMO neighbour -> Interface (branch A').
   Cube cube = uniformCube(GAS);
   cube[1][0][2] = FORMER_OBSTACLE; // (dx,dy,dz) = (+1,-1,0), same arbitrary non-centre cell
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Interface);
}

// ---------------------------------------------------------------------------------------------
// Stale-flag tests: MovingParticleMapping.h can leave a cell flagged OBSTACLE together with a stale phase bit from
// before it was covered (the mapping only rewrites flags of cells that were isDomain(), i.e. liquid|interface, so a
// gas cell keeps its gas bit). The obstacle check in classifyFormerObstacleCell must ignore such a neighbour whatever
// phase bit it carries.
// ---------------------------------------------------------------------------------------------

void testObstacleExcludesStaleGasEvidence()
{
   // All-liquid except one GAS|OBSTACLE neighbour: no gas evidence, so Liquid rather than Interface.
   Cube cube = uniformCube(LIQUID);
   cube[1][0][2] = GAS | OBSTACLE; // (dx,dy,dz) = (+1,-1,0)
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Liquid);
}

void testObstacleExcludesStaleLiquidEvidence()
{
   // All-gas except one LIQUID|OBSTACLE neighbour: hasLiquid stays false, so Gas.
   Cube cube = uniformCube(GAS);
   cube[1][0][2] = LIQUID | OBSTACLE;
   WALBERLA_CHECK(classifyCube(cube) == BognerNeighbourhoodState::Gas);
}

// ---------------------------------------------------------------------------------------------
// isPass3Candidate() truth table: is a cell that fell back to Interface in Pass 2 re-examined in Pass 3, once all
// same-pass FormerMO cells have a final state? True for branches A'/C and for branch E with an unresolved neighbour;
// false for branch D (a genuine surface) and for branch E without unresolved neighbours (nothing can change).
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
// Two-round (Pass 2, then Pass 3) tests: classify the same cube twice, resolving the FormerMO neighbour in between,
// as Pass 3 sees a queued candidate in PdfReconstructionManager.h's decide sub-loop.
// ---------------------------------------------------------------------------------------------

void testCubeTwoRoundLiquidPromotion()
{
   // Pass 2: all-liquid except one unresolved FormerMO neighbour -> Interface (branch C), a Pass-3 candidate.
   Cube cube = uniformCube(LIQUID);
   cube[1][0][2] = FORMER_OBSTACLE; // (dx,dy,dz) = (+1,-1,0)
   bool hasGas, hasLiquid, hasUnknown;
   const auto round1 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round1 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(isPass3Candidate(round1, hasGas, hasLiquid, hasUnknown));

   // Pass 3: the neighbour resolved to LIQUID -> Liquid.
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
   // Pass 2: all-interface except one unresolved FormerMO neighbour -> Interface (branch E); a Pass-3 candidate
   // because hasUnknown is true.
   Cube cube = uniformCube(INTERFACE);
   cube[1][0][2] = FORMER_OBSTACLE;
   bool hasGas, hasLiquid, hasUnknown;
   const auto round1 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round1 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(isPass3Candidate(round1, hasGas, hasLiquid, hasUnknown));

   // Pass 3: the neighbour resolved to INTERFACE. The cell must stay Interface (an interface neighbour is no phase
   // evidence), and hasGas/hasLiquid must both be false: the [PASS3-UNRESOLVED-ISLAND] condition in
   // PdfReconstructionManager.h.
   cube[1][0][2] = INTERFACE;
   const auto round2 = classifyCubeFull(cube, hasGas, hasLiquid, hasUnknown);
   WALBERLA_CHECK(round2 == BognerNeighbourhoodState::Interface);
   WALBERLA_CHECK(!hasGas);
   WALBERLA_CHECK(!hasLiquid);
   WALBERLA_CHECK(!hasUnknown);
}

// ---------------------------------------------------------------------------------------------
// seedInterfaceFillFromNeighbourhood (PdfReconstructionManager.h): fill level seeded into a freshly uncovered cell
// that became INTERFACE. Legacy rule: mean fill of the interface neighbours, 0.5 if none. All-neighbour rule (prm
// seedFillFromAllNeighbours): mean volume fraction of all resolved neighbours (gas 0, liquid 1, interface fill),
// excluding obstacle and still-unclassified (formerObstacle) cells; 0.5 if none.
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
   // branch C: liquid all around, one unclassified neighbour; legacy falls back to 0.5, all-neighbour gives 1
   Cube c = uniformCube(LIQUID);
   c[1][1][2] = FORMER_OBSTACLE;
   checkSeed("bulk liquid + unknown (C)", c, uniformFill(real_t(1)), real_t(0.5), true, real_t(1), false);
}

void testSeedAirWithUnknownNeighbour()
{
   // branch A': gas all around, one unclassified neighbour; legacy falls back to 0.5, all-neighbour gives 0
   Cube c = uniformCube(GAS);
   c[1][1][0] = FORMER_OBSTACLE;
   checkSeed("air + unknown (A')", c, uniformFill(real_t(0)), real_t(0.5), true, real_t(0), false);
}

void testSeedSurfaceLiquidBelowGasAbove()
{
   // branch D without interface neighbours: lower plane liquid, upper plane gas, middle plane obstacle;
   // all-neighbour: 9 liquid / 18 resolved = 0.5
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
   // An interface neighbour that still carries the formerObstacle flag: skipped by default (prm
   // seedSkipObstacleNeighbours = true) -> no evidence, 0.5 fallback; counted with skipObstacle = false -> its fill.
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

// Strict Bogner rule (classifyStrictBognerNeighbourhood / classifyFormerObstacleCellStrict): Bogner's three cases in
// his order; cells freed in the same step (formerObstacle) are not in B(x).
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
   // dry upper face: gas above, the rest of the freed sheet beside it, body below -> gas (fallback rule: interface)
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
