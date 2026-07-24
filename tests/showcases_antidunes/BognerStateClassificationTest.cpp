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
//! Unit tests for classifyFormerObstacleNeighbourhood() (PdfReconstructionManager.h), the pure
//! decision logic behind the Bogner & Rüde state-selection reconstruction's bognerTargetFlag.
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
//======================================================================================================================

#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <vector>

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
      LIQUID, GAS, FORMER_OBSTACLE, hasGas, hasLiquid, hasUnknown, unknownNbr);
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
// [WATERLINE-SHADOW] tests: computeRealPhaseComponents()/classifyViaWaterlineReachability() --
// DIAGNOSTIC-ONLY functions (see PdfReconstructionManager.h) that check whether a genuine real
// gas/liquid transition is reachable from an ambiguous cell via connectivity, not just "does the
// block contain any real gas/liquid anywhere." A permeable cell (formerObstacle, interface, or a
// particle-obstacle cell) never itself counts as evidence and never blocks the search -- only a
// real bognerGas cell blocks a real-liquid search, and vice versa. This is what makes a floating
// object sitting in open liquid NOT look like two separate liquid bodies, and a cluster of merely-
// uncovered cells NOT look like a waterline by itself.
//
// A flat 3D grid (row-major, indexed [z][y][x]) stands in for a FlagField region; getFlag(x,y,z)
// closes over it exactly like classifyCube's lambda above.
// ---------------------------------------------------------------------------------------------

constexpr MockFlag OBSTACLE = MockFlag(1) << 4; // a particle/wall cell -- mechanically just "not gas, not liquid"

using Grid = std::vector< std::vector< std::vector< MockFlag > > >; // [z][y][x]

Grid makeGrid(cell_idx_t xSize, cell_idx_t ySize, cell_idx_t zSize, MockFlag fill)
{
   return { size_t(zSize),
            std::vector< std::vector< MockFlag > >(size_t(ySize), std::vector< MockFlag >(size_t(xSize), fill)) };
}

PhaseReachability computePhases(const Grid& grid, cell_idx_t xSize, cell_idx_t ySize, cell_idx_t zSize)
{
   return computeRealPhaseComponents(
      [&](cell_idx_t x, cell_idx_t y, cell_idx_t z) { return grid[size_t(z)][size_t(y)][size_t(x)]; },
      cell_idx_t(0), xSize - cell_idx_t(1), cell_idx_t(0), ySize - cell_idx_t(1), cell_idx_t(0),
      zSize - cell_idx_t(1), LIQUID, GAS);
}

void testComponentsSingleLiquidSingleGasSimpleSurface()
{
   // z=0..2 real liquid, z=3 a single formerObstacle "surface" cell, z=4..6 real gas. The surface
   // cell genuinely borders both real phases (through the permeable layer, not directly) -> a
   // true local surface.
   const cell_idx_t xSize = 3, ySize = 3, zSize = 7;
   Grid grid = makeGrid(xSize, ySize, zSize, LIQUID);
   for (cell_idx_t y = 0; y < ySize; ++y)
      for (cell_idx_t x = 0; x < xSize; ++x)
      {
         grid[3][size_t(y)][size_t(x)] = FORMER_OBSTACLE;
         grid[4][size_t(y)][size_t(x)] = GAS;
         grid[5][size_t(y)][size_t(x)] = GAS;
         grid[6][size_t(y)][size_t(x)] = GAS;
      }

   const PhaseReachability pr = computePhases(grid, xSize, ySize, zSize);
   WALBERLA_CHECK_EQUAL(pr.totalGasComponents, uint_t(1));
   WALBERLA_CHECK_EQUAL(pr.totalLiquidComponents, uint_t(1));
   WALBERLA_CHECK(pr.gasConfidentlyReachable(1, 1, 3));
   WALBERLA_CHECK(pr.liquidConfidentlyReachable(1, 1, 3));
   WALBERLA_CHECK(!pr.gasInconclusive(1, 1, 3));
   WALBERLA_CHECK(!pr.liquidInconclusive(1, 1, 3));
   WALBERLA_CHECK(classifyViaWaterlineReachability(true, false, true, false) ==
                  WaterlineVerdict::StillAmbiguousSurfaceNearby);
}

void testComponentsObstacleDoesNotSplitLiquid()
{
   // A full-cross-section obstacle "wall" sitting inside an otherwise all-liquid domain, with one
   // formerObstacle cell in the middle of the wall (the cell under test). Without the permeability
   // rule this wall would look like a hard separator between the liquid on either side; with it,
   // both sides are recognised as the SAME connected body of liquid -- the direct regression test
   // for "an obstacle interrupting the waterline must not look like a second surface."
   const cell_idx_t xSize = 5, ySize = 5, zSize = 5;
   Grid grid = makeGrid(xSize, ySize, zSize, LIQUID);
   for (cell_idx_t z = 0; z < zSize; ++z)
      for (cell_idx_t y = 0; y < ySize; ++y)
         grid[size_t(z)][size_t(y)][2] = OBSTACLE; // full y-z wall at x=2
   grid[2][2][2] = FORMER_OBSTACLE;                 // the cell under test, in the middle of the wall

   const PhaseReachability pr = computePhases(grid, xSize, ySize, zSize);
   // The core claim: still exactly ONE liquid component, despite the wall spanning the full
   // cross-section -- the two "sides" are the same body of water, connected through the wall.
   WALBERLA_CHECK_EQUAL(pr.totalLiquidComponents, uint_t(1));
   WALBERLA_CHECK(pr.liquidConfidentlyReachable(2, 2, 2));
   // No real gas anywhere in this fixture at all.
   WALBERLA_CHECK(!pr.gasConfidentlyReachable(2, 2, 2));
}

void testComponentsFormerObstacleClusterIsNotEvidence()
{
   // A 3x3 cluster of formerObstacle cells sitting in the middle of an otherwise all-liquid plane
   // (inset from the domain boundary so the cluster's own gas-graph component doesn't touch the
   // scanned region's edge -- note the z-dimension needs thickness >1 too, or every cell would
   // trivially "touch" the degenerate single-layer z boundary). No real gas anywhere. Confirms the
   // cluster itself -- however large -- contributes zero evidence: every cell in it still correctly
   // sees the surrounding real liquid.
   const cell_idx_t xSize = 5, ySize = 5, zSize = 3;
   Grid grid = makeGrid(xSize, ySize, zSize, LIQUID);
   for (cell_idx_t y = 1; y <= 3; ++y)
      for (cell_idx_t x = 1; x <= 3; ++x)
         grid[1][size_t(y)][size_t(x)] = FORMER_OBSTACLE;

   const PhaseReachability pr = computePhases(grid, xSize, ySize, zSize);
   WALBERLA_CHECK_EQUAL(pr.totalGasComponents, uint_t(0));
   WALBERLA_CHECK_EQUAL(pr.totalLiquidComponents, uint_t(1));
   for (cell_idx_t y = 1; y <= 3; ++y)
      for (cell_idx_t x = 1; x <= 3; ++x)
      {
         WALBERLA_CHECK(pr.liquidConfidentlyReachable(x, y, 1));
         WALBERLA_CHECK(!pr.gasConfidentlyReachable(x, y, 1));
         WALBERLA_CHECK(!pr.gasInconclusive(x, y, 1));
         WALBERLA_CHECK(classifyViaWaterlineReachability(false, false, true, false) ==
                        WaterlineVerdict::WouldBeLiquid);
      }
}

void testComponentsTwoSeparateGasRegionsBothTracked()
{
   // Two real-gas patches far enough apart (and separated by real liquid) to be genuinely
   // disconnected components -- the "main surface plus an unrelated trapped cavity" case. A single
   // formerObstacle test cell sits adjacent to only ONE of the two patches. Confirms reachability
   // is tracked per-component: the cell correctly reaches the nearby patch and is NOT penalised by
   // the existence of the other, unreachable one elsewhere in the same scan.
   const cell_idx_t xSize = 9, ySize = 9, zSize = 2;
   Grid grid = makeGrid(xSize, ySize, zSize, LIQUID);
   // Patch A (near the test cell): z=1, x=1..2, y=1..2.
   for (cell_idx_t y = 1; y <= 2; ++y)
      for (cell_idx_t x = 1; x <= 2; ++x)
         grid[1][size_t(y)][size_t(x)] = GAS;
   grid[0][2][2] = FORMER_OBSTACLE; // test cell: directly below/adjacent to patch A
   // Patch B (far away, unreachable from the test cell without crossing real liquid): z=1, x=6..7, y=6..7.
   for (cell_idx_t y = 6; y <= 7; ++y)
      for (cell_idx_t x = 6; x <= 7; ++x)
         grid[1][size_t(y)][size_t(x)] = GAS;

   const PhaseReachability pr = computePhases(grid, xSize, ySize, zSize);
   WALBERLA_CHECK_EQUAL(pr.totalGasComponents, uint_t(2)); // two genuinely separate real-gas regions
   WALBERLA_CHECK(pr.gasConfidentlyReachable(2, 2, 0));    // reaches patch A
   WALBERLA_CHECK(pr.liquidConfidentlyReachable(2, 2, 0)); // plenty of real liquid nearby too
   WALBERLA_CHECK(classifyViaWaterlineReachability(true, false, true, false) ==
                  WaterlineVerdict::StillAmbiguousSurfaceNearby);
}

void testComponentsNoEvidenceAnywhereIsInconclusiveNotConfident()
{
   // A domain with NO real gas or liquid anywhere at all (every cell formerObstacle). The single
   // all-permeable component necessarily spans -- and therefore touches -- the scanned region's own
   // outer boundary (there is nothing else to stop it), so this resolves to "inconclusive" rather
   // than a confident "no evidence of either phase": with truly zero real evidence anywhere in
   // reach, the honest answer is "can't rule out more just past what we can see," not "confirmed
   // empty." (Enclosure by one real phase always also proves reachability of that same phase --
   // see class docs -- so a confident-and-non-edge-touching empty verdict never arises from a
   // partial enclosure either; it takes an actual real cell of the corresponding phase to resolve
   // confidently.)
   const cell_idx_t xSize = 5, ySize = 5, zSize = 5;
   const Grid grid = makeGrid(xSize, ySize, zSize, FORMER_OBSTACLE);

   const PhaseReachability pr = computePhases(grid, xSize, ySize, zSize);
   WALBERLA_CHECK_EQUAL(pr.totalGasComponents, uint_t(0));
   WALBERLA_CHECK_EQUAL(pr.totalLiquidComponents, uint_t(0));
   WALBERLA_CHECK(pr.gasInconclusive(2, 2, 2));
   WALBERLA_CHECK(pr.liquidInconclusive(2, 2, 2));
   WALBERLA_CHECK(classifyViaWaterlineReachability(false, true, false, true) ==
                  WaterlineVerdict::ExcludedNearUnresolvedEdge);
}

void testComponentsEdgeTruncation()
{
   // A formerObstacle region that runs all the way to the scanned window's outer boundary without
   // ever finding real gas -- the "we can't see far enough to be sure" case. Real liquid is present
   // elsewhere in the fixture (so this is not simply the no-evidence-anywhere case above), but the
   // specific cell under test only has a permeable path that reaches the edge on the gas side.
   const cell_idx_t xSize = 5, ySize = 5, zSize = 5;
   Grid grid = makeGrid(xSize, ySize, zSize, LIQUID);
   // Punch a formerObstacle corridor from the centre out to the x=0 boundary, never touching real
   // gas (there is none in this fixture at all) -- the centre cell's gas-graph component therefore
   // reaches the edge without ever becoming productive.
   for (cell_idx_t x = 0; x <= 2; ++x)
      grid[2][2][size_t(x)] = FORMER_OBSTACLE;

   const PhaseReachability pr = computePhases(grid, xSize, ySize, zSize);
   WALBERLA_CHECK(pr.gasInconclusive(2, 2, 2));
   WALBERLA_CHECK(!pr.gasConfidentlyReachable(2, 2, 2));
   // Liquid is directly adjacent (e.g. (2,1,2) or (1,2,2) are real liquid) and not edge-limited.
   WALBERLA_CHECK(pr.liquidConfidentlyReachable(2, 2, 2));
   WALBERLA_CHECK(!pr.liquidInconclusive(2, 2, 2));
   WALBERLA_CHECK(classifyViaWaterlineReachability(false, true, true, false) ==
                  WaterlineVerdict::ExcludedNearUnresolvedEdge);
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

   testComponentsSingleLiquidSingleGasSimpleSurface();
   testComponentsObstacleDoesNotSplitLiquid();
   testComponentsFormerObstacleClusterIsNotEvidence();
   testComponentsTwoSeparateGasRegionsBothTracked();
   testComponentsNoEvidenceAnywhereIsInconclusiveNotConfident();
   testComponentsEdgeTruncation();

   return EXIT_SUCCESS;
}
