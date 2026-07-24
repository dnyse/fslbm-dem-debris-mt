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

   return EXIT_SUCCESS;
}
