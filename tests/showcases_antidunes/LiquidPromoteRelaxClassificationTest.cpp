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
//! \file LiquidPromoteRelaxClassificationTest.cpp
//! Full-grid stress test of ModelParameters.strictBognerRelaxInterface (the Round 0/1/2 liquid-promote-then-relax
//! classification in PdfReconstructionManager.h). Like BognerStateClassificationTest.cpp it calls the pure decision
//! functions (classifyFormerObstacleCellStrict, classifyFormerObstacleCell, isPass3Candidate), here over hand-built 3D
//! flag grids in which a sphere, cylinder or cuboid floater moves and releases cells. Floater sizes follow the prm
//! files: sphere radius 5 LU (prms/clog_centred_pillar_startblock_cappid_gas.prm; also 4 LU and random 3-5 LU
//! variants), cylinder radius 3 / length 18 LU
//! (prms/hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm), box half-extents <9,18,9> LU
//! (prms/bogner3_r025_deep_resume_gas.prm). The ambient field is liquid below / interface / gas above, either a flat
//! waterline or a moving wave-front height field. Not covered: BoundaryHandling/StructuredBlockForest integration,
//! MPI block boundaries, ghost-layer communication timing and the real ParticleFunctions containment.
//!
//! Five classification strategies are run per scenario:
//!   (a) pure strict Bogner (Round 0)
//!   (b) + strictBognerNoLiquidAtGas (Round 0+1)
//!   (c) + strictBognerRelaxInterface (Round 0+1+2)
//!   (d) legacy non-strict Pass 2 + Pass-3 re-resolution
//!   (e) legacy Pass 2 only (usePass3Reresolution = false, the production default)
//! The only assertion: (c) never leaves a liquid cell touching a gas cell. The violation counts of the other strategies
//! are logged only.
//!
//! Also logged, for every strategy other than (c): the "out-of-place interface" count, i.e. released cells the
//! strategy classifies INTERFACE while (c) classifies the same cell LIQUID or GAS. It shows how much of a strategy's
//! safety comes from over-classifying cells as interface.
//======================================================================================================================

#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace walberla
{
namespace lbm_mesapd_coupling
{

using MockFlag                             = uint8_t;
constexpr MockFlag LIQUID                  = MockFlag(1) << 0;
constexpr MockFlag GAS                     = MockFlag(1) << 1;
constexpr MockFlag INTERFACE               = MockFlag(1) << 2;
constexpr MockFlag FORMER_OBSTACLE         = MockFlag(1) << 3;
constexpr MockFlag OBSTACLE                = MockFlag(1) << 4;

struct Vec3 { double x, y, z; };
struct IVec3 { int x, y, z; };

// Zero-evidence instrumentation: released cells whose real (non-same-step) neighbourhood has neither liquid nor gas
// (hasGas == hasLiquid == false in Round 0). classifyStrictBognerNeighbourhood (`if (!hasLiquid) return Gas; if
// (!hasGas) return Liquid; return Interface;`, Bogner's ordered rule) assigns them GAS rather than INTERFACE.
uint_t g_zeroEvidenceCount = 0;
uint_t g_zeroEvidenceMismatchVsRelax = 0; // Round-0 GAS differs from relax result
struct ZeroEvidenceExample { std::string scenario; IVec3 cell; MockFlag pureResult; MockFlag relaxResult; };
std::vector< ZeroEvidenceExample > g_zeroEvidenceExamples;
uint_t g_scenarioIndex = 0;
bool g_dumpedFullExample = false;
// Of the zero-evidence cells (set to GAS in Round 0), how many touch a LIQUID cell in Round 0's output, i.e. cause
// a liquid/gas contact.
uint_t g_zeroEvidenceCausesViolation = 0;

// ---------------------------------------------------------------------------------------------
// A plain 3D grid of MockFlags, large enough to hold the floater plus an ambient margin on every side.
// ---------------------------------------------------------------------------------------------
struct Grid
{
   int nx, ny, nz;
   std::vector< MockFlag > data;
   Grid(int nx_, int ny_, int nz_) : nx(nx_), ny(ny_), nz(nz_), data(size_t(nx_) * size_t(ny_) * size_t(nz_), MockFlag(0)) {}
   size_t idx(int x, int y, int z) const { return (size_t(z) * size_t(ny) + size_t(y)) * size_t(nx) + size_t(x); }
   bool inRange(int x, int y, int z) const { return x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz; }
   MockFlag& at(int x, int y, int z) { return data[idx(x, y, z)]; }
   MockFlag  get(int x, int y, int z) const { return data[idx(x, y, z)]; }
};

// Single flat ambient interface at z == waterlineZ, liquid below, gas above.
MockFlag ambientFlag(int z, int waterlineZ)
{
   if (z < waterlineZ) return LIQUID;
   if (z == waterlineZ) return INTERFACE;
   return GAS;
}

// ---------------------------------------------------------------------------------------------
// Wave-front height fields: the local waterline height is a continuous function of (x,y) (tilted, curved, or with
// several lobes), classified by ambientFlagWave below, to model a moving/deforming free-surface front.
// ---------------------------------------------------------------------------------------------

// A locally raised or lowered Gaussian lobe of the free-surface height field.
double gaussianLobe(double x, double y, double cx, double cy, double amp, double sigma)
{
   const double dx = x - cx, dy = y - cy;
   return amp * std::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
}

// Tilted-plane wave front, so a sweeping front intersects a body on more than one side at once. Returns the
// continuous height (see ambientFlagWave).
double tiltedWaterlineHeight(double x, double y, double wz0, double cx, double cy, double slopeX, double slopeY)
{
   return wz0 + slopeX * (x - cx) + slopeY * (y - cy);
}

// Curved (non-planar) wave front: a single Gaussian ridge/trough rather than a flat tilted plane.
double curvedWaterlineHeight(double x, double y, double wz0, double cx, double cy, double amp, double sigma)
{
   return wz0 + gaussianLobe(x, y, cx, cy, amp, sigma);
}

// Two-lobe wave front: a crest on one side of the body and a trough on the other in the same height field, so one
// step can raise the surface windward and lower it leeward (two disconnected release fronts).
double twoLobeWaterlineHeight(double x, double y, double wz0, double crestCx, double crestCy, double crestAmp,
                               double troughCx, double troughCy, double troughAmp, double sigma)
{
   return wz0 + gaussianLobe(x, y, crestCx, crestCy, crestAmp, sigma)
              + gaussianLobe(x, y, troughCx, troughCy, troughAmp, sigma);
}

// Interface half-width (LU) for the wave-front ambient classification only. A one-cell band suffices for the flat
// waterline, but on a tilted or curved surface diagonal neighbours in adjacent columns can then put ambient liquid
// directly against ambient gas, an artefact of the set-up rather than of the engines under test. For all slopes used
// below, 1.5 LU leaves no such ambient-only violations.
constexpr double WAVE_INTERFACE_HALF_WIDTH = 1.5;

// Wave-front ambient classification from the cell centre's continuous height relative to the local surface, so the
// interface band stays a connected sheet that separates liquid from gas (see WAVE_INTERFACE_HALF_WIDTH).
template< typename HeightFn >
MockFlag ambientFlagWave(double x, double y, int z, HeightFn heightAt)
{
   const double zc = z + 0.5;
   const double h  = heightAt(x, y);
   if (zc < h - WAVE_INTERFACE_HALF_WIDTH) return LIQUID;
   if (zc > h + WAVE_INTERFACE_HALF_WIDTH) return GAS;
   return INTERFACE;
}

// ---------------------------------------------------------------------------------------------
// Floater shapes with cell-centre-inside-shape membership (the MEM mapping rule, cf. src/mesa_pd/common/Contains.h:
// centre test, not voxel overlap).
// ---------------------------------------------------------------------------------------------
bool insideSphere(Vec3 c, Vec3 centre, double radius)
{
   const double dx = c.x - centre.x, dy = c.y - centre.y, dz = c.z - centre.z;
   return dx * dx + dy * dy + dz * dz < radius * radius;
}

// Cylinder, axis along +Y (FloatingCylinder axis <0,1,0> in the prms).
bool insideCylinderY(Vec3 c, Vec3 centre, double radius, double length)
{
   const double dx = c.x - centre.x, dz = c.z - centre.z;
   if (dx * dx + dz * dz >= radius * radius) return false;
   return std::abs(c.y - centre.y) < length * 0.5;
}

// Box rotated by `rotDeg` about the Y axis, enough to exercise non-axis-aligned releases.
bool insideBoxRotY(Vec3 c, Vec3 centre, Vec3 halfExtents, double rotDeg)
{
   const double rad = rotDeg * 3.14159265358979323846 / 180.0;
   const double cs = std::cos(-rad), sn = std::sin(-rad);
   const double lx = c.x - centre.x, ly = c.y - centre.y, lz = c.z - centre.z;
   const double rx = cs * lx + sn * lz;
   const double rz = -sn * lx + cs * lz;
   return std::abs(rx) < halfExtents.x && std::abs(ly) < halfExtents.y && std::abs(rz) < halfExtents.z;
}

// ---------------------------------------------------------------------------------------------
// Decision engines built on the pure functions of PdfReconstructionManager.h (classifyFormerObstacleCellStrict /
// classifyFormerObstacleCell / isPass3Candidate); only the decide-then-apply orchestration around them is
// re-implemented here, following PdfReconstructionManager.h.
// ---------------------------------------------------------------------------------------------

const char* flagName(MockFlag f)
{
   if (f == LIQUID) return "LIQUID";
   if (f == GAS) return "GAS";
   if (f == INTERFACE) return "INTERFACE";
   if (f == FORMER_OBSTACLE) return "FORMER_OBSTACLE(same-step)";
   if (f == OBSTACLE) return "OBSTACLE";
   return "NONE(0)";
}

// Logs the 26 neighbour flags of a released cell as roundZeroPureStrictBogner's `nbr` lambda sees them (same-step or
// out-of-range -> FORMER_OBSTACLE, else the ambient flag).
void dumpNeighbourEvidence(const char* label, const Grid& ambientPlusObstacle, const std::vector< bool >& isUncovered,
                            const IVec3& c)
{
   WALBERLA_LOG_INFO_ON_ROOT(label << ": neighbour evidence for released cell (" << c.x << "," << c.y << "," << c.z << ")");
   for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
         for (int dx = -1; dx <= 1; ++dx)
         {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
            MockFlag f;
            if (!ambientPlusObstacle.inRange(nx, ny, nz)) { f = FORMER_OBSTACLE; }
            else
            {
               const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
               f = isUncovered[ni] ? FORMER_OBSTACLE : ambientPlusObstacle.get(nx, ny, nz);
            }
            WALBERLA_LOG_INFO_ON_ROOT("    offset(" << dx << "," << dy << "," << dz << ") -> " << flagName(f));
         }
}

// Round 0: pure strict Bogner, decide-then-apply over every newly uncovered cell. Same-step neighbours are reported
// as FORMER_OBSTACLE, which classifyFormerObstacleCellStrict excludes from B(x), as in Bogner's rule.
std::vector< MockFlag > roundZeroPureStrictBogner(const Grid& ambientPlusObstacle, const std::vector< IVec3 >& uncovered,
                                                   std::vector< IVec3 >* zeroEvidenceCellsOut = nullptr)
{
   std::vector< MockFlag > decided(ambientPlusObstacle.data.size(), MockFlag(0));
   std::vector< bool > isUncovered(ambientPlusObstacle.data.size(), false);
   for (auto& c : uncovered) isUncovered[ambientPlusObstacle.idx(c.x, c.y, c.z)] = true;

   for (auto& c : uncovered)
   {
      bool hasGas, hasLiquid;
      auto nbr = [&](int dx, int dy, int dz) -> MockFlag {
         const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
         if (!ambientPlusObstacle.inRange(nx, ny, nz)) return FORMER_OBSTACLE;
         const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
         if (isUncovered[ni]) return FORMER_OBSTACLE; // same-step, not yet decided -- not in B(x), Bogner-strict
         return ambientPlusObstacle.get(nx, ny, nz);
      };
      const auto state = classifyFormerObstacleCellStrict< MockFlag >(nbr, LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas, hasLiquid);
      // Zero-evidence case: B(x) holds neither a liquid nor a gas neighbour. Recorded for the instrumentation.
      if (!hasGas && !hasLiquid && zeroEvidenceCellsOut != nullptr) zeroEvidenceCellsOut->push_back(c);
      decided[ambientPlusObstacle.idx(c.x, c.y, c.z)] =
         (state == BognerNeighbourhoodState::Gas) ? GAS : (state == BognerNeighbourhoodState::Liquid) ? LIQUID : INTERFACE;
   }
   return decided;
}

// Round 1: turn any Round-0 LIQUID cell touching GAS into INTERFACE (OBSTACLE neighbours do not count). Mirrors the
// strictBognerNoLiquidAtGas loop in PdfReconstructionManager.h.
void roundOnePromoteLiquidTouchingGas(const Grid& ambientPlusObstacle, const std::vector< IVec3 >& uncovered,
                                       std::vector< MockFlag >& decided)
{
   std::vector< MockFlag > frozen = decided;
   for (auto& c : uncovered)
   {
      const size_t ci = ambientPlusObstacle.idx(c.x, c.y, c.z);
      if (frozen[ci] != LIQUID) continue;
      bool touches = false;
      for (int dz = -1; dz <= 1 && !touches; ++dz)
         for (int dy = -1; dy <= 1 && !touches; ++dy)
            for (int dx = -1; dx <= 1 && !touches; ++dx)
            {
               if (dx == 0 && dy == 0 && dz == 0) continue;
               const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
               if (!ambientPlusObstacle.inRange(nx, ny, nz)) continue;
               const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
               const MockFlag nf = (frozen[ni] != MockFlag(0)) ? frozen[ni] : ambientPlusObstacle.get(nx, ny, nz);
               touches = (nf == GAS);
            }
      if (touches) decided[ci] = INTERFACE;
   }
}

// Round 2: relax an INTERFACE cell to GAS/LIQUID if its post-Round-1 neighbourhood shows only one phase. Mirrors
// the strictBognerRelaxInterface pass.
void roundTwoRelaxInterface(const Grid& ambientPlusObstacle, const std::vector< IVec3 >& uncovered,
                             std::vector< MockFlag >& decided, uint_t& relaxedToGas, uint_t& relaxedToLiquid)
{
   std::vector< MockFlag > frozen = decided;
   std::vector< MockFlag > result = decided;
   for (auto& c : uncovered)
   {
      const size_t ci = ambientPlusObstacle.idx(c.x, c.y, c.z);
      if (frozen[ci] != INTERFACE) continue;
      bool hasGas = false, hasLiquid = false;
      for (int dz = -1; dz <= 1; ++dz)
         for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
            {
               if (dx == 0 && dy == 0 && dz == 0) continue;
               const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
               if (!ambientPlusObstacle.inRange(nx, ny, nz)) continue;
               const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
               const MockFlag nf = (frozen[ni] != MockFlag(0)) ? frozen[ni] : ambientPlusObstacle.get(nx, ny, nz);
               if (nf == OBSTACLE) continue;
               if (nf == LIQUID) hasLiquid = true;
               if (nf == GAS) hasGas = true;
            }
      if (hasGas && !hasLiquid) { result[ci] = GAS; ++relaxedToGas; }
      else if (hasLiquid && !hasGas) { result[ci] = LIQUID; ++relaxedToLiquid; }
   }
   decided = result;
}

// Legacy non-strict Pass 2 (sequential, z outer / y / x inner: the non-OpenMP order of
// src/field/iterators/IteratorMacros.h), then Pass 3 (decide-then-apply against the frozen post-Pass-2 snapshot).
std::vector< MockFlag > legacyPass2Pass3(const Grid& ambientPlusObstacle, std::vector< IVec3 > uncovered)
{
   std::sort(uncovered.begin(), uncovered.end(), [](const IVec3& a, const IVec3& b) {
      return std::tie(a.z, a.y, a.x) < std::tie(b.z, b.y, b.x);
   });
   std::vector< MockFlag > decided(ambientPlusObstacle.data.size(), MockFlag(0));
   std::vector< bool > isUncovered(ambientPlusObstacle.data.size(), false);
   for (auto& c : uncovered) isUncovered[ambientPlusObstacle.idx(c.x, c.y, c.z)] = true;
   std::vector< bool > isDecided(ambientPlusObstacle.data.size(), false);

   struct Pass3Cand { IVec3 c; bool hasGas, hasLiquid, hasUnknown; };
   std::vector< Pass3Cand > candidates;

   for (auto& c : uncovered)
   {
      bool hasGas, hasLiquid, hasUnknown;
      uint_t unknownNbr;
      auto nbr = [&](int dx, int dy, int dz) -> MockFlag {
         const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
         if (!ambientPlusObstacle.inRange(nx, ny, nz)) return FORMER_OBSTACLE;
         const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
         if (isUncovered[ni]) return isDecided[ni] ? decided[ni] : FORMER_OBSTACLE;
         return ambientPlusObstacle.get(nx, ny, nz);
      };
      const auto state = classifyFormerObstacleCell< MockFlag >(nbr, LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas,
                                                                  hasLiquid, hasUnknown, unknownNbr);
      const size_t ci = ambientPlusObstacle.idx(c.x, c.y, c.z);
      decided[ci] = (state == BognerNeighbourhoodState::Gas) ? GAS : (state == BognerNeighbourhoodState::Liquid) ? LIQUID : INTERFACE;
      isDecided[ci] = true;
      if (isPass3Candidate(state, hasGas, hasLiquid, hasUnknown)) candidates.push_back({ c, hasGas, hasLiquid, hasUnknown });
   }

   // Pass 3: re-run the same classifier against the frozen post-Pass-2 snapshot, decide-then-apply.
   std::vector< MockFlag > frozen = decided;
   std::vector< MockFlag > afterPass3 = decided;
   for (auto& cand : candidates)
   {
      bool hasGas, hasLiquid, hasUnknown;
      uint_t unknownNbr;
      auto nbr = [&](int dx, int dy, int dz) -> MockFlag {
         const int nx = cand.c.x + dx, ny = cand.c.y + dy, nz = cand.c.z + dz;
         if (!ambientPlusObstacle.inRange(nx, ny, nz)) return FORMER_OBSTACLE;
         const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
         if (isUncovered[ni]) return frozen[ni]; // always a final Pass-2 flag by now
         return ambientPlusObstacle.get(nx, ny, nz);
      };
      const auto state = classifyFormerObstacleCell< MockFlag >(nbr, LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas,
                                                                  hasLiquid, hasUnknown, unknownNbr);
      if (state == BognerNeighbourhoodState::Gas) afterPass3[ambientPlusObstacle.idx(cand.c.x, cand.c.y, cand.c.z)] = GAS;
      else if (state == BognerNeighbourhoodState::Liquid) afterPass3[ambientPlusObstacle.idx(cand.c.x, cand.c.y, cand.c.z)] = LIQUID;
   }
   return afterPass3;
}

// Legacy non-strict Pass 2 only: as legacyPass2Pass3 but without Pass 3, i.e. usePass3Reresolution = false, the
// production default in PdfReconstructionManager.h.
std::vector< MockFlag > legacyPass2Only(const Grid& ambientPlusObstacle, std::vector< IVec3 > uncovered)
{
   std::sort(uncovered.begin(), uncovered.end(), [](const IVec3& a, const IVec3& b) {
      return std::tie(a.z, a.y, a.x) < std::tie(b.z, b.y, b.x);
   });
   std::vector< MockFlag > decided(ambientPlusObstacle.data.size(), MockFlag(0));
   std::vector< bool > isUncovered(ambientPlusObstacle.data.size(), false);
   for (auto& c : uncovered) isUncovered[ambientPlusObstacle.idx(c.x, c.y, c.z)] = true;
   std::vector< bool > isDecided(ambientPlusObstacle.data.size(), false);

   for (auto& c : uncovered)
   {
      bool hasGas, hasLiquid, hasUnknown;
      uint_t unknownNbr;
      auto nbr = [&](int dx, int dy, int dz) -> MockFlag {
         const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
         if (!ambientPlusObstacle.inRange(nx, ny, nz)) return FORMER_OBSTACLE;
         const size_t ni = ambientPlusObstacle.idx(nx, ny, nz);
         if (isUncovered[ni]) return isDecided[ni] ? decided[ni] : FORMER_OBSTACLE;
         return ambientPlusObstacle.get(nx, ny, nz);
      };
      const auto state = classifyFormerObstacleCell< MockFlag >(nbr, LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas,
                                                                  hasLiquid, hasUnknown, unknownNbr);
      const size_t ci = ambientPlusObstacle.idx(c.x, c.y, c.z);
      decided[ci] = (state == BognerNeighbourhoodState::Gas) ? GAS : (state == BognerNeighbourhoodState::Liquid) ? LIQUID : INTERFACE;
      isDecided[ci] = true;
      // No Pass-3 loop: the only difference from legacyPass2Pass3.
   }
   return decided;
}

// Count liquid-touching-gas violations (26-neighbourhood) over the whole grid, using `decided` for cells released
// this step and the grid's ambient/obstacle flags elsewhere. Counts (liquid cell, gas neighbour) pairs.
uint_t countViolations(const Grid& ambientPlusObstacle, const std::vector< MockFlag >& decided)
{
   auto flagAt = [&](int x, int y, int z) -> MockFlag {
      const MockFlag d = decided[ambientPlusObstacle.idx(x, y, z)];
      return (d != MockFlag(0)) ? d : ambientPlusObstacle.get(x, y, z);
   };
   uint_t violations = 0;
   for (int z = 0; z < ambientPlusObstacle.nz; ++z)
      for (int y = 0; y < ambientPlusObstacle.ny; ++y)
         for (int x = 0; x < ambientPlusObstacle.nx; ++x)
         {
            if (flagAt(x, y, z) != LIQUID) continue;
            for (int dz = -1; dz <= 1; ++dz)
               for (int dy = -1; dy <= 1; ++dy)
                  for (int dx = -1; dx <= 1; ++dx)
                  {
                     if (dx == 0 && dy == 0 && dz == 0) continue;
                     const int nx = x + dx, ny = y + dy, nz = z + dz;
                     if (!ambientPlusObstacle.inRange(nx, ny, nz)) continue;
                     if (flagAt(nx, ny, nz) == GAS) ++violations;
                  }
         }
   return violations;
}

// Does cell c have a LIQUID neighbour (same flagAt convention as countViolations)?
bool touchesLiquid(const Grid& ambientPlusObstacle, const std::vector< MockFlag >& decided, const IVec3& c)
{
   auto flagAt = [&](int x, int y, int z) -> MockFlag {
      const MockFlag d = decided[ambientPlusObstacle.idx(x, y, z)];
      return (d != MockFlag(0)) ? d : ambientPlusObstacle.get(x, y, z);
   };
   for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
         for (int dx = -1; dx <= 1; ++dx)
         {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const int nx = c.x + dx, ny = c.y + dy, nz = c.z + dz;
            if (!ambientPlusObstacle.inRange(nx, ny, nz)) continue;
            if (flagAt(nx, ny, nz) == LIQUID) return true;
         }
   return false;
}

uint_t countInterface(const std::vector< IVec3 >& uncovered, const Grid& g, const std::vector< MockFlag >& decided)
{
   uint_t n = 0;
   for (auto& c : uncovered)
      if (decided[g.idx(c.x, c.y, c.z)] == INTERFACE) ++n;
   return n;
}

// Out-of-place interface: released cells that `strategyDecided` classifies INTERFACE while `groundTruthDecided`
// (the strictBognerRelaxInterface result of the same scenario) classifies them LIQUID or GAS. Measured, not asserted.
uint_t countOutOfPlaceInterface(const std::vector< IVec3 >& uncovered, const Grid& g,
                                 const std::vector< MockFlag >& strategyDecided,
                                 const std::vector< MockFlag >& groundTruthDecided)
{
   uint_t n = 0;
   for (auto& c : uncovered)
   {
      const size_t ci = g.idx(c.x, c.y, c.z);
      if (strategyDecided[ci] != INTERFACE) continue;
      const MockFlag truth = groundTruthDecided[ci];
      if (truth == LIQUID || truth == GAS) ++n;
   }
   return n;
}

// ---------------------------------------------------------------------------------------------
// Scenario harness: build before/after obstacle membership for a shape and move, derive the newly uncovered cells,
// run all five strategies, report.
// ---------------------------------------------------------------------------------------------
struct ScenarioResult { uint_t nUncovered, pureViolations, noLiquidAtGasViolations, relaxViolations, legacyViolations,
                        pass2OnlyViolations, pureInterface, relaxInterface, legacyInterface, pass2OnlyInterface,
                        // Out-of-place interface counts (not computed for strictBognerRelaxInterface, the
                        // reference).
                        pureOutOfPlaceInterface, noLiquidAtGasOutOfPlaceInterface, legacyOutOfPlaceInterface,
                        pass2OnlyOutOfPlaceInterface,
                        // Zero-evidence cells in Round 0, and how many of them end in a different phase than
                        // under strictBognerRelaxInterface. They are set to GAS, never INTERFACE, so they cannot
                        // appear in pureOutOfPlaceInterface.
                        nZeroEvidence, nZeroEvidenceMismatch; };

template< typename BeforeFn, typename AfterFn >
ScenarioResult runScenario(int nx, int ny, int nz, int waterlineZ, BeforeFn insideBefore, AfterFn insideAfter)
{
   ++g_scenarioIndex;
   Grid grid(nx, ny, nz);
   std::vector< bool > before(size_t(nx) * ny * nz), after(size_t(nx) * ny * nz);
   for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
         for (int x = 0; x < nx; ++x)
         {
            const Vec3 c{ x + 0.5, y + 0.5, z + 0.5 };
            before[grid.idx(x, y, z)] = insideBefore(c);
            after[grid.idx(x, y, z)]  = insideAfter(c);
            grid.at(x, y, z) = before[grid.idx(x, y, z)] ? OBSTACLE : ambientFlag(z, waterlineZ);
         }
   // grid: BEFORE shape as OBSTACLE, ambient elsewhere. afterGrid: the same for the AFTER shape.
   Grid afterGrid(nx, ny, nz);
   std::vector< IVec3 > uncovered, covered;
   for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
         for (int x = 0; x < nx; ++x)
         {
            const bool b = before[grid.idx(x, y, z)], a = after[grid.idx(x, y, z)];
            afterGrid.at(x, y, z) = a ? OBSTACLE : ambientFlag(z, waterlineZ);
            if (b && !a) uncovered.push_back({ x, y, z });
            if (!b && a) covered.push_back({ x, y, z });
         }
   WALBERLA_UNUSED(covered);

   // Evidence grid read by every engine: OBSTACLE where the AFTER shape covers the cell, ambient elsewhere (the
   // ambient value of the released cells themselves is never read).
   Grid evidence(nx, ny, nz);
   for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
         for (int x = 0; x < nx; ++x)
         {
            const bool a = after[grid.idx(x, y, z)];
            evidence.at(x, y, z) = a ? OBSTACLE : ambientFlag(z, waterlineZ);
         }

   ScenarioResult r{};
   r.nUncovered = uint_t(uncovered.size());
   if (uncovered.empty()) return r;

   std::vector< IVec3 > zeroEvidenceCells;
   auto pure = roundZeroPureStrictBogner(evidence, uncovered, &zeroEvidenceCells);
   r.pureViolations = countViolations(evidence, pure);
   r.pureInterface  = countInterface(uncovered, evidence, pure);

   auto withPromote = pure;
   roundOnePromoteLiquidTouchingGas(evidence, uncovered, withPromote);
   r.noLiquidAtGasViolations = countViolations(evidence, withPromote);

   auto withRelax = withPromote;
   uint_t toGas = 0, toLiquid = 0;
   roundTwoRelaxInterface(evidence, uncovered, withRelax, toGas, toLiquid);
   r.relaxViolations = countViolations(evidence, withRelax);
   r.relaxInterface  = countInterface(uncovered, evidence, withRelax);

   auto legacy = legacyPass2Pass3(evidence, uncovered);
   r.legacyViolations = countViolations(evidence, legacy);
   r.legacyInterface  = countInterface(uncovered, evidence, legacy);

   auto pass2Only = legacyPass2Only(evidence, uncovered);
   r.pass2OnlyViolations = countViolations(evidence, pass2Only);
   r.pass2OnlyInterface  = countInterface(uncovered, evidence, pass2Only);

   // Out-of-place interface, measured against the strictBognerRelaxInterface result (withRelax).
   r.pureOutOfPlaceInterface         = countOutOfPlaceInterface(uncovered, evidence, pure, withRelax);
   r.noLiquidAtGasOutOfPlaceInterface = countOutOfPlaceInterface(uncovered, evidence, withPromote, withRelax);
   r.legacyOutOfPlaceInterface       = countOutOfPlaceInterface(uncovered, evidence, legacy, withRelax);
   r.pass2OnlyOutOfPlaceInterface    = countOutOfPlaceInterface(uncovered, evidence, pass2Only, withRelax);

   // Zero-evidence cells in Round 0 and their mismatches with the strictBognerRelaxInterface result.
   r.nZeroEvidence = uint_t(zeroEvidenceCells.size());
   uint_t mismatch = 0;
   for (auto& c : zeroEvidenceCells)
   {
      const size_t ci     = evidence.idx(c.x, c.y, c.z);
      const MockFlag pureFlag  = pure[ci];
      const MockFlag relaxFlag = withRelax[ci];
      if (pureFlag != relaxFlag)
      {
         ++mismatch;
         if (g_zeroEvidenceExamples.size() < 10)
            g_zeroEvidenceExamples.push_back({ "rigid#" + std::to_string(g_scenarioIndex), c, pureFlag, relaxFlag });
      }
      if (touchesLiquid(evidence, pure, c)) ++g_zeroEvidenceCausesViolation;
      if (!g_dumpedFullExample)
      {
         g_dumpedFullExample = true;
         std::vector< bool > isUncovered(evidence.data.size(), false);
         for (auto& u : uncovered) isUncovered[evidence.idx(u.x, u.y, u.z)] = true;
         dumpNeighbourEvidence("FIRST zero-evidence cell found (rigid-body scenario)", evidence, isUncovered, c);
         WALBERLA_LOG_INFO_ON_ROOT("    pure (Round 0) result = " << flagName(pureFlag)
            << " | strictBognerRelaxInterface ground truth = " << flagName(relaxFlag));
      }
   }
   r.nZeroEvidenceMismatch = mismatch;
   g_zeroEvidenceCount += r.nZeroEvidence;
   g_zeroEvidenceMismatchVsRelax += mismatch;

   return r;
}

// ---------------------------------------------------------------------------------------------
// Wave-driven scenario harness: the same five-strategy comparison as runScenario, but the ambient field comes from a
// moving free-surface height field instead of a flat waterline; the body can still move between before and after.
// beforeWL/afterWL return the continuous local waterline height at a cell's (x,y) centre before/after the step,
// classified by ambientFlagWave.
// ---------------------------------------------------------------------------------------------
template< typename BeforeFn, typename AfterFn, typename BeforeWLFn, typename AfterWLFn >
ScenarioResult runScenarioWave(int nx, int ny, int nz, BeforeWLFn beforeWL, AfterWLFn afterWL,
                                BeforeFn insideBefore, AfterFn insideAfter)
{
   ++g_scenarioIndex;
   Grid grid(nx, ny, nz);
   std::vector< bool > before(size_t(nx) * ny * nz), after(size_t(nx) * ny * nz);
   for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
         for (int x = 0; x < nx; ++x)
         {
            const Vec3 c{ x + 0.5, y + 0.5, z + 0.5 };
            before[grid.idx(x, y, z)] = insideBefore(c);
            after[grid.idx(x, y, z)]  = insideAfter(c);
            grid.at(x, y, z) = before[grid.idx(x, y, z)] ? OBSTACLE : ambientFlagWave(c.x, c.y, z, beforeWL);
         }

   std::vector< IVec3 > uncovered, covered;
   for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
         for (int x = 0; x < nx; ++x)
         {
            const bool b = before[grid.idx(x, y, z)], a = after[grid.idx(x, y, z)];
            if (b && !a) uncovered.push_back({ x, y, z });
            if (!b && a) covered.push_back({ x, y, z });
         }
   WALBERLA_UNUSED(covered);

   // Evidence grid: the AFTER body shape as OBSTACLE, the AFTER wave field elsewhere.
   Grid evidence(nx, ny, nz);
   for (int z = 0; z < nz; ++z)
      for (int y = 0; y < ny; ++y)
         for (int x = 0; x < nx; ++x)
         {
            const Vec3 c{ x + 0.5, y + 0.5, z + 0.5 };
            const bool a = after[grid.idx(x, y, z)];
            evidence.at(x, y, z) = a ? OBSTACLE : ambientFlagWave(c.x, c.y, z, afterWL);
         }

   ScenarioResult r{};
   r.nUncovered = uint_t(uncovered.size());
   if (uncovered.empty()) return r;

   std::vector< IVec3 > zeroEvidenceCells;
   auto pure = roundZeroPureStrictBogner(evidence, uncovered, &zeroEvidenceCells);
   r.pureViolations = countViolations(evidence, pure);
   r.pureInterface  = countInterface(uncovered, evidence, pure);

   auto withPromote = pure;
   roundOnePromoteLiquidTouchingGas(evidence, uncovered, withPromote);
   r.noLiquidAtGasViolations = countViolations(evidence, withPromote);

   auto withRelax = withPromote;
   uint_t toGas = 0, toLiquid = 0;
   roundTwoRelaxInterface(evidence, uncovered, withRelax, toGas, toLiquid);
   r.relaxViolations = countViolations(evidence, withRelax);
   r.relaxInterface  = countInterface(uncovered, evidence, withRelax);

   auto legacy = legacyPass2Pass3(evidence, uncovered);
   r.legacyViolations = countViolations(evidence, legacy);
   r.legacyInterface  = countInterface(uncovered, evidence, legacy);

   auto pass2Only = legacyPass2Only(evidence, uncovered);
   r.pass2OnlyViolations = countViolations(evidence, pass2Only);
   r.pass2OnlyInterface  = countInterface(uncovered, evidence, pass2Only);

   // Out-of-place interface, measured against the strictBognerRelaxInterface result (withRelax).
   r.pureOutOfPlaceInterface         = countOutOfPlaceInterface(uncovered, evidence, pure, withRelax);
   r.noLiquidAtGasOutOfPlaceInterface = countOutOfPlaceInterface(uncovered, evidence, withPromote, withRelax);
   r.legacyOutOfPlaceInterface       = countOutOfPlaceInterface(uncovered, evidence, legacy, withRelax);
   r.pass2OnlyOutOfPlaceInterface    = countOutOfPlaceInterface(uncovered, evidence, pass2Only, withRelax);

   // Zero-evidence tracking as in runScenario.
   r.nZeroEvidence = uint_t(zeroEvidenceCells.size());
   uint_t mismatch = 0;
   for (auto& c : zeroEvidenceCells)
   {
      const size_t ci     = evidence.idx(c.x, c.y, c.z);
      const MockFlag pureFlag  = pure[ci];
      const MockFlag relaxFlag = withRelax[ci];
      if (pureFlag != relaxFlag)
      {
         ++mismatch;
         if (g_zeroEvidenceExamples.size() < 10)
            g_zeroEvidenceExamples.push_back({ "wave#" + std::to_string(g_scenarioIndex), c, pureFlag, relaxFlag });
      }
      if (touchesLiquid(evidence, pure, c)) ++g_zeroEvidenceCausesViolation;
   }
   r.nZeroEvidenceMismatch = mismatch;
   g_zeroEvidenceCount += r.nZeroEvidence;
   g_zeroEvidenceMismatchVsRelax += mismatch;

   return r;
}

void reportAndCheck(const char* name, const ScenarioResult& r, uint_t& totalScenarios, uint_t& pureBad, uint_t& legacyBad,
                     uint_t& pass2OnlyBad, uint_t& pureOutOfPlace, uint_t& noLiquidAtGasOutOfPlace,
                     uint_t& legacyOutOfPlace, uint_t& pass2OnlyOutOfPlace)
{
   if (r.nUncovered == 0) return;
   ++totalScenarios;
   if (r.pureViolations > 0) ++pureBad;
   if (r.legacyViolations > 0) ++legacyBad;
   if (r.pass2OnlyViolations > 0) ++pass2OnlyBad;
   pureOutOfPlace += r.pureOutOfPlaceInterface;
   noLiquidAtGasOutOfPlace += r.noLiquidAtGasOutOfPlaceInterface;
   legacyOutOfPlace += r.legacyOutOfPlaceInterface;
   pass2OnlyOutOfPlace += r.pass2OnlyOutOfPlaceInterface;
   WALBERLA_LOG_INFO_ON_ROOT(
      name << ": uncovered=" << r.nUncovered << " | pure-Bogner violations=" << r.pureViolations
           << " (interface=" << r.pureInterface << ", outOfPlaceInterface=" << r.pureOutOfPlaceInterface << ")"
           << " | +noLiquidAtGas violations=" << r.noLiquidAtGasViolations
           << " (outOfPlaceInterface=" << r.noLiquidAtGasOutOfPlaceInterface << ")"
           << " | +relax violations=" << r.relaxViolations << " (interface=" << r.relaxInterface << ")"
           << " | legacy Pass2+Pass3 violations=" << r.legacyViolations << " (interface=" << r.legacyInterface
           << ", outOfPlaceInterface=" << r.legacyOutOfPlaceInterface << ")"
           << " | legacy Pass2-only (no Pass3) violations=" << r.pass2OnlyViolations
           << " (interface=" << r.pass2OnlyInterface << ", outOfPlaceInterface=" << r.pass2OnlyOutOfPlaceInterface << ")"
           << " | pure zero-real-evidence cells=" << r.nZeroEvidence << " (mismatch vs relax ground truth="
           << r.nZeroEvidenceMismatch << ")");
   // The only assertion of this test:
   WALBERLA_CHECK_EQUAL(r.relaxViolations, uint_t(0), name << ": strictBognerRelaxInterface must NEVER leave a liquid"
                                                            " cell directly touching a gas cell");
}

}
}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::lbm_mesapd_coupling;

   walberla::uint_t totalScenarios = 0, pureBad = 0, legacyBad = 0, pass2OnlyBad = 0;
   // Out-of-place interface running totals (summed across every rigid-body scenario below).
   walberla::uint_t pureOOP = 0, noLiquidAtGasOOP = 0, legacyOOP = 0, pass2OnlyOOP = 0;

   // ---- Single sphere, radius 5 LU (prms/clog_centred_pillar_startblock_cappid_gas.prm), ----
   // translating through the waterline along axis-aligned and diagonal directions, 1-4 LU/step.
   {
      const int N = 28, WZ = 14;
      const Vec3 centre{ 14, 14, 14 };
      const double R = 5.0;
      std::vector< Vec3 > moves = { {1,0,0}, {0,1,0}, {0,0,1}, {-1,0,0}, {0,-1,0}, {0,0,-1},
                                    {1,1,0}, {1,0,1}, {1,1,1}, {2,0,0}, {3,0,0}, {4,0,0}, {2,2,1},
                                    {0,0,2}, {0,0,3}, {-1,-1,-1} };
      for (auto& m : moves)
      {
         auto r = runScenario(N, N, N, WZ,
            [&](Vec3 c) { return insideSphere(c, centre, R); },
            [&](Vec3 c) { return insideSphere(c, Vec3{centre.x + m.x, centre.y + m.y, centre.z + m.z}, R); });
         reportAndCheck("sphere R=5 translate", r, totalScenarios, pureBad, legacyBad, pass2OnlyBad, pureOOP, noLiquidAtGasOOP, legacyOOP, pass2OnlyOOP);
      }
   }

   // ---- Single sphere, radius 4 LU, at several vertical offsets relative to the waterline. ----
   {
      const int N = 24, WZ = 12;
      const double R = 4.0;
      for (double zoff : { -1.5, -0.5, 0.0, 0.5, 1.5 })
      {
         const Vec3 centre{ 12, 12, 12 + zoff };
         auto r = runScenario(N, N, N, WZ,
            [&](Vec3 c) { return insideSphere(c, centre, R); },
            [&](Vec3 c) { return insideSphere(c, Vec3{centre.x + 1, centre.y, centre.z}, R); });
         reportAndCheck("sphere R=4 vertical-offset", r, totalScenarios, pureBad, legacyBad, pass2OnlyBad, pureOOP, noLiquidAtGasOOP, legacyOOP, pass2OnlyOOP);
      }
   }

   // ---- Single cylinder, radius 3, length 18 LU, axis Y (prms/hopper_clog_cyl...). ----
   {
      const int N = 30, WZ = 15;
      const Vec3 centre{ 15, 15, 15 };
      std::vector< Vec3 > moves = { {1,0,0}, {0,0,1}, {1,0,1}, {2,0,0}, {0,0,2}, {-1,0,1} };
      for (auto& m : moves)
      {
         auto r = runScenario(N, N, N, WZ,
            [&](Vec3 c) { return insideCylinderY(c, centre, 3.0, 18.0); },
            [&](Vec3 c) { return insideCylinderY(c, Vec3{centre.x + m.x, centre.y + m.y, centre.z + m.z}, 3.0, 18.0); });
         reportAndCheck("cylinder R=3 L=18 translate", r, totalScenarios, pureBad, legacyBad, pass2OnlyBad, pureOOP, noLiquidAtGasOOP, legacyOOP, pass2OnlyOOP);
      }
   }

   // ---- Single cuboid, half-extents <9,18,9> LU (prms/bogner3_r025...), including rotation. ----
   {
      const int N = 44, WZ = 22;
      const Vec3 centre{ 22, 22, 22 };
      const Vec3 he{ 9, 18, 9 };
      struct Case { Vec3 move; double rotBefore, rotAfter; };
      std::vector< Case > cases = { {{1,0,0}, 0, 0}, {{0,0,1}, 0, 0}, {{1,0,1}, 0, 0},
                                     {{1,0,0}, 0, 15}, {{0,0,1}, 10, 25}, {{1,0,1}, 5, 5} };
      for (auto& cs : cases)
      {
         auto r = runScenario(N, N, N, WZ,
            [&](Vec3 c) { return insideBoxRotY(c, centre, he, cs.rotBefore); },
            [&](Vec3 c) {
               return insideBoxRotY(c, Vec3{centre.x + cs.move.x, centre.y + cs.move.y, centre.z + cs.move.z}, he,
                                     cs.rotAfter);
            });
         reportAndCheck("cuboid half<9,18,9> translate+rotate", r, totalScenarios, pureBad, legacyBad, pass2OnlyBad, pureOOP, noLiquidAtGasOOP, legacyOOP, pass2OnlyOOP);
      }
   }

   // ---- Two spheres closing in on each other at the waterline (cf. FloatingSphere entries 12 LU apart in ----
   // clog_centred_pillar_startblock_cappid_gas.prm): released cells of both spheres can be same-step neighbours.
   {
      const int N = 30, WZ = 15;
      const double R = 5.0;
      std::vector< double > gaps = { 11.0, 10.2, 9.5, 9.0 }; // 2R=10 -> gaps < 10 overlap before the move
      for (double gap : gaps)
      {
         const Vec3 c1{ 15 - gap * 0.5, 15, 15 };
         const Vec3 c2{ 15 + gap * 0.5, 15, 15 };
         auto r = runScenario(N, N, N, WZ,
            [&](Vec3 c) { return insideSphere(c, c1, R) || insideSphere(c, c2, R); },
            [&](Vec3 c) {
               return insideSphere(c, Vec3{c1.x + 1, c1.y, c1.z}, R) || insideSphere(c, Vec3{c2.x - 1, c2.y, c2.z}, R);
            });
         reportAndCheck("two spheres R=5 closing gap", r, totalScenarios, pureBad, legacyBad, pass2OnlyBad, pureOOP, noLiquidAtGasOOP, legacyOOP, pass2OnlyOOP);
      }
   }

   // ---- Randomised jam patch: four spheres of random size moving in random directions near the waterline. ----
   {
      std::mt19937 rng(1337);
      std::uniform_real_distribution< double > posXY(6.0, 34.0), zOff(-2.0, 2.0), radiusD(3.0, 5.0),
         moveD(-2.0, 2.0);
      for (int trial = 0; trial < 6; ++trial)
      {
         const int N = 40, WZ = 20;
         struct S { Vec3 c0, c1; double r; };
         std::vector< S > spheres;
         for (int i = 0; i < 4; ++i)
         {
            Vec3 c0{ posXY(rng), posXY(rng), WZ + zOff(rng) };
            Vec3 c1{ c0.x + moveD(rng), c0.y + moveD(rng), c0.z + moveD(rng) * 0.3 };
            spheres.push_back({ c0, c1, radiusD(rng) });
         }
         auto r = runScenario(N, N, N, WZ,
            [&](Vec3 c) {
               for (auto& s : spheres) if (insideSphere(c, s.c0, s.r)) return true;
               return false;
            },
            [&](Vec3 c) {
               for (auto& s : spheres) if (insideSphere(c, s.c1, s.r)) return true;
               return false;
            });
         reportAndCheck("randomised 4-sphere jam patch", r, totalScenarios, pureBad, legacyBad, pass2OnlyBad, pureOOP, noLiquidAtGasOOP, legacyOOP, pass2OnlyOOP);
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("=== SUMMARY: " << totalScenarios << " scenarios with >=1 uncovered cell | "
      << "pure strict Bogner violated the invariant in " << pureBad << "/" << totalScenarios << " | "
      << "legacy Pass-2+Pass-3 violated the invariant in " << legacyBad << "/" << totalScenarios << " | "
      << "legacy Pass-2-only (no Pass-3) violated the invariant in " << pass2OnlyBad << "/" << totalScenarios << " | "
      << "strictBognerRelaxInterface violated it in 0/" << totalScenarios << " (every case asserted above) ===");
   WALBERLA_LOG_INFO_ON_ROOT("=== SUMMARY (out-of-place interface, summed over " << totalScenarios
      << " rigid-body scenarios): pure strict Bogner=" << pureOOP << " | +noLiquidAtGas=" << noLiquidAtGasOOP
      << " | legacy Pass-2+Pass-3=" << legacyOOP << " | legacy Pass-2-only (no Pass-3)=" << pass2OnlyOOP
      << " cells classified INTERFACE that strictBognerRelaxInterface itself resolved to LIQUID/GAS ===");

   // =============================================================================================
   // WAVE-FRONT scenarios: the ambient field is driven by a moving free-surface front (tilted, curved,
   // crest-then-trough, two-lobe) while the body moves slightly. Separate counters, plus a combined total.
   // =============================================================================================
   walberla::uint_t waveScenarios = 0, pureBadWave = 0, legacyBadWave = 0, pass2OnlyBadWave = 0;
   // Out-of-place interface running totals (summed across every wave-driven scenario below).
   walberla::uint_t pureOOPWave = 0, noLiquidAtGasOOPWave = 0, legacyOOPWave = 0, pass2OnlyOOPWave = 0;

   // ---- Tilted wave front sweeping across sphere R=5, intersecting it on several sides; body moves slightly. ----
   {
      const int N = 28, WZ = 14;
      const Vec3 centre{ 14, 14, 14 };
      const double R = 5.0;
      const Vec3 smallMove{ 0.4, 0.3, 0.2 };
      struct Tilt { double slopeX, slopeY, sweep; };
      std::vector< Tilt > tilts = { {0.30, 0.00, 1.0}, {0.00, 0.30, 1.0}, {0.25, 0.25, 1.0},
                                     {0.40, -0.20, 2.0}, {-0.30, 0.30, 1.0} };
      for (auto& t : tilts)
      {
         auto beforeWL = [&](double x, double y) { return tiltedWaterlineHeight(x, y, WZ, centre.x, centre.y, t.slopeX, t.slopeY); };
         auto afterWL  = [&](double x, double y) { return tiltedWaterlineHeight(x, y, WZ + t.sweep, centre.x, centre.y, t.slopeX, t.slopeY); };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideSphere(c, centre, R); },
            [&](Vec3 c) { return insideSphere(c, Vec3{centre.x + smallMove.x, centre.y + smallMove.y, centre.z + smallMove.z}, R); });
         reportAndCheck("wave: tilted front sweeping across sphere R=5", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
   }

   // ---- Curved (Gaussian-ridge) wave front across cylinder R=3 L=18, amplitude changing between before/after. ----
   {
      const int N = 30, WZ = 15;
      const Vec3 centre{ 15, 15, 15 };
      const Vec3 smallMove{ 0.3, 0.2, 0.0 };
      struct Curve { double cx, cy, amp, sigma, ampAfter; };
      std::vector< Curve > curves = { {8, 15, 2.0, 6.0, 3.0}, {15, 15, -1.5, 5.0, -2.5},
                                       {22, 15, 2.5, 7.0, 1.0}, {15, 8, 1.8, 6.0, 2.8} };
      for (auto& cv : curves)
      {
         auto beforeWL = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, cv.cx, cv.cy, cv.amp, cv.sigma); };
         auto afterWL  = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, cv.cx, cv.cy, cv.ampAfter, cv.sigma); };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideCylinderY(c, centre, 3.0, 18.0); },
            [&](Vec3 c) { return insideCylinderY(c, Vec3{centre.x + smallMove.x, centre.y + smallMove.y, centre.z + smallMove.z}, 3.0, 18.0); });
         reportAndCheck("wave: curved ridge front sweeping across cylinder R=3 L=18", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
   }

   // ---- A crest advances, then in the next step a trough recedes over the same side of the body ----
   // (covered, then uncovered, over two consecutive steps).
   // Sphere case:
   {
      const int N = 28, WZ = 14;
      const Vec3 centre{ 14, 14, 14 };
      const double R = 5.0;
      Vec3 centreStep1{ centre.x + 0.5, centre.y, centre.z };
      {
         auto beforeWL = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x + R, centre.y, 0.5, 5.0); };
         auto afterWL  = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x + R, centre.y, 3.0, 5.0); };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideSphere(c, centre, R); },
            [&](Vec3 c) { return insideSphere(c, centreStep1, R); });
         reportAndCheck("wave: crest advancing over sphere R=5 (step 1 of 2)", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
      {
         Vec3 centreStep2{ centreStep1.x + 0.5, centreStep1.y, centreStep1.z };
         auto beforeWL = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x + R, centre.y, 3.0, 5.0); };
         auto afterWL  = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x + R, centre.y, -3.0, 5.0); };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideSphere(c, centreStep1, R); },
            [&](Vec3 c) { return insideSphere(c, centreStep2, R); });
         reportAndCheck("wave: trough receding over sphere R=5, same side (step 2 of 2)", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
   }
   // Rotating cuboid case:
   {
      const int N = 44, WZ = 22;
      const Vec3 centre{ 22, 22, 22 };
      const Vec3 he{ 9, 18, 9 };
      Vec3 centreStep1{ centre.x, centre.y, centre.z + 0.4 };
      {
         auto beforeWL = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x - he.x, centre.y, 1.0, 10.0); };
         auto afterWL  = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x - he.x, centre.y, 4.0, 10.0); };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideBoxRotY(c, centre, he, 0.0); },
            [&](Vec3 c) { return insideBoxRotY(c, centreStep1, he, 6.0); });
         reportAndCheck("wave: crest advancing over rotating cuboid (step 1 of 2)", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
      {
         Vec3 centreStep2{ centreStep1.x, centreStep1.y, centreStep1.z + 0.4 };
         auto beforeWL = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x - he.x, centre.y, 4.0, 10.0); };
         auto afterWL  = [&](double x, double y) { return curvedWaterlineHeight(x, y, WZ, centre.x - he.x, centre.y, -4.0, 10.0); };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideBoxRotY(c, centreStep1, he, 6.0); },
            [&](Vec3 c) { return insideBoxRotY(c, centreStep2, he, 10.0); });
         reportAndCheck("wave: trough receding over rotating cuboid, same side (step 2 of 2)", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
   }

   // ---- Disconnected release patches on one body in one step: crest advancing windward, trough ----
   // receding leeward.
   // Sphere case:
   {
      const int N = 30, WZ = 15;
      const Vec3 centre{ 15, 15, 15 };
      const double R = 5.0;
      const Vec3 smallMove{ 0.3, 0.0, 0.1 };
      struct Lobes { double crestCx, crestCy, crestAmp, troughCx, troughCy, troughAmp, sigma; };
      std::vector< Lobes > cases = {
         { centre.x - R, centre.y, 2.5, centre.x + R, centre.y, -2.5, 4.0 },
         { centre.x, centre.y - R, 3.0, centre.x, centre.y + R, -2.0, 4.5 },
         { centre.x - R * 0.7, centre.y - R * 0.7, 2.0, centre.x + R * 0.7, centre.y + R * 0.7, -2.0, 4.0 }
      };
      for (auto& L : cases)
      {
         auto beforeWL = [&](double x, double y) {
            return twoLobeWaterlineHeight(x, y, WZ, L.crestCx, L.crestCy, L.crestAmp * 0.4, L.troughCx, L.troughCy, L.troughAmp * 0.4, L.sigma);
         };
         auto afterWL = [&](double x, double y) {
            return twoLobeWaterlineHeight(x, y, WZ, L.crestCx, L.crestCy, L.crestAmp, L.troughCx, L.troughCy, L.troughAmp, L.sigma);
         };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideSphere(c, centre, R); },
            [&](Vec3 c) { return insideSphere(c, Vec3{centre.x + smallMove.x, centre.y + smallMove.y, centre.z + smallMove.z}, R); });
         reportAndCheck("wave: simultaneous crest+trough disconnected patches on sphere R=5", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
   }
   // Rotating cuboid case:
   {
      const int N = 44, WZ = 22;
      const Vec3 centre{ 22, 22, 22 };
      const Vec3 he{ 9, 18, 9 };
      const Vec3 smallMove{ 0.3, 0.0, 0.2 };
      struct Lobes { double crestCx, crestCy, crestAmp, troughCx, troughCy, troughAmp, sigma, rotAfter; };
      std::vector< Lobes > cases = {
         { centre.x - he.x, centre.y, 2.0, centre.x + he.x, centre.y, -2.0, 6.0, 4.0 },
         { centre.x, centre.y - he.y * 0.5, 2.5, centre.x, centre.y + he.y * 0.5, -1.5, 7.0, 7.0 }
      };
      for (auto& L : cases)
      {
         auto beforeWL = [&](double x, double y) {
            return twoLobeWaterlineHeight(x, y, WZ, L.crestCx, L.crestCy, L.crestAmp * 0.3, L.troughCx, L.troughCy, L.troughAmp * 0.3, L.sigma);
         };
         auto afterWL = [&](double x, double y) {
            return twoLobeWaterlineHeight(x, y, WZ, L.crestCx, L.crestCy, L.crestAmp, L.troughCx, L.troughCy, L.troughAmp, L.sigma);
         };
         auto r = runScenarioWave(N, N, N, beforeWL, afterWL,
            [&](Vec3 c) { return insideBoxRotY(c, centre, he, 0.0); },
            [&](Vec3 c) {
               return insideBoxRotY(c, Vec3{centre.x + smallMove.x, centre.y + smallMove.y, centre.z + smallMove.z}, he, L.rotAfter);
            });
         reportAndCheck("wave: simultaneous crest+trough disconnected patches on rotating cuboid", r, waveScenarios, pureBadWave, legacyBadWave, pass2OnlyBadWave, pureOOPWave, noLiquidAtGasOOPWave, legacyOOPWave, pass2OnlyOOPWave);
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("=== WAVE-FRONT SUMMARY: " << waveScenarios << " new wave-driven scenarios with >=1"
      << " uncovered cell | pure strict Bogner violated the invariant in " << pureBadWave << "/" << waveScenarios
      << " | legacy Pass-2+Pass-3 violated the invariant in " << legacyBadWave << "/" << waveScenarios
      << " | legacy Pass-2-only (no Pass-3) violated the invariant in " << pass2OnlyBadWave << "/" << waveScenarios
      << " | strictBognerRelaxInterface violated it in 0/" << waveScenarios << " (every case asserted above) ===");
   WALBERLA_LOG_INFO_ON_ROOT("=== WAVE-FRONT SUMMARY (out-of-place interface, summed over " << waveScenarios
      << " wave-driven scenarios): pure strict Bogner=" << pureOOPWave << " | +noLiquidAtGas=" << noLiquidAtGasOOPWave
      << " | legacy Pass-2+Pass-3=" << legacyOOPWave << " | legacy Pass-2-only (no Pass-3)=" << pass2OnlyOOPWave
      << " cells classified INTERFACE that strictBognerRelaxInterface itself resolved to LIQUID/GAS ===");

   WALBERLA_LOG_INFO_ON_ROOT("=== COMBINED SUMMARY: " << (totalScenarios + waveScenarios) << " scenarios total ("
      << totalScenarios << " rigid-body + " << waveScenarios << " wave-driven) | pure strict Bogner violated the"
      << " invariant in " << (pureBad + pureBadWave) << "/" << (totalScenarios + waveScenarios)
      << " | legacy Pass-2+Pass-3 violated the invariant in " << (legacyBad + legacyBadWave) << "/"
      << (totalScenarios + waveScenarios)
      << " | legacy Pass-2-only (no Pass-3) violated the invariant in " << (pass2OnlyBad + pass2OnlyBadWave) << "/"
      << (totalScenarios + waveScenarios) << " | strictBognerRelaxInterface violated it in 0/"
      << (totalScenarios + waveScenarios) << " (every case asserted above) ===");
   WALBERLA_LOG_INFO_ON_ROOT("=== COMBINED SUMMARY (out-of-place interface, summed over " << (totalScenarios + waveScenarios)
      << " scenarios total): pure strict Bogner=" << (pureOOP + pureOOPWave) << " | +noLiquidAtGas="
      << (noLiquidAtGasOOP + noLiquidAtGasOOPWave) << " | legacy Pass-2+Pass-3=" << (legacyOOP + legacyOOPWave)
      << " | legacy Pass-2-only (no Pass-3)=" << (pass2OnlyOOP + pass2OnlyOOPWave)
      << " cells classified INTERFACE that strictBognerRelaxInterface itself resolved to LIQUID/GAS ===");

   // Zero-evidence report over all scenarios. These cells are set to GAS, not INTERFACE, so a wrong call shows up
   // only here or as a violation, never as out-of-place interface.
   WALBERLA_LOG_INFO_ON_ROOT("=== ZERO-REAL-EVIDENCE INSTRUMENTATION (author's (4,4) counterexample re-check): "
      << g_zeroEvidenceCount << " released cells, across all " << (totalScenarios + waveScenarios)
      << " scenarios, had zero real (non-same-step) liquid/gas evidence under pure Round-0 strict Bogner | "
      << g_zeroEvidenceMismatchVsRelax << " of those disagreed with strictBognerRelaxInterface's own ground"
      << " truth for the same cell | " << g_zeroEvidenceCausesViolation << " of those forced-GAS calls"
      << " directly touch a genuinely-evidenced LIQUID cell in Round 0's own output (a real liquid/gas"
      << " contact caused by this specific branch, not merely 'forced without consequence') ===");
   for (auto& ex : g_zeroEvidenceExamples)
      WALBERLA_LOG_INFO_ON_ROOT("    mismatch example [" << ex.scenario << "] cell(" << ex.cell.x << ","
         << ex.cell.y << "," << ex.cell.z << "): pure=" << flagName(ex.pureResult) << " vs relax ground truth="
         << flagName(ex.relaxResult));

   WALBERLA_LOG_INFO_ON_ROOT("[liquid-promote-relax] full-grid 3D stress test: ok");

   return EXIT_SUCCESS;
}
