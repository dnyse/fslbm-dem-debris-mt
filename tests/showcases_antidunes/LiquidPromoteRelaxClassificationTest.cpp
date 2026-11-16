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
//! [Claude, 2026-11-12] Full-grid (not single-cell) stress test of ModelParameters.strictBognerRelaxInterface
//! (the new "Round 0/1/2" liquid-promote-then-relax classification, PdfReconstructionManager.h), built on the
//! same pure decision functions BognerStateClassificationTest.cpp already unit-tests, now driven over REAL 3D
//! sphere/cylinder/cuboid floater releases at sizes matching this thesis's own prm files (not toy numbers):
//! sphere radius 4-5 LU (prms/clog_centred_pillar_startblock_cappid_gas.prm), cylinder radius 3 / length 18 LU
//! (prms/hopper_clog_cyl_startblock_skin0p5_cappid_gas.prm), box half-extents <9,18,9> LU
//! (prms/bogner3_r025_deep_resume_gas.prm).
//!
//! This is a pure-function-level test (same convention as BognerStateClassificationTest.cpp): it drives
//! classifyFormerObstacleCellStrict/classifyFormerObstacleCell/isPass3Candidate directly over a hand-built 3D
//! flag grid with a single flat ambient interface (liquid below/interface row/gas above, matching
//! chapters/02_background/float_obj/coupling.tex's own verified convention), not a full BoundaryHandling/
//! StructuredBlockForest integration sweep -- see the final report for what that leaves untested (MPI block
//! boundaries, ghost-layer communication timing, real ParticleFunctions containment).
//!
//! For every scenario, five classification strategies are run and compared:
//!   (a) PURE strict Bogner only (Round 0, no fix at all)            -- expected to VIOLATE the invariant
//!   (b) + strictBognerNoLiquidAtGas (Round 0+1)                     -- expected safe (not asserted, logged)
//!   (c) + strictBognerRelaxInterface (Round 0+1+2, the new toggle)  -- asserted safe, fewer interface cells
//!   (d) legacy non-strict Pass-2 + Pass-3 re-resolution (bundled on) -- expected to also sometimes VIOLATE
//!   (e) legacy Pass-2 only, Pass-3 explicitly off (usePass3Reresolution unset/false) -- the actual
//!       non-strict sequential rule every one of this thesis's reported clogging and equilibrium prm
//!       runs used (confirmed: usePass3Reresolution never set true in any of them) -- measured only,
//!       not asserted, like (b) and (d)
//!
//! For every strategy OTHER than (c) itself, an additional "out-of-place interface" count is also
//! measured per scenario and summed over all scenarios: the number of released cells that strategy
//! classifies INTERFACE while (c) strictBognerRelaxInterface's own result for the SAME scenario
//! classifies that same cell LIQUID or GAS instead. This generalizes the "avoidably interface"
//! concept already used for cell (4,4) in scripts/figures/coupling_classification.py
//! (round2_relax/caution_cells there) into a reusable per-strategy test metric: it makes visible,
//! with a real number, how much of a strategy's own "0 violations" (where that holds) is actually
//! just over-classifying cells as interface rather than genuinely resolving them. Measured only,
//! like the violation counts for (b), (d), (e) -- no new hard assertion.
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

// [Claude, 2026-11-12, instrumentation for the (4,4)-counterexample re-verification requested by the
// author] Global counters/examples for cells where Round 0 pure strict Bogner's REAL (non-same-step)
// neighbourhood shows neither real liquid nor real gas evidence at all (hasGas==false && hasLiquid==
// false going into classifyStrictBognerNeighbourhood). Under the current classifyStrictBognerNeighbourhood
// ordering (`if (!hasLiquid) return Gas; if (!hasGas) return Liquid; return Interface;`), this case takes
// the FIRST branch and is assigned GAS directly -- never INTERFACE -- which is structurally different from
// scripts/figures/coupling_classification.py's bogner_together() (`if hasG and not hasL: G; elif hasL and
// not hasG: L; else: I`), whose else-branch sends the SAME zero-evidence case to INTERFACE instead (exactly
// what happens to cell (4,4) in that figure's three-cell-move worked example).
uint_t g_zeroEvidenceCount = 0;
uint_t g_zeroEvidenceMismatchVsRelax = 0; // pure's Gas call disagrees with strictBognerRelaxInterface's own final answer
struct ZeroEvidenceExample { std::string scenario; IVec3 cell; MockFlag pureResult; MockFlag relaxResult; };
std::vector< ZeroEvidenceExample > g_zeroEvidenceExamples;
uint_t g_scenarioIndex = 0;
bool g_dumpedFullExample = false;
// [Claude, 2026-11-12] Extra instrumentation: of the zero-evidence cells (forced to GAS by Round 0's
// literal first-branch-wins rule), how many directly touch a genuinely-evidenced LIQUID cell in Round
// 0's own output -- i.e. actually contribute to a real liquid/gas-contact violation, as opposed to
// merely being "forced" without consequence.
uint_t g_zeroEvidenceCausesViolation = 0;

// ---------------------------------------------------------------------------------------------
// A plain 3D grid of MockFlags, big enough to hold the whole floater plus a comfortable ambient
// margin on every side -- no out-of-range handling needed anywhere below.
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

// Single flat ambient interface at z == waterlineZ, liquid below, gas above -- exactly the
// convention verified in chapters/02_background/float_obj/coupling.tex this session.
MockFlag ambientFlag(int z, int waterlineZ)
{
   if (z < waterlineZ) return LIQUID;
   if (z == waterlineZ) return INTERFACE;
   return GAS;
}

// ---------------------------------------------------------------------------------------------
// Wave-front height fields: instead of a single flat waterlineZ, the local ambient waterline
// index is now a function of (x,y) -- still a single-interface-row column (fed straight into
// ambientFlag above), but the column height itself can tilt, curve, and carry more than one
// independent lobe, modelling a moving/deforming free-surface front rather than rigid-body motion.
// ---------------------------------------------------------------------------------------------

// A locally raised or lowered Gaussian lobe of the free-surface height field -- the building
// block for composing wave fronts with more than one independent crest/trough.
double gaussianLobe(double x, double y, double cx, double cy, double amp, double sigma)
{
   const double dx = x - cx, dy = y - cy;
   return amp * std::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
}

// Tilted-plane wave front: an angled (not axis-aligned) water-surface height field, so a sweeping
// front intersects a body on more than one side at once instead of releasing cells only along a
// single directional trailing edge (cf. the pure-translation scenarios above). Returns the REAL
// (continuous) height, not a per-column rounded index -- see ambientFlagWave below for why.
double tiltedWaterlineHeight(double x, double y, double wz0, double cx, double cy, double slopeX, double slopeY)
{
   return wz0 + slopeX * (x - cx) + slopeY * (y - cy);
}

// Curved (non-planar) wave front: a single Gaussian ridge/trough rather than a flat tilted plane.
double curvedWaterlineHeight(double x, double y, double wz0, double cx, double cy, double amp, double sigma)
{
   return wz0 + gaussianLobe(x, y, cx, cy, amp, sigma);
}

// Two-lobe wave front: an advancing crest on one side of the body superposed with a receding
// trough on the other, in the SAME height field -- so a single step can raise the ambient
// surface on the windward side while lowering it on the leeward side, i.e. two disconnected
// release "fronts" active on the same body at once.
double twoLobeWaterlineHeight(double x, double y, double wz0, double crestCx, double crestCy, double crestAmp,
                               double troughCx, double troughCy, double troughAmp, double sigma)
{
   return wz0 + gaussianLobe(x, y, crestCx, crestCy, crestAmp, sigma)
              + gaussianLobe(x, y, troughCx, troughCy, troughAmp, sigma);
}

// Interface half-width (LU) used ONLY for the continuous wave-front ambient classification below.
// The flat-waterline ambientFlag() above gets away with a single-cell-wide (0.5) interface band
// because every column shares the exact same integer threshold. That band is NOT wide enough once
// the surface tilts or curves: a cell diagonally adjacent to a neighbouring column whose threshold
// differs can then land on the wrong side of a too-thin band and touch liquid directly against
// gas -- a pure artefact of the ambient-field construction, not a finding about the classification
// engines under test. Numerically verified (this session) against every tilt/curve/two-lobe slope
// used below that a half-width of 1.5 LU leaves zero such ambient-only violations with margin.
constexpr double WAVE_INTERFACE_HALF_WIDTH = 1.5;

// Continuous wave-front ambient classification: classify directly from the cell-centre's REAL
// (continuous) height above/below the local wave surface, instead of ambientFlag's per-column
// nearest-integer snap, so the interface band stays a connected sheet wide enough to keep liquid
// and gas apart everywhere the front is tilted or curved (see WAVE_INTERFACE_HALF_WIDTH above).
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
// Floater shapes, cell-centre-inside-shape membership (matches this thesis's own MEM mapping rule,
// sec:mem-mapping / src/mesa_pd/common/Contains.h convention: centre test, not voxel-overlap).
// ---------------------------------------------------------------------------------------------
bool insideSphere(Vec3 c, Vec3 centre, double radius)
{
   const double dx = c.x - centre.x, dy = c.y - centre.y, dz = c.z - centre.z;
   return dx * dx + dy * dy + dz * dz < radius * radius;
}

// Cylinder, axis along +Y (matches FloatingCylinder's own axis<0,1,0> convention in the real prms).
bool insideCylinderY(Vec3 c, Vec3 centre, double radius, double length)
{
   const double dx = c.x - centre.x, dz = c.z - centre.z;
   if (dx * dx + dz * dz >= radius * radius) return false;
   return std::abs(c.y - centre.y) < length * 0.5;
}

// Axis-aligned-in-its-own-frame box, rotated by `rotDeg` about the Y axis (a simple, representative
// rotation -- enough to exercise "not axis-aligned" releases without a full quaternion stack).
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
// Decision engines, all built directly on the SAME pure functions PdfReconstructionManager.h
// exposes (classifyFormerObstacleCellStrict / classifyFormerObstacleCell / isPass3Candidate) --
// no re-implementation of the classification RULES themselves, only of the decide-then-apply
// orchestration around them, faithfully mirroring PdfReconstructionManager.h's own structure.
// ---------------------------------------------------------------------------------------------

// Round 0: pure strict Bogner, decide-then-apply over every newly-uncovered cell. Same-step
// neighbours (still newlyUncovered, not yet in `decided`) are reported as FORMER_OBSTACLE, which
// classifyFormerObstacleCellStrict's own B(x) excludes entirely -- exactly Bogner's literal rule.
const char* flagName(MockFlag f)
{
   if (f == LIQUID) return "LIQUID";
   if (f == GAS) return "GAS";
   if (f == INTERFACE) return "INTERFACE";
   if (f == FORMER_OBSTACLE) return "FORMER_OBSTACLE(same-step)";
   if (f == OBSTACLE) return "OBSTACLE";
   return "NONE(0)";
}

// [Claude, 2026-11-12] Instrumentation: dump the literal evidence a zero-evidence cell's 26-neighbour
// D3Q27 scan returns under the SAME exclusion rule roundZeroPureStrictBogner's own `nbr` lambda uses
// (same-step/newly-uncovered -> FORMER_OBSTACLE, out-of-range -> FORMER_OBSTACLE, else the real ambient
// flag) -- concrete evidence for task (4) rather than assertion.
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
      // [Claude, 2026-11-12] Instrumentation: hasGas==false && hasLiquid==false is the "zero real
      // evidence at all" case -- B(x) contains neither a real liquid nor a real gas neighbour (every
      // informative neighbour was same-step/FORMER_OBSTACLE-excluded, or out of range). This is the
      // case the author's counterexample (cell (4,4), scripts/figures/coupling_classification.py's
      // three-cell move) puts into bogner_together()'s else-branch -> INTERFACE. Record every
      // occurrence so the caller can report real counts, not an assertion about whether it happens.
      if (!hasGas && !hasLiquid && zeroEvidenceCellsOut != nullptr) zeroEvidenceCellsOut->push_back(c);
      decided[ambientPlusObstacle.idx(c.x, c.y, c.z)] =
         (state == BognerNeighbourhoodState::Gas) ? GAS : (state == BognerNeighbourhoodState::Liquid) ? LIQUID : INTERFACE;
   }
   return decided;
}

// Round 1: promote any Round-0 LIQUID cell touching GAS (real neighbours, OBSTACLE-covered ones
// excluded) to INTERFACE. Mirrors PdfReconstructionManager.h's strictBognerNoLiquidAtGas loop.
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

// Round 2: relax any cell currently INTERFACE back to GAS/LIQUID if its real (post-Round-1)
// neighbourhood shows only one phase. Mirrors the new strictBognerRelaxInterface pass.
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

// Legacy non-strict Pass 2 (sequential, ascending z-outer/y/x-inner -- the real non-OMP build's
// deterministic order, confirmed this session against src/field/iterators/IteratorMacros.h) then
// Pass 3 (decide-then-apply against the frozen post-Pass-2 snapshot). Mirrors the existing,
// currently-default production path.
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

   // Pass 3: frozen post-Pass-2 snapshot (every cell in `decided` is now final from Pass 2's own
   // perspective); re-run the SAME classifier, decide-then-apply.
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
         if (isUncovered[ni]) return frozen[ni]; // ALWAYS a real, final Pass-2 flag by now
         return ambientPlusObstacle.get(nx, ny, nz);
      };
      const auto state = classifyFormerObstacleCell< MockFlag >(nbr, LIQUID, GAS, FORMER_OBSTACLE, OBSTACLE, hasGas,
                                                                  hasLiquid, hasUnknown, unknownNbr);
      if (state == BognerNeighbourhoodState::Gas) afterPass3[ambientPlusObstacle.idx(cand.c.x, cand.c.y, cand.c.z)] = GAS;
      else if (state == BognerNeighbourhoodState::Liquid) afterPass3[ambientPlusObstacle.idx(cand.c.x, cand.c.y, cand.c.z)] = LIQUID;
   }
   return afterPass3;
}

// Legacy non-strict Pass 2 ONLY (same sequential, ascending z-outer/y/x-inner scan order and
// undecided-same-step-neighbour-falls-back-to-FORMER_OBSTACLE evidence rule as legacyPass2Pass3 above),
// but WITHOUT the Pass-3 local re-resolution step -- i.e. usePass3Reresolution left/forced false, the
// production default (PdfReconstructionManager.h ~L345/L578). Confirmed this session: every one of this
// thesis's actual reported clogging and equilibrium prm runs left usePass3Reresolution unset (false), so
// THIS is the strategy those runs actually executed -- legacyPass2Pass3 above (which bundles Pass-3 on)
// does not represent them.
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
      // No Pass-3 loop here at all -- this IS the difference from legacyPass2Pass3 above.
   }
   return decided;
}

// Count genuine liquid-directly-touching-gas violations (26-neighbour sense) across the WHOLE
// grid, using `decided` for cells released this step and the grid's own ambient/obstacle flags
// everywhere else.
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

// [Claude, 2026-11-12] Instrumentation: does `decided`'s own call for cell c (using the same
// flagAt fallback-to-ambient convention as countViolations) directly touch a real LIQUID neighbour?
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

// Out-of-place interface: generalizes the "avoidably interface" idea used for cell (4,4) in
// scripts/figures/coupling_classification.py (round2_relax/caution_cells there) into a reusable
// per-strategy test metric. A released cell is "out of place" under `strategyDecided` if that
// strategy classifies it INTERFACE while the known-safe ground truth `groundTruthDecided` (always
// strictBognerRelaxInterface's own final result for the SAME scenario) classifies the SAME cell
// LIQUID or GAS instead -- i.e. the strategy was unnecessarily cautious on a cell the safe strategy
// itself proves could have been resolved. This is a measurement, not a new hard assertion: it
// exists to make visible how much of a strategy's own "0 violations" (if any) is actually just
// over-classifying cells as interface rather than genuinely resolving them.
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
// Scenario harness: build before/after obstacle membership for a shape + move, derive
// newly-uncovered/newly-covered, run all four strategies, report.
// ---------------------------------------------------------------------------------------------
struct ScenarioResult { uint_t nUncovered, pureViolations, noLiquidAtGasViolations, relaxViolations, legacyViolations,
                        pass2OnlyViolations, pureInterface, relaxInterface, legacyInterface, pass2OnlyInterface,
                        // Out-of-place interface counts (this scenario): released cells a strategy calls
                        // INTERFACE that strictBognerRelaxInterface's own result for the SAME scenario calls
                        // LIQUID or GAS instead -- not computed for strictBognerRelaxInterface itself (it IS
                        // the ground truth here).
                        pureOutOfPlaceInterface, noLiquidAtGasOutOfPlaceInterface, legacyOutOfPlaceInterface,
                        pass2OnlyOutOfPlaceInterface,
                        // [Claude, 2026-11-12] Instrumentation requested by the author re: the (4,4) counterexample:
                        // how many released cells this scenario had zero real (non-same-step) liquid/gas evidence
                        // under pure Round-0 strict Bogner, and of those, how many got a DIFFERENT final phase from
                        // strictBognerRelaxInterface's own ground truth for the same cell (a genuine, not just
                        // "out of place interface", mismatch -- since classifyStrictBognerNeighbourhood's zero-
                        // evidence case resolves straight to GAS, never to INTERFACE, it can never show up in the
                        // pureOutOfPlaceInterface metric above no matter what).
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
   // ambientPlusObstacle: ambient flags everywhere not covered by the BEFORE shape; the AFTER
   // shape's own interior (whether newly covered or still covered) reported as OBSTACLE for every
   // decision engine's "currently solid" exclusion, via a second grid snapshot.
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

   // Build the single "ambientPlusObstacle" grid every engine reads non-same-step evidence from:
   // ambient flag where neither before nor after covers the cell, OBSTACLE where `after` still/newly
   // covers it (so a true remaining-obstacle or newly-covered neighbour is correctly excluded), and
   // left at its ambient value (never queried directly) for the newly-uncovered cells themselves.
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

   // Out-of-place interface: measured against strictBognerRelaxInterface's own result (withRelax)
   // for this SAME scenario, for every strategy other than itself.
   r.pureOutOfPlaceInterface         = countOutOfPlaceInterface(uncovered, evidence, pure, withRelax);
   r.noLiquidAtGasOutOfPlaceInterface = countOutOfPlaceInterface(uncovered, evidence, withPromote, withRelax);
   r.legacyOutOfPlaceInterface       = countOutOfPlaceInterface(uncovered, evidence, legacy, withRelax);
   r.pass2OnlyOutOfPlaceInterface    = countOutOfPlaceInterface(uncovered, evidence, pass2Only, withRelax);

   // [Claude, 2026-11-12] Instrumentation (author request): real count of zero-real-evidence released
   // cells under pure Round-0 strict Bogner, and how many of those disagree with strictBognerRelaxInterface's
   // own ground truth for the SAME cell once the full picture is known.
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
// Wave-driven scenario harness: the SAME four-strategy comparison as runScenario above, but the
// ambient liquid/gas evidence a newly-uncovered cell sees comes from a MOVING/DEFORMING free-
// surface height field (tilted, curved, crest-then-trough, two-lobe) instead of a single flat
// waterline -- while the body itself can still translate/rotate between before and after, exactly
// as in runScenario (a wave-struck floater is not held perfectly fixed). beforeWL/afterWL give
// the local waterline Z index (an int, same single-interface-row convention as ambientFlag) as a
// function of a cell's (x,y) centre coordinates, evaluated BEFORE/AFTER the step respectively.
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

   // Evidence grid: the AFTER body shape as OBSTACLE, the AFTER wave field everywhere else --
   // the post-step ambient reality every decision engine reads non-same-step evidence from.
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

   // Out-of-place interface: measured against strictBognerRelaxInterface's own result (withRelax)
   // for this SAME scenario, for every strategy other than itself.
   r.pureOutOfPlaceInterface         = countOutOfPlaceInterface(uncovered, evidence, pure, withRelax);
   r.noLiquidAtGasOutOfPlaceInterface = countOutOfPlaceInterface(uncovered, evidence, withPromote, withRelax);
   r.legacyOutOfPlaceInterface       = countOutOfPlaceInterface(uncovered, evidence, legacy, withRelax);
   r.pass2OnlyOutOfPlaceInterface    = countOutOfPlaceInterface(uncovered, evidence, pass2Only, withRelax);

   // [Claude, 2026-11-12] Instrumentation (author request): same zero-real-evidence count/mismatch
   // tracking as runScenario above, for the wave-driven scenarios.
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
   // The one hard safety assertion this whole file exists to make:
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
   // translating through the waterline in every axis-aligned and diagonal direction, 1-4 LU/step.
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

   // ---- Two spheres approaching/overlapping the same patch of waterline simultaneously ----
   // (as in this thesis's own clogging/jam scenarios, e.g. FloatingSphere entries spaced 12 LU
   // apart in clog_centred_pillar_startblock_cappid_gas.prm) -- released cells from BOTH spheres
   // can be same-step neighbours of each other near the gap between them.
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

   // ---- Randomised multi-sphere jam patch: several floaters of varying size moving in random ----
   // directions near the same stretch of waterline, closer to an actual clogging snapshot than
   // any single-floater case above.
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
   // WAVE-FRONT scenarios: the before/after obstacle-liquid-interface-gas state is now driven by a
   // MOVING/DEFORMING free-surface front (tilted, curved, crest-then-trough, two-lobe) rather than
   // purely by the floater's own rigid translation/rotation -- while the body itself still nudges
   // slightly between steps, as a real wave-struck floater would. Separate counters so the new
   // category's own pass/fail is visible alongside the combined total with the 43 scenarios above.
   // =============================================================================================
   walberla::uint_t waveScenarios = 0, pureBadWave = 0, legacyBadWave = 0, pass2OnlyBadWave = 0;
   // Out-of-place interface running totals (summed across every wave-driven scenario below).
   walberla::uint_t pureOOPWave = 0, noLiquidAtGasOOPWave = 0, legacyOOPWave = 0, pass2OnlyOOPWave = 0;

   // ---- Tilted wave front (not axis-aligned) sweeping across sphere R=5 at an angle, so it ----
   // intersects the body on multiple sides at once rather than releasing cells only along one
   // directional trailing edge. Body also nudges slightly, as a struck floater would.
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

   // ---- Curved (Gaussian-ridge) wave front sweeping across cylinder R=3 L=18, a non-planar ----
   // front rather than a flat tilted plane, amplitude growing between before/after as it sweeps.
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

   // ---- Crest advances then, in the NEXT step, a trough recedes over the same side of the ----
   // body: covered-then-uncovered across two consecutive steps rather than monotonic release.
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

   // ---- Simultaneous, disconnected release patches on the SAME body: a windward side still ----
   // being covered by an advancing crest while the leeward side is uncovered by a receding
   // trough, both lobes active in the same single step -- more than one release "front" at once.
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

   // [Claude, 2026-11-12] Instrumentation (author request, re-verifying the (4,4) counterexample):
   // real count, across all 61 scenarios (43 rigid-body + 18 wave-driven), of released cells whose
   // REAL (non-same-step) D3Q27 neighbourhood shows neither a real liquid nor a real gas cell at all
   // under pure Round-0 strict Bogner -- the exact situation cell (4,4) is in, in the figure's three-
   // cell-move worked example. classifyStrictBognerNeighbourhood's `if (!hasLiquid) return Gas; ...`
   // ordering sends this case straight to GAS, never INTERFACE, so it structurally can never appear in
   // the "out-of-place interface" metric above (that metric only counts cells a strategy calls
   // INTERFACE) -- a genuine mismatch here would instead show up as a silent wrong-phase call, visible
   // only via this separate counter or as a knock-on violation count.
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
