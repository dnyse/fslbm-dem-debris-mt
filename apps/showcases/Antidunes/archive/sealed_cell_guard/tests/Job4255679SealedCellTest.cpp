//======================================================================================================================
//! \file Job4255679SealedCellTest.cpp
//! Proof, using the REAL geometry of job 4255679's own crash, that the sealed/pinched-cell fluid instability
//! (THESIS_PLAN.md class (B), root-caused in job 4252462) applies here too -- extends PinchedFluidCellsTest.cpp's
//! own voxelization method (same cell-centre-inside-body mapping as lbm_mesapd_coupling::ParticleMapping) from random
//! synthetic placements to the EXACT logged positions of the two culprit rods.
//!
//! Job 4255679 (antidunes-hopper-clog-cyl-halfrows-gjk-noguards-60k) aborted at step 3133 via STALE-PARTICLE-IDX.
//! Its own [CylForces] log gives the real world positions of rods uid 3530822107858468885/888 at step 3100 (radius
//! 3, length 18, axis close to world Y -- see this session's own TiltYOwnAxisSpinArtifactTest.cpp for why the
//! logged tiltY is NOT reliable evidence of a large axis tilt, so axis = world Y is used here, the simplest
//! defensible reading of the data): centre separation 5.931 LU, i.e. essentially touching (2R = 6). The run's own
//! WARNING-level log (not the fatal ERROR line) directly confirms what this predicts: at step 3132,
//! [MEM-GARBAGE-LINK] fires for a liquid cell bordering BOTH uids' obstacle surfaces, with populations of order
//! 1e9-1e12 and density collapsed to -0.066; by step 3133 both rods are launched to ~138,000 LU/step (vs ~0.02
//! LU/step normal cruising speed).
//!
//! This test shows the GEOMETRIC precondition is present at the real positions (not just in a random gap sweep):
//! pinched cells exist, strongly-pinched (>=3 opposite obstacle pairs, matching job 4251066's own exploded-cell
//! signature) cells exist, and -- the key point -- cells meeting closeSealedCellsMinNeighbours=24 (the threshold
//! actually used by SealedCellGuard.h in the guarded rerun, job 4259313) exist too, i.e. the guard would genuinely
//! have found and closed a cell here, not merely a threshold that happens not to trigger for this exact geometry.
//======================================================================================================================
#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
using V = std::array< double, 3 >;
V operator-(const V& a, const V& b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
double dot(const V& a, const V& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

struct Body { V c; V axis; double R, halfL; };
bool inside(const Body& b, const V& p)
{
   const V d = p - b.c;
   const double z = dot(d, b.axis);
   const V r = d - V{ b.axis[0] * z, b.axis[1] * z, b.axis[2] * z };
   return std::abs(z) <= b.halfL && dot(r, r) <= b.R * b.R;
}

constexpr int N = 44; // voxel window, world-cell-sized, centred on the pair's midpoint
struct Grid
{
   std::vector< char > obst;
   V origin; // world position of voxel (0,0,0)'s lower corner
   Grid() : obst(N * N * N, 0) {}
   char& at(int x, int y, int z) { return obst[(size_t(z) * N + y) * N + x]; }
};

Grid voxelize(const std::vector< Body >& bodies, const V& centre)
{
   Grid g;
   g.origin = { centre[0] - N / 2.0, centre[1] - N / 2.0, centre[2] - N / 2.0 };
   for (int z = 0; z < N; ++z)
      for (int y = 0; y < N; ++y)
         for (int x = 0; x < N; ++x)
         {
            const V p{ g.origin[0] + x + 0.5, g.origin[1] + y + 0.5, g.origin[2] + z + 0.5 };
            for (const auto& b : bodies)
               if (inside(b, p)) { g.at(x, y, z) = 1; break; }
         }
   return g;
}

struct Stats { int pinched = 0, strong3 = 0, sealed24 = 0; int worstObstacleNeighbours = 0; };

Stats scan(Grid& g)
{
   static std::vector< std::array< int, 3 > > half;
   if (half.empty())
      for (int dz = -1; dz <= 1; ++dz)
         for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
               if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0)) half.push_back({ dx, dy, dz });
   static std::vector< std::array< int, 3 > > full;
   if (full.empty())
      for (int dz = -1; dz <= 1; ++dz)
         for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
               if (!(dx == 0 && dy == 0 && dz == 0)) full.push_back({ dx, dy, dz });

   Stats s;
   for (int z = 1; z < N - 1; ++z)
      for (int y = 1; y < N - 1; ++y)
         for (int x = 1; x < N - 1; ++x)
         {
            if (g.at(x, y, z)) continue;
            int pairs = 0;
            for (const auto& d : half) pairs += (g.at(x + d[0], y + d[1], z + d[2]) && g.at(x - d[0], y - d[1], z - d[2]));
            int nObst = 0;
            for (const auto& d : full) nObst += g.at(x + d[0], y + d[1], z + d[2]);
            s.pinched += pairs > 0;
            s.strong3 += pairs >= 3;
            s.sealed24 += nObst >= 24;
            s.worstObstacleNeighbours = std::max(s.worstObstacleNeighbours, nObst);
         }
   return s;
}
}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };

   // Real logged step-3100 positions from job 4255679's own stdout ([CylForces], rank 49/50), radius/length from
   // its own prm (hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm: FloatingCylinder radius 3 length 18).
   const V pos885{ 287.687, 9.46882, 57.8675 };
   const V pos888{ 283.45, 10.1355, 61.9644 };
   const V axisY{ 0.0, 1.0, 0.0 };
   const double R = 3.0, halfL = 9.0;

   const double sep = std::sqrt(dot(pos888 - pos885, pos888 - pos885));
   std::cout << "real step-3100 centre separation: " << sep << " LU (2R = " << 2 * R << ")\n";
   check(sep < 2 * R + 1.0, "the real logged separation must be near-touching (< 2R+1) for this test to be meaningful");

   const V mid{ (pos885[0] + pos888[0]) / 2.0, (pos885[1] + pos888[1]) / 2.0, (pos885[2] + pos888[2]) / 2.0 };
   auto g = voxelize({ Body{ pos885, axisY, R, halfL }, Body{ pos888, axisY, R, halfL } }, mid);
   const Stats s = scan(g);

   std::cout << "at the REAL job-4255679 geometry (axis = world Y, the defensible reading -- see file header):\n"
             << "  pinched cells (>=1 opposite obstacle pair):        " << s.pinched << "\n"
             << "  strongly pinched cells (>=3 pairs, job 4251066's own signature): " << s.strong3 << "\n"
             << "  cells meeting closeSealedCellsMinNeighbours=24:    " << s.sealed24 << "\n"
             << "  worst single cell's obstacle-neighbour count:      " << s.worstObstacleNeighbours << " of 26\n";

   check(s.pinched > 0, "the real job-4255679 rod pair must produce at least one pinched fluid cell");
   check(s.strong3 > 0, "the real job-4255679 rod pair must produce at least one strongly-pinched (>=3 pair) cell, "
                          "matching job 4251066's own exploded-cell signature");
   check(s.sealed24 > 0, "the real job-4255679 rod pair must produce at least one cell meeting the guard's own "
                          "closeSealedCellsMinNeighbours=24 threshold -- otherwise the guard as configured would "
                          "not actually have found anything to close here");

   std::cout << (failures == 0 ?
      "\nALL CHECKS PASSED -- the sealed/pinched-cell precondition is present at the EXACT real geometry of job\n"
      "4255679's crash, not just in a random gap sweep, and the guard24 threshold used in the rerun (job 4259313)\n"
      "genuinely applies to it.\n"
      : "\nSOME CHECKS FAILED -- see FAIL lines above.\n");
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
