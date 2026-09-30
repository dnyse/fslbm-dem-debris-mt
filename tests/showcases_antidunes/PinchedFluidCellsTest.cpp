//======================================================================================================================
//! \file PinchedFluidCellsTest.cpp
//! Geometric precondition of the fluid-side instability seen in the cylinder hopper (THESIS_PLAN.md, class (B)): fluid cells that are SQUEEZED BETWEEN TWO OBSTACLES.
//! Cells are mapped exactly like lbm_mesapd_coupling::ParticleMapping does it: a lattice cell (dx = 1, centre at i+0.5) is obstacle iff its centre lies inside a body. A fluid cell is "pinched" when for
//! at least one of the 13 opposite direction pairs (d, -d) of the D3Q27 stencil BOTH neighbours are obstacle, i.e. it has a wall on two opposite sides at distance one cell (the fluid cell that exploded in job
//! 4251066 had such pairs, e.g. (1,-1,-1)/(-1,1,1) and (1,1,-1)/(-1,-1,1), and 11 of 26 neighbours were obstacle cells of two different rods).
//! A single convex body never creates a pinched cell; two bodies at a small gap do. We count them for two rods (r=3, L=18, side by side, parallel or slightly tilted) and for two spheres (same radius r=3, and
//! the same VOLUME r=4.95), over gaps 0..2 cells with random sub-cell registration and orientation. Lubrication (LubricationCorrectionKernel, cut-off 2/3 cell) exists for sphere-sphere pairs only.
//======================================================================================================================
#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace {
using V = std::array< double, 3 >;
V operator+(const V& a, const V& b) { return { a[0] + b[0], a[1] + b[1], a[2] + b[2] }; }
V operator-(const V& a, const V& b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
V operator*(const V& a, double s) { return { a[0] * s, a[1] * s, a[2] * s }; }
double dot(const V& a, const V& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V cross(const V& a, const V& b) { return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] }; }
V unit(const V& a) { const double l = std::sqrt(dot(a, a)); return a * (1.0 / l); }

struct Body { bool sphere; V c; V axis; double R, halfL; };
bool inside(const Body& b, const V& p)
{
   const V d = p - b.c;
   if (b.sphere) return dot(d, d) <= b.R * b.R;
   const double z = dot(d, b.axis); const V r = d - b.axis * z;
   return std::abs(z) <= b.halfL && dot(r, r) <= b.R * b.R;
}

constexpr int N = 40;                       // 40^3 voxel window around the pair
struct Grid { std::vector< char > obst; Grid() : obst(N * N * N, 0) {} char& at(int x, int y, int z) { return obst[(size_t(z) * N + y) * N + x]; } };

Grid voxelize(const std::vector< Body >& bodies)
{
   Grid g;
   for (int z = 0; z < N; ++z) for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x)
   {
      const V p{ x + 0.5, y + 0.5, z + 0.5 };
      for (const auto& b : bodies) if (inside(b, p)) { g.at(x, y, z) = 1; break; }
   }
   return g;
}

//! number of fluid cells with an obstacle on two opposite sides (over the 13 direction pairs of D3Q27); optionally histogram by the number of such pairs
int pinchedCells(Grid& g, int* atLeast3 = nullptr)
{
   static std::vector< std::array< int, 3 > > half;
   if (half.empty())
      for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
         if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0)) half.push_back({ dx, dy, dz });
   int n = 0, n3 = 0;
   for (int z = 1; z < N - 1; ++z) for (int y = 1; y < N - 1; ++y) for (int x = 1; x < N - 1; ++x)
   {
      if (g.at(x, y, z)) continue;
      int pairs = 0;
      for (const auto& d : half) pairs += (g.at(x + d[0], y + d[1], z + d[2]) && g.at(x - d[0], y - d[1], z - d[2]));
      n += pairs > 0; n3 += pairs >= 3;
   }
   if (atLeast3) *atLeast3 = n3;
   return n;
}
}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };
   std::mt19937 gen(99); std::uniform_real_distribution< double > ud(0.0, 1.0); std::normal_distribution< double > nd(0.0, 1.0);
   const double R = 3.0, L = 18.0, Rvol = std::cbrt(3.0 * (3.14159265358979323846 * R * R * L) / (4.0 * 3.14159265358979323846));   // sphere with the rod's volume

   // sanity: a single body never creates a pinched cell
   {
      int worst = 0;
      for (int t = 0; t < 30; ++t)
      {
         V ax = unit(V{ nd(gen), nd(gen), nd(gen) }); const V off{ ud(gen), ud(gen), ud(gen) }; const V c = V{ 20, 20, 20 } + off;
         for (bool sph : { false, true }) { auto g = voxelize({ Body{ sph, c, ax, sph ? Rvol : R, 0.5 * L } }); worst = std::max(worst, pinchedCells(g)); }
      }
      std::cout << "single body (rod or sphere), 60 random placements: max pinched cells = " << worst << "\n";
      check(worst == 0, "a single convex body must never create a pinched fluid cell");
   }

   std::cout << "\nmean number of pinched fluid cells per PAIR of bodies (300 random placements/orientations per row; in brackets: cells with >= 3 opposite obstacle pairs)\n"
             << "   gap [cells] | rod-rod (r=3, L=18) | sphere-sphere r=3 | sphere-sphere same volume (r=" << Rvol << ") | rod / sphere(r=3)\n";
   double rodSum = 0, sphSum = 0, sphVolSum = 0; int rows = 0;
   for (double gap : { 0.0, 0.25, 0.5, 1.0, 1.5, 2.0 })
   {
      double rod = 0, sph = 0, sphv = 0, rod3 = 0, sph3 = 0, sphv3 = 0; const int T = 300;
      for (int t = 0; t < T; ++t)
      {
         const V off{ ud(gen), ud(gen), ud(gen) }; const V mid = V{ 20, 20, 20 } + off;
         V a = unit(V{ nd(gen), nd(gen), nd(gen) });                       // rod axis
         V s = unit(cross(a, unit(V{ nd(gen), nd(gen), nd(gen) })));       // side-by-side direction, perpendicular to the axis
         const double shift = 8.0 * (ud(gen) - 0.5);                         // axial offset of the second rod
         const double tilt = 0.17 * ud(gen);                                // up to ~10 degrees
         V a2 = unit(a + s * 0.0 + unit(cross(a, s)) * std::tan(tilt));
         { int n3; auto g = voxelize({ Body{ false, mid - s * (R + 0.5 * gap), a, R, 0.5 * L }, Body{ false, mid + s * (R + 0.5 * gap) + a * shift, a2, R, 0.5 * L } }); rod += pinchedCells(g, &n3); rod3 += n3; }
         V d = unit(V{ nd(gen), nd(gen), nd(gen) });
         { int n3; auto g = voxelize({ Body{ true, mid - d * (R + 0.5 * gap), a, R, 0 }, Body{ true, mid + d * (R + 0.5 * gap), a, R, 0 } }); sph += pinchedCells(g, &n3); sph3 += n3; }
         { int n3; auto g = voxelize({ Body{ true, mid - d * (Rvol + 0.5 * gap), a, Rvol, 0 }, Body{ true, mid + d * (Rvol + 0.5 * gap), a, Rvol, 0 } }); sphv += pinchedCells(g, &n3); sphv3 += n3; }
      }
      rod /= T; sph /= T; sphv /= T; rod3 /= T; sph3 /= T; sphv3 /= T;
      std::printf("   %5.2f       | %6.1f  [%5.1f]      | %6.1f  [%5.1f]     | %6.1f  [%5.1f]                        | %5.1fx\n", gap, rod, rod3, sph, sph3, sphv, sphv3, sph > 0 ? rod / sph : 0.0);
      if (gap <= 1.0) { rodSum += rod; sphSum += sph; sphVolSum += sphv; ++rows; }
   }
   std::printf("average over gaps <= 1 cell: rods %.1f, spheres (r=3) %.1f, spheres (same volume) %.1f pinched cells per pair\n", rodSum / rows, sphSum / rows, sphVolSum / rows);
   check(rodSum > 2.0 * sphSum, "two rods at small gaps must produce at least 2x as many pinched fluid cells as two spheres of the same radius");
   check(rodSum > 1.5 * sphVolSum, "two rods at small gaps must produce more pinched fluid cells than two spheres of the same volume");
   std::cout << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failures)" << std::endl;
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
