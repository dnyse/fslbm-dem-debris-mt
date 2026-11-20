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
//! \file QuarticSolverComparisonTest.cpp
//! Time/accuracy comparison of the three quartic root finders in the `detail` namespace of
//! ClosedFormContactsCylinderPlaneGeometry.h: solveQuarticDurandKerner (default; at its 30-iteration default and at
//! 60), solveQuarticFerrari (closed form, no iteration) and solveQuarticAberthEhrlich (same simultaneous-refinement
//! family with a full Newton correction: cubic instead of quadratic convergence order; tested at 30 and at 10
//! iterations to see whether it gets away with fewer sweeps).
//!
//! Quartics are built from known roots, so each solver is compared with ground truth rather than with another
//! solver. Five root configurations, chosen because they stress the methods differently:
//!   0) four real roots, widely separated                     -- easy for every method.
//!   1) four real roots with one close pair                   -- stresses the iterative convergence rate.
//!   2) two real roots + one complex-conjugate pair           -- the common case in the geometric application
//!      (a quartic in tan(theta/2) generally has 0 or 2 real roots, not 4).
//!   3) two complex-conjugate pairs, no real roots            -- the "no valid theta" case.
//!   4) a near-repeated real root pair                        -- the hard case: Ferrari's resolvent-cubic
//!      discriminant is near its branch boundary, and the Durand-Kerner/Aberth-Ehrlich denominators (products/sums
//!      of root differences) get small and convergence slows.
//!
//! Timing uses one fixed category-2 quartic in a tight call loop, as in CylinderCylinderContactBenchmark.cpp
//! (volatile sink against dead-code elimination; fixed input, so only per-call cost is measured).
//! Informs whether solveQuarticDurandKerner should be replaced.
//======================================================================================================================

#include "ClosedFormContacts.h"

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{
namespace quartic_compare
{
using Complex = std::complex< real_t >;

//======================================================================================================================
// Ground truth: expand c4*(x-r0)(x-r1)(x-r2)(x-r3) from four known roots (real or conjugate pairs, so the
// coefficients are real by construction; .real() only strips floating-point noise).
//======================================================================================================================
struct QuarticProblem
{
   real_t  c4, c3, c2, c1, c0;
   Complex trueRoots[4];
};

QuarticProblem makeFromRoots(real_t c4, Complex r0, Complex r1, Complex r2, Complex r3)
{
   const Complex e1 = r0 + r1 + r2 + r3;
   const Complex e2 = r0 * r1 + r0 * r2 + r0 * r3 + r1 * r2 + r1 * r3 + r2 * r3;
   const Complex e3 = r0 * r1 * r2 + r0 * r1 * r3 + r0 * r2 * r3 + r1 * r2 * r3;
   const Complex e4 = r0 * r1 * r2 * r3;

   QuarticProblem p;
   p.c4        = c4;
   p.c3        = c4 * (-e1.real());
   p.c2        = c4 * e2.real();
   p.c1        = c4 * (-e3.real());
   p.c0        = c4 * e4.real();
   p.trueRoots[0] = r0;
   p.trueRoots[1] = r1;
   p.trueRoots[2] = r2;
   p.trueRoots[3] = r3;
   return p;
}

struct FuzzRng
{
   std::mt19937 rng;
   explicit FuzzRng(uint32_t seed) : rng(seed) {}
};

QuarticProblem randomProblem(FuzzRng& fuzz, int category)
{
   std::uniform_real_distribution< real_t > wide(real_t(-8), real_t(8));
   std::uniform_real_distribution< real_t > c4Dist(real_t(0.2), real_t(5.0));
   std::uniform_real_distribution< real_t > tiny(real_t(1e-4), real_t(1e-2));
   std::uniform_real_distribution< real_t > imagDist(real_t(0.3), real_t(6.0));
   const real_t                             c4 = c4Dist(fuzz.rng);

   switch (category)
   {
   case 0: {
      real_t r[4];
      for (auto& x : r)
         x = wide(fuzz.rng);
      return makeFromRoots(c4, r[0], r[1], r[2], r[3]);
   }
   case 1: {
      const real_t a = wide(fuzz.rng), b = wide(fuzz.rng), eps = tiny(fuzz.rng);
      return makeFromRoots(c4, a, a + eps, b, b - eps * real_t(0.5));
   }
   case 2: {
      const real_t a = wide(fuzz.rng), b = wide(fuzz.rng);
      const real_t re = wide(fuzz.rng), im = imagDist(fuzz.rng);
      return makeFromRoots(c4, Complex(a, 0), Complex(b, 0), Complex(re, im), Complex(re, -im));
   }
   case 3: {
      const real_t re1 = wide(fuzz.rng), im1 = imagDist(fuzz.rng);
      const real_t re2 = wide(fuzz.rng), im2 = imagDist(fuzz.rng);
      return makeFromRoots(c4, Complex(re1, im1), Complex(re1, -im1), Complex(re2, im2), Complex(re2, -im2));
   }
   default: {
      const real_t a = wide(fuzz.rng), b = wide(fuzz.rng);
      return makeFromRoots(c4, a, a, b, b + tiny(fuzz.rng));
   }
   }
}

//! Greedy nearest-neighbour matching of the 4 true roots to a solver's 4 roots; returns the worst (max) matched
//! distance, because a mean would dilute a single badly wrong root.
real_t maxMatchedRootError(const Complex trueRoots[4], const Complex computed[4])
{
   bool   used[4] = { false, false, false, false };
   real_t maxErr  = 0;
   for (int i = 0; i < 4; ++i)
   {
      int    best     = -1;
      real_t bestDist = std::numeric_limits< real_t >::max();
      for (int j = 0; j < 4; ++j)
      {
         if (used[j]) continue;
         const real_t d = std::abs(trueRoots[i] - computed[j]);
         if (d < bestDist)
         {
            bestDist = d;
            best     = j;
         }
      }
      used[best] = true;
      maxErr     = std::max(maxErr, bestDist);
   }
   return maxErr;
}

struct AccuracyStats
{
   int    trials = 0;
   real_t sumMaxErr = 0, worstMaxErr = 0;
};

void accumulate(AccuracyStats& s, real_t maxErr)
{
   ++s.trials;
   s.sumMaxErr += maxErr;
   s.worstMaxErr = std::max(s.worstMaxErr, maxErr);
}

void report(const char* categoryLabel, const char* solverLabel, const AccuracyStats& s)
{
   WALBERLA_LOG_INFO_ON_ROOT("[QuarticAccuracy] " << categoryLabel << " | " << solverLabel
                                                   << ": mean max-root-error=" << (s.sumMaxErr / real_t(s.trials))
                                                   << " | worst max-root-error=" << s.worstMaxErr);
}

//======================================================================================================================
// Timing harness as in CylinderCylinderContactBenchmark.cpp's timeIt: a volatile sink keeps the optimizer from
// deleting the loop.
//======================================================================================================================
volatile const void* g_optimizationSink = nullptr;

template< typename Fn >
double timeItNs(uint64_t iterations, Fn&& fn)
{
   for (uint64_t i = 0; i < std::min< uint64_t >(iterations / 10, uint64_t(1000)); ++i)
      fn();

   const auto start = std::chrono::steady_clock::now();
   for (uint64_t i = 0; i < iterations; ++i)
      fn();
   const auto end = std::chrono::steady_clock::now();

   return double(std::chrono::duration_cast< std::chrono::nanoseconds >(end - start).count()) / double(iterations);
}

} // namespace quartic_compare
} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla

//======================================================================================================================
// Main
//======================================================================================================================
int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla;
   using namespace walberla::antidunes::closed_form_contacts;
   using namespace walberla::antidunes::closed_form_contacts::quartic_compare;

   //===================================================================================================================
   // Phase 1: accuracy against known roots, 5 root-configuration categories.
   //===================================================================================================================
   constexpr int numTrialsPerCategory = 4000;
   const char*   categoryLabels[5]    = { "0) 4 real, well-separated", "1) 4 real, one CLOSE pair",
                                        "2) 2 real + 1 conjugate pair", "3) 2 conjugate pairs (no real roots)",
                                        "4) near-REPEATED real root" };

   FuzzRng fuzz(90210);
   for (int category = 0; category < 5; ++category)
   {
      AccuracyStats dk30, dk60, ferrari, ae30, ae10;
      for (int trial = 0; trial < numTrialsPerCategory; ++trial)
      {
         const QuarticProblem p = randomProblem(fuzz, category);

         Complex rootsDk30[4], rootsDk60[4], rootsFerrari[4], rootsAe30[4], rootsAe10[4];
         detail::solveQuarticDurandKerner(p.c4, p.c3, p.c2, p.c1, p.c0, rootsDk30, 30);
         detail::solveQuarticDurandKerner(p.c4, p.c3, p.c2, p.c1, p.c0, rootsDk60, 60);
         detail::solveQuarticFerrari(p.c4, p.c3, p.c2, p.c1, p.c0, rootsFerrari);
         detail::solveQuarticAberthEhrlich(p.c4, p.c3, p.c2, p.c1, p.c0, rootsAe30, 30);
         detail::solveQuarticAberthEhrlich(p.c4, p.c3, p.c2, p.c1, p.c0, rootsAe10, 10);

         accumulate(dk30, maxMatchedRootError(p.trueRoots, rootsDk30));
         accumulate(dk60, maxMatchedRootError(p.trueRoots, rootsDk60));
         accumulate(ferrari, maxMatchedRootError(p.trueRoots, rootsFerrari));
         accumulate(ae30, maxMatchedRootError(p.trueRoots, rootsAe30));
         accumulate(ae10, maxMatchedRootError(p.trueRoots, rootsAe10));
      }
      report(categoryLabels[category], "Durand-Kerner (30, current default)", dk30);
      report(categoryLabels[category], "Durand-Kerner (60, original)       ", dk60);
      report(categoryLabels[category], "Ferrari (closed-form)              ", ferrari);
      report(categoryLabels[category], "Aberth-Ehrlich (30)                ", ae30);
      report(categoryLabels[category], "Aberth-Ehrlich (10)                ", ae10);
      WALBERLA_LOG_INFO_ON_ROOT("");
   }

   //===================================================================================================================
   // Phase 2: timing, one fixed category-2 quartic (2 real roots + 1 conjugate pair), tight call loop.
   //===================================================================================================================
   constexpr uint64_t N = 2000000;
   const QuarticProblem timingProblem =
      makeFromRoots(real_t(1.7), Complex(real_t(-2.3), 0), Complex(real_t(1.1), 0), Complex(real_t(0.4), real_t(2.1)),
                    Complex(real_t(0.4), real_t(-2.1)));

   Complex sinkRoots[4];
   const double nsDk30 = timeItNs(N, [&]() {
      detail::solveQuarticDurandKerner(timingProblem.c4, timingProblem.c3, timingProblem.c2, timingProblem.c1,
                                        timingProblem.c0, sinkRoots, 30);
      g_optimizationSink = &sinkRoots[0];
   });
   const double nsDk60 = timeItNs(N, [&]() {
      detail::solveQuarticDurandKerner(timingProblem.c4, timingProblem.c3, timingProblem.c2, timingProblem.c1,
                                        timingProblem.c0, sinkRoots, 60);
      g_optimizationSink = &sinkRoots[0];
   });
   const double nsFerrari = timeItNs(N, [&]() {
      detail::solveQuarticFerrari(timingProblem.c4, timingProblem.c3, timingProblem.c2, timingProblem.c1,
                                   timingProblem.c0, sinkRoots);
      g_optimizationSink = &sinkRoots[0];
   });
   const double nsAe30 = timeItNs(N, [&]() {
      detail::solveQuarticAberthEhrlich(timingProblem.c4, timingProblem.c3, timingProblem.c2, timingProblem.c1,
                                         timingProblem.c0, sinkRoots, 30);
      g_optimizationSink = &sinkRoots[0];
   });
   const double nsAe10 = timeItNs(N, [&]() {
      detail::solveQuarticAberthEhrlich(timingProblem.c4, timingProblem.c3, timingProblem.c2, timingProblem.c1,
                                         timingProblem.c0, sinkRoots, 10);
      g_optimizationSink = &sinkRoots[0];
   });

   struct TimingRow
   {
      std::string label;
      double      ns;
   };
   std::vector< TimingRow > rows = { { "Durand-Kerner (30, current default)", nsDk30 },
                                      { "Durand-Kerner (60, original)", nsDk60 },
                                      { "Ferrari (closed-form)", nsFerrari },
                                      { "Aberth-Ehrlich (30)", nsAe30 },
                                      { "Aberth-Ehrlich (10)", nsAe10 } };
   std::sort(rows.begin(), rows.end(), [](const TimingRow& a, const TimingRow& b) { return a.ns > b.ns; });

   WALBERLA_LOG_INFO_ON_ROOT("[QuarticTiming] === Results, slowest to fastest (N=" << N << ") ===");
   const double fastest = rows.back().ns;
   for (const auto& row : rows)
      WALBERLA_LOG_INFO_ON_ROOT("[QuarticTiming]   " << row.label << ": " << row.ns << " ns/call  ("
                                                      << (row.ns / fastest) << "x the fastest entry)");

   return EXIT_SUCCESS;
}
