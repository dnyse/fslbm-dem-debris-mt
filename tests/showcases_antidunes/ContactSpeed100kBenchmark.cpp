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
//! \file ContactSpeed100kBenchmark.cpp
//! Microbenchmark, NOT a correctness test. Unlike CylinderCylinderContactBenchmark.cpp (one fixed,
//! hand-picked geometry per scenario, repeated identically 500000 times -- good for isolating a single
//! scenario's own cost, but not representative of a realistic mix of easy/hard cases, and vulnerable to
//! the branch predictor/cache learning that one fixed input), this one replays the SAME 100000-trial
//! random box_pairs.csv/cyl_pairs.csv dump already used for this thesis's own accuracy validation
//! (Section 5.2.1; ANTIDUNES_REFERENCE_DUMP_ONLY=1 ANTIDUNES_REFERENCE_TRIALS=100000, seed 13579) and
//! times the closed form against GJK/EPA over that whole realistic distribution (contacts and
//! near-misses alike, exactly the sampling band testIndependentReferenceDump already uses), reporting
//! one average ns/call per method.
//!
//! WHAT IS AND ISN'T TIMED: every row is parsed from the CSV and converted into plain Vec3/Mat3/Rot3
//! values BEFORE any timer starts -- file I/O and parsing are never in a timed region. The box-box SAT
//! and cylinder-cylinder closed form calls are pure, stateless free functions (boxBoxSAT,
//! cylinderCylinderContact with PRODUCTION DEFAULT arguments, i.e. the Ferrari quartic solver and
//! lazyThetaRef=true) taking the precomputed Vec3/Mat3 values directly -- nothing but the call itself is
//! timed. GJK/EPA needs mesa_pd particles, so (matching CylinderCylinderContactBenchmark.cpp's own
//! approach) the ParticleStorage and the two particles are created ONCE outside the timed region; per
//! trial, only updating their position/rotation/shape to the next row's values and calling
//! SolidCylinderContactDetection's operator() is timed -- this per-trial update is a handful of struct
//! field writes (nanoseconds), not a new allocation, and is the same marshaling any production call
//! into GJK/EPA with varying geometry would need, so it is honestly part of "running GJK/EPA on this
//! input", not avoidable overhead. Neither path's loop allocates per iteration.
//!
//! Usage: ANTIDUNES_PAIRS_DIR=<dir with box_pairs.csv/cyl_pairs.csv, 100000-trial dump>
//!   ./ANTIDUNES_ContactSpeed100kBenchmark
//======================================================================================================================

#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"
#include "ClosedFormContacts.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{
namespace benchmark100k
{
using Vec3 = mesa_pd::Vec3;
using Mat3 = mesa_pd::Mat3;
using Rot3 = mesa_pd::Rot3;

//! Duplicated from CylinderCylinderContactBenchmark.cpp's own copy (see that file for the derivation).
Rot3 rot3FromMatrix(const Mat3& m)
{
   const real_t w = std::sqrt(std::max(real_t(0), real_t(1) + m(0, 0) + m(1, 1) + m(2, 2))) / real_t(2);
   real_t       x = std::sqrt(std::max(real_t(0), real_t(1) + m(0, 0) - m(1, 1) - m(2, 2))) / real_t(2);
   real_t       y = std::sqrt(std::max(real_t(0), real_t(1) - m(0, 0) + m(1, 1) - m(2, 2))) / real_t(2);
   real_t       z = std::sqrt(std::max(real_t(0), real_t(1) - m(0, 0) - m(1, 1) + m(2, 2))) / real_t(2);
   x              = std::copysign(x, m(2, 1) - m(1, 2));
   y              = std::copysign(y, m(0, 2) - m(2, 0));
   z              = std::copysign(z, m(1, 0) - m(0, 1));
   return Rot3(mesa_pd::Quat(w, x, y, z));
}

//======================================================================================================================
// Minimal CSV reader: header-indexed, double-valued columns only (every column in box_pairs.csv /
// cyl_pairs.csv is numeric). Not timed -- runs entirely before any benchmark loop starts.
//======================================================================================================================
struct Csv
{
   std::unordered_map< std::string, size_t > col;
   std::vector< std::vector< double > >      rows;
   double operator()(size_t row, const std::string& name) const { return rows[row][col.at(name)]; }
};

Csv readCsv(const std::string& path)
{
   Csv           csv;
   std::ifstream f(path);
   WALBERLA_CHECK(f.is_open(), "ContactSpeed100kBenchmark: cannot open " << path);
   std::string header;
   std::getline(f, header);
   {
      std::stringstream ss(header);
      std::string       tok;
      size_t            i = 0;
      while (std::getline(ss, tok, ',')) csv.col[tok] = i++;
   }
   std::string line;
   while (std::getline(f, line))
   {
      std::vector< double > row;
      row.reserve(csv.col.size());
      std::stringstream ss(line);
      std::string       tok;
      while (std::getline(ss, tok, ','))
      {
         // cyl_pairs.csv's cf_scenario column is a string (e.g. "band-band-parallel", "none"), not a
         // number; every column this benchmark actually reads by name is numeric, so a non-numeric
         // token is always an unused column -- stored as NaN rather than aborting the parse.
         try
         {
            row.push_back(std::stod(tok));
         }
         catch (const std::invalid_argument&)
         {
            row.push_back(std::numeric_limits< double >::quiet_NaN());
         }
      }
      csv.rows.push_back(std::move(row));
   }
   return csv;
}

Mat3 matFrom(const Csv& csv, size_t row, const std::string& prefix)
{
   auto g = [&](int i, int j) { return real_t(csv(row, prefix + std::to_string(i) + std::to_string(j))); };
   return Mat3(g(0, 0), g(0, 1), g(0, 2), g(1, 0), g(1, 1), g(1, 2), g(2, 0), g(2, 1), g(2, 2));
}

//======================================================================================================================
// Timing harness (same style as CylinderCylinderContactBenchmark.cpp: steady_clock around the whole
// loop, a warm-up pass excluded from timing, and a volatile sink so the optimizer cannot elide the call).
//======================================================================================================================
volatile const void* g_optimizationSink = nullptr;

template< typename Fn >
double timeItVarying(uint64_t n, uint64_t warmup, Fn&& fn)
{
   for (uint64_t i = 0; i < warmup; ++i)
   {
      auto result = fn(i % n);
      g_optimizationSink = &result;
   }
   const auto start = std::chrono::steady_clock::now();
   for (uint64_t i = 0; i < n; ++i)
   {
      auto result = fn(i);
      g_optimizationSink = &result;
   }
   const auto end = std::chrono::steady_clock::now();
   return double(std::chrono::duration_cast< std::chrono::nanoseconds >(end - start).count()) / double(n);
}

} // namespace benchmark100k
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

   using namespace walberla::antidunes::closed_form_contacts;
   using namespace walberla::antidunes::closed_form_contacts::benchmark100k;

   const char* dirEnv = std::getenv("ANTIDUNES_PAIRS_DIR");
   WALBERLA_CHECK(dirEnv != nullptr, "set ANTIDUNES_PAIRS_DIR to the directory holding the 100000-trial "
                                        "box_pairs.csv/cyl_pairs.csv (ANTIDUNES_REFERENCE_DUMP_ONLY=1 rerun)");
   const std::string dir = dirEnv;

   constexpr uint64_t warmup = 1000;

   // ---- Box x Box ------------------------------------------------------------------------------------------------
   {
      const Csv    csv = readCsv(dir + "/box_pairs.csv");
      const size_t n   = csv.rows.size();
      WALBERLA_LOG_INFO_ON_ROOT("box_pairs.csv: " << n << " trials");

      std::vector< Mat3 > RA(n), RB(n);
      std::vector< Vec3 > halfA(n), halfB(n), posB(n);
      for (size_t i = 0; i < n; ++i)
      {
         RA[i]    = matFrom(csv, i, "Aa");
         RB[i]    = matFrom(csv, i, "Ab");
         halfA[i] = Vec3(real_t(csv(i, "hax")), real_t(csv(i, "hay")), real_t(csv(i, "haz")));
         halfB[i] = Vec3(real_t(csv(i, "hbx")), real_t(csv(i, "hby")), real_t(csv(i, "hbz")));
         posB[i]  = Vec3(real_t(csv(i, "xb")), real_t(csv(i, "yb")), real_t(csv(i, "zb")));
      }
      const Vec3 posA(0, 0, 0);

      const double satNs = timeItVarying(n, warmup, [&](size_t i) {
         return boxBoxSAT(posA, RA[i], halfA[i], posB[i], RB[i], halfB[i]);
      });

      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto body0 = ps->create();
      body0->setPosition(posA);
      auto body1 = ps->create();
      walberla::antidunes::SolidCylinderContactDetection gjk;
      const double gjkNs = timeItVarying(n, warmup, [&](size_t i) {
         body0->setRotation(rot3FromMatrix(RA[i]));
         body1->setPosition(posB[i]);
         body1->setRotation(rot3FromMatrix(RB[i]));
         const mesa_pd::data::Box boxA(real_t(2) * halfA[i]), boxB(real_t(2) * halfB[i]);
         const bool contactFound = gjk(0, 1, boxA, boxB, ac);
         return std::make_pair(contactFound, gjk.getPenetrationDepth());
      });

      WALBERLA_LOG_INFO_ON_ROOT("");
      WALBERLA_LOG_INFO_ON_ROOT("=== Box x Box, average over " << n << " trials of the 100000-trial reference sample ===");
      WALBERLA_LOG_INFO_ON_ROOT("  SAT (closed form): " << satNs << " ns/call");
      WALBERLA_LOG_INFO_ON_ROOT("  GJK/EPA:           " << gjkNs << " ns/call  (" << (gjkNs / satNs) << "x SAT)");
   }

   // ---- Cylinder x Cylinder ---------------------------------------------------------------------------------------
   {
      const Csv    csv = readCsv(dir + "/cyl_pairs.csv");
      const size_t n   = csv.rows.size();
      WALBERLA_LOG_INFO_ON_ROOT("cyl_pairs.csv: " << n << " trials");

      std::vector< Mat3 >   RI(n), RJ(n);
      std::vector< Vec3 >   posJ(n);
      std::vector< real_t > Ri(n), Li(n), Rj(n), Lj(n);
      for (size_t i = 0; i < n; ++i)
      {
         RI[i]   = matFrom(csv, i, "Ai");
         RJ[i]   = matFrom(csv, i, "Aj");
         posJ[i] = Vec3(real_t(csv(i, "xj")), real_t(csv(i, "yj")), real_t(csv(i, "zj")));
         Ri[i]   = real_t(csv(i, "Ri"));
         Li[i]   = real_t(csv(i, "Li"));
         Rj[i]   = real_t(csv(i, "Rj"));
         Lj[i]   = real_t(csv(i, "Lj"));
      }
      const Vec3 posI(0, 0, 0);

      const double cfNs = timeItVarying(n, warmup, [&](size_t i) {
         const CylinderState ci{ posI, RI[i], Ri[i], Li[i] };
         const CylinderState cj{ posJ[i], RJ[i], Rj[i], Lj[i] };
         return cylinderCylinderContact(ci, cj); // production defaults: Ferrari, lazyThetaRef=true
      });

      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto body0 = ps->create();
      body0->setPosition(posI);
      body0->setInteractionRadius(real_t(1e6)); // disable the cheap center-distance reject (see
                                                 // CylinderCylinderContactBenchmark.cpp's own comment on
                                                 // the same guard); it would otherwise reject every pair
                                                 // here since the radius defaults to 0.
      auto body1 = ps->create();
      body1->setInteractionRadius(real_t(1e6));
      walberla::antidunes::SolidCylinderContactDetection gjk;
      const double gjkNs = timeItVarying(n, warmup, [&](size_t i) {
         body0->setRotation(rot3FromMatrix(RI[i]));
         body1->setPosition(posJ[i]);
         body1->setRotation(rot3FromMatrix(RJ[i]));
         const mesa_pd::data::Cylinder cylI(Ri[i], Li[i]), cylJ(Rj[i], Lj[i]);
         const bool contactFound = gjk(0, 1, cylI, cylJ, ac);
         return std::make_pair(contactFound, gjk.getPenetrationDepth());
      });

      // ---- Hybrid dispatch (production path: PaperAnalyticalContactTest.cpp's own
      // cylinderCylinderContactHybrid logic, reimplemented here to reuse the SAME pre-allocated
      // particle storage as the pure-GJK/EPA benchmark above instead of that test's own gjkCylCyl,
      // which allocates a fresh ParticleStorage every call -- fine for a one-off dump, not for
      // timing, since it would add allocation cost to every fallback and make the average less
      // representative of the real per-pair dispatch cost) -- closed form first, GJK/EPA only if
      // cylinderCylinderResultIsUncertain flags it, exactly the condition the production dispatcher
      // uses (ParticleKernels.h).
      // Fallback rate, computed separately from (and before) the timed loop below: counting it
      // inside the timed lambda would also count the warm-up pass's repeated calls over the first
      // `warmup` trials, double-counting some of them.
      uint64_t fallbackCount = 0;
      for (size_t i = 0; i < n; ++i)
      {
         const CylinderState ci{ posI, RI[i], Ri[i], Li[i] };
         const CylinderState cj{ posJ[i], RJ[i], Rj[i], Lj[i] };
         if (cylinderCylinderResultIsUncertain(ci, cj, cylinderCylinderContact(ci, cj))) ++fallbackCount;
      }

      const double hybridNs = timeItVarying(n, warmup, [&](size_t i) {
         const CylinderState ci{ posI, RI[i], Ri[i], Li[i] };
         const CylinderState cj{ posJ[i], RJ[i], Rj[i], Lj[i] };
         const ContactResult analytical = cylinderCylinderContact(ci, cj);
         if (!cylinderCylinderResultIsUncertain(ci, cj, analytical))
            return std::make_pair(analytical.contact, analytical.overlap);
         body0->setRotation(rot3FromMatrix(RI[i]));
         body1->setPosition(posJ[i]);
         body1->setRotation(rot3FromMatrix(RJ[i]));
         const mesa_pd::data::Cylinder cylI(Ri[i], Li[i]), cylJ(Rj[i], Lj[i]);
         const bool contactFound = gjk(0, 1, cylI, cylJ, ac);
         return std::make_pair(contactFound, -gjk.getPenetrationDepth());
      });

      // ---- Hybrid dispatch, "no scenario" trusted directly as no contact -- the alternative
      // hybrid-dispatch.tex's own prose names as "also viable" but does not use: special contact
      // still falls back to GJK/EPA (not the option under discussion there), but a pair for which
      // none of Kodam et al.'s seven scenarios applies is reported as not touching immediately,
      // without ever calling GJK/EPA to check. Correctness cost is NOT recomputed here -- it is
      // already established in cylinder-cylinder.tex: the closed form alone (no fallback at all,
      // the same "no scenario" cases this variant also never corrects) misses 166 of the 33439
      // contacts of this exact sample.
      uint64_t noScenarioFallbackCount = 0;
      for (size_t i = 0; i < n; ++i)
      {
         const CylinderState ci{ posI, RI[i], Ri[i], Li[i] };
         const CylinderState cj{ posJ[i], RJ[i], Rj[i], Lj[i] };
         const ContactResult analytical = cylinderCylinderContact(ci, cj);
         if (analytical.contact && std::string(analytical.scenario) == "special-contact") ++noScenarioFallbackCount;
      }

      const double noScenarioNs = timeItVarying(n, warmup, [&](size_t i) {
         const CylinderState ci{ posI, RI[i], Ri[i], Li[i] };
         const CylinderState cj{ posJ[i], RJ[i], Rj[i], Lj[i] };
         const ContactResult analytical = cylinderCylinderContact(ci, cj);
         if (!analytical.contact) return std::make_pair(false, real_t(0)); // trusted directly, no GJK/EPA call
         if (std::string(analytical.scenario) != "special-contact") return std::make_pair(analytical.contact, analytical.overlap);
         body0->setRotation(rot3FromMatrix(RI[i]));
         body1->setPosition(posJ[i]);
         body1->setRotation(rot3FromMatrix(RJ[i]));
         const mesa_pd::data::Cylinder cylI(Ri[i], Li[i]), cylJ(Rj[i], Lj[i]);
         const bool contactFound = gjk(0, 1, cylI, cylJ, ac);
         return std::make_pair(contactFound, -gjk.getPenetrationDepth());
      });

      WALBERLA_LOG_INFO_ON_ROOT("");
      WALBERLA_LOG_INFO_ON_ROOT("=== Cylinder x Cylinder, average over " << n << " trials of the 100000-trial reference sample ===");
      WALBERLA_LOG_INFO_ON_ROOT("  Pure closed form (production defaults): " << cfNs << " ns/call");
      WALBERLA_LOG_INFO_ON_ROOT("  Pure GJK/EPA:                           " << gjkNs << " ns/call  (" << (gjkNs / cfNs) << "x pure closed form)");
      WALBERLA_LOG_INFO_ON_ROOT("  Hybrid dispatch (production path):      " << hybridNs << " ns/call  (" << (hybridNs / cfNs)
                                 << "x pure closed form), fell back to GJK/EPA on " << fallbackCount << " of " << n
                                 << " trials (" << (100.0 * double(fallbackCount) / double(n)) << " %)");
      WALBERLA_LOG_INFO_ON_ROOT("  Hybrid, no-scenario trusted as no contact: " << noScenarioNs << " ns/call  (" << (noScenarioNs / cfNs)
                                 << "x pure closed form), fell back to GJK/EPA (special contact only) on "
                                 << noScenarioFallbackCount << " of " << n << " trials ("
                                 << (100.0 * double(noScenarioFallbackCount) / double(n)) << " %)");
   }

   return EXIT_SUCCESS;
}
