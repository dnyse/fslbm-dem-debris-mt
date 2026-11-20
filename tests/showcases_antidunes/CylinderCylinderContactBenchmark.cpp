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
//! \file CylinderCylinderContactBenchmark.cpp
//! Microbenchmark, not a correctness test (see PaperAnalyticalContactTest.cpp): times each of the 8 named
//! cylinder-cylinder closed-form sub-cases of ClosedFormContactsCylinderCylinder.h (face-face, band-band-parallel,
//! face-band, band-band-skewed, face-edge, special-contact, band-edge, edge-edge), the GJK/EPA fallback
//! (SolidCylinderContactDetection, the production code the fallback calls) and cylinderCylinderResultIsUncertain(),
//! each in a tight loop with fixed geometry, so only per-call algorithmic cost is measured (no broad phase, MPI or
//! allocation; ParticleStorage/bodies for the GJK cases are built once outside the timed loop). Box-box SAT and
//! box-box GJK/EPA are timed at the widerstacks box size (12x2x3 LU) for comparison.
//!
//! Motivation: in the widerstacks shape ablation a box run forced through GJK/EPA for every box-box pair
//! (useAnalyticalContactDetection=false) ran ~2x faster per step (~30 steps/s) than a cylinder run on the hybrid
//! path (closed form first, GJK/EPA as fallback; ~15.9 steps/s). The circumscribed radii (~6.27 vs ~6.16 LU) and the
//! sync algorithm are nearly the same, so the benchmark asks whether the cylinder-cylinder closed-form dispatch
//! itself, even when it avoids GJK, costs more per call than box-box SAT and box-box GJK/EPA.
//!
//! Uses ClosedFormContacts*.h and SolidCylinderContactDetection unmodified; GjkResult, rot3FromMatrix and the GJK/EPA
//! wrappers (gjkCylCylOnce/gjkBoxBoxOnce, reduced from gjkCylCyl/gjkBoxBox) are duplicated from
//! PaperAnalyticalContactTest.cpp, which keeps such helpers out of the headers so that ClosedFormContacts*.h does not
//! depend on ParticleKernels.h.
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
#include <core/math/Constants.h>
#include <core/mpi/MPIManager.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{
namespace benchmark
{
using Vec3 = mesa_pd::Vec3;
using Mat3 = mesa_pd::Mat3;
using Rot3 = mesa_pd::Rot3;

//======================================================================================================================
// Duplicated from PaperAnalyticalContactTest.cpp (see its comment on cylinderCylinderContactHybrid for why these
// live outside ClosedFormContacts*.h).
//======================================================================================================================
struct GjkResult
{
   bool   contact;
   real_t penetrationDepth;
};

//! Duplicated from PaperAnalyticalContactTest.cpp's rot3FromMatrix (component magnitudes from the diagonal, signs from
//! the off-diagonal differences).
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

GjkResult gjkCylCylOnce(walberla::antidunes::SolidCylinderContactDetection& gjk, mesa_pd::data::ParticleAccessor& ac,
                        const mesa_pd::data::Cylinder& cylI, const mesa_pd::data::Cylinder& cylJ)
{
   const bool contactFound = gjk(0, 1, cylI, cylJ, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth() };
}

GjkResult gjkBoxBoxOnce(walberla::antidunes::SolidCylinderContactDetection& gjk, mesa_pd::data::ParticleAccessor& ac,
                        const mesa_pd::data::Box& boxA, const mesa_pd::data::Box& boxB)
{
   const bool contactFound = gjk(0, 1, boxA, boxB, ac);
   return GjkResult{ contactFound, gjk.getPenetrationDepth() };
}

//======================================================================================================================
// Timing harness
//======================================================================================================================
struct TimingResult
{
   std::string label;
   double      nsPerCall;
   uint64_t    iterations;
};

// Every fn() calls an inline header function (ClosedFormContacts*.h / SolidCylinderContactDetection) with fixed
// inputs; its unused result would let the optimizer delete the call. Storing the result's address through a
// volatile pointer forces it to be computed every iteration. Correctness (scenario/contact flag) is checked by
// WALBERLA_CHECK* calls outside each timed block.
volatile const void* g_optimizationSink = nullptr;

template< typename Fn >
TimingResult timeIt(const std::string& label, uint64_t iterations, Fn&& fn)
{
   // Warm-up (cache effects, branch predictor) -- not timed.
   for (uint64_t i = 0; i < std::min< uint64_t >(iterations / 10, uint64_t(1000)); ++i)
   {
      auto result = fn();
      g_optimizationSink = &result;
   }

   const auto start = std::chrono::steady_clock::now();
   for (uint64_t i = 0; i < iterations; ++i)
   {
      auto result = fn();
      g_optimizationSink = &result;
   }
   const auto end = std::chrono::steady_clock::now();

   const double totalNs = double(std::chrono::duration_cast< std::chrono::nanoseconds >(end - start).count());
   return TimingResult{ label, totalNs / double(iterations), iterations };
}

void printSortedResults(std::vector< TimingResult > results)
{
   std::sort(results.begin(), results.end(),
             [](const TimingResult& a, const TimingResult& b) { return a.nsPerCall > b.nsPerCall; });

   WALBERLA_LOG_INFO_ON_ROOT("");
   WALBERLA_LOG_INFO_ON_ROOT("=== Results, slowest to fastest ===");
   const double fastest = results.back().nsPerCall;
   for (const auto& r : results)
   {
      WALBERLA_LOG_INFO_ON_ROOT("  " << std::string(38 - std::min< size_t >(38, r.label.size()), ' ') << r.label
                                     << ": " << r.nsPerCall << " ns/call"
                                     << "  (" << (r.nsPerCall / fastest) << "x the fastest entry)"
                                     << "  [" << r.iterations << " iterations]");
   }
   WALBERLA_LOG_INFO_ON_ROOT("");
}

} // namespace benchmark
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
   using namespace walberla::antidunes::closed_form_contacts::benchmark;

   constexpr uint64_t N = 500000;
   std::vector< TimingResult > results;

   // ---- Cylinder x Cylinder, closed-form dispatch only (cylinderCylinderContact), one geometry per named sub-case,
   // copied from PaperAnalyticalContactTest.cpp's testCylinderCylinder* functions (validated to hit that branch).
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.5)), Mat3(), real_t(0.5), real_t(2) };
      const auto r = timeIt("face-face", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "face-face");
      results.push_back(r);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      CylinderState cj{ Vec3(real_t(0), real_t(0.8), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      const auto r = timeIt("band-band-parallel", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "band-band-parallel");
      results.push_back(r);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
      CylinderState cj{ Vec3(real_t(0), real_t(0.45), real_t(1.4)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
      const auto r = timeIt("face-band", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "face-band");
      results.push_back(r);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(4) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi * real_t(0.5));
      CylinderState cj{ Vec3(real_t(0), real_t(0.8), real_t(1.5)), rotJ.getMatrix(), real_t(0.5), real_t(4) };
      const auto r = timeIt("band-band-skewed", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "band-band-skewed");
      results.push_back(r);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(3.0), real_t(1.0) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(6));
      CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.0)), rotJ.getMatrix(), real_t(0.3), real_t(1.5) };
      const auto r = timeIt("face-edge", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "face-edge");
      results.push_back(r);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(1.0), real_t(3.0) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * (real_t(35) / real_t(180)));
      CylinderState cj{ Vec3(real_t(0.9), real_t(0), real_t(1.7)), rotJ.getMatrix(), real_t(1.0), real_t(0.6) };
      const auto r = timeIt("special-contact", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "special-contact");
      results.push_back(r);
   }
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(3.0), real_t(6.0) };
      Rot3          rotJ;
      rotJ.rotate(Vec3(real_t(0), real_t(1), real_t(0)), math::pi / real_t(3));
      CylinderState cj{ Vec3(real_t(3.3), real_t(0), real_t(1.0)), rotJ.getMatrix(), real_t(0.3), real_t(1.5) };
      const auto r = timeIt("band-edge", N, [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "band-edge");
      results.push_back(r);
   }
   {
      // geometry of testCylinderCylinderEdgeEdgeGenuine (see there for the matrices' provenance)
      const Mat3 rotICols(real_t(-0.390898), real_t(-0.198955), real_t(0.898675), real_t(0.196307),
                           real_t(0.935875), real_t(0.292579), real_t(-0.899257), real_t(0.290784),
                           real_t(-0.326775));
      const Mat3 rotJCols(real_t(0.131414), real_t(-0.370814), real_t(0.919362), real_t(-0.959446),
                           real_t(0.18573), real_t(0.212056), real_t(-0.249386), real_t(-0.909946),
                           real_t(-0.331369));
      const Rot3    rotI = rot3FromMatrix(rotICols);
      const Rot3    rotJ = rot3FromMatrix(rotJCols);
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), rotI.getMatrix(), real_t(1.44446), real_t(1.30996) };
      CylinderState cj{ Vec3(real_t(2.11131), real_t(0.251881), real_t(0.154264)), rotJ.getMatrix(), real_t(1.22407),
                         real_t(2.5082) };
      // baseline: Durand-Kerner, lazyThetaRef=false, passed explicitly because the defaults are Ferrari/true
      const auto r = timeIt("edge-edge (original: old ref, Durand-Kerner, eager)", N, [&]() {
         return cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::DurandKerner, false);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::DurandKerner, false).scenario),
         "edge-edge");
      results.push_back(r);

      // useFastEdgeEdgeReference=true with the baseline solver/laziness, isolating the fast-reference effect. Read
      // together with findRimCrossing's comment (ClosedFormContactsCylinderPlaneGeometry.h) and the accuracy
      // comparison in PaperAnalyticalContactTest.cpp's testEdgeEdgeThreeWayAccuracy10k.
      const auto rFast = timeIt("edge-edge (fast ref only, opt-in)", N, [&]() {
         return cylinderCylinderContact(ci, cj, true, 30, detail::QuarticSolverMethod::DurandKerner, false);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, true, 30, detail::QuarticSolverMethod::DurandKerner, false).scenario),
         "edge-edge");
      results.push_back(rFast);

      // current default (Ferrari + lazyThetaRef), called as production does: cylinderCylinderContact(ci, cj)
      const auto rProductionDefault =
         timeIt("edge-edge (CURRENT file default: Ferrari + lazyThetaRef)", N,
                [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "edge-edge");
      results.push_back(rProductionDefault);

      // Solver variants: same geometry, line-based reference and 30 iterations, quartic solver swapped. lazyThetaRef is
      // not passed and takes its default (true), so the Ferrari row is the same call as the production-default and
      // "Ferrari + lazyThetaRef" rows. QuarticSolverComparisonTest.cpp has the synthetic-quartic comparison.
      const auto rFerrari = timeIt("edge-edge (Ferrari solver, opt-in)", N, [&]() {
         return cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari).scenario),
         "edge-edge");
      results.push_back(rFerrari);

      const auto rAberthEhrlich = timeIt("edge-edge (Aberth-Ehrlich solver, opt-in)", N, [&]() {
         return cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::AberthEhrlich);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::AberthEhrlich).scenario),
         "edge-edge");
      results.push_back(rAberthEhrlich);

      // lazyThetaRef: Ferrari, but thetaRef's quartic is computed only if point A's or B's containment check gives
      // more than one candidate (never in testThetaRefTieBreakNecessity10k's fuzz corpus, 0/832 calls).
      // lazyThetaRef changes only when thetaRef is computed, not which reference.
      const auto rLazyThetaRef = timeIt("edge-edge (Ferrari + lazyThetaRef, opt-in)", N, [&]() {
         return cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, true);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, true).scenario),
         "edge-edge");
      results.push_back(rLazyThetaRef);
   }

   // ---- cylinderCylinderResultIsUncertain() alone: the gate that decides whether a closed-form result is trusted
   // or GJK/EPA is called.
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.5)), Mat3(), real_t(0.5), real_t(2) };
      const ContactResult res = cylinderCylinderContact(ci, cj);
      const auto r = timeIt("cylinderCylinderResultIsUncertain", N,
                             [&]() { return cylinderCylinderResultIsUncertain(ci, cj, res); });
      results.push_back(r);
   }

   // ---- Cylinder x Cylinder, GJK/EPA fallback (SolidCylinderContactDetection, as called by the hybrid dispatch),
   // geometry of PaperAnalyticalContactTest.cpp's testCylinderCylinderTiltedRimIntoBand (exactly touching, closed form
   // finds nothing). Bodies built once.
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto body0 = ps->create();
      body0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
      body0->setRotation(Rot3());
      body0->setInteractionRadius(std::sqrt(real_t(0.5) * real_t(0.5) + real_t(0.25) * real_t(4) * real_t(4)));
      auto body1 = ps->create();
      body1->setPosition(Vec3(real_t(0.75), real_t(0), real_t(1.2)));
      Rot3 rotJ;
      rotJ.rotate(Vec3(real_t(1), real_t(0), real_t(0)), math::pi * real_t(0.3));
      body1->setRotation(rotJ);
      body1->setInteractionRadius(std::sqrt(real_t(0.25) * real_t(0.25) + real_t(0.25) * real_t(1.5) * real_t(1.5)));
      const mesa_pd::data::Cylinder cylI(real_t(0.5), real_t(4)), cylJ(real_t(0.25), real_t(1.5));

      walberla::antidunes::SolidCylinderContactDetection gjk;
      const auto r = timeIt("cylinder-cylinder GJK/EPA fallback", N,
                             [&]() { return gjkCylCylOnce(gjk, ac, cylI, cylJ); });
      WALBERLA_CHECK(gjkCylCylOnce(gjk, ac, cylI, cylJ).contact);
      results.push_back(r);
   }

   // ---- Box x Box, SAT, widerstacks box size (halfExtents <6,1,1.5>, edges <12,2,3>), face-to-face overlap along x.
   {
      const Vec3 posA(real_t(0), real_t(0), real_t(0));
      const Vec3 posB(real_t(11), real_t(0), real_t(0)); // 1 LU overlap along x (edge length 12 each)
      const Mat3 R;
      const Vec3 half(real_t(6), real_t(1), real_t(1.5));
      const auto r = timeIt("box-box SAT", N, [&]() { return boxBoxSAT(posA, R, half, posB, R, half); });
      WALBERLA_CHECK(boxBoxSAT(posA, R, half, posB, R, half).contact);
      results.push_back(r);
   }

   // ---- Box x Box, GJK/EPA (Box overload of SolidCylinderContactDetection), same geometry as the SAT case.
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
      mesa_pd::data::ParticleAccessor ac(ps);
      auto body0 = ps->create();
      body0->setPosition(Vec3(real_t(0), real_t(0), real_t(0)));
      body0->setRotation(Rot3());
      auto body1 = ps->create();
      body1->setPosition(Vec3(real_t(11), real_t(0), real_t(0)));
      body1->setRotation(Rot3());
      const mesa_pd::data::Box boxA(Vec3(real_t(12), real_t(2), real_t(3))),
         boxB(Vec3(real_t(12), real_t(2), real_t(3)));

      walberla::antidunes::SolidCylinderContactDetection gjk;
      const auto r =
         timeIt("box-box GJK/EPA", N, [&]() { return gjkBoxBoxOnce(gjk, ac, boxA, boxB); });
      WALBERLA_CHECK(gjkBoxBoxOnce(gjk, ac, boxA, boxB).contact);
      results.push_back(r);
   }

   printSortedResults(results);

   return EXIT_SUCCESS;
}
