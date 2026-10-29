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
//! Microbenchmark, NOT a correctness test (see PaperAnalyticalContactTest.cpp for that): times each
//! of PaperAnalyticalContacts.h's 8 named cylinder-cylinder analytical sub-cases (face-face,
//! band-band-parallel, face-band, band-band-skewed, face-edge, special-contact, band-edge, edge-edge),
//! the GJK/EPA fallback path (SolidCylinderContactDetection, the same production code the fallback
//! actually calls), and cylinderCylinderResultIsUncertain() itself, in a tight repeated-call loop with
//! FIXED geometry per case (so what's measured is pure per-call algorithmic cost, not broad-phase,
//! MPI, or particle-storage/allocation overhead -- ParticleStorage/bodies for the GJK cases are built
//! ONCE outside the timed loop). Also benchmarks box-box SAT and box-box GJK/EPA at the widerstacks
//! scenario's own box size (12x2x3 LU) for a same-methodology comparison against the cylinder numbers.
//!
//! WHY THIS EXISTS: this session's own widerstacks-shape-ablation runs found something the user
//! called "crazy" -- a box-version run FORCED through GJK/EPA for every Box x Box pair (job 4192137,
//! useAnalyticalContactDetection=false) ran ~2x FASTER per step (steps 40000-45000: ~30 steps/sec)
//! than a cylinder-version run using the STANDARD hybrid path (analytical-first, GJK/EPA only as a
//! fallback; job 4192276, same step range: ~15.9 steps/sec) -- i.e. the path that's supposed to avoid
//! the "expensive" iterative solver most of the time was still slower than the path that always pays
//! for it. Geometry classification alone (SetupStages.h's blockSyncExtension calc) doesn't explain
//! this: box and cylinder have almost identical circumscribed radii for this scenario's piece sizes
//! (~6.27 vs ~6.16 LU) and use the same sync algorithm. This benchmark isolates the one remaining
//! candidate: is the cylinder-cylinder ANALYTICAL dispatch itself (not the fallback) the actual
//! bottleneck -- i.e. is classifying and closed-form-solving a cylinder-cylinder pair, even when it
//! successfully avoids GJK, more expensive per call than box-box SAT AND (surprisingly) more
//! expensive than box-box's own GJK/EPA path?
//!
//! NOT wired into the simulation and NOT a correctness check -- reuses (does not modify)
//! PaperAnalyticalContacts.h and ParticleKernels.h's SolidCylinderContactDetection exactly as
//! PaperAnalyticalContactTest.cpp does; gjkCylCyl/gjkBoxBox/HybridResult are duplicated here (not
//! shared with that file) for the same reason that file gives for defining them locally in the first
//! place: keeping PaperAnalyticalContacts.h free of any dependency on ParticleKernels.h.
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
// Duplicated from PaperAnalyticalContactTest.cpp (see that file's own comment on
// cylinderCylinderContactHybrid for why these live outside PaperAnalyticalContacts.h).
//======================================================================================================================
struct GjkResult
{
   bool   contact;
   real_t penetrationDepth;
};

//! Duplicated from PaperAnalyticalContactTest.cpp's own rot3FromMatrix -- see that file's copy for
//! the derivation (standard matrix-to-quaternion conversion, largest-diagonal-term-safe variant).
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

// Every fn() below calls an `inline` header function (PaperAnalyticalContacts.h /
// SolidCylinderContactDetection) with fixed inputs and returns a value that would otherwise go
// unused -- an optimizer is entitled to prove the whole call has no observable effect and delete it
// outright, which would silently turn a timed loop into an empty one. Storing each result's address
// through a `volatile` pointer forces the compiler to actually materialize it every iteration,
// without needing a type-erased result or a benchmarking library. Correctness (do the analytical
// scenario/contact-flag checks below still pass) is verified by the WALBERLA_CHECK_EQUAL/WALBERLA_CHECK
// calls placed OUTSIDE each timed block -- this sink exists purely to make the timing honest.
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

   // ---- Cylinder x Cylinder, analytical dispatch only (cylinderCylinderContact), one geometry per
   // named sub-case -- EXACT geometries copied from PaperAnalyticalContactTest.cpp's own
   // testCylinderCylinder* functions, so these are known-good, already-validated inputs for each case,
   // not hand-picked-for-this-benchmark geometry that might accidentally miss the intended branch.
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
      // testCylinderCylinderEdgeEdgeGenuine's geometry -- see that test for provenance of the matrices.
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
      // "original" pinned EXPLICITLY to the literal pre-this-session baseline (DurandKerner,
      // lazyThetaRef=false) -- the file's own defaults for quarticSolver/lazyThetaRef changed to
      // Ferrari/true later in this session, so leaving this call on cylinderCylinderContact(ci,cj,false)
      // (2-arg, relying on defaults) would have silently turned this into a duplicate of the Ferrari +
      // lazyThetaRef entry below instead of the baseline it's meant to be.
      const auto r = timeIt("edge-edge (original: old ref, Durand-Kerner, eager)", N, [&]() {
         return cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::DurandKerner, false);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::DurandKerner, false).scenario),
         "edge-edge");
      results.push_back(r);

      // Simplification #1 (this session): same geometry, useFastEdgeEdgeReference=true, solver/laziness
      // pinned to the ORIGINAL Durand-Kerner/eager settings so this isolates ONLY the fast-reference
      // effect -- see PaperAnalyticalContacts.h's findRimCrossing's own comment and
      // PaperAnalyticalContactTest.cpp's testEdgeEdgeFastReferenceAccuracy10k for the accuracy
      // comparison this speed number should always be read alongside.
      const auto rFast = timeIt("edge-edge (fast ref only, opt-in)", N, [&]() {
         return cylinderCylinderContact(ci, cj, true, 30, detail::QuarticSolverMethod::DurandKerner, false);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, true, 30, detail::QuarticSolverMethod::DurandKerner, false).scenario),
         "edge-edge");
      results.push_back(rFast);

      // The file's ACTUAL current default (this session's end state): Ferrari + lazyThetaRef, called
      // exactly as production would (no arguments beyond ci/cj) -- this is what a plain
      // cylinderCylinderContact(ci, cj) call anywhere else in this codebase actually runs today.
      const auto rProductionDefault =
         timeIt("edge-edge (CURRENT file default: Ferrari + lazyThetaRef)", N,
                [&]() { return cylinderCylinderContact(ci, cj); });
      WALBERLA_CHECK_EQUAL(std::string(cylinderCylinderContact(ci, cj).scenario), "edge-edge");
      results.push_back(rProductionDefault);

      // Solver-choice variants (this session's Ferrari/Aberth-Ehrlich addition) -- same edge-edge
      // geometry, old (line-based) reference and 30-iteration budget held fixed, only the quartic
      // solver itself swapped. See QuarticSolverComparisonTest.cpp for the synthetic-quartic
      // time/accuracy comparison that motivated checking real geometry too.
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

      // lazyThetaRef (this session's follow-up) -- Ferrari solver, but thetaRef's own quartic is only
      // computed if point A or B's containment check ever produces more than one candidate, which
      // testThetaRefTieBreakNecessity10k found is NEVER true on real geometry (0/832 calls). Same
      // geometry as production default (Ferrari, 30 iter, old useFastEdgeEdgeReference semantics --
      // lazyThetaRef only changes WHEN thetaRef is computed, not which reference it computes).
      const auto rLazyThetaRef = timeIt("edge-edge (Ferrari + lazyThetaRef, opt-in)", N, [&]() {
         return cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, true);
      });
      WALBERLA_CHECK_EQUAL(
         std::string(
            cylinderCylinderContact(ci, cj, false, 30, detail::QuarticSolverMethod::Ferrari, true).scenario),
         "edge-edge");
      results.push_back(rLazyThetaRef);
   }

   // ---- cylinderCylinderResultIsUncertain() alone -- the cheap classification gate every analytical
   // result passes through to decide whether to trust it or fall back to GJK/EPA.
   {
      CylinderState ci{ Vec3(real_t(0), real_t(0), real_t(0)), Mat3(), real_t(0.5), real_t(2) };
      CylinderState cj{ Vec3(real_t(0), real_t(0), real_t(1.5)), Mat3(), real_t(0.5), real_t(2) };
      const ContactResult res = cylinderCylinderContact(ci, cj);
      const auto r = timeIt("cylinderCylinderResultIsUncertain", N,
                             [&]() { return cylinderCylinderResultIsUncertain(ci, cj, res); });
      results.push_back(r);
   }

   // ---- Cylinder x Cylinder, GJK/EPA fallback path (production SolidCylinderContactDetection,
   // exactly what the hybrid dispatcher calls on a fallback) -- geometry from
   // testHybridTier1TriggersOnMissedContact, a real tier-1 (analytical-found-nothing) fallback case.
   // ParticleStorage/bodies built ONCE, outside the timed loop.
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

   // ---- Box x Box, SAT (analytical, same methodology as above) -- widerstacks scenario's own box
   // size (halfExtents<6,1,1.5> -> full edge <12,2,3>), two boxes overlapping face-to-face along x.
   {
      const Vec3 posA(real_t(0), real_t(0), real_t(0));
      const Vec3 posB(real_t(11), real_t(0), real_t(0)); // 1 LU overlap along x (edge length 12 each)
      const Mat3 R;
      const Vec3 half(real_t(6), real_t(1), real_t(1.5));
      const auto r = timeIt("box-box SAT", N, [&]() { return boxBoxSAT(posA, R, half, posB, R, half); });
      WALBERLA_CHECK(boxBoxSAT(posA, R, half, posB, R, half).contact);
      results.push_back(r);
   }

   // ---- Box x Box, GJK/EPA (production SolidCylinderContactDetection's Box overload) -- SAME
   // geometry as the SAT case above, so this is a direct apples-to-apples comparison of "GJK/EPA cost
   // for this scenario's own box size" against every cylinder-cylinder number above.
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
