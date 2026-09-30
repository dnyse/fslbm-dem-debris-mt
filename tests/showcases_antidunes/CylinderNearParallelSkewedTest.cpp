//======================================================================================================================
//! \file CylinderNearParallelSkewedTest.cpp
//! Regression test for the false band-band-skewed contacts that aborted the cylinder hopper ([RELEASE-OVERLAP], job 4248725).
//!
//! Kodam et al. (2010), Part I, Section 3.2.4.2 (Eqs. 55-66): for two cylinders with SKEWED axes the shortest line between the two infinite
//! centrelines has length d = |x_cyl_j^i . r_hat| with r_hat = (e_i x e_j)/|e_i x e_j| (Eqs. 55-56) and touches the centrelines at
//! l_i*, l_j* (Eqs. 59-60, Vega & Lago 1994). Contact requires ALL of  d < R_i+R_j,  |l_i*| < L_i/2,  |l_j*| < L_j/2  (Eq. 63).
//!
//! For NEARLY parallel axes (tilt theta ~ 1e-6 rad, i.e. ordinary floating-point rotation noise of a body that was placed axis-aligned)
//! the lines cross theta-close to each other only very FAR along the axes (|l*| ~ separation/theta ~ 1e6 LU), so Eq. 63's l* conditions
//! correctly REJECT. The implementation (PaperAnalyticalContacts.h, band-band-skewed) used to replace l* by 0 whenever
//! 1-(e_i.e_j)^2 < 1e-10 ("guard against division by ~0") -- which makes both l* conditions vacuously true, while d = |x . r_hat| with an
//! arbitrary noise direction r_hat can be as small as 0. Result: two rods clearly SEPARATED by a gap were reported as overlapping by up
//! to R_i+R_j = 6 LU. This test
//!   (1) implements Kodam Eqs. 55-66 literally (with the exact l* of Eq. 59, no guard) as an independent reference,
//!   (2) shows the reference rejects the separated near-parallel pairs (so the PAPER is fine, the guard was the bug),
//!   (3) requires paper_analytical::cylinderCylinderContact and the production dispatch (analytic on, tight layer-1 on/off) to agree,
//!   (4) requires genuinely overlapping near-parallel pairs to still be detected with the right overlap,
//!   (5) requires the dispatch's capsule bound (a contact is impossible if the axis SEGMENTS are farther apart than R_i+R_j).
//======================================================================================================================
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <cmath>
#include <iostream>
#include <random>

namespace walberla { namespace antidunes {
using Vec3 = mesa_pd::Vec3;

constexpr real_t kR = real_t(3), kL = real_t(18);

mesa_pd::Rot3 rotToAxis(const Vec3& axis)
{
   mesa_pd::Rot3 rot; const Vec3 ez(real_t(0), real_t(0), real_t(1));
   const real_t cosA = ez * axis;
   if (cosA < real_t(1) - real_t(1e-14))
   {
      const Vec3 rotAxis = (cosA > real_t(-1) + real_t(1e-14)) ? (ez % axis).getNormalized() : Vec3(real_t(1), real_t(0), real_t(0));
      rot.rotate(rotAxis, std::acos(std::clamp(cosA, real_t(-1), real_t(1))));
   }
   return rot;
}

struct Rod { Vec3 pos; Vec3 axis; };  // axis = unit vector (world) of the cylinder's symmetry axis

paper_analytical::CylinderState state(const Rod& r)
{
   return paper_analytical::CylinderState{ r.pos, rotToAxis(r.axis).getMatrix(), kR, kL };
}

//! Kodam et al. (2010) Eqs. 55-66, literal, in world vectors; NO guard on the 1-(e_i.e_j)^2 denominator.
//! Returns overlap (R_i+R_j-d, Eq. 64) if contact per Eq. 63, otherwise a negative number; also outputs the l* values.
double kodamBandBandSkewed(const Rod& i, const Rod& j, double& liStar, double& ljStar, double& d)
{
   const Vec3 ei = i.axis, ej = j.axis, x = j.pos - i.pos;        // x = x_cyl_j^i (cylinder j's position seen from i)
   liStar = ljStar = d = 0.0;
   const Vec3 cr = ei % ej; const double crLen = cr.length();
   if (crLen < 1e-15) return -1.0;                                  // exactly parallel: Eq. 55 undefined (the parallel branch, Eq. 52, applies instead)
   const Vec3 rHat = cr / real_t(crLen);                            // Eq. 55 (manual normalisation: getNormalized() asserts |v| > ~1e-8 in debug builds)
   d = std::abs(x * rHat);                                         // Eq. 56
   const double c = ei * ej;
   ljStar = ((c * (ei * x)) - (ej * x)) / (1.0 - c * c);           // Eq. 59
   liStar = ei * (x + ej * ljStar);                                // Eq. 60
   const bool contact = d < double(2 * kR) && std::abs(liStar) < 0.5 * double(kL) && std::abs(ljStar) < 0.5 * double(kL); // Eq. 63
   return contact ? double(2 * kR) - d : -1.0;                     // Eq. 64
}

//! Two Y-axis rods (as in the hopper), the second at (-centreDistance) along X, tilted by `tilt` rad about `tiltAxis`.
Rod tiltedRod(const Vec3& pos, const Vec3& tiltAxis, double tilt)
{
   mesa_pd::Rot3 r; r.rotate(tiltAxis, real_t(tilt));
   return Rod{ pos, r.getMatrix() * Vec3(real_t(0), real_t(1), real_t(0)) };
}

struct Outcome { bool kodam; double kodamOverlap; bool paper; double paperOverlap; bool dispatch; bool dispatchTight; double dispatchOverlap; const char* scenario; };

Outcome evaluate(const Rod& a, const Rod& b)
{
   Outcome o{};
   double li, lj, d;
   o.kodamOverlap = kodamBandBandSkewed(a, b, li, lj, d);
   o.kodam = o.kodamOverlap >= 0;
   const auto res = paper_analytical::cylinderCylinderContact(state(a), state(b));
   o.paper = res.contact; o.paperOverlap = res.contact ? double(res.overlap) : -1.0; o.scenario = res.scenario;

   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   for (const Rod* r : { &a, &b }) { auto p = ps->create(); p->setPosition(r->pos); p->setRotation(rotToAxis(r->axis)); }
   const mesa_pd::data::Cylinder shape(kR, kL);
   SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = true;
   o.dispatch = scd(0, 1, shape, shape, ac); o.dispatchOverlap = o.dispatch ? -double(scd.getPenetrationDepth()) : -1.0;
   SolidCylinderContactDetection scdT; scdT.useAnalyticalContactDetection = true; scdT.useTightLayer1 = true;
   o.dispatchTight = scdT(0, 1, shape, shape, ac);
   return o;
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::antidunes;
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };

   const Vec3 Xa(1, 0, 0), Za(0, 0, 1);
   const Vec3 base(real_t(249), real_t(30), real_t(49.13));

   // ---- (A) deterministic minimal reproducer: tilt about Z, separation along X => r_hat = +-Z is PERPENDICULAR to the separation, d = 0 ----
   for (double tilt : { 1e-7, 1e-6, 2e-6, 5e-6, 1e-5, 3e-5, 1e-4, 1e-3 })
      for (double gap : { 0.05, 0.1, 0.5, 2.0 })
      {
         const Rod a{ base, Vec3(0, 1, 0) };
         const Rod b = tiltedRod(base - Vec3(real_t(2 * kR + gap), 0, 0), Za, tilt);
         const auto o = evaluate(a, b);
         std::ostringstream tag; tag << "A tilt=" << tilt << " gap=" << gap;
         if (tilt == 2e-6 && gap == 0.1)
            std::cout << "reproducer (tilt 2e-6 about Z, gap 0.1): Kodam literal contact=" << o.kodam << "  paper_analytical contact=" << o.paper << " (" << o.scenario << ", overlap " << o.paperOverlap
                      << ")  dispatch=" << o.dispatch << " dispatch(tightLayer1)=" << o.dispatchTight << std::endl;
         check(!o.kodam,          tag.str() + ": Kodam Eq. 63 with the exact l* must REJECT (rods are separated)");
         check(!o.paper,          tag.str() + ": paper_analytical::cylinderCylinderContact reports a contact between SEPARATED rods (" + o.scenario + ", overlap " + std::to_string(o.paperOverlap) + ")");
         check(!o.dispatch,       tag.str() + ": production dispatch (analytic) reports a contact between SEPARATED rods");
         check(!o.dispatchTight,  tag.str() + ": production dispatch (analytic + useTightLayer1) reports a contact between SEPARATED rods");
      }

   // ---- (B) Monte Carlo: random tilt axis/magnitude, random position noise, separated by gap > 0 (any direction in the XZ plane) ----
   {
      std::mt19937 gen(4711); std::normal_distribution< double > nd(0.0, 1.0); std::uniform_real_distribution< double > ud(0.0, 6.283185307179586);
      int nFalse = 0, nTrials = 0;
      for (double tilt : { 1e-7, 1e-6, 1e-5, 1e-4 })
         for (double gap : { 0.02, 0.1, 1.0 })
            for (int t = 0; t < 400; ++t)
            {
               const double phi = ud(gen);                               // separation direction in the XZ plane
               const Vec3 sep(real_t(std::cos(phi)), 0, real_t(std::sin(phi)));
               const Vec3 tiltAxisA(real_t(nd(gen)), real_t(nd(gen)), real_t(nd(gen))), tiltAxisB(real_t(nd(gen)), real_t(nd(gen)), real_t(nd(gen)));
               const Rod a = tiltedRod(base, tiltAxisA.getNormalized(), tilt * nd(gen));
               const Rod b = tiltedRod(base + sep * real_t(2 * kR + gap + 1e-3), tiltAxisB.getNormalized(), tilt * nd(gen));
               // gap measured between the true (axis-aligned) rods; the sub-1e-3 tilt effect on the surface distance is < 9*1e-4 = 9e-4 < 1e-3 margin above
               const auto o = evaluate(a, b);
               ++nTrials;
               if (o.paper || o.dispatch || o.dispatchTight || o.kodam) ++nFalse;
            }
      std::cout << "Monte Carlo: separated near-parallel rod pairs " << nTrials << ", false contacts (any implementation incl. Kodam literal) " << nFalse << std::endl;
      check(nFalse == 0, "Monte Carlo: " + std::to_string(nFalse) + " of " + std::to_string(nTrials) + " separated near-parallel pairs reported as touching");
   }

   // ---- (C) genuinely overlapping near-parallel rods must STILL be detected, with the right overlap (= 2R - centre distance) ----
   for (double tilt : { 0.0, 1e-7, 1e-6, 1e-5, 1e-4 })
      for (double overlap : { 0.1, 0.5, 1.0 })
      {
         const Rod a{ base, Vec3(0, 1, 0) };
         const Rod b = tiltedRod(base - Vec3(real_t(2 * kR - overlap), 0, 0), Za, tilt);
         const auto o = evaluate(a, b);
         std::ostringstream tag; tag << "C tilt=" << tilt << " overlap=" << overlap;
         check(o.paper && o.dispatch && o.dispatchTight, tag.str() + ": overlapping rods NOT detected (paper=" + std::to_string(o.paper) + " dispatch=" + std::to_string(o.dispatch) + " tight=" + std::to_string(o.dispatchTight) + ")");
         check(o.dispatch && std::abs(o.dispatchOverlap - overlap) < 1e-3, tag.str() + ": dispatch overlap " + std::to_string(o.dispatchOverlap) + " != " + std::to_string(overlap));
      }

   // ---- (D) capsule bound: the two axis SEGMENTS are 40 LU apart, so NO analytic result may be accepted (defence in depth against ANY future
   //          dispatch bug), checked through the dispatch with several tilts/orientations ----
   for (double tilt : { 0.0, 1e-6, 1e-3, 0.3 })
   {
      const Rod a{ base, Vec3(0, 1, 0) };
      const Rod b = tiltedRod(base - Vec3(40, 0, 0), Za, tilt);
      const auto o = evaluate(a, b);
      check(!o.dispatch && !o.dispatchTight && !o.paper, "D: rods 34 LU apart reported as touching, tilt=" + std::to_string(tilt));
   }

   std::cout << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failures)" << std::endl;
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
