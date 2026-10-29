//======================================================================================================================
//! \file CylinderClosedFormRangeTest.cpp
//! Characterisation test: for which overlap depths is the closed-form Kodam+Guo Cylinder x Cylinder result (closed_form_contacts::cylinderCylinderContact) a correct penetration depth?
//!
//! Background (THESIS_PLAN.md 2026-09-20): in the cylinder hopper rods pile up and interpenetrate by 30-90 % of their volume. There the closed-form dispatch returns e.g. face-edge overlap 0.001 LU or edge-edge
//! overlap 18 LU for a pair whose true penetration depth is ~5.9 LU (GJK/EPA agrees with the truth). Claim to prove here:
//!   (1) the exact penetration depth of two convex bodies is  PD = min_n [ h_A(n) + h_B(n) - d.n ]  (h = support function; h(n) = (L/2)|n.a| + R sqrt(1-(n.a)^2)); this is the independent ground truth,
//!   (2) the closed form is accurate for SHALLOW overlaps and wrong for deep ones (the range is measured below, per depth class),
//!   (3) a cheap self-consistency check on the closed form -- "how far must B move ALONG THE REPORTED NORMAL to separate?"  needed(n) = h_A(n)+h_B(n)-d.n (min over the two normal orientations) --
//!       exposes the wrong results: needed(n) >= PD always, so an UNDER-estimate `overlap < needed - tol` is always detectable, and for shallow contacts needed(n) == overlap;
//!   (4) GJK/EPA (analytic path off) is accurate at every depth.
//! The check in (3) is NOT used by production (no fallback is implemented yet); this test only shows that it would separate right from wrong.
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
constexpr double kR = 3.0, kL = 18.0, kH = 9.0;

inline double support(const Vec3& n, const Vec3& a)
{
   const double c = std::abs(double(n * a));
   return kH * c + kR * std::sqrt(std::max(0.0, 1.0 - c * c));
}
inline double needed(const Vec3& n, const Vec3& a1, const Vec3& a2, const Vec3& d) { return support(n, a1) + support(n, a2) - double(d * n); }

//! cheap UPPER bound on the exact depth: PD <= needed(n) for EVERY direction n, so the minimum over a few canonical candidate directions bounds it from above
//! (axes, centre line, common perpendicular, radial directions of each rod). An analytic overlap above this bound (+tol) is certainly an over-estimate.
double upperBound(const Vec3& a1, const Vec3& a2, const Vec3& d, const Vec3& reported)
{
   auto unit = [](const Vec3& v, Vec3& out) { const double l = v.length(); if (l < 1e-9) return false; out = v / real_t(l); return true; };
   std::vector< Vec3 > cand{ reported, a1, a2 };
   Vec3 u;
   if (unit(d, u)) cand.push_back(u);
   if (unit(a1 % a2, u)) cand.push_back(u);
   if (unit(d - a1 * (d * a1), u)) cand.push_back(u);
   if (unit(d - a2 * (d * a2), u)) cand.push_back(u);
   double ub = 1e300;
   for (const Vec3& n : cand) ub = std::min(ub, std::min(needed(n, a1, a2, d), needed(-n, a1, a2, d)));
   return ub;
}

//! ground truth: minimum over unit directions (Fibonacci grid, then pattern-search refinement of the best few)
double exactPenetrationDepth(const Vec3& a1, const Vec3& a2, const Vec3& d)
{
   constexpr int N = 6000;
   std::vector< std::pair< double, Vec3 > > c; c.reserve(N);
   const double ga = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
   for (int i = 0; i < N; ++i)
   {
      const double z = 1.0 - 2.0 * (i + 0.5) / N, r = std::sqrt(1.0 - z * z), ph = ga * i;
      const Vec3 n(real_t(r * std::cos(ph)), real_t(r * std::sin(ph)), real_t(z));
      c.emplace_back(needed(n, a1, a2, d), n);
   }
   std::partial_sort(c.begin(), c.begin() + 8, c.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
   double best = 1e300;
   std::mt19937 gen(1); std::normal_distribution< double > nd(0.0, 1.0);
   for (int k = 0; k < 8; ++k)
   {
      Vec3 n = c[size_t(k)].second; double v = c[size_t(k)].first;
      for (double step : { 0.03, 0.01, 0.003, 0.001, 0.0003 })
         for (int it = 0; it < 80; ++it)
         {
            Vec3 m = n + Vec3(real_t(step * nd(gen)), real_t(step * nd(gen)), real_t(step * nd(gen)));
            m = m / real_t(m.length());
            const double w = needed(m, a1, a2, d);
            if (w < v) { v = w; n = m; }
         }
      best = std::min(best, v);
   }
   return best;
}

mesa_pd::Rot3 rotFromQuat(double w, double x, double y, double z)
{
   const double n = std::sqrt(w * w + x * x + y * y + z * z);
   return mesa_pd::Rot3(mesa_pd::Quat(real_t(w / n), real_t(x / n), real_t(y / n), real_t(z / n)));
}
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

struct Pair { Vec3 pa, pb; mesa_pd::Rot3 ra, rb; };
struct Eval { bool contact; double overlap; double need; double ub; double exact; double gjk; const char* scenario; };

Eval evaluate(const Pair& p, bool withGjk)
{
   Eval e{};
   const Vec3 a1 = p.ra.getMatrix() * Vec3(0, 0, 1), a2 = p.rb.getMatrix() * Vec3(0, 0, 1), d = p.pb - p.pa;
   e.exact = exactPenetrationDepth(a1, a2, d);
   const auto res = closed_form_contacts::cylinderCylinderContact({ p.pa, p.ra.getMatrix(), real_t(kR), real_t(kL) }, { p.pb, p.rb.getMatrix(), real_t(kR), real_t(kL) });
   e.contact = res.contact; e.overlap = res.contact ? double(res.overlap) : 0.0; e.scenario = res.scenario;
   if (res.contact)
   {
      const Vec3 n = res.normal;
      e.need = std::min(needed(n, a1, a2, d), needed(-n, a1, a2, d));
      e.ub   = upperBound(a1, a2, d, n);
   }
   if (withGjk)
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
      auto c0 = ps->create(); c0->setPosition(p.pa); c0->setRotation(p.ra);
      auto c1 = ps->create(); c1->setPosition(p.pb); c1->setRotation(p.rb);
      const mesa_pd::data::Cylinder shape{ real_t(kR), real_t(kL) };
      SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = false;
      e.gjk = scd(0, 1, shape, shape, ac) ? -double(scd.getPenetrationDepth()) : 0.0;
   }
   return e;
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::antidunes;
   int failures = 0;
   if (std::getenv("CYL_CENTRE_CAP") != nullptr) { closed_form_contacts::kodamCapSelectionByAxisDistance = false; std::cout << "[using the OLD 3D centre-distance cap selection instead of Kodam Eq. 68-70]\n"; }
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };
   const double tol = 0.05 * kR;   // verifier tolerance: 5 % of the radius

   // ---------------- (A) real pairs from the cylinder-hopper floater log (job 4250907, steps 5650/5850): deeply nested rods ----------------
   struct Real { const char* name; double x1, y1, z1, q1[4], x2, y2, z2, q2[4]; };
   const Real real[] = {
      { "step5650 rods 176/180", 353.428, 10.0374, 33.9364, { 0.675649, -0.632585, 0.378261, 0.0159054 }, 352.427, 10.0585, 34.7087, { 0.662445, -0.664601, -0.116257, -0.325508 } },
      { "step5850 rods 176/180", 353.161, 10.3119, 34.0849, { 0.675306, -0.640524, 0.364316, -0.031058 }, 352.46, 10.1747, 34.4936, { 0.618026, -0.676647, -0.127839, -0.379275 } },
      { "step5650 rods 161/165", 372.913, 10.2225, 4.5596, { 0.69175, -0.582072, 0.357119, -0.234819 }, 371.679, 12.7248, 3.88215, { 0.698659, -0.60131, 0.387673, 0.00345134 } },
      { "step5850 rods 782/786", 365.368, 30.0509, 9.80322, { 0.590926, -0.350542, -0.657838, -0.308507 }, 364.46, 29.9638, 8.91965, { 0.708915, -0.661348, -0.23094, 0.081992 } },
      { "step5850 rods 571/572", 297.125, 27.244, 21.3028, { 0.703483, -0.666895, -0.233038, 0.0778158 }, 297.231, 27.4026, 21.0605, { 0.654417, -0.71362, -0.0481814, 0.245281 } } };
   std::cout << "(A) real deeply-nested hopper pairs:\n";
   for (const auto& r : real)
   {
      const Pair p{ Vec3(real_t(r.x1), real_t(r.y1), real_t(r.z1)), Vec3(real_t(r.x2), real_t(r.y2), real_t(r.z2)), rotFromQuat(r.q1[0], r.q1[1], r.q1[2], r.q1[3]), rotFromQuat(r.q2[0], r.q2[1], r.q2[2], r.q2[3]) };
      const Eval e = evaluate(p, true);
      std::cout << "   " << r.name << ": exact PD " << e.exact << " | closed form " << e.scenario << " overlap " << e.overlap << " (needs " << e.need << " along its own normal) | GJK/EPA " << e.gjk << std::endl;
      check(e.exact > 0.9 * 2 * kR - 0.6, std::string(r.name) + ": exact depth should be ~2R for a nested pair");
      check(e.contact && e.overlap < 0.25 * e.exact, std::string(r.name) + ": characterisation -- closed form is expected to UNDER-report this deep overlap by >4x");
      check(e.need - e.overlap > tol || e.overlap - e.ub > tol, std::string(r.name) + ": self-consistency check must flag the wrong closed-form result");
      check(std::abs(e.gjk - e.exact) < 0.25, std::string(r.name) + ": GJK/EPA must match the exact depth (EPA is slightly off for near-coincident bodies, <= 0.25 LU here)");
   }

   // ---------------- (B) random pairs, per depth class ----------------
   std::mt19937 gen(20260920); std::normal_distribution< double > nd(0.0, 1.0); std::uniform_real_distribution< double > ud(0.0, 1.0);
   struct Bin { double lo, hi; int n = 0, detected = 0, under = 0, underFlagged = 0, over = 0, overFlagged = 0, right = 0, rightFlagged = 0, gjkN = 0, gjkBad = 0; double maxUnder = 0, maxOver = 0, sumErr = 0, gjkMax = 0; };
   std::vector< Bin > bins = { { 0.0, 0.05 * kR }, { 0.05 * kR, 0.15 * kR }, { 0.15 * kR, 0.3 * kR }, { 0.3 * kR, 0.6 * kR }, { 0.6 * kR, 1.0 * kR }, { 1.0 * kR, 2.0 * kR }, { 2.0 * kR, 2.5 * kR } };
   int accepted = 0;
   for (int trial = 0; trial < 400000 && accepted < 14000; ++trial)
   {
      Vec3 a1(real_t(nd(gen)), real_t(nd(gen)), real_t(nd(gen))), a2(real_t(nd(gen)), real_t(nd(gen)), real_t(nd(gen)));
      a1 = a1 / real_t(a1.length()); a2 = a2 / real_t(a2.length());
      Vec3 dir(real_t(nd(gen)), real_t(nd(gen)), real_t(nd(gen))); dir = dir / real_t(dir.length());
      const double dist = 12.0 * std::pow(ud(gen), 1.5);   // biased towards small centre distances (deep overlaps) but covers the whole range
      const Pair p{ Vec3(0, 0, 0), dir * real_t(dist), rotToAxis(a1), rotToAxis(a2) };
      const Vec3 ea1 = p.ra.getMatrix() * Vec3(0, 0, 1), ea2 = p.rb.getMatrix() * Vec3(0, 0, 1);
      const double exact = exactPenetrationDepth(ea1, ea2, p.pb - p.pa);
      if (exact <= 0.0) continue;                              // separated
      Bin* b = nullptr; for (auto& x : bins) if (exact >= x.lo && exact < x.hi) b = &x;
      if (!b) continue;
      const bool wantGjk = (b->gjkN < 60);
      const Eval e = evaluate(p, wantGjk);
      ++accepted; ++b->n;
      if (wantGjk) { ++b->gjkN; b->gjkMax = std::max(b->gjkMax, std::abs(e.gjk - exact)); if (std::abs(e.gjk - exact) > 0.25) ++b->gjkBad; }
      if (e.contact) ++b->detected;
      // a missed contact (closed form finds nothing) is NOT a silent error in production -- the dispatch falls back to GJK/EPA then -- so it is counted as "flagged" under-estimate
      const double err = e.contact ? e.overlap - exact : -exact;
      b->sumErr += std::abs(err);
      const bool flaggedUnder = !e.contact || (e.need - e.overlap > tol);
      const bool flaggedOver  = e.contact && (e.overlap - e.ub > tol);
      if (err < -tol)      { ++b->under; b->maxUnder = std::max(b->maxUnder, -err); if (flaggedUnder) ++b->underFlagged; }
      else if (err > tol)  { ++b->over;  b->maxOver  = std::max(b->maxOver, err);   if (flaggedOver || flaggedUnder) ++b->overFlagged; }
      else                 { ++b->right; if (flaggedUnder || flaggedOver) ++b->rightFlagged; }
   }
   std::cout << "(B) random rod pairs (r=3, L=18, random axes/offsets) by TRUE depth (fraction of R); wrong = |closed-form overlap - exact| > 5% R:\n";
   std::cout << "   depth class      pairs  found   UNDER-est.  caught   OVER-est.  caught   right: false alarms   max under  max over   GJK max err\n";
   for (const auto& b : bins)
   {
      if (b.n == 0) continue;
      std::printf("   %4.2f-%4.2f R  %7d  %5.1f%%   %6.2f%%  %6.1f%%   %6.2f%%  %6.1f%%       %6.2f%%          %7.3f   %7.3f      %6.3f\n", b.lo / kR, b.hi / kR, b.n, 100.0 * b.detected / b.n,
                  100.0 * b.under / b.n, b.under ? 100.0 * b.underFlagged / b.under : 100.0, 100.0 * b.over / b.n, b.over ? 100.0 * b.overFlagged / b.over : 100.0,
                  b.right ? 100.0 * b.rightFlagged / b.right : 0.0, b.maxUnder, b.maxOver, b.gjkMax);
      check(b.gjkBad == 0, "GJK/EPA deviates from the exact depth by more than 0.25 LU in some pair");
      check(b.underFlagged == b.under, "the self-consistency check must catch EVERY under-estimate (theorem: needed(n) >= exact depth)");
      check(b.right == 0 || double(b.rightFlagged) / b.right <= 0.08, "false-alarm rate on correct closed-form results should be <= 8 %");
   }
   int overAll = 0, overCaught = 0;
   for (const auto& b : bins) { overAll += b.over; overCaught += b.overFlagged; }
   std::cout << "   over-estimates in total: " << overAll << ", caught by the upper-bound/self-consistency check: " << overCaught << std::endl;
   check(overAll == 0 || double(overCaught) / overAll >= 0.95, "the upper-bound check should catch >= 95 % of all over-estimates");
   // range claim: (nearly) all shallow overlaps are right, a substantial fraction of the deep ones is wrong
   int shN = 0, shWrong = 0, deepN = 0, deepWrong = 0;
   for (size_t k = 0; k < bins.size(); ++k) { const Bin& b = bins[k]; if (k <= 2) { shN += b.n; shWrong += b.under + b.over; } if (k == 4 || k == 5) { deepN += b.n; deepWrong += b.under + b.over; } }
   std::cout << "   shallow (<0.3 R): " << shWrong << "/" << shN << " wrong; deep (0.6-2 R): " << deepWrong << "/" << deepN << " wrong" << std::endl;
   check(shN > 500 && double(shWrong) / shN < 0.01, "closed form should be right for >= 99 % of the shallow (< 0.3 R) overlaps");
   check(deepN > 500 && double(deepWrong) / deepN > 0.10, "closed form should be wrong for > 10 % of the deep (0.6-2 R) overlaps");

   std::cout << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failures)" << std::endl;
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
