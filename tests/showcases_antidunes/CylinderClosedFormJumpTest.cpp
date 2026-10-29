//======================================================================================================================
//! \file CylinderClosedFormJumpTest.cpp
//! Why does the closed-form Kodam+Guo cylinder dispatch report an overlap of 0.0034 LU and one sub-cycle later 2.5 LU for the SAME geometry?
//! Data: rods ...486/...489 of job 4251049 (analytic hybrid; launched at step 3493), poses of steps 3484..3489 as logged (CylinderClosedFormJumpData.h). We
//!   (1) recompute the closed form, GJK/EPA (the reference) and the exact penetration depth for EVERY logged sub-cycle,
//!   (2) walk between two consecutive logged poses in 50 tiny steps and locate the discontinuity,
//!   (3) show that it comes from the band-edge branch's choice of WHICH rim circle of the probe cylinder is examined: rimVsFaceOrBand() only looks at the probe's NEARER end-cap circle
//!       (cap centre closest to the host's centre) and returns the closest point of THAT circle to the host's axis; when the two caps are about equally far, the choice flips and the reported
//!       radial overlap jumps between the two circles' values (one circle just grazes the band, the other one dips deep into it).
//! The exact penetration depth stays ~3.5 LU throughout; GJK/EPA reports it (within EPA accuracy).
//======================================================================================================================
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"
#include "CylinderClosedFormJumpData.h"

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

struct Pose { Vec3 p; double q[4]; };
Pose fromRow(const double* r) { return Pose{ Vec3(real_t(r[0]), real_t(r[1]), real_t(r[2])), { r[3], r[4], r[5], r[6] } }; }
mesa_pd::Rot3 rotOf(const Pose& p)
{
   const double n = std::sqrt(p.q[0] * p.q[0] + p.q[1] * p.q[1] + p.q[2] * p.q[2] + p.q[3] * p.q[3]);
   return mesa_pd::Rot3(mesa_pd::Quat(real_t(p.q[0] / n), real_t(p.q[1] / n), real_t(p.q[2] / n), real_t(p.q[3] / n)));
}
closed_form_contacts::CylinderState stateOf(const Pose& p) { return closed_form_contacts::CylinderState{ p.p, rotOf(p).getMatrix(), real_t(kR), real_t(kL) }; }

double support(const Vec3& n, const Vec3& a) { const double c = std::abs(double(n * a)); return kH * c + kR * std::sqrt(std::max(0.0, 1.0 - c * c)); }
double exactPD(const Vec3& a1, const Vec3& a2, const Vec3& d)
{
   auto f = [&](const Vec3& n) { return support(n, a1) + support(n, a2) - double(d * n); };
   double best = 1e300; Vec3 bn;
   const int N = 8000; const double ga = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
   for (int i = 0; i < N; ++i) { const double z = 1.0 - 2.0 * (i + 0.5) / N, r = std::sqrt(1.0 - z * z), ph = ga * i; const Vec3 n(real_t(r * std::cos(ph)), real_t(r * std::sin(ph)), real_t(z)); const double v = f(n); if (v < best) { best = v; bn = n; } }
   std::mt19937 gen(3); std::normal_distribution< double > nd(0.0, 1.0);
   for (double step : { 0.03, 0.01, 0.003, 0.001, 0.0003 }) for (int it = 0; it < 100; ++it) { Vec3 m = bn + Vec3(real_t(step * nd(gen)), real_t(step * nd(gen)), real_t(step * nd(gen))); m = m / real_t(m.length()); const double v = f(m); if (v < best) { best = v; bn = m; } }
   return best;
}
double gjkDepth(const Pose& a, const Pose& b)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   auto c0 = ps->create(); c0->setPosition(a.p); c0->setRotation(rotOf(a));
   auto c1 = ps->create(); c1->setPosition(b.p); c1->setRotation(rotOf(b));
   const mesa_pd::data::Cylinder shape{ real_t(kR), real_t(kL) };
   SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = false;
   return scd(0, 1, shape, shape, ac) ? -double(scd.getPenetrationDepth()) : 0.0;
}
//! Radial distance from the HOST's axis line of the closest point of one of the PROBE's two rim circles (brute-force over the angle): overlap the band-edge branch would report for that circle.
double rimOverlap(const closed_form_contacts::CylinderState& host, const closed_form_contacts::CylinderState& probe, int cap)
{
   const Vec3 e = probe.axis(), eh = host.axis();
   Vec3 u, v; closed_form_contacts::detail::perpendicularBasis(e, u, v);
   const Vec3 centre = probe.pos + real_t(cap == 0 ? 0.5 : -0.5) * probe.length * e;
   double best = 1e300;
   for (int k = 0; k < 36000; ++k)
   {
      const double th = 2.0 * 3.14159265358979323846 * k / 36000.0;
      const Vec3 P = centre + real_t(kR * std::cos(th)) * u + real_t(kR * std::sin(th)) * v;
      const Vec3 rel = P - host.pos; const double rho = (rel - eh * (rel * eh)).length();
      best = std::min(best, rho);
   }
   return kR - best;
}
//! which of the probe's caps rimVsFaceOrBand examines: the cap centre closer to the host's centre (0 = +axis end, 1 = -axis end)
int chosenCap(const closed_form_contacts::CylinderState& host, const closed_form_contacts::CylinderState& probe)
{
   const Vec3 e = probe.axis(); const Vec3 top = probe.pos + real_t(0.5) * probe.length * e, bot = probe.pos - real_t(0.5) * probe.length * e;
   return ((top - host.pos).length() < (bot - host.pos).length()) ? 0 : 1;
}
Pose lerp(const Pose& a, const Pose& b, double t)
{
   Pose r; r.p = a.p * real_t(1 - t) + b.p * real_t(t);
   double dot = 0; for (int i = 0; i < 4; ++i) dot += a.q[i] * b.q[i];
   const double s = dot < 0 ? -1.0 : 1.0; for (int i = 0; i < 4; ++i) r.q[i] = (1 - t) * a.q[i] + t * s * b.q[i];
   return r;
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::antidunes;
   int failures = 0;
   closed_form_contacts::kodamCapSelectionByAxisDistance = false;   // parts (1)-(3) reproduce the behaviour at the time of the crash (3D centre-distance cap choice, deviation from Kodam Eq. 68-70)
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };
   constexpr size_t N = sizeof(kLoggedRows) / sizeof(kLoggedRows[0]);

   std::cout << "(1) every logged sub-cycle of the pair (A=...486, B=...489): logged closed-form result | recomputed closed form (from the 6-digit CSV poses) | GJK/EPA | exact depth | chosen rim caps (host=A/probe=B, host=B/probe=A)\n";
   int reproduced = 0;
   for (size_t k = 0; k < N; ++k)
   {
      const auto& r = kLoggedRows[k];
      const Pose a = fromRow(r.a), b = fromRow(r.b);
      const auto ca = stateOf(a), cb = stateOf(b);
      const auto res = closed_form_contacts::cylinderCylinderContact(ca, cb);
      const double gjk = gjkDepth(a, b), ex = exactPD(ca.axis(), cb.axis(), cb.pos - ca.pos);
      const bool same = res.contact && std::string(res.scenario) == r.scenario && std::abs(double(res.overlap) - r.overlap) < 0.05 * std::max(1.0, r.overlap);
      reproduced += same;
      if (r.sub == 0 || (r.step == 3485 && r.sub >= 5 && r.sub <= 8))
         std::printf("   %d.%d | logged %-10s %.4f | recomputed %-10s %.4f%s | GJK %.3f | exact %.3f | caps %d,%d\n", r.step, r.sub, r.scenario, r.overlap, res.contact ? res.scenario : "none", res.contact ? double(res.overlap) : 0.0,
                     same ? " (reproduced)" : " (differs)", gjk, ex, chosenCap(ca, cb), chosenCap(cb, ca));
   }
   if (std::getenv("CYL_JUMP_TABLE") != nullptr)
   {
      std::cout << "TABLE step.sub | closed form (recomputed) scenario | overlap | GJK/EPA depth | exact depth | closed form - GJK | jump vs previous row\n";
      double prev = -1;
      for (size_t k = 0; k < N; ++k)
      {
         const auto& r = kLoggedRows[k]; const Pose a = fromRow(r.a), b = fromRow(r.b); const auto ca = stateOf(a), cb = stateOf(b);
         const auto res = closed_form_contacts::cylinderCylinderContact(ca, cb); const double ov = res.contact ? double(res.overlap) : 0.0;
         const double gjk = gjkDepth(a, b), ex = exactPD(ca.axis(), cb.axis(), cb.pos - ca.pos);
         std::printf("TABLE %d.%d | %s | %.4f | %.4f | %.4f | %+.4f | %s\n", r.step, r.sub, res.contact ? res.scenario : "none", ov, gjk, ex, ov - gjk, (prev >= 0 && std::abs(ov - prev) > 0.5) ? "JUMP" : "");
         prev = ov;
      }
   }
   std::cout << "   recomputed == logged for " << reproduced << " of " << N << " sub-cycles (positions have only 6 significant digits in the CSV)\n";

   // (2)+(3): walk between logged poses 3485.6 and 3485.7 (the jump 0.0034 -> 2.53 in the log)
   const LoggedPoseRow *r0 = nullptr, *r1 = nullptr;
   for (size_t k = 0; k + 1 < N; ++k) if (kLoggedRows[k].step == 3485 && kLoggedRows[k].sub == 6) { r0 = &kLoggedRows[k]; r1 = &kLoggedRows[k + 1]; }
   check(r0 && r1, "logged rows 3485.6/3485.7 present");
   if (!(r0 && r1)) return EXIT_FAILURE;
   check(r0->overlap < 0.01 && r1->overlap > 2.0, "the log itself shows the jump: 3485.6 overlap < 0.01, 3485.7 overlap > 2");
   const Pose a0 = fromRow(r0->a), b0 = fromRow(r0->b), a1 = fromRow(r1->a), b1 = fromRow(r1->b);
   std::cout << "(2) walk between the logged poses 3485.6 -> 3485.7 (position change " << (a1.p - a0.p).length() + (b1.p - b0.p).length() << " LU in total):\n"
             << "    t     closed form (scenario overlap)   GJK    exact | caps host=A,B | rim-circle overlaps for host=A (cap0, cap1) | for host=B (cap0, cap1)\n";
   double prevOv = -1, jumpAt = -1, jumpSize = 0; int prevCap = -1, capFlipAt = -1; double maxGjkChange = 0, maxExactChange = 0, prevGjk = 0, prevEx = 0;
   for (int i = 0; i <= 50; ++i)
   {
      const double t = i / 50.0; const Pose a = lerp(a0, a1, t), b = lerp(b0, b1, t);
      const auto ca = stateOf(a), cb = stateOf(b);
      const auto res = closed_form_contacts::cylinderCylinderContact(ca, cb);
      const double ov = res.contact ? double(res.overlap) : 0.0;
      const int capA = chosenCap(ca, cb), capB = chosenCap(cb, ca);
      const double gjk = gjkDepth(a, b), ex = exactPD(ca.axis(), cb.axis(), cb.pos - ca.pos);
      if (i % 5 == 0 || (prevOv >= 0 && std::abs(ov - prevOv) > 0.5))
         std::printf("    %.2f  %-10s %8.4f            %6.3f %6.3f |   %d,%d       |   %7.4f, %7.4f                    |   %7.4f, %7.4f\n", t, res.contact ? res.scenario : "none", ov, gjk, ex, capA, capB,
                     rimOverlap(ca, cb, 0), rimOverlap(ca, cb, 1), rimOverlap(cb, ca, 0), rimOverlap(cb, ca, 1));
      if (prevOv >= 0) { if (std::abs(ov - prevOv) > jumpSize) { jumpSize = std::abs(ov - prevOv); jumpAt = t; } maxGjkChange = std::max(maxGjkChange, std::abs(gjk - prevGjk)); maxExactChange = std::max(maxExactChange, std::abs(ex - prevEx)); }
      if (prevCap >= 0 && capA != prevCap) capFlipAt = i;
      prevOv = ov; prevCap = capA; prevGjk = gjk; prevEx = ex;
   }
   std::cout << "    largest closed-form step between neighbouring t: " << jumpSize << " LU at t~" << jumpAt << "; largest exact-depth change per step " << maxExactChange << " LU, GJK " << maxGjkChange << " LU\n";
   check(jumpSize > 1.0, "the closed form must JUMP by > 1 LU somewhere in this tiny pose interval");
   check(maxExactChange < 0.05, "the exact depth changes smoothly (< 0.05 LU per step) over the same interval");
   check(capFlipAt >= 0, "the probe rim circle examined by the band-edge branch (host=A) flips within the interval");
   check(capFlipAt >= 0 && std::abs(jumpAt - capFlipAt / 50.0) < 0.041, "the jump coincides with the flip of the examined rim circle");

   // ---------------- (4) cross-check against the paper: Kodam Eq. 68-70 choose the cap closer to the host's AXIS (in the host's xy-plane), this implementation chooses the cap closer to the host's CENTRE (3D) ----------------
   if (std::getenv("CYL_JUMP_TABLE") != nullptr)
   {
      std::cout << "TABLE2 step.sub | 3D-centre rule: scenario overlap | Kodam axis rule: scenario overlap | GJK\n";
      for (size_t k = 0; k < N; ++k)
      {
         const auto& r = kLoggedRows[k]; const Pose a = fromRow(r.a), b = fromRow(r.b); const auto ca = stateOf(a), cb = stateOf(b);
         closed_form_contacts::kodamCapSelectionByAxisDistance = false; const auto r3 = closed_form_contacts::cylinderCylinderContact(ca, cb);
         closed_form_contacts::kodamCapSelectionByAxisDistance = true;  const auto rk = closed_form_contacts::cylinderCylinderContact(ca, cb);
         closed_form_contacts::kodamCapSelectionByAxisDistance = false;
         std::printf("TABLE2 %d.%d | %s %.4f | %s %.4f | %.4f\n", r.step, r.sub, r3.contact ? r3.scenario : "none", r3.contact ? double(r3.overlap) : 0.0, rk.contact ? rk.scenario : "none", rk.contact ? double(rk.overlap) : 0.0, gjkDepth(a, b));
      }
   }
   if (std::getenv("CYL_JUMP_TABLE") != nullptr)
   {
      std::cout << "TABLE3 step.sub | rule | scenario | reported overlap delta | distance needed along ITS OWN normal s(n) | gap s(n)-delta | exact depth | GJK\n";
      for (size_t k = 0; k < N; ++k)
      {
         const auto& r = kLoggedRows[k]; if (!(r.sub == 0 || (r.step == 3485 && r.sub == 7))) continue;
         const Pose a = fromRow(r.a), b = fromRow(r.b); const auto ca = stateOf(a), cb = stateOf(b);
         for (int rule = 0; rule < 2; ++rule)
         {
            closed_form_contacts::kodamCapSelectionByAxisDistance = (rule == 1);
            const auto res = closed_form_contacts::cylinderCylinderContact(ca, cb);
            closed_form_contacts::kodamCapSelectionByAxisDistance = false;
            const Vec3 d = cb.pos - ca.pos; const Vec3 n = res.normal;
            auto need = [&](const Vec3& m) { return support(m, ca.axis()) + support(m, cb.axis()) - double(d * m); };
            const double s_n = std::min(need(n), need(-n));
            std::printf("TABLE3 %d.%d | %s | %s | %.4f | %.4f | %.4f | %.4f | %.4f\n", r.step, r.sub, rule ? "Kodam axis rule" : "centre rule    ", res.contact ? res.scenario : "none", double(res.overlap), s_n, s_n - double(res.overlap),
                        exactPD(ca.axis(), cb.axis(), d), gjkDepth(a, b));
         }
      }
   }
   // ---------------- (5) with the paper's cap rule (now the default) the jump between the two logged poses is gone, and the early rows report the deep circle ----------------
   {
      closed_form_contacts::kodamCapSelectionByAxisDistance = true;
      double maxStep = 0, prevOv = -1;
      for (int i = 0; i <= 50; ++i)
      {
         const double t = i / 50.0; const auto ca = stateOf(lerp(a0, a1, t)), cb = stateOf(lerp(b0, b1, t));
         const auto res = closed_form_contacts::cylinderCylinderContact(ca, cb); const double ov = res.contact ? double(res.overlap) : 0.0;
         if (prevOv >= 0) maxStep = std::max(maxStep, std::abs(ov - prevOv));
         prevOv = ov;
      }
      std::cout << "(5) paper's cap rule (Eq. 68-70): largest closed-form change between neighbouring t over the same interval: " << maxStep << " LU\n";
      check(maxStep < 0.1, "with Kodam's cap selection the closed form is continuous over the interval that jumped by 2.5 LU with the centre-distance rule");
      int deepRows = 0;
      for (size_t k = 0; k < N; ++k)
      {
         const auto& r = kLoggedRows[k]; if (r.step == 3485 && r.sub > 6) break;
         const auto res = closed_form_contacts::cylinderCylinderContact(stateOf(fromRow(r.a)), stateOf(fromRow(r.b)));
         deepRows += (res.contact && double(res.overlap) > 2.0);
      }
      check(deepRows == 17, "with Kodam's cap selection all 17 early rows (3484.0-3485.6) report the deep rim circle (overlap > 2 LU) instead of ~0.003");
      closed_form_contacts::kodamCapSelectionByAxisDistance = false;
   }
   if (std::getenv("CYL_JUMP_TABLE") != nullptr)
   {
      std::cout << "TABLE4 step.sub | normal of the closed form . normal of GJK/EPA (+1 same direction, -1 exactly opposite: the closed-form force would push the rods DEEPER) | closed form scenario, overlap\n";
      closed_form_contacts::kodamCapSelectionByAxisDistance = true;
      for (size_t k = 0; k < N; ++k)
      {
         const auto& r = kLoggedRows[k]; if (r.sub % 3 != 0 && !(r.step == 3485 && r.sub == 7)) continue;
         const Pose a = fromRow(r.a), b = fromRow(r.b); const auto ca = stateOf(a), cb = stateOf(b);
         const auto res = closed_form_contacts::cylinderCylinderContact(ca, cb);
         auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
         auto c0 = ps->create(); c0->setPosition(a.p); c0->setRotation(rotOf(a)); auto c1 = ps->create(); c1->setPosition(b.p); c1->setRotation(rotOf(b));
         const mesa_pd::data::Cylinder shape{ real_t(kR), real_t(kL) };
         SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = false; scd(0, 1, shape, shape, ac);
         std::printf("TABLE4 %d.%d | %+.3f | %s %.3f\n", r.step, r.sub, double(res.normal * scd.getContactNormal()), res.contact ? res.scenario : "none", double(res.overlap));
      }
      closed_form_contacts::kodamCapSelectionByAxisDistance = false;
   }
   std::cout << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failures)" << std::endl;
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
