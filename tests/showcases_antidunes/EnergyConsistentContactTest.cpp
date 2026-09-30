//======================================================================================================================
//! \file EnergyConsistentContactTest.cpp
//! Feng-style ENERGY-CONSISTENT (ECC) normal contact for two cylinders vs the contact models we use today -- a reference implementation that lives ONLY in this test (nothing here is used by the simulation).
//!
//! Feng (2023, Acta Mech. Sin. 39, 722343), Sect. 6: the normal contact force must derive from a scalar contact-energy potential w(x1,theta1,x2,theta2) (F = -dw/dx, moments = -dw/dtheta); with the
//! MINKOWSKI OVERLAP delta_M (minimal translation that separates the bodies, ref. [132]) as the overlap and w = 1/2 k delta_M^2 the model is energy-conserving for any convex shape.
//! For two cylinders (axis a_i, half length H, radius R; support function h_i(n) = H|n.a_i| + R sqrt(1-(n.a_i)^2)):
//!     delta_M = min_n s(n),  s(n) = h_1(n) + h_2(n) - (x2-x1).n,   n* = argmin   (body 2 is pushed along +n*, body 1 along -n*)
//!     F_2 = +k delta n*,  F_1 = -k delta n*,   torque_i = -k delta g_i (a_i x n*),   g_i = H sgn(c_i) - R c_i / sqrt(1-c_i^2),  c_i = n*.a_i     (envelope theorem)
//! i.e. the force acts along the Minkowski normal at the two SUPPORT (witness) points; because the witness points differ by a vector parallel to n* the pair can equivalently act at their midpoint.
//! NOTE: for rods the generic optimum lies on a KINK of the support function (n perpendicular to an axis = line contact, n along an axis = cap contact); there the potential is not differentiable in the orientation
//! and the torque is a sub-gradient. We pick the element that keeps the two witness points on a common line parallel to n* and puts the resultant at the centre of the contact patch (closest points of the two
//! support sets in the plane perpendicular to n*; centre of the overlap for parallel lines/discs) -- the 'centroid of the contact area' choice of Feng's contact-area/volume models.
//! We (1) verify the analytic force/torque against finite-difference gradients of w (proof that the model is conservative), (2) compare normal / overlap / lever arm with what the simulation
//! uses today -- GJK/EPA (useAnalyticalContactDetection=false) and the closed-form Kodam+Guo hybrid (=true) -- on random pairs, (3) run elastic rod-rod impacts and measure the energy that each model
//! creates or destroys (an undamped spring must return exactly the kinetic energy it received), (4) measure normal jumps along a smooth relative motion.
//======================================================================================================================
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace walberla;
using walberla::antidunes::SolidCylinderContactDetection;
namespace {
using V = std::array< double, 3 >;
V operator+(const V& a, const V& b) { return { a[0] + b[0], a[1] + b[1], a[2] + b[2] }; }
V operator-(const V& a, const V& b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
V operator*(const V& a, double s) { return { a[0] * s, a[1] * s, a[2] * s }; }
V operator*(double s, const V& a) { return a * s; }
V& operator+=(V& a, const V& b) { a = a + b; return a; }
double dot(const V& a, const V& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V cross(const V& a, const V& b) { return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] }; }
double norm(const V& a) { return std::sqrt(dot(a, a)); }
V unit(const V& a) { return a * (1.0 / norm(a)); }
V rodrigues(const V& v, const V& w)   // rotate v about w by angle |w|
{
   const double th = norm(w); if (th < 1e-14) return v;
   const V k = w * (1.0 / th); return v * std::cos(th) + cross(k, v) * std::sin(th) + k * (dot(k, v) * (1.0 - std::cos(th)));
}

constexpr double kR = 3.0, kL = 18.0, kH = 9.0;
constexpr double kMass = 3.14159265358979323846 * kR * kR * kL;                    // density 1
constexpr double kIax = 0.5 * kMass * kR * kR, kIperp = kMass * (3.0 * kR * kR + kL * kL) / 12.0;

double support(const V& n, const V& a) { const double c = std::abs(dot(n, a)); return kH * c + kR * std::sqrt(std::max(0.0, 1.0 - c * c)); }
double sOf(const V& n, const V& a1, const V& a2, const V& d) { return support(n, a1) + support(n, a2) - dot(d, n); }

struct Mink { double depth = 0; V n{ 0, 0, 1 }; };
//! Minkowski overlap: global Fibonacci search (if no warm start) followed by pattern-search refinement
Mink minkowski(const V& a1, const V& a2, const V& d, const V* warm = nullptr)
{
   static std::mt19937 gen(12345); std::normal_distribution< double > nd(0.0, 1.0);
   const int N = 3000; const double ga = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
   std::vector< std::pair< double, V > > c; c.reserve(N + 1);
   for (int i = 0; i < N; ++i) { const double z = 1.0 - 2.0 * (i + 0.5) / N, rr = std::sqrt(1.0 - z * z), ph = ga * i; const V n{ rr * std::cos(ph), rr * std::sin(ph), z }; c.emplace_back(sOf(n, a1, a2, d), n); }
   std::partial_sort(c.begin(), c.begin() + 4, c.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
   std::vector< std::pair< double, V > > starts(c.begin(), c.begin() + 4);
   if (warm) starts.emplace_back(sOf(*warm, a1, a2, d), *warm);
   double bestAll = 1e300; V bestN = starts[0].second;
   for (auto& st : starts)
   {
      V n = st.second; double v = st.first;
      for (double step : { 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001, 3e-5, 1e-5, 3e-6, 1e-6 })
         for (int it = 0; it < 40; ++it) { V m = unit(n + V{ step * nd(gen), step * nd(gen), step * nd(gen) }); const double w = sOf(m, a1, a2, d); if (w < v) { v = w; n = m; } }
      if (v < bestAll) { bestAll = v; bestN = n; }
   }
   Mink r; r.depth = bestAll; r.n = bestN; return r;
}

struct Force { bool hit = false; V F1{ 0, 0, 0 }, F2{ 0, 0, 0 }, T1{ 0, 0, 0 }, T2{ 0, 0, 0 }; double depth = 0; V n{ 0, 0, 1 }; V p{ 0, 0, 0 }; double gap = 0; };   // n: direction body 2 is pushed along (from 1 to 2)

double gOf(const V& n, const V& a) { const double c = dot(n, a), s = std::sqrt(std::max(1e-12, 1.0 - c * c)); return kH * (c >= 0 ? 1.0 : -1.0) - kR * c / s; }

enum class SupType { Point, Line, Cap };
SupType supType(const V& n, const V& a) { const double c = std::abs(dot(n, a)); return c < 1e-3 ? SupType::Line : (c > 1.0 - 1e-3 ? SupType::Cap : SupType::Point); }
V pointSupport(const V& x, const V& a, const V& n) { const double c = dot(n, a), s = std::sqrt(std::max(1e-12, 1.0 - c * c)); return x + a * (kH * (c >= 0 ? 1.0 : -1.0)) + (n - a * c) * (kR / s); }
V perpUnit(const V& n, const V& a) { return unit(n - a * dot(n, a)); }
V capCentre(const V& x, const V& a, const V& n) { return x + a * (kH * (dot(n, a) >= 0 ? 1.0 : -1.0)); }
V clampDisc(const V& q, const V& c, const V& a) { V v = q - c; v = v - a * dot(v, a); const double l = norm(v); return l > kR ? c + v * (kR / l) : c + v; }
//! closest points of two segments (p1 + t d1, |t| <= h1) and (p2 + u d2, |u| <= h2); for (near-)parallel segments: centre of the overlap
void segSeg(const V& p1, const V& d1, const V& p2, const V& d2, V& q1, V& q2)
{
   const V r = p1 - p2; const double a = 1.0, b = dot(d1, d2), c = 1.0, d = dot(d1, r), e = dot(d2, r), den = a * c - b * b;
   double t, u;
   if (den < 1e-9)
   {   // parallel: axial overlap interval of segment 2 in the coordinate of segment 1
      const double s0 = dot(p2 - p1, d1), sgn = b >= 0 ? 1.0 : -1.0; (void)sgn;
      const double lo = std::max(-kH, s0 - kH), hi = std::min(kH, s0 + kH), mid = (lo <= hi) ? 0.5 * (lo + hi) : std::clamp(s0, -kH, kH);
      t = mid; u = dot(p1 + d1 * t - p2, d2);
   }
   else
   {
      t = std::clamp((b * e - c * d) / den, -kH, kH); u = std::clamp((a * e - b * d) / den, -kH, kH);
      t = std::clamp((b * u - d) / a, -kH, kH); u = std::clamp((b * t + e) / c, -kH, kH);
      for (int it = 0; it < 4; ++it) { t = std::clamp((b * u - d) / a, -kH, kH); u = std::clamp((b * t + e) / c, -kH, kH); }
   }
   q1 = p1 + d1 * t; q2 = p2 + d2 * u;
}
//! witness (contact-patch) points of the two bodies for the Minkowski normal n (body 1 supports along +n, body 2 along -n)
void witnessPoints(const V& x1, const V& a1, const V& x2, const V& a2, const V& n, V& p1, V& p2)
{
   const SupType t1 = supType(n, a1), t2 = supType(n * -1.0, a2);
   const V mn = n * -1.0;
   auto lineBase1 = [&]() { return x1 + perpUnit(n, a1) * kR; }; auto lineBase2 = [&]() { return x2 + perpUnit(mn, a2) * kR; };
   if (t1 == SupType::Point && t2 == SupType::Point) { p1 = pointSupport(x1, a1, n); p2 = pointSupport(x2, a2, mn); }
   else if (t1 == SupType::Line && t2 == SupType::Point) { p2 = pointSupport(x2, a2, mn); const double t = std::clamp(dot(p2 - x1, a1), -kH, kH); p1 = lineBase1() + a1 * t; }
   else if (t1 == SupType::Point && t2 == SupType::Line) { p1 = pointSupport(x1, a1, n); const double u = std::clamp(dot(p1 - x2, a2), -kH, kH); p2 = lineBase2() + a2 * u; }
   else if (t1 == SupType::Line && t2 == SupType::Line) segSeg(lineBase1() - a1 * 0.0, a1, lineBase2(), a2, p1, p2);
   else if (t1 == SupType::Cap && t2 == SupType::Point) { p2 = pointSupport(x2, a2, mn); const V c = capCentre(x1, a1, n); p1 = clampDisc(p2 - n * dot(p2 - c, n), c, a1); }
   else if (t1 == SupType::Point && t2 == SupType::Cap) { p1 = pointSupport(x1, a1, n); const V c = capCentre(x2, a2, mn); p2 = clampDisc(p1 - n * dot(p1 - c, n), c, a2); }
   else if (t1 == SupType::Cap && t2 == SupType::Cap) { const V c1 = capCentre(x1, a1, n), c2 = capCentre(x2, a2, mn); const V m = (c1 + c2) * 0.5; p1 = m - n * dot(m - c1, n); p2 = m - n * dot(m - c2, n); }
   else
   {   // cap vs line: scan the line, project along n onto the disc, take the closest projected pair (centroid of near-ties)
      const bool capFirst = (t1 == SupType::Cap); const V xc = capFirst ? x1 : x2, ac = capFirst ? a1 : a2, nc = capFirst ? n : mn; const V xl = capFirst ? x2 : x1, al = capFirst ? a2 : a1;
      const V c = capCentre(xc, ac, nc), lbase = xl + perpUnit(capFirst ? mn : n, al) * kR;
      double dmin = 1e300; std::vector< std::pair< double, double > > cand;
      for (int i = 0; i <= 90; ++i) { const double t = -kH + 2.0 * kH * i / 90.0; const V q = lbase + al * t; const V pc = clampDisc(q - nc * dot(q - c, nc), c, ac); const double dd = norm((q - pc) - nc * dot(q - pc, nc)); cand.emplace_back(dd, t); dmin = std::min(dmin, dd); }
      double ts = 0; int cnt = 0; for (const auto& cd : cand) if (cd.first <= dmin + 1e-3) { ts += cd.second; ++cnt; }
      const V q = lbase + al * (ts / std::max(1, cnt)); const V pc = clampDisc(q - nc * dot(q - c, nc), c, ac);
      if (capFirst) { p1 = pc; p2 = q; } else { p1 = q; p2 = pc; }
   }
}

Force eccForce(const V& x1, const V& a1, const V& x2, const V& a2, double k, V* warm)
{
   Force f; const V d = x2 - x1;
   if (norm(d) > 2.0 * std::hypot(kR, kH) + 1e-9) return f;
   const Mink m = minkowski(a1, a2, d, warm);
   if (warm) *warm = m.n;
   if (m.depth <= 0) return f;
   f.hit = true; f.depth = m.depth; f.n = m.n;
   f.F2 = m.n * (k * m.depth); f.F1 = f.F2 * -1.0;
   V p1, p2; witnessPoints(x1, a1, x2, a2, m.n, p1, p2);
   f.T1 = cross(p1 - x1, f.F1); f.T2 = cross(p2 - x2, f.F2); f.p = (p1 + p2) * 0.5; f.gap = norm(cross(p1 - p2, m.n));   // gap: witness points must lie on a common line parallel to n
   return f;
}

mesa_pd::Rot3 rotToAxis(const V& axis)
{
   mesa_pd::Rot3 rot; const mesa_pd::Vec3 ez{ real_t(0), real_t(0), real_t(1) }; const mesa_pd::Vec3 ax{ real_t(axis[0]), real_t(axis[1]), real_t(axis[2]) };
   const real_t cosA = ez * ax;
   if (cosA < real_t(1) - real_t(1e-14))
   {
      const mesa_pd::Vec3 rotAxis = (cosA > real_t(-1) + real_t(1e-14)) ? (ez % ax).getNormalized() : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
      rot.rotate(rotAxis, std::acos(std::clamp(cosA, real_t(-1), real_t(1))));
   }
   return rot;
}

//! contact of the simulation's own detection (GJK/EPA or the closed-form hybrid) + a linear spring F = k*depth along its normal at its contact point, applied to both bodies at that single point
Force simForce(const V& x1, const V& a1, const V& x2, const V& a2, double k, bool analytic)
{
   Force f; if (norm(x2 - x1) > 2.0 * std::hypot(kR, kH) + 1e-9) return f;
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   auto c0 = ps->create(); c0->setPosition(mesa_pd::Vec3(real_t(x1[0]), real_t(x1[1]), real_t(x1[2]))); c0->setRotation(rotToAxis(a1));
   auto c1 = ps->create(); c1->setPosition(mesa_pd::Vec3(real_t(x2[0]), real_t(x2[1]), real_t(x2[2]))); c1->setRotation(rotToAxis(a2));
   const mesa_pd::data::Cylinder shape{ real_t(kR), real_t(kL) };
   SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = analytic;
   if (!scd(0, 1, shape, shape, ac)) return f;
   const double depth = -double(scd.getPenetrationDepth()); if (depth <= 0) return f;
   const mesa_pd::Vec3 nn = scd.getContactNormal(), pp = scd.getContactPoint();
   const V n{ double(nn[0]), double(nn[1]), double(nn[2]) }, p{ double(pp[0]), double(pp[1]), double(pp[2]) };   // production convention: n points from body 2 toward body 1
   f.hit = true; f.depth = depth; f.n = n * -1.0;                                                                 // store as "from 1 to 2" like the ECC model
   f.F1 = n * (k * depth); f.F2 = f.F1 * -1.0; f.T1 = cross(p - x1, f.F1); f.T2 = cross(p - x2, f.F2); f.p = p;
   return f;
}

struct Body { V x, v, L, a; };
V omegaOf(const Body& b) { const double Ln = dot(b.L, b.a); const V Lax = b.a * Ln; const V Lperp = b.L - Lax; return Lax * (1.0 / kIax) + Lperp * (1.0 / kIperp); }
double energy(const Body& b) { return 0.5 * kMass * dot(b.v, b.v) + 0.5 * dot(b.L, omegaOf(b)); }
}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };
   std::mt19937 gen(2026); std::normal_distribution< double > nd(0.0, 1.0); std::uniform_real_distribution< double > ud(0.0, 1.0);
   const double meff = 0.5 * kMass, tc = 32.0, k = meff * 3.14159265358979323846 * 3.14159265358979323846 / (tc * tc);   // undamped spring with contact time tc (dt_DEM = 0.1)
   auto randAxis = [&]() { return unit(V{ nd(gen), nd(gen), nd(gen) }); };

   // ---------------------------------------------------------------- (1) analytic ECC force/torque == gradient of w = 1/2 k delta_M^2 (finite differences)
   //   generic optima sit on a KINK of the support function (line/cap contact), so the smooth formula is only checked on rim-rim (tip-to-tip) configurations, built explicitly
   {
      int n = 0, ok = 0; double worstF = 0, worstT = 0;
      for (int t = 0; t < 400 && n < 60; ++t)
      {
         const V n0 = randAxis(); V e1 = unit(cross(n0, V{ nd(gen), nd(gen), nd(gen) })), e2 = unit(cross(n0, V{ nd(gen), nd(gen), nd(gen) }));
         const V a1 = unit(n0 + e1 * (0.25 + 0.5 * ud(gen))), a2 = unit(n0 + e2 * (0.25 + 0.5 * ud(gen)));
         const double depth0 = 0.1 + 0.6 * ud(gen); const V d = n0 * (support(n0, a1) + support(n0, a2) - depth0);
         const Mink m = minkowski(a1, a2, d); if (m.depth < 0.02 || m.depth > 2.0) continue;
         const double c1 = std::abs(dot(m.n, a1)), c2 = std::abs(dot(m.n, a2)); if (c1 < 2e-2 || c1 > 1 - 2e-2 || c2 < 2e-2 || c2 > 1 - 2e-2) continue;
         V warm = m.n; const Force f = eccForce({ 0, 0, 0 }, a1, d, a2, k, &warm);
         auto w = [&](const V& x1, const V& b1, const V& x2, const V& b2) { V ww = m.n; const Mink q = minkowski(b1, b2, x2 - x1, &ww); return 0.5 * k * std::max(0.0, q.depth) * std::max(0.0, q.depth); };
         const double h = 2e-4; V Fnum{ 0, 0, 0 }, Tnum1{ 0, 0, 0 };
         for (int i = 0; i < 3; ++i)
         {
            V e{ 0, 0, 0 }; e[size_t(i)] = 1.0;
            Fnum[size_t(i)] = -(w({ 0, 0, 0 }, a1, d + e * h, a2) - w({ 0, 0, 0 }, a1, d - e * h, a2)) / (2 * h);                                   // -dw/dx2
            Tnum1[size_t(i)] = -(w({ 0, 0, 0 }, unit(rodrigues(a1, e * h)), d, a2) - w({ 0, 0, 0 }, unit(rodrigues(a1, e * -h)), d, a2)) / (2 * h);    // -dw/dtheta1
         }
         const V Tan1 = cross(a1, m.n) * (-k * m.depth * gOf(m.n, a1));                                                                              // analytic torque, closed form of the derivation above
         const double eF = norm(Fnum - f.F2) / std::max(1e-9, norm(f.F2)), eT = norm(Tnum1 - Tan1) / std::max(1e-9, norm(Tan1)), eW = norm(f.T1 - Tan1) / std::max(1e-9, norm(Tan1));
         ++n; worstF = std::max(worstF, eF); worstT = std::max(worstT, std::max(eT, eW)); ok += (eF < 2e-2 && eT < 5e-2 && eW < 1e-3);
      }
      std::printf("(1) analytic ECC force/torque vs finite-difference gradient of w = k delta_M^2 / 2 on smooth (tip-to-tip rim-rim) configurations: %d/%d pairs agree (force err <2%%, torque err <5%%, witness-point torque == closed form <0.1%%); worst force err %.2e, worst torque err %.2e\n", ok, n, worstF, worstT);
      check(n >= 20 && ok >= int(0.9 * n), "analytic ECC force/torque must equal the gradient of the potential (smooth rim-rim cases)");
   }

   // ---------------------------------------------------------------- (2) what the simulation reports today vs ECC, same random shallow/deep pairs
   {
      struct Acc { int n = 0; double angMean = 0, angMax = 0, depthErrMean = 0, depthErrMax = 0, torqueErrMean = 0, pointErr = 0, gapMax = 0; int flipped = 0, missed = 0; };
      const char* names[2] = { "GJK/EPA (analytic off)", "closed-form hybrid (analytic on)" };
      for (int model = 0; model < 2; ++model)
      {
         std::vector< Acc > acc(3); const char* cls[3] = { "shallow (depth < 0.3 R)", "mid (0.3-1 R)", "deep (1-2 R)" };
         for (int t = 0; t < 1200; ++t)
         {
            const V a1 = randAxis(), a2 = randAxis(); const V d = unit(V{ nd(gen), nd(gen), nd(gen) }) * (10.0 * std::pow(ud(gen), 1.6));
            const Mink m = minkowski(a1, a2, d); if (m.depth <= 0.02) continue;
            const int ci = m.depth < 0.3 * kR ? 0 : (m.depth < kR ? 1 : (m.depth < 2 * kR ? 2 : -1)); if (ci < 0) continue;
            V warm = m.n; const Force e = eccForce({ 0, 0, 0 }, a1, d, a2, k, &warm); const Force s = simForce({ 0, 0, 0 }, a1, d, a2, k, model == 1);
            Acc& A = acc[size_t(ci)]; ++A.n;
            if (!s.hit) { ++A.missed; continue; }
            const double cs = std::clamp(dot(s.n, e.n), -1.0, 1.0), ang = std::acos(cs) * 180.0 / 3.14159265358979323846;
            A.angMean += ang; A.angMax = std::max(A.angMax, ang); A.flipped += cs < 0;
            const double de = std::abs(s.depth - e.depth) / e.depth; A.depthErrMean += de; A.depthErrMax = std::max(A.depthErrMax, de);
            A.torqueErrMean += (norm(s.T1 - e.T1) + norm(s.T2 - e.T2)) / std::max(1e-9, 2.0 * k * e.depth * kH); A.pointErr += norm(s.p - e.p); A.gapMax = std::max(A.gapMax, e.gap);
         }
         std::cout << "(2) " << names[model] << " vs ECC (Minkowski normal/overlap/lever arm), 1200 random pairs:\n     class                       pairs | missed | normal angle mean/max [deg], flipped (>90) | rel. overlap error mean/max | torque error per body (units k*delta*H) | contact point distance [LU] | ECC witness-line residual max [LU]\n";
         for (int c = 0; c < 3; ++c)
         {
            const Acc& A = acc[size_t(c)]; const int nn = std::max(1, A.n - A.missed);
            std::printf("     %-26s %5d | %5d | %8.3f / %8.2f, %5d | %8.4f / %8.3f | %8.4f | %8.3f | %8.2e\n", cls[c], A.n, A.missed, A.angMean / nn, A.angMax, A.flipped, A.depthErrMean / nn, A.depthErrMax, A.torqueErrMean / nn, A.pointErr / nn, A.gapMax);
            if (model == 0) check(A.n > 30 && A.angMean / nn < 1.0 && A.depthErrMean / nn < 0.02 && A.torqueErrMean / nn < 0.02, std::string("GJK/EPA must agree with the Minkowski (ECC) normal, overlap and torque in class ") + cls[c]);
            if (model == 1 && c == 0) check(A.n > 30 && A.angMean / nn < 1.0 && A.depthErrMean / nn < 0.02, "the closed-form hybrid must agree with ECC for shallow contacts");
            if (model == 1 && c == 2) check(A.angMean / nn > 5.0 && A.depthErrMean / nn > 0.05 && A.flipped > 0, "the closed-form hybrid must visibly deviate from ECC (normal, overlap, flipped normals) for deep contacts");
         }
      }
   }

   // ---------------------------------------------------------------- (3) elastic rod-rod impacts: energy created/destroyed by each model (undamped spring, must be 0)
   {
      const double dt = 0.1;
      struct Res { double mean = 0, worst = 0; int n = 0; double meanAbs = 0; int contacts = 0; double maxDepth = 0; };
      const char* names[3] = { "ECC (Minkowski, energy-conserving)", "GJK/EPA + spring at its contact point", "closed-form hybrid + spring at its contact point" };
      std::cout << "(3) elastic impacts of two rods (undamped linear spring k = m_eff pi^2/tc^2, tc = 32, dt = 0.1, one rod at rest): relative change of the total (translational + rotational) kinetic energy after separation\n";
      for (double vrel : { 0.02, 0.06, 0.15 })
      {
         Res res[3];
         for (int t = 0; t < 16; ++t)
         {
            const V a1 = randAxis(), a2 = randAxis(); const V u = unit(V{ nd(gen), nd(gen), nd(gen) });
            // touching distance along u, then start 0.8 LU further out
            double lo = 0.0, hi = 20.0; for (int i = 0; i < 40; ++i) { const double mid = 0.5 * (lo + hi); if (minkowski(a1, a2, u * mid).depth > 0) lo = mid; else hi = mid; }
            const double dist0 = hi + 0.8;
            for (int model = 0; model < 3; ++model)
            {
               Body b1{ { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, a1 }, b2{ u * dist0, u * -vrel, { 0, 0, 0 }, a2 };
               const double E0 = energy(b1) + energy(b2); V warm = u; int quiet = 0; bool hadContact = false;
               for (int step = 0; step < 12000 && quiet < 40; ++step)
               {
                  const Force f = model == 0 ? eccForce(b1.x, b1.a, b2.x, b2.a, k, &warm) : simForce(b1.x, b1.a, b2.x, b2.a, k, model == 2);
                  const bool separating = dot(b2.x - b1.x, b2.v - b1.v) > 0; if (f.hit) { hadContact = true; res[model].maxDepth = std::max(res[model].maxDepth, f.depth); }
                  quiet = (!f.hit && separating) ? quiet + 1 : 0;
                  // velocity Verlet for translation; angular momentum kick, rigid rotation of the axis about omega
                  b1.v += f.F1 * (0.5 * dt / kMass); b2.v += f.F2 * (0.5 * dt / kMass); b1.L += f.T1 * (0.5 * dt); b2.L += f.T2 * (0.5 * dt);
                  b1.x += b1.v * dt; b2.x += b2.v * dt; b1.a = unit(rodrigues(b1.a, omegaOf(b1) * dt)); b2.a = unit(rodrigues(b2.a, omegaOf(b2) * dt));
                  const Force g = model == 0 ? eccForce(b1.x, b1.a, b2.x, b2.a, k, &warm) : simForce(b1.x, b1.a, b2.x, b2.a, k, model == 2);
                  b1.v += g.F1 * (0.5 * dt / kMass); b2.v += g.F2 * (0.5 * dt / kMass); b1.L += g.T1 * (0.5 * dt); b2.L += g.T2 * (0.5 * dt);
               }
               const double dE = (energy(b1) + energy(b2) - E0) / E0; res[model].contacts += hadContact;
               res[model].mean += dE; res[model].meanAbs += std::abs(dE); res[model].worst = std::max(res[model].worst, std::abs(dE)); ++res[model].n;
            }
         }
         std::printf("     |v_rel| = %.2f (peak overlap ~ %.1f LU):\n", vrel, vrel * std::sqrt(meff / k));
         for (int m = 0; m < 3; ++m) std::printf("       %-52s mean dE/E %+8.5f  mean|dE/E| %8.5f  worst %8.5f | impacts with contact %d/%d, deepest %.2f LU\n", names[m], res[m].mean / res[m].n, res[m].meanAbs / res[m].n, res[m].worst, res[m].contacts, res[m].n, res[m].maxDepth);
         check(res[0].meanAbs / res[0].n < 0.02, "the ECC model must conserve energy (mean |dE/E| < 2 %)");
         check(res[1].meanAbs / res[1].n < 1e-3, "GJK/EPA + spring must conserve energy (mean |dE/E| < 0.1 %): it realises the Minkowski/ECC model");
         if (vrel >= 0.15) check(res[2].meanAbs / res[2].n > 5.0 * res[1].meanAbs / res[1].n && res[2].mean / res[2].n > 0, "the closed-form hybrid must create (net positive) energy in deep impacts, far more than GJK/EPA");
      }
   }

   // ---------------------------------------------------------------- (4) shape-specific effective mass at the ECC contact point (what useRotationalEffectiveMass would use) vs the translational reduced mass the stiffness is calibrated with
   {
      auto invIw = [&](const V& a, const V& r) { const V ra = r - a * dot(r, a); return dot(r, r) < 1e-30 ? 0.0 : (dot(ra, ra) / kIperp + dot(r, a) * dot(r, a) / kIax); };   // r^T I_w^-1 r for a vector r
      std::vector< double > ratio; int used = 0;
      for (int t = 0; t < 1500 && used < 400; ++t)
      {
         const V a1 = randAxis(), a2 = randAxis(); const V d = unit(V{ nd(gen), nd(gen), nd(gen) }) * (2.0 + 6.0 * ud(gen));
         const Mink m = minkowski(a1, a2, d); if (m.depth < 0.02 || m.depth > 0.6) continue;
         V p1, p2; witnessPoints({ 0, 0, 0 }, a1, d, a2, m.n, p1, p2); const V pc = (p1 + p2) * 0.5;
         const double inv = 2.0 / kMass + invIw(a1, cross(pc, m.n)) + invIw(a2, cross(pc - d, m.n)); ratio.push_back((1.0 / inv) / meff); ++used;
      }
      std::sort(ratio.begin(), ratio.end());
      auto q = [&](double f) { return ratio[std::min(ratio.size() - 1, size_t(f * double(ratio.size())))]; };
      std::printf("(4) generalised effective mass m_eff = 1/(1/mA + 1/mB + (rA x n)^T IA^-1 (rA x n) + (rB x n)^T IB^-1 (rB x n)) at the ECC contact point / translational reduced mass m/2 (the value the stiffness is calibrated with), %zu shallow rod-rod contacts:\n"
                  "    ratio quantiles: min %.3f | 10%% %.3f | median %.3f | 90%% %.3f | max %.3f  -> the contact period for the same k is sqrt(ratio) x tc; median contact %.0f%% shorter than intended\n", ratio.size(), ratio.front(), q(0.1), q(0.5), q(0.9), ratio.back(), 100.0 * (1.0 - std::sqrt(q(0.5))));
      check(ratio.size() >= 40 && q(0.5) < 0.95 && ratio.front() < 0.5, "rod contacts must have a rotational contribution that lowers the effective mass noticeably (sphere-sphere would give exactly 1)");
   }

   std::cout << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failures)" << std::endl;
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
