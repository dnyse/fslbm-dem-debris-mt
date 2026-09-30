//======================================================================================================================
//! \file RodSpinSubstepResolutionTest.cpp
//! Follow-up to CylinderTippingRodJamTest.cpp for the job 4255679 crash: that test showed GJK/EPA correctly tracked a
//! genuine, large geometric overlap as rod uid 3530822107858468888 toppled past its neighbour uid ...885 -- so the
//! open question is why the rod rotated fast enough, within the window bracketed by the [CylForces] log, to sweep
//! through that overlap-generating arc at all. This does NOT recover the (unlogged) instantaneous angular velocity --
//! no FloaterHistoryLogger was enabled for this run, and the tiltY diagnostic is only printed once per
//! evaluationFrequency (100 steps here) -- instead it makes the DEM contact model's own substep-resolution assumption
//! explicit and checks the real, logged data against it:
//!
//!   - the model's own calibrated collision duration is tc = collisionTimeFactor * avgParticleDiameter [main LBM
//!     steps] (SimulationParams.h's own documented formula: LinearSpringDashpot's stiffness is picked so a normal
//!     head-on collision between two average-diameter particles resolves in ~tc steps, independent of numSubCycles).
//!     For hopper_clog_cyl_frontblock_pidon_2node_halfrows_gas.prm (job 4255679's own prm): collisionTimeFactor is
//!     not overridden (default 4), avgParticleDiameterOverride=8.001943372119138, numSubCycles=10.
//!   - the real log brackets uid 888's tiltY at -3.634916 deg (step 2700) and +167.240066 deg (step 2900): the only
//!     two logged points either side of the flip, net ~170.87 deg over 200 main steps.
//!   - re-runs CylinderTippingRodJamTest.cpp's own tilt sweep (same real logged centre separation, rod A at rest) to
//!     find the theta-interval where the true overlap exceeds 1 LU -- the part of the rotation that actually matters.
//!
//! What this computes: IF the observed net rotation happened at a perfectly UNIFORM rate across the whole 200-step
//! bracket (the best case for the model -- any real burstiness only makes the true local rate faster), how many
//! sub-steps would the rod have spent crossing the dangerous arc, versus tc? This is a LOWER BOUND on how fast the
//! true, local rotation must have been; it cannot identify WHAT made the rod spin up (a lost supporting contact, an
//! upstream collision, ...), but it can show whether the bracket's own average rate is already marginal against the
//! model's own calibration, or whether a substantially bursty rotation is required to explain the force spike via
//! genuine under-resolution.
//======================================================================================================================
#include <core/Environment.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace walberla { namespace antidunes {

struct Vec3d
{
   double x, y, z;
   Vec3d operator-(const Vec3d& o) const { return { x - o.x, y - o.y, z - o.z }; }
   Vec3d operator*(double s) const { return { x * s, y * s, z * s }; }
   Vec3d cross(const Vec3d& o) const { return { y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x }; }
   double dot(const Vec3d& o) const { return x * o.x + y * o.y + z * o.z; }
   double length() const { return std::sqrt(dot(*this)); }
   Vec3d normalized() const { const double l = length(); return { x / l, y / l, z / l }; }
};
inline Vec3d operator+(const Vec3d& a, const Vec3d& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }

constexpr double kR = 3.0, kH = 9.0; // radius, half-length (rod: radius 3, length 18)
inline double support(const Vec3d& n, const Vec3d& a)
{
   const double c = std::abs(n.dot(a));
   return kH * c + kR * std::sqrt(std::max(0.0, 1.0 - c * c));
}
inline double needed(const Vec3d& n, const Vec3d& a1, const Vec3d& a2, const Vec3d& d) { return support(n, a1) + support(n, a2) - d.dot(n); }

//! same ground-truth machinery as CylinderClosedFormRangeTest.cpp / CylinderTippingRodJamTest.cpp, re-derived here
//! (no header shared between these self-contained diagnostic tests, matching this test directory's own pattern).
double exactPenetrationDepth(const Vec3d& a1, const Vec3d& a2, const Vec3d& d)
{
   constexpr int N = 4000;
   std::vector< std::pair< double, Vec3d > > c; c.reserve(N);
   const double ga = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
   for (int i = 0; i < N; ++i)
   {
      const double z = 1.0 - 2.0 * (i + 0.5) / N, r = std::sqrt(1.0 - z * z), ph = ga * i;
      Vec3d n{ r * std::cos(ph), r * std::sin(ph), z };
      c.emplace_back(needed(n, a1, a2, d), n);
   }
   std::partial_sort(c.begin(), c.begin() + 8, c.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
   double best = 1e300; std::mt19937 gen(1); std::normal_distribution< double > nd(0.0, 1.0);
   for (int k = 0; k < 8; ++k)
   {
      Vec3d n = c[size_t(k)].second; double v = c[size_t(k)].first;
      for (double step : { 0.03, 0.01, 0.003, 0.001, 0.0003 })
         for (int it = 0; it < 60; ++it)
         {
            Vec3d m{ n.x + step * nd(gen), n.y + step * nd(gen), n.z + step * nd(gen) };
            m = m.normalized();
            const double w = needed(m, a1, a2, d);
            if (w < v) { v = w; n = m; }
         }
      best = std::min(best, v);
   }
   return best;
}
Vec3d rotateAbout(const Vec3d& v, const Vec3d& k, double rad)
{
   return v * std::cos(rad) + k.cross(v) * std::sin(rad) + k * (k.dot(v) * (1.0 - std::cos(rad)));
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   using namespace walberla::antidunes;
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };
   std::cout << std::fixed << std::setprecision(3);

   // ---- the DEM model's own calibrated collision duration, from job 4255679's own prm ----
   const double collisionTimeFactor = 4.0;               // default (SimulationParams.h), not overridden in this prm
   const double avgParticleDiameter = 8.001943372119138;  // avgParticleDiameterOverride
   const double numSubCycles        = 10.0;               // this prm's own numSubCycles
   const double tcMainSteps = collisionTimeFactor * avgParticleDiameter;
   const double tcSubSteps  = tcMainSteps * numSubCycles;
   std::cout << "calibrated DEM collision duration tc = " << tcMainSteps << " main steps = " << tcSubSteps << " sub-steps\n";
   std::cout << "  (collisionTimeFactor=" << collisionTimeFactor << " * avgParticleDiameter=" << avgParticleDiameter
             << ", numSubCycles=" << numSubCycles << ")\n\n";
   check(tcSubSteps > 0.0, "calibrated collision duration must be positive");

   // ---- real logged bracket for uid 888's tiltY flip (this run's own [CylForces] stdout, steps 2700/2900) ----
   const double thetaStart = -3.634916, thetaEnd = 167.240066;
   const double stepStart = 2700.0, stepEnd = 2900.0;
   const double dTheta = thetaEnd - thetaStart, dStep = stepEnd - stepStart;
   const double avgRateDegPerMainStep = dTheta / dStep;
   std::cout << "real logged bracket: tiltY " << thetaStart << " deg (step " << int(stepStart) << ") -> " << thetaEnd
             << " deg (step " << int(stepEnd) << "), net " << dTheta << " deg over " << dStep << " main steps\n";
   std::cout << "  implied UNIFORM-rate average: " << avgRateDegPerMainStep << " deg/main-step ("
             << avgRateDegPerMainStep / numSubCycles << " deg/sub-step)\n\n";
   check(dTheta > 90.0, "the logged bracket should show a large (>90 deg) net rotation -- sanity check on the pasted-in numbers");

   const Vec3d pa{ 287.687, 9.46882, 57.8675 }, pb{ 283.45, 10.1355, 61.9644 }; // real step-3100 positions, uids 885/888
   const Vec3d d = pb - pa;
   const Vec3d restAxis{ 0, 1, 0 };
   struct Sweep { const char* name; Vec3d rotAxis; };
   const Sweep sweeps[] = { { "about world X", Vec3d{ 1, 0, 0 } }, { "about world Z", Vec3d{ 0, 0, 1 } } };

   for (const auto& sw : sweeps)
   {
      double dangerLo = -1, dangerHi = -1;
      for (int i = 0; i <= 180; ++i)
      {
         const double theta = double(i), rad = theta * 3.14159265358979323846 / 180.0;
         const Vec3d axisB = rotateAbout(restAxis, sw.rotAxis, rad).normalized();
         const double ov = exactPenetrationDepth(restAxis, axisB, d);
         if (ov > 1.0) { if (dangerLo < 0) dangerLo = theta; dangerHi = theta; }
      }
      check(dangerLo >= 0.0, std::string(sw.name) + ": no theta in [0,180] reaches exact overlap > 1 LU -- sweep assumption failed");
      const double dangerWidth = dangerHi - dangerLo;
      std::cout << "=== rod B tilting " << sw.name << " ===\n";
      std::cout << "  dangerous arc (exact overlap > 1 LU): theta in [" << dangerLo << ", " << dangerHi << "] deg, width "
                << dangerWidth << " deg\n";

      const double crossingMainSteps = dangerWidth / avgRateDegPerMainStep;
      const double crossingSubSteps  = crossingMainSteps * numSubCycles;
      std::cout << "  best-case (uniform rate) time to cross it: " << crossingMainSteps << " main steps = "
                << crossingSubSteps << " sub-steps  (" << (crossingSubSteps / tcSubSteps) << "x tc)\n";

      const double thresholdRateDegPerMainStep = dangerWidth / tcSubSteps * numSubCycles;
      std::cout << "  uniform rate that would cross it in exactly tc: " << thresholdRateDegPerMainStep << " deg/main-step\n";
      std::cout << "  observed 200-step average / that threshold: " << (avgRateDegPerMainStep / thresholdRateDegPerMainStep)
                << "x\n\n";
   }

   std::cout <<
      "Interpretation: if \"best-case crossing sub-steps\" > tc, the 200-step-AVERAGE rate alone still leaves the\n"
      "model comfortably resolved -- under-resolution can only explain the force spike if the TRUE local rotation\n"
      "was concentrated well above this bracket's own average (by roughly the printed x-tc ratio) for at least part\n"
      "of the 200-step window, i.e. a burst, not a smooth topple. If instead \"best-case crossing sub-steps\" <= tc,\n"
      "the average rate is already fast enough on its own -- no extra burstiness is required to explain it.\n";

   std::cout << (failures == 0 ? "\nALL CHECKS PASSED.\n" : "\nSOME CHECKS FAILED -- see FAIL lines above.\n");
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
