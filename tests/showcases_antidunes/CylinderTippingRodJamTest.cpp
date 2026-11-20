//======================================================================================================================
//! \file CylinderTippingRodJamTest.cpp
//! Reconstruction of a recorded crash in the half-density rod hopper with GJK/EPA-only cylinder-cylinder dispatch
//! (useAnalyticalContactDetection=false): the [STALE-PARTICLE-IDX] guard in SimpleBB.h aborted after the
//! hydrodynamic force/torque of two neighbouring rods jumped by roughly 10-30x within 100 steps, while one rod's
//! tiltY diagnostic (ParticleKernels.h, atan2(R(0,2),R(0,0))) sat at ~167-171 deg, i.e. nearly axis-reversed.
//! The exact quaternions were not logged, so the test approximates the configuration from what was:
//!   - rod geometry: radius 3, length 18;
//!   - the logged relative position at step 3100, magnitude ~5.93 LU (just under 2R = 6, barely touching if parallel);
//!   - rod A at rest orientation (axis = world Y; its tiltY stayed below 10 deg);
//!   - rod B swept through tilt theta in [0,180] deg about world X or Z, covering the observed range.
//!
//! Question: is the force/torque explosion a GJK/EPA correctness bug, or a correctly detected genuine deep overlap
//! (the rod tipped fast enough between DEM sub-cycles to penetrate deeply), i.e. a substep-resolution question?
//! Pass criterion: |GJK/EPA depth - exact depth| < 0.25 LU over the whole sweep.
//!
//! The exact-depth support-function machinery is duplicated from CylinderClosedFormRangeTest.cpp.
//======================================================================================================================
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

namespace walberla { namespace antidunes {
using Vec3 = mesa_pd::Vec3;
constexpr double kR = 3.0, kL = 18.0, kH = 9.0;

inline double support(const Vec3& n, const Vec3& a)
{
   const double c = std::abs(double(n * a));
   return kH * c + kR * std::sqrt(std::max(0.0, 1.0 - c * c));
}
inline double needed(const Vec3& n, const Vec3& a1, const Vec3& a2, const Vec3& d) { return support(n, a1) + support(n, a2) - double(d * n); }

//! ground truth: minimum over unit directions (Fibonacci grid, then pattern-search refinement of the best few), as in
//! CylinderClosedFormRangeTest.cpp.
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

struct Row { double theta; double exact; double closedForm; bool closedContact; const char* scenario; double gjk; };

Row evaluate(const Vec3& axisA, const Vec3& axisB, const Vec3& pa, const Vec3& pb, double theta)
{
   Row row{}; row.theta = theta;
   const Vec3 d = pb - pa;
   row.exact = exactPenetrationDepth(axisA, axisB, d);

   const auto ra = rotToAxis(axisA), rb = rotToAxis(axisB);
   const auto res = closed_form_contacts::cylinderCylinderContact({ pa, ra.getMatrix(), real_t(kR), real_t(kL) },
                                                                { pb, rb.getMatrix(), real_t(kR), real_t(kL) });
   row.closedContact = res.contact; row.closedForm = res.contact ? double(res.overlap) : 0.0; row.scenario = res.scenario;

   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   auto c0 = ps->create(); c0->setPosition(pa); c0->setRotation(ra);
   auto c1 = ps->create(); c1->setPosition(pb); c1->setRotation(rb);
   const mesa_pd::data::Cylinder shape{ real_t(kR), real_t(kL) };
   // GJK/EPA-only dispatch with tight layer 1, as in the recorded run
   SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = false; scd.useTightLayer1 = true;
   row.gjk = scd(0, 1, shape, shape, ac) ? -double(scd.getPenetrationDepth()) : 0.0;
   return row;
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::antidunes;
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };

   // Logged step-3100 positions ([CylForces] output).
   const Vec3 pa(real_t(287.687), real_t(9.46882), real_t(57.8675));  // rod A
   const Vec3 pb(real_t(283.45), real_t(10.1355), real_t(61.9644));   // rod B
   const double sep = double((pb - pa).length());
   std::cout << "real step-3100 centre separation: " << sep << " LU (2R = " << 2 * kR << ")\n\n";

   const Vec3 restAxis(real_t(0), real_t(1), real_t(0));
   struct Sweep { const char* name; Vec3 rotAxis; };
   const Sweep sweeps[] = { { "tilting about world X", Vec3(real_t(1), real_t(0), real_t(0)) },
                             { "tilting about world Z", Vec3(real_t(0), real_t(0), real_t(1)) } };

   for (const auto& sw : sweeps)
   {
      std::cout << "=== rod B " << sw.name << ", rod A fixed at rest (axis = world Y), real logged offset ===\n";
      std::cout << std::fixed << std::setprecision(3);
      double worstGjkErr = 0.0, worstGjkTheta = 0.0;
      double firstSevereTheta = -1.0; // first angle where exact overlap exceeds 1 LU (genuine deep penetration)
      for (int i = 0; i <= 90; ++i)
      {
         const double theta = i * 2.0; // degrees, 0..180
         const double rad = theta * 3.14159265358979323846 / 180.0;
         // rotate restAxis about sw.rotAxis by `rad` (Rodrigues' formula)
         const Vec3& k = sw.rotAxis;
         const Vec3 axisB = restAxis * real_t(std::cos(rad)) + (k % restAxis) * real_t(std::sin(rad)) +
                             k * real_t((k * restAxis) * (1.0 - std::cos(rad)));
         const Row row = evaluate(restAxis, axisB.getNormalized(), pa, pb, theta);
         const double gjkErr = std::abs(row.gjk - row.exact);
         if (gjkErr > worstGjkErr) { worstGjkErr = gjkErr; worstGjkTheta = theta; }
         if (firstSevereTheta < 0.0 && row.exact > 1.0) firstSevereTheta = theta;

         const bool nearReal = theta >= 165.0 && theta <= 172.0; // the actually-observed tiltY range
         if (i % 10 == 0 || nearReal)
            std::cout << "  theta=" << std::setw(6) << theta << "  exact=" << std::setw(7) << row.exact
                       << "  closedForm(" << row.scenario << ")=" << std::setw(7) << row.closedForm
                       << "  gjk=" << std::setw(7) << row.gjk << (nearReal ? "   <-- observed tiltY range" : "") << "\n";
      }
      std::cout << "  worst |gjk - exact| = " << worstGjkErr << " LU at theta=" << worstGjkTheta << " deg\n";
      std::cout << "  exact overlap first exceeds 1 LU at theta="
                 << (firstSevereTheta < 0.0 ? std::string("never (<=180 deg)") : std::to_string(firstSevereTheta)) << " deg\n\n";
      check(worstGjkErr < 0.25, std::string(sw.name) + ": GJK/EPA (the actual dispatch job 4255679 used) must stay within 0.25 LU of the exact depth across the whole tilt sweep");
   }

   std::cout << (failures == 0 ? "\nALL CHECKS PASSED -- GJK/EPA tracked the true penetration depth correctly throughout the tilt sweep.\n"
                                  "The force/torque explosion is consistent with a GENUINE deep overlap developing (the rod topples\n"
                                  "through a large fraction of its own length arc in a handful of steps), not a contact-detection bug.\n"
                                : "\nSOME CHECKS FAILED -- see FAIL lines above.\n");
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
