//======================================================================================================================
//! \file CylinderParallelGapTest.cpp
//! Regression test (asserting) for the cylinder hopper release abort ([RELEASE-OVERLAP] at step 2001, job 4248725): two PARALLEL
//! cylinders (radius 3, length 18, axis rotated onto world Y exactly as SetupStages.h does) that are clearly SEPARATED (surface gaps
//! 1, 2, 8, 15 LU) must NOT report a contact, with the analytic Kodam+Guo hybrid on and off.
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

namespace walberla { namespace antidunes {
using Vec3 = mesa_pd::Vec3;

// returns true if a (wrong) contact was reported
bool probe(real_t centreDistance, bool analytic, const Vec3& axisOffsetYZ, bool secondRotationNoise)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10);
   mesa_pd::data::ParticleAccessor ac(ps);

   // rotation exactly as SetupStages.h: rotate about (ez x axis) by acos(ez . axis), axis = world Y
   mesa_pd::Rot3 rot;
   const Vec3 ez(real_t(0), real_t(0), real_t(1)), axis(real_t(0), real_t(1), real_t(0));
   rot.rotate((ez % axis).getNormalized(), std::acos(std::clamp(ez * axis, real_t(-1), real_t(1))));

   auto c0 = ps->create();
   c0->setPosition(Vec3(real_t(249), real_t(30), real_t(49.13)));
   c0->setRotation(rot);
   auto c1 = ps->create();
   c1->setPosition(Vec3(real_t(249) - centreDistance, real_t(30) + axisOffsetYZ[1], real_t(49.13) + axisOffsetYZ[2]));
   mesa_pd::Rot3 rot2 = rot;
   if (secondRotationNoise) rot2.rotate(Vec3(real_t(0), real_t(0), real_t(1)), real_t(1e-9));
   c1->setRotation(rot2);

   const mesa_pd::data::Cylinder shape(real_t(3), real_t(18));
   SolidCylinderContactDetection scd;
   scd.useAnalyticalContactDetection = analytic;
   const bool hit = scd(0, 1, shape, shape, ac);
   std::cout << "  centreDist " << centreDistance << " (surface gap " << centreDistance - 6 << ") analytic=" << analytic
             << " noise=" << secondRotationNoise << " -> " << (hit ? "CONTACT" : "no contact");
   if (hit) std::cout << "  penetration=" << scd.getPenetrationDepth() << " normal=" << scd.getContactNormal() << " point=" << scd.getContactPoint();
   std::cout << std::endl;
   return hit;
}
}}


namespace walberla { namespace antidunes {
mesa_pd::Rot3 rotToAxis(const Vec3& axis)
{
   mesa_pd::Rot3 rot; const Vec3 ez(real_t(0), real_t(0), real_t(1));
   const real_t cosA = ez * axis;
   if (cosA < real_t(1) - real_t(1e-10))
   {
      const Vec3 rotAxis = (cosA > real_t(-1) + real_t(1e-10)) ? (ez % axis).getNormalized() : Vec3(real_t(1), real_t(0), real_t(0));
      rot.rotate(rotAxis, std::acos(std::clamp(cosA, real_t(-1), real_t(1))));
   }
   return rot;
}
// rod A at (249,30,49.13) with axis a; rod B at distance d along -x with axis b; prints result
void probeAxes(const char* na, const Vec3& a, const char* nb, const Vec3& b, real_t d, bool analytic)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   auto c0 = ps->create(); c0->setPosition(Vec3(real_t(249), real_t(30), real_t(49.13))); c0->setRotation(rotToAxis(a));
   auto c1 = ps->create(); c1->setPosition(Vec3(real_t(249) - d, real_t(30), real_t(49.13))); c1->setRotation(rotToAxis(b));
   const mesa_pd::data::Cylinder shape(real_t(3), real_t(18));
   SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = analytic;
   const bool hit = scd(0, 1, shape, shape, ac);
   std::cout << "AXES " << na << "/" << nb << " d=" << d << " analytic=" << analytic << " -> " << (hit ? "CONTACT" : "no contact");
   if (hit) std::cout << " pen=" << scd.getPenetrationDepth() << " n=" << scd.getContactNormal();
   std::cout << std::endl;
}
}}

#include <random>
namespace walberla { namespace antidunes {
// Monte Carlo: two parallel Y-axis rods (r=3, L=18) whose orientations/positions carry tiny perturbations like the running simulation
// (quaternion components ~1e-7, positions ~1e-6). Reports how often a SEPARATED pair is reported as touching.
int monteCarlo(real_t centreDistance, bool analytic, double rotSigma, double posSigma, int trials, std::mt19937& gen)
{
   std::normal_distribution< double > nd(0.0, 1.0);
   int hits = 0; double firstPen = 0; Vec3 firstN; bool have = false;
   for (int t = 0; t < trials; ++t)
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
      auto mk = [&](real_t x) {
         auto c = ps->create();
         c->setPosition(Vec3(x + real_t(posSigma * nd(gen)), real_t(30) + real_t(posSigma * nd(gen)), real_t(49.13) + real_t(posSigma * nd(gen))));
         mesa_pd::Rot3 r = rotToAxis(Vec3(real_t(0), real_t(1), real_t(0)));
         r.rotate(Vec3(real_t(1), real_t(0), real_t(0)), real_t(rotSigma * nd(gen)));
         r.rotate(Vec3(real_t(0), real_t(1), real_t(0)), real_t(rotSigma * nd(gen)));
         r.rotate(Vec3(real_t(0), real_t(0), real_t(1)), real_t(rotSigma * nd(gen)));
         c->setRotation(r);
      };
      mk(real_t(249)); mk(real_t(249) - centreDistance);
      const mesa_pd::data::Cylinder shape(real_t(3), real_t(18));
      SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = analytic;
      if (scd(0, 1, shape, shape, ac)) { ++hits; if (!have) { have = true; firstPen = scd.getPenetrationDepth(); firstN = scd.getContactNormal(); } }
   }
   std::cout << "MC d=" << centreDistance << " (gap " << centreDistance - 6 << ") analytic=" << analytic << " rotSigma=" << rotSigma << " posSigma=" << posSigma
             << " -> wrong contacts " << hits << "/" << trials; if (have) std::cout << "  first pen=" << firstPen << " n=" << firstN; std::cout << std::endl;
   return hits;
}
}}

namespace walberla { namespace antidunes {
// Accuracy near touching: analytic (hybrid) vs the GJK/EPA reference for two near-parallel Y-axis rods whose TRUE surface gap is (centreDistance - 6).
void touchingAccuracy(real_t centreDistance, double rotSigma, bool tightLayer1, int trials, std::mt19937& gen)
{
   std::normal_distribution< double > nd(0.0, 1.0);
   int hitA = 0, hitG = 0, disagree = 0; double worstErr = 0, sumErr = 0; int nErr = 0;
   for (int t = 0; t < trials; ++t)
   {
      auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
      auto mk = [&](real_t x) {
         auto c = ps->create();
         c->setPosition(Vec3(x + real_t(1e-6 * nd(gen)), real_t(30) + real_t(1e-6 * nd(gen)), real_t(49.13) + real_t(1e-6 * nd(gen))));
         mesa_pd::Rot3 r = rotToAxis(Vec3(real_t(0), real_t(1), real_t(0)));
         r.rotate(Vec3(real_t(1), real_t(0), real_t(0)), real_t(rotSigma * nd(gen)));
         r.rotate(Vec3(real_t(0), real_t(1), real_t(0)), real_t(rotSigma * nd(gen)));
         r.rotate(Vec3(real_t(0), real_t(0), real_t(1)), real_t(rotSigma * nd(gen)));
         c->setRotation(r);
      };
      mk(real_t(249)); mk(real_t(249) - centreDistance);
      const mesa_pd::data::Cylinder shape(real_t(3), real_t(18));
      SolidCylinderContactDetection a, g; a.useAnalyticalContactDetection = true; a.useTightLayer1 = tightLayer1; g.useAnalyticalContactDetection = false;
      const bool ha = a(0, 1, shape, shape, ac), hg = g(0, 1, shape, shape, ac);
      hitA += ha; hitG += hg; if (ha != hg) ++disagree;
      if (ha && hg) { const double e = std::abs(double(a.getPenetrationDepth()) - double(g.getPenetrationDepth())); worstErr = std::max(worstErr, e); sumErr += e; ++nErr; }
   }
   std::cout << "TOUCH d=" << centreDistance << " (true overlap " << 6 - centreDistance << ") rotSigma=" << rotSigma << " tightLayer1=" << tightLayer1
             << " -> analytic contact " << hitA << "/" << trials << ", GJK contact " << hitG << "/" << trials << ", disagree " << disagree
             << "; where both agree on contact: mean |pen diff| " << (nErr ? sumErr / nErr : 0.0) << " worst " << worstErr << std::endl;
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::antidunes;
   int wrong = 0;
   for (bool analytic : { false, true })
      for (bool noise : { false, true })
         for (real_t d : { real_t(7), real_t(8), real_t(14), real_t(21) })
            wrong += probe(d, analytic, Vec3(0, 0, 0), noise) ? 1 : 0;
   {
      const Vec3 X(1, 0, 0), Y(0, 1, 0), Z(0, 0, 1);
      const std::pair< const char*, Vec3 > ax[] = { { "X", X }, { "Y", Y }, { "Z", Z } };
      for (bool analytic : { false, true })
         for (const auto& A : ax) for (const auto& B : ax)
            for (real_t d : { real_t(7), real_t(14) }) probeAxes(A.first, A.second, B.first, B.second, d, analytic);
   }
   {
      std::mt19937 gen(12345);
      for (bool analytic : { true, false })
         for (double rs : { 0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5 })
            for (real_t d : { real_t(7), real_t(14) }) wrong += monteCarlo(d, analytic, rs, 1e-6, 400, gen);
   }
   {
      std::mt19937 gen(777);
      for (bool tight : { false, true })
         for (double rs : { 0.0, 1e-6, 1e-5 })
            for (real_t d : { real_t(5.0), real_t(5.5), real_t(5.9), real_t(6.1) }) touchingAccuracy(d, rs, tight, 200, gen);
   }
   std::cout << "wrong contacts (separated rods reported as touching): " << wrong << std::endl;
   return wrong == 0 ? EXIT_SUCCESS : EXIT_FAILURE;   // regression: any contact between these SEPARATED rods is a failure
}
