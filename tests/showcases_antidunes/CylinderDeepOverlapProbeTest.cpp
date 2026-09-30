//======================================================================================================================
//! \file CylinderDeepOverlapProbeTest.cpp
//! Diagnostic: rod pairs (r=3, L=18) taken from the floater-history of the cylinder-hopper resume run (job 4250907, steps 5650/5850) that are
//! interpenetrated by 30-90 % of their volume (Monte-Carlo volume check in the plan) -- what do the analytic hybrid dispatch and plain GJK/EPA report?
//======================================================================================================================
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <iostream>

namespace walberla { namespace antidunes {
using Vec3 = mesa_pd::Vec3;
struct Pose { real_t x, y, z, qw, qx, qy, qz; };
mesa_pd::Rot3 rotFrom(const Pose& p)
{
   const real_t n = std::sqrt(p.qw * p.qw + p.qx * p.qx + p.qy * p.qy + p.qz * p.qz);
   return mesa_pd::Rot3(mesa_pd::Quat(p.qw / n, p.qx / n, p.qy / n, p.qz / n));
}
void probe(const char* name, const Pose& a, const Pose& b)
{
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   for (const Pose* p : { &a, &b }) { auto c = ps->create(); c->setPosition(Vec3(p->x, p->y, p->z)); c->setRotation(rotFrom(*p)); }
   const mesa_pd::data::Cylinder shape(real_t(3), real_t(18));
   const auto ci = paper_analytical::CylinderState{ Vec3(a.x, a.y, a.z), rotFrom(a).getMatrix(), real_t(3), real_t(18) };
   const auto cj = paper_analytical::CylinderState{ Vec3(b.x, b.y, b.z), rotFrom(b).getMatrix(), real_t(3), real_t(18) };
   const auto raw = paper_analytical::cylinderCylinderContact(ci, cj);
   std::cout << name << ": centre distance " << (Vec3(a.x, a.y, a.z) - Vec3(b.x, b.y, b.z)).length() << "\n   raw closed-form: contact=" << raw.contact << " scenario=" << raw.scenario << " overlap=" << raw.overlap << " normal=" << raw.normal << "\n";
   for (int mode = 0; mode < 2; ++mode)
   {
      SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = (mode == 1); scd.useTightLayer1 = true; scd.shadowAnalyticLog = (mode == 0);
      const bool hit = scd(0, 1, shape, shape, ac);
      if (mode == 0) for (size_t i = 0; i < SolidCylinderContactDetection::kShadowScenarios; ++i) { const auto& b = scd.shadowStats.byScenario[i]; if (b.n) std::cout << "   [shadow, GJK used] closed form would have said: " << SolidCylinderContactDetection::shadowScenarioName(i) << " wrong=" << b.wrong << " checkFlagged=" << b.verifierFlagged << " |d|=" << b.sumAbsErr << " normalFlipped=" << b.normalFlipped << "\n"; }
      std::cout << "   dispatch analytic=" << mode << ": " << (hit ? "CONTACT" : "no contact");
      if (hit) std::cout << " penetrationDepth=" << scd.getPenetrationDepth() << " normal=" << scd.getContactNormal();
      std::cout << "\n";
   }
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   using namespace walberla::antidunes;
   probe("5650 rods 176/180 (67% overlap)", { 353.428, 10.0374, 33.9364, 0.675649, -0.632585, 0.378261, 0.0159054 }, { 352.427, 10.0585, 34.7087, 0.662445, -0.664601, -0.116257, -0.325508 });
   probe("5850 rods 176/180 (75% overlap)", { 353.161, 10.3119, 34.0849, 0.675306, -0.640524, 0.364316, -0.031058 }, { 352.46, 10.1747, 34.4936, 0.618026, -0.676647, -0.127839, -0.379275 });
   probe("5650 rods 161/165 (45% overlap)", { 372.913, 10.2225, 4.5596, 0.69175, -0.582072, 0.357119, -0.234819 }, { 371.679, 12.7248, 3.88215, 0.698659, -0.60131, 0.387673, 0.00345134 });
   probe("5850 rods 782/786 (72% overlap)", { 365.368, 30.0509, 9.80322, 0.590926, -0.350542, -0.657838, -0.308507 }, { 364.46, 29.9638, 8.91965, 0.708915, -0.661348, -0.23094, 0.081992 });
   probe("5850 rods 571/572 (93% overlap)", { 297.125, 27.244, 21.3028, 0.703483, -0.666895, -0.233038, 0.0778158 }, { 297.231, 27.4026, 21.0605, 0.654417, -0.71362, -0.0481814, 0.245281 });
   return EXIT_SUCCESS;
}
