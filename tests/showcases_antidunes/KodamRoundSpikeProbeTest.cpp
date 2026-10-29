//======================================================================================================================
//! \file KodamRoundSpikeProbeTest.cpp
//! [Claude, 2026-10-10] Diagnostic for the energy spike of the round-container drop test (job 4292501, step 125000,
//! REV2-I-015): cylinders A (uid 2882303761517117441, upright on the floor at the container wall) and B (uid
//! 5764607523034234885, stacked on A) overlap by ~0.46 LU in the VTK output from step 118000 on. Poses from
//! vtk-out/floaterCylinder (VTK rows = R col0, R col2, -R col1; converted to quaternions by the scratch script poses.py).
//! For each output step: raw closed form, production dispatch (closed form + GJK/EPA fallback, useTightLayer1 as in the
//! prm) and plain GJK/EPA, plus the axial face-face overlap of the two (nearly coaxial) cylinders as a geometric check.
//======================================================================================================================
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

#include <cstdlib>
#include <iostream>

namespace walberla { namespace antidunes {
using Vec3 = mesa_pd::Vec3;
struct Pose { real_t x, y, z, qw, qx, qy, qz; };
struct Sample { int step; Pose a, b; };
const Sample samples[] = {
   { 118000, { 16.3973, 45.0097, 2.6765, 0.998129125, -0.000866025404, -0.00132287566, -0.0611207821 }, { 16.4584, 45.0455, 7.52879, 0.0126095202, 0.161458973, -0.986790251, 0.00412310563 } },
   { 118500, { 16.3971, 45.0094, 2.67713, 0.998131003, -0.001, -0.00141421356, -0.0610860049 }, { 16.4578, 45.0449, 7.52819, 0.0126392247, 0.161421033, -0.986795952, 0.00415331193 } },
   { 119000, { 16.3979, 45.0098, 2.67804, 0.998136639, -0.000866025404, -0.0015, -0.0609938521 }, { 16.4576, 45.0451, 7.52833, 0.0126095202, 0.161392379, -0.986800892, 0.00418330013 } },
   { 119500, { 16.3979, 45.0089, 2.67914, 0.998140396, -0.000866025404, -0.0016583124, -0.0609282365 }, { 16.4577, 45.045, 7.52835, 0.0126787223, 0.161479875, -0.986785564, 0.00421307489 } },
   { 120000, { 16.3974, 45.0093, 2.67891, 0.998144529, -0.001, -0.00158113883, -0.0608604962 }, { 16.4575, 45.0451, 7.52782, 0.0127964839, 0.161532504, -0.98677543, 0.00421307489 } },
   { 120500, { 16.3972, 45.0087, 2.67977, 0.998144403, -0.000866025404, -0.0016583124, -0.0608625501 }, { 16.4574, 45.0449, 7.52826, 0.0127671453, 0.161519349, -0.986777837, 0.00424264069 } },
   { 121000, { 16.3967, 45.0092, 2.6799, 0.998150415, -0.000866025404, -0.0016583124, -0.0607638873 }, { 16.4573, 45.0451, 7.52795, 0.012922848, 0.161618378, -0.986759596, 0.00424264069 } },
   { 121500, { 16.3957, 45.0088, 2.68005, 0.998148787, -0.001, -0.00158113883, -0.0607906243 }, { 16.4566, 45.0446, 7.52752, 0.0129711218, 0.161551076, -0.986769857, 0.00427200187 } },
   { 122000, { 16.3964, 45.0083, 2.68131, 0.998157052, -0.001, -0.00173205081, -0.0606506389 }, { 16.4566, 45.0446, 7.52854, 0.0128743932, 0.16155417, -0.986770363, 0.00433012702 } },
   { 122500, { 16.3948, 45.0077, 2.6815, 0.998155799, -0.000707106781, -0.00173205081, -0.0606753657 }, { 16.4561, 45.0442, 7.52776, 0.0130384048, 0.161588985, -0.986762636, 0.00430116263 } },
   { 123000, { 16.3941, 45.0074, 2.68229, 0.998161435, -0.000866025404, -0.00180277564, -0.0605784615 }, { 16.4551, 45.0439, 7.52814, 0.0130862523, 0.16154798, -0.986768337, 0.00438748219 } },
   { 123500, { 16.3938, 45.0067, 2.6836, 0.998162311, -0.001, -0.00187082869, -0.0605598877 }, { 16.4554, 45.044, 7.52703, 0.0133041347, 0.161672508, -0.986745155, 0.00435889894 } },
   { 124000, { 16.3928, 45.0069, 2.68313, 0.998165818, -0.001, -0.00173205081, -0.060506198 }, { 16.4553, 45.044, 7.52705, 0.0133322916, 0.161713481, -0.986737934, 0.00438748219 } },
   { 124500, { 16.3936, 45.0064, 2.68469, 0.998170702, -0.000866025404, -0.00193649167, -0.0604214366 }, { 16.4551, 45.0441, 7.52748, 0.013369742, 0.161745943, -0.986731853, 0.00444409721 } },
   { 125000, { 16.3919, 45.0063, 2.68563, 0.998178216, -0.000866025404, -0.00193649167, -0.0602971807 }, { 16.4542, 45.0436, 7.52746, 0.0134814688, 0.161744397, -0.986730333, 0.0045 } },
   { 125500, { 16.3906, 45.0251, 2.6457, 0.99838144, -0.00848528137, 0.00806225775, -0.0556551884 }, { 16.4558, 45.0224, 7.56399, 0.00339116499, 0.157332451, -0.987522911, 0.00578791845 } },
   { 126000, { 16.3455, 44.9703, 2.65396, 0.999309261, 0.00570087713, 0.000707106781, -0.0367151195 }, { 16.3663, 44.9879, 7.75197, 0.00304138127, 0.125732852, -0.99205834, -0.0015 } },
   { 126500, { 16.3852, 45.0103, 2.73577, 0.998615291, -0, -0.000707106781, -0.0526022813 }, { 16.4486, 45.0274, 7.59888, 0.0120208153, 0.126942901, -0.991836176, 0.00141421356 } },
   { 127000, { 16.3816, 45.0054, 2.6903, 0.998634943, -0.0005, -0.000866025404, -0.0522230792 }, { 16.3964, 45.0081, 7.57455, 0.00282842712, 0.126067442, -0.992017641, -5.26835606e-09 } },
};
mesa_pd::Rot3 rotFrom(const Pose& p)
{
   const real_t n = std::sqrt(p.qw * p.qw + p.qx * p.qx + p.qy * p.qy + p.qz * p.qz);
   return mesa_pd::Rot3(mesa_pd::Quat(p.qw / n, p.qx / n, p.qy / n, p.qz / n));
}
void probe(const Sample& s)
{
   const real_t r = real_t(4), L = real_t(5.3);
   auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
   for (const Pose* p : { &s.a, &s.b }) { auto c = ps->create(); c->setPosition(Vec3(p->x, p->y, p->z)); c->setRotation(rotFrom(*p)); }
   const mesa_pd::data::Cylinder shape(r, L);
   const auto Ra = rotFrom(s.a).getMatrix(), Rb = rotFrom(s.b).getMatrix();
   const Vec3 axA(Ra(0, 2), Ra(1, 2), Ra(2, 2)), axB(Rb(0, 2), Rb(1, 2), Rb(2, 2));
   const Vec3 d = Vec3(s.b.x, s.b.y, s.b.z) - Vec3(s.a.x, s.a.y, s.a.z);
   const real_t axial = std::abs(d * axA), lateral = (d - (d * axA) * axA).length();
   std::cout << "step " << s.step << ": axis angle " << std::acos(std::min(real_t(1), std::abs(axA * axB))) * real_t(180) / math::pi
             << " deg, centre distance along A's axis " << axial << " (face-face overlap ~ " << L - axial << "), lateral " << lateral << "\n";
   const auto ci = closed_form_contacts::CylinderState{ Vec3(s.a.x, s.a.y, s.a.z), Ra, r, L };
   const auto cj = closed_form_contacts::CylinderState{ Vec3(s.b.x, s.b.y, s.b.z), Rb, r, L };
   const auto raw = closed_form_contacts::cylinderCylinderContact(ci, cj);
   std::cout << "   raw closed form: contact=" << raw.contact << " scenario=" << raw.scenario << " overlap=" << raw.overlap << " normal=" << raw.normal << "\n";
   const auto ref = closed_form_contacts::cylinderCylinderContactRefined(ci, cj);
   std::cout << "   refined (refineCylinderContacts): contact=" << ref.contact << " scenario=" << ref.scenario << " overlap=" << ref.overlap << " normal=" << ref.normal << "\n";
   for (int mode = 0; mode < 2; ++mode)
   {
      SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = (mode == 1); scd.useTightLayer1 = true;
      const bool hit = scd(0, 1, shape, shape, ac);
      std::cout << "   " << (mode == 1 ? "production (closed form + fallback)" : "plain GJK/EPA                      ") << ": "
                << (hit ? "CONTACT" : "no contact");
      if (hit) std::cout << " depth=" << scd.getPenetrationDepth() << " normal=" << scd.getContactNormal() << " point=" << scd.getContactPoint();
      std::cout << "\n";
   }
}
// Sweep around the step-125000 pose: B's relative rotation to A scaled by f (0 = coaxial-parallel, 1 = as in the run) and
// B's lateral offset (perpendicular to A's axis) scaled by g; axial position unchanged. Prints scenario / depth / normal z
// of the production dispatch and the plain GJK/EPA depth.
void sweep()
{
   const Sample& s0 = samples[14]; // step 125000
   WALBERLA_CHECK_EQUAL(s0.step, 125000);
   const real_t r = real_t(4), L = real_t(5.3);
   const mesa_pd::Rot3 rotA = rotFrom(s0.a), rotB = rotFrom(s0.b);
   const mesa_pd::Quat qA = rotA.getQuaternion(), qB = rotB.getQuaternion();
   const mesa_pd::Quat qRel = qB * qA.getInverse();
   const real_t ang = qRel.getAngle();
   const mesa_pd::Vec3 relAxis = mesa_pd::Vec3(qRel[1], qRel[2], qRel[3]).getNormalized() * (qRel[0] < 0 ? real_t(-1) : real_t(1));
   const auto Ra = rotA.getMatrix();
   const Vec3 axA(Ra(0, 2), Ra(1, 2), Ra(2, 2));
   const Vec3 pA(s0.a.x, s0.a.y, s0.a.z), d = Vec3(s0.b.x, s0.b.y, s0.b.z) - pA;
   const Vec3 dAx = (d * axA) * axA, dLat = d - dAx;
   std::cout << "sweep at step 125000: relative rotation " << ang * real_t(180) / math::pi << " deg, lateral offset " << dLat.length() << "\n";
   const mesa_pd::data::Cylinder shape(r, L);
   // (1) B's own orientation (anti-parallel to A, own spin kept), its axis turned so that the angle to -axA is fDeg in the
   //     plane of the actual tilt; (2) the same with B flipped to be parallel to A (rotated 180 deg about a horizontal axis).
   const auto Rb0 = rotB.getMatrix();
   const Vec3 axB0(Rb0(0, 2), Rb0(1, 2), Rb0(2, 2));
   for (int flip = 0; flip < 2; ++flip)
      for (const real_t g : { real_t(0), real_t(1), real_t(3) })
         for (const real_t fDeg : { real_t(0), real_t(0.001), real_t(0.01), real_t(0.1), real_t(0.25), real_t(0.5), real_t(0.75), real_t(1), real_t(1.5), real_t(1.78), real_t(2.5), real_t(5) })
         {
            mesa_pd::Quat q0 = qB;
            Vec3 axB = axB0;
            if (flip) { q0 = mesa_pd::Quat(Vec3(1, 0, 0), math::pi) * qB; axB = -axB0; } // now parallel-ish to axA
            const Vec3 ref = (axB * axA > 0) ? axA : -axA;                               // the axis B should approach
            Vec3 perp = axB - (axB * ref) * ref;
            if (perp.length() < real_t(1e-12)) perp = Vec3(1, 0, 0) - (Vec3(1, 0, 0) * ref) * ref;
            perp = perp.getNormalized();
            const real_t th = fDeg * math::pi / real_t(180);
            const Vec3 target = std::cos(th) * ref + std::sin(th) * perp;
            // minimal rotation taking axB to target, applied on top of q0 (keeps the spin about the own axis)
            const Vec3 c = axB % target; const real_t sn = c.length(), cs = axB * target;
            const mesa_pd::Quat qCorr = (sn < real_t(1e-14)) ? mesa_pd::Quat() : mesa_pd::Quat(c / sn, std::atan2(sn, cs));
            const mesa_pd::Quat qBf = qCorr * q0;
            const Vec3 pB = pA + dAx + g * dLat;
            auto ps = std::make_shared< mesa_pd::data::ParticleStorage >(10); mesa_pd::data::ParticleAccessor ac(ps);
            auto c0 = ps->create(); c0->setPosition(pA); c0->setRotation(rotA);
            auto c1 = ps->create(); c1->setPosition(pB); c1->setRotation(mesa_pd::Rot3(qBf));
            const auto RbF = mesa_pd::Rot3(qBf).getMatrix();
            const auto raw = closed_form_contacts::cylinderCylinderContact(closed_form_contacts::CylinderState{ pA, Ra, r, L },
                                                                      closed_form_contacts::CylinderState{ pB, RbF, r, L });
            real_t depth[2] = { 0, 0 }; Vec3 n[2]; bool hit[2];
            for (int mode = 0; mode < 2; ++mode)
            {
               SolidCylinderContactDetection scd; scd.useAnalyticalContactDetection = (mode == 1); scd.useTightLayer1 = true;
               hit[mode] = scd(0, 1, shape, shape, ac);
               if (hit[mode]) { depth[mode] = scd.getPenetrationDepth(); n[mode] = scd.getContactNormal(); }
            }
            const Vec3 axBf(RbF(0, 2), RbF(1, 2), RbF(2, 2));
            std::cout << "   " << (flip ? "parallel     " : "anti-parallel") << " lat x" << g << " tilt " << fDeg << " deg (dot " << axA * axBf << "): "
                      << raw.scenario << " -> production depth " << (hit[1] ? depth[1] : real_t(0)) << " n_z " << (hit[1] ? n[1][2] : real_t(0))
                      << " | GJK depth " << (hit[0] ? depth[0] : real_t(0)) << "\n";
         }
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();
   walberla::antidunes::closed_form_contacts::kodamBandEdgeAxialCheck = (std::getenv("KODAM_AXIAL") != nullptr); // [Claude, 2026-10-10] the fix
   std::cout << "kodamBandEdgeAxialCheck = " << walberla::antidunes::closed_form_contacts::kodamBandEdgeAxialCheck << "\n";
   if (std::getenv("KODAM_SWEEP")) { walberla::antidunes::sweep(); return EXIT_SUCCESS; }
   for (const auto& s : walberla::antidunes::samples) walberla::antidunes::probe(s);
   return EXIT_SUCCESS;
}
