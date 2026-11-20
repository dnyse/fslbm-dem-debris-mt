//======================================================================================================================
//! \file TiltYOwnAxisSpinArtifactTest.cpp
//! Companion to CylinderTippingRodJamTest.cpp / RodSpinSubstepResolutionTest.cpp, which read a logged tiltY change
//! of -3.6 -> 167.2 deg as a physical topple. A fluid-free drop-test reproduction with the same rod geometry and
//! contact settings showed tiltY jumps of > 100 deg within one 10-step window while the rod's position, velocity and
//! contact force stayed smooth: the signature of a discontinuity in the diagnostic, not a real flip.
//!
//! Mechanism: tiltY = atan2(R(0,2), R(0,0)) (ParticleKernels.h) uses the world-X components of the local X and local Z
//! axes. For a rod with its axis (local Z) along world Y, a spin about that axis is a geometric no-op (body of
//! revolution: shape, overlaps and contact state are invariant) but sweeps local X around the plane perpendicular to
//! Y, the plane in which tiltY measures its angle. So an own-axis spin sweeps tiltY through the full +-180 deg range.
//======================================================================================================================
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/Cylinder.h>

#include <core/Environment.h>
#include <core/math/Constants.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace walberla { namespace antidunes {
using Mat3 = mesa_pd::Mat3;
using Vec3 = mesa_pd::Vec3;

//! same formula as the tiltY diagnostic in ParticleKernels.h
double tiltYDeg(const Mat3& R) { return std::atan2(double(R(0, 2)), double(R(0, 0))) * 180.0 / walberla::math::pi; }

//! rotation matrix for a rod whose axis (local Z) is world Y tilted by a small fixed `wobbleDeg` (the residual
//! transverse wobble of a settling rod), then spun by `spin` about that tilted axis. Columns = local X,Y,Z in world
//! coordinates. The rod's geometry depends only on the axis direction, never on `spin`.
Mat3 rodAtSpin(double wobbleDeg, double spin)
{
   const double w = wobbleDeg * walberla::math::pi / 180.0;
   const Vec3 axis(real_t(std::sin(w)), real_t(std::cos(w)), real_t(0)); // world Y tilted toward world X by `wobbleDeg`
   // any fixed perpendicular reference direction, e.g. world Z projected off the axis (Gram-Schmidt)
   Vec3 ref(real_t(0), real_t(0), real_t(1));
   ref = ref - axis * (ref * axis); ref = ref / ref.length();
   const Vec3 perp = axis % ref; // completes the orthonormal triad
   const double c = std::cos(spin), s = std::sin(spin);
   const Vec3 localX = ref * real_t(c) + perp * real_t(s);
   const Vec3 localY = perp * real_t(c) - ref * real_t(s);
   Mat3 R;
   R(0, 0) = localX[0]; R(1, 0) = localX[1]; R(2, 0) = localX[2];
   R(0, 1) = localY[0]; R(1, 1) = localY[1]; R(2, 1) = localY[2];
   R(0, 2) = axis[0];   R(1, 2) = axis[1];   R(2, 2) = axis[2];
   return R;
}
}}

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv); WALBERLA_UNUSED(env);
   using namespace walberla::antidunes;
   int failures = 0;
   auto check = [&](bool ok, const std::string& what) { if (!ok) { ++failures; std::cout << "FAIL: " << what << std::endl; } };
   std::cout << std::fixed << std::setprecision(3);

   // At fixed wobble (fixed geometry), sweep `spin` in 0.5 deg steps and record the largest tiltY change per step.
   // Repeated for several wobbles to show a genuine singularity (sensitivity diverges as wobble -> 0).
   double worstJumpOverall = 0.0;
   for (double wobbleDeg : { 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.02 })
   {
      double minTilt = 1e300, maxTilt = -1e300, prevTilt = 0.0, biggestJump = 0.0, biggestJumpAtSpin = 0.0;
      for (int i = 0; i <= 720; ++i) // 0.5 deg steps, two full turns
      {
         const double spinDeg = i * 0.5;
         const double t = tiltYDeg(rodAtSpin(wobbleDeg, spinDeg * walberla::math::pi / 180.0));
         minTilt = std::min(minTilt, t); maxTilt = std::max(maxTilt, t);
         if (i > 0)
         {
            double d = std::abs(t - prevTilt); if (d > 180.0) d = 360.0 - d;
            if (d > biggestJump) { biggestJump = d; biggestJumpAtSpin = spinDeg; }
         }
         prevTilt = t;
      }
      worstJumpOverall = std::max(worstJumpOverall, biggestJump);
      std::cout << "wobble=" << std::setw(6) << wobbleDeg << " deg (TRUE geometry fixed):  tiltY range=["
                << std::setw(7) << minTilt << ", " << std::setw(7) << maxTilt << "] deg,  largest jump for a single "
                   "0.5 deg spin step=" << std::setw(7) << biggestJump << " deg (at spin~" << biggestJumpAtSpin << " deg)\n";
   }
   // R(0,0) = sin(spin) cos(wobble) vanishes exactly on grid points (spin = 0, 180, 360 deg), where tiltY = 90 deg, so
   // the single-step jump approaches, but never reaches, 90 deg as wobble -> 0 (a step discontinuity in the limit).
   // Threshold 80 deg.
   check(worstJumpOverall > 80.0, "as the axis wobble shrinks, a single 0.5 deg change in the (physically "
                                    "meaningless) spin angle must approach a step discontinuity in tiltY, "
                                    "consistent with the real repro's >100 deg jumps within one 10-step (i.e. "
                                    "tiny true-rotation) logging window");

   std::cout << (failures == 0 ?
      "\nALL CHECKS PASSED -- tiltY's own-axis-spin sensitivity fully explains the large single-window jumps seen in\n"
      "the local repro (DropTest_hopper_cyl_repro_gjk_gas.prm), and therefore the earlier interpretation of job\n"
      "4255679's real tiltY=-3.6->167.2 deg reading (steps 2700->2900) as evidence of a physical TOPPLE needs\n"
      "revisiting: it is equally consistent with a benign spin about the rod's own axis, which changes nothing\n"
      "about its shape, overlap, or contact state. CylinderTippingRodJamTest.cpp's transverse-tilt sweep is still\n"
      "mathematically valid on its own terms, but its real-world premise (that uid 888 physically toppled by that\n"
      "range) is now in doubt; the real run's own logged POSITION for uid 888 (smooth, ~0.02 LU/step throughout\n"
      "steps 2700-3100) supports this -- a genuine topple sweeping the rod's far end through a wide arc should also\n"
      "show up as position noise, and it does not.\n"
      : "\nSOME CHECKS FAILED -- see FAIL lines above.\n");
   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
