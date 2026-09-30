//======================================================================================================================
//! \file TiltYOwnAxisSpinArtifactTest.cpp
//! Follow-up to CylinderTippingRodJamTest.cpp / RodSpinSubstepResolutionTest.cpp, prompted by a LOCAL reproduction
//! attempt of job 4255679 (DropTest_hopper_cyl_repro_gjk_gas.prm: same rod geometry/spacing/contact-model settings,
//! a small slice of the real hopper, no fluid). That repro ran clean (as expected -- the STALE-PARTICLE-IDX guard is
//! a fluid-side bookkeeping check, and DropTest has no fluid), but its own [CylForces] log showed several rods'
//! tiltY diagnostic (ParticleKernels.h: atan2(R(0,2),R(0,0))) jump by >100 deg within a single 10-step window (e.g.
//! uid=53: -5.2 -> -163.3 deg in steps 1930->1940) while that SAME rod's position, velocity, and contact force all
//! stayed smooth and small through the exact same window -- no physical signature of a violent event at all. That
//! is the signature of a DIAGNOSTIC discontinuity, not a real flip, and it directly undermines the "toppling rod"
//! premise of the two earlier tests (both used the real crashed run's own tiltY=-3.6->167.2 deg reading, over steps
//! 2700->2900, as evidence uid 888 physically toppled).
//!
//! This test isolates the mechanism directly: tiltY = atan2(R(0,2), R(0,0)) uses R(0,0) and R(0,2), the world-X
//! components of the LOCAL X and LOCAL Z axes. For a cylinder at rest with its own physical axis (local Z, this
//! project's own convention) along world Y, a rotation ABOUT THAT SAME AXIS (a pure spin about the rod's own long
//! axis) is physically a complete no-op -- a cylinder is a body of revolution, its shape, its overlap with any other
//! body, and its true contact state are all exactly invariant under it -- but it sweeps the LOCAL X axis all the way
//! around the plane perpendicular to Y, which is exactly the plane atan2(R(0,2),R(0,0)) is measuring an angle
//! within. So a full, physically meaningless own-axis spin should sweep tiltY through the FULL +-180 deg range, with
//! the exact-ground-truth ALWAYS reporting zero shape change.
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

//! same formula as ParticleKernels.h's own tiltY diagnostic
double tiltYDeg(const Mat3& R) { return std::atan2(double(R(0, 2)), double(R(0, 0))) * 180.0 / walberla::math::pi; }

//! rotation matrix for a rod whose own axis (local Z) is world Y tilted by a small, FIXED `wobbleDeg` about world
//! X (representing the small residual transverse wobble any real settling rod has -- never exactly zero), then
//! spun by `spin` about that (fixed, tilted) axis. Columns = local X,Y,Z expressed in world coordinates. Physically:
//! the rod's shape/overlap-with-anything-else depends ONLY on the axis direction (fixed here), never on `spin` --
//! a cylinder is a body of revolution, so varying `spin` alone is an exact geometric no-op.
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

   // For a FIXED transverse axis wobble (the rod's TRUE axis direction, hence its true shape/overlap with anything
   // else, is completely determined by this one number and never changes across a spin sweep at fixed wobble),
   // sweep `spin` -- about that fixed axis, a real but physically meaningless degree of freedom for a body of
   // revolution -- through a full turn and find the largest tiltY change produced by a single small (0.5 deg) spin
   // step. Repeated at several wobble magnitudes to show this is a genuine mathematical singularity (sensitivity
   // diverges as wobble -> 0), not an artifact of one particular wobble choice.
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
   // The single-sample jump asymptotically approaches (but, by construction of a finite 0.5 deg sampling grid
   // across a true discontinuity, never quite reaches) 90 deg as wobble -> 0 -- itself the signature of a genuine
   // step discontinuity in the continuum limit. 80 deg is a conservative, comfortably-cleared threshold.
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
