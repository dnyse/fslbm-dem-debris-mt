//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file AnalyticalSolutions.h
//! \brief Closed-form rigid-body-impact reference solutions for DEMValidation.h.
//!
//! ONE general solver (rigidBodyImpact, below) is applied to every shape rather than writing three
//! shape-specific formulas. This is the same standard single-point impulsive-contact-with-Coulomb-
//! friction method that the shape-specific literature results below each independently specialize:
//!
//!  - Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., & Wassgren, C. (2010). "Cylindrical object
//!    contact detection for use in discrete element method simulations. Part II - Experimental
//!    validation." Chemical Engineering Science, 65(22), 5863-5871. Their Test 1 (single cylinder
//!    dropped onto a flat wall at a prescribed tilt angle, compared to a closed-form analytical
//!    rebound) is exactly rigidBodyImpact() evaluated at the rim contact point their own Part I paper
//!    derives (reused here directly via paper_analytical::rimPointTowardPlane, PaperAnalyticalContacts.h).
//!  - Housner, G.W. (1963). "The behavior of inverted pendulum structures during earthquakes." Bulletin
//!    of the Seismological Society of America, 53(2), 403-417. The rocking-block restitution
//!    r = (omega_after/omega_before)^2 = (1 - 1.5*sin^2(alpha))^2 is the special case of rigidBodyImpact()
//!    with zero normal restitution (Housner's block does not bounce -- the impulsive event is the pivot
//!    transferring from one base corner to the other) and full-stick friction (no final slip at the new
//!    pivot corner), evaluated by conservation of angular momentum about that corner. selfCheckHousner()
//!    below asserts the general solver reproduces this closed form exactly, as a standing regression
//!    check on the general solver rather than a separate user-facing experiment.
//!  - The classic oblique sphere-wall impact (e.g. Walton, O.R. & Braun, R.L. (1986), or as used in
//!    Kruggel-Emden, H. et al. (2011), "Benchmark tests for verifying discrete element modelling codes
//!    at particle impact level," Granular Matter 13, 235-244) is rigidBodyImpact() with the trivial
//!    contact point r_c = -R*n (always directly below the sphere's centre along the wall normal).
//!
//! Restricted to PLANAR impacts (velocity, spin axis, and contact point all in one plane) throughout --
//! matching Kodam Part II's own cylinder test setup and Housner's rocking-block assumptions. Fully
//! general 3D oblique impacts are a harder, separately-studied problem (see e.g. Stronge, W.J. (2000),
//! "Impact Mechanics," Cambridge University Press) and out of scope here.
//======================================================================================================================

#pragma once

#include "PaperAnalyticalContacts.h"

#include <mesa_pd/data/DataTypes.h>

#include <core/math/Constants.h>
#include <core/math/Matrix3.h>

#include <algorithm>
#include <cmath>

namespace walberla
{
namespace antidunes
{
namespace analytical
{
using mesa_pd::Mat3;
using mesa_pd::Vec3;

//! Cross-product (skew-symmetric) matrix S such that S*x == r % x for any x.
inline Mat3 skew(const Vec3& r)
{
   return Mat3(real_t(0), -r[2], r[1], r[2], real_t(0), -r[0], -r[1], r[0], real_t(0));
}

//! Effective inverse-mass matrix K mapping an impulse J applied at contact point r (relative to this
//! body's centre of mass) to the resulting change in the CONTACT POINT's velocity: delta_v_c = K * J.
//! invInertiaWorld is the body's inverse inertia tensor already rotated into world frame (e.g. via
//! math::transformMatrixRART(R, invInertiaBF)). For a fixed/infinite-mass body (a wall), pass
//! invMass = 0 and invInertiaWorld = Mat3(0) -- K becomes the zero matrix, contributing nothing.
inline Mat3 effectiveInverseMassMatrix(real_t invMass, const Mat3& invInertiaWorld, const Vec3& r)
{
   const Mat3 S = skew(r);
   return Mat3::makeDiagonalMatrix(invMass) - S * invInertiaWorld * S;
}

//! One rigid body's kinematic + inertial state at the instant of contact.
struct BodyState
{
   real_t invMass          = real_t(0); // 0 for an infinite-mass (fixed) body, e.g. a wall
   Mat3   invInertiaWorld  = Mat3(real_t(0));
   Vec3   contactPoint     = Vec3(real_t(0)); // relative to this body's centre of mass, world frame
   Vec3   velocity         = Vec3(real_t(0)); // centre-of-mass linear velocity
   Vec3   angularVelocity  = Vec3(real_t(0));
};

struct ImpactResult
{
   Vec3 postVelocityA, postAngularVelocityA;
   Vec3 postVelocityB, postAngularVelocityB;
   Vec3 impulse;   // applied to A; -impulse applied to B
   bool stuck    = false; // true if the sticking (zero final slip) branch was used, else kinetic sliding
   bool contact  = false; // false if the bodies were already separating (vn >= 0) -- no impact computed
};

//! General single-point impulsive rigid-body collision with Newtonian normal restitution `e` and
//! Coulomb friction `mu`, following the standard two-branch impulse-momentum method (e.g. Baraff, D.
//! (2001), "Physically Based Modeling: Rigid Body Simulation," SIGGRAPH course notes, section on
//! collision response with friction; or Stronge (2000), ch. 6). `normal` points from B toward A (the
//! direction A is pushed to separate).
//!
//! CALIBRATION NOTE (found empirically running DEMValidationScenarios.h's wall-impact tests): this
//! model assumes kinetic friction saturates at mu*Jn for the ENTIRE duration sliding persists (the
//! textbook rigid-body-impact assumption). mesa_pd's actual contact law (LinearSpringDashpot) instead
//! uses a Cundall & Strack (1979)-style TANGENTIAL SPRING clamped by the Coulomb limit
//! (tangentialSpringDisplacement / fTLS in LinearSpringDashpot.h) -- over a short, finite, soft contact
//! duration that spring does not necessarily have time to ramp up to the full Coulomb-limited force for
//! the whole contact, so the simulated tangential impulse transferred (and hence induced spin) comes in
//! systematically BELOW this idealized-instantaneous-impulse prediction for fast oblique impacts --
//! observed ~15-25% low on angular velocity at tangential/normal speed ratios of 0.25-1.5, growing
//! roughly monotonically with the ratio, while the PURE NORMAL case (ratio=0, no friction involved)
//! matches to <1%. This is an expected, physically-understood difference between two different valid
//! tangential contact models, not a bug in either -- the wall-impact experiments' pass tolerance is set
//! accordingly (looser than the friction-free energy-conservation checks).
inline ImpactResult rigidBodyImpact(const BodyState& A, const BodyState& B, const Vec3& normal,
                                     real_t restitution, real_t friction)
{
   ImpactResult result;

   const Vec3 vcA = A.velocity + (A.angularVelocity % A.contactPoint);
   const Vec3 vcB = B.velocity + (B.angularVelocity % B.contactPoint);
   const Vec3 vrel = vcA - vcB;

   const Vec3   n  = normal.getNormalized();
   const real_t vn = vrel * n;
   if (vn >= real_t(0)) return result; // already separating -- no impact

   const Mat3 K = effectiveInverseMassMatrix(A.invMass, A.invInertiaWorld, A.contactPoint) +
                  effectiveInverseMassMatrix(B.invMass, B.invInertiaWorld, B.contactPoint);

   const Vec3   vtVec = vrel - vn * n;
   const real_t vt    = vtVec.length();

   Vec3 J;
   bool stuck = false;

   if (vt < real_t(1e-12) || friction <= real_t(0))
   {
      // No tangential slip (or frictionless): pure normal impulse.
      const real_t nKn = n * (K * n);
      const real_t Jn  = -(real_t(1) + restitution) * vn / nKn;
      J                = Jn * n;
   }
   else
   {
      const Vec3 t = vtVec / vt; // direction of pre-impact slip

      // --- Kinetic-sliding branch: friction impulse opposes slip at its Coulomb limit. ---
      const Vec3   dirSliding = n - friction * t;
      const real_t nKdir      = n * (K * dirSliding);

      // n*K*dirSliding <= 0 means the sliding-impulse direction has non-positive effective stiffness
      // along the normal -- the classic Coulomb-friction-impact indeterminacy/inconsistency for strong
      // rotational coupling + high friction (Painleve's paradox family; see Stronge (2000) ch. 6, or
      // Wang, Y. & Mason, M.T. (1992), ASME J. Applied Mechanics 59, on how planar frictional impact can
      // become ill-posed for some inertia/geometry/friction combinations -- large-lever-arm anisotropic
      // bodies like a box corner strike are exactly where this shows up; the sphere and cylinder tests
      // never hit it, r x n stays modest there). A naive sliding-branch solve would then produce a
      // negative or huge Jn (observed empirically: an UNPHYSICAL negative post-impact normal velocity,
      // i.e. penetrating further rather than separating). Skip straight to the sticking branch instead,
      // which stays well-posed as long as its own Knn (n*K*n, unconditionally positive for any real body)
      // is used as the leading coefficient.
      const bool slidingIllPosed = nKdir <= real_t(1e-10);

      Vec3   Jsliding(real_t(0));
      real_t vtFinalSliding = real_t(-1); // forces the sticking branch below when slidingIllPosed
      if (!slidingIllPosed)
      {
         const real_t JnSliding = -(real_t(1) + restitution) * vn / nKdir;
         Jsliding               = JnSliding * dirSliding;

         // Does the sliding-branch impulse over-correct the tangential velocity (i.e. would friction
         // actually have arrested slip before Jn finished ramping up)? If so, the sliding assumption is
         // inconsistent and the sticking branch applies instead.
         const Vec3 deltaVrelSliding = K * Jsliding;
         vtFinalSliding = t * (vtVec + deltaVrelSliding - vn * n); // component along original slip dir t
      }

      if (vtFinalSliding < real_t(0))
      {
         // --- Sticking branch: solve the 2x2 system enforcing v_n' = -e*vn AND v_t' = 0. ---
         const real_t Knn = n * (K * n);
         const real_t Knt = n * (K * t);
         const real_t Ktn = t * (K * n);
         const real_t Ktt = t * (K * t);
         const real_t rhsN = -(real_t(1) + restitution) * vn;
         const real_t rhsT = -vt;

         const real_t det = Knn * Ktt - Knt * Ktn;
         real_t Jn = real_t(0), Jt = real_t(0);
         if (std::abs(det) > real_t(1e-14))
         {
            Jn = (rhsN * Ktt - Knt * rhsT) / det;
            Jt = (Knn * rhsT - rhsN * Ktn) / det;
         }

         if (Jn > real_t(0) && std::abs(Jt) <= friction * Jn + real_t(1e-9) * std::max(real_t(1), Jn))
         {
            J     = Jn * n + Jt * t;
            stuck = true;
         }
         else if (!slidingIllPosed)
         {
            // Sticking would require more tangential impulse than Coulomb's limit allows -- the
            // sliding-branch result was actually consistent after all; fall back to it.
            J = Jsliding;
         }
         else
         {
            // Both branches degenerate (Painleve-paradox-adjacent geometry/friction combination): fall
            // back to a pure normal impulse (Knn = n*K*n is unconditionally positive for any physically
            // valid K, so this is always well-posed). Loses the frictional torque this step, but stays
            // physically sane (no spurious penetration) rather than propagating a nonsensical impulse.
            const real_t JnSafe = -(real_t(1) + restitution) * vn / Knn;
            J                   = JnSafe * n;
         }
      }
      else { J = Jsliding; }
   }

   const Vec3 deltaVrel = K * J;
   (void) deltaVrel;

   result.contact             = true;
   result.stuck                = stuck;
   result.impulse               = J;
   result.postVelocityA         = A.velocity + A.invMass * J;
   result.postAngularVelocityA  = A.angularVelocity + A.invInertiaWorld * (A.contactPoint % J);
   result.postVelocityB         = B.velocity - B.invMass * J;
   result.postAngularVelocityB  = B.angularVelocity - B.invInertiaWorld * (B.contactPoint % J);
   return result;
}

//! Convenience wrapper: single free body A against a fixed (infinite-mass) wall with normal `wallNormal`
//! (pointing away from the wall, into the domain). `contactPoint` is relative to A's centre of mass.
inline ImpactResult rigidBodyWallImpact(real_t massA, const Mat3& invInertiaWorldA, const Vec3& contactPoint,
                                         const Vec3& velocityA, const Vec3& angularVelocityA,
                                         const Vec3& wallNormal, real_t restitution, real_t friction)
{
   BodyState A;
   A.invMass         = real_t(1) / massA;
   A.invInertiaWorld = invInertiaWorldA;
   A.contactPoint    = contactPoint;
   A.velocity        = velocityA;
   A.angularVelocity = angularVelocityA;

   BodyState wall; // invMass = 0, invInertiaWorld = 0, at rest -- an infinite-mass fixed body
   return rigidBodyImpact(A, wall, wallNormal, restitution, friction);
}

//! Sphere-wall impact: the contact point is always -R*n (isotropic inertia, so invInertiaWorld is a
//! multiple of the identity regardless of orientation and needs no rotation tracking).
inline ImpactResult sphereWallImpact(real_t radius, real_t mass, const Vec3& velocity,
                                      const Vec3& angularVelocity, const Vec3& wallNormal,
                                      real_t restitution, real_t friction)
{
   const real_t I           = real_t(0.4) * mass * radius * radius; // 2/5 m R^2
   const Mat3   invInertia  = Mat3::makeDiagonalMatrix(real_t(1) / I);
   const Vec3   contactPoint = -radius * wallNormal.getNormalized();
   return rigidBodyWallImpact(mass, invInertia, contactPoint, velocity, angularVelocity, wallNormal,
                               restitution, friction);
}

//! Cylinder-wall impact, replicating Kodam et al. (2010) Part II Test 1: the contact point is the
//! point on the cylinder's rim closest to the wall plane, from paper_analytical::rimPointTowardPlane
//! (PaperAnalyticalContacts.h) -- reusing the not-yet-simulation-wired Kodam/Guo closed forms directly
//! as the reference geometry, exactly as that file's own header comment anticipates ("so a dedicated
//! test can compare these closed forms ... for identical geometry").
//!
//! FACE-PLANE / BAND-PLANE SPECIAL CASE: near axis-parallel (face-on) or axis-perpendicular (band-on)
//! orientations, the true contact is a whole flat face or the full band LENGTH touching the wall
//! simultaneously, not a single point -- rimPointTowardPlane still returns *a* correct-height point on
//! the rim, but an essentially ARBITRARY one off to one side, which would inject a spurious lever arm
//! (and hence spurious predicted spin) that the true symmetric, evenly-distributed contact does not
//! produce (confirmed empirically: DEMValidationScenarios.h's simulated angular velocity at tilt=0/90
//! is consistent with zero, matching the physically-expected symmetric case, not the nonzero value the
//! naive rim-point lever arm would predict). Mirrors paper_analytical::cylinderPlaneContact's own
//! face-plane/band-plane branch split, using the same axisDotN thresholds -- but here to pick the
//! physically-correct ON-AXIS effective contact point (zero net lever arm) instead.
inline ImpactResult cylinderWallImpact(const paper_analytical::CylinderState& cyl, real_t mass,
                                        const Mat3& invInertiaWorld, const Vec3& velocity,
                                        const Vec3& angularVelocity, const Vec3& wallPos,
                                        const Vec3& wallNormal, real_t restitution, real_t friction)
{
   const Vec3   nHat     = wallNormal.getNormalized();
   const Vec3   eZ       = cyl.axis();
   const real_t axisDotN = eZ * nHat;
   const real_t halfL    = real_t(0.5) * cyl.length;

   Vec3 contactPoint;
   if (std::abs(std::abs(axisDotN) - real_t(1)) < real_t(1e-6))
   {
      // face-plane: on-axis, offset by the half-length toward the wall
      contactPoint = -math::sign(axisDotN) * halfL * eZ;
   }
   else if (std::abs(axisDotN) < real_t(1e-6))
   {
      // band-plane: on the axis-perpendicular line through the centre, offset radially toward the wall
      contactPoint = -cyl.radius * nHat;
   }
   else
   {
      const Vec3 E = paper_analytical::detail::rimPointTowardPlane(cyl, wallPos, wallNormal);
      contactPoint = E - cyl.pos;
   }
   return rigidBodyWallImpact(mass, invInertiaWorld, contactPoint, velocity, angularVelocity, wallNormal,
                               restitution, friction);
}

//! Box-corner-first wall impact. `support` is the box's own support() function (its exact deepest-point
//! geometry, mesa_pd::data::Box::support) evaluated at the body-frame direction opposite the wall
//! normal -- reusing the identical corner-selection logic SolidCylinderContactDetection's Box x
//! HalfSpace overload already uses in production (ParticleKernels.h), so the analytical reference
//! targets the exact same contact point the DEM code under test will find.
inline ImpactResult boxCornerWallImpact(const Vec3& supportPointBodyFrame, const Mat3& R, real_t mass,
                                         const Mat3& invInertiaBF, const Vec3& velocity,
                                         const Vec3& angularVelocity, const Vec3& wallNormal,
                                         real_t restitution, real_t friction)
{
   const Mat3 invInertiaWorld = math::transformMatrixRART(R, invInertiaBF);
   const Vec3 contactPoint    = R * supportPointBodyFrame;
   return rigidBodyWallImpact(mass, invInertiaWorld, contactPoint, velocity, angularVelocity, wallNormal,
                               restitution, friction);
}

//! Housner (1963) rocking-block restitution: r = (omega_after / omega_before)^2, alpha = the block's
//! slenderness angle (angle between the base-corner-to-centroid diagonal and the vertical).
inline real_t housnerRockingRestitution(real_t slendernessAngleRad)
{
   const real_t s = std::sin(slendernessAngleRad);
   return (real_t(1) - real_t(1.5) * s * s) * (real_t(1) - real_t(1.5) * s * s);
}

//! Internal self-consistency check (not a user-facing experiment): Housner's pivot-transfer formula is
//! the special case of rigidBodyWallImpact with zero normal restitution (the block does not bounce --
//! the impulsive event is the pivot transferring from one base corner to the other) and full-stick
//! friction (mu large enough to force the sticking branch, so the new pivot corner's final velocity is
//! exactly zero), evaluated by angular-momentum conservation about the new pivot corner. Using
//! Housner's own reduced 2D parametrization directly (half-diagonal b, slenderness alpha, planar
//! rotation about Y) avoids needing a full 3D box just to exercise this identity. Returns the relative
//! error between the general solver's (omega_after/omega_before)^2 and housnerRockingRestitution(alpha);
//! should be ~0 (machine precision) if rigidBodyImpact() is implemented correctly.
inline real_t selfCheckHousnerRelativeError(real_t slendernessAngleRad)
{
   const real_t m = real_t(1), b = real_t(1);
   // Planar rectangular block, thin in the out-of-plane (Y) direction: I_O about a base corner O is the
   // standard rocking-block value I_O = (4/3) m b^2 (Housner 1963, eq. 1); about the centroid,
   // I_G = I_O - m*b^2 = (1/3) m b^2. Rotation axis is world Y throughout (planar problem).
   const real_t IG = (real_t(1) / real_t(3)) * m * b * b;
   const Mat3   invInertiaWorld = Mat3::makeDiagonalMatrix(real_t(1) / IG); // only the Y-component is used

   const real_t alpha = slendernessAngleRad;
   const real_t omegaBefore = real_t(1); // arbitrary reference magnitude -- the ratio is scale-invariant

   // Pre-impact: block rotating about the OLD corner O (position -b*(sin a, 0, -cos a) rel. to centroid,
   // i.e. down and to one side) with angular velocity omegaBefore about +Y. The NEW corner O' about to
   // touch down is at +b*(sin a, 0, -cos a) rel. to centroid. Centroid velocity from rigid rotation about
   // O: v_G = omega x (G - O) = omega*Y x (b*sin a, 0, -b*cos a).
   const Vec3 omegaVec(real_t(0), omegaBefore, real_t(0));
   const Vec3 rOldCorner(-b * std::sin(alpha), real_t(0), -b * std::cos(alpha));
   const Vec3 rNewCorner( b * std::sin(alpha), real_t(0), -b * std::cos(alpha));
   const Vec3 vG = omegaVec % (-rOldCorner); // omega x (G - O), G-O = -rOldCorner

   // Impulsive contact at the NEW corner, wall normal = +Z, e=0 (no bounce), mu effectively infinite
   // (guarantees the sticking branch is selected and enforced).
   const ImpactResult res = rigidBodyWallImpact(m, invInertiaWorld, rNewCorner, vG, omegaVec,
                                                 Vec3(real_t(0), real_t(0), real_t(1)), real_t(0),
                                                 real_t(1e6));

   const real_t omegaAfter = res.postAngularVelocityA[1];
   const real_t rGeneral   = (omegaAfter / omegaBefore) * (omegaAfter / omegaBefore);
   const real_t rHousner   = housnerRockingRestitution(alpha);
   return std::abs(rGeneral - rHousner) / std::max(real_t(1e-12), std::abs(rHousner));
}

//! Trivial two-body head-on restitution (r x n = 0 on both bodies -- central, face-aligned, no lever
//! arm, so rigidBodyImpact() collapses to the textbook v_rel_after = -e * v_rel_before regardless of
//! shape or inertia). Provided as the explicit reference value for the pairRestitution experiment.
inline real_t twoBodyHeadOnRestitution(real_t restitution) { return -restitution; }

} // namespace analytical
} // namespace antidunes
} // namespace walberla
