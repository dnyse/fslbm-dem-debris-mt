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
//! \file PaperAnalyticalContacts.h
//! \brief Closed-form (non-GJK/EPA) contact-detection algorithms taken directly from published
//!        papers, for the shape pairs SolidCylinderContactDetection (ParticleKernels.h) currently
//!        resolves via the generic GJK/EPA fallback: Box x Box and Cylinder x Cylinder.
//!
//! NOT used by the simulation. Nothing in BuildTimeloop.h / SetupStages.h / ParticleKernels.h
//! includes this header. It exists purely so a dedicated test (PaperAnalyticalContactTest.cpp) can
//! compare these closed forms against SolidCylinderContactDetection's existing GJK/EPA results for
//! identical geometry.
//!
//! -------------------------------------------------------------------------------------------------
//! Box x Box: Gottschalk, S., Lin, M. C., & Manocha, D. (1996). "OBBTree: A Hierarchical Structure
//! for Rapid Interference Detection." SIGGRAPH 1996 Proceedings, 171-180. Section 5, "Fast Overlap
//! Test for OBBs" (the 15-axis separating-axis-theorem test).
//!
//! The paper itself only derives a BOOLEAN overlap test (is there a separating axis, yes/no). The
//! minimum-translation-vector extension used here to also produce a penetration depth and contact
//! normal -- track the axis with the SMALLEST overlap instead of exiting on the first separating
//! axis found -- is the standard textbook extension of the same separating-axis theorem (see e.g.
//! Ericson, "Real-Time Collision Detection", 2004, ch. 4), not an additional result the 1996 paper
//! itself derives. A supplementary contact-POINT construction (also not from the paper) is provided
//! too, restricted to the one well-defined axis component -- the same restriction
//! ContactDetectionTest.cpp's testBoxBoxContact already applies to the GJK/EPA result, since a
//! flush face-face contact has no unique contact point in the first place (confirmed there by
//! cross-checking against pybullet, which picks a different, equally valid, point on the same
//! shared face).
//!
//! -------------------------------------------------------------------------------------------------
//! Cylinder x Cylinder / Cylinder x Plane: Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., &
//! Wassgren, C. (2010). "Cylindrical object contact detection for use in discrete element method
//! simulations. Part I - Contact detection algorithms." Chemical Engineering Science, 65(22),
//! 5852-5862. All "(eq. N)" comments below refer to this paper's equation numbers.
//!
//! Validation caveat (from the companion Part II paper, Chemical Engineering Science 65(22),
//! 5863-5871): Part II validates the cylinder-PLANE algorithm directly against a closed-form
//! analytical rebound formula (single cylinder bouncing off a wall). It does NOT validate
//! cylinder-CYLINDER contact geometry (overlap/normal/point) against any direct ground truth --
//! the only cylinder-cylinder-contact-bearing experiment (250 cylinders packed in a container) is
//! validated against an AGGREGATE bulk measurement (fill height, bulk solid fraction), not against
//! individual contact geometry. So "this is the published algorithm" is a stronger claim than "this
//! algorithm's cylinder-cylinder branch has been independently ground-truthed" -- keep that
//! distinction in mind when reading the test comparisons below.
//!
//! CORRECTIONS: Guo, Y., Wassgren, C., Ketterhagen, W., Hancock, B., & Curtis, J. (2012). "Some
//! computational considerations associated with discrete element modeling of cylindrical particles."
//! Powder Technology, 228, 193-198. An overlapping author team's own follow-up to Kodam et al. (2010),
//! documenting real classification bugs in the original face-edge/band-edge/band-band-skewed dispatch
//! logic ("undetected or wrongly detected contacts can lead to large overlaps") and one whole scenario
//! Kodam et al. never defined at all (their Fig. 8 "special contact"). This file now implements Guo et
//! al.'s corrected criteria (see rimVsFaceOrBand() and the band-band-skewed exclusion below) in place
//! of Kodam et al.'s original ones.
//!
//! HISTORY (kept for context, since this is exactly the kind of thing worth being able to point to):
//! an earlier version of this file used Kodam et al.'s ORIGINAL (uncorrected) face-edge/band-edge
//! dispatch logic and was demonstrated by PaperAnalyticalContactTest.cpp to produce false negatives
//! -- reporting no contact for genuinely overlapping, oblique-tilted cylinder pairs that GJK/EPA
//! correctly detects. Applying Guo et al.'s corrections (this version) targets exactly that class of
//! scenario-boundary misclassification; see PaperAnalyticalContactTest.cpp's own comments for which
//! specific previously-failing cases now pass and which, if any, still don't.
//!
//! Two sub-problems in the paper (band-edge, eq. 74-77; edge-edge point location, eq. 82-85) reduce
//! to a fourth-order polynomial in sin(theta). An earlier version of this file avoided hand-
//! transcribing that quartic (risk of an undetected sign/coefficient error with no compiler-
//! verified reference) and instead solved the equivalent 1D minimization numerically: a coarse
//! uniform scan over theta followed by golden-section refinement (up to ~1500 evaluations per
//! candidate contact). That numeric approach turned out to be a genuine performance problem at
//! scale -- dense scenarios with many simultaneous Cylinder x Cylinder contacts ran ~6x slower per
//! step than an equivalent Box x Box scenario purely because of this search.
//!
//! This version replaces it with the actual closed form: closestThetaOnCircleToLine derives the
//! quartic coefficients directly from the geometry (see that function's own comment for the
//! algebra) and solves them via the Durand-Kerner method (solveQuarticDurandKerner) -- chosen over
//! Ferrari's method specifically because it has no case-branching on a resolvent cubic's
//! discriminant to get subtly wrong; every root, real or complex, is refined by the same fixed
//! iteration loop. closestThetaOnCircleToPoint (used by the edge-edge alternating-projection loop)
//! turned out to have an even simpler exact answer -- a single atan2, no polynomial at all, since
//! its objective is a pure first-harmonic trig function of theta. Both closed forms were validated
//! against the previous golden-section implementation across thousands of randomized
//! configurations before replacing it (matching to floating-point precision in every trial), the
//! "compiler-verified reference" the original numeric approach didn't have.
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>

#include <core/math/Constants.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>

namespace walberla
{
namespace antidunes
{
namespace paper_analytical
{
using mesa_pd::Mat3;
using mesa_pd::Vec3;

//! Common return type for every function in this file. Overlap follows the PAPERS' OWN convention
//! (positive = penetrating), which is the OPPOSITE sign from SolidCylinderContactDetection's
//! penetrationDepth (negative = penetrating) -- callers/tests comparing the two must negate one of
//! them. Kept this way deliberately so the numbers here can be checked directly against the papers'
//! own equations without an extra sign translation layer while reading this file.
struct ContactResult
{
   bool        contact  = false;
   real_t      overlap  = real_t(0);
   Vec3        normal   = Vec3(real_t(0));
   Vec3        point    = Vec3(real_t(0));
   const char* scenario = "none";
};

//======================================================================================================================
// Box x Box -- Gottschalk, Lin & Manocha (1996), Section 5
//======================================================================================================================

//! \param posA,RA,halfA  box A: world position, body->world rotation matrix, half-extents (radii)
//! \param posB,RB,halfB  box B: same, for box B
//! Normal convention matches SolidCylinderContactDetection's pairwise overloads: points from B
//! (2nd argument) toward A (1st argument), e.g. testSphereSphereContact's contactNormal = -dir
//! where the 2nd sphere sits along +dir from the 1st.
inline ContactResult boxBoxSAT(const Vec3& posA, const Mat3& RA, const Vec3& halfA, const Vec3& posB,
                                const Mat3& RB, const Vec3& halfB)
{
   Vec3 axesA[3], axesB[3];
   for (uint_t i = 0; i < 3; ++i)
   {
      axesA[i] = Vec3(RA(0, i), RA(1, i), RA(2, i)); // columns of R = box's world-frame axes
      axesB[i] = Vec3(RB(0, i), RB(1, i), RB(2, i));
   }

   const Vec3 T = posB - posA;

   real_t minGap  = std::numeric_limits< real_t >::max();
   Vec3   bestAxis(real_t(1), real_t(0), real_t(0));
   real_t bestRa = real_t(0), bestRb = real_t(0);
   bool   found = false;

   // Returns false the instant a separating axis is found (boxes proven disjoint); otherwise folds
   // this axis's gap into the running minimum (the minimum-overlap axis is the contact normal /
   // penetration-depth choice -- see class-level comment on this being the standard MTD extension).
   auto testAxis = [&](Vec3 L) -> bool {
      const real_t len = L.length();
      if (len < real_t(1e-10)) return true; // degenerate (near-parallel edges) -- not a valid axis

      L = L / len;
      real_t ra = real_t(0), rb = real_t(0);
      for (uint_t i = 0; i < 3; ++i)
      {
         ra += halfA[i] * std::abs(axesA[i] * L);
         rb += halfB[i] * std::abs(axesB[i] * L);
      }
      const real_t dist = std::abs(T * L);
      const real_t gap  = ra + rb - dist; // > 0 while still overlapping along this axis

      if (gap < real_t(0)) return false; // separating axis -- disjoint

      if (gap < minGap)
      {
         minGap   = gap;
         bestAxis = L;
         bestRa   = ra;
         bestRb   = rb;
         found    = true;
      }
      return true;
   };

   for (int i = 0; i < 3; ++i)
      if (!testAxis(axesA[i])) return ContactResult{};
   for (int i = 0; i < 3; ++i)
      if (!testAxis(axesB[i])) return ContactResult{};
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         if (!testAxis(axesA[i] % axesB[j])) return ContactResult{};

   WALBERLA_ASSERT(found);

   const real_t sign   = (T * bestAxis >= real_t(0)) ? real_t(1) : real_t(-1);
   const Vec3   normal = -sign * bestAxis;

   ContactResult result;
   result.contact  = true;
   result.overlap  = minGap;
   result.normal   = normal;
   result.scenario = "SAT";
   // Supplementary contact point (NOT from the paper, see class-level comment): midpoint between
   // the two boxes' extreme points along the winning axis. Only the component along `bestAxis` is
   // geometrically well-defined for a flush face-face/edge-edge configuration.
   result.point = real_t(0.5) * (posA + posB) + real_t(0.5) * sign * (bestRa - bestRb) * bestAxis;
   return result;
}

//======================================================================================================================
// Cylinder x Cylinder / Cylinder x Plane -- Kodam, Bharadwaj, Curtis, Hancock & Wassgren (2010)
//======================================================================================================================

//! Minimal description of a solid, finite, right-circular cylinder needed by every Kodam et al.
//! formula: world position (centre of volume), body->world rotation (body z = cylinder axis),
//! radius R (paper's "band radius"), and full length L (paper's "face-to-face length").
struct CylinderState
{
   Vec3   pos;
   Mat3   R;
   real_t radius;
   real_t length;

   Vec3 axis() const { return Vec3(R(0, 2), R(1, 2), R(2, 2)); } // e_z in world = R * (0,0,1)
};

namespace detail
{
//! Closed-form global minimizer of theta -> |P(theta) - point|, where P(theta) = C + R*cos(theta)*u
//! + R*sin(theta)*v traces a circle of radius R centred at C in the plane spanned by (u, v).
//! |P(theta)-point|^2 = const + M*cos(theta) + N*sin(theta) (a pure first harmonic, since u,v are
//! orthonormal so the cross term vanishes) -- minimized in exactly one place, no search needed.
//! Used below both directly (closestThetaOnCircleToPoint) and as the point-case fallback inside
//! closestThetaOnCircleToLine when the line direction has no component in the (u,v) plane.
inline real_t closestThetaOnCircleToPointClosedForm(const Vec3& C, const Vec3& u, const Vec3& v,
                                                      real_t R, const Vec3& point)
{
   const Vec3   a = C - point;
   const real_t M = real_t(2) * R * (a * u);
   const real_t N = real_t(2) * R * (a * v);
   return std::atan2(-N, -M);
}

//! Real roots of c4*t^4 + c3*t^3 + c2*t^2 + c1*t + c0 = 0 via the Durand-Kerner method: all 4
//! (possibly complex) roots are refined simultaneously, so unlike Ferrari's method there is no
//! branching on the resolvent cubic's discriminant to get wrong. Fixed iteration count (matching
//! this file's existing fixed-iteration-count style, e.g. the old golden-section refinement) --
//! empirically converges to double precision within ~20-30 iterations for the geometries this is
//! used on; 60 leaves comfortable margin. Validated against the previous golden-section
//! implementation across thousands of random configurations (see git history / PR description) to
//! confirm this replacement produces the same answers, not just plausible-looking ones.
inline void solveQuarticDurandKerner(real_t c4, real_t c3, real_t c2, real_t c1, real_t c0,
                                      std::complex< real_t > roots[4])
{
   const real_t b3 = c3 / c4, b2 = c2 / c4, b1 = c1 / c4, b0 = c0 / c4;
   auto         poly = [&](std::complex< real_t > z) {
      return ((z + b3) * z + b2) * z * z + b1 * z + b0;
   };

   const real_t scale = real_t(1) + std::pow(std::max(std::max(std::abs(b3), std::abs(b2)),
                                                        std::max(std::abs(b1), std::abs(b0))),
                                              real_t(0.25));
   const std::complex< real_t > base(real_t(0.4), real_t(0.9)); // classic Durand-Kerner seed
   std::complex< real_t >       p = base;
   for (int k = 0; k < 4; ++k)
   {
      roots[k] = scale * p;
      p *= base;
   }

   for (int it = 0; it < 60; ++it)
   {
      std::complex< real_t > next[4];
      for (int i = 0; i < 4; ++i)
      {
         std::complex< real_t > denom(real_t(1), real_t(0));
         for (int j = 0; j < 4; ++j)
            if (j != i) denom *= (roots[i] - roots[j]);
         next[i] = roots[i] - poly(roots[i]) / denom;
      }
      for (int i = 0; i < 4; ++i)
         roots[i] = next[i];
   }
}

//! Closed-form global minimizer of theta -> |P(theta) - closest point on the given infinite line|,
//! where P(theta) = C + R*cos(theta)*u + R*sin(theta)*v traces a circle of radius R centred at C in
//! the plane spanned by (u, v). This is the paper's own quartic-in-sin(theta) (eq. 74-77 / 82-85,
//! see the file-level comment above): expanding |P(theta)-linePos|^2 - ((P(theta)-linePos).lineDir)^2
//! gives const + M'cos(theta) + N'sin(theta) - R^2*Amag^2*cos^2(theta-psi) (a mix of first- and
//! second-harmonic terms), whose stationary points satisfy R^2*Amag^2*sin(2*theta-2*psi) =
//! M'*sin(theta) - N'*cos(theta). Substituting t=tan(theta/2) (Weierstrass) and clearing
//! denominators turns this into a genuine quartic in t; solveQuarticDurandKerner finds its roots,
//! and each real root (plus theta=pi, the t->infinity point this substitution can't reach directly)
//! is evaluated against the true squared-distance to pick the global minimum among what are, in
//! general, 1-2 minima and 1-2 maxima of this function on the circle.
inline real_t closestThetaOnCircleToLine(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                          const Vec3& linePos, const Vec3& lineDir)
{
   const Vec3   a  = C - linePos;
   const real_t M  = real_t(2) * R * (a * u);
   const real_t N  = real_t(2) * R * (a * v);
   const real_t A0 = a * lineDir;
   const real_t A1 = u * lineDir;
   const real_t A2 = v * lineDir;

   // (The constant term K - A0^2 of f(theta) is omitted -- it doesn't affect where the minimum is,
   // only its value, and distSq() below recomputes the true squared distance directly anyway.)
   const real_t Mp = M - real_t(2) * A0 * R * A1;
   const real_t Np = N - real_t(2) * A0 * R * A2;

   auto distSq = [&](real_t theta) {
      const Vec3 P      = C + R * std::cos(theta) * u + R * std::sin(theta) * v;
      const Vec3 rel    = P - linePos;
      const Vec3 radial = rel - (rel * lineDir) * lineDir;
      return radial.sqrLength();
   };

   const real_t Amag2 = A1 * A1 + A2 * A2;
   real_t       bestTheta = math::pi;
   real_t       bestVal   = distSq(bestTheta);

   auto consider = [&](real_t theta) {
      const real_t val = distSq(theta);
      if (val < bestVal)
      {
         bestVal   = val;
         bestTheta = theta;
      }
   };

   if (Amag2 < real_t(1e-20))
   {
      // lineDir has (numerically) no component in the (u,v) plane -- the (rel.lineDir)^2 term is
      // then independent of theta, so this degenerates to the pure point-distance problem.
      consider(closestThetaOnCircleToPointClosedForm(C, u, v, R, linePos));
      return bestTheta;
   }

   const real_t psi = std::atan2(A2, A1);
   const real_t S2  = R * R * Amag2 * std::cos(real_t(2) * psi);
   const real_t S3  = R * R * Amag2 * std::sin(real_t(2) * psi);

   const real_t c4 = Np + S3;
   const real_t c3 = real_t(2) * Mp + real_t(4) * S2;
   const real_t c2 = real_t(-6) * S3;
   const real_t c1 = real_t(2) * Mp - real_t(4) * S2;
   const real_t c0 = S3 - Np;

   if (std::abs(c4) < real_t(1e-14) * std::max(real_t(1), std::abs(c0)))
   {
      // Leading coefficient (numerically) vanishes -- t=infinity (theta=pi, already in
      // bestTheta's initial candidate) is itself at least a double root; the remaining roots
      // solve the cubic-or-lower c3 t^3 + c2 t^2 + c1 t + c0 = 0. Reuse the quartic solver with
      // c4 nudged to an epsilon rather than special-casing a cubic solver separately: it just
      // means one of the four roots runs off towards a very large |t| (theta -> pi again), which
      // consider() naturally handles like any other candidate.
      const real_t            c4Safe = std::abs(c4) > real_t(1e-30) ? c4 : real_t(1e-14);
      std::complex< real_t >  roots[4];
      solveQuarticDurandKerner(c4Safe, c3, c2, c1, c0, roots);
      for (const auto& r : roots)
         if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
            consider(real_t(2) * std::atan(r.real()));
      return bestTheta;
   }

   std::complex< real_t > roots[4];
   solveQuarticDurandKerner(c4, c3, c2, c1, c0, roots);
   for (const auto& r : roots)
      if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
         consider(real_t(2) * std::atan(r.real()));

   return bestTheta;
}

//! MAXIMIZING counterpart to closestThetaOnCircleToLine -- identical derivation and quartic (setting
//! the derivative to zero finds ALL critical points of f(theta), minima and maxima alike; this just
//! keeps the largest instead of the smallest). Used by ComputeTunnelContacts (ParticleKernels.h) for
//! a Cylinder floater's exact farthest point from an arch tunnel's axis, per end-cap rim.
inline real_t farthestThetaOnCircleToLine(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                           const Vec3& linePos, const Vec3& lineDir)
{
   const Vec3   a  = C - linePos;
   const real_t M  = real_t(2) * R * (a * u);
   const real_t N  = real_t(2) * R * (a * v);
   const real_t A0 = a * lineDir;
   const real_t A1 = u * lineDir;
   const real_t A2 = v * lineDir;

   const real_t Mp = M - real_t(2) * A0 * R * A1;
   const real_t Np = N - real_t(2) * A0 * R * A2;

   auto distSq = [&](real_t theta) {
      const Vec3 P      = C + R * std::cos(theta) * u + R * std::sin(theta) * v;
      const Vec3 rel    = P - linePos;
      const Vec3 radial = rel - (rel * lineDir) * lineDir;
      return radial.sqrLength();
   };

   const real_t Amag2 = A1 * A1 + A2 * A2;
   real_t       bestTheta = math::pi;
   real_t       bestVal   = distSq(bestTheta);

   auto consider = [&](real_t theta) {
      const real_t val = distSq(theta);
      if (val > bestVal)
      {
         bestVal   = val;
         bestTheta = theta;
      }
   };

   if (Amag2 < real_t(1e-20))
   {
      // lineDir has no component in the (u,v) plane -- degenerates to point-distance; the FARTHEST
      // point from a point on a circle is diametrically opposite the CLOSEST one.
      consider(closestThetaOnCircleToPointClosedForm(C, u, v, R, linePos) + math::pi);
      return bestTheta;
   }

   const real_t psi = std::atan2(A2, A1);
   const real_t S2  = R * R * Amag2 * std::cos(real_t(2) * psi);
   const real_t S3  = R * R * Amag2 * std::sin(real_t(2) * psi);

   const real_t c4 = Np + S3;
   const real_t c3 = real_t(2) * Mp + real_t(4) * S2;
   const real_t c2 = real_t(-6) * S3;
   const real_t c1 = real_t(2) * Mp - real_t(4) * S2;
   const real_t c0 = S3 - Np;

   if (std::abs(c4) < real_t(1e-14) * std::max(real_t(1), std::abs(c0)))
   {
      const real_t            c4Safe = std::abs(c4) > real_t(1e-30) ? c4 : real_t(1e-14);
      std::complex< real_t >  roots[4];
      solveQuarticDurandKerner(c4Safe, c3, c2, c1, c0, roots);
      for (const auto& r : roots)
         if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
            consider(real_t(2) * std::atan(r.real()));
      return bestTheta;
   }

   std::complex< real_t > roots[4];
   solveQuarticDurandKerner(c4, c3, c2, c1, c0, roots);
   for (const auto& r : roots)
      if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
         consider(real_t(2) * std::atan(r.real()));

   return bestTheta;
}

//! Closed-form (no search) version of the same minimization as closestThetaOnCircleToLine, against
//! a fixed POINT instead of a line. Used by the edge-edge alternating-projection loop below.
inline real_t closestThetaOnCircleToPoint(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                           const Vec3& point)
{
   return closestThetaOnCircleToPointClosedForm(C, u, v, R, point);
}

//! Two arbitrary unit vectors spanning the plane perpendicular to `axis` (paper's u-hat/v-hat, eq.
//! 71-72), built more robustly than the paper's literal "cross with global z" (which degenerates
//! when axis IS global z) via a standard not-parallel helper vector.
inline void perpendicularBasis(const Vec3& axis, Vec3& u, Vec3& v)
{
   const Vec3 helper = (std::abs(axis[2]) < real_t(0.9)) ? Vec3(real_t(0), real_t(0), real_t(1))
                                                          : Vec3(real_t(1), real_t(0), real_t(0));
   u = (axis % helper).getNormalized();
   v = (axis % u).getNormalized();
}

//! Kodam et al. eq. 82-86: the two points (A, where `rim` crosses `target`'s infinite LATERAL
//! surface; B, where `rim` crosses `target`'s nearer end FACE plane) bounding the short arc of
//! `rim` that pokes into `target`'s finite capped-cylinder solid, if any. `rim` is a circle of
//! radius R centred at C, spanned by orthonormal (u, v) (NOT assumed tied to `target`'s axis the
//! way the paper's own eq. 71-72 basis is -- this file's `perpendicularBasis` picks an arbitrary
//! perpendicular pair, so both derivations below are the fully general versions of eq. 83 and
//! eq. 84-85 rather than the paper's axis-aligned-basis special case).
//!
//! Eq. 83 (Px^2+Py^2=R_target^2, "radial distance from rim point to target's axis equals
//! target's radius") expands -- exactly as the SAME expansion underlying closestThetaOnCircleToLine
//! above, just set equal to a target value instead of minimized -- into a genuine quartic in
//! t=tan(theta/2) after the Weierstrass substitution; solved via the same Durand-Kerner solver.
//! Verified by direct numerical brute-force root-scanning against this closed form (see PR
//! description / git history) to machine precision across thousands of random configurations.
//!
//! Eq. 84-85 (Pz=faceZ, "axial position of rim point equals target's face height") is
//! A1*cos(theta)+A2*sin(theta)=rhs, i.e. Amag*cos(theta-phi)=rhs -- solved directly via
//! phi+-acos(rhs/Amag), the general form of which the paper's own asin formula (eq. 85) is the
//! A1=0 special case (only true for the paper's own axis-tied basis choice, not this file's).
//!
//! Point A can have up to 4 candidate roots (plus theta=pi); point B always has exactly 2
//! (phi+-delta) when a solution exists at all. Neither the paper nor this file's predecessor
//! implementation specifies a tie-break when multiple candidates pass the containment check --
//! the choice here (nearest, angularly, to the rim's own closest point to target's axis) matches
//! the "near first contact, small segments" regime eq. 87's linearity assumption already commits
//! this whole scenario to.
struct RimCrossing
{
   bool contact = false;
   Vec3 A, B;
};

inline RimCrossing findRimCrossing(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                    const Vec3& targetPos, const Vec3& eTarget, real_t targetRadius,
                                    real_t targetHalfLength)
{
   RimCrossing out;

   const real_t thetaRef = closestThetaOnCircleToLine(C, u, v, R, targetPos, eTarget);

   auto pointAt = [&](real_t theta) { return C + R * std::cos(theta) * u + R * std::sin(theta) * v; };
   auto radialSqToAxis = [&](const Vec3& P) {
      const Vec3 rel    = P - targetPos;
      const Vec3 radial = rel - (rel * eTarget) * eTarget;
      return radial.sqrLength();
   };
   auto axialToAxis = [&](const Vec3& P) { return (P - targetPos) * eTarget; };
   auto angDiff = [](real_t a, real_t b) {
      real_t d = std::fmod(a - b + math::pi, real_t(2) * math::pi);
      if (d < real_t(0)) d += real_t(2) * math::pi;
      return std::abs(d - math::pi);
   };

   // --- Point A: rim crosses target's infinite lateral surface (eq. 83). ------------------------
   const Vec3   a  = C - targetPos;
   const real_t A0 = a * eTarget;
   const real_t A1 = u * eTarget;
   const real_t A2 = v * eTarget;
   const real_t Mp = real_t(2) * R * ((a * u) - A0 * A1);
   const real_t Np = real_t(2) * R * ((a * v) - A0 * A2);
   const real_t K2 = a.sqrLength() + R * R - A0 * A0 - targetRadius * targetRadius -
                      real_t(0.5) * R * R * (A1 * A1 + A2 * A2);
   const real_t P2 = real_t(-0.5) * R * R * (A1 * A1 - A2 * A2);
   const real_t Q2 = -R * R * A1 * A2;

   const real_t c4 = K2 - Mp + P2;
   const real_t c3 = real_t(2) * Np - real_t(4) * Q2;
   const real_t c2 = real_t(2) * K2 - real_t(6) * P2;
   const real_t c1 = real_t(2) * Np + real_t(4) * Q2;
   const real_t c0 = K2 + Mp + P2;

   bool   haveA = false;
   real_t bestDiffA = std::numeric_limits< real_t >::max();
   real_t thetaA     = real_t(0);

   auto considerA = [&](real_t theta) {
      const Vec3 P = pointAt(theta);
      if (std::abs(radialSqToAxis(P) - targetRadius * targetRadius) >
          real_t(1e-6) * std::max(real_t(1), targetRadius * targetRadius))
         return;
      if (std::abs(axialToAxis(P)) >= targetHalfLength) return;
      const real_t d = angDiff(theta, thetaRef);
      if (d < bestDiffA)
      {
         bestDiffA = d;
         thetaA    = theta;
         haveA     = true;
      }
   };

   if (std::abs(c4) > real_t(1e-12) * std::max(real_t(1), std::abs(c0)))
   {
      std::complex< real_t > roots[4];
      solveQuarticDurandKerner(c4, c3, c2, c1, c0, roots);
      for (const auto& r : roots)
         if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
            considerA(real_t(2) * std::atan(r.real()));
   }
   considerA(math::pi); // t -> infinity point the Weierstrass substitution can't reach directly

   if (!haveA) return out;
   out.A = pointAt(thetaA);

   // --- Point B: rim crosses target's nearer end face (eq. 84-85, general form). ------------------
   const real_t faceZ = math::sign(A0) * targetHalfLength; // nearer face, mirrors this file's near-cap picks elsewhere
   const real_t Amag  = std::sqrt(A1 * A1 + A2 * A2);
   if (Amag < real_t(1e-12)) return out; // rim lies in a plane parallel to target's axis -- no face crossing
   const real_t rhs = (faceZ - A0) / R;
   if (std::abs(rhs) > Amag) return out; // rim never reaches that face height

   const real_t phi   = std::atan2(A2, A1);
   const real_t delta = std::acos(std::clamp(rhs / Amag, real_t(-1), real_t(1)));

   bool   haveB = false;
   real_t bestDiffB = std::numeric_limits< real_t >::max();
   real_t thetaB     = real_t(0);

   auto considerB = [&](real_t theta) {
      const Vec3 P = pointAt(theta);
      if (radialSqToAxis(P) >= targetRadius * targetRadius) return;
      const real_t d = angDiff(theta, thetaRef);
      if (d < bestDiffB)
      {
         bestDiffB = d;
         thetaB    = theta;
         haveB     = true;
      }
   };
   considerB(phi + delta);
   considerB(phi - delta);

   if (!haveB) return out;
   out.B       = pointAt(thetaB);
   out.contact = true;
   return out;
}

//! The point on `cyl`'s rim closest to a plane through `planePos` with normal `planeNormal` --
//! shared closed-form building block for cylinder-plane edge-plane (eq. 15-17) AND, per Kodam et al.
//! Section 3.2.3 eq. 45-47, cylinder-cylinder face-edge (identical formula, with the "plane" being
//! the OTHER cylinder's flat face treated locally as a plane). CORRECTED per the derivation in
//! cylinderPlaneContact's edge-plane comment (an earlier version had the radial direction and a sign
//! wrong, caught by testCylinderPlaneEdge disagreeing with GJK/EPA by ~0.3).
inline Vec3 rimPointTowardPlane(const CylinderState& cyl, const Vec3& planePos, const Vec3& planeNormalIn)
{
   const Vec3   nHat     = planeNormalIn.getNormalized();
   const Vec3   eZ       = cyl.axis();
   const real_t axisDotN = eZ * nHat;
   const real_t zCyl     = (cyl.pos - planePos) * nHat;
   const real_t halfL    = real_t(0.5) * cyl.length;

   // e_r^P = e_z,cyl x (e_z,cyl x n) / |...| = -(component of n perpendicular to the axis, normalized)
   Vec3         radialAwayFromPlane = axisDotN * eZ - nHat;
   const real_t radialLen           = radialAwayFromPlane.length();
   radialAwayFromPlane = (radialLen > real_t(1e-10)) ? radialAwayFromPlane / radialLen : Vec3(real_t(0));

   return cyl.pos + math::sign(zCyl) * cyl.radius * radialAwayFromPlane -
          math::sign(zCyl) * math::sign(axisDotN) * halfL * eZ;
}
} // namespace detail

//! Cylinder-plane contact: face-plane, band-plane, edge-plane (eq. 2-21). Formulated directly in
//! world-space vectors rather than constructing the paper's separate plane-local rotation matrix
//! R^{PG} -- mathematically equivalent (R^{PG} is only ever used to take components along
//! `planeNormal` and within the plane), simpler to implement correctly without a compiler in the
//! loop to catch a transposed-matrix mistake.
inline ContactResult cylinderPlaneContact(const CylinderState& cyl, const Vec3& planePos,
                                           const Vec3& planeNormalIn)
{
   const Vec3   nHat      = planeNormalIn.getNormalized();
   const Vec3   eZ        = cyl.axis();
   const real_t axisDotN  = eZ * nHat;
   const real_t zCyl      = (cyl.pos - planePos) * nHat; // eq. 3, component along the plane normal
   const real_t halfL     = real_t(0.5) * cyl.length;
   const Vec3   centreOnPlane = cyl.pos - zCyl * nHat;     // cylinder centre projected onto the plane

   ContactResult result;

   // --- Face-plane (eq. 2): axis parallel to the plane normal -----------------------------------
   if (std::abs(std::abs(axisDotN) - real_t(1)) < real_t(1e-8) && std::abs(zCyl) < halfL)
   {
      result.contact  = true;
      result.overlap  = halfL - std::abs(zCyl);                                  // eq. 4
      result.normal   = math::sign(zCyl) * nHat;                                 // eq. 5
      result.point    = centreOnPlane - math::sign(zCyl) * real_t(0.5) * result.overlap * nHat; // eq. 6
      result.scenario = "face-plane";
      return result;
   }

   // --- Band-plane (eq. 9): axis perpendicular to the plane normal ------------------------------
   if (std::abs(axisDotN) < real_t(1e-8) && std::abs(zCyl) < cyl.radius)
   {
      result.contact  = true;
      result.overlap  = cyl.radius - std::abs(zCyl);                            // eq. 10
      result.normal   = math::sign(zCyl) * nHat;                                // eq. 11
      result.point    = centreOnPlane - math::sign(zCyl) * real_t(0.5) * result.overlap * nHat; // eq. 12
      result.scenario = "band-plane";
      return result;
   }

   // --- Edge-plane (eq. 13-21): general (tilted) orientation ------------------------------------
   // Point E, the point on the cylinder's rim closest to the plane (eq. 15-17), via the shared
   // detail::rimPointTowardPlane (see its comment for the corrected-formula derivation).
   const Vec3   E  = detail::rimPointTowardPlane(cyl, planePos, nHat);
   const real_t zE = (E - planePos) * nHat;

   if (zE * zCyl < real_t(0)) // eq. 18: opposite sides of the plane => contact
   {
      result.contact  = true;
      result.overlap  = std::abs(zE);                                          // eq. 19
      result.normal   = math::sign(zCyl) * nHat;                               // eq. 20
      const Vec3 eOnPlane = E - zE * nHat;
      result.point    = eOnPlane - math::sign(zE) * real_t(0.5) * result.overlap * nHat; // eq. 21
      result.scenario = "edge-plane";
   }
   return result;
}

//! Contact point shared by face-face (eq. 30-33) and band-band-parallel (paper: "same approach... as
//! Eqs. (30)-(33)", since band-band-parallel's overlap formula eq. 53 is literally identical to
//! face-face's radial overlap eq. 27). Only cylinder i's radius (Ri) and cylinder j's length (Li,
//! confusingly named after the paper's own i/j-swapped eq. 32) are needed -- deltaRadial and
//! deltaAxial already fold in the other cylinder's radius/length, so it is not a separate parameter.
inline Vec3 faceFaceStyleContactPoint(const Vec3& relBodyI, real_t Ri, real_t Li, real_t deltaAxial,
                                       real_t deltaRadial)
{
   const real_t x = relBodyI[0], y = relBodyI[1], z = relBodyI[2];
   const real_t rho = std::sqrt(x * x + y * y);
   const Vec3   rHat = (rho > real_t(1e-10)) ? Vec3(x, y, real_t(0)) / rho : Vec3(real_t(1), real_t(0), real_t(0));
   const Vec3   xyContact = (Ri - real_t(0.5) * deltaRadial) * rHat;                        // eq. 31
   const real_t zContact  = z - math::sign(z) * (real_t(0.5) * Li - real_t(0.5) * deltaAxial); // eq. 32
   return Vec3(xyContact[0], xyContact[1], zContact);                                        // body-i frame; caller transforms to world
}

//! Face-edge / special-contact / band-edge, all three tested against `host`'s face/band using a
//! candidate point on `probe`'s rim -- per Kodam et al. Section 3.2.3/3.2.5 AS CORRECTED by Guo,
//! Wassgren, Ketterhagen, Hancock & Curtis (2012), Powder Technology 228, 193-198, Table 1.
//!
//! Returned normal/point use the "natural" host->probe direction (like eq. 29/39/50/54/80 elsewhere
//! in this file); the caller negates when host corresponds to this module's 1st cylinderCylinderContact
//! argument, to match this file's toward-1st-argument convention (see the face-face branch's comment).
//!
//! Replaces an earlier version of this function that (a) used a single rim-circle search for BOTH
//! the face-edge and band-edge candidate points, when Kodam's face-edge actually has an exact closed
//! form (eq. 45-47, identical in structure to the now-corrected cylinder-plane edge-plane formula --
//! see detail::rimPointTowardPlane) rather than needing a numeric search at all; and (b) used Kodam's
//! ORIGINAL (pre-Guo-correction) face-edge/band-edge dispatch criteria, which PaperAnalyticalContactTest.cpp
//! demonstrated produce false negatives for oblique tilts. This version uses the closed form for the
//! face-edge candidate, the numeric search only for the band-edge/special-contact candidate (matching
//! Kodam's own eq. 67-77, which does need it), and Guo et al.'s corrected dispatch logic throughout.
inline ContactResult rimVsFaceOrBand(const CylinderState& host, const CylinderState& probe)
{
   ContactResult result;
   const real_t  halfLHost = real_t(0.5) * host.length;
   const Vec3    eHost     = host.axis();
   const Vec3    relProbeCentre = host.R.getTranspose() * (probe.pos - host.pos); // "O_j" in host frame

   // --- Face-edge candidate: E, the point on probe's rim closest to host's face (eq. 45-47), via the
   // same closed form as cylinder-plane edge-plane (detail::rimPointTowardPlane), with host's face
   // playing the role of the plane (position host.pos, normal eHost). -------------------------------
   const Vec3 E    = detail::rimPointTowardPlane(probe, host.pos, eHost);
   const Vec3 relE = host.R.getTranspose() * (E - host.pos);

   const bool eInFaceDisc = std::abs(relE[2]) < halfLHost &&
                             (relE[0] * relE[0] + relE[1] * relE[1]) < host.radius * host.radius;
   if (eInFaceDisc)
   {
      // Guo et al. Table 1, corrected face-edge criterion: E'.O_host' dot E'.O_probe' (projected
      // onto host's face plane, i.e. components 0/1 only) determines which of the two OR-branches
      // applies; the second branch additionally requires delta_a < delta_r (Kodam's original
      // condition, now only a fallback rather than the sole criterion).
      const real_t dot = (-relE[0]) * (relProbeCentre[0] - relE[0]) + (-relE[1]) * (relProbeCentre[1] - relE[1]);
      const real_t deltaA = halfLHost - std::abs(relE[2]);
      const real_t rhoE   = std::sqrt(relE[0] * relE[0] + relE[1] * relE[1]);
      const real_t deltaR = host.radius - rhoE;

      if (dot >= real_t(0) || deltaA < deltaR)
      {
         result.contact  = true;
         result.overlap  = deltaA; // eq. 49, unchanged by Guo's correction (only the dispatch changed)
         result.normal   = host.R * Vec3(real_t(0), real_t(0), math::sign(relE[2]));
         result.point    = E - result.normal * real_t(0.5) * deltaA;
         result.scenario = "face-edge";
         return result;
      }
   }

   // --- Band-edge / special-contact candidate: P, the point on probe's rim closest to host's BAND
   // (Kodam eq. 67-77, a quartic in sin(theta) solved in closed form by closestThetaOnCircleToLine
   // -- see that function's own comment for the derivation). --------------------------------------
   Vec3         uP, vP;
   detail::perpendicularBasis(probe.axis(), uP, vP);
   const Vec3   capTop      = probe.pos + real_t(0.5) * probe.length * probe.axis();
   const Vec3   capBottom   = probe.pos - real_t(0.5) * probe.length * probe.axis();
   const Vec3   capCentre   = ((capTop - host.pos).length() < (capBottom - host.pos).length()) ? capTop : capBottom;
   const real_t theta       = detail::closestThetaOnCircleToLine(capCentre, uP, vP, probe.radius, host.pos, eHost);
   const Vec3   P           = capCentre + probe.radius * std::cos(theta) * uP + probe.radius * std::sin(theta) * vP;
   const Vec3   relP        = host.R.getTranspose() * (P - host.pos);

   const bool pInBandRegion = std::abs(relP[2]) < halfLHost &&
                               (relP[0] * relP[0] + relP[1] * relP[1]) < host.radius * host.radius;
   if (!pInBandRegion) return result; // no contact from this host/probe pairing

   const real_t deltaAp = halfLHost - std::abs(relP[2]);
   const real_t rhoP    = std::sqrt(relP[0] * relP[0] + relP[1] * relP[1]);
   const real_t deltaRp = host.radius - rhoP;

   // Guo et al. Table 1, "Special contact" (their Fig. 8): NOT mentioned in Kodam et al.'s original
   // paper at all. Fires when P is on the SAME side of host's face as probe's own centre, but probe's
   // centre lies BEYOND that face (an axial pile-up band-edge would otherwise wrongly report with a
   // too-large radial overlap) -- treated as an axial-overlap contact with a face-normal-like normal.
   const bool sameSide            = relP[2] * relProbeCentre[2] > real_t(0);
   const bool probeCentreBeyondFace = std::abs(relProbeCentre[2]) > halfLHost;
   if (sameSide && probeCentreBeyondFace && deltaAp < deltaRp)
   {
      result.contact  = true;
      result.overlap  = deltaAp;
      result.normal   = host.R * Vec3(real_t(0), real_t(0), math::sign(relP[2]));
      result.point    = P - result.normal * real_t(0.5) * deltaAp;
      result.scenario = "special-contact";
      return result;
   }

   // Band-edge (eq. 79-81), CORRECTED per Guo et al.: Kodam's original delta_r < delta_a requirement
   // is dropped entirely (their Fig. 7 shows a genuine band-edge contact where it doesn't hold); the
   // "NOT face-edge" and "NOT special-contact" exclusions above (via the early returns) take its place.
   result.contact  = true;
   result.overlap  = deltaRp;
   const Vec3 nHat = host.R * ((rhoP > real_t(1e-10)) ? Vec3(relP[0], relP[1], real_t(0)) / rhoP
                                                       : Vec3(real_t(1), real_t(0), real_t(0)));
   result.normal   = nHat;
   result.point    = P + nHat * real_t(0.5) * deltaRp;
   result.scenario = "band-edge";
   return result;
}

//! Cylinder-cylinder contact: dispatches across all 7 scenarios in Fig. 9 (face-face, face-band,
//! face-edge, band-band parallel, band-band skewed, band-edge, edge-edge), in the order the paper
//! presents them. `i`/`j` naming mirrors the paper; `ci` is the 1st argument, `cj` the 2nd.
inline ContactResult cylinderCylinderContact(const CylinderState& ci, const CylinderState& cj)
{
   const Vec3   eI = ci.axis(), eJ = cj.axis();
   const Vec3   relI  = ci.R.getTranspose() * (cj.pos - ci.pos); // cylinder j's centre in i's body frame (eq. 23-25)
   const real_t zIJ   = relI[2];
   const real_t rhoIJ = std::sqrt(relI[0] * relI[0] + relI[1] * relI[1]);
   const real_t halfLi = real_t(0.5) * ci.length, halfLj = real_t(0.5) * cj.length;
   const real_t parallelDot = eI * eJ;

   ContactResult result;

   // --- Face-face / band-band (parallel) precondition (eq. 22) ----------------------------------
   if (std::abs(std::abs(parallelDot) - real_t(1)) < real_t(1e-8) && std::abs(zIJ) < (halfLi + halfLj) &&
       rhoIJ < (ci.radius + cj.radius))
   {
      const real_t deltaAxial  = (halfLi + halfLj) - std::abs(zIJ);        // eq. 26
      const real_t deltaRadial = (ci.radius + cj.radius) - rhoIJ;          // eq. 27 == eq. 53

      if (deltaAxial < deltaRadial)
      {
         // face-face (Section 3.2.1)
         result.contact  = true;
         result.overlap  = deltaAxial;                                     // eq. 28
         // eq. 29 defines this pointing i->j; negated to match this codebase's convention (normal
         // points from the 2nd argument toward the 1st, confirmed against production via
         // testSphereSphereContact et al. -- see the band-band-skewed branch below, which already
         // negates for the same reason).
         result.normal   = ci.R * Vec3(real_t(0), real_t(0), -math::sign(zIJ));

         result.point    = ci.pos + ci.R * faceFaceStyleContactPoint(relI, ci.radius, cj.length,
                                                                      deltaAxial, deltaRadial);
         result.scenario = "face-face";
      }
      else
      {
         // band-band, parallel axes (Section 3.2.4.1)
         result.contact  = true;
         result.overlap  = deltaRadial;                                    // eq. 53
         const Vec3 rHatWorld = ci.R * ((rhoIJ > real_t(1e-10)) ? Vec3(relI[0], relI[1], real_t(0)) / rhoIJ
                                                                 : Vec3(real_t(1), real_t(0), real_t(0)));
         result.normal   = -rHatWorld; // eq. 54 defines this i->j; negated, same reason as eq. 29 above
         result.point    = ci.pos + ci.R * faceFaceStyleContactPoint(relI, ci.radius, cj.length,
                                                                      deltaAxial, deltaRadial);
         result.scenario = "band-band-parallel";
      }
      return result;
   }

   // --- Face-band (Section 3.2.2, eq. 34-44): axes perpendicular ---------------------------------
   if (std::abs(parallelDot) < real_t(1e-8))
   {
      // eq. 35: S = |x_cyl_i^j . e_cyl_j,z|, i.e. projection of (ci.pos - cj.pos) onto cylinder j's
      // axis -- computed below as `Sj` directly from ci's position expressed in cylinder j's body
      // frame (equivalent to the paper's dot-product form, just taken via the frame transform
      // already needed for the rest of this scenario).
      const Vec3   relJ = cj.R.getTranspose() * (ci.pos - cj.pos);
      const real_t Sj   = std::abs(relJ[2]);
      const real_t t    = std::sqrt(relI[0] * relI[0] + relI[1] * relI[1] - Sj * Sj); // eq. 36 (guarded below)

      if (std::abs(zIJ) < (halfLi + cj.radius) && Sj < (halfLj + ci.radius) &&
          (relI[0] * relI[0] + relI[1] * relI[1]) >= Sj * Sj) // domain guard for the sqrt in eq. 36
      {
         const real_t tStar = (Sj < halfLj) ? ci.radius
                                             : std::sqrt(std::max(real_t(0), ci.radius * ci.radius -
                                                                                 (Sj - halfLj) * (Sj - halfLj)));
         if (t < tStar)
         {
            result.contact  = true;
            result.overlap  = halfLi + cj.radius - std::abs(zIJ);          // eq. 38
            // eq. 39, like eq. 29/54, is worded "points ... toward the cylinder's centre" -- i.e.
            // toward j (cj is the cylinder approaching ci's face here), so negated for the same
            // reason as face-face/band-band-parallel above.
            result.normal   = ci.R * Vec3(real_t(0), real_t(0), -math::sign(zIJ));
            // Contact point: the paper places it midway along line AB (cylinder j's band-edge
            // line, eq. 40-41), clamped into cylinder i's face -- distinguishing whether AB lies
            // fully (eq. 42-43) or only partially (eq. 44) inside the face. The "partially inside"
            // sub-case needs a line-circle intersection the paper explicitly omits ("the solution
            // is straightforward and will not be presented here"). Rather than silently guess at
            // that omitted derivation, both sub-cases are approximated here by the same
            // radial/axial split formula used for face-face and band-band-parallel
            // (faceFaceStyleContactPoint) -- exact for the "fully inside" case, an approximation
            // for "partially inside" (flagged, not hidden).
            result.point    = ci.pos + ci.R * faceFaceStyleContactPoint(relI, ci.radius, cj.length,
                                                                         result.overlap, real_t(2) * cj.radius);
            result.scenario = "face-band";
            return result;
         }
      }
   }

   // --- Band-band, skewed axes (Section 3.2.4.2, eq. 55-66) -------------------------------------
   {
      const Vec3 rHat = (eI % eJ);
      const real_t rHatLen = rHat.length();
      if (rHatLen > real_t(1e-6)) // non-parallel (parallel already handled above)
      {
         const Vec3   rHatN = rHat / rHatLen;
         // eq. 56. NOTE: must use the WORLD-frame relative position here, not the body-frame `relI`
         // computed above for the parallel-axis branches -- rHatN (built from ci.axis()/cj.axis(),
         // both world-frame) dotted against a body-frame vector is a frame mismatch that silently
         // produces a wrong (generally too-small) `d` whenever ci itself is rotated away from
         // identity. Every hand-picked test in this file only ever rotates cj, never ci, so this
         // stayed invisible until PaperAnalyticalContactTest.cpp's randomized fuzz test (which
         // rotates both bodies) found it as a source of false-positive band-band-skewed contacts.
         const real_t d = std::abs((cj.pos - ci.pos) * rHatN);

         // Closest points on the two (infinite) centrelines. The paper derives this via Vega & Lago
         // (1994)'s eq. 59-60, expressed in cylinder i's body frame. Implemented here instead as the
         // standard closest-point-between-two-lines formula directly in WORLD vectors (mathematically
         // the same computation Vega & Lago's l_i*/l_j* solve for; this form avoids an extra
         // frame transform that would otherwise need independent double-checking).
         {
            const Vec3   w0     = ci.pos - cj.pos;
            const real_t a_     = eI * eI, b_ = eI * eJ, c_ = eJ * eJ;
            const real_t d_     = eI * w0, e_ = eJ * w0;
            const real_t denom2 = a_ * c_ - b_ * b_;
            real_t       liW = real_t(0), ljW = real_t(0);
            if (std::abs(denom2) > real_t(1e-10))
            {
               liW = (b_ * e_ - c_ * d_) / denom2;
               ljW = (a_ * e_ - b_ * d_) / denom2;
            }

            // Guo, Wassgren, Ketterhagen, Hancock & Curtis (2012), Powder Technology 228, 193-198,
            // Table 1 "Band-band contact (skewed)": Kodam et al.'s original criterion (the first four
            // conditions below) can fire on a near-face-face contact perturbed slightly off-parallel,
            // reporting an unphysically large overlap (their Fig. 9). Their fix adds the fifth
            // condition (eq. 6), excluding cases where the reported overlap would exceed how much
            // room is actually left on EITHER cylinder's axial extent -- such cases correctly fall
            // through to face-edge or the "special contact" scenario below instead.
            const bool guoExclusion =
               (ci.radius + cj.radius - d) <
               std::min(halfLi - std::abs(liW), halfLj - std::abs(ljW));

            if (std::abs(liW) < halfLi && std::abs(ljW) < halfLj && d < (ci.radius + cj.radius) &&
                guoExclusion)
            {
               result.contact  = true;
               result.overlap  = (ci.radius + cj.radius) - d;              // eq. 64
               const Vec3 Pi = ci.pos + liW * eI;
               const Vec3 Pj = cj.pos + ljW * eJ;
               result.point    = real_t(0.5) * (Pi + Pj);                  // eq. 65
               const Vec3 diff = Pj - Pi;
               result.normal   = (diff.length() > real_t(1e-10)) ? -diff.getNormalized() : rHatN; // eq. 66, i<-j flipped to match this codebase's idx2->idx1 convention
               result.scenario = "band-band-skewed";
               return result;
            }
         }
      }
   }

   // --- Face-edge / special-contact / band-edge: capped surfaces facing each other ---------------
   // Checked in BOTH directions -- j's rim against i's face/band, then i's rim against j's face/band
   // ("if the edge of cylinder i contacts with the face of cylinder j, the same algorithm may be
   // used, but with superscripts i and j interchanged" -- Kodam et al., end of Section 3.2.3) -- using
   // rimVsFaceOrBand's Guo-et-al.-corrected criteria. Only entered if none of the simpler scenarios
   // above matched. This replaces an earlier, uncorrected version of this section that
   // PaperAnalyticalContactTest.cpp demonstrated produced false negatives for oblique tilts (see the
   // git history / that test file's comments for the specific failing geometries); Guo et al.'s
   // corrections target exactly this class of scenario-boundary misclassification, though whether
   // they close EVERY oblique case those tests constructed, versus only the near-degenerate
   // perturbation cases Guo et al.'s own figures diagram, is an empirical question the retested suite
   // answers below rather than one this comment can promise in advance.
   {
      // j's rim against i's face/band: rimVsFaceOrBand returns i->j (host->probe); negate for this
      // file's toward-1st-argument convention (see the face-face branch above).
      ContactResult jRimVsI = rimVsFaceOrBand(ci, cj);
      if (jRimVsI.contact)
      {
         jRimVsI.normal = -jRimVsI.normal;
         return jRimVsI;
      }

      // i's rim against j's face/band: host=cj, probe=ci, so rimVsFaceOrBand's natural host->probe
      // direction (j->i) IS ALREADY this file's toward-1st-argument (i) convention -- no negation.
      const ContactResult iRimVsJ = rimVsFaceOrBand(cj, ci);
      if (iRimVsJ.contact) { return iRimVsJ; }
   }

   // --- Edge-edge (Section 3.2.6, eq. 82-101): both rims (circular edges) facing each other. ------
   //
   // An earlier version of this branch computed the closest-point PAIR between the two FULL rims
   // via alternating projection, gated purely by "is that distance small" -- a proximity test, not
   // Kodam's actual existence test, and it reported `overlap = -dist` with `dist` a Euclidean norm
   // (always >= 0), which is structurally incompatible with this file's "positive = penetrating"
   // convention no matter how deep the true interpenetration was. Fuzz-tested against GJK/EPA that
   // version had a 657/10000 false-positive rate plus 310/10000 magnitude mismatches, 100% traced to
   // this one branch.
   //
   // This version implements Kodam's actual eq. 82-101: findRimCrossing finds points A, B on rim j
   // (where it crosses cylinder i's lateral surface and near face respectively, eq. 83-86) -- their
   // combined success IS Kodam's own existence test, not a proximity threshold -- and symmetrically
   // C, D on rim i vs cylinder j. Given both chords exist, eq. 87-98's exact skew-line closest
   // distance between them replaces the old alternating-projection approximation, and since contact
   // is now established by genuine geometric containment rather than "rims happen to be close," the
   // resulting distance is used directly as a small positive overlap (matching this file's
   // convention) rather than negated -- unlike the old branch, existence and a valid sign are no
   // longer decoupled. Re-tested against the same 10000-trial fuzz harness: see
   // PaperAnalyticalContactTest.cpp for the current pass rate; cylinderCylinderResultIsUncertain's
   // Tier 1/Tier 2 gates remain in place regardless as the safety net for whatever this doesn't
   // resolve (near-degenerate configurations, multiple-root ambiguity Kodam's own paper doesn't
   // fully disambiguate for the general skewed case) -- GJK/EPA fallback is unconditionally
   // available, so an imperfect analytical answer here costs an extra fallback call, not a wrong one.
   {
      Vec3 uJ, vJ, uI, vI;
      detail::perpendicularBasis(eJ, uJ, vJ);
      detail::perpendicularBasis(eI, uI, vI);
      const Vec3 capJ = cj.pos + math::sign((ci.pos - cj.pos) * eJ) * halfLj * eJ;
      const Vec3 capI = ci.pos + math::sign((cj.pos - ci.pos) * eI) * halfLi * eI;

      const detail::RimCrossing jVsI = detail::findRimCrossing(capJ, uJ, vJ, cj.radius, ci.pos, eI, ci.radius, halfLi);
      if (jVsI.contact)
      {
         const detail::RimCrossing iVsJ =
            detail::findRimCrossing(capI, uI, vI, ci.radius, cj.pos, eJ, cj.radius, halfLj);
         if (iVsJ.contact)
         {
            // eq. 87-92: lines through chord AB (rim j) and chord CD (rim i).
            const Vec3   rJ = real_t(0.5) * (jVsI.A + jVsI.B);
            const Vec3   rI = real_t(0.5) * (iVsJ.A + iVsJ.B);
            const Vec3   dJ = jVsI.B - jVsI.A;
            const Vec3   dI = iVsJ.B - iVsJ.A;
            const real_t lenJ = dJ.length(), lenI = dI.length();

            if (lenJ > real_t(1e-12) && lenI > real_t(1e-12))
            {
               const Vec3   eJHat = dJ / lenJ;
               const Vec3   eIHat = dI / lenI;
               const real_t dot   = eIHat * eJHat;
               const real_t denom = real_t(1) - dot * dot;

               if (std::abs(denom) > real_t(1e-9)) // eq. 93-96, guarded against near-parallel chords
               {
                  const Vec3   rij    = rI - rJ;
                  const real_t muStar = -((rij * eIHat) - (rij * eJHat) * dot) / denom;
                  const real_t lamStar = ((rij * eJHat) - (rij * eIHat) * dot) / denom;

                  const Vec3   Pj   = rJ + lamStar * eJHat; // eq. 97
                  const Vec3   Pi   = rI + muStar * eIHat;  // eq. 98
                  const real_t dist = (Pi - Pj).length();   // eq. 93/99: delta_contact = d_perp

                  result.contact  = true;
                  result.overlap  = dist; // contact already established geometrically above (unlike
                                           // the old proximity-only test), so d_perp here is a
                                           // genuine (small, near-first-contact) penetration estimate
                  result.normal   = (Pi - Pj).getNormalizedOrZero();
                  result.point    = real_t(0.5) * (Pi + Pj); // eq. 101
                  result.scenario = "edge-edge";
               }
            }
         }
      }
   }

   return result;
}

//! Conservative "should this analytical result be trusted, or should the caller fall back to a
//! general algorithm (GJK/EPA) instead?" check. Deliberately biased toward over-triggering: a false
//! "needs fallback" only costs an extra GJK/EPA call; a false "safe to trust" silently returns a
//! wrong answer, which is what this function exists to prevent. See PaperAnalyticalContactTest.cpp
//! for the specific geometries this was calibrated and tested against (cylinderCylinderContactHybrid
//! and its accompanying tests).
//!
//! Two independent triggers:
//!
//! TIER 1 (dispatch found nothing, but the bodies are close enough that a real contact is
//! plausible): a cheap bounding-sum gate -- centre-to-centre distance vs. the sum of both cylinders'
//! "reach" (radius + half-length each, i.e. the distance at which even a fully-extended,
//! perfectly-aimed pair could just barely touch tip-to-tip). This is the SAME style of coarse
//! bounding check DEM codes already use before ever calling fine contact detection (cf. Guo et al.
//! 2012 Section 2.1's circumscribing-sphere/cube coarse detection) -- deliberately loose, so it will
//! flag plenty of genuinely-separated pairs too; that's fine here since triggering it only means
//! "double-check with GJK/EPA", not "report a contact".
//!
//! TIER 2 (dispatch found a contact of a scenario type demonstrated -- by PaperAnalyticalContactTest
//! .cpp's randomized fuzz test, testCylinderCylinderRandomizedFuzz -- to sometimes report a WRONG
//! answer, including reporting contact where GJK/EPA finds none at all): cross-check the reported
//! overlap against each cylinder's CIRCUMSCRIBING SPHERE overlapping the other's face plane. A
//! cylinder of radius r and half-length h fits entirely inside a sphere of radius sqrt(r^2+h^2)
//! centred at its own centre -- ALWAYS true regardless of orientation, so "circumscribing-sphere
//! overlap of the other cylinder's face plane" is a valid (if loose) UPPER BOUND on how much that
//! cylinder could possibly be penetrating, independent of tilt angle. If this bound exceeds what the
//! scenario actually reported, more penetration than the analytical answer accounts for is
//! geometrically possible -- so the reported number (or its very existence) should not be trusted.
//!
//! Originally this only covered face-edge/special-contact/band-edge (the three scenarios where a
//! WRONG MAGNITUDE was already known, from hand-picked test geometries, because the true closest
//! feature can be on a curved band neither Kodam et al. nor Guo et al.'s corrections ever search).
//! band-band-skewed and edge-edge were excluded, on the belief they were "fully verified against
//! GJK/EPA with no open discrepancy" -- the fuzz test disproved that directly: with separation
//! sampled near the true GJK/EPA contact boundary (5000 random cylinder pairs, realistic overlap
//! depth), band-band-skewed and edge-edge produced ~8% FALSE POSITIVES -- reporting contact for
//! pairs GJK/EPA confirms are genuinely separated -- every single one silently trusted by this
//! function exactly because they were excluded from Tier 2. They are no longer excluded.
//!
//! face-face and band-band-parallel remain excluded: those two apply only under the near-exact-
//! parallel-axis precondition (eq. 22, checked to 1e-8), where the closed form is a direct algebraic
//! solution with no curved-band/rim search involved at all -- the fuzz test's own results confirm
//! they're never the source of a mismatch (all logged mismatches are band-band-skewed or edge-edge).
inline bool cylinderCylinderResultIsUncertain(const CylinderState& ci, const CylinderState& cj,
                                               const ContactResult& analytical)
{
   if (!analytical.contact)
   {
      const real_t centreDist = (cj.pos - ci.pos).length();
      const real_t reach = ci.radius + cj.radius + real_t(0.5) * ci.length + real_t(0.5) * cj.length;
      return centreDist < reach; // Tier 1
   }

   const bool isFullyVerifiedScenario =
      std::string(analytical.scenario) == "face-face" || std::string(analytical.scenario) == "band-band-parallel";
   if (isFullyVerifiedScenario) return false;

   auto circumscribingSphereOverlapOfOtherFace = [](const CylinderState& probe, const CylinderState& host) {
      const real_t circumRadius = std::sqrt(probe.radius * probe.radius +
                                             real_t(0.25) * probe.length * probe.length);
      // Use whichever of host's two face planes probe's centre is nearer to -- the relevant one.
      const Vec3   eHost  = host.axis();
      const real_t halfLh = real_t(0.5) * host.length;
      const Vec3   topFace = host.pos + halfLh * eHost;
      const Vec3   botFace = host.pos - halfLh * eHost;
      const real_t zTop = std::abs((probe.pos - topFace) * eHost);
      const real_t zBot = std::abs((probe.pos - botFace) * eHost);
      const real_t zFace = std::min(zTop, zBot);
      return circumRadius - zFace; // Tier 2, per-direction
   };

   const real_t maxPossibleJ = circumscribingSphereOverlapOfOtherFace(cj, ci); // cj's reach into ci's face
   const real_t maxPossibleI = circumscribingSphereOverlapOfOtherFace(ci, cj); // ci's reach into cj's face

   return (maxPossibleJ > analytical.overlap) || (maxPossibleI > analytical.overlap); // Tier 2
}

} // namespace paper_analytical
} // namespace antidunes
} // namespace walberla
