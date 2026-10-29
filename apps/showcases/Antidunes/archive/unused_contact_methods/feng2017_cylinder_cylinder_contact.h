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
//! \file feng2017_cylinder_cylinder_contact.h
//! \brief ARCHIVED, unused alternative cylinder-cylinder contact method (Feng, Han & Owen 2017).
//! Zero callers anywhere in this showcase -- never wired into production or any DEMValidation
//! scenario. Moved out of PaperAnalyticalContacts.h 2026-10-29; see archive/ARCHIVE_README.md.
#pragma once

#include "../../PaperAnalyticalContacts.h"

namespace walberla
{
namespace antidunes
{
//======================================================================================================================
// Cylinder x Cylinder, ALTERNATIVE method -- Feng, Y.T., Han, K., & Owen, D.R.J. (2017). "A generic
// contact detection framework for cylindrical particles in discrete element modelling." Computer
// Methods in Applied Mechanics and Engineering, 315, 632-651.
//
// Unlike paper_analytical::cylinderCylinderContact (Kodam et al. 2010 / Guo et al. 2012), which
// enumerates 7 named scenarios and searches specific rim/edge candidate points for each, this method
// never enumerates scenarios at all: it builds a 4-parameter local frame from the shortest line
// segment between the two cylinder AXES (Sec. 2), then reduces "does cylinder j's rim get within
// radius of cylinder i's axis" to a closed-form closest-point-on-an-ellipse-to-the-origin problem
// (Sec. 3.1, Appendix A) -- the paper's own motivation for this is explicitly the same rim-vs-band
// gap this codebase's own testing independently found in Kodam/Guo: "the computations of the contact
// overlap and contact normal for most of the identified cases are often of an ad-hoc nature... the
// inappropriate choice of a contact normal may introduce superficial energy into the system."
//
// SCOPE, deliberately limited relative to the full paper:
//  - Implements Sec. 3.1's non-overlap check (Step 1 interior contact, Step 3 base-disc-vs-axis
//    quartic check) as a genuine existence + approximate overlap test. Step 2 (the profile-rectangle
//    early-reject) is skipped -- the paper itself positions it purely as a speed optimisation ahead
//    of Step 3, not something Step 3's own correctness depends on.
//  - Does NOT implement Sec. 3.3/3.4's full layered/edge overlap-REGION representation. Instead, of
//    the (up to 4) base-disc-vs-axis candidates Step 3 considers, whichever has the largest
//    penetration is used directly for overlap/point -- an approximation, not the paper's own
//    region-based quantity, flagged here exactly as this file flags its other approximations.
//  - The contact normal uses Theorem 1's own PHYSICAL statement (the normal must pass through both
//    cylinder axes) re-derived directly as the cross product of the two axis-plane normals, rather
//    than transcribing the paper's eq. (23): that equation's two claimed-equal forms use n_{p1} in
//    one and an apparently-mistyped "n_p" (missing its "2" subscript, or a rendering artefact) in the
//    other, and independently deriving eq. (33) from the paper's OWN clearly-stated eq. (29)/(31)
//    produces the opposite sign from what eq. (33) itself prints. Rather than propagate whichever
//    side of that inconsistency is the actual typo, the normal below is derived fresh from the
//    unambiguous physical statement instead (see fengContactNormal's own comment).
//  - Exact axis-parallel configurations (cos(alpha) ~ 1) are explicitly out of scope, matching the
//    paper's own Sec. 2 scoping ("cos alpha = 1... is a special case that will be considered
//    [separately]") -- paper_analytical's face-face/band-band-parallel branches already handle that
//    case exactly, or GJK/EPA.
//======================================================================================================================
namespace feng2017
{
using mesa_pd::Mat3;
using mesa_pd::Vec3;
using paper_analytical::ContactResult;
using paper_analytical::CylinderState;
using paper_analytical::detail::closestThetaOnEllipseToOrigin;
using paper_analytical::detail::closestThetaOnEllipseWithinAxialRange;
using paper_analytical::detail::intersectEllipseWithCircle;

//! Theorem 1 (Contact Normal Theorem), re-derived from the paper's own physical statement rather
//! than its eq. (23) (see the namespace-level comment for why): the normal line through contact
//! point `p` must lie in the plane containing `p` and cylinder 1's axis (normal to that plane:
//! np1 = n1 x (c1-p), exactly the paper's own eq. (30)/(36) definition) AND in the analogous plane
//! for cylinder 2's axis (normal np2 = n2 x (c2-p)). A direction lying in both planes must be
//! perpendicular to both planes' normals, i.e. n = np1 x np2 -- symmetric in the two cylinders,
//! unlike eq. (23)'s asymmetric t-parametrised construction, and needing no separate "prove the
//! first part similarly" step the paper's own proof leaves to the reader.
inline Vec3 fengContactNormal(const Vec3& p, const Vec3& n1, const Vec3& c1, const Vec3& n2, const Vec3& c2)
{
   const Vec3 np1 = n1 % (c1 - p);
   const Vec3 np2 = n2 % (c2 - p);
   Vec3       n   = (np1 % np2).getNormalizedIfNotZero(); // zero iff p, axis 1 and axis 2 are
                                                            // co-planar -- Theorem 1's own documented
                                                            // exceptional case
   // n x (-n) are both valid solutions of "perpendicular to both plane normals" -- orient toward
   // cylinder 1, matching this file's own "points toward the 1st argument" convention (established
   // for every paper_analytical scenario and sphere-sphere alike).
   if (n * (c1 - c2) < real_t(0)) n = -n;
   return n;
}

//! `debugClosestMissGap`, when non-null, is set to the smallest (distToAxis - hostRadius) found
//! across every Step-3 candidate that had a genuinely axial-range-valid closest point (i.e.
//! closestThetaOnEllipseWithinAxialRange succeeded), or left at +infinity if NONE of the 4
//! candidates ever had a valid one at all. Diagnostic only, purely to distinguish "missed by a
//! hair" (a small positive gap: a legitimately-in-range candidate came close but not close enough)
//! from "no candidate ever qualified" (a structural gap -- +infinity) when this disagrees with an
//! independent reference. Does not affect `result`.
inline ContactResult fengCylinderCylinderContact(const CylinderState& c1, const CylinderState& c2,
                                                  real_t* debugClosestMissGap = nullptr)
{
   ContactResult result;
   real_t        closestMissGap = std::numeric_limits< real_t >::max();

   const Vec3 n1 = c1.axis();
   Vec3       n2 = c2.axis();
   real_t     cosAlpha = n1 * n2;
   if (cosAlpha < real_t(0))
   {
      n2 = -n2;
      cosAlpha = -cosAlpha;
   }
   if (cosAlpha > real_t(1) - real_t(1e-6))
   {
      result.scenario = "feng-parallel-unsupported"; // out of scope, see namespace comment
      return result;
   }

   const real_t sinAlpha = std::sqrt(std::max(real_t(0), real_t(1) - cosAlpha * cosAlpha));
   const Vec3   ex  = (n2 % n1).getNormalized();  // eq. (2): common-perpendicular direction
   const Vec3   ey1 = n1 % ex;                    // eq. (3)/(6): cylinder 1's own local y-axis
   const Vec3   ey2 = n2 % ex;                    // cylinder 2's own local y-axis

   const Vec3   c2c1 = c2.pos - c1.pos;
   const real_t d  = ex * c2c1;                   // eq. (4): signed shortest distance between axes
   const real_t t1 = n1 * c2c1;
   const real_t t2 = n2 * c2c1;
   const real_t s1 = (t1 - t2 * cosAlpha) / (sinAlpha * sinAlpha); // eq. (7)
   const real_t s2 = -(t2 + t1 * cosAlpha) / (sinAlpha * sinAlpha);

   const real_t h1 = real_t(0.5) * c1.length, h2 = real_t(0.5) * c2.length;

   // --- Step 1: interior contact (Sec. 3.1) -- the closest-approach points o1, o2 between the two
   // infinite axes both fall within the cylinders' own finite half-lengths. -----------------------
   if (std::abs(d) < c1.radius + c2.radius && std::abs(s1) < h1 && std::abs(s2) < h2)
   {
      const Vec3 o1 = c1.pos + s1 * n1;
      const Vec3 o2 = c2.pos + s2 * n2;
      result.contact  = true;
      result.overlap  = (c1.radius + c2.radius) - std::abs(d);
      result.point    = real_t(0.5) * (o1 + o2);
      result.normal   = fengContactNormal(result.point, n1, c1.pos, n2, c2.pos);
      if (result.normal.sqrLength() < real_t(1e-12)) result.normal = math::sign(-d) * ex; // co-planar
                                                                                            // fallback
      result.scenario = "feng-interior";
      return result;
   }
   if (std::abs(d) >= c1.radius + c2.radius) return result; // axes too far apart, definitely separated

   // --- Step 3: base-disc-vs-axis quartic check (Sec. 3.1, Appendix A) -----------------------------
   // Up to 4 candidates: each of the two cylinders' two base-disc rims, projected as an ellipse into
   // the OTHER cylinder's cross-sectional plane, checked against that other cylinder's radius. See
   // this function's own derivation (namespace comment) for why the ellipse's semi-axes are exactly
   // (probeRadius, probeRadius*cosAlpha) and why the ellipse parameter equals the true 3D rim angle
   // directly (ex, ey_host, and ey_probe share a common perpendicular decomposition).
   struct Candidate
   {
      bool        valid   = false;
      real_t      overlap = real_t(0);
      Vec3        point, normal;
      const char* kind = "feng-rim-edge"; // vs. "feng-rim-face" -- see the two tryDisc branches below;
                                            // kept distinct so their reliability can be measured
                                            // separately (they have different failure profiles).
   };
   Candidate best;

   auto tryDisc = [&](const Vec3& hostPos, const Vec3& hostN, const Vec3& hostEx, const Vec3& hostEy,
                       real_t hostRadius, real_t hostHalfLen, const Vec3& probeEy, const Vec3& baseCentre,
                       real_t probeRadius) {
      const real_t a = probeRadius, b = probeRadius * cosAlpha;
      const real_t x0 = hostEx * (baseCentre - hostPos);
      const real_t y0 = hostEy * (baseCentre - hostPos);
      const real_t zC = hostN * (baseCentre - hostPos);

      // Face-type candidate (Kodam's "face-band" mechanism: probe's flat disc INTERIOR, not its rim,
      // is the closest feature) -- third attempt, using Appendix B to detect the case properly rather
      // than inferring it from Appendix A alone.
      //   Attempt 1: "origin inside the ellipse", disc-centre axial reference (zC). Net regression:
      //     309->467 total mismatches (245 new false positives).
      //   Attempt 2: added the "closest boundary point clears hostRadius" guard, and corrected the
      //     axial reference to zHostAtOrigin (the disc PLANE's crossing point at host's axis, not the
      //     disc's own centre). Net regression: 309->339 total mismatches (53 new false positives).
      // Both attempts used a SINGLE axial reference point, when the disc is a tilted plane whose
      // axial position along host varies continuously across it. This attempt instead uses
      // intersectEllipseWithCircle (Appendix B) as the authoritative containment DETECTOR (zero real
      // roots at radius hostRadius means the boundary genuinely never crosses -- either fully inside
      // or fully outside, with no ambiguity left to a separate closed-form guess), and computes the
      // full RANGE of axial positions the disc's plane spans across host's entire cross-section
      // (not one point) for the overlap magnitude.
      {
         real_t crossings[4];
         int    crossingCount;
         intersectEllipseWithCircle(a, b, x0, y0, hostRadius, crossings, crossingCount);
         const real_t insideParam = (a > real_t(1e-12) && b > real_t(1e-12))
                                        ? (x0 * x0) / (a * a) + (y0 * y0) / (b * b)
                                        : real_t(2);
         if (crossingCount == 0 && insideParam < real_t(1))
         {
            // No crossing anywhere on the rim, AND host's own axis point lies inside the ellipse ->
            // the rim never dips inside hostRadius (checked next) means host's ENTIRE cross-section,
            // not just its axis, is swallowed by the disc.
            const real_t thetaU          = closestThetaOnEllipseToOrigin(a, b, x0, y0);
            const real_t xU              = x0 + a * std::cos(thetaU);
            const real_t yU              = y0 + b * std::sin(thetaU);
            const real_t boundaryMinDist = std::sqrt(xU * xU + yU * yU);
            if (boundaryMinDist > hostRadius)
            {
               // z1(x1,y1) on the disc's plane is LINEAR in (x1,y1) (see the derivation in the
               // conversation this was worked out in): z1 = zHostAtOrigin + x1*uA + y1*(uB/cosAlpha).
               // Its extreme values over host's entire cross-sectional disk (radius hostRadius,
               // centred at the origin) occur at the disk boundary in the gradient direction:
               // zHostAtOrigin +- hostRadius*|gradient|. Using the RANGE this spans (not a single
               // point) is the correction over both earlier attempts.
               const real_t uA = hostN * hostEx, uB = hostN * probeEy;
               const real_t sAtOrigin = -x0, tAtOrigin = -y0 / cosAlpha;
               const real_t zHostAtOrigin = zC + sAtOrigin * uA + tAtOrigin * uB;
               const real_t gradMag = std::sqrt(uA * uA + (uB / cosAlpha) * (uB / cosAlpha));
               const real_t zMin = zHostAtOrigin - hostRadius * gradMag;
               const real_t zMax = zHostAtOrigin + hostRadius * gradMag;
               // Point in [zMin, zMax] closest to host's own centre (z1=0) -- the most-overlapping
               // slice of the disc's coverage relative to host's finite extent.
               const real_t zBest = std::clamp(real_t(0), zMin, zMax);

               const real_t faceOverlap = hostHalfLen - std::abs(zBest);
               if (faceOverlap > real_t(0) && (!best.valid || faceOverlap > best.overlap))
               {
                  const real_t zSigned   = math::sign(zBest) * (hostHalfLen - real_t(0.5) * faceOverlap);
                  const Vec3   facePoint = hostPos + zSigned * hostN;
                  best.valid   = true;
                  best.overlap = faceOverlap;
                  best.point   = facePoint;
                  best.kind    = "feng-rim-face";
                  best.normal  = fengContactNormal(facePoint, n1, c1.pos, n2, c2.pos);
                  if (best.normal.sqrLength() < real_t(1e-12))
                     best.normal = (-math::sign(zBest) * hostN).getNormalizedIfNotZero();
               }
            }
         }
      }

      // zHost(theta) = hostN . (rimPoint(theta) - hostPos): the SAME rim point's axial position along
      // the host's own axis, as a function of the SAME theta -- also just a first-harmonic function,
      // since rimPoint(theta) is itself a first-harmonic function of theta. Constraining the ellipse
      // search to where |zHost(theta)| < hostHalfLen (rather than solving unconstrained and rejecting
      // the answer if it lands outside that range) is the fix for the false-negative gap fuzz-testing
      // found in the unconstrained version -- see closestThetaOnEllipseWithinAxialRange's own comment.
      const real_t zA = probeRadius * (hostN * hostEx);
      const real_t zB = probeRadius * (hostN * probeEy);

      real_t theta;
      if (!closestThetaOnEllipseWithinAxialRange(a, b, x0, y0, zC, zA, zB, hostHalfLen, theta))
         return; // the ENTIRE rim lies outside the host's finite length -- genuinely no candidate at
                  // all here (NOT tracked in closestMissGap: the unconstrained minimum's own gap is
                  // meaningless here precisely because it's known to lie outside the valid range)

      // theta is simultaneously the ellipse parameter AND the true 3D angle around the probe's own
      // base-disc rim (both bases share the same (ex, ey_probe) orthonormal in-plane frame), so the
      // actual 3D rim point is reconstructed directly, with no separate inverse mapping needed.
      const Vec3   rimPoint = baseCentre + probeRadius * (std::cos(theta) * hostEx + std::sin(theta) * probeEy);
      const real_t zHost    = hostN * (rimPoint - hostPos);

      const Vec3   radial     = (rimPoint - (hostPos + zHost * hostN));
      const real_t distToAxis = radial.length();
      closestMissGap = std::min(closestMissGap, distToAxis - hostRadius);
      if (distToAxis >= hostRadius) return;

      const real_t overlap = hostRadius - distToAxis;
      if (!best.valid || overlap > best.overlap)
      {
         Vec3 radialHat = radial.getNormalizedIfNotZero(); // points from host axis toward the rim
                                                             // point (i.e. away from host)
         if (radialHat * (c1.pos - c2.pos) < real_t(0)) radialHat = -radialHat; // orient toward c1,
                                                                                  // same convention as
                                                                                  // fengContactNormal
         best.valid   = true;
         best.overlap = overlap;
         best.point   = rimPoint - real_t(0.5) * overlap * radialHat;
         best.normal  = fengContactNormal(best.point, n1, c1.pos, n2, c2.pos);
         if (best.normal.sqrLength() < real_t(1e-12)) best.normal = radialHat; // co-planar fallback
      }
   };

   tryDisc(c1.pos, n1, ex, ey1, c1.radius, h1, ey2, c2.pos + h2 * n2, c2.radius);
   tryDisc(c1.pos, n1, ex, ey1, c1.radius, h1, ey2, c2.pos - h2 * n2, c2.radius);
   tryDisc(c2.pos, n2, ex, ey2, c2.radius, h2, ey1, c1.pos + h1 * n1, c1.radius);
   tryDisc(c2.pos, n2, ex, ey2, c2.radius, h2, ey1, c1.pos - h1 * n1, c1.radius);

   if (best.valid)
   {
      result.contact  = true;
      result.overlap  = best.overlap;
      result.point    = best.point;
      result.normal   = best.normal;
      result.scenario = best.kind; // "feng-rim-edge" or "feng-rim-face" -- see Candidate's own comment
   }
   if (debugClosestMissGap) *debugClosestMissGap = closestMissGap;
   return result;
}

} // namespace feng2017
} // namespace antidunes
} // namespace walberla
