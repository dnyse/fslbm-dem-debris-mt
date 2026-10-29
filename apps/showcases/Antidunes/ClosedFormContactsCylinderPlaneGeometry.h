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
//! \file ClosedFormContactsCylinderPlaneGeometry.h
//! \brief Cylinder x Plane / Cylinder x Sphere closed-form contact: Kodam, Bharadwaj, Curtis, Hancock
//! & Wassgren (2010), Chemical Engineering Science 65(22), 5852-5862 (Part I). CylinderState, the
//! quartic/ellipse/rim-crossing solvers, cylinderPlaneContact, sphereCylinderContact,
//! faceFaceStyleContactPoint. Split out of PaperAnalyticalContacts.h (namespace paper_analytical)
//! 2026-10-29; see ClosedFormContactsInfra.h for the shared ContactResult type.

#pragma once

#include "ClosedFormContactsInfra.h"

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{

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
// maxIterations: default 30 (see below) preserves current behaviour for every one of this function's
// existing call sites (none pass this argument). Exposed as a parameter (not just a local constant)
// specifically so PaperAnalyticalContactTest.cpp can A/B/C-compare 30 against the original 60 directly
// on identical trials -- see testEdgeEdgeFastReferenceAccuracy10k / any 3-way iteration-count
// comparison built on top of it.
inline void solveQuarticDurandKerner(real_t c4, real_t c3, real_t c2, real_t c1, real_t c0,
                                      std::complex< real_t > roots[4], int maxIterations = 30)
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

   // Default 30 (down from the original 60): this file's own comment above already documented that
   // convergence to double precision was empirically observed within ~20-30 iterations, with 60 kept
   // purely as "comfortable margin" -- i.e. the original count was never load-bearing for correctness,
   // just headroom. Halving it (to the upper end of the already-observed-sufficient range, not the
   // aggressive low end) roughly halves the cost of every quartic solve in this file
   // (closestThetaOnCircleToLine, allCriticalThetasOnEllipseToOrigin, and everything built on them:
   // band-edge, special-contact, edge-edge, face-edge). Re-validated against the full
   // PaperAnalyticalContactTest.cpp suite (fuzz tests included) after this change, not just asserted
   // from the original comment's claim -- see also the direct 30-vs-60 A/B comparison enabled by this
   // being a parameter now, not a hardcoded constant.
   for (int it = 0; it < maxIterations; ++it)
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

//! One real root of the depressed-then-un-depressed real cubic t^3 + a*t^2 + b*t + c = 0. A real cubic
//! always has at least one real root, so this never needs to return a complex value -- used below only
//! to pick a resolvent-cubic root for solveQuarticFerrari, which needs exactly one (any one) real root,
//! not all three. Branches between Cardano's formula (one-real-root case, safe/direct) and the
//! trigonometric form (three-real-roots case, avoids the cube root of a near-zero or catastrophically
//! cancelling complex quantity that Cardano's formula would otherwise hit) -- precisely the kind of
//! discriminant-sign branch that solveQuarticDurandKerner's own comment calls out as Ferrari's method's
//! weak point; isolating it to a single real root of a cubic (rather than needing all 4 roots of the
//! quartic to come out of the same branch) is the most that can be done to contain it.
inline real_t solveRealCubicOneRealRoot(real_t a, real_t b, real_t c)
{
   const real_t p      = b - a * a / real_t(3);
   const real_t q      = real_t(2) * a * a * a / real_t(27) - a * b / real_t(3) + c;
   const real_t offset = -a / real_t(3);

   const real_t discriminant = q * q / real_t(4) + p * p * p / real_t(27);
   if (discriminant >= real_t(0))
   {
      const real_t sqrtDisc = std::sqrt(discriminant);
      const real_t u1       = std::cbrt(-q / real_t(2) + sqrtDisc);
      const real_t u2       = std::cbrt(-q / real_t(2) - sqrtDisc);
      return u1 + u2 + offset;
   }
   // Three real roots: p is necessarily < 0 here (discriminant < 0 forces it), so -p/3 > 0 and the
   // sqrt below is safe.
   const real_t r   = std::sqrt(-p * p * p / real_t(27));
   const real_t phi = std::acos(std::clamp(-q / (real_t(2) * r), real_t(-1), real_t(1)));
   const real_t u    = real_t(2) * std::sqrt(-p / real_t(3)) * std::cos(phi / real_t(3));
   return u + offset;
}

//! Real roots of c4*t^4 + c3*t^3 + c2*t^2 + c1*t + c0 = 0 via Ferrari's closed-form method: depress to
//! y^4 + p*y^2 + q*y + r (x = y - c3/(4*c4)), solve the resolvent cubic
//! m^3 + p*m^2 + (p^2/4 - r)*m - q^2/8 = 0 for one real root m (solveRealCubicOneRealRoot), then factor
//! the quartic into two quadratics via that m and solve each in closed form. No iteration at all --
//! O(1) instead of solveQuarticDurandKerner's fixed 30 sweeps.
//! [Claude, 2026-09-29] This IS the production solver: every quartic solve in the cylinder-cylinder
//! dispatch defaults to QuarticSolverMethod::Ferrari (closestThetaOnCircleToLine, findRimCrossing,
//! cylinderCylinderContact). The previous comment here ("not wired into any production call site",
//! "Durand-Kerner in production") was outdated. No case has been found where it fails. Known
//! theoretical weak points, kept here as documentation: (1) the resolvent cubic's discriminant sign is a
//! branch (handled in solveRealCubicOneRealRoot); (2) m = 0 (q = 0, biquadratic quartic) needs the
//! special handling below to avoid dividing by sqrt(2m); (3) like any closed-form quartic formula it can
//! lose accuracy through cancellation when roots nearly coincide, where iterative solvers degrade more
//! gracefully. Durand-Kerner / Aberth-Ehrlich remain available for comparison only.
inline void solveQuarticFerrari(real_t c4, real_t c3, real_t c2, real_t c1, real_t c0,
                                 std::complex< real_t > roots[4])
{
   const real_t b3 = c3 / c4, b2 = c2 / c4, b1 = c1 / c4, b0 = c0 / c4;

   const real_t p = b2 - real_t(3) * b3 * b3 / real_t(8);
   const real_t q = b3 * b3 * b3 / real_t(8) - b3 * b2 / real_t(2) + b1;
   const real_t r = -real_t(3) * b3 * b3 * b3 * b3 / real_t(256) + b3 * b3 * b2 / real_t(16) -
                     b3 * b1 / real_t(4) + b0;
   const real_t shift = -b3 / real_t(4);

   // Biquadratic special case (q == 0, i.e. the depressed quartic has no odd-degree term): solving the
   // resolvent cubic here would hit m=0 (or a spurious nonzero root that still works, but needlessly),
   // since the general construction below divides by q. Solving y^4+p*y^2+r=0 directly as a quadratic
   // in y^2 sidesteps that degeneracy entirely rather than special-casing sqrt2m~=0 further down.
   if (std::abs(q) < real_t(1e-13) * std::max(real_t(1), std::abs(p)))
   {
      const std::complex< real_t > disc =
         std::sqrt(std::complex< real_t >(p * p - real_t(4) * r, real_t(0)));
      const std::complex< real_t > ySq1 = (-p + disc) / real_t(2);
      const std::complex< real_t > ySq2 = (-p - disc) / real_t(2);
      const std::complex< real_t > y1 = std::sqrt(ySq1), y2 = std::sqrt(ySq2);
      roots[0] = y1 + shift;
      roots[1] = -y1 + shift;
      roots[2] = y2 + shift;
      roots[3] = -y2 + shift;
      return;
   }

   const real_t m = solveRealCubicOneRealRoot(p, p * p / real_t(4) - r, -q * q / real_t(8));

   const std::complex< real_t > sqrt2m = std::sqrt(std::complex< real_t >(real_t(2) * m, real_t(0)));
   const std::complex< real_t > qTerm  = q / (real_t(2) * sqrt2m);

   auto solveQuadratic = [](std::complex< real_t > B, std::complex< real_t > C,
                             std::complex< real_t >& r1, std::complex< real_t >& r2) {
      const std::complex< real_t > disc = std::sqrt(B * B - real_t(4) * C);
      r1                                = (-B + disc) / real_t(2);
      r2                                = (-B - disc) / real_t(2);
   };

   std::complex< real_t > y1, y2, y3, y4;
   solveQuadratic(-sqrt2m, std::complex< real_t >(p / real_t(2) + m) + qTerm, y1, y2);
   solveQuadratic(sqrt2m, std::complex< real_t >(p / real_t(2) + m) - qTerm, y3, y4);

   roots[0] = y1 + shift;
   roots[1] = y2 + shift;
   roots[2] = y3 + shift;
   roots[3] = y4 + shift;
}

//! Real roots of c4*t^4 + c3*t^3 + c2*t^2 + c1*t + c0 = 0 via the Aberth-Ehrlich (Ehrlich-Aberth)
//! method: same simultaneous-refinement family as solveQuarticDurandKerner (all 4 roots updated
//! together each sweep, no resolvent-cubic branching), but each step applies a full Newton correction
//! w_i = P(x_i)/P'(x_i) with an Aberth correction term (the sum of reciprocal distances to the other
//! current root estimates) rather than Durand-Kerner's plain Weierstrass product-of-differences
//! correction. This gives locally cubic convergence per root (vs Durand-Kerner's quadratic), so it
//! should in principle reach the same accuracy in fewer sweeps -- exactly what
//! QuarticSolverComparisonTest.cpp checks directly, both at matched iteration count (30) and at a
//! reduced count, rather than assuming the textbook convergence-order claim carries over unchanged to
//! this file's specific coefficient scales. Not wired into any production call site.
inline void solveQuarticAberthEhrlich(real_t c4, real_t c3, real_t c2, real_t c1, real_t c0,
                                       std::complex< real_t > roots[4], int maxIterations = 30)
{
   const real_t b3 = c3 / c4, b2 = c2 / c4, b1 = c1 / c4, b0 = c0 / c4;
   auto         poly = [&](std::complex< real_t > z) {
      return ((z + b3) * z + b2) * z * z + b1 * z + b0;
   };
   auto dpoly = [&](std::complex< real_t > z) {
      return ((real_t(4) * z + real_t(3) * b3) * z + real_t(2) * b2) * z + b1;
   };

   const real_t scale = real_t(1) + std::pow(std::max(std::max(std::abs(b3), std::abs(b2)),
                                                        std::max(std::abs(b1), std::abs(b0))),
                                              real_t(0.25));
   const std::complex< real_t > base(real_t(0.4), real_t(0.9)); // same seed as Durand-Kerner
   std::complex< real_t >       p = base;
   for (int k = 0; k < 4; ++k)
   {
      roots[k] = scale * p;
      p *= base;
   }

   for (int it = 0; it < maxIterations; ++it)
   {
      std::complex< real_t > w[4];
      for (int i = 0; i < 4; ++i)
      {
         const std::complex< real_t > dz = dpoly(roots[i]);
         w[i] = (std::abs(dz) > real_t(1e-300)) ? poly(roots[i]) / dz : std::complex< real_t >(0);
      }
      std::complex< real_t > next[4];
      for (int i = 0; i < 4; ++i)
      {
         std::complex< real_t > sumRecip(real_t(0), real_t(0));
         for (int j = 0; j < 4; ++j)
            if (j != i) sumRecip += real_t(1) / (roots[i] - roots[j]);
         const std::complex< real_t > denom = real_t(1) - w[i] * sumRecip;
         next[i] = (std::abs(denom) > real_t(1e-300)) ? roots[i] - w[i] / denom : roots[i] - w[i];
      }
      for (int i = 0; i < 4; ++i)
         roots[i] = next[i];
   }
}

//! Selects which of the three quartic solvers above (see QuarticSolverComparisonTest.cpp for the
//! time/accuracy comparison). [Claude, 2026-09-29] Ferrari is the production default everywhere
//! (closestThetaOnCircleToLine / findRimCrossing / cylinderCylinderContact); DurandKerner and
//! AberthEhrlich are kept only as comparison alternatives. (This comment previously said the opposite,
//! which was outdated.)
enum class QuarticSolverMethod
{
   DurandKerner,
   Ferrari,
   AberthEhrlich
};

//! Dispatches to one of the three quartic solvers by QuarticSolverMethod tag, so every call site needs
//! only one solver-choice parameter instead of duplicating itself per method. maxIterations is ignored
//! for Ferrari (closed-form, nothing to iterate).
inline void solveQuarticGeneric(real_t c4, real_t c3, real_t c2, real_t c1, real_t c0,
                                 std::complex< real_t > roots[4], QuarticSolverMethod method,
                                 int maxIterations = 30)
{
   switch (method)
   {
   case QuarticSolverMethod::Ferrari:
      solveQuarticFerrari(c4, c3, c2, c1, c0, roots);
      return;
   case QuarticSolverMethod::AberthEhrlich:
      solveQuarticAberthEhrlich(c4, c3, c2, c1, c0, roots, maxIterations);
      return;
   default:
      solveQuarticDurandKerner(c4, c3, c2, c1, c0, roots, maxIterations);
      return;
   }
}

//! Closed-form global minimizer of theta -> x(theta)^2+y(theta)^2, where x(theta)=x0+a*cos(theta),
//! y(theta)=y0+b*sin(theta) traces an ellipse (semi-axes a,b, centre (x0,y0)) -- i.e. the closest
//! point on an ellipse to the ORIGIN. This is the same problem as Feng, Han & Owen (2017)'s Appendix
//! A ("A generic contact detection framework for cylindrical particles in discrete element
//! modelling," Comput. Methods Appl. Mech. Engrg. 315, 632-651), used there to test a cylinder's
//! base-disc edge (projected as an ellipse into the OTHER cylinder's cross-sectional plane) against
//! that other cylinder's radius.
//!
//! Re-derived here from scratch rather than transcribing that paper's own Appendix A coefficients:
//! the intermediate step there ("F(theta) = -a*x0*sin(theta)+b*y0*cos(theta) = -a*x*x'+b*y*y'", then
//! "F(theta)=0 gives a*x*x'=b*y*y'") is inconsistent with its own stated substitution x'=cos(theta),
//! y'=sin(theta) two lines above it -- substituting directly gives a*x0*y'=b*y0*x', not a*x*x'=b*y*y'.
//! Rather than risk propagating whichever side of that inconsistency is the actual typo, this
//! rederives the quartic independently (exactly the same Weierstrass-substitution technique as
//! closestThetaOnCircleToLine below) and validates it against a brute-force numeric sweep in
//! PaperAnalyticalContactTest.cpp before trusting it anywhere else.
//!
//! f(theta) = (x0+a*cos theta)^2 + (y0+b*sin theta)^2; F(theta) = f'(theta)/2 =
//! -a*x0*sin(theta) + b*y0*cos(theta) + (b^2-a^2)*sin(theta)*cos(theta) -- a mix of first- and
//! second-harmonic terms. Substituting t=tan(theta/2) and clearing denominators gives a quartic in t
//! with no cubic-resolvent branching needed (Durand-Kerner, same as every other quartic in this file).
//! All (up to 5, including the t->infinity/theta=pi point the Weierstrass substitution can't reach
//! directly) real critical thetas of theta -> x(theta)^2+y(theta)^2 -- i.e. every local min AND max
//! on the ellipse, not just the global minimum. Needed by the constrained version of this search
//! (closestThetaOnEllipseWithinAxialRange below): when the unconstrained global minimum falls outside
//! a required range, the CONSTRAINED minimum is generically at another critical point or at the
//! range boundary, not simply "no valid point at all."
inline void allCriticalThetasOnEllipseToOrigin(real_t a, real_t b, real_t x0, real_t y0,
                                                real_t (&thetasOut)[5], int& countOut)
{
   countOut          = 1;
   thetasOut[0]      = math::pi;

   const real_t c4 = -b * y0;
   const real_t c3 = real_t(2) * (a * a - b * b - a * x0);
   const real_t c1 = real_t(-2) * (a * a - b * b + a * x0);
   const real_t c0 = b * y0;

   const bool   c4Vanishes = std::abs(c4) < real_t(1e-14) * std::max(real_t(1), std::abs(c0));
   const real_t c4Safe     = c4Vanishes ? (std::abs(c4) > real_t(1e-30) ? c4 : real_t(1e-14)) : c4;

   std::complex< real_t > roots[4];
   solveQuarticDurandKerner(c4Safe, c3, real_t(0), c1, c0, roots);
   for (const auto& r : roots)
      if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
         thetasOut[countOut++] = real_t(2) * std::atan(r.real());
}

//! See allCriticalThetasOnEllipseToOrigin's own comment: this is just its global-minimum reduction,
//! kept as the simple entry point for callers (and the brute-force cross-check test) that don't need
//! the constrained search.
inline real_t closestThetaOnEllipseToOrigin(real_t a, real_t b, real_t x0, real_t y0)
{
   auto distSq = [&](real_t theta) {
      const real_t x = x0 + a * std::cos(theta);
      const real_t y = y0 + b * std::sin(theta);
      return x * x + y * y;
   };

   real_t thetas[5];
   int    count;
   allCriticalThetasOnEllipseToOrigin(a, b, x0, y0, thetas, count);

   real_t bestTheta = thetas[0];
   real_t bestVal   = distSq(bestTheta);
   for (int i = 1; i < count; ++i)
   {
      const real_t val = distSq(thetas[i]);
      if (val < bestVal)
      {
         bestVal   = val;
         bestTheta = thetas[i];
      }
   }
   return bestTheta;
}

//! Closest point on the ellipse to the origin, CONSTRAINED to the sub-arc where a second quantity,
//! zHost(theta) = zC + zA*cos(theta) + zB*sin(theta) (itself just a first-harmonic function of theta,
//! since it's the axial coordinate of the same rim point along a fixed third direction), stays within
//! [-zHalfRange, zHalfRange]. This is the exact fix for the gap Feng et al. (2017) themselves flag in
//! their own unconstrained version of this check ("the master cylinder is implicitly assumed
//! infinitely long, and no consideration is given to a possible overlap between the base discs") --
//! but cutting the OTHER way from how they describe it: their concern is a false POSITIVE (finding a
//! close point that's actually beyond the host's finite end); fuzz-testing this implementation found
//! the reverse is at least as important in practice -- a false NEGATIVE, where the unconstrained
//! global minimum falls just outside the valid range while a nearby, still-overlapping point on the
//! same ellipse remains within it. Solving the constrained problem in closed form (rather than simply
//! rejecting the unconstrained answer) closes both directions at once.
//!
//! Returns false if the valid sub-arc is empty (the entire rim lies outside the host's finite range).
inline bool closestThetaOnEllipseWithinAxialRange(real_t a, real_t b, real_t x0, real_t y0, real_t zC,
                                                   real_t zA, real_t zB, real_t zHalfRange, real_t& thetaOut)
{
   auto distSq = [&](real_t theta) {
      const real_t x = x0 + a * std::cos(theta);
      const real_t y = y0 + b * std::sin(theta);
      return x * x + y * y;
   };
   auto zHost = [&](real_t theta) { return zC + zA * std::cos(theta) + zB * std::sin(theta); };
   auto inRange = [&](real_t theta) { return std::abs(zHost(theta)) < zHalfRange; };

   const real_t zAmag = std::sqrt(zA * zA + zB * zB);

   bool   have = false;
   real_t bestVal = real_t(0);
   auto   consider = [&](real_t theta) {
      if (!inRange(theta)) return;
      const real_t val = distSq(theta);
      if (!have || val < bestVal)
      {
         have     = true;
         bestVal  = val;
         thetaOut = theta;
      }
   };

   // Candidate set 1: the unconstrained critical points (global/local minima AND maxima) that
   // happen to already satisfy the range constraint.
   real_t thetas[5];
   int    count;
   allCriticalThetasOnEllipseToOrigin(a, b, x0, y0, thetas, count);
   for (int i = 0; i < count; ++i)
      consider(thetas[i]);

   // Candidate set 2: the range boundary itself. zHost(theta) = +-zHalfRange is a plain
   // R*cos(theta-phi) = k first-harmonic equation (no quartic needed) -- if |k/R| > 1 the entire rim
   // lies on one side (either always in range, already covered by set 1, or never in range, and this
   // loop simply adds nothing).
   if (zAmag > real_t(1e-12))
   {
      const real_t phi = std::atan2(zB, zA);
      for (real_t target : { zHalfRange, -zHalfRange })
      {
         const real_t k = (target - zC) / zAmag;
         if (std::abs(k) <= real_t(1))
         {
            const real_t delta = std::acos(std::clamp(k, real_t(-1), real_t(1)));
            consider(phi + delta);
            consider(phi - delta);
         }
      }
   }

   return have;
}

//! Feng, Han & Owen (2017) Appendix B: intersection points between a circle of radius R at the
//! origin and an ellipse (semi-axes a, b, centre (x0,y0)) -- i.e. where the ellipse's rim actually
//! CROSSES the circle boundary, as opposed to closestThetaOnEllipseToOrigin's closest-POINT problem.
//! Returns up to 4 real thetas via thetasOut/countOut (countOut may be 0, 2, or 4; a tangency counts
//! as a repeated root).
//!
//! Re-derived from scratch rather than transcribed from the paper's own printed equations, which
//! have two internal inconsistencies here: (1) the circle is introduced as "radius r" but then
//! written as x^2+y^2=1 (not r^2), while the resulting quartic's own a0=r^2(1-tau^2) coefficient only
//! makes sense if the circle equation actually was x^2+y^2=r^2; (2) Appendix C's eq. (33), a
//! one-line derivation from its own eq. (31)-(32) just above it, drops a minus sign (redoing that
//! one line independently does not reproduce what eq. (33) prints). Given two independent
//! transcription errors already found in this paper (see closestThetaOnEllipseToOrigin's own
//! comment for a third, in Appendix A), this quartic is derived fresh here and cross-checked at
//! theta=0, pi, and pi/2 (each independently re-deriving the corresponding coefficient) before ever
//! being trusted, then validated against a brute-force numeric sweep in
//! PaperAnalyticalContactTest.cpp, the same policy as every other quartic in this file.
//!
//! Derivation: substituting the ellipse's own parametrisation into x^2+y^2=R^2 and expanding
//! cos^2(theta), sin^2(theta) via the double-angle identities gives a mix of first- and second-
//! harmonic terms in theta, exactly the same shape as every other trig equation in this file;
//! substituting t=tan(theta/2) and clearing denominators gives:
inline void intersectEllipseWithCircle(real_t a, real_t b, real_t x0, real_t y0, real_t R,
                                        real_t (&thetasOut)[4], int& countOut)
{
   countOut = 0;

   const real_t c4 = (x0 - a) * (x0 - a) + y0 * y0 - R * R; // = value at theta=pi (t -> infinity)
   const real_t c0 = (x0 + a) * (x0 + a) + y0 * y0 - R * R; // = value at theta=0 (t=0)
   const real_t c3 = real_t(4) * b * y0;
   const real_t c1 = c3;
   const real_t c2 = real_t(2) * x0 * x0 + real_t(2) * y0 * y0 - real_t(2) * R * R - real_t(2) * a * a +
                      real_t(4) * b * b;

   // theta=pi is itself an exact root exactly when c4=0 (the point at theta=pi lies on the circle);
   // the Weierstrass substitution can't reach it as a finite t, so it must be added explicitly rather
   // than relying on solveQuarticDurandKerner to find it as t->infinity.
   if (std::abs(c4) < real_t(1e-12) * std::max(real_t(1), std::abs(c0)))
   {
      thetasOut[countOut++] = math::pi;
      // Remaining roots solve the cubic c3 t^3 + c2 t^2 + c1 t + c0 = 0 -- reuse the quartic solver
      // with c4 nudged to a small epsilon (same technique used throughout this file), which just
      // sends one of the four roots towards a very large |t| (theta -> pi again, already added).
      const real_t c4Safe = real_t(1e-14);
      std::complex< real_t > roots[4];
      solveQuarticDurandKerner(c4Safe, c3, c2, c1, c0, roots);
      for (const auto& r : roots)
         if (countOut < 4 && std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())) &&
             std::abs(r.real()) < real_t(1e6)) // exclude the spurious large-|t| root standing in for pi
            thetasOut[countOut++] = real_t(2) * std::atan(r.real());
      return;
   }

   std::complex< real_t > roots[4];
   solveQuarticDurandKerner(c4, c3, c2, c1, c0, roots);
   for (const auto& r : roots)
      if (countOut < 4 && std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
         thetasOut[countOut++] = real_t(2) * std::atan(r.real());
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
// maxIterations forwarded to solveQuarticDurandKerner (default 30, unchanged behaviour for every
// existing caller that doesn't pass it) -- see that function's own comment.
// quarticSolver: default changed to Ferrari (this session, after QuarticSolverComparisonTest.cpp's
// synthetic-quartic comparison and PaperAnalyticalContactTest.cpp's edge-edge N-way comparison on real
// geometry both showed it ~2.3x faster than Durand-Kerner with bit-identical-or-indistinguishable
// accuracy) -- see QuarticSolverMethod's own comment. Durand-Kerner/Aberth-Ehrlich remain fully
// available via this parameter, not deleted.
inline real_t closestThetaOnCircleToLine(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                          const Vec3& linePos, const Vec3& lineDir, int maxIterations = 30,
                                          QuarticSolverMethod quarticSolver = QuarticSolverMethod::Ferrari)
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
      solveQuarticGeneric(c4Safe, c3, c2, c1, c0, roots, quarticSolver, maxIterations);
      for (const auto& r : roots)
         if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
            consider(real_t(2) * std::atan(r.real()));
      return bestTheta;
   }

   std::complex< real_t > roots[4];
   solveQuarticGeneric(c4, c3, c2, c1, c0, roots, quarticSolver, maxIterations);
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

// useFastReference: OPTIONAL simplification #1 from this session's own performance investigation
// (see CylinderCylinderContactBenchmark.cpp's file header) -- thetaRef exists PURELY as a tie-break
// among point A's (up to 4) candidate roots (see this function's own RimCrossing comment above: "the
// choice here... matches the... regime", i.e. this file's own predecessor's heuristic, not something
// Kodam's or Guo's paper specifies at all). false (DEFAULT, unchanged behaviour) computes it via
// closestThetaOnCircleToLine -- the exact closest point on the rim to target's full axis LINE, itself
// a second 60(30)-iteration Durand-Kerner quartic solve, on top of point A's own. true replaces that
// with closestThetaOnCircleToPoint (a single atan2, already validated closed-form -- see that
// function's own comment: "turned out to have an even simpler exact answer, no polynomial at all"),
// using target's CENTRE POSITION as a cheap stand-in reference point instead of its full axis line --
// geometrically a different (cheaper, approximate) reference, not merely a faster way to compute the
// SAME reference, so this is genuinely an accuracy/speed trade-off, not a free simplification. Old
// path kept fully intact (not deleted) specifically so the two can be A/B compared for accuracy --
// see PaperAnalyticalContactTest.cpp's fuzz suite, run against both settings.
// maxIterations forwarded to every Durand-Kerner/Aberth-Ehrlich solve this function performs (thetaRef
// when NOT useFastReference, and point A's own quartic always) -- default 30, unchanged behaviour for
// both of this function's existing callers (cylinderCylinderContact's edge-edge branch, which now
// forwards its own maxIterations parameter -- see that function's own comment).
// quarticSolver: forwarded to both of this function's own quartic solves (thetaRef's, when not using
// the fast point-based reference, and point A's own) -- default Ferrari (see closestThetaOnCircleToLine's
// own comment for why); see QuarticSolverMethod's own comment.
// candidateCountA/candidateCountB: OPTIONAL diagnostic out-parameters (default nullptr, zero effect on
// production behaviour) -- count how many candidate thetas pass EACH of point A's/B's own containment
// checks, independent of thetaRef's tie-break, to answer directly whether that tie-break (and by
// extension thetaRef's own quartic, when not using the fast reference) is ever actually load-bearing
// for real geometry, or whether at most one candidate ever survives in practice. See
// PaperAnalyticalContactTest.cpp's testThetaRefTieBreakNecessity10k, which found 0/832 real calls ever
// had more than one candidate for EITHER point.
// lazyThetaRef: OPTIONAL simplification (this session, on top of #1 above) -- defers thetaRef's
// computation until AFTER point A's (and, if needed, point B's) candidates are actually collected, and
// only computes it AT ALL if either list ends up with more than one entry (the collect-then-select
// restructuring below is mathematically identical to the old inline considerA/considerB running-best
// tracking either way -- same traversal order, same strict "<" tie-break, so ties resolve to the same
// candidate). Default changed to true this session, after testThetaRefTieBreakNecessity10k found this
// is NEVER needed on real geometry (0/832 calls ever had >1 candidate for either point) and
// testLazyThetaRefAccuracy10k confirmed bit-identical results with it enabled. false (the original
// eager behaviour, unconditional thetaRef computation up front) remains fully available.
inline RimCrossing findRimCrossing(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                    const Vec3& targetPos, const Vec3& eTarget, real_t targetRadius,
                                    real_t targetHalfLength, bool useFastReference = false,
                                    int maxIterations = 30,
                                    // [Claude, 2026-09-29] default changed DurandKerner -> Ferrari to match
                                    // every other entry point (production always passed Ferrari explicitly
                                    // via cylinderCylinderContact; the only direct test call passes a solver).
                                    QuarticSolverMethod quarticSolver = QuarticSolverMethod::Ferrari,
                                    int* candidateCountA = nullptr, int* candidateCountB = nullptr,
                                    bool lazyThetaRef = true)
{
   if (candidateCountA) *candidateCountA = 0;
   if (candidateCountB) *candidateCountB = 0;
   RimCrossing out;

   real_t thetaRef         = real_t(0);
   bool   thetaRefComputed = false;
   auto   ensureThetaRef   = [&]() {
      if (thetaRefComputed) return;
      thetaRef = useFastReference
         ? closestThetaOnCircleToPoint(C, u, v, R, targetPos)
         : closestThetaOnCircleToLine(C, u, v, R, targetPos, eTarget, maxIterations, quarticSolver);
      thetaRefComputed = true;
   };
   if (!lazyThetaRef) ensureThetaRef(); // EAGER (default): unconditional, matches original cost exactly.

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
   // The quartic below is Kodam's own eq. 83 -- unavoidable, unaffected by lazyThetaRef. What changes
   // is only WHEN (if ever) thetaRef is computed to pick among the candidates this produces.
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

   // Every candidate theta passing containment is collected here FIRST (independent of thetaRef, which
   // this containment check never references) -- at most 5: up to 4 quartic roots plus theta=pi.
   real_t candidatesA[5];
   int    numA = 0;
   auto   collectA = [&](real_t theta) {
      const Vec3 P = pointAt(theta);
      if (std::abs(radialSqToAxis(P) - targetRadius * targetRadius) >
          real_t(1e-6) * std::max(real_t(1), targetRadius * targetRadius))
         return;
      if (std::abs(axialToAxis(P)) >= targetHalfLength) return;
      if (numA < 5) candidatesA[numA++] = theta;
   };

   if (std::abs(c4) > real_t(1e-12) * std::max(real_t(1), std::abs(c0)))
   {
      std::complex< real_t > roots[4];
      solveQuarticGeneric(c4, c3, c2, c1, c0, roots, quarticSolver, maxIterations);
      for (const auto& r : roots)
         if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
            collectA(real_t(2) * std::atan(r.real()));
   }
   collectA(math::pi); // t -> infinity point the Weierstrass substitution can't reach directly

   if (candidateCountA) *candidateCountA = numA;

   // Selection: 0 candidates -> no point A at all; exactly 1 -> thetaRef is irrelevant, skip it
   // entirely (this is the actual savings when lazyThetaRef); >1 -> thetaRef genuinely needed to break
   // the tie (computed here, lazily, the FIRST time this ever actually happens for either point).
   bool   haveA  = false;
   real_t thetaA = real_t(0);
   if (numA == 1) { thetaA = candidatesA[0]; haveA = true; }
   else if (numA > 1)
   {
      ensureThetaRef();
      real_t bestDiffA = std::numeric_limits< real_t >::max();
      for (int i = 0; i < numA; ++i)
      {
         const real_t d = angDiff(candidatesA[i], thetaRef);
         if (d < bestDiffA)
         {
            bestDiffA = d;
            thetaA    = candidatesA[i];
            haveA     = true;
         }
      }
   }

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

   real_t candidatesB[2];
   int    numB = 0;
   auto   collectB = [&](real_t theta) {
      const Vec3 P = pointAt(theta);
      if (radialSqToAxis(P) >= targetRadius * targetRadius) return;
      if (numB < 2) candidatesB[numB++] = theta;
   };
   collectB(phi + delta);
   collectB(phi - delta);

   if (candidateCountB) *candidateCountB = numB;

   bool   haveB  = false;
   real_t thetaB = real_t(0);
   if (numB == 1) { thetaB = candidatesB[0]; haveB = true; }
   else if (numB > 1)
   {
      ensureThetaRef(); // no-op if point A's own selection already computed it above
      real_t bestDiffB = std::numeric_limits< real_t >::max();
      for (int i = 0; i < numB; ++i)
      {
         const real_t d = angDiff(candidatesB[i], thetaRef);
         if (d < bestDiffB)
         {
            bestDiffB = d;
            thetaB    = candidatesB[i];
            haveB     = true;
         }
      }
   }

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

//! [Claude, 2026-09-29] ALL crossings of a rim circle (centre C, in-plane basis u/v, radius R) with the
//! INFINITE lateral surface of a target cylinder (axis through targetPos along eTarget, radius
//! targetRadius) -- Kodam et al. (2010) eq. 83, the same quartic findRimCrossing solves for its point A.
//! findRimCrossing keeps only ONE crossing (tie-broken by a reference angle, as its edge-edge chord
//! construction needs); the corrected special-contact depth below needs every crossing, so this returns
//! all real roots that satisfy the crossing condition (up to 4, plus theta = pi, which the Weierstrass
//! substitution t = tan(theta/2) cannot reach). Deliberately a separate function: findRimCrossing is left
//! byte-for-byte unchanged so the edge-edge results cannot shift. Returns the number of thetas written.
inline int allRimCrossingsWithLateralSurface(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                             const Vec3& targetPos, const Vec3& eTarget, real_t targetRadius,
                                             real_t thetasOut[5])
{
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

   int  n = 0;
   auto accept = [&](real_t theta) {
      const Vec3   P      = C + R * std::cos(theta) * u + R * std::sin(theta) * v;
      const Vec3   rel    = P - targetPos;
      const Vec3   radial = rel - (rel * eTarget) * eTarget;
      if (std::abs(radial.sqrLength() - targetRadius * targetRadius) >
          real_t(1e-6) * std::max(real_t(1), targetRadius * targetRadius))
         return;
      if (n < 5) thetasOut[n++] = theta;
   };
   if (std::abs(c4) > real_t(1e-12) * std::max(real_t(1), std::abs(c0)))
   {
      std::complex< real_t > roots[4];
      solveQuarticGeneric(c4, c3, c2, c1, c0, roots, QuarticSolverMethod::Ferrari);
      for (const auto& r : roots)
         if (std::abs(r.imag()) < real_t(1e-6) * std::max(real_t(1), std::abs(r.real())))
            accept(real_t(2) * std::atan(r.real()));
   }
   accept(math::pi);
   return n;
}

//! Kodam et al. (2010) eq. 40-44: exact face-band contact point, replacing an earlier version of
//! this file that approximated it by reusing the face-face/band-band-parallel point construction
//! (faceFaceStyleContactPoint) with 2*cj.radius substituted for the radial overlap -- verified by
//! hand against the paper's own PDF (Kodam, Bharadwaj, Curtis, Hancock & Wassgren (2010), Chemical
//! Engineering Science 65(22), 5852-5862) to be a materially different point in general, not merely
//! a stylistic simplification: for PaperAnalyticalContactTest.cpp's own testCylinderCylinderFaceBand
//! geometry, the approximation returns cylinder i's own axis (0,0) purely because that test happens
//! to use equal radii (ci.radius == cj.radius, cancelling eq. 43's (Ri-Rj) term to zero); the exact
//! construction below returns (0, 0.45) instead -- offset toward cylinder j's actual centre, as
//! physically expected.
//!
//! `relI` is cylinder j's centre in cylinder i's body frame (so cylinder i's near face lies at
//! z = sign(relI[2])*halfLi); `axisJInI` is cylinder j's axis expressed in the SAME frame; `halfLj`
//! is cylinder j's half-length; `Ri` is cylinder i's radius; `overlap` is the already-computed axial
//! overlap (eq. 38). Returns the contact point in cylinder i's body frame (caller transforms to
//! world via ci.pos + ci.R * result).
//!
//! Eq. 40-41 construct line AB by taking cylinder j's own centreline and translating it (purely
//! along cylinder i's z-axis, NOT radially -- there is no cj.radius term in eq. 40-41 at all) onto
//! cylinder i's near face plane. Eq. 42/43 place the contact point at AB's own midpoint when AB lies
//! fully inside i's face disc (radius Ri); eq. 44 places it at the midpoint of B and the disc-
//! boundary crossing point D when AB pokes out through one end, a line-circle intersection the paper
//! declines to spell out ("the solution is straightforward and will not be presented here"). Solved
//! here directly via the standard finite line-segment/circle intersection quadratic, which reduces
//! exactly to eq. 43 (both endpoints inside -> tLo=-halfLj, tHi=halfLj -> midpoint = AB's own centre)
//! and eq. 44 (one endpoint inside -> tLo or tHi = the surviving endpoint, the other = the crossing
//! point) in their respective regimes, and additionally covers a case neither eq. 42 nor eq. 44
//! explicitly names -- AB poking out through BOTH ends (a long, narrow cylinder j straddling a much
//! smaller face i, as in the geometry above) -- by midpointing the two crossing points instead.
inline Vec3 faceBandContactPoint(const Vec3& relI, const Vec3& axisJInI, real_t halfLj, real_t Ri,
                                  real_t halfLi, real_t overlap)
{
   const real_t faceZ = math::sign(relI[2]) * halfLi;
   const real_t cx = relI[0], cy = relI[1];
   const real_t ux = axisJInI[0], uy = axisJInI[1];

   const real_t k  = ux * ux + uy * uy;
   const real_t m  = cx * ux + cy * uy;
   const real_t c2 = cx * cx + cy * cy - Ri * Ri;

   real_t tLo = -halfLj, tHi = halfLj; // fallback assumption: fully inside (eq. 43)
   if (k > real_t(1e-12))
   {
      const real_t disc = m * m - k * c2;
      bool         clampToNearest = (disc < real_t(0));
      if (!clampToNearest)
      {
         const real_t sq = std::sqrt(disc);
         const real_t t1 = (-m - sq) / k, t2 = (-m + sq) / k; // t1 <= t2: infinite-line/disc crossings
         tLo = std::max(t1, -halfLj);
         tHi = std::min(t2, halfLj);
         clampToNearest = (tLo >= tHi); // no surviving overlap on the finite segment
      }
      if (clampToNearest)
      {
         // Shouldn't occur given the caller's own containment gate (eq. 34/36/37) already established
         // a real contact; if floating-point edge cases ever land here anyway, degrade gracefully to
         // the point on the finite segment nearest cylinder i's axis rather than an undefined chord.
         const real_t tNear = std::clamp(-m / k, -halfLj, halfLj);
         tLo = tHi = tNear;
      }
   }

   const real_t tMid = real_t(0.5) * (tLo + tHi);
   const real_t pz   = faceZ - math::sign(relI[2]) * real_t(0.5) * overlap;
   return Vec3(cx + tMid * ux, cy + tMid * uy, pz);
}

//! Closest squared distance between two finite line segments P1->Q1 and P2->Q2 -- thin wrapper around
//! detail::closestPointsBetweenSegments (defined near the top of this file, before boxBoxSAT, since
//! that function needs it too -- kept as its own entry point since most callers, e.g. the
//! cylinder-cylinder capsule bound below, only need the distance, not the witness points).
inline real_t closestDistanceSqBetweenSegments(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2)
{
   Vec3 c1, c2;
   closestPointsBetweenSegments(p1, q1, p2, q2, c1, c2);
   return (c1 - c2).sqrLength();
}

//! The two endpoints of a cylinder's own axis segment (its centreline), in world coordinates.
inline void cylinderAxisEndpoints(const CylinderState& c, Vec3& p0, Vec3& p1)
{
   const Vec3 halfAxis = (real_t(0.5) * c.length) * c.axis();
   p0                  = c.pos - halfAxis;
   p1                  = c.pos + halfAxis;
}

//! An oriented box's own support-function value (half-width of its projection) along a given WORLD
//! direction -- dirWorld MUST be unit length. Standard OBB-projection formula (the same one full SAT
//! uses internally for each of its 15 candidate axes, here applied to just ONE direction): the box's
//! projected half-extent is the sum of each local half-extent weighted by how aligned that local axis
//! is with dirWorld.
inline real_t boxSupportAlongDirection(const Vec3& halfExtents, const Mat3& R, const Vec3& dirWorld)
{
   return halfExtents[0] * std::abs(dirWorld * Vec3(R(0, 0), R(1, 0), R(2, 0))) +
          halfExtents[1] * std::abs(dirWorld * Vec3(R(0, 1), R(1, 1), R(2, 1))) +
          halfExtents[2] * std::abs(dirWorld * Vec3(R(0, 2), R(1, 2), R(2, 2)));
}

//! CONSERVATIVE box-box separation test using only ONE candidate separating axis (the line between
//! the two centres) instead of full SAT's 15 -- cheap (a handful of dot products, no iteration), at
//! the cost of completeness: unlike full SAT, this can produce false NEGATIVES (says "can't rule out
//! contact" for a pair that's actually separated along some other axis, e.g. a true edge-edge
//! separating axis this never tests), but per the separating-axis theorem it can NEVER produce a
//! false positive rejection -- if the two boxes' projections onto this one legitimate candidate axis
//! don't overlap, the boxes themselves provably don't overlap, regardless of what any other axis
//! would have shown. Intended purely as a cheap Layer 1 pre-filter ahead of the real (exact) boxBoxSAT
//! test, exactly mirroring closestDistanceSqBetweenSegments's role for cylinders -- NOT a replacement
//! for SAT, which remains the sole source of truth for actual contact/overlap/normal.
inline bool boxesDefinitelySeparatedAlongCenterLine(const Vec3& posA, const Mat3& RA, const Vec3& halfA,
                                                     const Vec3& posB, const Mat3& RB, const Vec3& halfB)
{
   const Vec3   delta = posB - posA;
   const real_t dist  = delta.length();
   if (dist < real_t(1e-12)) return false; // coincident centres -- can't determine anything from this axis
   const Vec3 dir = delta / dist;
   const real_t reach = boxSupportAlongDirection(halfA, RA, dir) + boxSupportAlongDirection(halfB, RB, dir);
   return dist > reach;
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
   // overlap = halfL - |zCyl|  (Kodam et al. 2010 Part I, eq. 4-6)
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
   // overlap = R - |zCyl|  (Kodam et al. 2010 Part I, eq. 10-12)
   if (std::abs(axisDotN) < real_t(1e-8) && std::abs(zCyl) < cyl.radius)
   {
      result.contact  = true;
      result.overlap  = cyl.radius - std::abs(zCyl);                            // eq. 10
      result.normal   = math::sign(zCyl) * nHat;                                // eq. 11
      result.point    = centreOnPlane - math::sign(zCyl) * real_t(0.5) * result.overlap * nHat; // eq. 12
      result.scenario = "band-plane";
      return result;
   }

   // --- Edge-plane (eq. 13-21): general (tilted) orientation, Kodam et al. 2010 Part I -----------
   // Point E, the point on the cylinder's rim closest to the plane (eq. 15-17), via the shared
   // detail::rimPointTowardPlane (see its comment for the corrected-formula derivation).
   const Vec3   E  = detail::rimPointTowardPlane(cyl, planePos, nHat);
   const real_t zE = (E - planePos) * nHat;

   if (zE * zCyl < real_t(0)) // eq. 18: opposite sides of the plane => contact
   {
      result.contact  = true;
      result.overlap  = std::abs(zE);                                          // eq. 19
      result.normal   = math::sign(zCyl) * nHat;                               // eq. 20
      // eq. 21 is z_E^P - sign(z_E^P)*0.5*delta, i.e. measured from E's OWN height, not from the
      // plane (z=0). An earlier version first projected E onto the plane (eOnPlane = E - zE*nHat,
      // height forced to exactly 0) before applying this offset, which silently shifts the result by
      // an extra zE and places the point on the wrong side of the plane entirely -- e.g. zCyl=+2,
      // zE=-1 (a genuine overlap, midpoint physically at z=-0.5) previously produced z=+0.5, outside
      // the overlap region altogether. rimVsFaceOrBand's structurally identical face-edge branch
      // (below) never had this bug, since it applies the offset to E directly; this now matches it.
      result.point    = E - math::sign(zE) * real_t(0.5) * result.overlap * nHat; // eq. 21
      result.scenario = "edge-plane";
   }
   return result;
}

//! Sphere-Cylinder contact -- NOT part of Kodam et al. or Guo et al. (neither paper addresses mixed
//! shape pairs at all); developed for this project, unlike sphere-sphere and sphere-box, which predate
//! it in mesa_pd's own AnalyticCollisionFunctions.h (detectSphereSphereCollision/detectSphereBoxCollision,
//! S. Eibl 2019). Exact and closed-form, and -- like boxBoxSAT -- needs no GJK/EPA fallback at all: the
//! cylinder's own closest-point-to-an-external-point mapping (the natural "clamp" counterpart to its
//! own support mapping, exactly the relationship detectSphereBoxCollision's box clamp has to the box's
//! own support mapping) is a COMPLETE, exact solution to the sphere-cylinder minimum-translation-
//! distance problem, not a named-scenario heuristic the way Kodam et al.'s cylinder-cylinder dispatch
//! is -- there is no configuration this construction fails to cover, so unlike cylinderCylinderContact
//! it carries no companion "is this result uncertain" gate.
//!
//! Convention matches this file throughout: normal points from the 2nd argument (the cylinder) toward
//! the 1st (the sphere); overlap positive = penetrating.
inline ContactResult sphereCylinderContact(const Vec3& spherePos, real_t sphereRadius, const CylinderState& cyl)
{
   ContactResult result;
   result.scenario = "sphere-cylinder";

   const Vec3   eZ    = cyl.axis();
   const real_t halfL = real_t(0.5) * cyl.length;
   const Vec3   rel   = spherePos - cyl.pos;      // sphere centre relative to cylinder centre, world frame
   const real_t z     = rel * eZ;                  // axial component
   const Vec3   radialVec = rel - z * eZ;           // radial component, world frame
   const real_t rho   = radialVec.length();

   // Arbitrary but well-defined radial direction for the degenerate case (sphere centre exactly ON
   // the cylinder's own axis, rho=0) -- same "pick any consistent perpendicular" resolution already
   // used throughout this file (detail::perpendicularBasis) rather than leaving a 0/0 direction.
   Vec3 uArb, vArb;
   detail::perpendicularBasis(eZ, uArb, vArb);
   const Vec3 radialDir = (rho > real_t(1e-12)) ? radialVec / rho : uArb;

   if (rho > cyl.radius || std::abs(z) > halfL)
   {
      // Sphere centre outside the (solid, capped) cylinder: closest surface point is the independent
      // radial/axial clamp -- covers band, cap, and rim proximity uniformly, with no separate case
      // count the way a box's 3 independent per-axis clamps need (a cylinder has only two surface
      // families, band and cap, instead of three pairs of faces).
      const real_t rhoC = std::min(rho, cyl.radius);
      const real_t zC   = std::clamp(z, -halfL, halfL);
      const Vec3   q    = rhoC * radialDir + zC * eZ; // closest point, relative to cylinder centre

      const Vec3   diff = rel - q;
      const real_t dist = diff.length();

      result.contact = dist < sphereRadius;
      result.overlap = sphereRadius - dist;
      result.normal  = (dist > real_t(1e-9)) ? diff / dist : radialDir;
      result.point   = cyl.pos + q;
   }
   else
   {
      // Sphere centre embedded inside the cylinder: push out through whichever surface (band or
      // nearer cap) has the least remaining clearance, the same "closest way out" logic
      // detectSphereBoxCollision uses for an embedded sphere centre, generalized from 3 face pairs
      // to this shape's 2 surface families.
      const real_t radialClearance = cyl.radius - rho;
      const real_t axialClearance  = halfL - std::abs(z);

      result.contact = true;
      if (radialClearance < axialClearance)
      {
         result.normal  = radialDir;
         result.overlap = radialClearance + sphereRadius;
      }
      else
      {
         result.normal  = math::sign(z) * eZ;
         result.overlap = axialClearance + sphereRadius;
      }
      result.point = spherePos;
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

} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla
