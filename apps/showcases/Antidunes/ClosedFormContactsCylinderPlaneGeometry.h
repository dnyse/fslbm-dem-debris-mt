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
//! faceFaceStyleContactPoint; shared ContactResult type in ClosedFormContactsInfra.h.

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
//! branching on the resolvent cubic's discriminant to get wrong. Fixed iteration count; empirically
//! converges to double precision within ~20-30 iterations for the geometries this is used on.
// maxIterations: default 30, validated against the full PaperAnalyticalContactTest.cpp suite (fuzz tests included);
// exposed so PaperAnalyticalContactTest.cpp can compare iteration counts on identical trials
// (testEdgeEdgeThreeWayAccuracy10k).
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
//! not all three. Branches between Cardano's formula (one-real-root case) and the trigonometric form
//! (three-real-roots case, avoids the cube root of a near-zero or catastrophically cancelling complex
//! quantity). This discriminant-sign branch is Ferrari's method's weak point; isolating it to a single
//! real root of a cubic (rather than needing all 4 quartic roots from the same branch) contains it.
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
//! This is the production solver: every quartic solve in the cylinder-cylinder dispatch defaults to
//! QuarticSolverMethod::Ferrari (closestThetaOnCircleToLine, findRimCrossing, cylinderCylinderContact); no
//! failing case is known. Theoretical weak points: (1) the resolvent cubic's discriminant sign is a branch
//! (handled in solveRealCubicOneRealRoot); (2) m = 0 (q = 0, biquadratic quartic) needs the special handling
//! below to avoid dividing by sqrt(2m); (3) like any closed-form quartic formula it can lose accuracy through
//! cancellation when roots nearly coincide, where iterative solvers degrade more gracefully. Durand-Kerner /
//! Aberth-Ehrlich remain for comparison only.
inline void solveQuarticFerrari(real_t c4, real_t c3, real_t c2, real_t c1, real_t c0,
                                 std::complex< real_t > roots[4])
{
   const real_t b3 = c3 / c4, b2 = c2 / c4, b1 = c1 / c4, b0 = c0 / c4;

   const real_t p = b2 - real_t(3) * b3 * b3 / real_t(8);
   const real_t q = b3 * b3 * b3 / real_t(8) - b3 * b2 / real_t(2) + b1;
   const real_t r = -real_t(3) * b3 * b3 * b3 * b3 / real_t(256) + b3 * b3 * b2 / real_t(16) -
                     b3 * b1 / real_t(4) + b0;
   const real_t shift = -b3 / real_t(4);

   // Biquadratic special case (q == 0, no odd-degree term in the depressed quartic): the general construction
   // below divides by sqrt(2m) with m = 0, so y^4+p*y^2+r=0 is solved directly as a quadratic in y^2.
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
//! correction. This gives locally cubic convergence per root (vs Durand-Kerner's quadratic);
//! QuarticSolverComparisonTest.cpp checks whether that holds at these coefficient scales, at matched (30)
//! and reduced iteration counts. Not wired into any production call site.
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

//! Selects one of the three quartic solvers above (time/accuracy comparison in QuarticSolverComparisonTest.cpp).
//! Ferrari is the production default (closestThetaOnCircleToLine / findRimCrossing / cylinderCylinderContact);
//! DurandKerner and AberthEhrlich are comparison alternatives.
enum class QuarticSolverMethod
{
   DurandKerner,
   Ferrari,
   AberthEhrlich
};

//! Dispatches to one of the three quartic solvers by QuarticSolverMethod tag. maxIterations is ignored for
//! Ferrari (closed form).
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
//! Re-derived here rather than transcribed from that paper's Appendix A: its intermediate step
//! ("F(theta) = -a*x*sin(theta)+b*y*cos(theta) = -a*x*x'+b*y*y'", then "F(theta)=0 gives a*x*x'=b*y*y'") is
//! inconsistent with its own substitution x'=cos(theta), y'=sin(theta) two lines above -- substituting directly
//! gives a*x*y'=b*y*x', not a*x*x'=b*y*y'. The quartic is therefore derived independently (same Weierstrass-
//! substitution technique as closestThetaOnCircleToLine below) and validated against a brute-force numeric sweep
//! (archive/unused_contact_methods/feng2017_validation_tests.cpp.excerpt). The ellipse helpers here are used only by
//! the archived Feng et al. (2017) cylinder-cylinder method (archive/unused_contact_methods/).
//!
//! f(theta) = (x0+a*cos theta)^2 + (y0+b*sin theta)^2; F(theta) = f'(theta)/2 =
//! -a*x0*sin(theta) + b*y0*cos(theta) + (b^2-a^2)*sin(theta)*cos(theta) -- a mix of first- and
//! second-harmonic terms. Substituting t=tan(theta/2) and clearing denominators gives a quartic in t
//! (c2 = 0), solved here with Durand-Kerner (no resolvent-cubic branching).
//! Returns all (up to 5, including the t->infinity/theta=pi point the Weierstrass substitution cannot reach)
//! real critical thetas -- every local min AND max on the ellipse, not just the global minimum. The
//! constrained search (closestThetaOnEllipseWithinAxialRange below) needs them: when the unconstrained global
//! minimum falls outside a required range, the CONSTRAINED minimum is generically at another critical point or
//! at the range boundary.
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

//! Global-minimum reduction of allCriticalThetasOnEllipseToOrigin, for callers (and the archived brute-force
//! cross-check) that do not need the constrained search.
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
//! [-zHalfRange, zHalfRange]. This closes the gap Feng et al. (2017) flag in their unconstrained version
//! ("the master cylinder is implicitly assumed infinitely long, and no consideration is given to a possible
//! overlap between the base discs") in both directions: their concern is a false POSITIVE (a close point beyond
//! the host's finite end); fuzz-testing shows the false NEGATIVE matters at least as much -- the unconstrained
//! global minimum falls just outside the valid range while a nearby, still-overlapping point on the same ellipse
//! remains within it. Solving the constrained problem in closed form (rather than rejecting the unconstrained
//! answer) handles both.
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

   // Candidate set 2: the range boundary. zHost(theta) = +-zHalfRange is a first-harmonic equation
   // R*cos(theta-phi) = k (no quartic); if |k/R| > 1 the entire rim lies on one side (always in range, covered
   // by set 1, or never in range) and this adds nothing.
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
//! Re-derived rather than transcribed from the paper's printed equations, which are internally inconsistent here:
//! (1) the circle is introduced as "radius r" but written as x^2+y^2=1 (not r^2), while the resulting quartic's
//! a0=r^2(1-tau^2) coefficient only makes sense for x^2+y^2=r^2; (2) Appendix C's eq. (33), a one-line derivation
//! from its eq. (31)-(32), drops a minus sign (redoing that line does not reproduce eq. (33)). (A third
//! inconsistency, in Appendix A, is noted at closestThetaOnEllipseToOrigin.) This quartic is derived fresh,
//! cross-checked at theta=0, pi and pi/2 (each re-deriving the corresponding coefficient), and validated against a
//! brute-force numeric sweep (archive/unused_contact_methods/feng2017_validation_tests.cpp.excerpt).
//!
//! Derivation: substituting the ellipse's parametrisation into x^2+y^2=R^2 and expanding cos^2(theta),
//! sin^2(theta) via the double-angle identities gives a mix of first- and second-harmonic terms in theta;
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

   // theta=pi is an exact root exactly when c4=0 (the point at theta=pi lies on the circle); the Weierstrass
   // substitution cannot reach it as a finite t, so it is added explicitly.
   if (std::abs(c4) < real_t(1e-12) * std::max(real_t(1), std::abs(c0)))
   {
      thetasOut[countOut++] = math::pi;
      // Remaining roots solve the cubic c3 t^3 + c2 t^2 + c1 t + c0 = 0: the quartic solver with c4 nudged to a
      // small epsilon sends one of the four roots towards a very large |t| (theta -> pi again, already added).
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
//! the plane spanned by (u, v). The band-edge minimisation the paper reduces to a quartic in sin(theta) (eq. 74-77),
//! solved here via a quartic in t = tan(theta/2): expanding |P(theta)-linePos|^2 - ((P(theta)-linePos).lineDir)^2
//! gives const + M'cos(theta) + N'sin(theta) - R^2*Amag^2*cos^2(theta-psi) (a mix of first- and
//! second-harmonic terms), whose stationary points satisfy R^2*Amag^2*sin(2*theta-2*psi) =
//! M'*sin(theta) - N'*cos(theta). Substituting t=tan(theta/2) (Weierstrass) and clearing
//! denominators turns this into a genuine quartic in t; solveQuarticGeneric finds its roots,
//! and each real root (plus theta=pi, the t->infinity point this substitution can't reach directly)
//! is evaluated against the true squared-distance to pick the global minimum among what are, in
//! general, 1-2 minima and 1-2 maxima of this function on the circle.
// maxIterations: forwarded to the iterative solvers (default 30, see solveQuarticDurandKerner).
// quarticSolver: default Ferrari; QuarticSolverComparisonTest.cpp (synthetic quartics) and
// PaperAnalyticalContactTest.cpp (edge-edge comparison on real geometry) show it ~2.3x faster than Durand-Kerner
// with bit-identical-or-indistinguishable accuracy (see QuarticSolverMethod).
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

   // (The constant term K - A0^2 of f(theta) is omitted: it does not move the minimum, and distSq() below
   // recomputes the true squared distance.)
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
      // Leading coefficient (numerically) vanishes: t=infinity (theta=pi, bestTheta's initial candidate) is a
      // root; the remaining roots solve the cubic-or-lower c3 t^3 + c2 t^2 + c1 t + c0 = 0. The quartic solver
      // with c4 nudged to an epsilon sends one root towards a very large |t| (theta -> pi again), which
      // consider() handles like any other candidate.
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
//! keeps the largest instead of the smallest). Used by containment_detail::circleFarthestFromAxisClipped
//! (ParticleKernels.h) for a cylinder floater's exact farthest rim point from a containment axis.
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
//! a fixed POINT instead of a line. Used by findRimCrossing's useFastReference option below.
inline real_t closestThetaOnCircleToPoint(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                           const Vec3& point)
{
   return closestThetaOnCircleToPointClosedForm(C, u, v, R, point);
}

//! Two arbitrary unit vectors spanning the plane perpendicular to `axis` (paper's u-hat/v-hat, eq.
//! 71-72), built more robustly than the paper's literal "cross with cylinder i's z-axis" (eq. 71, which degenerates
//! when the two axes are parallel) via a standard not-parallel helper vector.
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
//! Eq. 83 (Px^2+Py^2=R_target^2, "radial distance from rim point to target's axis equals target's radius")
//! expands -- the same expansion as in closestThetaOnCircleToLine, set equal to a target value instead of
//! minimised -- into a genuine quartic in t=tan(theta/2) after the Weierstrass substitution. Matches direct
//! numerical brute-force root-scanning to machine precision across thousands of random configurations.
//!
//! Eq. 84-85 (Pz=faceZ, "axial position of rim point equals target's face height") is
//! A1*cos(theta)+A2*sin(theta)=rhs, i.e. Amag*cos(theta-phi)=rhs -- solved directly via
//! phi+-acos(rhs/Amag), the general form of which the paper's own asin formula (eq. 85) is the
//! A1=0 special case (only true for the paper's own axis-tied basis choice, not this file's).
//!
//! Point A can have up to 4 candidate roots (plus theta=pi); point B always has exactly 2
//! (phi+-delta) when a solution exists at all. The paper specifies no tie-break when multiple candidates pass
//! the containment check; the choice here (nearest, angularly, to the rim's own closest point to target's axis) matches
//! the "near first contact, small segments" regime eq. 87's linearity assumption already commits
//! this whole scenario to.
struct RimCrossing
{
   bool contact = false;
   Vec3 A, B;
};

// useFastReference: thetaRef serves only as the tie-break among point A's (up to 4) candidate roots (see
// RimCrossing; a heuristic of this file, not specified by Kodam or Guo). false (default) computes it via
// closestThetaOnCircleToLine -- the exact closest rim point to target's axis LINE, a second quartic solve on top
// of point A's. true uses closestThetaOnCircleToPoint (a single atan2) with target's CENTRE as reference point:
// a different, approximate reference, i.e. an accuracy/speed trade-off (compared in PaperAnalyticalContactTest.cpp's
// fuzz suite; see also CylinderCylinderContactBenchmark.cpp).
// maxIterations: forwarded to the iterative solvers of thetaRef (when not useFastReference) and point A.
// quarticSolver: forwarded to both quartic solves; default Ferrari (see closestThetaOnCircleToLine).
// candidateCountA/candidateCountB: optional diagnostic out-parameters (default nullptr, no effect) counting the
// candidate thetas passing each of point A's/B's containment checks, independent of the tie-break.
// testThetaRefTieBreakNecessity10k found 0/832 real calls with more than one candidate for either point.
// lazyThetaRef (default true): compute thetaRef only after the candidates are collected, and only if either list
// has more than one entry. Collect-then-select keeps the same traversal order and strict "<" tie-break, so ties
// resolve identically (bit-identical results, testEdgeEdgeThreeWayAccuracy10k); false computes it eagerly up front.
inline RimCrossing findRimCrossing(const Vec3& C, const Vec3& u, const Vec3& v, real_t R,
                                    const Vec3& targetPos, const Vec3& eTarget, real_t targetRadius,
                                    real_t targetHalfLength, bool useFastReference = false,
                                    int maxIterations = 30,
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
   if (!lazyThetaRef) ensureThetaRef(); // eager: unconditional

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
   // Kodam's eq. 83; lazyThetaRef only changes when (if ever) thetaRef is computed to pick among its roots.
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

   // Candidates passing containment (independent of thetaRef): at most 5, up to 4 quartic roots plus theta=pi.
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

   // Selection: 0 candidates -> no point A; exactly 1 -> thetaRef not needed (the saving of lazyThetaRef);
   // >1 -> thetaRef breaks the tie (computed here on first need for either point).
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
//! the OTHER cylinder's flat face treated locally as a plane). Validated against GJK/EPA by
//! testCylinderPlaneEdge.
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

//! ALL crossings of a rim circle (centre C, in-plane basis u/v, radius R) with the
//! INFINITE lateral surface of a target cylinder (axis through targetPos along eTarget, radius
//! targetRadius) -- Kodam et al. (2010) eq. 83, the same quartic findRimCrossing solves for its point A.
//! findRimCrossing keeps only ONE crossing (tie-broken by a reference angle, as its edge-edge chord
//! construction needs); the corrected special-contact depth in rimVsFaceOrBand (ClosedFormContactsCylinderCylinder.h)
//! needs every crossing, so this returns
//! all real roots that satisfy the crossing condition (up to 4, plus theta = pi, which the Weierstrass
//! substitution t = tan(theta/2) cannot reach). Kept separate from findRimCrossing so the edge-edge results
//! are unaffected. Returns the number of thetas written.
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

//! Kodam et al. (2010) eq. 40-44: exact face-band contact point. Reusing the face-face/band-band-parallel point
//! construction (faceFaceStyleContactPoint) with 2*cj.radius substituted for the radial overlap gives a materially
//! different point in general: for PaperAnalyticalContactTest.cpp's testCylinderCylinderFaceBand geometry it
//! returns cylinder i's own axis (0,0) only because the test uses equal radii (cancelling eq. 43's (Ri-Rj) term),
//! whereas the exact construction below returns (0, 0.45), offset toward cylinder j's centre as expected.
//!
//! `relI` is cylinder j's centre in cylinder i's body frame (so cylinder i's near face lies at
//! z = sign(relI[2])*halfLi); `axisJInI` is cylinder j's axis expressed in the SAME frame; `halfLj`
//! is cylinder j's half-length; `Ri` is cylinder i's radius; `overlap` is the already-computed axial
//! overlap (eq. 38). Returns the contact point in cylinder i's body frame (caller transforms to
//! world via ci.pos + ci.R * result).
//!
//! Eq. 40-41 construct line AB by translating cylinder j's centreline (purely along cylinder i's z-axis, NOT
//! radially -- eq. 40-41 have no cj.radius term) onto cylinder i's near face plane. Eq. 42/43 place the contact
//! point at AB's midpoint when AB lies fully inside i's face disc (radius Ri); eq. 44 places it at the midpoint of
//! B and the disc-boundary crossing point D when AB pokes out through one end, a line-circle intersection the
//! paper does not spell out ("the solution is straightforward and will not be presented here"). Solved here via
//! the finite segment/circle intersection quadratic, which reduces exactly to eq. 43 (both endpoints inside ->
//! tLo=-halfLj, tHi=halfLj -> AB's centre) and eq. 44 (one endpoint inside -> tLo or tHi = the surviving endpoint,
//! the other = the crossing point), and also covers AB poking out through BOTH ends (a long, narrow cylinder j
//! straddling a much smaller face i, as in the geometry above), which neither eq. 42 nor eq. 44 names, by
//! midpointing the two crossing points.
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
         // Should not occur after the caller's containment gate (eq. 34/36/37); on floating-point edge cases,
         // use the point on the finite segment nearest cylinder i's axis rather than an undefined chord.
         const real_t tNear = std::clamp(-m / k, -halfLj, halfLj);
         tLo = tHi = tNear;
      }
   }

   const real_t tMid = real_t(0.5) * (tLo + tHi);
   const real_t pz   = faceZ - math::sign(relI[2]) * real_t(0.5) * overlap;
   return Vec3(cx + tMid * ux, cy + tMid * uy, pz);
}

//! Closest squared distance between two finite line segments P1->Q1 and P2->Q2: wrapper around
//! detail::closestPointsBetweenSegments (ClosedFormContactsBoxPoint.h) for callers that need only the distance,
//! e.g. the cylinder-cylinder capsule bound.
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
//! direction -- dirWorld MUST be unit length. Standard OBB-projection formula (as used by full SAT for each of
//! its 15 candidate axes): the sum of each local half-extent weighted by |local axis . dirWorld|.
inline real_t boxSupportAlongDirection(const Vec3& halfExtents, const Mat3& R, const Vec3& dirWorld)
{
   return halfExtents[0] * std::abs(dirWorld * Vec3(R(0, 0), R(1, 0), R(2, 0))) +
          halfExtents[1] * std::abs(dirWorld * Vec3(R(0, 1), R(1, 1), R(2, 1))) +
          halfExtents[2] * std::abs(dirWorld * Vec3(R(0, 2), R(1, 2), R(2, 2)));
}

//! CONSERVATIVE box-box separation test using only ONE candidate separating axis (the line between the two
//! centres) instead of full SAT's 15: a handful of dot products, but incomplete -- it can fail to reject a pair
//! separated along some other axis (e.g. an edge-edge axis), yet by the separating-axis theorem it never rejects
//! an overlapping pair: if the projections onto this legitimate candidate axis do not overlap, the boxes do not
//! overlap. A cheap pre-filter ahead of the exact boxBoxSAT (as closestDistanceSqBetweenSegments is for
//! cylinders), not a replacement for it.
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
//! R^{PG} -- mathematically equivalent (R^{PG} is only used to take components along `planeNormal` and within
//! the plane) and less prone to a transposed-matrix mistake.
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
   // detail::rimPointTowardPlane (see its comment for the derivation).
   const Vec3   E  = detail::rimPointTowardPlane(cyl, planePos, nHat);
   const real_t zE = (E - planePos) * nHat;

   if (zE * zCyl < real_t(0)) // eq. 18: opposite sides of the plane => contact
   {
      result.contact  = true;
      result.overlap  = std::abs(zE);                                          // eq. 19
      result.normal   = math::sign(zCyl) * nHat;                               // eq. 20
      // eq. 21 is z_E^P - sign(z_E^P)*0.5*delta, i.e. measured from E's OWN height, not from the plane (z=0).
      // Applying the offset to E projected onto the plane would shift the result by an extra zE onto the wrong
      // side of the plane (e.g. zCyl=+2, zE=-1: midpoint at z=-0.5, not +0.5). rimVsFaceOrBand's face-edge branch
      // applies the offset to E the same way.
      result.point    = E - math::sign(zE) * real_t(0.5) * result.overlap * nHat; // eq. 21
      result.scenario = "edge-plane";
   }
   return result;
}

//! Sphere-Cylinder contact -- NOT part of Kodam et al. or Guo et al. (neither paper addresses mixed
//! shape pairs at all); developed for this project, unlike sphere-sphere and sphere-box, which predate
//! it in mesa_pd's own AnalyticCollisionFunctions.h (detectSphereSphereCollision/detectSphereBoxCollision,
//! S. Eibl 2019). Exact and closed-form and, like boxBoxSAT, needs no GJK/EPA fallback: the cylinder's
//! closest-point-to-an-external-point mapping (the "clamp" counterpart to its support mapping, as
//! detectSphereBoxCollision's clamp is for the box) is a COMPLETE, exact solution to the sphere-cylinder
//! minimum-translation-distance problem, not a named-scenario heuristic like Kodam et al.'s cylinder-cylinder
//! dispatch. It covers every configuration, so unlike cylinderCylinderContact it has no "is this result
//! uncertain" gate.
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

   // Arbitrary but well-defined radial direction for the degenerate case (sphere centre exactly ON the axis,
   // rho=0), via detail::perpendicularBasis, instead of a 0/0 direction.
   Vec3 uArb, vArb;
   detail::perpendicularBasis(eZ, uArb, vArb);
   const Vec3 radialDir = (rho > real_t(1e-12)) ? radialVec / rho : uArb;

   if (rho > cyl.radius || std::abs(z) > halfL)
   {
      // Sphere centre outside the (solid, capped) cylinder: the closest surface point is the independent
      // radial/axial clamp, covering band, cap and rim proximity uniformly (two surface families, band and cap).
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
      // Sphere centre embedded inside the cylinder: push out through the surface (band or nearer cap) with the
      // least clearance, the "closest way out" logic of detectSphereBoxCollision for an embedded centre.
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
//! face-face's radial overlap eq. 27). Only cylinder i's radius (Ri) and cylinder j's length (eq. 32's L_j, passed
//! as the parameter `Li`; the caller passes cj.length) are needed -- deltaRadial and deltaAxial already fold in the
//! other cylinder's radius/length, so it is not a separate parameter.
inline Vec3 faceFaceStyleContactPoint(const Vec3& relBodyI, real_t Ri, real_t Li, real_t deltaAxial,
                                       real_t deltaRadial)
{
   const real_t x = relBodyI[0], y = relBodyI[1], z = relBodyI[2];
   const real_t rho = std::sqrt(x * x + y * y);
   const Vec3   rHat = (rho > real_t(1e-10)) ? Vec3(x, y, real_t(0)) / rho : Vec3(real_t(1), real_t(0), real_t(0));
   const Vec3   xyContact = (Ri - real_t(0.5) * deltaRadial) * rHat;                        // eq. 31
   const real_t zContact  = z - math::sign(z) * (real_t(0.5) * Li - real_t(0.5) * deltaAxial); // eq. 32
   // body-i frame; caller transforms to world
   return Vec3(xyContact[0], xyContact[1], zContact);
}

} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla
