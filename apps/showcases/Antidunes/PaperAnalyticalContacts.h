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
//! USED BY THE SIMULATION, opt-in: ParticleKernels.h's ComputeContacts includes this header and calls
//! cylinderCylinderContact/cylinderCylinderResultIsUncertain directly, gated behind the runtime prm
//! flag ModelParameters.useAnalyticalContactDetection (default false -- GJK/EPA-only unless a scenario
//! explicitly turns this on). When on, a result cylinderCylinderResultIsUncertain flags still falls
//! back to GJK/EPA UNLESS ModelParameters.usePureAnalyticalCylinder is ALSO set, which trusts the flagged
//! analytical answer with no safety net -- see ComputeContacts::operator()'s own comment for the exact
//! tiering. PaperAnalyticalContactTest.cpp additionally compares these closed forms against
//! SolidCylinderContactDetection's GJK/EPA results directly, standalone from the production dispatch.
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
#include <vector>

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

namespace detail
{
//! Closest points between two finite line segments P1->Q1 and P2->Q2 -- the standard closed-form
//! algorithm (Ericson, "Real-Time Collision Detection" 2004, Section 5.1.9, ClosestPtSegmentSegment),
//! no iteration. Defined here (before boxBoxSAT, which needs it for its own useRoundedCorners option)
//! rather than down in the cylinder-cylinder section further below, where the squared-distance-only
//! wrapper closestDistanceSqBetweenSegments (used by the cylinder-cylinder capsule bound) still lives.
inline void closestPointsBetweenSegments(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2,
                                          Vec3& outC1, Vec3& outC2)
{
   const Vec3   d1  = q1 - p1;
   const Vec3   d2  = q2 - p2;
   const Vec3   r   = p1 - p2;
   const real_t a   = d1 * d1;
   const real_t e   = d2 * d2;
   const real_t f   = d2 * r;
   const real_t eps = real_t(1e-14);

   real_t s, t;
   if (a <= eps && e <= eps) { s = t = real_t(0); }
   else if (a <= eps)
   {
      s = real_t(0);
      t = std::clamp(f / e, real_t(0), real_t(1));
   }
   else
   {
      const real_t c = d1 * r;
      if (e <= eps)
      {
         t = real_t(0);
         s = std::clamp(-c / a, real_t(0), real_t(1));
      }
      else
      {
         const real_t b     = d1 * d2;
         const real_t denom = a * e - b * b; // always >= 0 (Cauchy-Schwarz); 0 exactly when parallel
         s = (denom > eps) ? std::clamp((b * f - c * e) / denom, real_t(0), real_t(1)) : real_t(0);
         t = (b * s + f) / e;

         if (t < real_t(0))
         {
            t = real_t(0);
            s = std::clamp(-c / a, real_t(0), real_t(1));
         }
         else if (t > real_t(1))
         {
            t = real_t(1);
            s = std::clamp((b - c) / a, real_t(0), real_t(1));
         }
      }
   }

   outC1 = p1 + d1 * s;
   outC2 = p2 + d2 * t;
}
} // namespace detail

//======================================================================================================================
// Box x Box -- Gottschalk, Lin & Manocha (1996), Section 5
//======================================================================================================================

//! \param posA,RA,halfA  box A: world position, body->world rotation matrix, half-extents (radii)
//! \param posB,RB,halfB  box B: same, for box B
//! Normal convention matches SolidCylinderContactDetection's pairwise overloads: points from B
//! (2nd argument) toward A (1st argument), e.g. testSphereSphereContact's contactNormal = -dir
//! where the 2nd sphere sits along +dir from the 1st.
//!
//! useRoundedCorners/roundingRadius (both default OFF/0, exact pre-existing behaviour otherwise):
//! mitigation for the corner/edge contact-normal discontinuity documented for polyhedral DEM (Lu,
//! Third & Müller 2015, Chem. Eng. Sci. 127, Fig. 3/Section 2.1.3, citing Fraige et al. 2008a) and
//! reproduced directly against this project's own production parameters by
//! DEMValidationScenarios.h's runCornerNormalDiscontinuity (47.7% spurious kinetic-energy gain with
//! this OFF). SAT's 15 candidate axes give a genuinely discontinuous winning normal near a
//! near-degenerate edge-edge/vertex configuration -- a tiny relative rotation can flip which of two
//! nearly-tied axes has the smaller gap, snapping the reported normal to a different direction between
//! consecutive sub-cycles. Treating each box as the Minkowski sum of the sharp box with a small ball of
//! radius `roundingRadius` (a "rounded box"/spherosimplex, Lu et al.'s Table 2/3) removes this
//! discontinuity at its source for the specific case that actually triggers it -- an EDGE-vs-EDGE
//! winning axis: instead of trusting SAT's raw cross-product axis direction (fixed, box-orientation-
//! relative, and what flips discontinuously), the normal is recomputed as the vector between the
//! ACTUAL closest points on the two finite edge SEGMENTS involved (closestPointsBetweenSegments above)
//! -- a quantity that varies smoothly as either box translates/rotates, by the standard continuity
//! property of the distance between two convex sets away from measure-zero symmetric configurations.
//! Face-axis wins are left untouched (already smooth: a box has only 3 fixed face normals, 90 degrees
//! apart, essentially never near-tied against each other in a way that matters physically). The
//! overlap is inflated by 2*roundingRadius for EVERY winning-axis case (not just edge-edge), matching
//! the actual Minkowski-sum geometry of two rounded boxes each dilated by roundingRadius.
//! maxOverlapCap (default: no cap, exact pre-existing behaviour): clamps the reported overlap to at
//! most this value. See boxBoxMultiPointContacts' own header comment for the full rationale -- SAT's
//! overlap is bounded only by whichever half-extent happens to win as the separating axis, which for
//! two DEEPLY interpenetrating boxes (a jam, or a floater released directly on top of a neighbour) can
//! legitimately reach that entire half-extent, producing a spring force far stronger than any shallow,
//! near-degenerate contact this project's contact model was ever tuned against.
//!
//! previousNormal/hysteresisMargin (default: nullptr/no hysteresis, exact pre-existing behaviour):
//! the OTHER candidate fix for the same axis-flip discontinuity useRoundedCorners targets, requested
//! directly after this session's own diagnosis of the defect -- a "sticky" axis selection instead of
//! geometric rounding. Every candidate axis with gap within hysteresisMargin (relative) of the TRUE
//! minimum gap is treated as "competitive"; among those, the one whose properly-signed normal aligns
//! best with previousNormal is used, even if it is not the strict argmin. This directly targets the
//! chatter itself (two nearly-tied axes flip-flopping every sub-cycle) without inflating the overlap
//! the way useRoundedCorners' Minkowski-sum correction does -- see this project's own
//! multiPointDeepOverlapStability/EdgeCases finding that useRoundedCorners, despite fixing the
//! axis-tie case cleanly, is CATASTROPHICALLY worse than doing nothing on a shape-mismatched deep
//! overlap (+884% on "plank-vs-near-cube-deep") -- a side effect hysteresis has no mechanism to
//! reproduce, since it never changes the reported overlap magnitude, only WHICH already-valid axis is
//! reported. Caller owns the persistence (previousNormal is read-only here) -- see
//! SolidCylinderContactDetection's own boxBoxHysteresisNormal_ map (ParticleKernels.h) for how
//! production keys this per uid-pair across sub-cycles.
inline ContactResult boxBoxSAT(const Vec3& posA, const Mat3& RA, const Vec3& halfA, const Vec3& posB,
                                const Mat3& RB, const Vec3& halfB, bool useRoundedCorners = false,
                                real_t roundingRadius = real_t(0),
                                real_t maxOverlapCap = std::numeric_limits< real_t >::max(),
                                const Vec3* previousNormal = nullptr,
                                real_t hysteresisMargin = real_t(0.15))
{
   Vec3 axesA[3], axesB[3];
   for (uint_t i = 0; i < 3; ++i)
   {
      axesA[i] = Vec3(RA(0, i), RA(1, i), RA(2, i)); // columns of R = box's world-frame axes
      axesB[i] = Vec3(RB(0, i), RB(1, i), RB(2, i));
   }

   const Vec3 T = posB - posA;

   // Every valid (non-separating) candidate axis is kept, not just the running minimum -- hysteresis
   // below needs to compare ALL competitive candidates against previousNormal, not only the argmin.
   struct Candidate
   {
      Vec3   axis;
      real_t gap;
      real_t ra, rb;
      int    edgeIA, edgeJB; // -1,-1 = face-of-A or face-of-B axis, see useRoundedCorners' own comment
   };
   std::vector< Candidate > candidates;
   candidates.reserve(15);

   // Returns false the instant a separating axis is found (boxes proven disjoint); otherwise records
   // this axis as a valid candidate (see Candidate's own comment on why ALL are kept, not just the
   // running minimum).
   auto testAxis = [&](Vec3 L, int ia, int jb) -> bool {
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

      candidates.push_back(Candidate{ L, gap, ra, rb, ia, jb });
      return true;
   };

   for (int i = 0; i < 3; ++i)
      if (!testAxis(axesA[i], -1, -1)) return ContactResult{};
   for (int i = 0; i < 3; ++i)
      if (!testAxis(axesB[i], -1, -1)) return ContactResult{};
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         if (!testAxis(axesA[i] % axesB[j], i, j)) return ContactResult{};

   WALBERLA_ASSERT(!candidates.empty());

   size_t bestIdx = 0;
   for (size_t i = 1; i < candidates.size(); ++i)
      if (candidates[i].gap < candidates[bestIdx].gap) bestIdx = i;
   const real_t minGap = candidates[bestIdx].gap;

   if (previousNormal != nullptr && previousNormal->sqrLength() > real_t(1e-12))
   {
      const Vec3   prevN          = previousNormal->getNormalizedIfNotZero();
      real_t       bestAlignment  = real_t(-2); // below any real cosine, so the loop always picks something
      size_t       hysteresisIdx  = bestIdx;
      const real_t competitiveGap = minGap * (real_t(1) + hysteresisMargin);
      for (size_t i = 0; i < candidates.size(); ++i)
      {
         if (candidates[i].gap > competitiveGap) continue; // not close enough to the true minimum to matter
         const real_t signI = (T * candidates[i].axis >= real_t(0)) ? real_t(1) : real_t(-1);
         const Vec3   nI    = -signI * candidates[i].axis;
         const real_t align = nI * prevN;
         if (align > bestAlignment)
         {
            bestAlignment = align;
            hysteresisIdx = i;
         }
      }
      bestIdx = hysteresisIdx;
   }

   const Candidate& best    = candidates[bestIdx]; // NOT necessarily the global minimum once hysteresis has acted
   const Vec3&      bestAxis = best.axis;
   const real_t     bestRa = best.ra, bestRb = best.rb;
   const int         edgeIA = best.edgeIA, edgeJB = best.edgeJB;

   const real_t sign   = (T * bestAxis >= real_t(0)) ? real_t(1) : real_t(-1);
   Vec3         normal = -sign * bestAxis;

   ContactResult result;
   result.contact  = true;
   result.overlap  = best.gap;
   result.scenario = "SAT";
   // Supplementary contact point (NOT from the paper, see class-level comment): midpoint between
   // the two boxes' extreme points along the winning axis. Only the component along `bestAxis` is
   // geometrically well-defined for a flush face-face/edge-edge configuration.
   result.point = real_t(0.5) * (posA + posB) + real_t(0.5) * sign * (bestRa - bestRb) * bestAxis;

   if (useRoundedCorners && edgeIA >= 0)
   {
      // Edge-vs-edge winning axis: rebuild the two specific finite edge SEGMENTS involved (the edge of
      // A running along axesA[edgeIA], the edge of B running along axesB[edgeJB]) and use the vector
      // between their actual closest points as the normal instead of the raw (discontinuous) SAT axis
      // -- see this function's own header comment for the derivation. Each edge's OTHER two axes are
      // offset by the sign of that axis's projection onto T (A) / -T (B), which selects the specific
      // one of the four axis-parallel edges that faces the other box.
      const uint_t edgeIAu = uint_t(edgeIA), edgeJBu = uint_t(edgeJB);
      Vec3         edgeCenterA = posA;
      for (uint_t o = 0; o < 3; ++o)
      {
         if (o == edgeIAu) continue;
         edgeCenterA += (axesA[o] * T >= real_t(0) ? real_t(1) : real_t(-1)) * halfA[o] * axesA[o];
      }
      Vec3 edgeCenterB = posB;
      for (uint_t p = 0; p < 3; ++p)
      {
         if (p == edgeJBu) continue;
         edgeCenterB += (axesB[p] * T >= real_t(0) ? real_t(-1) : real_t(1)) * halfB[p] * axesB[p];
      }
      const Vec3 halfEdgeA = halfA[edgeIAu] * axesA[edgeIAu];
      const Vec3 halfEdgeB = halfB[edgeJBu] * axesB[edgeJBu];

      Vec3 closestA, closestB;
      detail::closestPointsBetweenSegments(edgeCenterA - halfEdgeA, edgeCenterA + halfEdgeA,
                                            edgeCenterB - halfEdgeB, edgeCenterB + halfEdgeB, closestA,
                                            closestB);
      const Vec3   diff = closestA - closestB; // points from B's edge toward A's edge
      const real_t len  = diff.length();
      if (len > real_t(1e-10))
      {
         normal       = diff / len;
         result.point = real_t(0.5) * (closestA + closestB);
      }
   }

   result.normal = normal;
   if (useRoundedCorners) result.overlap += real_t(2) * roundingRadius;
   result.overlap = std::min(result.overlap, maxOverlapCap);
   return result;
}

//! Closed-form (non-iterative) distance between two solid oriented boxes -- the analytical
//! alternative to running mesa_pd's GJK/EPA with a margin for useSpheroBoxContact (see
//! ModelParameters::useSpheroBoxAnalytical's own doc comment, SimulationParams.h). Ported from
//! David Eberly's GeometricTools library (github.com/davideberly/GeometricTools,
//! GTE/Mathematics/{DistPointCanonicalBox,DistLine3CanonicalBox3,DistSegment3CanonicalBox3,
//! DistPlane3CanonicalBox3,DistRectangle3CanonicalBox3,DistRectangle3OrientedBox3,
//! DistOrientedBox3OrientedBox3}.h) -- the line-vs-box case dispatch (faceCase/doQuery3D/2D/1D/0D
//! below) is documented there as being from Section 10.9.4 "Linear Component to Oriented
//! Bounding Box" of Schneider, P.J. & Eberly, D.H. (2002) "Geometric Tools for Computer
//! Graphics", Morgan Kaufmann. Every function here is a DETERMINISTIC, finite case dispatch --
//! no simplex refinement, no convergence tolerance, no epsilon-loop -- unlike GJK/EPA, which is
//! why this exists as a toggle: it targets the exact same (point,normal,overlap) triplet GJK/EPA
//! produces for a spherosimplex contact, but via closed-form geometry instead of iteration.
namespace eberly_box_distance
{

inline Vec3 clampToCanonicalBox(const Vec3& point, const Vec3& extent)
{
   Vec3 c;
   for (int i = 0; i < 3; ++i) c[i] = std::max(-extent[i], std::min(extent[i], point[i]));
   return c;
}

// ---- Line vs canonical (axis-aligned, origin-centred) box: Schneider & Eberly Section 10.9.4 ----
struct LineBoxResult
{
   real_t parameter   = real_t(0);
   Vec3   closestBox;  // in the SAME (possibly axis-reflected) frame `origin` is passed in
   real_t sqrDistance = real_t(0);
};

// direction[i0] is strictly positive here (the line crosses the x[i0]=+extent[i0] face plane);
// resolves which of the up-to-9 sub-regions of that face's exterior the closest point falls in.
inline void faceCase(int i0, int i1, int i2, Vec3& origin, const Vec3& direction, const Vec3& PmE,
                      const Vec3& extent, LineBoxResult& result)
{
   const Vec3 PpE = origin + extent;
   if (direction[i0] * PpE[i1] >= direction[i1] * PmE[i0])
   {
      if (direction[i0] * PpE[i2] >= direction[i2] * PmE[i0])
      {
         origin[i0] = extent[i0];
         origin[i1] -= direction[i1] * PmE[i0] / direction[i0];
         origin[i2] -= direction[i2] * PmE[i0] / direction[i0];
         result.parameter = -PmE[i0] / direction[i0];
      }
      else
      {
         real_t lenSqr = direction[i0] * direction[i0] + direction[i2] * direction[i2];
         real_t tmp = lenSqr * PpE[i1] - direction[i1] * (direction[i0] * PmE[i0] + direction[i2] * PpE[i2]);
         if (tmp <= real_t(2) * lenSqr * extent[i1])
         {
            real_t t = tmp / lenSqr;
            lenSqr += direction[i1] * direction[i1];
            real_t tmp2 = PpE[i1] - t;
            real_t delta = direction[i0] * PmE[i0] + direction[i1] * tmp2 + direction[i2] * PpE[i2];
            result.parameter = -delta / lenSqr;
            result.sqrDistance += PmE[i0] * PmE[i0] + tmp2 * tmp2 + PpE[i2] * PpE[i2] + delta * result.parameter;
            origin[i0] = extent[i0]; origin[i1] = t - extent[i1]; origin[i2] = -extent[i2];
         }
         else
         {
            lenSqr += direction[i1] * direction[i1];
            real_t delta = direction[i0] * PmE[i0] + direction[i1] * PmE[i1] + direction[i2] * PpE[i2];
            result.parameter = -delta / lenSqr;
            result.sqrDistance += PmE[i0] * PmE[i0] + PmE[i1] * PmE[i1] + PpE[i2] * PpE[i2] + delta * result.parameter;
            origin[i0] = extent[i0]; origin[i1] = extent[i1]; origin[i2] = -extent[i2];
         }
      }
   }
   else
   {
      if (direction[i0] * PpE[i2] >= direction[i2] * PmE[i0])
      {
         real_t lenSqr = direction[i0] * direction[i0] + direction[i1] * direction[i1];
         real_t tmp = lenSqr * PpE[i2] - direction[i2] * (direction[i0] * PmE[i0] + direction[i1] * PpE[i1]);
         if (tmp <= real_t(2) * lenSqr * extent[i2])
         {
            real_t t = tmp / lenSqr;
            lenSqr += direction[i2] * direction[i2];
            real_t tmp2 = PpE[i2] - t;
            real_t delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] + direction[i2] * tmp2;
            result.parameter = -delta / lenSqr;
            result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] + tmp2 * tmp2 + delta * result.parameter;
            origin[i0] = extent[i0]; origin[i1] = -extent[i1]; origin[i2] = t - extent[i2];
         }
         else
         {
            lenSqr += direction[i2] * direction[i2];
            real_t delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] + direction[i2] * PmE[i2];
            result.parameter = -delta / lenSqr;
            result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] + PmE[i2] * PmE[i2] + delta * result.parameter;
            origin[i0] = extent[i0]; origin[i1] = -extent[i1]; origin[i2] = extent[i2];
         }
      }
      else
      {
         real_t lenSqr = direction[i0] * direction[i0] + direction[i2] * direction[i2];
         real_t tmp = lenSqr * PpE[i1] - direction[i1] * (direction[i0] * PmE[i0] + direction[i2] * PpE[i2]);
         if (tmp >= real_t(0))
         {
            if (tmp <= real_t(2) * lenSqr * extent[i1])
            {
               real_t t = tmp / lenSqr;
               lenSqr += direction[i1] * direction[i1];
               real_t tmp2 = PpE[i1] - t;
               real_t delta = direction[i0] * PmE[i0] + direction[i1] * tmp2 + direction[i2] * PpE[i2];
               result.parameter = -delta / lenSqr;
               result.sqrDistance += PmE[i0] * PmE[i0] + tmp2 * tmp2 + PpE[i2] * PpE[i2] + delta * result.parameter;
               origin[i0] = extent[i0]; origin[i1] = t - extent[i1]; origin[i2] = -extent[i2];
            }
            else
            {
               lenSqr += direction[i1] * direction[i1];
               real_t delta = direction[i0] * PmE[i0] + direction[i1] * PmE[i1] + direction[i2] * PpE[i2];
               result.parameter = -delta / lenSqr;
               result.sqrDistance += PmE[i0] * PmE[i0] + PmE[i1] * PmE[i1] + PpE[i2] * PpE[i2] + delta * result.parameter;
               origin[i0] = extent[i0]; origin[i1] = extent[i1]; origin[i2] = -extent[i2];
            }
            return;
         }
         lenSqr = direction[i0] * direction[i0] + direction[i1] * direction[i1];
         tmp = lenSqr * PpE[i2] - direction[i2] * (direction[i0] * PmE[i0] + direction[i1] * PpE[i1]);
         if (tmp >= real_t(0))
         {
            if (tmp <= real_t(2) * lenSqr * extent[i2])
            {
               real_t t = tmp / lenSqr;
               lenSqr += direction[i2] * direction[i2];
               real_t tmp2 = PpE[i2] - t;
               real_t delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] + direction[i2] * tmp2;
               result.parameter = -delta / lenSqr;
               result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] + tmp2 * tmp2 + delta * result.parameter;
               origin[i0] = extent[i0]; origin[i1] = -extent[i1]; origin[i2] = t - extent[i2];
            }
            else
            {
               lenSqr += direction[i2] * direction[i2];
               real_t delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] + direction[i2] * PmE[i2];
               result.parameter = -delta / lenSqr;
               result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] + PmE[i2] * PmE[i2] + delta * result.parameter;
               origin[i0] = extent[i0]; origin[i1] = -extent[i1]; origin[i2] = extent[i2];
            }
            return;
         }
         lenSqr += direction[i2] * direction[i2];
         real_t delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] + direction[i2] * PpE[i2];
         result.parameter = -delta / lenSqr;
         result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] + PpE[i2] * PpE[i2] + delta * result.parameter;
         origin[i0] = extent[i0]; origin[i1] = -extent[i1]; origin[i2] = -extent[i2];
      }
   }
}

inline void lineDoQuery3D(Vec3& origin, const Vec3& direction, const Vec3& extent, LineBoxResult& result)
{
   const Vec3 PmE = origin - extent;
   const real_t prodDxPy = direction[0] * PmE[1], prodDyPx = direction[1] * PmE[0];
   if (prodDyPx >= prodDxPy)
   {
      const real_t prodDzPx = direction[2] * PmE[0], prodDxPz = direction[0] * PmE[2];
      if (prodDzPx >= prodDxPz) faceCase(0, 1, 2, origin, direction, PmE, extent, result);
      else faceCase(2, 0, 1, origin, direction, PmE, extent, result);
   }
   else
   {
      const real_t prodDzPy = direction[2] * PmE[1], prodDyPz = direction[1] * PmE[2];
      if (prodDzPy >= prodDyPz) faceCase(1, 2, 0, origin, direction, PmE, extent, result);
      else faceCase(2, 0, 1, origin, direction, PmE, extent, result);
   }
}

inline void lineDoQuery2D(int i0, int i1, int i2, Vec3& origin, const Vec3& direction,
                           const Vec3& extent, LineBoxResult& result)
{
   const real_t PmE0 = origin[i0] - extent[i0], PmE1 = origin[i1] - extent[i1];
   const real_t prod0 = direction[i1] * PmE0, prod1 = direction[i0] * PmE1;
   if (prod0 >= prod1)
   {
      origin[i0] = extent[i0];
      const real_t PpE1 = origin[i1] + extent[i1];
      const real_t delta = prod0 - direction[i0] * PpE1;
      if (delta >= real_t(0))
      {
         const real_t lenSqr = direction[i0] * direction[i0] + direction[i1] * direction[i1];
         result.sqrDistance += delta * delta / lenSqr;
         origin[i1] = -extent[i1];
         result.parameter = -(direction[i0] * PmE0 + direction[i1] * PpE1) / lenSqr;
      }
      else
      {
         origin[i1] -= prod0 / direction[i0];
         result.parameter = -PmE0 / direction[i0];
      }
   }
   else
   {
      origin[i1] = extent[i1];
      const real_t PpE0 = origin[i0] + extent[i0];
      const real_t delta = prod1 - direction[i1] * PpE0;
      if (delta >= real_t(0))
      {
         const real_t lenSqr = direction[i0] * direction[i0] + direction[i1] * direction[i1];
         result.sqrDistance += delta * delta / lenSqr;
         origin[i0] = -extent[i0];
         result.parameter = -(direction[i0] * PpE0 + direction[i1] * PmE1) / lenSqr;
      }
      else
      {
         origin[i0] -= prod1 / direction[i1];
         result.parameter = -PmE1 / direction[i1];
      }
   }
   if (origin[i2] < -extent[i2])
   {
      const real_t delta = origin[i2] + extent[i2];
      result.sqrDistance += delta * delta;
      origin[i2] = -extent[i2];
   }
   else if (origin[i2] > extent[i2])
   {
      const real_t delta = origin[i2] - extent[i2];
      result.sqrDistance += delta * delta;
      origin[i2] = extent[i2];
   }
}

inline void lineDoQuery1D(int i0, int i1, int i2, Vec3& origin, const Vec3& direction,
                           const Vec3& extent, LineBoxResult& result)
{
   result.parameter = (extent[i0] - origin[i0]) / direction[i0];
   origin[i0] = extent[i0];
   for (int i : { i1, i2 })
   {
      if (origin[i] < -extent[i]) { const real_t d = origin[i] + extent[i]; result.sqrDistance += d * d; origin[i] = -extent[i]; }
      else if (origin[i] > extent[i]) { const real_t d = origin[i] - extent[i]; result.sqrDistance += d * d; origin[i] = extent[i]; }
   }
}

inline void lineDoQuery0D(Vec3& origin, const Vec3& extent, LineBoxResult& result)
{
   for (int i = 0; i < 3; ++i)
   {
      if (origin[i] < -extent[i]) { const real_t d = origin[i] + extent[i]; result.sqrDistance += d * d; origin[i] = -extent[i]; }
      else if (origin[i] > extent[i]) { const real_t d = origin[i] - extent[i]; result.sqrDistance += d * d; origin[i] = extent[i]; }
   }
}

//! Distance between a line (origin + t*direction) and a canonical (origin-centred,
//! axis-aligned) box. Returns the line parameter t and the closest box point; `origin` is
//! reflected into the first octant internally and un-reflected before returning.
inline LineBoxResult lineCanonicalBox(const Vec3& lineOrigin, const Vec3& lineDirection, const Vec3& extent)
{
   LineBoxResult result;
   Vec3 origin = lineOrigin, direction = lineDirection;
   bool reflect[3] = { false, false, false };
   for (int i = 0; i < 3; ++i)
      if (direction[i] < real_t(0)) { origin[i] = -origin[i]; direction[i] = -direction[i]; reflect[i] = true; }

   if (direction[0] > real_t(0))
   {
      if (direction[1] > real_t(0))
      {
         if (direction[2] > real_t(0)) lineDoQuery3D(origin, direction, extent, result);
         else lineDoQuery2D(0, 1, 2, origin, direction, extent, result);
      }
      else
      {
         if (direction[2] > real_t(0)) lineDoQuery2D(0, 2, 1, origin, direction, extent, result);
         else lineDoQuery1D(0, 1, 2, origin, direction, extent, result);
      }
   }
   else
   {
      if (direction[1] > real_t(0))
      {
         if (direction[2] > real_t(0)) lineDoQuery2D(1, 2, 0, origin, direction, extent, result);
         else lineDoQuery1D(1, 0, 2, origin, direction, extent, result);
      }
      else
      {
         if (direction[2] > real_t(0)) lineDoQuery1D(2, 0, 1, origin, direction, extent, result);
         else lineDoQuery0D(origin, extent, result);
      }
   }

   for (int i = 0; i < 3; ++i) if (reflect[i]) origin[i] = -origin[i];
   result.closestBox = origin;
   return result;
}

// ---- Segment vs canonical box: clamp the line result's parameter into [0,1] ----
struct SegmentBoxResult
{
   real_t parameter   = real_t(0);
   Vec3   closestSeg;
   Vec3   closestBox;
   real_t sqrDistance = real_t(0);
};

inline SegmentBoxResult segmentCanonicalBox(const Vec3& p0, const Vec3& p1, const Vec3& extent)
{
   SegmentBoxResult result;
   const Vec3 dir = p1 - p0;
   const auto lb = lineCanonicalBox(p0, dir, extent);
   if (lb.parameter >= real_t(0) && lb.parameter <= real_t(1))
   {
      result.parameter = lb.parameter;
      result.closestSeg = p0 + lb.parameter * dir;
      result.closestBox = lb.closestBox;
      result.sqrDistance = lb.sqrDistance;
   }
   else
   {
      const Vec3& endpoint = (lb.parameter < real_t(0)) ? p0 : p1;
      result.parameter = (lb.parameter < real_t(0)) ? real_t(0) : real_t(1);
      result.closestSeg = endpoint;
      result.closestBox = clampToCanonicalBox(endpoint, extent);
      const Vec3 diff = endpoint - result.closestBox;
      result.sqrDistance = diff * diff;
   }
   return result;
}

// ---- Plane vs canonical box ----
struct PlaneBoxResult
{
   Vec3 closestPlane;
   Vec3 closestBox;
};

inline void planeDoQuery3D(const Vec3& origin, const Vec3& normal, const Vec3& extent, PlaneBoxResult& result)
{
   const real_t dmin = -(normal * (extent + origin));
   if (dmin >= real_t(0)) { result.closestPlane = -extent - dmin * normal; result.closestBox = -extent; return; }
   const real_t dmax = normal * (extent - origin);
   if (dmax <= real_t(0)) { result.closestPlane = extent - dmax * normal; result.closestBox = extent; return; }
   const real_t s = real_t(2) * dmin / (dmin - dmax) - real_t(1);
   result.closestPlane = s * extent;
   result.closestBox = result.closestPlane;
}

inline void planeDoQuery2D(int i0, int i1, int i2, const Vec3& origin, const Vec3& normal,
                            const Vec3& extent, PlaneBoxResult& result)
{
   const real_t dmin = -(normal[i0] * (extent[i0] + origin[i0]) + normal[i1] * (extent[i1] + origin[i1]));
   if (dmin >= real_t(0))
   {
      result.closestPlane[i0] = -extent[i0] - dmin * normal[i0];
      result.closestPlane[i1] = -extent[i1] - dmin * normal[i1];
      result.closestPlane[i2] = extent[i2];
      result.closestBox[i0] = -extent[i0]; result.closestBox[i1] = -extent[i1]; result.closestBox[i2] = extent[i2];
      return;
   }
   const real_t dmax = normal[i0] * (extent[i0] - origin[i0]) + normal[i1] * (extent[i1] - origin[i1]);
   if (dmax <= real_t(0))
   {
      result.closestPlane[i0] = extent[i0] - dmax * normal[i0];
      result.closestPlane[i1] = extent[i1] - dmax * normal[i1];
      result.closestPlane[i2] = extent[i2];
      result.closestBox[i0] = extent[i0]; result.closestBox[i1] = extent[i1]; result.closestBox[i2] = extent[i2];
      return;
   }
   const real_t s = real_t(2) * dmin / (dmin - dmax) - real_t(1);
   result.closestPlane[i0] = s * extent[i0]; result.closestPlane[i1] = s * extent[i1]; result.closestPlane[i2] = extent[i2];
   result.closestBox = result.closestPlane;
}

inline void planeDoQuery1D(int i0, int i1, int i2, const Vec3& origin, const Vec3& extent, PlaneBoxResult& result)
{
   result.closestPlane[i0] = origin[i0]; result.closestPlane[i1] = extent[i1]; result.closestPlane[i2] = extent[i2];
   result.closestBox[i0] = std::max(-extent[i0], std::min(extent[i0], origin[i0]));
   result.closestBox[i1] = extent[i1]; result.closestBox[i2] = extent[i2];
}

//! Distance between a plane (through `pointOnPlane`, unit normal `planeNormal`) and a canonical
//! box; only the two closest points are needed here (the case dispatch mirrors planeDoQuery3D
//! reflected into the first octant of normal-space).
inline PlaneBoxResult planeCanonicalBox(const Vec3& planeNormal, const Vec3& pointOnPlane, const Vec3& extent)
{
   PlaneBoxResult result;
   const real_t d = planeNormal * pointOnPlane;
   Vec3 origin = d * planeNormal;
   Vec3 normal = planeNormal;
   bool reflect[3] = { false, false, false };
   for (int i = 0; i < 3; ++i)
      if (normal[i] < real_t(0)) { origin[i] = -origin[i]; normal[i] = -normal[i]; reflect[i] = true; }

   if (normal[0] > real_t(0))
   {
      if (normal[1] > real_t(0))
      {
         if (normal[2] > real_t(0)) planeDoQuery3D(origin, normal, extent, result);
         else planeDoQuery2D(0, 1, 2, origin, normal, extent, result);
      }
      else
      {
         if (normal[2] > real_t(0)) planeDoQuery2D(0, 2, 1, origin, normal, extent, result);
         else planeDoQuery1D(0, 1, 2, origin, extent, result);
      }
   }
   else
   {
      if (normal[1] > real_t(0))
      {
         if (normal[2] > real_t(0)) planeDoQuery2D(1, 2, 0, origin, normal, extent, result);
         else planeDoQuery1D(1, 2, 0, origin, extent, result);
      }
      else
      {
         if (normal[2] > real_t(0)) planeDoQuery1D(2, 0, 1, origin, extent, result);
         else { result.closestPlane = pointOnPlane; result.closestBox = clampToCanonicalBox(pointOnPlane, extent); }
      }
   }

   for (int i = 0; i < 3; ++i)
      if (reflect[i]) { result.closestPlane[i] = -result.closestPlane[i]; result.closestBox[i] = -result.closestBox[i]; }
   return result;
}

// ---- Rectangle (center, two unit in-plane axes, two extents) vs canonical box ----
struct RectangleBoxResult
{
   Vec3   closestRect;
   Vec3   closestBox;
   real_t sqrDistance = real_t(0);
};

inline RectangleBoxResult rectangleCanonicalBox(const Vec3& rectCenter, const Vec3& axis0, const Vec3& axis1,
                                                 real_t ext0, real_t ext1, const Vec3& extent)
{
   RectangleBoxResult result;
   const Vec3 normal = (axis0 % axis1).getNormalizedIfNotZero();
   const auto pb = planeCanonicalBox(normal, rectCenter, extent);
   const Vec3 delta = pb.closestPlane - rectCenter;
   const real_t s0 = axis0 * delta, s1 = axis1 * delta;

   if (std::fabs(s0) <= ext0 && std::fabs(s1) <= ext1)
   {
      result.closestRect = pb.closestPlane;
      result.closestBox = pb.closestBox;
      const Vec3 diff = result.closestRect - result.closestBox;
      result.sqrDistance = diff * diff;
      return result;
   }

   // The closest plane point falls outside the rectangle -- the true closest point must then
   // lie on one of the rectangle's 4 edges (checked exhaustively; matches the reference
   // implementation's own TODO note that a cheaper clamp is not known to always be correct).
   const Vec3 verts[4] = { rectCenter - ext0 * axis0 - ext1 * axis1, rectCenter + ext0 * axis0 - ext1 * axis1,
                            rectCenter - ext0 * axis0 + ext1 * axis1, rectCenter + ext0 * axis0 + ext1 * axis1 };
   const int edges[4][2] = { { 0, 1 }, { 2, 3 }, { 0, 2 }, { 1, 3 } };
   bool have = false;
   for (const auto& e : edges)
   {
      const auto sb = segmentCanonicalBox(verts[e[0]], verts[e[1]], extent);
      if (!have || sb.sqrDistance < result.sqrDistance)
      {
         have = true;
         result.sqrDistance = sb.sqrDistance;
         result.closestRect = sb.closestSeg;
         result.closestBox = sb.closestBox;
      }
   }
   return result;
}

//! Rectangle (world-space center/axes/extents) vs a general oriented box: transform into the
//! box's own local frame (making it canonical), solve there, transform the result back.
inline RectangleBoxResult rectangleOrientedBox(const Vec3& rectCenter, const Vec3& rectAxis0,
                                                const Vec3& rectAxis1, real_t ext0, real_t ext1,
                                                const Vec3& boxPos, const Mat3& boxR, const Vec3& boxHalf)
{
   const Mat3 boxRT = boxR.getTranspose();
   const Vec3 localCenter = boxRT * (rectCenter - boxPos);
   const Vec3 localAxis0 = boxRT * rectAxis0;
   const Vec3 localAxis1 = boxRT * rectAxis1;

   auto local = rectangleCanonicalBox(localCenter, localAxis0, localAxis1, ext0, ext1, boxHalf);

   RectangleBoxResult result;
   result.closestRect = boxPos + boxR * local.closestRect;
   result.closestBox  = boxPos + boxR * local.closestBox;
   result.sqrDistance = local.sqrDistance;
   return result;
}

struct BoxBoxResult
{
   Vec3   closestA, closestB;
   real_t distance = real_t(0);
};

//! Closed-form minimum distance (and closest point pair) between two solid oriented boxes:
//! decompose each box's boundary into its 6 faces (3 axis pairs x 2 signs), run
//! rectangleOrientedBox for each of those 12 face-vs-other-box queries, keep the global
//! minimum. Complete and exact for two SEPARATED (or exactly touching) convex boxes -- every
//! possible closest-feature pair (vertex-vertex, vertex-edge, vertex-face, edge-edge, edge-face,
//! face-face) is covered because the decomposition spans each box's entire boundary.
inline BoxBoxResult boxBoxDistance(const Vec3& posA, const Mat3& RA, const Vec3& halfA,
                                    const Vec3& posB, const Mat3& RB, const Vec3& halfB)
{
   BoxBoxResult best;
   bool have = false;
   Vec3 axesA[3], axesB[3];
   for (int i = 0; i < 3; ++i)
   {
      axesA[i] = Vec3(RA(0, i), RA(1, i), RA(2, i));
      axesB[i] = Vec3(RB(0, i), RB(1, i), RB(2, i));
   }

   auto consider = [&](const RectangleBoxResult& r, bool swapped) {
      if (!have || r.sqrDistance < best.distance * best.distance)
      {
         have = true;
         best.distance = std::sqrt(r.sqrDistance);
         best.closestA = swapped ? r.closestBox : r.closestRect;
         best.closestB = swapped ? r.closestRect : r.closestBox;
      }
   };

   // Faces of A vs box B.
   for (int i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++)
   {
      const Vec3 scaled = halfA[i2] * axesA[i2];
      auto rPlus  = rectangleOrientedBox(posA + scaled, axesA[i0], axesA[i1], halfA[i0], halfA[i1], posB, RB, halfB);
      auto rMinus = rectangleOrientedBox(posA - scaled, axesA[i0], axesA[i1], halfA[i0], halfA[i1], posB, RB, halfB);
      consider(rPlus, false);
      consider(rMinus, false);
   }
   // Faces of B vs box A.
   for (int i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++)
   {
      const Vec3 scaled = halfB[i2] * axesB[i2];
      auto rPlus  = rectangleOrientedBox(posB + scaled, axesB[i0], axesB[i1], halfB[i0], halfB[i1], posA, RA, halfA);
      auto rMinus = rectangleOrientedBox(posB - scaled, axesB[i0], axesB[i1], halfB[i0], halfB[i1], posA, RA, halfA);
      consider(rPlus, true);
      consider(rMinus, true);
   }
   return best;
}

} // namespace eberly_box_distance

//! Analytical (closed-form, non-GJK) spherosimplex contact for two boxes rounded by `radius`
//! each (combined dilation 2*radius) -- the analytical counterpart to the GJK/EPA-with-margin
//! path useSpheroBoxContact normally uses. See ModelParameters::useSpheroBoxAnalytical's own doc
//! comment (SimulationParams.h). Normal convention matches boxBoxSAT: points from B (2nd
//! argument) toward A (1st argument).
inline ContactResult boxBoxSpheroAnalytical(const Vec3& posA, const Mat3& RA, const Vec3& halfA,
                                             const Vec3& posB, const Mat3& RB, const Vec3& halfB,
                                             real_t radius)
{
   const auto d = eberly_box_distance::boxBoxDistance(posA, RA, halfA, posB, RB, halfB);
   ContactResult result; // default-constructed: contact=false
   // Defensive: the ported case dispatch has measure-zero configurations (exact axis-aligned
   // ties, e.g. two perfectly coincident boxes) where floating-point round-off can push a
   // boundary-equality check to the wrong side and enter a division branch with a near-zero
   // denominator. Rather than propagate a NaN/Inf into the simulation, treat a non-finite
   // result as "detection failed this sub-cycle" (no contact) -- safe because a genuine deep
   // contact will be detected again next sub-cycle once the exact tie has moved.
   if (!std::isfinite(d.distance)) return result;
   const real_t overlap = real_t(2) * radius - d.distance;
   if (overlap <= real_t(0)) return result;

   const Vec3   diff = d.closestA - d.closestB;
   const real_t len  = diff.length();
   if (!std::isfinite(len)) return ContactResult{};
   result.contact  = true;
   result.overlap  = overlap;
   result.normal   = (len > real_t(1e-10)) ? diff / len : Vec3(0, 0, 1);
   result.scenario = "EberlyBoxBoxSphero";
   result.point    = real_t(0.5) * (d.closestA + d.closestB);
   return result;
}

//! One vertex-in-the-other-box overlap, for boxBoxMultiPointContacts below.
struct BoxVertexContact
{
   Vec3   point;
   Vec3   normal; // points OUT of the box the vertex is embedded in, toward the vertex's own box
   real_t overlap;
   // Stable identity for the SAME physical vertex across sub-cycles: 0-7 = one of A's 8 corners
   // (encoded as 4*(sx>0)+2*(sy>0)+1*(sz>0) in A's own body frame), 8-15 = the equivalent for B's
   // corners. A vertex's IDENTITY never changes between sub-cycles (it's always "A's +x+y+z corner"),
   // only whether it happens to be embedded and how deep -- this lets a persistent contact manifold
   // (ParticleKernels.h's useContactManifold) match this sub-cycle's points to last sub-cycle's by
   // exact identity instead of a fragile world-distance threshold.
   int identity = -1;
};

//! Multi-point Box x Box contact: tests all 8 vertices of A against B and all 8 vertices of B against
//! A (the "simple particle geometry" polyhedron test described in Lu et al. 2015, Section 2.2.3 --
//! "it is first assessed whether any vertex from one cube falls inside the other cube and vice versa"),
//! returning one BoxVertexContact per vertex actually found embedded in the other box. Mitigation for
//! the SAME literature-documented defect as boxBoxSAT's useRoundedCorners (see that function's own
//! comment), via the OTHER remedy the review names specifically for this failure mode (Fraige et al.
//! 2008a/2008b): "a multi-point contact concept... multiple simultaneous contacts are modelled. This
//! 'stabilizes' the collision... and effectively prevents a sudden change of the contact force (both
//! direction and amplitude) over the contact region" -- resolving several simultaneously-overlapping
//! vertices each with their own local normal/overlap, rather than the single SAT-chosen deepest point/
//! axis, means no single sub-cycle-to-sub-cycle axis flip can dominate the total force by itself.
//! Deliberately independent of boxBoxSAT/useRoundedCorners -- either mitigation can be enabled alone or
//! together (see ModelParameters.useMultiPointBoxContact/useRoundedBoxCorners, both togglable).
//!
//! maxOverlapCap (default: no cap, exact pre-existing behaviour): clamps EACH vertex's own overlap to
//! at most this value before it is returned. Unlike boxBoxSAT (one overlap number, one force), this
//! function can return several simultaneously-embedded vertices, each contributing its OWN independent
//! spring force in ComputeContacts (ParticleKernels.h) at its OWN world point -- and each is bounded
//! only by min(halfBox) with no shared cap across points, so under deep interpenetration (a real jam,
//! or a floater released directly on top of a neighbour -- see the RELEASE-OVERLAP diagnostic) a single
//! vertex can carry an overlap approaching the box's own half-extent, applied at a point that may sit
//! far from the centroid -- exactly the ingredient for a large, badly-placed spurious torque.
//! DEMValidationScenarios.h's multiPointDeepOverlapStability/Fuzz/EdgeCases scenarios confirm this
//! directly: useMultiPointBoxContact, uncapped, can inject MORE energy than plain boxBoxSAT at
//! identical deep-overlap geometry that a real production jam produces (e.g. the "grazing-plank-stack"
//! edge case, +618% vs baseline's -0.9%) -- precisely the geometry this mitigation was meant to help,
//! not hurt. Capping brings deep-overlap behaviour back in line with what the underlying contact model
//! was actually tuned for (shallow, near-degenerate contact), at the cost of under-representing the
//! true geometric overlap once a vertex is embedded deeper than the cap -- an intentional trade-off:
//! a bounded, physically-plausible force beats an unbounded, geometrically-exact one that can inject
//! spurious energy.
inline std::vector< BoxVertexContact > boxBoxMultiPointContacts(
   const Vec3& posA, const Mat3& RA, const Vec3& halfA, const Vec3& posB, const Mat3& RB, const Vec3& halfB,
   real_t maxOverlapCap = std::numeric_limits< real_t >::max())
{
   std::vector< BoxVertexContact > out;

   auto testVerticesAgainstBox = [&](const Vec3& posV, const Mat3& RV, const Vec3& halfV, const Vec3& posBox,
                                      const Mat3& RBox, const Vec3& halfBox, real_t sign) {
      Vec3 axesV[3], axesBox[3];
      for (uint_t i = 0; i < 3; ++i)
      {
         axesV[i]   = Vec3(RV(0, i), RV(1, i), RV(2, i));
         axesBox[i] = Vec3(RBox(0, i), RBox(1, i), RBox(2, i));
      }
      for (int sx = -1; sx <= 1; sx += 2)
         for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2)
            {
               const Vec3 vertex = posV + real_t(sx) * halfV[0] * axesV[0] + real_t(sy) * halfV[1] * axesV[1] +
                                   real_t(sz) * halfV[2] * axesV[2];
               const Vec3 rel = vertex - posBox;
               real_t     local[3];
               for (uint_t k = 0; k < 3; ++k)
                  local[k] = rel * axesBox[k];
               if (std::abs(local[0]) > halfBox[0] || std::abs(local[1]) > halfBox[1] ||
                   std::abs(local[2]) > halfBox[2])
                  continue; // this vertex is outside the box -- no contact from it

               // Embedded: penetration is the SMALLEST distance to any of the box's 3 pairs of faces
               // (the same "deepest point" logic Lu et al. Section 2.1.3 describes for polygon vertex-
               // in-polygon contacts), and the outward normal is that face's own normal, signed to
               // point away from the box (toward the embedded vertex's own body).
               real_t bestOverlap = std::numeric_limits< real_t >::max();
               uint_t bestAxis    = 0;
               real_t bestSign    = real_t(1);
               for (uint_t k = 0; k < 3; ++k)
               {
                  const real_t overlapK = halfBox[k] - std::abs(local[k]);
                  if (overlapK < bestOverlap)
                  {
                     bestOverlap = overlapK;
                     bestAxis    = k;
                     bestSign    = (local[k] >= real_t(0)) ? real_t(1) : real_t(-1);
                  }
               }
               BoxVertexContact c;
               c.point    = vertex;
               c.normal   = sign * bestSign * axesBox[bestAxis];
               c.overlap  = std::min(bestOverlap, maxOverlapCap);
               c.identity = (sign > real_t(0) ? 0 : 8) + (sx > 0 ? 4 : 0) + (sy > 0 ? 2 : 0) + (sz > 0 ? 1 : 0);
               out.push_back(c);
            }
   };

   // A's vertices embedded in B: normal points from B toward A (this file's convention, see boxBoxSAT's
   // own header comment) -- i.e. away from B, toward the vertex's own box (A), sign=+1 with bestSign
   // already pointing away from B's centre along the box's own local axis.
   testVerticesAgainstBox(posA, RA, halfA, posB, RB, halfB, real_t(1));
   // B's vertices embedded in A: normal must still point from B toward A, i.e. TOWARD A's centre here
   // (away from A would be backwards for this pairing) -- sign=-1 flips bestSign's "away from the
   // embedding box" convention to match.
   testVerticesAgainstBox(posB, RB, halfB, posA, RA, halfA, real_t(-1));

   return out;
}

namespace detail
{
//! Clips a world-space segment (p,q) against a box's own local extent [-half,half] (standard slab/
//! Liang-Barsky parametric clip). Returns false if the segment misses the box entirely; otherwise
//! outTMin/outTMax give the surviving parameter range along p+t*(q-p), t in [0,1].
inline bool clipSegmentToBoxLocal(const Vec3& p, const Vec3& q, const Vec3& pos, const Mat3& R,
                                    const Vec3& half, real_t& outTMin, real_t& outTMax)
{
   Vec3 axes[3];
   for (uint_t i = 0; i < 3; ++i)
      axes[i] = Vec3(R(0, i), R(1, i), R(2, i));
   const Vec3 relP(( p - pos) * axes[0], (p - pos) * axes[1], (p - pos) * axes[2]);
   const Vec3 relQ(( q - pos) * axes[0], (q - pos) * axes[1], (q - pos) * axes[2]);

   real_t tmin = real_t(0), tmax = real_t(1);
   for (uint_t k = 0; k < 3; ++k)
   {
      const real_t d = relQ[k] - relP[k];
      if (std::abs(d) < real_t(1e-14))
      {
         if (relP[k] < -half[k] || relP[k] > half[k]) return false; // parallel to slab, outside it
         continue;
      }
      real_t t1 = (-half[k] - relP[k]) / d;
      real_t t2 = (half[k] - relP[k]) / d;
      if (t1 > t2) std::swap(t1, t2);
      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);
      if (tmin > tmax) return false;
   }
   outTMin = tmin;
   outTMax = tmax;
   return true;
}
} // namespace detail

//! Literal replica of Fraige, Langston & Chen (2008), Powder Technology 186(3):224-240, Section 3.3
//! ("Modelling cube geometry") / Figs. 3-5 -- a DIFFERENT, standalone option from
//! boxBoxMultiPointContacts above (which only implements the paper's vertex-in-box step, sourced
//! secondhand via Lu, Third & Muller 2015's literature review). This function additionally implements
//! the two contact modes the paper describes but boxBoxMultiPointContacts does not (Fig. 5's edge-edge
//! and edge-face cases -- geometries where two boxes overlap with NEITHER box having a vertex embedded
//! in the other, e.g. two elongated boxes crossing like an "X"), plus the paper's own near-corner
//! normal blend (Fig. 4b: "vertex near corner or edge of I" switches the normal from the nearest-face
//! rule to the direction from I's CENTRE to the vertex, specifically to avoid a discontinuity as the
//! nearest-face choice flips right at a corner -- the same normal-discontinuity defect
//! boxBoxSAT::useRoundedCorners and this file's useAxisHysteresis mitigation target via completely
//! different mechanisms).
//!
//! Kept entirely independent of boxBoxMultiPointContacts (own vertex test, own identity numbering) so
//! enabling this new option can never change useMultiPointBoxContact's already-validated behaviour.
//!
//! Contact-point/normal/overlap rules, exactly as the paper describes them:
//!  - Vertex-in-box (Fig. 4): unchanged from boxBoxMultiPointContacts EXCEPT for the near-corner blend
//!    (see above). "Near a corner/edge" is interpreted as: at least 2 of the vertex's 3 per-axis
//!    clearances from the embedding box's own faces are within 5% of that axis's half-extent (a vertex
//!    can only be simultaneously close to 2+ faces near where those faces meet, i.e. an edge or
//!    corner) -- matching the paper's own "within 0.05 of this transformed scale" threshold (Fig. 4's
//!    own caption), generalised from the 2D illustration (where "near corner" is the only such case) to
//!    3D (where 2 faces meet at an edge, 3 at a corner).
//!  - Edge-face (Fig. 5b): for each edge of one box NEITHER of whose own endpoints is already embedded
//!    in the other box (avoiding double-counting what the vertex step above already resolves), clip the
//!    edge segment against the other box's local extent; a nonzero-length surviving sub-segment is a
//!    genuine edge-face contact, point = its midpoint ("centre of edge-length overlap part", the
//!    paper's own phrase), normal/overlap = the nearest-face rule evaluated at that midpoint.
//!  - Edge-edge (Fig. 5a): for every pair of edges (one from each box) where NEITHER endpoint of either
//!    edge is already embedded in the other box, find the closest points between the two finite
//!    segments (closestPointsBetweenSegments); if each closest point actually lies inside the OTHER
//!    box (mutual embedding -- true interpenetration, not just proximity), record a contact at their
//!    midpoint with normal = the two edges' own cross product ("normal to both edges", the paper's own
//!    phrase), signed toward A.
inline std::vector< BoxVertexContact > boxBoxFraigeMultiPointContacts(
   const Vec3& posA, const Mat3& RA, const Vec3& halfA, const Vec3& posB, const Mat3& RB, const Vec3& halfB,
   real_t maxOverlapCap = std::numeric_limits< real_t >::max())
{
   std::vector< BoxVertexContact > out;

   Vec3 axesA[3], axesB[3];
   for (uint_t i = 0; i < 3; ++i)
   {
      axesA[i] = Vec3(RA(0, i), RA(1, i), RA(2, i));
      axesB[i] = Vec3(RB(0, i), RB(1, i), RB(2, i));
   }

   //! World-space endpoints of edge `edgeIdx` (0-11) of a box: edgeIdx = axis*4 + 2*(su>0) + (sv>0),
   //! where `axis` is the axis the edge runs ALONG and (u,v) are the other two axes in cyclic order.
   auto edgeSegment = [](const Vec3& pos, const Vec3 axes[3], const Vec3& half, int edgeIdx, Vec3& outP,
                          Vec3& outQ) {
      const uint_t axis  = uint_t(edgeIdx / 4);
      const int    combo = edgeIdx % 4;
      const uint_t u     = (axis + 1) % 3;
      const uint_t v     = (axis + 2) % 3;
      const real_t su    = (combo & 2) ? real_t(1) : real_t(-1);
      const real_t sv    = (combo & 1) ? real_t(1) : real_t(-1);
      outP = pos + su * half[u] * axes[u] + sv * half[v] * axes[v] - half[axis] * axes[axis];
      outQ = outP + real_t(2) * half[axis] * axes[axis];
   };

   //! Is world point `p` inside box (pos,axes,half)? Same test boxBoxMultiPointContacts' own vertex
   //! step uses -- shared here to decide whether an edge's endpoint is already vertex-resolved.
   auto isInsideBox = [](const Vec3& p, const Vec3& pos, const Vec3 axes[3], const Vec3& half) {
      const Vec3 rel = p - pos;
      for (uint_t k = 0; k < 3; ++k)
         if (std::abs(rel * axes[k]) > half[k]) return false;
      return true;
   };

   // --- Vertex-in-box, both directions (Fig. 4), with the near-corner normal blend added -----------
   auto testVertices = [&](const Vec3& posV, const Mat3& RV, const Vec3& halfV, const Vec3& posBox,
                            const Mat3& RBox, const Vec3& halfBox, real_t sign, int identityBase) {
      Vec3 axesV[3], axesBox[3];
      for (uint_t i = 0; i < 3; ++i)
      {
         axesV[i]   = Vec3(RV(0, i), RV(1, i), RV(2, i));
         axesBox[i] = Vec3(RBox(0, i), RBox(1, i), RBox(2, i));
      }
      for (int sx = -1; sx <= 1; sx += 2)
         for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2)
            {
               const Vec3 vertex = posV + real_t(sx) * halfV[0] * axesV[0] + real_t(sy) * halfV[1] * axesV[1] +
                                   real_t(sz) * halfV[2] * axesV[2];
               const Vec3 rel = vertex - posBox;
               real_t     local[3];
               for (uint_t k = 0; k < 3; ++k)
                  local[k] = rel * axesBox[k];
               if (std::abs(local[0]) > halfBox[0] || std::abs(local[1]) > halfBox[1] ||
                   std::abs(local[2]) > halfBox[2])
                  continue;

               real_t bestOverlap = std::numeric_limits< real_t >::max();
               uint_t bestAxis    = 0;
               real_t bestSign    = real_t(1);
               int    nearCount   = 0;
               for (uint_t k = 0; k < 3; ++k)
               {
                  const real_t clearanceK = halfBox[k] - std::abs(local[k]);
                  if (clearanceK < real_t(0.05) * halfBox[k]) ++nearCount;
                  if (clearanceK < bestOverlap)
                  {
                     bestOverlap = clearanceK;
                     bestAxis    = k;
                     bestSign    = (local[k] >= real_t(0)) ? real_t(1) : real_t(-1);
                  }
               }
               BoxVertexContact c;
               c.point   = vertex;
               c.overlap = std::min(bestOverlap, maxOverlapCap);
               if (nearCount >= 2)
               {
                  // Near a corner/edge of the box (2+ faces simultaneously close, Fig. 4b): blend to
                  // the radial direction from the box's own centre to the vertex instead of the single
                  // nearest-face axis, which is what actually flips discontinuously right at a corner.
                  const real_t relLen = rel.length();
                  c.normal = (relLen > real_t(1e-10)) ? sign * (rel / relLen) : sign * bestSign * axesBox[bestAxis];
               }
               else { c.normal = sign * bestSign * axesBox[bestAxis]; }
               c.identity = identityBase + (sx > 0 ? 4 : 0) + (sy > 0 ? 2 : 0) + (sz > 0 ? 1 : 0);
               out.push_back(c);
            }
   };
   testVertices(posA, RA, halfA, posB, RB, halfB, real_t(1), 0);
   testVertices(posB, RB, halfB, posA, RA, halfA, real_t(-1), 8);

   // --- Edge-face (Fig. 5b): edges of one box whose own endpoints are NOT already vertex-embedded ---
   auto testEdgeFace = [&](const Vec3& posEdgeOwner, const Vec3 axesEdgeOwner[3], const Vec3& halfEdgeOwner,
                            const Vec3& posBox, const Mat3& RBox, const Vec3& halfBox, const Vec3 axesBox[3],
                            real_t sign, int identityBase) {
      for (int e = 0; e < 12; ++e)
      {
         Vec3 p, q;
         edgeSegment(posEdgeOwner, axesEdgeOwner, halfEdgeOwner, e, p, q);
         if (isInsideBox(p, posBox, axesBox, halfBox) || isInsideBox(q, posBox, axesBox, halfBox))
            continue; // already resolved by the vertex step above

         real_t tmin, tmax;
         if (!detail::clipSegmentToBoxLocal(p, q, posBox, RBox, halfBox, tmin, tmax)) continue;
         if (tmax - tmin < real_t(1e-9)) continue; // touches at a single point only -- not a real overlap

         const Vec3 mid = p + real_t(0.5) * (tmin + tmax) * (q - p);
         const Vec3 rel = mid - posBox;
         real_t     local[3];
         for (uint_t k = 0; k < 3; ++k)
            local[k] = rel * axesBox[k];

         real_t bestOverlap = std::numeric_limits< real_t >::max();
         uint_t bestAxis    = 0;
         real_t bestSign    = real_t(1);
         for (uint_t k = 0; k < 3; ++k)
         {
            const real_t clearanceK = halfBox[k] - std::abs(local[k]);
            if (clearanceK < bestOverlap)
            {
               bestOverlap = clearanceK;
               bestAxis    = k;
               bestSign    = (local[k] >= real_t(0)) ? real_t(1) : real_t(-1);
            }
         }
         BoxVertexContact c;
         c.point    = mid;
         c.normal   = sign * bestSign * axesBox[bestAxis];
         c.overlap  = std::min(bestOverlap, maxOverlapCap);
         c.identity = identityBase + e;
         out.push_back(c);
      }
   };
   // A's edges against B: normal points from B toward A, i.e. away from B (sign=+1, matching the
   // vertex step's own convention above).
   testEdgeFace(posA, axesA, halfA, posB, RB, halfB, axesB, real_t(1), 16);
   // B's edges against A: normal must still point from B toward A, i.e. TOWARD A here (sign=-1).
   testEdgeFace(posB, axesB, halfB, posA, RA, halfA, axesA, real_t(-1), 28);

   // --- Edge-edge (Fig. 5a): genuine crossing, neither edge's endpoints already vertex-embedded -----
   for (int ea = 0; ea < 12; ++ea)
   {
      Vec3 pa, qa;
      edgeSegment(posA, axesA, halfA, ea, pa, qa);
      if (isInsideBox(pa, posB, axesB, halfB) || isInsideBox(qa, posB, axesB, halfB)) continue;

      for (int eb = 0; eb < 12; ++eb)
      {
         Vec3 pb, qb;
         edgeSegment(posB, axesB, halfB, eb, pb, qb);
         if (isInsideBox(pb, posA, axesA, halfA) || isInsideBox(qb, posA, axesA, halfA)) continue;

         Vec3 c1, c2;
         detail::closestPointsBetweenSegments(pa, qa, pb, qb, c1, c2);
         // Mutual embedding of the two closest points is what distinguishes a genuine interpenetrating
         // crossing from two edges that merely pass near each other without overlapping.
         if (!isInsideBox(c1, posB, axesB, halfB) || !isInsideBox(c2, posA, axesA, halfA)) continue;

         const Vec3   dirA = (qa - pa);
         const Vec3   dirB = (qb - pb);
         Vec3         n    = dirA % dirB; // cross product -- "normal to both edges" (paper's own phrase)
         const real_t nLen = n.length();
         if (nLen < real_t(1e-10)) continue; // near-parallel edges: cross product degenerate, skip
         n /= nLen;
         if (n * (posA - posB) < real_t(0)) n = -n; // sign toward A, this file's own convention

         const Vec3   relB = c2 - posB;
         const Vec3   relA = c1 - posA;
         real_t       clearB = std::numeric_limits< real_t >::max();
         real_t       clearA = std::numeric_limits< real_t >::max();
         for (uint_t k = 0; k < 3; ++k)
         {
            clearB = std::min(clearB, halfB[k] - std::abs(relB * axesB[k]));
            clearA = std::min(clearA, halfA[k] - std::abs(relA * axesA[k]));
         }

         BoxVertexContact c;
         c.point    = real_t(0.5) * (c1 + c2);
         c.normal   = n;
         c.overlap  = std::min(std::min(clearA, clearB), maxOverlapCap);
         c.identity = 40 + ea * 12 + eb;
         out.push_back(c);
      }
   }

   return out;
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
//! Which end-cap circle of the probe the band-edge branch examines. true (DEFAULT since 2026-09-20) = Kodam et al. (2010) Eq. 68-70: S = sign(d_bottom - d_top) with d = the cap centre's distance from
//! the HOST'S AXIS measured in the host's xy-plane, i.e. the cap circle closer to the host's axis. false = the previous behaviour of this file: the cap whose centre is closer to the host's CENTRE
//! (3D distance). The old rule was a deviation from the paper: for a rod pair that is nested it picked the circle that merely grazes the band (overlap ~0.003 LU, wrong normal) instead of the one
//! that dips into it (job 4251049, rods ...486/...489, steps 3484-3485.6; see THESIS_PLAN.md). The switch is kept so tests can reproduce the old behaviour.
inline bool kodamCapSelectionByAxisDistance = true;

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
         // [Claude, 2026-09-29] contact point sign fix: E is the probe's rim point that has penetrated
         // deepest into the host (|relE_z| < halfLHost), deltaA its distance below the host face at
         // z = sign(relE_z)*halfLHost, and result.normal (host->probe here) points from E towards that
         // face. The overlap along the normal therefore spans [E, E + deltaA*normal], whose midpoint is
         // E + 0.5*deltaA*normal. The previous "E - ..." placed the point 0.5*deltaA deeper inside the
         // host, outside the probe (i.e. outside the overlap region) -- inconsistent with band-edge
         // below, which already used P + 0.5*deltaRp*nHat. Checked by testCylinderCylinderFaceEdgeClean.
         result.point    = E + result.normal * real_t(0.5) * deltaA;
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
   Vec3 capCentre;
   if (kodamCapSelectionByAxisDistance)
   {
      const Vec3 relTop = host.R.getTranspose() * (capTop - host.pos), relBottom = host.R.getTranspose() * (capBottom - host.pos);
      const real_t dTop = std::sqrt(relTop[0] * relTop[0] + relTop[1] * relTop[1]), dBottom = std::sqrt(relBottom[0] * relBottom[0] + relBottom[1] * relBottom[1]);
      capCentre = (dTop < dBottom) ? capTop : capBottom;
   }
   else
      capCentre = ((capTop - host.pos).length() < (capBottom - host.pos).length()) ? capTop : capBottom;
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
   // [Claude, 2026-09-29] Guo's criterion compared delta_a AT P with delta_r. P (closest to the host axis) is the
   // probe's deepest point RADIALLY, but not axially, so delta_a(P) underestimates the axial penetration and the
   // comparison picks the axial ("special") contact too often: in the fuzz suite, most special-contact
   // disagreements with GJK/EPA had GJK normals ~90 degrees from the host axis, i.e. the true escape direction
   // was radial. The comparison now uses the corrected axial depth computed below (deepest rim point inside the
   // host); if it is not smaller than delta_r, the contact falls through to band-edge (radial depth at P).
   if (sameSide && probeCentreBeyondFace && deltaAp < deltaRp)
   {
      // [Claude, 2026-09-29] CORRECTED special-contact depth. Guo et al. (2012, Table 1 / Fig. 8) place this
      // contact "at the point P with overlap delta_a", P being the rim point closest to the host AXIS. For a
      // tilted probe, P is in general NOT the probe's deepest point inside the host: the rim is a tilted
      // circle whose height varies along it, and its deepest point inside the host lies where the rim leaves
      // the host's face disc, i.e. where it crosses the host's lateral surface (or at the rim's lowest point
      // itself, if that lies inside the disc). On testCylinderCylinderSpecialContact's geometry Guo's value is
      // 0.181 while GJK/EPA and an independent brute-force MTD give 0.610. Candidates: (a) the rim's own
      // deepest point towards the host face, if inside the host; (b) every crossing of the rim with the
      // host's lateral surface inside the host's length (Kodam eq. 83 quartic). The deepest candidate
      // defines the overlap; normal stays axial (Guo). P itself is kept as a last resort.
      const real_t s        = math::sign(relP[2]); // side of the host face the probe comes from
      const Vec3   eHostW   = eHost;
      real_t       bestSZ   = relP[2] * s;          // s*z of the deepest candidate so far (smaller = deeper)
      Vec3         bestPt   = P;
      auto consider = [&](const Vec3& Q) {
         const Vec3   relQ = host.R.getTranspose() * (Q - host.pos);
         const real_t rho2 = relQ[0] * relQ[0] + relQ[1] * relQ[1];
         if (std::abs(relQ[2]) >= halfLHost) return;                                   // outside host length
         if (rho2 > host.radius * host.radius * (real_t(1) + real_t(1e-9))) return;    // outside host disc
         if (relQ[2] * s < bestSZ) { bestSZ = relQ[2] * s; bestPt = Q; }
      };
      // (a) lowest rim point towards the host face: minimises s*(Q - host.pos).eHost along the rim.
      {
         const real_t A1 = uP * eHostW, A2 = vP * eHostW;
         const real_t thetaLow = std::atan2(-s * A2, -s * A1);
         consider(capCentre + probe.radius * std::cos(thetaLow) * uP + probe.radius * std::sin(thetaLow) * vP);
      }
      // (b) crossings of the rim with the host's lateral surface.
      real_t thetas[5];
      const int nCross = detail::allRimCrossingsWithLateralSurface(capCentre, uP, vP, probe.radius, host.pos, eHostW,
                                                                   host.radius, thetas);
      for (int k = 0; k < nCross; ++k)
         consider(capCentre + probe.radius * std::cos(thetas[k]) * uP + probe.radius * std::sin(thetas[k]) * vP);

      const real_t deltaCorrected = halfLHost - bestSZ;
      if (deltaCorrected < deltaRp)
      {
      result.contact  = true;
      result.overlap  = deltaCorrected;
      result.normal   = host.R * Vec3(real_t(0), real_t(0), s);
      result.point    = bestPt + result.normal * real_t(0.5) * deltaCorrected;
      result.scenario = "special-contact";
      return result;
      }
      // else: axial penetration is the larger one -> band-edge below (radial depth at P).
   }
   // [Claude, 2026-09-29] the original Guo et al. special-contact result (overlap = deltaAp at P, point
   // P + 0.5*deltaAp*n) was replaced by the corrected depth above.

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
// useFastEdgeEdgeReference: forwarded to findRimCrossing's own useFastReference (see that function's
// own comment) for the edge-edge branch's two rim-crossing calls -- default false, i.e. UNCHANGED
// behaviour for every existing caller. Has no effect on any scenario other than edge-edge.
// maxIterations: forwarded ONLY to the edge-edge branch's two findRimCrossing calls below (NOT to
// band-edge/special-contact/face-edge's own quartic solves inside rimVsFaceOrBand, which still use
// that function's own unparameterized default of 30 regardless of what's passed here). Default 30
// (unchanged behaviour); exposed so PaperAnalyticalContactTest.cpp can run a genuine 3-way comparison
// on edge-edge specifically (original 60 vs current-default 30-with-old-reference vs
// 30-with-new-reference) on IDENTICAL trials in one process, instead of separate builds/runs.
// quarticSolver: this PARAMETER is forwarded ONLY to the edge-edge branch's two findRimCrossing calls
// (same scoping as maxIterations above -- rimVsFaceOrBand's own quartic solves never see whatever is
// passed HERE). Default changed to Ferrari (from DurandKerner) this session -- since that is also
// closestThetaOnCircleToLine's own new default, rimVsFaceOrBand's unparameterized calls (band-edge,
// special-contact) pick up Ferrari too, transitively, not via this parameter's forwarding but via that
// function's own default. Net effect: EVERY quartic solve in this file's cylinder-cylinder dispatch now
// defaults to Ferrari, validated across all scenarios (not just edge-edge) by
// testCylinderCylinderFullComparisonAndHybridFeasibility, run at 10000 trials for this change. See
// QuarticSolverComparisonTest.cpp for the standalone time/accuracy comparison (synthetic quartics,
// known ground-truth roots) that motivated exposing this here, and
// PaperAnalyticalContactTest.cpp's edge-edge N-way comparison for the same check against real
// cylinder-cylinder geometry and GJK/EPA.
// lazyThetaRef: forwarded to the edge-edge branch's two findRimCrossing calls (same scoping as
// maxIterations/quarticSolver above -- rimVsFaceOrBand is unaffected). Default changed to true this
// session (same as findRimCrossing's own default -- see that function's own comment for the empirical
// finding, 0/832 real calls ever needed thetaRef at all, and testLazyThetaRefAccuracy10k's confirmation
// of bit-identical results). false (the original eager behaviour) remains fully available.
inline ContactResult cylinderCylinderContact(const CylinderState& ci, const CylinderState& cj,
                                              bool useFastEdgeEdgeReference = false, int maxIterations = 30,
                                              detail::QuarticSolverMethod quarticSolver =
                                                 detail::QuarticSolverMethod::Ferrari,
                                              bool lazyThetaRef = true)
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
            // Contact point: eq. 40-44, via detail::faceBandContactPoint (see its own comment for
            // the derivation, including the line-circle intersection eq. 44 declines to spell out).
            const Vec3 axisJInI = ci.R.getTranspose() * eJ;
            result.point    = ci.pos + ci.R * detail::faceBandContactPoint(relI, axisJInI, halfLj,
                                                                            ci.radius, halfLi, result.overlap);
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
         //
         // NEAR-PARALLEL axes (1 - (e_i.e_j)^2 = sin^2(theta) < 1e-10, i.e. theta < 1e-5 rad -- ordinary
         // floating-point rotation noise of a body that was placed axis-aligned): the two infinite
         // centrelines only approach each other about (separation / theta) ~ 1e5-1e6 LU away along the
         // axes, so eq. 63's |l*| < L/2 conditions REJECT the band-band contact -- and that is the
         // right answer, because any near-parallel pair that really touches band-to-band or face-to-face
         // is caught by the parallel-axes branch at the top of this function (|zIJ| < (L_i+L_j)/2 AND
         // rhoIJ < R_i+R_j). This branch used to set l* = 0 in that regime ("guard against division by
         // ~0"), which made BOTH l* conditions vacuously true while r_hat = (e_i x e_j)/|e_i x e_j| is a
         // pure noise direction (any unit vector perpendicular to both axes, e.g. perpendicular to the
         // actual separation => d = |x.r_hat| = 0): two rods separated by a clear gap were reported as
         // overlapping by R_i+R_j (a 6 LU false contact -> cylinder hopper [RELEASE-OVERLAP] abort, job
         // 4248725). Kodam's equations are fine; only that guard deviated from them. Now the branch is
         // simply not applicable in that regime (see CylinderNearParallelSkewedTest.cpp, which checks
         // against a literal implementation of eq. 55-66 with the exact l*).
         const Vec3   w0     = ci.pos - cj.pos;
         const real_t a_     = eI * eI, b_ = eI * eJ, c_ = eJ * eJ;
         const real_t d_     = eI * w0, e_ = eJ * w0;
         const real_t denom2 = a_ * c_ - b_ * b_;
         if (denom2 > real_t(1e-10))
         {
            const real_t liW = (b_ * e_ - c_ * d_) / denom2;
            const real_t ljW = (a_ * e_ - b_ * d_) / denom2;

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

      const detail::RimCrossing jVsI =
         detail::findRimCrossing(capJ, uJ, vJ, cj.radius, ci.pos, eI, ci.radius, halfLi, useFastEdgeEdgeReference,
                                  maxIterations, quarticSolver, nullptr, nullptr, lazyThetaRef);
      if (jVsI.contact)
      {
         const detail::RimCrossing iVsJ =
            detail::findRimCrossing(capI, uI, vI, ci.radius, cj.pos, eJ, cj.radius, halfLj, useFastEdgeEdgeReference,
                                     maxIterations, quarticSolver, nullptr, nullptr, lazyThetaRef);
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
                  // Pi is on rim i, positioned where rim i crosses INTO cylinder j (so physically near
                  // cj); Pj is on rim j, positioned where rim j crosses INTO cylinder i (so physically
                  // near ci). This file's convention is "normal points toward the 1st argument (ci)"
                  // (see the face-face branch's own comment) -- i.e. from the cj-side point toward the
                  // ci-side point, which is Pj-Pi, NOT Pi-Pj. Trusted directly (no Tier 2 substitution --
                  // that gate was removed, see cylinderCylinderResultIsUncertain's own comment): a fuzz
                  // run against GJK/EPA confirms this direction agrees to a few degrees on trials where
                  // the reported magnitude is also correct (see PaperAnalyticalContactTest.cpp's
                  // edge-edge normal-angle breakdown, PASSING rows). The residual edge-edge disagreement
                  // documented elsewhere in this file is a wrong-MAGNITUDE issue (the curvature-coupling
                  // feature gap), not a sign/direction error in this formula.
                  result.normal   = (Pj - Pi).getNormalizedOrZero();
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
//! TIER 1 (dispatch found nothing, but the bodies are close enough that a real contact is
//! plausible): a cheap bounding-sum gate -- centre-to-centre distance vs. the sum of both cylinders'
//! "reach" (radius + half-length each, i.e. the distance at which even a fully-extended,
//! perfectly-aimed pair could just barely touch tip-to-tip). This is the SAME style of coarse
//! bounding check DEM codes already use before ever calling fine contact detection (cf. Guo et al.
//! 2012 Section 2.1's circumscribing-sphere/cube coarse detection) -- deliberately loose, so it will
//! flag plenty of genuinely-separated pairs too; that's fine here since triggering it only means
//! "double-check with GJK/EPA", not "report a contact". A 10000-trial fuzz run found this catches
//! 99.2% of the "no scenario matched at all" gap correctly (the residual 0.8% are real contacts this
//! bound still misses -- see PaperAnalyticalContactTest.cpp's dispatch-gap characterisation -- but
//! those are always resolved correctly downstream anyway, since GJK/EPA IS the fallback).
//!
//! NOTE: Kodam et al. 2010's own ten equations (Part I) only ever test for POSITIVE contact -- none
//! of them is a "definitely separated" fast path. Their own Section 4 conclusion recommends adding
//! "an intermediate detection scheme such as a bounding sphere or rectangular prism" ON TOP of their
//! algorithms for efficiency, without specifying one -- this Tier 1 gate, and ParticleKernels.h's own
//! separate, TIGHTER per-pair reject (circumscribing-sphere sum, sqrt(radius^2+halfLength^2) each,
//! run BEFORE cylinderCylinderContact is even called) are this codebase's own answer to that
//! recommendation, not anything derived from the paper itself. Because that outer reject is strictly
//! tighter than this Tier 1 bound (r+L/2 per side is always >= sqrt(r^2+(L/2)^2) per side), any pair
//! that ever reaches this Tier 1 check has, by construction, already passed the outer guard -- so this
//! check can never actually reject anything in production; every "no scenario matched" case
//! unconditionally falls back to GJK/EPA. Most of those calls are expected to be genuine non-contacts
//! (pairs within possible-touch range but not actually touching -- normal DEM broad-phase population,
//! not a sign anything is wrong), not real contacts the dispatch is silently missing -- see
//! ParticleKernels.h's ContactDetectionStats::cylFallbackTier1GjkFoundContact/GjkFoundNoContact
//! counters, added specifically to measure that split from real runs instead of assuming it.
//!
//! special-contact: unconditional distrust. This scenario's own overlap formula (Guo et al. 2012,
//! Table 1: the axial distance from a single rim point to the host's face, delta_a) is Guo's fix for
//! this configuration's DETECTION -- ensuring it isn't silently dropped or misclassified as
//! band-edge -- not a claim about magnitude accuracy, and Guo's own paper never validates it as one.
//! Verified independently three ways for the hand-picked test geometry
//! (testCylinderCylinderSpecialContact): Kodam/Guo reports 0.181, GJK/EPA reports ~0.610, and a
//! THIRD, fully independent brute-force minimum-translation-distance search directly over each
//! cylinder's own support function (no dependency on either implementation -- see
//! testCylinderCylinderSpecialContactBruteForceMTD) finds 0.6095, matching GJK/EPA almost exactly and
//! confirming Guo's own formula genuinely underestimates by >3x. A 10000-trial fuzz run found this
//! scenario wrong 93% of the time it fires, the single largest concentrated error source in the raw
//! analytical dispatch -- and it is rare enough (under 2% of all real contacts) that routing it to
//! GJK/EPA unconditionally costs almost nothing.
//!
//! EVERY OTHER SCENARIO IS TRUSTED DIRECTLY, WITH NO FURTHER CHECK. This function used to also run a
//! Tier 2: a circumscribing-sphere magnitude cross-check applied to every non-exempt scenario
//! (face-edge/band-edge/edge-edge, plus band-band-skewed/face-band whenever the fuzz suite found
//! those unreliable too), falling back whenever a looser upper bound suggested more overlap was
//! geometrically possible than what was reported. It was REMOVED (see PaperAnalyticalContactTest.cpp,
//! testCylinderCylinderFullComparisonAndHybridFeasibility's history) once a dedicated 5000/10000-trial
//! master comparison showed:
//!   (a) it pushed the GJK/EPA fallback rate to ~90% -- for most trials, the closed-form dispatch was
//!       being computed and then thrown away, defeating most of the point of having one;
//!   (b) it is structurally ONE-SIDED: "is MORE overlap possible than what was reported?" reliably
//!       catches UNDERESTIMATES (band-edge/edge-edge/face-edge/special-contact's typical failure
//!       mode), but has no mechanism at all for OVERESTIMATES. A dedicated deep/degenerate-overlap
//!       stress test (testCylinderCylinderLargeOverlapDegenerateStress, deliberately pushed to 30%-
//!       150% of the smaller radius -- well past the 15% "realistic" cap every other fuzz test in
//!       this file uses) found a case where the FULL hybrid, Tier 2 included, still returned an
//!       overlap ~4.7x too large (372.6% relative error) because the bad value itself already
//!       exceeded Tier 2's own bound, making "is more possible?" trivially false. Tier 2 was
//!       therefore providing a false sense of safety it could not actually deliver in exactly the
//!       regime (large errors) that mattered most.
//! Dropping Tier 2 and keeping only Tier 1 plus special-contact's unconditional routing reduces the
//! GJK/EPA fallback rate to ~3.3% (from ~90%), at the cost of a residual wrong-magnitude rate of
//! ~2.3-2.5% of all genuine contacts -- entirely concentrated in band-edge/edge-edge/face-edge, where
//! the errors are modest (median ~2.7% relative, max ~13.5% absolute in a 10000-trial run) rather
//! than the severe (>90%) errors special-contact alone used to produce. face-face, band-band-parallel
//! and band-band-skewed remain provably exact under their own preconditions regardless (see their own
//! branch comments above), so trusting them directly was never in question.
// useTightTier1: OPTIONAL simplification (this session) -- default false reproduces the ORIGINAL
// Tier 1 bound exactly (centre-to-centre distance vs. radius+halfLength sum per side). true replaces
// it with detail::closestDistanceSqBetweenSegments's capsule bound: the two cylinders' own AXIS
// SEGMENTS vs. the sum of their radii, instead of collapsing each cylinder to a sphere first. Still
// fully conservative (never wrongly declares "definitely separated" for a real contact -- see that
// function's own comment for why), but strictly tighter for any cylinder with length > 0, which is
// every cylinder this file is ever used on. The original bound is ALSO provably redundant in
// production regardless of this flag: ParticleKernels.h's own outer per-pair reject (circumscribing-
// sphere sum) always runs first and is strictly tighter than the original center-distance Tier 1
// bound, so any pair that ever reaches this check has already passed a stricter version of it -- the
// original Tier 1 can never actually reject anything. useTightTier1=true is the first version of this
// check that can. Old bound kept fully intact (not deleted) so the two can be A/B compared -- see
// PaperAnalyticalContactTest.cpp's testTightTier1Accuracy10k.
// skipTier1: OPTIONAL (this session) -- default false runs Tier 1 as normal (either bound above).
// true skips the dispatch-gap check ENTIRELY -- no distance computation of any kind, "no contact" is
// trusted immediately. This is NOT a general-purpose speedup: it is only correct when a sufficiently
// tight Layer 1 pre-filter (ParticleKernels.h's own cheap reject, ModelParameters.useTightLayer1) has
// ALREADY run upstream and rejected every pair that could possibly need this check -- proven via the
// triangle inequality (see useTightTier1's own comment: any pair passing the tight Layer 1 bound also
// passes -- i.e. is never rejected by -- Tier 1's own bound, tight OR original), and confirmed
// empirically (job 4194491: cylTier1CheaplyProvedSeparated=0 at every interval with both flags on).
// Using skipTier1=true WITHOUT a tight Layer 1 filter upstream would silently trust "no contact" for
// pairs Tier 1 would otherwise have caught -- a real correctness risk, not just a missed optimization.
inline bool cylinderCylinderResultIsUncertain(const CylinderState& ci, const CylinderState& cj,
                                               const ContactResult& analytical, bool useTightTier1 = false,
                                               bool skipTier1 = false)
{
   if (!analytical.contact)
   {
      if (skipTier1) return false; // Tier 1 disabled entirely -- see this function's own comment above.
      if (useTightTier1)
      {
         Vec3 p0i, p1i, p0j, p1j;
         detail::cylinderAxisEndpoints(ci, p0i, p1i);
         detail::cylinderAxisEndpoints(cj, p0j, p1j);
         const real_t reach = ci.radius + cj.radius;
         return detail::closestDistanceSqBetweenSegments(p0i, p1i, p0j, p1j) < reach * reach;
      }
      const real_t centreDist = (cj.pos - ci.pos).length();
      const real_t reach = ci.radius + cj.radius + real_t(0.5) * ci.length + real_t(0.5) * cj.length;
      return centreDist < reach; // Tier 1 (original, loose center-distance bound)
   }

   return std::string(analytical.scenario) == "special-contact";
}

//======================================================================================================================
// [Claude, 2026-09-30] Support-function refinement of cylinder-cylinder contacts (thesis Section 5.2).
//
// For two convex, centrally symmetric bodies, the overlap of their projections onto a unit direction L (pointing from
// ci towards cj) is
//    D(L) = r_i(L) + r_j(L) - L.(x_j - x_i),       r_k(L) = (length_k/2) |L.a_k| + R_k sqrt(1 - (L.a_k)^2),
// the bodies are separated iff D(L) < 0 for SOME L (separating-axis theorem), and the penetration depth is min_L D(L),
// attained along the contact normal. Seen as a positively homogeneous function of L in R^3, D has the gradient
//    grad D(L) = S_i(L) - S_j(-L),
// the difference of the two support points (the point of ci farthest along L and the point of cj farthest along -L), and
// grad D(L).L = D(L). The minimisation over the unit sphere is therefore a projected-gradient descent whose gradient
// costs two support-point evaluations -- a few dozen flops, no simplex, no polytope.
//
// Used for two things:
//  (a) refineCylinderContact: the closed-form scenarios that are not exact (edge-edge, band-edge, face-edge, special
//      contact) get their depth and normal from a descent started at their own normal, which is typically close to the
//      optimum. The exact scenarios (face-face, band-band parallel/skewed, face-band) are returned unchanged.
//  (b) the "no scenario matched" case: a descent from a few seed directions. Any direction with D < 0 PROVES
//      separation (early exit). If none is found, the pair is reported as a contact with the refined depth.
// A local descent can end in a local minimum of D; the result is therefore an UPPER bound of the true depth, and in
// case (b) a pair could in principle be reported as touching although a separating direction exists elsewhere. Both are
// measured against GJK/EPA and an independent brute-force reference in PaperAnalyticalContactTest.cpp.
//======================================================================================================================
namespace detail
{
inline real_t cylinderHalfWidth(const CylinderState& c, const Vec3& L)
{
   const real_t ca = std::abs(L * c.axis());
   return real_t(0.5) * c.length * ca + c.radius * std::sqrt(std::max(real_t(0), real_t(1) - ca * ca));
}

//! Point of cylinder c farthest along the unit direction L (on an edge where the support set is a face or band segment,
//! its midpoint along the degenerate direction).
inline Vec3 cylinderSupportPoint(const CylinderState& c, const Vec3& L)
{
   const Vec3   a   = c.axis();
   const real_t d   = L * a;
   Vec3         p   = c.pos;
   if (d > real_t(1e-12)) p += (real_t(0.5) * c.length) * a;
   else if (d < -real_t(1e-12)) p -= (real_t(0.5) * c.length) * a;
   const Vec3   rad = L - d * a;
   const real_t rl  = rad.length();
   if (rl > real_t(1e-12)) p += (c.radius / rl) * rad;
   return p;
}

//! Projection overlap D(L) of ci and cj along the unit direction L (from ci towards cj).
inline real_t cylinderProjectionOverlap(const CylinderState& ci, const CylinderState& cj, const Vec3& L)
{
   return cylinderHalfWidth(ci, L) + cylinderHalfWidth(cj, L) - L * (cj.pos - ci.pos);
}

struct OverlapMinimum
{
   real_t value;     //!< smallest D found
   Vec3   direction; //!< L at which it was found (unit, from ci towards cj)
   int    evaluations;
};

//! Projected-gradient descent of D on the unit sphere from the given seed directions (need not be normalised; zero
//! seeds are skipped). stopIfSeparated: return as soon as some D < 0 is found (existence only).
inline OverlapMinimum minimiseCylinderProjectionOverlap(const CylinderState& ci, const CylinderState& cj,
                                                        const Vec3* seeds, int numSeeds, bool stopIfSeparated,
                                                        int maxIterations = 40, real_t tolerance = real_t(1e-10))
{
   OverlapMinimum best{ std::numeric_limits< real_t >::max(), Vec3(real_t(0)), 0 };
   const real_t   scale = ci.radius + cj.radius + real_t(0.5) * (ci.length + cj.length);
   for (int s = 0; s < numSeeds; ++s)
   {
      const real_t sl = seeds[s].length();
      if (sl < real_t(1e-12)) continue;
      Vec3   L = seeds[s] / sl;
      real_t D = cylinderProjectionOverlap(ci, cj, L);
      ++best.evaluations;
      real_t step = real_t(0.5) / scale; // descent step in 1/length units: L - step * grad, grad has length units
      for (int it = 0; it < maxIterations; ++it)
      {
         if (stopIfSeparated && D < real_t(0)) break;
         const Vec3   g  = cylinderSupportPoint(ci, L) - cylinderSupportPoint(cj, -L);
         const Vec3   gt = g - (g * L) * L; // tangential part
         const real_t gn = gt.length();
         if (gn * scale < tolerance) break;
         bool improved = false;
         for (int bt = 0; bt < 30; ++bt)
         {
            Vec3 Ln = L - step * gt;
            Ln /= Ln.length();
            const real_t Dn = cylinderProjectionOverlap(ci, cj, Ln);
            ++best.evaluations;
            if (Dn < D)
            {
               L = Ln; D = Dn; improved = true;
               step *= real_t(2);
               break;
            }
            step *= real_t(0.5);
         }
         if (!improved) break;
      }
      if (D < best.value) { best.value = D; best.direction = L; }
      if (stopIfSeparated && best.value < real_t(0)) break;
   }
   return best;
}

//! Seed directions (from ci towards cj) for the descent: the closest points of the two axis segments, the centre line,
//! the common perpendicular of the axes, and the two axes, each oriented towards cj. Returns the number written (5).
inline int cylinderSeedDirections(const CylinderState& ci, const CylinderState& cj, Vec3* seeds)
{
   Vec3 p0i, p1i, p0j, p1j, c1, c2;
   cylinderAxisEndpoints(ci, p0i, p1i);
   cylinderAxisEndpoints(cj, p0j, p1j);
   closestPointsBetweenSegments(p0i, p1i, p0j, p1j, c1, c2);
   const Vec3 w = cj.pos - ci.pos;
   int        n = 0;
   seeds[n++]   = c2 - c1;
   seeds[n++]   = w;
   Vec3 perp    = ci.axis() % cj.axis();
   if (perp * w < real_t(0)) perp = -perp;
   seeds[n++] = perp;
   // the two axes themselves (face normals): separating directions when one body lies beyond the other's end face
   seeds[n++] = (ci.axis() * w >= real_t(0)) ? ci.axis() : -ci.axis();
   seeds[n++] = (cj.axis() * w >= real_t(0)) ? cj.axis() : -cj.axis();
   return n;
}
//! Minimum of D on the great circle of directions perpendicular to the axis a (the band plane of a cylinder): D is not
//! differentiable across this circle (the axial term |L.a| has a kink there), and the true minimum of an edge-on-band
//! contact lies exactly on it, where a smooth descent stalls. 1D problem: coarse scan, then golden-section refinement.
inline OverlapMinimum minimiseOverlapOnBandCircle(const CylinderState& ci, const CylinderState& cj, const Vec3& a)
{
   Vec3 u = a % Vec3(real_t(1), real_t(0), real_t(0));
   if (u.length() < real_t(0.5)) u = a % Vec3(real_t(0), real_t(1), real_t(0));
   u /= u.length();
   const Vec3 v  = a % u;
   auto       Lp = [&](real_t phi) { return std::cos(phi) * u + std::sin(phi) * v; };
   auto       Dp = [&](real_t phi) { return cylinderProjectionOverlap(ci, cj, Lp(phi)); };
   constexpr int N = 48;
   const real_t  h = real_t(2) * math::pi / real_t(N);
   int           kBest = 0;
   real_t        dBest = Dp(0);
   for (int k = 1; k < N; ++k)
   {
      const real_t d = Dp(real_t(k) * h);
      if (d < dBest) { dBest = d; kBest = k; }
   }
   real_t lo = (real_t(kBest) - real_t(1)) * h, hi = (real_t(kBest) + real_t(1)) * h;
   const real_t gr = real_t(0.6180339887498949);
   real_t x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo), f1 = Dp(x1), f2 = Dp(x2);
   int evals = N + 2;
   for (int it = 0; it < 60 && hi - lo > real_t(1e-11); ++it, ++evals)
   {
      if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = Dp(x1); }
      else { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = Dp(x2); }
   }
   const real_t phi = (f1 < f2) ? x1 : x2;
   const real_t d   = std::min(f1, f2);
   return (d < dBest) ? OverlapMinimum{ d, Lp(phi), evals } : OverlapMinimum{ dBest, Lp(real_t(kBest) * h), evals };
}

//! The best of the non-smooth candidates of D: the two band circles and the four axis directions (end faces).
inline OverlapMinimum minimiseOverlapOnKinks(const CylinderState& ci, const CylinderState& cj)
{
   OverlapMinimum best = minimiseOverlapOnBandCircle(ci, cj, ci.axis());
   const OverlapMinimum bj = minimiseOverlapOnBandCircle(ci, cj, cj.axis());
   if (bj.value < best.value) best = bj;
   for (const Vec3& L : { ci.axis(), -ci.axis(), cj.axis(), -cj.axis() })
   {
      const real_t d = cylinderProjectionOverlap(ci, cj, L);
      ++best.evaluations;
      if (d < best.value) { best.value = d; best.direction = L; }
   }
   return best;
}

//! Global check for the rare case in which the descent on the sphere ends at a POSITIVE local minimum for a pair without
//! any closed-form scenario. Extended to the unit BALL, D(x) = r_i(x) + r_j(x) - x.(x_j - x_i) is convex (a support
//! function of the Minkowski difference), so a projected subgradient method converges to its global minimum: negative
//! (at |x| = 1) iff the bodies are separated, 0 (at x = 0) iff they overlap. The subgradient is the support-point
//! difference at x/|x|. Returns the smallest D found on the sphere along the iterates (and its direction); stops as
//! soon as D < 0 (separation proven). The result replaces the local minimum when it is smaller.
inline OverlapMinimum cylindersConvexOverlapCheck(const CylinderState& ci, const CylinderState& cj, const Vec3& start,
                                                  real_t startValue, int maxIterations = 400)
{
   OverlapMinimum best{ startValue, start, 0 };
   Vec3           x = start;
   for (int k = 1; k <= maxIterations; ++k)
   {
      const real_t xl = x.length();
      if (xl < real_t(1e-9)) break; // converged to the origin: overlapping
      const Vec3   L = x / xl;
      const real_t D = cylinderProjectionOverlap(ci, cj, L);
      ++best.evaluations;
      if (D < best.value) { best.value = D; best.direction = L; }
      if (D < real_t(0)) break;
      const Vec3   g  = cylinderSupportPoint(ci, L) - cylinderSupportPoint(cj, -L);
      const real_t gl = g.length();
      if (gl < real_t(1e-14)) break;
      x -= (real_t(0.2) / std::sqrt(real_t(k))) * (g / gl);
      const real_t nl = x.length();
      if (nl > real_t(1)) x /= nl;
   }
   return best;
}

inline bool cylindersSeparatedConvexCheck(const CylinderState& ci, const CylinderState& cj, const Vec3& start,
                                          int maxIterations = 400)
{
   return cylindersConvexOverlapCheck(ci, cj, start, cylinderProjectionOverlap(ci, cj, start), maxIterations).value < real_t(0);
}
} // namespace detail

//! Refined cylinder-cylinder contact (see the block comment above). `raw` is the closed-form result for (ci, cj).
//! Exact scenarios pass through unchanged; inexact scenarios are refined from their own normal (plus the standard
//! seeds); "no scenario" pairs are decided by the descent (separation proven if D < 0 is found).
inline ContactResult refineCylinderContact(const CylinderState& ci, const CylinderState& cj, const ContactResult& raw)
{
   const std::string sc(raw.scenario);
   if (raw.contact && (sc == "face-face" || sc == "band-band-parallel" || sc == "band-band-skewed" || sc == "face-band"))
      return raw;

   Vec3 seeds[6];
   int  n = 0;
   if (raw.contact)
   {
      // fast path: if the closed-form normal already satisfies the optimality condition (no tangential component of the
      // gradient S_i(L) - S_j(-L)) and D along it equals the closed-form depth, the closed form found the minimum
      const Vec3   L0    = -raw.normal; // raw.normal points towards ci; L points towards cj
      const real_t scale = ci.radius + cj.radius + real_t(0.5) * (ci.length + cj.length);
      const real_t D0    = detail::cylinderProjectionOverlap(ci, cj, L0);
      if (std::abs(D0 - raw.overlap) < real_t(1e-6) * scale)
      {
         // local optimality test by direct comparison (the minima of the inexact scenarios often lie ON a kink of D,
         // where the gradient does not vanish): no direction tilted by a small angle may give a smaller overlap
         Vec3 t1 = L0 % Vec3(real_t(1), real_t(0), real_t(0));
         if (t1.length() < real_t(0.5)) t1 = L0 % Vec3(real_t(0), real_t(1), real_t(0));
         t1 /= t1.length();
         const Vec3   t2  = L0 % t1;
         const real_t eps = real_t(1e-4);
         bool         optimal = true;
         for (int k = 0; k < 8 && optimal; ++k)
         {
            const real_t phi = real_t(k) * real_t(0.25) * math::pi;
            Vec3         Lk  = L0 + eps * (std::cos(phi) * t1 + std::sin(phi) * t2);
            Lk /= Lk.length();
            optimal = detail::cylinderProjectionOverlap(ci, cj, Lk) >= D0 - real_t(1e-12) * scale;
         }
         if (optimal)
         {
            // a local minimum need not be the global one (a misclassified scenario can sit in another basin of D):
            // accept only if a short descent from the closest-axis-points direction does not find a smaller overlap.
            // Known limitation: for nearly parallel rods D is almost flat, two minima can differ by a few percent, and the
            // accepted normal can then be far off (reference pair 2031: depth +3.8 %, normal 54 deg; GJK/EPA 27 deg there)
            Vec3 std5[5];
            detail::cylinderSeedDirections(ci, cj, std5);
            const detail::OverlapMinimum q = detail::minimiseCylinderProjectionOverlap(ci, cj, std5, 1, false, 20, real_t(1e-7));
            if (q.value >= D0 - real_t(1e-9) * scale) return raw;
         }
      }
      // otherwise descend from the closed-form normal and the closest-axis-points direction only
      Vec3 std5[5];
      detail::cylinderSeedDirections(ci, cj, std5);
      seeds[n++] = L0;
      seeds[n++] = std5[0];
   }
   else
      n += detail::cylinderSeedDirections(ci, cj, seeds + n);

   detail::OverlapMinimum m = detail::minimiseCylinderProjectionOverlap(ci, cj, seeds, n, !raw.contact, 40, real_t(1e-9));
   ContactResult r;
   if (m.value <= real_t(0)) return r; // separated (or exactly touching)
   // the non-smooth candidates (band circles, end-face normals), where the descent stalls: searched when the descent
   // ended with a remaining tangential gradient (it stopped on a kink, not at a smooth minimum) and always for pairs
   // without a scenario
   const real_t lenScale = ci.radius + cj.radius + real_t(0.5) * (ci.length + cj.length);
   const Vec3   gm = detail::cylinderSupportPoint(ci, m.direction) - detail::cylinderSupportPoint(cj, -m.direction);
   const bool   stalled = (gm - (gm * m.direction) * m.direction).length() > real_t(1e-7) * lenScale;
   if (!raw.contact || stalled)
   {
      const detail::OverlapMinimum k = detail::minimiseOverlapOnKinks(ci, cj);
      if (k.value < m.value) m = k;
      if (m.value <= real_t(0)) return r;
   }
   // no scenario and a positive local minimum: run the convex (global) check, which either proves separation or finds a
   // smaller overlap than the local minimum; in the latter case polish it with the descent from the new direction
   if (!raw.contact)
   {
      const detail::OverlapMinimum c = detail::cylindersConvexOverlapCheck(ci, cj, m.direction, m.value);
      if (c.value <= real_t(0)) return r;
      if (c.value < m.value)
      {
         const detail::OverlapMinimum p = detail::minimiseCylinderProjectionOverlap(ci, cj, &c.direction, 1, true);
         m = (p.value < c.value) ? p : c;
         if (m.value <= real_t(0)) return r;
      }
   }
   r.contact  = true;
   r.overlap  = m.value;
   r.normal   = -m.direction;
   // Contact point: the overlap along L spans [S_j(-L).L, S_i(L).L]; its midpoint is S_j(-L) + D/2 L = S_i(L) - D/2 L. The
   // support point of a cylinder is not unique when L is (numerically) its axis (a whole end face) or perpendicular to it
   // (a whole band line); cylinderSupportPoint then returns the face centre / band middle, which can lie far from the
   // actual touching region and would give a wrong torque. So the point is taken from the body whose support point is
   // unique; if neither is (flat on flat), the closed-form point is kept when there is one.
   {
      const Vec3   L  = m.direction;
      auto degenerate = [&](const CylinderState& c, const Vec3& dir) {
         const real_t d = dir * c.axis();
         return std::abs(d) < real_t(1e-6) || (dir - d * c.axis()).length() < real_t(1e-6);
      };
      const bool   di = degenerate(ci, L), dj = degenerate(cj, -L);
      const Vec3   si = detail::cylinderSupportPoint(ci, L), sj = detail::cylinderSupportPoint(cj, -L);
      if (di && !dj) r.point = sj + (real_t(0.5) * m.value) * L;
      else if (dj && !di) r.point = si - (real_t(0.5) * m.value) * L;
      else if (di && dj && raw.contact) r.point = raw.point;
      else r.point = real_t(0.5) * (si + sj);
   }
   r.scenario = raw.contact ? raw.scenario : "refined-no-scenario";
   return r;
}

//! Closed form followed by the refinement -- the complete detection without any GJK/EPA call.
inline ContactResult cylinderCylinderContactRefined(const CylinderState& ci, const CylinderState& cj)
{
   return refineCylinderContact(ci, cj, cylinderCylinderContact(ci, cj));
}

} // namespace paper_analytical

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
