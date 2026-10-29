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
//! \file ClosedFormContactsBoxPoint.h
//! \brief Box x Box closed-form contact: Gottschalk, Lin & Manocha (1996) SAT (boxBoxSAT) plus this
//! codebase's own contact-point construction (box_contact_point, pullIntoBoxIntersection,
//! placeOnNormalLine, boxBoxContactPoint) and vertex-only multi-point contact (BoxVertexContact,
//! boxBoxMultiPointContacts). Split out of PaperAnalyticalContacts.h (namespace paper_analytical)
//! 2026-10-29; see ClosedFormContactsInfra.h for the shared ContactResult type.

#pragma once

#include "ClosedFormContactsInfra.h"

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{

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

// [Claude, 2026-10-01] Contact point of two overlapping boxes for a given contact normal and penetration depth.
// boxBoxSAT used to place the point on the line between the centres ("midpoint between the extreme points along the
// winning axis"), which for elongated, rotated boxes can lie several units from the actual contact and outside both
// bodies; the force then acts with a wrong lever arm, the contact does net work, and elastic collisions gained up to
// 98x their kinetic energy (DEMValidation cuboidPairEnergyFuzz; launched cuboids in the randomised drop tests,
// THESIS_PLAN.md 2026-10-01). The point is now taken from the touching features:
//  1. Contact patch of each box towards the other along the normal: the body-frame bounding box of the part of the box
//     within `depth` of its deepest point (see supportFeature) -- a vertex for a just-touching corner, growing
//     continuously into an edge strip or a whole face.
//  2. The two features are projected onto the plane perpendicular to the normal and intersected (point, segment, or
//     convex polygon); the contact point is the centroid of the intersection, placed at mid-depth between the two
//     features along the normal. Edge-edge -> crossing point; face-face -> centre of the common area; vertex-face ->
//     the vertex. If the projections do not intersect (round-off at a feature boundary), the closest pair of points of
//     the two features is used instead. The patches are bounding boxes and slightly larger than the penetrating parts,
//     so the point is finally pulled into the intersection by alternating projection (pullIntoBoxIntersection).
//  [Claude, 2026-10-02] Step 3 is now closed form (placeOnNormalLine): the line through the in-plane point along the
//     normal is clipped against the slabs of both boxes (clipLineToBox) and the point is the midpoint of the common
//     interval; this moves the point only along the normal, the force's line of action, so the torque is unchanged.
//     Where the line misses the intersection (wedge-shaped overlaps), the point is first moved in-plane to the nearest
//     point of the intersection's projection (boxIntersectionVertices: edges of each box clipped by the other). The
//     alternating projection converges only linearly; in the thin wedges of a settled pile (overlaps ~1e-5) it needed up
//     to its 2000-iteration cap, and in a 16-rank replay of the random_cuboids drop-test pile (4284396, step 70000) it
//     took 70-77 % of the busy ranks' time. It remains only as a last resort against round-off.
namespace box_contact_point
{
struct P2
{
   real_t x, y;
};
inline real_t cross2(const P2& a, const P2& b, const P2& c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

// convex hull (counter-clockwise) of up to 4 points; collinear/duplicate points removed
inline std::vector< P2 > hull(std::vector< P2 > pts, real_t eps)
{
   std::sort(pts.begin(), pts.end(), [](const P2& a, const P2& b) { return a.x < b.x || (!(b.x < a.x) && a.y < b.y); });
   std::vector< P2 > h;
   for (int pass = 0; pass < 2; ++pass)
   {
      const size_t start = h.size();
      for (const auto& p : pts)
      {
         while (h.size() >= start + 2 && cross2(h[h.size() - 2], h[h.size() - 1], p) <= eps) h.pop_back();
         h.push_back(p);
      }
      h.pop_back();
      std::reverse(pts.begin(), pts.end());
   }
   if (h.empty() && !pts.empty()) h.push_back(pts[0]);
   return h;
}

// Sutherland-Hodgman: clip polygon `subj` (any size, may be degenerate) by convex CCW polygon `clip` (>= 3 vertices)
inline std::vector< P2 > clipPolygon(const std::vector< P2 >& subj, const std::vector< P2 >& clip, real_t eps)
{
   std::vector< P2 > out = subj;
   for (size_t e = 0; e < clip.size() && !out.empty(); ++e)
   {
      const P2 a = clip[e], b = clip[(e + 1) % clip.size()];
      std::vector< P2 > in;
      in.swap(out);
      for (size_t i = 0; i < in.size(); ++i)
      {
         const P2     cur = in[i], prv = in[(i + in.size() - 1) % in.size()];
         const real_t dc = cross2(a, b, cur), dp = cross2(a, b, prv);
         const bool   cin = dc >= -eps, pin = dp >= -eps;
         if (cin)
         {
            if (!pin && in.size() > 1) out.push_back(P2{ prv.x + (cur.x - prv.x) * dp / (dp - dc), prv.y + (cur.y - prv.y) * dp / (dp - dc) });
            out.push_back(cur);
         }
         else if (pin && in.size() > 1)
            out.push_back(P2{ prv.x + (cur.x - prv.x) * dp / (dp - dc), prv.y + (cur.y - prv.y) * dp / (dp - dc) });
      }
   }
   return out;
}

// centroid of a convex polygon (area-weighted), of a segment (midpoint), or of a point
inline P2 centroid(const std::vector< P2 >& poly, real_t eps)
{
   real_t a = 0, cx = 0, cy = 0;
   for (size_t i = 0; i < poly.size(); ++i)
   {
      const P2&    p = poly[i];
      const P2&    q = poly[(i + 1) % poly.size()];
      const real_t c = p.x * q.y - q.x * p.y;
      a += c;
      cx += (p.x + q.x) * c;
      cy += (p.y + q.y) * c;
   }
   if (std::abs(a) > eps)
      return P2{ cx / (real_t(3) * a), cy / (real_t(3) * a) };
   // degenerate (segment or point): midpoint of the extreme points along the longest direction
   P2 lo = poly[0], hi = poly[0];
   real_t best = -1;
   for (const auto& p : poly)
      for (const auto& q : poly)
      {
         const real_t d = (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
         if (d > best) { best = d; lo = p; hi = q; }
      }
   return P2{ real_t(0.5) * (lo.x + hi.x), real_t(0.5) * (lo.y + hi.y) };
}

// [Claude, 2026-10-01, revised] contact patch of a box in direction d (unit): the bounding box (in the body frame) of the
// part of the box that lies within `depth` of its deepest point along d -- along axis i it reaches from the deepest end
// over L_i = min(2 h_i, f depth / |d.a_i|), f = kPatchDepthFraction. Grows continuously from a vertex (depth -> 0) via an edge strip to a whole
// face; the earlier yes/no classification (axis free if 2 h_i |d.a_i| <= depth) made the point jump by up to half an
// edge where the classification switched (testBoxBoxContactPointContinuity). Returns the 8 corners of the patch.
// Patch depth fraction: the patch is the part within kPatchDepthFraction * depth of the deepest point. A single-point
// force k * depth conserves energy only if it acts where the depth is measured (the deepest point, as for GJK/EPA); a
// patch over the full depth moved the point by up to depth / (2 |d.a_i|) along nearly parallel features and lost up to
// 3 % in single elastic collisions (cuboidPairEnergyFuzz, 200 collisions). A small fraction keeps the point at the
// deepest feature unless the features are parallel to within ~ this fraction of the depth, and still varies continuously.
constexpr real_t kPatchDepthFraction = real_t(0.1);
inline std::vector< Vec3 > supportFeature(const Vec3& pos, const Vec3 axes[3], const Vec3& half, const Vec3& d, real_t depth)
{
   depth *= kPatchDepthFraction;
   std::vector< Vec3 > pts{ pos };
   for (uint_t i = 0; i < 3; ++i)
   {
      const real_t c    = d * axes[i];
      const real_t sgn  = (c >= real_t(0)) ? real_t(1) : real_t(-1);
      const real_t len  = (std::abs(c) * real_t(2) * half[i] <= depth) ? real_t(2) * half[i] : depth / std::abs(c);
      const real_t hi   = sgn * half[i], lo = sgn * (half[i] - len);
      std::vector< Vec3 > next;
      for (const auto& p : pts)
      {
         next.push_back(p + hi * axes[i]);
         if (len > real_t(0)) next.push_back(p + lo * axes[i]);
      }
      pts.swap(next);
   }
   return pts;
}
} // namespace box_contact_point

//! [Claude, 2026-10-02] iteration statistics of pullIntoBoxIntersection (profiling; per thread-unsafe process counter)
struct BoxPointIterationStats
{
   uint64_t calls = 0, iterations = 0, maxIterations = 0, capped = 0;
   uint64_t lineClipHits = 0;   // [Claude, 2026-10-02] the line through the planned point hit the intersection
   uint64_t polytopeShifts = 0; // missed; moved in-plane onto the projected intersection polytope, then hit
   uint64_t projectionFallbacks = 0; // neither: alternating projection as the last resort
   void record(uint64_t n) { ++calls; iterations += n; maxIterations = std::max(maxIterations, n); if (n >= 2000) ++capped; }
};
inline BoxPointIterationStats& boxPointIterationStats()
{
   static BoxPointIterationStats st;
   return st;
}

//! [Claude, 2026-10-01] moves a point into the intersection of two boxes by alternating projection (clamping in each
//! body frame); continuous in the point and the configuration, and exact once the point is inside both
inline Vec3 pullIntoBoxIntersection(Vec3 x, const Vec3& posA, const Vec3 axesA[3], const Vec3& halfA, const Vec3& posB,
                                    const Vec3 axesB[3], const Vec3& halfB)
{
   auto clampInto = [](const Vec3& p, const Vec3& c, const Vec3 ax[3], const Vec3& h) {
      Vec3 r = c;
      for (uint_t i = 0; i < 3; ++i) r += std::clamp((p - c) * ax[i], -h[i], h[i]) * ax[i];
      return r;
   };
   for (int it = 0; it < 2000; ++it) // converges linearly; slowly only for thin wedge-shaped intersections
   {
      const Vec3 y = clampInto(clampInto(x, posA, axesA, halfA), posB, axesB, halfB);
      if ((y - x).sqrLength() < real_t(1e-28))
      {
         boxPointIterationStats().record(it + 1);
         return y;
      }
      x = y;
   }
   boxPointIterationStats().record(2000);
   return x;
}

//! [Claude, 2026-10-02] clips the line x0 + t*dir against a box (slab test per body axis); narrows [tLo, tHi], which
//! is empty (tLo > tHi) afterwards if the line misses the box
inline void clipLineToBox(const Vec3& x0, const Vec3& dir, const Vec3& pos, const Vec3 axes[3], const Vec3& half,
                          real_t& tLo, real_t& tHi)
{
   const Vec3 rel = x0 - pos;
   for (uint_t i = 0; i < 3; ++i)
   {
      const real_t s = rel * axes[i], d = dir * axes[i];
      if (std::abs(d) < real_t(1e-14))
      {
         if (std::abs(s) > half[i]) { tLo = real_t(1); tHi = real_t(0); return; }
         continue;
      }
      real_t t0 = (-half[i] - s) / d, t1 = (half[i] - s) / d;
      if (t0 > t1) std::swap(t0, t1);
      tLo = std::max(tLo, t0);
      tHi = std::min(tHi, t1);
   }
}

//! [Claude, 2026-10-02] vertices of the intersection polytope of two boxes: the parts of every edge of one box that lie
//! inside the other (clipLineToBox on the segment), for both boxes. These are all vertex-in-box points and all
//! edge-face crossings, i.e. every vertex of the intersection (plus points on its edges, which do not change its hull).
inline void boxIntersectionVertices(const Vec3& posA, const Vec3 axesA[3], const Vec3& halfA, const Vec3& posB,
                                    const Vec3 axesB[3], const Vec3& halfB, std::vector< Vec3 >& out)
{
   auto clipEdges = [&](const Vec3& pos, const Vec3 axes[3], const Vec3& half, const Vec3& oPos, const Vec3 oAxes[3],
                        const Vec3& oHalf) {
      for (uint_t i = 0; i < 3; ++i)
      {
         const uint_t j = (i + 1) % 3, k = (i + 2) % 3;
         const Vec3   d = real_t(2) * half[i] * axes[i];
         for (int sj = -1; sj <= 1; sj += 2)
            for (int sk = -1; sk <= 1; sk += 2)
            {
               const Vec3 p0 = pos - half[i] * axes[i] + real_t(sj) * half[j] * axes[j] + real_t(sk) * half[k] * axes[k];
               real_t     tLo = 0, tHi = 1;
               clipLineToBox(p0, d, oPos, oAxes, oHalf, tLo, tHi);
               if (tLo <= tHi)
               {
                  out.push_back(p0 + tLo * d);
                  out.push_back(p0 + tHi * d);
               }
            }
      }
   };
   clipEdges(posA, axesA, halfA, posB, axesB, halfB);
   clipEdges(posB, axesB, halfB, posA, axesA, halfA);
}

//! [Claude, 2026-10-02] point of the intersection of both boxes on the line through x0 along `normal`: midpoint of the
//! common interval of the two line-box clips. If the line misses the intersection, x0 is first moved perpendicular to
//! the normal to the closest point of the intersection's projection onto the contact plane (the convex hull of the
//! projected intersection vertices); the move is continuous across the hit/miss boundary, and as it is perpendicular to
//! the normal it changes only the lever arm, by the least amount that puts the point on the body overlap. The
//! alternating projection is the last resort if even that line misses (round-off).
inline Vec3 placeOnNormalLine(const Vec3& x0, const Vec3& normal, const Vec3& posA, const Vec3 axesA[3], const Vec3& halfA,
                              const Vec3& posB, const Vec3 axesB[3], const Vec3& halfB, bool legacyProjection = false)
{
   if (legacyProjection) return pullIntoBoxIntersection(x0, posA, axesA, halfA, posB, axesB, halfB);
   auto onLine = [&](const Vec3& x, Vec3& result) {
      real_t tLo = -std::numeric_limits< real_t >::max(), tHi = std::numeric_limits< real_t >::max();
      clipLineToBox(x, normal, posA, axesA, halfA, tLo, tHi);
      if (tLo <= tHi) clipLineToBox(x, normal, posB, axesB, halfB, tLo, tHi);
      if (tLo > tHi) return false;
      result = x + (real_t(0.5) * (tLo + tHi)) * normal;
      return true;
   };
   Vec3 result;
   if (onLine(x0, result))
   {
      ++boxPointIterationStats().lineClipHits;
      return result;
   }
   using namespace box_contact_point;
   std::vector< Vec3 > verts;
   boxIntersectionVertices(posA, axesA, halfA, posB, axesB, halfB, verts);
   if (!verts.empty())
   {
      Vec3 u = normal % Vec3(real_t(1), real_t(0), real_t(0));
      if (u.length() < real_t(0.5)) u = normal % Vec3(real_t(0), real_t(1), real_t(0));
      u /= u.length();
      const Vec3   v     = normal % u;
      const real_t scale = std::max(halfA.max(), halfB.max());
      std::vector< P2 > pts;
      for (const auto& p : verts) pts.push_back(P2{ p * u, p * v });
      const auto poly = hull(pts, real_t(1e-12) * scale * scale);
      const P2   q{ x0 * u, x0 * v };
      // closest point of the (convex) polygon to q; q itself lies outside, otherwise the line would have hit
      P2     best = poly[0];
      real_t bestD = std::numeric_limits< real_t >::max();
      for (size_t i = 0; i < poly.size(); ++i)
      {
         const P2     a = poly[i], b = poly[(i + 1) % poly.size()];
         const real_t ex = b.x - a.x, ey = b.y - a.y, L2 = ex * ex + ey * ey;
         const real_t t  = L2 > real_t(0) ? std::clamp(((q.x - a.x) * ex + (q.y - a.y) * ey) / L2, real_t(0), real_t(1)) : real_t(0);
         const P2     c{ a.x + t * ex, a.y + t * ey };
         const real_t d = (c.x - q.x) * (c.x - q.x) + (c.y - q.y) * (c.y - q.y);
         if (d < bestD) { bestD = d; best = c; }
      }
      // the closest point lies on the boundary of the projection, where the line touches the intersection only at its
      // rim; a step of 1e-9 scale towards the polygon's vertex mean keeps it inside against round-off
      P2 mean{ 0, 0 };
      for (const auto& p : poly) { mean.x += p.x; mean.y += p.y; }
      mean.x /= real_t(poly.size());
      mean.y /= real_t(poly.size());
      const real_t mx = mean.x - best.x, my = mean.y - best.y, ml = std::sqrt(mx * mx + my * my);
      if (ml > real_t(0))
      {
         const real_t step = std::min(ml, real_t(1e-9) * scale);
         best.x += step * mx / ml;
         best.y += step * my / ml;
      }
      if (onLine(x0 + (best.x - q.x) * u + (best.y - q.y) * v, result))
      {
         ++boxPointIterationStats().polytopeShifts;
         return result;
      }
   }
   ++boxPointIterationStats().projectionFallbacks;
   return pullIntoBoxIntersection(x0, posA, axesA, halfA, posB, axesB, halfB);
}

//! normal points from B towards A (boxBoxSAT convention); overlap > 0
//! legacyProjection: the 2026-10-01 point (alternating projection only), kept for the comparison in
//! PaperAnalyticalContactTest.cpp; boxBoxSAT never sets it
inline Vec3 boxBoxContactPoint(const Vec3& posA, const Vec3 axesA[3], const Vec3& halfA, const Vec3& posB,
                               const Vec3 axesB[3], const Vec3& halfB, const Vec3& normal, real_t overlap,
                               bool legacyProjection = false)
{
   using namespace box_contact_point;
   const real_t depth = std::max(overlap, real_t(0));
   // A's deepest feature points towards B (-normal), B's towards A (+normal)
   const auto fA = supportFeature(posA, axesA, halfA, -normal, depth);
   const auto fB = supportFeature(posB, axesB, halfB, normal, depth);
   // in-plane basis
   Vec3 u = normal % Vec3(real_t(1), real_t(0), real_t(0));
   if (u.length() < real_t(0.5)) u = normal % Vec3(real_t(0), real_t(1), real_t(0));
   u /= u.length();
   const Vec3 v = normal % u;
   const real_t scale = std::max(halfA.max(), halfB.max());
   const real_t eps   = real_t(1e-12) * scale * scale;
   auto proj = [&](const std::vector< Vec3 >& f) {
      std::vector< P2 > r;
      for (const auto& p : f) r.push_back(P2{ p * u, p * v });
      return hull(r, eps);
   };
   const auto hA = proj(fA), hB = proj(fB);
   std::vector< P2 > inter;
   if (hA.size() >= 3) inter = clipPolygon(hB, hA, real_t(1e-9) * scale);
   else if (hB.size() >= 3) inter = clipPolygon(hA, hB, real_t(1e-9) * scale);
   else if (hA.size() == 1) inter = hA;
   else if (hB.size() == 1) inter = hB;
   else
   {
      // segment vs segment: crossing point if not parallel, else the overlap of the two collinear intervals
      const P2 a0 = hA[0], a1 = hA[1], b0 = hB[0], b1 = hB[1];
      const real_t rx = a1.x - a0.x, ry = a1.y - a0.y, sx = b1.x - b0.x, sy = b1.y - b0.y;
      const real_t den = rx * sy - ry * sx;
      if (std::abs(den) > real_t(1e-9) * (rx * rx + ry * ry + sx * sx + sy * sy))
      {
         const real_t t = ((b0.x - a0.x) * sy - (b0.y - a0.y) * sx) / den;
         const real_t w = ((b0.x - a0.x) * ry - (b0.y - a0.y) * rx) / den;
         if (t >= real_t(-1e-9) && t <= real_t(1) + real_t(1e-9) && w >= real_t(-1e-9) && w <= real_t(1) + real_t(1e-9))
            inter.push_back(P2{ a0.x + t * rx, a0.y + t * ry });
      }
      else
      {
         const real_t L2 = rx * rx + ry * ry;
         if (L2 > eps)
         {
            auto par = [&](const P2& p) { return ((p.x - a0.x) * rx + (p.y - a0.y) * ry) / L2; };
            const real_t lo = std::max(real_t(0), std::min(par(b0), par(b1)));
            const real_t hi = std::min(real_t(1), std::max(par(b0), par(b1)));
            if (lo <= hi + real_t(1e-9))
            {
               const real_t m = real_t(0.5) * (lo + hi);
               inter.push_back(P2{ a0.x + m * rx, a0.y + m * ry });
            }
         }
      }
   }
   // mid-depth level along the normal between the two features
   // mid-depth between the deepest points of the two boxes along the normal (continuous in the configuration)
   real_t levelA = std::numeric_limits< real_t >::max(), levelB = -std::numeric_limits< real_t >::max();
   for (const auto& p : fA) levelA = std::min(levelA, p * normal);
   for (const auto& p : fB) levelB = std::max(levelB, p * normal);
   const real_t level  = real_t(0.5) * (levelA + levelB);
   if (!inter.empty())
   {
      const P2 c = centroid(inter, eps);
      return placeOnNormalLine(c.x * u + c.y * v + level * normal, normal, posA, axesA, halfA, posB, axesB, halfB,
                               legacyProjection);
   }
   // fallback: closest points between the two features (as segments / points)
   Vec3 best = real_t(0.5) * (fA[0] + fB[0]);
   real_t bestD = std::numeric_limits< real_t >::max();
   for (size_t i = 0; i < fA.size(); ++i)
      for (size_t k = 0; k < fB.size(); ++k)
      {
         Vec3 cA, cB;
         closestPointsBetweenSegments(fA[i], fA[(i + 1) % fA.size()], fB[k], fB[(k + 1) % fB.size()], cA, cB);
         const Vec3 dAB = cA - cB;
         const real_t dperp = (dAB - (dAB * normal) * normal).sqrLength();
         if (dperp < bestD) { bestD = dperp; best = real_t(0.5) * (cA + cB); }
      }
   return placeOnNormalLine(best - (best * normal) * normal + level * normal, normal, posA, axesA, halfA, posB, axesB,
                            halfB, legacyProjection);
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
   // [Claude, 2026-10-01] point from the touching features (detail::boxBoxContactPoint), no longer on the line between
   // the centres (former: 0.5 (posA + posB) + 0.5 sign (ra - rb) axis); see that function's comment.
   // ANTIDUNES_SAT_FORMER_POINT (env, any value): the former point, for before/after comparisons in the validation runs only
   static const bool formerPoint = std::getenv("ANTIDUNES_SAT_FORMER_POINT") != nullptr;
   result.point = formerPoint ? Vec3(real_t(0.5) * (posA + posB) + real_t(0.5) * sign * (bestRa - bestRb) * bestAxis)
                              : detail::boxBoxContactPoint(posA, axesA, halfA, posB, axesB, halfB, normal, best.gap);

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

// [Claude, 2026-10-29] eberly_box_distance (Eberly's closed-form box-box distance, ported for
// useSpheroBoxAnalytical) and boxBoxSpheroAnalytical removed entirely -- unused by any thesis-cited
// result. Archived at archive/unused_contact_methods/spherosimplex_box_contact.h (see ARCHIVE_README.md).

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

} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla
