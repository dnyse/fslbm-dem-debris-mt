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
// ARCHIVED fraige_multipoint_box_contact.h -- unused Fraige et al. (2008) multi-point box-box
// contact variant (edge-edge/edge-face modes, superset of the live vertex-only boxBoxMultiPointContacts).
// Unused by any thesis-cited result; removed from the live build 2026-10-29 (ParticleKernels.h's
// useFraigeBoxContactModel dispatch and SimulationParams.h's flag removed too). See archive/ARCHIVE_README.md.
#pragma once

#include "../../ClosedFormContacts.h"

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{
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
} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla
