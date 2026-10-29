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
//! \file spherosimplex_box_contact.h
//! \brief ARCHIVED, unused spherosimplex box-box contact method (GJK/EPA-with-margin + Eberly
//! closed-form distance variants). Unused by any thesis-cited result; removed from the live build
//! 2026-10-29 (ParticleKernels.h's useSpheroBoxContact/useSpheroBoxAnalytical dispatch and
//! SimulationParams.h's flags removed too). See archive/ARCHIVE_README.md.
#pragma once

#include "../../PaperAnalyticalContacts.h"

namespace walberla
{
namespace antidunes
{
namespace paper_analytical
{
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
} // namespace paper_analytical
} // namespace antidunes
} // namespace walberla
