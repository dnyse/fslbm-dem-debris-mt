//======================================================================================================================
//! \file multipoint_box_contact.h
//! ARCHIVED 2026-11-20 from ClosedFormContactsBoxPoint.h (useMultiPointBoxContact): one contact point per vertex of
//! one cuboid embedded in the other, as an alternative to the single SAT contact point. Not used by any thesis result;
//! kept for reference only and not compiled. See ../ARCHIVE_README.md.
//======================================================================================================================
#pragma once

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
