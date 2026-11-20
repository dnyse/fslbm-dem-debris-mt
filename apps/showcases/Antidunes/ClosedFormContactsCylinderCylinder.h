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
//! \file ClosedFormContactsCylinderCylinder.h
//! \brief Cylinder x Cylinder closed-form contact: Kodam et al. (2010) Part I, corrected per Guo,
//! Wassgren, Ketterhagen, Hancock & Curtis (2012), Powder Technology 228, 193-198.
//! rimVsFaceOrBand, cylinderCylinderContact, cylinderCylinderResultIsUncertain,
//! refineCylinderContact, cylinderCylinderContactRefined. Depends on ClosedFormContactsCylinderPlaneGeometry.h
//! (CylinderState, faceFaceStyleContactPoint) and ClosedFormContactsInfra.h.

#pragma once

#include "ClosedFormContactsCylinderPlaneGeometry.h"

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{

//! Face-edge / special-contact / band-edge, all three tested against `host`'s face/band using a
//! candidate point on `probe`'s rim -- per Kodam et al. Section 3.2.3/3.2.5 AS CORRECTED by Guo,
//! Wassgren, Ketterhagen, Hancock & Curtis (2012), Powder Technology 228, 193-198, Table 1.
//!
//! Returned normal/point use the "natural" host->probe direction (like eq. 29/39/50/54/80 elsewhere
//! in this file); the caller negates when host corresponds to this module's 1st cylinderCylinderContact
//! argument, to match this file's toward-1st-argument convention (see the face-face branch's comment).
//!
//! The face-edge candidate uses Kodam's exact closed form (eq. 45-47, same structure as the cylinder-plane
//! edge-plane formula, detail::rimPointTowardPlane); the band-edge/special-contact candidate uses the rim-circle
//! quartic (Kodam eq. 67-77); the dispatch follows Guo et al.'s corrected criteria throughout (Kodam's original
//! face-edge/band-edge criteria give false negatives for oblique tilts, PaperAnalyticalContactTest.cpp).
//!
//! kodamCapSelectionByAxisDistance: which end-cap circle of the probe the band-edge branch examines. true (default)
//! = Kodam et al. (2010) eq. 68-70: S = sign(d_bottom - d_top), d = the cap centre's distance from the HOST'S AXIS
//! measured in the host's xy-plane, i.e. the cap circle closer to the host's axis. false = the cap whose centre is
//! closer to the host's CENTRE (3D distance), a deviation from the paper: for a nested rod pair it picks the circle
//! that merely grazes the band (overlap ~0.003 LU, wrong normal) instead of the one that dips into it. Kept so tests
//! can reproduce it.
inline bool kodamCapSelectionByAxisDistance = true;

//! ModelParameters.kodamBandEdgeAxialCheck (default false). A band-edge result of cylinderCylinderContact is checked
//! against the projection overlap D(n) along its own normal (detail::cylinderProjectionOverlap: the distance the
//! bodies must move along n until their projections onto n separate, an upper bound of the depth that equals it at
//! the true contact normal). If D(n) > kBandEdgeConsistencyFactor * overlap, the closed form is demonstrably wrong for
//! this pair, and the result is replaced by the support-function refinement (refineCylinderContact, thesis appendix
//! app:cc-refinement), started from it. Failure mode: two equal-radius cylinders stacked nearly coaxially overlap
//! face-to-face (0.46 in the round-container drop test); face-edge is rejected because the probe rim's lowest point
//! lies just outside the host's face disc, special contact because delta_a(P) > delta_r(P), so band-edge reports the
//! radial depth at the single point P (0.005, horizontal normal) although D along that normal is 7.9; when the pose
//! leaves that window the depth jumps to ~0.5 and injects an energy spike. Returning face-edge with D along the host
//! face normal instead over-estimates oblique contacts (the face normal is not their contact normal).
//! Threshold: in the 1019-contact reference sample, D(n)/overlap is at most 2.06 for the correct band-edge depths and
//! 6.8-149 for four of the five wrong ones (the fifth has D(n) = overlap and is not detectable this way).
inline bool   kodamBandEdgeAxialCheck    = false;
inline real_t kBandEdgeConsistencyFactor = real_t(3);
namespace detail
{
inline real_t cylinderProjectionOverlap(const CylinderState& ci, const CylinderState& cj, const Vec3& L);
} // namespace detail
struct RefineTrace;
inline ContactResult refineCylinderContact(const CylinderState& ci, const CylinderState& cj, const ContactResult& raw,
                                           RefineTrace* trace);

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
         // contact point: E is the probe's rim point that has penetrated deepest into the host
         // (|relE_z| < halfLHost), deltaA its distance below the host face at z = sign(relE_z)*halfLHost, and
         // result.normal (host->probe here) points from E towards that face. The overlap along the normal spans
         // [E, E + deltaA*normal], whose midpoint is E + 0.5*deltaA*normal (as band-edge's P + 0.5*deltaRp*nHat
         // below). Checked by testCylinderCylinderFaceEdgeClean.
         result.point    = E + result.normal * real_t(0.5) * deltaA;
         result.scenario = "face-edge";
         return result;
      }
   }

   // --- Band-edge / special-contact candidate: P, the point on probe's rim closest to host's BAND
   // (Kodam eq. 67-77, a quartic in sin(theta); solved here by closestThetaOnCircleToLine as a closed-form quartic in
   // tan(theta/2) -- see that function's own comment for the derivation). ---------------------------
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
   // Guo's criterion compared delta_a AT P with delta_r. P (closest to the host axis) is the
   // probe's deepest point RADIALLY, but not axially, so delta_a(P) underestimates the axial penetration and the
   // comparison picks the axial ("special") contact too often: in the fuzz suite, most special-contact
   // disagreements with GJK/EPA had GJK normals ~90 degrees from the host axis, i.e. the true escape direction
   // was radial. The comparison therefore uses the corrected axial depth computed below (deepest rim point inside
   // the host); if it is not smaller than delta_r, the contact falls through to band-edge (radial depth at P).
   if (sameSide && probeCentreBeyondFace && deltaAp < deltaRp)
   {
      // Corrected special-contact depth. Guo et al. (2012, Table 1 / Fig. 8) place this
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

   // Band-edge (eq. 79-81), CORRECTED per Guo et al.: Kodam's original delta_r < delta_a requirement
   // is dropped entirely (Guo et al.'s Fig. 7 shows a genuine band-edge contact where it doesn't hold); the
   // "NOT face-edge" and "NOT special-contact" exclusions above (via the early returns) take its place.
   result.contact  = true;
   // overlap = delta_r(P), radial depth at the rim point closest to the host axis (Kodam et al. 2010
   // Part I, eq. 79-81), reached via Guo et al. (2012) Powder Technology 228:193-198's corrected dispatch
   result.overlap  = deltaRp;
   const Vec3 nHat = host.R * ((rhoP > real_t(1e-10)) ? Vec3(relP[0], relP[1], real_t(0)) / rhoP
                                                       : Vec3(real_t(1), real_t(0), real_t(0)));
   result.normal   = nHat;
   result.point    = P + nHat * real_t(0.5) * deltaRp;
   result.scenario = "band-edge";
   return result;
}

//! Cylinder-cylinder contact: dispatches across all 7 scenarios in Fig. 9 (face-face, face-band,
//! face-edge, band-band parallel, band-band skewed, band-edge, edge-edge) plus Guo et al.'s special contact,
//! tested in the order face-face / band-band parallel, face-band, band-band skewed, face-edge / special contact /
//! band-edge (both directions), edge-edge. `i`/`j` naming mirrors the paper; `ci` is the 1st argument, `cj` the 2nd.
// useFastEdgeEdgeReference: forwarded to findRimCrossing's useFastReference for the edge-edge branch's two
// rim-crossing calls (default false); affects only edge-edge.
// maxIterations, quarticSolver, lazyThetaRef: forwarded ONLY to the edge-edge branch's two findRimCrossing calls, not
// to rimVsFaceOrBand's quartic solves (band-edge/special-contact/face-edge), which use closestThetaOnCircleToLine's own
// defaults (30 iterations, Ferrari). Every quartic solve of the dispatch therefore defaults to Ferrari, validated
// across all scenarios by testCylinderCylinderFullComparisonAndHybridFeasibility (5000 trials); see also
// QuarticSolverComparisonTest.cpp (synthetic quartics with known roots) and PaperAnalyticalContactTest.cpp's edge-edge
// comparison (real geometry vs GJK/EPA). They are exposed so PaperAnalyticalContactTest.cpp can compare iteration
// counts, solvers and reference choices on identical trials in one process; lazyThetaRef: see findRimCrossing.
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
   // deltaAxial = (Li/2+Lj/2) - |zIJ|, deltaRadial = (Ri+Rj) - rhoIJ  (Kodam et al. 2010 Part I, eq. 26-27/53)
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
         // testSphereSphereContact et al. -- see the band-band-skewed branch below, which negates for
         // the same reason).
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
      // eq. 35: S = |x_cyl_i^j . e_cyl_j,z|, the projection of (ci.pos - cj.pos) onto cylinder j's axis, computed
      // as `Sj` from ci's position in cylinder j's body frame (equivalent to the paper's dot-product form).
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
            // eq. 39, like eq. 29/54, is worded "points ... toward the cylinder's centre", i.e. toward j (cj
            // approaches ci's face here), so negated as for face-face/band-band-parallel above.
            result.normal   = ci.R * Vec3(real_t(0), real_t(0), -math::sign(zIJ));
            // Contact point: eq. 40-44, via detail::faceBandContactPoint (derivation there, including the
            // line-circle intersection eq. 44 does not spell out).
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
         // eq. 56. Must use the WORLD-frame relative position, not the body-frame `relI` of the parallel-axis
         // branches: rHatN is built from world-frame axes, and dotting it with a body-frame vector gives a wrong
         // (generally too small) `d`, i.e. false-positive band-band-skewed contacts, whenever ci is rotated
         // (covered by PaperAnalyticalContactTest.cpp's randomised fuzz test, which rotates both bodies).
         const real_t d = std::abs((cj.pos - ci.pos) * rHatN);

         // Closest points on the two (infinite) centrelines. The paper derives this via Vega & Lago (1994)'s
         // eq. 59-60 in cylinder i's body frame; implemented here as the standard closest-point-between-two-lines
         // formula in WORLD vectors (the same quantities as Vega & Lago's l_i*/l_j*, without an extra frame
         // transform).
         //
         // NEAR-PARALLEL axes (1 - (e_i.e_j)^2 = sin^2(theta) < 1e-10, i.e. theta < 1e-5 rad -- ordinary
         // floating-point rotation noise of a body placed axis-aligned): the two infinite centrelines approach
         // each other only about (separation / theta) ~ 1e5-1e6 LU away along the axes, so eq. 63's |l*| < L/2
         // conditions REJECT the band-band contact. That is the right answer: any near-parallel pair that really
         // touches band-to-band or face-to-face is caught by the parallel-axes branch at the top of this function
         // (|zIJ| < (L_i+L_j)/2 AND rhoIJ < R_i+R_j). Setting l* = 0 in this regime ("guard against division by ~0")
         // would make both l* conditions vacuously true while r_hat = (e_i x e_j)/|e_i x e_j| is a pure noise
         // direction (e.g. perpendicular to the actual separation => d = |x.r_hat| = 0), reporting two rods with a
         // clear gap as overlapping by R_i+R_j. As in Kodam's equations, the branch is simply not applicable here
         // (CylinderNearParallelSkewedTest.cpp checks against a literal implementation of eq. 55-66 with the exact
         // l*).
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
               const Vec3 diff = Pj - Pi;
               // contact point: Kodam et al.'s eq. 65 takes the midpoint of the two closest AXIS points P*, Q*,
               // which lies in the overlap only for equal radii; otherwise it is off by |R_i - R_j|/2 along the
               // normal (171 of 183 unequal-radius band-band-skewed points outside the intersection in
               // PaperAnalyticalContactTest.cpp's contact-point check). The midpoint of the two SURFACE points,
               // P* + R_i u and Q* - R_j u with u = (Q* - P*)/|Q* - P*|, is used instead; identical to eq. 65 for
               // R_i = R_j.
               {
                  const real_t dl = diff.length();
                  const Vec3   u  = (dl > real_t(1e-10)) ? diff / dl : -rHatN;
                  result.point    = real_t(0.5) * (Pi + Pj) + (real_t(0.5) * (ci.radius - cj.radius)) * u;
               }
               // eq. 66, i<-j flipped to match this codebase's idx2->idx1 convention
               result.normal   = (diff.length() > real_t(1e-10)) ? -diff.getNormalized() : rHatN;
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
   // above matched.
   {
      // j's rim against i's face/band: rimVsFaceOrBand returns i->j (host->probe); negate for this
      // file's toward-1st-argument convention (see the face-face branch above).
      // kodamBandEdgeAxialCheck: an inconsistent band-edge result is refined (see the flag's comment)
      const auto checked = [&](const ContactResult& res) {
         if (!kodamBandEdgeAxialCheck || std::string(res.scenario) != "band-edge") return res;
         if (detail::cylinderProjectionOverlap(ci, cj, -res.normal) <= kBandEdgeConsistencyFactor * res.overlap) return res;
         return refineCylinderContact(ci, cj, res, nullptr);
      };
      ContactResult jRimVsI = rimVsFaceOrBand(ci, cj);
      if (jRimVsI.contact)
      {
         jRimVsI.normal = -jRimVsI.normal;
         return checked(jRimVsI);
      }

      // i's rim against j's face/band: host=cj, probe=ci, so rimVsFaceOrBand's natural host->probe
      // direction (j->i) IS ALREADY this file's toward-1st-argument (i) convention -- no negation.
      const ContactResult iRimVsJ = rimVsFaceOrBand(cj, ci);
      if (iRimVsJ.contact) { return checked(iRimVsJ); }
   }

   // --- Edge-edge (Section 3.2.6, eq. 82-101): both rims (circular edges) facing each other. ------
   //
   // Implements Kodam's eq. 82-101: findRimCrossing finds points A, B on rim j (where it crosses cylinder i's
   // lateral surface and near face respectively, eq. 83-86) -- their combined success IS Kodam's existence test,
   // not a proximity threshold -- and symmetrically C, D on rim i vs cylinder j. Given both chords, eq. 87-98's
   // exact skew-line closest distance between them is used directly as a small positive overlap; since contact is
   // established by geometric containment, existence and sign are consistent. Pairs for which no chord pair
   // exists end as "no scenario" and fall back via Tier 1 (cylinderCylinderResultIsUncertain); a found edge-edge
   // contact is trusted directly.
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
                  result.overlap  = dist; // contact already established geometrically above, so d_perp is a
                                           // (small, near-first-contact) penetration estimate
                  // Pi is on rim i where it crosses INTO cylinder j (near cj); Pj is on rim j where it crosses
                  // INTO cylinder i (near ci). The normal points toward the 1st argument (ci), i.e. from the
                  // cj-side point to the ci-side point: Pj-Pi. Against GJK/EPA this direction agrees to a few
                  // degrees wherever the magnitude is also correct (PaperAnalyticalContactTest.cpp's edge-edge
                  // normal-angle breakdown); the residual edge-edge disagreement is a wrong-MAGNITUDE issue (the
                  // curvature-coupling feature gap), not a direction error.
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
//! wrong answer. Calibrated and tested in PaperAnalyticalContactTest.cpp (cylinderCylinderContactHybrid and its
//! accompanying tests).
//!
//! TIER 1 (dispatch found nothing, but a real contact is plausible): a cheap bounding-sum gate -- centre distance vs.
//! the sum of both cylinders' "reach" (radius + half-length each, the distance at which a fully-extended,
//! perfectly-aimed pair could just touch tip-to-tip). The same style of coarse bounding check DEM codes use before
//! fine contact detection (cf. Guo et al. 2012 Section 2.1's circumscribing-sphere/cube coarse detection);
//! deliberately loose, since triggering it only means "double-check with GJK/EPA". A 10000-trial fuzz run found it
//! catches 99.2% of the "no scenario matched" gap (the residual 0.8% are real contacts this bound misses, see
//! PaperAnalyticalContactTest.cpp's dispatch-gap characterisation, resolved downstream by GJK/EPA).
//!
//! Kodam et al. 2010's ten contact scenarios (Part I) only test for POSITIVE contact; none is a "definitely
//! separated" fast path. Their Section 4 conclusion recommends adding "an intermediate detection scheme such as a
//! bounding sphere or rectangular prism" for efficiency, without specifying one. This Tier 1 gate and
//! ParticleKernels.h's TIGHTER per-pair reject run before cylinderCylinderContact (circumscribing-sphere sum,
//! sqrt(radius^2+halfLength^2) each, or the axis-segment capsule bound with useTightLayer1) are this codebase's answer
//! to that recommendation, not derived from the paper. Since that outer reject is strictly tighter than Tier 1 (r+L/2
//! per side >= sqrt(r^2+(L/2)^2) per side), any pair reaching Tier 1 has already passed it, so Tier 1 never rejects
//! in production: every "no scenario matched" case falls back to GJK/EPA. Most such calls are expected to be genuine
//! non-contacts (normal broad-phase population); ParticleKernels.h's ContactDetectionStats::
//! cylFallbackTier1GjkFoundContact/GjkFoundNoContact counters measure that split.
//!
//! special-contact: unconditional distrust. Guo et al.'s overlap formula (2012, Table 1: axial distance from a
//! single rim point to the host's face, delta_a) is their fix for this configuration's DETECTION (not silently
//! dropped or misclassified as band-edge), not a validated magnitude. For testCylinderCylinderSpecialContact's
//! geometry it gives 0.181, GJK/EPA ~0.610, and a fully independent brute-force minimum-translation-distance
//! search over each cylinder's support function (testCylinderCylinderSpecialContactBruteForceMTD) 0.6095, i.e. Guo's
//! formula underestimates by >3x. rimVsFaceOrBand replaces it by the deepest rim point inside the host (0.614 on that
//! geometry), but the scenario is still routed to GJK/EPA. With Guo's formula a 10000-trial fuzz run found this
//! scenario wrong 93% of the time it fires; it is rare (under 2% of real contacts), so routing it to GJK/EPA costs
//! almost nothing.
//!
//! EVERY OTHER SCENARIO IS TRUSTED DIRECTLY. A magnitude cross-check against a looser circumscribing-sphere upper
//! bound ("is MORE overlap possible than reported?") is not used: it pushes the GJK/EPA fallback rate to ~90%, and it
//! is structurally ONE-SIDED -- it catches UNDERESTIMATES (the typical failure of band-edge/edge-edge/face-edge/
//! special-contact) but not OVERESTIMATES; the deep/degenerate-overlap stress test
//! (testCylinderCylinderLargeOverlapDegenerateStress, 30%-150% of the smaller radius, beyond the 15% cap of the other
//! fuzz tests) has a case ~4.7x too large (372.6% relative error) that already exceeds the bound itself.
//! With Tier 1 plus special-contact's routing the GJK/EPA fallback rate is ~3.3%, with a residual wrong-magnitude rate
//! of ~2.3-2.5% of genuine contacts, concentrated in band-edge/edge-edge/face-edge with modest errors (median ~2.7%
//! relative, max ~13.5% absolute in a 10000-trial run). face-face, band-band-parallel and band-band-skewed are exact
//! under their own preconditions (see their branch comments).
inline bool cylinderCylinderResultIsUncertain(const CylinderState& ci, const CylinderState& cj,
                                               const ContactResult& analytical)
{
   if (!analytical.contact)
   {
      const real_t centreDist = (cj.pos - ci.pos).length();
      const real_t reach = ci.radius + cj.radius + real_t(0.5) * ci.length + real_t(0.5) * cj.length;
      return centreDist < reach; // Tier 1 (loose centre-distance bound)
   }

   return std::string(analytical.scenario) == "special-contact";
}

//======================================================================================================================
// Support-function refinement of cylinder-cylinder contacts (thesis appendix, app:cc-refinement).
//
// For two convex, centrally symmetric bodies, the overlap of their projections onto a unit direction L (pointing from
// ci towards cj) is
//    D(L) = r_i(L) + r_j(L) - L.(x_j - x_i),       r_k(L) = (length_k/2) |L.a_k| + R_k sqrt(1 - (L.a_k)^2),
// the bodies are separated iff D(L) < 0 for SOME L (separating-axis theorem), and the penetration depth is min_L D(L),
// attained along the contact normal. Seen as a positively homogeneous function of L in R^3, D has the gradient
//    grad D(L) = S_i(L) - S_j(-L),
// the difference of the two support points (the point of ci farthest along L, the point of cj farthest along -L), and
// grad D(L).L = D(L). The minimisation over the unit sphere is therefore a projected-gradient descent whose gradient
// costs two support-point evaluations -- a few dozen flops, no simplex, no polytope.
//
// Used for two things:
//  (a) refineCylinderContact: the closed-form scenarios that are not exact (edge-edge, band-edge, face-edge, special
//      contact) get their depth and normal from a descent started at their own normal, which is typically close to the
//      optimum. The exact scenarios (face-face, band-band parallel/skewed, face-band) are returned unchanged.
//  (b) the "no scenario matched" case: a descent from the five seed directions (cylinderSeedDirections), then the
//      non-smooth candidates (minimiseOverlapOnKinks) and the convex check (cylindersConvexOverlapCheck). Any
//      direction with D < 0 PROVES separation (early exit). If none is found, the pair is reported as a contact with
//      the refined depth.
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

//! Point of cylinder c farthest along the unit direction L (where the support set is an end face or a band segment
//! rather than a single point, its centre).
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

//! Global check for the rare case in which the descent on the sphere ends at a POSITIVE local minimum for a pair
//! without any closed-form scenario. Extended to the unit BALL, D(x) = r_i(x) + r_j(x) - x.(x_j - x_i) is convex
//! (a support function of the Minkowski difference), so a projected subgradient method converges to its global
//! minimum: negative (at |x| = 1) iff the bodies are separated, 0 (at x = 0) iff they overlap. The subgradient is the
//! support-point difference at x/|x|. Returns the smallest D found on the sphere along the iterates (and its
//! direction); stops as soon as D < 0 (separation proven). The result replaces the local minimum when it is smaller.
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

namespace detail
{
//! Closest point of the solid cylinder c to p (identity if p is inside): clamp the axial and radial coordinates.
inline Vec3 projectIntoCylinder(const CylinderState& c, const Vec3& p)
{
   const Vec3   a  = c.axis();
   const Vec3   d  = p - c.pos;
   real_t       z  = d * a;
   Vec3         rv = d - z * a;
   const real_t r  = rv.length();
   z               = std::clamp(z, -real_t(0.5) * c.length, real_t(0.5) * c.length);
   if (r > c.radius) rv *= c.radius / r;
   return c.pos + z * a + rv;
}

//! Moves p into the intersection of the two cylinders by alternating projection (converges for two convex bodies that
//! overlap); used because at minima on a kink of D the two support points are not exactly aligned along L and the
//! midpoint of the overlap can lie slightly outside one body.
inline Vec3 pullIntoIntersection(const CylinderState& ci, const CylinderState& cj, Vec3 p, int iterations = 8)
{
   for (int k = 0; k < iterations; ++k)
   {
      const Vec3 q = projectIntoCylinder(cj, projectIntoCylinder(ci, p));
      if ((q - p).sqrLength() < real_t(1e-24)) break;
      p = q;
   }
   return p;
}

//! Contact point of a refined contact with direction L (from ci towards cj) and overlap D: the midpoint of the overlap
//! along L, S_j(-L) + D/2 L = S_i(L) - D/2 L, taken from the body whose support point is unique; `fallback` if neither
//! is unique (flat on flat).
inline Vec3 refinedContactPointRaw(const CylinderState& ci, const CylinderState& cj, const Vec3& L, real_t D, const Vec3& fallback)
{
   auto degenerate = [&](const CylinderState& c, const Vec3& dir) {
      const real_t d = dir * c.axis();
      return std::abs(d) < real_t(1e-6) || (dir - d * c.axis()).length() < real_t(1e-6);
   };
   const bool di = degenerate(ci, L), dj = degenerate(cj, -L);
   const Vec3 si = cylinderSupportPoint(ci, L), sj = cylinderSupportPoint(cj, -L);
   if (di && !dj) return sj + (real_t(0.5) * D) * L;
   if (dj && !di) return si - (real_t(0.5) * D) * L;
   if (di && dj) return fallback;
   return real_t(0.5) * (si + sj);
}
inline Vec3 refinedContactPoint(const CylinderState& ci, const CylinderState& cj, const Vec3& L, real_t D, const Vec3& fallback)
{
   return pullIntoIntersection(ci, cj, refinedContactPointRaw(ci, cj, L, D, fallback));
}

} // namespace detail

//! Which exit refineCylinderContact took, for profiling (ParticleKernels.h ContactDetectionStats::cylRefine*):
//! exact scenario passed through, closed form confirmed by the fast path, closed-form contact re-descended (contact or
//! separated), pair without a scenario decided separated or in contact. `fastReject` says why the fast path failed;
//! `kinks` / `convex` whether the kink search / the convex check ran.
enum class RefinePath : int { Exact = 0, FastAccepted, ContactDescended, ContactDescendedSeparated, NoScenarioSeparated,
                              NoScenarioContact, Count };
enum class FastReject : int { None = 0, DepthMismatch, NotLocallyOptimal, SmallerElsewhere, Count };
struct RefineTrace
{
   RefinePath path       = RefinePath::Exact;
   FastReject fastReject = FastReject::None;
   bool       kinks      = false;
   bool       convex     = false;
   real_t     fastD0     = real_t(-1); // overlap along the closed-form normal (fast path only)
   real_t     fastQ      = real_t(-1); // the fast path's short descent minimum (when it ran)
};

//! Refined cylinder-cylinder contact (see the block comment above). `raw` is the closed-form result for (ci, cj).
//! Exact scenarios pass through unchanged; inexact scenarios are accepted if a fast-path check finds their normal
//! optimal, else refined by a descent from their own normal and the closest-axis-points direction; "no scenario"
//! pairs are decided by the descent (separation proven if D < 0 is found).
inline ContactResult refineCylinderContact(const CylinderState& ci, const CylinderState& cj, const ContactResult& raw,
                                           RefineTrace* trace = nullptr)
{
   RefineTrace tr;
   const auto  done = [&](const ContactResult& res, RefinePath path) {
      tr.path = path;
      if (trace) *trace = tr;
      return res;
   };
   const std::string sc(raw.scenario);
   if (raw.contact && (sc == "face-face" || sc == "band-band-parallel" || sc == "band-band-skewed" || sc == "face-band"))
      return done(raw, RefinePath::Exact);

   Vec3 seeds[6];
   int  n = 0;
   // relative tolerance of the fast path (see below), and the short descent's minimum, which the full
   // refinement keeps as a candidate when the fast path rejects
   const real_t           kFastRelTol = real_t(1e-3);
   detail::OverlapMinimum fastQ{ std::numeric_limits< real_t >::max(), Vec3(real_t(0)), 0 };
   if (raw.contact)
   {
      // fast path: if the closed-form normal already satisfies the optimality condition (no tangential component of the
      // gradient S_i(L) - S_j(-L)) and D along it equals the closed-form depth, the closed form found the minimum
      const Vec3   L0    = -raw.normal; // raw.normal points towards ci; L points towards cj
      const real_t scale = ci.radius + cj.radius + real_t(0.5) * (ci.length + cj.length);
      const real_t D0    = detail::cylinderProjectionOverlap(ci, cj, L0);
      tr.fastD0          = D0;
      // relative tolerance: the face-edge closed form differs from D along its own normal by a median of 0.13 % in a
      // dense packing; an absolute 1e-6 * scale alone would send these contacts to the full refinement for nothing
      if (std::abs(D0 - raw.overlap) < std::max(real_t(1e-6) * scale, kFastRelTol * D0))
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
            // Known limitation: for nearly parallel rods D is almost flat, two minima can differ by a few percent, and
            // the accepted normal can be far off (reference pair 2031: depth +3.8 %, normal 54 deg; GJK/EPA 27 deg)
            Vec3 std5[5];
            detail::cylinderSeedDirections(ci, cj, std5);
            const detail::OverlapMinimum q = detail::minimiseCylinderProjectionOverlap(ci, cj, std5, 1, false, 20, real_t(1e-7));
            tr.fastQ                       = q.value;
            fastQ                          = q;
            // tolerance relative to the depth: in a dense packing (DropTest Kodam packing, ANTIDUNES_REFINE_DUMP), the
            // short descent improves on 60 % of the face-edge results, by a median of 3e-4 D0 and by less than 1e-2 D0
            // in 99 % of them; each rejection costs the full refinement, ~5x the fast path. An improvement below
            // kFastRelTol of the depth is accepted, with the smaller of the two overlaps as the depth.
            if (q.value >= D0 - std::max(real_t(1e-9) * scale, kFastRelTol * D0))
            {
               ContactResult accepted = raw;
               accepted.overlap       = std::min(D0, q.value);
               accepted.point         = detail::refinedContactPoint(ci, cj, L0, accepted.overlap, raw.point);
               return done(accepted, RefinePath::FastAccepted);
            }
            tr.fastReject = FastReject::SmallerElsewhere;
         }
         else
            tr.fastReject = FastReject::NotLocallyOptimal;
      }
      else
         tr.fastReject = FastReject::DepthMismatch;
      // otherwise descend from the closed-form normal and the closest-axis-points direction only
      Vec3 std5[5];
      detail::cylinderSeedDirections(ci, cj, std5);
      seeds[n++] = L0;
      seeds[n++] = std5[0];
   }
   else
      n += detail::cylinderSeedDirections(ci, cj, seeds + n);

   detail::OverlapMinimum m = detail::minimiseCylinderProjectionOverlap(ci, cj, seeds, n, !raw.contact, 40, real_t(1e-9));
   if (fastQ.value < m.value) m = fastQ; // every D(L) bounds the depth from above: keep the smallest one found
   ContactResult    r;
   const RefinePath sepPath = raw.contact ? RefinePath::ContactDescendedSeparated : RefinePath::NoScenarioSeparated;
   if (m.value <= real_t(0)) return done(r, sepPath); // separated (or exactly touching)
   // the non-smooth candidates (band circles, end-face normals), where the descent stalls: searched when the descent
   // ended with a remaining tangential gradient (it stopped on a kink, not at a smooth minimum) and always for pairs
   // without a scenario
   const real_t lenScale = ci.radius + cj.radius + real_t(0.5) * (ci.length + cj.length);
   const Vec3   gm = detail::cylinderSupportPoint(ci, m.direction) - detail::cylinderSupportPoint(cj, -m.direction);
   const bool   stalled = (gm - (gm * m.direction) * m.direction).length() > real_t(1e-7) * lenScale;
   if (!raw.contact || stalled)
   {
      tr.kinks                       = true;
      const detail::OverlapMinimum k = detail::minimiseOverlapOnKinks(ci, cj);
      if (k.value < m.value) m = k;
      if (m.value <= real_t(0)) return done(r, sepPath);
   }
   // no scenario and a positive local minimum: run the convex (global) check, which either proves separation or finds a
   // smaller overlap than the local minimum; in the latter case polish it with the descent from the new direction
   if (!raw.contact)
   {
      tr.convex                      = true;
      const detail::OverlapMinimum c = detail::cylindersConvexOverlapCheck(ci, cj, m.direction, m.value);
      if (c.value <= real_t(0)) return done(r, sepPath);
      if (c.value < m.value)
      {
         const detail::OverlapMinimum p = detail::minimiseCylinderProjectionOverlap(ci, cj, &c.direction, 1, true);
         m = (p.value < c.value) ? p : c;
         if (m.value <= real_t(0)) return done(r, sepPath);
      }
   }
   r.contact  = true;
   r.overlap  = m.value;
   r.normal   = -m.direction;
   // Contact point: midpoint of the overlap along L from the body with the unique support point
   // (detail::refinedContactPoint)
   r.point = detail::refinedContactPoint(ci, cj, m.direction, m.value, raw.contact ? raw.point : r.point);
   r.scenario = raw.contact ? raw.scenario : "refined-no-scenario";
   return done(r, raw.contact ? RefinePath::ContactDescended : RefinePath::NoScenarioContact);
}

//! Closed form followed by the refinement -- the complete detection without any GJK/EPA call.
inline ContactResult cylinderCylinderContactRefined(const CylinderState& ci, const CylinderState& cj)
{
   return refineCylinderContact(ci, cj, cylinderCylinderContact(ci, cj));
}

} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla
