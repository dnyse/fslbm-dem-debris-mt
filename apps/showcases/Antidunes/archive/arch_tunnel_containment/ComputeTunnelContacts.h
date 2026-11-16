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
//! \file ComputeTunnelContacts.h
//! \brief ARCHIVED: the full arch-tunnel/bridge containment feature (ArchTunnelParams and its
//! geometry helpers, originally SimulationParams.h; createArchTunnel, originally Utility.h;
//! ComputeTunnelContacts, originally ParticleKernels.h) for the real-scale Pont Ferroviere bridge
//! case study (chapters/03_implementation/case-study/, chapters/04_results/case-study/ --
//! intentionally excluded from the thesis, see CLAUDE.md). Archived in full 2026-11-16 at the
//! author's request, completing the fuller cleanup flagged (but deliberately deferred) by the
//! original 2026-10-29 arch_tunnel_containment/ archival -- see archive/ARCHIVE_README.md.
//!
//! NOT standalone-compilable in isolation, but stays compilable by hand alongside the live
//! showcase: includes ParticleKernels.h (for FloaterHistoryLogger and
//! closed_form_contacts::detail::farthestThetaOnCircleToLine, transitively via ClosedFormContacts.h)
//! and Utility.h (for the mesa_pd body-creation helpers createArchTunnel itself follows the pattern
//! of, e.g. createCylinder).
#pragma once

#include "../../ParticleKernels.h"
#include "../../Utility.h"

namespace walberla
{
namespace antidunes
{

// An arch's tunnel opening: a finite hollow-cylinder interior a floater can pass through, bounded to
// the arch's own axial (X) span so it never affects anything outside that span. Unlike CylindricalBoundary
// (infinite along its axis, used for solid pillars), this is not a mesa_pd body at all -- it's consumed
// directly by ParticleKernels.h's ComputeTunnelContacts. center's X is ignored (the tunnel's axis is
// always world-X); only Y/Z (the tunnel's cross-section centre, matching the arch's BodyBC cylinder)
// matter.
struct ArchTunnelParams {
   Vector3< real_t > center;
   real_t            radius;
   real_t            xMin;
   real_t            xMax;

   // Optional solid-geometry extensions for useArchTunnelAsMesaPDObject's deck-fill/haunch DEM
   // bodies (SetupStages.h::createArchSolidElements, called from setupParticles) -- NOT used by
   // ComputeTunnelContacts or the curved containment body, only by those two additional solid
   // pieces. All three default to reproducing "nothing extra" (no deck-fill, no haunch) so every
   // existing ArchTunnels block that doesn't set them keeps its exact prior behaviour:
   //   - yMin/yMax default to center.y -+ radius, i.e. the arch's own box-span exactly matches the
   //     circular opening's own diameter (true for every arch in this project's bridge scenarios
   //     today -- confirmed directly against their own BodyBC Y-spans) -- so the haunch boxes
   //     Utility.h builds from these (constant-width = radius, safely never encroaching into the
   //     true bore -- see that function's own comment for the derivation) come out zero-width and
   //     are skipped. Set yMin < center.y-radius / yMax > center.y+radius for a scenario whose
   //     abutment is wider than its own opening, to get real, correctly-sized haunch boxes.
   //   - deckTopZ defaults to center.z+radius (the crown height itself), giving a zero-height (and
   //     therefore skipped) deck-fill box. Set it to the real deck top (matching the fluid-side
   //     BodyBC arch box's own max Z) to get a real solid slab from the crown up to the deck.
   real_t            yMin;
   real_t            yMax;
   real_t            deckTopZ;

   ArchTunnelParams(const config::Config::BlockHandle& b)
   {
      center   = b.getParameter< Vector3< real_t > >("center");
      radius   = b.getParameter< real_t >           ("radius");
      xMin     = b.getParameter< real_t >           ("xMin");
      xMax     = b.getParameter< real_t >           ("xMax");
      yMin     = b.getParameter< real_t >           ("yMin", center[1] - radius);
      yMax     = b.getParameter< real_t >           ("yMax", center[1] + radius);
      deckTopZ = b.getParameter< real_t >           ("deckTopZ", center[2] + radius);
   }

   // Direct-value constructor for callers with no prm block to parse from (e.g. DEMValidationScenarios.h's
   // fluid-free local validation harness, which builds scenarios entirely in code). yMin/yMax/deckTopZ
   // default the same way the prm-parsing constructor does -- see that constructor's own comment.
   ArchTunnelParams(const Vector3< real_t >& c, real_t r, real_t xmn, real_t xmx,
      real_t ymn = std::numeric_limits< real_t >::quiet_NaN(),
      real_t ymx = std::numeric_limits< real_t >::quiet_NaN(),
      real_t deckZ = std::numeric_limits< real_t >::quiet_NaN())
      : center(c), radius(r), xMin(xmn), xMax(xmx)
      , yMin(std::isnan(ymn) ? c[1] - r : ymn)
      , yMax(std::isnan(ymx) ? c[1] + r : ymx)
      , deckTopZ(std::isnan(deckZ) ? c[2] + r : deckZ)
   {}
};

// Pure center/edgeLength geometry for the deck-fill and haunch solid pieces described in
// ArchTunnelParams's own comment above -- factored out of SetupStages.h::createArchSolidElements
// specifically so this placement math is independently unit-testable (ContactDetectionTest.cpp)
// without needing a live mesa_pd::data::ParticleStorage. Returns std::nullopt for a degenerate
// (zero-or-negative-size) piece -- exactly the "nothing extra" backward-compatible default
// described on ArchTunnelParams, so callers can skip createBox() outright rather than build a
// zero-volume body.
struct ArchSolidBox {
   Vector3< real_t > center;
   Vector3< real_t > edgeLength;
};

inline std::optional< ArchSolidBox > archDeckFillBoxGeometry(const ArchTunnelParams& tunnel)
{
   const real_t crownZ = tunnel.center[2] + tunnel.radius;
   if (tunnel.deckTopZ <= crownZ + real_t(1e-6)) return std::nullopt;

   const real_t midX = real_t(0.5) * (tunnel.xMin + tunnel.xMax);
   const real_t lenX = tunnel.xMax - tunnel.xMin;
   return ArchSolidBox{ Vector3< real_t >(midX, tunnel.center[1], real_t(0.5) * (crownZ + tunnel.deckTopZ)),
                         Vector3< real_t >(lenX, tunnel.yMax - tunnel.yMin, tunnel.deckTopZ - crownZ) };
}

// Pure dims->ArchTunnelParams conversion for the dimension-based createArchBlock builder
// (SetupStages.h): X/Y/Z extents + bore radius, rather than requiring the caller to pre-compute
// xMin/xMax/yMin/yMax/deckTopZ by hand. `base` is the bottom CENTER of the block -- (X center, tube
// cross-section Y center, springing-line Z) -- matching ArchTunnelParams.center's own convention
// (its X is ignored; the tube axis is always world-X). Kept separate from particle creation (same
// reasoning as archDeckFillBoxGeometry/archHaunchBoxGeometry above) so the dims->params arithmetic is
// unit-testable without a live mesa_pd::data::ParticleStorage. Does NOT validate yDim>=2*radius /
// zDim>=radius -- SetupStages.h::createArchBlock does that before calling this, since a pure function
// like this one has no good way to fail loudly (WALBERLA_CHECK aborts the whole run, which a test
// exploring invalid inputs on purpose would rather avoid).
inline ArchTunnelParams archBlockParams(const Vector3< real_t >& base, real_t xDim, real_t yDim, real_t zDim, real_t radius)
{
   const real_t xHalf = real_t(0.5) * xDim;
   const real_t yHalf = real_t(0.5) * yDim;
   return ArchTunnelParams(base, radius, base[0] - xHalf, base[0] + xHalf,
                            base[1] - yHalf, base[1] + yHalf, base[2] + zDim);
}

enum class ArchHaunchSide { Left, Right };

// Left = the -Y side of the tunnel's own opening (toward yMin); Right = the +Y side (toward yMax).
// Both sit beside the curve, spanning the same [center.z, center.z+radius] band as the curve's own
// quarter-height (a conservative constant-width slab -- see ArchTunnelParams's own comment for why
// this never encroaches into the true bore).
inline std::optional< ArchSolidBox > archHaunchBoxGeometry(const ArchTunnelParams& tunnel, ArchHaunchSide side)
{
   const real_t midX = real_t(0.5) * (tunnel.xMin + tunnel.xMax);
   const real_t lenX = tunnel.xMax - tunnel.xMin;
   const real_t midZ = tunnel.center[2] + real_t(0.5) * tunnel.radius;

   if (side == ArchHaunchSide::Left)
   {
      const real_t width = (tunnel.center[1] - tunnel.radius) - tunnel.yMin;
      if (width <= real_t(1e-6)) return std::nullopt;
      return ArchSolidBox{ Vector3< real_t >(midX, tunnel.yMin + real_t(0.5) * width, midZ),
                            Vector3< real_t >(lenX, width, tunnel.radius) };
   }
   else
   {
      const real_t width = tunnel.yMax - (tunnel.center[1] + tunnel.radius);
      if (width <= real_t(1e-6)) return std::nullopt;
      return ArchSolidBox{ Vector3< real_t >(midX, tunnel.yMax - real_t(0.5) * width, midZ),
                            Vector3< real_t >(lenX, width, tunnel.radius) };
   }
}

inline std::vector< ArchTunnelParams > parseArchTunnels(const shared_ptr< config::Config >& cfg)
{
   std::vector< ArchTunnelParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("ArchTunnels", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks tunnelBlocks;
   containers[0].getBlocks("ArchTunnel", tunnelBlocks, 0, 100);
   for (const auto& b : tunnelBlocks)
      result.emplace_back(b);
   return result;
}

// Arch tunnel as a REAL mesa_pd body (ModelParams.useArchTunnelAsMesaPDObject -- see that flag's own
// doc comment): a CylindricalBoundary in CONTAINMENT mode (SolidCylinderContactDetection's new
// containment overloads, ParticleKernels.h), axially bounded to the tunnel's own [xMin,xMax] span and
// half-clipped to only the semicircle above its own springing line (a real arch is solid material only
// from the springing line up to the crown -- see ComputeTunnelContacts' own class comment for the same
// physical fact, ported here as a general CylindricalBoundary field instead of Antidunes-specific code).
// Same FIXED/INFINITE/GLOBAL/NON_COMMUNICATING pattern as createCylinder (piers) above -- infinite mass,
// always a broad-phase candidate, narrow-phase math (axial + half-clip gates) does the real cutoff.
inline auto createArchTunnel(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& center,
                             real_t radius, real_t xMin, real_t xMax)
{
   const real_t axisX = real_t(0.5) * (xMin + xMax);
   auto p0 = ps.create(true);
   p0->setPosition(mesa_pd::Vec3(axisX, center[1], center[2]));
   auto shape = std::make_shared< mesa_pd::data::CylindricalBoundary >(radius, mesa_pd::Vec3(1, 0, 0));
   shape->setContainment(true);
   shape->setAxialBounds(xMin - axisX, xMax - axisX);
   shape->setHalfClipNormal(mesa_pd::Vec3(0, 0, 1)); // solid only at/above the springing line (center[2])
   p0->setBaseShape(shape);
   p0->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(2);
   p0->setInteractionRadius(std::numeric_limits< real_t >::infinity());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}

// An arch tunnel is a finite hollow-cylinder INTERIOR a floater can pass through -- the opposite of
// every other contact in this file, which is a solid a floater collides with from the OUTSIDE
// (CylindricalBoundary pillars, HalfSpace walls). Not a mesa_pd body at all: DoubleCast is a
// code-generated closed switch over 7 hardcoded shape types, and CylindricalBoundary (SHAPE_TYPE==2)
// is already claimed by the pillars' outside-solid overloads above, so the same C++ shape type cannot
// get a second, inside-out interpretation through overload resolution. Instead this iterates active
// floaters directly (forEachParticle, not forEachParticlePairHalf/DoubleCast) against a small
// in-memory list of tunnel definitions (ArchTunnelParams — center Y/Z, radius, axial X-range).
//
// mesa_pd::kernel::LinearSpringDashpot::operator() cannot be called one-sided (it dereferences a
// SECOND real particle index for type/uid/contact-history/velocity) -- see plan Correction A. This
// class uses LinearSpringDashpot only for its public stiffness/damping GETTERS (never operator()),
// and applies the resulting force directly via mesa_pd::addForceAtWFPosAtomic -- no second particle,
// no dummy/sentinel body.
//
// Sphere and Box floaters use their TRUE farthest point from the tunnel axis (see below), not an
// isotropic circumscribing-sphere bound. This used to use ac.getInteractionRadius(idx) (the
// project-wide circumscribing-sphere bound, e.g. a Box floater's half-diagonal) for every shape alike
// -- exact for a sphere, but a real bug for an elongated box: a log with halfExtents<6,1,1.5> has
// interactionRadius = sqrt(6²+1²+1.5²) = 6.26, dominated entirely by its 12-LU LENGTH, even though its
// true radial (Y-Z, perpendicular-to-the-tunnel-axis) reach is only sqrt(1²+1.5²) = 1.8. Since
// 6.26 > every arch's own radius (5/10), that meant ANY log whose nose entered a tunnel's axial
// [xMin,xMax] span registered a near-maximal "overlap" (dist-to-axis + 6.26 - radius) REGARDLESS of
// its actual Y/Z position -- confirmed against a real run (job 4134552): floater
// uid=144115188075855873 hit [TUNNEL-ENVELOP] with dist-to-axis=9.88 (true gap to the radius-5 wall
// ~3 LU) purely because of this, producing a single-subcycle contactForce_z of ~0.64 against a normal
// buoyancy force of ~0.0003 -- the "log flies out of the water" event. Cylinder floaters (unused by
// any current scenario) still fall back to the isotropic bound below: an exact farthest-point-from-a-
// line calc for a solid cylinder needs a continuous optimisation over its rim circles (not a discrete
// vertex set), not worth the complexity until a scenario actually exercises it.
//
// The fix is the MIRROR IMAGE of the pillars' Box/CylindricalBoundary CLOSEST-point clamp above, but
// easier: containment needs the box's FARTHEST point from the axis (does any part of it poke outside
// the tunnel wall?), and radial-distance-to-a-line is a CONVEX function of the box's body-frame
// point -- a convex function maximized over a convex polytope (the box) always attains its max at a
// VERTEX, so the exact answer is just "check all 8 corners, keep the largest," no iteration needed
// (unlike the pillars' MINIMUM-distance search, which genuinely can land in the box's interior and
// needs the iterative project-and-clamp).
//
// Friction is similarly simplified relative to LinearSpringDashpot's full stick-slip tangential
// spring (which needs a persistent per-contact history keyed against a second real particle): plain
// Coulomb-capped viscous drag instead (proportional to tangential slip velocity, capped at
// mu_dynamic * fN — the "always slipping" limit of the full model). The physically essential claim
// here is containment (fN keeps the floater inside the tunnel); faithful stick/slip fidelity is not.
class ComputeTunnelContacts {
public:
   ComputeTunnelContacts(std::vector< ArchTunnelParams > tunnels, real_t restitutionCoefficient,
      real_t kappa, real_t collisionTime, real_t frictionCoefficient)
      : tunnels_(std::move(tunnels))
      , stiffnessCalc_(1)
      , restitutionCoefficient_(restitutionCoefficient)
      , kappa_(kappa)
      , collisionTime_(collisionTime)
   {
      stiffnessCalc_.setFrictionCoefficientDynamic(0, 0, frictionCoefficient);
   }

   // Wires ModelParameters.logFloaterHistory's CSV writer in -- same rationale/pattern as
   // ComputeContacts::setHistoryLogger/setStepInfo, except this is a one-sided kernel (no second real
   // particle), so it calls FloaterHistoryLogger::logTunnelContact instead of log().
   void setHistoryLogger(FloaterHistoryLogger* logger) { historyLogger_ = logger; }
   void setStepInfo(uint_t mainStep, uint_t subCycle) { mainStep_ = mainStep; subCycle_ = subCycle; }

   template< typename ParticleAccessor_T >
   void operator()(const size_t idx, ParticleAccessor_T& ac)
   {
      const uint_t t = ac.getType(idx);
      if (t != uint_t(1) && t != uint_t(3) && t != uint_t(4)) return; // floaters only

      const mesa_pd::Vec3 pos  = ac.getPosition(idx);
      const real_t        effR = ac.getInteractionRadius(idx);

      // Exact world-frame corner set for a Box floater -- see class comment for the "farthest point
      // from a line is a convex-over-polytope maximum, so it's exactly one of the 8 corners"
      // derivation. Computed once per floater, reused for every tunnel below (the corners themselves
      // don't depend on which tunnel we're testing against).
      mesa_pd::Vec3 corners[8];
      bool          haveCorners = false;
      if (t == uint_t(3) && ac.getBaseShape(idx) &&
          ac.getBaseShape(idx)->getShapeType() == mesa_pd::data::Box::SHAPE_TYPE)
      {
         const auto&          box = *static_cast< const mesa_pd::data::Box* >(ac.getBaseShape(idx).get());
         const mesa_pd::Vec3  h   = box.getEdgeLength() * real_t(0.5); // getEdgeLength() is FULL length
         const mesa_pd::Mat3& R   = ac.getRotation(idx).getMatrix();
         uint_t c = 0;
         for (real_t sx : { real_t(-1), real_t(1) })
            for (real_t sy : { real_t(-1), real_t(1) })
               for (real_t sz : { real_t(-1), real_t(1) })
                  corners[c++] = pos + R * mesa_pd::Vec3(sx * h[0], sy * h[1], sz * h[2]);
         haveCorners = true;
      }

      for (size_t ti = 0; ti < tunnels_.size(); ++ti)
      {
         const auto& tunnel = tunnels_[ti];

         real_t        farthestRadial; // exact (box) or conservative (sphere/cylinder) outward reach
         real_t        xLo, xHi;       // world-X extent, for the axial gate below
         mesa_pd::Vec3 worstPoint;     // the actual floater-surface point driving farthestRadial

         // SPRINGING-LINE CONSTRAINT: a real arch is only a SEMICIRCLE of solid material -- the LBM
         // BodyBC builds it as LogicalAND(deck box, LogicalNOT(cylinder)), which only carves/bounds the
         // opening from the springing line (tunnel.center[2], by construction — this Antidunes bridge
         // always centers its ArchTunnel on the springing line, e.g. Z=56) up to the crown
         // (center[2]+radius). BELOW the springing line there is no solid at all -- that's the open
         // water passage between the piers. ArchTunnelParams itself has no such restriction (it's a
         // full circle, center[2]±radius both ways), so without this constraint ComputeTunnelContacts
         // was treating that entire lower, nonexistent half as solid tunnel wall -- confirmed against
         // real data: 171226 of 176345 logged tunnel-contact rows (97%) from job 4134746 had
         // contactZ < 56, i.e. were firing against open water, not real structure. Below, a Box corner
         // or a Sphere/Cylinder's own top point below the springing line simply can't have hit anything
         // real, so it's excluded from the farthest-point search entirely (not merely a coarser bound).
         if (haveCorners)
         {
            farthestRadial = real_t(-1);
            xLo = corners[0][0];
            xHi = corners[0][0];
            for (const auto& corner : corners)
            {
               xLo = std::min(xLo, corner[0]);
               xHi = std::max(xHi, corner[0]);
               if (corner[2] < tunnel.center[2]) continue; // below the springing line -- no real wall here
               const real_t cdy = corner[1] - tunnel.center[1];
               const real_t cdz = corner[2] - tunnel.center[2];
               const real_t r   = std::sqrt(cdy * cdy + cdz * cdz);
               if (r > farthestRadial) { farthestRadial = r; worstPoint = corner; }
            }
            // No corner reached above the springing line at all -- the whole floater is in the open
            // water passage under this tunnel, nothing to contain it against. farthestRadial stays -1,
            // which safely fails the overlap check below regardless of tunnel.radius, but skip explicitly
            // for clarity.
            if (farthestRadial < real_t(0)) continue;
         }
         else
         {
            if (pos[2] + effR < tunnel.center[2]) continue; // wholly below the springing line
            xLo = pos[0] - effR;
            xHi = pos[0] + effR;

            if (t == uint_t(1))
            {
               // Sphere: exact and springing-line-restricted, matching the box's per-corner
               // rigor. CONFIRMED NECESSARY, not just theoretical: job 4135905 (raised-water
               // sphere run) showed 87.3% of logged tunnel-contact rows below the springing
               // line before this fix -- the naive "unconstrained centre+effR outward point"
               // frequently lands below Z=springing-line for a sphere whose CENTRE is below it
               // even when the sphere's top pokes above (the only thing the top-point gate
               // above actually rejects).
               //
               // Distance-to-the-axis is a strictly convex function of angle around the
               // sphere's own Y-Z great circle (center=pos, radius=effR): it has exactly one
               // global max and one global min (diametrically opposite), and is monotonic
               // between them each way around. Restricting to the springing-line arc
               // (Z >= tunnel.center[2], itself a single contiguous arc since Z = pos.z +
               // effR*sin(theta) >= threshold defines one interval of theta per period): if the
               // UNCONSTRAINED global max already lies in that arc, it's still the answer (no
               // change from before). Otherwise the arc lies within a single monotonic piece
               // (whether or not it also contains the global min -- either way the function is
               // "high at both ends, monotonic or valley-shaped in between"), so the
               // CONSTRAINED max is provably at one of the arc's own two endpoints: the two
               // points where the sphere's great circle crosses the springing line.
               const real_t dy   = pos[1] - tunnel.center[1];
               const real_t dz   = pos[2] - tunnel.center[2];
               const real_t dist = std::sqrt(dy * dy + dz * dz);
               mesa_pd::Vec3 candidate = pos;
               if (dist > real_t(1e-10))
                  candidate = pos + (mesa_pd::Vec3(real_t(0), dy, dz) / dist) * effR;
               else
                  candidate = pos + mesa_pd::Vec3(real_t(0), effR, real_t(0));

               if (candidate[2] < tunnel.center[2])
               {
                  const real_t dz0   = tunnel.center[2] - pos[2]; // pos.z -> springing line
                  const real_t under = effR * effR - dz0 * dz0;   // >= 0 guaranteed by the top-point gate
                  const real_t halfChord = std::sqrt(std::max(under, real_t(0)));
                  const mesa_pd::Vec3 p1(pos[0], pos[1] + halfChord, tunnel.center[2]);
                  const mesa_pd::Vec3 p2(pos[0], pos[1] - halfChord, tunnel.center[2]);
                  const real_t d1 = std::abs(p1[1] - tunnel.center[1]);
                  const real_t d2 = std::abs(p2[1] - tunnel.center[1]);
                  candidate = (d1 > d2) ? p1 : p2;
               }
               worstPoint = candidate;
               const real_t wdy = worstPoint[1] - tunnel.center[1];
               const real_t wdz = worstPoint[2] - tunnel.center[2];
               farthestRadial    = std::sqrt(wdy * wdy + wdz * wdz);
            }
            else
            {
               // Cylinder: exact farthest point from the tunnel axis, springing-line-restricted,
               // matching the box's per-corner and the sphere's arc-endpoint rigor. A solid
               // cylinder's farthest point from ANY line is always on one of its two end-cap RIMS,
               // never on the interior of its lateral surface or its length: for a fixed angle theta
               // around the rim, the squared radial distance to the tunnel axis is a convex quadratic
               // in the axial parameter s (Py,Pz are affine in s, and squared-distance-to-a-line is a
               // positive-semidefinite quadratic form of an affine function), so its max over
               // s in [-h,h] is always at an endpoint (s=+h or s=-h) -- i.e. one of the two rims --
               // exactly the same "convex function maximized over a convex body attains its max at an
               // extreme point" argument used for the box's corners above, just applied per-angle
               // before extremizing over angle too.
               //
               // Each rim is then its OWN circle-vs-line farthest-point problem
               // (closed_form_contacts::detail::farthestThetaOnCircleToLine, the maximizing counterpart
               // of the closest-point search already used elsewhere in this file), restricted to the
               // springing-line arc the same way the sphere fix is: if the rim's unconstrained
               // farthest point already satisfies Z >= springing line, it's exact and unchanged; if
               // not, fall back to the two points where the rim crosses the springing line itself
               // (same simplification precedent as the sphere fix above -- provably exact for a
               // single-extremum function; here f(theta) can in principle have a second local
               // maximum, so this is a slight, deliberately-conservative approximation rather than a
               // proof, but it is already a large improvement over the isotropic dist+effR bound this
               // replaces, which had no springing-line awareness or radial precision at all).
               if (!ac.getBaseShape(idx) ||
                   ac.getBaseShape(idx)->getShapeType() != mesa_pd::data::Cylinder::SHAPE_TYPE)
               {
                  // Defensive fallback only -- t==4 should always carry a real Cylinder shape; this
                  // just avoids an invalid cast if that invariant is ever violated, using the same
                  // isotropic bound this whole branch otherwise replaces.
                  const real_t dy   = pos[1] - tunnel.center[1];
                  const real_t dz   = pos[2] - tunnel.center[2];
                  const real_t dist = std::sqrt(dy * dy + dz * dz);
                  farthestRadial    = dist + effR;
                  const mesa_pd::Vec3 outward = (dist > real_t(1e-10))
                                                   ? mesa_pd::Vec3(real_t(0), dy, dz) / dist
                                                   : mesa_pd::Vec3(real_t(0), real_t(1), real_t(0));
                  worstPoint = pos + outward * effR;
               }
               else
               {
               const auto& cyl = *static_cast< const mesa_pd::data::Cylinder* >(ac.getBaseShape(idx).get());
               const real_t         cylR    = cyl.getRadius();
               const real_t         cylH    = real_t(0.5) * cyl.getHeight();
               const mesa_pd::Vec3  axis    = ac.getRotation(idx).getMatrix() * mesa_pd::Vec3(real_t(0), real_t(0), real_t(1));
               const mesa_pd::Vec3  tLinePos(real_t(0), tunnel.center[1], tunnel.center[2]);
               const mesa_pd::Vec3  tLineDir(real_t(1), real_t(0), real_t(0));

               auto perpBasis = [](const mesa_pd::Vec3& ax, mesa_pd::Vec3& u, mesa_pd::Vec3& v) {
                  const mesa_pd::Vec3 helper = (std::abs(ax[2]) < real_t(0.9))
                                                   ? mesa_pd::Vec3(real_t(0), real_t(0), real_t(1))
                                                   : mesa_pd::Vec3(real_t(1), real_t(0), real_t(0));
                  u = (ax % helper).getNormalized();
                  v = (ax % u).getNormalized();
               };
               mesa_pd::Vec3 uB, vB;
               perpBasis(axis, uB, vB);

               farthestRadial = real_t(-1);
               for (const real_t sgn : { real_t(-1), real_t(1) })
               {
                  const mesa_pd::Vec3 rimC = pos + sgn * cylH * axis;

                  auto pointAt = [&](real_t theta) {
                     return rimC + cylR * std::cos(theta) * uB + cylR * std::sin(theta) * vB;
                  };
                  auto radialSq = [&](const mesa_pd::Vec3& P) {
                     const real_t ddy = P[1] - tunnel.center[1];
                     const real_t ddz = P[2] - tunnel.center[2];
                     return ddy * ddy + ddz * ddz;
                  };

                  real_t thetaBest = closed_form_contacts::detail::farthestThetaOnCircleToLine(
                     rimC, uB, vB, cylR, tLinePos, tLineDir);
                  mesa_pd::Vec3 candidate = pointAt(thetaBest);

                  if (candidate[2] < tunnel.center[2])
                  {
                     // Unconstrained farthest point falls below the springing line -- restrict to
                     // the arc's own two endpoints (where the rim crosses Z = springing line),
                     // same closed-form amplitude/phase solve used elsewhere in this file for a
                     // circle-crosses-a-plane-height problem: A1*cos(theta)+A2*sin(theta) = rhs.
                     const real_t A1 = uB[2], A2 = vB[2];
                     const real_t Amag = std::sqrt(A1 * A1 + A2 * A2);
                     if (Amag > real_t(1e-12))
                     {
                        const real_t rhs = (tunnel.center[2] - rimC[2]) / cylR;
                        if (std::abs(rhs) <= Amag)
                        {
                           const real_t phi   = std::atan2(A2, A1);
                           const real_t delta = std::acos(std::clamp(rhs / Amag, real_t(-1), real_t(1)));
                           const mesa_pd::Vec3 e1 = pointAt(phi + delta);
                           const mesa_pd::Vec3 e2 = pointAt(phi - delta);
                           candidate = (radialSq(e1) > radialSq(e2)) ? e1 : e2;
                        }
                        else { continue; } // whole rim below the springing line -- no real wall here
                     }
                     else { continue; } // rim's plane parallel to Z -- degenerate, skip (no crossing)
                  }

                  const real_t r = std::sqrt(radialSq(candidate));
                  if (r > farthestRadial) { farthestRadial = r; worstPoint = candidate; }
               }
               if (farthestRadial < real_t(0)) continue; // both rims entirely below the springing line
               }
            }
         }

         // Axial (world-X) gate: the tunnel's cross-section is constant along its whole length, so a
         // floater whose own X-extent doesn't overlap [xMin,xMax] at all is unaffected by this tunnel,
         // full stop -- no phantom wall outside the arch's own footprint. Now uses the floater's TRUE
         // X-extent (box corners) where available, instead of the isotropic effR bound.
         if (xHi < tunnel.xMin || xLo > tunnel.xMax) continue;

         // Radial gate: mirrors the axial (X) gate above. All arch tunnels in this bridge share the same
         // [xMin,xMax] span (they're all under the same deck), so the X-gate alone lets a floater near
         // ONE tunnel also pass the gate for the OTHER tunnels it isn't anywhere near. Without this check,
         // `overlap` below is unbounded for any floater merely passing through the shared X-span far from
         // this particular tunnel's axis, producing a runaway spring force toward the wrong tunnel's wall
         // (this caused a real crash: the resulting single-substep velocity spike flung a floater far
         // enough to outrun ghost-particle sync, producing a stale-UID heap-buffer-overflow in
         // SimpleBB::treatDirection, caught via ASan on job 4065016). REGRESSION FIX: the exact-corner
         // rewrite above dropped the actual gate check (only kept this comment!) -- confirmed as a real
         // bug via job 4134671, where uid=144115188075855873 sitting at Y=23.6,Z=46.1 (nowhere near any
         // of the three tunnel axes) still tripped [TUNNEL-ENVELOP] against all three simultaneously,
         // since farthestRadial grows without bound as the floater moves away and nothing rejected it
         // first. dist-effR is a valid LOWER bound on the floater's true minimum distance to the axis
         // for ANY shape (effR is always an upper bound on the floater's own reach from its centre), so
         // this reuses the original pre-fix gate unchanged -- cheap, and correct regardless of anisotropy.
         const real_t dyCenter   = pos[1] - tunnel.center[1];
         const real_t dzCenter   = pos[2] - tunnel.center[2];
         const real_t distCenter = std::sqrt(dyCenter * dyCenter + dzCenter * dzCenter);
         if (distCenter - effR > tunnel.radius) continue; // nowhere near this tunnel's wall

         real_t overlap = farthestRadial - tunnel.radius; // > 0 → floater's true outer reach crosses the wall
         if (overlap <= real_t(0)) continue;              // safely inside the tunnel bore

         // Absolute safety clamp regardless of floater size -- mirrors the pillar overloads' own
         // envelope clamp above. Now backed by an exact (box) or still-conservative (sphere/cylinder)
         // overlap, so this should trigger far less often than before; kept as a genuine safety net for
         // real fast single-substep impacts, not a routine crutch for the geometry bug it used to mask.
         const real_t maxOverlap = real_t(1.5);
         if (overlap > maxOverlap)
         {
            static uint_t envelopLogCount = uint_t(0);
            if (envelopLogCount < uint_t(20))
            {
               WALBERLA_LOG_WARNING("[TUNNEL-ENVELOP] floater uid=" << ac.getUid(idx)
                  << " deeply overlaps arch tunnel wall (raw overlap=" << overlap
                  << ", farthest-point radial dist=" << farthestRadial << ") — clamping overlap to "
                  << maxOverlap << "; likely a fast single-substep impact rather than a converged contact.");
               ++envelopLogCount;
            }
            overlap = maxOverlap;
         }

         const real_t worstDy = worstPoint[1] - tunnel.center[1];
         const real_t worstDz = worstPoint[2] - tunnel.center[2];
         const real_t worstDist = std::sqrt(worstDy * worstDy + worstDz * worstDz);
         const mesa_pd::Vec3 outward = (worstDist > real_t(1e-10))
                                          ? mesa_pd::Vec3(real_t(0), worstDy, worstDz) / worstDist
                                          : mesa_pd::Vec3(real_t(0), real_t(1), real_t(0));
         const mesa_pd::Vec3 normal       = -outward; // points inward, toward the axis (opposite of a
                                                       // solid pillar's outward-pointing normal)
         const mesa_pd::Vec3 contactPoint = worstPoint; // the actual (box corner, or sphere-surface) point

         // Effective mass: wall side is infinite mass (same physical assumption CylindricalBoundary/
         // HalfSpace pillar contacts already make), so 1/(invMass_floater + 0) = the floater's own mass.
         const real_t meff = real_t(1) / ac.getInvMass(idx);
         stiffnessCalc_.setStiffnessAndDamping(0, 0, restitutionCoefficient_, collisionTime_, kappa_, meff);
         const real_t stiffnessN = stiffnessCalc_.getStiffnessN(0, 0);
         const real_t dampingN   = stiffnessCalc_.getDampingN(0, 0);
         const real_t dampingT   = stiffnessCalc_.getDampingT(0, 0);
         const real_t muDynamic  = stiffnessCalc_.getFrictionCoefficientDynamic(0, 0);

         const mesa_pd::Vec3 velAtContact = mesa_pd::getVelocityAtWFPoint(idx, ac, contactPoint);
         const real_t        velN         = velAtContact * normal;
         const mesa_pd::Vec3 velT         = velAtContact - velN * normal;

         // Normal force: same linear spring-dashpot formula LinearSpringDashpot::operator() uses
         // (stiffnessN*overlap + dampingN*closing-velocity), one-sided.
         real_t fN = stiffnessN * overlap - dampingN * velN;
         if (fN < real_t(0)) fN = real_t(0); // a spring-dashpot must never pull the floater toward the wall

         mesa_pd::Vec3 fT = -dampingT * velT;
         const real_t  fTmag = fT.length();
         const real_t  fTmax = muDynamic * fN;
         if (fTmag > fTmax && fTmag > real_t(1e-12)) fT *= (fTmax / fTmag);

         const bool           doLog       = historyLogger_ && historyLogger_->enabled();
         const mesa_pd::Vec3  forceBefore = doLog ? ac.getForceRef(idx) : mesa_pd::Vec3(real_t(0));

         mesa_pd::addForceAtWFPosAtomic(idx, ac, fN * normal + fT, contactPoint);

         if (doLog)
         {
            // penetrationDepth convention matches ComputeContacts/ComputeLubrication's CSV column:
            // negative = penetrating (this class's own `overlap` is the opposite sign, positive =
            // penetrating, so negate it here).
            historyLogger_->logTunnelContact(mainStep_, subCycle_, ac, idx, ti, contactPoint, normal,
               -overlap, ac.getForceRef(idx) - forceBefore);
         }
      }
   }

private:
   std::vector< ArchTunnelParams >      tunnels_;
   mesa_pd::kernel::LinearSpringDashpot stiffnessCalc_; // getters only -- operator() is never called
   real_t restitutionCoefficient_;
   real_t kappa_;
   real_t collisionTime_;

   FloaterHistoryLogger* historyLogger_ = nullptr;
   uint_t mainStep_ = uint_t(0);
   uint_t subCycle_ = uint_t(0);
};

} // namespace antidunes
} // namespace walberla
