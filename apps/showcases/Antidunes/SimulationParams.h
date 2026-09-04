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
//! \file SimulationParams.h
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
// Plain-data parameter structs. Nothing in here owns any waLBerla resources.
//======================================================================================================================

#pragma once

#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "Utility.h"

#include <cmath>
#include <limits>
#include <optional>


namespace walberla {
namespace antidunes {

struct BlockForestParams {
   Vector3< uint_t > cellsPerBlock;
   Vector3< bool >   periodicity;
   bool              loadSnapshot;
   bool              storeSnapshot;
   uint_t            snapshotFrequency;
   std::string       snapshotBaseFolder;

   BlockForestParams(const config::Config::BlockHandle& b)
      : loadSnapshot       { b.getParameter< bool >               ("loadSnapshot") }
      , storeSnapshot      { b.getParameter< bool >               ("storeSnapshot") }
      , snapshotFrequency  { b.getParameter< uint_t >             ("snapshotFrequency") }
      , snapshotBaseFolder { b.getParameter< std::string >        ("snapshotBaseFolder") }
   {
      cellsPerBlock =        b.getParameter< Vector3< uint_t > >  ("cellsPerBlock");
      periodicity =          b.getParameter< Vector3< bool > >    ("periodicity");
   }
};

struct DomainParams {
   Vector3< uint_t > domainSize;
   uint_t            wavePeriods;
   real_t            liquidHeightFactor;
   real_t            floorHeightFactor;
   real_t            initialAmplitude;

   DomainParams(const config::Config::BlockHandle& b)
      : wavePeriods         { b.getParameter< uint_t >            ("wavePeriods") }
      , liquidHeightFactor  { b.getParameter< real_t >            ("liquidHeightFactor") }
      , floorHeightFactor   { b.getParameter< real_t >            ("floorHeightFactor") }
      , initialAmplitude    { b.getParameter< real_t >            ("initialAmplitude") }
   {
      domainSize =            b.getParameter< Vector3< uint_t > > ("domainSize");
   }
};

struct PIDParams {
   real_t targetMeanVelocityMagnitude;
   real_t proportionalGain;
   real_t derivativeGain;
   real_t integralGain;
   real_t maxRamp;
   real_t minActuatingVariable;
   real_t maxActuatingVariable;

   PIDParams(const config::Config::BlockHandle& b)
      : targetMeanVelocityMagnitude { b.getParameter< real_t >("targetMeanVelocityMagnitude") }
      , proportionalGain            { b.getParameter< real_t >("proportionalGain") }
      , derivativeGain              { b.getParameter< real_t >("derivativeGain") }
      , integralGain                { b.getParameter< real_t >("integralGain") }
      , maxRamp                     { b.getParameter< real_t >("maxRamp") }
      , minActuatingVariable        { b.getParameter< real_t >("minActuatingVariable") }
      , maxActuatingVariable        { b.getParameter< real_t >("maxActuatingVariable") }
   {}
};

struct PhysicsParams {
   bool   enableWetting;
   uint_t timesteps;
   real_t Re;
   real_t Fr;
   real_t We;
   // When true, the U-derived LBM units (viscosity, relaxationRate, gravity, surfaceTension) use the
   // fixed reference values instead of being computed from targetMeanVelocityMagnitude. This lets the
   // mean velocity be set to 0 (still water / no flow forcing) while keeping valid physical constants.
   bool   useHardcodedUnits;
   // Optional damping knob: when > 0, sets the LBM relaxation rate directly (and derives the matching
   // viscosity), overriding the value above. Lets a still-water settle test (U=0) be made more viscous
   // — e.g. so a tumbling cube/log damps to rest instead of coasting in near-inviscid water. 0 = off.
   real_t relaxationRateOverride;
   // Optional gravity knob: when > 0, sets the lattice-unit gravity directly, overriding both the
   // hardcoded and U-derived values above. Lets a still-water validation run be matched against a
   // specific reference value from the literature (e.g. Bogner & Rude 2013's g=7.5e-4, used for their
   // Fekken-cube floating-density validation) independent of this project's own Fr-derived antidune-flow
   // scaling, which for the default config lands ~100x smaller (7.02618e-06). 0 = off.
   real_t gravityOverride;

   // Independent of Antidunes.cpp's own WATCH_FREQ/WATCH_DENSE_FROM debug watchers (trackWatchCells,
   // checkClosedInterfaceLayer, logFloaterLinkBalance) -- these three ONLY gate checkPdfNaN's cadence.
   // A NaN blowup is normally only caught once every WATCH_FREQ=1000 steps, by which point the fluid
   // field may already be corrupted across much of the domain (streaming propagates a NaN outward by
   // ~1 cell/step) -- not enough to localize the true origin step/cell. Resume from the last checkpoint
   // before a known crash (loadSnapshot=true) and set denseNanScanEnabled=true with a
   // from/toStep window bracketing the crash to scan every step in just that window, without paying
   // per-step overhead for the rest of the run. Off by default (matches the pre-existing behaviour
   // exactly: checkPdfNaN only runs on the normal watchNow cadence).
   bool   denseNanScanEnabled;
   uint_t denseNanScanFromStep;
   uint_t denseNanScanToStep;

   PhysicsParams(const config::Config::BlockHandle& b)
      : enableWetting          { b.getParameter< bool >  ("enableWetting") }
      , timesteps              { b.getParameter< uint_t >("timesteps") }
      , Re                     { b.getParameter< real_t >("Re") }
      , Fr                     { b.getParameter< real_t >("Fr") }
      , We                     { b.getParameter< real_t >("We") }
      , useHardcodedUnits      { b.getParameter< bool >  ("useHardcodedUnits", false) }
      , relaxationRateOverride { b.getParameter< real_t >("relaxationRateOverride", real_c(0)) }
      , gravityOverride        { b.getParameter< real_t >("gravityOverride", real_c(0)) }
      , denseNanScanEnabled    { b.getParameter< bool >  ("denseNanScanEnabled", false) }
      , denseNanScanFromStep   { b.getParameter< uint_t >("denseNanScanFromStep", uint_t(0)) }
      , denseNanScanToStep     { b.getParameter< uint_t >("denseNanScanToStep", timesteps) }
   {}
};

struct ParticleParams {
   std::string inFileName;
   uint_t      bedCopiesInX;
   uint_t      bedCopiesInY;
   real_t      densityRatio;
   real_t      fixingHeightFactor;
   real_t      frictionCoefficient;
   real_t      restitutionCoefficient;
   uint_t      numSubCycles;
   bool        useLubricationCorrection;
   bool        useNoSlipParticles;
   bool        useParticleBed;
   real_t      floaterRingBand;   // sphere V_sub plane-fit: keep interface cells with |dist-r|<band.
                                  // 1.5 = contact-ring restriction (default). A large value (e.g. 1e6)
                                  // reproduces the old full-box sampling, for clean A/B comparison.
   // Optional pair of overrides that together bypass reading inFileName's header for SI-unit scaling
   // (getAvgDiameterScalingFactor, Utility.h) entirely -- avgParticleDiameter and particleScalingFactor
   // are normally derived from an (often unrelated) sediment bed's own particle sizes, even in
   // floater-only scenarios (useParticleBed=false) where that bed is never placed. For those scenarios,
   // liquidHeight/floorHeight/particleFixingHeight (via avgParticleDiameterOverride, in LU) and
   // dx_SI/dt_SI/the gravityOverride SI calibration (via dxSIOverride, in metres/LU) can instead be set
   // directly, independent of whatever inFileName happens to contain. Both must be > 0 to take effect
   // (either alone is ignored and the bed file is read as before) -- see ComputedParams's constructor.
   real_t      avgParticleDiameterOverride;
   real_t      dxSIOverride;

   ParticleParams(const config::Config::BlockHandle& b)
      : inFileName               { b.getParameter< std::string >("inFileName") }
      , bedCopiesInX             { b.getParameter< uint_t >     ("bedCopiesInX") }
      , bedCopiesInY             { b.getParameter< uint_t >     ("bedCopiesInY") }
      , densityRatio             { b.getParameter< real_t >     ("densityRatio") }
      , fixingHeightFactor       { b.getParameter< real_t >     ("fixingHeightFactor") }
      , frictionCoefficient      { b.getParameter< real_t >     ("frictionCoefficient") }
      , restitutionCoefficient   { b.getParameter< real_t >     ("restitutionCoefficient") }
      , numSubCycles             { b.getParameter< uint_t >     ("numSubCycles") }
      , useLubricationCorrection { b.getParameter< bool >       ("useLubricationCorrection") }
      , useNoSlipParticles       { b.getParameter< bool >       ("useNoSlipParticles") }
      , useParticleBed           { b.getParameter< bool >       ("useParticleBed", true) }
      , floaterRingBand          { b.getParameter< real_t >     ("floaterRingBand", real_t(1.5)) }
      , avgParticleDiameterOverride { b.getParameter< real_t >  ("avgParticleDiameterOverride", real_c(0)) }
      , dxSIOverride                { b.getParameter< real_t >  ("dxSIOverride", real_c(0)) }
   {}
};

struct PillarParams {
   Vector3< real_t > center;
   Vector3< real_t > halfExtents;

   PillarParams(const config::Config::BlockHandle& b)
   {
      center      = b.getParameter< Vector3< real_t > >("center");
      halfExtents = b.getParameter< Vector3< real_t > >("halfExtents");
   }
};

struct CylinderPillarParams {
   Vector3< real_t > center;
   real_t            radius;

   CylinderPillarParams(const config::Config::BlockHandle& b)
   {
      center = b.getParameter< Vector3< real_t > >("center");
      radius = b.getParameter< real_t >           ("radius");
   }
};

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

struct FloatingSphereParams {
   Vector3< real_t > position;
   real_t            radius;
   real_t            densityRatio;
   // Number of LBM timesteps the sphere is held fixed before being released to float freely.
   // Use this when starting with the sphere at or above the free surface: the FSLBM needs time
   // to establish a stable fill-level distribution around the contact line before the sphere moves.
   // 0 = release immediately (default, correct for fully submerged start).
   uint_t            settlingSteps;
   // Initial linear velocity (LU/timestep) imparted at creation, e.g. to drive two floaters
   // together for a collision test. Default <0,0,0> = at rest (legacy behaviour). Only meaningful
   // with settlingSteps 0 (a FIXED settling particle can't move regardless of velocity).
   Vector3< real_t > initialVelocity;
   // Same pin semantics as FloatingBoxParams::pinPosition/pinAxes: holds selected TRANSLATIONAL DOFs
   // at `position` every DEM sub-cycle (rotation, if any, stays free — though a sphere's rotation has
   // no visible effect on buoyancy). Use <1,1,0> to let a sphere find its own buoyant depth in place
   // without drifting laterally, e.g. for a side-by-side different-density buoyancy comparison.
   bool              pinPosition;
   Vector3< real_t > pinAxes;
   // Holds selected ROTATIONAL DOFs at zero angular velocity every DEM sub-cycle (world-frame axes,
   // same mask convention as pinAxes): a nonzero pinRotationAxes component zeroes that world axis's
   // angular velocity component every sub-cycle, so the body can never accumulate rotation about that
   // axis; a zero component leaves it genuinely free. Independent of pinPosition/pinAxes. E.g.
   // pinRotationAxes <1,0,1> with pinRotation true allows rotation about Y only (roll — rotating
   // within the world XZ plane) while suppressing pitch (about X) and yaw (about Z). A sphere's own
   // rotation has no visible effect on buoyancy, so this mainly exists for API symmetry with
   // FloatingBoxParams/FloatingCylinderParams.
   bool              pinRotation;
   Vector3< real_t > pinRotationAxes;
   // true → the floater is not "in the simulation" at all until settlingSteps elapses: excluded from
   // fluid mapping AND DEM contact (dormant), not merely held FIXED. Lets several floaters share the
   // exact same spawn position with different settlingSteps, dropped one at a time with no step-0
   // collision/mapping pile-up between overlapping siblings. Requires settlingSteps > 0 (otherwise
   // there is no trigger to ever activate it -- a warning is logged if this is violated).
   bool              startInactive;

   FloatingSphereParams(const config::Config::BlockHandle& b)
   {
      position        = b.getParameter< Vector3< real_t > >("position");
      radius          = b.getParameter< real_t >           ("radius");
      densityRatio    = b.getParameter< real_t >           ("densityRatio");
      settlingSteps   = b.getParameter< uint_t >           ("settlingSteps", uint_t(0));
      initialVelocity = b.getParameter< Vector3< real_t > >("initialVelocity", Vector3< real_t >(real_t(0)));
      pinPosition     = b.getParameter< bool >             ("pinPosition", false);
      pinAxes         = b.getParameter< Vector3< real_t > >("pinAxes", Vector3< real_t >(real_t(1)));
      pinRotation     = b.getParameter< bool >             ("pinRotation", false);
      pinRotationAxes = b.getParameter< Vector3< real_t > >("pinRotationAxes", Vector3< real_t >(real_t(1)));
      startInactive   = b.getParameter< bool >             ("startInactive", false);
      if (startInactive && settlingSteps == uint_t(0))
         WALBERLA_LOG_WARNING("FloatingSphere startInactive=true but settlingSteps=0 -- this floater "
                              "will never be activated (releaseSettledFloatingSpheres no-ops when "
                              "settlingSteps==0).");
   }
};

inline std::vector< PillarParams > parsePillars(const shared_ptr< config::Config >& cfg)
{
   std::vector< PillarParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("Pillars", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks pillarBlocks;
   containers[0].getBlocks("Pillar", pillarBlocks, 0, 100);
   for (const auto& b : pillarBlocks)
      result.emplace_back(b);
   return result;
}

inline std::vector< CylinderPillarParams > parseCylinderPillars(const shared_ptr< config::Config >& cfg)
{
   std::vector< CylinderPillarParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("CylinderPillars", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks pillarBlocks;
   containers[0].getBlocks("CylinderPillar", pillarBlocks, 0, 100);
   for (const auto& b : pillarBlocks)
      result.emplace_back(b);
   return result;
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

inline std::vector< FloatingSphereParams > parseFloatingSpheres(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingSphereParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingSpheres", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks sphereBlocks;
   containers[0].getBlocks("FloatingSphere", sphereBlocks, 0, 10000);
   for (const auto& b : sphereBlocks)
      result.emplace_back(b);
   return result;
}

struct FloatingBoxParams {
   Vector3< real_t > position;
   Vector3< real_t > halfExtents;   // half-size in each dimension; edgeLength = 2 * halfExtents
   real_t            densityRatio;
   uint_t            settlingSteps;      // same semantics as FloatingSphereParams::settlingSteps
   real_t            initialTiltAngleDeg; // initial rotation about the Y-axis (degrees); 0 = upright
   real_t            comOffsetZ;          // signed body-frame z-offset of the centre of mass from the
                                          // geometric centroid (LU). Negative = COM BELOW centroid →
                                          // pendulum righting torque (raises GM, prevents tumbling).
                                          // 0 = uniform-density box (COM at centroid).
   bool              pinPosition;         // true → the box's TRANSLATION is held at `position` for the
                                          // whole run while ROTATION stays free (translational DOFs
                                          // pinned each DEM sub-cycle; not the FIXED flag, which would
                                          // also freeze rotation). Use to observe the free-floating
                                          // equilibrium tilt in isolation. false → normal dynamics.
   Vector3< real_t > pinAxes;             // per-axis pin mask (only used when pinPosition true): a nonzero
                                          // component pins that world axis to `position`, a zero leaves it
                                          // FREE. Default <1,1,1> = pin all three (rigid hold). <1,1,0>
                                          // pins x,y but lets z float, so the box finds its own buoyant
                                          // depth (≈50% submerged) while not drifting laterally.
   bool              pinRotation;         // true → selected ROTATIONAL DOFs are held at zero angular
                                          // velocity every DEM sub-cycle (see pinRotationAxes). false →
                                          // rotation fully free (normal dynamics).
   Vector3< real_t > pinRotationAxes;     // per-axis rotation pin mask (only used when pinRotation true):
                                          // a nonzero component zeroes that world axis's angular velocity
                                          // every sub-cycle, so the box can never accumulate rotation
                                          // about that axis; a zero component leaves it genuinely free.
                                          // Default <1,1,1> = pin all rotation (rigid hold, like the old
                                          // capture-and-freeze behaviour). <1,0,1> allows rotation about Y
                                          // only (roll, within the world XZ plane) while suppressing pitch
                                          // (X) and yaw (Z) — e.g. to let a 1x2x1 box roll toward its
                                          // corner-down equilibrium without also picking up an incidental
                                          // yaw from asymmetric contact/splash noise.
   bool              startInactive;      // see FloatingSphereParams::startInactive for the full comment.

   FloatingBoxParams(const config::Config::BlockHandle& b)
   {
      position             = b.getParameter< Vector3< real_t > >("position");
      halfExtents          = b.getParameter< Vector3< real_t > >("halfExtents");
      densityRatio         = b.getParameter< real_t >           ("densityRatio");
      settlingSteps        = b.getParameter< uint_t >           ("settlingSteps", uint_t(0));
      initialTiltAngleDeg  = b.getParameter< real_t >           ("initialTiltAngleDeg", real_t(0));
      comOffsetZ           = b.getParameter< real_t >           ("comOffsetZ", real_t(0));
      pinPosition          = b.getParameter< bool >             ("pinPosition", false);
      pinAxes              = b.getParameter< Vector3< real_t > >("pinAxes", Vector3< real_t >(real_t(1)));
      pinRotation          = b.getParameter< bool >             ("pinRotation", false);
      pinRotationAxes      = b.getParameter< Vector3< real_t > >("pinRotationAxes", Vector3< real_t >(real_t(1)));
      startInactive        = b.getParameter< bool >             ("startInactive", false);
      if (startInactive && settlingSteps == uint_t(0))
         WALBERLA_LOG_WARNING("FloatingBox startInactive=true but settlingSteps=0 -- this floater will "
                              "never be activated (releaseSettledFloatingSpheres no-ops when "
                              "settlingSteps==0).");
   }
};

inline std::vector< FloatingBoxParams > parseFloatingBoxes(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingBoxParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingBoxes", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks boxBlocks;
   containers[0].getBlocks("FloatingBox", boxBlocks, 0, 10000); // matches FloatingSphere's own cap
   for (const auto& b : boxBlocks)
      result.emplace_back(b);
   return result;
}

// Floating cylinder ("log", type 4): a finite solid mesa_pd Cylinder, fully two-way coupled like
// the floating box. The shape's axis is the BODY z-axis; `axis` gives its initial WORLD direction
// (default <0,1,0> = spanwise log lying across the flow). All box knobs carry over 1:1.
struct FloatingCylinderParams {
   Vector3< real_t > position;
   real_t            radius;
   real_t            length;              // full axis length; halfLength = length / 2
   Vector3< real_t > axis;                // initial world direction of the cylinder axis (normalised on parse)
   real_t            densityRatio;
   uint_t            settlingSteps;       // same semantics as FloatingSphereParams::settlingSteps
   real_t            initialTiltAngleDeg; // extra rotation about the world Y-axis, applied AFTER axis alignment
   real_t            comOffsetZ;          // signed body-frame offset of the COM below the centroid (same as box)
   bool              pinPosition;         // hold TRANSLATION at `position`, rotation free (same as box)
   Vector3< real_t > pinAxes;             // per-axis pin mask (same as box)
   bool              pinRotation;         // per-axis rotation pin, zero angular velocity (same as box)
   Vector3< real_t > pinRotationAxes;     // per-axis rotation pin mask (same as box)
   bool              startInactive;       // see FloatingSphereParams::startInactive for the full comment.

   FloatingCylinderParams(const config::Config::BlockHandle& b)
   {
      position             = b.getParameter< Vector3< real_t > >("position");
      radius               = b.getParameter< real_t >           ("radius");
      length               = b.getParameter< real_t >           ("length");
      axis                 = b.getParameter< Vector3< real_t > >("axis", Vector3< real_t >(real_t(0), real_t(1), real_t(0)));
      densityRatio         = b.getParameter< real_t >           ("densityRatio");
      settlingSteps        = b.getParameter< uint_t >           ("settlingSteps", uint_t(0));
      initialTiltAngleDeg  = b.getParameter< real_t >           ("initialTiltAngleDeg", real_t(0));
      comOffsetZ           = b.getParameter< real_t >           ("comOffsetZ", real_t(0));
      pinPosition          = b.getParameter< bool >             ("pinPosition", false);
      pinAxes              = b.getParameter< Vector3< real_t > >("pinAxes", Vector3< real_t >(real_t(1)));
      pinRotation          = b.getParameter< bool >             ("pinRotation", false);
      pinRotationAxes      = b.getParameter< Vector3< real_t > >("pinRotationAxes", Vector3< real_t >(real_t(1)));
      startInactive        = b.getParameter< bool >             ("startInactive", false);
      if (startInactive && settlingSteps == uint_t(0))
         WALBERLA_LOG_WARNING("FloatingCylinder startInactive=true but settlingSteps=0 -- this floater "
                              "will never be activated (releaseSettledFloatingSpheres no-ops when "
                              "settlingSteps==0).");
      WALBERLA_CHECK_GREATER(axis.length(), real_t(1e-10), "FloatingCylinder axis must be nonzero!");
      axis = axis.getNormalized();
   }
};

inline std::vector< FloatingCylinderParams > parseFloatingCylinders(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingCylinderParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingCylinders", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks cylinderBlocks;
   containers[0].getBlocks("FloatingCylinder", cylinderBlocks, 0, 10000); // matches FloatingSphere's own cap
   for (const auto& b : cylinderBlocks)
      result.emplace_back(b);
   return result;
}

struct ModelParams {
   std::string pdfReconstructionModel;
   std::string pdfRefillingModel;
   std::string excessMassDistributionModel;
   std::string curvatureModel;
   bool        useSimpleMassExchange;
   real_t      cellConversionThreshold;
   real_t      cellConversionForceThreshold;
   // true  → Bogner & Rüde (2013) per-cell gas/liquid/interface state selection for FormerMO cells.
   // false → legacy uniform reconstruction: every uncovered cell becomes an interface cell
   //         (no fill=1 over-promotion band; the all-interface path that never drifted/NaN'd).
   bool        useBognerReconstruction;

   // true  → φ(x_b) fill-weighting of the MovingObstacle (SimpleBB) MEM force (Bogner & Rüde 2013,
   //         Eq. 5.9): interface cells at a floater's waterline contribute proportionally less force,
   //         suppressing the phantom full-liquid buoyant lift.
   // false → pass a null fill accessor so the raw two-phase MEM force is used (reference / ablation).
   bool        useFillWeighting;

   // true  → (DEFAULT) disable the static/dynamic MEM reference height H entirely: SimpleBB is
   //         constructed with a null hydrostaticDensityFct, so no rhoHydStat is ever subtracted from the
   //         boundary-link PDFs. Combined with useFillWeighting (Eq. 5.9's φ-weighted blend, already the
   //         default), this reproduces Bogner & Rüde (2013)'s original three-phase MEM force EXACTLY as
   //         published — no static-H ablation machinery, no separate hydrostatic reference at all.
   //         staticReferenceHeightOverride/useDynamicReferenceHeight are ignored while this is true
   //         (there is nothing left to feed them into). Paired with useMemOnlyBuoyancy=true (also the
   //         default) below, this is now the project default; NOT pairing them (i.e. leaving
   //         useMemOnlyBuoyancy=true while overriding this back to false) reintroduces gravity-only body
   //         force with rhoHydStat still subtracted from the MEM force -- no buoyant support at all,
   //         which is exactly what silently sank every floater in the bridge/1000-sphere scenarios
   //         before this flag's default was flipped (both prms predated this default and never set it
   //         explicitly, so they inherited the old, inconsistent false/true pairing).
   // false → subtract rhoHydStat as before (Antidunes' extended, non-thesis behaviour; the default before
   //         this flag was flipped to true).
   bool        useBognerThreePhaseMEM;

   // true  → (DEFAULT) drop the separate analytic geometric buoyancy term (vSub·g, computeSubmergedVolumes)
   //         from AddExternalForces' vertical body force entirely -- gravity only (rho*vFull*a_z). Buoyant
   //         support then comes ONLY from the MEM force (hydro_force), matching Bogner & Rüde (2013): no
   //         separate buoyancy channel, the hydrostatic pressure already carried by the (uncancelled,
   //         fill-weighted) boundary-link PDFs supplies it. computeSubmergedVolumes()/computeSubmergedVolumesBox()
   //         are SKIPPED entirely in this mode (ParticleKernels::doDEMSubCycles) -- not just excluded from the
   //         force term -- since Bogner's own righting-moment derivation (Eq. 18) uses the raw MEM torque alone,
   //         with no separate analytic buoyancy-lever addition on top; this also matches the torque path,
   //         not just the force. Most physically consistent when paired with useBognerThreePhaseMEM=true (no
   //         double subtraction nor double-counted buoyancy, and now both this project's defaults); pairing
   //         this alone with rhoHydStat subtraction (useBognerThreePhaseMEM=false) removes ALL buoyant
   //         support, not the phantom-lift problem that subtraction was originally added to avoid -- see
   //         useBognerThreePhaseMEM's own comment above for the concrete failure this caused.
   // false → gravity minus vSub*g, plus the analytic vSub-derived righting torque (Antidunes' extended,
   //         non-thesis split; the pre-existing behaviour before this flag defaulted to true).
   bool        useMemOnlyBuoyancy;

   // true  → re-seed each FormerMO cell to the depth-consistent hydrostatic density
   //         ρ_hyd(z) = 1 + 3·a_z·(z − H) instead of the EquilibriumReconstructor's fill-blind
   //         unweighted mean (interface cells report ρ≈1 however empty they are). Kills the ρ≈1
   //         over-density at a moving floater's waterline that the MEM otherwise reads as phantom
   //         upward lift. Robust to negative contact-line fill drift (uses no φ). Off by default.
   bool        useHydrostaticReconstruction;

   // true  → after Pass 2 of FormerMO reconstruction (PdfReconstructionManager) has given every
   //         freshly-uncovered cell in the block a final state, re-examine cells that fell back to
   //         interface purely because of same-pass ambiguity (an unresolved FormerMO neighbour at
   //         the moment they were classified) and resolve them directly against the now-fully-known
   //         neighbourhood. A cheap, purely local second sweep -- no flood fill, no block-wide scan.
   // false → plain Pass-2-only fallback (default; every such cell stays interface, as before).
   bool        usePass3Reresolution;

   // true  → floating-box submerged volume & centre of buoyancy computed ANALYTICALLY by clipping the
   //         oriented box against the fitted free-surface plane (exact 6-tetrahedron OBB∩half-space
   //         clip; no discretisation error, rotation-exact).
   // false → legacy 12³ point-sampling clip of the box against the same plane (default; kept as the
   //         reference and for cross-checking the analytical result).
   bool        useAnalyticBoxBuoyancy;

   // true  → subtract the analytic phantom force F_pred_z = ΔP·A_proj·taper (the double-counted pile-up
   //         head, validated live to ratio ≈ 1) from the VERTICAL MEM force on each floating sphere,
   //         downstream of the fluid (particle side only, no PDF feedback). Removes the spurious waterline
   //         levitation while leaving horizontal drag and the dynamic vertical fluctuations intact
   //         (Option A). Self-disables where there is no pile-up (Δh→0 ⇒ F_pred→0).
   // false → no correction (default; reproduces the pre-fix behaviour).
   bool        usePhantomForceCorrection;

   // true  → the MEM hydrostatic reference height H (used by SimpleBB's rhoHydStat, force-side only —
   //         never touches the fluid/PDF fields) is tracked LIVE instead of held at the flat
   //         absoluteLiquidHeight: each step, the local waterline (localSurfaceZ) of every floater with a
   //         clean fitted ring is averaged and EMA-blended (α=0.2) into a single shared reference height.
   //         Cancels the uncancelled pile-up head at its source (exact, no A_proj/taper geometry needed),
   //         but is a purely STATIC/offset correction — it cannot remove a genuinely velocity-dependent
   //         (dynamic suction) residual, and moves only the WHOLE force curve, not its slope.
   // false → static flat H = absoluteLiquidHeight (default; legacy behaviour).
   bool        useDynamicReferenceHeight;

   // >= 0 → replace the STATIC MEM reference height H (SimpleBB's rhoHydStat, the hydrostatic
   //        reconstruction target, and the PhantomDecomp flatRef) with this absolute z value [LU]
   //        instead of the flat absoluteLiquidHeight. Diagnostic ablation only: shifting H adds a
   //        CONSTANT density offset Δρ = 3·|a_z|·ΔH to the subtracted reference — closed (fully
   //        submerged) bodies should feel nothing, while a floater straddling the surface takes an
   //        unbalanced Δp·A_proj on its open wetted cap. Does NOT touch the fluid initialisation.
   // <  0 → disabled (default; H = absoluteLiquidHeight).
   real_t      staticReferenceHeightOverride;

   // Lateral reach of the free-surface plane-fit used for a sphere floater's own localSurfaceZ AND
   // (when useDynamicReferenceHeight is on) the dynamic-H ambient-waterline target, as a MULTIPLE of
   // the floater's radius: planeLatR = floaterPlaneFitRadiusFactor * r. Default 1.8 reaches past a
   // calmly-floating body's own crater/gas-moat, but NOT past the splash/wave thrown up by a fast
   // impact — the fitted "ambient" surface is still contaminated by the floater's own disturbance
   // during a drop-through-surface transient (see changes.md, dynamic-H self-pollution). Raise this
   // for scenarios with a violent entry.
   real_t      floaterPlaneFitRadiusFactor;

   // Lateral reach of the free-surface plane-fit shell used for a BOX/CYLINDER floater's own
   // submerged-volume clip (computeSubmergedVolumesBox()), as a MULTIPLE of the floater's own largest
   // half-extent: band = boxPlaneFitBandFactor * max(halfExtents) + 2. Same purpose as
   // floaterPlaneFitRadiusFactor above, but boxes/cylinders had this hardcoded to 0.8 (no prm knob)
   // until this field was added. Default 0.8 keeps existing behaviour unchanged; raise it if the
   // fitted plane still looks contaminated by the floater's own local disturbance (crater/meniscus)
   // once the body has stopped moving — see changes.md for the box-stillwater-showcase investigation
   // that motivated this.
   real_t      boxPlaneFitBandFactor;

   // true  → Box-Box and Cylinder-Cylinder DEM contact detection uses the closed-form algorithms in
   //         PaperAnalyticalContacts.h (Gottschalk, Lin & Manocha 1996 OBBTree SAT for Box-Box; Kodam
   //         et al. 2010 + Guo et al. 2012's corrections for Cylinder-Cylinder) instead of the generic
   //         GJK/EPA fallback. Box-Box is trusted directly (fully verified against GJK/EPA, no known
   //         discrepancy). Cylinder-Cylinder is a HYBRID: cylinderCylinderResultIsUncertain() gates
   //         every analytical result and falls back to GJK/EPA automatically for the oblique-tilt cases
   //         the papers' own scenario taxonomy does not cover cleanly — see PaperAnalyticalContacts.h
   //         and PaperAnalyticalContactTest.cpp (33 cross-checked cases) for the full derivation.
   // false → GJK/EPA for both pairs, exactly the pre-existing behaviour (default; flip this back to
   //         false to return to the known-working path if the analytical route misbehaves).
   bool        useAnalyticalContactDetection;

   // true  → (only meaningful together with useAnalyticalContactDetection=true) removes the
   //         Cylinder-Cylinder GJK/EPA safety net entirely: cylinderCylinderResultIsUncertain()'s
   //         verdict is IGNORED, and the closed-form Kodam et al. (2010) + Guo et al. (2012) dispatch
   //         is trusted as the sole answer -- exactly "pure Kodam+Guo", nothing added, matching how
   //         those two papers' own algorithm stands on its own (neither paper uses GJK/EPA or any
   //         other iterative solver anywhere; this project's hybrid fallback is entirely this
   //         codebase's own addition on top of both papers, not part of either's methodology). A
   //         tier-1 case (analytical found no contact at all) is accepted as "no contact", no
   //         fallback attempted; a tier-2 case (analytical found a contact but flagged its own
   //         magnitude untrustworthy) is still applied as-is. The independent raw-center-distance
   //         guard in SolidCylinderContactDetection (rejects pairs geometrically too far apart to
   //         possibly touch, added after a false-positive band-band-skewed crash -- see ParticleKernels.h)
   //         stays active regardless of this flag: it is not part of either paper's scenario dispatch,
   //         just a cheap pre-filter, so disabling it isn't needed to run "pure" Kodam+Guo honestly.
   // false → (DEFAULT) the existing hybrid behaviour: fall back to GJK/EPA whenever
   //         cylinderCylinderResultIsUncertain() flags the analytical result.
   bool        usePureAnalyticalCylinder;

   // true  → every rank writes its own CSV file under floaterHistoryBaseFolder/ (rank_<N>.csv,
   //         truncated at simulation start), one row per contact/lubrication pair it actually computed
   //         and accepted (i.e. that passed mesa_pd::mpi::ContactFilter), for EVERY DEM sub-cycle of
   //         EVERY main LBM step (not gated behind the usual WATCH_FREQ/evaluationFrequency cadence
   //         used by the coarser [CylForces]/[BoxForces] logs -- this is meant to catch the exact
   //         sub-step transient those miss). Each row records: which rank computed it, both particles'
   //         uid/isGhost, and -- critically -- ac.getOwner(idx) for both particles, i.e. mesa_pd's own
   //         synced-to-every-copy "canonical owner" field (see mesa_pd::mpi::ContactFilter's `*`/`#`
   //         rules, which is the actual arbitration rule this project's still-unimplemented "Option 2"
   //         tiebreak would replace). Comparing the logging rank against getOwner() across all ranks'
   //         CSVs for the same (uidA,uidB,mainStep,subCycle) is exactly how to catch a double-accept:
   //         two different ranks each believing they should compute this contact in the same cycle.
   //         Also logs the contact point/normal/penetrationDepth and the FORCE DELTA this specific
   //         contact/lubrication pair contributed (captured via force-ref before/after readback, since
   //         LinearSpringDashpot/LubricationCorrectionKernel add atomically onto whatever force other
   //         contacts already accumulated this sub-cycle). Every row is flushed immediately (not
   //         buffered) so data survives an abrupt WALBERLA_ABORT/MPI_ABORT mid-run, mirroring the
   //         WALBERLA_LOG_WARNING-not-_INFO lesson from ParticleKernels' own logCylComputation. OFF by
   //         default: this is a heavy, opt-in diagnostic (one row per pair per sub-cycle -- easily
   //         millions of rows for a long dense run) meant for actively debugging the ContactFilter race,
   //         not routine use.
   bool        logFloaterHistory;

   // Base folder for logFloaterHistory's per-rank CSVs, sibling to vtk-out/eval/snapshots (created by
   // Antidunes.cpp alongside those, root-only, only when logFloaterHistory is true).
   std::string floaterHistoryBaseFolder;

   // 0 (DEFAULT) → logFloaterHistory (both the per-contact log() and the per-floater logState()
   // world-state snapshot) covers every floater, as before. Nonzero → restrict BOTH logs to only the
   // floaters whose settlingSteps equals this value (i.e. one specific staggered-release batch), by
   // UID -- everyone else is skipped entirely, not just filtered after writing. Exists because
   // logFloaterHistory's per-subcycle, local+ghost-copy world-state snapshot is what blew a 126GB
   // quota when it covered all 1000 floaters in the collapsed-spawn scenario; scoping it to one
   // ~100-floater release batch keeps the same per-substep position/velocity resolution at a
   // fraction of the size. Matched by FloatingBox/FloatingSphere/FloatingCylinder settlingSteps in
   // ParticleKernels.h's constructor (built once into a UID set, wired into FloaterHistoryLogger via
   // setUidFilter()) -- has no effect on which floaters actually get released/settle, only on logging.
   uint_t      logFloaterHistoryReleaseStepFilter;

   // 0 (DEFAULT) -> logFloaterHistory writes from step 0, as before. Nonzero -> both log() (per-
   // contact) and logState() (per-floater world-state) skip every row until mainStep reaches this
   // value, same zero-cost-when-disabled spirit as logFloaterHistoryReleaseStepFilter above but
   // scoping by TIME instead of by floater identity. For a targeted rerun of a crash whose approximate
   // step is already known (e.g. from a prior run's own [LOW-RHO]/[MIGRATION-NAN] step number), this
   // lets logFloaterHistory run at full (unfiltered-by-UID) resolution for just the steps that matter
   // instead of paying its per-subcycle cost for the entire run up to that point.
   uint_t      logFloaterHistoryStartStep;

   // Diagnostic for the "how did water get up there" question: a job 4132582 [LOW-RHO] crash was
   // traced (VTK velocity field, decoded directly) to a region that was completely dry gas 592 steps
   // before the corruption appeared, with no floater anywhere nearby -- meaning some cell above the
   // normal waterline was converted from gas to liquid/interface with NO log trail at all (neither
   // waLBerla's own CellConversionSweep, src/lbm/free_surface/dynamics/CellConversionSweep.h, nor this
   // project's wrapping of it in BuildTimeloop.h ever logs anything; only the separate, DEM-specific
   // FormerMO reconstruction path has its own [FORMERMO-CONVERT] counters). This scans every LOCAL
   // cell above highAltitudeWaterThresholdZ, every step, and warns the instant one is liquid or
   // interface -- catching the conversion itself, not just its aftermath. OFF by default (0 threshold
   // is meaningless/inactive unless explicitly enabled): opt-in for a targeted rerun once a crash's
   // approximate altitude is already known (e.g. from a prior run's own [LOW-RHO] cell coordinates),
   // not something to leave on for routine production runs -- a healthy run should trigger it zero
   // times, but a run that's actively corrupting the free surface at height could log heavily.
   bool        highAltitudeWaterScanEnabled;
   real_t      highAltitudeWaterThresholdZ;

   // true  → each ArchTunnelParams entry (SimulationParams.h) is instantiated as a REAL mesa_pd body
   //         (a CylindricalBoundary in containment mode -- see that shape's own doc comment and
   //         SolidCylinderContactDetection's new containment overloads in ParticleKernels.h) at
   //         setupParticles time, going through the exact same forEachParticlePairHalf +
   //         LinearSpringDashpot two-particle contact/force pipeline piers already use. The old
   //         one-sided ComputeTunnelContacts kernel is NOT constructed/called when this is true
   //         (both active at once would double-apply the containment force for the same tunnel).
   // false → (DEFAULT) exact pre-existing behaviour: arch tunnels are the one-sided
   //         ComputeTunnelContacts kernel, no real mesa_pd body is created for them. Flip this back
   //         to false to return to the known-working path if the real-object route misbehaves.
   bool        useArchTunnelAsMesaPDObject;

   // true  -> mesa_pd::mpi::ContactFilter::useOwnerRankTiebreak is set on every ComputeContacts/
   //         ComputeLubrication instance: the local-remote/remote-remote ownership-arbitration rule
   //         switches from the pre-existing per-rank geometric test (domain.isContainedInProcessSubdomain
   //         on each rank's own independently-computed contactPoint -- provably able to disagree between
   //         ranks, see ContactFilter.h's own doc comment) to a deterministic owner-rank-number tiebreak
   //         (accept iff myRank == min(ac.getOwner(idx1), ac.getOwner(idx2))). Root-caused specifically
   //         for job 4166652's ContactHistoryNotification duplicate-insert crash (Cylinder x Box pair
   //         straddling a rank boundary).
   // false -> (DEFAULT) exact pre-existing ContactFilter behaviour, unchanged.
   bool        useContactFilterOwnerTiebreak;

   ModelParams(const config::Config::BlockHandle& b)
      : pdfReconstructionModel      { b.getParameter< std::string >("pdfReconstructionModel") }
      , pdfRefillingModel           { b.getParameter< std::string >("pdfRefillingModel") }
      , excessMassDistributionModel { b.getParameter< std::string >("excessMassDistributionModel") }
      , curvatureModel              { b.getParameter< std::string >("curvatureModel") }
      , useSimpleMassExchange       { b.getParameter< bool >       ("useSimpleMassExchange") }
      , cellConversionThreshold     { b.getParameter< real_t >     ("cellConversionThreshold") }
      , cellConversionForceThreshold{ b.getParameter< real_t >     ("cellConversionForceThreshold") }
      , useBognerReconstruction     { b.getParameter< bool >       ("useBognerReconstruction", true) }
      , useFillWeighting            { b.getParameter< bool >       ("useFillWeighting", true) }
      , useBognerThreePhaseMEM      { b.getParameter< bool >       ("useBognerThreePhaseMEM", true) }
      , useMemOnlyBuoyancy          { b.getParameter< bool >       ("useMemOnlyBuoyancy", true) }
      , useHydrostaticReconstruction { b.getParameter< bool >      ("useHydrostaticReconstruction", false) }
      , usePass3Reresolution        { b.getParameter< bool >       ("usePass3Reresolution", false) }
      , useAnalyticBoxBuoyancy      { b.getParameter< bool >       ("useAnalyticBoxBuoyancy", false) }
      , usePhantomForceCorrection   { b.getParameter< bool >       ("usePhantomForceCorrection", false) }
      , useDynamicReferenceHeight   { b.getParameter< bool >       ("useDynamicReferenceHeight", false) }
      , staticReferenceHeightOverride{ b.getParameter< real_t >    ("staticReferenceHeightOverride", real_t(-1)) }
      , floaterPlaneFitRadiusFactor  { b.getParameter< real_t >    ("floaterPlaneFitRadiusFactor", real_t(1.8)) }
      , boxPlaneFitBandFactor        { b.getParameter< real_t >    ("boxPlaneFitBandFactor", real_t(0.8)) }
      , useAnalyticalContactDetection{ b.getParameter< bool >      ("useAnalyticalContactDetection", false) }
      , usePureAnalyticalCylinder   { b.getParameter< bool >       ("usePureAnalyticalCylinder", false) }
      , logFloaterHistory           { b.getParameter< bool >       ("logFloaterHistory", false) }
      , floaterHistoryBaseFolder    { b.getParameter< std::string >("floaterHistoryBaseFolder", "floater-history") }
      , logFloaterHistoryReleaseStepFilter{ b.getParameter< uint_t >("logFloaterHistoryReleaseStepFilter", uint_t(0)) }
      , logFloaterHistoryStartStep   { b.getParameter< uint_t >     ("logFloaterHistoryStartStep", uint_t(0)) }
      , highAltitudeWaterScanEnabled { b.getParameter< bool >       ("highAltitudeWaterScanEnabled", false) }
      , highAltitudeWaterThresholdZ  { b.getParameter< real_t >     ("highAltitudeWaterThresholdZ", real_t(0)) }
      , useArchTunnelAsMesaPDObject  { b.getParameter< bool >       ("useArchTunnelAsMesaPDObject", false) }
      , useContactFilterOwnerTiebreak{ b.getParameter< bool >       ("useContactFilterOwnerTiebreak", false) }
   {
      if (usePureAnalyticalCylinder && !useAnalyticalContactDetection)
         WALBERLA_LOG_WARNING("usePureAnalyticalCylinder=true but useAnalyticalContactDetection=false -- "
            "this flag only has an effect inside the Cylinder x Cylinder analytical dispatch, which is "
            "off, so it does nothing here.");
   }
};

struct EvaluationParams {
   uint_t      performanceLogFrequency;
   uint_t      evaluationFrequency;
   std::string baseFolderName;

   EvaluationParams(const config::Config::BlockHandle& b)
      : performanceLogFrequency { b.getParameter< uint_t >     ("performanceLogFrequency") }
      , evaluationFrequency     { b.getParameter< uint_t >     ("evaluationFrequency") }
      , baseFolderName          { b.getParameter< std::string >("baseFolderName") }
   {}
};

struct VTKParams {
   uint_t            writeFrequency;
   uint_t            ghostLayers;
   std::string       baseFolder;
   real_t            samplingResolution;

   // Windowed higher-frequency fluid_field VTK output around an expected crash step -- same
   // rationale/pattern as PhysicsParams::denseNanScanEnabled (see its own doc comment), applied to
   // VTK dumps instead of the NaN scanner: a NaN blowup found by the periodic (every writeFrequency
   // steps) checkPdfNaN/VTK cadence is normally already spread across much of a block by the time
   // it's caught -- not enough to see the actual onset. Parsed from an OPTIONAL second named VTK
   // output block "fluid_field_dense" (see VtkWriter.h's addVTKOutput, which is what actually gates
   // its execution to [denseFromStep, denseToStep) -- outside that window it is never invoked at
   // all, so it costs nothing for the rest of the run). If that block is absent, denseEnabled is
   // false and every existing prm behaves exactly as before.
   bool              denseEnabled;
   uint_t            denseFromStep;
   uint_t            denseToStep;

   VTKParams(const config::Config::BlockHandle& b)
      : writeFrequency    { b.getOneBlock("fluid_field").getParameter< uint_t >     ("writeFrequency") }
      , ghostLayers       { b.getOneBlock("fluid_field").getParameter< uint_t >     ("ghostLayers") }
      , baseFolder        { b.getOneBlock("fluid_field").getParameter< std::string >("baseFolder") }
      , samplingResolution{ b.getOneBlock("fluid_field").getParameter< real_t >     ("samplingResolution") }
      , denseEnabled      { static_cast< bool >(b.getBlock("fluid_field_dense")) }
      , denseFromStep     { denseEnabled ? b.getOneBlock("fluid_field_dense").getParameter< uint_t >("denseFromStep") : uint_t(0) }
      , denseToStep       { denseEnabled ? b.getOneBlock("fluid_field_dense").getParameter< uint_t >("denseToStep") : uint_t(0) }
   {}
};

struct SimulationParams;

// Parameters that are computed once and for all (derived from config params)
//-----------------------------------------------------------------------
struct ComputedParams {
   real_t avgParticleDiameter;
   real_t particleScalingFactor;
   real_t particleKappa;
   real_t particleCollisionTime;
   real_t particleTimestepSize;
   real_t liquidHeight;
   real_t floorHeight;
   real_t absoluteLiquidHeight;
   real_t particleFixingHeight;
   real_t viscosity;
   real_t relaxationRate;
   real_t gravity;
   real_t surfaceTension;
   real_t dx_SI;
   real_t dt_SI;
   real_t sinusAmplitude;
   real_t sinusOffset;
   real_t sinusWavelength;

   ComputedParams(
      const ParticleParams& particle,
      const DomainParams& domain,
      const PIDParams& pid,
      const PhysicsParams& physics)
   {
      // avgParticleDiameterOverride/dxSIOverride (both > 0): skip inFileName's header entirely --
      // see ParticleParams's doc comment. Otherwise fall back to the bed-file-derived scaling, exactly
      // as before (also covers useParticleBed=true, which always needs the real bed-derived values).
      if (particle.avgParticleDiameterOverride > real_c(0) && particle.dxSIOverride > real_c(0))
      {
         avgParticleDiameter   = particle.avgParticleDiameterOverride;
         particleScalingFactor = real_c(1) / particle.dxSIOverride;
         WALBERLA_LOG_INFO_ON_ROOT("avgParticleDiameter/particleScalingFactor overridden directly "
            "(avgParticleDiameterOverride=" << particle.avgParticleDiameterOverride
            << ", dxSIOverride=" << particle.dxSIOverride << ") -- " << particle.inFileName
            << "'s header is NOT read for scaling purposes.")
      }
      else
      {
         getAvgDiameterScalingFactor(
            particle.inFileName, domain.domainSize,
            particle.bedCopiesInX, particle.bedCopiesInY,
            avgParticleDiameter, particleScalingFactor);
      }

      const real_t particlePoissonsRatio = real_c(0.22);
      particleKappa = real_c(2) * (real_c(1) - particlePoissonsRatio) / (real_c(2) - particlePoissonsRatio);
      const real_t particleCollisionTimeNonDim = real_c(4);
      particleCollisionTime = particleCollisionTimeNonDim * avgParticleDiameter;
      particleTimestepSize = real_c(1) / real_c(particle.numSubCycles);
      
      liquidHeight         = avgParticleDiameter * domain.liquidHeightFactor;
      floorHeight          = avgParticleDiameter * domain.floorHeightFactor;
      WALBERLA_CHECK_FLOAT_UNEQUAL(liquidHeight, 0.0)
      WALBERLA_CHECK_FLOAT_UNEQUAL(floorHeight, 0.0)
      absoluteLiquidHeight = liquidHeight + floorHeight;
      particleFixingHeight = particle.fixingHeightFactor * avgParticleDiameter;

      const real_t U  = pid.targetMeanVelocityMagnitude;
      if (physics.useHardcodedUnits)
      {
         // Fixed reference units (the values the dynamic formulas below produce at the default config:
         // U=0.02, Re=3100, Fr=1.31, We=15.6188). Use these to run with targetMeanVelocityMagnitude=0
         // (still water) without zeroing viscosity/gravity/surface tension.
         viscosity      = real_c(0.000214026);
         relaxationRate = real_c(1.99743);
         gravity        = real_c(7.02618e-06);
         surfaceTension = real_c(0.000849592);
      }
      else
      {
         viscosity      = U * liquidHeight / physics.Re;                         // 0.000214026
         relaxationRate = real_c(1) / (real_c(3) * viscosity + real_c(0.5));     // 1.99743
         gravity        = (U / physics.Fr) * (U / physics.Fr) / liquidHeight;    // 7.02618e-06
         surfaceTension = real_c(1) * U * U * liquidHeight / physics.We;         // 0.000849592
      }

      // Damping override: set the relaxation rate directly and derive the matching viscosity
      // (ν = (1/ω − 0.5)/3). Lower ω → more viscous → more damping. 0 = keep the value above.
      if (physics.relaxationRateOverride > real_c(0))
      {
         relaxationRate = physics.relaxationRateOverride;
         viscosity      = (real_c(1) / relaxationRate - real_c(0.5)) / real_c(3);
      }

      // Gravity override: see PhysicsParams::gravityOverride doc comment. 0 = keep the value above.
      if (physics.gravityOverride > real_c(0)) { gravity = physics.gravityOverride; }

      const real_t viscosity_SI = real_c(1.0016e-6);  // kinemtic viscosity of water at 20°C at 1 bar
      dx_SI = real_c(1) / particleScalingFactor;
      dt_SI = viscosity / viscosity_SI * dx_SI * dx_SI;

      sinusAmplitude  = real_c(0.5) * domain.initialAmplitude;
      sinusOffset     = floorHeight;
      sinusWavelength = real_c(domain.domainSize[0]) / real_c(domain.wavePeriods);
   }

   void logOnRoot() const
   {
      WALBERLA_LOG_INFO_ON_ROOT("Computed parameters --------------------")
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(avgParticleDiameter);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleScalingFactor);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleKappa);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleCollisionTime);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleTimestepSize);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(liquidHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(floorHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(absoluteLiquidHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleFixingHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(gravity);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(surfaceTension);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dx_SI);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dt_SI);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(sinusAmplitude);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(sinusOffset);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(sinusWavelength);
      WALBERLA_LOG_INFO_ON_ROOT("------------------------------------------")
   }

   void writeToFile(const std::string& baseFolderName) const
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file(baseFolderName + "/computedParams.txt");
         file << "avgParticleDiameter:   " << avgParticleDiameter    << "\n";
         file << "particleScalingFactor: " << particleScalingFactor  << "\n";
         file << "particleKappa:         " << particleKappa          << "\n";
         file << "particleCollisionTime: " << particleCollisionTime  << "\n";
         file << "particleTimestepSize:  " << particleTimestepSize   << "\n";
         file << "liquidHeight:          " << liquidHeight           << "\n";
         file << "floorHeight:           " << floorHeight            << "\n";
         file << "absoluteLiquidHeight:  " << absoluteLiquidHeight   << "\n";
         file << "particleFixingHeight:  " << particleFixingHeight   << "\n";
         file << "viscosity:             " << viscosity              << "\n";
         file << "relaxationRate:        " << relaxationRate         << "\n";
         file << "gravity:               " << gravity                << "\n";
         file << "surfaceTension:        " << surfaceTension         << "\n";
         file << "dx_SI:                 " << dx_SI                  << "\n";
         file << "dt_SI:                 " << dt_SI                  << "\n";
         file << "sinusAmplitude:        " << sinusAmplitude         << "\n";
         file << "sinusOffset:           " << sinusOffset            << "\n";
         file << "sinusWavelength:       " << sinusWavelength        << "\n";
         file.close();
      }
   }
};

// Container for all simulation parameters 
//-----------------------------------------------------------------------
struct SimulationParams {
   // config ----------
   shared_ptr< config::Config > cfgPtr;
   
   // config params ----------
   BlockForestParams blockForest;
   DomainParams      domain;
   PIDParams         pid;
   PhysicsParams     physics;
   ParticleParams                       particle;
   std::vector< PillarParams >          pillars;
   std::vector< CylinderPillarParams >  cylinderPillars;
   std::vector< ArchTunnelParams >      archTunnels;
   std::vector< FloatingSphereParams >  floatingSpheres;
   std::vector< FloatingBoxParams >     floatingBoxes;
   std::vector< FloatingCylinderParams > floatingCylinders;
   ModelParams       model;
   EvaluationParams  evaluation;
   VTKParams         vtk;

   // computed params ----------
   ComputedParams    computed;

   // Effective STATIC MEM reference height H: the flat still-water level, unless
   // ModelParameters.staticReferenceHeightOverride (>= 0) replaces it (diagnostic ablation).
   real_t staticReferenceHeight() const
   {
      return model.staticReferenceHeightOverride >= real_t(0) ? model.staticReferenceHeightOverride
                                                              : computed.absoluteLiquidHeight;
   }

   SimulationParams(const shared_ptr< config::Config >& cfg)
      : cfgPtr       { cfg }
      , blockForest  { cfg->getOneBlock("BlockForestParameters") }
      , domain       { cfg->getOneBlock("DomainParameters") }
      , pid          { cfg->getOneBlock("PIDParameters") }
      , physics      { cfg->getOneBlock("PhysicsParameters") }
      , particle     { cfg->getOneBlock("ParticleParameters") }
      , pillars          { parsePillars(cfg) }
      , cylinderPillars  { parseCylinderPillars(cfg) }
      , archTunnels      { parseArchTunnels(cfg) }
      , floatingSpheres  { parseFloatingSpheres(cfg) }
      , floatingBoxes    { parseFloatingBoxes(cfg) }
      , floatingCylinders{ parseFloatingCylinders(cfg) }
      , model        { cfg->getOneBlock("ModelParameters") }
      , evaluation   { cfg->getOneBlock("EvaluationParameters") }
      , vtk          { cfg->getOneBlock("VTK") }
      , computed     { particle, domain, pid, physics }
   {}
};


} // namespace antidunes
} // namespace walberla