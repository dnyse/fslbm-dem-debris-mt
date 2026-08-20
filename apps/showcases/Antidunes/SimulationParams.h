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
#include "core/math/Vector3.h"

#include "Utility.h"


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

   PhysicsParams(const config::Config::BlockHandle& b)
      : enableWetting          { b.getParameter< bool >  ("enableWetting") }
      , timesteps              { b.getParameter< uint_t >("timesteps") }
      , Re                     { b.getParameter< real_t >("Re") }
      , Fr                     { b.getParameter< real_t >("Fr") }
      , We                     { b.getParameter< real_t >("We") }
      , useHardcodedUnits      { b.getParameter< bool >  ("useHardcodedUnits", false) }
      , relaxationRateOverride { b.getParameter< real_t >("relaxationRateOverride", real_c(0)) }
      , gravityOverride        { b.getParameter< real_t >("gravityOverride", real_c(0)) }
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
   }
};

inline std::vector< FloatingBoxParams > parseFloatingBoxes(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingBoxParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingBoxes", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks boxBlocks;
   containers[0].getBlocks("FloatingBox", boxBlocks, 0, 100);
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
   containers[0].getBlocks("FloatingCylinder", cylinderBlocks, 0, 100);
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

   // true  → disable the static/dynamic MEM reference height H entirely: SimpleBB is constructed with
   //         a null hydrostaticDensityFct, so no rhoHydStat is ever subtracted from the boundary-link
   //         PDFs. Combined with useFillWeighting (Eq. 5.9's φ-weighted blend, already the default),
   //         this reproduces Bogner & Rüde (2013)'s original three-phase MEM force EXACTLY as
   //         published — no static-H ablation machinery, no separate hydrostatic reference at all.
   //         staticReferenceHeightOverride/useDynamicReferenceHeight are ignored while this is true
   //         (there is nothing left to feed them into).
   // false → subtract rhoHydStat as before (default; Antidunes' extended, non-thesis behaviour).
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
   //         double subtraction nor double-counted buoyancy); pairing this alone with the default rhoHydStat
   //         subtraction re-introduces the original phantom-lift problem that subtraction was added to avoid.
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
      , useBognerThreePhaseMEM      { b.getParameter< bool >       ("useBognerThreePhaseMEM", false) }
      , useMemOnlyBuoyancy          { b.getParameter< bool >       ("useMemOnlyBuoyancy", true) }
      , useHydrostaticReconstruction { b.getParameter< bool >      ("useHydrostaticReconstruction", false) }
      , usePass3Reresolution        { b.getParameter< bool >       ("usePass3Reresolution", false) }
      , useAnalyticBoxBuoyancy      { b.getParameter< bool >       ("useAnalyticBoxBuoyancy", false) }
      , usePhantomForceCorrection   { b.getParameter< bool >       ("usePhantomForceCorrection", false) }
      , useDynamicReferenceHeight   { b.getParameter< bool >       ("useDynamicReferenceHeight", false) }
      , staticReferenceHeightOverride{ b.getParameter< real_t >    ("staticReferenceHeightOverride", real_t(-1)) }
      , floaterPlaneFitRadiusFactor  { b.getParameter< real_t >    ("floaterPlaneFitRadiusFactor", real_t(1.8)) }
      , boxPlaneFitBandFactor        { b.getParameter< real_t >    ("boxPlaneFitBandFactor", real_t(0.8)) }
   {}
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

   VTKParams(const config::Config::BlockHandle& b)
      : writeFrequency    { b.getOneBlock("fluid_field").getParameter< uint_t >     ("writeFrequency") }
      , ghostLayers       { b.getOneBlock("fluid_field").getParameter< uint_t >     ("ghostLayers") }
      , baseFolder        { b.getOneBlock("fluid_field").getParameter< std::string >("baseFolder") }
      , samplingResolution{ b.getOneBlock("fluid_field").getParameter< real_t >     ("samplingResolution") }
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