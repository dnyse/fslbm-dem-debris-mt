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
#include <set>
#include <sstream>
#include <string>


namespace walberla {
namespace antidunes {

struct BlockForestParams {
   Vector3< uint_t > cellsPerBlock;
   Vector3< bool >   periodicity;
   bool              loadSnapshot;
   bool              storeSnapshot;
   uint_t            snapshotFrequency;
   std::string       snapshotBaseFolder;

   // Diagnostic (needs storeSnapshot): from step keepAllSnapshotsFromStep on, every step that is not a regular
   // checkpoint step is checkpointed into its own step_<N> subfolder. Default false / 0.
   bool              keepAllSnapshots;
   uint_t            keepAllSnapshotsFromStep;

   // true -> the regular checkpoint (every snapshotFrequency steps) alternates between snapshotBaseFolder/gen_0 and
   // gen_1, so the two most recent checkpoints survive. Loading is unchanged: point a resume prm's snapshotBaseFolder
   // at the newer gen_<N> (compare beginTimeStep in simState.file). Default false (one rolling checkpoint).
   bool              keepLastTwoSnapshots;

   BlockForestParams(const config::Config::BlockHandle& b)
      : loadSnapshot       { b.getParameter< bool >               ("loadSnapshot") }
      , storeSnapshot      { b.getParameter< bool >               ("storeSnapshot") }
      , snapshotFrequency  { b.getParameter< uint_t >             ("snapshotFrequency") }
      , snapshotBaseFolder { b.getParameter< std::string >        ("snapshotBaseFolder") }
      , keepAllSnapshots         { b.getParameter< bool >   ("keepAllSnapshots", false) }
      , keepAllSnapshotsFromStep { b.getParameter< uint_t > ("keepAllSnapshotsFromStep", uint_c(0)) }
      , keepLastTwoSnapshots     { b.getParameter< bool >   ("keepLastTwoSnapshots", false) }
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

   // false -> the PID update is disabled and acceleration[0] stays at its initial value from setupSimState (fixed
   // hydraulic gradient instead of a controlled flow; ForcingAdjuster, Forcing.h). Default true.
   bool   enabled;

   // > 0 -> replaces the initial accelerationX [LU] (also the PID's initial actuating variable). The default formula
   // 3*U*nu/liquidHeight^2 assumes an unobstructed channel and under-drives flows past a pillar; intended for PID-off
   // runs, e.g. with the converged actuatingVariable of a matching PID-on run. 0 = off (default).
   real_t accelerationXOverride;

   // switchOffStep > 0: PID control until that absolute timestep, then acceleration_x is frozen at the mean of the
   // PID's commands over the last switchOffAveragingWindow steps (the mean, because the command oscillates around its
   // converged value) and never adjusted again; window 0 freezes the last command. Lets the flow spin up under PID
   // control before debris arrives. switchOffStep 0 = off (default); window default 5000. Ignored if enabled=false.
   uint_t switchOffStep;
   uint_t switchOffAveragingWindow;

   PIDParams(const config::Config::BlockHandle& b)
      : targetMeanVelocityMagnitude { b.getParameter< real_t >("targetMeanVelocityMagnitude") }
      , proportionalGain            { b.getParameter< real_t >("proportionalGain") }
      , derivativeGain              { b.getParameter< real_t >("derivativeGain") }
      , integralGain                { b.getParameter< real_t >("integralGain") }
      , maxRamp                     { b.getParameter< real_t >("maxRamp") }
      , minActuatingVariable        { b.getParameter< real_t >("minActuatingVariable") }
      , maxActuatingVariable        { b.getParameter< real_t >("maxActuatingVariable") }
      , enabled                     { b.getParameter< bool >  ("enabled", true) }
      , accelerationXOverride       { b.getParameter< real_t >("accelerationXOverride", real_c(0)) }
      , switchOffStep               { b.getParameter< uint_t >("switchOffStep", uint_c(0)) }
      , switchOffAveragingWindow    { b.getParameter< uint_t >("switchOffAveragingWindow", uint_c(5000)) }
   {}
};

struct PhysicsParams {
   bool   enableWetting;
   uint_t timesteps;
   real_t Re;
   real_t Fr;
   real_t We;
   // true -> the U-derived LBM units (viscosity, relaxationRate, gravity, surfaceTension) use fixed reference values
   // (ComputedParams), so the mean velocity can be 0 (still water) with valid physical constants. Default false.
   bool   useHardcodedUnits;
   // > 0 -> sets the LBM relaxation rate directly and derives the matching viscosity, e.g. to damp a still-water
   // settle test. 0 = off (default).
   real_t relaxationRateOverride;
   // > 0 -> sets the lattice gravity directly, overriding both values above, e.g. g = 7.5e-4 of Bogner & Rude (2013)
   // for their floating-cube validation (the default Fr-derived value is ~7.0e-6). 0 = off (default).
   real_t gravityOverride;
   // Surface tension [LU]: >= 0 replaces the derived/hardcoded value (0 = no surface tension); < 0 (default -1)
   // keeps it.
   real_t surfaceTensionOverride;

   // Dense NaN scan: checkPdfNaN normally runs every WATCH_FREQ steps (Antidunes.cpp), by which time a NaN has spread
   // (~1 cell/step). denseNanScanEnabled=true runs it every step in [denseNanScanFromStep, denseNanScanToStep)
   // (defaults 0 and timesteps), typically in a rerun from a checkpoint shortly before a known crash. Default false.
   // The precursor thresholds below are evaluated on the same cadence.
   bool   denseNanScanEnabled;
   uint_t denseNanScanFromStep;
   uint_t denseNanScanToStep;

   // > 0 -> checkPdfNaN (SetupStages.h) also warns [LOW-RHO-PRECURSOR] (first 20 per rank) for liquid/interface cells
   // whose density is finite but below this value (e.g. 0.1; nominal rho ~ 1). 0 = off (default).
   real_t lowRhoPrecursorThreshold;

   // > 0 -> also warns [HIGH-VELOCITY-PRECURSOR] (first 20 per rank) for liquid/interface cells whose finite velocity
   // magnitude exceeds this value (e.g. 10x targetMeanVelocityMagnitude). 0 = off (default).
   real_t highVelocityPrecursorThreshold;

   // > 0 -> also warns [POP-PRECURSOR] (first 40 per rank) for liquid/interface cells with some |f_i| above this value
   // (physical populations are <~ 0.3 in LU), usually several steps before |u| grows. The line reports rho, fill and
   // the obstacle geometry: pinched opposite D3Q27 direction pairs (of 13) and obstacle neighbours (of 26).
   // 0 = off (default).
   real_t populationPrecursorThreshold;

   // Density window, each side enabled by a value > 0 (default 0 = off). checkPdfNaN reports liquid/interface cells
   // (liquid only if densityPrecursorLiquidOnly) with rho < densityPrecursorLow or rho > densityPrecursorHigh (typical
   // 0.1 / 1.9): per rank one [DENSITY-PRECURSOR] line per step (first 60) with counts and details of the first,
   // lowest- and highest-rho cell, plus a root [DENSITY-WINDOW] summary (first 400). Same cadence as the dense NaN
   // scan. Drained interface cells at the free surface legitimately have rho -> 0; densityPrecursorLiquidOnly skips
   // them.
   real_t densityPrecursorLow;
   real_t densityPrecursorHigh;
   bool   densityPrecursorLiquidOnly;

   // Fluid sanity sweep (default off): extra sweep at the end of every main step (after the evaluators, Antidunes.cpp),
   // every fluidSanityEveryNSteps steps in [fluidSanityFromStep, fluidSanityToStep), checking every liquid/interface
   // cell (liquid only if fluidSanityLiquidOnly) for: non-finite f_i/rho/u; rho outside [DensityMin, DensityMax]
   // (defaults 0.1 / 1.9); |u| > fluidSanityMaxVelocity; max|f_i| > fluidSanityMaxPopulation (default 2); fill outside
   // [-tol, 1+tol] (interface) or [1-tol, 1+tol] (liquid), tol = fluidSanityFillTolerance (default 0.05).
   // fluidSanityMaxVelocity default 0.17 = lattice Mach number 0.3 (Kruger et al. 2017, Sect. 7.2, p. 279): an accuracy
   // warning, not a collapse threshold. Output: per rank [FLUID-SANITY] (first 60 lines), root [FLUID-SANITY-SUMMARY]
   // (first 400). fluidSanityAbortAbove > 0 aborts once one step has more violating cells (0 = report only).
   // A threshold of 0 (fill tolerance < 0) disables that category.
   bool   fluidSanityEnabled;
   uint_t fluidSanityFromStep;
   uint_t fluidSanityToStep;
   uint_t fluidSanityEveryNSteps;
   real_t fluidSanityDensityMin;
   real_t fluidSanityDensityMax;
   real_t fluidSanityMaxVelocity;
   real_t fluidSanityMaxPopulation;
   real_t fluidSanityFillTolerance;
   bool   fluidSanityLiquidOnly;
   uint_t fluidSanityAbortAbove;
   // Additional opt-in checks (0 / false = off): |excess mass|; largest |u| / |rho| jump to a face neighbour in the
   // same block (spike detector); liquid cell touching gas (broken interface layer). fluidSanityStatsEveryNSteps > 0
   // prints [FLUID-STATS] (global extrema) every N steps regardless of violations.
   real_t fluidSanityMaxExcessMass;
   real_t fluidSanityMaxNeighbourVelocityJump;
   real_t fluidSanityMaxNeighbourDensityJump;
   bool   fluidSanityCheckTopology;
   uint_t fluidSanityStatsEveryNSteps;

   PhysicsParams(const config::Config::BlockHandle& b)
      : enableWetting          { b.getParameter< bool >  ("enableWetting") }
      , timesteps              { b.getParameter< uint_t >("timesteps") }
      , Re                     { b.getParameter< real_t >("Re") }
      , Fr                     { b.getParameter< real_t >("Fr") }
      , We                     { b.getParameter< real_t >("We") }
      , useHardcodedUnits      { b.getParameter< bool >  ("useHardcodedUnits", false) }
      , relaxationRateOverride { b.getParameter< real_t >("relaxationRateOverride", real_c(0)) }
      , gravityOverride        { b.getParameter< real_t >("gravityOverride", real_c(0)) }
      , surfaceTensionOverride { b.getParameter< real_t >("surfaceTensionOverride", real_c(-1)) }
      , denseNanScanEnabled    { b.getParameter< bool >  ("denseNanScanEnabled", false) }
      , denseNanScanFromStep   { b.getParameter< uint_t >("denseNanScanFromStep", uint_t(0)) }
      , denseNanScanToStep     { b.getParameter< uint_t >("denseNanScanToStep", timesteps) }
      , lowRhoPrecursorThreshold { b.getParameter< real_t >("lowRhoPrecursorThreshold", real_c(0)) }
      , highVelocityPrecursorThreshold { b.getParameter< real_t >("highVelocityPrecursorThreshold", real_c(0)) }
      , populationPrecursorThreshold { b.getParameter< real_t >("populationPrecursorThreshold", real_c(0)) }
      , densityPrecursorLow          { b.getParameter< real_t >("densityPrecursorLow", real_c(0)) }
      , densityPrecursorHigh         { b.getParameter< real_t >("densityPrecursorHigh", real_c(0)) }
      , densityPrecursorLiquidOnly   { b.getParameter< bool >  ("densityPrecursorLiquidOnly", false) }
      , fluidSanityEnabled           { b.getParameter< bool >  ("fluidSanityEnabled", false) }
      , fluidSanityFromStep          { b.getParameter< uint_t >("fluidSanityFromStep", uint_t(0)) }
      , fluidSanityToStep            { b.getParameter< uint_t >("fluidSanityToStep", timesteps) }
      , fluidSanityEveryNSteps       { b.getParameter< uint_t >("fluidSanityEveryNSteps", uint_t(1)) }
      , fluidSanityDensityMin        { b.getParameter< real_t >("fluidSanityDensityMin", real_t(0.1)) }
      , fluidSanityDensityMax        { b.getParameter< real_t >("fluidSanityDensityMax", real_t(1.9)) }
      , fluidSanityMaxVelocity       { b.getParameter< real_t >("fluidSanityMaxVelocity", real_t(0.17)) }
      , fluidSanityMaxPopulation     { b.getParameter< real_t >("fluidSanityMaxPopulation", real_t(2.0)) }
      , fluidSanityFillTolerance     { b.getParameter< real_t >("fluidSanityFillTolerance", real_t(0.05)) }
      , fluidSanityLiquidOnly        { b.getParameter< bool >  ("fluidSanityLiquidOnly", false) }
      , fluidSanityAbortAbove        { b.getParameter< uint_t >("fluidSanityAbortAbove", uint_t(0)) }
      , fluidSanityMaxExcessMass     { b.getParameter< real_t >("fluidSanityMaxExcessMass", real_t(0)) }
      , fluidSanityMaxNeighbourVelocityJump { b.getParameter< real_t >("fluidSanityMaxNeighbourVelocityJump", real_t(0)) }
      , fluidSanityMaxNeighbourDensityJump  { b.getParameter< real_t >("fluidSanityMaxNeighbourDensityJump", real_t(0)) }
      , fluidSanityCheckTopology     { b.getParameter< bool >  ("fluidSanityCheckTopology", false) }
      , fluidSanityStatsEveryNSteps  { b.getParameter< uint_t >("fluidSanityStatsEveryNSteps", uint_t(0)) }
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
   real_t      floaterRingBand;   // sphere V_sub plane fit: keep interface cells with |dist-r| < band. Default 1.5
                                  // (contact ring); a huge value (e.g. 1e6) samples the full box.
   // Both > 0: the SI scaling skips inFileName's header (getAvgDiameterScalingFactor, Utility.h) and uses
   // avgParticleDiameterOverride [LU] (sets liquidHeight, floorHeight, particleFixingHeight) and dxSIOverride [m/LU]
   // (sets dx_SI, dt_SI) directly, e.g. for floater-only scenarios. Either alone is ignored. Default 0.
   real_t      avgParticleDiameterOverride;
   real_t      dxSIOverride;

   // DEM collision time tc = collisionTimeFactor * avgParticleDiameter [LBM steps], default 4. The spring-dashpot
   // stiffness k = m_eff (pi^2 + ln^2 e) / tc^2 scales as 1/collisionTimeFactor^2 and is independent of numSubCycles.
   // Smaller factor = stiffer contacts, less overlap; keep tc * numSubCycles >~ 20-40 DEM steps for stability.
   real_t      collisionTimeFactor;

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
      , collisionTimeFactor         { b.getParameter< real_t >  ("collisionTimeFactor", real_c(4)) }
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

// Vertical round container (DropTest only): a CylindricalBoundary in containment mode, axis along z through `center`
// (x, y used), bodies stay inside `radius`. Reproduces the container of Kodam et al. (2010b, Sec. 3.3, 50.6 mm inner
// diameter). Optional prm block: ContainerCylinders { ContainerCylinder { center < x, y, z >; radius r; } }
struct ContainerCylinderParams {
   Vector3< real_t > center;
   real_t            radius;

   ContainerCylinderParams(const config::Config::BlockHandle& b)
   {
      center = b.getParameter< Vector3< real_t > >("center");
      radius = b.getParameter< real_t >           ("radius");
   }
};

struct FloatingSphereParams {
   Vector3< real_t > position;
   real_t            radius;
   real_t            densityRatio;
   // LBM steps the sphere is held fixed before release, so the FSLBM can settle the interface around a body started at
   // or above the free surface. 0 = release immediately (default).
   uint_t            settlingSteps;
   // Initial linear velocity [LU/step], e.g. for collision tests. Default <0,0,0>; only effective with settlingSteps 0.
   Vector3< real_t > initialVelocity;
   // pinPosition/pinAxes: hold the translational DOFs selected by pinAxes (nonzero = pinned, default <1,1,1>) at
   // `position` every DEM sub-cycle; rotation stays free. <1,1,0> lets a body find its buoyant depth without lateral
   // drift. Default false.
   bool              pinPosition;
   Vector3< real_t > pinAxes;
   // pinRotation/pinRotationAxes: zero the world-frame angular velocity components selected by pinRotationAxes
   // (default <1,1,1>) every DEM sub-cycle; independent of pinPosition. For symmetry with boxes/cylinders.
   // Default false.
   bool              pinRotation;
   Vector3< real_t > pinRotationAxes;
   // true -> dormant (no fluid mapping, no DEM contact) until settlingSteps elapse, so several floaters can share one
   // spawn position and be released one at a time. Requires settlingSteps > 0 (warned otherwise). Default false.
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

inline std::vector< ContainerCylinderParams > parseContainerCylinders(const shared_ptr< config::Config >& cfg)
{
   std::vector< ContainerCylinderParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("ContainerCylinders", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks blocks;
   containers[0].getBlocks("ContainerCylinder", blocks, 0, 10);
   for (const auto& b : blocks)
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
   // Axis of the initial rotation (default world Y); with a random axis and angle the box starts in any orientation.
   Vector3< real_t > initialTiltAxis;
   real_t            comOffsetZ;          // body-frame z-offset of the COM from the centroid [LU]; negative = COM
                                          // below the centroid (righting torque). Default 0 (uniform density).
   bool              pinPosition;         // true -> translation held at `position` every DEM sub-cycle, rotation
                                          // free (unlike the FIXED flag). Default false.
   Vector3< real_t > pinAxes;             // pin mask for pinPosition, nonzero = pinned. Default <1,1,1>; <1,1,0>
                                          // lets the box find its buoyant depth without lateral drift.
   bool              pinRotation;         // true -> the angular velocity components selected by pinRotationAxes
                                          // are zeroed every DEM sub-cycle. Default false.
   Vector3< real_t > pinRotationAxes;     // world-axis mask for pinRotation. Default <1,1,1> (no rotation);
                                          // <1,0,1> allows rotation about Y only.
   bool              startInactive;      // see FloatingSphereParams::startInactive.

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
      initialTiltAxis      = b.getParameter< Vector3< real_t > >("initialTiltAxis", Vector3< real_t >(real_t(0), real_t(1), real_t(0)));
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
   bool              startInactive;       // see FloatingSphereParams::startInactive.

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

// Floating capsule (mesa_pd::data::Capsule, SHAPE_TYPE 7): a cylinder with hemispherical caps. `length` is the full
// end-to-end length including the caps, as in FloatingCylinderParams. There is no closed-form contact routine: every
// capsule pair goes through GeneralContactDetection's GJK/EPA.
struct FloatingCapsuleParams {
   Vector3< real_t > position;
   real_t            radius;
   real_t            length;              // full end-to-end length INCLUDING the hemispherical caps
   Vector3< real_t > axis;                // initial world direction of the capsule axis (normalised on parse)
   real_t            densityRatio;
   uint_t            settlingSteps;       // same semantics as FloatingSphereParams::settlingSteps
   real_t            initialTiltAngleDeg; // extra rotation about the world Y-axis, applied AFTER axis alignment
   bool              pinPosition;         // hold TRANSLATION at `position`, rotation free (same as box/cylinder)
   Vector3< real_t > pinAxes;             // per-axis pin mask (same as box/cylinder)
   bool              pinRotation;         // per-axis rotation pin, zero angular velocity (same as box/cylinder)
   Vector3< real_t > pinRotationAxes;     // per-axis rotation pin mask (same as box/cylinder)
   bool              startInactive;       // see FloatingSphereParams::startInactive.

   FloatingCapsuleParams(const config::Config::BlockHandle& b)
   {
      position             = b.getParameter< Vector3< real_t > >("position");
      radius               = b.getParameter< real_t >           ("radius");
      length               = b.getParameter< real_t >           ("length");
      axis                 = b.getParameter< Vector3< real_t > >("axis", Vector3< real_t >(real_t(0), real_t(1), real_t(0)));
      densityRatio         = b.getParameter< real_t >           ("densityRatio");
      settlingSteps        = b.getParameter< uint_t >           ("settlingSteps", uint_t(0));
      initialTiltAngleDeg  = b.getParameter< real_t >           ("initialTiltAngleDeg", real_t(0));
      pinPosition          = b.getParameter< bool >             ("pinPosition", false);
      pinAxes              = b.getParameter< Vector3< real_t > >("pinAxes", Vector3< real_t >(real_t(1)));
      pinRotation          = b.getParameter< bool >             ("pinRotation", false);
      pinRotationAxes      = b.getParameter< Vector3< real_t > >("pinRotationAxes", Vector3< real_t >(real_t(1)));
      startInactive        = b.getParameter< bool >             ("startInactive", false);
      if (startInactive && settlingSteps == uint_t(0))
         WALBERLA_LOG_WARNING("FloatingCapsule startInactive=true but settlingSteps=0 -- this floater "
                              "will never be activated (releaseSettledFloatingSpheres no-ops when "
                              "settlingSteps==0).");
      WALBERLA_CHECK_GREATER(axis.length(), real_t(1e-10), "FloatingCapsule axis must be nonzero!");
      axis = axis.getNormalized();
   }
};

inline std::vector< FloatingCapsuleParams > parseFloatingCapsules(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingCapsuleParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingCapsules", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks capsuleBlocks;
   containers[0].getBlocks("FloatingCapsule", capsuleBlocks, 0, 10000); // matches FloatingSphere's own cap
   for (const auto& b : capsuleBlocks)
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
   // true (default) -> Bogner & Rude (2013) per-cell gas/liquid/interface selection for freshly uncovered (FormerMO)
   //         cells. false -> every uncovered cell becomes an interface cell.
   bool        useBognerReconstruction;

   // true (default) -> fill-weighted MovingObstacle (SimpleBB) MEM force (Bogner 2017, Eq. 5.9): interface cells at a
   //         floater's waterline contribute in proportion to their fill, suppressing phantom buoyant lift.
   // false -> null fill accessor, raw two-phase MEM force (ablation).
   bool        useFillWeighting;

   // true -> interpolated bounce-back (CLI, Ginzburg et al. 2008; lbm_mesapd_coupling::CurvedLinear) on MovingObstacle
   //         links whose fluid neighbour is liquid/interface, with the exact wall distance from a ray-particle
   //         intersection; other links use simple bounce-back (SimpleBB.h).
   // false (default) -> simple bounce-back on every link (published Antidunes/Bogner behaviour).
   bool        useInterpolatedBounceBack;

   // Experimental (MovingParticleMapping.h, cellCenterTestOffset): the cell-centre containment test uses the centre
   // shifted by mappingTieOffset*(dx, 1.3 dy, 1.7 dz), so a body face lying exactly on cell centres is mapped the same
   // way every step. 0 = off (default); 1e-4 recommended.
   real_t      mappingTieOffset;

   // true (default) -> no hydrostatic reference density is subtracted in the MEM force (SimpleBB gets a null
   //         hydrostaticDensityFct). With useFillWeighting this is the three-phase MEM force of Bogner & Rude (2013);
   //         the reference height H then has no effect on the MEM force. Pair with useMemOnlyBuoyancy=true: false here
   //         together with useMemOnlyBuoyancy=true leaves floaters without buoyant support.
   // false -> subtract rhoHydStat relative to H (staticReferenceHeight / useDynamicReferenceHeight; extended,
   //         non-thesis behaviour).
   bool        useBognerThreePhaseMEM;

   // true (default) -> the vertical body force is gravity only (rho*vFull*a_z); buoyancy comes solely from the MEM
   //         force, as in Bogner & Rude (2013). computeSubmergedVolumes*() are skipped entirely, so there is no
   //         analytic righting torque either. Pair with useBognerThreePhaseMEM=true (see there).
   // false -> gravity minus vSub*g plus the analytic vSub-based righting torque (extended, non-thesis split).
   bool        useMemOnlyBuoyancy;

   // Streamwise/lateral drive of floaters (types 1/3/4) in AddExternalForces: the periodic channel is driven by a body
   // force a_x, a_y, which the MEM does not see as a pressure gradient.
   // false (default) -> F_xy = vFull * a_xy (original Antidunes term, derived for fully submerged sediment).
   // true -> F_xy = rho_p * vFull * a_xy = m_p * a_xy (slope reading of the body force, Schwarzmeier et al. 2023),
   //         independent of submergence. Non-floaters are not affected.
   bool        driveFloatersByMass;

   // Only read when useMemOnlyBuoyancy=false. true -> vSub = vFull and no righting lever for every floater, skipping
   //         the free-surface plane fit (valid for fully submerged floaters; MEM vs. analytic buoyancy comparison).
   // false (default) -> plane fit / clip / sampling.
   bool        useFullVolumeAnalyticBuoyancy;

   // true -> re-seed each FormerMO cell to the hydrostatic density rho_hyd(z) = 1 + 3*a_z*(z - H) instead of the
   //         EquilibriumReconstructor's fill-blind mean (interface cells report rho ~ 1 however empty), removing the
   //         over-density at a moving floater's waterline. Does not use the fill. Default false.
   bool        useHydrostaticReconstruction;

   // true -> after Pass 2 of the FormerMO reconstruction (PdfReconstructionManager), cells that fell back to interface
   //         only because a same-pass FormerMO neighbour was unresolved are re-resolved against the now-known
   //         neighbourhood (cheap local sweep). false (default) -> they stay interface.
   bool        usePass3Reresolution;

   // true -> a freshly uncovered cell that becomes interface gets the mean volume fraction of all resolved neighbours
   //         (gas 0, liquid 1, interface its fill) as fill (PdfReconstructionManager.h,
   //         seedInterfaceFillFromNeighbourhood). false (default) -> mean fill of the interface neighbours only, 0.5
   //         if none; this seeds half-full cells inside the water or in the air when a cell is interface only because
   //         of an unclassified same-pass neighbour. Changes only the fill, never the cell state.
   bool        seedFillFromAllNeighbours;

   // true (default) -> the seeding above ignores neighbours flagged obstacle / formerObstacle. false -> an
   //         interface-flagged neighbour counts even if it also carries such a flag.
   bool        seedSkipObstacleNeighbours;

   // true -> freshly uncovered cells follow Bogner's rule literally (Bogner 2017, pp. 113-114): B(x) = non-obstacle
   //         neighbours, cells freed in the same step counted as obstacle; no liquid -> gas, else no gas -> liquid,
   //         else interface (fill = mean of the interface cells in B(x), 0.5 if none). No same-step interface fallback
   //         (PdfReconstructionManager.h, classifyStrictBognerNeighbourhood). false (default) -> the fallback rule.
   bool        strictBognerClassification;
   // true -> (with strictBognerClassification) a freed cell decided liquid that touches gas after the step's decisions
   //         are applied becomes interface with fill 1 (PdfReconstructionManager.h, Pass 2). false (default) -> only
   //         counted ([FORMERMO-STRICT]).
   bool        strictBognerNoLiquidAtGas;
   // true -> (needs strictBognerClassification and strictBognerNoLiquidAtGas) a third local pass re-scans every cell
   //         that is now interface against the fully decided flag field and demotes it to gas or liquid if its D3Q27
   //         neighbourhood shows only one phase (the promotions of the previous pass can remove a neighbour's only
   //         liquid evidence). Never leaves liquid touching gas (LiquidPromoteRelaxClassificationTest.cpp), unlike
   //         usePass3Reresolution, which sees only frozen Pass-2 values. No-op otherwise. Default false.
   bool        strictBognerRelaxInterface;
   // true -> log [MOMENTUM-BUDGET] every evaluationFrequency steps: momentum carried into the fluid by freed cells,
   //         fluid momentum in newly covered cells, total fluid momentum (liquid + fill-weighted interface).
   //         Diagnostic only. Default false.
   bool        logMomentumBudget;

   // true -> finite-difference curvature (CurvatureSweepFiniteDifferences) skips neighbour directions whose flag is in
   //         no group (upstream reads an uninitialised normal there) or whose normal is non-finite, sets curvature 0 if
   //         none remain, and logs [CURVATURE-GUARD]. Identical results for valid inputs. false (default) -> upstream.
   bool        guardCurvatureInputs;
   // true -> after flag initialisation, the outer ghost layer (layer 2) beyond every non-periodic domain border also
   //         becomes wall (same NoSlip/FreeSlip flag as the layer-1 wall, fill 1), i.e. walls are two cells thick.
   //         Otherwise layer 2 keeps the liquid/interface flags set by fillWithDomain(nrOfGhostLayers), is never reset
   //         or communicated, and ObstacleFillLevelSweep averages its fill into the wall cells' smoothed fill, which
   //         can give NaN normals at side walls ([GHOST2-FLUID] check, BuildTimeloop.h). false (default) -> unchanged.
   bool        wallOuterGhostLayer;

   // true -> after the normal sweep, scan for non-finite interface normals, smoothed fill levels and obstacle normals;
   //         log [NAN-NORMAL] / [NAN-SMOOTHFILL] / [NAN-OBSTACLENORMAL] with their inputs and a global count every 1000
   //         steps (BuildTimeloop.h, diagnoseNaNNormals). Diagnostic only, full field scan per step. Default false.
   bool        diagnoseNaNNormals;

   // true -> subtract the analytic phantom force F_pred_z = dP * A_proj * taper (pile-up head) from the vertical MEM
   //         force on each floating sphere (particle side only, no PDF feedback); vanishes without pile-up. Applied
   //         only when useMemOnlyBuoyancy=false (AddExternalForces). false (default) -> no correction.
   bool        usePhantomForceCorrection;

   // true -> the MEM hydrostatic reference height H (SimpleBB's rhoHydStat, force side only) is tracked live: each step
   //         the local waterlines of all floaters with a clean fitted ring are averaged and EMA-blended (alpha 0.2)
   //         into one shared H. An offset correction only; it cannot remove velocity-dependent residuals. No effect on
   //         the MEM force while useBognerThreePhaseMEM=true. false (default) -> H = absoluteLiquidHeight.
   bool        useDynamicReferenceHeight;

   // >= 0 -> replaces the static reference height H (SimpleBB's rhoHydStat, the hydrostatic reconstruction target, the
   //         PhantomDecomp flatRef) by this absolute z [LU] instead of absoluteLiquidHeight. Diagnostic ablation:
   //         shifts the subtracted reference by drho = 3*|a_z|*dH, which a fully submerged body does not feel but a
   //         floater does. Fluid initialisation unaffected. < 0 (default -1) -> disabled.
   real_t      staticReferenceHeightOverride;

   // Lateral reach of the free-surface plane fit for a sphere floater's localSurfaceZ (and the dynamic-H target), as a
   // multiple of its radius. Default 1.8 reaches past a calmly floating body's own crater but not past the splash of a
   // fast impact; raise it for violent entries.
   real_t      floaterPlaneFitRadiusFactor;

   // Same for box/cylinder floaters (computeSubmergedVolumesBox): band = boxPlaneFitBandFactor * max(halfExtents) + 2.
   // Default 0.8; raise it if the fitted plane is contaminated by the body's own crater/meniscus.
   real_t      boxPlaneFitBandFactor;

   // true -> Box-Box and Cylinder-Cylinder contacts use the closed forms in ClosedFormContacts*.h (OBB separating-axis
   //         test of Gottschalk, Lin & Manocha 1996 for Box-Box; Kodam et al. 2010 with the corrections of Guo et al.
   //         2012 for Cylinder-Cylinder) instead of GJK/EPA. Cylinder-Cylinder is hybrid: when
   //         cylinderCylinderResultIsUncertain() flags a result (configurations the closed forms do not cover
   //         cleanly) it falls back to GJK/EPA (see PaperAnalyticalContactTest.cpp).
   // false (default) -> GJK/EPA for both pairs.
   bool        useAnalyticalContactDetection;

   // true -> (only with useAnalyticalContactDetection) no Cylinder-Cylinder GJK/EPA fallback: the Kodam + Guo dispatch
   //         is the sole answer, as in the papers (tier 1 "no contact" accepted, tier 2 contacts applied as-is). The
   //         raw centre-distance pre-filter in SolidCylinderContactDetection stays active.
   // false (default) -> fall back to GJK/EPA whenever cylinderCylinderResultIsUncertain() flags the result.
   bool        usePureAnalyticalCylinder;

   // true -> Cylinder-Cylinder contacts (only with useAnalyticalContactDetection) are refined by minimising the
   //         projection overlap r_i(L) + r_j(L) - L.(x_j - x_i) over directions L
   //         (closed_form_contacts::refineCylinderContact): exact closed-form scenarios are unchanged, inexact ones get
   //         the minimising depth/normal, "no scenario" pairs are decided by the same descent; no GJK/EPA call. Also
   //         used for rod vs. solid cylindrical pillar. Validated in PaperAnalyticalContactTest.cpp (depth error
   //         <= 5e-4; separated pairs closer than ~0.01 may be reported as grazing contacts, overlap <= 0.006).
   // false (default) -> GJK/EPA fallback for "no scenario" and special contacts.
   bool        refineCylinderContacts;

   // true -> cylinderCylinderContact (ClosedFormContactsCylinderCylinder.h) passes a band-edge result to
   //         refineCylinderContact when the projection overlap D(n) along its normal exceeds
   //         kBandEdgeConsistencyFactor (3) times its depth (closed_form_contacts::kodamBandEdgeAxialCheck); fixes
   //         stacked equal-radius cylinders in the round-container drop test. false (default) -> unchanged dispatch.
   bool        kodamBandEdgeAxialCheck;

   // true (default) -> the cheap pre-filter before contact detection uses tight, orientation-aware bounds instead of
   //         circumscribing spheres. Cylinder x Cylinder: distance between the axis segments vs. the sum of the radii
   //         (closestDistanceSqBetweenSegments). Box x Box: projection of both boxes onto the line between the
   //         centres (boxesDefinitelySeparatedAlongCenterLine, one SAT axis instead of 15) ahead of boxBoxSAT,
   //         effective only with useAnalyticalContactDetection. Both are conservative (PaperAnalyticalContactTest.cpp,
   //         testTightLayer1BoxAccuracy10k); derivations in ClosedFormContacts*.h.
   // false -> circumscribing-sphere reject for Cylinder x Cylinder, no pre-filter for Box x Box.
   bool        useTightLayer1;

   // Experimental: build cylinder floaters (FloatingCylinder entries) as Capsules (SHAPE_TYPE 7, hemispherical caps)
   // instead of flat-capped Cylinders, avoiding flat-cap/rim contact edge cases. Capsule pairs always use GJK/EPA,
   // regardless of useAnalyticalContactDetection. Geometry parameters unchanged; only the shape built in SetupStages.h
   // differs. Default false.
   bool        useCapsuleShapeForCylinders;

   // Diagnostic caps for cylinder floaters (type 4), all 0 = off (default). They damp the rod dynamics, change the
   // physics and are not for production. First hits per rank are logged ([ROD-CAP ...]); counts and the largest
   // uncapped value appear in the [RodCaps] line of [ContactDetectionStats].
   //  rodCapHydroForce/rodCapHydroTorque: |F_hydro| / |T_hydro| scaled down to this value (non-finite values zeroed);
   //  rodCapVelocity/rodCapAngularVelocity: |v| / |omega| scaled down after every DEM sub-cycle's velocity update;
   //  rodCapOverlap: penetration depth fed into the cylinder contact force law limited to this value [LU].
   real_t      rodCapHydroForce;
   real_t      rodCapHydroTorque;
   real_t      rodCapVelocity;
   real_t      rodCapAngularVelocity;
   real_t      rodCapOverlap;

   // Test option (all 0 = off, default): prescribed rigid motion of all rods (type 4). After every DEM sub-cycle's
   // velocity update the linear velocity is set to this vector and the angular velocity to 0 (forces and torques
   // ignored), so moving rods cover/uncover cells without DEM dynamics. Changes the physics; not for production.
   real_t      rodPrescribedVelocityX;   // three scalars: the config proxy cannot brace-initialise a Vector3 member
   real_t      rodPrescribedVelocityY;
   real_t      rodPrescribedVelocityZ;

   // true -> shadow mode for Cylinder x Cylinder, meaningful only with useAnalyticalContactDetection=false: every rod
   //         pair in contact also gets the closed-form result computed (never used for forces): chosen scenario,
   //         overlap/normal, the GJK/EPA result and the self-consistency check of CylinderClosedFormRangeTest.cpp.
   //         Per-scenario statistics go to [ShadowAnalytic] (end of run and with [ContactDetectionStats]); a per-rank
   //         CSV (shadow-analytic/rank_<N>.csv, one row per pair and sub-cycle) is written for main steps in
   //         [shadowAnalyticStartStep, shadowAnalyticEndStep] (EndStep 0 = no CSV). false (default) -> off.
   bool        shadowAnalyticCylinderLog;
   uint_t      shadowAnalyticStartStep;
   uint_t      shadowAnalyticEndStep;

   // true -> the normal stiffness/damping of LinearSpringDashpot::setStiffnessAndDamping use a generalised effective
   //         mass including the rotational inertia at the contact point:
   //             1/meff = 1/mA + 1/mB + (rA x n)^T IA_world^-1 (rA x n) + (rB x n)^T IB_world^-1 (rB x n)
   //         (ParticleKernels.h, computeRotationalEffectiveMass; Baraff, "Physically Based Modeling: Rigid Body
   //         Simulation"). Exact for sphere-sphere, where the rotational terms vanish. With the translational meff the
   //         requested restitution holds only for the closing velocity; corner/edge contacts of boxes/cylinders put
   //         energy into rotation that the damping was not sized for.
   // false (default) -> translational reduced mass 1/(1/mA + 1/mB).
   bool        useRotationalEffectiveMass;

   // Contact skin m [LU] of the cylinder floaters (type 4), default 0 = off. For contact detection only, every
   // floating cylinder is inflated to radius R + m and half-length h + m, so the contact force starts at a surface gap
   // of 2m between two rods and of m between a rod and a wall plane or pillar (at the rims the sharp inflated edge
   // reaches up to (sqrt 2 - 1) m further than a rounded skin). Purpose: keep at least one resolved fluid cell between
   // touching rods, whose line contacts otherwise seal sub-cell liquid channels. Same idea as the "effective diameter"
   // of Boutt et al. (2007, Sec. 2.3), who shrink the fluid body instead (cited for this by Rettinger 2023, p. 7).
   // The fluid mapping is not dilated. Implementation (ParticleKernels.h, SolidCylinderContactDetection): the
   // Cylinder x Cylinder, Cylinder x HalfSpace and Cylinder x CylindricalBoundary overloads call themselves with the
   // inflated shapes, reusing the validated detection (a GJK/EPA margin would not work: mesa_pd's EPA gives up when
   // only the margins overlap, EPA.cpp "numPoints == 1"). The interaction radius (cheap reject, ghost sync range)
   // grows by 2m (SetupStages.h), which bounds the inflated circumscribing radius sqrt((R+m)^2 + (h+m)^2).
   // Not applied to Cylinder x Sphere / Cylinder x Box.
   real_t      cylinderContactSkin;

   // Box analogue of cylinderContactSkin for FloatingBox floaters (type 3), default 0 = off: for contact detection
   // only, each half-extent is inflated by m (contact at a gap of 2m box-box, m box-wall/pillar). The Box x Box,
   // Box x HalfSpace and Box x CylindricalBoundary overloads reuse the unchanged detection (SAT or GJK/EPA) on the
   // inflated boxes; the interaction radius grows by 2m (SetupStages.h). Not applied to Box x Sphere / Box x Cylinder.
   real_t      boxContactSkin;

   // true -> every rank writes floaterHistoryBaseFolder/rank_<N>.csv (truncated at start), one row per contact or
   //         lubrication pair it computed and accepted (passed mesa_pd::mpi::ContactFilter), for every DEM sub-cycle.
   //         Rows hold the computing rank, both uids/isGhost and ac.getOwner(idx) of both particles (comparing the
   //         ranks' CSVs for one (uidA, uidB, mainStep, subCycle) reveals double-accepted contacts), the contact
   //         point/normal/penetration depth and the force delta of that pair. Rows are flushed immediately so they
   //         survive an abort. Heavy (one row per pair per sub-cycle); default false.
   bool        logFloaterHistory;

   // Folder of logFloaterHistory's CSVs (default "floater-history"), created by Antidunes.cpp when enabled.
   std::string floaterHistoryBaseFolder;

   // 0 (default) -> logFloaterHistory (per-contact log() and per-floater logState()) covers every floater. Nonzero ->
   // only floaters whose settlingSteps equals this value (one staggered-release batch), filtered by uid
   // (FloaterHistoryLogger::setUidFilter, ParticleKernels.h). Affects logging only.
   uint_t      logFloaterHistoryReleaseStepFilter;

   // 0 (default) -> logFloaterHistory writes from step 0; nonzero -> rows before this main step are skipped.
   uint_t      logFloaterHistoryStartStep;

   // 0 (default) -> logFloaterHistory writes to the end of the run; nonzero -> rows after this main step are skipped.
   uint_t      logFloaterHistoryEndStep;

   // true -> every rank writes fluidBoxHistoryBaseFolder/rank_<N>.csv with one row per obstacle or interface cell
   //         inside the world box [fluidBoxMin, fluidBoxMax], once per main step in [logFluidBoxHistoryStartStep,
   //         logFluidBoxHistoryEndStep]. Columns: mainStep, rank, world x/y/z, flag, owning particle uid, density,
   //         fill, max|f_i|. Runs inside the diagnostic sweep, so it also needs diagnosticSweepEnabled (and that
   //         sweep's step window). Default false (zero cost).
   bool        logFluidBoxHistory;

   // Folder of logFluidBoxHistory's CSVs (default "fluid-box-history"), created by Antidunes.cpp when enabled.
   std::string fluidBoxHistoryBaseFolder;

   // World-space corners of the logFluidBoxHistory box (whole cells, inclusive). Default (0,0,0) for both, an empty
   // box: nothing is logged until set.
   Vector3< real_t > fluidBoxMin;
   Vector3< real_t > fluidBoxMax;

   // 0 (default) -> from step 0; nonzero -> main steps before this value are skipped.
   uint_t      logFluidBoxHistoryStartStep;

   // 0 (default) -> to the end of the run; nonzero -> main steps after this value are skipped.
   uint_t      logFluidBoxHistoryEndStep;

   // true -> the ranks owning the world column (densityColumnX, densityColumnY) write densityColumnBaseFolder/
   //         rank_<N>.csv with one row per cell of that column (all flags), once per main step in
   //         [logDensityColumnStartStep, logDensityColumnEndStep], to check the hydrostatic stratification of the
   //         settled fluid. Columns: mainStep, rank, worldX, worldY, worldZ, flag, density, fill. Default false.
   bool        logDensityColumn;

   // Folder of logDensityColumn's CSVs (default "density-column"), created by Antidunes.cpp when enabled.
   std::string densityColumnBaseFolder;

   // World x/y of the logged column; choose a point far from any floater. Default 0/0 (must be set).
   real_t      densityColumnX;
   real_t      densityColumnY;

   // 0 (default) -> from step 0; nonzero -> main steps before this value are skipped.
   uint_t      logDensityColumnStartStep;

   // 0 (default) -> to the end of the run; nonzero -> main steps after this value are skipped.
   uint_t      logDensityColumnEndStep;

   // true -> free-surface height from the fill integral (sum of fill over liquid + interface cells along z) on a
   // waterHeightGridX x waterHeightGridY grid of columns (default 0 x 0), every waterHeightEveryNSteps (default 250)
   // main steps over the whole run. Per-rank CSV rows (mainStep, rank, gx, gy, partialFillSum, nObstacleCells) in
   // waterHeightBaseFolder (default "water-height"); a column's height is the sum of partialFillSum over the ranks,
   // nObstacleCells > 0 marks a column under or inside a moving body. Checks fluid-mass conservation. Default false.
   bool        logWaterHeight;

   // true -> after every collide, obstacle cells whose PDF sum is <= 0 or non-finite are reset to rest equilibrium
   // (SetupStages.h, repairDegenerateObstacleCells). Default false; required by the hopper/rod-bundle clogging prms.
   bool        repairDegenerateObstaclePdfs;

   // true -> seed the scratch PDF clone of the stream-reconstruct-advect sweep
   // (src/lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h) from the real PDF field on its first creation (it
   // is otherwise uninitialised memory). Affects only obstacle-cell PDFs, never fill/mass. Default false (upstream).
   bool        seedStreamPdfCloneOnFirstUse;
   std::string waterHeightBaseFolder;
   uint_t      waterHeightGridX;
   uint_t      waterHeightGridY;
   uint_t      waterHeightEveryNSteps;

   // true -> one sweep per main step on timeloopEval (BuildTimeloop.h, addAfterTimeStepDiagnosticSweep, after DEM and
   //         everything else) that (1) drains and logs SimulationObjects::diagnosticEvents, filled by per-sub-cycle
   //         checks in ParticleKernels.h; (2) logs every liquid/interface cell with density outside
   //         [densityPrecursorLow, densityPrecursorHigh] (same thresholds and liquidOnly flag as the dense NaN scan);
   //         (3) logs every local floater's position, quaternion and hydrodynamic/contact force and torque.
   //         Default false.
   bool        diagnosticSweepEnabled;

   // true -> MovingParticleMapping gets BuildTimeloop.h's onOwnershipHandoff callback: right before a cell's obstacle
   //         ownership passes from one particle to another, its PDFs are reset to rest equilibrium (rho 1, u 0), so
   //         the cell does not inherit the previous owner's stale PDF history (a source of growing checkerboard
   //         instabilities between touching rods). false (default) -> no callback, the cell keeps its PDFs.
   bool        useOwnershipHandoffPdfReset;

   // 0 (default) -> diagnosticSweepEnabled runs from step 0; nonzero -> main steps before this value are skipped.
   uint_t      diagnosticSweepFromStep;

   // 0 (default) -> to the end of the run; nonzero -> main steps at or after this value are skipped.
   uint_t      diagnosticSweepToStep;

   // true -> three WALBERLA_LOG_WARNING watchers in ComputeContacts, active in [contactAnomalyStartStep,
   //         contactAnomalyEndStep]:
   //           [NORMAL-FLIP]           the same uid pair's contact normal turns by more than ~45 degrees between
   //                                   sub-cycles (SAT/Kodam scenario switch; impossible for sphere-sphere);
   //           [PROLONGED-CONTACT]     the same pair stays in contact for many consecutive sub-cycles (episode
   //                                   length, peak overlap);
   //           [SIMULTANEOUS-CONTACTS] one floater touches 3+ bodies in the same sub-cycle (the summed stiffness
   //                                   may exceed what the explicit integrator step was tuned for).
   // false (default) -> off, no cost.
   bool        logContactAnomalies;
   // 0 (DEFAULT) -> from step 0, if logContactAnomalies is on.
   uint_t      contactAnomalyStartStep;
   // 0 (DEFAULT) -> unbounded (through end of run), if logContactAnomalies is on.
   uint_t      contactAnomalyEndStep;

   // true -> every step, warn as soon as a local cell above highAltitudeWaterThresholdZ [LU] is liquid or interface,
   // catching gas-to-liquid conversions above the waterline that no other log reports (CellConversionSweep logs
   // nothing). A healthy run triggers it zero times. Default false / 0.
   bool        highAltitudeWaterScanEnabled;
   real_t      highAltitudeWaterThresholdZ;

   // true -> sets mesa_pd::mpi::ContactFilter::useOwnerRankTiebreak on every ComputeContacts/ComputeLubrication: a
   //         contact is accepted iff myRank == min(owner(idx1), owner(idx2)) instead of the per-rank geometric test on
   //         each rank's own contact point, which can disagree between ranks (see ContactFilter.h) and cause duplicate
   //         contact-history inserts for pairs straddling a rank boundary. false (default) -> geometric test.
   bool        useContactFilterOwnerTiebreak;

   // true (default) -> reduceAndAvgHydroForces() (ParticleKernels.h) averages the hydrodynamic force/torque of the
   //         current and the previous step, F = 0.5 (F(t) + F(t - dt)) (AverageHydrodynamicForceTorqueKernel; Ladd
   //         1994; Watanabe et al. 2023, Eqs. 48-49, for floating objects).
   // false -> DEM sees the raw force of the current step (A/B test).
   bool        useHydrodynamicForceTimeAveraging;

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
      , useInterpolatedBounceBack   { b.getParameter< bool >       ("useInterpolatedBounceBack", false) }
      , mappingTieOffset          { b.getParameter< real_t >    ("mappingTieOffset", real_c(0)) }
      , useBognerThreePhaseMEM      { b.getParameter< bool >       ("useBognerThreePhaseMEM", true) }
      , useMemOnlyBuoyancy          { b.getParameter< bool >       ("useMemOnlyBuoyancy", true) }
      , driveFloatersByMass         { b.getParameter< bool >       ("driveFloatersByMass", false) }
      , useFullVolumeAnalyticBuoyancy { b.getParameter< bool >     ("useFullVolumeAnalyticBuoyancy", false) }
      , useHydrostaticReconstruction { b.getParameter< bool >      ("useHydrostaticReconstruction", false) }
      , usePass3Reresolution        { b.getParameter< bool >       ("usePass3Reresolution", false) }
      , seedFillFromAllNeighbours   { b.getParameter< bool >       ("seedFillFromAllNeighbours", false) }
      , seedSkipObstacleNeighbours  { b.getParameter< bool >       ("seedSkipObstacleNeighbours", true) }
      , strictBognerClassification  { b.getParameter< bool >       ("strictBognerClassification", false) }
      , strictBognerNoLiquidAtGas   { b.getParameter< bool >       ("strictBognerNoLiquidAtGas", false) }
      , strictBognerRelaxInterface  { b.getParameter< bool >       ("strictBognerRelaxInterface", false) }
      , logMomentumBudget           { b.getParameter< bool >       ("logMomentumBudget", false) }
      , guardCurvatureInputs        { b.getParameter< bool >       ("guardCurvatureInputs", false) }
      , wallOuterGhostLayer         { b.getParameter< bool >       ("wallOuterGhostLayer", false) }
      , diagnoseNaNNormals          { b.getParameter< bool >       ("diagnoseNaNNormals", false) }
      , usePhantomForceCorrection   { b.getParameter< bool >       ("usePhantomForceCorrection", false) }
      , useDynamicReferenceHeight   { b.getParameter< bool >       ("useDynamicReferenceHeight", false) }
      , staticReferenceHeightOverride{ b.getParameter< real_t >    ("staticReferenceHeightOverride", real_t(-1)) }
      , floaterPlaneFitRadiusFactor  { b.getParameter< real_t >    ("floaterPlaneFitRadiusFactor", real_t(1.8)) }
      , boxPlaneFitBandFactor        { b.getParameter< real_t >    ("boxPlaneFitBandFactor", real_t(0.8)) }
      , useAnalyticalContactDetection{ b.getParameter< bool >      ("useAnalyticalContactDetection", false) }
      , usePureAnalyticalCylinder   { b.getParameter< bool >       ("usePureAnalyticalCylinder", false) }
      , refineCylinderContacts      { b.getParameter< bool >       ("refineCylinderContacts", false) }
      , kodamBandEdgeAxialCheck     { b.getParameter< bool >       ("kodamBandEdgeAxialCheck", false) }
      , useTightLayer1              { b.getParameter< bool >       ("useTightLayer1", true) }
      , useCapsuleShapeForCylinders { b.getParameter< bool >       ("useCapsuleShapeForCylinders", false) }
      , rodCapHydroForce            { b.getParameter< real_t >     ("rodCapHydroForce", real_t(0)) }
      , rodCapHydroTorque           { b.getParameter< real_t >     ("rodCapHydroTorque", real_t(0)) }
      , rodCapVelocity              { b.getParameter< real_t >     ("rodCapVelocity", real_t(0)) }
      , rodCapAngularVelocity       { b.getParameter< real_t >     ("rodCapAngularVelocity", real_t(0)) }
      , rodCapOverlap               { b.getParameter< real_t >     ("rodCapOverlap", real_t(0)) }
      , rodPrescribedVelocityX      { b.getParameter< real_t >     ("rodPrescribedVelocityX", real_t(0)) }
      , rodPrescribedVelocityY      { b.getParameter< real_t >     ("rodPrescribedVelocityY", real_t(0)) }
      , rodPrescribedVelocityZ      { b.getParameter< real_t >     ("rodPrescribedVelocityZ", real_t(0)) }
      , shadowAnalyticCylinderLog   { b.getParameter< bool >       ("shadowAnalyticCylinderLog", false) }
      , shadowAnalyticStartStep     { b.getParameter< uint_t >     ("shadowAnalyticStartStep", uint_t(0)) }
      , shadowAnalyticEndStep       { b.getParameter< uint_t >     ("shadowAnalyticEndStep", uint_t(0)) }
      , useRotationalEffectiveMass  { b.getParameter< bool >       ("useRotationalEffectiveMass", false) }
      , cylinderContactSkin         { b.getParameter< real_t >     ("cylinderContactSkin", real_t(0)) }
      , boxContactSkin               { b.getParameter< real_t >     ("boxContactSkin", real_t(0)) }
      , logFloaterHistory           { b.getParameter< bool >       ("logFloaterHistory", false) }
      , floaterHistoryBaseFolder    { b.getParameter< std::string >("floaterHistoryBaseFolder", "floater-history") }
      , logFloaterHistoryReleaseStepFilter{ b.getParameter< uint_t >("logFloaterHistoryReleaseStepFilter", uint_t(0)) }
      , logFloaterHistoryStartStep   { b.getParameter< uint_t >     ("logFloaterHistoryStartStep", uint_t(0)) }
      , logFloaterHistoryEndStep     { b.getParameter< uint_t >     ("logFloaterHistoryEndStep", uint_t(0)) }
      , logFluidBoxHistory           { b.getParameter< bool >       ("logFluidBoxHistory", false) }
      , fluidBoxHistoryBaseFolder    { b.getParameter< std::string >("fluidBoxHistoryBaseFolder", "fluid-box-history") }
      // getParameter< Vector3 >(name, default) does not compile for this config binding, hence isDefined() + ternary.
      , fluidBoxMin                  { b.isDefined("fluidBoxMin") ? b.getParameter< Vector3< real_t > >("fluidBoxMin") : Vector3< real_t >(real_t(0)) }
      , fluidBoxMax                  { b.isDefined("fluidBoxMax") ? b.getParameter< Vector3< real_t > >("fluidBoxMax") : Vector3< real_t >(real_t(0)) }
      , logFluidBoxHistoryStartStep  { b.getParameter< uint_t >     ("logFluidBoxHistoryStartStep", uint_t(0)) }
      , logFluidBoxHistoryEndStep    { b.getParameter< uint_t >     ("logFluidBoxHistoryEndStep", uint_t(0)) }
      , logDensityColumn             { b.getParameter< bool >       ("logDensityColumn", false) }
      , densityColumnBaseFolder      { b.getParameter< std::string >("densityColumnBaseFolder", "density-column") }
      , densityColumnX               { b.getParameter< real_t >     ("densityColumnX", real_t(0)) }
      , densityColumnY               { b.getParameter< real_t >     ("densityColumnY", real_t(0)) }
      , logDensityColumnStartStep    { b.getParameter< uint_t >     ("logDensityColumnStartStep", uint_t(0)) }
      , logDensityColumnEndStep      { b.getParameter< uint_t >     ("logDensityColumnEndStep", uint_t(0)) }
      , logWaterHeight               { b.getParameter< bool >       ("logWaterHeight", false) }
      , repairDegenerateObstaclePdfs { b.getParameter< bool >       ("repairDegenerateObstaclePdfs", false) }
      , seedStreamPdfCloneOnFirstUse { b.getParameter< bool >       ("seedStreamPdfCloneOnFirstUse", false) }
      , waterHeightBaseFolder        { b.getParameter< std::string >("waterHeightBaseFolder", "water-height") }
      , waterHeightGridX             { b.getParameter< uint_t >     ("waterHeightGridX", uint_t(0)) }
      , waterHeightGridY             { b.getParameter< uint_t >     ("waterHeightGridY", uint_t(0)) }
      , waterHeightEveryNSteps       { b.getParameter< uint_t >     ("waterHeightEveryNSteps", uint_t(250)) }
      , diagnosticSweepEnabled       { b.getParameter< bool >       ("diagnosticSweepEnabled", false) }
      , useOwnershipHandoffPdfReset  { b.getParameter< bool >       ("useOwnershipHandoffPdfReset", false) }
      , diagnosticSweepFromStep      { b.getParameter< uint_t >     ("diagnosticSweepFromStep", uint_t(0)) }
      , diagnosticSweepToStep        { b.getParameter< uint_t >     ("diagnosticSweepToStep", uint_t(0)) }
      , logContactAnomalies          { b.getParameter< bool >       ("logContactAnomalies", false) }
      , contactAnomalyStartStep      { b.getParameter< uint_t >     ("contactAnomalyStartStep", uint_t(0)) }
      , contactAnomalyEndStep        { b.getParameter< uint_t >     ("contactAnomalyEndStep", uint_t(0)) }
      , highAltitudeWaterScanEnabled { b.getParameter< bool >       ("highAltitudeWaterScanEnabled", false) }
      , highAltitudeWaterThresholdZ  { b.getParameter< real_t >     ("highAltitudeWaterThresholdZ", real_t(0)) }
      , useContactFilterOwnerTiebreak{ b.getParameter< bool >       ("useContactFilterOwnerTiebreak", false) }
      , useHydrodynamicForceTimeAveraging{ b.getParameter< bool >   ("useHydrodynamicForceTimeAveraging", true) }
   {
      if (usePureAnalyticalCylinder && !useAnalyticalContactDetection)
         WALBERLA_LOG_WARNING("usePureAnalyticalCylinder=true but useAnalyticalContactDetection=false -- "
            "this flag only has an effect inside the Cylinder x Cylinder analytical dispatch, which is "
            "off, so it does nothing here.");
      if (logFloaterHistoryEndStep != uint_t(0) && logFloaterHistoryEndStep <= logFloaterHistoryStartStep)
         WALBERLA_LOG_WARNING("logFloaterHistoryEndStep=" << logFloaterHistoryEndStep << " <= "
            "logFloaterHistoryStartStep=" << logFloaterHistoryStartStep << " -- this window is empty, "
            "logFloaterHistory will never write a row regardless of its own enable flag.");
      if (logContactAnomalies && contactAnomalyEndStep != uint_t(0) && contactAnomalyEndStep <= contactAnomalyStartStep)
         WALBERLA_LOG_WARNING("contactAnomalyEndStep=" << contactAnomalyEndStep << " <= "
            "contactAnomalyStartStep=" << contactAnomalyStartStep << " -- this window is empty, "
            "logContactAnomalies' watchers will never fire.");
   }
};

struct EvaluationParams {
   uint_t      performanceLogFrequency;
   uint_t      evaluationFrequency;
   // Frequency of AverageDataSliceEvaluator (BuildTimeloop.h), which writes three text files per call to
   // eval/avgDataSlices/. 0 / unset = evaluationFrequency. A larger value throttles the slice output when
   // evaluationFrequency is small (e.g. per-step force logging).
   uint_t      dataSliceFrequency;
   std::string baseFolderName;
   // false (default) -> ParticleKernels::logContactDetectionStatistics() runs once at the end of the run. true ->
   // also every performanceLogFrequency steps, followed by a reset, so each [ContactDetectionStats] line covers one
   // interval.
   bool        logContactDetectionStatsPeriodically;

   // Optional "ClogMeasurement" block: probe positions for ClogMeasurementEvaluator, which records the water elevation
   // and pressure upstream and downstream of an obstruction (backwater rise). Keys: upstreamPos, downstreamPos,
   // probeHalfWidthXY (default 2). Absent block -> clogMeasurementEnabled false.
   bool               clogMeasurementEnabled;
   Vector3< real_t >  clogUpstreamPos;
   Vector3< real_t >  clogDownstreamPos;
   real_t             clogProbeHalfWidthXY;

   // Optional second probe pair at the pillar's y (keys pillarUpstreamPos / pillarDownstreamPos, clogPillar* below):
   // Schalko, Follett & Nepf (2023, WRR) found the water depth to differ between the blocked region and the open gap
   // of a partial-width obstruction and measured separate profiles at the jam and gap centres.
   // Key flowRateEnabled (default false): also write the volumetric flow rate Q and the wet cross-section area A of
   // the y-z plane at the probes' x (clog{Upstream,Downstream}{FlowRate,WetArea}; ClogMeasurementEvaluator::
   // flowRateAt). First probe pair only.
   bool               clogFlowRateEnabled;
   // Key blockageRatioEnabled (default false): also write the fraction of the reference (undisturbed) wetted cross-
   // section occupied by the pillar(s) and the floater jam at the probes' x (clog{Upstream,Downstream}BlockageRatio
   // and the clog*PillarArea / clog*JamArea components; ClogMeasurementEvaluator::blockageAt), for afflux vs.
   // blockage ratio as in Panici & Kripakaran (2018). First probe pair only; the pillar geometry and channel width
   // come from simParams at the BuildTimeloop.h call site.
   bool               clogBlockageRatioEnabled;
   bool               clogPillarMeasurementEnabled;
   Vector3< real_t >  clogPillarUpstreamPos;
   Vector3< real_t >  clogPillarDownstreamPos;

   EvaluationParams(const config::Config::BlockHandle& b)
      : performanceLogFrequency { b.getParameter< uint_t >     ("performanceLogFrequency") }
      , evaluationFrequency     { b.getParameter< uint_t >     ("evaluationFrequency") }
      , dataSliceFrequency      { b.getParameter< uint_t >     ("dataSliceFrequency", uint_c(0)) == uint_c(0)
                                     ? b.getParameter< uint_t >("evaluationFrequency")
                                     : b.getParameter< uint_t >("dataSliceFrequency") }
      , baseFolderName          { b.getParameter< std::string >("baseFolderName") }
      , logContactDetectionStatsPeriodically{ b.getParameter< bool >("logContactDetectionStatsPeriodically", false) }
      , clogMeasurementEnabled  { static_cast< bool >(b.getBlock("ClogMeasurement")) }
      , clogUpstreamPos         { clogMeasurementEnabled ? b.getOneBlock("ClogMeasurement").getParameter< Vector3< real_t > >("upstreamPos")   : Vector3< real_t >(real_t(0)) }
      , clogDownstreamPos       { clogMeasurementEnabled ? b.getOneBlock("ClogMeasurement").getParameter< Vector3< real_t > >("downstreamPos") : Vector3< real_t >(real_t(0)) }
      , clogProbeHalfWidthXY    { clogMeasurementEnabled ? b.getOneBlock("ClogMeasurement").getParameter< real_t >("probeHalfWidthXY", real_t(2)) : real_t(0) }
      , clogFlowRateEnabled     { clogMeasurementEnabled && b.getOneBlock("ClogMeasurement").getParameter< bool >("flowRateEnabled", false) }
      , clogBlockageRatioEnabled { clogMeasurementEnabled && b.getOneBlock("ClogMeasurement").getParameter< bool >("blockageRatioEnabled", false) }
      , clogPillarMeasurementEnabled { clogMeasurementEnabled && b.getOneBlock("ClogMeasurement").isDefined("pillarUpstreamPos") }
      , clogPillarUpstreamPos   { clogPillarMeasurementEnabled ? b.getOneBlock("ClogMeasurement").getParameter< Vector3< real_t > >("pillarUpstreamPos")   : Vector3< real_t >(real_t(0)) }
      , clogPillarDownstreamPos { clogPillarMeasurementEnabled ? b.getOneBlock("ClogMeasurement").getParameter< Vector3< real_t > >("pillarDownstreamPos") : Vector3< real_t >(real_t(0)) }
   {}
};

// Parses a comma-separated list of unsigned integers, e.g. "1000, 2000, 8000,16000" -> {1000, 2000, 8000, 16000}
// (whitespace ignored; std::set sorts and de-duplicates). Used by VTKParams::snapshotSteps.
// getParameter<std::string>() returns a quoted prm value with its quote characters, which would make the first entry
// fail to parse, so they are stripped first.
inline std::set< uint_t > parseUintList(const std::string& csvIn)
{
   std::string csv = csvIn;
   if (!csv.empty() && csv.front() == '"') csv.erase(0, 1);
   if (!csv.empty() && csv.back() == '"') csv.pop_back();

   std::set< uint_t > result;
   std::istringstream iss(csv);
   std::string token;
   while (std::getline(iss, token, ','))
   {
      std::istringstream tokenStream(token);
      uint_t value;
      if (tokenStream >> value) result.insert(value);
   }
   return result;
}

struct VTKParams {
   uint_t            writeFrequency;
   uint_t            ghostLayers;
   std::string       baseFolder;
   real_t            samplingResolution;

   // Optional VTK block "fluid_field_dense": extra fluid_field output executed only in [denseFromStep, denseToStep)
   // (VtkWriter.h, addVTKOutput), e.g. at a higher frequency around an expected crash. Absent -> denseEnabled false.
   bool              denseEnabled;
   uint_t            denseFromStep;
   uint_t            denseToStep;

   // Optional VTK block "fluid_field_snapshots": fluid_field output only at the listed steps (key "steps", e.g.
   // "0,1000,2000,4000,8000"; VtkWriter.h, addVTKOutput). Its writeFrequency is unused (vtk::initializeVTKOutput
   // requires it; leave it at 1). Absent -> snapshotsEnabled false.
   bool               snapshotsEnabled;
   std::set< uint_t > snapshotSteps;

   VTKParams(const config::Config::BlockHandle& b)
      : writeFrequency    { b.getOneBlock("fluid_field").getParameter< uint_t >     ("writeFrequency") }
      , ghostLayers       { b.getOneBlock("fluid_field").getParameter< uint_t >     ("ghostLayers") }
      , baseFolder        { b.getOneBlock("fluid_field").getParameter< std::string >("baseFolder") }
      , samplingResolution{ b.getOneBlock("fluid_field").getParameter< real_t >     ("samplingResolution") }
      , denseEnabled      { static_cast< bool >(b.getBlock("fluid_field_dense")) }
      , denseFromStep     { denseEnabled ? b.getOneBlock("fluid_field_dense").getParameter< uint_t >("denseFromStep") : uint_t(0) }
      , denseToStep       { denseEnabled ? b.getOneBlock("fluid_field_dense").getParameter< uint_t >("denseToStep") : uint_t(0) }
      , snapshotsEnabled  { static_cast< bool >(b.getBlock("fluid_field_snapshots")) }
      , snapshotSteps     { snapshotsEnabled
                                ? parseUintList(b.getOneBlock("fluid_field_snapshots").getParameter< std::string >("steps"))
                                : std::set< uint_t >{} }
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
      // avgParticleDiameterOverride and dxSIOverride both > 0: skip inFileName's header (see ParticleParams);
      // otherwise use the bed-file-derived scaling.
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
      const real_t particleCollisionTimeNonDim = particle.collisionTimeFactor;
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
      // Surface tension override: >= 0 replaces the value above (0 = no surface tension), < 0 keeps it.
      if (physics.surfaceTensionOverride >= real_c(0)) { surfaceTension = physics.surfaceTensionOverride; }

      const real_t viscosity_SI = real_c(1.0016e-6);  // kinematic viscosity of water at 20°C at 1 bar
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
   std::vector< ContainerCylinderParams > containerCylinders; // DropTest round container
   std::vector< FloatingSphereParams >  floatingSpheres;
   std::vector< FloatingBoxParams >     floatingBoxes;
   std::vector< FloatingCylinderParams > floatingCylinders;
   std::vector< FloatingCapsuleParams >  floatingCapsules;
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
      , containerCylinders { parseContainerCylinders(cfg) }
      , floatingSpheres  { parseFloatingSpheres(cfg) }
      , floatingBoxes    { parseFloatingBoxes(cfg) }
      , floatingCylinders{ parseFloatingCylinders(cfg) }
      , floatingCapsules { parseFloatingCapsules(cfg) }
      , model        { cfg->getOneBlock("ModelParameters") }
      , evaluation   { cfg->getOneBlock("EvaluationParameters") }
      , vtk          { cfg->getOneBlock("VTK") }
      , computed     { particle, domain, pid, physics }
   {}
};


} // namespace antidunes
} // namespace walberla