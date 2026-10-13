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

   // [Claude, 2026-09-24] Diagnostic-only: when true, EVERY step from keepAllSnapshotsFromStep onward
   // is checkpointed into its own step_<N> subfolder instead of the normal single rolling checkpoint
   // (which overwrites itself each time). Used to capture a fine-grained per-step history leading up
   // to a reproducible crash for post-mortem cell-by-cell inspection. Default false/0 reproduces every
   // existing prm's exact prior behavior.
   bool              keepAllSnapshots;
   uint_t            keepAllSnapshotsFromStep;

   // [Claude, 2026-09-27] Per explicit request ("save always the last two checkpoint not only the
   // last one"): when true, the regular-cadence checkpoint (every snapshotFrequency steps) alternates
   // between two subfolders (snapshotBaseFolder/gen_0, gen_1) instead of overwriting one single
   // rolling checkpoint in place -- so at any point in time BOTH the most recent AND the
   // second-most-recent checkpoint survive on disk (whichever one wasn't just written), giving one
   // extra step of recovery margin if a write is interrupted or a crash happens to land exactly on a
   // checkpoint step. Mutually exclusive in practice with keepAllSnapshots (that mode is for a
   // different, per-step diagnostic use case); if both are set, keepAllSnapshots's own `else if`
   // branch (Antidunes.cpp) never fires when the regular cadence also matches, same as before this
   // change. To resume from whichever generation is newest, compare gen_0/simState.file vs
   // gen_1/simState.file's beginTimeStep and point a resume prm's own snapshotBaseFolder at that
   // specific gen_<N> subfolder directly -- loading itself is unchanged, only writing rotates.
   // Default false reproduces every existing prm's exact prior (single rolling checkpoint) behavior.
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

   // false disables the PID controller's own forcing update entirely (Forcing.h's ForcingAdjuster
   // then leaves `acceleration[0]` pinned at its setupSimState-computed initial value for the whole
   // run -- a "fixed hydraulic gradient" rather than a controlled flume; see THESIS_PLAN.md item 3).
   // Default true reproduces every existing prm's exact prior behavior (unconditional PID control).
   bool   enabled;

   // Optional knob for a PID-off (enabled=false) run: setupSimState's default fixed accelerationX
   // is the plain-channel formula 3*U*viscosity/liquidHeight^2, which assumes unobstructed flow and
   // has no knowledge of any BodyBC obstacle in the domain. For a scenario with a real flow-blocking
   // pillar, that formula under-drives the flow (measured ~8x too low against this project's own
   // pillar geometry -- see THESIS_PLAN.md, 2026-09-18). When > 0, this overrides accelerationX
   // directly with a caller-supplied value (e.g. the converged actuatingVariable read back from a
   // matching PID-on run's own eval database, once it has settled) instead of the analytical formula.
   // 0 = off (default; reproduces the pre-existing analytical-formula behavior exactly).
   real_t accelerationXOverride;

   // Optional mid-run switch from PID control to a constant force (ForcingAdjuster in Forcing.h).
   // switchOffStep > 0: PID control runs normally until that absolute timestep, then acceleration_x is
   // frozen at the MEAN of the PID's own commanded values over the last switchOffAveragingWindow steps
   // and never adjusted again. Lets a flow be spun up self-consistently under PID control (correct force
   // for the target velocity, no hand-set value or startup transient), then switched to constant forcing
   // right before debris arrives. The mean (not the instantaneous value) is frozen because the PID's
   // command oscillates around its converged value; freezing an arbitrary phase would bias the force by
   // the oscillation amplitude. 0 = off (default: PID control, or none if enabled=false, for the whole
   // run). Ignored when enabled=false. Window 0 freezes at the last commanded value instead of a mean.
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
   // [Claude, 2026-09-21] surface tension override in lattice units: >= 0 replaces the derived/hardcoded value (0 = NO surface tension); < 0 (DEFAULT -1) = keep the value from the units above
   real_t surfaceTensionOverride;

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

   // 0 (DEFAULT, off) -> checkPdfNaN behaves exactly as before (isfinite only). > 0 -> checkPdfNaN
   // (SetupStages.h) ALSO logs a throttled (first 20 per rank) [LOW-RHO-PRECURSOR] warning for any
   // liquid/interface cell whose aggregate density is still finite but has drained below this value --
   // catching the step(s) BEFORE a cell goes fully non-finite, which the plain isfinite scan can never
   // see by construction. Same intended usage as denseNanScanEnabled above (a targeted, windowed
   // rerun from a checkpoint near an already-known crash, not routine production use) -- gated on the
   // SAME nanScanNow cadence in Antidunes.cpp, no separate window needed. Pick a threshold well below
   // the fluid's nominal rho (~1 in lattice units); e.g. 0.1.
   real_t lowRhoPrecursorThreshold;

   // 0 (DEFAULT, off) -> no effect. > 0 -> checkPdfNaN ALSO logs a throttled (first 20 per rank)
   // [HIGH-VELOCITY-PRECURSOR] warning for any liquid/interface cell whose velocity magnitude is
   // still finite but exceeds this value -- the symmetric counterpart to lowRhoPrecursorThreshold
   // above: a blowup does not always drain toward rho=0 first, it can equally show up as an
   // unrealistically large (but not yet Inf/NaN) velocity for a step or two beforehand. Same intended
   // usage (targeted windowed rerun, not routine production use). Pick a value well above the flow's
   // own target velocity (e.g. 10x targetMeanVelocityMagnitude).
   real_t highVelocityPrecursorThreshold;

   // 0 (DEFAULT, off) -> no effect. > 0 -> checkPdfNaN ALSO logs a throttled (first 40 per rank) [POP-PRECURSOR] warning for any liquid/interface cell in which some population |f_i| exceeds this
   // value (physical populations are <~ 0.3 in lattice units): the earliest visible sign of the exploding-pore mode (zero-mass f_i ~ prod(3c^2-2) pattern, THESIS_PLAN.md class (B)), typically several
   // steps before |u| becomes large. The log line also reports the density, the fill level and the cell's OBSTACLE GEOMETRY (number of the 13 opposite D3Q27 direction pairs with obstacle on both sides
   // = "pinched" pairs, and the number of obstacle neighbours of 26), so the pinched-cell hypothesis can be checked directly. Same cadence/window as the other precursors (denseNanScan*).
   real_t populationPrecursorThreshold;

   // Density WINDOW watcher (0 = that side off, DEFAULT off on both sides). > 0 -> checkPdfNaN (SetupStages.h) checks the density rho = sum_i f_i (the 0th moment, exactly what PdfField::getDensity returns) of every
   // liquid/interface cell (liquid cells only if densityPrecursorLiquidOnly) and reports cells with rho < densityPrecursorLow or rho > densityPrecursorHigh (physical rho ~ 1; typical values 0.1 and 1.9):
   //  * per rank and scanned step ONE [DENSITY-PRECURSOR t=.. after=..] line (first 60 per rank): number of cells below/above the window and, for the FIRST offending cell in scan order, the LOWEST-rho and the
   //    HIGHEST-rho cell: cell class (LIQ/INTF/fMO), world position, rho, fill, max|f|, pinched opposite obstacle pairs (of 13), obstacle neighbours (of 26), former-moving-obstacle / interface / gas neighbour
   //    counts and the uids of the rods owning the obstacle neighbours;
   //  * on the root ONE [DENSITY-WINDOW t=.. after=..] summary line per scanned step with the global counts (first 400 lines).
   // Same cadence/window as the other precursors (denseNanScanEnabled / denseNanScanFromStep / denseNanScanToStep). Benign hits exist: drained interface cells at the free surface have rho -> 0
   // (see SimpleBB's dry-link comment) -- use densityPrecursorLiquidOnly=true to look only at the bulk fluid.
   real_t densityPrecursorLow;
   real_t densityPrecursorHigh;
   bool   densityPrecursorLiquidOnly;

   // FLUID SANITY SWEEP (DEFAULT off): an extra per-step sweep at the END of every main step (after the evaluators, see Antidunes.cpp) that checks every liquid/interface cell of every rank against configurable
   // ranges and reports (or aborts). Independent of the dense NaN scan above, and with an additional fill-level check. Per checked cell: rho = sum_i f_i, max_i |f_i|, |u| (PdfField::getVelocity), fill.
   //   fluidSanityMaxVelocity default 0.17 = lattice Mach number 0.3 (Ma = |u|/cs, cs = 1/sqrt(3)); Kruger et al. 2017 Sect. 7.2 p. 279: Ma < 0.3 is 'small', above it compressibility errors grow (Bogner 2017 Eq. 2.11 similar).
   //   It is an accuracy warning, NOT a collapse threshold: cells at 0.21-0.25 survived in the surviving runs. Earlier sanity prms (4251480..4252463) explicitly set 0.3 (Ma 0.52).
   //   violation categories: non-finite population/rho/u; rho < fluidSanityDensityMin; rho > fluidSanityDensityMax; |u| > fluidSanityMaxVelocity; max|f_i| > fluidSanityMaxPopulation;
   //                         fill outside [-tol, 1+tol] for interface cells or outside [1-tol, 1+tol] for liquid cells (tol = fluidSanityFillTolerance).
   // Output: per rank and checked step ONE [FLUID-SANITY t=..] line if it has violations (first offending cell in scan order, lowest-rho, highest-|u|, highest-|f| cell with class, world position, rho, fill, max|f|,
   // |u|, obstacle geometry, neighbour classes, owners of the obstacle neighbours; first 60 lines per rank) and ONE root line [FLUID-SANITY-SUMMARY t=..] with the global per-category counts (first 400 lines).
   // fluidSanityAbortAbove > 0 -> the run aborts (WALBERLA_ABORT) as soon as the global number of violating cells in one step exceeds this value (0 = never abort, report only).
   // A category is switched off by setting its threshold to 0 (density min/max, velocity, population) or the fill tolerance to a negative value.
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
   // additional opt-in checks (0 / false = off): |excess mass| of a checked cell; largest |u| resp. |rho| difference to a checked face neighbour of the same block (a shock / single-cell spike detector that triggers long before
   // the absolute limits); liquid cell touching a gas cell (broken interface layer); and fluidSanityStatsEveryNSteps > 0 additionally prints [FLUID-STATS t=..] every N steps (global rho range, max|u|, max|f_i|, fill range and the
   // fastest cell) regardless of any violation, so the growth of the extrema can be read off before a threshold is crossed.
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

   // Non-dimensional DEM collision time: tc = collisionTimeFactor * avgParticleDiameter [LBM steps] (default 4, the previously hard-coded value). The linear spring-dashpot stiffness is
   // k = m_eff (pi^2 + ln^2(e)) / tc^2 (LinearSpringDashpot::setStiffnessAndDamping), so k ~ 1/collisionTimeFactor^2 and is INDEPENDENT of numSubCycles (more sub-cycles only resolve the same
   // spring better). A smaller factor = stiffer contacts = smaller steady-state overlap under a given load (overlap = F/k); keep tc >~ 20-40 DEM steps (tc*numSubCycles) for explicit-integration stability.
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

// [Claude, 2026-10-03] A vertical round container (DropTest only): a CylindricalBoundary in containment mode with its
// axis along z through `center` (only x, y used), infinite in z -- the bodies stay inside radius `radius`. Added to
// reproduce the round container of Kodam et al. (2010b, Sec. 3.3: 50.6 mm inner diameter) instead of an area-matched
// square box. Optional prm block: ContainerCylinders { ContainerCylinder { center < x, y, z >; radius r; } }
struct ContainerCylinderParams {
   Vector3< real_t > center;
   real_t            radius;

   ContainerCylinderParams(const config::Config::BlockHandle& b)
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
   // [Claude, 2026-10-01] axis of the initial rotation (default world Y, i.e. unchanged for existing prms); with a random
   // axis and angle a box can be released in an arbitrary orientation (randomised drop tests).
   Vector3< real_t > initialTiltAxis;
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

// [Claude, 2026-09-23] "Pill" shape (mesa_pd::data::Capsule, SHAPE_TYPE 7): a cylindrical core with two
// hemispherical caps. `length` has the SAME "full end-to-end extent including the caps" meaning as
// FloatingCylinderParams::length, so a scenario can be ported between the two shapes at identical
// radius/length numbers. No closed-form contact algorithm exists for this shape anywhere in the
// project (deliberately, for now -- see Capsule::support()'s own doc comment): every capsule pair
// (capsule-capsule, capsule-anything) goes through GeneralContactDetection's generic GJK/EPA path,
// the same way Cylinder x Cuboid already does, with no dispatch code of its own required.
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
   bool              startInactive;       // see FloatingSphereParams::startInactive for the full comment.

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
   // true  → Bogner & Rüde (2013) per-cell gas/liquid/interface state selection for FormerMO cells.
   // false → legacy uniform reconstruction: every uncovered cell becomes an interface cell
   //         (no fill=1 over-promotion band; the all-interface path that never drifted/NaN'd).
   bool        useBognerReconstruction;

   // true  → φ(x_b) fill-weighting of the MovingObstacle (SimpleBB) MEM force (Bogner & Rüde 2013,
   //         Eq. 5.9): interface cells at a floater's waterline contribute proportionally less force,
   //         suppressing the phantom full-liquid buoyant lift.
   // false → pass a null fill accessor so the raw two-phase MEM force is used (reference / ablation).
   bool        useFillWeighting;

   // [Claude, 2026-09-21] true  → INTERPOLATED bounce-back (CLI, Ginzburg et al. 2008 = lbm_mesapd_coupling::CurvedLinear) on the MovingObstacle links whose further fluid cell is liquid/interface
   //         (exact wall distance from a ray-particle intersection; second-order wall placement for curved/tilted bodies); other links fall back to simple bounce-back (see SimpleBB.h).
   // false → (DEFAULT) plain simple bounce-back = the published Antidunes/Bogner behaviour (wall halfway on every cut link, stair-case surface).
   bool        useInterpolatedBounceBack;
   // Experimental (SealedCellGuard.h): fluid cells with at least this many of their 26 D3Q27 neighbours covered by a
   // moving obstacle / no-slip wall (and >= 1 moving obstacle) are turned into obstacle cells right after the particle
   // mapping ("close sub-cell channels between touching bodies"). 0 = off (default), typical 24.
   uint_t      closeSealedCellsMinNeighbours;

   // Experimental (MovingParticleMapping.h, 'cellCenterTestOffset'): the cell-centre containment test is done at the cell
   // centre plus mappingTieOffset*(dx, 1.3 dy, 1.7 dz) so that a body face plane through cell centres (axis-aligned box on
   // a half-integer position, 'tie') is mapped the same way every step. 0 = off (default, bit-identical to before); 1e-4 recommended.
   real_t      mappingTieOffset;

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

   // Streamwise/lateral drive of FLOATERS (types 1/3/4) in AddExternalForces. The periodic channel is
   // driven by a body force a_x, a_y; the MEM sees no streamwise pressure gradient, so the particle
   // needs this explicit term.
   // false → (DEFAULT) unchanged: F_xy = vFull * a_xy (fluid density 1 times full volume; the
   //         original Antidunes term, derived for fully submerged sediment).
   // true  → F_xy = rho_p * vFull * a_xy = m_p * a_xy ("tilted gravity" / slope reading of the body
   //         force, Schwarzmeier et al. 2023): every body is pushed with its own mass, independent
   //         of how much of it is submerged. Non-floaters are not affected.
   bool        driveFloatersByMass;

   // Only consulted when useMemOnlyBuoyancy=false (computeSubmergedVolumes()/computeSubmergedVolumesBox()
   // are only ever called from that branch, doDEMSubCycles) -- meaningless and never read otherwise.
   // true  → force vSub = vFull (and buoyancyCenterOffset_ = 0, no righting lever) for every floater
   //         unconditionally, skipping the entire local-free-surface plane-fit/clip/sampling machinery
   //         in computeSubmergedVolumes()/computeSubmergedVolumesBox() entirely. This is "the original,
   //         simple" full-volume-submerged buoyancy computation (correct as long as the floater never
   //         approaches the free surface) -- added for thesis item 1b's MEM-vs-analytical comparison,
   //         which is scoped to the fully-submerged case only and must never depend on incidental
   //         floater-to-surface proximity engaging the (separately validated, partial-submersion-only)
   //         plane-fit machinery mid-run.
   // false → (DEFAULT) unchanged existing behaviour: plane-fit/clip/sampling as before.
   bool        useFullVolumeAnalyticBuoyancy;

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

   // [Claude, 2026-10-01] true  -> a freshly uncovered cell that becomes INTERFACE gets the mean volume fraction of ALL its
   //         resolved neighbours (gas 0, liquid 1, interface its fill) as fill level (PdfReconstructionManager.h,
   //         seedInterfaceFillFromNeighbourhood). false -> legacy: mean fill of the INTERFACE neighbours only, 0.5 if there
   //         is none (default). The legacy rule seeds 0.5 into cells that are interface only because a same-pass
   //         neighbour was still unclassified -- half-empty cells inside the water (branch C, ~1300 per 100 steps in the
   //         hopper runs) or half a cell of water in the air (A'). Changes only the fill, never the cell state.
   bool        seedFillFromAllNeighbours;

   // [Claude, 2026-10-08] true (default) -> the seeding above ignores neighbours flagged obstacle / formerObstacle, as since
   //         2026-10-01; false -> legacy rule as before 2026-10-01 (an interface-flagged neighbour counts even if it also
   //         carries such a flag). See seedInterfaceFillFromNeighbourhood; false avoids the NaN of the free r075 cylinder.
   bool        seedSkipObstacleNeighbours;

   // [Claude, 2026-10-08] true -> freshly uncovered cells follow Bogner's rule literally (Bogner 2017 pp. 113-114): B(x) =
   //         non-obstacle neighbours, with cells freed in the same step counted as obstacle; no liquid -> gas, else no gas ->
   //         liquid, else interface (fill = mean of the interface cells in B(x), 0.5 if none). No same-step interface fallback.
   //         See PdfReconstructionManager.h, classifyStrictBognerNeighbourhood. false (default) -> the fallback rule.
   bool        strictBognerClassification;
   // [Claude, 2026-10-10] true -> (with strictBognerClassification) a freed cell decided liquid that touches gas after the step's
   //         decisions are applied becomes interface with fill 1 (PdfReconstructionManager.h Pass 2). Fixes the waterline NaN of the
   //         r025 cuboid (THESIS_PLAN.md 2026-10-10). false (default) -> counted only ([FORMERMO-STRICT]).
   bool        strictBognerNoLiquidAtGas;
   // [Claude, 2026-10-09] true -> log [MOMENTUM-BUDGET] every evaluationFrequency steps: momentum the freed cells carry into the
   //         fluid, momentum of the fluid in newly covered cells, and the total fluid momentum (liquid + fill-weighted interface).
   //         Diagnostic only. false (default) -> off.
   bool        logMomentumBudget;

   // [Claude, 2026-10-02] true -> finite-difference curvature (CurvatureSweepFiniteDifferences) skips neighbour directions
   //         whose flag is in no group (upstream reads an UNINITIALISED neighbour normal there) or whose normal is not
   //         finite, sets curvature 0 if no direction remains, and logs each case as [CURVATURE-GUARD] (cell, neighbour
   //         flag, neighbour normal). false -> unchanged upstream code (default). Identical results for valid inputs.
   //         Motivation: run 4284164 (192 rods) broke down at step 9046 with curvature=nan in healthy interface cells
   //         (fill 0.994, rho 1.005, finite normal, no obstacle neighbour) on both sides of the block boundary x = 448.
   bool        guardCurvatureInputs;
   // [Claude, 2026-10-03] true -> after the flag initialisation, the outer ghost layer (layer 2) beyond every non-periodic
   //         domain border becomes wall as well (same NoSlip/FreeSlip flag as the layer-1 wall cell next to it, fill 1),
   //         i.e. the wall is two cells thick. false -> unchanged (default). Reason: a Border wall with walldistance -1
   //         occupies only layer 1; layer 2 behind it is set to liquid|interface by fillWithDomain(nrOfGhostLayers)
   //         (mapParticlesOntoCells), never reset (initFlagsFromFillLevels: interior only) and never communicated, and
   //         ObstacleFillLevelSweep (ghost layer 1 included) averages its fill into the wall cells' smoothed fill
   //         ([GHOST2-FLUID] check, BuildTimeloop.h: 11520 of 12000 such cells flagged fluid in a small side-wall case).
   //         The NaN normals of the hopper runs start in the smoothed fill of exactly these side-wall cells (job 4292672).
   bool        wallOuterGhostLayer;

   // [Claude, 2026-10-02] true -> after the normal sweep, scan for non-finite interface normals, smoothed fill levels and
   //         obstacle normals and log [NAN-NORMAL] / [NAN-SMOOTHFILL] / [NAN-OBSTACLENORMAL] with the inputs of the normal
   //         computation, plus a global count every 1000 steps (BuildTimeloop.h, diagnoseNaNNormals). Diagnostic only, no
   //         effect on the simulation; a full field scan per step. Default false.
   bool        diagnoseNaNNormals;

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

   // [Claude, 2026-09-30] true -> (only with useAnalyticalContactDetection=true) Cylinder x Cylinder contacts are
   //         refined by minimising the projection overlap r_i(L) + r_j(L) - L.(x_j - x_i) over directions L
   //         (paper_analytical::refineCylinderContact): exact closed-form scenarios unchanged, inexact ones get the
   //         minimising depth/normal, "no scenario" pairs are decided by the same descent. No GJK/EPA call at all.
   //         Validated in PaperAnalyticalContactTest.cpp (testRefinedCylinderContactVsGjk and the independent
   //         brute-force reference): no missed contacts, depth error <= 5e-4; separated pairs closer than ~0.01
   //         can be reported as grazing contacts (overlap <= 0.006).
   // false -> (DEFAULT) unchanged behaviour (GJK/EPA fallback for "no scenario" and special contacts).
   bool        refineCylinderContacts;

   // [Claude, 2026-10-10] true -> cylinderCylinderContact (PaperAnalyticalContacts.h) passes a band-edge result to
   //         refineCylinderContact when the projection overlap D(n) along its normal exceeds
   //         kBandEdgeConsistencyFactor (3) times its depth (see paper_analytical::kodamBandEdgeAxialCheck). Fixes the
   //         stacked equal-radius cylinders of the round-container drop test (THESIS_PLAN.md 2026-10-10, "FIX B REVISED").
   // false -> (DEFAULT) unchanged closed-form dispatch.
   bool        kodamBandEdgeAxialCheck;

   // true  → (only meaningful together with useAnalyticalContactDetection=true) cylinderCylinderResultIsUncertain's
   //         Tier 1 check (fired when the analytical dispatch found no scenario at all) uses the tight,
   //         capsule-based bound (detail::closestDistanceSqBetweenSegments: the two cylinders' own axis
   //         SEGMENTS vs. the sum of their radii) instead of the original loose center-distance bound
   //         (radius+halfLength sum per side). Both are fully conservative -- neither can ever wrongly
   //         declare "definitely separated" for a pair that's actually touching (validated directly
   //         against GJK/EPA in PaperAnalyticalContactTest.cpp's testTightTier1Accuracy10k, 0 unsafe
   //         cases across 10000 trials) -- but the tight bound is strictly tighter for any cylinder
   //         with length > 0, so it resolves far more dispatch-gap pairs as "definitely separated,
   //         skip GJK/EPA entirely" instead of falling back on every one of them. See
   //         PaperAnalyticalContacts.h's own comment on cylinderCylinderResultIsUncertain's
   //         useTightTier1 parameter for the full derivation.
   // false → (DEFAULT) the original loose center-distance Tier 1 bound, unchanged behaviour.
   bool        useTightTier1;

   // true  → the "cheap raw-center-distance reject" that runs BEFORE any real contact detection at
   //         all (both Cylinder x Cylinder and, newly, Box x Box) uses a tight, ORIENTATION-AWARE bound
   //         instead of collapsing each body to its circumscribing sphere. Cylinder x Cylinder: the
   //         same capsule/segment-distance bound as useTightTier1 above, applied one level earlier (so
   //         a pair it resolves never even reaches the analytical dispatch, not just skips GJK/EPA).
   //         Box x Box: a NEW cheap single-candidate-axis test (detail::boxesDefinitelySeparatedAlongCenterLine
   //         -- projects each box onto the line between the two centres and compares to the sum of
   //         radii, i.e. one axis of full SAT instead of all 15) run ahead of boxBoxSAT, which
   //         previously had no pre-filter at all. Both remain fully conservative (validated directly
   //         against GJK/EPA: 0 unsafe cases across 10000 trials each, see
   //         PaperAnalyticalContactTest.cpp's testTightTier1Accuracy10k and testTightLayer1BoxAccuracy10k)
   //         -- see PaperAnalyticalContacts.h's own comments on closestDistanceSqBetweenSegments and
   //         boxesDefinitelySeparatedAlongCenterLine for the full derivation of each.
   // false → the original circumscribing-sphere-sum reject for Cylinder x Cylinder, unchanged; NO
   //         pre-filter at all for Box x Box (boxBoxSAT runs directly), unchanged.
   // DEFAULT true since 2026-09-20 (user request; a prm can still set it false explicitly).
   bool        useTightLayer1;

   // [Claude, 2026-09-26] EXPERIMENTAL: construct cylinder floaters (type 4, "FloatingCylinder" prm
   // entries) as a Capsule (mesa_pd::data::Capsule, SHAPE_TYPE 7 -- a cylinder with hemispherical end
   // caps, i.e. a "pill" shape / spherocylinder) instead of a flat-capped Cylinder. Tries to sidestep
   // the flat-cap/rim contact-resolution edge cases (e.g. the "band-band-skewed" scenario) that this
   // session traced as the trigger for a 90% rod-rod overlap -> MEM-GARBAGE-LINK -> NaN cascade
   // (THESIS_PLAN.md, job 4267519). No analytic Capsule-Capsule contact routine exists yet, so every
   // capsule pair falls through to the generic GJK/EPA path (mesa_pd::collision_detection::
   // GeneralContactDetection) regardless of useAnalyticalContactDetection -- slower per-pair, but a
   // continuous closest-point-between-segments calculation with no flat-cap/rim special cases.
   // Geometry (radius/length/axis/position/densityRatio/settlingSteps) is UNCHANGED -- only the
   // mesa_pd shape constructed in SetupStages.h differs. DEFAULT false (existing scenarios unaffected).
   bool        useCapsuleShapeForCylinders;

   // DIAGNOSTIC ARTIFICIAL CAPS for cylinder floaters (rods, particle type 4) -- all 0 = OFF (DEFAULT, no effect). They deliberately DAMP everything to answer "are excessive rod motion/force and launches
   // the cause of the crashes, or only their symptom?" (THESIS_PLAN.md). They change the physics and are NOT for production runs. Every cap logs its first hits per rank ([ROD-CAP ...]) and the
   // cumulative counts + the largest value seen before capping appear in the [RodCaps] line of the periodic/end-of-run [ContactDetectionStats] logging.
   //  rodCapHydroForce/rodCapHydroTorque: |F_hydro| / |T_hydro| of the momentum-exchange force applied to the rod are scaled down to this value (also non-finite values are zeroed when a cap is on);
   //  rodCapVelocity/rodCapAngularVelocity: |v| / |omega| of the rod are scaled down to this value after every DEM sub-cycle's velocity update;
   //  rodCapOverlap: the penetration depth fed into the Cylinder contact force law is limited to this value (LU).
   real_t      rodCapHydroForce;
   real_t      rodCapHydroTorque;
   real_t      rodCapVelocity;
   real_t      rodCapAngularVelocity;
   real_t      rodCapOverlap;

   // [Claude, 2026-09-21] TEST OPTION (rodPrescribedVelocityX/Y/Z, default all 0 = off): PRESCRIBED rigid motion of all rods (type 4): after every DEM sub-cycle's velocity update the rods' linear velocity is set to this constant vector and their angular
   // velocity to 0 (forces/torques are ignored, positions integrate with this velocity). Used by the local 'rod bundle channel' test (three touching parallel rods dragged through still water at ~0.012 like in the
   // PID-off crash 4252482) so that lattice cells are covered/uncovered by the moving bundle without any DEM/contact dynamics. Changes the physics -- NOT for production.
   real_t      rodPrescribedVelocityX;   // x/y/z components (three scalars: the config proxy cannot brace-initialise a Vector3 member)
   real_t      rodPrescribedVelocityY;
   real_t      rodPrescribedVelocityZ;

   // true  → SHADOW MODE for Cylinder x Cylinder, meaningful only with useAnalyticalContactDetection=false (GJK/EPA is what is USED): every rod pair for which GJK/EPA or the closed-form
   //         Kodam+Guo dispatch reports a contact ALSO gets the closed-form result computed (never used for forces): which scenario it picked, its overlap/normal, what GJK/EPA reported
   //         instead, and whether the "needed distance along the reported normal" self-consistency check (see CylinderClosedFormRangeTest.cpp) would flag it. Cumulative per-scenario
   //         statistics go into the [ShadowAnalytic] log line (at the end of the run and at every periodic [ContactDetectionStats]); a per-rank CSV (shadow-analytic/rank_<N>.csv, one row per
   //         pair per sub-cycle, both ranks' duplicate computations included) is written for main steps in [shadowAnalyticStartStep, shadowAnalyticEndStep] (EndStep 0 = no CSV).
   // false → (DEFAULT) nothing is computed or logged.
   bool        shadowAnalyticCylinderLog;
   uint_t      shadowAnalyticStartStep;
   uint_t      shadowAnalyticEndStep;

   // true  → (only correct together with useTightLayer1=true -- see below) skips
   //         cylinderCylinderResultIsUncertain's Tier 1 check ENTIRELY for dispatch-gap pairs (no
   //         distance computation of any kind; "no contact" is trusted immediately) instead of merely
   //         swapping which bound it uses. Justified by a proof, not just a benchmark: any pair that
   //         survives the tight Layer 1 pre-filter is GUARANTEED (via the triangle inequality) to also
   //         survive Tier 1's own bound, tight or original -- so once Layer 1 is tight, Tier 1's
   //         dispatch-gap check can provably never reject anything, and running it is pure wasted
   //         work. Confirmed empirically too (job 4194491: cylTier1CheaplyProvedSeparated=0 at every
   //         interval with both useTightLayer1 and useTightTier1 on together). See
   //         PaperAnalyticalContacts.h's own comment on cylinderCylinderResultIsUncertain's skipTier1
   //         parameter for the full derivation.
   // false → (DEFAULT) Tier 1 always runs (using whichever bound useTightTier1 selects).
   // WARNING: skipTier1=true with useTightLayer1=false is NOT just a missed optimization -- it is
   //          incorrect, and will silently trust "no contact" for real dispatch-gap contacts the loose
   //          Layer 1 reject lets through. A startup warning fires if this combination is detected.
   bool        skipTier1;

   // true  → Box x Box's boxBoxSAT (only meaningful together with useAnalyticalContactDetection=true)
   //         replaces its normal with a closest-point-between-edge-segments vector whenever the
   //         winning SAT axis is edge-vs-edge, and inflates the reported overlap by
   //         2*boxRoundingRadius -- see PaperAnalyticalContacts.h's own comment on boxBoxSAT's
   //         useRoundedCorners parameter for the full derivation and the corner/edge contact-normal
   //         DISCONTINUITY this fixes (Lu, Third & Müller 2015, Chem. Eng. Sci. 127, Section 2.1.3,
   //         citing Fraige et al. 2008a) -- proven present in this project's own production Box x Box
   //         dispatch by DEMValidationScenarios.h's runCornerNormalDiscontinuity (47.7% spurious
   //         kinetic-energy gain with this OFF, using this project's own real FloatingBox parameters).
   // false → (DEFAULT) exact pre-existing boxBoxSAT behaviour, unchanged.
   bool        useRoundedBoxCorners;
   // Rounding radius used only when useRoundedBoxCorners is true. 0.05 (DEFAULT) is small relative to
   // a typical FloatingBox's own half-extents (this scenario family's own logs use halfExtents as
   // small as <6,1,1.5>) -- enough to smooth the corner/edge transition without materially changing
   // the box's effective size.
   real_t      boxRoundingRadius;

   // true  → Box x Box additionally resolves via paper_analytical::boxBoxMultiPointContacts (the
   //         OTHER literature-documented remedy for the same defect useRoundedBoxCorners targets --
   //         Fraige et al. 2008a/2008b's "multi-point contact concept", see that function's own doc
   //         comment in PaperAnalyticalContacts.h): every one of both boxes' 8 vertices found embedded
   //         in the other box gets its own local normal/overlap and its own (history-free, same
   //         formula as the existing arch/tunnel secondary-contact code) spring-dashpot force, each
   //         scaled by 1/N (N = total active contact points this sub-cycle) so multiple simultaneous
   //         points don't multiply the total force -- rather than one single SAT-chosen deepest
   //         point/axis whose identity can flip discontinuously between sub-cycles. Independent of
   //         useRoundedBoxCorners -- either, both, or neither can be enabled.
   // false → (DEFAULT) exact pre-existing single-point behaviour, unchanged.
   bool        useMultiPointBoxContact;

   // true  → both paper_analytical::boxBoxSAT's single overlap AND every individual vertex overlap
   //         paper_analytical::boxBoxMultiPointContacts reports are clamped to at most
   //         boxOverlapCapFraction * min(halfExtents over BOTH boxes' 3 components) before being used
   //         for force computation -- see boxBoxMultiPointContacts' own doc comment (PaperAnalyticalContacts.h)
   //         for the full rationale. Confirmed necessary by DEMValidationScenarios.h's
   //         multiPointDeepOverlapStability/Fuzz/EdgeCases scenarios: UNCAPPED, useMultiPointBoxContact
   //         can inject dramatically MORE spurious energy than doing nothing at all under deep
   //         interpenetration (e.g. +618% vs baseline's -0.9% on the "grazing-plank-stack" edge case,
   //         the everyday two-slightly-askew-planks geometry a real jam actually produces) -- exactly
   //         the opposite of what the mitigation is for. Independent of useRoundedBoxCorners/
   //         useMultiPointBoxContact -- meaningful with either, both, or (via boxBoxSAT's own cap)
   //         even neither enabled.
   // false → (DEFAULT) exact pre-existing unbounded overlap behaviour, unchanged.
   bool        useBoxOverlapCap;
   // Fraction of min(halfExtents) used as the cap when useBoxOverlapCap is true. 0.3 (DEFAULT) sits
   // comfortably above the shallow, near-degenerate overlaps this contact model was actually tuned
   // against (typically a few percent of the smaller half-extent, per this project's own
   // STRONG-CONTACT-FORCE/DEM-WATCH log evidence) while still well below the geometric maximum
   // (1.0 = a vertex sitting exactly at the other box's centre).
   real_t      boxOverlapCapFraction;

   // true  → Box x Box's boxBoxSAT axis selection becomes "sticky": among all candidate axes within
   //         axisHysteresisMargin (relative) of the true minimum-gap axis, the one whose normal aligns
   //         best with the PREVIOUS sub-cycle's reported normal for this exact uid pair is used, even
   //         when it is not the strict minimum -- see paper_analytical::boxBoxSAT's own doc comment
   //         (PaperAnalyticalContacts.h) for the full derivation. The OTHER candidate fix for the same
   //         axis-flip discontinuity useRoundedBoxCorners targets, without useRoundedBoxCorners' own
   //         deep-overlap side effect (confirmed by this project's own multiPointDeepOverlapStability/
   //         EdgeCases scenarios: useRoundedBoxCorners is +884% on "plank-vs-near-cube-deep" despite
   //         fixing the axis-tie case cleanly). Independent of useRoundedBoxCorners/
   //         useMultiPointBoxContact/useBoxOverlapCap -- any combination is meaningful.
   // false → (DEFAULT) exact pre-existing behaviour, unchanged.
   bool        useAxisHysteresis;
   // Relative margin (DEFAULT 0.15 = 15%) within which a competing axis is considered "as good as" the
   // true minimum for hysteresis purposes -- see useAxisHysteresis's own comment.
   real_t      axisHysteresisMargin;

   // true  → boxBoxMultiPointContacts' own vertex SET is smoothed across sub-cycles into a persistent
   //         manifold (one entry per possible vertex identity, up to 16, per uid pair) instead of being
   //         used fresh every sub-cycle -- a newly-embedded vertex fades its overlap contribution IN
   //         over several sub-cycles, and one that just separated fades OUT the same way, rather than
   //         appearing/vanishing in one discontinuous step. See ParticleKernels.h's
   //         ComputeContacts::boxBoxManifold_ own doc comment for the full rationale -- a simplified,
   //         same-spirit version of the persistent contact manifold + warm-starting real physics
   //         engines (Bullet/PhysX) use, targeting a THIRD distinct symptom from useAxisHysteresis
   //         (which axis wins) and useBoxOverlapCap (how deep one point's force can be): the vertex SET
   //         itself changing abruptly. Only meaningful together with useMultiPointBoxContact=true.
   // false → (DEFAULT) exact pre-existing behaviour, unchanged.
   bool        useContactManifold;
   // Fraction of full strength gained/lost per sub-cycle a manifold entry fades by (DEFAULT 0.25 = a
   // vertex reaches full/zero strength over 4 sub-cycles).
   real_t      manifoldFadeRate;

   // true  → the normal-direction stiffness/damping fed into LinearSpringDashpot::setStiffnessAndDamping
   //         is computed from a GENERALIZED effective mass that includes each body's rotational inertia
   //         at the actual contact point, not just 1/(1/mA+1/mB) (pure translational reduced mass):
   //             1/meff = 1/mA + 1/mB + (rA×n)ᵀ IA_world⁻¹ (rA×n) + (rB×n)ᵀ IB_world⁻¹ (rB×n)
   //         where rA/rB are the contact point's offset from each body's centre of mass and IA_world⁻¹/
   //         IB_world⁻¹ are each body's inverse inertia tensor rotated into the WORLD frame (see
   //         ParticleKernels.h's computeRotationalEffectiveMass, and mesa_pd::kernel::VelocityVerlet's own
   //         use of math::transformMatrixRART for the same body→world inertia transform). This is the
   //         standard generalized/rigid-body-collision effective mass (Baraff, "Physically Based
   //         Modeling: Rigid Body Simulation"), and is EXACT for sphere-sphere contact (there rA/rB are
   //         always parallel to n, so the rotational terms are identically zero and meff reduces to the
   //         plain translational formula) -- consistent with sphere-sphere floaters never showing the
   //         energy-injection defect this whole family of options targets.
   //
   //         Root-cause rationale: LinearSpringDashpot's stiffness/damping are calibrated ONCE per
   //         particle-type pair from the plain translational effectiveMass, so the requested coefficient
   //         of restitution is only actually delivered to the linear closing-velocity DOF. Any contact
   //         whose point has a nonzero moment arm off the line of centres (any Box/Cylinder corner or
   //         edge contact) bleeds kinetic energy into rotation that this damping was never sized to
   //         control -- under-damping that combination and injecting energy on its own, BEFORE and
   //         INDEPENDENT of any SAT axis-flip discontinuity. useAxisHysteresis/useContactManifold/
   //         useRoundedBoxCorners all instead fix WHICH (point,normal,overlap) triplet feeds the force
   //         law; this option is the only one of the five that changes the force law's own calibration,
   //         so it is meaningful together with any/all of the others, or entirely on its own.
   // false → (DEFAULT) exact pre-existing behaviour: plain translational effectiveMass, unchanged.
   bool        useRotationalEffectiveMass;

   // true  → Box x Box contact is resolved as a genuine "rounded box" (spherosimplex): the Minkowski
   //         sum of the SHARP box (its true faces/edges/vertices) with a sphere of radius
   //         spheroBoxRadius. Unlike useRoundedBoxCorners (which patches boxBoxSAT's discrete 15-axis
   //         dispatch, only fixing the normal for the ONE case where the winning axis is edge x edge)
   //         or useMultiPointBoxContact (which resolves the SHARP corner discontinuity via multiple
   //         simultaneous points), this replaces the whole detection path with mesa_pd's own GJK/EPA
   //         (the SAME general-purpose convex-convex algorithm this project already trusts as its
   //         fallback for every other shape pair), called with margin=spheroBoxRadius on each of the
   //         two SHARP boxes (mesa_pd::collision_detection::GJK::doGJKmargin/EPA::doEPAmargin already
   //         implement exactly this "core shape + margin" construction -- the identical mechanism real
   //         physics engines use for collision margins). Because GJK/EPA computes the TRUE closest
   //         points (or true penetration via EPA) between the two sharp skeletons across EVERY feature
   //         pair (vertex-vertex, vertex-edge, vertex-face, edge-edge, edge-face, face-face) through
   //         ONE unified algorithm, the reported normal is continuous through every possible feature
   //         transition, not just the edge-edge one useRoundedBoxCorners special-cases -- the literature
   //         analogue is Galindo-Torres & Pedroso (2010, Phys. Rev. E 81, 061303)'s Voronoi-based
   //         spheropolyhedra, cited in Lu, Third & Müller (2015) Section 2.2.5 as the general-polyhedron
   //         version of this same "minimum distance between skeletons, dilated by spheroradius" idea --
   //         our case is the box-specialised instance, using mesa_pd's own already-tested GJK/EPA
   //         instead of re-deriving a closed-form OBB-OBB distance routine from scratch.
   //         MUTUALLY EXCLUSIVE with useRoundedBoxCorners/useMultiPointBoxContact/useBoxOverlapCap/
   //         useAxisHysteresis/useContactManifold -- when true, this option's detection path is used
   //         INSTEAD of all of those (they only apply to the boxBoxSAT dispatch this bypasses
   //         entirely). Still composable with useRotationalEffectiveMass, which only changes the force
   //         law's calibration, independent of how the (point,normal,overlap) triplet was obtained.
   // false → (DEFAULT) exact pre-existing behaviour, unchanged.
   bool        useSpheroBoxContact;
   // Rounding radius for useSpheroBoxContact, applied identically to both boxes in a pair (combined
   // Minkowski dilation = 2*spheroBoxRadius, matching two equally-rounded bodies). DEFAULT 0.05 matches
   // boxRoundingRadius's own default for the SAT-patch approach, for a like-for-like comparison.
   real_t      spheroBoxRadius;

   // [Claude, 2026-09-30] cylinderContactSkin (DEFAULT 0 = off): contact skin of the cylindrical floaters (type 4),
   // in lattice units. For contact detection only, every floating cylinder is replaced by an inflated cylinder with
   // radius R + m and half-length h + m: the DEM contact force starts when the SURFACE GAP between two rods drops below
   // 2m (face-face, band-band), and between a rod and a wall plane or a pillar when it drops below m. At the rims the
   // sharp inflated edge reaches up to (sqrt 2 - 1) m further than a rounded (Minkowski) skin would. Purpose: keep at least one resolved fluid
   // cell between touching rods, whose line contacts otherwise seal sub-cell liquid channels (THESIS_PLAN.md,
   // 2026-09-21 rod-bundle test: sealed channel cells collapse in density). [Claude, 2026-10-01, checked against the
   // PDFs in cld_kbase] This is the "effective diameter" remedy of Boutt et al. (2007, JGR 112 B10209, Sec. 2.3,
   // Fig. 1): in their 2-D LBM-DEM, densely packed discs leave no connected flow path, so the fluid sees a smaller
   // "hydraulic" diameter (multiplier 0.8) than the mechanical contact diameter -- the same separation of contact
   // size and fluid size, applied there by shrinking the fluid body instead of enlarging the contact body. Rettinger
   // (2023, p. 7) cites it for this purpose and notes that it lowers the predictive power. Related but with a
   // different aim: Glowinski et al. (2001, J. Comput. Phys. 169, Sec. 5) add a short-range repulsion between bodies
   // within a range of about one mesh spacing to prevent them from overlapping. The fluid mapping (MEM) is NOT
   // dilated. Implementation (ParticleKernels.h, SolidCylinderContactDetection): the Cylinder x Cylinder,
   // Cylinder x HalfSpace, and Cylinder x CylindricalBoundary (solid pillar) overloads call themselves once with the
   // inflated shapes, so the unchanged, validated detection (closed forms, GJK/EPA fallback) is used. A GJK/EPA
   // "margin" formulation was tried first and rejected: mesa_pd's EPA gives up when only the margins overlap
   // (EPA.cpp, "numPoints == 1"), which is exactly the parallel-rod line contact the skin is meant for.
   // The rod interaction radius (cheap reject, ghost sync range) is enlarged by 2m (SetupStages.h), which bounds the
   // circumscribing radius of the inflated cylinder, sqrt((R+m)^2 + (h+m)^2) <= sqrt(R^2+h^2) + sqrt(2) m.
   // Not applied to Cylinder x Sphere / Cylinder x Box (no such pairs in the rod scenarios).
   real_t      cylinderContactSkin;

   // [Claude, 2026-10-01] boxContactSkin (DEFAULT 0 = off): box analog of cylinderContactSkin above
   // -- same purpose, same mechanism, same formula, just for axis-aligned-at-rest FloatingBox
   // floaters (type 3) instead of cylinders. Every floating box is replaced, for contact detection
   // only, by an inflated box with each half-extent + m: the DEM contact force starts when the
   // SURFACE GAP between two boxes drops below 2m (face-face), and between a box and a wall plane or
   // a pillar when it drops below m. Implementation (ParticleKernels.h, SolidCylinderContactDetection):
   // the Box x Box, Box x HalfSpace, and Box x CylindricalBoundary (solid pillar) overloads call
   // themselves/a shared core once with the inflated box(es), so the unchanged, validated detection
   // (SAT, sphero-box, or GJK/EPA fallback, whichever this call would otherwise have taken) is used.
   // The box interaction radius (cheap reject, ghost sync range) is enlarged by 2m (SetupStages.h),
   // the same simple safe over-approximation used for cylinderContactSkin's own circumR. Not applied
   // to Box x Sphere / Box x Cylinder (no such pairs in the box-hopper scenarios this was built for).
   real_t      boxContactSkin;

   // Only meaningful when useSpheroBoxContact is true -- selects WHICH algorithm computes the
   // spherosimplex (point,normal,overlap) triplet:
   // false → (DEFAULT) mesa_pd's own GJK/EPA run with margin=spheroBoxRadius on the two sharp
   //         boxes (paper_analytical::boxBoxSAT's sibling comment describes the mechanism).
   // true  → a closed-form, non-iterative analytical distance between the two sharp boxes
   //         instead (paper_analytical::boxBoxSpheroAnalytical / eberly_box_distance namespace,
   //         PaperAnalyticalContacts.h) -- ported from David Eberly's GeometricTools library
   //         (github.com/davideberly/GeometricTools), whose line-vs-box case dispatch is
   //         documented there as Section 10.9.4 "Linear Component to Oriented Bounding Box" of
   //         Schneider, P.J. & Eberly, D.H. (2002) "Geometric Tools for Computer Graphics",
   //         Morgan Kaufmann. Unlike GJK/EPA (an iterative simplex refinement with its own
   //         convergence tolerance), this is a deterministic finite case dispatch -- added
   //         specifically to A/B-test against the GJK/EPA path on the exact same test batteries.
   bool        useSpheroBoxAnalytical;

   // true  → Box x Box's multi-point contact list (see useMultiPointBoxContact above) is computed by
   //         paper_analytical::boxBoxFraigeMultiPointContacts INSTEAD OF boxBoxMultiPointContacts --
   //         a fuller, more literal replica of Fraige, Langston & Chen (2008, Powder Technology
   //         186(3):224-240, Section 3.3)'s own method, adding the two contact modes
   //         boxBoxMultiPointContacts omits (their Fig. 5 edge-edge and edge-face cases -- two boxes
   //         overlapping with NEITHER having a vertex embedded in the other, e.g. crossing like an
   //         "X") plus their own near-corner normal blend (Fig. 4b). Takes precedence over
   //         useMultiPointBoxContact if both are set (see ComputeContacts' own Box x Box operator()
   //         for the exact dispatch) -- kept as a fully separate flag/function pair specifically so
   //         enabling it can never change useMultiPointBoxContact's own already-validated behaviour.
   //         Composable with useContactManifold/useBoxOverlapCap/useRotationalEffectiveMass exactly
   //         like useMultiPointBoxContact is.
   // false → (DEFAULT) exact pre-existing behaviour, unchanged.
   bool        useFraigeBoxContactModel;

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

   // 0 (DEFAULT) -> logFloaterHistory writes through to the end of the run, as before. Nonzero -> both
   // log() and logState() skip every row once mainStep exceeds this value -- the upper-bound
   // counterpart to logFloaterHistoryStartStep above, so a targeted rerun can scope full-resolution
   // logging to a closed [start, end] window (e.g. "90000 to 100000") instead of paying the per-
   // subcycle cost from the start step through to whenever the run happens to end/crash.
   uint_t      logFloaterHistoryEndStep;

   // [Claude, 2026-09-27] Per explicit request ("create a box around all the cylinders in the crash
   // area, log out their values and see the crash evolution"): true -> every rank writes its own CSV
   // (fluidBoxHistoryBaseFolder/rank_<N>.csv) with ONE ROW PER OBSTACLE-OR-INTERFACE CELL inside the
   // fixed world-space box [fluidBoxMin, fluidBoxMax], once per MAIN step (not sub-cycle -- the fluid
   // field itself only changes once per main step, see CLAUDE.md's per-timestep execution order),
   // within [logFluidBoxHistoryStartStep, logFluidBoxHistoryEndStep]. Columns: mainStep, rank,
   // world x/y/z, flag, owning particleField_ uid, density, fill, max|PDF component|. Complements
   // logFloaterHistory (which only sees the FLOATERS' own state) and the single-cell trackWatchCells
   // probes (BuildTimeloop.h, hardcoded to a handful of specific world coordinates from earlier crash
   // instances) with the FULL spatial evolution of every obstacle/interface cell in a region known to
   // contain a recurring crash's rods -- restricted to OBST/INTF cells (not every cell in the box) to
   // keep row count bounded: a 3-rod-sized box is a few thousand cells, not tens of thousands.
   // OFF by default (zero cost when disabled -- the caller checks this before doing anything).
   bool        logFluidBoxHistory;

   // Base folder for logFluidBoxHistory's per-rank CSVs, sibling to vtk-out/eval/snapshots (created by
   // Antidunes.cpp alongside those, root-only, only when logFluidBoxHistory is true).
   std::string fluidBoxHistoryBaseFolder;

   // World-space corners of the fixed box logFluidBoxHistory scans every main step (inclusive of both
   // corners' whole cells). Default (0,0,0) both -- an empty box, so logFluidBoxHistory silently logs
   // nothing if left unset even when enabled; the caller (BuildTimeloop.h) does not warn about this
   // the way logFloaterHistoryEndStep<=StartStep does, since an empty box is a normal "not configured
   // yet" state, not a likely mistake in an otherwise-populated prm.
   Vector3< real_t > fluidBoxMin;
   Vector3< real_t > fluidBoxMax;

   // 0 (DEFAULT) -> logFluidBoxHistory writes from step 0, as before. Nonzero -> skip every main step
   // before this value -- same windowed-diagnostic convention as logFloaterHistoryStartStep, so a
   // targeted rerun of an already-approximately-located crash doesn't pay the per-step box-scan cost
   // (and doesn't fill the CSV) for steps far from where the crash actually happens.
   uint_t      logFluidBoxHistoryStartStep;

   // 0 (DEFAULT) -> logFluidBoxHistory writes through to the end of the run, as before. Nonzero ->
   // skip every main step once mainStep exceeds this value -- upper-bound counterpart to
   // logFluidBoxHistoryStartStep above.
   uint_t      logFluidBoxHistoryEndStep;

   // [Claude, 2026-09-28] Per explicit request to investigate the sphere/cylinder/cuboid equilibrium
   // position bias (still-water buoyancy sweeps consistently settle HIGHER -- less submerged -- than
   // the naive constant-density Archimedes target, growing with density ratio; see THESIS_PLAN.md):
   // true -> every rank that owns the world column [densityColumnX, densityColumnY] (all z) writes
   // its own CSV (densityColumnBaseFolder/rank_<N>.csv) with ONE ROW PER CELL in that column --
   // EVERY flag type, unlike logFluidBoxHistory's OBST/INTF-only restriction, since the whole point
   // here is the LIQUID density profile far from the floater (to test whether the settled fluid is
   // hydrostatically stratified under the gravity body force, rather than the uniform rho=1 the
   // simple Archimedes target assumes). Once per MAIN step, within [logDensityColumnStartStep,
   // logDensityColumnEndStep]. Columns: mainStep, rank, worldX, worldY, worldZ, flag, density, fill.
   // OFF by default (zero cost when disabled).
   bool        logDensityColumn;

   // Base folder for logDensityColumn's per-rank CSVs, sibling to vtk-out/eval/snapshots (created by
   // Antidunes.cpp alongside those, root-only, only when logDensityColumn is true).
   std::string densityColumnBaseFolder;

   // World-space X/Y of the column logDensityColumn scans every main step (every z, full domain
   // height); pick a point far from any floater's footprint so the logged density reflects the
   // ambient/far-field fluid, not the floater's own local disturbance. Default 0/0 -- a real prm
   // must set both explicitly, same "silently logs nothing if left unset" convention as fluidBoxMin/
   // fluidBoxMax above.
   real_t      densityColumnX;
   real_t      densityColumnY;

   // 0 (DEFAULT) -> logDensityColumn writes from step 0. Nonzero -> skip every main step before this
   // value -- same windowed-diagnostic convention as logFluidBoxHistoryStartStep.
   uint_t      logDensityColumnStartStep;

   // 0 (DEFAULT) -> logDensityColumn writes through to the end of the run. Nonzero -> skip every main
   // step once mainStep exceeds this value.
   uint_t      logDensityColumnEndStep;

   // [Claude, 2026-09-28] logWaterHeight (true -> ON, DEFAULT false = zero cost): free-surface height from the FILL
   // integral (sum of fill over liquid+interface cells along z) of a regular grid of waterHeightGridX x
   // waterHeightGridY columns, every waterHeightEveryNSteps main steps, over the WHOLE run (unlike logDensityColumn's
   // single column/window). Per-rank CSV rows (mainStep, rank, gx, gy, partialFillSum, nObstacleCells): a column spans
   // several blocks in z, so the height of column (gx, gy) is the SUM of partialFillSum over ranks; nObstacleCells > 0
   // marks a column that currently lies under/inside a moving body (its height is NOT the free level). Purpose: track
   // the basin's water level at many points to see whether/where the level drifts (fluid-mass non-conservation test).
   bool        logWaterHeight;

   // [Claude, 2026-09-28] Option B (SetupStages.h repairDegenerateObstacleCells, THESIS_PLAN.md 2026-09-25): after every
   // collide, obstacle cells whose PDF sum is <= 0 or non-finite are reset to rest equilibrium. DEFAULT FALSE (off).
   // The hopper/rod-bundle clogging prms validated with the guard on (they crashed within ~3000 steps without it) must set
   // repairDegenerateObstaclePdfs true explicitly.
   bool        repairDegenerateObstaclePdfs;

   // [Claude, 2026-09-28] Option A (src/lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h): true -> seed the scratch PDF
   // clone of the stream-reconstruct-advect sweep from the real PDF field once, on its first creation (it is otherwise
   // uninitialised memory). DEFAULT FALSE (off = upstream behaviour). Only affects obstacle-cell PDFs, never fill/mass.
   bool        seedStreamPdfCloneOnFirstUse;
   std::string waterHeightBaseFolder;
   uint_t      waterHeightGridX;
   uint_t      waterHeightGridY;
   uint_t      waterHeightEveryNSteps;

   // 0 (DEFAULT) -> disabled, zero cost. Nonzero -> at exactly this main step (checked once,
   // subCycle 0), every floater this rank LOCALLY owns has every currently-mapped obstacle cell
   // forcibly cleared: moving-obstacle flag removed, PDFs reset to quiescent equilibrium (rho=1,
   // u=0), fill level reset to 1 (liquid), particleField_ reset to the invalid-UID sentinel -- the
   // exact same clear-and-reset sequence a since-reverted STALE-PARTICLE-IDX guard used on
   // genuinely orphaned cells (see THESIS_PLAN.md, 2026-09-26). TEST-ONLY: lets a small,
   // deterministic scenario reproduce that guard's own side effect (an abrupt, geometry-blind
   // fill-level discontinuity destabilising the free-surface curvature/normal reconstruction) on a
   // REAL, currently-valid floater, without needing the rare cross-rank ghost-ownership desync that
   // originally triggered it. Not meant for production scenarios.
   uint_t      debugForceObstacleClearAtStep;

   // false (DEFAULT) -> disabled, zero cost. true -> activates a single toggleable sweep, added once
   // per MAIN step (BuildTimeloop.h's addAfterTimeStepDiagnosticSweep, on timeloopEval so it runs
   // after hydro forces/DEM subcycles/everything else this step has settled), that: (1) drains and
   // logs SimulationObjects::diagnosticEvents -- the shared buffer several per-sub-cycle checks in
   // ParticleKernels.h push a lightweight record into instead of logging directly, so their output is
   // centralised here instead of scattered inline; (2) scans every liquid/interface cell in the
   // domain for density outside [densityPrecursorLow, densityPrecursorHigh] (same thresholds/
   // liquidOnly flag as the existing denseNanScan machinery -- no duplicate config surface) and logs
   // each one's world position; (3) logs every local floater's position, quaternion and hydrodynamic/
   // contact force+torque. Per explicit request ("introduce a toggleable sweep... checking the
   // density of every cell and flagging cells that have weird densities... log the floaters
   // position quaternion and forces"), 2026-09-27.
   bool        diagnosticSweepEnabled;

   // false (DEFAULT) -> no callback is passed to MovingParticleMapping (empty default): the original
   // pre-fix behaviour, a handed-off cell silently keeps whatever PDF history its previous owner left
   // it. true -> BuildTimeloop.h's onOwnershipHandoff callback is wired in instead: exactly once,
   // right before a genuine cross-particle obstacle-ownership handoff is committed, the handed-off
   // cell's PDF is reset to rest equilibrium (rho=1, u=0), breaking the silent-inheritance-of-stale-
   // history chain traced as the origin of a recurring growing-checkerboard STALE-PARTICLE-IDX crash
   // family (two rods' surfaces repeatedly handing a shared cell off with no reset -- THESIS_PLAN.md,
   // 2026-09-27, confirmed via real crashed-run data, a code read, and an independent local test
   // reproduction). Kept as a prm toggle (like the codebase's other A/B-testable mitigations) rather
   // than made permanently unconditional, so the fix can be switched on for comparison without a
   // rebuild; left off by default per explicit request.
   bool        useOwnershipHandoffPdfReset;

   // 0 (DEFAULT) -> diagnosticSweepEnabled (if true) runs from step 0. Nonzero -> skip every main
   // step before this value, same windowed-diagnostic convention as logFloaterHistoryStartStep.
   uint_t      diagnosticSweepFromStep;

   // 0 (DEFAULT) -> diagnosticSweepEnabled (if true) runs through to the end of the run. Nonzero ->
   // skip every main step at or after this value, same convention as logFloaterHistoryEndStep.
   uint_t      diagnosticSweepToStep;

   // true  → enables three lightweight WALBERLA_LOG_WARNING-based watchers inside ComputeContacts,
   //         active only within [contactAnomalyStartStep, contactAnomalyEndStep] below (same windowed-
   //         diagnostic spirit as logFloaterHistoryStartStep/EndStep, but independent of whether the
   //         heavier per-contact CSV logger is also on). Added investigating why Box/Cylinder floaters
   //         ("logs") occasionally build up 50%+ geometric overlap over many sub-cycles before any
   //         corrective force engages, then "explode" apart once contact is finally resolved -- a
   //         failure mode Sphere floaters in the identical jam never show:
   //           [NORMAL-FLIP]          -- the SAME uid pair's contact normal swings by more than ~45
   //                                     degrees between the sub-cycles it's seen in contact (the
   //                                     discrete-mode-switch symptom of SAT/Kodam-scenario dispatch
   //                                     near a near-symmetric configuration -- impossible for
   //                                     sphere-sphere, which has only one possible contact feature).
   //           [PROLONGED-CONTACT]    -- the SAME uid pair stays in contact for many consecutive
   //                                     sub-cycles; reports episode length and peak overlap seen so
   //                                     far, to distinguish a slow overlap build-up from a sudden one.
   //           [SIMULTANEOUS-CONTACTS] -- a single floater is resolved against 3+ other bodies within
   //                                     the SAME sub-cycle (the "many floaters piled on top of each
   //                                     other" packing case) -- each contact's spring-dashpot force is
   //                                     computed independently and summed, so the effective combined
   //                                     stiffness on that body may exceed what the explicit
   //                                     integrator's timestep was tuned to handle for one contact.
   // false → (DEFAULT) none of the three run; zero added cost.
   bool        logContactAnomalies;
   // 0 (DEFAULT) -> from step 0, if logContactAnomalies is on.
   uint_t      contactAnomalyStartStep;
   // 0 (DEFAULT) -> unbounded (through end of run), if logContactAnomalies is on.
   uint_t      contactAnomalyEndStep;

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

   // true  → (DEFAULT) exact pre-existing behaviour: reduceAndAvgHydroForces() (ParticleKernels.h)
   //         blends this step's freshly-reduced hydrodynamic force/torque with the PREVIOUS step's via
   //         AverageHydrodynamicForceTorqueKernel -- F = 0.5*(F(t)+F(t-Δt)), a rolling average (Ladd
   //         1994, "Numerical simulations of particulate suspensions via a discretized Boltzmann
   //         equation"). Watanabe et al. 2023 (Eq. 48-49, citing Aidun et al. 1998) use the identical
   //         mechanism "to improve the stability of the fluid-solid coupling with floating objects" --
   //         this codebase already had it wired in unconditionally before this flag existed.
   // false → skip the averaging kernel: DEM sees the raw just-reduced hydrodynamic force for this step,
   //         undamped by the previous step's value. A/B-testing knob for whether the averaging is
   //         actually load-bearing for this codebase's STALE-PARTICLE-IDX/runaway-force crash family, or
   //         masking a different problem (see ParticleKernels.h's [TOTAL-FORCE-THRESHOLD]/[VELOCITY-JUMP]
   //         diagnostics).
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
      , closeSealedCellsMinNeighbours { b.getParameter< uint_t >    ("closeSealedCellsMinNeighbours", uint_c(0)) }
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
      , logMomentumBudget           { b.getParameter< bool >       ("logMomentumBudget", false) }
      , guardCurvatureInputs        { b.getParameter< bool >       ("guardCurvatureInputs", false) }
      , wallOuterGhostLayer         { b.getParameter< bool >       ("wallOuterGhostLayer", false) }
      , diagnoseNaNNormals          { b.getParameter< bool >       ("diagnoseNaNNormals", false) }
      , useAnalyticBoxBuoyancy      { b.getParameter< bool >       ("useAnalyticBoxBuoyancy", false) }
      , usePhantomForceCorrection   { b.getParameter< bool >       ("usePhantomForceCorrection", false) }
      , useDynamicReferenceHeight   { b.getParameter< bool >       ("useDynamicReferenceHeight", false) }
      , staticReferenceHeightOverride{ b.getParameter< real_t >    ("staticReferenceHeightOverride", real_t(-1)) }
      , floaterPlaneFitRadiusFactor  { b.getParameter< real_t >    ("floaterPlaneFitRadiusFactor", real_t(1.8)) }
      , boxPlaneFitBandFactor        { b.getParameter< real_t >    ("boxPlaneFitBandFactor", real_t(0.8)) }
      , useAnalyticalContactDetection{ b.getParameter< bool >      ("useAnalyticalContactDetection", false) }
      , usePureAnalyticalCylinder   { b.getParameter< bool >       ("usePureAnalyticalCylinder", false) }
      , refineCylinderContacts      { b.getParameter< bool >       ("refineCylinderContacts", false) }
      , kodamBandEdgeAxialCheck     { b.getParameter< bool >       ("kodamBandEdgeAxialCheck", false) }
      , useTightTier1               { b.getParameter< bool >       ("useTightTier1", false) }
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
      , skipTier1                   { b.getParameter< bool >       ("skipTier1", false) }
      , useRoundedBoxCorners        { b.getParameter< bool >       ("useRoundedBoxCorners", false) }
      , boxRoundingRadius           { b.getParameter< real_t >     ("boxRoundingRadius", real_t(0.05)) }
      , useMultiPointBoxContact     { b.getParameter< bool >       ("useMultiPointBoxContact", false) }
      , useBoxOverlapCap            { b.getParameter< bool >       ("useBoxOverlapCap", false) }
      , boxOverlapCapFraction       { b.getParameter< real_t >     ("boxOverlapCapFraction", real_t(0.3)) }
      , useAxisHysteresis           { b.getParameter< bool >       ("useAxisHysteresis", false) }
      , axisHysteresisMargin        { b.getParameter< real_t >     ("axisHysteresisMargin", real_t(0.15)) }
      , useContactManifold          { b.getParameter< bool >       ("useContactManifold", false) }
      , manifoldFadeRate            { b.getParameter< real_t >     ("manifoldFadeRate", real_t(0.25)) }
      , useRotationalEffectiveMass  { b.getParameter< bool >       ("useRotationalEffectiveMass", false) }
      , useSpheroBoxContact         { b.getParameter< bool >       ("useSpheroBoxContact", false) }
      , spheroBoxRadius             { b.getParameter< real_t >     ("spheroBoxRadius", real_t(0.05)) }
      , cylinderContactSkin         { b.getParameter< real_t >     ("cylinderContactSkin", real_t(0)) }
      , boxContactSkin               { b.getParameter< real_t >     ("boxContactSkin", real_t(0)) }
      , useSpheroBoxAnalytical      { b.getParameter< bool >       ("useSpheroBoxAnalytical", false) }
      , useFraigeBoxContactModel    { b.getParameter< bool >       ("useFraigeBoxContactModel", false) }
      , logFloaterHistory           { b.getParameter< bool >       ("logFloaterHistory", false) }
      , floaterHistoryBaseFolder    { b.getParameter< std::string >("floaterHistoryBaseFolder", "floater-history") }
      , logFloaterHistoryReleaseStepFilter{ b.getParameter< uint_t >("logFloaterHistoryReleaseStepFilter", uint_t(0)) }
      , logFloaterHistoryStartStep   { b.getParameter< uint_t >     ("logFloaterHistoryStartStep", uint_t(0)) }
      , logFloaterHistoryEndStep     { b.getParameter< uint_t >     ("logFloaterHistoryEndStep", uint_t(0)) }
      , logFluidBoxHistory           { b.getParameter< bool >       ("logFluidBoxHistory", false) }
      , fluidBoxHistoryBaseFolder    { b.getParameter< std::string >("fluidBoxHistoryBaseFolder", "fluid-box-history") }
      // Note: getParameter< Vector3<real_t> >(name, defaultVal) (the two-argument overload) fails a
      // static_assert deep in Vector3's own generic single-arg constructor for this config binding --
      // sidestepped the same way clogUpstreamPos/clogPillarUpstreamPos above do, via isDefined() + a
      // ternary instead of relying on getParameter's own default mechanism for a Vector3.
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
      , debugForceObstacleClearAtStep{ b.getParameter< uint_t >     ("debugForceObstacleClearAtStep", uint_t(0)) }
      , diagnosticSweepEnabled       { b.getParameter< bool >       ("diagnosticSweepEnabled", false) }
      , useOwnershipHandoffPdfReset  { b.getParameter< bool >       ("useOwnershipHandoffPdfReset", false) }
      , diagnosticSweepFromStep      { b.getParameter< uint_t >     ("diagnosticSweepFromStep", uint_t(0)) }
      , diagnosticSweepToStep        { b.getParameter< uint_t >     ("diagnosticSweepToStep", uint_t(0)) }
      , logContactAnomalies          { b.getParameter< bool >       ("logContactAnomalies", false) }
      , contactAnomalyStartStep      { b.getParameter< uint_t >     ("contactAnomalyStartStep", uint_t(0)) }
      , contactAnomalyEndStep        { b.getParameter< uint_t >     ("contactAnomalyEndStep", uint_t(0)) }
      , highAltitudeWaterScanEnabled { b.getParameter< bool >       ("highAltitudeWaterScanEnabled", false) }
      , highAltitudeWaterThresholdZ  { b.getParameter< real_t >     ("highAltitudeWaterThresholdZ", real_t(0)) }
      , useArchTunnelAsMesaPDObject  { b.getParameter< bool >       ("useArchTunnelAsMesaPDObject", false) }
      , useContactFilterOwnerTiebreak{ b.getParameter< bool >       ("useContactFilterOwnerTiebreak", false) }
      , useHydrodynamicForceTimeAveraging{ b.getParameter< bool >   ("useHydrodynamicForceTimeAveraging", true) }
   {
      if (usePureAnalyticalCylinder && !useAnalyticalContactDetection)
         WALBERLA_LOG_WARNING("usePureAnalyticalCylinder=true but useAnalyticalContactDetection=false -- "
            "this flag only has an effect inside the Cylinder x Cylinder analytical dispatch, which is "
            "off, so it does nothing here.");
      if (useTightTier1 && !useAnalyticalContactDetection)
         WALBERLA_LOG_WARNING("useTightTier1=true but useAnalyticalContactDetection=false -- "
            "this flag only has an effect inside the Cylinder x Cylinder analytical dispatch's Tier 1 "
            "check, which is off, so it does nothing here.");
      if (skipTier1 && !useAnalyticalContactDetection)
         WALBERLA_LOG_WARNING("skipTier1=true but useAnalyticalContactDetection=false -- "
            "this flag only has an effect inside the Cylinder x Cylinder analytical dispatch's Tier 1 "
            "check, which is off, so it does nothing here.");
      if (useRoundedBoxCorners && !useAnalyticalContactDetection)
         WALBERLA_LOG_WARNING("useRoundedBoxCorners=true but useAnalyticalContactDetection=false -- "
            "this flag only has an effect inside boxBoxSAT, which is off (Box x Box uses GJK/EPA), so "
            "it does nothing here.");
      if (skipTier1 && !useTightLayer1)
         WALBERLA_LOG_WARNING("skipTier1=true but useTightLayer1=false -- this is NOT just a missed "
            "optimization, it is INCORRECT: without a sufficiently tight Layer 1 pre-filter upstream, "
            "skipping Tier 1 entirely will silently trust \"no contact\" for real dispatch-gap "
            "contacts the loose Layer 1 reject lets through. Set useTightLayer1=true as well, or turn "
            "skipTier1 back off.");
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
   // [Claude, 2026-09-22] Separate, OPTIONAL frequency for AverageDataSliceEvaluator (BuildTimeloop.h) -- that evaluator
   // writes THREE per-call text files (svfSlice_N.txt/fillSlice_N.txt/velXSlice_N.txt to eval/avgDataSlices/), a
   // completely different cost profile from the sqlite-row/log-line evaluators that also share evaluationFrequency.
   // Defaults to evaluationFrequency itself (0 = "not set", same behaviour as before this option existed) so every
   // existing prm is unaffected; set LOWER than evaluationFrequency to log e.g. [BoxForces]/[BoxTorque] every
   // timestep while keeping the slice dump throttled (a job that set evaluationFrequency=1 for per-step torque
   // logging otherwise floods eval/avgDataSlices/ with ~3 files/step -- 300000 files over 100000 steps blew a vault
   // disk-quota mid-archive and killed job 4255656 before its 2nd of 6 sub-runs even started).
   uint_t      dataSliceFrequency;
   std::string baseFolderName;
   // false (DEFAULT) -> ParticleKernels::logContactDetectionStatistics() is called exactly once, at
   // end-of-run, as before (no behaviour change for any existing scenario). true -> ALSO called
   // (paired with resetContactDetectionStatistics()) every performanceLogFrequency steps, so each
   // printed [ContactDetectionStats] line reports counts for just that interval instead of a
   // cumulative running total -- lets the cylinder-cylinder sub-case mix (and GJK fallback rate) be
   // tracked as a jam develops over a run, not just seen as one aggregate total at the end.
   bool        logContactDetectionStatsPeriodically;

   // Optional block (THESIS_PLAN.md item 3): probe positions for ClogMeasurementEvaluator, tracking
   // local water elevation and pressure at an "upstream" and a "downstream"/reference point, e.g.
   // either side of a pier/obstruction, to measure a developing clog's backwater effect directly
   // rather than only observing it visually. Parsed from an OPTIONAL "ClogMeasurement" block (same
   // pattern as VTKParams::denseEnabled's "fluid_field_dense" block below) -- absent block means
   // clogMeasurementEnabled is false and every existing prm is completely unaffected.
   bool               clogMeasurementEnabled;
   Vector3< real_t >  clogUpstreamPos;
   Vector3< real_t >  clogDownstreamPos;
   real_t             clogProbeHalfWidthXY;

   // Second, optional probe pair placed at the pillar's own Y (rather than the gap centerline) --
   // added 2026-09-14 per Schalko, Follett & Nepf (2023), Water Resources Research, "Impact of Lateral
   // Gap on Flow Distribution, Backwater Rise, and Turbulence Generated by a Logjam": that study found
   // real lateral variation in water depth (~7.5%+-1.8%) between a partial-width obstruction's own
   // blocked-flow region and the open gap beside it, and took SEPARATE longitudinal profiles at the
   // jam centre and the gap centre rather than assuming one location represents both -- the pillar-face
   // location is where the strongest backwater pile-up is actually expected, not the gap centreline
   // (which still has open flow through it). Optional sibling params `pillarUpstreamPos`/
   // `pillarDownstreamPos` inside the SAME "ClogMeasurement" block; absent means disabled (existing
   // prms using only upstreamPos/downstreamPos are completely unaffected).
   // Optional "flowRateEnabled" key inside the same "ClogMeasurement" block (default false): additionally writes the
   // volumetric flow rate Q and the wet cross-section area A of the y-z plane at the probes' x positions
   // (clogUpstreamFlowRate, clogDownstreamFlowRate, clogUpstreamWetArea, clogDownstreamWetArea; see
   // ClogMeasurementEvaluator::flowRateAt). Only the first probe pair (prefix "clog") writes them.
   bool               clogFlowRateEnabled;
   // Optional "blockageRatioEnabled" key inside the same "ClogMeasurement" block (default false):
   // additionally writes the fraction of the REFERENCE (undisturbed) wetted cross-section occupied
   // by the pillar(s) and the jam of floaters at the probes' x positions (clogUpstreamBlockageRatio,
   // clogDownstreamBlockageRatio, plus the raw clog*PillarArea/clog*JamArea components; see
   // ClogMeasurementEvaluator::blockageAt). Per review.md item 4 ("measure the blockage ratio instead
   // of only defining it ... plot afflux against blockage ratio, the relation Panici & Kripakaran
   // (2018) use"). Only the first probe pair (prefix "clog") writes them. Needs the pillar geometry
   // and channel width, threaded through from `simParams->cylinderPillars` and
   // `simParams->domain.domainSize[1]` at the BuildTimeloop.h call site, not read from this block.
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

// Parses a comma-separated list of unsigned integers, e.g. "1000, 2000, 8000,16000" -> {1000, 2000,
// 8000, 16000} (whitespace around each entry is ignored; a std::set both de-duplicates and sorts,
// neither of which the caller needs to care about). Used by VTKParams::snapshotSteps below.
//
// [Claude, 2026-09-29] FIX: config::Config::BlockHandle::getParameter<std::string>() returns a
// quoted prm value (e.g. steps "0,250,...";) WITH the surrounding quote characters still attached
// to the string -- confirmed by direct test against the real parser, not documented anywhere
// obvious. Left unstripped, the token adjacent to each quote mark fails/corrupts on extraction: a
// leading '"' before a digit makes operator>> fail immediately (silently dropping that whole
// entry -- this is exactly how "0" in "0,250,...,100000" went missing, while "100000" at the far
// end still parsed fine because operator>> just stops at trailing garbage instead of failing).
// Stripped here so the function behaves the same whether the prm value happens to be quoted or not.
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

   // A DISCRETE, typically non-uniform set of steps to dump fluid_field VTK at, instead of (or
   // alongside) the uniform every-writeFrequency-steps cadence above -- e.g. "0,1000,2000,4000,
   // 8000,16000,32000,64000,100000" to see a settling transient at increasing spacing plus one
   // exact end-of-run frame, without paying for a uniform fine cadence through the whole run (per
   // "implement the option so we can have a discrete set of vtk that we can print instead of every
   // thousand, like at 1k 2k 8k 16k..."). Parsed from an OPTIONAL third named VTK output block
   // "fluid_field_snapshots" (see VtkWriter.h's addVTKOutput, which gates its execution to exactly
   // this set of steps -- every other step it is never invoked at all, so it costs nothing outside
   // the listed steps). Same "absent block -> disabled, no change to existing prms" contract as
   // fluid_field_dense above. That block's own "writeFrequency" field is unused by this mechanism
   // (kept only because vtk::initializeVTKOutput requires it to construct the underlying VTKOutput
   // object) and should just be left at 1.
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
      // Surface tension override (>= 0 replaces the value from the units above; 0 = off)
      if (physics.surfaceTensionOverride >= real_c(0)) { surfaceTension = physics.surfaceTensionOverride; }

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
   std::vector< ContainerCylinderParams > containerCylinders; // [Claude, 2026-10-03] DropTest round container
   std::vector< ArchTunnelParams >      archTunnels;
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
      , archTunnels      { parseArchTunnels(cfg) }
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