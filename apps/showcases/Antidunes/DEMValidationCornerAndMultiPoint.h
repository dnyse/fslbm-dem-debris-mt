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
//! \file DEMValidationCornerAndMultiPoint.h
//! \brief Corner-contact-normal-discontinuity and multi-point box-contact DEM validation scenarios
//! -- split out of DEMValidationScenarios.h 2026-10-29. runCornerNormalDiscontinuity,
//! runCornerFixComparison, runMultiPointDeepOverlapStability/Fuzz, runMultiPointEdgeCases,
//! runContactStiffnessCalibrationSweep, runFallOntoPileStability, runMultiPointBoxContactFuzz.

#pragma once

#include "DEMValidationInfra.h"

namespace walberla
{
namespace antidunes
{
namespace dem_validation
{
using namespace mesa_pd;

//======================================================================================================================
// 5b: corner/edge contact-normal discontinuity (Box x Box analytical dispatch) -- REPRODUCTION test
//======================================================================================================================

// Investigates a real defect spotted visually in production job 4196224 (bridge "widerstacks" box
// scenario, ParaView output around step 26000-27000): a batch of staggered-release plank floaters
// lands together and jams, and one plank is subsequently ejected/"flies" out of the pile with no
// single log line ever showing a force large enough to explain it. Literature on polyhedral DEM (e.g.
// the "improved common-plane concept" and "energy-conserving contact" papers for arbitrarily-shaped
// elements) identifies the mechanism directly: near a corner-to-corner or corner-to-edge
// configuration, a tiny relative positional/rotational change can flip which face/edge/vertex pair the
// contact solver reports as the separating axis, snapping the contact normal to a different direction
// between consecutive steps. Because the spring-dashpot force is applied along whatever normal is
// reported THIS step with no continuity constraint against the PREVIOUS step's normal, that snap can
// inject kinetic energy with no corresponding physical source -- exactly what turns an ordinary
// pile-up into a "flying" launch.
//
// This reproduces the geometry realistically, not as a contrived worst case: both bodies use this
// project's own production FloatingBox size/density/Poisson's ratio (test_bridge_v1_..._boxes_..._gas.
// prm's own FloatingBox halfExtents <6,1,1.5>, densityRatio 0.6, and the scenario's own derived
// particleCollisionTime/particleTimestepSize/particleKappa from computedParams.txt), the SAME
// analytical Box x Box dispatch (boxBoxSAT) that scenario actually runs with
// useAnalyticalContactDetection=true -- NOT the generic GJK/EPA fallback every OTHER DEMValidation
// scenario in this file defaults to, since GJK/EPA and boxBoxSAT are different code with different
// failure modes and this defect is specific to the analytical path -- and a realistic approach speed
// (~0.025 LU/step, matching the ~0.01-0.03 LU/step range actually observed in job 4196224's own
// STRONG-CONTACT-FORCE log lines, i.e. what free-fall under this scenario's own tiny gravity,
// 7.03e-6 LU/step^2, produces after a plausible multi-thousand-step drop -- simulating the drop itself
// would just be a long uneventful prelude, so this starts right at that realistic pre-impact speed).
//
// Two identical planks approach roughly along their own long (X) axis -- exactly how two staggered-
// release logs converging on the same jam location would meet -- but offset slightly in Y/Z and tilted
// a few degrees (a plausible natural tumble; ParticleKernels.h's own FloaterHistoryLogger comment notes
// a floater's angular velocity "is often nonzero well before any 'expected' collision"), so their
// SHORT end-cap CORNERS meet first rather than a clean face-to-face impact -- precisely the
// near-degenerate configuration the literature identifies as the trigger. e=1 (perfectly elastic,
// frictionless, no gravity) is used for the PASS/FAIL bound itself -- the strictest, least ambiguous
// bar: total system kinetic energy can dip during compression (converted to non-negative spring
// potential energy) but can never legitimately exceed its own pre-contact value at any instant, since
// there is no other energy source in this system. This is exactly the "ECC" criterion
// runEnergyConservation above already established for this file, applied here as a running bound
// across every step instead of only checked once before/after. restitutionCoefficient/
// frictionCoefficient remain overridable via the prm to rerun under the scenario's literal production
// damping (0.3 "wet/waterlogged wood", 0.5) instead, if a looser reproduction is wanted -- with e<1,
// mu>0 the bound is no longer "exactly conserved" but "monotonically non-increasing once separated",
// so keGainRel is still meaningful as a smoking gun, just not as mathematically clean.
namespace detail
{
//! Shared parameters for one run of the corner/edge contact-normal-discontinuity scenario -- see
//! runCornerNormalDiscontinuity's own (still-present, unchanged) header comment for the full
//! literature/production-incident background and the realism rationale behind every default value
//! here. Factored out so runCornerFixComparison below can run the exact SAME geometry/physics several
//! times, varying only useRoundedBoxCorners/useMultiPointBoxContact.
struct CornerScenarioParams
{
   Vec3   halfExtents      = Vec3(6, 1, 1.5);
   real_t density          = real_t(0.6);
   real_t restitution      = real_t(1);
   real_t friction         = real_t(0);
   real_t poisson          = real_t(0.22);
   real_t collisionTime    = real_t(32.0078);
   real_t dt               = real_t(0.05);
   real_t approachSpeed    = real_t(0.025);
   real_t lateralOffsetY   = real_t(0.6);
   real_t lateralOffsetZ   = real_t(0.9);
   real_t tiltDegreesA     = real_t(8);
   real_t tiltDegreesB     = real_t(-11);
   real_t spinRateB        = real_t(0.004);
   real_t separation       = real_t(20);
   real_t tolerance        = real_t(0.02);
   real_t flipCosThreshold = real_t(0.7); // ~45 degrees, matching ParticleKernels.h's own [NORMAL-FLIP]
   uint_t maxSteps         = uint_t(20000);
   uint_t vtkSpacing       = uint_t(0);
   std::string outputBaseName = "corner_normal_discontinuity";
   // false skips openCsv/setupVtk entirely -- for a large fuzz sweep (runMultiPointBoxContactFuzz)
   // calling this thousands of times, per-trial CSV file I/O would dominate the runtime for no benefit
   // (nobody inspects 2000 near-identical single-trial CSVs).
   bool enableOutput = true;

   // The two togglable mitigations under test -- see SimulationParams.h's
   // useRoundedBoxCorners/useMultiPointBoxContact doc comments. Both default OFF (exact pre-existing
   // behaviour, i.e. the original reproduction case).
   bool   useRoundedBoxCorners   = false;
   real_t boxRoundingRadius      = real_t(0.05);
   bool   useMultiPointBoxContact = false;
   bool   useBoxOverlapCap       = false;
   real_t boxOverlapCapFraction  = real_t(0.3);
   bool   useAxisHysteresis      = false;
   real_t axisHysteresisMargin   = real_t(0.15);
   bool   useContactManifold     = false;
   real_t manifoldFadeRate       = real_t(0.25);
   bool   useRotationalEffectiveMass = false;
};

struct CornerScenarioResult
{
   real_t keBefore        = real_t(0);
   real_t keAfter         = real_t(0);
   real_t keMax           = real_t(0);
   uint_t keMaxStep       = uint_t(0);
   real_t peInitial       = real_t(0); // elastic PE already stored in any t=0 overlap, see
                                       // computeInitialElasticPE's own comment -- 0 for every scenario
                                       // that starts from clean separation (the common case), nonzero
                                       // only for the deep-overlap scenarios that start pre-embedded.
   real_t keGainRel       = real_t(0); // (keMax - (keBefore+peInitial)) / max(1e-12, keBefore+peInitial)
   uint_t normalFlipCount = uint_t(0);
   // The normalFlipCount>0 diagnostic only ever catches a SUDDEN >45-degree single-step jump in the
   // raw SAT normal -- this session's own multiPointDeepOverlapStability data confirms that's real but
   // INCOMPLETE: baseline reproducibly injects energy at several depths with normalFlipCount==0. Two
   // candidate mechanisms this doesn't see: (a) the normal rotating a LOT but gradually (never crossing
   // 45 degrees in any single step, so never flagged, yet still non-conservative -- a linear spring is
   // only exactly conservative along a FIXED axis), and (b) the reported contact POINT jumping (torque
   // arm changes discontinuously) while the normal itself stays stable. Both tracked directly below.
   real_t cumulativeNormalRotationDeg = real_t(0); // sum of EVERY step's angle change, no threshold
   real_t maxPointJumpFrac            = real_t(0); // largest single-step contact-point jump / halfExtents.length()
   uint_t pointJumpCount              = uint_t(0);  // steps where that jump exceeded pointJumpThreshold
   bool   everContacted   = false;
   bool   pass            = false; // keGainRel < tolerance (no anomalous energy gain observed)
   // [Claude, 2026-10-29] the four "eberly*" diagnostics that used to live here (tracking
   // paper_analytical::boxBoxSpheroAnalytical's own reported normal/point) were removed along with
   // useSpheroBoxContact/useSpheroBoxAnalytical -- see archive/unused_contact_methods/.
};

// The "no energy from nowhere" bound every scenario in this file checks is: total mechanical energy at
// any instant <= total mechanical energy at t=0. For every scenario that starts from clean separation
// (approachSpeed>0, zero initial overlap -- runCornerNormalDiscontinuity, runMultiPointBoxContactFuzz,
// etc.) t=0's mechanical energy IS just KE_before, since there is no overlap yet to store any elastic
// potential energy in. But the deep-overlap scenarios (runMultiPointDeepOverlapStability/Fuzz/
// EdgeCases) deliberately start the two bodies ALREADY overlapping at rest -- KE_before is exactly 0
// there, but t=0's true mechanical energy is NOT 0, it is the elastic PE already locked into that
// overlap by construction. Comparing keMax against KE_before alone (as an earlier version of these
// deep-overlap tests did) makes every contacted trial report an appraoch to infinite relative gain
// regardless of whether the config actually behaves correctly, since it divides by ~0 -- this function
// computes the missing term so the pass/fail bound is physically meaningful even from a static release.
//
// Uses the SAME formula ComputeContacts/LinearSpringDashpot itself uses (LinearSpringDashpot.h's own
// setStiffnessAndDamping: stiffnessN = effectiveMass*(pi^2 + ln(e)^2)/collisionTime^2, PE = 1/2*k*x^2 for
// a linear spring compressed by x) and the SAME contact-point set the dynamics will actually use this
// run (boxBoxMultiPointContacts's vertex list -- or its own SAT-fallback single point when that list is
// empty -- if useMultiPointBoxContact, otherwise boxBoxSAT's own single overlap): each point's force is
// scaled by 1/N in ComputeContacts's multi-point loop, so its effective per-point stiffness is
// stiffnessN/N, giving PE = (stiffnessN/(2N)) * sum(overlap_i^2) -- reduces to the ordinary 1/2*k*x^2
// for the N=1 single-point case (baseline/rounded, or multipoint with no vertex embedded).
// pointPE: potential energy stored by ONE point of effective stiffness k compressed by the TRUE
// (uncapped) overlap x, under the SAME piecewise force law ComputeContacts actually applies when a cap
// is active: F(x) = k*x for x<=cap (ordinary linear spring), F(x) = k*cap (constant, clipped) for
// x>cap. This is still a perfectly conservative force (a well-defined potential exists for any
// force that depends only on x), so integrating it gives the CORRECT stored energy at depth x --
// 1/2*k*cap^2 (the ordinary spring energy up to the cap) PLUS k*cap*(x-cap) (the extra distance
// travelled at the now-constant force) for x>cap. Using the naive 1/2*k*min(x,cap)^2 instead (as an
// earlier version of this function did) silently THROWS AWAY that second term, understating the true
// starting energy at any overlap deeper than the cap and making a correctly-behaving capped
// contact look like it "gained" energy from nowhere when it did not.
inline real_t pointPE(real_t k, real_t x, real_t cap)
{
   if (x <= cap) return real_t(0.5) * k * x * x;
   return real_t(0.5) * k * cap * cap + k * cap * (x - cap);
}

inline real_t computeInitialElasticPE(const Vec3& posA, const Rot3& rotA, const Vec3& halfA, real_t densityA,
                                       const Vec3& posB, const Rot3& rotB, const Vec3& halfB, real_t densityB,
                                       real_t restitution, real_t collisionTime, bool useRoundedBoxCorners,
                                       real_t boxRoundingRadius, bool useMultiPointBoxContact,
                                       bool useBoxOverlapCap = false, real_t boxOverlapCapFraction = real_t(0.3),
                                       bool useRotationalEffectiveMass = false)
{
   const real_t massA = densityA * real_t(8) * halfA[0] * halfA[1] * halfA[2];
   const real_t massB = densityB * real_t(8) * halfB[0] * halfB[1] * halfB[2];
   const real_t meffTranslational = real_t(1) / (real_t(1) / massA + real_t(1) / massB);
   const real_t lnE   = std::log(restitution);

   // Same generalized effective-mass formula as ParticleKernels.h's own computeRotationalEffectiveMass
   // (see ModelParameters::useRotationalEffectiveMass's doc comment) -- re-derived here from raw geometry
   // since this function has no live mesa_pd Accessor/ParticleStorage to query. Each box's world-frame
   // inverse inertia only depends on its (static, at t=0) rotation, so it's built once per box and reused
   // for whichever single (point,normal) ends up chosen as PRIMARY below -- exactly mirroring
   // ComputeContacts' own behaviour, where ONE meff (from the primary SAT point/normal, computed BEFORE
   // branching into the multi-point loop) sets stiffnessN/dampingN for every point applied that sub-cycle.
   Mat3 invInertiaWF_A, invInertiaWF_B;
   if (useRotationalEffectiveMass)
   {
      FloatingBox boxA(real_t(2) * halfA);
      boxA.updateMassAndInertia(densityA);
      FloatingBox boxB(real_t(2) * halfB);
      boxB.updateMassAndInertia(densityB);
      invInertiaWF_A = math::transformMatrixRART(rotA.getMatrix(), boxA.getInvInertiaBF());
      invInertiaWF_B = math::transformMatrixRART(rotB.getMatrix(), boxB.getInvInertiaBF());
   }
   const auto stiffnessAtPrimary = [&](const Vec3& point, const Vec3& normal) -> real_t {
      real_t meff = meffTranslational;
      if (useRotationalEffectiveMass)
      {
         const Vec3 rA = point - posA;
         const Vec3 rB = point - posB;
         const Vec3 crossA = rA % normal;
         const Vec3 crossB = rB % normal;
         real_t invMeff = real_t(1) / massA + real_t(1) / massB;
         if (crossA.sqrLength() > real_t(0)) invMeff += crossA * (invInertiaWF_A * crossA);
         if (crossB.sqrLength() > real_t(0)) invMeff += crossB * (invInertiaWF_B * crossB);
         meff = real_t(1) / invMeff;
      }
      return meff * (math::pi * math::pi + lnE * lnE) / (collisionTime * collisionTime);
   };
   // Must match ComputeContacts' own per-pair cap exactly (ParticleKernels.h) so this ceiling reflects
   // the SAME cap value the actual dynamics will use. Always computed from the TRUE (uncapped)
   // geometric overlap below -- see pointPE's own comment for why the cap changes HOW the energy is
   // computed (piecewise), not whether the true starting overlap/energy is smaller.
   const real_t overlapCap = useBoxOverlapCap
      ? boxOverlapCapFraction * std::min({ halfA[0], halfA[1], halfA[2], halfB[0], halfB[1], halfB[2] })
      : std::numeric_limits< real_t >::max();

   // [Claude, 2026-10-29] useSpheroBoxContact branch removed -- see archive/unused_contact_methods/.

   if (useMultiPointBoxContact)
   {
      // Uncapped on purpose -- need the TRUE overlap at each point for pointPE below, not whatever a
      // cap would have clipped it to.
      auto points = paper_analytical::boxBoxMultiPointContacts(posA, rotA.getMatrix(), halfA, posB,
                                                                 rotB.getMatrix(), halfB);
      const auto sat = paper_analytical::boxBoxSAT(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(),
                                                     halfB, useRoundedBoxCorners, boxRoundingRadius);
      if (points.empty())
      {
         if (!sat.contact) return real_t(0);
         return pointPE(stiffnessAtPrimary(sat.point, sat.normal), sat.overlap, overlapCap);
      }
      // PRIMARY point/normal: sat's own result if it found contact, else the deepest vertex -- exactly
      // matching ParticleKernels.h's Box x Box operator() (see its own comment right after
      // "SAT itself found no contact, but useMultiPointBoxContact's independent vertex test did").
      Vec3 primaryPoint, primaryNormal;
      if (sat.contact) { primaryPoint = sat.point; primaryNormal = sat.normal; }
      else
      {
         const auto deepest = std::max_element(points.begin(), points.end(),
            [](const auto& a, const auto& b) { return a.overlap < b.overlap; });
         primaryPoint = deepest->point; primaryNormal = deepest->normal;
      }
      const real_t stiffnessN = stiffnessAtPrimary(primaryPoint, primaryNormal);
      const real_t kEff = stiffnessN / real_t(points.size()); // matches ComputeContacts' own 1/N force scaling
      real_t sum = real_t(0);
      for (const auto& pt : points) sum += pointPE(kEff, pt.overlap, overlapCap);
      return sum;
   }

   const auto sat = paper_analytical::boxBoxSAT(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(), halfB,
                                                 useRoundedBoxCorners, boxRoundingRadius);
   if (!sat.contact) return real_t(0);
   return pointPE(stiffnessAtPrimary(sat.point, sat.normal), sat.overlap, overlapCap);
}

inline CornerScenarioResult runCornerScenarioOnce(const CornerScenarioParams& p)
{
   // useAnalyticalContactDetection=true: the whole point -- reproduce against the ACTUAL production
   // Box x Box dispatch (boxBoxSAT), not the GJK/EPA fallback every other scenario in this file uses.
   DEMEngine engine(p.dt, p.restitution, p.friction, p.poisson, p.collisionTime, real_t(0),
                     p.separation + real_t(40), /*useAnalyticalContactDetection=*/true,
                     p.useRoundedBoxCorners, p.boxRoundingRadius, p.useMultiPointBoxContact,
                     p.useBoxOverlapCap, p.boxOverlapCapFraction, p.useAxisHysteresis,
                     p.axisHysteresisMargin, p.useContactManifold, p.manifoldFadeRate,
                     p.useRotationalEffectiveMass);

   const Rot3 rotA(Vec3(0, 0, 1), p.tiltDegreesA * math::pi / real_t(180));
   const Rot3 rotB(Vec3(0, 0, 1), p.tiltDegreesB * math::pi / real_t(180));

   const size_t idxA = engine.addBox(Vec3(-p.separation / 2, 0, 0), Vec3(p.approachSpeed / 2, 0, 0), rotA,
                                      real_t(2) * p.halfExtents, p.density);
   const size_t idxB = engine.addBox(Vec3(p.separation / 2, p.lateralOffsetY, p.lateralOffsetZ),
                                      Vec3(-p.approachSpeed / 2, 0, 0), rotB, real_t(2) * p.halfExtents,
                                      p.density);
   // Small residual tumble on B ONLY, about an axis (Y) that keeps sweeping its end-cap corner through
   // A's corner region as the two close in -- NOT an artificial imposed spin to force the bug, but the
   // same order of magnitude as the pre-contact angular velocities this scenario's own floaters
   // routinely carry (see runCornerNormalDiscontinuity's own header comment).
   engine.ac().setAngularVelocity(idxB, Vec3(0, p.spinRateB, 0));

   engine.finalizeSetup();

   if (p.enableOutput)
   {
      engine.openCsv(p.outputBaseName + ".csv",
                      { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { idxA, idxB });
      engine.setupVtk(p.outputBaseName, p.vtkSpacing);
   }

   CornerScenarioResult result;
   result.keBefore  = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
   result.peInitial = computeInitialElasticPE(
      Vec3(-p.separation / 2, 0, 0), rotA, p.halfExtents, p.density,
      Vec3(p.separation / 2, p.lateralOffsetY, p.lateralOffsetZ), rotB, p.halfExtents, p.density,
      p.restitution, p.collisionTime, p.useRoundedBoxCorners, p.boxRoundingRadius, p.useMultiPointBoxContact,
      p.useBoxOverlapCap, p.boxOverlapCapFraction, p.useRotationalEffectiveMass);
   result.keMax    = result.keBefore;

   real_t penetration = real_t(0);
   Vec3   lastNormal;
   Vec3   lastPoint;
   bool   haveLastNormal = false;
   const real_t pointJumpThreshold = real_t(0.15); // fraction of halfExtents.length(), see doc comment above

   for (uint_t s = 0; s < p.maxSteps; ++s)
   {
      const bool contactNow = engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0);
      if (contactNow)
      {
         result.everContacted = true;
         // Diagnostic-only: reads back the SAME closed-form dispatch that is already driving the
         // dynamics this step, purely to record its own reported normal/point -- not a second/
         // different code path. Deliberately does NOT pass p.useRoundedBoxCorners here: this
         // diagnostic exists to characterize the RAW SAT axis's own discontinuity regardless of which
         // mitigation (if any) is actually resolving the dynamics this run.
         const auto res = paper_analytical::boxBoxSAT(engine.ac().getPosition(idxA),
                                                        engine.ac().getRotation(idxA).getMatrix(), p.halfExtents,
                                                        engine.ac().getPosition(idxB),
                                                        engine.ac().getRotation(idxB).getMatrix(), p.halfExtents);
         if (res.contact)
         {
            if (haveLastNormal)
            {
               const real_t cosAngle = std::max(real_t(-1), std::min(real_t(1), res.normal * lastNormal));
               if (cosAngle < p.flipCosThreshold) ++result.normalFlipCount;
               // Unconditional -- unlike normalFlipCount, this accumulates EVERY step's angle change,
               // so a normal that rotates a lot gradually (never crossing 45 degrees in one step) still
               // shows up here, distinguishing "genuinely fixed axis" from "slowly but continuously
               // rotating axis" -- both look identical to normalFlipCount, only the latter can still
               // inject energy (a linear spring is only exactly conservative along a FIXED direction).
               result.cumulativeNormalRotationDeg += std::acos(cosAngle) * real_t(180) / math::pi;

               const real_t pointJumpFrac = (res.point - lastPoint).length() / p.halfExtents.length();
               if (pointJumpFrac > result.maxPointJumpFrac) result.maxPointJumpFrac = pointJumpFrac;
               if (pointJumpFrac > pointJumpThreshold) ++result.pointJumpCount;
            }
            lastNormal     = res.normal;
            lastPoint      = res.point;
            haveLastNormal = true;
         }
         // [Claude, 2026-10-29] the useSpheroBoxContact+useSpheroBoxAnalytical Eberly-diagnostic
         // block that used to live here was removed -- see archive/unused_contact_methods/.
      }
      engine.step();
      engine.writeCsvRow();
      engine.writeVtkIfDue();

      const real_t keNow = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      if (keNow > result.keMax)
      {
         result.keMax     = keNow;
         result.keMaxStep = s;
      }

      if (result.everContacted && !contactNow && s > 0) break;
   }

   result.keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
   // Ceiling = KE_before + the elastic PE already stored in any t=0 overlap (0 for every scenario that
   // starts from clean separation -- see computeInitialElasticPE's own comment for why this term
   // matters for the deep-overlap scenarios, which start pre-embedded at rest with KE_before=0).
   const real_t ceiling = result.keBefore + result.peInitial;
   result.keGainRel = (result.keMax - ceiling) / std::max(real_t(1e-12), ceiling);
   result.pass      = result.everContacted && result.keGainRel < p.tolerance;
   return result;
}
} // namespace detail

// Investigates a real defect spotted visually in production job 4196224 (bridge "widerstacks" box
// scenario, ParaView output around step 26000-27000): a batch of staggered-release plank floaters
// lands together and jams, and one plank is subsequently ejected/"flies" out of the pile with no
// single log line ever showing a force large enough to explain it. Literature on polyhedral DEM (e.g.
// the "improved common-plane concept" and "energy-conserving contact" papers for arbitrarily-shaped
// elements) identifies the mechanism directly: near a corner-to-corner or corner-to-edge
// configuration, a tiny relative positional/rotational change can flip which face/edge/vertex pair the
// contact solver reports as the separating axis, snapping the contact normal to a different direction
// between consecutive steps. Because the spring-dashpot force is applied along whatever normal is
// reported THIS step with no continuity constraint against the PREVIOUS step's normal, that snap can
// inject kinetic energy with no corresponding physical source -- exactly what turns an ordinary
// pile-up into a "flying" launch.
//
// This reproduces the geometry realistically, not as a contrived worst case: both bodies use this
// project's own production FloatingBox size/density/Poisson's ratio (test_bridge_v1_..._boxes_..._gas.
// prm's own FloatingBox halfExtents <6,1,1.5>, densityRatio 0.6, and the scenario's own derived
// particleCollisionTime/particleTimestepSize/particleKappa from computedParams.txt), the SAME
// analytical Box x Box dispatch (boxBoxSAT) that scenario actually runs with
// useAnalyticalContactDetection=true -- NOT the generic GJK/EPA fallback every OTHER DEMValidation
// scenario in this file defaults to, since GJK/EPA and boxBoxSAT are different code with different
// failure modes and this defect is specific to the analytical path -- and a realistic approach speed
// (~0.025 LU/step, matching the ~0.01-0.03 LU/step range actually observed in job 4196224's own
// STRONG-CONTACT-FORCE log lines, i.e. what free-fall under this scenario's own tiny gravity,
// 7.03e-6 LU/step^2, produces after a plausible multi-thousand-step drop -- simulating the drop itself
// would just be a long uneventful prelude, so this starts right at that realistic pre-impact speed).
//
// Two identical planks approach roughly along their own long (X) axis -- exactly how two staggered-
// release logs converging on the same jam location would meet -- but offset slightly in Y/Z and tilted
// a few degrees (a plausible natural tumble; ParticleKernels.h's own FloaterHistoryLogger comment notes
// a floater's angular velocity "is often nonzero well before any 'expected' collision"), so their
// SHORT end-cap CORNERS meet first rather than a clean face-to-face impact -- precisely the
// near-degenerate configuration the literature identifies as the trigger. e=1 (perfectly elastic,
// frictionless, no gravity) is used for the PASS/FAIL bound itself -- the strictest, least ambiguous
// bar: total system kinetic energy can dip during compression (converted to non-negative spring
// potential energy) but can never legitimately exceed its own pre-contact value at any instant, since
// there is no other energy source in this system. useRoundedBoxCorners/useMultiPointBoxContact
// (default OFF, matching the original reproduction) can be set true from the prm to instead verify one
// or both mitigations clear this bound -- see runCornerFixComparison below for a single run that
// checks all four combinations side by side.
inline bool runCornerNormalDiscontinuity(const config::Config::BlockHandle& b)
{
   detail::CornerScenarioParams p;
   p.halfExtents             = b.getParameter< Vec3 >("halfExtents", p.halfExtents);
   p.density                 = b.getParameter< real_t >("density", p.density);
   p.restitution             = b.getParameter< real_t >("restitutionCoefficient", p.restitution);
   p.friction                = b.getParameter< real_t >("frictionCoefficient", p.friction);
   p.poisson                 = b.getParameter< real_t >("poissonsRatio", p.poisson);
   p.collisionTime           = b.getParameter< real_t >("collisionTime", p.collisionTime);
   p.dt                      = b.getParameter< real_t >("dt", p.dt);
   p.approachSpeed           = b.getParameter< real_t >("approachSpeed", p.approachSpeed);
   p.lateralOffsetY          = b.getParameter< real_t >("lateralOffsetY", p.lateralOffsetY);
   p.lateralOffsetZ          = b.getParameter< real_t >("lateralOffsetZ", p.lateralOffsetZ);
   p.tiltDegreesA            = b.getParameter< real_t >("tiltDegreesA", p.tiltDegreesA);
   p.tiltDegreesB            = b.getParameter< real_t >("tiltDegreesB", p.tiltDegreesB);
   p.spinRateB               = b.getParameter< real_t >("spinRateB", p.spinRateB);
   p.separation              = b.getParameter< real_t >("initialSeparation", p.separation);
   p.tolerance               = b.getParameter< real_t >("toleranceRelError", p.tolerance);
   p.flipCosThreshold        = b.getParameter< real_t >("flipCosThreshold", p.flipCosThreshold);
   p.maxSteps                = b.getParameter< uint_t >("maxSteps", p.maxSteps);
   p.vtkSpacing              = b.getParameter< uint_t >("vtkSpacing", p.vtkSpacing);
   p.outputBaseName          = b.getParameter< std::string >("outputBaseName", p.outputBaseName);
   p.useRoundedBoxCorners    = b.getParameter< bool >("useRoundedBoxCorners", p.useRoundedBoxCorners);
   p.boxRoundingRadius       = b.getParameter< real_t >("boxRoundingRadius", p.boxRoundingRadius);
   p.useMultiPointBoxContact = b.getParameter< bool >("useMultiPointBoxContact", p.useMultiPointBoxContact);

   const auto r = detail::runCornerScenarioOnce(p);

   detail::PassFailAccumulator acc;
   acc.record(r.pass);

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerNormalDiscontinuity  KE_before=" << r.keBefore
      << " KE_after=" << r.keAfter << " KE_max=" << r.keMax << " (at step " << r.keMaxStep << ")"
      << " maxRelGain=" << r.keGainRel << " normalFlipCount=" << r.normalFlipCount
      << " everContacted=" << r.everContacted << "  " << (r.pass ? "PASS" : "FAIL (energy injected)"));

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerNormalDiscontinuity OVERALL: "
      << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total << ")");
   return acc.allPassed();
}

// Runs the EXACT SAME corner/tilt/offset geometry as runCornerNormalDiscontinuity above four times,
// varying only useRoundedBoxCorners/useMultiPointBoxContact, and reports a side-by-side comparison --
// the direct answer to "do the two togglable mitigations actually fix the reproduced defect, and does
// combining them help further or just cost more". Geometry/physics parameters are read from the SAME
// prm keys as runCornerNormalDiscontinuity (minus useRoundedBoxCorners/useMultiPointBoxContact
// themselves, which this scenario controls); outputBaseName is suffixed per configuration so all four
// CSVs/VTK series are kept side by side rather than overwriting each other.
inline bool runCornerFixComparison(const config::Config::BlockHandle& b)
{
   detail::CornerScenarioParams p;
   p.halfExtents      = b.getParameter< Vec3 >("halfExtents", p.halfExtents);
   p.density          = b.getParameter< real_t >("density", p.density);
   p.restitution      = b.getParameter< real_t >("restitutionCoefficient", p.restitution);
   p.friction         = b.getParameter< real_t >("frictionCoefficient", p.friction);
   p.poisson          = b.getParameter< real_t >("poissonsRatio", p.poisson);
   p.collisionTime    = b.getParameter< real_t >("collisionTime", p.collisionTime);
   p.dt               = b.getParameter< real_t >("dt", p.dt);
   p.approachSpeed    = b.getParameter< real_t >("approachSpeed", p.approachSpeed);
   p.lateralOffsetY   = b.getParameter< real_t >("lateralOffsetY", p.lateralOffsetY);
   p.lateralOffsetZ   = b.getParameter< real_t >("lateralOffsetZ", p.lateralOffsetZ);
   p.tiltDegreesA     = b.getParameter< real_t >("tiltDegreesA", p.tiltDegreesA);
   p.tiltDegreesB     = b.getParameter< real_t >("tiltDegreesB", p.tiltDegreesB);
   p.spinRateB        = b.getParameter< real_t >("spinRateB", p.spinRateB);
   p.separation       = b.getParameter< real_t >("initialSeparation", p.separation);
   p.tolerance        = b.getParameter< real_t >("toleranceRelError", p.tolerance);
   p.flipCosThreshold = b.getParameter< real_t >("flipCosThreshold", p.flipCosThreshold);
   p.maxSteps         = b.getParameter< uint_t >("maxSteps", p.maxSteps);
   p.vtkSpacing       = b.getParameter< uint_t >("vtkSpacing", p.vtkSpacing);
   const real_t boxRoundingRadius = b.getParameter< real_t >("boxRoundingRadius", p.boxRoundingRadius);
   p.boxRoundingRadius = boxRoundingRadius;

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
   };
   const Config configs[] = {
      { "baseline (both OFF)", false, false },
      { "useRoundedBoxCorners", true, false },
      { "useMultiPointBoxContact", false, true },
      { "both", true, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison -- same corner/tilt/offset geometry "
      "(the exact case runCornerNormalDiscontinuity reproduces), 4 configurations:");

   uint_t configsWithoutInjection = 0;
   bool   baselineFailed          = false;
   bool   anyFixPassed            = false;
   for (const auto& cfg : configs)
   {
      p.useRoundedBoxCorners    = cfg.rounded;
      p.useMultiPointBoxContact = cfg.multiPoint;
      p.outputBaseName          = std::string("corner_fix_comparison_") + cfg.name;

      const auto r = detail::runCornerScenarioOnce(p);
      if (r.pass) ++configsWithoutInjection;
      if (std::string(cfg.name) == "baseline (both OFF)" && !r.pass) baselineFailed = true;
      if (std::string(cfg.name) != "baseline (both OFF)" && r.pass) anyFixPassed = true;

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison  config=\"" << cfg.name << "\""
         << "  KE_before=" << r.keBefore << " KE_max=" << r.keMax << " (at step " << r.keMaxStep << ")"
         << " maxRelGain=" << r.keGainRel << " normalFlipCount=" << r.normalFlipCount
         << " everContacted=" << r.everContacted
         << "  " << (r.pass ? "PASS (no injection)" : "FAIL (energy injected)"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison OVERALL: " << configsWithoutInjection
      << "/" << (sizeof(configs) / sizeof(configs[0])) << " configurations avoided energy injection");

   // This is a comparison/characterization tool, not a strict per-scenario regression gate -- the
   // meaningful pass criterion is "the defect is still confirmed present at baseline (otherwise this
   // whole comparison is moot) AND at least one mitigation actually clears it" -- both being true is
   // the direct, minimal evidence that (a) the reproduction still holds and (b) a real fix exists.
   // [Claude, 2026-10-01] since the SAT contact-point fix (detail::boxBoxContactPoint) the baseline no longer injects energy
   // here, so "baseline reproduces the defect" no longer holds; the criterion is now that NO configuration injects
   // energy (with ANTIDUNES_SAT_FORMER_POINT set, the former point and the former outcome are reproduced).
   const bool pass = configsWithoutInjection == uint_t(sizeof(configs) / sizeof(configs[0]));
   (void) anyFixPassed;
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison: baseline reproduces the defect="
      << baselineFailed << ", at least one mitigation fixes it=" << anyFixPassed << "  "
      << (pass ? "PASS" : "FAIL"));
   return pass;
}

// Targets a DIFFERENT regime than runCornerNormalDiscontinuity/runCornerFixComparison/
// runMultiPointBoxContactFuzz above: all three of those start the two boxes SEPARATED and let them
// approach at a slow, realistic speed, so contact only ever begins at a SHALLOW, near-degenerate
// overlap -- exactly the regime Fraige et al. 2008a's multi-point concept targets. None of them ever
// exercise what boxBoxMultiPointContacts (PaperAnalyticalContacts.h) does when a vertex is embedded
// DEEP inside the other box -- e.g. a jam where several floaters have already compressed hard against
// each other, or a body released directly into an existing pile (see the RELEASE-OVERLAP diagnostic
// in ParticleKernels.h). That function's per-vertex overlap is bounded only by
// min(halfExtents) (up to ~1 LU for this project's own <6,1,1.5> floating boxes, i.e. NOT a small
// near-surface number) and has no shared cap across points the way a single SAT overlap implicitly
// does -- each embedded vertex contributes its OWN independent spring force at its OWN world point,
// summed. The concern (raised directly against real production VTK evidence: uid 3314649325744685089
// tipping to a near-vertical pose over ~9000 steps in job 4201127, the useMultiPointBoxContact rerun of
// the crashed 200k "boxes/layer1only" scenario) is that under deep interpenetration this can inject
// MORE spurious energy/torque than the single-point boxBoxSAT path it's meant to improve on, since N
// independently-deep vertices each store their own "spring energy" rather than one shared overlap
// value.
//
// Starts the two boxes ALREADY OVERLAPPING by a controlled amount (approachSpeed=0 -- no closing
// momentum at all, so ANY kinetic energy that appears after release comes purely from resolving the
// static initial overlap, not from a real collision) at several depths from shallow (a bare touch) to
// severe (well past what boxBoxSAT's own single-axis result could plausibly represent), each run
// against baseline/rounded/multiPoint/both -- mirroring runCornerFixComparison's own 4-way structure,
// but swept over depth instead of compared at one fixed shallow contact. Same strict e=1, mu=0, no-
// gravity ECC bound as every other scenario in this file: total system KE must never exceed its own
// pre-release value, since there is no other energy source. A tilt is kept on both bodies (same
// defaults as CornerScenarioParams) specifically so MULTIPLE vertices, not just one clean face-face
// overlap, are actually embedded -- a flat face-face deep overlap wouldn't exercise the multi-point-
// vs-single-point comparison this test exists for at all.
inline bool runMultiPointDeepOverlapStability(const config::Config::BlockHandle& b)
{
   detail::CornerScenarioParams p;
   p.halfExtents      = b.getParameter< Vec3 >("halfExtents", p.halfExtents);
   p.density          = b.getParameter< real_t >("density", p.density);
   p.restitution      = b.getParameter< real_t >("restitutionCoefficient", p.restitution);
   p.friction         = b.getParameter< real_t >("frictionCoefficient", p.friction);
   p.poisson          = b.getParameter< real_t >("poissonsRatio", p.poisson);
   p.collisionTime    = b.getParameter< real_t >("collisionTime", p.collisionTime);
   p.dt               = b.getParameter< real_t >("dt", p.dt);
   p.approachSpeed    = real_t(0); // static overlap only -- no closing momentum, see header comment
   p.lateralOffsetY   = b.getParameter< real_t >("lateralOffsetY", p.lateralOffsetY);
   p.lateralOffsetZ   = b.getParameter< real_t >("lateralOffsetZ", p.lateralOffsetZ);
   p.tiltDegreesA     = b.getParameter< real_t >("tiltDegreesA", p.tiltDegreesA);
   p.tiltDegreesB     = b.getParameter< real_t >("tiltDegreesB", p.tiltDegreesB);
   p.spinRateB        = real_t(0); // isolate overlap-release energy from any imposed spin too
   p.tolerance        = b.getParameter< real_t >("toleranceRelError", p.tolerance);
   p.maxSteps         = b.getParameter< uint_t >("maxSteps", uint_t(4000));
   p.vtkSpacing       = b.getParameter< uint_t >("vtkSpacing", p.vtkSpacing);
   p.boxRoundingRadius = b.getParameter< real_t >("boxRoundingRadius", p.boxRoundingRadius);

   // Separation values relative to the box's own full length (2*halfExtents[0] = 12 by default):
   // 10 = shallow corner touch (same order as runCornerNormalDiscontinuity's own 20 minus its own
   // approach distance), down to -4 = centres nearly coincident, well past anything a real single
   // sub-cycle's worth of approach velocity could produce -- deliberately including an extreme point
   // so a latent "unbounded with depth" defect cannot hide just outside a narrower sweep.
   const std::vector< real_t > depths = { real_t(10), real_t(6), real_t(3), real_t(0), real_t(-4) };

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
   };
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false },
      { "both", true, true, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false },
      { "useAxisHysteresis", false, false, false, true, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapStability -- same tilt/offset "
      "geometry as runCornerNormalDiscontinuity, started ALREADY OVERLAPPING (approachSpeed=0) at "
      << depths.size() << " depths x " << (sizeof(configs) / sizeof(configs[0])) << " configurations:");

   bool multiPointEverWorseThanBaseline = false;
   bool anyConfigFailed                 = false;
   for (const real_t sep : depths)
   {
      p.separation = sep;
      real_t baselineGain = real_t(0);
      for (const auto& cfg : configs)
      {
         p.useRoundedBoxCorners    = cfg.rounded;
         p.useMultiPointBoxContact = cfg.multiPoint;
         p.useBoxOverlapCap        = cfg.overlapCap;
         p.useAxisHysteresis       = cfg.hysteresis;
         p.useContactManifold      = cfg.manifold;
         p.outputBaseName = std::string("multipoint_deep_overlap_sep") + std::to_string(sep) + "_" + cfg.name;

         const auto r = detail::runCornerScenarioOnce(p);
         if (std::string(cfg.name) == "baseline (both OFF)") baselineGain = r.keGainRel;
         if (std::string(cfg.name) == "useMultiPointBoxContact" && r.keGainRel > baselineGain + p.tolerance)
            multiPointEverWorseThanBaseline = true;
         if (!r.pass) anyConfigFailed = true;

         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapStability  separation=" << sep
            << " config=\"" << cfg.name << "\"  KE_before=" << r.keBefore << " PE_initial=" << r.peInitial
            << " KE_max=" << r.keMax
            << " (step " << r.keMaxStep << ") maxRelGain=" << r.keGainRel
            << " rawSatNormalFlipCount=" << r.normalFlipCount
            << " cumulativeNormalRotationDeg=" << r.cumulativeNormalRotationDeg
            << " maxPointJumpFrac=" << r.maxPointJumpFrac << " pointJumpCount=" << r.pointJumpCount
            << " everContacted=" << r.everContacted << "  " << (r.pass ? "PASS" : "FAIL"));
      }
   }

   // Two independent findings reported separately: (a) does EVERY configuration stay inside the
   // strict ECC bound at every depth (the same bar every other scenario in this file uses), and (b)
   // even where multiPoint itself still passes that bound, does it ever inject MORE relative energy
   // than the plain single-point baseline at the SAME depth -- the specific "multi-point is not an
   // improvement, it behaves more like the original defect" concern this test was written to check,
   // which a bare PASS/FAIL against the tolerance alone would not surface.
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapStability OVERALL: "
      "anyConfigFailedEccBound=" << anyConfigFailed << "  multiPointEverWorseThanBaselineAtSameDepth="
      << multiPointEverWorseThanBaseline << "  "
      << (!anyConfigFailed && !multiPointEverWorseThanBaseline ? "PASS" : "FAIL"));
   return !anyConfigFailed && !multiPointEverWorseThanBaseline;
}

// Randomized counterpart to runMultiPointDeepOverlapStability, the same relationship
// runMultiPointBoxContactFuzz already has to runCornerNormalDiscontinuity/runCornerFixComparison: the
// hand-picked depths/tilt/offset in that scenario are FIVE specific points on what is actually a
// continuous, multi-dimensional space (box aspect ratio, tilt on each body, lateral offset, starting
// overlap depth) -- this sweeps that space randomly, looking for the worst-case combination rather
// than trusting five samples to be representative. Every trial:
//   (1) GEOMETRY, direct on paper_analytical::boxBoxMultiPointContacts (same sanity checks as
//       runMultiPointBoxContactFuzz's own geometry pass: finite, unit normals, non-negative overlap)
//       PLUS the new check this whole investigation exists for: records the WORST (largest) single
//       vertex overlap seen across all trials, and separately counts how many trials produce at least
//       one vertex overlap deep enough to trip [LARGE-OVERLAP-MULTIPOINT]'s own threshold (half the
//       smaller body's bounding-sphere radius, ParticleKernels.h) -- quantifying how often a "normal"
//       random configuration (not a hand-crafted worst case) actually reaches that regime.
//   (2) DYNAMICS, all 4 configurations (baseline/rounded/multiPoint/both) per trial, approachSpeed=0
//       (static release only, isolating overlap-driven injection from momentum transfer -- same
//       rationale as runMultiPointDeepOverlapStability), strict e=1/mu=0/no-gravity ECC bound. Tracks
//       whether multiPoint (or both) EVER injects more relative energy than the plain baseline at the
//       SAME randomly-drawn trial -- the specific "not an improvement, more like the original defect"
//       question this whole line of testing exists to answer, which a bare pass/fail vs. a fixed
//       tolerance would not surface on its own.
inline bool runMultiPointDeepOverlapFuzz(const config::Config::BlockHandle& b)
{
   const uint_t numTrials        = b.getParameter< uint_t >("numTrials", uint_t(300));
   const real_t tolerance        = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t maxStepsPerTrial = b.getParameter< uint_t >("maxStepsPerTrial", uint_t(4000));
   const uint_t randomSeed       = b.getParameter< uint_t >("randomSeed", uint_t(1337));

   math::seedRandomGenerator(randomSeed);

   uint_t geometryChecked        = uint_t(0);
   uint_t geometryNaNOrInfCount  = uint_t(0);
   uint_t geometryBadNormalCount = uint_t(0);
   uint_t geometryNegativeOverlapCount = uint_t(0);
   real_t worstVertexOverlap     = real_t(0);
   real_t worstVertexOverlapHalfExtentMin = real_t(0); // min(halfExtents) for the worst trial, for context
   uint_t trialsExceedingLargeOverlapThreshold = uint_t(0); // mirrors [LARGE-OVERLAP-MULTIPOINT]'s own 0.5x bound

   uint_t contactedCount           = uint_t(0);
   uint_t baselineFailCount        = uint_t(0);
   uint_t multiPointFailCount      = uint_t(0);
   uint_t multiPointWorseThanBaselineCount = uint_t(0);
   real_t worstMultiPointGain      = real_t(0);
   uint_t worstMultiPointGainTrial = uint_t(0);
   uint_t cappedFailCount               = uint_t(0);
   uint_t cappedWorseThanBaselineCount  = uint_t(0);
   uint_t hysteresisFailCount               = uint_t(0);
   uint_t hysteresisWorseThanBaselineCount  = uint_t(0);
   uint_t manifoldFailCount                 = uint_t(0);
   uint_t manifoldWorseThanBaselineCount    = uint_t(0);

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
   };
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false },
      { "both", true, true, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false },
      { "useAxisHysteresis", false, false, false, true, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true },
   };

   for (uint_t trial = 0; trial < numTrials; ++trial)
   {
      const Vec3   halfExtents(math::realRandom< real_t >(real_t(0.5), real_t(8)),
                                math::realRandom< real_t >(real_t(0.5), real_t(4)),
                                math::realRandom< real_t >(real_t(0.5), real_t(4)));
      const real_t tiltA = math::realRandom< real_t >(real_t(-90), real_t(90));
      const real_t tiltB = math::realRandom< real_t >(real_t(-90), real_t(90));
      const real_t offY  = math::realRandom< real_t >(-halfExtents[1], halfExtents[1]);
      const real_t offZ  = math::realRandom< real_t >(-halfExtents[2], halfExtents[2]);
      // Full range from clearly separated-but-close through severely embedded (centres closer than
      // either body's own smallest half-extent) -- the added dimension vs. runMultiPointBoxContactFuzz,
      // which never starts anywhere but a fresh approach from clear separation.
      const real_t maxHalfLength = halfExtents[0];
      const real_t separation    = math::realRandom< real_t >(-real_t(1.5) * maxHalfLength, real_t(1.2) * maxHalfLength);

      // (1) Geometry sanity + worst-vertex-overlap tracking, direct on the vertex-test function, at
      // the exact pose this trial's dynamics run will start from (posA/posB/rotA/rotB match
      // detail::runCornerScenarioOnce's own construction below).
      {
         const Vec3 posA(-separation / real_t(2), real_t(0), real_t(0));
         const Vec3 posB(separation / real_t(2), offY, offZ);
         const Rot3 rotA(Vec3(0, 0, 1), tiltA * math::pi / real_t(180));
         const Rot3 rotB(Vec3(0, 0, 1), tiltB * math::pi / real_t(180));

         ++geometryChecked;
         const auto points = paper_analytical::boxBoxMultiPointContacts(posA, rotA.getMatrix(), halfExtents,
                                                                          posB, rotB.getMatrix(), halfExtents);
         const real_t smallerRad = real_t(2) * halfExtents.length() * real_t(0.5); // bounding-sphere radius (both bodies identical this trial)
         bool anyExceedsThreshold = false;
         for (const auto& pt : points)
         {
            const bool finite = std::isfinite(pt.point[0]) && std::isfinite(pt.point[1]) &&
                                 std::isfinite(pt.point[2]) && std::isfinite(pt.normal[0]) &&
                                 std::isfinite(pt.normal[1]) && std::isfinite(pt.normal[2]) &&
                                 std::isfinite(pt.overlap);
            if (!finite) { ++geometryNaNOrInfCount; continue; }
            if (pt.overlap < real_t(-1e-9)) ++geometryNegativeOverlapCount;
            if (std::abs(pt.normal.length() - real_t(1)) > real_t(1e-6)) ++geometryBadNormalCount;
            if (pt.overlap > worstVertexOverlap)
            {
               worstVertexOverlap = pt.overlap;
               worstVertexOverlapHalfExtentMin = std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
            }
            if (pt.overlap > real_t(0.5) * smallerRad) anyExceedsThreshold = true;
         }
         if (anyExceedsThreshold) ++trialsExceedingLargeOverlapThreshold;
      }

      // (2) Full dynamics, all 4 configurations, same randomly-drawn geometry.
      detail::CornerScenarioParams p;
      p.halfExtents             = halfExtents;
      p.density                 = real_t(0.6);
      p.restitution             = real_t(1);
      p.friction                = real_t(0);
      p.poisson                 = real_t(0.22);
      p.collisionTime           = real_t(32.0078);
      p.dt                      = real_t(0.05);
      p.approachSpeed           = real_t(0); // static release only, see header comment
      p.lateralOffsetY          = offY;
      p.lateralOffsetZ          = offZ;
      p.tiltDegreesA            = tiltA;
      p.tiltDegreesB            = tiltB;
      p.spinRateB               = real_t(0);
      p.separation              = separation;
      p.tolerance               = tolerance;
      p.maxSteps                = maxStepsPerTrial;
      p.enableOutput            = false;

      real_t baselineGain = real_t(0);
      bool   everContactedThisTrial = false;
      for (const auto& cfg : configs)
      {
         p.useRoundedBoxCorners    = cfg.rounded;
         p.useMultiPointBoxContact = cfg.multiPoint;
         p.useBoxOverlapCap        = cfg.overlapCap;
         p.useAxisHysteresis       = cfg.hysteresis;
         p.useContactManifold      = cfg.manifold;

         const auto r = detail::runCornerScenarioOnce(p);
         everContactedThisTrial |= r.everContacted;
         if (!r.everContacted) continue;

         if (std::string(cfg.name) == "baseline (both OFF)")
         {
            baselineGain = r.keGainRel;
            if (!r.pass) ++baselineFailCount;
         }
         if (std::string(cfg.name) == "useMultiPointBoxContact")
         {
            if (!r.pass) ++multiPointFailCount;
            if (r.keGainRel > baselineGain + tolerance)
            {
               ++multiPointWorseThanBaselineCount;
               if (r.keGainRel > worstMultiPointGain)
               {
                  worstMultiPointGain      = r.keGainRel;
                  worstMultiPointGainTrial = trial;
               }
               static uint_t worseLogCount = uint_t(0);
               if (worseLogCount < uint_t(20))
               {
                  WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz WORSE-THAN-BASELINE "
                     "TRIAL " << trial << ": halfExtents=" << halfExtents << " tiltA=" << tiltA
                     << " tiltB=" << tiltB << " offY=" << offY << " offZ=" << offZ
                     << " separation=" << separation << " | baselineGain=" << baselineGain
                     << " multiPointGain=" << r.keGainRel);
                  ++worseLogCount;
               }
            }
         }
         if (std::string(cfg.name) == "useMultiPointBoxContact+overlapCap")
         {
            if (!r.pass) ++cappedFailCount;
            if (r.keGainRel > baselineGain + tolerance) ++cappedWorseThanBaselineCount;
         }
         if (std::string(cfg.name) == "useAxisHysteresis")
         {
            if (!r.pass) ++hysteresisFailCount;
            if (r.keGainRel > baselineGain + tolerance) ++hysteresisWorseThanBaselineCount;
         }
         if (std::string(cfg.name) == "useMultiPointBoxContact+manifold")
         {
            if (!r.pass) ++manifoldFailCount;
            if (r.keGainRel > baselineGain + tolerance) ++manifoldWorseThanBaselineCount;
         }
      }
      if (everContactedThisTrial) ++contactedCount;
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  GEOMETRY: " << geometryChecked
      << " random configurations | NaN/Inf=" << geometryNaNOrInfCount << " | negative overlap="
      << geometryNegativeOverlapCount << " | non-unit normal=" << geometryBadNormalCount
      << " | worst single-vertex overlap=" << worstVertexOverlap << " (min(halfExtents) that trial="
      << worstVertexOverlapHalfExtentMin << ") | trials with >=1 vertex past [LARGE-OVERLAP-MULTIPOINT]'s "
         "own threshold=" << trialsExceedingLargeOverlapThreshold << "/" << numTrials);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  DYNAMICS: " << numTrials
      << " randomized trials | contacted=" << contactedCount << " | baseline ECC-bound failures="
      << baselineFailCount << " | multiPoint ECC-bound failures=" << multiPointFailCount
      << " | multiPoint WORSE than baseline at same trial=" << multiPointWorseThanBaselineCount
      << " | worst multiPoint relGain=" << worstMultiPointGain << " (trial " << worstMultiPointGainTrial << ")");
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  CAPPED: multiPoint+overlapCap ECC-bound "
      "failures=" << cappedFailCount << " | WORSE than baseline at same trial=" << cappedWorseThanBaselineCount);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  HYSTERESIS: useAxisHysteresis ECC-bound "
      "failures=" << hysteresisFailCount << " | WORSE than baseline at same trial=" << hysteresisWorseThanBaselineCount);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  MANIFOLD: multiPoint+manifold ECC-bound "
      "failures=" << manifoldFailCount << " | WORSE than baseline at same trial=" << manifoldWorseThanBaselineCount);

   // A geometry violation (NaN/non-unit-normal/negative overlap) is a hard bug regardless of depth --
   // always fails. Deep-overlap ECC-bound failures and the worse-than-baseline count are NOT gated
   // into pass/fail here (unlike runMultiPointBoxContactFuzz's shallow-regime equivalent): this test
   // exists specifically to CHARACTERIZE how often and how badly the known deep-overlap gap occurs, not
   // to assert it doesn't -- see multiPointDeepOverlapStability (which already established the defect
   // is real) for the strict pass/fail gate. Report-only for the dynamics half.
   const bool geometryOk = geometryNaNOrInfCount == uint_t(0) && geometryNegativeOverlapCount == uint_t(0) &&
                            geometryBadNormalCount == uint_t(0);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz OVERALL: geometryOk=" << geometryOk
      << "  " << (geometryOk ? "PASS" : "FAIL"));
   return geometryOk;
}

// Hand-picked, named counterpart to the two randomized sweeps above: specific geometric configurations
// chosen because each is individually a plausible worst case for a DIFFERENT reason, not because a
// random search happened to land there. A fuzz sweep tells you a bad case exists somewhere in the
// space; these tell you WHICH kind of geometry to worry about, by name, every time this test runs --
// useful as a fixed regression set (adding a new mitigation later, its effect on each named case here
// is directly legible, unlike a fuzz seed that may or may not still land near a similar case) and
// because a couple of these (coincident, engulfed) sit at exact geometric singularities a random draw
// would almost never hit precisely.
inline bool runMultiPointEdgeCases(const config::Config::BlockHandle& b)
{
   // DEFAULT 1 = the strict ECC diagnostic bound every OTHER scenario in this file uses. Overridable
   // to check whether production's own REAL damping (e=0.3) masks a failure this strict bound shows --
   // see runContactStiffnessCalibrationSweep's own finding that e=0.3 rescues the generic corner-tilt
   // geometry completely; this lets the SAME check run against each of this function's own named
   // cases, including "grazing-plank-stack" (the everyday geometry, and the worst failure found this
   // session under e=1: +618% for useMultiPointBoxContact alone).
   const real_t restitutionOverride = b.getParameter< real_t >("restitutionCoefficient", real_t(1));
   struct EdgeCase
   {
      const char* name;
      const char* why;
      Vec3        halfExtentsA;
      Vec3        halfExtentsB;
      real_t      tiltA, tiltB;   // degrees about Z, same convention as CornerScenarioParams
      real_t      offY, offZ;
      real_t      separation;     // X-axis centre-to-centre offset parameter (see runCornerScenarioOnce)
   };

   const std::vector< EdgeCase > cases = {
      { "exact-coincident", "two identical boxes at the SAME position/orientation -- overlap is maximal "
        "and every SAT axis is simultaneously tied (zero gap on all 3 face-pairs), the purest possible "
        "normal-direction singularity.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(0), real_t(0), real_t(0), real_t(0), real_t(0) },

      { "fully-engulfed", "a small box entirely inside a much larger one, centred -- every one of the "
        "small box's 8 vertices is embedded simultaneously, each with a DIFFERENT nearest-face axis, "
        "exercising the full N=8 branch of the per-point force loop at once.",
        Vec3(6, 4, 4), Vec3(1, 1, 1), real_t(0), real_t(0), real_t(0), real_t(0), real_t(0) },

      { "perpendicular-cross", "two identical planks with perpendicular long axes crossing through each "
        "other's centres (a 3D '+' shape) -- the most symmetric possible tie between two candidate "
        "separating axes, both equally valid.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(0), real_t(90), real_t(0), real_t(0), real_t(0) },

      { "corner-into-face-center", "one box's corner driven exactly into the geometric centre of the "
        "other's largest face -- a single vertex deeply embedded with no others nearby, isolating the "
        "single-deep-vertex failure mode directly (no averaging-out from multiple simultaneous points).",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(45), real_t(0), real_t(0), real_t(0), real_t(-5) },

      { "45deg-symmetric-tiebreak", "both bodies tilted +-45 degrees symmetrically about the line "
        "connecting their centres -- SAT's 15-axis test has no single best-separated axis by "
        "construction, a textbook near-degenerate configuration.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(45), real_t(-45), real_t(0), real_t(0), real_t(-2) },

      { "plank-vs-near-cube-deep", "this project's own production floating-box aspect ratio (12x2x3) "
        "deeply overlapping a nearly-cubic body of comparable volume -- the actual shape mismatch a real "
        "jam produces (planks piling against blockier debris), not two identical bodies.",
        Vec3(6, 1, 1.5), Vec3(2.5, 2.3, 2.2), real_t(20), real_t(-15), real_t(0.3), real_t(-0.2), real_t(-3) },

      { "flush-stack-axis-aligned", "both boxes perfectly axis-aligned, offset along X only, well past "
        "full-face overlap -- the SIMPLEST possible deep overlap (no tilt at all), a sanity baseline: if "
        "even this triggers the defect, tilt/offset were never the real trigger, overlap depth alone is.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(0), real_t(0), real_t(0), real_t(0), real_t(-4) },

      { "extreme-size-ratio-deep", "a tiny box (comparable to boxRoundingRadius's own scale) deeply "
        "embedded in a very large one -- checks the per-vertex overlap bound (or lack of one) doesn't "
        "behave pathologically differently when the two bodies' own sizes are very different, unlike "
        "every other case here which uses comparable-sized bodies.",
        Vec3(6, 4, 4), Vec3(0.5, 0.5, 0.5), real_t(30), real_t(0), real_t(0), real_t(0), real_t(0) },

      { "grazing-plank-stack", "two production-shaped planks nearly parallel (5 degree relative tilt), "
        "offset laterally by almost a full half-extent -- mimics two real floaters stacked slightly "
        "askew in a pile, the everyday case rather than a constructed extreme.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(2), real_t(-3), real_t(0.9), real_t(1.3), real_t(-1) },
   };

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
      bool        rotationalEffMass = false;
   };
   // [Claude, 2026-10-29] two "useSpheroBoxContact (GJK/Eberly analytical)" rows removed -- see
   // archive/unused_contact_methods/.
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false, false },
      { "both", true, true, false, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false, false },
      { "useAxisHysteresis", false, false, false, true, false, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true, false },
      { "useRotationalEffectiveMass", false, false, false, false, false, true },
      { "useMultiPointBoxContact+rotationalEffMass", false, true, false, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases -- " << cases.size()
      << " hand-picked geometric configurations x " << (sizeof(configs) / sizeof(configs[0]))
      << " configurations, approachSpeed=0 (static release, same rationale as "
         "multiPointDeepOverlapStability):");

   uint_t multiPointWorseThanBaselineCount = uint_t(0);
   for (const auto& ec : cases)
   {
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases  case=\"" << ec.name << "\" -- " << ec.why);

      detail::CornerScenarioParams p;
      // halfExtentsA/B differ in general -- CornerScenarioParams only has one halfExtents field (both
      // bodies share it in every other scenario in this file), so the smaller/asymmetric-size cases
      // here are built by scaling B's EFFECTIVE size via a direct engine build instead of going through
      // detail::runCornerScenarioOnce (which assumes equal sizes). Kept self-contained rather than
      // extending CornerScenarioParams for a one-off need used only by 3 of these 9 cases.
      const bool sameSize = (ec.halfExtentsA == ec.halfExtentsB);

      real_t baselineGain = real_t(0);
      for (const auto& cfg : configs)
      {
         detail::CornerScenarioResult r;
         if (sameSize)
         {
            p.halfExtents             = ec.halfExtentsA;
            p.density                 = real_t(0.6);
            p.restitution             = restitutionOverride;
            p.friction                = real_t(0);
            p.poisson                 = real_t(0.22);
            p.collisionTime           = real_t(32.0078);
            p.dt                      = real_t(0.05);
            p.approachSpeed           = real_t(0);
            p.lateralOffsetY          = ec.offY;
            p.lateralOffsetZ          = ec.offZ;
            p.tiltDegreesA            = ec.tiltA;
            p.tiltDegreesB            = ec.tiltB;
            p.spinRateB               = real_t(0);
            p.separation              = ec.separation;
            p.tolerance               = real_t(0.02);
            p.maxSteps                = uint_t(6000);
            p.enableOutput            = false;
            p.useRoundedBoxCorners    = cfg.rounded;
            p.boxRoundingRadius       = real_t(0.05);
            p.useMultiPointBoxContact = cfg.multiPoint;
            p.useBoxOverlapCap        = cfg.overlapCap;
            p.useAxisHysteresis       = cfg.hysteresis;
            p.useContactManifold      = cfg.manifold;
            p.useRotationalEffectiveMass = cfg.rotationalEffMass;
            r = detail::runCornerScenarioOnce(p);
         }
         else
         {
            // Asymmetric-size cases: build the engine directly rather than through
            // runCornerScenarioOnce, which only accepts one shared halfExtents for both bodies.
            DEMEngine engine(real_t(0.05), restitutionOverride, real_t(0), real_t(0.22), real_t(32.0078), real_t(0),
                              real_t(200), /*useAnalyticalContactDetection=*/true, cfg.rounded, real_t(0.05),
                              cfg.multiPoint, cfg.overlapCap, real_t(0.3), cfg.hysteresis, real_t(0.15),
                              cfg.manifold, real_t(0.25), cfg.rotationalEffMass);
            const Rot3 rotA(Vec3(0, 0, 1), ec.tiltA * math::pi / real_t(180));
            const Rot3 rotB(Vec3(0, 0, 1), ec.tiltB * math::pi / real_t(180));
            const size_t idxA = engine.addBox(Vec3(-ec.separation / 2, 0, 0), Vec3(0, 0, 0), rotA,
                                               real_t(2) * ec.halfExtentsA, real_t(0.6));
            const size_t idxB = engine.addBox(Vec3(ec.separation / 2, ec.offY, ec.offZ), Vec3(0, 0, 0), rotB,
                                               real_t(2) * ec.halfExtentsB, real_t(0.6));
            engine.finalizeSetup();

            r.keBefore  = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
            r.peInitial = detail::computeInitialElasticPE(
               Vec3(-ec.separation / 2, 0, 0), rotA, ec.halfExtentsA, real_t(0.6),
               Vec3(ec.separation / 2, ec.offY, ec.offZ), rotB, ec.halfExtentsB, real_t(0.6),
               restitutionOverride, real_t(32.0078), cfg.rounded, real_t(0.05), cfg.multiPoint, cfg.overlapCap,
               real_t(0.3), cfg.rotationalEffMass);
            r.keMax    = r.keBefore;
            for (uint_t s = 0; s < uint_t(6000); ++s)
            {
               real_t penetration = real_t(0);
               if (engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0)) r.everContacted = true;
               engine.step();
               const real_t keNow = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
               if (keNow > r.keMax) { r.keMax = keNow; r.keMaxStep = s; }
            }
            r.keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
            const real_t ceiling = r.keBefore + r.peInitial;
            r.keGainRel = (r.keMax - ceiling) / std::max(real_t(1e-12), ceiling);
            r.pass      = r.everContacted && r.keGainRel < real_t(0.02);
         }

         if (std::string(cfg.name) == "baseline (both OFF)") baselineGain = r.keGainRel;
         if (std::string(cfg.name) == "useMultiPointBoxContact" && r.keGainRel > baselineGain + real_t(0.02))
            ++multiPointWorseThanBaselineCount;

         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases    config=\"" << cfg.name << "\""
            << "  KE_before=" << r.keBefore << " PE_initial=" << r.peInitial << " KE_max=" << r.keMax
            << " (step " << r.keMaxStep << ")"
            << " maxRelGain=" << r.keGainRel
            << " rawSatNormalFlipCount=" << r.normalFlipCount
            << " cumulativeNormalRotationDeg=" << r.cumulativeNormalRotationDeg
            << " maxPointJumpFrac=" << r.maxPointJumpFrac << " pointJumpCount=" << r.pointJumpCount
            << (sameSize ? "" : " (N/A, asymmetric-size branch)")
            << " everContacted=" << r.everContacted
            << "  " << (r.pass ? "PASS" : "FAIL"));
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases OVERALL: multiPointWorseThanBaseline in "
      << multiPointWorseThanBaselineCount << "/" << cases.size() << " named cases");
   // Characterization tool like multiPointDeepOverlapFuzz -- report-only, not a hard gate (some of
   // these cases are deliberately singular/degenerate and not expected to cleanly pass any config).
   return true;
}

// multiPointDeepOverlapStability's own dt-independence result (16x timestep refinement, <0.01%
// change in baseline's gain at every depth) rules out numerical/integration error as the cause of
// baseline's OWN deep-overlap energy injection (separations 6/3/0/-4 all fail with ZERO detected
// axis flips, near-zero cumulative normal rotation, and a contact-point jump that itself shrinks to
// noise as dt refines -- yet the energy gain does not move at all). None of useRoundedBoxCorners/
// useMultiPointBoxContact/useBoxOverlapCap/useAxisHysteresis/useContactManifold change this, because
// all five still compute their force from the SAME mesa_pd::kernel::LinearSpringDashpot::
// setStiffnessAndDamping calibration (stiffnessN/dampingN), which uses ONLY the translational
// effective mass (meff = 1/(1/mA + 1/mB)) -- no correction for the ROTATIONAL inertia excited by a
// contact point offset from either body's centroid, which every case in this file has by
// construction (that offset is exactly what makes the near-degenerate/multi-point geometry
// interesting in the first place). At restitution=1 (this file's own strict diagnostic bound) that
// calibration gives EXACTLY ZERO damping -- correct for a head-on, non-rotating collision, but with
// no mechanism to correctly balance a contact that also does rotational work.
//
// This sweeps the EXISTING parameters that calibration already exposes, to see whether the effect is
// masked (not necessarily fixed) by realistic production values, before concluding a genuine formula
// change (e.g. an effective mass that also accounts for rotational inertia at the contact point,
// 1/meff = 1/mA + 1/mB + (r x n)^T I^-1 (r x n) summed appropriately) is actually required:
//   (1) restitutionCoefficient: e<1 introduces REAL damping (dampingN != 0) -- does the resulting
//       dissipation happen to offset the rotational-energy leak enough to mask it? Production itself
//       uses e=0.3 ("wet/waterlogged wood"), never e=1 -- this checks whether the strict e=1 bound
//       this whole file uses is somehow overstating a problem production's own real damping already
//       absorbs.
//   (2) collisionTime: controls the ABSOLUTE stiffness magnitude for a fixed restitution -- a softer
//       spring (larger collisionTime) produces less peak force for the same overlap, hence less
//       torque at a fixed lever arm; does simply detuning the spring change the picture?
//   (3) frictionCoefficient: acts on the TANGENTIAL component only (this test's own geometry has
//       nonzero relative tangential velocity at the contact point from p.spinRateB) -- tangential
//       dissipation could not fix a normal-direction energy leak directly, but is included as a
//       direct empirical check rather than an assumption.
inline bool runContactStiffnessCalibrationSweep(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents = b.getParameter< Vec3 >("halfExtents", Vec3(6, 1, 1.5));
   const real_t tolerance   = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t maxSteps    = b.getParameter< uint_t >("maxSteps", uint_t(4000));

   const std::vector< real_t > separations = { real_t(6), real_t(3), real_t(0), real_t(-4) };
   struct RestitutionCase
   {
      const char* name;
      real_t      value;
   };
   const std::vector< RestitutionCase > restitutions = {
      { "e=1.00 (strict diagnostic bound)", real_t(1.0) },
      { "e=0.99", real_t(0.99) },
      { "e=0.90", real_t(0.90) },
      { "e=0.50", real_t(0.50) },
      { "e=0.30 (PRODUCTION value)", real_t(0.30) },
   };
   struct CollisionTimeCase
   {
      const char* name;
      real_t      value;
   };
   const std::vector< CollisionTimeCase > collisionTimes = {
      { "collisionTime=8 (stiffer)", real_t(8) },
      { "collisionTime=32.0078 (PRODUCTION value)", real_t(32.0078) },
      { "collisionTime=128 (softer)", real_t(128) },
   };
   const std::vector< real_t > frictions = { real_t(0), real_t(0.5) }; // 0.5 = PRODUCTION value
   // Baseline alone isolates the stiffness calibration itself; ALSO sweeping useMultiPointBoxContact
   // answers the natural follow-up once baseline's own mystery gain turned out to vanish completely at
   // e=0.30 (production's own real value): does that SAME real damping also rescue multipoint's own,
   // much larger, "grazing-plank-stack" failure (+618% at e=1, DEMValidationScenarios.h's own
   // multiPointEdgeCases), or is multipoint's defect a genuinely different, damping-independent one?
   const std::vector< bool > multiPointOptions = { false, true };
   // useRotationalEffectiveMass: the third candidate fix researched this session (see
   // ModelParameters::useRotationalEffectiveMass's own doc comment, SimulationParams.h) -- unlike
   // useMultiPointBoxContact/useRoundedBoxCorners/useAxisHysteresis/useContactManifold (which all change
   // WHICH (point,normal,overlap) triplet feeds LinearSpringDashpot), this changes the force law's OWN
   // calibration: stiffnessN/dampingN are set from the generalized (translational+rotational) effective
   // mass at the actual contact point instead of the plain 1/(1/mA+1/mB). This scenario's own tiltDegreesA/
   // tiltDegreesB/lateralOffsetY/Z/spinRateB (see CornerScenarioParams defaults) give every contact here a
   // genuine moment arm off the line of centres -- precisely the geometry the rotational correction targets
   // -- so sweeping it crossed with useMultiPointBoxContact directly answers: does correcting the force
   // law's calibration alone (no change to which point/normal is used) already fix baseline's own mystery
   // gain, independent of and possibly BEFORE reaching for a different (point,normal) mitigation at all?
   const std::vector< bool > rotationalEffMassOptions = { false, true };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] contactStiffnessCalibrationSweep -- " << separations.size()
      << " deep-overlap depths x " << restitutions.size() << " restitution values x "
      << collisionTimes.size() << " collision times x " << frictions.size() << " friction values x "
      << multiPointOptions.size() << " useMultiPointBoxContact settings x "
      << rotationalEffMassOptions.size() << " useRotationalEffectiveMass settings, testing whether the "
         "LinearSpringDashpot stiffness/damping calibration itself (translational effective mass only, "
         "no rotational-inertia correction) is masked or fixed by any EXISTING parameter combination, "
         "for baseline AND for multipoint, with and without the rotational-effective-mass correction:");

   uint_t anyPassCount   = uint_t(0);
   uint_t totalRuns      = uint_t(0);
   real_t worstGainSeen  = real_t(0);
   std::string worstGainConfig;
   // Per-restitution/collisionTime/friction/multiPoint/rotationalEffMass combo: does EVERY depth pass
   // under it? This is the meaningful question -- a combo that only happens to pass at one depth isn't a
   // real fix.
   for (const bool useRotationalEffMass : rotationalEffMassOptions)
   for (const bool useMultiPoint : multiPointOptions)
   for (const auto& restCase : restitutions)
   {
      for (const auto& ctCase : collisionTimes)
      {
         for (const real_t friction : frictions)
         {
            bool allDepthsPassed = true;
            real_t worstGainThisCombo = real_t(0);
            for (const real_t sep : separations)
            {
               detail::CornerScenarioParams p;
               p.halfExtents   = halfExtents;
               p.density       = real_t(0.6);
               p.restitution   = restCase.value;
               p.friction      = friction;
               p.poisson       = real_t(0.22);
               p.collisionTime = ctCase.value;
               p.dt            = real_t(0.05);
               p.approachSpeed = real_t(0); // static release, same rationale as multiPointDeepOverlapStability
               p.separation    = sep;
               p.tolerance     = tolerance;
               p.maxSteps      = maxSteps;
               p.enableOutput  = false;
               p.useRoundedBoxCorners        = false;
               p.useMultiPointBoxContact     = useMultiPoint;
               p.useBoxOverlapCap            = false;
               p.useAxisHysteresis           = false;
               p.useContactManifold          = false;
               p.useRotationalEffectiveMass  = useRotationalEffMass;

               const auto r = detail::runCornerScenarioOnce(p);
               ++totalRuns;
               if (!r.pass) allDepthsPassed = false;
               if (r.keGainRel > worstGainThisCombo) worstGainThisCombo = r.keGainRel;
               if (r.keGainRel > worstGainSeen)
               {
                  worstGainSeen = r.keGainRel;
                  worstGainConfig = std::string("useRotationalEffectiveMass=") +
                                    (useRotationalEffMass ? "true" : "false") +
                                    " / useMultiPointBoxContact=" + (useMultiPoint ? "true" : "false") +
                                    " / " + restCase.name + " / " + ctCase.name + " / friction=" +
                                    std::to_string(friction) + " / separation=" + std::to_string(sep);
               }
            }
            if (allDepthsPassed) ++anyPassCount;
            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] contactStiffnessCalibrationSweep  "
               << "useRotationalEffectiveMass=" << useRotationalEffMass
               << " / useMultiPointBoxContact=" << useMultiPoint << " / " << restCase.name
               << " / " << ctCase.name << " / friction=" << friction << "  worstGainAcrossDepths="
               << worstGainThisCombo << "  " << (allDepthsPassed ? "ALL 4 DEPTHS PASS" : "at least one depth FAILS"));
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] contactStiffnessCalibrationSweep OVERALL: " << anyPassCount
      << "/" << (restitutions.size() * collisionTimes.size() * frictions.size() * multiPointOptions.size() *
                 rotationalEffMassOptions.size())
      << " parameter combinations pass at EVERY depth (" << totalRuns << " total runs) | worst gain seen="
      << worstGainSeen << " at [" << worstGainConfig << "]");
   // Report-only, like the other characterization scenarios in this file -- the point is to learn
   // WHICH (if any) existing parameter combination masks the effect, not to assert a pass/fail verdict.
   return anyPassCount > uint_t(0);
}

// Every scenario above starts the two bodies ALREADY overlapping at rest -- a controlled but somewhat
// artificial static-release probe of the contact FORMULA itself, isolated from how such an overlap
// would ever actually arise in production. This scenario instead reproduces the literal PHYSICAL
// process suspected of producing them: this project's own FloatingBox (halfExtents <6,1,1.5>,
// densityRatio 0.6) free-falling under REAL production gravity onto a small pile of IDENTICAL boxes
// already at rest on a floor -- exactly "a floater falling onto already-floating debris", the
// mechanism asked about directly. No water/buoyancy here (this whole file is a dry DEM-only harness,
// same as every other scenario) -- this isolates whether DEM contact resolution alone, given a
// realistic impact velocity from a real fall height, can launch something out of a pile, independent
// of any fluid effect.
//
// Two things are swept, both directly targeting "how/why do these get flying, especially with more
// floaters already present":
//   (1) pileSize: how many boxes are already resting side-by-side before the drop (1, 3, 6) -- tests
//       whether MORE simultaneous neighbours makes the outcome worse, i.e. whether this is fundamentally
//       a [SIMULTANEOUS-CONTACTS]-style effect (ParticleKernels.h's own diagnostic for exactly this
//       concern: "the EFFECTIVE combined stiffness ... may exceed what the explicit integrator's
//       timestep was tuned to handle for a single dominant contact").
//   (2) landingOffset: WHERE lateral the falling box lands relative to the pile -- dead-centre on one
//       pile box (clean, symmetric single contact), on the SEAM between two adjacent boxes (the
//       corner-into-two-neighbours case a real staggered/jittered pile produces most often), or
//       overhanging the pile's own edge (asymmetric, one side supported one side not).
//
// Energy bound here is NOT the elastic-PE-in-overlap accounting the scenarios above use (which is
// entangled with the contact model's own formula) -- it is the much more fundamental, model-
// independent one: total system mechanical energy E = sum(KE_i) + sum(m_i * g * z_i) can only
// DECREASE over time (restitution=0.3, "wet/waterlogged wood", matches production exactly -- real
// dissipation is expected, not just "no worse than start"). A rise in E at any step, by any amount
// beyond floating-point tolerance, is unambiguous evidence of energy injected from nowhere -- gravity
// is the only OTHER force in this dry harness and its own contribution is exact and trivial (mgz), so
// there is no formula-specific approximation error to account for here at all.
inline bool runFallOntoPileStability(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("halfExtents", Vec3(6, 1, 1.5));
   const real_t density       = b.getParameter< real_t >("density", real_t(0.6));
   const real_t restitution   = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t friction      = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(32.0078));
   const real_t dt            = b.getParameter< real_t >("dt", real_t(0.05));
   const real_t gravity       = b.getParameter< real_t >("gravity", real_t(7.02618e-6)); // production value
   const uint_t maxSteps      = b.getParameter< uint_t >("maxSteps", uint_t(8000));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const bool   writeVtk      = b.getParameter< bool >("writeVtk", false);
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(4));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "fall_onto_pile_stability");

   // CORRECTED MECHANISM (per direct feedback on an earlier version of this test): the RELEASE-OVERLAP
   // penetration is not a symptom of "a newly-released floater drops onto an already-settled pile" --
   // the just-released floater is FIXED (frozen, zero velocity) at its exact spawn point right up until
   // the instant of release, so at the very first contact check it hasn't moved AT ALL. The overlap
   // exists because a DIFFERENT, EARLIER-launched floater -- already flying due to the instability
   // this whole session investigates -- has drifted/climbed into that exact spot in the meantime (this
   // session's own VTK evidence: uid 1585267068834414616 climbed +16.7 LU in the 1000 steps before
   // detection, uid 5044031582654955533 climbed +10.6 LU). So the STATIONARY body here (the "pile")
   // is the correct stand-in for the just-released floater, and the ONE body given real velocity
   // (flyingSpeed below) is the stand-in for the already-launched neighbour arriving at it -- the
   // previous version of this test had these two roles backwards (both bodies started at rest).
   //
   // dropHeight is still the CLEAR GAP at t=0 (negative = starting already embedded) -- -0.71/-0.5
   // directly replay the two ACTUAL measured RELEASE-OVERLAP penetrations (job 4201127, step 122001);
   // 0 is the boundary case (already touching, zero gap).
   // quickMode (default false, unchanged behaviour): the full 3x2x2x1 geometry sweep x 9 configs = 108
   // runs is prohibitively slow when a config's deep initial overlap keeps [DEM-WATCH]/
   // [STRONG-CONTACT-FORCE] firing every single sub-cycle for the whole maxSteps window (each combo can
   // take several minutes just formatting/writing that log volume) -- true collapses the geometry sweep
   // to the SINGLE most representative combo (the actual measured release-overlap magnitude and speed,
   // dropHeight=-0.71, flyingSpeed=0.015, pileSize=6) so all 9 mitigation configs can still be compared
   // quickly, at the cost of not covering the other 11 geometry combinations this run.
   const bool quickMode = b.getParameter< bool >("quickMode", false);
   const std::vector< real_t > dropHeights =
      quickMode ? std::vector< real_t >{ real_t(-0.71) } : std::vector< real_t >{ real_t(-0.71), real_t(-0.5), real_t(0) };
   // Initial speed of the "flying" body (toward the stationary one, see the mechanism note above).
   // 0 isolates pure static-overlap release (no momentum at all, matching this file's OTHER deep-
   // overlap scenarios); 0.015 matches the ORDER OF MAGNITUDE this session directly measured for the
   // two real flying neighbours (+10.6/+16.7 LU over 1000 steps => ~0.011-0.017 LU/step average).
   const std::vector< real_t > flyingSpeeds =
      quickMode ? std::vector< real_t >{ real_t(0.015) } : std::vector< real_t >{ real_t(0), real_t(0.015) };

   const std::vector< uint_t > pileSizes = quickMode ? std::vector< uint_t >{ uint_t(6) } : std::vector< uint_t >{ uint_t(1), uint_t(6) };
   enum class Landing
   {
      OnCentre,
      OnSeam
   };
   const std::vector< std::pair< Landing, const char* > > landings = { { Landing::OnCentre, "on-centre" } };

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
      bool        rotationalEffMass = false;
   };
   // [Claude, 2026-10-29] "useSpheroBoxContact" row removed -- see archive/unused_contact_methods/.
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false, false },
      { "both", true, true, false, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false, false },
      { "useAxisHysteresis", false, false, false, true, false, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true, false },
      { "useRotationalEffectiveMass", false, false, false, false, false, true },
      { "useMultiPointBoxContact+rotationalEffMass", false, true, false, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] fallOntoPileStability -- a STATIONARY just-released "
      "FloatingBox vs. an ALREADY-FLYING neighbour under real production gravity (" << gravity << "), "
      << dropHeights.size() << " release gaps (negative = already overlapping) x " << flyingSpeeds.size()
      << " flying speeds x " << pileSizes.size() << " pile sizes x " << landings.size()
      << " landing offsets x " << (sizeof(configs) / sizeof(configs[0])) << " configurations:");

   bool anyConfigFailed              = false;
   bool multiPointEverWorseThanBase  = false;
   real_t worstRelGainOverall        = real_t(0);
   for (const real_t dropHeight : dropHeights)
   for (const real_t flyingSpeed : flyingSpeeds)
   for (const uint_t pileSize : pileSizes)
   {
      for (const auto& landing : landings)
      {
         real_t baselineGain = real_t(0);
         for (const auto& cfg : configs)
         {
            // Floor at z=0. Pile boxes rest side-by-side along Y, long axis (halfExtents[0]) along X,
            // flat face down (identity rotation) -- the resting orientation this project's own VTK
            // evidence shows these boxes actually settle into. Touching (zero gap) between neighbours,
            // matching an already-packed/jammed pile rather than a loosely-spaced one.
            const real_t pileZ    = halfExtents[2]; // box centre sits exactly halfExtents.z above the floor
            const real_t pileSpan = real_t(2) * halfExtents[1] * real_t(pileSize); // total Y extent of the pile
            const real_t pileY0   = -pileSpan / real_t(2) + halfExtents[1]; // first box's centre Y

            DEMEngine engine(dt, restitution, friction, poisson, collisionTime, gravity,
                              /*domainHalfExtent=*/real_t(500), /*useAnalyticalContactDetection=*/true,
                              cfg.rounded, real_t(0.05), cfg.multiPoint, cfg.overlapCap, real_t(0.3),
                              cfg.hysteresis, real_t(0.15), cfg.manifold, real_t(0.25),
                              cfg.rotationalEffMass);
            engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1)); // floor

            std::vector< size_t > pileIdx;
            for (uint_t i = 0; i < pileSize; ++i)
            {
               const real_t y = pileY0 + real_t(i) * real_t(2) * halfExtents[1];
               pileIdx.push_back(engine.addBox(Vec3(0, y, pileZ), Vec3(0, 0, 0), Rot3(),
                                                real_t(2) * halfExtents, density));
            }

            // Landing Y offset relative to the pile: dead centre of pile box 0, or the seam between
            // box 0 and box 1 (only meaningful for pileSize>=2 -- falls back to on-centre for
            // pileSize==1, still a valid single-neighbour baseline for that pile size).
            real_t landingY = pileY0;
            if (landing.first == Landing::OnSeam && pileSize >= 2) landingY = pileY0 + halfExtents[1];

            // This body is the STAND-IN FOR THE ALREADY-FLYING NEIGHBOUR (see the mechanism note
            // above) -- moving at flyingSpeed TOWARD the stationary pile, not the newly-released
            // floater itself (which is correctly represented by the stationary pile boxes below).
            const real_t dropZ = pileZ + halfExtents[2] + dropHeight + halfExtents[2];
            const size_t fallIdx = engine.addBox(Vec3(0, landingY, dropZ), Vec3(0, 0, -flyingSpeed), Rot3(),
                                                  real_t(2) * halfExtents, density);

            engine.finalizeSetup();

            if (writeVtk)
            {
               std::string safeName = cfg.name;
               for (char& c : safeName)
                  if (!std::isalnum(static_cast< unsigned char >(c))) c = '_';
               engine.setupVtk(outBase + "_" + safeName, vtkSpacing);
            }

            auto totalMechanicalEnergy = [&]() {
               real_t e = real_t(0);
               for (size_t idx : pileIdx) e += engine.kineticEnergy(idx) + engine.ac().getMass(idx) * gravity * engine.ac().getPosition(idx)[2];
               e += engine.kineticEnergy(fallIdx) + engine.ac().getMass(fallIdx) * gravity * engine.ac().getPosition(fallIdx)[2];
               return e;
            };

            const real_t e0     = totalMechanicalEnergy();
            real_t       eMax   = e0;
            uint_t       eMaxStep = uint_t(0);
            real_t       peakZAfterImpact = real_t(0);
            bool         everContacted    = false;
            for (uint_t s = 0; s < maxSteps; ++s)
            {
               real_t pen = real_t(0);
               if (engine.queryContact(fallIdx, pileIdx[0], pen) && pen < real_t(0)) everContacted = true;
               engine.step();
               if (writeVtk) engine.writeVtkIfDue();
               const real_t eNow = totalMechanicalEnergy();
               if (eNow > eMax) { eMax = eNow; eMaxStep = s; }
               if (everContacted)
               {
                  const real_t z = engine.ac().getPosition(fallIdx)[2];
                  if (z > peakZAfterImpact) peakZAfterImpact = z;
               }
            }

            const real_t relGain = (eMax - e0) / std::max(real_t(1e-12), std::abs(e0));
            const bool   pass    = relGain < tolerance;
            if (!pass) anyConfigFailed = true;
            if (std::string(cfg.name) == "baseline (both OFF)") baselineGain = relGain;
            if (std::string(cfg.name) == "useMultiPointBoxContact" && relGain > baselineGain + tolerance)
               multiPointEverWorseThanBase = true;
            if (relGain > worstRelGainOverall) worstRelGainOverall = relGain;

            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] fallOntoPileStability  dropHeight=" << dropHeight
               << " flyingSpeed=" << flyingSpeed << " pileSize=" << pileSize
               << " landing=" << landing.second << " config=\"" << cfg.name << "\"  E0=" << e0
               << " Emax=" << eMax << " (step " << eMaxStep << ") relGain=" << relGain
               << " peakZAfterImpact=" << peakZAfterImpact << " (dropZ was " << dropZ << ", pile top z="
               << (pileZ + halfExtents[2]) << ")  everContacted=" << everContacted
               << "  " << (pass ? "PASS" : "FAIL"));
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] fallOntoPileStability OVERALL: anyConfigFailed="
      << anyConfigFailed << "  multiPointEverWorseThanBaseline=" << multiPointEverWorseThanBase
      << "  worstRelGainSeen=" << worstRelGainOverall << "  " << (!anyConfigFailed ? "PASS" : "FAIL"));
   return !anyConfigFailed;
}

// Stress test for useMultiPointBoxContact before trusting it on a full production run: the hand-tuned
// geometry in runCornerNormalDiscontinuity/runCornerFixComparison is ONE specific corner/tilt/offset
// configuration -- this scenario sweeps many randomized ones (box aspect ratio, tilt about a random
// axis on EACH body, lateral offset, approach speed, spin) to check the mitigation holds broadly, not
// just for the one case it was validated against. Two independent checks per trial:
//   (1) GEOMETRY sanity, direct on paper_analytical::boxBoxMultiPointContacts: no NaN/Inf in any
//       returned point/normal/overlap, every normal is a unit vector, every overlap is >= 0, and --
//       the key cross-check against boxBoxSAT -- whenever SAT proves the pair definitely SEPARATED
//       (result.contact == false), the vertex-embedding test must ALSO find nothing (a genuinely
//       separated pair can never have one body's vertex sitting inside the other; a violation here
//       would mean a sign/axis bug in the vertex test itself, independent of any dynamics).
//   (2) DYNAMICS: the exact same strict ECC energy bound (e=1, mu=0, no gravity: total system KE must
//       never exceed its own pre-contact value) as runCornerNormalDiscontinuity, but run under
//       useMultiPointBoxContact=true across every randomized trial that actually made contact.
// A trial's box sizes/speeds/offsets are drawn from ranges that comfortably bracket this project's own
// production FloatingBox family (halfExtents as small as <6,1,1.5> and as chunky as roughly cube-
// shaped) and the realistic approach speeds this session's own forensics established (see
// runCornerNormalDiscontinuity's own header comment) -- not an arbitrary/unbounded fuzz domain.
inline bool runMultiPointBoxContactFuzz(const config::Config::BlockHandle& b)
{
   const uint_t numTrials         = b.getParameter< uint_t >("numTrials", uint_t(500));
   const real_t tolerance         = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t maxStepsPerTrial  = b.getParameter< uint_t >("maxStepsPerTrial", uint_t(8000));
   const uint_t randomSeed        = b.getParameter< uint_t >("randomSeed", uint_t(42));

   math::seedRandomGenerator(randomSeed);

   uint_t geometryChecked         = uint_t(0);
   uint_t geometryNaNOrInfCount   = uint_t(0);
   uint_t geometryBadNormalCount  = uint_t(0);
   uint_t geometryNegativeOverlapCount = uint_t(0);
   uint_t geometrySatDisagreementCount = uint_t(0); // SAT says separated, vertex test says embedded

   uint_t contactedCount = uint_t(0);
   uint_t energyPassCount = uint_t(0);
   real_t worstGainSeen   = real_t(0);
   uint_t worstGainTrial  = uint_t(0);

   for (uint_t trial = 0; trial < numTrials; ++trial)
   {
      const Vec3 halfExtents(math::realRandom< real_t >(real_t(0.5), real_t(8)),
                              math::realRandom< real_t >(real_t(0.5), real_t(4)),
                              math::realRandom< real_t >(real_t(0.5), real_t(4)));

      // (1) Pure-geometry sanity, several random relative poses per trial's box size -- independent of
      // any dynamics, so cheap to sample densely.
      for (uint_t g = 0; g < 5; ++g)
      {
         const Vec3 posA(real_t(0), real_t(0), real_t(0));
         const Rot3 rotA(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                               math::realRandom< real_t >(-1, 1))
                            .getNormalizedIfNotZero(),
                         math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
         // Offset chosen so roughly half the samples land in genuine overlap and half in clear
         // separation -- a pure "always separated" or "always overlapping" sweep would never exercise
         // the boundary this check cares about.
         const real_t maxReach = halfExtents.length() * real_t(2.2);
         const Vec3   posB(math::realRandom< real_t >(-maxReach, maxReach),
                            math::realRandom< real_t >(-maxReach, maxReach),
                            math::realRandom< real_t >(-maxReach, maxReach));
         const Rot3 rotB(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                               math::realRandom< real_t >(-1, 1))
                            .getNormalizedIfNotZero(),
                         math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));

         ++geometryChecked;
         const auto satResult = paper_analytical::boxBoxSAT(posA, rotA.getMatrix(), halfExtents, posB,
                                                              rotB.getMatrix(), halfExtents);
         const auto points = paper_analytical::boxBoxMultiPointContacts(posA, rotA.getMatrix(), halfExtents,
                                                                          posB, rotB.getMatrix(), halfExtents);

         if (!satResult.contact && !points.empty()) ++geometrySatDisagreementCount;

         for (const auto& pt : points)
         {
            const bool finite = std::isfinite(pt.point[0]) && std::isfinite(pt.point[1]) &&
                                 std::isfinite(pt.point[2]) && std::isfinite(pt.normal[0]) &&
                                 std::isfinite(pt.normal[1]) && std::isfinite(pt.normal[2]) &&
                                 std::isfinite(pt.overlap);
            if (!finite) { ++geometryNaNOrInfCount; continue; }
            if (pt.overlap < real_t(-1e-9)) ++geometryNegativeOverlapCount;
            const real_t normalLen = pt.normal.length();
            if (std::abs(normalLen - real_t(1)) > real_t(1e-6)) ++geometryBadNormalCount;
         }
      }

      // (2) Full dynamic energy-conservation check, ONE randomized corner/tilt/offset/speed
      // configuration per trial, useMultiPointBoxContact=true throughout.
      detail::CornerScenarioParams p;
      p.halfExtents             = halfExtents;
      p.density                 = real_t(0.6);
      p.restitution             = real_t(1); // strict ECC bound, see this function's own header comment
      p.friction                = real_t(0);
      p.poisson                 = real_t(0.22);
      p.collisionTime           = real_t(32.0078);
      p.dt                      = real_t(0.05);
      p.approachSpeed           = math::realRandom< real_t >(real_t(0.005), real_t(0.05));
      p.lateralOffsetY          = math::realRandom< real_t >(-halfExtents[1], halfExtents[1]);
      p.lateralOffsetZ          = math::realRandom< real_t >(-halfExtents[2], halfExtents[2]);
      p.tiltDegreesA            = math::realRandom< real_t >(real_t(-60), real_t(60));
      p.tiltDegreesB            = math::realRandom< real_t >(real_t(-60), real_t(60));
      p.spinRateB               = math::realRandom< real_t >(real_t(-0.01), real_t(0.01));
      p.separation              = real_t(12);
      p.tolerance               = tolerance;
      p.maxSteps                = maxStepsPerTrial;
      p.enableOutput            = false;
      p.useRoundedBoxCorners    = false;
      p.useMultiPointBoxContact = true;

      const auto r = detail::runCornerScenarioOnce(p);
      if (r.everContacted)
      {
         ++contactedCount;
         if (r.pass) ++energyPassCount;
         else
         {
            static uint_t failureLogCount = uint_t(0);
            if (failureLogCount < uint_t(20))
            {
               WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz FAILING TRIAL "
                  << trial << ": halfExtents=" << p.halfExtents << " approachSpeed=" << p.approachSpeed
                  << " lateralOffsetY=" << p.lateralOffsetY << " lateralOffsetZ=" << p.lateralOffsetZ
                  << " tiltDegreesA=" << p.tiltDegreesA << " tiltDegreesB=" << p.tiltDegreesB
                  << " spinRateB=" << p.spinRateB << " | KE_before=" << r.keBefore
                  << " KE_max=" << r.keMax << " (step " << r.keMaxStep << ") maxRelGain=" << r.keGainRel);
               ++failureLogCount;
            }
         }
         if (r.keGainRel > worstGainSeen)
         {
            worstGainSeen  = r.keGainRel;
            worstGainTrial = trial;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz  GEOMETRY: " << geometryChecked
      << " random pose pairs checked | NaN/Inf=" << geometryNaNOrInfCount
      << " | negative overlap=" << geometryNegativeOverlapCount
      << " | non-unit normal=" << geometryBadNormalCount
      << " | SAT-says-separated-but-vertex-embedded=" << geometrySatDisagreementCount);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz  DYNAMICS: " << numTrials
      << " randomized trials | contacted=" << contactedCount << " | energy-bound PASS="
      << energyPassCount << "/" << contactedCount << " | worst maxRelGain=" << worstGainSeen
      << " (trial " << worstGainTrial << ")");

   const bool geometryOk = geometryNaNOrInfCount == uint_t(0) && geometryNegativeOverlapCount == uint_t(0) &&
                            geometryBadNormalCount == uint_t(0) && geometrySatDisagreementCount == uint_t(0);
   const bool dynamicsOk = contactedCount > uint_t(0) && energyPassCount == contactedCount;
   const bool pass       = geometryOk && dynamicsOk;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz OVERALL: geometryOk=" << geometryOk
      << " dynamicsOk=" << dynamicsOk << "  " << (pass ? "PASS" : "FAIL"));
   return pass;
}


} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
