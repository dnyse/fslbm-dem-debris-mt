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
//! \file DEMValidationWallImpact.h
//! \brief Single-body and GJK-vs-analytical wall-impact DEM validation scenarios -- split out of
//! DEMValidationScenarios.h 2026-10-29. runWallImpactSphere/Cylinder/Box,
//! runGjkVsAnalyticalCylinderWallCheck, runBoxFaceDropRestitution, runBoxInclineSlideStraightness.

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
// 1-3: single-body wall-impact tests
//======================================================================================================================

inline bool runWallImpactSphere(const config::Config::BlockHandle& b)
{
   const real_t radius        = b.getParameter< real_t >("sphereRadius", real_t(1));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.85));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.15));
   const real_t speed         = b.getParameter< real_t >("impactSpeed", real_t(1));
   const real_t ratioStart    = b.getParameter< real_t >("tangentialRatioStart", real_t(0));
   const real_t ratioEnd      = b.getParameter< real_t >("tangentialRatioEnd", real_t(1.5));
   const real_t ratioStep     = b.getParameter< real_t >("tangentialRatioStep", real_t(0.25));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "wall_impact_sphere");

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);

   for (real_t ratio = ratioStart; ratio <= ratioEnd + real_t(1e-9); ratio += ratioStep)
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), real_t(50) * radius);
      const real_t gap = real_t(0.01) * radius;
      const size_t sphereIdx =
         engine.addSphere(Vec3(0, 0, radius + gap), Vec3(ratio * speed, 0, -speed), radius, density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_ratio" + std::to_string(ratio);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { sphereIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(sphereIdx);
      const real_t keBefore = engine.kineticEnergy(sphereIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, sphereIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(sphereIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(sphereIdx);
      const real_t keAfter   = engine.kineticEnergy(sphereIdx);

      const auto analytic = analytical::sphereWallImpact(radius, engine.ac().getMass(sphereIdx), preVel, Vec3(0),
                                                           wallNormal, e, mu);

      const std::string caseTag = "wallImpactSphere ratio=" + std::to_string(ratio);
      detail::logCaseResult(caseTag, simVel, analytic.postVelocityA, simAngVel, analytic.postAngularVelocityA,
                             tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const real_t sphereI = real_t(0.4) * engine.ac().getMass(sphereIdx) * radius * radius; // isotropic
      const real_t keExpected = detail::kineticEnergyOf(engine.ac().getMass(sphereIdx), analytic.postVelocityA,
                                                          analytic.postAngularVelocityA,
                                                          Mat3::makeDiagonalMatrix(sphereI));
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactSphere OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

inline bool runWallImpactCylinder(const config::Config::BlockHandle& b)
{
   const real_t radius        = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t length        = b.getParameter< real_t >("cylinderLength", real_t(4));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.85));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.15));
   const real_t speed         = b.getParameter< real_t >("impactSpeed", real_t(1));
   // Default 1.0 (equal tangential/normal, this app's own general oblique-impact sweep). Kodam et al.
   // (2010) Part II's actual Test 1 protocol is PURELY NORMAL (zero tangential component, zero initial
   // spin, frictionless) -- see dem_validation_wall_cylinder_kodam_test1.prm, which sets this to 0.
   const real_t tangentialFraction = b.getParameter< real_t >("tangentialSpeedFraction", real_t(1));
   const real_t tiltStartDeg  = b.getParameter< real_t >("tiltStartDeg", real_t(0));
   const real_t tiltEndDeg    = b.getParameter< real_t >("tiltEndDeg", real_t(90));
   const real_t tiltStepDeg   = b.getParameter< real_t >("tiltStepDeg", real_t(15));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "wall_impact_cylinder");
   // useRotationalEffectiveMass (see ModelParameters' own doc comment, SimulationParams.h, and
   // computeRotationalEffectiveMass, ParticleKernels.h): LinearSpringDashpot's stiffness/damping are
   // calibrated ONCE from the plain TRANSLATIONAL effective mass, which is only exact when the contact
   // point has zero lever arm about the line of centres -- exactly true at tilt=0/90 (on-axis contact)
   // and increasingly false at intermediate tilts, matching this test's own observed error pattern
   // (clean at 0/90, worst in the middle). Default false reproduces every existing prm's exact prior
   // behaviour; the Kodam-replica diagnostic run turns this on to test whether it explains the gap.
   const bool useRotationalEffectiveMass = b.getParameter< bool >("useRotationalEffectiveMass", false);

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);

   data::Cylinder refShape(radius, length);
   refShape.updateMassAndInertia(density);

   for (real_t tiltDeg = tiltStartDeg; tiltDeg <= tiltEndDeg + real_t(1e-9); tiltDeg += tiltStepDeg)
   {
      const detail::YRotation yrot(detail::degToRad(tiltDeg));
      // Cylinder's own support() gives the exact "reach" toward the wall for this tilt, matching the
      // formula production contact detection uses -- see AnalyticalSolutions.h's boxCornerWallImpact
      // comment for the same reasoning applied to the box case.
      const Vec3   bfDir  = yrot.mat.getTranspose() * (-wallNormal);
      const Vec3   lowest = yrot.mat * refShape.support(bfDir);
      const real_t reach  = -lowest[2];

      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), real_t(20) * (radius + length),
                        /*useAnalyticalContactDetection=*/false, /*useRoundedBoxCorners=*/false,
                        /*boxRoundingRadius=*/real_t(0.05), /*useMultiPointBoxContact=*/false,
                        /*useBoxOverlapCap=*/false, /*boxOverlapCapFraction=*/real_t(0.3),
                        /*useAxisHysteresis=*/false, /*axisHysteresisMargin=*/real_t(0.15),
                        /*useContactManifold=*/false, /*manifoldFadeRate=*/real_t(0.25),
                        useRotationalEffectiveMass);
      const real_t gap    = real_t(0.01) * radius;
      const size_t cylIdx = engine.addCylinder(Vec3(0, 0, reach + gap),
                                                Vec3(tangentialFraction * speed, 0, -speed), yrot.rot, radius,
                                                length, density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_tilt" + std::to_string(tiltDeg);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { cylIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(cylIdx);
      const real_t keBefore = engine.kineticEnergy(cylIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, cylIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(cylIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(cylIdx);
      const real_t keAfter   = engine.kineticEnergy(cylIdx);

      const closed_form_contacts::CylinderState cylState{ Vec3(0, 0, reach + gap), yrot.mat, radius, length };
      const Mat3 invInertiaWorld = math::transformMatrixRART(yrot.mat, refShape.getInvInertiaBF());
      const auto analytic = analytical::cylinderWallImpact(cylState, refShape.getMass(), invInertiaWorld, preVel,
                                                             Vec3(0), wallPos, wallNormal, e, mu);

      const std::string caseTag = "wallImpactCylinder tilt=" + std::to_string(tiltDeg) + "deg";
      detail::logCaseResult(caseTag, simVel, analytic.postVelocityA, simAngVel, analytic.postAngularVelocityA,
                             tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const Mat3   inertiaWorld = math::transformMatrixRART(yrot.mat, refShape.getInertiaBF());
      const real_t keExpected   = detail::kineticEnergyOf(refShape.getMass(), analytic.postVelocityA,
                                                            analytic.postAngularVelocityA, inertiaWorld);
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);

      // DIRECT comparison against Kodam et al. (2010) Part II's own stated numeric result (p. 5866): at
      // face-on (0deg) and band-on (90deg), the dimensionless rebound speed Vz+/Vz- equals EXACTLY the
      // restitution coefficient (their reported value: 0.85). This checks the simulated ratio against
      // that PAPER-STATED number directly -- independent of this file's own analytical solver -- rather
      // than only cross-checking against a formula we derived ourselves. Only meaningful at frictionless
      // 0/90deg (see dem_validation_wall_cylinder_kodam_test1.prm); silently skipped otherwise since
      // Kodam's paper gives no other explicit numeric target across the sweep.
      if (std::abs(mu) < real_t(1e-9) &&
          (std::abs(tiltDeg) < real_t(1e-6) || std::abs(tiltDeg - real_t(90)) < real_t(1e-6)))
      {
         const real_t simRebound = std::abs(simVel[2]) / std::abs(preVel[2]);
         const real_t relErr     = std::abs(simRebound - e) / e;
         const bool   pass       = relErr < real_t(0.02);
         acc.record(pass);
         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactCylinder tilt="
                                    << tiltDeg << "deg  [vs. Kodam et al. 2010 Part II, p. 5866, stated result] "
                                    << "simRebound=" << simRebound << " paperReported=" << e
                                    << " (=configured restitution, their calibrated value) relErr=" << relErr
                                    << "  " << (pass ? "PASS" : "FAIL"));
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactCylinder OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

//! Diagnostic-only cross-check (not itself a pass/fail validation experiment against a literature
//! result -- a self-consistency check between two of OUR OWN geometric methods): does mesa_pd's
//! generic iterative GJK/EPA solver agree with the exact closed-form Cylinder x HalfSpace contact
//! geometry that production actually uses (SolidCylinderContactDetection's cyl.support()-based
//! formula, ParticleKernels.h:338 -- algebraically identical to closed_form_contacts::rimPointTowardPlane,
//! PaperAnalyticalContacts.h:1980, the formula AnalyticalSolutions.h::cylinderWallImpact uses as
//! wallImpactCylinder's own analytical reference)?
//!
//! A true HalfSpace has NO support() function at all (mesa_pd::data::HalfSpace.h) -- GJK/EPA cannot
//! run on a Cylinder x HalfSpace pair in this codebase AT ALL, by construction, not by configuration
//! choice. This test approximates the wall as a very large, thin FloatingBox (which DOES have
//! support(), so GJK/EPA is well-defined) purely to get an independent third answer for the SAME
//! cylinder geometry, at the SAME tilt angles wallImpactCylinder's own Kodam-replica sweep uses.
inline bool runGjkVsAnalyticalCylinderWallCheck(const config::Config::BlockHandle& b)
{
   const real_t radius       = b.getParameter< real_t >("cylinderRadius", real_t(4));
   const real_t length       = b.getParameter< real_t >("cylinderLength", real_t(5.3));
   const real_t overlapDepth = b.getParameter< real_t >("overlapDepth", real_t(0.05));
   const real_t tiltStartDeg = b.getParameter< real_t >("tiltStartDeg", real_t(0));
   const real_t tiltEndDeg   = b.getParameter< real_t >("tiltEndDeg", real_t(90));
   const real_t tiltStepDeg  = b.getParameter< real_t >("tiltStepDeg", real_t(10));
   const real_t tolerance    = b.getParameter< real_t >("toleranceRelError", real_t(0.02));

   data::Cylinder cyl(radius, length);
   FloatingBox    wallProxy(Vec3(real_t(2000), real_t(2000), real_t(100)));
   wallProxy.updateMassAndInertia(real_t(1));
   const Vec3 wallPos(0, 0, -real_t(50)); // top face exactly at world z=0, matching HalfSpace's own plane
   const Rot3 wallRot;                    // identity
   const Vec3 wallNormal(0, 0, 1);

   detail::PassFailAccumulator acc;

   for (real_t tiltDeg = tiltStartDeg; tiltDeg <= tiltEndDeg + real_t(1e-9); tiltDeg += tiltStepDeg)
   {
      const detail::YRotation yrot(detail::degToRad(tiltDeg));

      // Exact support-point reach (same formula production contact detection + the analytical
      // reference both use), positioned to overlap the wall by exactly `overlapDepth`.
      const Vec3   bfDir  = yrot.mat.getTranspose() * (-wallNormal);
      const Vec3   lowest = yrot.mat * cyl.support(bfDir);
      const real_t reach  = -lowest[2];
      const Vec3   cylPos(0, 0, reach - overlapDepth);

      // (a) exact closed form -- identical to ParticleKernels.h's SolidCylinderContactDetection
      const Vec3   exactContact          = cylPos + yrot.mat * cyl.support(bfDir);
      const real_t exactPenetrationDepth = wallNormal * exactContact; // wall plane passes through the origin

      // (b) GJK/EPA against the large-box wall proxy
      mesa_pd::collision_detection::Support supportCyl(cylPos, yrot.rot, cyl);
      mesa_pd::collision_detection::Support supportWall(wallPos, wallRot, wallProxy);
      mesa_pd::collision_detection::GJK     gjk;
      bool   gjkFound            = false;
      Vec3   gjkPoint(0), gjkNormal(0);
      real_t gjkPenetrationDepth = 0;
      if (gjk.doGJKmargin(supportCyl, supportWall, real_t(0)))
      {
         mesa_pd::collision_detection::EPA epa;
         epa.useSphereOptimization(false);
         gjkFound =
            epa.doEPAmargin(supportCyl, supportWall, gjk, gjkNormal, gjkPoint, gjkPenetrationDepth, real_t(0));
      }

      const real_t overlapDiff =
         gjkFound ? std::abs(std::abs(gjkPenetrationDepth) - std::abs(exactPenetrationDepth)) : real_t(-1);
      const real_t pointDist = gjkFound ? (gjkPoint - exactContact).length() : real_t(-1);

      const bool pass = gjkFound && overlapDiff < tolerance * radius;
      acc.record(pass);
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] gjkVsAnalyticalCylinderWall tilt="
                                 << tiltDeg << "deg  exactPoint=" << exactContact
                                 << " exactOverlap=" << exactPenetrationDepth << "  gjkFound=" << gjkFound
                                 << " gjkPoint=" << gjkPoint << " gjkOverlap=" << gjkPenetrationDepth
                                 << "  overlapDiff=" << overlapDiff << " pointDist=" << pointDist << "  "
                                 << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] gjkVsAnalyticalCylinderWall OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

inline bool runWallImpactBox(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 2, 2));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.6));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.2));
   const real_t speed         = b.getParameter< real_t >("impactSpeed", real_t(1));
   // [Claude, 2026-10-29] Mirrors wallImpactCylinder's own tangentialSpeedFraction (above): default 1.0 keeps
   // every existing prm's exact prior behaviour (equal tangential/normal, this app's own general oblique-impact
   // sweep); set to 0.0 for a Kodam-et-al.-Test-1-style purely normal impact, now available for the box the same
   // way it already was for the cylinder (author: "can you recreate the kodam version for the cuboids").
   const real_t tangentialFraction = b.getParameter< real_t >("tangentialSpeedFraction", real_t(1));
   const real_t tiltStartDeg  = b.getParameter< real_t >("tiltStartDeg", real_t(5));
   const real_t tiltEndDeg    = b.getParameter< real_t >("tiltEndDeg", real_t(40));
   const real_t tiltStepDeg   = b.getParameter< real_t >("tiltStepDeg", real_t(10));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   // [Claude, 2026-09-30] was not wired for the box (always the plain translational effective mass); default false keeps
   // the previous behaviour.
   const bool   useRotationalEffectiveMass = b.getParameter< bool >("useRotationalEffectiveMass", false);
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "wall_impact_box");

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);
   const Vec3 edgeLength = real_t(2) * halfExtents;

   FloatingBox refShape(edgeLength);
   refShape.updateMassAndInertia(density);

   for (real_t tiltDeg = tiltStartDeg; tiltDeg <= tiltEndDeg + real_t(1e-9); tiltDeg += tiltStepDeg)
   {
      const detail::YRotation yrot(detail::degToRad(tiltDeg));
      const Vec3   bfDir  = yrot.mat.getTranspose() * (-wallNormal);
      const Vec3   corner = refShape.support(bfDir); // body-frame corner, used for POSITIONING (height) only
      const Vec3   lowest = yrot.mat * corner;
      const real_t reach  = -lowest[2];
      // A pure Y-axis tilt never changes a point's world Y (rotation about Y leaves Y invariant) or its
      // contribution to world Z (which depends only on body X/Z), so the "lowest" support corner is
      // always TIED between its +Y and -Y twin -- the true first contact is the full bottom EDGE
      // spanning Y, matching Housner's 2D rocking-block corner exactly in the X-Z cross-section, not a
      // single 3D corner. support()'s arbitrary tie-break otherwise injects a spurious Y lever arm (and,
      // via the K matrix's rotational coupling, corrupts the X-Z answer too) that the true symmetric
      // edge contact does not have -- confirmed empirically (a non-zero-Y support point predicted a
      // physically impossible negative post-bounce Z-velocity). Analytical contact point is the same
      // corner projected onto the Y=0 symmetry plane.
      const Vec3 analyticCornerBF(corner[0], real_t(0), corner[2]);

      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0),
                       real_t(30) * std::max({ halfExtents[0], halfExtents[1], halfExtents[2] }),
                       /*useAnalyticalContactDetection=*/false, /*useRoundedBoxCorners=*/false,
                       /*boxRoundingRadius=*/real_t(0.05), /*useMultiPointBoxContact=*/false,
                       /*useBoxOverlapCap=*/false, /*boxOverlapCapFraction=*/real_t(0.3),
                       /*useAxisHysteresis=*/false, /*axisHysteresisMargin=*/real_t(0.15),
                       /*useContactManifold=*/false, /*manifoldFadeRate=*/real_t(0.25), useRotationalEffectiveMass);
      const real_t gap = real_t(0.01) * std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
      // Horizontal velocity component aimed so the leading (lowest) corner strikes squarely, matching
      // Housner's planar rocking-block assumption (velocity, spin axis and contact point all in the
      // X-Z plane -- tilt is about world Y, so no Y motion is introduced).
      const size_t boxIdx = engine.addBox(Vec3(0, 0, reach + gap),
                                            Vec3(tangentialFraction * speed, 0, -speed), yrot.rot, edgeLength,
                                            density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_tilt" + std::to_string(tiltDeg);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { boxIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(boxIdx);
      const real_t keBefore = engine.kineticEnergy(boxIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, boxIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(boxIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(boxIdx);
      const real_t keAfter   = engine.kineticEnergy(boxIdx);

      const auto analytic = analytical::boxCornerWallImpact(analyticCornerBF, yrot.mat, refShape.getMass(),
                                                              refShape.getInvInertiaBF(), preVel, Vec3(0), wallNormal,
                                                              e, mu);

      const std::string caseTag = "wallImpactBox tilt=" + std::to_string(tiltDeg) + "deg";
      detail::logCaseResult(caseTag, simVel, analytic.postVelocityA, simAngVel, analytic.postAngularVelocityA,
                             tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const Mat3   inertiaWorld = math::transformMatrixRART(yrot.mat, refShape.getInertiaBF());
      const real_t keExpected   = detail::kineticEnergyOf(refShape.getMass(), analytic.postVelocityA,
                                                            analytic.postAngularVelocityA, inertiaWorld);
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactBox OVERALL: " << (acc.allPassed() ? "PASS" : "FAIL") << " ("
                                                                        << acc.passed << "/" << acc.total << ")");
   return acc.allPassed();
}

//======================================================================================================================
// 3b: literature-grounded cuboid single-particle regression checks (Fraige et al. 2008 / Mack et al. 2011)
//======================================================================================================================

//! Literal replica of Mack et al. (2011)'s single-cube coefficient-of-restitution drop test methodology
//! (Powder Technology 214:431-442, Section 2.4: e = sqrt(h1/h0), measured for real hexahedra as e=0.54,
//! the default below) -- but FACE-first (flush, zero tilt), not corner-first like wallImpactBox's own
//! sweep (which starts at tiltStartDeg=5 and never exercises the flush case). This is deliberately the
//! EXACT degenerate configuration that Cylinder/Box::support()'s tie-break bug (fixed 2026-09-13, see
//! DEM_VALIDATION_MATH_REFERENCE.md's box-support() note) corrupted: pre-fix, a flush box's support()
//! picked an arbitrary corner instead of the true face centre, injecting a spurious lever arm/torque even
//! though the true physics here is a trivial 1-DOF normal bounce by symmetry -- analyticVel =
//! (0,0,+e*impactSpeed), analyticAngVel = (0,0,0) exactly, no AnalyticalSolutions.h call needed. Any
//! nonzero simulated angular velocity, or a normal-velocity ratio deviating from the input e, is
//! therefore a direct regression signal for that bug, not a friction/calibration artifact -- there is no
//! tangential velocity here to mobilise frictionCoefficient at all.
inline bool runBoxFaceDropRestitution(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 1, 1));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   // Mack et al. (2011), Powder Technology 214:431-442, Section 2.4 / Table 2: measured hexahedron
   // (cube) particle-wall coefficient of restitution e=0.54 and sliding friction mu=0.28 (Section 2.3).
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.54));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.28));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const real_t speedStart    = b.getParameter< real_t >("impactSpeedStart", real_t(0.5));
   const real_t speedEnd      = b.getParameter< real_t >("impactSpeedEnd", real_t(2.0));
   const real_t speedStep     = b.getParameter< real_t >("impactSpeedStep", real_t(0.5));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "box_face_drop_restitution");

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);
   const Vec3 edgeLength = real_t(2) * halfExtents;

   for (real_t speed = speedStart; speed <= speedEnd + real_t(1e-9); speed += speedStep)
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0),
                        real_t(30) * std::max({ halfExtents[0], halfExtents[1], halfExtents[2] }));
      const real_t gap = real_t(0.01) * std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
      const size_t boxIdx =
         engine.addBox(Vec3(0, 0, halfExtents[2] + gap), Vec3(0, 0, -speed), Rot3(), edgeLength, density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_v" + std::to_string(speed);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { boxIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(boxIdx);
      const real_t keBefore = engine.kineticEnergy(boxIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, boxIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(boxIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(boxIdx);
      const real_t keAfter   = engine.kineticEnergy(boxIdx);

      // Trivial by symmetry (flush, zero tangential velocity): normal velocity flips and scales by e,
      // angular velocity stays exactly zero -- no AnalyticalSolutions.h call needed for this case.
      // v_z^+ = -e*v_z^-, e=0.54 (Mack et al. 2011, Table 2)
      const Vec3 analyticVel(0, 0, -e * preVel[2]);
      const Vec3 analyticAngVel(0, 0, 0);

      const std::string caseTag = "boxFaceDropRestitution v=" + std::to_string(speed);
      detail::logCaseResult(caseTag, simVel, analyticVel, simAngVel, analyticAngVel, tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const real_t keExpected = real_t(0.5) * engine.ac().getMass(boxIdx) * analyticVel.sqrLength();
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] boxFaceDropRestitution OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")  [Mack et al. 2011, Powder Technology 214:431-442, measured hexahedron "
                              << "e=0.54, mu=0.28]");
   return acc.allPassed();
}

//! Fraige et al. (2008)'s "slide test" methodology (Powder Technology 186:224-240, Section 3.6 / Figs.
//! 7-8): a cube resting FACE-flush (not corner-first) on a plane inclined steeper than the friction angle
//! should slide straight down-slope with no lateral deviation -- their own single-contact-point DEM-Cube
//! version showed the contact force "jumping" between corners as the box translates, injecting a
//! spurious lateral force/torque, fixed there by a multiple-contact-point model. This is the SAME family
//! of degenerate flush-contact tie-break bug fixed in Cylinder/Box::support() this session
//! (2026-09-13) -- a box resting flush on a tilted plane is exactly the configuration where the old
//! support() implementation picked an arbitrary corner instead of the true face centre. The plane is
//! tilted about world X, so by symmetry the box's world-X position AND the two OUT-OF-PLANE angular
//! velocity components (wY, wZ, which would only appear from an asymmetric/"jumping" contact point) must
//! stay EXACTLY zero for all time -- no closed-form dynamics needed, just mirror symmetry -- so any
//! nonzero drift there is a direct regression signal for that bug. The IN-PLANE pitch component (wX,
//! rotation about the incline's own tilt axis) is NOT expected to be zero throughout: the box is released
//! through a small gap and free-falls under world gravity (not purely along the tilted wall normal)
//! before first contact, so a real, physically legitimate settling/rocking transient in wX is expected
//! immediately after impact (confirmed empirically: it decays from a transient ~0.07 rad/s down to
//! ~0.007 rad/s by the last 20% of a 2000-step run) -- only checked in the last `steadyStateFraction` of
//! the run, where it should have damped down close to zero for genuine steady sliding.
inline bool runBoxInclineSlideStraightness(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents  = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 1, 1));
   const real_t density      = b.getParameter< real_t >("density", real_t(1));
   const real_t e            = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   // Fraige et al. (2008), Powder Technology 186:224-240, Table 2: measured cube-wall sliding friction
   // coefficient mu=0.4.
   const real_t mu           = b.getParameter< real_t >("frictionCoefficient", real_t(0.4));
   const real_t poisson      = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt           = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t gravity      = b.getParameter< real_t >("gravityMagnitude", real_t(9.81));
   // Must exceed atan(mu) (~21.8deg for mu=0.4) for the box to actually slide rather than stick.
   const real_t inclineDeg   = b.getParameter< real_t >("inclineAngleDeg", real_t(30));
   const uint_t slideSteps   = b.getParameter< uint_t >("slideSteps", uint_t(2000));
   const real_t steadyStateFraction     = b.getParameter< real_t >("steadyStateFraction", real_t(0.2));
   const real_t toleranceLateralAbs     = b.getParameter< real_t >("toleranceLateralAbs", real_t(1e-6));
   const real_t toleranceOutOfPlaneAbs  = b.getParameter< real_t >("toleranceOutOfPlaneAngVelAbs", real_t(1e-6));
   const real_t toleranceSteadyPitchAbs = b.getParameter< real_t >("toleranceSteadyPitchAbs", real_t(0.01));
   const uint_t vtkSpacing   = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase = b.getParameter< std::string >("outputBaseName", "box_incline_slide");

   const real_t angle      = detail::degToRad(inclineDeg);
   const Rot3   rot(Vec3(1, 0, 0), angle);
   const Vec3   wallNormal = rot.getMatrix() * Vec3(0, 0, 1);
   const Vec3   wallPos(0, 0, 0);
   const Vec3   edgeLength = real_t(2) * halfExtents;

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity,
                     real_t(30) * std::max({ halfExtents[0], halfExtents[1], halfExtents[2] }));
   // Box's local frame shares the wall's own rotation, so in the wall's local frame this is exactly an
   // axis-aligned box resting flush on a z=0 plane -- distance from centre to wall along the normal is
   // exactly halfExtents[2], regardless of inclineDeg. A tiny gap lets gravity pull it into full contact
   // over the first few steps rather than starting at an ambiguous exactly-zero overlap.
   const real_t gap      = real_t(0.01) * std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
   const Vec3   startPos = wallNormal * (halfExtents[2] + gap);
   const size_t boxIdx   = engine.addBox(startPos, Vec3(0, 0, 0), rot, edgeLength, density);
   const size_t wallIdx  = engine.addWall(wallPos, wallNormal);
   engine.finalizeSetup();
   static_cast< void >(wallIdx);

   const std::string tag = outBase + "_" + std::to_string(inclineDeg) + "deg";
   engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                   { boxIdx });
   engine.setupVtk(tag, vtkSpacing);

   const uint_t steadyStateStartStep = uint_t(real_t(slideSteps) * (real_t(1) - steadyStateFraction));
   real_t maxLateralX    = real_t(0);
   real_t maxOutOfPlaneW = real_t(0);
   real_t maxSteadyPitch = real_t(0);
   for (uint_t s = 0; s < slideSteps; ++s)
   {
      engine.step();
      engine.writeCsvRow();
      engine.writeVtkIfDue();
      const Vec3 w = engine.ac().getAngularVelocity(boxIdx);
      maxLateralX    = std::max(maxLateralX, std::abs(engine.ac().getPosition(boxIdx)[0]));
      maxOutOfPlaneW = std::max({ maxOutOfPlaneW, std::abs(w[1]), std::abs(w[2]) });
      if (s >= steadyStateStartStep) maxSteadyPitch = std::max(maxSteadyPitch, std::abs(w[0]));
   }

   detail::PassFailAccumulator acc;
   const bool lateralPass = maxLateralX < toleranceLateralAbs;
   const bool outOfPlanePass = maxOutOfPlaneW < toleranceOutOfPlaneAbs;
   const bool steadyPitchPass = maxSteadyPitch < toleranceSteadyPitchAbs;
   acc.record(lateralPass);
   acc.record(outOfPlanePass);
   acc.record(steadyPitchPass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] boxInclineSlideStraightness incline=" << inclineDeg << "deg  "
                              << "maxLateralX=" << maxLateralX << " (tol=" << toleranceLateralAbs << ")  "
                              << (lateralPass ? "PASS" : "FAIL") << "   maxOutOfPlaneAngVel(wY,wZ)="
                              << maxOutOfPlaneW << " (tol=" << toleranceOutOfPlaneAbs << ")  "
                              << (outOfPlanePass ? "PASS" : "FAIL") << "   maxSteadyStatePitch(wX, last "
                              << (steadyStateFraction * real_t(100)) << "%)=" << maxSteadyPitch
                              << " (tol=" << toleranceSteadyPitchAbs << ")  " << (steadyPitchPass ? "PASS" : "FAIL"));
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] boxInclineSlideStraightness OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")  [Fraige et al. 2008, Powder Technology 186:224-240, mu=0.4]");
   return acc.allPassed();
}


} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
