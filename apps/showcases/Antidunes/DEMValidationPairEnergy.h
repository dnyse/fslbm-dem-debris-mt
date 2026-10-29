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
//! \file DEMValidationPairEnergy.h
//! \brief Pairwise central head-on restitution/energy-conservation DEM validation scenarios --
//! split out of DEMValidationScenarios.h 2026-10-29. runPairRestitution,
//! runCylinderPairEnergyFuzz, runCuboidPairEnergyFuzz, runEnergyConservation.

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
// 4-5: pairwise central head-on tests (all 6 shape-pair combinations)
//======================================================================================================================

namespace detail
{
//! Creates a shape by name ("sphere"/"box"/"cylinder") centred at `pos`, moving at `vel`, oriented so
//! that a pair of these shapes approaching head-on along world X makes CENTRAL, face-aligned contact
//! (r x n = 0 on both sides -- no lever arm): box/cylinder are sized and oriented so the contact point
//! sits exactly on the line of centres (box: flat face perpendicular to X, unrotated; cylinder: axis
//! along world Y, so its band's nearest surface point to a partner on the X-axis is exactly (radius,0,0)).
inline size_t createPairMember(DEMEngine& engine, const std::string& shape, const Vec3& pos, const Vec3& vel,
                                real_t sphereRadius, const Vec3& boxHalfExtents, real_t cylRadius, real_t cylLength,
                                real_t density)
{
   if (shape == "sphere") return engine.addSphere(pos, vel, sphereRadius, density);
   if (shape == "box") return engine.addBox(pos, vel, Rot3(), real_t(2) * boxHalfExtents, density);
   // cylinder: body-z aligned with world Y (band faces the partner along X)
   Rot3 rot(Vec3(1, 0, 0), math::pi / real_t(2));
   return engine.addCylinder(pos, vel, rot, cylRadius, cylLength, density);
}

inline const std::vector< std::pair< std::string, std::string > >& allShapePairs()
{
   static const std::vector< std::pair< std::string, std::string > > pairs = {
      { "sphere", "sphere" }, { "sphere", "box" },   { "sphere", "cylinder" },
      { "box", "box" },       { "box", "cylinder" }, { "cylinder", "cylinder" }
   };
   return pairs;
}

//! Runs `engine` until the tracked dynamic-body pair (idxA, idxB) has made and then broken contact, or
//! maxSteps is reached. Returns true if a contact was actually observed (false = the bodies never
//! touched -- misconfigured scenario, not a physical result).
inline bool runUntilSeparated(DEMEngine& engine, size_t idxA, size_t idxB, uint_t maxSteps, bool writeOutput)
{
   bool   everContacted = false;
   real_t penetration   = real_t(0);
   for (uint_t s = 0; s < maxSteps; ++s)
   {
      const bool contactNow = engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0);
      if (contactNow) everContacted = true;
      engine.step();
      if (writeOutput)
      {
         engine.writeCsvRow();
         engine.writeVtkIfDue();
      }
      if (everContacted && !contactNow && s > 0) break; // just separated this step
   }
   return everContacted;
}
} // namespace detail

inline bool runPairRestitution(const config::Config::BlockHandle& b)
{
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 2, 2));
   const real_t cylRadius      = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t cylLength      = b.getParameter< real_t >("cylinderLength", real_t(6));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t e              = b.getParameter< real_t >("restitutionCoefficient", real_t(0.7));
   const real_t mu             = b.getParameter< real_t >("frictionCoefficient", real_t(0.0));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance      = b.getParameter< real_t >("toleranceRelError", real_t(0.1));
   const real_t approachSpeed  = b.getParameter< real_t >("approachSpeed", real_t(1));
   const real_t separation     = b.getParameter< real_t >("initialSeparation", real_t(20));
   const uint_t vtkSpacing     = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase   = b.getParameter< std::string >("outputBaseName", "pair_restitution");

   detail::PassFailAccumulator acc;
   const real_t analyticRatio = analytical::twoBodyHeadOnRestitution(e);

   for (const auto& pr : detail::allShapePairs())
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), separation + real_t(20));
      const size_t idxA = detail::createPairMember(engine, pr.first, Vec3(-separation / 2, 0, 0),
                                                     Vec3(approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      const size_t idxB = detail::createPairMember(engine, pr.second, Vec3(separation / 2, 0, 0),
                                                     Vec3(-approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      engine.finalizeSetup();

      const std::string tag = outBase + "_" + pr.first + "_" + pr.second;
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { idxA, idxB });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   pBefore  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                              engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);
      const real_t keBefore = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);

      const uint_t maxSteps = uint_t(2) * uint_t(separation / (approachSpeed * dt)) +
                               uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runUntilSeparated(engine, idxA, idxB, maxSteps, true);

      const Vec3   pAfter  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                             engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);
      const real_t keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);

      const real_t vRelBefore = approachSpeed; // by construction: closing speed along X at t=0
      const real_t vRelAfter  = (engine.ac().getLinearVelocity(idxA) - engine.ac().getLinearVelocity(idxB))[0];
      const real_t simRatio   = vRelAfter / vRelBefore;
      const real_t relErr     = std::abs(simRatio - analyticRatio) / std::max(real_t(1e-8), std::abs(analyticRatio));
      const bool   pass       = relErr < tolerance;
      acc.record(pass);

      const std::string pairTag = std::string("pairRestitution ") + pr.first + "-" + pr.second;
      detail::logMomentumConservation(pairTag, pBefore, pAfter, acc);
      detail::logConservation(pairTag, keBefore, keAfter, acc);
      // Exact closed form for a CENTRAL (r x n = 0, guaranteed by createPairMember's geometry),
      // FRICTIONLESS (mu == 0, this test's default) head-on collision: total KE splits into COM-frame
      // translation (unchanged -- momentum conservation fixes the COM velocity) plus relative-motion KE,
      // which scales by e^2. Only claimed exact when mu is actually 0 -- with friction on, rotational
      // coupling makes the energy budget depend on the same friction-model details already discussed for
      // the wall-impact tests, so the check is skipped rather than silently claiming false precision.
      if (std::abs(mu) < real_t(1e-9))
      {
         const real_t mA = engine.ac().getMass(idxA), mB = engine.ac().getMass(idxB);
         const real_t reducedMass = mA * mB / (mA + mB);
         const real_t keExpected  = keBefore - real_t(0.5) * reducedMass * vRelBefore * vRelBefore *
                                                  (real_t(1) - e * e);
         detail::logEnergyBudget(pairTag, keAfter, keExpected, acc, tolerance);
      }

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] pairRestitution "
                                 << pr.first << "-" << pr.second << "  simRatio=" << simRatio
                                 << " analyticRatio=" << analyticRatio << " relErr=" << relErr << "  "
                                 << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] pairRestitution OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

//! [Claude, 2026-09-30] Energy conservation of OBLIQUE cylinder-cylinder collisions under each cylinder detection mode
//! (thesis Section 5.2, "continuity"). Two free cylinders with uniformly random orientations collide with a random
//! lateral offset (so that rims, bands and faces all meet), e = 1, mu = 0, no gravity: the kinetic energy after the
//! collision must equal the energy before. A detection that jumps in depth or normal while the bodies move injects or
//! removes energy. The SAME set of initial conditions is run with each mode:
//!   gjk      GJK/EPA only (useAnalyticalContactDetection false)
//!   hybrid   closed form + GJK/EPA fallback (default analytical mode)
//!   pure     closed form without fallback (usePureAnalyticalCylinder)
//!   refined  closed form + support-function refinement (refineCylinderContacts)
//! Reported per mode: mean and max |KE_after/KE_before - 1| over the collisions that made contact, max relative momentum
//! error, max penetration / radius, and collisions that never made contact although the other modes did.
inline bool runCylinderPairEnergyFuzz(const config::Config::BlockHandle& b)
{
   const real_t radius        = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t length        = b.getParameter< real_t >("cylinderLength", real_t(3));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t speed         = b.getParameter< real_t >("approachSpeed", real_t(1));
   const uint_t numCollisions = b.getParameter< uint_t >("numCollisions", uint_t(200));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t seed          = b.getParameter< uint_t >("seed", uint_t(4242));

   struct Init { Rot3 rotA, rotB; Vec3 offset; };
   std::mt19937                             rng(seed);
   std::uniform_real_distribution< real_t > u(real_t(0), real_t(1));
   auto uniformRot = [&]() {
      const real_t u1 = u(rng), u2 = u(rng), u3 = u(rng);
      const real_t a = std::sqrt(real_t(1) - u1), bb = std::sqrt(u1);
      const real_t w = a * std::sin(real_t(2) * math::pi * u2), x = a * std::cos(real_t(2) * math::pi * u2);
      const real_t y = bb * std::sin(real_t(2) * math::pi * u3), z = bb * std::cos(real_t(2) * math::pi * u3);
      Rot3 r;
      const Vec3 axis = (w < 0 ? real_t(-1) : real_t(1)) * Vec3(x, y, z);
      if (axis.length() > real_t(1e-12)) r.rotate(axis.getNormalized(), real_t(2) * std::acos(std::min(std::abs(w), real_t(1))));
      return r;
   };
   const real_t reach = radius + real_t(0.5) * length;
   std::vector< Init > inits;
   for (uint_t c = 0; c < numCollisions; ++c)
      inits.push_back(Init{ uniformRot(), uniformRot(),
                            Vec3(0, (u(rng) - real_t(0.5)) * real_t(2) * reach, (u(rng) - real_t(0.5)) * real_t(2) * reach) });

   std::vector< std::string > modes;
   {
      std::stringstream ss(b.getParameter< std::string >("modes", "gjk,hybrid,pure,refined"));
      for (std::string item; std::getline(ss, item, ',');) modes.push_back(item);
   }
   std::map< std::string, std::vector< real_t > > keErr, pErr, pen;
   std::map< std::string, std::vector< bool > >   contacted;
   const real_t sep0 = real_t(2.5) * reach;
   for (const auto& mode : modes)
   {
      for (const auto& in : inits)
      {
         DEMEngine engine(dt, real_t(1), real_t(0), poisson, collisionTime, real_t(0), real_t(10) * reach,
                          /*useAnalyticalContactDetection=*/mode != "gjk");
         engine.computeContacts().setUsePureAnalyticalCylinder(mode == "pure");
         engine.computeContacts().setRefineCylinderContacts(mode == "refined");
         const size_t A = engine.addCylinder(Vec3(-sep0 / 2, 0, 0) - in.offset * real_t(0.5), Vec3(speed / 2, 0, 0), in.rotA,
                                             radius, length, density);
         const size_t B = engine.addCylinder(Vec3(sep0 / 2, 0, 0) + in.offset * real_t(0.5), Vec3(-speed / 2, 0, 0), in.rotB,
                                             radius, length, density);
         engine.finalizeSetup();
         const real_t ke0 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p0  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         bool   touched = false;
         real_t maxPen = 0, penetration = 0;
         uint_t quiet = 0;
         const uint_t maxSteps = uint_t(real_t(3) * sep0 / (speed * dt));
         for (uint_t s = 0; s < maxSteps; ++s)
         {
            engine.step();
            const bool c = engine.queryContact(A, B, penetration) && penetration < real_t(0);
            if (mode == "refined" && std::getenv("ANTIDUNES_ENERGY_DEBUG") &&
                int(&in - &inits[0]) == std::atoi(std::getenv("ANTIDUNES_ENERGY_DEBUG")))
            {
               const paper_analytical::CylinderState ci{ engine.ac().getPosition(A), engine.ac().getRotation(A).getMatrix(), radius, length };
               const paper_analytical::CylinderState cj{ engine.ac().getPosition(B), engine.ac().getRotation(B).getMatrix(), radius, length };
               const auto raw = paper_analytical::cylinderCylinderContact(ci, cj);
               const auto rf  = paper_analytical::refineCylinderContact(ci, cj, raw);
               if (c || rf.contact || raw.contact)
               WALBERLA_LOG_INFO_ON_ROOT("[energy-debug] step " << s << " gjkContact=" << c << " gjkPen=" << -penetration << " raw=" << raw.scenario
                  << " " << raw.overlap << " refined=" << rf.scenario << " " << rf.overlap << " n=" << rf.normal << " p=" << rf.point
                  << " rawP=" << raw.point << " KE=" << engine.kineticEnergy(A) + engine.kineticEnergy(B));
            }
            if (c) { touched = true; quiet = 0; maxPen = std::max(maxPen, -penetration); }
            else if (touched && ++quiet > uint_t(20) * uint_t(collisionTime / dt + real_t(1))) break;
         }
         const real_t ke1 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p1  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         contacted[mode].push_back(touched);
         keErr[mode].push_back(touched ? std::abs(ke1 / ke0 - real_t(1)) : real_t(0));
         pErr[mode].push_back((p1 - p0).length() / (engine.ac().getMass(A) * speed));
         pen[mode].push_back(maxPen / radius);
      }
   }

   detail::PassFailAccumulator acc;
   for (const auto& mode : modes)
   {
      real_t sum = 0, mx = 0, pmx = 0, penmx = 0;
      int    n = 0, over = 0, missedVsGjk = 0;
      for (size_t c = 0; c < inits.size(); ++c)
      {
         if (contacted.count("gjk") && contacted["gjk"][c] && !contacted[mode][c]) ++missedVsGjk;
         if (!contacted[mode][c]) continue;
         ++n; sum += keErr[mode][c]; mx = std::max(mx, keErr[mode][c]);
         pmx = std::max(pmx, pErr[mode][c]); penmx = std::max(penmx, pen[mode][c]);
         if (keErr[mode][c] > tolerance)
         {
            ++over;
            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cylinderPairEnergyFuzz mode=" << mode << " collision " << c
                                       << " |KE ratio-1|=" << keErr[mode][c] << " maxPenetration/R=" << pen[mode][c]);
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cylinderPairEnergyFuzz mode=" << mode << " collisions=" << n << "/"
         << inits.size() << " mean|KE ratio-1|=" << (n ? sum / real_t(n) : real_t(0)) << " max=" << mx
         << " over tolerance(" << tolerance << ")=" << over << " maxMomentumErr=" << pmx << " maxPenetration/R=" << penmx
         << " noContactWhereGjkHad=" << missedVsGjk);
      if (mode == "gjk" || mode == "refined") acc.record(over == 0 && missedVsGjk == 0);
   }
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cylinderPairEnergyFuzz OVERALL: " << (acc.allPassed() ? "PASS" : "FAIL") << " ("
                             << acc.passed << "/" << acc.total << ", gjk and refined modes checked)");
   return acc.allPassed();
}

//! [Claude, 2026-10-01] Cuboid version of runCylinderPairEnergyFuzz: oblique elastic (e = 1, mu = 0) collisions of two
//! freely rotating cuboids with random orientations and lateral offsets, no gravity; kinetic energy must be conserved.
//! Modes: "gjk" (GJK/EPA) and "sat" (closed-form separating-axis test, as in the drop tests). Written to test whether
//! the energy gain seen in the randomised cuboid drop tests comes from the SAT contact point (placed on the line between
//! the centres, not at the contact).
inline bool runCuboidPairEnergyFuzz(const config::Config::BlockHandle& b)
{
   const Vec3   halfExt       = b.getParameter< Vec3 >("cuboidHalfExtents", Vec3(6, 1, real_t(1.5)));
   const real_t radius        = std::min(halfExt[0], std::min(halfExt[1], halfExt[2])); // penetration reported / this
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t speed         = b.getParameter< real_t >("approachSpeed", real_t(1));
   const uint_t numCollisions = b.getParameter< uint_t >("numCollisions", uint_t(200));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t seed          = b.getParameter< uint_t >("seed", uint_t(4242));

   struct Init { Rot3 rotA, rotB; Vec3 offset; };
   std::mt19937                             rng(seed);
   std::uniform_real_distribution< real_t > u(real_t(0), real_t(1));
   auto uniformRot = [&]() {
      const real_t u1 = u(rng), u2 = u(rng), u3 = u(rng);
      const real_t a = std::sqrt(real_t(1) - u1), bb = std::sqrt(u1);
      const real_t w = a * std::sin(real_t(2) * math::pi * u2), x = a * std::cos(real_t(2) * math::pi * u2);
      const real_t y = bb * std::sin(real_t(2) * math::pi * u3), z = bb * std::cos(real_t(2) * math::pi * u3);
      Rot3 r;
      const Vec3 axis = (w < 0 ? real_t(-1) : real_t(1)) * Vec3(x, y, z);
      if (axis.length() > real_t(1e-12)) r.rotate(axis.getNormalized(), real_t(2) * std::acos(std::min(std::abs(w), real_t(1))));
      return r;
   };
   const real_t reach = halfExt.length();
   std::vector< Init > inits;
   for (uint_t c = 0; c < numCollisions; ++c)
      inits.push_back(Init{ uniformRot(), uniformRot(),
                            Vec3(0, (u(rng) - real_t(0.5)) * real_t(2) * reach, (u(rng) - real_t(0.5)) * real_t(2) * reach) });

   std::vector< std::string > modes;
   {
      std::stringstream ss(b.getParameter< std::string >("modes", "gjk,sat"));
      for (std::string item; std::getline(ss, item, ',');) modes.push_back(item);
   }
   std::map< std::string, std::vector< real_t > > keErr, pErr, pen;
   std::map< std::string, std::vector< bool > >   contacted;
   const real_t sep0 = real_t(2.5) * reach;
   for (const auto& mode : modes)
   {
      for (const auto& in : inits)
      {
         // modes: gjk (GJK/EPA), sat (plain SAT, as in the drop tests), multipoint (useMultiPointBoxContact).
         // [Claude, 2026-10-29] "fraige" mode (useFraigeBoxContactModel) removed -- see
         // archive/unused_contact_methods/ (fraige_multipoint_box_contact.h + fraige_validation_tests.cpp.excerpt).
         DEMEngine engine(dt, real_t(1), real_t(0), poisson, collisionTime, real_t(0), real_t(10) * reach,
                          /*useAnalyticalContactDetection=*/mode != "gjk", /*useRoundedBoxCorners=*/false,
                          /*boxRoundingRadius=*/real_t(0.05), /*useMultiPointBoxContact=*/mode == "multipoint",
                          false, real_t(0.3), false, real_t(0.15), false, real_t(0.25), false);
         const size_t A = engine.addBox(Vec3(-sep0 / 2, 0, 0) - in.offset * real_t(0.5), Vec3(speed / 2, 0, 0), in.rotA,
                                        real_t(2) * halfExt, density);
         const size_t B = engine.addBox(Vec3(sep0 / 2, 0, 0) + in.offset * real_t(0.5), Vec3(-speed / 2, 0, 0), in.rotB,
                                        real_t(2) * halfExt, density);
         engine.finalizeSetup();
         const real_t ke0 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p0  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         bool   touched = false;
         real_t maxPen = 0, penetration = 0;
         uint_t quiet = 0;
         const uint_t maxSteps = uint_t(real_t(3) * sep0 / (speed * dt));
         for (uint_t s = 0; s < maxSteps; ++s)
         {
            engine.step();
            const bool c = engine.queryContact(A, B, penetration) && penetration < real_t(0);
            if (std::getenv("ANTIDUNES_ENERGY_DEBUG") && int(&in - &inits[0]) == std::atoi(std::getenv("ANTIDUNES_ENERGY_DEBUG")))
            {
               const auto r = paper_analytical::boxBoxSAT(engine.ac().getPosition(A), engine.ac().getRotation(A).getMatrix(),
                                                          halfExt, engine.ac().getPosition(B),
                                                          engine.ac().getRotation(B).getMatrix(), halfExt);
               if (c || r.contact)
                  WALBERLA_LOG_INFO_ON_ROOT("[energy-debug] mode=" << mode << " step " << s << " contact=" << c << " pen="
                     << -penetration << " sat=" << r.overlap << " n=" << r.normal << " p=" << r.point
                     << " KE=" << engine.kineticEnergy(A) + engine.kineticEnergy(B));
            }
            if (c) { touched = true; quiet = 0; maxPen = std::max(maxPen, -penetration); }
            else if (touched && ++quiet > uint_t(20) * uint_t(collisionTime / dt + real_t(1))) break;
         }
         const real_t ke1 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p1  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         contacted[mode].push_back(touched);
         keErr[mode].push_back(touched ? std::abs(ke1 / ke0 - real_t(1)) : real_t(0));
         pErr[mode].push_back((p1 - p0).length() / (engine.ac().getMass(A) * speed));
         pen[mode].push_back(maxPen / radius);
      }
   }

   detail::PassFailAccumulator acc;
   for (const auto& mode : modes)
   {
      real_t sum = 0, mx = 0, pmx = 0, penmx = 0;
      int    n = 0, over = 0, missedVsGjk = 0;
      for (size_t c = 0; c < inits.size(); ++c)
      {
         if (contacted.count("gjk") && contacted["gjk"][c] && !contacted[mode][c]) ++missedVsGjk;
         if (!contacted[mode][c]) continue;
         ++n; sum += keErr[mode][c]; mx = std::max(mx, keErr[mode][c]);
         pmx = std::max(pmx, pErr[mode][c]); penmx = std::max(penmx, pen[mode][c]);
         if (keErr[mode][c] > tolerance)
         {
            ++over;
            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cuboidPairEnergyFuzz mode=" << mode << " collision " << c
                                       << " |KE ratio-1|=" << keErr[mode][c] << " maxPenetration/R=" << pen[mode][c]);
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cuboidPairEnergyFuzz mode=" << mode << " collisions=" << n << "/"
         << inits.size() << " mean|KE ratio-1|=" << (n ? sum / real_t(n) : real_t(0)) << " max=" << mx
         << " over tolerance(" << tolerance << ")=" << over << " maxMomentumErr=" << pmx << " maxPenetration/R=" << penmx
         << " noContactWhereGjkHad=" << missedVsGjk);
      acc.record(over == 0 && missedVsGjk == 0);
   }
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cuboidPairEnergyFuzz OVERALL: " << (acc.allPassed() ? "PASS" : "FAIL") << " ("
                             << acc.passed << "/" << acc.total << ", all modes checked)");
   return acc.allPassed();
}

inline bool runEnergyConservation(const config::Config::BlockHandle& b)
{
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 2, 2));
   const real_t cylRadius      = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t cylLength      = b.getParameter< real_t >("cylinderLength", real_t(6));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance      = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const real_t approachSpeed  = b.getParameter< real_t >("approachSpeed", real_t(1));
   const real_t separation     = b.getParameter< real_t >("initialSeparation", real_t(20));
   const uint_t vtkSpacing     = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase   = b.getParameter< std::string >("outputBaseName", "energy_conservation");

   detail::PassFailAccumulator acc;
   // Feng (2023)'s ECC criterion: e=1, mu=0, no gravity -- total KE must be conserved across an
   // elastic collision. This is the theoretically rigorous minimum bar for a contact model, and
   // specifically stresses Box x Box / Cylinder x Cylinder, the two pairs ComputeContacts resolves via
   // the generic GJK/EPA fallback (no analytic overload) rather than a closed form -- see
   // ParticleKernels.h's own comment on SolidCylinderContactDetection.
   const real_t e = real_t(1), mu = real_t(0);

   for (const auto& pr : detail::allShapePairs())
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), separation + real_t(20));
      const size_t idxA = detail::createPairMember(engine, pr.first, Vec3(-separation / 2, 0, 0),
                                                     Vec3(approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      const size_t idxB = detail::createPairMember(engine, pr.second, Vec3(separation / 2, 0, 0),
                                                     Vec3(-approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      engine.finalizeSetup();

      const std::string tag = outBase + "_" + pr.first + "_" + pr.second;
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { idxA, idxB });
      engine.setupVtk(tag, vtkSpacing);

      const real_t keBefore = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      const Vec3   pBefore  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                              engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);

      // Live cross-check of the not-yet-simulation-wired closed forms (PaperAnalyticalContacts.h)
      // against the production GJK/EPA geometry actually driving the dynamics this step, for the two
      // pairs that have a closed form there (box-box, cylinder-cylinder). Diagnostic only -- does not
      // change what drives the physics, per that file's own stated intent.
      const bool crossCheckable =
         (pr.first == "box" && pr.second == "box") || (pr.first == "cylinder" && pr.second == "cylinder");
      real_t maxOverlapDeltaSeen = real_t(0);

      const uint_t maxSteps = uint_t(2) * uint_t(separation / (approachSpeed * dt)) +
                               uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      bool   everContacted = false;
      real_t penetration   = real_t(0);
      for (uint_t s = 0; s < maxSteps; ++s)
      {
         const bool contactNow = engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0);
         if (contactNow)
         {
            everContacted = true;
            if (crossCheckable)
            {
               real_t paperOverlap = real_t(0);
               if (pr.first == "box")
               {
                  const auto res = paper_analytical::boxBoxSAT(
                     engine.ac().getPosition(idxA), engine.ac().getRotation(idxA).getMatrix(), boxHalfExtents,
                     engine.ac().getPosition(idxB), engine.ac().getRotation(idxB).getMatrix(), boxHalfExtents);
                  if (res.contact) paperOverlap = res.overlap;
               }
               else
               {
                  const paper_analytical::CylinderState csA{ engine.ac().getPosition(idxA),
                                                               engine.ac().getRotation(idxA).getMatrix(), cylRadius,
                                                               cylLength };
                  const paper_analytical::CylinderState csB{ engine.ac().getPosition(idxB),
                                                               engine.ac().getRotation(idxB).getMatrix(), cylRadius,
                                                               cylLength };
                  const auto res = paper_analytical::cylinderCylinderContact(csA, csB);
                  if (res.contact) paperOverlap = res.overlap;
               }
               // paper_analytical uses the OPPOSITE overlap sign convention (positive = penetrating);
               // SolidCylinderContactDetection's penetrationDepth is negative = penetrating.
               maxOverlapDeltaSeen = std::max(maxOverlapDeltaSeen, std::abs(paperOverlap - (-penetration)));
            }
         }
         engine.step();
         engine.writeCsvRow();
         engine.writeVtkIfDue();
         if (everContacted && !contactNow && s > 0) break;
      }

      const real_t keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      const Vec3   pAfter  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                             engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);
      const real_t relErr  = std::abs(keAfter - keBefore) / std::max(real_t(1e-12), keBefore);
      const bool   pass    = everContacted && relErr < tolerance;
      acc.record(pass);
      detail::logMomentumConservation(std::string("energyConservation ") + pr.first + "-" + pr.second, pBefore,
                                       pAfter, acc);

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] energyConservation "
                                 << pr.first << "-" << pr.second << "  KE_before=" << keBefore
                                 << " KE_after=" << keAfter << " relErr=" << relErr
                                 << (crossCheckable ? "  [paperAnalyticalMaxOverlapDelta=" +
                                                          std::to_string(maxOverlapDeltaSeen) + "]"
                                                    : "")
                                 << "  " << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] energyConservation OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}


} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
