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
//! \file DEMValidationPackingAndDrum.h
//! \brief Packing, angle-of-repose, and rotating-drum DEM validation scenarios: runPacking, runAngleOfRepose,
//! runPackingSingleShape, runRotatingDrum.

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
// 6-7: bulk mixed-shape tests
//======================================================================================================================

inline bool runPacking(const config::Config::BlockHandle& b)
{
   const uint_t numPerShape    = b.getParameter< uint_t >("numPerShape", uint_t(30));
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(0.9, 0.9, 0.9));
   const real_t cylRadius      = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength      = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t e              = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu             = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity        = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t containerX     = b.getParameter< real_t >("containerSizeX", real_t(20));
   const real_t containerY     = b.getParameter< real_t >("containerSizeY", real_t(20));
   const uint_t maxSteps       = b.getParameter< uint_t >("maxSteps", uint_t(50000));
   const uint_t vtkSpacing     = b.getParameter< uint_t >("vtkSpacing", uint_t(500));
   const std::string outBase   = b.getParameter< std::string >("outputBaseName", "packing");

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity, real_t(2) * (containerX + containerY));
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1));
   engine.addWall(Vec3(containerX, 0, 0), Vec3(-1, 0, 0));
   engine.addWall(Vec3(0, 0, 0), Vec3(1, 0, 0));
   engine.addWall(Vec3(0, containerY, 0), Vec3(0, -1, 0));
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 1, 0));

   real_t       totalVolume = real_t(0);
   uint_t       created     = 0;
   const real_t dropHeight  = real_t(1);
   const real_t cellSize    = real_t(2.2) * std::max({ sphereRadius, boxHalfExtents[0], boxHalfExtents[1],
                                                        boxHalfExtents[2],
                                                        std::sqrt(cylRadius * cylRadius +
                                                                   real_t(0.25) * cylLength * cylLength) });
   const uint_t perRow      = std::max< uint_t >(1, uint_t(containerX / cellSize));

   for (uint_t shapeType = 0; shapeType < 3; ++shapeType)
      for (uint_t i = 0; i < numPerShape; ++i)
      {
         const uint_t row = created % perRow, col = (created / perRow) % perRow;
         const real_t x =
            cellSize * (real_t(row) + real_t(0.5)) + math::realRandom< real_t >(-real_t(0.1), real_t(0.1));
         const real_t y =
            cellSize * (real_t(col) + real_t(0.5)) + math::realRandom< real_t >(-real_t(0.1), real_t(0.1));
         const real_t z = dropHeight + (real_t(created) / real_t(perRow * perRow)) * cellSize;
         const Vec3   pos(std::min(std::max(x, real_t(0.5)), containerX - real_t(0.5)),
                           std::min(std::max(y, real_t(0.5)), containerY - real_t(0.5)), z);

         if (shapeType == 0)
         {
            engine.addSphere(pos, Vec3(0), sphereRadius, density);
            totalVolume += (real_t(4) / real_t(3)) * math::pi * sphereRadius * sphereRadius * sphereRadius;
         }
         else if (shapeType == 1)
         {
            const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                                 math::realRandom< real_t >(-1, 1))
                              .getNormalizedIfNotZero(),
                           math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
            engine.addBox(pos, Vec3(0), rot, real_t(2) * boxHalfExtents, density);
            totalVolume += real_t(8) * boxHalfExtents[0] * boxHalfExtents[1] * boxHalfExtents[2];
         }
         else
         {
            const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                                 math::realRandom< real_t >(-1, 1))
                              .getNormalizedIfNotZero(),
                           math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
            engine.addCylinder(pos, Vec3(0), rot, cylRadius, cylLength, density);
            totalVolume += math::pi * cylRadius * cylRadius * cylLength;
         }
         ++created;
      }

   engine.finalizeSetup();
   engine.setupVtk(outBase, vtkSpacing);

   real_t avgSpeed          = std::numeric_limits< real_t >::max();
   real_t prevEnergy        = std::numeric_limits< real_t >::max();
   real_t worstEnergyIncrease = real_t(0); // largest checkpoint-to-checkpoint relative increase seen
   uint_t step              = 0;
   for (; step < maxSteps && avgSpeed > real_t(1e-3) * std::sqrt(gravity * sphereRadius); ++step)
   {
      engine.step();
      engine.writeVtkIfDue();
      // step > 0 guard: VelocityVerletPostForceUpdate averages oldForce and newForce, and oldForce is zero on
      // the first step, so avgSpeed at step==0 reflects about HALF a gravity kick and could read as "settled".
      if (step > 0 && step % 200 == 0)
      {
         real_t sumSpeed = real_t(0);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            sumSpeed += engine.ac().getLinearVelocity(i).length();
         }
         avgSpeed = sumSpeed / real_t(created);

         const real_t energy = detail::totalMechanicalEnergy(engine, gravity);
         if (prevEnergy < std::numeric_limits< real_t >::max())
         {
            const real_t increase = (energy - prevEnergy) / std::max(real_t(1e-9), std::abs(prevEnergy));
            worstEnergyIncrease = std::max(worstEnergyIncrease, increase);
         }
         prevEnergy = energy;

         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packing step=" << step << " avgSpeed=" << avgSpeed
                                    << " totalMechanicalEnergy=" << energy);
      }
   }

   // Checkpoint-to-checkpoint tolerance (not step-to-step): a sub-cycled spring-dashpot contact is not exactly
   // energy-monotonic on every timestep, but should show no sustained upward trend sampled every 200 steps.
   const bool energyMonotonic = worstEnergyIncrease < real_t(0.01);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packing [energy non-increase] worst checkpoint-to-checkpoint "
                              << "relative increase=" << worstEnergyIncrease << "  "
                              << (energyMonotonic ? "PASS" : "FAIL (energy created -- contact/integration unstable)"));

   real_t maxTop = real_t(0);
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      maxTop = std::max(maxTop, engine.ac().getPosition(i)[2] + engine.ac().getInteractionRadius(i));
   }

   // solid fraction = totalVolume / (containerX*containerY*maxTop), checked against the band below
   const real_t bulkSolidFraction = totalVolume / (containerX * containerY * std::max(maxTop, real_t(1e-6)));
   // Literature ordering for random-packed basic 3D shapes (Li, S., Zhao, J., Lu, P., & Xie, Y. (2010),
   // "Maximum packing densities of basic 3D objects," Chinese Science Bulletin 55(2), 114-119,
   // doi:10.1007/s11434-009-0650-0): the MAXIMALLY-JAMMED upper bounds are sphere~0.64, cylinder~0.72,
   // cube~0.78. A gravity-poured bed settles BELOW these maxima, so the pass band below is a broad plausibility
   // check ("not obviously broken"), not an exact-match test.
   // NOT gated on step < maxSteps: a nonzero-restitution mixed-shape bed settles into a persistent low-level dynamic
   // equilibrium (particles micro-resettling against irregular neighbours) rather than zero velocity -- avgSpeed
   // plateaus at a small, non-decaying value well before maxSteps -- so the strict threshold would only measure
   // "ran out of budget". maxSteps reaches a representative quasi-steady state; the bulk-fraction band is the check.
   const bool pass = bulkSolidFraction > real_t(0.30) && bulkSolidFraction < real_t(0.80) && energyMonotonic;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packing steps="
                              << step << " settledHeight=" << maxTop << " bulkSolidFraction=" << bulkSolidFraction
                              << "  (literature RCP upper bounds: sphere~0.64, cylinder~0.72, cube~0.78)"
                              << "  " << (pass ? "PASS" : "FAIL"));
   return pass;
}

inline bool runAngleOfRepose(const config::Config::BlockHandle& b)
{
   const uint_t numParticles    = b.getParameter< uint_t >("numParticles", uint_t(60));
   const real_t sphereRadius    = b.getParameter< real_t >("sphereRadius", real_t(1));
   const real_t cylRadius       = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength       = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const bool   useCylinders    = b.getParameter< bool >("useCylinders", true);
   const real_t density         = b.getParameter< real_t >("density", real_t(1));
   const real_t e               = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu              = b.getParameter< real_t >("frictionCoefficient", real_t(0.6));
   const real_t poisson         = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime   = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt              = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity         = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t floorHalfExtent = b.getParameter< real_t >("floorHalfExtent", real_t(30));
   const uint_t maxSteps        = b.getParameter< uint_t >("maxSteps", uint_t(60000));
   const uint_t vtkSpacing      = b.getParameter< uint_t >("vtkSpacing", uint_t(500));
   const std::string outBase    = b.getParameter< std::string >("outputBaseName", "angle_of_repose");

   // Domain must be tall enough for the INITIAL POUR COLUMN (see the placement loop below), not just wide enough
   // for the floor footprint; a pour column poking outside the domain trips AssocToBlock's block-assignment
   // assertion.
   //
   // A random (r,phi) pour column seats particles overlapping at creation; LinearSpringDashpot's force scales with
   // penetration depth, so that overlap is an explosive first-step force (a numerical blow-up to r,z in the
   // thousands, which averaging can make read as "flat"). Hence a regular grid (cell size with a safety margin
   // over the largest particle's circumscribing radius, as in runPacking), so no two particles start in contact.
   const real_t cellSize =
      real_t(2.2) * std::max(sphereRadius, std::sqrt(cylRadius * cylRadius + real_t(0.25) * cylLength * cylLength));
   const uint_t gridSide          = uint_t(4); // 4x4 = 16 particles per layer
   const real_t pourColumnHeight  = real_t(2) + cellSize * real_t(numParticles / (gridSide * gridSide) + 1) + real_t(2);
   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity,
                     std::max(real_t(3) * floorHalfExtent, pourColumnHeight + real_t(5)));
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1));

   // Even with a clean pour, the pile stays essentially flat on a fully open floor: mesa_pd's LinearSpringDashpot
   // models pure Coulomb sliding friction with no rolling-resistance torque, so a sphere or rolling cylinder forms a
   // bounded pile only if it cannot roll away (a well-documented limitation of plain-Coulomb DEM contact models).
   // Containing the pour in a bounded pit (the wall pattern of runPacking) sidesteps it, closer to the "pour into a
   // hopper, measure the settled slope" method of the DEM calibration literature than Kodam's free-standing cone.
   const real_t pitHalfExtent = b.getParameter< real_t >("pitHalfExtent", real_t(10));
   engine.addWall(Vec3(pitHalfExtent, 0, 0), Vec3(-1, 0, 0));
   engine.addWall(Vec3(-pitHalfExtent, 0, 0), Vec3(1, 0, 0));
   engine.addWall(Vec3(0, pitHalfExtent, 0), Vec3(0, -1, 0));
   engine.addWall(Vec3(0, -pitHalfExtent, 0), Vec3(0, 1, 0));

   // Pour from a small column above the centre -- the standard static-pile method of DEM calibration studies
   // for measuring the angle of repose.
   uint_t created = 0;
   for (uint_t i = 0; i < numParticles; ++i)
   {
      const uint_t row = (i / gridSide) % gridSide, col = i % gridSide, layer = i / (gridSide * gridSide);
      const real_t gx  = cellSize * (real_t(row) - real_t(gridSide - 1) * real_t(0.5)) +
                         math::realRandom< real_t >(-real_t(0.05), real_t(0.05));
      const real_t gy = cellSize * (real_t(col) - real_t(gridSide - 1) * real_t(0.5)) +
                        math::realRandom< real_t >(-real_t(0.05), real_t(0.05));
      const Vec3 pos(gx, gy, real_t(2) + cellSize * real_t(layer));
      // ALL cylinders when useCylinders (not a 50/50 mix with spheres): without a rolling-resistance torque (see
      // above), a sphere (or a short, disk-like cylinder) can roll indefinitely instead of settling into a bounded
      // pile -- spherical/near-spherical particles need an explicit rolling-friction model for a realistic angle of
      // repose. An elongated rod (see cylLength's aspect ratio) resists end-over-end rolling far better.
      if (useCylinders)
      {
         const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                              math::realRandom< real_t >(-1, 1))
                           .getNormalizedIfNotZero(),
                        math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
         engine.addCylinder(pos, Vec3(0), rot, cylRadius, cylLength, density);
      }
      else { engine.addSphere(pos, Vec3(0), sphereRadius, density); }
      ++created;
   }

   engine.finalizeSetup();
   engine.setupVtk(outBase, vtkSpacing);

   real_t avgSpeed            = std::numeric_limits< real_t >::max();
   real_t prevEnergy          = std::numeric_limits< real_t >::max();
   real_t worstEnergyIncrease = real_t(0);
   uint_t step                = 0;
   for (; step < maxSteps && avgSpeed > real_t(1e-3) * std::sqrt(gravity * sphereRadius); ++step)
   {
      engine.step();
      engine.writeVtkIfDue();
      // step > 0 guard: VelocityVerletPostForceUpdate averages oldForce and newForce, and oldForce is zero on
      // the first step, so avgSpeed at step==0 reflects about HALF a gravity kick and could read as "settled".
      if (step > 0 && step % 200 == 0)
      {
         real_t sumSpeed = real_t(0);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            sumSpeed += engine.ac().getLinearVelocity(i).length();
         }
         avgSpeed = sumSpeed / real_t(created);

         const real_t energy = detail::totalMechanicalEnergy(engine, gravity);
         if (prevEnergy < std::numeric_limits< real_t >::max())
         {
            const real_t increase = (energy - prevEnergy) / std::max(real_t(1e-9), std::abs(prevEnergy));
            worstEnergyIncrease   = std::max(worstEnergyIncrease, increase);
         }
         prevEnergy = energy;
      }
   }

   const bool energyMonotonic = worstEnergyIncrease < real_t(0.01);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] angleOfRepose [energy non-increase] worst checkpoint-to-checkpoint "
                              << "relative increase=" << worstEnergyIncrease << "  "
                              << (energyMonotonic ? "PASS" : "FAIL (energy created -- contact/integration unstable)"));

   // Least-squares fit of height vs. radial distance from the pour axis over the settled pile's outer
   // shell -- slope magnitude approximates tan(angle of repose).
   real_t sumR = 0, sumZ = 0, sumRZ = 0, sumRR = 0;
   real_t minR = std::numeric_limits< real_t >::max(), maxR = real_t(0);
   real_t minZ = std::numeric_limits< real_t >::max(), maxZ = real_t(0);
   uint_t n = 0;
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      const Vec3   pos = engine.ac().getPosition(i);
      const real_t r   = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      sumR += r;
      sumZ += pos[2];
      sumRZ += r * pos[2];
      sumRR += r * r;
      minR = std::min(minR, r);
      maxR = std::max(maxR, r);
      minZ = std::min(minZ, pos[2]);
      maxZ = std::max(maxZ, pos[2]);
      ++n;
   }
   // Sanity-check range: a genuinely settled pile should have maxR/maxZ of the same order as the pour
   // footprint/height, not orders of magnitude larger (the latter would indicate a numerical blow-up,
   // e.g. particles seeded overlapping at creation, rather than a physically settled result).
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] angleOfRepose pile extent: n=" << n << " rRange=[" << minR << ","
                              << maxR << "] zRange=[" << minZ << "," << maxZ << "]");
   const real_t denom = real_t(n) * sumRR - sumR * sumR;
   const real_t slope = (std::abs(denom) > real_t(1e-9)) ? (real_t(n) * sumRZ - sumR * sumZ) / denom : real_t(0);
   // repose angle = atan(|slope|) of the pile's r-vs-z least-squares fit; checked against the
   // plausibility band below, anchored by Kodam et al. 2010 Part II's 39.6deg sim / 43.9deg experiment
   const real_t reposeAngleDeg = std::atan(std::abs(slope)) * real_t(180) / math::pi;

   // Plausibility band only (~15-55 deg spans essentially all real granular materials); Kodam et al.
   // (2010) Part II's cylinder result (39.6 deg simulated vs. 43.9 deg experimental, aspect ratio
   // ~0.66) is one concrete literature anchor within this band. NOT gated on step < maxSteps -- see
   // runPacking's identical reasoning: a nonzero-restitution mixed-shape pile settles into a persistent
   // low-level dynamic equilibrium, not exact rest.
   const bool pass = reposeAngleDeg > real_t(15) && reposeAngleDeg < real_t(55) && energyMonotonic;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] angleOfRepose steps="
                              << step << " reposeAngle=" << reposeAngleDeg
                              << "deg  (plausible range 15-55deg; Kodam et al. 2010 cylinder reference: "
                              << "39.6deg sim / 43.9deg experiment)  " << (pass ? "PASS" : "FAIL"));
   if (!pass && reposeAngleDeg < real_t(5))
      WALBERLA_LOG_INFO_ON_ROOT(
         "[DEMValidation] angleOfRepose DIAGNOSIS: a near-zero angle after settling (not just below "
         "tolerance) points to the known mesa_pd LinearSpringDashpot limitation documented above this "
         "function -- pure Coulomb sliding friction with no rolling-resistance torque lets particles roll/"
         "slide out to a flat, low layer rather than piling up, even when contained. This is a genuine "
         "contact-model gap this validation run has surfaced, not a broken test.");
   return pass;
}

//======================================================================================================================
// 8-9: single-shape Kodam et al. (2010) Part II Fig. 10/11 methodology recreations -- unlike runPacking/
// runAngleOfRepose above (generalised/simplified variants), these follow Kodam's protocol as closely as a
// fluid-free DEM dry-run can: ONE shape at a time (`shape` prm key: "sphere"/"box"/"cylinder"), a circular
// container/drum wall (polygon-of-planes approximation, the "curved wall from flat plates" idea of their Fig. 6(a)
// friability tester), particles with a random initial orientation, and -- for the repose test -- a ROTATING drum
// (DEMEngine::setWallPose every tumbling step) rather than a static pour-pile.
//======================================================================================================================

namespace detail
{
//! Creates one particle of `shape` ("sphere"/"box"/"cylinder") at `pos`, at rest, with a uniformly random
//! 3D orientation (Kodam et al.'s "random initial orientation" drop protocol); skipped for orientation-invariant
//! spheres.
inline size_t addRandomOrientedShape(DEMEngine& engine, const std::string& shape, const Vec3& pos,
                                      real_t sphereRadius, const Vec3& boxHalfExtents, real_t cylRadius,
                                      real_t cylLength, real_t density)
{
   if (shape == "sphere") return engine.addSphere(pos, Vec3(0), sphereRadius, density);
   const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                        math::realRandom< real_t >(-1, 1))
                     .getNormalizedIfNotZero(),
                  math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
   if (shape == "box") return engine.addBox(pos, Vec3(0), rot, real_t(2) * boxHalfExtents, density);
   if (shape == "cylinder") return engine.addCylinder(pos, Vec3(0), rot, cylRadius, cylLength, density);
   WALBERLA_ABORT("Unknown shape '" << shape << "'. Valid values: sphere, box, cylinder.");
}

//! Circumscribing-sphere radius for `shape` -- used for non-overlapping grid spacing and fill-height
//! measurement (same convention as DEMEngine::addSphere/addBox/addCylinder's setInteractionRadius).
inline real_t shapeCircumRadius(const std::string& shape, real_t sphereRadius, const Vec3& boxHalfExtents,
                                 real_t cylRadius, real_t cylLength)
{
   if (shape == "sphere") return sphereRadius;
   if (shape == "box") return boxHalfExtents.length();
   if (shape == "cylinder") return std::sqrt(cylRadius * cylRadius + real_t(0.25) * cylLength * cylLength);
   WALBERLA_ABORT("Unknown shape '" << shape << "'. Valid values: sphere, box, cylinder.");
}

//! Exact solid volume for `shape` (used for bulk solid fraction, Kodam et al. Eq. (7)).
inline real_t shapeVolume(const std::string& shape, real_t sphereRadius, const Vec3& boxHalfExtents,
                           real_t cylRadius, real_t cylLength)
{
   if (shape == "sphere") return (real_t(4) / real_t(3)) * math::pi * sphereRadius * sphereRadius * sphereRadius;
   if (shape == "box") return real_t(8) * boxHalfExtents[0] * boxHalfExtents[1] * boxHalfExtents[2];
   if (shape == "cylinder") return math::pi * cylRadius * cylRadius * cylLength;
   WALBERLA_ABORT("Unknown shape '" << shape << "'. Valid values: sphere, box, cylinder.");
}

//! Builds a regular `numSegments`-sided polygon of flat HalfSpace walls approximating a cylindrical
//! container/drum wall of interior (apothem) radius `radius`, axis along world Y, centred at
//! (axisX, *, axisZ) in the X-Z plane. `angleOffset` rotates the whole polygon; passing a time-varying
//! offset (see updateDrumWalls) animates a rotating drum. Each wall is tangent to the target circle at
//! its segment's angle, so particles are contained within `radius` at the flats and bulge very slightly
//! beyond it at the corners (by a factor 1/cos(pi/numSegments), negligible for numSegments>=16).
inline std::vector< size_t > addDrumWalls(DEMEngine& engine, real_t axisX, real_t axisZ, real_t radius,
                                           uint_t numSegments, real_t angleOffset = real_t(0))
{
   std::vector< size_t > walls;
   walls.reserve(numSegments);
   for (uint_t i = 0; i < numSegments; ++i)
   {
      const real_t theta = real_t(2) * math::pi * real_t(i) / real_t(numSegments) + angleOffset;
      const Vec3   normal(-std::cos(theta), real_t(0), -std::sin(theta)); // inward
      const Vec3   pos(axisX + radius * std::cos(theta), real_t(0), axisZ + radius * std::sin(theta));
      walls.push_back(engine.addWall(pos, normal));
   }
   return walls;
}

//! Re-poses an existing set of drum walls (from addDrumWalls) to a new `angleOffset`, keeping the same
//! axis/radius/segment count -- called every step during the tumbling phase of runRotatingDrum.
inline void updateDrumWalls(DEMEngine& engine, const std::vector< size_t >& walls, real_t axisX, real_t axisZ,
                             real_t radius, real_t angleOffset)
{
   const uint_t numSegments = uint_t(walls.size());
   for (uint_t i = 0; i < numSegments; ++i)
   {
      const real_t theta = real_t(2) * math::pi * real_t(i) / real_t(numSegments) + angleOffset;
      const Vec3   normal(-std::cos(theta), real_t(0), -std::sin(theta));
      const Vec3   pos(axisX + radius * std::cos(theta), real_t(0), axisZ + radius * std::sin(theta));
      engine.setWallPose(walls[i], pos, normal);
   }
}
} // namespace detail

//! Kodam et al. (2010) Fig. 10 / Section 3.3 recreation: ONE shape (swap via `shape`), released one at a
//! time from a random point near the top with negligible initial KE and a random orientation, into a
//! genuinely circular container (detail::addDrumWalls, static), each particle waited out (settle-speed
//! threshold or dropTimeoutSteps) before the next is released -- Kodam's protocol keeps the drop rate
//! far below 1000/s so that the drop rate does not bias the packing (Macrae & Gray, 1961; Zhang
//! et al., 2001, both cited in the paper); one-at-a-time is the limiting case of that. Measures fill
//! height and bulk solid fraction (their Eq. (7)) against shape-specific plausibility bands anchored to
//! the random-close-packing literature (Li, Zhao, Lu & Xie, 2010, Chinese Science Bulletin 55(2):114-119)
//! and Kodam's measured cylinder values (Table 2: 0.617 experiment / 0.589 true-cylinder-DEM sim).
inline bool runPackingSingleShape(const config::Config::BlockHandle& b)
{
   const std::string shape        = b.getParameter< std::string >("shape", "cylinder");
   const uint_t       numParticles = b.getParameter< uint_t >("numParticles", uint_t(250));
   const real_t sphereRadius      = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents    = b.getParameter< Vec3 >("boxHalfExtents", Vec3(0.9, 0.9, 0.9));
   const real_t cylRadius         = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength         = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const real_t density           = b.getParameter< real_t >("density", real_t(1));
   const real_t e                 = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu                = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson           = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime     = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt                = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity           = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t containerRadius   = b.getParameter< real_t >("containerRadius", real_t(6));
   const uint_t numWallSegments   = b.getParameter< uint_t >("numWallSegments", uint_t(24));
   const uint_t settleCheckEvery  = b.getParameter< uint_t >("settleCheckEvery", uint_t(200));
   const uint_t dropTimeoutSteps  = b.getParameter< uint_t >("dropTimeoutSteps", uint_t(6000));
   const uint_t maxSteps          = b.getParameter< uint_t >("maxSteps", uint_t(600000));
   const uint_t vtkSpacing        = b.getParameter< uint_t >("vtkSpacing", uint_t(200));
   const std::string outBase      = b.getParameter< std::string >("outputBaseName", "packing_single_shape");

   const real_t circumRadius   = detail::shapeCircumRadius(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);
   const real_t particleVolume = detail::shapeVolume(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);

   // Domain must contain the tallest plausible settled pile PLUS the drop height above it, sized from the actual
   // particle count/volume so it holds across shape/count combinations (an under-sized domain fails as noted in
   // runAngleOfRepose).
   const real_t estPileHeight = real_t(1.6) * particleVolume * real_t(numParticles) /
                                 (math::pi * containerRadius * containerRadius);
   const real_t dropHeight     = real_t(2) * containerRadius + real_t(2) * circumRadius + estPileHeight;
   // Margin scaled to particle size, not a fixed constant: at Kodam's SI scale (containerRadius ~2.5e-2 m) a fixed
   // "+10" margin would exceed the physical domain ~1000x, and step()'s LinkedCells (cell width from the mm-scale
   // particle diameter) would allocate on the order of 1e10 cells.
   const real_t domainHalfExtent = dropHeight + real_t(2) * containerRadius + real_t(4) * circumRadius;

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity, domainHalfExtent);
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1)); // floor
   detail::addDrumWalls(engine, real_t(0), real_t(0), containerRadius, numWallSegments);

   const real_t dropRadius = real_t(0.6) * containerRadius; // random release point, away from the wall

   // No dynamic particles exist yet (only GLOBAL walls): pass the known max diameter explicitly (see
   // DEMEngine::finalizeSetup) so the linked-cell grid fits the particles the drop loop creates.
   engine.finalizeSetup(real_t(2) * circumRadius);
   engine.setupVtk(outBase, vtkSpacing);

   uint_t step    = 0;
   uint_t created = 0;
   for (; created < numParticles && step < maxSteps; ++created)
   {
      // Release height tracks the current pile top (with margin): a fixed drop height would eventually spawn a
      // particle INSIDE the growing pile, an explosive initial-overlap contact force (see runAngleOfRepose).
      real_t currentTop = real_t(0);
      for (size_t i = 0; i < engine.ac().size(); ++i)
      {
         if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
         currentTop = std::max(currentTop, engine.ac().getPosition(i)[2] + engine.ac().getInteractionRadius(i));
      }
      const real_t dropZ = std::max(real_t(2) * containerRadius, currentTop + real_t(3) * circumRadius);

      const real_t r   = dropRadius * std::sqrt(math::realRandom< real_t >(real_t(0), real_t(1)));
      const real_t phi = math::realRandom< real_t >(real_t(0), real_t(2) * math::pi);
      detail::addRandomOrientedShape(engine, shape, Vec3(r * std::cos(phi), real_t(0), dropZ), sphereRadius,
                                      boxHalfExtents, cylRadius, cylLength, density);

      // Wait for this particle (and whatever it disturbs) to settle before releasing the next one (Kodam et al.'s
      // one-at-a-time protocol). dropTimeoutSteps is a safety valve: a particle in a stable-but-nonzero
      // micro-equilibrium (rocking against neighbours) must not stall the sequence.
      for (uint_t waited = 0; waited < dropTimeoutSteps && step < maxSteps; ++waited, ++step)
      {
         engine.step();
         engine.writeVtkIfDue();
         if (waited > 0 && waited % settleCheckEvery == 0)
         {
            real_t sumSpeed = real_t(0);
            uint_t n         = uint_t(0);
            for (size_t i = 0; i < engine.ac().size(); ++i)
            {
               if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
               sumSpeed += engine.ac().getLinearVelocity(i).length();
               ++n;
            }
            if (n > 0 && sumSpeed / real_t(n) < real_t(1e-2) * std::sqrt(gravity * circumRadius)) break;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packingSingleShape shape=" << shape << " created=" << created << "/"
                              << numParticles << " steps=" << step);

   real_t maxTop = real_t(0);
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      maxTop = std::max(maxTop, engine.ac().getPosition(i)[2] + engine.ac().getInteractionRadius(i));
   }

   // solid fraction = N*V_particle / (pi*R_container^2*maxTop), checked against the band below
   const real_t bulkSolidFraction = real_t(created) * particleVolume /
                                     (math::pi * containerRadius * containerRadius * std::max(maxTop, real_t(1e-6)));

   // Shape-specific plausibility bands, anchored to the literature RCP upper bounds (Li, Zhao, Lu & Xie, 2010,
   // Chinese Science Bulletin 55(2):114-119, doi:10.1007/s11434-009-0650-0 -- sphere~0.64, cylinder~0.72,
   // cube~0.78 maximally jammed; not verified against a local copy, so a soft, non-gating ceiling only) and Kodam et
   // al.'s measured cylinder result (0.617 experiment / 0.589 true-cylinder-DEM sim, Table 2, verified against
   // kodam_cylindrical_2010b.txt) as a concrete anchor. A gravity-poured bed sits below the jammed maxima, so these
   // are broad "not obviously broken" bands: this dry run's contact parameters are not calibrated to Kodam's
   // Plexiglas-on-Plexiglas values, so exact numeric reproduction is not the goal.
   real_t bandHi = real_t(0.64);
   if (shape == "cylinder") bandHi = real_t(0.72);
   else if (shape == "box") bandHi = real_t(0.78);
   const bool pass = created == numParticles && bulkSolidFraction > real_t(0.30) && bulkSolidFraction < bandHi;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packingSingleShape shape="
                              << shape << " settledHeight=" << maxTop << " bulkSolidFraction=" << bulkSolidFraction
                              << " (plausible band [0.30," << bandHi << "])"
                              << "  (Kodam et al. 2010 Table 2 cylinder reference: 0.617 experiment / 0.589 "
                                 "true-cylinder-DEM sim)  "
                              << (pass ? "PASS" : "FAIL"));
   return pass;
}

//! Kodam et al. (2010) Fig. 11 / Section 3.4 recreation: ONE shape (swap via `shape`) tumbling inside an
//! ACTUALLY ROTATING drum (detail::addDrumWalls + updateDrumWalls, animated every tumbling step at `rpm`
//! -- not runAngleOfRepose's static pour-pile), axis along world Y, bounded by two flat endcap walls at
//! +-drumHalfLength. Particles fill the drum's interior directly (a 3D grid clipped to the circular cross-section,
//! non-overlapping by construction; gravity sorts them into a resting bed during the settle phase). After settling
//! with the drum stationary, the drum rotates for `tumbleRevolutions`
//! revolutions; `numSnapshots` free-surface readings are taken across that window (Kodam et al.: 10
//! photos at 1s intervals during their real 25rpm test) -- each reading bins particles by X into
//! `numRadialBins` bins and records the topmost occupied Z per occupied bin, the same "profile from
//! images" idea as their Fig. 12(a) pixel extraction, sampled directly from particle geometry instead of
//! pixels. All snapshot points are pooled and least-squares line-fit (their Fig. 12(b) method) to get one
//! dynamic angle of repose.
inline bool runRotatingDrum(const config::Config::BlockHandle& b)
{
   const std::string shape         = b.getParameter< std::string >("shape", "cylinder");
   const uint_t       numParticles  = b.getParameter< uint_t >("numParticles", uint_t(100));
   const real_t sphereRadius       = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents     = b.getParameter< Vec3 >("boxHalfExtents", Vec3(0.9, 0.9, 0.9));
   const real_t cylRadius          = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength          = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const real_t density            = b.getParameter< real_t >("density", real_t(1));
   const real_t e                  = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu                 = b.getParameter< real_t >("frictionCoefficient", real_t(0.6));
   const real_t poisson            = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime      = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt                 = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity            = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t drumRadius         = b.getParameter< real_t >("drumRadius", real_t(10));
   const real_t drumHalfLength     = b.getParameter< real_t >("drumHalfLength", real_t(5.6));
   const uint_t numWallSegments    = b.getParameter< uint_t >("numWallSegments", uint_t(24));
   const real_t rpm                = b.getParameter< real_t >("rpm", real_t(25));       // Kodam et al.'s value
   const uint_t settleSteps        = b.getParameter< uint_t >("settleSteps", uint_t(30000));
   const uint_t tumbleRevolutions  = b.getParameter< uint_t >("tumbleRevolutions", uint_t(3));
   const uint_t numSnapshots       = b.getParameter< uint_t >("numSnapshots", uint_t(10)); // Kodam: 10 photos
   const uint_t numRadialBins      = b.getParameter< uint_t >("numRadialBins", uint_t(20));
   const uint_t vtkSpacing         = b.getParameter< uint_t >("vtkSpacing", uint_t(200));
   const std::string outBase       = b.getParameter< std::string >("outputBaseName", "rotating_drum");

   const real_t circumRadius = detail::shapeCircumRadius(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);
   const real_t omega        = rpm * real_t(2) * math::pi / real_t(60); // rad/s

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity, real_t(3) * (drumRadius + drumHalfLength));
   auto wallIndices = detail::addDrumWalls(engine, real_t(0), real_t(0), drumRadius, numWallSegments);
   engine.addWall(Vec3(0, drumHalfLength, 0), Vec3(0, -1, 0));  // +Y endcap
   engine.addWall(Vec3(0, -drumHalfLength, 0), Vec3(0, 1, 0));  // -Y endcap

   // Fill the drum's interior directly with a 3D grid clipped to the circular cross-section (a real drum is
   // loaded closed). Non-overlapping by construction (cellSize has a safety margin over circumRadius), so any
   // subset works as a starting condition; gravity settles it into a resting bed during the stationary phase.
   const real_t cellSize       = real_t(2.1) * circumRadius;
   const real_t interiorRadius = std::max(real_t(0.1), drumRadius - real_t(1.05) * circumRadius);
   const real_t yMax           = std::max(real_t(0.1), drumHalfLength - real_t(1.05) * circumRadius);
   const uint_t nx              = uint_t(real_t(2) * interiorRadius / cellSize) + uint_t(1);
   const uint_t ny              = uint_t(real_t(2) * yMax / cellSize) + uint_t(1);

   uint_t created = 0;
   for (uint_t iz = 0; iz < nx && created < numParticles; ++iz)
   {
      const real_t gz = -interiorRadius + cellSize * (real_t(iz) + real_t(0.5));
      for (uint_t ix = 0; ix < nx && created < numParticles; ++ix)
      {
         const real_t gx = -interiorRadius + cellSize * (real_t(ix) + real_t(0.5));
         if (gx * gx + gz * gz > interiorRadius * interiorRadius) continue;
         for (uint_t iy = 0; iy < ny && created < numParticles; ++iy)
         {
            const real_t gy = (ny > 1) ? (-yMax + cellSize * (real_t(iy) + real_t(0.5))) : real_t(0);
            detail::addRandomOrientedShape(engine, shape, Vec3(gx, gy, gz), sphereRadius, boxHalfExtents, cylRadius,
                                            cylLength, density);
            ++created;
         }
      }
   }
   if (created < numParticles)
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum WARNING: only fit "
                                 << created << "/" << numParticles
                                 << " particles inside the drum at this grid resolution -- increase drumRadius/"
                                    "drumHalfLength or reduce numParticles/particle size for a larger bed.");

   // DIAGNOSTIC: verify empirically, before any physics runs, that the fill grid placed particles
   // non-overlapping (cellSize=2.1*circumRadius > 2*circumRadius, so no two grid points can be closer than
   // cellSize regardless of shape/orientation), to rule out spawn overlap as the cause of a settle-phase energy
   // spike. O(N^2) is cheap at N=250 (~31000 pairs), done once.
   {
      real_t minPairDist = std::numeric_limits< real_t >::max();
      size_t minI = 0, minJ = 0;
      for (size_t i = 0; i < engine.ac().size(); ++i)
      {
         if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
         for (size_t j = i + 1; j < engine.ac().size(); ++j)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(j), data::particle_flags::GLOBAL)) continue;
            const real_t d = (engine.ac().getPosition(i) - engine.ac().getPosition(j)).length();
            if (d < minPairDist) { minPairDist = d; minI = i; minJ = j; }
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum fill-grid check: minPairwiseCenterDistance="
                                 << minPairDist << " (expected >= cellSize=" << cellSize
                                 << ", i.e. >= 2*circumRadius=" << (real_t(2) * circumRadius) << ") "
                                 << (minPairDist >= real_t(2) * circumRadius ? "OK" : "*** OVERLAP AT SPAWN ***")
                                 << " pair=(" << minI << "," << minJ << ") pos_i=" << engine.ac().getPosition(minI)
                                 << " pos_j=" << engine.ac().getPosition(minJ));
   }

   engine.finalizeSetup();
   engine.setupVtk(outBase, vtkSpacing);

   // Phase 1: settle under gravity with the drum stationary (same avgSpeed/energy-monotonicity pattern as
   // runPacking/runAngleOfRepose above). This phase's energy check IS meaningful (no external driving
   // yet); the tumbling phase's is not (the rotating wall does net positive work on the bed by design),
   // so it is only checked here.
   real_t avgSpeed            = std::numeric_limits< real_t >::max();
   real_t prevEnergy          = std::numeric_limits< real_t >::max();
   real_t worstEnergyIncrease = real_t(0);
   uint_t step                = 0;
   for (; step < settleSteps && avgSpeed > real_t(1e-3) * std::sqrt(gravity * circumRadius); ++step)
   {
      engine.step();
      engine.writeVtkIfDue();
      if (step > 0 && step % 200 == 0)
      {
         real_t sumSpeed = real_t(0);
         real_t maxSpeed  = real_t(0);
         uint_t n         = uint_t(0);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            const real_t speed = engine.ac().getLinearVelocity(i).length();
            sumSpeed += speed;
            maxSpeed = std::max(maxSpeed, speed);
            ++n;
         }
         avgSpeed = (n > 0) ? sumSpeed / real_t(n) : real_t(0);

         const real_t energy = detail::totalMechanicalEnergy(engine, gravity);
         // Diagnostic trace for the settle-phase energy-monotonicity check (as runPacking's per-checkpoint print):
         // logs every checkpoint at INFO so a spike's step/energy/maxSpeed can be localised.
         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum settle step=" << step << " avgSpeed=" << avgSpeed
                                    << " maxSpeed=" << maxSpeed << " totalMechanicalEnergy=" << energy);
         if (prevEnergy < std::numeric_limits< real_t >::max())
         {
            const real_t increase = (energy - prevEnergy) / std::max(real_t(1e-9), std::abs(prevEnergy));
            if (increase > worstEnergyIncrease)
               WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum settle step=" << step
                                          << " NEW WORST relative energy increase=" << increase
                                          << " (prevEnergy=" << prevEnergy << " energy=" << energy << ")");
            worstEnergyIncrease = std::max(worstEnergyIncrease, increase);
         }
         prevEnergy = energy;
      }
   }
   const bool settlePhaseStable = worstEnergyIncrease < real_t(0.01);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum settle phase: steps=" << step << " [energy non-increase] "
                              << "worst checkpoint-to-checkpoint relative increase=" << worstEnergyIncrease << "  "
                              << (settlePhaseStable ? "PASS" : "FAIL (energy created -- contact/integration unstable)"));

   // Phase 2: rotate the drum at `rpm` for `tumbleRevolutions` full turns, taking `numSnapshots` free-
   // surface readings evenly spaced across that window.
   const real_t stepsPerRevolution = (real_t(2) * math::pi / omega) / dt;
   const uint_t tumbleSteps        = uint_t(real_t(tumbleRevolutions) * stepsPerRevolution);
   const uint_t snapshotEvery      = std::max< uint_t >(uint_t(1), tumbleSteps / std::max< uint_t >(1, numSnapshots));

   std::vector< std::pair< real_t, real_t > > surfacePoints;
   for (uint_t t = 0; t < tumbleSteps; ++t)
   {
      detail::updateDrumWalls(engine, wallIndices, real_t(0), real_t(0), drumRadius, omega * (real_t(t) * dt));
      engine.step();
      engine.writeVtkIfDue();

      if (t % snapshotEvery == 0)
      {
         std::vector< real_t > binMaxZ(numRadialBins, -std::numeric_limits< real_t >::max());
         std::vector< bool >   binHas(numRadialBins, false);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            const Vec3 pos = engine.ac().getPosition(i);
            if (std::abs(pos[0]) > drumRadius) continue;
            const uint_t bin = std::min< uint_t >(
               numRadialBins - 1,
               uint_t((pos[0] + drumRadius) / (real_t(2) * drumRadius) * real_t(numRadialBins)));
            const real_t topZ = pos[2] + engine.ac().getInteractionRadius(i);
            if (topZ > binMaxZ[bin])
            {
               binMaxZ[bin] = topZ;
               binHas[bin]  = true;
            }
         }
         for (uint_t bi = 0; bi < numRadialBins; ++bi)
            if (binHas[bi])
            {
               const real_t xCenter =
                  -drumRadius + (real_t(bi) + real_t(0.5)) * (real_t(2) * drumRadius / real_t(numRadialBins));
               surfacePoints.emplace_back(xCenter, binMaxZ[bi]);
            }
      }
   }

   // Least-squares fit of the pooled free-surface points -- Kodam et al.'s Fig. 12(b) method (average
   // the free surface over the ten images, then fit one line; here pooling plays the role of averaging).
   real_t sumX = 0, sumZ = 0, sumXZ = 0, sumXX = 0;
   uint_t n = 0;
   for (const auto& p : surfacePoints)
   {
      sumX += p.first;
      sumZ += p.second;
      sumXZ += p.first * p.second;
      sumXX += p.first * p.first;
      ++n;
   }
   const real_t denom          = real_t(n) * sumXX - sumX * sumX;
   const real_t slope          = (std::abs(denom) > real_t(1e-9)) ? (real_t(n) * sumXZ - sumX * sumZ) / denom : real_t(0);
   // repose angle = atan(|slope|) of the pooled free-surface points' least-squares fit; checked against
   // the plausibility band below, anchored by Kodam et al. 2010 Part II's 39.6deg sim / 43.9deg experiment
   const real_t reposeAngleDeg = std::atan(std::abs(slope)) * real_t(180) / math::pi;

   // Plausibility band only (~5-60deg spans essentially all real granular materials in a rotating drum);
   // Kodam et al. (2010) Part II's cylinder result (39.6deg simulated vs. 43.9deg experimental) is one
   // concrete literature anchor within this band.
   const bool pass = n > 0 && reposeAngleDeg > real_t(5) && reposeAngleDeg < real_t(60) && settlePhaseStable;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum shape="
                              << shape << " created=" << created << " tumbleSteps=" << tumbleSteps
                              << " snapshotPoints=" << n << " reposeAngle=" << reposeAngleDeg
                              << "deg  (plausible range 5-60deg; Kodam et al. 2010 cylinder reference: "
                              << "39.6deg sim / 43.9deg experiment)  " << (pass ? "PASS" : "FAIL"));
   if (!pass && n > 0 && reposeAngleDeg < real_t(5))
      WALBERLA_LOG_INFO_ON_ROOT(
         "[DEMValidation] rotatingDrum DIAGNOSIS: a near-zero angle points to the known mesa_pd "
         "LinearSpringDashpot limitation documented in runAngleOfRepose above -- pure Coulomb sliding "
         "friction with no rolling-resistance torque lets particles (especially spheres) roll/slide into a "
         "flat, low layer rather than a sloped bed, even inside a closed, rotating drum. This is a genuine "
         "contact-model gap this validation run has surfaced, not a broken test.");
   return pass;
}


} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
