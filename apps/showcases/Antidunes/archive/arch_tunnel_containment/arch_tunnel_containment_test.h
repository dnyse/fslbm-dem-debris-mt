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
//! \file arch_tunnel_containment_test.h
//! \brief ARCHIVED DEM-validation scenario for ComputeTunnelContacts (ParticleKernels.h) -- the
//! arch/tunnel pillar-containment check for the real-scale Pont Ferroviere bridge case study
//! (chapters/03_implementation/case-study/, chapters/04_results/case-study/ -- intentionally
//! excluded from the thesis, see CLAUDE.md). Moved out of DEMValidationScenarios.h 2026-10-29;
//! see archive/ARCHIVE_README.md. NOTE: ComputeTunnelContacts itself is NOT archived -- it is
//! still live production code in ParticleKernels.h, used by the much larger test_bridge_v1_*
//! prm family (also not archived); only this one DEM-validation test moved.
#pragma once

#include "../../DEMValidationScenarios.h"

namespace walberla
{
namespace antidunes
{
namespace dem_validation
{
//======================================================================================================================
// 10: arch tunnel containment -- fluid-free local validation for ComputeTunnelContacts (ParticleKernels.h),
// isolated from the fluid-coupling/dormant-floater changes so a containment bug surfaces here first
// rather than in a full bridge run. Two claims per shape ("sphere"/"box"/"cylinder", swap via the
// `shape` prm key):
//   1. CONTAINMENT: dropped from rest near the tunnel axis under gravity, the body's radial distance
//      from the axis (plus its own circumscribing radius) never exceeds the tunnel radius -- it settles
//      resting against the wall rather than clipping through it.
//   2. NO PHANTOM WALL: the SAME drop, but positioned outside the tunnel's axial [xMin,xMax] span,
//      shows ZERO lateral (Y) deflection -- confirms the X-range gate doesn't leak a force outside the
//      arch's own footprint (there's nothing else in this bare scenario to deflect it sideways, so any
//      Y motion at all means the gate leaked).
//======================================================================================================================
inline bool runArchTunnelContainment(const config::Config::BlockHandle& b)
{
   const std::string shape     = b.getParameter< std::string >("shape", "sphere");
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(real_t(0.8)));
   const real_t cylRadius      = b.getParameter< real_t >("cylRadius", real_t(0.8));
   const real_t cylLength      = b.getParameter< real_t >("cylLength", real_t(2));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t e              = b.getParameter< real_t >("restitutionCoefficient", real_t(0.5));
   const real_t mu             = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t gravityMag     = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t tunnelRadius   = b.getParameter< real_t >("tunnelRadius", real_t(5));
   const real_t tunnelHalfSpan = b.getParameter< real_t >("tunnelHalfSpan", real_t(5)); // xMin=-half, xMax=+half
   const uint_t settleSteps    = b.getParameter< uint_t >("settleSteps", uint_t(4000));
   const real_t tolerance      = b.getParameter< real_t >("containmentTolerance", real_t(0.05)); // LU slack

   detail::PassFailAccumulator acc;
   const real_t effR = detail::shapeCircumRadius(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);

   // --- Case 1: containment -- start just ABOVE where the body will come to rest against the tunnel
   // wall (a small gap, not at the axis), well inside [xMin,xMax], so most of the run is spent actually
   // resolving the wall contact and settling rather than in free-fall approaching it. -----------------
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, gravityMag, tunnelHalfSpan + tunnelRadius + real_t(20));
      engine.addArchTunnel(Vec3(0, 0, 0), tunnelRadius, -tunnelHalfSpan, tunnelHalfSpan, mu);
      const Vec3   startPos = Vec3(0, 0, -(tunnelRadius - effR) + real_t(0.3) * effR);
      const size_t idx      = detail::addRandomOrientedShape(engine, shape, startPos, sphereRadius,
                                                              boxHalfExtents, cylRadius, cylLength, density);
      engine.finalizeSetup(real_t(2) * effR);

      real_t maxRadialReach = real_t(0);
      const uint_t maxSteps = settleSteps;
      for (uint_t s = 0; s < maxSteps; ++s)
      {
         engine.step();
         const Vec3   pos    = engine.ac().getPosition(idx);
         const real_t radial = std::sqrt(pos[1] * pos[1] + pos[2] * pos[2]) + effR;
         maxRadialReach       = std::max(maxRadialReach, radial);
      }
      const bool pass = maxRadialReach <= tunnelRadius + tolerance;
      acc.record(pass);
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] archTunnelContainment shape=" << shape
                                 << " maxRadialReach=" << maxRadialReach << " tunnelRadius=" << tunnelRadius
                                 << " (tolerance=" << tolerance << ")  " << (pass ? "PASS" : "FAIL"));
   }

   // --- Case 2: no phantom wall -- same drop, positioned outside the tunnel's axial span. -------------
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, gravityMag, tunnelHalfSpan + tunnelRadius + real_t(20));
      engine.addArchTunnel(Vec3(0, 0, 0), tunnelRadius, -tunnelHalfSpan, tunnelHalfSpan, mu);
      const real_t startY    = real_t(0.3) * (tunnelRadius - effR); // same lateral offset as case 1
      const Vec3   startPos  = Vec3(tunnelHalfSpan + real_t(5) + effR, startY, 0);
      const size_t idx       = detail::addRandomOrientedShape(engine, shape, startPos, sphereRadius,
                                                               boxHalfExtents, cylRadius, cylLength, density);
      engine.finalizeSetup(real_t(2) * effR);

      for (uint_t s = 0; s < settleSteps; ++s)
         engine.step();

      const Vec3   pos      = engine.ac().getPosition(idx);
      const real_t yDrift   = std::abs(pos[1] - startY);
      const bool   pass     = yDrift < real_t(1e-6); // pure free-fall: Y must be untouched
      acc.record(pass);
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] archTunnelContainment (outside xRange) shape=" << shape
                                 << " yDrift=" << yDrift << " (expect ~0, no phantom wall)  "
                                 << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] archTunnelContainment OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}
} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
