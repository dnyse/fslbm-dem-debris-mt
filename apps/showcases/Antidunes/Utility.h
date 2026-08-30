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
//! \file   Utility.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Box.h"
#include "mesa_pd/data/shape/CylindricalBoundary.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"

#include <algorithm>
#include <core/mpi/Broadcast.h>
#include <core/mpi/MPITextFile.h>
#include <core/mpi/Reduce.h>
#include <functional>
#include <iterator>

namespace walberla
{
namespace antidunes
{

inline void renameFile(const std::string& oldName, const std::string& newName)
{
   int result = std::rename(oldName.c_str(), newName.c_str());
   if (result != 0)
      WALBERLA_LOG_WARNING_ON_ROOT("Could not rename file " << oldName << " to " << newName << " with error code "
                                                            << result);
}

inline void write2DVectorToFile(const std::vector< real_t >& vec, uint_t len1, uint_t len2, std::string filename)
{
   std::ofstream file;
   file.open(filename.c_str());

   if (!file.is_open()) {
      WALBERLA_ABORT("Could not open " << filename);
   }

   file.precision(5);

   file << "# " << len1 << " " << len2 << "\n";

   for (uint_t j = uint_t(0); j < len2; ++j)
   {
      for (uint_t i = uint_t(0); i < len1; ++i)
      {
         file << vec[i + j * len1] << "\n";
      }
   }
   file.close();
}

inline void writeSphereInformationToFile(const std::string& filename, walberla::mesa_pd::data::ParticleStorage& ps,
                                  Vector3< real_t >& domainSize, int precision = 12)
{
   std::ostringstream ossData;
   ossData << std::setprecision(precision);

   WALBERLA_ROOT_SECTION() { ossData << domainSize[0] << " " << domainSize[1] << " " << domainSize[2] << "\n"; }

   for (auto pIt : ps)
   {
      using namespace walberla::mesa_pd::data;
      if (pIt->getBaseShape()->getShapeType() != Sphere::SHAPE_TYPE) continue;
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(pIt->getFlags(), GHOST)) continue;
      auto sp = static_cast< Sphere* >(pIt->getBaseShape().get());

      auto position = pIt->getPosition();

      ossData << pIt->getUid() << " " << position[0] << " " << position[1] << " " << position[2] << " "
              << sp->getRadius() << '\n';
   }

   walberla::mpi::writeMPITextFile(filename, ossData.str());
}

inline void getAvgDiameterScalingFactor(const std::string& filename, const Vector3< uint_t >& domainSize,
                                 const uint_t bedCopiesInX, const uint_t bedCopiesInY, real_t& avgDiameter,
                                 real_t& scalingFactor)
{
   using namespace walberla;

   std::string textFile;

   WALBERLA_ROOT_SECTION()
   {
      std::ifstream t(filename.c_str());
      if (!t) { WALBERLA_ABORT("Invalid input file " << filename << "\n"); }
      std::stringstream buffer;
      buffer << t.rdbuf();
      textFile = buffer.str();

      std::istringstream fileIss(textFile);
      std::string line;

      // first line contains generation domain sizes
      std::getline(fileIss, line);
      Vector3< real_t > generationDomainSize_SI(0_r);
      std::istringstream firstLine(line);
      firstLine >> generationDomainSize_SI[0] >> generationDomainSize_SI[1] >> generationDomainSize_SI[2];

      real_t diameter_SI  = 0.0;
      uint_t numParticles = 0;
      while (std::getline(fileIss, line))
      {
         std::istringstream iss(line);

         mesa_pd::data::ParticleStorage::uid_type uID;
         mesa_pd::data::ParticleStorage::position_type pos(0_r);
         walberla::real_t radius = 0;
         iss >> uID >> pos[0] >> pos[1] >> pos[2] >> radius;
         WALBERLA_CHECK_GREATER(radius, 0_r, "Invalid radius of " << radius << " found in input file!")

         diameter_SI += 2_r * radius;

         numParticles++;
      }
      diameter_SI /= real_t(numParticles);

      scalingFactor = real_c(domainSize[0]) / (generationDomainSize_SI[0] * real_c(bedCopiesInX));
      avgDiameter   = diameter_SI * scalingFactor;

      WALBERLA_CHECK_EQUAL(uint_c(scalingFactor * generationDomainSize_SI[1] * real_c(bedCopiesInY)), domainSize[1],
                           "Error: Generated bed with copies and simulation domain do not match!")
   }

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(scalingFactor)
   walberla::mpi::broadcastObject(scalingFactor);
   walberla::mpi::broadcastObject(avgDiameter);
}

template< typename ParticleAccessor_T>
inline void getAverageVelocity(const ParticleAccessor_T& ac, real_t& averageVelocity,
                        real_t& maxVelocity, uint_t& numParticles, real_t& maxHeight)
{
   averageVelocity = real_t(0);
   maxVelocity     = real_t(0);
   numParticles    = uint_t(0);
   maxHeight       = real_t(0);
   for (uint_t i = 0; i < ac.size(); ++i)
   {
      if (isSet(ac.getFlags(i), walberla::mesa_pd::data::particle_flags::GHOST)) continue;
      if (isSet(ac.getFlags(i), walberla::mesa_pd::data::particle_flags::GLOBAL)) continue;

      ++numParticles;
      real_t velMagnitude = ac.getLinearVelocity(i).length();
      averageVelocity += velMagnitude;
      maxVelocity = std::max(maxVelocity, velMagnitude);
      maxHeight   = std::max(maxHeight, ac.getPosition(i)[2]);
   }

   walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
   walberla::mpi::allReduceInplace(averageVelocity, walberla::mpi::SUM);
   walberla::mpi::allReduceInplace(maxVelocity, walberla::mpi::MAX);
   walberla::mpi::allReduceInplace(maxHeight, walberla::mpi::MAX);

   averageVelocity /= real_t(numParticles);
}

// Box subclass that behaves like a wall: infinite mass so meff = sphere mass (identical to HalfSpace behavior)
class FixedBox : public mesa_pd::data::Box {
public:
   explicit FixedBox(const mesa_pd::Vec3& edgeLength) : Box(edgeLength) {}
   void updateMassAndInertia(const real_t /*density*/) override {
      mass_         = std::numeric_limits<real_t>::infinity();
      invMass_      = real_t(0);
      inertiaBF_    = mesa_pd::Mat3(std::numeric_limits<real_t>::infinity());
      invInertiaBF_ = mesa_pd::Mat3(real_t(0));
   }
};

// FloatingBox is a freely rotating mesa_pd Box — translation AND rotation are fully free.
// The base Box::updateMassAndInertia computes the correct uniform-box inertia tensor:
//   I_x = m/12 * (Ly² + Lz²),  I_y = m/12 * (Lx² + Lz²),  I_z = m/12 * (Lx² + Ly²).
// Contact detection (Box×HalfSpace, Box×CylindricalBoundary) accounts for the current rotation
// via the R^T body-frame transform. The submerged-volume estimate in computeSubmergedVolumesBox()
// integrates the fill field over the box transformed into its body frame, so it is EXACT for any
// tilt angle (not an axis-aligned / small-angle approximation). Orientation is written to VTK
// ("floaterBox" output, "rotation" 3x3 tensor) for rendering the oriented solid in ParaView.
class FloatingBox : public mesa_pd::data::Box {
public:
   explicit FloatingBox(const mesa_pd::Vec3& edgeLength) : Box(edgeLength) {}
   // Inherits Box::updateMassAndInertia — no rotation lock.
};

inline auto createBox(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& center, const mesa_pd::Vec3& edgeLength)
{
   auto p0 = ps.create(true);
   p0->setPosition(center);
   p0->setBaseShape(std::make_shared< FixedBox >(edgeLength));
   p0->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(0);
   p0->setInteractionRadius(std::numeric_limits< real_t >::infinity());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}

// Fixed-mass cylinder using CylindricalBoundary shape — repurposed for solid exterior contact via
// SolidCylinderContactDetection. axis should be normalised (e.g. Z-axis for a vertical pillar).
// End caps are handled by the existing floor/ceiling HalfSpace planes.
inline auto createCylinder(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& center,
                           real_t radius, const mesa_pd::Vec3& axis = mesa_pd::Vec3(0, 0, 1))
{
   auto p0 = ps.create(true);
   p0->setPosition(center);
   p0->setBaseShape(std::make_shared< mesa_pd::data::CylindricalBoundary >(radius, axis));
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

inline auto createPlane(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& pos, const mesa_pd::Vec3& normal)
{
   auto p0 = ps.create(true);
   p0->setPosition(pos);
   p0->setBaseShape(std::make_shared< mesa_pd::data::HalfSpace >(normal));
   p0->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(0);
   p0->setInteractionRadius(std::numeric_limits< real_t >::infinity());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}

} // namespace antidunes
} // namespace walberla