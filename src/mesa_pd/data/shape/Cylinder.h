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
//! \file
//! \author Dennys Huber <dennys.huber@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

#include <core/math/Constants.h>

namespace walberla {
namespace mesa_pd {
namespace data {

/**
 * Finite solid cylinder. The cylinder axis is the BODY-FRAME Z-AXIS; the world-frame axis
 * direction follows from the particle's rotation. In contrast to CylindricalBoundary (an
 * infinite hollow container boundary), this is a finite convex solid, so it works with the
 * GJK/EPA-based GeneralContactDetection via support().
 */
class Cylinder : public BaseShape
{
public:
   explicit Cylinder(const real_t radius = real_t(1), const real_t height = real_t(1))
      : BaseShape(Cylinder::SHAPE_TYPE), radius_(radius), height_(height)
   {}

   const real_t& getRadius() const { return radius_; }
   const real_t& getHeight() const { return height_; }

   real_t getVolume() const override { return math::pi * radius_ * radius_ * height_; }
   void   updateMassAndInertia(const real_t density) override;

   Vec3 support( const Vec3& bfD ) const override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 6; ///< Unique shape type identifier for cylinders.\ingroup mesa_pd_shape

private:
   real_t radius_; ///< radius of the cylinder
   real_t height_; ///< full length of the cylinder along its (body-frame z) axis
};

inline
void Cylinder::updateMassAndInertia(const real_t density)
{
   const real_t m  = getVolume() * density;
   // Solid cylinder about its centroid: I_axis = m r²/2, I_perp = m (3r² + h²)/12.
   const Mat3   I  = Mat3::makeDiagonalMatrix(
         m * (real_t(3) * radius_ * radius_ + height_ * height_) / real_t(12),
         m * (real_t(3) * radius_ * radius_ + height_ * height_) / real_t(12),
         m * radius_ * radius_ * real_t(0.5) );

   mass_         = m;
   invMass_      = real_t(1.0) / m;

   inertiaBF_    = I;
   invInertiaBF_ = I.getInverse();
}

inline
Vec3 Cylinder::support( const Vec3& bfD ) const
{
   // bfD is normalized (see collision_detection::Support). The farthest point is on the rim of
   // the end cap facing bfD: radial part follows the lateral direction, axial part is ±h/2.
   const real_t lateralLen = std::sqrt(bfD[0] * bfD[0] + bfD[1] * bfD[1]);
   Vec3 supportPoint;
   if (lateralLen > real_t(1e-12))
   {
      const real_t scale = radius_ / lateralLen;
      supportPoint[0] = bfD[0] * scale;
      supportPoint[1] = bfD[1] * scale;
   }
   else
   {
      supportPoint[0] = real_t(0);
      supportPoint[1] = real_t(0);
   }
   supportPoint[2] = math::sign(bfD[2]) * height_ * real_t(0.5);
   return supportPoint;
}

inline
void Cylinder::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << radius_;
   buf << height_;
}
inline
void Cylinder::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> radius_;
   buf >> height_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
