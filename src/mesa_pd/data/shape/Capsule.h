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
 * Finite solid capsule ("pill" shape): a cylindrical core with two hemispherical caps, i.e. the
 * Minkowski sum of a line segment (the body-frame Z-axis, from -coreHalfLength to +coreHalfLength)
 * and a ball of the given radius. `height` is the FULL end-to-end extent INCLUDING the rounded caps,
 * like Cylinder's `height`, so both shapes take the same (radius, height); hence
 * coreHalfLength = max(0, height/2 - radius). A capsule with height <= 2*radius degenerates to a
 * sphere of that radius (coreHalfLength = 0).
 *
 * mesa_pd provides no analytic capsule contact functions; GeneralContactDetection handles capsules
 * through support() alone (GJK/EPA, or a single support query against a HalfSpace or
 * CylindricalBoundary); support() is the segment's support point plus radius * direction.
 */
class Capsule : public BaseShape
{
public:
   explicit Capsule(const real_t radius = real_t(1), const real_t height = real_t(1))
      : BaseShape(Capsule::SHAPE_TYPE), radius_(radius), height_(height)
   {}

   const real_t& getRadius() const { return radius_; }
   const real_t& getHeight() const { return height_; }
   real_t getCoreHalfLength() const { return std::max(real_t(0), height_ * real_t(0.5) - radius_); }

   real_t getVolume() const override
   {
      const real_t coreLen = real_t(2) * getCoreHalfLength();
      // cylindrical core + a full sphere's worth of volume from the two hemispherical caps combined
      return math::pi * radius_ * radius_ * coreLen + (real_t(4) / real_t(3)) * math::pi * radius_ * radius_ * radius_;
   }
   void updateMassAndInertia(const real_t density) override;

   Vec3 support( const Vec3& bfD ) const override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 7; ///< Unique shape type identifier for capsules.\ingroup mesa_pd_shape

private:
   real_t radius_; ///< radius of the capsule (core cylinder radius = cap sphere radius)
   real_t height_; ///< full end-to-end length of the capsule INCLUDING the rounded caps
};

inline
void Capsule::updateMassAndInertia(const real_t density)
{
   const real_t r  = radius_;
   const real_t a  = getCoreHalfLength();      // core cylinder half-length
   const real_t Lc = real_t(2) * a;             // core cylinder full length

   const real_t mCyl  = math::pi * r * r * Lc * density;                              // cylindrical core mass
   const real_t mCaps = (real_t(4) / real_t(3)) * math::pi * r * r * r * density;      // BOTH hemispheres combined (= one full sphere)
   const real_t mHemi = real_t(0.5) * mCaps;                                          // ONE hemisphere's mass

   // Axial (spin) moment: core cylinder m r^2/2 plus the caps. A hemisphere's moment about its symmetry
   // axis equals a full sphere's, (2/5) m r^2, independent of its axial position, so both caps add
   // (2/5) mCaps r^2.
   const real_t Iaxis = mCyl * r * r * real_t(0.5) + mCaps * (real_t(2) / real_t(5)) * r * r;

   // Perpendicular moment about the capsule centre (parallel-axis theorem for the two caps):
   // core cylinder about its own centre (= capsule centre): m(3r^2+Lc^2)/12.
   // Hemisphere about a diameter of its flat face: by mirror symmetry half of a full sphere's moment,
   // (2/5) m_hemi r^2. Its centroid lies 3r/8 from the flat face, so about the centroid:
   // (2/5 - 9/64) m_hemi r^2 = (83/320) m_hemi r^2. The centroid sits at a + 3r/8 from the capsule
   // centre (hemiOffset). For a = 0 this gives the solid-sphere value (2/5) m r^2.
   const real_t IcylPerp   = mCyl * (real_t(3) * r * r + Lc * Lc) / real_t(12);
   const real_t hemiOffset = a + (real_t(3) / real_t(8)) * r;
   const real_t IhemiPerp  = (real_t(83) / real_t(320)) * mHemi * r * r + mHemi * hemiOffset * hemiOffset;
   const real_t Iperp      = IcylPerp + real_t(2) * IhemiPerp;

   const real_t m = mCyl + mCaps;
   const Mat3   I = Mat3::makeDiagonalMatrix(Iperp, Iperp, Iaxis);

   mass_         = m;
   invMass_      = real_t(1.0) / m;

   inertiaBF_    = I;
   invInertiaBF_ = I.getInverse();
}

inline
Vec3 Capsule::support( const Vec3& bfD ) const
{
   // Support of the core segment plus the ball's support radius_ * bfD (bfD is normalised). The segment
   // end is chosen as in Cylinder::support()'s axial component, with the same near-perpendicular guard
   // (segment centre when |bfD[2]| <= 1e-12).
   const real_t coreHalfLen = getCoreHalfLength();
   Vec3 segPoint(real_t(0), real_t(0),
                 (std::abs(bfD[2]) > real_t(1e-12)) ? math::sign(bfD[2]) * coreHalfLen : real_t(0));
   return segPoint + bfD * radius_;
}

inline
void Capsule::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << radius_;
   buf << height_;
}
inline
void Capsule::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> radius_;
   buf >> height_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
