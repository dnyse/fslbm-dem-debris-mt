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
 * and a ball of the given radius. `height` is the FULL end-to-end extent INCLUDING the rounded caps
 * (matching Cylinder's own `height` = full flat-capped length, so the two are directly comparable /
 * swappable at the same radius+height numbers); the straight core segment's half-length is therefore
 * coreHalfLength = max(0, height/2 - radius). A capsule with height <= 2*radius degenerates to a
 * sphere of that radius (coreHalfLength = 0).
 *
 * Added specifically so GJK/EPA-based GeneralContactDetection (support()-only, no closed form) can be
 * exercised without writing any analytic capsule contact algorithm: see this shape's own support()
 * below, which is a 3-line reduction to Cylinder's own analytic support() plus a radius offset.
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

   // Axial (spin) moment: the core cylinder's own m r^2/2, plus the caps' contribution -- by symmetry a
   // hemisphere's moment about the axis THROUGH ITS OWN POLE is identical to a full sphere's, (2/5) m r^2,
   // regardless of where along that axis it sits, so both hemispheres combined contribute (2/5) mCaps r^2.
   const real_t Iaxis = mCyl * r * r * real_t(0.5) + mCaps * (real_t(2) / real_t(5)) * r * r;

   // Perpendicular moment about the capsule's own centre (parallel-axis theorem for the two caps):
   // core cylinder about its own centre (= capsule centre): m(3r^2+Lc^2)/12.
   // Each hemisphere's own-centroid perpendicular moment is (83/320) m_hemi r^2 -- derived here
   // rather than quoted, since an earlier version of this derivation had the two standard values
   // swapped (caught by ANTIDUNES_CapsuleShape's testCapsuleDegenerateSphereLimit, which checks the
   // exact height<=2r sphere limit): split a full sphere of mass M=2*m_hemi through its own centre;
   // the flat-face-diameter axis of either resulting hemisphere lies EXACTLY IN that splitting
   // plane, so by mirror symmetry (z -> -z swaps the two hemispheres and fixes that axis) each
   // hemisphere contributes EQUALLY to the full sphere's own (2/5) M r^2 about it -- i.e. (1/5) M r^2
   // = (2/5) m_hemi r^2 per hemisphere, about the FLAT-FACE diameter (through the flat face's own
   // centre, not the hemisphere's centroid). Shifting via parallel axis IN REVERSE to the
   // hemisphere's own centroid, a distance 3r/8 further along the symmetry axis, SUBTRACTS
   // m_hemi*(3r/8)^2: 2/5 - 9/64 = 128/320 - 45/320 = 83/320. (Sanity check: at a=0, degenerating to
   // a single sphere, hemiOffset=3r/8 and adding that same m_hemi*(3r/8)^2 back via the forward
   // parallel-axis step below exactly restores 2/5 m_hemi r^2 per hemisphere, i.e. (2/5) M r^2 total
   // -- the correct, exact solid-sphere moment.)
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
   // Minkowski sum of the core segment's own support (identical construction to Cylinder::support(),
   // with height_/2 replaced by coreHalfLength -- including the same near-perpendicular guard: without
   // it, floating-point noise in bfD[2] around an intended-exactly-zero axial component still resolves
   // sign(bfD[2]) to a definite +-1, picking an arbitrary end instead of the physically-symmetric
   // segment-centre point) and the cap sphere's own support (radius_ * bfD, bfD already normalised).
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
