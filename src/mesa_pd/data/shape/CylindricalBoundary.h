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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/shape/BaseShape.h>

#include <core/math/Constants.h>

namespace walberla {
namespace mesa_pd {
namespace data {

class CylindricalBoundary : public BaseShape
{
public:
   explicit CylindricalBoundary(const real_t& radius = real_t(1),
                                const Vec3& axis = Vec3(real_t(1), real_t(0), real_t(0)))
      : BaseShape(CylindricalBoundary::SHAPE_TYPE)
      , radius_(radius)
      , axis_(axis)
   {
      WALBERLA_CHECK_FLOAT_EQUAL(axis.sqrLength(), real_t(1.0), "Axis has to be normalized!");
   }

   const real_t& getRadius() const { return radius_; }
   const Vec3&   getAxis() const   { return axis_; }

   real_t getVolume() const override { return std::numeric_limits<real_t>::infinity(); }
   void   updateMassAndInertia(const real_t density) override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   // --- Optional containment-mode fields. Default (isContainment=false, no axial bounds, no half
   // clip) is byte-for-byte the original shape: a solid, axially-infinite pillar a body collides with
   // from OUTSIDE. Setting isContainment=true instead makes this a HOLLOW TUBE INTERIOR -- a body is
   // contained by / pushed back toward the axis instead of repelled from it. axialMin_/axialMax_
   // (optional) bound the tube to a finite span along its own axis, measured as an offset from this
   // shape's own body position (0 = at the body's own centre); a halfClipNormal_ (optional) further
   // restricts the solid material to the half of the tube on its positive side of a plane THROUGH the
   // axis (e.g. a bridge arch is only a semicircle of solid material above its own springing line,
   // open below -- clip normal = world Z, plane through the arch body's own position). All three are
   // additive/backward-compatible: any existing CylindricalBoundary user (pillars, domain-boundary
   // containment in AnalyticCollisionFunctions.h, ParticlePacking/Mixer) never sets them and sees no
   // behavior change.
   bool   isContainment() const { return isContainment_; }
   void   setContainment(bool v) { isContainment_ = v; }

   bool   hasAxialBounds() const
   {
      return axialMin_ > -std::numeric_limits<real_t>::infinity() ||
             axialMax_ <  std::numeric_limits<real_t>::infinity();
   }
   real_t getAxialMin() const { return axialMin_; }
   real_t getAxialMax() const { return axialMax_; }
   void   setAxialBounds(real_t lo, real_t hi) { axialMin_ = lo; axialMax_ = hi; }

   bool        hasHalfClip() const { return halfClipNormal_.sqrLength() > real_t(0); }
   const Vec3& getHalfClipNormal() const { return halfClipNormal_; }
   void        setHalfClipNormal(const Vec3& n) { halfClipNormal_ = n; }

   constexpr static int SHAPE_TYPE = 2; ///< Unique shape type identifier for cylindrical boundaries.\ingroup mesa_pd_shape

private:
      real_t radius_; ///< radius of the cylinder
      Vec3   axis_;   ///< axis of the cylinder

      bool   isContainment_   = false;
      real_t axialMin_        = -std::numeric_limits<real_t>::infinity();
      real_t axialMax_        =  std::numeric_limits<real_t>::infinity();
      Vec3   halfClipNormal_  = Vec3(real_t(0), real_t(0), real_t(0)); // zero vector = disabled
};

inline
void CylindricalBoundary::updateMassAndInertia(const real_t /*density*/)
{
   mass_         = std::numeric_limits<real_t>::infinity();
   invMass_      = real_t(0);

   inertiaBF_    = Mat3(std::numeric_limits<real_t>::infinity());
   invInertiaBF_ = Mat3(real_t(0));
}

inline
void CylindricalBoundary::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << radius_;
   buf << axis_;
   buf << isContainment_;
   buf << axialMin_;
   buf << axialMax_;
   buf << halfClipNormal_;
}
inline
void CylindricalBoundary::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> radius_;
   buf >> axis_;
   buf >> isContainment_;
   buf >> axialMin_;
   buf >> axialMax_;
   buf >> halfClipNormal_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
