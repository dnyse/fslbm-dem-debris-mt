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

class Box : public BaseShape
{
public:
   explicit Box(const Vec3& edgeLength = Vec3(real_t(1)))
      : BaseShape(Box::SHAPE_TYPE)
      , edgeLength_(edgeLength)
   {}

   const Vec3&   getEdgeLength() const   { return edgeLength_; }

   real_t getVolume() const override { return edgeLength_[0] * edgeLength_[1] * edgeLength_[2]; }
   void   updateMassAndInertia(const real_t density) override;

   Vec3 support( const Vec3& bfD ) const override;

   void pack(walberla::mpi::SendBuffer& buf) override;
   void unpack(walberla::mpi::RecvBuffer& buf) override;

   constexpr static int SHAPE_TYPE = 3; ///< Unique shape type identifier for boxes.\ingroup mesa_pd_shape

private:
   Vec3   edgeLength_;   ///< edge length of the box
};

inline
void Box::updateMassAndInertia(const real_t density)
{
   const real_t m = getVolume() * density;
   const Mat3   I = Mat3::makeDiagonalMatrix(
         edgeLength_[1]*edgeLength_[1] + edgeLength_[2]*edgeLength_[2] ,
         edgeLength_[0]*edgeLength_[0] + edgeLength_[2]*edgeLength_[2] ,
         edgeLength_[0]*edgeLength_[0] + edgeLength_[1]*edgeLength_[1] ) * (m / static_cast<real_t>( 12 ));

   mass_ = m;
   invMass_      = real_t(1.0) / m;

   inertiaBF_ = I;
   invInertiaBF_ = I.getInverse();
}

inline
Vec3 Box::support( const Vec3& bfD ) const
{
   // Per axis: zero offset unless bfD has a genuinely nonzero component along that axis. Without this
   // guard, sign(bfD[i]) still resolves to a definite +-1 even when bfD[i] is floating-point noise
   // around an intended-exactly-zero component (direction lying in, or very close to, the plane
   // perpendicular to that axis -- e.g. a box face flush against a wall) -- the true support region
   // there is the WHOLE corresponding face/edge (every point on it is an equally valid maximizer), and
   // jumping to a full +-edgeLength/2 offset on that axis picks one arbitrary corner among several
   // equally-valid ones instead of the physically-symmetric face/edge-centre point. Worse than merely
   // wrong: since the resulting spurious lever arm's torque then rotates the body, the NEXT support()
   // call can flip which near-zero component's sign resolves which way, flickering between different
   // corners step-to-step and injecting an oscillating torque -- found directly via DEMValidation's
   // wallImpactBox at tilt=0 degrees (flush face-on), where consecutive [STRONG-CONTACT-FORCE] log
   // lines within the SAME contact event showed the contact point jumping between near-opposite
   // corners of the same flush face, with torque flipping sign correspondingly. Same fix as
   // Cylinder::support() for the same class of degenerate-direction ambiguity.
   constexpr real_t epsilon = real_t(1e-12);
   return Vec3( (std::abs(bfD[0]) > epsilon) ? math::sign(bfD[0]) * edgeLength_[0] * real_t(0.5) : real_t(0),
                (std::abs(bfD[1]) > epsilon) ? math::sign(bfD[1]) * edgeLength_[1] * real_t(0.5) : real_t(0),
                (std::abs(bfD[2]) > epsilon) ? math::sign(bfD[2]) * edgeLength_[2] * real_t(0.5) : real_t(0) );
}

inline
void Box::pack(walberla::mpi::SendBuffer& buf)
{
   BaseShape::pack(buf);
   buf << edgeLength_;
}
inline
void Box::unpack(walberla::mpi::RecvBuffer& buf)
{
   BaseShape::unpack(buf);
   buf >> edgeLength_;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
