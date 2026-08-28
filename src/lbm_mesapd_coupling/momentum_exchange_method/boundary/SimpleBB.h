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
//! \file SimpleBB.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/logging/all.h"

#include "field/FlagField.h"

#include "lbm/field/PdfField.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/ParticleFunctions.h"
#include "mesa_pd/data/IAccessor.h"

#include "stencil/Directions.h"

#include <algorithm>
#include <limits>

namespace walberla {
namespace lbm_mesapd_coupling {

//**************************************************************************************************************************************
/*!
*   \brief Bounce back boundary handling for moving particles
*
*   This boundary condition implements the bounce back scheme to model a no-slip boundary condition for moving particles
*
*/
//**************************************************************************************************************************************

template< typename LatticeModel_T, typename FlagField_T, typename ParticleAccessor_T >
class SimpleBB : public Boundary< typename FlagField_T::flag_t >
{
   using PDFField_T = lbm::PdfField< LatticeModel_T >;
   using Stencil_T = typename LatticeModel_T::Stencil;
   using flag_t = typename FlagField_T::flag_t;

   static_assert(std::is_base_of_v<mesa_pd::data::IAccessor, ParticleAccessor_T>, "Provide a valid accessor as template");

public:

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& )
   {
      WALBERLA_ABORT( "A SimpleBB boundary cannot be created from a config file" );
      return make_shared<BoundaryConfiguration>();
   }

   inline SimpleBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                    ParticleField_T * const particleField, const shared_ptr<ParticleAccessor_T>& ac,
                    const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block,
                    std::function<real_t(const Vector3<real_t>&)> hydrostaticDensityFct = nullptr,
                    std::function<real_t(const Cell&)> fillLevelFct = nullptr );

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T >
   void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) {}

   void registerCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t, const BoundaryConfiguration& ) {}
   void registerCells( const flag_t, const CellInterval&, const BoundaryConfiguration& ) {}
   template< typename CellIterator >
   void registerCells( const flag_t, const CellIterator&, const CellIterator&, const BoundaryConfiguration& ) {}

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );

private:

   const FlagUID uid_;


   PDFField_T  * const  pdfField_;
   const FlagField_T * const flagField_;
   ParticleField_T * const particleField_;
   shared_ptr<ParticleAccessor_T> ac_;

   flag_t domainMask_;

   const StructuredBlockStorage & blockStorage_;
   const IBlock & block_;

   std::function<real_t(const Vector3<real_t>&)> hydrostaticDensityFct_;

   // Optional fill-level (φ) accessor for the three-phase MEM fill-weighting (Bogner & Rüde 2013, Eq. 5.9).
   // Returns the free-surface fill level of a given cell (1 = full liquid, 0 = gas). When set, the
   // hydrodynamic force from each obstacle link is scaled by φ of the near-boundary cell, so interface
   // cells at a floating particle's waterline do not inject a phantom (full-liquid) buoyant force. Left
   // null for non-free-surface couplings, in which case the standard (Eq. 5.6) force is used unchanged.
   std::function<real_t(const Cell&)> fillLevelFct_;

   real_t lengthScalingFactor_;
   real_t forceScalingFactor_;

}; // class SimpleBB


template< typename LatticeModel_T, typename FlagField_T, typename ParticleAccessor_T >
inline SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >::SimpleBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                                                                              ParticleField_T * const particleField, const shared_ptr<ParticleAccessor_T>& ac,
                                                                              const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block,
                                                                              std::function<real_t(const Vector3<real_t>&)> hydrostaticDensityFct,
                                                                              std::function<real_t(const Cell&)> fillLevelFct):
Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ), flagField_( flagField ), particleField_( particleField ), ac_( ac ), domainMask_(domain),
blockStorage_( blockStorage ), block_( block ), hydrostaticDensityFct_(hydrostaticDensityFct), fillLevelFct_(fillLevelFct)
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField_ );
   WALBERLA_ASSERT_NOT_NULLPTR( particleField_ );
   WALBERLA_ASSERT( flagField_->isRegistered( domainMask_ )  );

   // force scaling factor to account for different dx on this block
   const real_t dxCurrentLevel = blockStorage_.dx( blockStorage_.getLevel(block) );
   lengthScalingFactor_ = dxCurrentLevel;
   forceScalingFactor_ = lengthScalingFactor_ * lengthScalingFactor_;
}

// requires: getVelocityAtPosition (was velFromWF), addForceAtPosition (was addForceAtPos)
template< typename LatticeModel_T, typename FlagField_T, typename ParticleAccessor_T >
#ifndef NDEBUG
inline void SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                         const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
{
   WALBERLA_ASSERT_EQUAL  ( nx, x + cell_idx_t( stencil::cx[ dir ] ) );
   WALBERLA_ASSERT_EQUAL  ( ny, y + cell_idx_t( stencil::cy[ dir ] ) );
   WALBERLA_ASSERT_EQUAL  ( nz, z + cell_idx_t( stencil::cz[ dir ] ) );
   WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
   WALBERLA_ASSERT_EQUAL  ( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                               // current implementation of this boundary condition
   WALBERLA_ASSERT_UNEQUAL( particleField_->get(nx,ny,nz), ac_->getInvalidUid(), "UID of particle is invalid!" );

   real_t pdf_old = pdfField_->get( x, y, z, Stencil_T::idx[dir] );
   const real_t pdf_old_raw = pdf_old; // pre-hydrostatic-subtraction value, kept only for [NAN-LINK] logging below

   // Dry-link guard for free-surface coupling: a cell on the fluid side of this obstacle link may be
   // a GAS cell (e.g. at the waterline of a floating particle), which carries no valid PDF — its value
   // can be NaN/Inf. Performing momentum exchange across such a link would (a) write the garbage back
   // into the obstacle cell and (b) inject it into the hydrodynamic force on the particle (forceMEM
   // below), producing a NaN particle force / invalid quaternion. Physically, an air cell exerts no
   // momentum on the particle, so we skip the link entirely. NOTE: in a fully submerged setup this can
   // never trigger (all obstacle-adjacent cells are liquid/interface); it is specific to particles at
   // the free surface. Genuine corruption of a *liquid/interface* PDF is caught separately by the
   // checkPdfNaN scanner, which is not silenced by this guard.
   if( !std::isfinite(pdf_old) ) {
      // TEMP DIAGNOSTIC (remove after verifying): log the flag of the dry-link fluid cell. The domain
      // mask is liquid|interface, so this cell can never be gas — the expectation is that it is a
      // drained INTERFACE cell (ρ→0). Throttled to the first 30 hits per rank to avoid flooding.
      static uint_t dryLinkLogCount = uint_t(0);
      if( dryLinkLogCount < uint_t(30) ) {
         const auto fl = flagField_->get(x, y, z);
         std::string flagStr;
         for( const char* name : { "liquid", "interface", "gas" } ) {
            const field::FlagUID uid(name);
            if( flagField_->flagExists(uid) && field::isFlagSet(fl, flagField_->getFlag(uid)) )
               flagStr += std::string(name) + " ";
         }
         WALBERLA_LOG_WARNING("[DRY-LINK] non-finite pdf_old at cell (" << x << "," << y << "," << z
            << ") flags=[ " << flagStr << "] rawFlag=" << uint_t(fl) << " pdf_old=" << pdf_old);
         ++dryLinkLogCount;
      }
      // Expected for a free-surface particle; skip the link (air / drained cell exerts no momentum).
      return;
   }

   // apply bounce back scheme
   const real_t delta = real_c( 0.5 );
   real_t pdf_new = pdf_old;
   const real_t alpha = real_c(2);

   // get coordinates of fluid cell center
   Cell nearBoundaryCell(x,y,z);
   Vector3< real_t > cellCenter = blockStorage_.getBlockLocalCellCenter(block_, nearBoundaryCell );

   // get vector from fluid cell center to boundary cell center
   Vector3< real_t > direction( lengthScalingFactor_ * real_c( stencil::cx[ dir ] ),
                                lengthScalingFactor_ * real_c( stencil::cy[ dir ] ),
                                lengthScalingFactor_ * real_c( stencil::cz[ dir ] ) );

   //get particle index
   auto particleIdx = ac_->uidToIdx(particleField_->get( nx, ny, nz ));
   WALBERLA_ASSERT_UNEQUAL( particleIdx, ac_->getInvalidIdx(), "Index of particle is invalid!" );

   // Stale-mapping guard, same reasoning as the dry-link guard above: particleField_'s UID -> idx
   // lookup can go stale (the owning particle's storage slot was reused, migrated, or removed
   // without this cell's obstacle-ownership entry being refreshed yet) between when a cell was last
   // mapped as this particle's obstacle and when this sweep actually runs. The WALBERLA_ASSERT above
   // only catches this in a Debug build -- in Release/RelWithDebInfo (i.e. every production run, and
   // the previous ASan run that WAS built RelWithDebInfo) it's compiled out, so an invalid index
   // silently proceeded into getVelocityAtWFPoint() below, reading position/velocity through a
   // garbage index -- an AddressSanitizer-confirmed heap-buffer-overflow (8 bytes before a
   // ParticleStorage-backed vector's allocation), consistent with this project's own history of
   // heap-corruption / invalid-quaternion crashes whose symptom varied run to run depending on what
   // happened to occupy the over-read memory. Treated the same as the dry-link case: this specific
   // obstacle link contributes no momentum exchange this sub-step rather than reading garbage.
   if( particleIdx == ac_->getInvalidIdx() ) {
      static uint_t staleParticleIdxLogCount = uint_t(0);
      if( staleParticleIdxLogCount < uint_t(30) ) {
         WALBERLA_LOG_WARNING("[STALE-PARTICLE-IDX] cell (" << nx << "," << ny << "," << nz
            << ") is flagged as owned by particle uid=" << particleField_->get(nx,ny,nz)
            << ", but that uid no longer maps to a valid local particle index -- skipping this "
               "obstacle link's momentum exchange instead of reading through an invalid index.");
         ++staleParticleIdxLogCount;
      }
      return;
   }

   // assumed boundary position
   const Vector3<real_t> boundaryPosition( cellCenter + delta * direction );

   auto boundaryVelocity = mesa_pd::getVelocityAtWFPoint(particleIdx, *ac_, boundaryPosition );

   // include effect of boundary velocity
   // Hoisted out of the compressible branch (was block-scoped to the `if` below) so [NAN-LINK] logging
   // further down can report it -- this is the near-boundary cell's FULL 27-PDF density (same quantity
   // ReconstructInterfaceCellABB.h's [LOW-RHO] guard checks), not just the single pdf_old already
   // finite-checked above. If OTHER PDF components at this cell have gone non-finite while pdf_old
   // (this one link's direction) is still finite, this density read silently carries that corruption
   // into pdf_new/forceMEM/force below without ever tripping the dry-link guard above.
   real_t density = std::numeric_limits< real_t >::quiet_NaN(); // stays NaN (marker for "n/a") on the incompressible path
   if( LatticeModel_T::compressible )
   {
       density  = pdfField_->getDensity(x,y,z);
       pdf_new -= real_c(3.0) * alpha * density * LatticeModel_T::w[ Stencil_T::idx[dir] ] *
                     ( real_c( stencil::cx[ dir ] ) * boundaryVelocity[0] +
                       real_c( stencil::cy[ dir ] ) * boundaryVelocity[1] +
                       real_c( stencil::cz[ dir ] ) * boundaryVelocity[2] );
   }
   else
   {
       pdf_new -= real_c(3.0) * alpha * LatticeModel_T::w[ Stencil_T::idx[dir] ] *
                     ( real_c( stencil::cx[ dir ] ) * boundaryVelocity[0] +
                       real_c( stencil::cy[ dir ] ) * boundaryVelocity[1] +
                       real_c( stencil::cz[ dir ] ) * boundaryVelocity[2] );
   }

   // carry out the boundary handling
   pdfField_->get( nx, ny, nz, Stencil_T::invDirIdx(dir) ) = pdf_new;

   // check if fluid cell is not inside ghost layer ( = is in inner part), since then no force is allowed to be added
   if( !pdfField_->isInInnerPart( nearBoundaryCell ) )
   {
      return;
   }

   // calculate the force on the obstacle
   // original work (MEM):      Ladd - Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1. Theoretical foundation (1994)
   // improved version (GIMEM): Wen et al. - Galilean invariant fluid-solid interfacial dynamics in lattice Boltzmann simulations (2014)

   if(LatticeModel_T::compressible)
   {
      // this artificially introduces the zero-centering of the PDF values here (see also /lbm/field/PdfField.h)
      // it is important to have correct force values when two or more particles are so close such that there is no fluid cell between
      // as a consequence, some (non-zero) PDF contributions would be missing after summing up the force contributions
      // those would need to be added artificially, see e.g. Ernst, Dietzel, Sommerfeld - A lattice Boltzmann method for simulating transport and agglomeration of resolved particles, Acta Mech, 2013
      // instead, we use the trick there that we just require the deviations from the equilibrium to get the correct force as it is already used for the incompressible case
       if (hydrostaticDensityFct_ == nullptr)
       {
          pdf_old -= LatticeModel_T::w[Stencil_T::idx[dir]];
          pdf_new -= LatticeModel_T::w[Stencil_T::idx[dir]];
       }
       else
       {
          const real_t rhoHydStat = hydrostaticDensityFct_(cellCenter);
          pdf_old -= rhoHydStat * LatticeModel_T::w[Stencil_T::idx[dir]];
          pdf_new -= rhoHydStat * LatticeModel_T::w[Stencil_T::idx[dir]];
       }
   }

   // MEM: F = pdf_old + pdf_new - common
   const real_t forceMEM = pdf_old + pdf_new;

   // correction from Wen
   const real_t correction = pdf_old - pdf_new;


   // force consists of the MEM part and the galilean invariance correction including the boundary velocity
   Vector3<real_t> force( real_c( stencil::cx[dir] ) * forceMEM - correction * boundaryVelocity[0],
                          real_c( stencil::cy[dir] ) * forceMEM - correction * boundaryVelocity[1],
                          real_c( stencil::cz[dir] ) * forceMEM - correction * boundaryVelocity[2] );

   // TEMP DIAGNOSTIC (remove after investigation): dumps every intermediate quantity that feeds into
   // this link's contribution to the particle's hydrodynamic force, the first time that contribution
   // comes out non-finite despite pdf_old (this link's own PDF) having passed the dry-link finite-check
   // above. This is the gap that check leaves open: pdf_old is only ONE of 27 populations at cell
   // (x,y,z); `density` (used above, compressible path) is getDensity(x,y,z) -- the SUM of all 27 --
   // so another direction's PDF can be corrupted (a [LOW-RHO] cell) while this link's own pdf_old is
   // still perfectly finite, silently carrying the corruption through pdf_new/forceMEM/correction into
   // `force` here without ever tripping [DRY-LINK]. Throttled to the first 50 hits per rank.
   if( !(std::isfinite(force[0]) && std::isfinite(force[1]) && std::isfinite(force[2])) )
   {
      static uint_t nanLinkLogCount = uint_t(0);
      if( nanLinkLogCount < uint_t(50) )
      {
         WALBERLA_LOG_WARNING("[NAN-LINK] non-finite force on particle uid="
            << ac_->getUid(particleIdx) << " from link cell (" << x << "," << y << "," << z
            << ") -> obstacle cell (" << nx << "," << ny << "," << nz << ") dir=" << uint_t(dir)
            << "  pdf_old_raw=" << pdf_old_raw << " pdf_old(post-hydrostatic)=" << pdf_old
            << " density=" << density << " boundaryVelocity=" << boundaryVelocity
            << " pdf_new=" << pdf_new << " forceMEM=" << forceMEM << " correction=" << correction
            << " force(pre-fill-weight,pre-scale)=" << force
            << " forceScalingFactor=" << forceScalingFactor_);
         ++nanLinkLogCount;
      }
   }

   // Three-phase MEM fill-weighting (Bogner & Rüde 2013, Eq. 5.9). The bounce-back force above (Eq. 5.6)
   // treats the near-boundary cell as full liquid (ρ≈1) regardless of its actual fill — but at a floating
   // particle's waterline that cell is an INTERFACE cell that is mostly gas. Counting it at full weight
   // injects a phantom buoyant force that holds floaters too high. Scaling the link's momentum exchange by
   // the boundary cell's fill level φ(x_b) makes an almost-empty interface cell contribute proportionally
   // less, recovering Eq. 5.9's liquid term c_q·φ·[f̃+f]. φ is clamped to [0,1] (free-surface dynamics can
   // briefly push it slightly outside). IMPORTANT: only the force returned to the particle is weighted — the
   // bounce-back PDF written back into the fluid (above) is untouched, so the free-surface LBM is unaffected.
   if( fillLevelFct_ )
   {
      const real_t fillRaw    = fillLevelFct_( nearBoundaryCell );
      const real_t fillWeight = std::clamp( fillRaw, real_c(0), real_c(1) );

      // TEMP DIAGNOSTIC (remove after verifying the fill-weighting): this branch is ONLY reached when a
      // fill-level accessor was wired in (i.e. the Eq. 5.9 fix is compiled in), so any line below is itself
      // proof the fix is live. We log the cells whose force is actually changed by the weighting (fillWeight
      // noticeably below 1), reporting the cell type, the raw/clamped fill, and the z-force before vs after
      // weighting (pre-dx-scaling). This shows which cell types + φ values carry the phantom buoyant force.
      // Throttled to the first 60 hits per rank to avoid flooding.
      if( fillWeight < real_c(0.99) )
      {
         static uint_t fillWeightLogCount = uint_t(0);
         if( fillWeightLogCount < uint_t(60) )
         {
            const auto fl = flagField_->get(x, y, z);
            std::string flagStr;
            for( const char* name : { "liquid", "interface", "gas" } ) {
               const field::FlagUID flagUID(name);
               if( flagField_->flagExists(flagUID) && field::isFlagSet(fl, flagField_->getFlag(flagUID)) )
                  flagStr += std::string(name) + " ";
            }
            WALBERLA_LOG_INFO("[FILL-WEIGHT] cell (" << x << "," << y << "," << z << ") flags=[ " << flagStr
               << "] phiRaw=" << fillRaw << " weight=" << fillWeight << " dir=" << uint_t(dir)
               << " forceZ: " << force[2] << " -> " << (force[2] * fillWeight)
               << " (forceMEM=" << forceMEM << ")");
            ++fillWeightLogCount;
         }
      }

      force *= fillWeight;
   }

   force *= forceScalingFactor_;

   // add the force onto the particle at the obstacle boundary
   lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic( particleIdx, *ac_, force, boundaryPosition );

}

} // namespace lbm_mesapd_coupling
} // namespace walberla
