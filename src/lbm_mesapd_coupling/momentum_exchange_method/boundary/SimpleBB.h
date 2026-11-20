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

#include <limits>
#include <sstream>
#include <string>

#include "field/FlagField.h"

#include "lbm/field/PdfField.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/ParticleFunctions.h"
#include "mesa_pd/common/RayParticleIntersection.h"
#include "mesa_pd/kernel/SingleCast.h"
#include "mesa_pd/data/IAccessor.h"

#include "core/mpi/MPIManager.h"
#include "stencil/Directions.h"

#include <algorithm>
#include <array>
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

// Momentum budget diagnostic (Antidunes ModelParameters.logMomentumBudget): rank-local sums over the
// inner links treated since the last reset, x and z components each, pre dx scaling. With phi the clamped fill of the
// near-boundary cell (1 without fill weighting):
//   [0,1] c (f_old + f_new)          momentum the fluid cell loses through the link
//   [2,3] phi c (f_old + f_new)      the same, fill weighted like the fluid mass
//   [4,5] phi c 2 rho_ref w_q        reference subtraction (rho_ref = 1, or rho_hs with a hydrostatic reference)
//   [6,7] -phi (f_old - f_new) u_b   Wen correction
//   [8,9] force given to the body    = [2,3] - [4,5] + [6,7]
// Only written while enabled() is true; not thread safe (Antidunes builds without OpenMP).
struct MemLinkMomentumBudget
{
   static std::array< real_t, 10 >& sums() { static std::array< real_t, 10 > s{}; return s; }
   static bool& enabled() { static bool e = false; return e; }
};

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
                    std::function<real_t(const Cell&)> fillLevelFct = nullptr,
                    bool useInterpolatedBounceBack = false );

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

   // Optional fill-level (φ) accessor for the three-phase MEM fill weighting (Bogner 2017, Eq. 5.9): returns a
   // cell's free-surface fill level (1 = liquid, 0 = gas). When set, each link's force is scaled by φ of the
   // near-boundary cell, so waterline interface cells do not act as full liquid. Null (non-free-surface couplings, or
   // Antidunes with useFillWeighting = false): the unweighted momentum exchange (Bogner 2017, Eq. 5.6).
   std::function<real_t(const Cell&)> fillLevelFct_;

   real_t lengthScalingFactor_;
   real_t forceScalingFactor_;

   // Optional interpolated bounce-back (default false: simple bounce-back). true: a link whose further fluid cell
   // x_ff = x - c_i is in the domain mask (Antidunes: liquid/interface) uses the CLI scheme of Ginzburg et al. (2008),
   // Eq. 4.1 / Tables 3-4 (as lbm_mesapd_coupling::CurvedLinear), with the wall distance delta from a ray-particle
   // intersection (analytic for spheres, half-spaces and ellipsoids, bisection for other shapes such as boxes and
   // cylinders). Other links (x_ff outside the field or not in the domain mask, non-finite PDFs, delta not in [0,1])
   // fall back to simple bounce-back; delta = 1/2 reproduces it exactly (kappa0 = 0, alpha = 2). The force below uses
   // the resulting pdf_new.
   bool useInterpolatedBounceBack_;

}; // class SimpleBB


template< typename LatticeModel_T, typename FlagField_T, typename ParticleAccessor_T >
inline SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >::SimpleBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField_T * const pdfField, const FlagField_T * const flagField,
                                                                              ParticleField_T * const particleField, const shared_ptr<ParticleAccessor_T>& ac,
                                                                              const flag_t domain, const StructuredBlockStorage & blockStorage, const IBlock & block,
                                                                              std::function<real_t(const Vector3<real_t>&)> hydrostaticDensityFct,
                                                                              std::function<real_t(const Cell&)> fillLevelFct,
                                                                              bool useInterpolatedBounceBack):
Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ), flagField_( flagField ), particleField_( particleField ), ac_( ac ), domainMask_(domain),
blockStorage_( blockStorage ), block_( block ), hydrostaticDensityFct_(hydrostaticDensityFct), fillLevelFct_(fillLevelFct), useInterpolatedBounceBack_(useInterpolatedBounceBack)
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
   const real_t pdf_old_raw = pdf_old; // value before the reference subtraction below (garbage-link guard, diagnostics)

   // Dry-link guard for free-surface coupling: a non-finite PDF in the near-boundary cell (liquid or interface, as the
   // domain mask excludes gas; e.g. a drained interface cell at a floater's waterline) would be written into the
   // obstacle cell and into the particle force (NaN force, invalid quaternion). Such a cell exerts no momentum, so the
   // link is skipped.
   // Corrupted liquid/interface PDFs elsewhere are still caught by Antidunes' checkPdfNaN scan.
   if( !std::isfinite(pdf_old) ) {
      // Logs the near-boundary cell's flags, first 30 hits per rank.
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
      return;
   }

   // Garbage-link guard: a finite but absurd population (|pdf_old| > 10; ordinary PDFs are about w_q rho, i.e.
   // 0.005-0.3) is skipped like a dry link. forceMEM uses pdf_old directly, so a single exploded population would
   // launch the particle, and the cell density (a sum) can look sane while individual populations have diverged, e.g.
   // in a liquid cell pinched between two bodies.
   if( std::abs(pdf_old_raw) > real_c(10) )
   {
      static uint_t garbageLinkLogCount = uint_t(0);
      if( garbageLinkLogCount < uint_t(40) )
      {
         ++garbageLinkLogCount;
         auto flagNames = [&](const cell_idx_t cx, const cell_idx_t cy, const cell_idx_t cz) {
            std::string out;
            const auto fl = flagField_->get(cx, cy, cz);
            for( const char* name : { "liquid", "interface", "gas", "former moving obstacle", "convert to liquid", "convert to gas", "convert to interface for gas neighbors",
                                      "convert to interface for liquid neighbors", "keep interface", "moving obstacle", "obstacle", "noSlip", "freeSlip", "ubb", "pressure" } ) {
               const field::FlagUID fuid(name);
               if( flagField_->flagExists(fuid) && field::isFlagSet(fl, flagField_->getFlag(fuid)) ) { if( !out.empty() ) out += "|"; out += name; }
            }
            return out.empty() ? std::string("rawFlag=") + std::to_string(uint_t(fl)) : out;
         };
         std::ostringstream pdfs;
         for( uint_t q = 0; q < Stencil_T::Size; ++q ) pdfs << (q ? "," : "") << pdfField_->get(x, y, z, q);
         std::ostringstream nb;
         for( cell_idx_t dz = -1; dz <= 1; ++dz ) for( cell_idx_t dy = -1; dy <= 1; ++dy ) for( cell_idx_t dx = -1; dx <= 1; ++dx )
            if( dx || dy || dz ) { const cell_idx_t cx = x + dx, cy = y + dy, cz = z + dz;
               if( pdfField_->xyzSizeWithGhostLayer().contains(cx, cy, cz) ) nb << " (" << dx << "," << dy << "," << dz << "):" << flagNames(cx, cy, cz); }
         const real_t fillHere = fillLevelFct_ ? fillLevelFct_(Cell(x, y, z)) : std::numeric_limits< real_t >::quiet_NaN();
         WALBERLA_LOG_WARNING("[MEM-GARBAGE-LINK] particle uid=" << particleField_->get(nx, ny, nz) << " fluid cell block-local (" << x << "," << y << "," << z << ") flags=[" << flagNames(x, y, z)
            << "] fill=" << fillHere << " density=" << pdfField_->getDensity(x, y, z) << " | obstacle cell (" << nx << "," << ny << "," << nz << ") flags=[" << flagNames(nx, ny, nz)
            << "] dir=" << uint_t(dir) << " pdf_old(dir)=" << pdf_old_raw << " | all populations of the fluid cell: " << pdfs.str() << " | neighbour flags:" << nb.str()
            << " -- SKIPPING this link (no momentum exchange contributed) rather than feeding this value into forceMEM.");
      }
      return;
   }

   // apply bounce back scheme (simple bounce-back by default; optionally overwritten by the interpolated variant below)
   real_t delta = real_c( 0.5 );
   real_t pdf_new = pdf_old;
   real_t alpha = real_c(2);

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

   // Stale-mapping guard: the owner UID stored in particleField_ may no longer map to a local particle (the
   // obstacle mapping has desynced from mesa_pd's storage). The assert above is compiled out in release builds, where
   // the invalid index would read out of bounds in getVelocityAtWFPoint() below; skipping the link instead would leave
   // reduceAndAvgHydroForces() with a partial force sum for the particle (which can cause a force spike and then an
   // invalid quaternion). Hence a hard abort here, reporting the near-boundary cell's fill, density and populations.
   if( particleIdx == ac_->getInvalidIdx() ) {
      std::ostringstream pdfsAtAbort;
      for( uint_t q = 0; q < Stencil_T::Size; ++q ) pdfsAtAbort << (q ? "," : "") << pdfField_->get(x, y, z, q);
      const real_t fillAtAbort = fillLevelFct_ ? fillLevelFct_(Cell(x, y, z)) : std::numeric_limits< real_t >::quiet_NaN();
      WALBERLA_ABORT("[STALE-PARTICLE-IDX] rank=" << mpi::MPIManager::instance()->rank()
         << " cell (" << nx << "," << ny << "," << nz << ") [near-boundary fluid cell ("
         << x << "," << y << "," << z << "), dir=" << dir
         << "] is flagged as owned by particle uid=" << particleField_->get(nx,ny,nz)
         << ", but that uid no longer maps to a valid local particle index. The fluid-side "
            "obstacle-cell mapping (particleField_) has desynced from mesa_pd's local particle "
            "storage for this UID -- its storage slot was reused, migrated, or removed without "
            "this cell's obstacle-ownership entry being refreshed. Silently skipping this link "
            "would leave reduceAndAvgHydroForces() with a partial/inconsistent momentum-exchange "
            "sum for this particle, which has been observed to produce a multi-order-of-magnitude "
            "hydrodynamic force spike and a subsequent invalid-quaternion crash several steps "
            "later -- aborting here instead, at the actual point of desync."
         << " | near-boundary fluid cell state AT THIS EXACT MOMENT: fill=" << fillAtAbort
         << " density=" << pdfField_->getDensity(x, y, z)
         << " | all populations of the fluid cell: " << pdfsAtAbort.str());
   }

   // Optional interpolated bounce-back, see the member comment of useInterpolatedBounceBack_.
   if( useInterpolatedBounceBack_ )
   {
      static uint_t ibbAppliedCount = uint_t(0), ibbFallbackCount = uint_t(0);
      bool applied = false;
      const cell_idx_t xff = x - cell_idx_c( stencil::cx[ dir ] );
      const cell_idx_t yff = y - cell_idx_c( stencil::cy[ dir ] );
      const cell_idx_t zff = z - cell_idx_c( stencil::cz[ dir ] );
      if( pdfField_->xyzSizeWithGhostLayer().contains( xff, yff, zff ) && flagField_->isPartOfMaskSet( xff, yff, zff, domainMask_ ) )
      {
         const real_t pdf_far = pdfField_->get( xff, yff, zff, Stencil_T::idx[dir] );
         const real_t pdf_inv = pdfField_->get( x, y, z, Stencil_T::invDirIdx(dir) );
         if( std::isfinite(pdf_far) && std::isfinite(pdf_inv) )
         {
            const real_t tolerance = real_c(1e-4) * lengthScalingFactor_;
            mesa_pd::kernel::SingleCast singleCast;
            mesa_pd::RayParticleIntersectionRatioFunctor intersectionRatioFctr;
            const real_t deltaRay = singleCast( particleIdx, *ac_, intersectionRatioFctr, *ac_, cellCenter, direction, tolerance );
            if( std::isfinite(deltaRay) && deltaRay >= real_c(0) && deltaRay <= real_c(1) )
            {
               delta = deltaRay;
               const real_t kappa0 = ( real_c(1) - real_c(2) * delta ) / ( real_c(1) + real_c(2) * delta ); // Table 3
               pdf_new = pdf_old + kappa0 * ( pdf_far - pdf_inv );                                           // Eq. (4.1)
               alpha   = real_c(4) / ( real_c(1) + real_c(2) * delta );                                      // Table 4
               applied = true;
            }
         }
      }
      if( applied ) { ++ibbAppliedCount; } else { ++ibbFallbackCount; }
      static uint_t ibbLogCount = uint_t(0);
      if( ibbLogCount < uint_t(6) && ( applied ? ( ibbAppliedCount == uint_t(1) || ibbAppliedCount == uint_t(1000) || ibbAppliedCount == uint_t(100000) ) : ( ibbFallbackCount == uint_t(1) || ibbFallbackCount == uint_t(1000) ) ) )
      {
         ++ibbLogCount;
         WALBERLA_LOG_INFO("[IBB] interpolated bounce-back active on this rank: links so far interpolated=" << ibbAppliedCount << ", fell back to simple bounce-back=" << ibbFallbackCount
            << " | this link (" << x << "," << y << "," << z << ") dir=" << uint_t(dir) << " " << ( applied ? "delta=" : "fallback (delta=" ) << delta << ( applied ? "" : ")" ) );
      }
   }

   // assumed boundary position
   const Vector3<real_t> boundaryPosition( cellCenter + delta * direction );

   auto boundaryVelocity = mesa_pd::getVelocityAtWFPoint(particleIdx, *ac_, boundaryPosition );

   // include effect of boundary velocity
   // density is declared outside the branch so [NAN-LINK] below can report it (NaN on the incompressible path). It
   // sums all populations of the cell, so it can carry a corrupted population that the pdf_old checks above miss.
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

   const real_t rawExchange = pdf_old + pdf_new; // momentum budget: before the reference subtraction below

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

   // Diagnostic: logs all intermediate quantities when this link's force is non-finite although pdf_old passed the
   // checks above (e.g. via density), first 50 hits per rank.
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

   // Three-phase MEM fill weighting (Bogner 2017, Eq. 5.9): the force above treats the near-boundary cell as full
   // liquid, but at a floater's waterline it is an interface cell that is partly gas, and full weight overstates the
   // momentum exchanged there. The link force (Wen correction included) is scaled by the clamped fill level φ(x_b) of
   // that cell. With the rho_ref = 1 subtraction above (compressible, no hydrostaticDensityFct_) the momentum-exchange
   // part becomes c_q φ [f̃ + f - 2 w_q], which is Eq. 5.9 with gas density 1, c_q φ [f̃ + f] + 2 c_q (1 - φ) w_q,
   // minus the constant 2 c_q w_q; gas cells, whose links are not treated (domain mask), give 0, the same expression
   // at φ = 0. Only the force on the particle is weighted; the PDF written back into the fluid is not.
   if( fillLevelFct_ )
   {
      const real_t fillRaw    = fillLevelFct_( nearBoundaryCell );
      const real_t fillWeight = std::clamp( fillRaw, real_c(0), real_c(1) );

      // Diagnostic: logs links whose force the weighting reduces (weight < 0.99), first 60 hits per rank.
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

   if( MemLinkMomentumBudget::enabled() )
   {
      const real_t phi = fillLevelFct_ ? std::clamp( fillLevelFct_( nearBoundaryCell ), real_c(0), real_c(1) ) : real_c(1);
      const real_t cx  = real_c( stencil::cx[dir] ), cz = real_c( stencil::cz[dir] );
      const real_t ref = rawExchange - forceMEM; // 2 rho_ref w_q (0 on the incompressible path)
      auto& s = MemLinkMomentumBudget::sums();
      s[0] += cx * rawExchange;                         s[1] += cz * rawExchange;
      s[2] += phi * cx * rawExchange;                   s[3] += phi * cz * rawExchange;
      s[4] += phi * cx * ref;                           s[5] += phi * cz * ref;
      s[6] -= phi * correction * boundaryVelocity[0];   s[7] -= phi * correction * boundaryVelocity[2];
      s[8] += force[0];                                 s[9] += force[2];
   }

   force *= forceScalingFactor_;

   // add the force onto the particle at the obstacle boundary
   lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic( particleIdx, *ac_, force, boundaryPosition );

}

} // namespace lbm_mesapd_coupling
} // namespace walberla
