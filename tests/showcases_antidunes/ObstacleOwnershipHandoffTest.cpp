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
//! \file ObstacleOwnershipHandoffTest.cpp
//! Regression test for the obstacle-to-obstacle ownership handoff in lbm_mesapd_coupling::MovingParticleMapping,
//! where an obstacle cell passes directly from one particle to another (no liquid/interface state in between). In the
//! "cell is already an obstacle (maybe from another particle)" branch, MovingParticleMapping itself only updates
//! (*particleField)(x,y,z) = particleUid and leaves the PDF field alone, so without a callback the new owner inherits
//! whatever the previous owner's obstacle dynamics left in the cell.
//!
//! On a genuine handoff (formerParticleUid != particleUid, not the routine same-owner re-confirmation of every step)
//! of an inner cell, MovingParticleMapping calls the optional onOwnershipHandoff_ callback (default empty). Antidunes'
//! callback in BuildTimeloop.h (prm switch useOwnershipHandoffPdfReset) resets the cell's PDF to D3Q27 rest
//! equilibrium.
//!
//! The test drives the real MovingParticleMapping with two mesa_pd spheres: A maps a cell, a non-equilibrium PDF is
//! written into it, then A moves away and B moves in before the next mapping call (a direct, same-step handoff).
//! Checked: ownership switches to B, the cell is still obstacle, the PDF equals rest equilibrium (rho=1, u=0) and
//! not the injected values, and the callback fired. A second check: the routine case (B re-confirming its own cell)
//! must not reset the PDF, otherwise every covered cell would be frozen at equilibrium.
//!
//! Scope: only the handoff itself is tested, not any instability growing from an inherited PDF; isolated collide
//! iteration (CollideIsolatedObstacle*Test) does not reproduce such growth.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition, same workaround as ClipAnalyticGeometryTest.cpp.
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h" // pulls in SimulationObjects.h -> LatticeModel_T, PdfField_T, VectorFieldFlattened_T

#include "blockforest/Initialization.h"

#include "boundary/BoundaryHandling.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"

#include "lbm/field/AddToStorage.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/kernel/ParticleSelector.h"

#include <core/Environment.h>
#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

#include <cmath>
#include <iomanip>

namespace walberla { namespace antidunes {

// flag_t, FlagField_T, ParticleAccessor_T, ParticleField_T come from SimulationObjects.h (via ParticleKernels.h):
// the production types (flag_t=uint32_t, ParticleAccessor_T=ParticleAccessorWithBaseShape), not redeclared.

const FlagUID Fluid_Flag( "fluid" );
const FlagUID Obstacle_Flag( "moving obstacle" );
const FlagUID FormerObstacle_Flag( "former moving obstacle" );

using SBB_T             = lbm_mesapd_coupling::SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
using BoundaryHandling_T = BoundaryHandling< FlagField_T, LatticeModel_T::Stencil, SBB_T >;

// Minimal single-block factory, same construction as IbbEquivalenceDragTest.cpp's MyBoundaryHandling.
class TestBoundaryHandling
{
public:
   TestBoundaryHandling( const BlockDataID& flagFieldID, const BlockDataID& pdfFieldID,
                         const BlockDataID& particleFieldID, const shared_ptr< ParticleAccessor_T >& ac )
      : flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac )
   {}

   BoundaryHandling_T* operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      auto* flagField     = block->getData< FlagField_T >( flagFieldID_ );
      auto* pdfField       = block->getData< PdfField_T >( pdfFieldID_ );
      auto* particleField = block->getData< lbm_mesapd_coupling::ParticleField_T >( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag )
                                                              : flagField->registerFlag( Fluid_Flag );

      auto* handling =
         new BoundaryHandling_T( "boundary handling", flagField, fluid,
            SBB_T( "SBB", Obstacle_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ) );
      handling->fillWithDomain( flagField->nrOfGhostLayers() );
      return handling;
   }

private:
   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;
   shared_ptr< ParticleAccessor_T > ac_;
};

void runTest()
{
   int failures = 0;
   auto check   = [&]( bool ok, const std::string& what ) {
      if (!ok)
      {
         ++failures;
         WALBERLA_LOG_INFO_ON_ROOT( "FAIL: " << what );
      }
   };

   // Single-block domain, just large enough for both spheres plus a margin.
   auto blocks = blockforest::createUniformBlockGrid( uint_t( 1 ), uint_t( 1 ), uint_t( 1 ), uint_t( 12 ),
                                                       uint_t( 12 ), uint_t( 12 ), real_t( 1 ), uint_t( 1 ),
                                                       uint_t( 1 ), uint_t( 1 ), false, false, false );

   auto forceFieldID =
      field::addToStorage< VectorFieldFlattened_T >( blocks, "force field", real_c( 0 ), field::fzyx, uint_c( 1 ) );
   const real_t omega = real_t( 1.9 );
   LatticeModel_T latticeModel( forceFieldID, omega );
   auto pdfFieldID = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, Vector3< real_t >( real_t( 0 ) ),
                                                 real_t( 1.0 ), field::fzyx );
   auto flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   auto ps       = std::make_shared< mesa_pd::data::ParticleStorage >( 2 );
   auto accessor = std::make_shared< ParticleAccessor_T >( ps );
   const real_t radius = real_t( 2.5 );

   auto particleFieldID = field::addToStorage< lbm_mesapd_coupling::ParticleField_T >(
      blocks, "particle field", accessor->getInvalidUid(), field::fzyx, uint_t( 1 ) );

   // Both particles exist from the start; B starts far outside the domain, so only A maps the test cell in the
   // first mapping call.
   const Vector3< real_t > testCellCentre( real_t( 5.5 ), real_t( 5.5 ), real_t( 5.5 ) );
   const Vector3< real_t > farAway( real_t( 1000 ), real_t( 1000 ), real_t( 1000 ) );

   walberla::id_t uidA, uidB;
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition( testCellCentre );
      p.setInteractionRadius( radius );
      p.setOwner( 0 );
      p.setBaseShape( std::make_shared< mesa_pd::data::Sphere >( radius ) ); // as in SetupStages.h
      uidA = p.getUid();
   }
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition( farAway );
      p.setInteractionRadius( radius );
      p.setOwner( 0 );
      p.setBaseShape( std::make_shared< mesa_pd::data::Sphere >( radius ) );
      uidB = p.getUid();
   }

   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
      TestBoundaryHandling( flagFieldID, pdfFieldID, particleFieldID, accessor ), "boundary handling" );

   IBlock* block             = &( *blocks->begin() );
   auto* pdfField           = block->getData< PdfField_T >( pdfFieldID );
   auto* particleField      = block->getData< lbm_mesapd_coupling::ParticleField_T >( particleFieldID );
   auto* boundaryHandling   = block->getData< BoundaryHandling_T >( boundaryHandlingID );

   const cell_idx_t ix = cell_idx_t( 5 ), iy = cell_idx_t( 5 ), iz = cell_idx_t( 5 ); // == testCellCentre

   // The reset is an opt-in callback (onOwnershipHandoff_), so the generic class makes no PDF/stencil assumption.
   // This callback mirrors BuildTimeloop.h's: manual D3Q27 weights (as in repairDegenerateObstacleCells) instead of
   // the generic setToEquilibrium() (see MovingParticleMapping.h).
   uint_t handoffCallbackInvocations = uint_t( 0 );
   auto onOwnershipHandoff = [&]( IBlock&, cell_idx_t hx, cell_idx_t hy, cell_idx_t hz,
                                   walberla::id_t, walberla::id_t ) {
      ++handoffCallbackInvocations;
      pdfField->get( hx, hy, hz, uint_t( 0 ) ) = real_t( 8.0 / 27.0 );
      for (uint_t q = uint_t( 1 ); q <= uint_t( 6 ); ++q) pdfField->get( hx, hy, hz, q ) = real_t( 2.0 / 27.0 );
      for (uint_t q = uint_t( 7 ); q <= uint_t( 18 ); ++q) pdfField->get( hx, hy, hz, q ) = real_t( 1.0 / 54.0 );
      for (uint_t q = uint_t( 19 ); q <= uint_t( 26 ); ++q) pdfField->get( hx, hy, hz, q ) = real_t( 1.0 / 216.0 );
   };

   lbm_mesapd_coupling::MovingParticleMapping< PdfField_T, BoundaryHandling_T, ParticleAccessor_T,
                                                mesa_pd::kernel::SelectAll >
      mapping( blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, Obstacle_Flag,
               FormerObstacle_Flag, mesa_pd::kernel::SelectAll(), false, real_t( 0 ), onOwnershipHandoff );

   // Step 1: map particle A only (B is far away): the test cell becomes obstacle, owned by A.
   mapping( block );

   check( boundaryHandling->isBoundary( ix, iy, iz ), "the test cell must be an obstacle after mapping particle A" );
   check( ( *particleField )( ix, iy, iz ) == uidA, "the test cell must be owned by particle A after the first mapping" );

   // Inject a non-equilibrium PDF, standing in for what A's obstacle dynamics may leave behind: D3Q27 rest
   // equilibrium (w0=8/27, w[1-6]=2/27, w[7-18]=1/54, w[19-26]=1/216) with one asymmetric perturbation. Any nonzero
   // deviation would do; the question is whether the handoff keeps it unchanged.
   std::array< real_t, 27 > injectedPdf{};
   injectedPdf[0] = real_t( 8.0 / 27.0 );
   for (uint_t q = uint_t( 1 ); q <= uint_t( 6 ); ++q) injectedPdf[q] = real_t( 2.0 / 27.0 );
   for (uint_t q = uint_t( 7 ); q <= uint_t( 18 ); ++q) injectedPdf[q] = real_t( 1.0 / 54.0 );
   for (uint_t q = uint_t( 19 ); q <= uint_t( 26 ); ++q) injectedPdf[q] = real_t( 1.0 / 216.0 );
   injectedPdf[1] += real_t( 1e-3 );  // arbitrary asymmetric perturbation, one opposite-direction pair
   injectedPdf[2] -= real_t( 1e-3 );
   injectedPdf[19] += real_t( 5e-4 );
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q)
      pdfField->get( ix, iy, iz, q ) = injectedPdf[q];

   // Step 2: move A away, move B to the position A vacated, and call the mapping once more: a direct handoff. The cell
   // stays obstacle; only B's pass visits it, hitting the "cell is already an obstacle (maybe from another particle)"
   // branch of MovingParticleMapping.h.
   {
      const size_t idxA = accessor->uidToIdx( uidA );
      const size_t idxB = accessor->uidToIdx( uidB );
      accessor->setPosition( idxA, farAway );
      accessor->setPosition( idxB, testCellCentre );
   }
   mapping( block );

   check( ( *particleField )( ix, iy, iz ) == uidB,
          "after the handoff, the test cell must be owned by particle B (the SAME cell, direct "
          "obstacle-to-obstacle ownership transfer, no intervening liquid/interface state)" );
   check( boundaryHandling->isBoundary( ix, iy, iz ),
          "the test cell must STILL be flagged obstacle after the handoff (never became fluid in between)" );

   // Expected D3Q27 rest equilibrium (rho=1, u=0), set manually as in the callback (not via setToEquilibrium()).
   std::array< real_t, 27 > expectedEquilibrium{};
   expectedEquilibrium[0] = real_t( 8.0 / 27.0 );
   for (uint_t q = uint_t( 1 ); q <= uint_t( 6 ); ++q) expectedEquilibrium[q] = real_t( 2.0 / 27.0 );
   for (uint_t q = uint_t( 7 ); q <= uint_t( 18 ); ++q) expectedEquilibrium[q] = real_t( 1.0 / 54.0 );
   for (uint_t q = uint_t( 19 ); q <= uint_t( 26 ); ++q) expectedEquilibrium[q] = real_t( 1.0 / 216.0 );

   real_t maxAbsDiffFromInjected    = real_t( 0 );
   real_t maxAbsDiffFromEquilibrium = real_t( 0 );
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q)
   {
      const real_t after = pdfField->get( ix, iy, iz, q );
      maxAbsDiffFromInjected    = std::max( maxAbsDiffFromInjected, std::abs( after - injectedPdf[q] ) );
      maxAbsDiffFromEquilibrium = std::max( maxAbsDiffFromEquilibrium, std::abs( after - expectedEquilibrium[q] ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "after the handoff: max|PDF - injected(stale)| = " << std::setprecision( 12 )
                              << maxAbsDiffFromInjected << ", max|PDF - rest equilibrium| = "
                              << maxAbsDiffFromEquilibrium << ", handoff callback invocations = "
                              << handoffCallbackInvocations );
   check( maxAbsDiffFromEquilibrium <= real_t( 0 ),
          "FIX CHECK: on a genuine cross-particle handoff, the onOwnershipHandoff callback must reset "
          "the cell's PDF to rest equilibrium -- if this fails, the callback is missing, wrong, or "
          "not being invoked (MovingParticleMapping.h:180-193 / BuildTimeloop.h's own callback)." );
   check( maxAbsDiffFromInjected > real_t( 0 ),
          "FIX CHECK: the handoff must NOT leave the stale injected PDF in place -- if this fails, "
          "the fix did not actually change anything (the pre-fix bug this test previously proved is "
          "still present)." );
   // B's sphere (radius 2.5) covers a volume of cells that A also covered in step 1, so every cell in it is handed
   // off; hence ">= 1" invocations, not "== 1".
   const uint_t handoffInvocationsAfterGenuineHandoff = handoffCallbackInvocations;
   check( handoffInvocationsAfterGenuineHandoff >= uint_t( 1 ),
          "FIX CHECK: the handoff callback must fire at least once for this genuine handoff -- zero "
          "would mean the callback is not wired up or not being invoked." );

   // Second scenario, the routine case: B re-confirms a cell it already owns (formerParticleUid == particleUid).
   // Inject a fresh non-equilibrium PDF and map again without moving B; the PDF must not be reset.
   std::array< real_t, 27 > secondInjection{};
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q) secondInjection[q] = injectedPdf[q] + real_t( 1e-4 ) * real_c( q );
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q) pdfField->get( ix, iy, iz, q ) = secondInjection[q];
   mapping( block ); // B has not moved: routine re-confirmation pass

   check( ( *particleField )( ix, iy, iz ) == uidB, "particle B must still own the cell after re-confirming it" );
   real_t maxAbsDiffAfterReconfirm = real_t( 0 );
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q)
      maxAbsDiffAfterReconfirm = std::max( maxAbsDiffAfterReconfirm, std::abs( pdfField->get( ix, iy, iz, q ) - secondInjection[q] ) );
   WALBERLA_LOG_INFO_ON_ROOT( "after the SAME owner re-confirming its own cell: max|PDF - injected| = "
                              << std::setprecision( 12 ) << maxAbsDiffAfterReconfirm );
   check( maxAbsDiffAfterReconfirm <= real_t( 0 ),
          "FIX CHECK: the fix must be scoped to GENUINE handoffs only (formerParticleUid != "
          "particleUid) -- a particle re-confirming a cell it already owns every routine step must "
          "NOT have its PDF reset, or every covered cell would be permanently frozen at rest "
          "equilibrium forever, which is a much more drastic (and wrong) change than intended." );
   check( handoffCallbackInvocations == handoffInvocationsAfterGenuineHandoff,
          "FIX CHECK: the routine same-owner re-confirmation must NOT invoke the handoff callback "
          "again -- invocation count must stay exactly where it was after the earlier genuine "
          "handoff, with zero further increments." );

   WALBERLA_LOG_INFO_ON_ROOT(
      ( failures == 0
           ? "\nALL CHECKS PASSED -- the fix works exactly as intended: a GENUINE cross-particle "
             "handoff now resets the cell's PDF to rest equilibrium (breaking the silent-inheritance "
             "chain), while the SAME particle routinely re-confirming a cell it already owns leaves "
             "the PDF untouched, exactly as before.\n"
           : "\nSOME CHECKS FAILED -- see FAIL lines above.\n" ) );

   WALBERLA_CHECK( failures == 0 );
}

}} // namespace walberla::antidunes

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::debug::enterTestMode();
   walberla::mpi::MPIManager::instance()->useWorldComm();
   walberla::antidunes::runTest();
   return EXIT_SUCCESS;
}
