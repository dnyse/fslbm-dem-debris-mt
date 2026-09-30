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
//! Small, isolated regression test for the 2026-09-27 THESIS_PLAN.md origin finding and its fix:
//! tracing a real crashed run's own box-logger data (job 4271714) found that cell
//! (603.5,20.5,27.5) was owned by rod uid...162 through step 21936 (rock-steady, density~1.003),
//! that ownership switched DIRECTLY to rod uid...780 at step 21938 (the cell stayed flagged obstacle
//! throughout -- no intervening liquid/interface state), and a growing odd-even (checkerboard)
//! oscillation began in that SAME step and never stopped until the eventual blow-up and crash.
//! Reading lbm_mesapd_coupling::MovingParticleMapping's own "cell is already an obstacle (maybe from
//! another particle)" branch found why: when a second particle's containment test finds a cell
//! already flagged obstacle (owned by a first, non-fixed particle), the code ONLY updated
//! `(*particleField)(x,y,z) = particleUid` -- it never touched the PDF field at all, so the second
//! particle silently inherited whatever the first particle's own obstacle dynamics left behind, with
//! zero re-initialization.
//!
//! **Fix** (MovingParticleMapping.h, same branch): on a GENUINE handoff (formerParticleUid !=
//! particleUid -- NOT the routine same-owner re-confirmation that fires every step a particle keeps
//! a cell it already owns), the cell's PDF is now also reset to rest equilibrium via PdfField's own
//! stencil-agnostic setToEquilibrium(), breaking the silent-inheritance chain at its source.
//!
//! This test exercises the REAL, patched MovingParticleMapping class (not a reimplementation) with
//! two real mesa_pd Sphere particles to prove the fix directly: particle A maps a cell, an arbitrary
//! (deliberately non-equilibrium, i.e. "whatever A's own imperfect obstacle dynamics might have left
//! behind") PDF is written into it, then A moves away and B moves in during the SAME mapping sweep
//! call (matching the real per-main-step cadence, where MovingParticleMapping's operator() iterates
//! every particle once per step) -- reproducing a direct, same-step handoff. Checked: ownership
//! (particleField) switches to B's uid, AND the cell's PDF is now rest equilibrium (rho=1, u=0),
//! NOT the stale injected values -- confirming the fix actually breaks the inheritance. A second
//! check confirms the OTHER, routine case (same particle re-confirming a cell it already owns) does
//! NOT reset the PDF -- the fix must not turn every covered cell into a permanently-frozen one.
//!
//! Scope note: this test deliberately does NOT attempt to also reproduce the subsequent exponential
//! growth phase in equal isolation. `CollideIsolatedObstacleIterationTest.cpp` (2026-09-25) already
//! tested repeated ISOLATED collide iteration (zero force, zero streaming input) on two different
//! real near-equilibrium seeds and found NEITHER diverges -- a genuine fixed point in both cases. So
//! whatever drives the real, observed growth needs more than "collide alone, forever" -- most likely
//! the interplay with SimpleBB's own per-step partial PDF rewrite (2026-09-25 finding: SimpleBB
//! writes directly into an obstacle cell's own PDF array every step, for the links facing fluid
//! neighbours) against a genuinely time-varying real neighbour flow field, which a single isolated
//! cell with no real neighbours cannot reproduce. That would need a heavier, less "isolated" two-rod
//! real-time reproduction (a small full scenario, not a unit test) -- not attempted here.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition itself -- match ClipAnalyticGeometryTest.cpp's own workaround.
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

// flag_t, FlagField_T, ParticleAccessor_T, ParticleField_T are already defined by SimulationObjects.h
// (pulled in transitively via ParticleKernels.h above) -- reused directly rather than redeclared,
// both to avoid a conflicting-declaration error and because these ARE exactly the real production
// types (flag_t=uint32_t, ParticleAccessor_T=ParticleAccessorWithBaseShape).

const FlagUID Fluid_Flag( "fluid" );
const FlagUID Obstacle_Flag( "moving obstacle" );
const FlagUID FormerObstacle_Flag( "former moving obstacle" );

using SBB_T             = lbm_mesapd_coupling::SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
using BoundaryHandling_T = BoundaryHandling< FlagField_T, LatticeModel_T::Stencil, SBB_T >;

// Minimal single-block factory, mirroring IbbEquivalenceDragTest.cpp's own MyBoundaryHandling (that
// test already exercises MovingParticleMappingKernel + SimpleBB together; this reuses the identical
// construction pattern for the per-step MovingParticleMapping class instead).
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

   // Tiny single-block domain -- only large enough to hold both spheres and a comfortable margin.
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

   // Both particles exist from the start; B starts far outside the domain so only A maps the test
   // cell in the FIRST mapping call -- matching the real scenario where both rods already exist and
   // simply move relative to each other over time.
   const Vector3< real_t > testCellCentre( real_t( 5.5 ), real_t( 5.5 ), real_t( 5.5 ) );
   const Vector3< real_t > farAway( real_t( 1000 ), real_t( 1000 ), real_t( 1000 ) );

   walberla::id_t uidA, uidB;
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition( testCellCentre );
      p.setInteractionRadius( radius );
      p.setOwner( 0 );
      p.setBaseShape( std::make_shared< mesa_pd::data::Sphere >( radius ) ); // matches SetupStages.h's own pattern
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

   // [Claude, 2026-09-27] The reset itself now lives in an OPT-IN callback (MovingParticleMapping's
   // onOwnershipHandoff_, default no-op) rather than being unconditional inside the generic class --
   // see MovingParticleMapping.h's own doc comment for why (the generic lbm::Equilibrium<>::set()/
   // PdfField::setToEquilibrium() API turned out to be unsafe to call from that generic file for
   // Antidunes' generated Cumulant model in the real multi-file production link, confirmed by direct
   // bisection on the cluster: jobs 4271763/4271773 segfaulted with it inside the generic class,
   // job 4271777 -- same code, callback removed -- got the ordinary STALE-PARTICLE-IDX abort
   // instead). This test's own callback mirrors BuildTimeloop.h's real one exactly (manual D3Q27
   // weight assignment, the same convention repairDegenerateObstacleCells already uses), so it
   // verifies the SAME reset Antidunes itself now performs, not a simplified stand-in.
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

   // Step 1: map particle A only (B is far away) -- the test cell should become obstacle, owned by A.
   mapping( block );

   check( boundaryHandling->isBoundary( ix, iy, iz ), "the test cell must be an obstacle after mapping particle A" );
   check( ( *particleField )( ix, iy, iz ) == uidA, "the test cell must be owned by particle A after the first mapping" );

   // Inject a deliberately NON-equilibrium PDF -- standing in for "whatever A's own imperfect
   // obstacle dynamics (unmasked collide, SimpleBB's own per-step partial rewrite) might have left
   // behind after some real time under A's ownership". Base state: D3Q27 rest equilibrium (rho=1,
   // u=0; same weights used throughout this session's own repairDegenerateObstacleCells/
   // trackWatchCells code: w0=8/27, w[1-6]=2/27, w[7-18]=1/54, w[19-26]=1/216), perturbed on one
   // arbitrary component pair -- the specific direction/magnitude is not the point; ANY nonzero,
   // asymmetric deviation from equilibrium is what a real, slightly-imperfect obstacle-cell history
   // would leave behind, and the claim under test is that the handoff preserves it UNCHANGED,
   // whatever it is.
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

   // Step 2: move A far away, move B into the SAME position A just vacated, then call the SAME
   // mapping sweep ONCE MORE -- reproducing a direct, same-step handoff (the cell is still flagged
   // obstacle throughout; A's own pass this step no longer visits it at all since A moved away, so
   // only B's pass touches it, hitting the "cell is already an obstacle (maybe from another
   // particle)" branch, MovingParticleMapping.h:180-193).
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

   // Expected D3Q27 rest equilibrium (rho=1, u=0) -- computed the SAME manual way the callback
   // itself does (and repairDegenerateObstacleCells elsewhere in this codebase), deliberately NOT
   // via PdfField::setToEquilibrium() -- see this file's own header comment for why that generic
   // API is not used here at all any more.
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
   // NOTE: particle B's own sphere (radius 2.5) covers a whole VOLUME of cells centred on
   // testCellCentre, not just the single tracked (ix,iy,iz) -- and A occupied that SAME volume in
   // step 1, so EVERY cell in it undergoes a genuine handoff here, not only the one this test
   // otherwise tracks. That is correct, expected behaviour (matching two real overlapping rods
   // handing off their whole shared region, not one point) -- so ">= 1", not "== 1", is the right
   // check; the important thing is that it fired AT LEAST once (the callback is wired up and
   // reachable), not an exact count for an incidental sphere volume.
   const uint_t handoffInvocationsAfterGenuineHandoff = handoffCallbackInvocations;
   check( handoffInvocationsAfterGenuineHandoff >= uint_t( 1 ),
          "FIX CHECK: the handoff callback must fire at least once for this genuine handoff -- zero "
          "would mean the callback is not wired up or not being invoked." );

   // Second scenario: the ROUTINE case -- particle B re-confirms a cell it ALREADY owns on a later
   // step (formerParticleUid == particleUid, no real handoff). Inject a fresh non-equilibrium PDF
   // and re-run the mapping sweep with B still in place (no position change) -- the fix must NOT
   // reset this, or every covered cell would be permanently frozen at rest equilibrium forever.
   std::array< real_t, 27 > secondInjection{};
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q) secondInjection[q] = injectedPdf[q] + real_t( 1e-4 ) * real_c( q );
   for (uint_t q = uint_t( 0 ); q < uint_t( 27 ); ++q) pdfField->get( ix, iy, iz, q ) = secondInjection[q];
   mapping( block ); // B has not moved -- this must be the routine "re-confirm my own cell" pass

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
