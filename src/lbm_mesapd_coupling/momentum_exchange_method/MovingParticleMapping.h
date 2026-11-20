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
//! \file MovingParticleMapping.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/cell/Cell.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagField.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/mapping/ParticleBoundingBox.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/IAccessor.h"
#include "mesa_pd/kernel/SingleCast.h"
#include "mesa_pd/data/Flags.h"

#include <functional>
#include <limits>

namespace walberla {
namespace lbm_mesapd_coupling {

/*!\brief Maps the moving particles into the simulation domain and updates the mapping
 *
 * Cells that are inside a particle, will be marked with the 'obstacle' flag.
 * 'Inside' means that the cell center is contained in the particle.
 * Thus, a containsPoint function has to be available for all particles/shapes.
 *
 * Cells that in the last time step were inside the particle but are now outside of it, i.e. the particle has moved,
 * will be marked with the 'formerObstacle' flag.
 * The 'formerObstacle' flag is used in a second step by the PDFReconstruction class (see PDFReconstruction.h) to
 * re-initialize the missing PDFs. Afterwards, the 'formerObstacle' flag is removed and the 'fluid' flag is set.
 *
 * It is not strictly required that the mapping has been initialized with the mapping kernel from below.
 *
 * The 'mappingParticleSelector' can be used to include/exclude certain particles from the mapping.
 *
 * cellCenterTestOffset (fraction of a cell, default 0 = off): the containment test uses the cell centre shifted by
 * cellCenterTestOffset * (dx, 1.3 dy, 1.7 dz). 'Inside' is inclusive, so a body face lying exactly on cell centres maps
 * an extra cell layer, which floating-point noise in the body's pose can flip from step to step; a small fixed offset
 * (e.g. 1e-4) resolves every such tie the same way. Use the same offset in MovingParticleMappingKernel (initial
 * mapping). Cell centres used for the momentum bookkeeping are not shifted.
 */
template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename ParticleSelector_T>
class MovingParticleMapping
{
   static_assert(std::is_base_of_v<mesa_pd::data::IAccessor, ParticleAccessor_T>, "Provide a valid accessor as template");

public:

   using FlagField_T = typename BoundaryHandling_T::FlagField;
   using flag_t = typename BoundaryHandling_T::flag_t;

   // Called right before a cell's obstacle ownership passes from one particle to another (not on the per-step
   // re-confirmation by its current owner), only for cells in the PDF field's inner part:
   // (block, x, y, z, formerUid, newUid). No-op by default. This class never writes PDFs itself, since it is used with
   // other lattice models too; a caller that wants the handed-off cell reset does so with its own model (Antidunes:
   // BuildTimeloop.h assigns the D3Q27 rest weights by hand, as SetupStages.h's repairDegenerateObstacleCells does).
   using OwnershipHandoffCallback_T =
      std::function< void( IBlock&, cell_idx_t, cell_idx_t, cell_idx_t, walberla::id_t, walberla::id_t ) >;
   // Diagnostic: called for an inner domain cell right before the particle covers it (block, x, y, z, particleIdx);
   // the cell's PDFs and flags are still the fluid's. No-op by default.
   using NewlyCoveredCallback_T = std::function< void( IBlock&, cell_idx_t, cell_idx_t, cell_idx_t, size_t ) >;

   MovingParticleMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                          const BlockDataID & pdfFieldID,
                          const BlockDataID & boundaryHandlingID,
                          const BlockDataID & particleFieldID,
                          const shared_ptr<ParticleAccessor_T>& ac,
                          const FlagUID & obstacle, const FlagUID & formerObstacle,
                          const ParticleSelector_T& mappingParticleSelector,
                          bool conserveMomentumWhenMapping, real_t cellCenterTestOffset = real_t(0),
                          OwnershipHandoffCallback_T onOwnershipHandoff = OwnershipHandoffCallback_T(),
                          NewlyCoveredCallback_T onNewlyCovered = NewlyCoveredCallback_T())
   : blockStorage_( blockStorage ), pdfFieldID_( pdfFieldID ), boundaryHandlingID_( boundaryHandlingID ),
     particleFieldID_( particleFieldID ), ac_( ac ),
     obstacle_( obstacle ), formerObstacle_( formerObstacle ), mappingParticleSelector_( mappingParticleSelector ),
     conserveMomentumWhenMapping_( conserveMomentumWhenMapping ), cellCenterTestOffset_( cellCenterTestOffset ),
     onOwnershipHandoff_( onOwnershipHandoff ), onNewlyCovered_( onNewlyCovered )
   {}

   void inline operator()( IBlock * const block )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      PdfField_T * pdfField = block->getData< PdfField_T >( pdfFieldID_ );
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );
      auto * flagField = boundaryHandling->getFlagField();
      ParticleField_T * particleField = block->getData< ParticleField_T >( particleFieldID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );
      WALBERLA_ASSERT_NOT_NULLPTR( particleField );
      WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), particleField->xyzSize() );
      WALBERLA_ASSERT( flagField->flagExists( obstacle_ ) );

      const flag_t       obstacle = flagField->getFlag( obstacle_ );
      const flag_t formerObstacle = flagField->flagExists( formerObstacle_ ) ? flagField->getFlag( formerObstacle_ ) :
                                    flagField->registerFlag( formerObstacle_ );

      const real_t dx = blockStorage_->dx( blockStorage_->getLevel(*block) );
      const real_t dy = blockStorage_->dy( blockStorage_->getLevel(*block) );
      const real_t dz = blockStorage_->dz( blockStorage_->getLevel(*block) );

      for( size_t idx = 0; idx < ac_->size(); ++idx )
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            mapParticleAndUpdateMapping( idx, block, pdfField, boundaryHandling, flagField, particleField, obstacle, formerObstacle, dx, dy, dz );
         }
      }
   }

private:

   void inline mapParticleAndUpdateMapping(const size_t particleIdx, IBlock * const block,
                                           PdfField_T * pdfField, BoundaryHandling_T * boundaryHandling, FlagField_T * flagField, ParticleField_T * particleField,
                                           const flag_t & obstacle, const flag_t & formerObstacle,
                                           real_t dx, real_t dy, real_t dz)
   {
      // policy: every particle manages only its own flags

      auto cellBB = getParticleCellBB( particleIdx, *ac_, *block, *blockStorage_, flagField->nrOfGhostLayers() );

      Vector3<real_t> startCellCenter = blockStorage_->getBlockLocalCellCenter( *block, cellBB.min() );

      mesa_pd::kernel::SingleCast singleCast;
      mesa_pd::ContainsPointFunctor containsPointFctr;

      auto particleUid = ac_->getUid(particleIdx);

      Vector3<real_t> currentCellCenter = startCellCenter;
      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
      {
         currentCellCenter[1] = startCellCenter[1];
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
         {
            currentCellCenter[0] = startCellCenter[0];
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {

               flag_t & cellFlagPtr = flagField->get(x,y,z);

               const Vector3<real_t> testPoint = currentCellCenter +
                  Vector3<real_t>( cellCenterTestOffset_ * dx, real_t(1.3) * cellCenterTestOffset_ * dy, real_t(1.7) * cellCenterTestOffset_ * dz );
               if( singleCast(particleIdx, *ac_, containsPointFctr, *ac_, testPoint) )
               {
                  // cell is inside particle, now check the flag of the cell
                  // and let this cell refer to the present particle (via Uid)

                  // cell is currently fluid and now newly occupied by the particle
                  if( boundaryHandling->isDomain(x,y,z))
                  {
                     if( onNewlyCovered_ && pdfField->isInInnerPart(Cell(x,y,z)) ) onNewlyCovered_( *block, x, y, z, particleIdx );
                     // set obstacle flag
                     boundaryHandling->forceBoundary( obstacle, x, y, z );
                     (*particleField)(x, y, z) = particleUid;

                     if( conserveMomentumWhenMapping_ && pdfField->isInInnerPart(Cell(x,y,z)) ) {
                        // this removes the fluid from the simulation, together with the containing momentum
                        // to ensure momentum conservation, the momentum of this fluid cell has to be added to the particle
                        // see Aidun, C. K., Lu, Y., & Ding, E. J. (1998). Direct analysis of particulate suspensions with inertia using the discrete Boltzmann equation. Journal of Fluid Mechanics, 373, 287-311.

                        // force = momentum / dt, with dt = 1
                        Vector3<real_t> force = pdfField->getMomentumDensity(x, y, z);
                        lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic( particleIdx, *ac_, force, currentCellCenter );
                     }

                  }
                  else if( isFlagSet( cellFlagPtr, obstacle ) )
                  {
                     // cell is already an obstacle (maybe from another particle)
                     auto formerParticleUid = (*particleField)(x, y, z);
                     auto formerParticleIdx = ac_->uidToIdx(formerParticleUid);

                     // uidToIdx() returns SIZE_MAX if the former owner is no longer known locally (its ghost copy
                     // was dropped). Such an owner cannot be fixed (fixed particles are never dropped), so it is
                     // treated as not fixed and the cell is claimed, instead of indexing getFlags() out of bounds.
                     const bool formerParticleKnownLocally = formerParticleIdx != std::numeric_limits< size_t >::max();
                     if( !formerParticleKnownLocally ||
                         !isSet(ac_->getFlags(formerParticleIdx), mesa_pd::data::particle_flags::FIXED) )
                     {
                        // only claim this cell if it not from a fixed particle, i.e. the fixed particle keeps the cell in all cases
                        // this is IMPORTANT to not mess up the mapping of cells that are contained inside fixed cells but are sporadically
                        // occupied by a moving particle as well (during collision and the there present overlap),
                        // since the fixed particle usually do not get remapped in every timestep

                        // The current owner re-enters this branch every step (formerParticleUid == particleUid);
                        // only a genuine handoff to another particle calls onOwnershipHandoff_, so the caller can
                        // reset the PDF history the cell would otherwise carry over from its former owner (see
                        // tests/showcases_antidunes/ObstacleOwnershipHandoffTest.cpp). cellBB spans the flag
                        // field's ghost layers, which can outnumber the PDF field's (Antidunes: 2 vs. 1), hence
                        // the same isInInnerPart guard as for the getMomentumDensity read above.
                        if( formerParticleUid != particleUid && onOwnershipHandoff_ && pdfField->isInInnerPart( Cell( x, y, z ) ) )
                        {
                           onOwnershipHandoff_( *block, x, y, z, formerParticleUid, particleUid );
                        }
                        (*particleField)(x, y, z) = particleUid;
                     }
                  }
                  else if( isFlagSet( cellFlagPtr, formerObstacle ) )
                  {
                     // cell is a former obstacle cell (maybe from another particle that has moved away)
                     boundaryHandling->setBoundary( obstacle, x, y, z );
                     removeFlag( cellFlagPtr, formerObstacle );
                     (*particleField)(x, y, z) = particleUid;
                  }

                  // else: IMPORTANT in all other cases, another boundary flag is already present in the cell
                  // this could e.g. be an inflow boundary condition
                  // we do NOT want to overwrite this boundary condition and thus skip it, even though the particle overlaps

               }
               else
               {
                  // cell is outside particle
                  if( isFlagSet( cellFlagPtr, obstacle ) && ((*particleField)(x,y,z) == particleUid) )
                  {
                     // cell was previously occupied by this particle
                     boundaryHandling->removeBoundary( obstacle, x, y, z );
                     addFlag( cellFlagPtr, formerObstacle );
                     // entry at (*particleField)(x,y,z) should still point to the previous particle.
                     // If during initialization the overlap between neighboring blocks
                     // was chosen correctly/large enough, the particle should still be on this block.
                     // The particle information is needed in the PDF reconstruction step.
                     // There, the flag will be removed and replaced by a domain flag after reconstruction.
                  }
               }
               currentCellCenter[0] += dx;
            }
            currentCellCenter[1] += dy;
         }
         currentCellCenter[2] += dz;
      }
   }

   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;

   const FlagUID obstacle_;
   const FlagUID formerObstacle_;

   ParticleSelector_T mappingParticleSelector_;

   bool conserveMomentumWhenMapping_;
   real_t cellCenterTestOffset_;
   OwnershipHandoffCallback_T onOwnershipHandoff_;
   NewlyCoveredCallback_T onNewlyCovered_;

};


template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename ParticleSelector_T>
MovingParticleMapping< PdfField_T, BoundaryHandling_T, ParticleAccessor_T, ParticleSelector_T>
makeMovingParticleMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                           const BlockDataID & pdfFieldID,
                           const BlockDataID & boundaryHandlingID,
                           const BlockDataID & particleFieldID,
                           const shared_ptr<ParticleAccessor_T>& ac,
                           const FlagUID & obstacle, const FlagUID & formerObstacle,
                           const ParticleSelector_T& mappingParticleSelector,
                           bool conserveMomentumWhenMapping, real_t cellCenterTestOffset = real_t(0),
                           typename MovingParticleMapping<PdfField_T, BoundaryHandling_T, ParticleAccessor_T, ParticleSelector_T>::OwnershipHandoffCallback_T
                              onOwnershipHandoff = typename MovingParticleMapping<PdfField_T, BoundaryHandling_T, ParticleAccessor_T, ParticleSelector_T>::OwnershipHandoffCallback_T(),
                           typename MovingParticleMapping<PdfField_T, BoundaryHandling_T, ParticleAccessor_T, ParticleSelector_T>::NewlyCoveredCallback_T
                              onNewlyCovered = typename MovingParticleMapping<PdfField_T, BoundaryHandling_T, ParticleAccessor_T, ParticleSelector_T>::NewlyCoveredCallback_T())
{
   return MovingParticleMapping<PdfField_T, BoundaryHandling_T, ParticleAccessor_T, ParticleSelector_T>
         (blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac, obstacle, formerObstacle, mappingParticleSelector, conserveMomentumWhenMapping, cellCenterTestOffset, onOwnershipHandoff, onNewlyCovered);
}

////////////////////
// MAPPING KERNEL //
////////////////////


/*!\brief Kernel that maps all particles onto all blocks to the "moving obstacle" boundary condition, that requires the additional particleField
 *
 * Cells that are inside the particles are set to 'obstacle' and the UID of the particle is stored in the particleField for later reference.
 *
 * Example:
 * lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);
 * ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);
 *
 * This kernel is usually applied before the start of the simulation, after the particle and fluid data structures have been created.
 * As a result, a consistent mapping of the particles into the domain is obtained.
 * For the update of the mapping during time stepping, use the class MovingParticleMapping (above).
 *
 * This mapping can also be used for stationary particle if the force and torque on this particle are to be evaluated.
 */
template< typename BoundaryHandling_T>
class MovingParticleMappingKernel
{
public:
   MovingParticleMappingKernel(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID, const BlockDataID & particleFieldID,
                               real_t cellCenterTestOffset = real_t(0))
         : blockStorage_(blockStorage), boundaryHandlingID_(boundaryHandlingID), particleFieldID_(particleFieldID), cellCenterTestOffset_(cellCenterTestOffset){}

   template<typename ParticleAccessor_T>
   void operator()(const size_t particleIdx, const ParticleAccessor_T& ac, const FlagUID & obstacle )
   {
      static_assert(std::is_base_of_v<mesa_pd::data::IAccessor, ParticleAccessor_T>, "Provide a valid accessor as template");

      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         mapMovingParticleOnBlock( particleIdx, ac, *blockIt, obstacle );
      }
   }

private:

   template< typename ParticleAccessor_T >
   void mapMovingParticleOnBlock( const size_t particleIdx, const ParticleAccessor_T& ac,
                                  IBlock & block, const FlagUID & obstacle )
   {
      static_assert(std::is_base_of_v<mesa_pd::data::IAccessor, ParticleAccessor_T>, "Provide a valid accessor as template");

      WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage_->getBlockStorage()) );

      BoundaryHandling_T * boundaryHandling = block.getData< BoundaryHandling_T >( boundaryHandlingID_);
      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

      auto * flagField = boundaryHandling->getFlagField();
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );
      WALBERLA_ASSERT( flagField->flagExists( obstacle ) );

      const auto obstacleFlag = flagField->getFlag( obstacle );

      auto * particleField = block.getData< ParticleField_T >( particleFieldID_ );
      WALBERLA_ASSERT_NOT_NULLPTR( particleField );
      WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), particleField->xyzSize() );

      auto cellBB = getParticleCellBB(particleIdx, ac, block, *blockStorage_, flagField->nrOfGhostLayers() );

      if( cellBB.empty() ) return;

      mesa_pd::kernel::SingleCast singleCast;
      mesa_pd::ContainsPointFunctor containsPointFctr;

      Vector3<real_t> startCellCenter = blockStorage_->getBlockLocalCellCenter( block, cellBB.min() );
      auto blockLevel = blockStorage_->getLevel(block);
      const real_t dx = blockStorage_->dx( blockLevel );
      const real_t dy = blockStorage_->dy( blockLevel );
      const real_t dz = blockStorage_->dz( blockLevel );

      real_t cz = startCellCenter[2];
      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
      {
         real_t cy = startCellCenter[1];
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
         {
            real_t cx = startCellCenter[0];
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {
               if( singleCast(particleIdx, ac, containsPointFctr, ac,
                              Vector3<real_t>(cx + cellCenterTestOffset_ * dx, cy + real_t(1.3) * cellCenterTestOffset_ * dy, cz + real_t(1.7) * cellCenterTestOffset_ * dz)) )
               {
                  if(boundaryHandling->isDomain(x,y,z))
                  {
                     boundaryHandling->forceBoundary(obstacleFlag, x, y, z);
                     (*particleField)(x,y,z) = ac.getUid(particleIdx);
                  }
                  // no else -> will not overwrite pre-existing boundary cells, e.g. of previously applied applied mapping (e.g. when first mapping global particles)
               }
               cx += dx;
            }
            cy += dy;
         }
         cz += dz;
      }
   }

   shared_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID boundaryHandlingID_;
   BlockDataID particleFieldID_;
   real_t cellCenterTestOffset_;
};


} // namespace lbm_rpd_coupling
} // namespace walberla
