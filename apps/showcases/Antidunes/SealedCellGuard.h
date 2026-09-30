//======================================================================================================================
//
//  Experimental guard against sub-cell fluid "channels" between touching bodies.
//
//  [Claude, 2026-09-21] Three touching parallel rods of radius r enclose a channel of inscribed radius 0.1547 r
//  (0.46 LU for r = 3), i.e. less than one cell: the cell-centre mapping (MovingParticleMapping: a cell is obstacle if
//  its centre lies inside the body) leaves a string of fluid cells whose 26 D3Q27 neighbours are (almost) all obstacle.
//  Such a cell cannot exchange information with any fluid neighbour, its populations just bounce back and forth
//  (period-2 sign flip / density excursions, see THESIS_PLAN.md). This sweep closes them: every non-gas fluid cell
//  that has at least `minObstacleNeighbours` of its 26 neighbours covered by a moving obstacle or a no-slip wall (and
//  at least one moving obstacle) is turned into an obstacle cell of the neighbouring body that owns most of its
//  obstacle neighbours (particle field UID), exactly like the mapping does.
//
//  [Claude, 2026-09-22] CORRECTNESS FIX: this sweep now runs AFTER PdfReconstructionManager (BuildTimeloop.h), not
//  before it -- the earlier version intercepted 'former moving obstacle' cells that were STILL sealed and slammed
//  them straight back to obstacle, skipping reconstruction entirely; that cell's PDF/fill/excess-mass state was then
//  whatever it happened to be before it was ever mapped, frozen and never made consistent again. A cell doing this
//  repeatedly for dozens of steps (as happens in a real, densely packed rod pile, not in the small local test that
//  first validated this sweep) fed stale/diverging values into a neighbour's hydrodynamic-force link computation --
//  traced to the [MEM-GARBAGE-LINK]/NaN failure of job 4253834 (THESIS_PLAN.md 2026-09-22). Running AFTER
//  reconstruction means every cell this sweep can see has ALREADY been resolved to a fully consistent domain state
//  this step (freshly reconstructed just now if it was vacated, or unchanged from a previous step otherwise) -- so
//  the special 'formerObstacle, still sealed' case no longer exists; every closure below acts on a cell that is
//  simply an ordinary, valid liquid/interface cell at the moment it is closed.
//
//  Prm: ModelParameters.closeSealedCellsMinNeighbours (0 = off, default; typical 24). OFF by default -> no change for
//  any existing prm. Decisions are taken in two phases (collect, then apply) so the result does not depend on the
//  sweep order, and on the inner cells plus the first ghost layer (flag/particle fields have two ghost layers) so
//  that neighbouring blocks agree on the flags of their shared cells.
//
//======================================================================================================================
#pragma once

#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/cell/CellInterval.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/IBlock.h"

#include "field/FlagUID.h"

#include "lbm_mesapd_coupling/DataTypes.h"

#include <algorithm>
#include <map>
#include <memory>
#include <vector>

namespace walberla
{
namespace antidunes
{

template< typename BoundaryHandling_T >
class CloseSealedCellsSweep
{
 public:
   using FlagField_T = typename BoundaryHandling_T::FlagField;

   CloseSealedCellsSweep(const BlockDataID& boundaryHandlingID, const BlockDataID& particleFieldID,
                         const FlagUID& obstacle, const FlagUID& noSlip, const FlagUID& gas,
                         uint_t minObstacleNeighbours)
      : boundaryHandlingID_(boundaryHandlingID), particleFieldID_(particleFieldID), obstacle_(obstacle),
        noSlip_(noSlip), gas_(gas), minObstacleNeighbours_(minObstacleNeighbours),
        closedTotal_(std::make_shared< uint_t >(uint_t(0)))
   {}

   void operator()(IBlock* const block)
   {
      if (minObstacleNeighbours_ == uint_t(0)) return;

      BoundaryHandling_T* boundaryHandling = block->getData< BoundaryHandling_T >(boundaryHandlingID_);
      FlagField_T* flagField               = boundaryHandling->getFlagField();
      auto* particleField                  = block->getData< lbm_mesapd_coupling::ParticleField_T >(particleFieldID_);

      if (!flagField->flagExists(obstacle_)) return; // nothing has ever been mapped on this block
      const auto obstacle = flagField->getFlag(obstacle_);
      const auto noSlip   = flagField->flagExists(noSlip_) ? flagField->getFlag(noSlip_) : obstacle;
      const auto gas      = flagField->flagExists(gas_) ? flagField->getFlag(gas_) : obstacle;

      struct Closure
      {
         Cell cell;
         walberla::id_t uid;
      };
      std::vector< Closure > closures;

      CellInterval range = flagField->xyzSizeWithGhostLayer();
      range.expand(cell_idx_t(-1)); // inner cells + first ghost layer: every 26-neighbourhood is inside the allocation

      for (cell_idx_t z = range.zMin(); z <= range.zMax(); ++z)
         for (cell_idx_t y = range.yMin(); y <= range.yMax(); ++y)
            for (cell_idx_t x = range.xMin(); x <= range.xMax(); ++x)
            {
               // PdfReconstructionManager has already run this step (see BuildTimeloop.h), so every cell here is
               // already an ordinary, fully-resolved domain cell -- no 'formerObstacle' case can remain to intercept.
               if (!boundaryHandling->isDomain(x, y, z)) continue;
               if (flagField->isFlagSet(x, y, z, gas)) continue; // gas cells are harmless

               // count the non-solid neighbours and give up as soon as the threshold cannot be reached any more
               // (almost every cell is left after ~3 flag reads, so the sweep is cheap)
               const uint_t maxOpen = uint_t(26) - std::min(minObstacleNeighbours_, uint_t(26));
               uint_t solid = uint_t(0), open = uint_t(0);
               std::map< walberla::id_t, uint_t > owners;
               bool hopeless = false;
               for (cell_idx_t dz = -1; dz <= 1 && !hopeless; ++dz)
                  for (cell_idx_t dy = -1; dy <= 1 && !hopeless; ++dy)
                     for (cell_idx_t dx = -1; dx <= 1; ++dx)
                     {
                        if (dx == 0 && dy == 0 && dz == 0) continue;
                        const cell_idx_t nx = x + dx, ny = y + dy, nz = z + dz;
                        if (flagField->isFlagSet(nx, ny, nz, obstacle))
                        {
                           ++solid;
                           ++owners[(*particleField)(nx, ny, nz)];
                        }
                        else if (flagField->isFlagSet(nx, ny, nz, noSlip)) { ++solid; }
                        else if (++open > maxOpen) { hopeless = true; break; }
                     }

               if (hopeless || solid < minObstacleNeighbours_ || owners.empty()) continue;

               walberla::id_t bestUid = owners.begin()->first;
               uint_t bestCount       = uint_t(0);
               for (const auto& o : owners)
                  if (o.second > bestCount) { bestCount = o.second; bestUid = o.first; }

               closures.push_back(Closure{ Cell(x, y, z), bestUid });
            }

      for (const auto& c : closures)
      {
         boundaryHandling->forceBoundary(obstacle, c.cell.x(), c.cell.y(), c.cell.z());
         (*particleField)(c.cell.x(), c.cell.y(), c.cell.z()) = c.uid;

         ++(*closedTotal_);
         if (*closedTotal_ <= uint_t(20) || *closedTotal_ % uint_t(2000) == uint_t(0))
         {
            WALBERLA_LOG_WARNING("[SEALED-CLOSE] rank " << mpi::MPIManager::instance()->rank() << " block-local cell ("
                                 << c.cell.x() << "," << c.cell.y() << "," << c.cell.z()
                                 << ") fluid cell turned into obstacle of particle uid " << c.uid
                                 << " (closures on this rank so far: " << *closedTotal_ << ")");
         }
      }
   }

 private:
   BlockDataID boundaryHandlingID_;
   BlockDataID particleFieldID_;
   FlagUID obstacle_;
   FlagUID noSlip_;
   FlagUID gas_;
   uint_t minObstacleNeighbours_;
   std::shared_ptr< uint_t > closedTotal_;
};

} // namespace antidunes
} // namespace walberla
