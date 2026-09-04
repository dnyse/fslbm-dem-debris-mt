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

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/domain/IDomain.h>

#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

#include <algorithm>

namespace walberla {
namespace mesa_pd {
namespace mpi {

/**
 * Contact filter for parallel execution of collision detection.
 *
 * This contact filter decides if a contact should be treated on this process.
 * Contact filtering has to be applied to avoid double treatment by multiple processes.
 */
class ContactFilter
{
public:
   ContactFilter() = default;
   ContactFilter(const ContactFilter& other) = default;
   ContactFilter(ContactFilter&& other) = default;
   ContactFilter& operator=(const ContactFilter& other) = default;
   ContactFilter& operator=(ContactFilter&& other) = default;

   /**
    * \param idx1 index of the first collision partner
    * \param idx2 index of the second collision partner
    * \param ac accessor data structure to access particle properties
    * \param contactPoint contact point of the two particles
    * \param domain domain datastructure
    * \return true if the contact should be treated by this process
    */
   template <typename Accessor>
   bool operator()( const size_t idx1,
                    const size_t idx2,
                    Accessor& ac,
                    const Vec3& contactPoint,
                    const domain::IDomain& domain ) const;

   // OFF (default) -> exact pre-existing behaviour: the local-remote/remote-remote (different-owner)
   // rule (marked `*` in the doc comment below) accepts iff domain.isContainedInProcessSubdomain(myRank_,
   // contactPoint) -- a PURELY LOCAL, PER-RANK geometric test using each rank's own independently
   // computed contactPoint. Two ranks holding slightly different local views of the same pair (ghost-sync
   // staleness, amplified by GJK/EPA's combinatorial sensitivity near ambiguous geometry -- see the
   // investigation this flag came out of) can each compute a different contactPoint that independently
   // satisfies its OWN rank's subdomain test, so both accept the same contact in the same cycle: a
   // ContactHistoryNotification duplicate-insert crash (confirmed: job 4166652, Cylinder x Box pair
   // straddling a rank boundary) and, more subtly, a silently doubled-then-partially-cancelled force
   // contribution even when it doesn't crash (see ForceTorqueNotification's sum-based reduce).
   // ON -> implements the doc comment's own documented-but-previously-unimplemented Option 2 `§` rule
   // (primary case only -- see below) for exactly this branch: accept iff myRank_ == the SMALLER of
   // ac.getOwner(idx1)/ac.getOwner(idx2). getOwner() is a stored, MPI-synced field (set once at particle
   // creation/migration, not recomputed from position every sub-cycle), so every rank witnessing this
   // pair computes the IDENTICAL answer -- no geometry, no possible disagreement, and (per
   // SyncNextNeighbors' own blockSyncExtension margin, sized to the largest possible particle in a
   // scenario) both true owners of any genuinely-touching pair are GUARANTEED to witness it via ghosting,
   // so the smaller-owner rank always actually runs this filter for that pair. Does NOT implement the
   // doc comment's secondary fallback ("if none of the owners witness the contact, the smallest rank
   // witnessing it") -- that requires knowing the full witnessing-rank set, which needs communication
   // this simple per-rank check doesn't have; a contact witnessed only by a third rank that owns neither
   // particle (a genuine multi-rank domain-corner configuration) is not covered. Believed rare/absent in
   // practice for real (non-degenerate) contacts given the ghosting guarantee above; flagged here as a
   // known, documented limitation rather than silently assumed away.
   bool useOwnerRankTiebreak = false;
private:
   uint_t myRank_ = uint_c( walberla::mpi::MPIManager::instance()->rank() );
};

template <typename Accessor>
inline bool ContactFilter::operator()( const size_t idx1,
                                       const size_t idx2,
                                       Accessor& ac,
                                       const Vec3& contactPoint,
                                       const domain::IDomain& domain ) const
{
   /* Contact filtering rules
    *
    * L: Local body
    * G: Global body
    * R: Remote body
    *
    * Option 1:              Option 2:
    * +---+---+---+---+      +---+---+---+---+
    * |   | L | G | R |      |   | L | G | R |
    * +---+---+---+---+      +---+---+---+---+
    * | L | + | + | * |      | L |§/+| + | § |
    * +---+---+---+---+      +---+---+---+---+
    * | G | + | ~ | - |      | G | + | ~ | - |
    * +---+---+---+---+      +---+---+---+---+
    * | R | * | - | # |      | R | § | - | § |
    * +---+---+---+---+      +---+---+---+---+
    *
    *  + Accept contact unconditionally
    *  - Reject contact unconditionally
    *  * Accept contact if we own the contact point
    *  # Accept contact if we own the contact point and the owners of the involved bodies are not the same
    *  ~ Accept contact only on root process
    *  § Accept contact if we are the owner with the smallest rank witnessing the contact or if none of the owners witness the contact and we are the process with the smallest rank witnessing the contact
    *
    * Note: Local-global contacts actually require a reduction of the contact reactions applied to the global body (unless it is fixed).
    * => MPI_Allreduce for all non-fixed global bodies before time-integration.
    */

   if( !data::particle_flags::isSet( ac.getFlags(idx1), data::particle_flags::GHOST) &&
       !data::particle_flags::isSet( ac.getFlags(idx2), data::particle_flags::GHOST) )
   {
      // local-local, local-global, global-global contacts

      if( data::particle_flags::isSet( ac.getFlags(idx1), data::particle_flags::GLOBAL) &&
          data::particle_flags::isSet( ac.getFlags(idx2), data::particle_flags::GLOBAL) )
      {
         // Resolve global-global contacts only on root process
         if( myRank_ != 0 )
         {
            WALBERLA_LOG_DETAIL( "Rejecting global-global contact " << contactPoint << " on non-root process." );
            return false;
         }
      } else
      {
         // Always resolve local-local and local-global contacts even if they are outside of our domain
         return true;
      }
   } else
   {
      // local-remote, global-remote or remote-remote contact

      if( data::particle_flags::isSet( ac.getFlags(idx1), data::particle_flags::GLOBAL) ||
          data::particle_flags::isSet( ac.getFlags(idx2), data::particle_flags::GLOBAL) )
      {
         // Never resolve remote-global contacts
         WALBERLA_LOG_DETAIL( "Rejecting global-remote contact " << contactPoint << "." );
         return false;
      } else if( data::particle_flags::isSet( ac.getFlags(idx1), data::particle_flags::GHOST) &&
                 data::particle_flags::isSet( ac.getFlags(idx2), data::particle_flags::GHOST) &&
                 ac.getOwner(idx1) == ac.getOwner(idx2) )
      {
         WALBERLA_LOG_DETAIL( "Rejecting remote-remote contact since it will be a local-local contact at the owner process: " << contactPoint << "." );
         return false;
      } else if( useOwnerRankTiebreak )
      {
         // See useOwnerRankTiebreak's own doc comment above for the full rationale.
         const uint_t owner1 = ac.getOwner(idx1);
         const uint_t owner2 = ac.getOwner(idx2);
         if( myRank_ != std::min( owner1, owner2 ) )
         {
            WALBERLA_LOG_DETAIL( "Rejecting contact " << contactPoint << " -- not the smallest-rank owner "
               "(owner1=" << owner1 << " owner2=" << owner2 << " myRank=" << myRank_ << ")." );
            return false;
         }
      } else
      {
         if( !domain.isContainedInProcessSubdomain(myRank_, contactPoint ) )
         {
            if( data::particle_flags::isSet( ac.getFlags(idx1), data::particle_flags::GHOST) &&
                data::particle_flags::isSet( ac.getFlags(idx2), data::particle_flags::GHOST) )
            {
               WALBERLA_LOG_DETAIL( "Rejecting remote-remote contact " << contactPoint << " since we don't own it." );
               return false;
            } else
            {
               WALBERLA_LOG_DETAIL( "Rejecting remote-local contact " << contactPoint << " since we don't own it." );
               return false;
            }
         }
      }
   }
   return true;
}

} //namespace mpi
} //namespace mesa_pd
} //namespace walberla
