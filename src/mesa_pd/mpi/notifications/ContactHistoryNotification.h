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

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ContactHistory.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/mpi/notifications/reset.h>

#include <core/logging/Logging.h>
#include <core/mpi/BufferDataTypeExtensions.h>
#include <core/mpi/Datatype.h>
#include <core/mpi/MPIManager.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <functional>
#include <map>

namespace walberla {
namespace mesa_pd {

// Optional diagnostic hook: given the two uids of a contact whose history entry is a duplicate in
// reduce() below, return a description of how this rank detected that contact (e.g. closed form vs.
// GJK/EPA fallback); printed with the duplicate warning. nullptr (default) disables it. Set by the
// Antidunes showcase (apps/showcases/Antidunes/ParticleKernels.h, ComputeContacts constructor).
inline std::function<std::string(id_t, id_t)> contactHistoryDuplicateDiagnosticHook = nullptr;

/**
 * Transmits the contact history
 */
class ContactHistoryNotification
{
public:
   struct Parameters
   {
      id_t uid_;
      std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> contactHistory_;
   };

   inline explicit ContactHistoryNotification( const data::Particle& p ) : p_(p) {}

   const data::Particle& p_;
};

template <>
inline void reset<ContactHistoryNotification>(data::Particle& p)
{
   p.setNewContactHistory(std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>());
}

inline void reduce(data::Particle&& p, const ContactHistoryNotification::Parameters& objparam)
{
   auto& ch = p.getNewContactHistoryRef();
   for (auto& entry : objparam.contactHistory_)
   {
      auto ret = ch.insert(entry);
      if (!ret.second)
      {
         // Diagnostic output before the WALBERLA_CHECK below aborts. Two ranks submitted a new
         // contact-history entry for the same (particle, partner) pair in one reduction, i.e.
         // ContactFilter accepted the contact on both (typically near a rank or periodic boundary).
         // Both entries are printed so that near-identical copies of one contact can be told apart
         // from genuinely different ones.
         WALBERLA_LOG_WARNING(
            "[ContactHistory DUPLICATE] rank=" << walberla::mpi::MPIManager::instance()->rank()
            << " -- two ranks each submitted a NEW contact-history entry for the same pair in the "
               "same reduction cycle (likely mesa_pd::mpi::ContactFilter double-accepting a contact "
               "near a rank or periodic boundary -- check both particles' positions against the "
               "block decomposition and any periodic wrap).\n"
            << "  owner particle   uid=" << p.getUid() << " pos=" << p.getPosition() << "\n"
            << "  contact partner  uid=" << entry.first << "\n"
            << "  ALREADY PRESENT (first arrival, kept): " << ret.first->second << "\n"
            << "  INCOMING (this one, rejected):         " << entry.second);
         if (contactHistoryDuplicateDiagnosticHook)
         {
            WALBERLA_LOG_WARNING(
               "[ContactHistory DUPLICATE] detection path (this rank, this contact, this cycle): "
               << contactHistoryDuplicateDiagnosticHook(p.getUid(), entry.first));
         }
      }
      WALBERLA_CHECK(ret.second, "entry already present: " << entry.first << " / " << entry.second);
   }
}

}  // namespace mesa_pd
}  // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ContactHistoryNotification& obj )
{
   buf.addDebugMarker( "ft" );
   buf << obj.p_.getUid();
   buf << obj.p_.getNewContactHistory();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ContactHistoryNotification::Parameters& objparam )
{
   buf.readDebugMarker( "ft" );
   buf >> objparam.uid_;
   buf >> objparam.contactHistory_;
   return buf;
}

} // mpi
} // walberla