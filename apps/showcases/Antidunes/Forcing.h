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
//! \file   Forcing.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "PIDController.h"
#include "Utility.h"

#include <deque>
#include <utility>

namespace walberla
{
namespace antidunes
{
// sweep for computing the force density from the fluid density and fill level
template< typename PdfField_T, typename FlagField_T, typename VectorFieldFlattened_T, typename ScalarField_T >
class ForceDensityCodegenSweep
{
 public:
   ForceDensityCodegenSweep(BlockDataID forceDensityFieldID, ConstBlockDataID pdfFieldID, ConstBlockDataID flagFieldID,
                            ConstBlockDataID fillFieldID,
                            const walberla::free_surface::FlagInfo< FlagField_T >& flagInfo,
                            std::shared_ptr< Vector3< real_t > > globalAcceleration)
      : forceDensityFieldID_(forceDensityFieldID), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID),
        fillFieldID_(fillFieldID), flagInfo_(flagInfo), globalAcceleration_(globalAcceleration)
   {}

   void operator()(IBlock* const block)
   {
      VectorFieldFlattened_T* const forceDensityField = block->getData< VectorFieldFlattened_T >(forceDensityFieldID_);
      const PdfField_T* const pdfField                = block->getData< const PdfField_T >(pdfFieldID_);
      const FlagField_T* const flagField              = block->getData< const FlagField_T >(flagFieldID_);
      const ScalarField_T* const fillField            = block->getData< const ScalarField_T >(fillFieldID_);

      WALBERLA_FOR_ALL_CELLS(forceDensityFieldIt, forceDensityField, pdfFieldIt, pdfField, flagFieldIt, flagField,
                             fillFieldIt, fillField, {
                                flag_t flag = *flagFieldIt;

                                // set force density in cells to acceleration * density * fillLevel (see equation 15
                                // in Koerner et al., 2005);
                                if (flagInfo_.isInterface(flag))
                                {
                                   const real_t density   = pdfField->getDensity(pdfFieldIt.cell());
                                   forceDensityFieldIt[0] = (*globalAcceleration_)[0] * *fillFieldIt * density;
                                   forceDensityFieldIt[1] = (*globalAcceleration_)[1] * *fillFieldIt * density;
                                   forceDensityFieldIt[2] = (*globalAcceleration_)[2] * *fillFieldIt * density;
                                }
                                else
                                {
                                   if (flagInfo_.isLiquid(flag))
                                   {
                                      const real_t density   = pdfField->getDensity(pdfFieldIt.cell());
                                      forceDensityFieldIt[0] = (*globalAcceleration_)[0] * density;
                                      forceDensityFieldIt[1] = (*globalAcceleration_)[1] * density;
                                      forceDensityFieldIt[2] = (*globalAcceleration_)[2] * density;
                                   }
                                }
                             }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   using flag_t = typename FlagField_T::flag_t;

   BlockDataID forceDensityFieldID_;
   ConstBlockDataID pdfFieldID_;
   ConstBlockDataID flagFieldID_;
   ConstBlockDataID fillFieldID_;
   walberla::free_surface::FlagInfo< FlagField_T > flagInfo_;
   std::shared_ptr< Vector3< real_t > > globalAcceleration_;
};


template< typename FreeSurfaceBoundaryHandling_T, typename PdfField_T, typename FlagField_T >
class MeanVelocityComputer
{
 public:
   MeanVelocityComputer(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                        const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                        const ConstBlockDataID& pdfFieldID, const std::shared_ptr< Vector3< real_t > >& meanVelocity,
                        real_t averagingFactor, uint_t frequency)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), pdfFieldID_(pdfFieldID),
        meanVelocity_(meanVelocity), averagingFactor_(averagingFactor), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {
         auto blockForest = blockForest_.lock();
         WALBERLA_CHECK_NOT_NULLPTR(blockForest);

         auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
         WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

         getMeanVelocity(blockForest, freeSurfaceBoundaryHandling);
      }

      ++executionCounter_;
   }

   void getMeanVelocity(const std::shared_ptr< const StructuredBlockForest >& blockForest,
                        const std::shared_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling)
   {
      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      // use separate variables for the velocity in each direction; syntax meanVelocity[0] does not work in OMP-macro
      real_t meanVelocityX = real_c(0);
      real_t meanVelocityY = real_c(0);
      real_t meanVelocityZ = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID);
         const PdfField_T* const pdfField   = blockIt->template getData< const PdfField_T >(pdfFieldID_);

         WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, pdfFieldIt, pdfField,
                                    omp parallel for schedule(static) reduction(+:meanVelocityX)
                                       reduction(+:meanVelocityY) reduction(+:meanVelocityZ),
                                    {
            if (flagInfo.isLiquid(flagFieldIt) || flagInfo.isInterface(flagFieldIt))
            {
               const Vector3< real_t > velocity = pdfField->getVelocity(pdfFieldIt.cell());

               meanVelocityX += velocity[0];
               meanVelocityY += velocity[1];
               meanVelocityZ += velocity[2];

               //++cellCount;
            }
                                    }) // WALBERLA_FOR_ALL_CELLS_OMP
      }

      Vector3< real_t > meanVelocity(meanVelocityX, meanVelocityY, meanVelocityZ);
      mpi::allReduceInplace< real_t >(meanVelocity, mpi::SUM);
      // mpi::allReduceInplace< uint_t >(cellCount, mpi::SUM);

      meanVelocity *= averagingFactor_;
      *meanVelocity_ = meanVelocity;
   };

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   const ConstBlockDataID pdfFieldID_;

   std::shared_ptr< Vector3< real_t > > meanVelocity_;
   real_t averagingFactor_;

   uint_t frequency_;
   uint_t executionCounter_;
};


class ForcingAdjuster
{
 public:
   // switchOffStep > 0 (PIDParameters.switchOffStep): PID control runs normally until that ABSOLUTE
   // timestep (beginTimeStep + calls so far, so it stays correct across a checkpoint resume), then
   // acceleration_[0] is frozen at the mean of the PID's own commanded values over the last
   // switchOffAveragingWindow steps and never touched again -- a "constant force" phase that starts
   // from a self-consistently spun-up flow instead of a hand-set force. 0 = feature off.
   ForcingAdjuster(const shared_ptr< Vector3<real_t> >& meanVelocity,
                   const shared_ptr< Vector3<real_t> >& acceleration,
                   const shared_ptr< PIDController >& pidController,
                   uint_t frequency,
                   bool enabled = true,
                   uint_t beginTimeStep = uint_c(0),
                   uint_t switchOffStep = uint_c(0),
                   uint_t switchOffAveragingWindow = uint_c(0))
      : meanVelocity_(meanVelocity)
      , acceleration_(acceleration)
      , pidController_(pidController)
      , frequency_(frequency)
      , executionCounter_(uint_c(0))
      , enabled_(enabled)
      , beginTimeStep_(beginTimeStep)
      , switchOffStep_(switchOffStep)
      , switchOffAveragingWindow_(switchOffAveragingWindow)
      , frozen_(false)
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // PIDParameters.enabled=false (THESIS_PLAN.md item 3, "constant forcing" mode): never touch
      // acceleration_ or the PID controller -- it stays pinned at its setupSimState-computed initial
      // value ("fixed hydraulic gradient") for the whole run. MeanVelocityComputer keeps running
      // regardless (wired up unconditionally in BuildTimeloop.h), so meanVelocity is still tracked.
      if (!enabled_ || frozen_) { ++executionCounter_; return; }

      const uint_t step = beginTimeStep_ + executionCounter_;

      // PIDParameters.switchOffStep: freeze at the mean of the PID's recent commanded values.
      // history_ is identical on every rank (each entry is the broadcast value), so every rank
      // computes the identical frozen value without another reduction. On a resume at/after the
      // switch-off step history_ is empty; acceleration_[0] then already holds the checkpointed
      // (frozen) value, so it is left untouched.
      if (switchOffStep_ > uint_c(0) && step >= switchOffStep_)
      {
         if (!history_.empty())
         {
            real_t sum = real_c(0);
            for (const auto& entry : history_) { sum += entry.second; }
            const real_t mean          = sum / real_c(history_.size());
            const real_t instantaneous = (*acceleration_)[0];
            (*acceleration_)[0]        = mean;
            WALBERLA_LOG_INFO_ON_ROOT("[PID-SWITCHOFF] step=" << step << ": PID control ended, acceleration_x frozen at "
                                      << mean << " (mean of " << history_.size() << " PID updates over the last "
                                      << switchOffAveragingWindow_ << " steps; instantaneous value was " << instantaneous << ")");
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT("[PID-SWITCHOFF] step=" << step << ": PID control ended, no PID history in this run "
                                      "(resumed at/after the switch-off step) -- acceleration_x stays at the loaded value "
                                      << (*acceleration_)[0]);
         }
         frozen_ = true;
         ++executionCounter_;
         return;
      }

      real_t newForceX = real_c(0);
      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {
         // compute new forcing value on root (since flow rate only known on root)
         WALBERLA_ROOT_SECTION()
         {
            newForceX = pidController_->update(meanVelocity_->length());
         }

         // send updated external forcing to all other processes
         mpi::broadcastObject(newForceX);

         // update the acceleration
         (*acceleration_)[0] = newForceX;

         if (switchOffStep_ > uint_c(0))
         {
            history_.emplace_back(step, newForceX);
            // keep only the last switchOffAveragingWindow_ steps (always at least the newest entry;
            // window 0 therefore degenerates to the instantaneous value)
            const uint_t window = switchOffAveragingWindow_ > uint_c(0) ? switchOffAveragingWindow_ : uint_c(1);
            while (history_.size() > uint_c(1) && history_.front().first + window <= step) { history_.pop_front(); }
         }
      }

      ++executionCounter_;
   }


 private:
   shared_ptr< Vector3< real_t > > meanVelocity_;
   shared_ptr< Vector3< real_t > > acceleration_;
   shared_ptr< PIDController > pidController_;

   uint_t frequency_;
   uint_t executionCounter_;
   bool enabled_;

   uint_t beginTimeStep_;
   uint_t switchOffStep_;
   uint_t switchOffAveragingWindow_;
   bool   frozen_;
   std::deque< std::pair< uint_t, real_t > > history_; // (absolute step, PID-commanded acceleration_x)
};


} // namespace antidunes
} // namespace walberla