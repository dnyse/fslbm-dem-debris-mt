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
//! \file   Evaluators.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "Utility.h"

namespace walberla
{
namespace antidunes
{

// evaluates X-Z slices of solid volume fraction, fill level and X-velocity
template< typename AntidunesBoundaryHandling_T, typename PdfField_T, typename FlagField_T, typename ScalarField_T >
class AverageDataSliceEvaluator
{
 public:
   AverageDataSliceEvaluator(const shared_ptr< StructuredBlockStorage >& blocks, const ConstBlockDataID& flagFieldID,
                             const ConstBlockDataID& fillFieldID, const ConstBlockDataID& pdfFieldID, const std::string& baseFolderName, uint_t frequency)
      : blocks_(blocks), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID), pdfFieldID_(pdfFieldID), 
      baseFolderName_(baseFolderName), frequency_(frequency), executionCounter_(uint_c(0))
   {
      xlen_ = blocks_->getNumberOfXCells();
      ylen_ = blocks_->getNumberOfYCells();
      zlen_ = blocks_->getNumberOfZCells();

      x_z_SolidVolumeFraction_ = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_FillLevel_           = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_VelocityX_           = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_FluidCellCount_      = std::vector< uint_t >(xlen_ * zlen_, 0);
      maxFluidZPos_            = uint_t(0);
   }

   void operator()(){
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {  
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   void compute()
   {
      // erase data
      std::fill(x_z_SolidVolumeFraction_.begin(), x_z_SolidVolumeFraction_.end(), real_t(0));
      std::fill(x_z_FillLevel_.begin(), x_z_FillLevel_.end(), real_t(0));
      std::fill(x_z_VelocityX_.begin(), x_z_VelocityX_.end(), real_t(0));
      std::fill(x_z_FluidCellCount_.begin(), x_z_FluidCellCount_.end(), 0);

      // fill contributions
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const PdfField_T* const pdfField     = block->getData< const PdfField_T >(pdfFieldID_);
         const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);
         const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);

         const auto solidMO     = flagField->getFlag(AntidunesBoundaryHandling_T::movingObstacleFlagID);
         const auto solidNoSlip = flagField->getFlag(AntidunesBoundaryHandling_T::noSlipFlagID);

         CellInterval xyz = flagField->xyzSize();
         Cell globalCell;

         maxFluidZPos_ = uint_t(0);

         // iterate all (inner) cells in the field
         for (auto cell : xyz)
         {
            blocks_->transformBlockLocalToGlobalCell(globalCell, *block, cell);
            auto entryIdx = uint_c(globalCell.x()) + uint_c(globalCell.z()) * xlen_;
            if (flagField->isFlagSet(cell, solidMO) || flagField->isFlagSet(cell, solidNoSlip))
            {
               x_z_SolidVolumeFraction_[entryIdx] += real_t(1);
               x_z_FillLevel_[entryIdx] += real_t(1);
            }
            else
            {
               auto fillLevel = fillField->get(cell);
               x_z_FillLevel_[entryIdx] += fillLevel;
               if (fillLevel > 0_r)
               {
                  x_z_VelocityX_[entryIdx] += pdfField->getVelocity(cell)[0];
                  ++x_z_FluidCellCount_[entryIdx];

                  maxFluidZPos_ = std::max(uint_t(globalCell.z()), maxFluidZPos_);
               }
            }
         }
      }

      // reduce this information to the root process
      mpi::reduceInplace(x_z_SolidVolumeFraction_, mpi::SUM);
      mpi::reduceInplace(x_z_FillLevel_, mpi::SUM);
      mpi::reduceInplace(x_z_VelocityX_, mpi::SUM);
      mpi::reduceInplace(x_z_FluidCellCount_, mpi::SUM);
      mpi::reduceInplace(maxFluidZPos_, mpi::MAX);

      WALBERLA_ROOT_SECTION()
      {
         // normalize
         for (uint_t i = 0; i < x_z_VelocityX_.size(); ++i)
         {
            if (x_z_FluidCellCount_[i] > 0) x_z_VelocityX_[i] /= real_c(x_z_FluidCellCount_[i]);
         }
         real_t invNumYCells = 1_r / real_c(ylen_);
         std::for_each(x_z_SolidVolumeFraction_.begin(), x_z_SolidVolumeFraction_.end(),
                     [invNumYCells](real_t& n) { n *= invNumYCells; });
         std::for_each(x_z_FillLevel_.begin(), x_z_FillLevel_.end(), [invNumYCells](real_t& n) { n *= invNumYCells; });
      }
      // note: only root process has the complete information!
   }

   void writeToFile(){
      WALBERLA_ROOT_SECTION()
      {
         const std::string base = baseFolderName_ + "/avgDataSlices";
         std::filesystem::create_directories(base);
         write2DVectorToFile(x_z_SolidVolumeFraction_, xlen_, zlen_,
                              base + "/svfSlice_" + std::to_string(executionCounter_) + ".txt");
         write2DVectorToFile(x_z_FillLevel_, xlen_, zlen_,
                              base + "/fillSlice_" + std::to_string(executionCounter_) + ".txt");
         write2DVectorToFile(x_z_VelocityX_, xlen_, zlen_,
                              base + "/velXSlice_" + std::to_string(executionCounter_) + ".txt");
         // WALBERLA_LOG_DEVEL(" Max fluid z-position = " << maxFluidZPos_);
      }
   }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const ConstBlockDataID flagFieldID_;
   const ConstBlockDataID fillFieldID_;
   const ConstBlockDataID pdfFieldID_;

   uint_t xlen_;
   uint_t ylen_;
   uint_t zlen_;
   uint_t maxFluidZPos_;
   std::vector< real_t > x_z_SolidVolumeFraction_;
   std::vector< real_t > x_z_FillLevel_;
   std::vector< real_t > x_z_VelocityX_;
   std::vector< uint_t > x_z_FluidCellCount_;

   std::string baseFolderName_;
   uint_t frequency_;
   uint_t executionCounter_;
};


template< typename ParticleAccessor_T >
class BedloadTransportEvaluator
{
 public:
   BedloadTransportEvaluator(const shared_ptr< ParticleAccessor_T >& ac, real_t normalizationFactor,
                             uint_t numParticles, const std::string& baseFolderName, uint_t frequency)
      : ac_(ac), normalizationFactor_(normalizationFactor), numParticles_(numParticles),
        baseFolderName_(baseFolderName), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {  
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   void compute(){
      real_t transportRate(real_t(0));
      real_t velocity(real_t(0));

      for (uint_t i = uint_t(0); i < ac_->size(); ++i)
      {
         if (!isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GHOST) &&
             !isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GLOBAL))
         {
            auto velX = ac_->getLinearVelocity(i)[0];
            transportRate += velX * ac_->getVolume(i);
            velocity += velX;
         }
      }

      // only reduce to root
      mpi::reduceInplace(transportRate, mpi::SUM);
      mpi::reduceInplace(velocity, mpi::SUM);

      avgTransportRate_ = transportRate * normalizationFactor_;
      averageVelocity_  = velocity / real_c(numParticles_);
      // note: only root process has the complete information!
   }

   void writeToFile(){
      WALBERLA_ROOT_SECTION()
      {
         // append values to the same file
         std::ofstream bedloadFile(baseFolderName_ + "/bedload.txt", std::ofstream::app);
         bedloadFile << executionCounter_ << " " << avgTransportRate_ << " " << averageVelocity_ << std::endl;
         bedloadFile.close();
      }
   }

   // sum_i V_p,i * u_p,i / (L*W)
   real_t getTransportRate() const { return avgTransportRate_; }
   real_t getAverageVelocity() const { return averageVelocity_; }

 private:
   shared_ptr< ParticleAccessor_T > ac_;
   real_t normalizationFactor_;
   uint_t numParticles_;
   real_t averageVelocity_;
   real_t avgTransportRate_;

   std::string baseFolderName_;
   uint_t frequency_;
   uint_t executionCounter_;
};


template< typename ParticleAccessor_T >
class TotalHydrodynamicForceEvaluator
{
 public:
   TotalHydrodynamicForceEvaluator(const shared_ptr< ParticleAccessor_T >& ac,
       const std::string& baseFolderName, uint_t frequency)
      : ac_(ac), baseFolderName_(baseFolderName),
      frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0))
      {  
         compute();
         writeToFile();
      }

      ++executionCounter_;
   }

   void compute(){
      totalHydrodynamicForce_ = Vector3< real_t >(real_c(0));
      maxParticlePosZ_ = real_c(0);

      for (uint_t i = uint_t(0); i < ac_->size(); ++i)
      {
         if (!isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GHOST) &&
            !isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GLOBAL))
         {
            totalHydrodynamicForce_ += ac_->getHydrodynamicForce(i);
            maxParticlePosZ_ = std::max(ac_->getPosition(i)[2], maxParticlePosZ_);
         }
      }

      // only reduce to root
      mpi::reduceInplace(totalHydrodynamicForce_, mpi::SUM);
      mpi::reduceInplace(maxParticlePosZ_, mpi::MAX);
   }

   void writeToFile(){
      WALBERLA_ROOT_SECTION()
      {
         // append values to the same file
         std::ofstream bedloadFile(baseFolderName_ + "/totalParticleHydroAndMaxPosZ.txt", std::ofstream::app);
         bedloadFile << executionCounter_ << " " << totalHydrodynamicForce_ << " " << maxParticlePosZ_ << std::endl;
         bedloadFile.close();
      }
   }

   Vector3< real_t > getTotalHydrodynamicForce() const { return totalHydrodynamicForce_; }

 private:
   shared_ptr< ParticleAccessor_T > ac_;
   Vector3< real_t > totalHydrodynamicForce_;
   real_t maxParticlePosZ_;

   std::string baseFolderName_;
   uint_t frequency_;
   uint_t executionCounter_;
};


} // namespace antidunes
} // namespace walberla
