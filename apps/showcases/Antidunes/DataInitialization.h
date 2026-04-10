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
//! \file   Utility.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once


#include "core/mpi/all.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"


namespace walberla
{
namespace antidunes
{

// data handling for loading a field of type ScalarField_T from file
template< typename ScalarField_T >
class ScalarFieldHandling : public field::BlockDataHandling< ScalarField_T >
{
 public:
   ScalarFieldHandling(const weak_ptr< StructuredBlockStorage >& blocks, uint_t numberGhostLayer)
      : blocks_(blocks), numberGhostLayer_(numberGhostLayer)
   {}

 protected:
   ScalarField_T* allocate(IBlock* const block) override { return allocateDispatch(block); }

   ScalarField_T* reallocate(IBlock* const block) override { return allocateDispatch(block); }

 private:
   weak_ptr< StructuredBlockStorage > blocks_;
   uint_t numberGhostLayer_;

   ScalarField_T* allocateDispatch(IBlock* const block)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks);

      return new ScalarField_T(blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block),
                               blocks->getNumberOfZCells(*block), numberGhostLayer_, real_c(0), field::fzyx);
   }
}; // class ScalarFieldHandling

// data handling for loading a field of type VectorFieldFlattened_T from file
template< typename VectorFieldFlattened_T >
class VectorFieldFlattenedHandling : public field::BlockDataHandling< VectorFieldFlattened_T >
{
 public:
   VectorFieldFlattenedHandling(const weak_ptr< StructuredBlockStorage >& blocks, uint_t numberGhostLayer)
      : blocks_(blocks), numberGhostLayer_(numberGhostLayer)
   {}

 protected:
   VectorFieldFlattened_T* allocate(IBlock* const block) override { return allocateDispatch(block); }

   VectorFieldFlattened_T* reallocate(IBlock* const block) override { return allocateDispatch(block); }

 private:
   weak_ptr< StructuredBlockStorage > blocks_;
   uint_t numberGhostLayer_;

   VectorFieldFlattened_T* allocateDispatch(IBlock* const block)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks);

      return new VectorFieldFlattened_T(blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block),
                                        blocks->getNumberOfZCells(*block), numberGhostLayer_, field::fzyx);
   }
}; // class VectorFieldFlattenedHandling


class HydrostaticDensityFtor {
public:
   HydrostaticDensityFtor(const std::shared_ptr< Vector3 < real_t > >& acceleration, const real_t absoluteLiquidHeight) :
   acceleration_(acceleration), absoluteLiquidHeight_(absoluteLiquidHeight) {}

   real_t operator()(const Vector3<real_t>& pos){
      return real_c(1) + real_c(3) * (*acceleration_)[2] * (pos[2] - absoluteLiquidHeight_);
   }
   
private:
   const std::shared_ptr< Vector3 < real_t > > acceleration_;
   const real_t absoluteLiquidHeight_ ;
};

// function describing the global initialization profile
inline real_t initializationProfile(real_t x, real_t amplitude, real_t offset, real_t wavelength)
{
   return amplitude * std::cos(x / wavelength * real_c(2) * math::pi + math::pi) + offset;
}


void initFillFieldSinus(const std::shared_ptr<const SimulationParams>& simParams, const std::shared_ptr<SimulationObjects>& simObjs)
{
   // samples used in the Monte-Carlo like estimation of the fill level
   const uint_t fillLevelInitSamples = uint_c(101);
   const uint_t numTotalPoints = fillLevelInitSamples * fillLevelInitSamples;
   const real_t stepsize       = real_c(1) / real_c(fillLevelInitSamples - 1);

   auto& bf             = simObjs->blockForest;
   // initialize free-surface sine profile
   for (auto blockIt = bf->begin(); blockIt != bf->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(simObjs->fillFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         // cell in block-local coordinates
         // TODO: check if it's cell center. If yes, sampling is wrong cuz it's shifted by 0.5
         const Cell localCell = fillFieldIt.cell();

         // get cell in global coordinates
         Cell globalCell = localCell;
         bf->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

         // Monte-Carlo like estimation of the fill level:
         // create uniformly-distributed sample points in each cell and count the number of points below the sine
         // profile; this fraction of points is used as the fill level to initialize the profile
         uint_t numPointsBelow = uint_c(0);

         for (uint_t xSample = uint_c(0); xSample < fillLevelInitSamples; ++xSample)
         {
            // Pascal et al. (2021) defined the amplitude to span from minimum peak to maximum peak; in
            // initializationProfile(), the amplitude is defined to range from the average to the maximum peak
            const real_t functionValue = initializationProfile(real_c(globalCell[0]) + real_c(xSample) * stepsize, 
                  simParams->computed.sinusAmplitude, simParams->computed.absoluteLiquidHeight + real_c(0.5), simParams->computed.sinusWavelength);

            for (uint_t zSample = uint_c(0); zSample < fillLevelInitSamples; ++zSample)
            {
               const real_t zPoint = real_c(globalCell[2]) + real_c(zSample) * stepsize;
               // with operator <, a fill level of 1 can not be reached when the line is equal to the cell's top
               // border; with operator <=, a fill level of 0 can not be reached when the line is equal to the cell's
               // bottom border
               if (zPoint < functionValue) { ++numPointsBelow; }
            }
         }

         // fill level is fraction of points below sine profile
         fillField->get(localCell) = real_c(numPointsBelow) / real_c(numTotalPoints);
      }) // WALBERLA_FOR_ALL_CELLS
   }
}

void initPdfFieldPoiseuille(const std::shared_ptr<const SimulationParams>& simParams, const std::shared_ptr<SimulationObjects>& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing Poiseuille velocity profile");
   auto bf = simObjs->blockForest;
   const real_t rho = real_c(1);

   for (auto blockIt = bf->begin(); blockIt != bf->end(); ++blockIt)
   {
      PdfField_T* const pdfField           = blockIt->getData< PdfField_T >(simObjs->pdfFieldID);
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(simObjs->fillFieldID);

      WALBERLA_FOR_ALL_CELLS_XYZ(
         pdfField, const Vector3< real_t > coord = bf->getBlockLocalCellCenter(*blockIt, Cell(x, y, z));

         Vector3< real_t > velocity(real_c(0));

         auto localBedHeight = initializationProfile(
            coord[0], simParams->computed.sinusAmplitude, simParams->computed.floorHeight, simParams->computed.sinusWavelength);
         auto heightAboveBed = coord[2] - localBedHeight;

         const real_t fillLevel = fillField->get(x, y, z);

         if (heightAboveBed >= real_c(0) && fillLevel > real_c(0)) {
            velocity[0] = (*simObjs->acceleration)[0] / (real_c(2) * simParams->computed.viscosity) * heightAboveBed *
                          (real_c(2) * (simParams->computed.liquidHeight + real_c(0.5)) - heightAboveBed);
         } 

         pdfField->setToEquilibrium(x, y, z, velocity, rho);
      ) // WALBERLA_FOR_ALL_CELLS_XYZ
   }
}

/***********************************************************************************************************************
 * Initialize the hydrostatic pressure in the direction in which a force is acting in ALL cells (regardless of a cell's
 * flag). The velocity remains unchanged.
 *
 * The force vector must have only one component, i.e., the direction of the force can only be in x-, y- or z-axis.
 * The variable fluidHeight determines the height at which the density is equal to reference density (=1).
 **********************************************************************************************************************/
void initPdfFieldHydrostaticPressure(const std::shared_ptr<const SimulationParams>& simParams, const std::shared_ptr<SimulationObjects>& simObjs)
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing hydrostatic pressure");

   auto bf = simObjs->blockForest;
   HydrostaticDensityFtor hydrostaticDensityFtor(simObjs->acceleration, simParams->computed.absoluteLiquidHeight);

   for (auto blockIt = bf->begin(); blockIt != bf->end(); ++blockIt)
   {
      PdfField_T* const pdfField = blockIt->getData< PdfField_T >(simObjs->pdfFieldID);

      CellInterval local = pdfField->xyzSizeWithGhostLayer(); // block-, i.e., process-local cell interval

      for (const auto& cell : local)
      {
         // initialize the (hydrostatic) pressure, i.e., LBM density
         // Bernoulli: p = p0 + density * gravity * height
         // => LBM (density=1): rho = rho0 + gravity * height = 1 + 1/cs^2 * g * h = 1 + 3 * g * h
         // shift global cell by 0.5 since density is set for cell center

         Vector3< real_t > cellCenter = bf->getBlockLocalCellCenter(*blockIt, cell);
         const real_t rho             = hydrostaticDensityFtor(cellCenter);

         const Vector3< real_t > velocity = pdfField->getVelocity(cell);

         pdfField->setDensityAndVelocity(cell, velocity, rho);
      }
   }
}

void initSpheresFromFile(const std::string& filename, walberla::mesa_pd::data::ParticleStorage& ps,
                         const walberla::mesa_pd::domain::IDomain& domain, walberla::real_t density,
                         const Vector3< uint_t >& domainSize,
                         const std::function< bool(walberla::Vector3< real_t >) >& particleCreateFunction,
                         math::AABB simulationDomain, uint_t bedCopiesInX, uint_t bedCopiesInY, uint_t& numParticles,
                         real_t& maxParticleHeight, const real_t& scalingFactor)
{
   using namespace walberla::mesa_pd;
   using namespace walberla::mesa_pd::data;

   auto rank = walberla::mpi::MPIManager::instance()->rank();

   std::string textFile;

   WALBERLA_ROOT_SECTION()
   {
      std::ifstream t(filename.c_str());
      if (!t) { WALBERLA_ABORT("Invalid input file " << filename << "\n"); }
      std::stringstream buffer;
      buffer << t.rdbuf();
      textFile = buffer.str();
   }

   walberla::mpi::broadcastObject(textFile);

   std::istringstream fileIss(textFile);
   std::string line;

   // first line contains generation domain sizes
   std::getline(fileIss, line);
   Vector3< real_t > generationDomainSize_SI(0_r);
   std::istringstream firstLine(line);
   firstLine >> generationDomainSize_SI[0] >> generationDomainSize_SI[1] >> generationDomainSize_SI[2];
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(generationDomainSize_SI)

   WALBERLA_CHECK_EQUAL(uint_c(scalingFactor * generationDomainSize_SI[0] * real_c(bedCopiesInX)), domainSize[0],
                        "Error: Generated bed with copies and simulation domain do not match in x!")
   WALBERLA_CHECK_EQUAL(uint_c(scalingFactor * generationDomainSize_SI[1] * real_c(bedCopiesInY)), domainSize[1],
                        "Error: Generated bed with copies and simulation domain do not match in y!")

   numParticles      = 0;
   maxParticleHeight = 0_r;

   while (std::getline(fileIss, line))
   {
      std::istringstream iss(line);

      data::ParticleStorage::uid_type uID;
      data::ParticleStorage::position_type pos;
      walberla::real_t radius;
      iss >> uID >> pos[0] >> pos[1] >> pos[2] >> radius;
      radius *= scalingFactor;

      for (uint_t copyInYDir = 0; copyInYDir < bedCopiesInY; ++copyInYDir)
      {
         for (uint_t copyInXDir = 0; copyInXDir < bedCopiesInX; ++copyInXDir)
         {
            auto particlePos = pos;

            particlePos[0] += real_c(copyInXDir) * generationDomainSize_SI[0];
            particlePos[1] += real_c(copyInYDir) * generationDomainSize_SI[1];

            particlePos *= scalingFactor;

            maxParticleHeight = std::max(maxParticleHeight, particlePos[2] + radius);

            if (!particleCreateFunction(particlePos)) continue;

            WALBERLA_CHECK(simulationDomain.contains(particlePos),
                           "Particle read from file is not contained in simulation domain");

            if (!domain.isContainedInProcessSubdomain(uint_c(rank), particlePos)) continue;

            auto pIt = ps.create();
            pIt->setPosition(particlePos);
            pIt->getBaseShapeRef() = std::make_shared< data::Sphere >(radius);
            pIt->getBaseShapeRef()->updateMassAndInertia(density);
            pIt->setInteractionRadius(radius);
            pIt->setOwner(rank);
            pIt->setType(0);

            numParticles++;
         }
      }

      WALBERLA_CHECK_EQUAL(iss.tellg(), -1);
   }
   walberla::mpi::allReduceInplace(maxParticleHeight, walberla::mpi::MAX);
   walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
}



} // namespace antidunes
} // namespace walberla