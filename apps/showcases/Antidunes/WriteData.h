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
//! \file WriteData.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "SimulationParams.h" 
#include "SimulationObjects.h" 

namespace walberla {
namespace antidunes {

void writeCheckpoint(const std::shared_ptr<const SimulationParams>& simParams, const std::shared_ptr<SimulationObjects>& simObjs, uint_t timestep)
{
   // Write to tmp_ files first, then atomically rename. This prevents corrupt snapshots if the job is killed mid-write.
   const std::string& baseTmp = simParams->blockForest.snapshotBaseFolder + "/tmp_";
   const std::string& baseActual = simParams->blockForest.snapshotBaseFolder + "/";

   // save fields ----------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT("Storing snapshots in " << baseTmp << "*");
   simObjs->blockForest->saveBlockData(baseTmp + simObjs->pdfFieldFile,        simObjs->pdfFieldID);
   simObjs->blockForest->saveBlockData(baseTmp + simObjs->fillFieldFile,       simObjs->fillFieldID);
   simObjs->blockForest->saveBlockData(baseTmp + simObjs->excessMassFieldFile, simObjs->excessMassFieldID);
   simObjs->blockForest->saveBlockData(baseTmp + simObjs->particleStorageFile, simObjs->particleStorageID);
   
   // save simState ----------------------------------------
   WALBERLA_ROOT_SECTION()
   {
      std::string tmpSimStateFile = baseTmp + simObjs->simStateFile;
      std::ofstream file;
      file.open(tmpSimStateFile.c_str());
      file << std::setprecision(16);
      file << timestep + 1 << "\n";                // beginTimestep
      file << (*simObjs->acceleration)[0] << "\n"; // x-acceleration
      file.close();
   }

   // save pidState ----------------------------------------
   simObjs->pidController->writeStateToFile(baseTmp + simObjs->pidStateFile);

   WALBERLA_MPI_BARRIER();

   WALBERLA_LOG_INFO_ON_ROOT("Stored snapshots in " << baseTmp << "*");

   // rename files ----------------------------------------
   WALBERLA_ROOT_SECTION() {
      for (const std::string& f : { simObjs->pdfFieldFile, 
                                    simObjs->fillFieldFile,
                                    simObjs->excessMassFieldFile, 
                                    simObjs->particleStorageFile,
                                    simObjs->simStateFile,
                                    simObjs->pidStateFile
                                 })
      {
         renameFile(baseTmp + f, baseActual + f);
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT(baseTmp << "_* files were renamed to " << baseActual << "* (old snapshots are replaced). Checkpointing is completed.");
}   

void writeFluidInfo(const std::shared_ptr<const SimulationParams>& simParams, const std::shared_ptr<SimulationObjects>& simObjs, uint_t timestep)
{
   WALBERLA_ROOT_SECTION()
   {
      const real_t meanVelocityX = (*simObjs->meanVelocity)[0];      
      const real_t froudeNumber = meanVelocityX / real_c(std::sqrt(simParams->computed.liquidHeight * std::abs((*simObjs->acceleration)[2])));
      const real_t reynoldsNumber = meanVelocityX * simParams->computed.liquidHeight / simParams->computed.viscosity;
      const real_t weberNumber = real_c(1.0) * meanVelocityX * meanVelocityX * simParams->computed.liquidHeight / simParams->computed.surfaceTension;

      std::ofstream fluidInfoFile(simParams->evaluation.baseFolderName + "/fluidInfo.txt", std::ofstream::app);
      fluidInfoFile << timestep << " " << (*simObjs->acceleration)[0] << " " << meanVelocityX << " " << (*simObjs->totalFluidMass) 
                        << " " << froudeNumber << " " << reynoldsNumber << " " << weberNumber << "\n";
      fluidInfoFile.close();
      
      if (std::isnan(reynoldsNumber)) WALBERLA_ABORT("reynoldsNumber is nan!")
   }
}

} // namespace antidunes
} // namespace walberla
