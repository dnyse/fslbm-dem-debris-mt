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

#include <filesystem> 

namespace walberla {
namespace antidunes {

void writeCheckpoint(const std::shared_ptr<const SimulationParams>& simParams, const std::shared_ptr<SimulationObjects>& simObjs, uint_t timestep, bool keepUnique = false, bool alternateGeneration = false)
{
   // keepUnique=true (BlockForestParams::keepAllSnapshots): write into a dedicated step_<timestep>
   // subfolder instead of the shared rolling checkpoint, so every call is retained on disk rather
   // than overwriting the previous one. Diagnostic use only (see SimulationParams.h).
   //
   // alternateGeneration=true (BlockForestParams::keepLastTwoSnapshots): write into
   // snapshotBaseFolder/gen_0 or gen_1, alternating by checkpoint index parity, so the write that
   // just completed and the one before it both survive on disk, instead of a single rolling
   // checkpoint where only the latest write survives. Mutually exclusive with keepUnique in practice.
   const std::string& snapshotDir = keepUnique
      ? simParams->blockForest.snapshotBaseFolder + "/step_" + std::to_string(timestep)
      : alternateGeneration
         ? simParams->blockForest.snapshotBaseFolder + "/gen_" +
              std::to_string((timestep / simParams->blockForest.snapshotFrequency) % uint_t(2))
         : simParams->blockForest.snapshotBaseFolder;
   if (keepUnique || alternateGeneration)
   {
      WALBERLA_ROOT_SECTION()
      {
         std::filesystem::create_directories(snapshotDir);
         // setupBlockForest (SetupStages.h) only writes blockForestFile ("blockForest.file") on a
         // fresh start that also stores -- during a resume, that file is never (re)written, so a
         // keepUnique/alternateGeneration subfolder would otherwise lack its own copy and fail to
         // load on its own. The block decomposition never changes during a run, so a plain copy from
         // the shared base folder (where it always exists) is valid and cheap; copied on every
         // checkpoint call rather than once, for simplicity, since the source never changes.
         const std::string& baseBlockForestFile = simParams->blockForest.snapshotBaseFolder + "/" + simObjs->blockForestFile;
         const std::string& subfolderBlockForestFile = snapshotDir + "/" + simObjs->blockForestFile;
         if (std::filesystem::exists(baseBlockForestFile))
            std::filesystem::copy_file(baseBlockForestFile, subfolderBlockForestFile,
               std::filesystem::copy_options::overwrite_existing);
         else
            WALBERLA_LOG_WARNING("[writeCheckpoint] " << baseBlockForestFile << " does not exist -- "
               << subfolderBlockForestFile << " will be missing too, and that subfolder will not be "
               "loadable on its own. This should only happen if snapshotBaseFolder itself was never "
               "properly seeded with a blockForest.file before this run started.");
      }
      WALBERLA_MPI_BARRIER();
   }

   // Write to tmp_ files first, then atomically rename. This prevents corrupt snapshots if the job is killed mid-write.
   const std::string& baseTmp = snapshotDir + "/tmp_";
   const std::string& baseActual = snapshotDir + "/";

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

      simObjs->dataCols->add("timestep", int(timestep));
      simObjs->dataCols->add("acceleration_x", (*simObjs->acceleration)[0]);
      simObjs->dataCols->add("meanVelocity_x", (*simObjs->meanVelocity)[0]);
      simObjs->dataCols->add("meanVelocity_y", (*simObjs->meanVelocity)[1]);
      simObjs->dataCols->add("meanVelocity_z", (*simObjs->meanVelocity)[2]);
      simObjs->dataCols->add("totalFluidMass", *simObjs->totalFluidMass);
      simObjs->dataCols->add("froudeNumber", froudeNumber);
      simObjs->dataCols->add("reynoldsNumber", reynoldsNumber);
      simObjs->dataCols->add("weberNumber", weberNumber);
      
      if (std::isnan(reynoldsNumber)) WALBERLA_ABORT("reynoldsNumber is nan!")
   }
}

} // namespace antidunes
} // namespace walberla
