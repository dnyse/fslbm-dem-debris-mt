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
//! \file DEMValidation.cpp
//! \brief Fluid-free DEM validation dry-run: runs a battery of physically-grounded validation
//!        experiments against the production DEM code (contact detection + LinearSpringDashpot +
//!        VelocityVerlet, exactly as used by Antidunes.cpp's doDEMSubCycles), with NO fluid coupling.
//!        Existing tests only check static contact geometry; this checks dynamics against closed-form
//!        analytical solutions (Kodam et al. 2010, Housner 1963) and literature bulk-behaviour ranges.
//!        See AnalyticalSolutions.h / DEMValidationScenarios.h for the physics and experiment design.
//!
//! Usage: ./DEMValidation <config>.prm  -- the [DEMValidation] block's `testCase` key selects which
//! experiment to run. One .prm file per scenario (see prms/dem_validation_*.prm).
//======================================================================================================================

#include "DEMValidationScenarios.h"

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include <cstdlib>
#include <string>

namespace walberla
{
namespace antidunes
{

int main(int argc, char** argv)
{
   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified! Usage: DEMValidation <config>.prm");

   const config::Config::BlockHandle block = cfg->getBlock("DEMValidation");
   const std::string testCase              = block.getParameter< std::string >("testCase");

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] running testCase=" << testCase);

   bool pass = false;
   if (testCase == "wallImpactSphere") pass = dem_validation::runWallImpactSphere(block);
   else if (testCase == "wallImpactCylinder") pass = dem_validation::runWallImpactCylinder(block);
   else if (testCase == "wallImpactBox") pass = dem_validation::runWallImpactBox(block);
   else if (testCase == "boxFaceDropRestitution") pass = dem_validation::runBoxFaceDropRestitution(block);
   else if (testCase == "boxInclineSlideStraightness") pass = dem_validation::runBoxInclineSlideStraightness(block);
   else if (testCase == "gjkVsAnalyticalCylinderWall") pass = dem_validation::runGjkVsAnalyticalCylinderWallCheck(block);
   else if (testCase == "pairRestitution") pass = dem_validation::runPairRestitution(block);
   else if (testCase == "energyConservation") pass = dem_validation::runEnergyConservation(block);
   else if (testCase == "cornerNormalDiscontinuity") pass = dem_validation::runCornerNormalDiscontinuity(block);
   else if (testCase == "cornerFixComparison") pass = dem_validation::runCornerFixComparison(block);
   else if (testCase == "multiPointDeepOverlapStability") pass = dem_validation::runMultiPointDeepOverlapStability(block);
   else if (testCase == "multiPointDeepOverlapFuzz") pass = dem_validation::runMultiPointDeepOverlapFuzz(block);
   else if (testCase == "multiPointEdgeCases") pass = dem_validation::runMultiPointEdgeCases(block);
   else if (testCase == "fallOntoPileStability") pass = dem_validation::runFallOntoPileStability(block);
   else if (testCase == "contactStiffnessCalibrationSweep") pass = dem_validation::runContactStiffnessCalibrationSweep(block);
   else if (testCase == "multiPointBoxContactFuzz") pass = dem_validation::runMultiPointBoxContactFuzz(block);
   else if (testCase == "productionCollisionReplay") pass = dem_validation::runProductionCollisionReplay(block);
   else if (testCase == "packing") pass = dem_validation::runPacking(block);
   else if (testCase == "angleOfRepose") pass = dem_validation::runAngleOfRepose(block);
   else if (testCase == "packingSingleShape") pass = dem_validation::runPackingSingleShape(block);
   else if (testCase == "rotatingDrum") pass = dem_validation::runRotatingDrum(block);
   else if (testCase == "archTunnelContainment") pass = dem_validation::runArchTunnelContainment(block);
   else WALBERLA_ABORT("Unknown testCase '" << testCase
                                             << "'. Valid values: wallImpactSphere, wallImpactCylinder, "
                                                "wallImpactBox, pairRestitution, energyConservation, "
                                                "cornerNormalDiscontinuity, cornerFixComparison, "
                                                "multiPointDeepOverlapStability, multiPointDeepOverlapFuzz, "
                                                "multiPointEdgeCases, fallOntoPileStability, "
                                                "contactStiffnessCalibrationSweep, "
                                                "multiPointBoxContactFuzz, productionCollisionReplay, packing, "
                                                "angleOfRepose, packingSingleShape, rotatingDrum, "
                                                "archTunnelContainment.");

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] FINAL RESULT for " << testCase << ": " << (pass ? "PASS" : "FAIL"));
   return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv) { return walberla::antidunes::main(argc, argv); }
