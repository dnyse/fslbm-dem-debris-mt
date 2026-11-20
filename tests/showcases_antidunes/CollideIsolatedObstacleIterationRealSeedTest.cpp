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
//! \file CollideIsolatedObstacleIterationRealSeedTest.cpp
//! Variant of CollideIsolatedObstacleIterationTest.cpp seeded with a recorded post-collide state of an
//! obstacle-flagged cell (world (267.5,18.5,19.5), step 2700, rho = 1.007356, near-equilibrium PDFs, force 0). In the
//! recorded run this cell stayed OBST while its rho drifted by ~1e-4 to ~5e-4 per step until its PDFs became
//! non-finite at step 2828. If isolated collide iteration (zero force, no other input, omega = 1.99753) is the whole
//! mechanism, iterating from this seed should reproduce (a) that per-iteration drift and (b) eventual divergence to
//! non-finite values. Runs up to 200 iterations and logs rho and drho per iteration; diagnostic only, no assertion.
//======================================================================================================================

#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h" // pulls in SimulationObjects.h -> LatticeModel_T, PdfField_T, VectorFieldFlattened_T

#include "blockforest/Initialization.h"

#include "field/AddToStorage.h"

#include "lbm/field/AddToStorage.h"

#include <core/Environment.h>
#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

#include <cmath>
#include <iomanip>

namespace walberla { namespace antidunes {

void runTest()
{
   auto blocks = blockforest::createUniformBlockGrid(uint_t(1), uint_t(1), uint_t(1), uint_t(64), uint_t(20),
                                                       uint_t(30), real_t(1), uint_t(1), uint_t(1), uint_t(1),
                                                       false, false, false);

   auto forceFieldID =
      field::addToStorage< VectorFieldFlattened_T >(blocks, "Force field", real_c(0), field::fzyx, uint_c(1));

   const real_t omega = real_t(1.99753); // relaxationRate of the scenario
   LatticeModel_T latticeModel(forceFieldID, omega);

   auto pdfFieldID =
      lbm::addPdfFieldToStorage(blocks, "PDF field", latticeModel, Vector3< real_t >(real_t(0)), real_t(1.0), field::fzyx);

   IBlock* block = &(*blocks->begin());
   auto* pdfField   = block->getData< PdfField_T >(pdfFieldID);
   auto* forceField = block->getData< VectorFieldFlattened_T >(forceFieldID);

   const cell_idx_t ix = cell_idx_t(11);
   const cell_idx_t iy = cell_idx_t(18);
   const cell_idx_t iz = cell_idx_t(19);

   // Recorded state [WATCH t=2700 HOPPER-after-collide]: world=(267.5,18.5,19.5) local=(11,18,19) flag=OBST
   // rho=1.007356e+00 fill=1.0 force=(0,0,0).
   const real_t f[27] = { real_t(2.9844e-01), real_t(7.4683e-02), real_t(7.4412e-02), real_t(7.6309e-02),
                          real_t(7.2937e-02), real_t(7.3678e-02), real_t(7.5695e-02), real_t(1.9029e-02),
                          real_t(1.8316e-02), real_t(1.9093e-02), real_t(1.8122e-02), real_t(1.8451e-02),
                          real_t(1.8357e-02), real_t(1.8883e-02), real_t(1.7964e-02), real_t(1.8928e-02),
                          real_t(1.8888e-02), real_t(1.9309e-02), real_t(1.8543e-02), real_t(4.5143e-03),
                          real_t(4.7125e-03), real_t(4.4601e-03), real_t(4.7212e-03), real_t(4.6530e-03),
                          real_t(4.8115e-03), real_t(4.6109e-03), real_t(4.8347e-03) };
   for (uint_t q = uint_t(0); q < uint_t(27); ++q)
      pdfField->get(ix, iy, iz, q) = f[q];

   // Force was exactly zero for this cell throughout its recorded OBST phase.
   forceField->get(ix, iy, iz, uint_t(0)) = real_t(0);
   forceField->get(ix, iy, iz, uint_t(1)) = real_t(0);
   forceField->get(ix, iy, iz, uint_t(2)) = real_t(0);

   LatticeModel_T::Sweep sweep(pdfFieldID);

   const uint_t maxIters = uint_t(200);
   uint_t divergedAt      = uint_t(0);
   real_t prevRho         = real_t(0);
   for (uint_t q = uint_t(0); q < uint_t(27); ++q) prevRho += f[q];

   for (uint_t iter = uint_t(1); iter <= maxIters; ++iter)
   {
      sweep.collide(block, uint_t(0));

      real_t rho = real_t(0);
      bool allFinite = true;
      for (uint_t q = uint_t(0); q < uint_t(27); ++q)
      {
         const real_t v = pdfField->get(ix, iy, iz, q);
         rho += v;
         if (!std::isfinite(v)) allFinite = false;
      }

      WALBERLA_LOG_INFO_ON_ROOT("iter=" << iter << " rho=" << std::setprecision(12) << rho
                                         << " drho=" << (rho - prevRho)
                                         << " allFinite=" << (allFinite ? "true" : "false"));
      prevRho = rho;

      if (!allFinite)
      {
         divergedAt = iter;
         break;
      }
   }

   if (divergedAt > uint_t(0))
   {
      WALBERLA_LOG_INFO_ON_ROOT("CONFIRMED: repeated isolated collide iteration from this run's REAL "
                                 "OBST-phase seed (zero force, zero streaming input, omega=" << omega
                                 << ") diverges to non-finite PDFs after " << divergedAt << " iterations.");
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Did NOT diverge within " << maxIters << " iterations from the real seed either -- "
                                 "check the printed per-iteration drho trend above: if it is not monotonic/decaying "
                                 "toward zero, isolated collide iteration alone still does not explain the real "
                                 "multi-step drift, and something else must be perturbing this cell every step.");
   }
}

}} // namespace walberla::antidunes

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   walberla::debug::enterTestMode();
   walberla::antidunes::runTest();
   return EXIT_SUCCESS;
}
