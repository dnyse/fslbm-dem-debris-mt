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
//! \file CollideIsolatedObstacleIterationTest.cpp
//! Hypothesis test for a collide-stage NaN in the hopper: the "treat non-interface cells" branch of
//! StreamReconstructAdvectSweep.h never writes PDFs for an obstacle-flagged cell (only the fill level is copied, no
//! streaming or reconstruction), while the generated collision kernel runs on every cell including obstacles (the
//! CollideSweep in BuildTimeloop.h has no flag gate). PDFs of an obstacle cell that no bounce-back link refreshes are
//! therefore changed only by collide, with no streaming input and no force update (ForceDensityCodegenSweep skips
//! obstacle cells; the force is held at zero here): a repeated, isolated iteration of the collision map at
//! relaxationRate 1.99753, close to the omega = 2 stability limit. A single application on near-equilibrium input is
//! finite (CollideKernelCrashReplayTest.cpp); this test iterates the same map up to 20000 times from that input and
//! reports whether it diverges. Diagnostic only: both outcomes are logged, no assertion.
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

   // Same starting PDFs as CollideKernelCrashReplayTest.cpp (a recorded pre-collide state, an ordinary
   // near-equilibrium point, not a specially chosen unstable one).
   const real_t f[27] = { real_t(2.9828e-01), real_t(7.4444e-02), real_t(7.4698e-02), real_t(7.5976e-02),
                          real_t(7.3192e-02), real_t(7.2307e-02), real_t(7.6906e-02), real_t(1.8962e-02),
                          real_t(1.8267e-02), real_t(1.9026e-02), real_t(1.8329e-02), real_t(1.8046e-02),
                          real_t(1.8107e-02), real_t(1.8417e-02), real_t(1.7742e-02), real_t(1.9194e-02),
                          real_t(1.9259e-02), real_t(1.9589e-02), real_t(1.8871e-02), real_t(4.4280e-03),
                          real_t(4.5965e-03), real_t(4.4431e-03), real_t(4.6122e-03), real_t(4.7097e-03),
                          real_t(4.8889e-03), real_t(4.7258e-03), real_t(4.9056e-03) };
   for (uint_t q = uint_t(0); q < uint_t(27); ++q)
      pdfField->get(ix, iy, iz, q) = f[q];

   // Zero force throughout: ForceDensityCodegenSweep does not update obstacle-flagged cells (the isolated-iteration
   // condition under test).
   forceField->get(ix, iy, iz, uint_t(0)) = real_t(0);
   forceField->get(ix, iy, iz, uint_t(1)) = real_t(0);
   forceField->get(ix, iy, iz, uint_t(2)) = real_t(0);

   LatticeModel_T::Sweep sweep(pdfFieldID);

   const uint_t maxIters = uint_t(20000);
   uint_t divergedAt      = uint_t(0);
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

      if (iter <= uint_t(10) || iter % uint_t(100) == uint_t(0) || !allFinite)
      {
         WALBERLA_LOG_INFO_ON_ROOT("iter=" << iter << " rho=" << std::setprecision(12) << rho
                                            << " allFinite=" << (allFinite ? "true" : "false"));
      }

      if (!allFinite)
      {
         divergedAt = iter;
         break;
      }
   }

   if (divergedAt > uint_t(0))
   {
      WALBERLA_LOG_INFO_ON_ROOT("CONFIRMED: repeated isolated collide iteration (zero force, zero streaming "
                                 "input, omega=" << omega << ") diverges to non-finite PDFs after " << divergedAt
                                 << " iterations, starting from an ordinary near-equilibrium point.");
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Did NOT diverge within " << maxIters
                                 << " iterations -- this specific hypothesis (isolated repeated collide "
                                    "iteration as the divergence mechanism) is NOT supported by this test.");
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
