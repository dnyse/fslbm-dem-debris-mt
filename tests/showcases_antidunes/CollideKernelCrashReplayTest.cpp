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
//! \file CollideKernelCrashReplayTest.cpp
//! Replays a recorded hopper-clogging collide-stage crash through the compiled generated Cumulant collision kernel
//! (AntidunesLatticeModel::Sweep::collide), not a re-implementation of its formula.
//!
//! The crash cell (world (267.5,18.5,19.5), local block index (11,18,19) for cellsPerBlock <64,20,30>) went from an
//! ordinary near-equilibrium PDF/force state, recorded in the [WATCH ... HOPPER-after-streamReconstructAdvect] line,
//! to all 27 PDFs NaN after the next stage ([WATCH ... HOPPER-after-collide]). The recorded PDFs are close to the
//! D3Q27 equilibrium weights, the force density is small, omega is never modified per cell or block, and the
//! CollideSweep in BuildTimeloop.h uses numberOfGhostLayersToInclude=0, consistent with the field-size/stride
//! computation of Sweep::collide for the PDF and force fields; a Python port of the generated kernel's arithmetic
//! gives finite output for the same input.
//!
//! This test feeds the recorded input into the compiled collide() in a single block of the same cellsPerBlock, at the
//! same cell position relative to the block edge (Y=18, 2 cells from the Y=20 edge), so the real pointer/stride
//! address computation is exercised. Pass criterion: all 27 post-collide PDFs finite. A pass means the kernel's own
//! arithmetic and addressing are not the cause; the NaN must then come from memory corruption elsewhere, which only
//! a memory checker (ASan/Valgrind) on the full binary can catch.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition, same workaround as ClipAnalyticGeometryTest.cpp.
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
   // Same cellsPerBlock as the recorded run, single block/single process: only the relative cell position within
   // the block matters for address-computation edge effects.
   auto blocks = blockforest::createUniformBlockGrid(uint_t(1), uint_t(1), uint_t(1), uint_t(64), uint_t(20),
                                                       uint_t(30), real_t(1), uint_t(1), uint_t(1), uint_t(1),
                                                       false, false, false);

   auto forceFieldID =
      field::addToStorage< VectorFieldFlattened_T >(blocks, "Force field", real_c(0), field::fzyx, uint_c(1));

   // relaxationRate of the recorded run
   const real_t omega = real_t(1.99753);
   LatticeModel_T latticeModel(forceFieldID, omega);

   auto pdfFieldID =
      lbm::addPdfFieldToStorage(blocks, "PDF field", latticeModel, Vector3< real_t >(real_t(0)), real_t(1.0), field::fzyx);

   IBlock* block = &(*blocks->begin());
   auto* pdfField   = block->getData< PdfField_T >(pdfFieldID);
   auto* forceField = block->getData< VectorFieldFlattened_T >(forceFieldID);

   // Local block index of the crash cell (world (267.5,18.5,19.5)).
   const cell_idx_t ix = cell_idx_t(11);
   const cell_idx_t iy = cell_idx_t(18);
   const cell_idx_t iz = cell_idx_t(19);

   // PDFs from [WATCH t=2838 HOPPER-after-streamReconstructAdvect] world=(267.5,18.5,19.5): the state immediately
   // before the collide that produced NaN.
   const real_t f[27] = { real_t(2.9828e-01), real_t(7.4444e-02), real_t(7.4698e-02), real_t(7.5976e-02),
                          real_t(7.3192e-02), real_t(7.2307e-02), real_t(7.6906e-02), real_t(1.8962e-02),
                          real_t(1.8267e-02), real_t(1.9026e-02), real_t(1.8329e-02), real_t(1.8046e-02),
                          real_t(1.8107e-02), real_t(1.8417e-02), real_t(1.7742e-02), real_t(1.9194e-02),
                          real_t(1.9259e-02), real_t(1.9589e-02), real_t(1.8871e-02), real_t(4.4280e-03),
                          real_t(4.5965e-03), real_t(4.4431e-03), real_t(4.6122e-03), real_t(4.7097e-03),
                          real_t(4.8889e-03), real_t(4.7258e-03), real_t(4.9056e-03) };

   for (uint_t q = uint_t(0); q < uint_t(27); ++q)
      pdfField->get(ix, iy, iz, q) = f[q];

   // Same log line's force=(4.531218e-07,0,-2.933050e-05).
   forceField->get(ix, iy, iz, uint_t(0)) = real_t(4.531218e-07);
   forceField->get(ix, iy, iz, uint_t(1)) = real_t(0.0);
   forceField->get(ix, iy, iz, uint_t(2)) = real_t(-2.933050e-05);

   real_t rhoBefore = real_t(0);
   for (uint_t q = uint_t(0); q < uint_t(27); ++q)
      rhoBefore += pdfField->get(ix, iy, iz, q);
   WALBERLA_LOG_INFO_ON_ROOT("rho before collide = " << std::setprecision(9) << rhoBefore);
   WALBERLA_CHECK(std::isfinite(rhoBefore), "Sanity check: the injected input itself must be finite.");

   LatticeModel_T::Sweep sweep(pdfFieldID);
   sweep.collide(block, uint_t(0)); // numberOfGhostLayersToInclude=0, as in CollideSweep

   bool allFinite  = true;
   real_t rhoAfter = real_t(0);
   for (uint_t q = uint_t(0); q < uint_t(27); ++q)
   {
      const real_t v = pdfField->get(ix, iy, iz, q);
      rhoAfter += v;
      if (!std::isfinite(v))
      {
         allFinite = false;
         WALBERLA_LOG_INFO_ON_ROOT("  f[" << q << "] = " << v << "  <-- NON-FINITE");
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("rho after collide = " << std::setprecision(9) << rhoAfter
                                                     << ", allFinite=" << (allFinite ? "true" : "false"));

   WALBERLA_CHECK(allFinite,
                  "The REAL compiled generated collide kernel produced a non-finite PDF from the exact "
                  "captured pre-collide state of job 4264997's crash cell -- this would mean the kernel's "
                  "own arithmetic/address computation (not external memory corruption) is the actual root "
                  "cause, contradicting the standalone Python replay done earlier this session.");
}

}} // namespace walberla::antidunes

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   walberla::debug::enterTestMode();
   walberla::antidunes::runTest();
   return EXIT_SUCCESS;
}
