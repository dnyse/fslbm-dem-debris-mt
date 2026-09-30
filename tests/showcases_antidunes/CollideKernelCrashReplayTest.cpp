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
//! Direct replay of job 4264997's hopper-clogging collide-stage crash through the REAL, compiled generated
//! Cumulant collision kernel (not a hand-ported reimplementation -- see the 2026-09-25 THESIS_PLAN.md entries
//! for that separate, already-completed Python-side check, which came out clean too).
//!
//! Context: the crash cell (world (267.5,18.5,19.5), local block index (11,18,19) for this scenario's
//! cellsPerBlock=<64,20,30>) went from a clean, ordinary near-equilibrium PDF/force state (captured verbatim
//! from the [WATCH ... HOPPER-after-streamReconstructAdvect] log line immediately preceding the crash) to
//! all-27-components NaN in the very next sweep stage, [WATCH ... HOPPER-after-collide]. Prior checks this
//! session ruled out: corrupted PDFs (raw values matched D3Q27 equilibrium weights closely), a corrupted
//! force-density field (captured value was small and ordinary), a corrupted/mutated omega (never modified
//! per-cell/per-block anywhere in this codebase), and a ghost-layer/stride mismatch at the collide Sweep call
//! site (BuildTimeloop.h's CollideSweep uses numberOfGhostLayersToInclude=0, matching AntidunesLatticeModel
//! ::Sweep::collide's own field-size/stride computation exactly for both the PDF and force fields). A
//! standalone Python port of the generated kernel's arithmetic (AntidunesLatticeModel.cpp:352-682) fed this
//! exact captured input and got fully finite output.
//!
//! This test closes the last gap in that verification: it feeds the SAME captured input into the real,
//! compiled C++ AntidunesLatticeModel::Sweep::collide() -- exercising the actual pointer/stride/pystencils
//! address computation, not a hand re-typed formula -- inside a single-block field allocated with the SAME
//! cellsPerBlock<64,20,30> and the SAME relative-to-block-edge cell position as the real crash (Y=18 sits
//! 2 cells from this block's own Y=20 edge), in case some address computation bug is specific to that exact
//! offset. If this ALSO comes out finite, the true root cause is not in the collision kernel's own address
//! arithmetic at all -- it must be genuine memory corruption from elsewhere landing in this exact field slot
//! between the log's read and collide's own read, which only a memory checker (ASan/Valgrind) on the full
//! simulation binary can catch directly.
//======================================================================================================================

// SimulationObjects.h (pulled in transitively by ParticleKernels.h) uses blockforest::SimpleCommunication
// without including its definition itself -- match ClipAnalyticGeometryTest.cpp's own workaround.
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
   // Same cellsPerBlock as hopper_clog_cyl_frontblock_pidon_2node_halfrows_pass3_precursor_gas.prm (job
   // 4264997), single block/single process -- irrelevant to the physics here, only the relative cell
   // position within the block matters for testing address-computation edge effects.
   auto blocks = blockforest::createUniformBlockGrid(uint_t(1), uint_t(1), uint_t(1), uint_t(64), uint_t(20),
                                                       uint_t(30), real_t(1), uint_t(1), uint_t(1), uint_t(1),
                                                       false, false, false);

   auto forceFieldID =
      field::addToStorage< VectorFieldFlattened_T >(blocks, "Force field", real_c(0), field::fzyx, uint_c(1));

   // relaxationRate printed at startup of job 4264997: "relaxationRate = 1.99753"
   const real_t omega = real_t(1.99753);
   LatticeModel_T latticeModel(forceFieldID, omega);

   auto pdfFieldID =
      lbm::addPdfFieldToStorage(blocks, "PDF field", latticeModel, Vector3< real_t >(real_t(0)), real_t(1.0), field::fzyx);

   IBlock* block = &(*blocks->begin());
   auto* pdfField   = block->getData< PdfField_T >(pdfFieldID);
   auto* forceField = block->getData< VectorFieldFlattened_T >(forceFieldID);

   // Local block-index of the real crash cell (world (267.5,18.5,19.5) in job 4264997's scenario).
   const cell_idx_t ix = cell_idx_t(11);
   const cell_idx_t iy = cell_idx_t(18);
   const cell_idx_t iz = cell_idx_t(19);

   // Exact PDF values from [WATCH t=2838 HOPPER-after-streamReconstructAdvect] world=(267.5,18.5,19.5)
   // (job 4264997's own log), i.e. the state immediately before collide ran and produced NaN.
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
   sweep.collide(block, uint_t(0)); // numberOfGhostLayersToInclude=0, matching CollideSweep's real default

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
