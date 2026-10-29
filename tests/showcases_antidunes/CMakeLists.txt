###################################################################################################
#
# Tests for the Antidunes showcase (apps/showcases/Antidunes) -- analytic geometry / contact math
# that can be unit-tested without a full blockforest/MPI simulation.
#
###################################################################################################

# No codegen dependency (doesn't touch ParticleKernels.h) -- exercises the Bogner state-selection
# classifier used by the free-surface/mesa_pd coupling directly, so it builds unconditionally.
waLBerla_add_test_executable ( ANTIDUNES_BognerStateClassification BognerStateClassificationTest.cpp )
target_link_libraries ( ANTIDUNES_BognerStateClassification
      walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::mesa_pd walberla::lbm_mesapd_coupling )
waLBerla_execute_test( NAME ANTIDUNES_BognerStateClassification )

# Follow-up to ANTIDUNES_CylinderTippingRodJam for job 4255679: is the ~171 deg rod flip bracketed by the real
# [CylForces] log (steps 2700/2900) fast enough, relative to the DEM contact model's own calibrated collision
# duration (tc = collisionTimeFactor * avgParticleDiameter), to plausibly under-resolve the resulting deep overlap?
# Pure arithmetic over real logged numbers -- no particle storage / mesa_pd needed, builds unconditionally.
waLBerla_add_test_executable ( ANTIDUNES_RodSpinSubstepResolution RodSpinSubstepResolutionTest.cpp )
target_link_libraries ( ANTIDUNES_RodSpinSubstepResolution walberla::core )
waLBerla_execute_test( NAME ANTIDUNES_RodSpinSubstepResolution )

# Prompted by DropTest_hopper_cyl_repro_gjk_gas.prm (local reproduction attempt of job 4255679): several rods'
# tiltY diagnostic jumped >100 deg in a single 10-step window with smooth position/velocity/force throughout --
# the signature of tiltY's own atan2 formula being sensitive to a physically-meaningless spin about the rod's own
# axis, not a real topple. Isolates and quantifies that mechanism directly.
waLBerla_add_test_executable ( ANTIDUNES_TiltYOwnAxisSpinArtifact TiltYOwnAxisSpinArtifactTest.cpp )
target_link_libraries ( ANTIDUNES_TiltYOwnAxisSpinArtifact walberla::core walberla::mesa_pd )
waLBerla_execute_test( NAME ANTIDUNES_TiltYOwnAxisSpinArtifact )

if ( WALBERLA_BUILD_WITH_CODEGEN )

   set( ANTIDUNES_SRC_DIR ${walberla_SOURCE_DIR}/apps/showcases/Antidunes )

   waLBerla_add_test_executable ( ANTIDUNES_ClipAnalyticGeometry ClipAnalyticGeometryTest.cpp )
   target_include_directories ( ANTIDUNES_ClipAnalyticGeometry PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_ClipAnalyticGeometry
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_ClipAnalyticGeometry )

   # Replays job 4264997's hopper-clogging collide-stage crash (world (267.5,18.5,19.5), step 2838) through
   # the REAL compiled generated Cumulant collision kernel with the exact captured pre-collide PDF/force
   # state -- checks whether the kernel's own address/arithmetic path (not just a hand-ported formula, see
   # THESIS_PLAN.md's 2026-09-25 entries) can reproduce the NaN from genuinely clean input.
   waLBerla_add_test_executable ( ANTIDUNES_CollideKernelCrashReplay CollideKernelCrashReplayTest.cpp )
   target_include_directories ( ANTIDUNES_CollideKernelCrashReplay PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CollideKernelCrashReplay
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CollideKernelCrashReplay )

   # Tests the "repeated isolated collide iteration" hypothesis for the hopper-clogging crash
   # (THESIS_PLAN.md, 2026-09-25): StreamReconstructAdvectSweep.h never writes PDFs for obstacle
   # cells, and collide runs unconditionally on every cell -- so a permanently-covered cell's PDFs
   # are updated ONLY by repeatedly re-applying collide to its own last output, with zero streaming
   # input and zero force, at this scenario's near-unstable omega=1.99753. Checks whether that
   # feedback loop diverges on its own from an ordinary near-equilibrium starting point.
   waLBerla_add_test_executable ( ANTIDUNES_CollideIsolatedObstacleIteration CollideIsolatedObstacleIterationTest.cpp )
   target_include_directories ( ANTIDUNES_CollideIsolatedObstacleIteration PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CollideIsolatedObstacleIteration
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CollideIsolatedObstacleIteration )

   # Follow-up: re-seeds the same isolated-iteration check with THIS run's own REAL captured OBST-phase
   # state (step 2700, rho=1.007356) instead of a different occurrence's post-reseed value, to check
   # whether repeated collide alone reproduces the real observed per-step drift (~1e-4-5e-4) and eventual
   # collapse ~128 real steps later.
   waLBerla_add_test_executable ( ANTIDUNES_CollideIsolatedObstacleIterationRealSeed CollideIsolatedObstacleIterationRealSeedTest.cpp )
   target_include_directories ( ANTIDUNES_CollideIsolatedObstacleIterationRealSeed PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CollideIsolatedObstacleIterationRealSeed
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CollideIsolatedObstacleIterationRealSeed )

   # [Claude, 2026-09-27] Per explicit request ("can you reproduce this in a small isolated test"):
   # small, isolated reproduction of the 2026-09-27 origin finding -- tracing a real crashed run's
   # box-logger data found that an obstacle cell's ownership handed off directly between two rods
   # (uid...162 -> uid...780) with no intervening liquid/interface state, and a growing checkerboard
   # oscillation began in that SAME step. Reading MovingParticleMapping.h's own cross-particle handoff
   # branch (lines 180-193) found it only updates particleField's owning uid, never the PDF field --
   # this test exercises the REAL MovingParticleMapping class with two real mesa_pd Sphere particles
   # to confirm that directly: a cell handed off between them carries its PDF forward bit-for-bit
   # unchanged.
   waLBerla_add_test_executable ( ANTIDUNES_ObstacleOwnershipHandoff ObstacleOwnershipHandoffTest.cpp )
   target_include_directories ( ANTIDUNES_ObstacleOwnershipHandoff PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_ObstacleOwnershipHandoff
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_ObstacleOwnershipHandoff )

   waLBerla_add_test_executable ( ANTIDUNES_ContactDetection ContactDetectionTest.cpp )
   target_include_directories ( ANTIDUNES_ContactDetection PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_ContactDetection
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_ContactDetection )

   # [Claude, 2026-09-26] Unit tests for mesa_pd::data::Capsule ("pill" shape), per explicit request
   # ("Create a test suite for the pills"). Volume/inertia/support() correctness, ContainsPointFunctor
   # via the exact SingleCast dispatch MovingParticleMapping.h uses (regression coverage for a real
   # bug found this session: Capsule had no ContainsPointFunctor overload in Contains.h, which would
   # have hard-aborted the fluid-obstacle mapping), and contact detection parity with Cylinder
   # (HalfSpace/wall, CylindricalBoundary/pillar, Capsule x Capsule) through the same production
   # SolidCylinderContactDetection dispatch real scenarios use.
   waLBerla_add_test_executable ( ANTIDUNES_CapsuleShape CapsuleShapeTest.cpp )
   target_include_directories ( ANTIDUNES_CapsuleShape PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CapsuleShape
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CapsuleShape )

   # Diagnostic for the cylinder-hopper [RELEASE-OVERLAP] abort (job 4248725): separated parallel Y-axis rods must not report a contact.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderParallelGap CylinderParallelGapTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderParallelGap PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderParallelGap
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderParallelGap )

   # Kodam et al. (2010) Eqs. 55-66 near-parallel (theta ~ 1e-6) band-band-skewed false contacts: literal-paper reference vs implementation.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderNearParallelSkewed CylinderNearParallelSkewedTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderNearParallelSkewed PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderNearParallelSkewed
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderNearParallelSkewed )

   # Closed-form Kodam+Guo Cylinder x Cylinder: for which overlap depths is it correct? (vs exact support-function depth and GJK/EPA)
   waLBerla_add_test_executable ( ANTIDUNES_CylinderClosedFormRange CylinderClosedFormRangeTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderClosedFormRange PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderClosedFormRange
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderClosedFormRange )

   # Why does the closed-form band-edge result jump 0.0034 -> 2.5 LU for (almost) the same geometry? Real logged poses of the launched rods of job 4251049.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderClosedFormJump CylinderClosedFormJumpTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderClosedFormJump PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderClosedFormJump
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderClosedFormJump )

   # Reproduction attempt for job 4255679's [STALE-PARTICLE-IDX] abort (half-density GJK-only hopper, no guards):
   # sweeps a rod's tilt angle through the ~170 deg toppled orientation observed in that run's own logs, checking
   # whether GJK/EPA (the dispatch that run actually used) tracks the true penetration depth or is itself the bug.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderTippingRodJam CylinderTippingRodJamTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderTippingRodJam PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderTippingRodJam
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderTippingRodJam )

   # Geometric precondition of the fluid-side instability in the cylinder hopper: fluid cells squeezed between two obstacles, rods vs spheres.
   waLBerla_add_test_executable ( ANTIDUNES_PinchedFluidCells PinchedFluidCellsTest.cpp )
   target_include_directories ( ANTIDUNES_PinchedFluidCells PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_PinchedFluidCells walberla::core )
   waLBerla_execute_test( NAME ANTIDUNES_PinchedFluidCells )

   # [Claude, 2026-10-29] ANTIDUNES_Job4255679SealedCell removed -- tested the guard24 threshold of
   # SealedCellGuard.h, which was superseded by cylinderContactSkin/boxContactSkin (THESIS_PLAN.md
   # 2026-09-25/09-30: the guard "delays but does NOT eliminate the sealed-cell mechanism" and
   # "checks the wrong signature"; contact skin, implemented next, became the actual production
   # approach and is the one used in every thesis-cited run since). Test moved alongside the guard
   # to apps/showcases/Antidunes/archive/sealed_cell_guard/tests/.

   # Feng-style energy-consistent (Minkowski-overlap) cylinder contact: reference implementation vs GJK/EPA and the closed-form hybrid (test only, not used by the simulation).
   waLBerla_add_test_executable ( ANTIDUNES_EnergyConsistentContact EnergyConsistentContactTest.cpp )
   target_include_directories ( ANTIDUNES_EnergyConsistentContact PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_EnergyConsistentContact
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_EnergyConsistentContact )

   # Diagnostic: deeply interpenetrated rod pairs from the cylinder-hopper resume run -- analytic hybrid vs GJK/EPA.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderDeepOverlapProbe CylinderDeepOverlapProbeTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderDeepOverlapProbe PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderDeepOverlapProbe
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )

   # [Claude, 2026-10-10] Diagnostic: A/B cylinder pair of the round-container drop-test energy spike (job 4292501).
   waLBerla_add_test_executable ( ANTIDUNES_KodamRoundSpikeProbe KodamRoundSpikeProbeTest.cpp )
   target_include_directories ( ANTIDUNES_KodamRoundSpikeProbe PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_KodamRoundSpikeProbe
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )

   # ForcingAdjuster's mid-run PID switch-off (PIDParameters.switchOffStep): freezes acceleration_x at the
   # mean of the PID's recent commanded values, including across checkpoint resumes. Needs PIDController.cpp
   # compiled in (it is a separate translation unit in the app, not a library).
   waLBerla_add_test_executable ( ANTIDUNES_ForcingAdjuster ForcingAdjusterTest.cpp ${ANTIDUNES_SRC_DIR}/PIDController.cpp )
   target_include_directories ( ANTIDUNES_ForcingAdjuster PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_ForcingAdjuster
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_ForcingAdjuster )

   # Closed-form Box x Box (OBBTree SAT) and Cylinder x Cylinder / Cylinder x Plane (Kodam et al.
   # 2010) algorithms from PaperAnalyticalContacts.h, compared here against
   # SolidCylinderContactDetection's existing GJK/EPA-fallback results for identical geometry. Used by the
   # simulation only when ModelParameters.useAnalyticalContactDetection is set (opt-in, default false) -- see
   # PaperAnalyticalContacts.h's file header. [Claude, 2026-09-30] comment corrected (said "NOT used by the simulation").
   # ANTIDUNES_REFERENCE_DUMP=<dir> additionally dumps random pairs for the thesis's independent brute-force reference.
   waLBerla_add_test_executable ( ANTIDUNES_PaperAnalyticalContacts PaperAnalyticalContactTest.cpp )
   target_include_directories ( ANTIDUNES_PaperAnalyticalContacts PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_PaperAnalyticalContacts
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_PaperAnalyticalContacts )

   # Microbenchmark (not a correctness test -- see PaperAnalyticalContactTest.cpp for that): times
   # each PaperAnalyticalContacts.h cylinder-cylinder analytical sub-case, the GJK/EPA fallback path,
   # and box-box SAT/GJK for comparison, to find which part of the cylinder-cylinder contact pipeline
   # actually dominates a run's cost -- see CylinderCylinderContactBenchmark.cpp's own file header for
   # the widerstacks-shape-ablation finding (cylinder+hybrid ran slower than box+forced-GJK) that
   # motivated writing this.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderCylinderContactBenchmark CylinderCylinderContactBenchmark.cpp )
   target_include_directories ( ANTIDUNES_CylinderCylinderContactBenchmark PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderCylinderContactBenchmark
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderCylinderContactBenchmark )

   # Standalone time/accuracy comparison of the three quartic-root solvers available in
   # PaperAnalyticalContacts.h::detail (Durand-Kerner, Ferrari closed-form, Aberth-Ehrlich) against
   # KNOWN synthetic roots -- not a correctness test of any production path (this header is unused by
   # the simulation), purely informs whether solveQuarticDurandKerner should ever be replaced.
   waLBerla_add_test_executable ( ANTIDUNES_QuarticSolverComparison QuarticSolverComparisonTest.cpp )
   target_include_directories ( ANTIDUNES_QuarticSolverComparison PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_QuarticSolverComparison
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_QuarticSolverComparison )

   # [Claude, 2026-10-29] ANTIDUNES_TsunamiDamBreak removed -- tested TsunamiDataInitialization.h,
   # which was archived (unused by any thesis-cited result; see
   # apps/showcases/Antidunes/archive/ARCHIVE_README.md). Test moved alongside it to
   # apps/showcases/Antidunes/archive/tsunami_pole_forest_driftwood/tests/.

endif ()

# [Claude, 2026-09-21] Equivalence of the optional interpolated bounce-back in SimpleBB (useInterpolatedBounceBack) with lbm_mesapd_coupling::CurvedLinear on the sphere drag-force test (copy of
# LBM_MESAPD_COUPLING_MEM_DragForceSphere with a third variant IBB). Build only; compare the printed 'IBB-EQUIV-DRAG' lines of `--MEMVariant CLI` and `--MEMVariant IBB` (run: --shortrun or full, 2 processes).
waLBerla_add_test_executable ( ANTIDUNES_IbbEquivalenceDrag IbbEquivalenceDragTest.cpp )
target_link_libraries ( ANTIDUNES_IbbEquivalenceDrag walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
