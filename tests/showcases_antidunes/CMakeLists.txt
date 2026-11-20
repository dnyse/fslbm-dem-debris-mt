###################################################################################################
#
# Tests for the Antidunes showcase (apps/showcases/Antidunes): free-surface reconstruction classifiers,
# contact detection and geometry, PID forcing, MEM coupling pieces, diagnostics and microbenchmarks, each
# without running the full simulation.
#
###################################################################################################

# No codegen dependency (doesn't touch ParticleKernels.h) -- exercises the Bogner state-selection
# classifier used by the free-surface/mesa_pd coupling directly, so it builds unconditionally.
waLBerla_add_test_executable ( ANTIDUNES_BognerStateClassification BognerStateClassificationTest.cpp )
target_link_libraries ( ANTIDUNES_BognerStateClassification
      walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::mesa_pd walberla::lbm_mesapd_coupling )
waLBerla_execute_test( NAME ANTIDUNES_BognerStateClassification )

# Full-3D-grid stress test (sphere/cylinder/cuboid floater sizes from the prm files) of the
# strictBognerRelaxInterface pass ("Round 2") in PdfReconstructionManager.h, combined with
# strictBognerClassification ("Round 0") and strictBognerNoLiquidAtGas ("Round 1"), compared with pure
# strict Bogner and the non-strict Pass-2+Pass-3 path over the same releases: the Round 0+1+2 combination
# never leaves a liquid cell touching a gas cell. Pure classifier functions only, so it builds unconditionally.
waLBerla_add_test_executable ( ANTIDUNES_LiquidPromoteRelaxClassification LiquidPromoteRelaxClassificationTest.cpp )
target_link_libraries ( ANTIDUNES_LiquidPromoteRelaxClassification
      walberla::core walberla::domain_decomposition walberla::field walberla::lbm
      walberla::mesa_pd walberla::lbm_mesapd_coupling )
waLBerla_execute_test( NAME ANTIDUNES_LiquidPromoteRelaxClassification )

# Companion to ANTIDUNES_CylinderTippingRodJam: compares the average rotation rate of a logged ~171 deg rod
# flip with the DEM collision duration tc = collisionTimeFactor * avgParticleDiameter, i.e. whether the
# overlap could be under-resolved by the sub-cycles. Pure arithmetic, no mesa_pd, builds unconditionally.
waLBerla_add_test_executable ( ANTIDUNES_RodSpinSubstepResolution RodSpinSubstepResolutionTest.cpp )
target_link_libraries ( ANTIDUNES_RodSpinSubstepResolution walberla::core )
waLBerla_execute_test( NAME ANTIDUNES_RodSpinSubstepResolution )

# Shows that the tiltY diagnostic (atan2 formula in ParticleKernels.h) jumps under a physically meaningless
# spin of a rod about its own axis, so a large tiltY change need not be a real topple.
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

   # Replays a recorded hopper-clogging collide-stage crash cell through the compiled generated Cumulant
   # collision kernel with the recorded pre-collide PDF/force state: can the kernel itself turn this
   # near-equilibrium input into NaN?
   waLBerla_add_test_executable ( ANTIDUNES_CollideKernelCrashReplay CollideKernelCrashReplayTest.cpp )
   target_include_directories ( ANTIDUNES_CollideKernelCrashReplay PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CollideKernelCrashReplay
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CollideKernelCrashReplay )

   # Repeated isolated collide iteration: StreamReconstructAdvectSweep.h never writes PDFs for obstacle
   # cells and collide runs on every cell, so a permanently covered cell's PDFs are updated only by
   # re-applying collide to its own output (zero streaming input, zero force) at omega = 1.99753. Checks
   # whether that iteration diverges from a near-equilibrium start. Diagnostic only.
   waLBerla_add_test_executable ( ANTIDUNES_CollideIsolatedObstacleIteration CollideIsolatedObstacleIterationTest.cpp )
   target_include_directories ( ANTIDUNES_CollideIsolatedObstacleIteration PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CollideIsolatedObstacleIteration
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CollideIsolatedObstacleIteration )

   # Same isolated iteration, seeded with the recorded state of the collapsed cell at the start of its
   # obstacle phase: does repeated collide alone reproduce the recorded per-step drift and collapse?
   waLBerla_add_test_executable ( ANTIDUNES_CollideIsolatedObstacleIterationRealSeed CollideIsolatedObstacleIterationRealSeedTest.cpp )
   target_include_directories ( ANTIDUNES_CollideIsolatedObstacleIterationRealSeed PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CollideIsolatedObstacleIterationRealSeed
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CollideIsolatedObstacleIterationRealSeed )

   # Obstacle-to-obstacle ownership handoff in MovingParticleMapping with two mesa_pd spheres: on a
   # genuine handoff the onOwnershipHandoff callback (useOwnershipHandoffPdfReset) resets the cell's PDF
   # to rest equilibrium, while a particle re-confirming its own cell leaves the PDF unchanged.
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

   # Unit tests for mesa_pd::data::Capsule: volume/inertia, support(), the ContainsPointFunctor overload
   # (Contains.h) through the SingleCast dispatch MovingParticleMapping.h uses, and contact detection
   # (Capsule x HalfSpace, x CylindricalBoundary, x Capsule) through SolidCylinderContactDetection.
   waLBerla_add_test_executable ( ANTIDUNES_CapsuleShape CapsuleShapeTest.cpp )
   target_include_directories ( ANTIDUNES_CapsuleShape PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CapsuleShape
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CapsuleShape )

   # Separated parallel Y-axis rods must not report a contact (analytic hybrid on and off).
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

   # Closed-form band-edge overlap jump (0.0034 -> 2.5 LU) for nearly the same geometry, on logged rod poses (CylinderClosedFormJumpData.h).
   waLBerla_add_test_executable ( ANTIDUNES_CylinderClosedFormJump CylinderClosedFormJumpTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderClosedFormJump PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderClosedFormJump
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderClosedFormJump )

   # Approximated recorded rod-hopper configuration (GJK/EPA-only dispatch): sweeps one rod's tilt through the
   # ~170 deg toppled orientation and checks whether GJK/EPA tracks the true penetration depth.
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


   # Feng-style energy-consistent (Minkowski-overlap) cylinder contact: reference implementation vs GJK/EPA and the closed-form hybrid (test only, not used by the simulation).
   waLBerla_add_test_executable ( ANTIDUNES_EnergyConsistentContact EnergyConsistentContactTest.cpp )
   target_include_directories ( ANTIDUNES_EnergyConsistentContact PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_EnergyConsistentContact
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_EnergyConsistentContact )

   # Diagnostic (build only): deeply interpenetrated rod pairs from a recorded cylinder-hopper run -- analytic hybrid vs GJK/EPA.
   waLBerla_add_test_executable ( ANTIDUNES_CylinderDeepOverlapProbe CylinderDeepOverlapProbeTest.cpp )
   target_include_directories ( ANTIDUNES_CylinderDeepOverlapProbe PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderDeepOverlapProbe
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )

   # Diagnostic (build only): stacked cylinder pair of a round-container drop-test energy spike -- closed form, refined, dispatch, GJK/EPA.
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

   # Closed-form Box x Box (OBB SAT) and Cylinder x Cylinder / Cylinder x Plane (Kodam et al. 2010, Guo et al.
   # 2012) contacts of ClosedFormContacts*.h against GJK/EPA (SolidCylinderContactDetection) and hand-derived
   # cases. The simulation uses them only with ModelParameters.useAnalyticalContactDetection (default false).
   # ANTIDUNES_REFERENCE_DUMP=<dir> additionally dumps random pairs for an independent brute-force reference.
   waLBerla_add_test_executable ( ANTIDUNES_PaperAnalyticalContacts PaperAnalyticalContactTest.cpp )
   target_include_directories ( ANTIDUNES_PaperAnalyticalContacts PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_PaperAnalyticalContacts
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_PaperAnalyticalContacts )

   # Microbenchmark (not a correctness test -- see PaperAnalyticalContactTest.cpp for that): times each
   # cylinder-cylinder closed-form sub-case of ClosedFormContactsCylinderCylinder.h, the GJK/EPA fallback
   # and box-box SAT/GJK on fixed geometries (motivation in the file header).
   waLBerla_add_test_executable ( ANTIDUNES_CylinderCylinderContactBenchmark CylinderCylinderContactBenchmark.cpp )
   target_include_directories ( ANTIDUNES_CylinderCylinderContactBenchmark PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_CylinderCylinderContactBenchmark
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_CylinderCylinderContactBenchmark )

   # Microbenchmark (build only, not a correctness test): closed form vs. GJK/EPA averaged over a 100000-trial
   # random box_pairs.csv/cyl_pairs.csv dump (ANTIDUNES_PAIRS_DIR) of the accuracy validation, a realistic mix
   # of easy and hard cases. The file header states what is and isn't timed.
   waLBerla_add_test_executable ( ANTIDUNES_ContactSpeed100kBenchmark ContactSpeed100kBenchmark.cpp )
   target_include_directories ( ANTIDUNES_ContactSpeed100kBenchmark PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_ContactSpeed100kBenchmark
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )

   # Time/accuracy comparison of the three quartic root finders in the `detail` namespace of
   # ClosedFormContactsCylinderPlaneGeometry.h (Durand-Kerner, Ferrari closed form -- the production default --,
   # Aberth-Ehrlich) against known synthetic roots. Reports only, no assertions.
   waLBerla_add_test_executable ( ANTIDUNES_QuarticSolverComparison QuarticSolverComparisonTest.cpp )
   target_include_directories ( ANTIDUNES_QuarticSolverComparison PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_QuarticSolverComparison
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_QuarticSolverComparison )


endif ()

# Equivalence of the optional interpolated bounce-back in SimpleBB (useInterpolatedBounceBack) with lbm_mesapd_coupling::CurvedLinear on the sphere drag-force test (adapted from
# LBM_MESAPD_COUPLING_MEM_DragForceSphere with a third variant IBB). Build only; compare the printed 'IBB-EQUIV-DRAG' lines of `--MEMVariant CLI` and `--MEMVariant IBB` (run: --shortrun or full, 2 processes).
waLBerla_add_test_executable ( ANTIDUNES_IbbEquivalenceDrag IbbEquivalenceDragTest.cpp )
target_link_libraries ( ANTIDUNES_IbbEquivalenceDrag walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
