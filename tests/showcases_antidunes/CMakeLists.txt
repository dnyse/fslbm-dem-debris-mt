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

   waLBerla_add_test_executable ( ANTIDUNES_ContactDetection ContactDetectionTest.cpp )
   target_include_directories ( ANTIDUNES_ContactDetection PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_ContactDetection
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_ContactDetection )

   # Closed-form Box x Box (OBBTree SAT) and Cylinder x Cylinder / Cylinder x Plane (Kodam et al.
   # 2010) algorithms from PaperAnalyticalContacts.h, compared here against
   # SolidCylinderContactDetection's existing GJK/EPA-fallback results for identical geometry. NOT
   # used by the simulation -- see PaperAnalyticalContacts.h's file header.
   waLBerla_add_test_executable ( ANTIDUNES_PaperAnalyticalContacts PaperAnalyticalContactTest.cpp )
   target_include_directories ( ANTIDUNES_PaperAnalyticalContacts PRIVATE ${ANTIDUNES_SRC_DIR} )
   target_link_libraries ( ANTIDUNES_PaperAnalyticalContacts
         walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
         walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
         walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk
         AntidunesLatticeModelGeneration )
   waLBerla_execute_test( NAME ANTIDUNES_PaperAnalyticalContacts )

endif ()
