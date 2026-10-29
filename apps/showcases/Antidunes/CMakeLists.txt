waLBerla_link_files_to_builddir ( *.prm )

if ( WALBERLA_BUILD_WITH_CODEGEN )
   walberla_generate_target_from_python ( NAME AntidunesLatticeModelGeneration
         FILE AntidunesLatticeModelGeneration.py
         OUT_FILES AntidunesLatticeModel.cpp AntidunesLatticeModel.h )

   add_executable ( Antidunes Antidunes.cpp PIDController.cpp )
   target_link_libraries ( Antidunes walberla::blockforest walberla::boundary walberla::core
           walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
           walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk AntidunesLatticeModelGeneration )

endif ()

add_executable ( BedGeneration BedGeneration.cpp )
target_link_libraries ( BedGeneration walberla::blockforest walberla::core walberla::domain_decomposition
        walberla::mesa_pd walberla::vtk )

if ( WALBERLA_BUILD_WITH_CODEGEN )
   # Fluid-free DEM validation dry-run (DEMValidationScenarios.h, AnalyticalSolutions.h). Reuses
   # SolidCylinderContactDetection/ComputeContacts from ParticleKernels.h verbatim, which transitively
   # pulls in the generated LBM lattice model -- same reason ContactDetectionTest.cpp and
   # PaperAnalyticalContactTest.cpp require this flag despite not touching LBM directly.
   add_executable ( DEMValidation DEMValidation.cpp )
   target_link_libraries ( DEMValidation walberla::blockforest walberla::boundary walberla::core
           walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
           walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk AntidunesLatticeModelGeneration )

   # Fluid-free, DEM-only cuboid drop test (DropTest.cpp, SetupStages_drop_test.h). Reuses
   # ParticleKernels.h/doDEMSubCycles verbatim, which transitively pulls in the generated LBM lattice
   # model -- same reason DEMValidation above requires this flag despite not touching LBM directly.
   # PIDController.cpp is linked (though DropTest never constructs a PIDController at runtime) because
   # SetupStages_drop_test.h's setupSimState() -- an ordinary function DropTest.cpp never calls, but
   # which still compiles into the object file since it's not a template -- references it.
   add_executable ( DropTest DropTest.cpp PIDController.cpp )
   target_link_libraries ( DropTest walberla::blockforest walberla::boundary walberla::core
           walberla::domain_decomposition walberla::field walberla::lbm walberla::mesa_pd walberla::lbm_mesapd_coupling
           walberla::postprocessing walberla::sqlite walberla::timeloop walberla::vtk AntidunesLatticeModelGeneration )
endif ()
