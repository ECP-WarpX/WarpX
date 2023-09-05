foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        ElectrostaticSolver.cpp
        WarpXPushFieldsEM.cpp
        WarpXPushFieldsHybridPIC.cpp
        WarpX_QED_Field_Pushers.cpp
    )
endforeach()

add_subdirectory(FiniteDifferenceSolver)
add_subdirectory(MagnetostaticSolver)
if(WarpX_PSATD)
    add_subdirectory(SpectralSolver)
endif()
