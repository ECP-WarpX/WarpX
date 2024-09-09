foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        ElectrostaticBase.cpp
        ElectrostaticSolver.cpp
        LabFrameExplicitES.cpp
        PoissonBoundaryHandler.cpp
        RelativisticExplicitES.cpp
    )
endforeach()
