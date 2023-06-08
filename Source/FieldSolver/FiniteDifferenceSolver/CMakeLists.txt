foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        ComputeDivE.cpp
        EvolveB.cpp
        EvolveBPML.cpp
        EvolveE.cpp
        EvolveEPML.cpp
        EvolveF.cpp
        EvolveFPML.cpp
        EvolveG.cpp
        EvolveECTRho.cpp
        FiniteDifferenceSolver.cpp
        MacroscopicEvolveE.cpp
        ApplySilverMuellerBoundary.cpp
    )
endforeach()

add_subdirectory(MacroscopicProperties)
