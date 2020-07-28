target_sources(WarpX
  PRIVATE
    SpectralFieldData.cpp
    SpectralKSpace.cpp
    SpectralSolver.cpp
)

if(ENABLE_CUDA)
    target_sources(WarpX
      PRIVATE
        WrapCuFFT.cpp
    )
else()
    target_sources(WarpX
      PRIVATE
        WrapFFTW.cpp
    )
endif()

if(WarpX_DIMS STREQUAL RZ)
    target_sources(WarpX
      PRIVATE
        SpectralSolverRZ.cpp
        SpectralFieldDataRZ.cpp
        SpectralKSpaceRZ.cpp
        SpectralBinomialFilter.cpp
    )
    add_subdirectory(SpectralHankelTransform)
endif()

add_subdirectory(SpectralAlgorithms)
