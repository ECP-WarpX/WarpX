target_sources(WarpX
  PRIVATE
    SpectralFieldData.cpp
    SpectralKSpace.cpp
    SpectralSolver.cpp
)

add_subdirectory(SpectralAlgorithms)
