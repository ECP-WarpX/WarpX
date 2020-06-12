target_sources(WarpX
  PRIVATE
    GalileanAlgorithm.cpp
    PMLPsatdAlgorithm.cpp
    PsatdAlgorithm.cpp
    SpectralBaseAlgorithm.cpp
)

if(WarpX_DIMS STREQUAL RZ)
    target_sources(WarpX
      PRIVATE
        SpectralBaseAlgorithmRZ.cpp
        PsatdAlgorithmRZ.cpp
    )
endif()
