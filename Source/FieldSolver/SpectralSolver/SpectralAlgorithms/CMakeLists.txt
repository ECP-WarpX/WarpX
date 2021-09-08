target_sources(WarpX
  PRIVATE
    PsatdAlgorithm.cpp
    PMLPsatdAlgorithm.cpp
    SpectralBaseAlgorithm.cpp
    ComovingPsatdAlgorithm.cpp
)

if(WarpX_DIMS STREQUAL RZ)
    target_sources(WarpX
      PRIVATE
        SpectralBaseAlgorithmRZ.cpp
        PsatdAlgorithmRZ.cpp
        GalileanPsatdAlgorithmRZ.cpp
    )
endif()
