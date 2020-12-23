target_sources(WarpX
  PRIVATE
    GalileanAlgorithm.cpp
    PMLPsatdAlgorithm.cpp
    PsatdAlgorithm.cpp
    SpectralBaseAlgorithm.cpp
    AvgGalileanAlgorithm.cpp
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
