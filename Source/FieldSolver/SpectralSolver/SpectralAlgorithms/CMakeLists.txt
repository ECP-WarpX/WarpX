target_sources(WarpX
  PRIVATE
    PsatdAlgorithmJConstantInTime.cpp
    PsatdAlgorithmJLinearInTime.cpp
    PsatdAlgorithmPml.cpp
    SpectralBaseAlgorithm.cpp
    PsatdAlgorithmComoving.cpp
)

if(WarpX_DIMS STREQUAL RZ)
    target_sources(WarpX
      PRIVATE
        SpectralBaseAlgorithmRZ.cpp
        PsatdAlgorithmRZ.cpp
        PsatdAlgorithmGalileanRZ.cpp
        PsatdAlgorithmPmlRZ.cpp
    )
endif()
