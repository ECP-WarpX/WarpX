target_sources(WarpX
  PRIVATE
    Ionization.cpp
    QEDPairGeneration.cpp
    QEDPhotonEmission.cpp
)

add_subdirectory(QEDInternals)

