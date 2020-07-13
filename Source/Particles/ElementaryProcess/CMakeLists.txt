target_sources(WarpX
  PRIVATE
    Ionization.cpp
)

if(WarpX_HAVE_QED)
    target_sources(WarpX
      PRIVATE
        QEDPairGeneration.cpp
        QEDPhotonEmission.cpp
    )
    add_subdirectory(QEDInternals)
endif()
