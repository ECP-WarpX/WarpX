target_sources(WarpX
  PRIVATE
    BackTransformedDiagnostic.cpp
    Diagnostics.cpp
    FieldIO.cpp
    FullDiagnostics.cpp
    MultiDiagnostics.cpp
    ParticleIO.cpp
    SliceDiagnostic.cpp
    WarpXIO.cpp
    WarpXOpenPMD.cpp
    BTDiagnostics.cpp
)

add_subdirectory(ComputeDiagFunctors)
add_subdirectory(FlushFormats)
add_subdirectory(ParticleDiag)
add_subdirectory(ReducedDiags)
