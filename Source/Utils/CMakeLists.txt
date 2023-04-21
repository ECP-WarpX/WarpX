target_sources(WarpX
  PRIVATE
    Interpolate.cpp
    MPIInitHelpers.cpp
    ParticleUtils.cpp
    RelativeCellPosition.cpp
    WarpXAlgorithmSelection.cpp
    WarpXMovingWindow.cpp
    WarpXTagging.cpp
    WarpXUtil.cpp
    WarpXrocfftUtil.cpp
    WarpXVersion.cpp
)

add_subdirectory(Logo)
add_subdirectory(Parser)
