target_sources(WarpX
  PRIVATE
    CoarsenIO.cpp
    CoarsenMR.cpp
    Interpolate.cpp
    IntervalsParser.cpp
    MPIInitHelpers.cpp
    ParticleUtils.cpp
    RelativeCellPosition.cpp
    WarnManager.cpp
    WarpXAlgorithmSelection.cpp
    WarpXMovingWindow.cpp
    WarpXTagging.cpp
    WarpXUtil.cpp
)

add_subdirectory(MsgLogger)
