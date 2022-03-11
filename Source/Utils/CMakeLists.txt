target_sources(WarpX
  PRIVATE
    CoarsenIO.cpp
    CoarsenMR.cpp
    Interpolate.cpp
    IntervalsParser.cpp
    ParticleUtils.cpp
    RelativeCellPosition.cpp
    WarnManager.cpp
    WarpXAlgorithmSelection.cpp
    WarpXMovingWindow.cpp
    WarpXTagging.cpp
    WarpXUtil.cpp
    WarpXVersion.cpp
)

add_subdirectory(MsgLogger)

target_sources(ablastr
  PRIVATE
    MPIInitHelpers.cpp
)
