foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        Interpolate.cpp
        ParticleUtils.cpp
        RelativeCellPosition.cpp
        WarpXAlgorithmSelection.cpp
        WarpXMovingWindow.cpp
        WarpXTagging.cpp
        WarpXUtil.cpp
        WarpXrocfftUtil.cpp
        WarpXVersion.cpp
    )
endforeach()

add_subdirectory(Logo)
add_subdirectory(Parser)
