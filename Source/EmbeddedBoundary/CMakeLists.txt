foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        WarpXInitEB.cpp
        WarpXFaceExtensions.cpp
        WarpXFaceInfoBox.H
    )
endforeach()
