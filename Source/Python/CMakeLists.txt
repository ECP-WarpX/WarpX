foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        WarpX_py.cpp
        WarpXWrappers.cpp
    )
endforeach()
