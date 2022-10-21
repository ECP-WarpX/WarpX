###############################################################################
# These are the files equivalent to the primary C++ implementation - but here
# we define how they will appear in our Python module (aka Python bindings).
#
foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        WarpX_py.cpp
    )
    if(WarpX_PYTHON)
        target_sources(pyWarpX_${SD}
          PRIVATE
            WarpX.cpp
        )
    endif()
endforeach()
