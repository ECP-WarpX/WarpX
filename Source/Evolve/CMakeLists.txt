foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        WarpXEvolve.cpp
        WarpXComputeDt.cpp
        WarpXOneStepImplicitPicard.cpp
    )
endforeach()
