foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    if(ABLASTR_FFT AND D EQUAL 3)
        target_sources(ablastr_${SD}
          PRIVATE
            IntegratedGreenFunctionSolver.cpp
        )
    endif()
endforeach()
