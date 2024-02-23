foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        FlushFormatAscent.cpp
        FlushFormatCheckpoint.cpp
        FlushFormatPlotfile.cpp
        FlushFormatSensei.cpp
        # may need to put a flag and ifdef for PlotPlus
        FlushFormatPlotPlus.cpp
    )

    if(WarpX_HAVE_OPENPMD)
        target_sources(lib_${SD}
          PRIVATE
            FlushFormatOpenPMD.cpp
        )
    endif()
endforeach()
