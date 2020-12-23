target_sources(WarpX
  PRIVATE
    FlushFormatAscent.cpp
    FlushFormatCheckpoint.cpp
    FlushFormatPlotfile.cpp
)

if(WarpX_HAVE_OPENPMD)
    target_sources(WarpX
      PRIVATE
        FlushFormatOpenPMD.cpp
    )
endif()
