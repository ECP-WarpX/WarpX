target_sources(WarpX
  PRIVATE
    PML.cpp
    WarpXEvolvePML.cpp
    WarpXFieldBoundaries.cpp
    WarpX_PerfectConductivity.cpp
)

if(WarpX_DIMS STREQUAL RZ)
    target_sources(WarpX
      PRIVATE
        PML_RZ.cpp
    )
endif()
