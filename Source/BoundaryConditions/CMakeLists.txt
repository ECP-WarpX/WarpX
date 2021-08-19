target_sources(WarpX
  PRIVATE
    PML.cpp
    WarpXEvolvePML.cpp
    WarpXFieldBoundaries.cpp
    WarpX_PEC.cpp
)

if(WarpX_DIMS STREQUAL RZ)
    target_sources(WarpX
      PRIVATE
        PML_RZ.cpp
    )
endif()
