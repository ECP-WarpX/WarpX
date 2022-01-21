target_sources(WarpX
  PRIVATE
    PushPX.cpp
)

if(WarpX_QED)
    target_sources(WarpX
      PRIVATE
      PushPXQED.cpp
    )
endif()
