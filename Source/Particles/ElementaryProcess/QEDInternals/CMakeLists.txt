target_sources(WarpX
  PRIVATE
    BreitWheelerEngineWrapper.cpp
    QuantumSyncEngineWrapper.cpp
)

#if(WARPX_QED_TABLE_GEN)
#    target_sources(WarpX
#      PRIVATE
#        BreitWheelerEngineTableBuilder.cpp
#        QuantumSyncEngineTableBuilder.cpp
#endif()
