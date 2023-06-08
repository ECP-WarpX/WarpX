foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        CellCenterFunctor.cpp
        DivBFunctor.cpp
        DivEFunctor.cpp
        RhoFunctor.cpp
        PartPerCellFunctor.cpp
        PartPerGridFunctor.cpp
        BackTransformFunctor.cpp
        BackTransformParticleFunctor.cpp
        ParticleReductionFunctor.cpp
    )
endforeach()
