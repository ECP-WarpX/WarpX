foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        BeamRelevant.cpp
        ColliderRelevant.cpp
        DifferentialLuminosity.cpp
        FieldEnergy.cpp
        FieldProbe.cpp
        FieldProbeParticleContainer.cpp
        FieldMomentum.cpp
        LoadBalanceCosts.cpp
        LoadBalanceEfficiency.cpp
        MultiReducedDiags.cpp
        ParticleEnergy.cpp
        ParticleMomentum.cpp
        ParticleHistogram.cpp
        ParticleHistogram2D.cpp
        ReducedDiags.cpp
        FieldMaximum.cpp
        ParticleExtrema.cpp
        RhoMaximum.cpp
        ParticleNumber.cpp
        FieldReduction.cpp
        FieldProbe.cpp
        ChargeOnEB.cpp
    )
endforeach()
