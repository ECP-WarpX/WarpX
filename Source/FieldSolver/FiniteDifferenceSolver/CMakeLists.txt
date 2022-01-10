target_sources(WarpX
  PRIVATE
    ComputeDivE.cpp
    EvolveB.cpp
    EvolveBPML.cpp
    EvolveE.cpp
    EvolveEPML.cpp
    EvolveF.cpp
    EvolveFPML.cpp
    EvolveG.cpp
    FiniteDifferenceSolver.cpp
    MacroscopicEvolveE.cpp
    MacroscopicEvolveEPML.cpp
    ApplySilverMuellerBoundary.cpp
)

add_subdirectory(MacroscopicProperties)
