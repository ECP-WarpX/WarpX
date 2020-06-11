target_sources(WarpX
  PRIVATE
    ComputeDivE.cpp
    EvolveB.cpp
    EvolveBPML.cpp
    EvolveE.cpp
    EvolveEPML.cpp
    EvolveF.cpp
    EvolveFPML.cpp
    FiniteDifferenceSolver.cpp
    MacroscopicEvolveE.cpp
)

add_subdirectory(MacroscopicProperties)
