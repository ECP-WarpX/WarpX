target_sources(WarpX
  PRIVATE
    MultiParticleContainer.cpp
    PhotonParticleContainer.cpp
    PhysicalParticleContainer.cpp
    RigidInjectedParticleContainer.cpp
    WarpXParticleContainer.cpp
)

add_subdirectory(Collision)
#add_subdirectory(Deposition)
add_subdirectory(ElementaryProcess)
add_subdirectory(Gather)
add_subdirectory(ParticleCreation)
#add_subdirectory(Pusher)
add_subdirectory(Resampling)
add_subdirectory(Sorting)
