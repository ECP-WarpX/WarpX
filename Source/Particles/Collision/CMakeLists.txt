target_sources(WarpX
  PRIVATE
    CollisionHandler.cpp
    CollisionBase.cpp
    BackgroundMCCCollision.cpp
    MCCProcess.cpp
)

add_subdirectory(BinaryCollision)
