target_sources(WarpX
  PRIVATE
    CollisionHandler.cpp
    CollisionBase.cpp
)

add_subdirectory(BinaryCollision)
add_subdirectory(BackgroundMCC)
add_subdirectory(BackgroundStopping)
