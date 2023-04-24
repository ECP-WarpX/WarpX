target_sources(ablastr
  PRIVATE
    Communication.cpp
    SignalHandling.cpp
    TextMsg.cpp
    UsedInputsFile.cpp
)

add_subdirectory(msg_logger)
add_subdirectory(text)
add_subdirectory(timer)
