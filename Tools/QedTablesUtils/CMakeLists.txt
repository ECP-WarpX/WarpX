########
# Common argument parser library
########
set(arg_parser_cpp "Source/ArgParser/QedTablesArgParser.cpp")
set(arg_parser_h   "Source/ArgParser/QedTablesArgParser.H")

add_library( QedToolsArgParser ${arg_parser_cpp} ${arg_parser_h} )

########
# Build QED lookup tables generator
########
set(qed_generator_name "qed_table_generator")
set(qed_generator_src  "Source/QedTableGenerator.cpp")

add_executable(${qed_generator_name} ${qed_generator_src})

target_link_libraries(${qed_generator_name} PRIVATE PXRMP_QED)
target_link_libraries(${qed_generator_name} PRIVATE QedToolsArgParser)

target_compile_features(${qed_generator_name} PUBLIC cxx_std_17)
set_target_properties(${qed_generator_name} PROPERTIES CXX_EXTENSIONS OFF)

set_target_properties(${qed_generator_name} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/qed_tools)

message(STATUS "QED lookup tables generator: ${qed_generator_name} (${qed_generator_src})")

########
# Build QED lookup tables reader
########
set(qed_reader_name "qed_table_reader")
set(qed_reader_src  "Source/QedTableReader.cpp")

add_executable(${qed_reader_name} ${qed_reader_src})

target_link_libraries(${qed_reader_name} PRIVATE PXRMP_QED)
target_link_libraries(${qed_reader_name} PRIVATE QedToolsArgParser)

target_compile_features(${qed_reader_name} PUBLIC cxx_std_17)
set_target_properties(${qed_reader_name} PROPERTIES CXX_EXTENSIONS OFF)

set_target_properties(${qed_reader_name} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/qed_tools)

message(STATUS "QED lookup tables generator: ${qed_reader_name} (${qed_reader_src})")
