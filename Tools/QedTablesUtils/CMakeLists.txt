########
# Header
########

cmake_minimum_required(VERSION 3.12)

project(
  QedTableTools
  DESCRIPTION
    "Tools to generate or parse QED lookup tables for WarpX+PICSAR")

#____________________________________________________________________________________________________________


########
# Options
########

option(TOOLS_QED_GENERATOR        "Builds the QED lookup tables generator (needs Boost)" ON)
option(TOOLS_QED_READER           "Builds the QED lookup tables reader" ON)
option(TOOLS_QED_OMP              "Enable OpenMP support for lookup tables generator" ON)
option(TOOLS_QED_PICSAR_internal  "Download & build PICSAR" ON)

#____________________________________________________________________________________________________________


########
# Download PICSAR library
########

# This macro downloads PICSAR from git if a TOOLS_QED_picsar_src is not provided
macro(find_picsar)

    if(TOOLS_QED_PICSAR_src)
        message(STATUS "Compiling local PICSAR ...")
        message(STATUS "PICSAR source path: ${TOOLS_QED_PICSAR_src}")
    elseif(TOOLS_QED_PICSAR_internal)
        message(STATUS "Downloading PICSAR ...")
        message(STATUS "PICSAR repository: ${TOOLS_QED_PICSAR_repo} (${TOOLS_QED_PICSAR_branch})")
        include(FetchContent)
    else()
        message(FATAL_ERROR "Either TOOLS_QED_PICSAR_src or TOOLS_QED_PICSAR_internal must be defined!")
    endif()

    set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

    # Enable or disable QED lookup tables generation
    # If table generation is enabled, enable or disable
    # openMP support depending on TOOLS_QED_OMP
    if(TOOLS_QED_GENERATOR)
        set(PXRMP_QED_TABLEGEN ON CACHE INTERNAL "")
        if(TOOLS_QED_OMP)
            set(PXRMP_QED_OMP ON CACHE INTERNAL "")
        else()
            set(PXRMP_QED_OMP OFF CACHE INTERNAL "")
        endif()
    else()
        set(PXRMP_QED_TABLEGEN OFF CACHE INTERNAL "")
        set(PXRMP_QED_OMP OFF CACHE INTERNAL "")
    endif()

    # Always disable tests
    set (PXRMP_QED_TEST OFF CACHE INTERNAL "")

    if(TOOLS_QED_PICSAR_src)
        add_subdirectory(
            ${TOOLS_QED_PICSAR_src}/multi_physics/QED
            _deps/localpicsar-build/
        )
    else()
        FetchContent_Declare(fetchedpicsar
            GIT_REPOSITORY ${TOOLS_QED_PICSAR_repo}
            GIT_TAG        ${TOOLS_QED_PICSAR_branch}
            BUILD_IN_SOURCE 0
        )
        FetchContent_GetProperties(fetchedpicsar)

        if(NOT fetchedpicsar_POPULATED)
            FetchContent_Populate(fetchedpicsar)
            add_subdirectory(
                ${fetchedpicsar_SOURCE_DIR}/multi_physics/QED
                ${fetchedpicsar_BINARY_DIR}
            )
        endif()
        if(NOT PXRMP_QED_GIT_VERSION)
            set(PXRMP_QED_GIT_VERSION "${TOOLS_QED_PICSAR_branch}" CACHE INTERNAL "")
        endif()

        # advanced fetch options
        mark_as_advanced(FETCHCONTENT_BASE_DIR)
        mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_QUIET)
        mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDPICSAR)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDPICSAR)
    endif()

    # advance PICSAR options not relevant for the user of the table generator/table reader
    mark_as_advanced(PXRMP_QED_TABLEGEN)
    mark_as_advanced(PXRMP_QED_OMP)
    mark_as_advanced(PXRMP_QED_TEST)
    mark_as_advanced(PXRMP_BOOST_TEST_DYN_LINK)
    mark_as_advanced(PXRMP_DPCPP_FIX)
    mark_as_advanced(PXRMP_KOKKOS_EXAMPLE)
    mark_as_advanced(PXRMP_QED_PYTHON_BINDINGS)
    mark_as_advanced(PXRMP_QED_TOOLS)

endmacro()

# local source-tree
set(TOOLS_QED_PICSAR_src ""
    CACHE PATH
    "Local path to PICSAR source directory (preferred if set)")

# Git fetcher
set(TOOLS_QED_PICSAR_repo "https://github.com/ECP-WarpX/picsar.git"
    CACHE STRING
    "Repository URI to pull and build PICSAR from if(TOOLS_QED_PICSAR_internal)")
set(TOOLS_QED_PICSAR_branch "2becfe0665590063fdfca10abaa05c4c1bd3ec22"
    CACHE STRING
    "commit for TOOLS_QED_PICSAR_repo if(TOOLS_QED_PICSAR_internal)")

find_picsar()

#____________________________________________________________________________________________________________

########
# Build QED lookup tables generator
########

set(arg_parser_cpp "Source/ArgParser/QedTablesArgParser.cpp")
set(arg_parser_h   "Source/ArgParser/QedTablesArgParser.H")

add_library( QedToolsArgParser ${arg_parser_cpp} ${arg_parser_h} )

if(TOOLS_QED_GENERATOR)

    set(qed_generator_name "qed_table_generator")
    set(qed_generator_src  "Source/QedTableGenerator.cpp")

    add_executable(${qed_generator_name} ${qed_generator_src})

    target_link_libraries(${qed_generator_name} PRIVATE PXRMP_QED)
    target_link_libraries(${qed_generator_name} PRIVATE QedToolsArgParser)

    target_compile_features(${qed_generator_name} PUBLIC cxx_std_17)
    set_target_properties(${qed_generator_name} PROPERTIES CXX_EXTENSIONS OFF)

    set_target_properties(${qed_generator_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

    message(STATUS "QED lookup tables generator: ${qed_generator_name} (${qed_generator_src})")

endif()

if(TOOLS_QED_GENERATOR)

    set(qed_reader_name "qed_table_reader")
    set(qed_reader_src  "Source/QedTableReader.cpp")

    add_executable(${qed_reader_name} ${qed_reader_src})

    target_link_libraries(${qed_reader_name} PRIVATE PXRMP_QED)
    target_link_libraries(${qed_reader_name} PRIVATE QedToolsArgParser)

    target_compile_features(${qed_reader_name} PUBLIC cxx_std_17)
    set_target_properties(${qed_reader_name} PROPERTIES CXX_EXTENSIONS OFF)

    set_target_properties(${qed_reader_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

    message(STATUS "QED lookup tables generator: ${qed_reader_name} (${qed_reader_src})")

endif()

#____________________________________________________________________________________________________________
