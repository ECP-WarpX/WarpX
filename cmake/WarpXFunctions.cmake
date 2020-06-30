# find the CCache tool and use it if found
#
macro(set_ccache)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
        if(ENABLE_CUDA)
            set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
        endif()
    endif()
    mark_as_advanced(CCACHE_PROGRAM)
endmacro()


# set names and paths of temporary build directories
# the defaults in CMake are sub-ideal for historic reasons, lets make them more
# Unix-ish and portable.
#
macro(set_default_build_dirs)
    if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                CACHE PATH "Build directory for archives")
        mark_as_advanced(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    endif()
    if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                CACHE PATH "Build directory for libraries")
        mark_as_advanced(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    endif()
    if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                CACHE PATH "Build directory for binaries")
        mark_as_advanced(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    endif()
endmacro()


# set names and paths of install directories
# the defaults in CMake are sub-ideal for historic reasons, lets make them more
# Unix-ish and portable.
#
macro(set_default_install_dirs)
    include(GNUInstallDirs)
    if(NOT CMAKE_INSTALL_CMAKEDIR)
        set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/WarpX"
                CACHE PATH "CMake config package location for installed targets")
        if(WIN32)
            set(CMAKE_INSTALL_LIBDIR Lib
                    CACHE PATH "Object code libraries")
            set_property(CACHE CMAKE_INSTALL_CMAKEDIR PROPERTY VALUE "cmake")
        endif()
        mark_as_advanced(CMAKE_INSTALL_CMAKEDIR)
    endif()
endmacro()


# change the default CMAKE_BUILD_TYPE
# the default in CMake is Debug for historic reasons
#
macro(set_default_build_type default_build_type)
    set(CMAKE_CONFIGURATION_TYPES "Release;Debug;MinSizeRel;RelWithDebInfo")
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE ${default_build_type}
            CACHE STRING
            "Choose the build type, e.g. Release, Debug, or RelWithDebInfo." FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE
            PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES})
    endif()
endmacro()

# Set CXX
# Note: this is a bit legacy and one should use CMake TOOLCHAINS instead.
#
macro(set_cxx_warnings)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # list(APPEND CMAKE_CXX_FLAGS "-fsanitize=address") # address, memory, undefined
        # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
        # set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
        # set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fsanitize=address")

        # note: might still need a
        #   export LD_PRELOAD=libclang_rt.asan.so
        # or on Debian 9 with Clang 6.0
        #   export LD_PRELOAD=/usr/lib/llvm-6.0/lib/clang/6.0.0/lib/linux/libclang_rt.asan-x86_64.so:
        #                     /usr/lib/llvm-6.0/lib/clang/6.0.0/lib/linux/libclang_rt.ubsan_minimal-x86_64.so
        # at runtime when used with symbol-hidden code (e.g. pybind11 module)

        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wshadow -Woverloaded-virtual -Wextra-semi -Wunreachable-code")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wshadow -Woverloaded-virtual -Wunreachable-code")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        # Warning C4503: "decorated name length exceeded, name was truncated"
        # Symbols longer than 4096 chars are truncated (and hashed instead)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4503")
        # Yes, you should build against the same C++ runtime and with same
        # configuration (Debug/Release). MSVC does inconvenient choices for their
        # developers, so be it. (Our Windows-users use conda-forge builds, which
        # are consistent.)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4251")
    endif ()
endmacro()

# Take an <imported_target> and expose it as INTERFACE target with
# WarpX::thirdparty::<propagated_name> naming and SYSTEM includes.
#
function(make_third_party_includes_system imported_target propagated_name)
    add_library(WarpX::thirdparty::${propagated_name} INTERFACE IMPORTED)
    target_link_libraries(WarpX::thirdparty::${propagated_name} INTERFACE ${imported_target})
    get_target_property(ALL_INCLUDES ${imported_target} INCLUDE_DIRECTORIES)
    set_target_properties(WarpX::thirdparty::${propagated_name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
    target_include_directories(WarpX::thirdparty::${propagated_name} SYSTEM INTERFACE ${ALL_INCLUDES})
endfunction()


# Set an MPI_TEST_EXE variable for test runs which runs num_ranks
# ranks. On some systems, you might need to use the a specific
# mpiexec wrapper, e.g. on Summit (ORNL) pass the hint
# -DMPIEXEC_EXECUTABLE=$(which jsrun) to run ctest.
#
function(configure_mpiexec num_ranks)
    # OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
    if("$ENV{USER}" STREQUAL "root")
        # calling even --help as root will abort and warn on stderr
        execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --help
            ERROR_VARIABLE MPIEXEC_HELP_TEXT
            OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(${MPIEXEC_HELP_TEXT} MATCHES "^.*allow-run-as-root.*$")
                set(MPI_ALLOW_ROOT --allow-run-as-root)
            endif()
    endif()
    set(MPI_TEST_EXE
        ${MPIEXEC_EXECUTABLE}
        ${MPI_ALLOW_ROOT}
        ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
        PARENT_SCOPE
    )
endfunction()


# Prints a summary of WarpX options at the end of the CMake configuration
#
function(warpx_print_summary)
    message("")
    message("WarpX build configuration:")
    message("  Version: ${WarpX_VERSION}")
    message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
                            "${CMAKE_CXX_COMPILER_VERSION} "
                            "${CMAKE_CXX_COMPILER_WRAPPER}")
    message("    ${CMAKE_CXX_COMPILER}")
    message("")
    message("  Installation prefix: ${CMAKE_INSTALL_PREFIX}")
    message("        bin: ${CMAKE_INSTALL_BINDIR}")
    message("        lib: ${CMAKE_INSTALL_LIBDIR}")
    message("    include: ${CMAKE_INSTALL_INCLUDEDIR}")
    message("      cmake: ${CMAKE_INSTALL_CMAKEDIR}")
    if(WarpX_HAVE_PYTHON)
        message("     python: ${CMAKE_INSTALL_PYTHONDIR}")
    endif()
    message("")
    message("  Build type: ${CMAKE_BUILD_TYPE}")
    #if(BUILD_SHARED_LIBS)
    #    message("  Library: shared")
    #else()
    #    message("  Library: static")
    #endif()
    message("  Testing: ${BUILD_TESTING}")
    #message("  Invasive Tests: ${WarpX_USE_INVASIVE_TESTS}")
    #message("  Internal VERIFY: ${WarpX_USE_VERIFY}")
    message("  Build options:")
    message("    ASCENT: ${WarpX_ASCENT}")
    message("    COMPUTE: ${WarpX_COMPUTE}")
    message("    DIMS: ${WarpX_DIMS}")
    message("    MPI: ${WarpX_MPI}")
    message("    PSATD: ${WarpX_PSATD}")
    message("    PRECISION: ${WarpX_PRECISION}")
    message("    OPENPMD: ${WarpX_OPENPMD}")
    message("    QED: ${WarpX_QED}")
    message("")
endfunction()
