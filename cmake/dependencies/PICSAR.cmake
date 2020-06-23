function(find_picsar)
    if(WarpX_picsar_internal)
        include(FetchContent)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        # FIXME no option to control WarpX_QED_TABLE_GEN / Boost trigger yet

        FetchContent_Declare(fetchedpicsar
            GIT_REPOSITORY ${WarpX_picsar_repo}
            GIT_TAG        ${WarpX_picsar_branch}
            BUILD_IN_SOURCE 0
        )
        FetchContent_GetProperties(fetchedpicsar)

        if(NOT fetchedpicsar_POPULATED)
            FetchContent_Populate(fetchedpicsar)
            add_subdirectory(${fetchedpicsar_SOURCE_DIR}/src/multi_physics ${fetchedpicsar_BINARY_DIR})
        endif()

        # advanced fetch options
        mark_as_advanced(FETCHCONTENT_BASE_DIR)
        mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_QUIET)
        mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDPICSAR)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDPICSAR)

        # PICSAR options not relevant to WarpX users
        mark_as_advanced(DIM)
        mark_as_advanced(USE_XSDK_DEFAULTS)

        message(STATUS "PICSAR: Using INTERNAL version '${PICSAR_VERSION}'")
    else()
    # not supported by PICSAR
    #    find_package(PICSAR 20.05 CONFIG REQUIRED QED)
    #    message(STATUS "PICSAR: Found version '${PICSAR_VERSION}'")
        message(FATAL_ERROR "PICSAR: Cannot be used as externally installed library yet.")
    endif()
endfunction()

if(WarpX_QED)
    option(WarpX_picsar_internal   "Download & build PICSAR" ON)
    set(WarpX_picsar_repo "https://github.com/ECP-WarpX/picsar.git"
        CACHE STRING
        "Repository URI to pull and build PICSAR from if(WarpX_picsar_internal)")
    set(WarpX_picsar_branch "master"
        CACHE STRING
        "Repository branch for WarpX_picsar_repo if(WarpX_picsar_internal)")

    cmake_dependent_option(WarpX_QED_TABLE_GEN "generate QED lookup tables (requires boost)"
        ON "WarpX_QED" OFF)

    find_picsar()
endif()
