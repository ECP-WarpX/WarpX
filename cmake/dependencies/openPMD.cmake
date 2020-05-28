function(find_openpmd)
    if(WarpX_openpmd_internal)
        include(FetchContent)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        # see https://openpmd-api.readthedocs.io/en/0.11.1-alpha/dev/buildoptions.html
        set(openPMD_USE_MPI    ${ENABLE_MPI} CACHE INTERNAL "")
        set(openPMD_USE_PYTHON OFF           CACHE INTERNAL "")
        set(BUILD_CLI_TOOLS    OFF           CACHE INTERNAL "")  # FIXME
        set(BUILD_EXAMPLES     OFF           CACHE INTERNAL "")  # FIXME
        set(BUILD_TESTING      OFF           CACHE INTERNAL "")  # FIXME
        set(openPMD_INSTALL    ${BUILD_SHARED_LIBS} CACHE INTERNAL "")

        FetchContent_Declare(fetchedopenpmd
            GIT_REPOSITORY ${WarpX_openpmd_repo}
            GIT_TAG        ${WarpX_openpmd_branch}
            BUILD_IN_SOURCE 0
        )
        FetchContent_GetProperties(fetchedopenpmd)

        if(NOT fetchedopenpmd_POPULATED)
            FetchContent_Populate(fetchedopenpmd)
            add_subdirectory(${fetchedopenpmd_SOURCE_DIR} ${fetchedopenpmd_BINARY_DIR})
        endif()

        # advanced fetch options
        mark_as_advanced(FETCHCONTENT_BASE_DIR)
        mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_QUIET)
        mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDOPENPMD)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDOPENPMD)

        # openPMD options not relevant to WarpX users
        mark_as_advanced(openPMD_USE_INTERNAL_VARIANT)
        mark_as_advanced(openPMD_USE_INTERNAL_CATCH)
        mark_as_advanced(openPMD_USE_INTERNAL_PYBIND11)
        mark_as_advanced(openPMD_USE_INTERNAL_JSON)

        message(STATUS "openPMD-api: Using INTERNAL version '${WarpX_openpmd_branch}'")
    else()
        find_package(openPMD 0.11.1 CONFIG REQUIRED)  # FIXME: COMPONENTS MPI ...
        message(STATUS "openPMD-api: Found version '${openPMD_VERSION}'")
    endif()
endfunction()

if(WarpX_USE_OPENPMD)
    set(WarpX_HAVE_OPENPMD TRUE)
    find_openpmd()
endif()
