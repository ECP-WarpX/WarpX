function(find_picsar)
    if(WarpX_picsar_src)
        message(STATUS "Compiling local PICSAR ...")
        message(STATUS "PICSAR source path: ${WarpX_picsar_src}")
    elseif(WarpX_picsar_internal)
        message(STATUS "Downloading PICSAR ...")
        message(STATUS "PICSAR repository: ${WarpX_picsar_repo} (${WarpX_picsar_branch})")
        include(FetchContent)
    endif()
    if(WarpX_picsar_internal OR WarpX_picsar_src)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        # Enable or disable QED lookup tables generation

        # If table generation is enabled, enable or disable
        # openMP support depending on WarpX_COMPUTE
        if(WarpX_QED_TABLE_GEN)
            set(PXRMP_QED_TABLEGEN ON CACHE INTERNAL "")
            if(WarpX_COMPUTE STREQUAL OMP)
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

        if(WarpX_COMPUTE STREQUAL SYCL)
            set (PXRMP_DPCPP_FIX ON CACHE INTERNAL "")
        endif()

        if(WarpX_picsar_src)
            add_subdirectory(
                ${WarpX_picsar_src}/multi_physics/QED
                _deps/localpicsar-build/
            )
            get_source_version(PXRMP_QED ${WarpX_picsar_src})
        else()
            FetchContent_Declare(fetchedpicsar
                GIT_REPOSITORY ${WarpX_picsar_repo}
                GIT_TAG        ${WarpX_picsar_branch}
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
            get_source_version(PXRMP_QED ${fetchedpicsar_SOURCE_DIR})
            if(NOT PXRMP_QED_GIT_VERSION)
                set(PXRMP_QED_GIT_VERSION "${WarpX_picsar_branch}" CACHE INTERNAL "")
            endif()

            # advanced fetch options
            mark_as_advanced(FETCHCONTENT_BASE_DIR)
            mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_QUIET)
            mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDPICSAR)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDPICSAR)
        endif()

        # PICSAR options not relevant to WarpX users
        mark_as_advanced(DIM)
        mark_as_advanced(USE_XSDK_DEFAULTS)
        mark_as_advanced(PXRMP_QED_TABLEGEN)
        mark_as_advanced(PXRMP_QED_OMP)
        mark_as_advanced(PXRMP_QED_TEST)
        mark_as_advanced(PXRMP_BOOST_TEST_DYN_LINK)
        mark_as_advanced(PXRMP_DPCPP_FIX)


        # PICSAR_VERSION: not yet defined
        #message(STATUS "PICSAR: Using version '${PICSAR_VERSION}'")
    else()
        # not supported by PICSAR (yet)
        #find_package(PICSAR 22.01 CONFIG REQUIRED QED)
        #message(STATUS "PICSAR: Found version '${PICSAR_VERSION}'")
        message(FATAL_ERROR "PICSAR: Cannot be used as externally installed "
            "library yet. "
            "Note that you can point to an external source tree with "
            "-DWarpX_picsar_src=<path>"
        )
    endif()
endfunction()

if(WarpX_QED)
    # local source-tree
    set(WarpX_picsar_src ""
        CACHE PATH
        "Local path to PICSAR source directory (preferred if set)")

    # Git fetcher
    option(WarpX_picsar_internal   "Download & build PICSAR" ON)
    set(WarpX_picsar_repo "https://github.com/ECP-WarpX/picsar.git"
        CACHE STRING
        "Repository URI to pull and build PICSAR from if(WarpX_picsar_internal)")
    set(WarpX_picsar_branch "7b5449f92a4b30a095cc4a67f0a8b1fc69680e15"
        CACHE STRING
        "Repository branch for WarpX_picsar_repo if(WarpX_picsar_internal)")

    cmake_dependent_option(WarpX_QED_TABLE_GEN "generate QED lookup tables (requires boost)"
        ON "WarpX_QED" OFF)

    find_picsar()
endif()
