macro(find_amrex)
    if(WarpX_amrex_src)
        message(STATUS "Compiling local AMReX ...")
        message(STATUS "AMReX source path: ${WarpX_amrex_src}")
    elseif(WarpX_amrex_internal)
        message(STATUS "Downloading AMReX ...")
        message(STATUS "AMReX repository: ${WarpX_amrex_repo} (${WarpX_amrex_branch})")
        include(FetchContent)
    endif()
    if(WarpX_amrex_internal OR WarpX_amrex_src)
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        # see https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options
        if(WarpX_ASCENT)
            set(AMReX_ASCENT ON CACHE INTERNAL "")
            set(AMReX_CONDUIT ON CACHE INTERNAL "")
        endif()

        if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
            set(AMReX_ASSERTIONS ON CACHE BOOL "")
            # note: floating-point exceptions can slow down debug runs a lot
            set(AMReX_FPE ON CACHE BOOL "")
        else()
            set(AMReX_ASSERTIONS OFF CACHE BOOL "")
            set(AMReX_FPE OFF CACHE BOOL "")
        endif()

        if(WarpX_COMPUTE STREQUAL OMP)
            set(AMReX_GPU_BACKEND  "NONE" CACHE INTERNAL "")
            set(AMReX_OMP          ON     CACHE INTERNAL "")
        elseif(WarpX_COMPUTE STREQUAL NOACC)
            set(AMReX_GPU_BACKEND  "NONE" CACHE INTERNAL "")
            set(AMReX_OMP          OFF    CACHE INTERNAL "")
        else()
            set(AMReX_GPU_BACKEND  "${WarpX_COMPUTE}" CACHE INTERNAL "")
            set(AMReX_OMP          OFF    CACHE INTERNAL "")
        endif()

        if(WarpX_EB)
            set(AMReX_EB ON CACHE INTERNAL "")
        else()
            set(AMReX_EB OFF CACHE INTERNAL "")
        endif()

        if(WarpX_MPI)
            set(AMReX_MPI ON CACHE INTERNAL "")
            if(WarpX_MPI_THREAD_MULTIPLE)
                set(AMReX_MPI_THREAD_MULTIPLE ON CACHE INTERNAL "")
            else()
                set(AMReX_MPI_THREAD_MULTIPLE OFF CACHE INTERNAL "")
            endif()
        else()
            set(AMReX_MPI OFF CACHE INTERNAL "")
        endif()

        if(WarpX_PRECISION STREQUAL "DOUBLE")
            set(AMReX_PRECISION "DOUBLE" CACHE INTERNAL "")
        else()
            set(AMReX_PRECISION "SINGLE" CACHE INTERNAL "")
        endif()

        if(WarpX_PARTICLE_PRECISION STREQUAL "DOUBLE")
            set(AMReX_PARTICLES_PRECISION "DOUBLE" CACHE INTERNAL "")
        else()
            set(AMReX_PARTICLES_PRECISION "SINGLE" CACHE INTERNAL "")
        endif()

        if(WarpX_SENSEI)
            set(AMReX_SENSEI ON CACHE INTERNAL "")
        endif()

        if(DEFINED AMReX_BUILD_SHARED_LIBS)
            set(AMReX_INSTALL ${AMReX_BUILD_SHARED_LIBS} CACHE INTERNAL "")
        else()
            set(AMReX_INSTALL ${BUILD_SHARED_LIBS} CACHE INTERNAL "")
        endif()
        set(AMReX_AMRLEVEL OFF CACHE INTERNAL "")
        set(AMReX_ENABLE_TESTS OFF CACHE INTERNAL "")
        set(AMReX_FORTRAN OFF CACHE INTERNAL "")
        set(AMReX_FORTRAN_INTERFACES OFF CACHE INTERNAL "")
        set(AMReX_BUILD_TUTORIALS OFF CACHE INTERNAL "")
        set(AMReX_PARTICLES ON CACHE INTERNAL "")
        set(AMReX_PROBINIT OFF CACHE INTERNAL "")
        set(AMReX_TINY_PROFILE ON CACHE BOOL "")

        if(WarpX_ASCENT OR WarpX_SENSEI)
            set(AMReX_GPU_RDC ON CACHE BOOL "")
        else()
            # we don't need RDC and disabling it simplifies the build
            # complexity and potentially improves code optimization
            set(AMReX_GPU_RDC OFF CACHE BOOL "")
        endif()

        # shared libs, i.e. for Python bindings, need relocatable code
        #   openPMD: currently triggers shared libs (TODO)
        if(WarpX_LIB OR ABLASTR_POSITION_INDEPENDENT_CODE OR BUILD_SHARED_LIBS OR WarpX_OPENPMD)
            set(AMReX_PIC ON CACHE INTERNAL "")
        endif()

        # IPO/LTO
        if(WarpX_IPO)
            set(AMReX_IPO ON CACHE INTERNAL "")
            if(WarpX_COMPUTE STREQUAL CUDA)
                set(AMReX_CUDA_LTO ON CACHE BOOL "")
            endif()
        endif()

        # RZ is AMReX 2D
        set(WarpX_amrex_dim ${WarpX_DIMS})
        list(TRANSFORM WarpX_amrex_dim REPLACE RZ 2)
        list(REMOVE_DUPLICATES WarpX_amrex_dim)
        set(AMReX_SPACEDIM ${WarpX_amrex_dim} CACHE INTERNAL "")

        if(WarpX_amrex_src)
            list(APPEND CMAKE_MODULE_PATH "${WarpX_amrex_src}/Tools/CMake")
            if(WarpX_COMPUTE STREQUAL CUDA)
                enable_language(CUDA)
                # AMReX 21.06+ supports CUDA_ARCHITECTURES
            endif()
            add_subdirectory(${WarpX_amrex_src} _deps/localamrex-build/)
        else()
            FetchContent_Declare(fetchedamrex
                GIT_REPOSITORY ${WarpX_amrex_repo}
                GIT_TAG        ${WarpX_amrex_branch}
                BUILD_IN_SOURCE 0
            )
            FetchContent_GetProperties(fetchedamrex)

            if(NOT fetchedamrex_POPULATED)
                FetchContent_Populate(fetchedamrex)
                list(APPEND CMAKE_MODULE_PATH "${fetchedamrex_SOURCE_DIR}/Tools/CMake")
                if(WarpX_COMPUTE STREQUAL CUDA)
                    enable_language(CUDA)
                    # AMReX 21.06+ supports CUDA_ARCHITECTURES
                endif()
                add_subdirectory(${fetchedamrex_SOURCE_DIR} ${fetchedamrex_BINARY_DIR})
            endif()

            # advanced fetch options
            mark_as_advanced(FETCHCONTENT_BASE_DIR)
            mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_QUIET)
            mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDAMREX)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
            mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDAMREX)
        endif()

        # AMReX options not relevant to most WarpX users
        mark_as_advanced(AMREX_BUILD_DATETIME)
        mark_as_advanced(AMReX_DIFFERENT_COMPILER)
        mark_as_advanced(AMReX_ENABLE_TESTS)
        mark_as_advanced(AMReX_SPACEDIM)
        mark_as_advanced(AMReX_AMRDATA)
        mark_as_advanced(AMReX_BASE_PROFILE) # mutually exclusive to tiny profile
        mark_as_advanced(AMReX_CONDUIT)
        mark_as_advanced(AMReX_CUDA)
        mark_as_advanced(AMReX_CUDA_COMPILATION_TIMER)
        mark_as_advanced(AMReX_CUDA_ERROR_CAPTURE_THIS)
        mark_as_advanced(AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL)
        mark_as_advanced(AMReX_CUDA_FASTMATH)
        mark_as_advanced(AMReX_CUDA_KEEP_FILES)
        mark_as_advanced(AMReX_CUDA_LTO)
        mark_as_advanced(AMReX_CUDA_MAXREGCOUNT)
        mark_as_advanced(AMReX_CUDA_MAX_THREADS)
        mark_as_advanced(AMReX_CUDA_PTX_VERBOSE)
        mark_as_advanced(AMReX_CUDA_SHOW_CODELINES)
        mark_as_advanced(AMReX_CUDA_SHOW_LINENUMBERS)
        mark_as_advanced(AMReX_CUDA_WARN_CAPTURE_THIS)
        mark_as_advanced(AMReX_GPU_RDC)
        mark_as_advanced(AMReX_PARTICLES)
        mark_as_advanced(AMReX_PARTICLES_PRECISION)
        mark_as_advanced(AMReX_DPCPP)
        mark_as_advanced(AMReX_EB)
        mark_as_advanced(AMReX_FPE)
        mark_as_advanced(AMReX_FORTRAN)
        mark_as_advanced(AMReX_FORTRAN_INTERFACES)
        mark_as_advanced(AMReX_HDF5)  # we do HDF5 I/O (and more) via openPMD-api
        mark_as_advanced(AMReX_HIP)
        mark_as_advanced(AMReX_HYPRE)
        mark_as_advanced(AMReX_IPO)
        mark_as_advanced(AMReX_LINEAR_SOLVERS)
        mark_as_advanced(AMReX_MEM_PROFILE)
        mark_as_advanced(AMReX_MPI)
        mark_as_advanced(AMReX_MPI_THREAD_MULTIPLE)
        mark_as_advanced(AMReX_OMP)
        mark_as_advanced(AMReX_PROBINIT)
        mark_as_advanced(AMReX_PETSC)
        mark_as_advanced(AMReX_PIC)
        mark_as_advanced(AMReX_SENSEI)
        mark_as_advanced(AMReX_SUNDIALS)
        mark_as_advanced(AMReX_TINY_PROFILE)
        mark_as_advanced(AMReX_TP_PROFILE)
        mark_as_advanced(Boost_INCLUDE_DIR)
        mark_as_advanced(LIBNVTOOLSEXT)
        mark_as_advanced(PXRMP_QED_PYTHON_BINDINGS)
        mark_as_advanced(USE_XSDK_DEFAULTS)

        message(STATUS "AMReX: Using version '${AMREX_PKG_VERSION}' (${AMREX_GIT_VERSION})")
    else()
        message(STATUS "Searching for pre-installed AMReX ...")
        # https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#importing-amrex-into-your-cmake-project
        if(WarpX_ASCENT)
            set(COMPONENT_ASCENT ASCENT CONDUIT)
        else()
            set(COMPONENT_ASCENT)
        endif()

        set(WarpX_amrex_dim ${WarpX_DIMS})  # RZ is AMReX 2D
        list(TRANSFORM WarpX_amrex_dim REPLACE RZ 2)
        list(REMOVE_DUPLICATES WarpX_amrex_dim)
        set(COMPONENT_DIMS)
        foreach(D IN LISTS WarpX_amrex_dim)
            set(COMPONENT_DIMS ${COMPONENT_DIMS} ${D}D)
        endforeach()
        if(WarpX_EB)
            set(COMPONENT_EB EB)
        else()
            set(COMPONENT_EB)
        endif()
        if(WarpX_LIB)
            set(COMPONENT_PIC PIC)
        else()
            set(COMPONENT_PIC)
        endif()
        if(WarpX_SENSEI)
            set(COMPONENT_SENSEI SENSEI)
        else()
            set(COMPONENT_SENSEI)
        endif()
        set(COMPONENT_PRECISION ${WarpX_PRECISION} P${WarpX_PARTICLE_PRECISION})

        find_package(AMReX 23.05 CONFIG REQUIRED COMPONENTS ${COMPONENT_ASCENT} ${COMPONENT_DIMS} ${COMPONENT_EB} PARTICLES ${COMPONENT_PIC} ${COMPONENT_PRECISION} ${COMPONENT_SENSEI} TINYP LSOLVERS)
        message(STATUS "AMReX: Found version '${AMReX_VERSION}'")
    endif()
endmacro()

# local source-tree
set(WarpX_amrex_src ""
    CACHE PATH
    "Local path to AMReX source directory (preferred if set)")

# Git fetcher
set(WarpX_amrex_repo "https://github.com/AMReX-Codes/amrex.git"
    CACHE STRING
    "Repository URI to pull and build AMReX from if(WarpX_amrex_internal)")
set(WarpX_amrex_branch "342799dd78ec14b689507ede31d2acf2e7c3b217"
    CACHE STRING
    "Repository branch for WarpX_amrex_repo if(WarpX_amrex_internal)")

find_amrex()
