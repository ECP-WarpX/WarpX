# Configuration ###############################################################
#
if(WarpX_MPI)
    # OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
    if("$ENV{USER}" STREQUAL "root")
        # calling even --help as root will abort and warn on stderr
        execute_process(
            COMMAND ${MPIEXEC_EXECUTABLE} --help
            ERROR_VARIABLE MPIEXEC_HELP_TEXT
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(${MPIEXEC_HELP_TEXT} MATCHES "^.*allow-run-as-root.*$")
            set(MPI_ALLOW_ROOT --allow-run-as-root)
        endif()
    endif()
endif()

# Add a WarpX test set (with sub-tests)
#
function(add_warpx_test
    name
    dims
    mpi
    mpi_np
    inputs
    analysis
    output
)
    # cannot run MPI tests w/o MPI support
    if(mpi AND NOT WarpX_MPI)
        return()
    endif()

    if(NOT dims IN_LIST WarpX_DIMS)
        return()
    endif()

    # set dimension suffix
    warpx_set_suffix_dims(SD ${dims})

    # make a unique run directory
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})
    set(THIS_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})

    # copy local input files to run directory
    get_filename_component(INPUTS_PATH ${inputs} REALPATH BASE_DIR "${CMAKE_BINARYDIR}")
    get_filename_component(INPUTS_DIR ${INPUTS_PATH} DIRECTORY)
    file(GLOB TEST_FILES "${INPUTS_DIR}/inputs*")
    file(COPY ${TEST_FILES} DESTINATION ${THIS_WORKING_DIR})

    set(THIS_MPI_TEST_EXE)
    if(mpi)
        set(THIS_MPI_TEST_EXE
            ${MPIEXEC_EXECUTABLE}
            ${MPI_ALLOW_ROOT}
            ${MPIEXEC_NUMPROC_FLAG} ${mpi_np}
            ${MPIEXEC_POSTFLAGS}
            ${MPIEXEC_PREFLAGS}
        )
    endif()

    # FIXME
    set(runtime_params
        "amrex.abort_on_unused_inputs = 1"
        "amrex.fpe_trap_invalid = 1"
        "amrex.fpe_trap_overflow = 1"
        "amrex.fpe_trap_zero = 1"
        "amrex.throw_exception = 1"
        "warpx.always_warn_immediately = 1"
        #"warpx.abort_on_warning_threshold = low"
        "warpx.do_dynamic_scheduling = 0"
        "warpx.serialize_initial_conditions = 1"
    )

    # test run
    add_test(
        NAME ${name}.run
        COMMAND
            ${THIS_MPI_TEST_EXE}
            $<TARGET_FILE:app_${SD}>
            ${WarpX_SOURCE_DIR}/${inputs}
            ${runtime_params}
        WORKING_DIRECTORY ${THIS_WORKING_DIR}
    )

    if(mpi)
        set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
    endif()

    # test analysis
    if(analysis)
        add_test(
            NAME ${name}.analysis
            COMMAND
                ${WarpX_SOURCE_DIR}/${analysis}
                ${output}
            WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        set_property(TEST ${name}.analysis APPEND PROPERTY DEPENDS "${name}.run")
    endif()

endfunction()

# Add tests (alphabetical order) ##############################################
#

add_subdirectory(Tests)
add_subdirectory(Physics_applications)
