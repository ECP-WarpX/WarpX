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
    nprocs
    eb
    inputs
    analysis
    output
    # ${ARGV7} optional command line arguments for Python tests
)
    # cannot run MPI tests w/o MPI build
    if(NOT WarpX_MPI)
            message(WARNING "cannot run MPI tests without MPI build")
        return()
    endif()

    ## do not run no-MPI tests w/ MPI build
    #if(NOT mpi AND WarpX_MPI)
    #    return()
    #endif()

    # cannot run EB tests w/o EB build
    if(eb AND NOT WarpX_EB)
            message(WARNING "cannot run EB tests without EB build")
        return()
    endif()

    # do not run no-EB tests w/ EB build
    if(NOT eb AND WarpX_EB)
        return()
    endif()

    # cannot run tests with unsupported geometry
    if(NOT dims IN_LIST WarpX_DIMS)
        return()
    endif()

    # Python test?
    set(python OFF)
    if(inputs MATCHES ".*\.py$")
        set(python ON)
    endif()

    # cannot run Python tests w/o Python support
    if(python AND NOT WarpX_PYTHON)
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
    file(GLOB TEST_FILES "${INPUTS_DIR}/inputs*" "${INPUTS_DIR}/*.py")
    file(COPY ${TEST_FILES} DESTINATION ${THIS_WORKING_DIR})

    set(THIS_MPI_TEST_EXE
        ${MPIEXEC_EXECUTABLE}
        ${MPI_ALLOW_ROOT}
        ${MPIEXEC_NUMPROC_FLAG} ${nprocs}
        ${MPIEXEC_POSTFLAGS}
        ${MPIEXEC_PREFLAGS}
    )

    set(THIS_Python_EXE)
    if(python)
        set(THIS_Python_EXE ${Python_EXECUTABLE})
    endif()

    # test run
    if(python)
        set(command_line_args "")
        if(${ARGC} GREATER 7)
            set(command_line_args ${ARGV7})
        endif()
        add_test(
            NAME ${name}.run
            COMMAND
                ${THIS_MPI_TEST_EXE}
                ${THIS_Python_EXE}
                ${WarpX_SOURCE_DIR}/${inputs}
                ${command_line_args}
            WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        # FIXME Use helper function to handle Windows exceptions
        set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${PYTHONPATH}:${CMAKE_PYTHON_OUTPUT_DIRECTORY}")
    else()
        set(runtime_params
            "amrex.abort_on_unused_inputs = 1"
            "amrex.fpe_trap_invalid = 1"
            "amrex.fpe_trap_overflow = 1"
            "amrex.fpe_trap_zero = 1"
            "amrex.throw_exception = 1"
            "warpx.always_warn_immediately = 1"
            # FIXME should come before input file
            #"warpx.abort_on_warning_threshold = low"
            "warpx.do_dynamic_scheduling = 0"
            "warpx.serialize_initial_conditions = 1"
        )
        add_test(
            NAME ${name}.run
            COMMAND
                ${THIS_MPI_TEST_EXE}
                $<TARGET_FILE:app_${SD}>
                ${WarpX_SOURCE_DIR}/${inputs}
                ${runtime_params}
            WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
    endif()

    set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")

    # test analysis
    if(analysis)
        add_test(
            NAME ${name}.analysis
            COMMAND
                ${WarpX_SOURCE_DIR}/${analysis}
                ${output}
            WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        # test analysis depends on test run
        set_property(TEST ${name}.analysis APPEND PROPERTY DEPENDS "${name}.run")
        # FIXME Use helper function to handle Windows exceptions
        set(PYTHONPATH "${PYTHONPATH}:${CMAKE_PYTHON_OUTPUT_DIRECTORY}")
        # add paths for custom Python modules
        set(PYTHONPATH "${PYTHONPATH}:${WarpX_SOURCE_DIR}/Regression/Checksum")
        set(PYTHONPATH "${PYTHONPATH}:${WarpX_SOURCE_DIR}/Regression/PostProcessingUtils")
        set(PYTHONPATH "${PYTHONPATH}:${WarpX_SOURCE_DIR}/Tools/Parser")
        set(PYTHONPATH "${PYTHONPATH}:${WarpX_SOURCE_DIR}/Tools/PostProcessing")
        set_property(TEST ${name}.analysis APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${PYTHONPATH}")
    endif()

endfunction()

# Add tests (alphabetical order) ##############################################
#

add_subdirectory(Tests)
add_subdirectory(Physics_applications)
