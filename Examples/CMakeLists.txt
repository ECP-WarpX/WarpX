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
# name: unique name of this test
# dims: 1,2,RZ,3
# nprocs: 1 or 2 (maybe refactor later on to just depend on WarpX_MPI)
# eb: needs EB support? (temporary until handled as runtime parameter)
# inputs: inputs file or PICMI script, WarpX_MPI decides w/ or w/o MPI
# analysis: analysis script, always run without MPI
# output: output file (dir) to analyze
# dependency: name of base test that must run first
#
function(add_warpx_test
    name
    dims
    nprocs
    eb
    inputs
    analysis
    output
    dependency
)
    # cannot run MPI tests w/o MPI build
    if(nprocs GREATER_EQUAL 2 AND NOT WarpX_MPI)
        message(WARNING "${name}: cannot run MPI tests without MPI build")
        return()
    endif()

    # cannot run EB tests w/o EB build
    if(eb AND NOT WarpX_EB)
        message(WARNING "${name}: cannot run EB tests without EB build")
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

    # cannot run tests with unfulfilled dependencies
    if(dependency AND NOT TEST ${dependency}.run)
        return()
    endif()

    # set dimension suffix
    warpx_set_suffix_dims(SD ${dims})

    # make a unique run directory
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})
    set(THIS_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})

    get_filename_component(INPUTS_ALL ${WarpX_SOURCE_DIR}/${inputs} REALPATH)

    # separate input file (full path) and optional command-line arguments
    separate_arguments(INPUTS_FILE_ARGS UNIX_COMMAND PROGRAM "${INPUTS_ALL}")

    # get input file (full path)
    list(GET INPUTS_FILE_ARGS 0 INPUTS_FILE)

    # get optional command-line arguments
    list(SUBLIST INPUTS_FILE_ARGS 1 -1 INPUTS_ARGS)
    list(JOIN INPUTS_FILE_ARGS " " INPUTS_FILE_ARGS)

    # get input file's directory
    get_filename_component(INPUTS_DIR ${INPUTS_FILE} DIRECTORY)

    # copy base input files to run directory
    file(GLOB INPUTS_BASE "${INPUTS_DIR}/inputs_base_${SD}*")
    file(COPY ${INPUTS_BASE} DESTINATION ${THIS_WORKING_DIR})

    # copy test input file to run directory
    file(COPY ${INPUTS_FILE} DESTINATION ${THIS_WORKING_DIR})

    # for restart tests, copy input file of original test (without "_restart" suffix)
    string(REPLACE "_restart" "" INPUTS_FILE_NORESTART ${INPUTS_FILE})
    file(COPY ${INPUTS_FILE_NORESTART} DESTINATION ${THIS_WORKING_DIR})

    # Python test?
    set(python OFF)
    if(${INPUTS_FILE} MATCHES ".*\.py$")
        set(python ON)
    endif()

    # cannot run Python tests w/o Python support
    if(python AND NOT WarpX_PYTHON)
        return()
    endif()

    # set MPI executable
    set(THIS_MPI_TEST_EXE
        ${MPIEXEC_EXECUTABLE}
        ${MPI_ALLOW_ROOT}
        ${MPIEXEC_NUMPROC_FLAG} ${nprocs}
        ${MPIEXEC_POSTFLAGS}
        ${MPIEXEC_PREFLAGS}
    )

    # set Python executable
    set(THIS_Python_EXE)
    if(python)
        set(THIS_Python_EXE ${Python_EXECUTABLE})
    endif()

    # test run
    if(python)
        # for argparse, do not pass command-line arguments as one quoted string
        separate_arguments(INPUTS_ARGS UNIX_COMMAND "${INPUTS_ARGS}")
        add_test(
            NAME ${name}.run
            COMMAND
                ${THIS_MPI_TEST_EXE}
                ${THIS_Python_EXE}
                ${INPUTS_FILE}
                ${INPUTS_ARGS}
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
                ${INPUTS_FILE}
                ${runtime_params}
                ${INPUTS_ARGS}
            WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
    endif()

    # run all tests with 1 OpenMP thread by default
    set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")

    set(THIS_Python_SCRIPT_EXE)
    if(python OR WIN32)
        set(THIS_Python_SCRIPT_EXE ${Python_EXECUTABLE})
    endif()

    # test analysis
    if(analysis)
        add_test(
            NAME ${name}.analysis
            COMMAND
                ${THIS_Python_SCRIPT_EXE} ${WarpX_SOURCE_DIR}/${analysis}
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

    # CI: remove test directory after run
    if(WarpX_TEST_CLEANUP)
        add_test(
            NAME ${name}.cleanup
            COMMAND ${CMAKE_COMMAND} -E rm -rf ${THIS_WORKING_DIR}
        )
        # test cleanup depends on test run
        set_property(TEST ${name}.cleanup APPEND PROPERTY DEPENDS "${name}.run")
        if(analysis)
            # test cleanup depends on test analysis
            set_property(TEST ${name}.cleanup APPEND PROPERTY DEPENDS "${name}.analysis")
        endif()
    endif()

    # Do we depend on another test?
    if(dependency)
        # current test depends on dependency test run (and analysis)
        set_property(TEST ${name}.run APPEND PROPERTY DEPENDS "${dependency}.run")
        if(analysis)
            set_property(TEST ${name}.run APPEND PROPERTY DEPENDS "${dependency}.analysis")
        endif()
        if(WarpX_TEST_CLEANUP)
            # do not clean up dependency test before current test is completed
            set_property(TEST ${dependency}.cleanup APPEND PROPERTY DEPENDS "${name}.cleanup")
        endif()
    endif()
endfunction()

# Add tests (alphabetical order) ##############################################
#

add_subdirectory(Tests)
add_subdirectory(Physics_applications)
