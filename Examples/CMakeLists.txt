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
    inputs_base
    inputs_test
    analysis
    output
)
    # cannot run MPI tests w/o MPI support
    if(mpi AND NOT WarpX_MPI)
        return()
    endif()

    # make a unique run directory
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})
    set(THIS_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})

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

    # test run
    add_test(
        NAME ${name}.run
        COMMAND
            ${THIS_MPI_TEST_EXE}
            $<TARGET_FILE:app_${dims}d>
            ${WarpX_SOURCE_DIR}/${inputs_base}
            amrex.abort_on_unused_inputs = 1
            amrex.fpe_trap_invalid = 1
            amrex.fpe_trap_overflow = 1
            amrex.fpe_trap_zero = 1
            amrex.throw_exception = 1
            warpx.always_warn_immediately = 1
            warpx.abort_on_warning_threshold = low
            warpx.do_dynamic_scheduling = 0
            warpx.serialize_initial_conditions = 1
            # FILE parameter must be at the end,
            # to possibly overwrite previous inputs
            FILE = ${WarpX_SOURCE_DIR}/${inputs_test}
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

# Langmuir_multi_2d_MR ########################################################
#
if(2 IN_LIST WarpX_DIMS)
    set(name Langmuir_multi_2d_MR)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 2)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()

# Langmuir_multi_2d_MR_anisotropic ############################################
#
if(2 IN_LIST WarpX_DIMS)
    set(name Langmuir_multi_2d_MR_anisotropic)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 2)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_anisotropic)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()

# Langmuir_multi_2d_MR_momentum_conserving ####################################
#
if(2 IN_LIST WarpX_DIMS)
    set(name Langmuir_multi_2d_MR_momentum_conserving)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 2)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_momentum_conserving)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()

# Langmuir_multi_2d_MR_psatd ##################################################
#
if(2 IN_LIST WarpX_DIMS AND WarpX_FFT)
    set(name Langmuir_multi_2d_MR_psatd)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 2)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_psatd)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()

# Langmuir_multi_2d_nodal #####################################################
#
if(2 IN_LIST WarpX_DIMS)
    set(name Langmuir_multi_2d_nodal)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 2)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_nodal)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()

# Langmuir_multi_2d_psatd #####################################################
#
if(2 IN_LIST WarpX_DIMS AND WarpX_FFT)
    set(name Langmuir_multi_2d_psatd)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 2)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()

# Langmuir_multi_2d_psatd_current_correction ##################################
#
if(2 IN_LIST WarpX_DIMS AND WarpX_FFT)
    set(name Langmuir_multi_2d_psatd_current_correction)
    set(dims 2)
    set(mpi ON)
    set(mpi_np 1)
    set(inputs_base Examples/Tests/langmuir/inputs_2d_base)
    set(inputs_test Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_current_correction)
    set(analysis Examples/Tests/langmuir/analysis_2d.py)
    set(output diags/diag1000080)
    add_warpx_test(
        ${name}
        ${dims}
        ${mpi}
        ${mpi_np}
        ${inputs_base}
        ${inputs_test}
        ${analysis}
        ${output}
    )
endif()
