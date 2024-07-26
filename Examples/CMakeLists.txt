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
    analysis
    output
    #inputs_test (optional)
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

    set(runtime_params
        "amrex.abort_on_unused_inputs = 1"
        "amrex.fpe_trap_invalid = 1"
        "amrex.fpe_trap_overflow = 1"
        "amrex.fpe_trap_zero = 1"
        "amrex.throw_exception = 1"
        "warpx.always_warn_immediately = 1"
        "warpx.abort_on_warning_threshold = low"
        "warpx.do_dynamic_scheduling = 0"
        "warpx.serialize_initial_conditions = 1"
    )

    # FILE parameter must be at the end,
    # to possibly overwrite previous inputs
    # FIXME Avoid magic numbers
    if(${ARGC} GREATER 7)
        list(APPEND runtime_params "FILE = ${WarpX_SOURCE_DIR}/${ARGV7}")
    endif()

    # test run
    add_test(
        NAME ${name}.run
        COMMAND
            ${THIS_MPI_TEST_EXE}
            $<TARGET_FILE:app_${SD}>
            ${WarpX_SOURCE_DIR}/${inputs_base}
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

# averaged_galilean_2d_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_2d_psatd # name
        2  # dims
        ON # mpi
        1  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_2d_base # inputs_base
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
        Examples/Tests/nci_psatd_stability/inputs_avg_2d_test_averaged_galilean_2d_psatd # inputs_test
    )
endif()

# averaged_galilean_2d_psatd_hybrid ###########################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_2d_psatd_hybrid # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_2d_base # inputs_base
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
        Examples/Tests/nci_psatd_stability/inputs_avg_2d_test_averaged_galilean_2d_psatd_hybrid # inputs_test
    )
endif()

# averaged_galilean_3d_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_3d_psatd # name
        3  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_base # inputs_base
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000160 # output
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_test_averaged_galilean_3d_psatd # inputs_test
    )
endif()

# averaged_galilean_3d_psatd_hybrid ###########################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_3d_psatd_hybrid # name
        3  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_base # inputs_base
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000160 # output
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_test_averaged_galilean_3d_psatd_hybrid # inputs_test
    )
endif()

# background_mcc ##############################################################
#
add_warpx_test(
    background_mcc # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Physics_applications/capacitive_discharge/inputs_2d_base # inputs_base
    Examples/analysis_default_regression.py # analysis
    diags/diag1000050 # output
    Examples/Physics_applications/capacitive_discharge/inputs_2d_test_background_mcc # inputs_test
)

# background_mcc_dp_psp #######################################################
#
# FIXME
#add_warpx_test(
#    background_mcc_dp_psp # name
#    2  # dims
#    ON # mpi
#    2  # mpi_np
#    Examples/Physics_applications/capacitive_discharge/inputs_2d_base # inputs_base
#    Examples/analysis_default_regression.py # analysis
#    diags/diag1000050 # output
#    Examples/Physics_applications/capacitive_discharge/inputs_2d_test_background_mcc_dp_psp # inputs_test
#)

# bilinear_filter #############################################################
#
add_warpx_test(
    bilinear_filter # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/single_particle/inputs_2d_base # inputs_base
    Examples/Tests/single_particle/analysis_bilinear_filter.py # analysis
    diags/diag1000001 # output
    Examples/Tests/single_particle/inputs_2d_test_bilinear_filter # inputs_test
)

# BTD_rz ######################################################################
#
add_warpx_test(
    BTD_rz # name
    RZ # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/btd_rz/inputs_rz_z_boosted_BTD # inputs_base
    Examples/Tests/btd_rz/analysis_BTD_laser_antenna.py # analysis
    diags/diag1000289 # output
)

# Langmuir_multi_2d_MR ########################################################
#
add_warpx_test(
    Langmuir_multi_2d_MR # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/langmuir/inputs_2d_base # inputs_base
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR # inputs_test
)

# Langmuir_multi_2d_MR_anisotropic ############################################
#
add_warpx_test(
    Langmuir_multi_2d_MR_anisotropic # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/langmuir/inputs_2d_base # inputs_base
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_anisotropic # inputs_test
)

# Langmuir_multi_2d_MR_momentum_conserving ####################################
#
add_warpx_test(
    Langmuir_multi_2d_MR_momentum_conserving # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/langmuir/inputs_2d_base # inputs_base
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_momentum_conserving # inputs_test
)

# Langmuir_multi_2d_MR_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_MR_psatd # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_psatd # inputs_test
    )
endif()

# Langmuir_multi_2d_nodal #####################################################
#
add_warpx_test(
    Langmuir_multi_2d_nodal # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/langmuir/inputs_2d_base # inputs_base
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_nodal # inputs_test
)

# Langmuir_multi_2d_psatd #####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_current_correction ##################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_current_correction # name
        2  # dims
        ON # mpi
        1  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_current_correction # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_current_correction_nodal ############################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_current_correction_nodal # name
        2  # dims
        ON # mpi
        1  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_current_correction_nodal # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_momentum_conserving #################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_momentum_conserving # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_momentum_conserving # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_multiJ ##############################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_multiJ # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_multiJ # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_multiJ_nodal ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_multiJ_nodal # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_multiJ_nodal # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_nodal ###############################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_nodal # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_nodal # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_Vay_deposition ######################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_Vay_deposition # name
        2  # dims
        ON # mpi
        1  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_Vay_deposition # inputs_test
    )
endif()

# Langmuir_multi_2d_psatd_Vay_deposition_particle_shape_4 #####################
#
# FIXME
#if(WarpX_FFT)
#    add_warpx_test(
#        Langmuir_multi_2d_psatd_Vay_deposition_particle_shape_4 # name
#        2  # dims
#        ON # mpi
#        1  # mpi_np
#        Examples/Tests/langmuir/inputs_2d_base # inputs_base
#        Examples/Tests/langmuir/analysis_2d.py # analysis
#        diags/diag1000080 # output
#        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_Vay_deposition_particle_shape_4 # inputs_test
#    )
#endif()

# Langmuir_multi_2d_psatd_Vay_deposition_nodal ################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_Vay_deposition_nodal # name
        2  # dims
        ON # mpi
        1  # mpi_np
        Examples/Tests/langmuir/inputs_2d_base # inputs_base
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_Vay_deposition_nodal # inputs_test
    )
endif()
