# Add tests (alphabetical order) ##############################################
#

# Langmuir_multi ##############################################################
#
add_warpx_test(
    Langmuir_multi # name
    3 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi # inputs
    Examples/Tests/langmuir/analysis_3d.py # analysis
    diags/diag1000040 # output
)

# Langmuir_multi_1d ###########################################################
#
add_warpx_test(
    Langmuir_multi_1d # name
    1 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_1d_test_Langmuir_multi_1d # inputs
    Examples/Tests/langmuir/analysis_1d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_multi_2d_MR ########################################################
#
add_warpx_test(
    Langmuir_multi_2d_MR # name
    2 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_multi_2d_MR_anisotropic ############################################
#
add_warpx_test(
    Langmuir_multi_2d_MR_anisotropic # name
    2 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_anisotropic # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_multi_2d_MR_momentum_conserving ####################################
#
add_warpx_test(
    Langmuir_multi_2d_MR_momentum_conserving # name
    2 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_momentum_conserving # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_multi_2d_MR_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_MR_psatd # name
        2 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_MR_psatd # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_nodal #####################################################
#
add_warpx_test(
    Langmuir_multi_2d_nodal # name
    2 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_nodal # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_multi_2d_picmi #####################################################
#
add_warpx_test(
    Langmuir_multi_2d_picmi # name
    2 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

# Langmuir_multi_2d_psatd #####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd # name
        2 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_current_correction ##################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_current_correction # name
        2 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_current_correction_nodal ############################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_current_correction_nodal # name
        2 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_current_correction_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_momentum_conserving #################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_momentum_conserving # name
        2 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_momentum_conserving # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_multiJ ##############################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_multiJ # name
        2 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_multiJ_nodal ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_multiJ_nodal # name
        2 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_multiJ_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_nodal ###############################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_nodal # name
        2 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_2d_psatd_Vay_deposition ######################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_Vay_deposition # name
        2 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_Vay_deposition # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# FIXME
# Langmuir_multi_2d_psatd_Vay_deposition_particle_shape_4 #####################
#
#if(WarpX_FFT)
#    add_warpx_test(
#        Langmuir_multi_2d_psatd_Vay_deposition_particle_shape_4 # name
#        2 # dims
#        1 # nprocs
#        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_Vay_deposition_particle_shape_4 # inputs
#        Examples/Tests/langmuir/analysis_2d.py # analysis
#        diags/diag1000080 # output
#    )
#endif()

# Langmuir_multi_2d_psatd_Vay_deposition_nodal ################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_2d_psatd_Vay_deposition_nodal # name
        2 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_2d_test_Langmuir_multi_2d_psatd_Vay_deposition_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_nodal ########################################################
#
add_warpx_test(
    Langmuir_multi_nodal # name
    3 # dims
    2 # nprocs
    Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_nodal # inputs
    Examples/Tests/langmuir/analysis_3d.py # analysis
    diags/diag1000040 # output
)

# Langmuir_multi_picmi ########################################################
#
add_warpx_test(
    Langmuir_multi_picmi # name
    3 # dims
    1 # nprocs
    Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

# Langmuir_multi_psatd ########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd # name
        3 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_current_correction #####################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_current_correction # name
        3 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_current_correction_nodal ###############################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_current_correction_nodal # name
        3 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_current_correction_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_div_cleaning ###########################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_div_cleaning # name
        3 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_div_cleaning # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_momentum_conserving ####################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_momentum_conserving # name
        3 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_momentum_conserving # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_multiJ #################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_multiJ # name
        3 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_multiJ_nodal ###########################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_multiJ_nodal # name
        3 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_multiJ_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_nodal ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_nodal # name
        3 # dims
        2 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_Vay_deposition #########################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_Vay_deposition # name
        3 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_Vay_deposition # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_psatd_Vay_deposition_nodal ###################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_psatd_Vay_deposition_nodal # name
        3 # dims
        1 # nprocs
        Examples/Tests/langmuir/inputs_3d_test_Langmuir_multi_psatd_Vay_deposition_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# Langmuir_multi_rz ###########################################################
#
add_warpx_test(
    Langmuir_multi_rz # name
    RZ # dims
    2  # nprocs
    Examples/Tests/langmuir/inputs_rz_test_Langmuir_multi_rz # inputs
    Examples/Tests/langmuir/analysis_rz.py # analysis
    diags/diag1000080 # output
)

# Langmuir_multi_rz_picmi #####################################################
#
add_warpx_test(
    Langmuir_multi_rz_picmi # name
    RZ # dims
    2  # nprocs
    Examples/Tests/langmuir/inputs_rz_test_Langmuir_multi_rz_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

# Langmuir_multi_rz_psatd #####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_rz_psatd # name
        RZ # dims
        2  # nprocs
        Examples/Tests/langmuir/inputs_rz_test_Langmuir_multi_rz_psatd # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_rz_psatd_current_correction ##################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_rz_psatd_correction # name
        RZ # dims
        1  # nprocs
        Examples/Tests/langmuir/inputs_rz_test_Langmuir_multi_rz_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()

# Langmuir_multi_rz_psatd_multiJ ##############################################
#
if(WarpX_FFT)
    add_warpx_test(
        Langmuir_multi_rz_psatd_multiJ # name
        RZ # dims
        2  # nprocs
        Examples/Tests/langmuir/inputs_rz_test_Langmuir_multi_rz_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()
