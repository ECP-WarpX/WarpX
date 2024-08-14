# Add tests (alphabetical order) ##############################################
#

# langmuir_multi_1d ###########################################################
#
add_warpx_test(
    langmuir_multi_1d # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_1d_test_langmuir_multi_1d # inputs
    Examples/Tests/langmuir/analysis_1d.py # analysis
    diags/diag1000080 # output
)

# langmuir_multi_2d_mr ########################################################
#
add_warpx_test(
    langmuir_multi_2d_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_mr # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# langmuir_multi_2d_mr_anisotropic ############################################
#
add_warpx_test(
    langmuir_multi_2d_mr_anisotropic # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_mr_anisotropic # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# langmuir_multi_2d_mr_momentum_conserving ####################################
#
add_warpx_test(
    langmuir_multi_2d_mr_momentum_conserving # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_mr_momentum_conserving # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# langmuir_multi_2d_mr_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_mr_psatd # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_mr_psatd # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_nodal #####################################################
#
add_warpx_test(
    langmuir_multi_2d_nodal # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_nodal # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# langmuir_multi_2d_picmi #####################################################
#
add_warpx_test(
    langmuir_multi_2d_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

# langmuir_multi_2d_psatd #####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_current_correction ##################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_current_correction # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_current_correction_nodal ############################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_current_correction_nodal # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_current_correction_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_momentum_conserving #################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_momentum_conserving # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_momentum_conserving # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_multiJ ##############################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_multiJ # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_multiJ_nodal ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_multiJ_nodal # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_multiJ_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_nodal ###############################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_nodal # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_2d_psatd_vay_deposition ######################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_vay_deposition # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_vay_deposition # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# FIXME
# langmuir_multi_2d_psatd_vay_deposition_particle_shape_4 #####################
#
#if(WarpX_FFT)
#    add_warpx_test(
#        langmuir_multi_2d_psatd_vay_deposition_particle_shape_4 # name
#        2 # dims
#        1 # nprocs
#        OFF # eb
#        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_vay_deposition_particle_shape_4 # inputs
#        Examples/Tests/langmuir/analysis_2d.py # analysis
#        diags/diag1000080 # output
#    )
#endif()

# langmuir_multi_2d_psatd_vay_deposition_nodal ################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_2d_psatd_vay_deposition_nodal # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_2d_test_langmuir_multi_2d_psatd_vay_deposition_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_3d ###########################################################
#
add_warpx_test(
    langmuir_multi_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d # inputs
    Examples/Tests/langmuir/analysis_3d.py # analysis
    diags/diag1000040 # output
)

# langmuir_multi_3d_nodal #####################################################
#
add_warpx_test(
    langmuir_multi_3d_nodal # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_nodal # inputs
    Examples/Tests/langmuir/analysis_3d.py # analysis
    diags/diag1000040 # output
)

# langmuir_multi_3d_picmi #####################################################
#
add_warpx_test(
    langmuir_multi_3d_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

# langmuir_multi_3d_psatd #####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_current_correction ##################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_current_correction # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_current_correction_nodal ############################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_current_correction_nodal # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_current_correction_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_div_cleaning ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_div_cleaning # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_div_cleaning # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_momentum_conserving #################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_momentum_conserving # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_momentum_conserving # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_multiJ ##############################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_multiJ # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_multiJ_nodal ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_multiJ_nodal # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_multiJ_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_nodal ###############################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_nodal # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_vay_deposition ######################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_vay_deposition # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_vay_deposition # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_3d_psatd_vay_deposition_nodal ################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_3d_psatd_vay_deposition_nodal # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_3d_test_langmuir_multi_3d_psatd_vay_deposition_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

# langmuir_multi_rz ###########################################################
#
add_warpx_test(
    langmuir_multi_rz # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_rz_test_langmuir_multi_rz # inputs
    Examples/Tests/langmuir/analysis_rz.py # analysis
    diags/diag1000080 # output
)

# langmuir_multi_rz_picmi #####################################################
#
add_warpx_test(
    langmuir_multi_rz_picmi # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_rz_test_langmuir_multi_rz_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

# langmuir_multi_rz_psatd #####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_rz_psatd # name
        RZ # dims
        2  # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_rz_test_langmuir_multi_rz_psatd # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_rz_psatd_current_correction ##################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_rz_psatd_current_correction # name
        RZ # dims
        1  # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_rz_test_langmuir_multi_rz_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()

# langmuir_multi_rz_psatd_multiJ ##############################################
#
if(WarpX_FFT)
    add_warpx_test(
        langmuir_multi_rz_psatd_multiJ # name
        RZ # dims
        2  # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_rz_test_langmuir_multi_rz_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()
