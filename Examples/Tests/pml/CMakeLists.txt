# Add tests (alphabetical order) ##############################################
#

# pml_psatd_dive_divb_cleaning ################################################
#
if(WarpX_FFT)
    add_warpx_test(
        pml_psatd_dive_divb_cleaning # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/pml/inputs_3d_test_pml_psatd_dive_divb_cleaning # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000100 # output
    )
endif()

# pml_psatd_rz ################################################################
#
if(WarpX_FFT)
    add_warpx_test(
        pml_psatd_rz # name
        RZ # dims
        2  # nprocs
        OFF # eb
        Examples/Tests/pml/inputs_rz_test_pml_psatd_rz # inputs
        Examples/Tests/pml/analysis_pml_psatd_rz.py # analysis
        diags/diag1000500 # output
    )
endif()

# pml_x_ckc ###################################################################
#
add_warpx_test(
    pml_x_ckc # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pml/inputs_2d_test_pml_x_ckc # inputs
    Examples/Tests/pml/analysis_pml_ckc.py # analysis
    diags/diag1000300 # output
)

# pml_x_galilean ##############################################################
#
add_warpx_test(
    pml_x_galilean # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pml/inputs_2d_test_pml_x_galilean # inputs
    Examples/Tests/pml/analysis_pml_psatd.py # analysis
    diags/diag1000300 # output
)

# pml_x_psatd #################################################################
#
add_warpx_test(
    pml_x_psatd # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pml/inputs_2d_test_pml_x_psatd # inputs
    Examples/Tests/pml/analysis_pml_psatd.py # analysis
    diags/diag1000300 # output
)

# pml_x_psatd_restart #########################################################
#
add_warpx_test(
    pml_x_psatd_restart # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pml/inputs_2d_test_pml_x_psatd_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)

# pml_x_yee ###################################################################
#
add_warpx_test(
    pml_x_yee # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pml/inputs_2d_test_pml_x_yee # inputs
    Examples/Tests/pml/analysis_pml_yee.py # analysis
    diags/diag1000300 # output
)

# pml_x_yee_restart ###########################################################
#
add_warpx_test(
    pml_x_yee_restart # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pml/inputs_2d_test_pml_x_yee_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)
