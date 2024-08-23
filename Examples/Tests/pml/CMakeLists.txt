# Add tests (alphabetical order) ##############################################
#

if(WarpX_FFT)
    add_warpx_test(
        test_3d_pml_psatd_dive_divb_cleaning # name
        3 # dims
        2 # nprocs
        OFF # eb
    OFF # preparation script
        Examples/Tests/pml/inputs_test_3d_pml_psatd_dive_divb_cleaning # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000100 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_rz_pml_psatd # name
        RZ # dims
        2  # nprocs
        OFF # eb
    OFF # preparation script
        Examples/Tests/pml/inputs_test_rz_pml_psatd # inputs
        Examples/Tests/pml/analysis_pml_psatd_rz.py # analysis
        diags/diag1000500 # output
    )
endif()

add_warpx_test(
    test_2d_pml_x_ckc # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pml/inputs_test_2d_pml_x_ckc # inputs
    Examples/Tests/pml/analysis_pml_ckc.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_2d_pml_x_galilean # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pml/inputs_test_2d_pml_x_galilean # inputs
    Examples/Tests/pml/analysis_pml_psatd.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_2d_pml_x_psatd # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pml/inputs_test_2d_pml_x_psatd # inputs
    Examples/Tests/pml/analysis_pml_psatd.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_2d_pml_x_psatd_restart # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pml/inputs_test_2d_pml_x_psatd_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_2d_pml_x_yee # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pml/inputs_test_2d_pml_x_yee # inputs
    Examples/Tests/pml/analysis_pml_yee.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_2d_pml_x_yee_restart # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pml/inputs_test_2d_pml_x_yee_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)
