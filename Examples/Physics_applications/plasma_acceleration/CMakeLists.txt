# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_plasma_acceleration_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_1d_plasma_acceleration_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1001000 # output
)

add_warpx_test(
    test_2d_plasma_acceleration_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_2d_plasma_acceleration_mr # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000400 # output
)

add_warpx_test(
    test_2d_plasma_acceleration_mr_momentum_conserving # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_2d_plasma_acceleration_mr_momentum_conserving # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000400 # output
)

add_warpx_test(
    test_2d_plasma_acceleration_boosted # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_2d_plasma_acceleration_boosted # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000020 # output
)

add_warpx_test(
    test_3d_plasma_acceleration_boosted # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_3d_plasma_acceleration_boosted # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000005 # output
)

add_warpx_test(
    test_3d_plasma_acceleration_boosted_hybrid # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_3d_plasma_acceleration_boosted_hybrid # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000025 # output
)

add_warpx_test(
    test_3d_plasma_acceleration_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_3d_plasma_acceleration_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

add_warpx_test(
    test_3d_plasma_acceleration_mr_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/plasma_acceleration/inputs_test_3d_plasma_acceleration_mr_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000002 # output
)
