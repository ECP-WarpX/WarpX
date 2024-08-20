# Add tests (alphabetical order) ##############################################
#

# plasma_acceleration_1d_picmi ################################################
#
add_warpx_test(
    plasma_acceleration_1d_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_1d_test_plasma_acceleration_1d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1001000 # output
)

# plasma_acceleration_2d_mr ###################################################
#
add_warpx_test(
    plasma_acceleration_2d_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_2d_test_plasma_acceleration_2d_mr # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000400 # output
)

# plasma_acceleration_2d_mr_momentum_conserving ###############################
#
add_warpx_test(
    plasma_acceleration_2d_mr_momentum_conserving # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_2d_test_plasma_acceleration_2d_mr_momentum_conserving # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000400 # output
)

# plasma_acceleration_boosted_2d ##############################################
#
add_warpx_test(
    plasma_acceleration_boosted_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_2d_test_plasma_acceleration_boosted_2d # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000020 # output
)

# plasma_acceleration_boosted_3d ##############################################
#
add_warpx_test(
    plasma_acceleration_boosted_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_3d_test_plasma_acceleration_boosted_3d # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000005 # output
)

# plasma_acceleration_boosted_3d_hybrid #######################################
#
add_warpx_test(
    plasma_acceleration_boosted_3d_hybrid # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_3d_test_plasma_acceleration_boosted_3d_hybrid # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000025 # output
)

# plasma_acceleration_3d_picmi ################################################
#
add_warpx_test(
    plasma_acceleration_3d_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_3d_test_plasma_acceleration_3d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# plasma_acceleration_3d_mr_picmi #############################################
#
add_warpx_test(
    plasma_acceleration_3d_mr_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_acceleration/inputs_3d_test_plasma_acceleration_3d_mr_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000002 # output
)
