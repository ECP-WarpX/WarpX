# Add tests (alphabetical order) ##############################################
#

# uniform_plasma_2d ###########################################################
#
add_warpx_test(
    uniform_plasma_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/uniform_plasma/inputs_2d_test_uniform_plasma_2d # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# uniform_plasma_3d ###########################################################
#
add_warpx_test(
    uniform_plasma_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/uniform_plasma/inputs_3d_test_uniform_plasma_3d # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# uniform_plasma_3d_restart ###################################################
#
add_warpx_test(
    uniform_plasma_3d_restart # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/uniform_plasma/inputs_3d_test_uniform_plasma_3d_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)
