# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_uniform_plasma # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/uniform_plasma/inputs_test_2d_uniform_plasma # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # depend on this test
)

add_warpx_test(
    test_3d_uniform_plasma # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/uniform_plasma/inputs_test_3d_uniform_plasma # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # depend on this test
)

add_warpx_test(
    test_3d_uniform_plasma_restart # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/uniform_plasma/inputs_test_3d_uniform_plasma_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # depend on this test
)
