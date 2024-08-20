# Add tests (alphabetical order) ##############################################
#

# plasma_mirror ###############################################################
#
add_warpx_test(
    plasma_mirror # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/plasma_mirror/inputs_2d_test_plasma_mirror # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000020 # output
)
