# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_laser_on_fine # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/laser_on_fine/inputs_test_2d_laser_on_fine # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000050 # output
)
