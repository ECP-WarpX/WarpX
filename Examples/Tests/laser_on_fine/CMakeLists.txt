# Add tests (alphabetical order) ##############################################
#

# laser_on_fine ###############################################################
#
add_warpx_test(
    laser_on_fine # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/laser_on_fine/inputs_2d_test_laser_on_fine # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000050 # output
)
