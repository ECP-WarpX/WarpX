# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_laser_acceleration_btd # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/boosted_diags/inputs_test_3d_laser_acceleration_btd # inputs
    Examples/Tests/boosted_diags/analysis.py # analysis
    diags/diag1000003 # output
)
