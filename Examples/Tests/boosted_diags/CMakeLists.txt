# Add tests (alphabetical order) ##############################################
#

# laser_acceleration_btd ######################################################
#
add_warpx_test(
    laser_acceleration_btd # name
    3  # dims
    ON # mpi
    Examples/Tests/boosted_diags/inputs_3d_test_laser_acceleration_btd # inputs
    Examples/Tests/boosted_diags/analysis.py # analysis
    diags/diag1000003 # output
)
