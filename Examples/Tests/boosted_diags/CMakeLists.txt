# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_laser_acceleration_btd  # name
    3  # dims
    2  # nprocs
    OFF  # eb
    inputs_test_3d_laser_acceleration_btd  # inputs
    analysis.py  # analysis
    diags/diag1000003  # output
    OFF  # dependency
)
