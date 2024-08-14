# Add tests (alphabetical order) ##############################################
#

# laser_acceleration_1d_picmi #################################################
#
add_warpx_test(
    laser_acceleration_1d_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_1d_test_laser_acceleration_1d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
)
