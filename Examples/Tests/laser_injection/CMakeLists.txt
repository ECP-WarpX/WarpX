# Add tests (alphabetical order) ##############################################
#

# laser_injection_1d ##########################################################
#
add_warpx_test(
    laser_injection_1d # name
    1  # dims
    ON # mpi
    Examples/Tests/laser_injection/inputs_1d_test_laser_injection_1d # inputs
    Examples/Tests/laser_injection/analysis_1d.py # analysis
    diags/diag1000240 # output
)

# laser_injection_2d ##########################################################
#
add_warpx_test(
    laser_injection_2d # name
    2  # dims
    ON # mpi
    Examples/Tests/laser_injection/inputs_2d_test_laser_injection_2d # inputs
    Examples/Tests/laser_injection/analysis_2d.py # analysis
    diags/diag1000240 # output
)

# laser_injection_3d ##########################################################
#
add_warpx_test(
    laser_injection_3d # name
    3  # dims
    ON # mpi
    Examples/Tests/laser_injection/inputs_3d_test_laser_injection_3d # inputs
    Examples/Tests/laser_injection/analysis_3d.py # analysis
    diags/diag1001000 # output
)
