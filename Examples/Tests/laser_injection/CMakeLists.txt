# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_laser_injection # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/laser_injection/inputs_test_1d_laser_injection # inputs
    Examples/Tests/laser_injection/analysis_1d.py # analysis
    diags/diag1000240 # output
)

add_warpx_test(
    test_2d_laser_injection # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/laser_injection/inputs_test_2d_laser_injection # inputs
    Examples/Tests/laser_injection/analysis_2d.py # analysis
    diags/diag1000240 # output
)

add_warpx_test(
    test_3d_laser_injection # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/laser_injection/inputs_test_3d_laser_injection # inputs
    Examples/Tests/laser_injection/analysis_3d.py # analysis
    diags/diag1000020 # output
)
