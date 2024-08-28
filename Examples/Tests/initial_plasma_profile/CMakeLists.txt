# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_parabolic_channel_initialization # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/initial_plasma_profile/inputs_test_2d_parabolic_channel_initialization # inputs
    Examples/Tests/initial_plasma_profile/analysis.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)
