# Add tests (alphabetical order) ##############################################
#

# rigid_injection_btd #########################################################
#
add_warpx_test(
    rigid_injection_btd # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/rigid_injection/inputs_2d_test_rigid_injection_btd # inputs
    Examples/Tests/rigid_injection/analysis_rigid_injection_btd.py # analysis
    diags/diag1000001 # output
)

# rigid_injection_lab #########################################################
#
add_warpx_test(
    rigid_injection_lab # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/rigid_injection/inputs_2d_test_rigid_injection_lab # inputs
    Examples/Tests/rigid_injection/analysis_rigid_injection_lab.py # analysis
    diags/diag1000289 # output
)
