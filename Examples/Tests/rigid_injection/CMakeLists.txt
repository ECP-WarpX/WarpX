# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_rigid_injection_btd # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/rigid_injection/inputs_test_2d_rigid_injection_btd # inputs
    Examples/Tests/rigid_injection/analysis_rigid_injection_btd.py # analysis
    diags/diag1000001 # output
    OFF # depend on this test
)

add_warpx_test(
    test_2d_rigid_injection_lab # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/rigid_injection/inputs_test_2d_rigid_injection_lab # inputs
    Examples/Tests/rigid_injection/analysis_rigid_injection_lab.py # analysis
    diags/diag1000289 # output
    OFF # depend on this test
)
