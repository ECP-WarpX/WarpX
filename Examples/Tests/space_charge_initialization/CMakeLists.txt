# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_space_charge_initialization # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/space_charge_initialization/inputs_test_2d_space_charge_initialization # inputs
    Examples/Tests/space_charge_initialization/analysis.py # analysis
    diags/diag1000001 # output
    OFF # depend on this test
)

add_warpx_test(
    test_3d_space_charge_initialization # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/space_charge_initialization/inputs_test_3d_space_charge_initialization # inputs
    Examples/Tests/space_charge_initialization/analysis.py # analysis
    diags/diag1000001 # output
    OFF # depend on this test
)
