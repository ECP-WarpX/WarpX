# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_relativistic_space_charge_initialization # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/relativistic_space_charge_initialization/inputs_test_3d_relativistic_space_charge_initialization # inputs
    Examples/Tests/relativistic_space_charge_initialization/analysis.py # analysis
    diags/diag1000001 # output
)
