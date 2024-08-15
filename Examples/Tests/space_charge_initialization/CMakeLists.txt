# Add tests (alphabetical order) ##############################################
#

# space_charge_initialization_2d ##############################################
#
add_warpx_test(
    space_charge_initialization_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/space_charge_initialization/inputs_3d_test_space_charge_initialization_2d # inputs
    Examples/Tests/space_charge_initialization/analysis.py # analysis
    diags/diag1000001 # output
)

# space_charge_initialization_3d ##############################################
#
add_warpx_test(
    space_charge_initialization_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/space_charge_initialization/inputs_3d_test_space_charge_initialization_3d # inputs
    Examples/Tests/space_charge_initialization/analysis.py # analysis
    diags/diag1000001 # output
)
