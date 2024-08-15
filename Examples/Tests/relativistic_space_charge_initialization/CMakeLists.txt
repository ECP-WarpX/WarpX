# Add tests (alphabetical order) ##############################################
#

# relativistic_space_charge_initialization ####################################
#
add_warpx_test(
    relativistic_space_charge_initialization # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/relativistic_space_charge_initialization/inputs_3d_test_relativistic_space_charge_initialization # inputs
    Examples/Tests/relativistic_space_charge_initialization/analysis.py # analysis
    diags/diag1000001 # output
)
