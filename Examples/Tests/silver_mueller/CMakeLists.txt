# Add tests (alphabetical order) ##############################################
#

# silver_mueller_1d ###########################################################
#
add_warpx_test(
    silver_mueller_1d # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/silver_mueller/inputs_1d_test_silver_mueller_1d # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)

# silver_mueller_2d_x #########################################################
#
add_warpx_test(
    silver_mueller_2d_x # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/silver_mueller/inputs_2d_test_silver_mueller_2d_x # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)

# silver_mueller_2d_z #########################################################
#
add_warpx_test(
    silver_mueller_2d_z # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/silver_mueller/inputs_2d_test_silver_mueller_2d_z # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)

# silver_mueller_rz_z #########################################################
#
add_warpx_test(
    silver_mueller_rz_z # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/silver_mueller/inputs_rz_test_silver_mueller_rz_z # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)
