# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_silver_mueller # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/silver_mueller/inputs_test_1d_silver_mueller # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)

add_warpx_test(
    test_2d_silver_mueller_x # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/silver_mueller/inputs_test_2d_silver_mueller_x # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)

add_warpx_test(
    test_2d_silver_mueller_z # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/silver_mueller/inputs_test_2d_silver_mueller_z # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)

add_warpx_test(
    test_rz_silver_mueller_z # name
    RZ # dims
    2  # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/silver_mueller/inputs_test_rz_silver_mueller_z # inputs
    Examples/Tests/silver_mueller/analysis.py # analysis
    diags/diag1000500 # output
)
