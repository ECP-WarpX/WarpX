# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_dive_cleaning # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/dive_cleaning/inputs_test_2d_dive_cleaning # inputs
    Examples/Tests/dive_cleaning/analysis.py # analysis
    diags/diag1000128 # output
)

add_warpx_test(
    test_3d_dive_cleaning # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/dive_cleaning/inputs_test_3d_dive_cleaning # inputs
    Examples/Tests/dive_cleaning/analysis.py # analysis
    diags/diag1000128 # output
)
