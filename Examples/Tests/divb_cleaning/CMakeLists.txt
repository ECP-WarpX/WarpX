# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_divb_cleaning # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/divb_cleaning/inputs_test_3d_divb_cleaning # inputs
    Examples/Tests/divb_cleaning/analysis.py # analysis
    diags/diag1000400 # output
    OFF # depend on this test
)
