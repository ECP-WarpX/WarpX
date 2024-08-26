# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_larmor # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/larmor/inputs_test_2d_larmor # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # depend on this test
)
