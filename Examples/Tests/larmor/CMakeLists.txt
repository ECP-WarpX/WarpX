# Add tests (alphabetical order) ##############################################
#

# larmor ######################################################################
#
add_warpx_test(
    larmor # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/larmor/inputs_2d_test_larmor # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)
