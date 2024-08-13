# Add tests (alphabetical order) ##############################################
#

# Larmor ######################################################################
#
add_warpx_test(
    Larmor # name
    2 # dims
    2 # nprocs
    Examples/Tests/larmor/inputs_2d_test_Larmor # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)
