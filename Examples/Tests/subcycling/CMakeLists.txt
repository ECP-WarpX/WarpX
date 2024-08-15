# Add tests (alphabetical order) ##############################################
#

# subcycling_mr ###############################################################
#
add_warpx_test(
    subcycling_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/subcycling/inputs_2d_test_subcycling_mr # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000250 # output
)
