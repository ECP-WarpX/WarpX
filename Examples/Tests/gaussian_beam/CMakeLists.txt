# Add tests (alphabetical order) ##############################################
#

# gaussian_beam_picmi #########################################################
#
add_warpx_test(
    gaussian_beam_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/gaussian_beam/inputs_2d_test_gaussian_beam_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)
