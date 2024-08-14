# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    Python_gaussian_beam # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/gaussian_beam/PICMI_inputs_gaussian_beam.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)
