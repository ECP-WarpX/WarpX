# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_projection_divb_cleaner_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/projection_divb_cleaner/inputs_test_3d_projection_divb_cleaner_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_projection_divb_cleaner_callback_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/projection_divb_cleaner/inputs_test_3d_projection_divb_cleaner_callback_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)

add_warpx_test(
    test_rz_projection_divb_cleaner # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/projection_divb_cleaner/inputs_test_rz_projection_divb_cleaner # inputs
    Examples/Tests/projection_divb_cleaner/analysis.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)
