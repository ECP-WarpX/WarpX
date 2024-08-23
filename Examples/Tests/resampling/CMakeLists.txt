# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_leveling_thinning # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/resampling/inputs_test_2d_leveling_thinning # inputs
    Examples/Tests/resampling/analysis.py # analysis
    diags/diag1000008 # output
)

add_warpx_test(
    test_1d_resample_velocity_coincidence_thinning # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/resampling/inputs_test_1d_resample_velocity_coincidence_thinning # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000004 # output
)

add_warpx_test(
    test_1d_resample_velocity_coincidence_thinning_cartesian # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/resampling/inputs_test_1d_resample_velocity_coincidence_thinning_cartesian # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000004 # output
)
