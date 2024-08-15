# Add tests (alphabetical order) ##############################################
#

# leveling_thinning ###########################################################
#
add_warpx_test(
    leveling_thinning # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/resampling/inputs_2d_test_leveling_thinning # inputs
    Examples/Tests/resampling/analysis.py # analysis
    diags/diag1000008 # output
)

# resample_velocity_coincidence_thinning ######################################
#
add_warpx_test(
    resample_velocity_coincidence_thinning # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/resampling/inputs_2d_test_resample_velocity_coincidence_thinning # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000004 # output
)

# resample_velocity_coincidence_thinning_cartesian ############################
#
add_warpx_test(
    resample_velocity_coincidence_thinning_cartesian # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/resampling/inputs_2d_test_resample_velocity_coincidence_thinning_cartesian # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000004 # output
)
