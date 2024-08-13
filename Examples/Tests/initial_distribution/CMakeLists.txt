# Add tests (alphabetical order) ##############################################
#

# initial_distribution ########################################################
#
add_warpx_test(
    initial_distribution # name
    3 # dims
    1 # nprocs
    Examples/Tests/initial_distribution/inputs_3d_test_initial_distribution # inputs
    Examples/Tests/initial_distribution/analysis_distribution.py # analysis
    diags/diag1000001 # output
)
