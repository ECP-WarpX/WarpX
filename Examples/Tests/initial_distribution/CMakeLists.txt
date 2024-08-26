# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_initial_distribution # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/initial_distribution/inputs_test_3d_initial_distribution # inputs
    Examples/Tests/initial_distribution/analysis.py # analysis
    diags/diag1000001 # output
    OFF # depend on this test
)
