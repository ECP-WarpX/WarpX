# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_bilinear_filter # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/single_particle/inputs_test_2d_bilinear_filter # inputs
    Examples/Tests/single_particle/analysis.py # analysis
    diags/diag1000001 # output
    OFF # depend on this test
)
