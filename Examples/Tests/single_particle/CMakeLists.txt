# Add tests (alphabetical order) ##############################################
#

# bilinear_filter #############################################################
#
add_warpx_test(
    bilinear_filter # name
    2  # dims
    ON # mpi
    Examples/Tests/single_particle/inputs_2d_base # inputs
    Examples/Tests/single_particle/analysis_bilinear_filter.py # analysis
    diags/diag1000001 # output
)
