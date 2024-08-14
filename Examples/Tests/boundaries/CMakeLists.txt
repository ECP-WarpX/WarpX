# Add tests (alphabetical order) ##############################################
#

# particle_boundaries_3d ######################################################
#
add_warpx_test(
    particle_boundaries_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/boundaries/inputs_3d_test_particle_boundaries_3d # inputs
    Examples/Tests/boundaries/analysis.py # analysis
    diags/diag1000008 # output
)
