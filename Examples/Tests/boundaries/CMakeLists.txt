# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_particle_boundaries # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/boundaries/inputs_test_3d_particle_boundaries # inputs
    Examples/Tests/boundaries/analysis.py # analysis
    diags/diag1000008 # output
    OFF # dependency
)
