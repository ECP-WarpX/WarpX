# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_particle_pusher # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/particle_pusher/inputs_test_3d_particle_pusher # inputs
    Examples/Tests/particle_pusher/analysis.py # analysis
    diags/diag1010000 # output
    OFF # dependency
)
