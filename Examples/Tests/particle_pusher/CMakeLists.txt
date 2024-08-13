# Add tests (alphabetical order) ##############################################
#

# particle_pusher #############################################################
#
add_warpx_test(
    particle_pusher # name
    3 # dims
    1 # nprocs
    Examples/Tests/particle_pusher/inputs_3d_test_particle_pusher # inputs
    Examples/Tests/particle_pusher/analysis_pusher.py # analysis
    diags/diag1010000 # output
)
