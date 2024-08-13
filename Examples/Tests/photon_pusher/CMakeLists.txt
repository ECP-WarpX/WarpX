# Add tests (alphabetical order) ##############################################
#

# photon_pusher ###############################################################
#
add_warpx_test(
    photon_pusher # name
    3 # dims
    2 # nprocs
    Examples/Tests/photon_pusher/inputs_3d_test_photon_pusher # inputs
    Examples/Tests/photon_pusher/analysis_photon_pusher.py # analysis
    diags/diag1000050 # output
)
