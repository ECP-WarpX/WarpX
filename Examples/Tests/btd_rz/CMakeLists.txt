# Add tests (alphabetical order) ##############################################
#

# BTD_rz ######################################################################
#
add_warpx_test(
    BTD_rz # name
    RZ # dims
    ON # mpi
    Examples/Tests/btd_rz/inputs_rz_test_BTD_rz # inputs
    Examples/Tests/btd_rz/analysis_BTD_laser_antenna.py # analysis
    diags/diag1000289 # output
)
