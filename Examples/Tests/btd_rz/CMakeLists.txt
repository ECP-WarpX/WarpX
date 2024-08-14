# Add tests (alphabetical order) ##############################################
#

# btd_rz ######################################################################
#
add_warpx_test(
    btd_rz # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/btd_rz/inputs_rz_test_btd_rz # inputs
    Examples/Tests/btd_rz/analysis_btd_laser_antenna.py # analysis
    diags/diag1000289 # output
)
